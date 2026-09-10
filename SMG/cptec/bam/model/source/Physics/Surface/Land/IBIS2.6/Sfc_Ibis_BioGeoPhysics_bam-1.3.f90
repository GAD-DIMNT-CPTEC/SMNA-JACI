MODULE Sfc_Ibis_LsxMain

  ! lsxmain_______setsoi
  !            |
  !            |__fwetcal
  !            |
  !            |__solset
  !            |
  !            |__solsur
  !            |
  !            |__solalb___twostr__twoset
  !            |        |
  !            |        |__twostr__twoset
  !            |
  !            |__solarf
  !            |
  !            |__irrad
  !            |
  !            |__cascade___mix
  !            |         | 
  !            |         |__mix
  !            |         |
  !            |         |__steph2o__mix
  !            |         |        |
  !            |         |        |__mix
  !            |         |
  !            |         |__steph2o__mix
  !            |         |        |
  !            |         |        |__mix
  !            |         |
  !            |         |__mix
  !            |         |
  !            |         |__mix
  !            |         |
  !            |         |__steph2o__mix
  !            |         |        |
  !            |         |        |__mix  
  !            |         |
  !            |         |
  !            |         |__mix
  !            |         |
  !            |         |__mix
  !            | 
  !            |__fwetcal
  !            |
  !            |__canopy__canini
  !            |       |
  !            |       |__drystress
  !            |       |
  !            |       |__turcof__fstrat
  !            |       |       |
  !            |       |       |__fstrat
  !            |       |
  !            |       |__stomata
  !            |       |
  !            |       |__turvap___impexp
  !            |       |       | 
  !            |       |       |__impexp
  !            |       |       |
  !            |       |       |__impexp
  !            |       |       |
  !            |       |       |__impexp
  !            |       |       |
  !            |       |       |__impexp2
  !            |       |       |
  !            |       |       |__linsolve
  !            |       |
  !            |       |__tscreen
  !            |       |
  !            |       |__tscreen
  !            |
  !            |__cascad2__steph2o2
  !            |        |
  !            |        |__steph2o2
  !            |        |
  !            |        |__steph2o2
  !            |
  !            |__noveg
  !            |
  !            |__snow__snowheat___tridia
  !            |     |
  !            |     |__vadapt
  !            |     |
  !            |     |__MixSnow
  !            |
  !            |__soilctl__soilh2o__tridia
  !                     |
  !                     |__soilheat__tridia
  !                     |
  !                     |__wadjust
  !                     |
  !                     |__wadjust
  USE PhysicalFunctions, Only:fpvs2es5,ftdp
  IMPLICIT NONE
SAVE

  PRIVATE
  INTEGER, PARAMETER :: r4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers 
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers

  !
  ! ---------------------------------------------------------------------
  ! statement functions and associated parameters
  ! ---------------------------------------------------------------------
  !
  ! polynomials for svp(t), d(svp)/dt over water and ice are from
  ! lowe(1977),jam,16,101-103.
  !
  !
  REAL(KIND=r8), PARAMETER :: hvap = 2.5104e+6_r8  ! latent heat of vaporization of water (J kg-1)
  REAL(KIND=r8), PARAMETER :: hfus = 0.3336e+6_r8! latent heat of fusion of water (J kg-1)
  REAL(KIND=r8), PARAMETER :: hsub = hvap + hfus ! latent heat of sublimation of ice (J kg-1)
  REAL(KIND=r8), PARAMETER :: cice = 2.106e+3_r8  ! specific heat of ice (J deg-1 kg-1)

  REAL(KIND=r8), PARAMETER :: ch2o = 4.218e+3_r8  ! specific heat of liquid water (J deg-1 kg-1)
  REAL(KIND=r8), PARAMETER :: cvap  = 1.81e+3_r8   ! specific heat of water vapor at constant pressure (J deg-1 kg-1)


  REAL(KIND=r8), PARAMETER :: asat0 =  6.1078000_r8
  REAL(KIND=r8), PARAMETER :: asat1 =  4.4365185e-1_r8
  REAL(KIND=r8), PARAMETER :: asat2 =  1.4289458e-2_r8
  REAL(KIND=r8), PARAMETER :: asat3 =  2.6506485e-4_r8
  REAL(KIND=r8), PARAMETER :: asat4 =  3.0312404e-6_r8
  REAL(KIND=r8), PARAMETER :: asat5 =  2.0340809e-8_r8
  REAL(KIND=r8), PARAMETER :: asat6 =  6.1368209e-11_r8 
  !
  !
  REAL(KIND=r8), PARAMETER :: bsat0 =  6.1091780_r8
  REAL(KIND=r8), PARAMETER :: bsat1 =  5.0346990e-1_r8
  REAL(KIND=r8), PARAMETER :: bsat2 =  1.8860134e-2_r8
  REAL(KIND=r8), PARAMETER :: bsat3 =  4.1762237e-4_r8
  REAL(KIND=r8), PARAMETER :: bsat4 =  5.8247203e-6_r8
  REAL(KIND=r8), PARAMETER :: bsat5 =  4.8388032e-8_r8
  REAL(KIND=r8), PARAMETER :: bsat6 =  1.8388269e-10_r8 
  !
  !
  REAL(KIND=r8), PARAMETER :: csat0 =  4.4381000e-1_r8
  REAL(KIND=r8), PARAMETER :: csat1 =  2.8570026e-2_r8
  REAL(KIND=r8), PARAMETER :: csat2 =  7.9380540e-4_r8
  REAL(KIND=r8), PARAMETER :: csat3 =  1.2152151e-5_r8
  REAL(KIND=r8), PARAMETER :: csat4 =  1.0365614e-7_r8
  REAL(KIND=r8), PARAMETER :: csat5 =  3.5324218e-10_r8
  REAL(KIND=r8), PARAMETER :: csat6 = -7.0902448e-13_r8
  !
  !
  REAL(KIND=r8), PARAMETER :: dsat0 =  5.0303052e-1_r8
  REAL(KIND=r8), PARAMETER :: dsat1 =  3.7732550e-2_r8
  REAL(KIND=r8), PARAMETER :: dsat2 =  1.2679954e-3_r8
  REAL(KIND=r8), PARAMETER :: dsat3 =  2.4775631e-5_r8
  REAL(KIND=r8), PARAMETER :: dsat4 =  3.0056931e-7_r8
  REAL(KIND=r8), PARAMETER :: dsat5 =  2.1585425e-9_r8
  REAL(KIND=r8), PARAMETER :: dsat6 =  7.1310977e-12_r8
  REAL(KIND=r8), PARAMETER, PUBLIC :: MAPL_AIRMW  = 28.97_r8               ! kg/Kmole
  REAL(KIND=r8), PARAMETER, PUBLIC :: MAPL_H2OMW  = 18.01_r8                  ! kg/Kmole

  REAL(KIND=r8), PARAMETER, PUBLIC :: MAPL_VIREPS = MAPL_AIRMW/MAPL_H2OMW-1.0_r8   ! --
  ! use uniform value 1.0 for average diffuse optical depth
  ! (although an array for solar, all values are set to 1 in twoset).
  ! The typical values of emissivity are 0.80-0.95 for bare soil,
  ! 0.95-0.97 for vegetated areas and 0.99 for snow 
  ! (Wilber et al. 1999; Jin 2004; Jin and Liang 2004, manuscript submitted to J. Climate)


  REAL(KIND=r8), PARAMETER :: avmuir = 1.0_r8      ! average diffuse optical depth

  PUBLIC :: lsxmain
  PUBLIC :: solset
  PUBLIC :: solsur
  PUBLIC :: solalb
  PUBLIC :: solarf
  PUBLIC :: irrad
  PUBLIC :: fwetcal
  PUBLIC :: co2
CONTAINS
  ! ---------------------------------------------------------------
  SUBROUTINE lsxmain(ginvap       , &! INTENT(OUT  )!local
       gsuvap       , &! INTENT(OUT  )!local
       gtrans       , &! INTENT(OUT  )!local
       gtransu      , &! INTENT(OUT  )!local
       gtransl      , &! INTENT(OUT  )!local
       grunof       , &! INTENT(OUT  )!local
       gdrain       , &! INTENT(OUT  )!local
       gadjust      , &! INTENT(INOUT) !local
       a10scalparamu, &! INTENT(INOUT) !global
       a10daylightu , &! INTENT(INOUT) !global
       a10scalparaml, &! INTENT(INOUT) !global
       a10daylightl , &! INTENT(INOUT) !global
       vmax_pft     , &! INTENT(IN) !global
       tau15        , &! INTENT(IN) !global
       kc15         , &! INTENT(IN) !global
       ko15         , &! INTENT(IN) !global
       cimax        , &! INTENT(IN) !global
       gammaub      , &! INTENT(IN) !global
       alpha3       , &! INTENT(IN) !global
       theta3       , &! INTENT(IN) !global
       beta3        , &! INTENT(IN) !global
       coefmub      , &! INTENT(IN) !global
       coefbub      , &! INTENT(IN) !global
       gsubmin      , &! INTENT(IN) !global
       gammauc      , &! INTENT(IN) !global
       coefmuc      , &! INTENT(IN) !global
       coefbuc      , &! INTENT(IN) !global
       gsucmin      , &! INTENT(IN) !global
       gammals      ,&! INTENT(IN) !global
       coefmls      , &! INTENT(IN) !global
       coefbls      , &! INTENT(IN) !global
       gslsmin      , &! INTENT(IN) !global
       gammal3      , &! INTENT(IN) !global
       coefml3      , &! INTENT(IN) !global
       coefbl3      , &! INTENT(IN) !global
       gsl3min      , &! INTENT(IN    ) !global
       gammal4      , &! INTENT(IN) !global
       alpha4       , &! INTENT(IN) !global
       theta4       ,&! INTENT(IN) !global
       beta4        , &! INTENT(IN) !global
       coefml4      , &! INTENT(IN) !global
       coefbl4      , &! INTENT(IN) !global
       gsl4min      , &! INTENT(IN) !global
       wliqu        , &! INTENT(INOUT ) !global
       wliqumax     , &! INTENT(IN) !global
       wsnou        , &! INTENT(INOUT ) !global
       wsnoumax     , &! INTENT(IN) !global
       tu           , &! INTENT(INOUT ) !global
       wliqs        , &! INTENT(INOUT ) !global
       wliqsmax     , &! INTENT(IN) !global
       wsnos        , &! INTENT(INOUT ) !global
       wsnosmax     , &! INTENT(IN) !global
       ts           , &! INTENT(INOUT ) !global
       wliql        , &! INTENT(INOUT ) !global
       wliqlmax     , &! INTENT(IN) !global
       wsnol        , &! INTENT(INOUT ) !global
       wsnolmax     , &! INTENT(IN) !global
       tl           , &! INTENT(INOUT ) !global
       topparu      , &! INTENT(INOUT ) !local
       topparl      , &! INTENT(INOUT ) !local
       fl           , &! INTENT(IN) !global
       fu              , &! INTENT(IN) !global
       lai          , &! INTENT(IN) !global
       sai              , &! INTENT(IN) !global
       rhoveg       , &! INTENT(IN) !global
       tauveg       , &! INTENT(IN) !global
       orieh        , &! INTENT(IN) !global
       oriev        , &! INTENT(IN) !global
       wliqmin      , &! INTENT(INOUT ) !local
       wsnomin      , &! INTENT(INOUT ) !local
       t12          , &! INTENT(INOUT ) !global
       tdripu       , &! INTENT(IN) !global
       tblowu       , &! INTENT(IN) !global
       tdrips       , &! INTENT(IN) !global
       tblows       , &! INTENT(IN) !global
       t34          , &! INTENT(INOUT ) !global
       tdripl       , &! INTENT(IN) !global
       tblowl       , &! INTENT(IN) !global
       ztop         , &! INTENT(IN) !global
       alaiml       , &! INTENT(IN) !global
       zbot              , &! INTENT(IN) !global
       alaimu       , &! INTENT(IN) !global
       froot        , &! INTENT(IN) !global
       q34          , &! INTENT(INOUT ) !global
       q12          , &! INTENT(INOUT ) !global
       su           , &! INTENT(INOUT)  !local
       cleaf        , &! INTENT(IN) !global
       dleaf        , &! INTENT(IN) !global
       ss           , &! INTENT(INOUT)  !local
       cstem        , &! INTENT(IN) !global
       dstem        , &! INTENT(IN) !global
       sl           , &! INTENT(INOUT) !local
       cgrass       , &! INTENT(IN) !global
       ciub         , &! INTENT(INOUT) !global
       ciuc         , &! INTENT(INOUT) !global
       exist        , &! INTENT(IN) !global
       csub              , &! INTENT(INOUT) !global
       gsub         , &! INTENT(INOUT) !global
       csuc         , &! INTENT(INOUT) !global
       gsuc              , &! INTENT(INOUT) !global
       agcub        , &! INTENT(OUT  ) !local
       agcuc        , &! INTENT(OUT  ) !local
       ancub        , &! INTENT(OUT  ) !local
       ancuc        , &! INTENT(OUT  ) !local
       totcondub    , &! INTENT(INOUT) !local
       totconduc    , &! INTENT(INOUT) !local
       cils         , &! INTENT(INOUT) !global
       cil3         , &! INTENT(INOUT) !global
       cil4         , &! INTENT(INOUT) !global
       csls              , &! INTENT(INOUT) !global
       gsls         , &! INTENT(INOUT) !global
       csl3         , &! INTENT(INOUT) !global
       gsl3         , &! INTENT(INOUT) !global
       csl4         , &! INTENT(INOUT) !global
       gsl4         , &! INTENT(INOUT) !global
       agcls        , &! INTENT(OUT  ) !local
       agcl4        , &! INTENT(OUT  ) !local
       agcl3        , &! INTENT(OUT  ) !local
       ancls        , &! INTENT(OUT  ) !local
       ancl4        , &! INTENT(OUT  ) !local
       ancl3        , &! INTENT(OUT  ) !local
       totcondls    , &! INTENT(INOUT) !local
       totcondl3    ,&! INTENT(INOUT) !local
       totcondl4    , &! INTENT(INOUT) !local
       chu          , &! INTENT(IN) !global
       chs          , &! INTENT(IN) !global
       chl          , &! INTENT(IN) !global
       frac         , &! INTENT(IN) !global
       tlsub              , &! INTENT(INOUT)!global
       z0sno        , &! INTENT(IN) !global
       rhos              , &! INTENT(IN) !global
       consno       , &! INTENT(IN) !global
       hsnotop      , &! INTENT(IN) !global
       hsnomin      , &! INTENT(IN) !global
       fimin        , &! INTENT(IN) !global
       fimax        , &! INTENT(IN) !global
       fi              , &! INTENT(INOUT)!global
       tsno              , &! INTENT(INOUT)!global
       hsno         , &! INTENT(INOUT)!global
       sand              , &! INTENT(IN) !global
       clay         , &! INTENT(IN) !global
       poros        , &! INTENT(IN) !global
       wsoi         , &! INTENT(INOUT)!global
       wisoi        , &! INTENT(INOUT)!global
       consoi       , &! INTENT(INOUT) !local
       zwpmax       , &! INTENT(IN) !global
       wpud         , &! INTENT(INOUT)!global
       wipud        , &! INTENT(INOUT)!global
       wpudmax      , &! INTENT(IN) !global
       qglif        , &! INTENT(INOUT) !local
       tsoi         , &! INTENT(INOUT)!global
       hvasug       , &! INTENT(INOUT) !local
       hvasui       , &! INTENT(INOUT) !local
       albsav       , &! INTENT(IN) !global
       albsan       , &! INTENT(IN) !global
       tg           , &! INTENT(INOUT) !global
       ti           , &! INTENT(INOUT) !global
       z0soi        , &! INTENT(IN) !global
       swilt        , &! INTENT(IN) !global
       sfield       , &! INTENT(IN) !global
       stressl      , &! INTENT(INOUT) !local
       stressu      , &! INTENT(INOUT) !local
       stresstl     , &! INTENT(INOUT) !local
       stresstu     , &! INTENT(INOUT) !local
       csoi         , &! INTENT(IN) !global
       rhosoi       , &! INTENT(IN) !global
       hsoi         , &! INTENT(IN) !global
       suction      , &! INTENT(IN) !global
       bex          , &! INTENT(IN) !global
       upsoiu       , &! INTENT(INOUT) !local
       upsoil       , &! INTENT(INOUT) !local
       heatg        , &! INTENT(INOUT) !local
       heati        , &! INTENT(INOUT) !local
       hydraul      , &! INTENT(IN) !global
       porosflo     , &! INTENT(INOUT)!global
       ibex         , &! INTENT(IN) !global
       bperm        , &! INTENT(IN) !global
       hflo         , &! INTENT(INOUT)!global
       ta              , &! INTENT(IN) !global
       asurd        , &! INTENT(INOUT)! local
       asuri        , &! INTENT(INOUT)! local
       coszen       , &! INTENT(IN) !global
       solad        , &! INTENT(IN) !global
       solai        , &! INTENT(IN) !global
       fira         , &! INTENT(IN) !global
       raina        , &! INTENT(IN) !global
       qa           , &! INTENT(IN) !global
       psurf        , &! INTENT(IN) !global
       snowa        , &! INTENT(IN) !global
       ua           , &! INTENT(IN) !global
       o2conc       , &! INTENT(IN) !global
       co2conc      ,&! INTENT(IN) !global
       fwetu        , &! INTENT(INOUT)! local
       rliqu        , &! INTENT(INOUT)! local
       fwets        , &! INTENT(INOUT)! local
       rliqs        , &! INTENT(INOUT)! local
       fwetl        , &! INTENT(INOUT)! local
       rliql        , &! INTENT(INOUT)! local
                                !nsol         , &! INTENT(INOUT) !local
       solu         , &! INTENT(INOUT)! local
       sols         , &! INTENT(INOUT)! local
       soll         ,&! INTENT(INOUT)! local
       solg              , &! INTENT(INOUT)! local
       soli         , &! INTENT(INOUT)! local
       scalcoefl    , &! INTENT(INOUT)! local
       scalcoefu    , &! INTENT(INOUT)! local
       indsol       , &! INTENT(INOUT)! local
       albsod       , &! INTENT(INOUT)! local
       albsoi       , &! INTENT(INOUT)! local
       albsnd       , &! INTENT(INOUT)! local
       albsni       , &! INTENT(INOUT)! local
       relod        , &! INTENT(OUT  )! local
       reloi              , &! INTENT(OUT  ) !local
       reupd        , &! INTENT(OUT  ) !local
       reupi        , &! INTENT(OUT  ) !local
       ablod        , &! INTENT(INOUT)! local
       abloi        , &! INTENT(INOUT)! local
       flodd        , &! INTENT(INOUT)! local
       dummy        , &! INTENT(OUT  ) !local
       flodi        , &! INTENT(INOUT)! local
       floii        , &! INTENT(INOUT)! local
       terml        , &! INTENT(INOUT)! local
       termu              , &! INTENT(INOUT)! local
       abupd        , &! INTENT(INOUT)! local
       abupi        , &! INTENT(INOUT)! local
       fupdd        , &! INTENT(INOUT)! local
       fupdi        , &! INTENT(INOUT)! local
       fupii        , &! INTENT(INOUT)! local
       sol2d        , &! INTENT(OUT  ) !local
       sol2i        , &! INTENT(OUT  ) !local
       sol3d        , &! INTENT(OUT  ) !local
       sol3i        , &! INTENT(OUT  ) !local
       firb         , &! INTENT(INOUT)! local
       firs         , &! INTENT(INOUT)! local
       firu         , &! INTENT(INOUT)! local
       firl         , &! INTENT(INOUT)! local
       firg         , &! INTENT(INOUT)! local
       firi         , &! INTENT(INOUT)! local
       snowg        , &! INTENT(INOUT)! local
       tsnowg       , &! INTENT(INOUT)! local
       tsnowl       , &! INTENT(INOUT)! local
       pfluxl       , &! INTENT(INOUT)! local
       raing        , &! INTENT(INOUT)! local
       traing       , &! INTENT(INOUT)! local
       trainl       , &! INTENT(INOUT)! local
       snowl        , &! INTENT(INOUT)! local
       tsnowu       , &! INTENT(INOUT)! local
       pfluxu       , &! INTENT(INOUT)! local
       rainu        , &! INTENT(INOUT)! local
       trainu       , &! INTENT(INOUT)! local
       snowu        , &! INTENT(INOUT)! local
       pfluxs       , &! INTENT(INOUT)! local
       rainl        , &! INTENT(INOUT)! local
       bps          , &! INTENT(INOUT)! local
       rhoa         , &! INTENT(INOUT) !local
       cp           , &! INTENT(INOUT) !local
       za           , &! INTENT(INOUT) !local
       bdl          , &! INTENT(INOUT) !local
       dil          , &! INTENT(INOUT) !local
       z3           , &! INTENT(INOUT) !local
       z4           , &! INTENT(INOUT) !local
       z34          , &! INTENT(INOUT) !local
       exphl        , &! INTENT(INOUT) !local
       expl         , &! INTENT(INOUT) !local
       displ        , &! INTENT(OUT  ) !local
       bdu          , &! INTENT(INOUT) !local
       diu          , &! INTENT(INOUT) !local
       z1           , &! INTENT(INOUT) !local
       z2           , &! INTENT(INOUT) !local
       z12          , &! INTENT(INOUT) !local
       exphu        , &! INTENT(INOUT) !local
       expu         , &! INTENT(INOUT) !local
       dispu        , &! INTENT(OUT  ) !local
       alogg        , &! INTENT(OUT  ) !local
       alogi        , &! INTENT(OUT  ) !local
       alogav       , &! INTENT(INOUT) !local
       alog4        , &! INTENT(INOUT) !local
       alog3        , &! INTENT(INOUT) !local
       alog2        , &! INTENT(OUT  ) !local
       alog1        , &! INTENT(INOUT) !local
       aloga        , &! INTENT(INOUT) !local
       u2           , &! INTENT(INOUT) !local
       alogu        , &! INTENT(INOUT) !local
       alogl        , &! INTENT(INOUT) !local
       richl        , &! INTENT(OUT  ) !local
       straml       , &! INTENT(OUT  ) !local
       strahl       , &! INTENT(OUT  ) !local
       richu        , &! INTENT(OUT  ) !local
       stramu       , &! INTENT(OUT  ) !local
       strahu       , &! INTENT(OUT  ) !local
       u1           , &! INTENT(OUT  ) !local
       u12          , &! INTENT(OUT  ) !local
       u3           , &! INTENT(OUT  ) !local
       u34          , &! INTENT(OUT  ) !local
       u4           , &! INTENT(OUT  ) !local
       cu           , &! INTENT(INOUT) !local
       cl           , &! INTENT(INOUT) !local
       sg           , &! INTENT(INOUT) !local
       si           , &! INTENT(INOUT) !local
       fwetux       , &! INTENT(OUT  ) !local
       fwetsx       , &! INTENT(OUT  ) !local
       fwetlx       , &! INTENT(OUT  ) !local
       fsena        , &! INTENT(INOUT) !local
       fseng        , &! INTENT(OUT  ) !local
       fseni        , &! INTENT(OUT  ) !local
       fsenu        , &! INTENT(OUT  ) !local
       fsens        , &! INTENT(OUT  ) !local
       fsenl              , &! INTENT(OUT  ) !local
       fvapa        , &! INTENT(INOUT) !local
       fvaput       , &! INTENT(OUT  ) !local
       fvaps        , &! INTENT(INOUT) !local
       fvaplw       , &! INTENT(INOUT) !local
       fvaplt       , &! INTENT(OUT  ) !local
       fvapg        , &! INTENT(INOUT) !local
       fvapi        , &! INTENT(INOUT) !local
       fvapuw       , &! INTENT(INOUT) !local
       npoi         , &! INTENT(IN) !global
       nband              , &! INTENT(IN) !global
       nsoilay      , &! INTENT(IN) !global
       nsnolay      , &! INTENT(IN) !global
       npft              , &! INTENT(IN) !global
       epsilon      , &! INTENT(IN) !global
       dtime              , &! INTENT(IN) !global
       stef         , &! INTENT(IN) !global
       vonk         , &! INTENT(IN) !global
       grav         , &! INTENT(IN) !global
       tmelt        , &! INTENT(IN) !global
       hfus         , &! INTENT(IN) !global
       hvap         , &! INTENT(IN) !global
       hsub         , &! INTENT(IN) !global
       ch2o         , &! INTENT(IN) !global
       cice         , &! INTENT(IN) !global
       cair         , &! INTENT(IN) !global
       cvap         , &! INTENT(IN) !global
       rair         , &! INTENT(IN) !global
       rvap         , &! INTENT(IN) !global
       cappa        , &! INTENT(IN) !global
       rhow         , &! INTENT(IN) !global
       vzero        , &! INTENT(IN) !global
       pi           , &! INTENT(IN) !global
       ux           , &! INTENT(IN    ) !global
       uy           , &! INTENT(IN    ) !global
       taux         , &! INTENT(OUT   ) !local
       tauy         , &! INTENT(OUT   ) !local
       bstar        , &! INTENT(OUT   ) !local
       ts2          , &! INTENT(OUT   ) !local
       qs2          , &! INTENT(OUT   ) !local
       vegtype0     , &! INTENT(IN   ) !local
       stressfac    , &
       avmuir_factor, &
       nVegClass     )




    ! ---------------------------------------------------------------
    !
    ! common blocks
    !
    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN   ) :: pi 
    INTEGER, INTENT(IN   )  :: nVegClass


    INTEGER, INTENT(IN   )  :: npoi                ! total number of land points
    INTEGER, INTENT(IN   )  :: nband               ! number of solar radiation wavebands
    INTEGER, INTENT(IN   )  :: nsoilay             ! number of soil layers
    INTEGER, INTENT(IN   )  :: nsnolay             ! number of snow layers
    INTEGER, INTENT(IN   )  :: npft                     ! number of plant functional types
    REAL(KIND=r8), INTENT(IN   )  :: epsilon             ! small quantity to avoid zero-divides and other
    ! truncation or machine-limit troubles with small
    ! values. should be slightly greater than o(1)
    ! machine precision
    REAL(KIND=r8), INTENT(IN   )  :: dtime               ! model timestep (seconds)
    REAL(KIND=r8), INTENT(IN   )  :: stef                ! stefan-boltzmann constant (W m-2 K-4)
    REAL(KIND=r8), INTENT(IN   )  :: vonk                ! von karman constant (dimensionless)
    REAL(KIND=r8), INTENT(IN   )  :: grav                ! gravitational acceleration (m s-2)
    REAL(KIND=r8), INTENT(IN   )  :: tmelt               ! freezing point of water (K)
    REAL(KIND=r8), INTENT(IN   )  :: hfus                ! latent heat of fusion of water (J kg-1)
    REAL(KIND=r8), INTENT(IN   )  :: hvap                ! latent heat of vaporization of water (J kg-1)
    REAL(KIND=r8), INTENT(IN   )  :: hsub                ! latent heat of sublimation of ice (J kg-1)
    REAL(KIND=r8), INTENT(IN   )  :: ch2o             ! specific heat of liquid water (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN   )  :: cice             ! specific heat of ice (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN   )  :: cair             ! specific heat of dry air at constant pressure (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN   )  :: cvap             ! specific heat of water vapor at constant pressure (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN   )  :: rair             ! gas constant for dry air (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN   )  :: rvap             ! gas constant for water vapor (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN   )  :: cappa               ! rair/cair
    REAL(KIND=r8), INTENT(IN   )  :: rhow             ! density of liquid water (all types) (kg m-3)
    REAL(KIND=r8), INTENT(IN   )  :: vzero   (npoi)      ! a real array of zeros, of length npoi

    !      INCLUDE 'com1d.h'
    REAL(KIND=r8) , INTENT(INOUT) :: fwetu   (npoi)      ! fraction of upper canopy leaf area wetted by intercepted liquid and/or snow
    REAL(KIND=r8) , INTENT(INOUT) :: rliqu   (npoi)      ! proportion of fwetu due to liquid
    REAL(KIND=r8) , INTENT(INOUT) :: fwets   (npoi)      ! fraction of upper canopy stem area wetted by intercepted liquid and/or snow
    REAL(KIND=r8) , INTENT(INOUT) :: rliqs   (npoi)      ! proportion of fwets due to liquid
    REAL(KIND=r8) , INTENT(INOUT) :: fwetl   (npoi)      ! fraction of lower canopy stem & leaf area wetted by intercepted liquid and/or snow
    REAL(KIND=r8) , INTENT(INOUT) :: rliql   (npoi)      ! proportion of fwetl due to liquid
    INTEGER  :: nsol                ! number of points in indsol
    REAL(KIND=r8) , INTENT(INOUT) :: solu   (npoi)       ! solar flux (direct + diffuse) absorbed by upper canopy leaves per unit canopy area (W m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: sols   (npoi)       ! solar flux (direct + diffuse) absorbed by upper canopy stems per unit canopy area (W m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: soll   (npoi)       ! solar flux (direct + diffuse) absorbed by lower canopy leaves and stems per unit canopy area (W m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: solg   (npoi)       ! solar flux (direct + diffuse) absorbed by unit snow-free soil (W m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: soli   (npoi)       ! solar flux (direct + diffuse) absorbed by unit snow surface (W m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: scalcoefl(npoi,4)   ! term needed in lower canopy scaling
    REAL(KIND=r8) , INTENT(INOUT) :: scalcoefu(npoi,4)   ! term needed in upper canopy scaling
    INTEGER , INTENT(INOUT) :: indsol (npoi)         ! index of current strip for points with positive coszen
    REAL(KIND=r8) , INTENT(INOUT) :: albsod (npoi)       ! direct  albedo for soil surface (visible or IR)
    REAL(KIND=r8) , INTENT(INOUT) :: albsoi (npoi)       ! diffuse albedo for soil surface (visible or IR)
    REAL(KIND=r8) , INTENT(INOUT) :: albsnd (npoi)       ! direct  albedo for snow surface (visible or IR)
    REAL(KIND=r8) , INTENT(INOUT) :: albsni (npoi)       ! diffuse albedo for snow surface (visible or IR)
    REAL(KIND=r8) , INTENT(OUT  ) :: relod  (npoi)       ! upward direct radiation per unit icident direct beam on lower canopy (W m-2)
    REAL(KIND=r8) , INTENT(OUT  ) :: reloi  (npoi)       ! upward diffuse radiation per unit incident diffuse radiation on lower canopy (W m-2)
    REAL(KIND=r8) , INTENT(OUT  ) :: reupd  (npoi)       ! upward direct radiation per unit incident direct radiation on upper canopy (W m-2)
    REAL(KIND=r8) , INTENT(OUT  ) :: reupi  (npoi)       ! upward diffuse radiation per unit incident diffuse radiation on upper canopy (W m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: ablod  (npoi)       ! fraction of direct  radiation absorbed by lower canopy
    REAL(KIND=r8) , INTENT(INOUT) :: abloi  (npoi)       ! fraction of diffuse radiation absorbed by lower canopy
    REAL(KIND=r8) , INTENT(INOUT) :: flodd  (npoi)       ! downward direct radiation per unit incident direct radiation on lower canopy (W m-2)
    REAL(KIND=r8) , INTENT(OUT  ) :: dummy  (npoi)       ! placeholder, always = 0: no direct flux produced for diffuse incident
    REAL(KIND=r8) , INTENT(INOUT) :: flodi  (npoi)       ! downward diffuse radiation per unit incident direct radiation on lower canopy (W m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: floii  (npoi)       ! downward diffuse radiation per unit incident diffuse radiation on lower canopy
    REAL(KIND=r8) , INTENT(INOUT) :: terml  (npoi,7)     ! term needed in lower canopy scaling
    REAL(KIND=r8) , INTENT(INOUT) :: termu  (npoi,7)     ! term needed in upper canopy scaling
    REAL(KIND=r8) , INTENT(INOUT) :: abupd  (npoi)       ! fraction of direct  radiation absorbed by upper canopy
    REAL(KIND=r8) , INTENT(INOUT) :: abupi  (npoi)       ! fraction of diffuse radiation absorbed by upper canopy
    REAL(KIND=r8) , INTENT(INOUT) :: fupdd  (npoi)       ! downward direct radiation per unit incident direct beam on upper canopy (W m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: fupdi  (npoi)       ! downward diffuse radiation per unit icident direct radiation on upper canopy (W m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: fupii  (npoi)       ! downward diffuse radiation per unit incident diffuse radiation on upper canopy (W m-2)
    REAL(KIND=r8) , INTENT(OUT  ) :: sol2d  (npoi)       ! direct downward radiation  out of upper canopy per unit vegetated (upper) area (W m-2)
    REAL(KIND=r8) , INTENT(OUT  ) :: sol2i  (npoi)       ! diffuse downward radiation out of upper canopy per unit vegetated (upper) area(W m-2)
    REAL(KIND=r8) , INTENT(OUT  ) :: sol3d  (npoi)       ! direct downward radiation  out of upper canopy + gaps per unit grid cell area (W m-2)
    REAL(KIND=r8) , INTENT(OUT  ) :: sol3i  (npoi)       ! diffuse downward radiation out of upper canopy + gaps per unit grid cell area (W m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: firb   (npoi)       ! net upward ir radiation at reference atmospheric level za (W m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: firs   (npoi)       ! ir radiation absorbed by upper canopy stems (W m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: firu   (npoi)       ! ir raditaion absorbed by upper canopy leaves (W m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: firl   (npoi)       ! ir radiation absorbed by lower canopy leaves and stems (W m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: firg   (npoi)       ! ir radiation absorbed by soil/ice (W m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: firi   (npoi)       ! ir radiation absorbed by snow (W m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: snowg  (npoi)       ! snowfall rate at soil level (kg h2o m-2 s-1)
    REAL(KIND=r8) , INTENT(INOUT) :: tsnowg (npoi)       ! snowfall temperature at soil level (K) 
    REAL(KIND=r8) , INTENT(INOUT) :: tsnowl (npoi)       ! snowfall temperature below upper canopy (K)
    REAL(KIND=r8) , INTENT(INOUT) :: pfluxl (npoi)       ! heat flux on lower canopy leaves & stems due to intercepted h2o (W m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: raing  (npoi)       ! rainfall rate at soil level (kg m-2 s-1)
    REAL(KIND=r8) , INTENT(INOUT) :: traing (npoi)       ! rainfall temperature at soil level (K)
    REAL(KIND=r8) , INTENT(INOUT) :: trainl (npoi)       ! rainfall temperature below upper canopy (K)
    REAL(KIND=r8) , INTENT(INOUT) :: snowl  (npoi)       ! snowfall rate below upper canopy (kg h2o m-2 s-1)
    REAL(KIND=r8) , INTENT(INOUT) :: tsnowu (npoi)       ! snowfall temperature above upper canopy (K)
    REAL(KIND=r8) , INTENT(INOUT) :: pfluxu (npoi)       ! heat flux on upper canopy leaves due to intercepted h2o (W m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: rainu  (npoi)       ! rainfall rate above upper canopy (kg m-2 s-1)
    REAL(KIND=r8) , INTENT(INOUT) :: trainu (npoi)       ! rainfall temperature above upper canopy (K)
    REAL(KIND=r8) , INTENT(INOUT) :: snowu  (npoi)       ! snowfall rate above upper canopy (kg h2o m-2 s-1)
    REAL(KIND=r8) , INTENT(INOUT) :: pfluxs (npoi)       ! heat flux on upper canopy stems due to intercepted h2o (W m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: rainl  (npoi)       ! rainfall rate below upper canopy (kg m-2 s-1)
    REAL(KIND=r8) , INTENT(INOUT) :: bps    (npoi)             ! (ps/p) ** (rair/cair) for atmospheric level  (const)
    REAL(KIND=r8) , INTENT(INOUT) :: rhoa   (npoi)       ! air density at za (allowing for h2o vapor) (kg m-3)
    REAL(KIND=r8) , INTENT(INOUT) :: cp     (npoi)       ! specific heat of air at za (allowing for h2o vapor) (J kg-1 K-1)
    REAL(KIND=r8) , INTENT(INOUT) :: za     (npoi)       ! height above the surface of atmospheric forcing (m)
    REAL(KIND=r8) , INTENT(INOUT) :: bdl    (npoi)       ! aerodynamic coefficient ([(tau/rho)/u**2] for laower canopy (A31/A30 Pollard & Thompson 1995)
    REAL(KIND=r8) , INTENT(INOUT) :: dil    (npoi)       ! inverse of momentum diffusion coefficient within lower canopy (m)
    REAL(KIND=r8) , INTENT(INOUT) :: z3     (npoi)       ! effective top of the lower canopy (for momentum) (m)
    REAL(KIND=r8) , INTENT(INOUT) :: z4     (npoi)       ! effective bottom of the lower canopy (for momentum) (m)
    REAL(KIND=r8) , INTENT(INOUT) :: z34    (npoi)       ! effective middle of the lower canopy (for momentum) (m)
    REAL(KIND=r8) , INTENT(INOUT) :: exphl  (npoi)       ! exp(lamda/2*(z3-z4)) for lower canopy (A30 Pollard & Thompson)
    REAL(KIND=r8) , INTENT(INOUT) :: expl   (npoi)       ! exphl**2
    REAL(KIND=r8) , INTENT(OUT  ) :: displ  (npoi)       ! zero-plane displacement height for lower canopy (m)
    REAL(KIND=r8) , INTENT(INOUT) :: bdu    (npoi)       ! aerodynamic coefficient ([(tau/rho)/u**2] for upper canopy (A31/A30 Pollard & Thompson 1995)
    REAL(KIND=r8) , INTENT(INOUT) :: diu    (npoi)       ! inverse of momentum diffusion coefficient within upper canopy (m)
    REAL(KIND=r8) , INTENT(INOUT) :: z1     (npoi)       ! effective top of upper canopy (for momentum) (m)
    REAL(KIND=r8) , INTENT(INOUT) :: z2     (npoi)       ! effective bottom of the upper canopy (for momentum) (m)
    REAL(KIND=r8) , INTENT(INOUT) :: z12    (npoi)       ! effective middle of the upper canopy (for momentum) (m)
    REAL(KIND=r8) , INTENT(INOUT) :: exphu  (npoi)       ! exp(lamda/2*(z3-z4)) for upper canopy (A30 Pollard & Thompson)
    REAL(KIND=r8) , INTENT(INOUT) :: expu   (npoi)       ! exphu**2
    REAL(KIND=r8) , INTENT(OUT  ) :: dispu  (npoi)       ! zero-plane displacement height for upper canopy (m)
    REAL(KIND=r8) , INTENT(OUT  ) :: alogg  (npoi)       ! log of soil roughness
    REAL(KIND=r8) , INTENT(OUT  ) :: alogi  (npoi)       ! log of snow roughness
    REAL(KIND=r8) , INTENT(INOUT) :: alogav (npoi)       ! average of alogi and alogg 
    REAL(KIND=r8) , INTENT(INOUT) :: alog4  (npoi)       ! log (max(z4, 1.1*z0sno, 1.1*z0soi)) 
    REAL(KIND=r8) , INTENT(INOUT) :: alog3  (npoi)       ! log (z3 - displ)
    REAL(KIND=r8) , INTENT(OUT  ) :: alog2  (npoi)       ! log (z2 - displ)
    REAL(KIND=r8) , INTENT(INOUT) :: alog1  (npoi)       ! log (z1 - dispu) 
    REAL(KIND=r8) , INTENT(INOUT) :: aloga  (npoi)       ! log (za - dispu) 
    REAL(KIND=r8) , INTENT(INOUT) :: u2     (npoi)       ! wind speed at level z2 (m s-1)
    REAL(KIND=r8) , INTENT(INOUT) :: alogu  (npoi)       ! log (roughness length of upper canopy)
    REAL(KIND=r8) , INTENT(INOUT) :: alogl  (npoi)       ! log (roughness length of lower canopy)
    REAL(KIND=r8) , INTENT(OUT  ) :: richl  (npoi)       ! richardson number for air above upper canopy (z3 to z2)
    REAL(KIND=r8) , INTENT(OUT  ) :: straml (npoi)       ! momentum correction factor for stratif between upper & lower canopy (z3 to z2) (louis et al.)
    REAL(KIND=r8) , INTENT(OUT  ) :: strahl (npoi)       ! heat/vap correction factor for stratif between upper & lower canopy (z3 to z2) (louis et al.)
    REAL(KIND=r8) , INTENT(OUT  ) :: richu  (npoi)       ! richardson number for air between upper & lower canopy (z1 to za)
    REAL(KIND=r8) , INTENT(OUT  ) :: stramu (npoi)       ! momentum correction factor for stratif above upper canopy (z1 to za) (louis et al.)
    REAL(KIND=r8) , INTENT(OUT  ) :: strahu (npoi)       ! heat/vap correction factor for stratif above upper canopy (z1 to za) (louis et al.)
    REAL(KIND=r8) , INTENT(OUT  ) :: u1     (npoi)       ! wind speed at level z1 (m s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: u12    (npoi)       ! wind speed at level z12 (m s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: u3     (npoi)       ! wind speed at level z3 (m s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: u34    (npoi)       ! wind speed at level z34 (m s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: u4     (npoi)       ! wind speed at level z4 (m s-1)
    REAL(KIND=r8) , INTENT(INOUT) :: cu     (npoi)       ! air transfer coefficient (*rhoa) (m s-1 kg m-3) for upper air region (z12 --> za) (A35 Pollard & Thompson 1995)
    REAL(KIND=r8) , INTENT(INOUT) :: cl     (npoi)       ! air transfer coefficient (*rhoa) (m s-1 kg m-3) between the 2 canopies (z34 --> z12) (A36 Pollard & Thompson 1995)
    REAL(KIND=r8) , INTENT(INOUT) :: sg     (npoi)       ! air-soil transfer coefficient
    REAL(KIND=r8) , INTENT(INOUT) :: si     (npoi)       ! air-snow transfer coefficient
    REAL(KIND=r8) , INTENT(OUT  ) :: fwetux (npoi)       ! fraction of upper canopy leaf area wetted if dew forms
    REAL(KIND=r8) , INTENT(OUT  ) :: fwetsx (npoi)       ! fraction of upper canopy stem area wetted if dew forms
    REAL(KIND=r8) , INTENT(OUT  ) :: fwetlx (npoi)       ! fraction of lower canopy leaf and stem area wetted if dew forms
    REAL(KIND=r8) , INTENT(INOUT) :: fsena  (npoi)       ! downward sensible heat flux between za & z12 at za (W m-2)
    REAL(KIND=r8) , INTENT(OUT  ) :: fseng  (npoi)       ! upward sensible heat flux between soil surface & air at z34 (W m-2)
    REAL(KIND=r8) , INTENT(OUT  ) :: fseni  (npoi)       ! upward sensible heat flux between snow surface & air at z34 (W m-2)
    REAL(KIND=r8) , INTENT(OUT  ) :: fsenu  (npoi)       ! sensible heat flux from upper canopy leaves to air (W m-2)
    REAL(KIND=r8) , INTENT(OUT  ) :: fsens  (npoi)       ! sensible heat flux from upper canopy stems to air (W m-2)
    REAL(KIND=r8) , INTENT(OUT  ) :: fsenl  (npoi)       ! sensible heat flux from lower canopy to air (W m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: fvapa  (npoi)       ! downward h2o vapor flux between za & z12 at za (kg m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: fvaput (npoi)       ! h2o vapor flux (transpiration from dry parts) between upper canopy leaves and air at z12 (kg m-2 s-1/ LAI upper canopy/ fu)
    REAL(KIND=r8) , INTENT(INOUT) :: fvaps  (npoi)       ! h2o vapor flux (evaporation from wet surface) between upper canopy stems and air at z12 (kg m-2 s-1 / SAI lower canopy / fu)
    REAL(KIND=r8) , INTENT(INOUT) :: fvaplw (npoi)       ! h2o vapor flux (evaporation from wet surface) between lower canopy leaves & stems and air at z34 (kg m-2 s-1/ LAI lower canopy/ fl)
    REAL(KIND=r8) , INTENT(OUT  ) :: fvaplt (npoi)       ! h2o vapor flux (transpiration) between lower canopy & air at z34 (kg m-2 s-1 / LAI lower canopy / fl)
    REAL(KIND=r8) , INTENT(INOUT) :: fvapg  (npoi)       ! h2o vapor flux (evaporation) between soil & air at z34 (kg m-2 s-1/bare ground fraction)
    REAL(KIND=r8) , INTENT(INOUT) :: fvapi  (npoi)       ! h2o vapor flux (evaporation) between snow & air at z34 (kg m-2 s-1 / fi )
    REAL(KIND=r8) , INTENT(INOUT) :: fvapuw (npoi)       ! h2o vapor flux (evaporation from wet parts) between upper canopy leaves and air at z12 (kg m-2 s-1/ LAI upper canopy/ fu)

    !      INCLUDE 'comatm.h'
    REAL(KIND=r8) , INTENT(IN   ) :: ta     (npoi)       ! air temperature (K)
    REAL(KIND=r8) , INTENT(INOUT) :: asurd  (npoi,nband) ! direct albedo of surface system
    REAL(KIND=r8) , INTENT(INOUT) :: asuri  (npoi,nband) ! diffuse albedo of surface system 
    REAL(KIND=r8) , INTENT(IN   ) :: coszen (npoi)       ! cosine of solar zenith angle
    REAL(KIND=r8) , INTENT(IN   ) :: solad  (npoi,nband) ! direct downward solar flux (W m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: solai  (npoi,nband) ! diffuse downward solar flux (W m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: fira   (npoi)       ! incoming ir flux (W m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: raina  (npoi)       ! rainfall rate (mm/s or kg m-2 s-1)
    REAL(KIND=r8) , INTENT(IN   ) :: qa     (npoi)       ! specific humidity (kg_h2o/kg_air)
    REAL(KIND=r8) , INTENT(IN   ) :: psurf  (npoi)       ! surface pressure (Pa)
    REAL(KIND=r8) , INTENT(IN   ) :: snowa  (npoi)       ! snowfall rate (mm/s or kg m-2 s-1 of water)
    REAL(KIND=r8) , INTENT(IN   ) :: ua     (npoi)       ! wind speed (m s-1)
    REAL(KIND=r8) , INTENT(IN   ) :: o2conc              ! o2 concentration (mol/mol)
    REAL(KIND=r8) , INTENT(IN   ) :: co2conc(npoi)             ! co2 concentration (mol/mol)
    !      INCLUDE 'comsoi.h'
    REAL(KIND=r8) , INTENT(IN   ) :: sand    (npoi,nsoilay)    ! percent sand of soil
    REAL(KIND=r8) , INTENT(IN   ) :: clay    (npoi,nsoilay)    ! percent clay of soil
    REAL(KIND=r8) , INTENT(IN   ) :: poros   (npoi,nsoilay)    ! porosity (mass of h2o per unit vol at sat / rhow)
    REAL(KIND=r8) , INTENT(INOUT) :: wsoi    (npoi,nsoilay)    ! fraction of soil pore space containing liquid water
    REAL(KIND=r8) , INTENT(INOUT) :: wisoi   (npoi,nsoilay)    ! fraction of soil pore space containing ice
    REAL(KIND=r8) , INTENT(INOUT) :: consoi  (npoi,nsoilay)    ! thermal conductivity of each soil layer (W m-1 K-1)
    REAL(KIND=r8) , INTENT(IN   ) :: zwpmax                    ! assumed maximum fraction of soil surface 
    ! covered by puddles (dimensionless)
    REAL(KIND=r8) , INTENT(INOUT) :: wpud    (npoi)            ! liquid content of puddles per soil area (kg m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: wipud   (npoi)            ! ice content of puddles per soil area (kg m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: wpudmax                   ! normalization constant for puddles (kg m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: qglif   (npoi,4)          ! 1: fraction of soil evap (fvapg) from soil liquid
    ! 2: fraction of soil evap (fvapg) from soil ice
    ! 3: fraction of soil evap (fvapg) from puddle liquid
    ! 4: fraction of soil evap (fvapg) from puddle ice
    REAL(KIND=r8) , INTENT(INOUT) :: tsoi    (npoi,nsoilay)    ! soil temperature for each layer (K)
    REAL(KIND=r8) , INTENT(INOUT) :: hvasug  (npoi)            ! latent heat of vap/subl, for soil surface (J kg-1)
    REAL(KIND=r8) , INTENT(INOUT) :: hvasui  (npoi)            ! latent heat of vap/subl, for snow surface (J kg-1)
    REAL(KIND=r8) , INTENT(IN   ) :: albsav  (npoi)            ! saturated soil surface albedo (visible waveband)
    REAL(KIND=r8) , INTENT(IN   ) :: albsan  (npoi)            ! saturated soil surface albedo (near-ir waveband)
    REAL(KIND=r8) , INTENT(INOUT) :: tg      (npoi)            ! soil skin temperature (K)
    REAL(KIND=r8) , INTENT(INOUT) :: ti      (npoi)            ! snow skin temperature (K)
    REAL(KIND=r8) , INTENT(IN   ) :: z0soi   (npoi)            ! roughness length of soil surface (m)
    REAL(KIND=r8) , INTENT(IN   ) :: swilt   (npoi,nsoilay)    ! wilting soil moisture value (fraction of pore space)
    REAL(KIND=r8) , INTENT(IN   ) :: sfield  (npoi,nsoilay)    ! field capacity soil moisture value (fraction of pore space)
    REAL(KIND=r8) , INTENT(INOUT) :: stressl (npoi,nsoilay)    ! soil moisture stress factor for the lower canopy (dimensionless)
    REAL(KIND=r8) , INTENT(INOUT) :: stressu (npoi,nsoilay)    ! soil moisture stress factor for the upper canopy (dimensionless)
    REAL(KIND=r8) , INTENT(INOUT) :: stresstl(npoi)            ! sum of stressl over all 6 soil layers (dimensionless)
    REAL(KIND=r8) , INTENT(INOUT) :: stresstu(npoi)            ! sum of stressu over all 6 soil layers (dimensionless)
    REAL(KIND=r8) , INTENT(IN   ) :: csoi    (npoi,nsoilay)    ! specific heat of soil, no pore spaces (J kg-1 deg-1)
    REAL(KIND=r8) , INTENT(IN   ) :: rhosoi  (npoi,nsoilay)    ! soil density (without pores, not bulk) (kg m-3)
    REAL(KIND=r8) , INTENT(IN   ) :: hsoi    (npoi,nsoilay+1)       ! soil layer thickness (m)
    REAL(KIND=r8) , INTENT(IN   ) :: suction (npoi,nsoilay)    ! saturated matric potential (m-h2o)
    REAL(KIND=r8) , INTENT(IN   ) :: bex     (npoi,nsoilay)    ! exponent "b" in soil water potential
    REAL(KIND=r8) , INTENT(INOUT) :: upsoiu  (npoi,nsoilay)    ! soil water uptake from transpiration (kg_h2o m-2 s-1)
    REAL(KIND=r8) , INTENT(INOUT) :: upsoil  (npoi,nsoilay)    ! soil water uptake from transpiration (kg_h2o m-2 s-1)
    REAL(KIND=r8) , INTENT(INOUT) :: heatg   (npoi)            ! net heat flux into soil surface (W m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: heati   (npoi)            ! net heat flux into snow surface (W m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: hydraul (npoi,nsoilay)    ! saturated hydraulic conductivity (m/s)
    REAL(KIND=r8) , INTENT(INOUT) :: porosflo(npoi,nsoilay)    ! porosity after reduction by ice content
    INTEGER , INTENT(IN   ) :: ibex    (npoi,nsoilay)    ! nint(bex), used for cpu speed
    REAL(KIND=r8) , INTENT(IN   ) :: bperm (npoi)                     ! lower b.c. for soil profile drainage 
    ! (0.0 = impermeable; 1.0 = fully permeable)
    REAL(KIND=r8) , INTENT(INOUT) :: hflo    (npoi,nsoilay+1)  ! downward heat transport through soil layers (W m-2)

    !      INCLUDE 'comsat.h'    
    !      include 'comsno.h'
    REAL(KIND=r8) , INTENT(IN   ) :: z0sno                        ! roughness length of snow surface (m)
    REAL(KIND=r8) , INTENT(IN   ) :: rhos                         ! density of snow (kg m-3)
    REAL(KIND=r8) , INTENT(IN   ) :: consno                       ! thermal conductivity of snow (W m-1 K-1)
    REAL(KIND=r8) , INTENT(IN   ) :: hsnotop               ! thickness of top snow layer (m)
    REAL(KIND=r8) , INTENT(IN   ) :: hsnomin               ! minimum total thickness of snow (m)
    REAL(KIND=r8) , INTENT(IN   ) :: fimin                        ! minimum fractional snow cover
    REAL(KIND=r8) , INTENT(IN   ) :: fimax                 ! maximum fractional snow cover
    REAL(KIND=r8) , INTENT(INOUT) :: fi     (npoi)               ! fractional snow cover
    REAL(KIND=r8) , INTENT(INOUT) :: tsno   (npoi,nsnolay) ! temperature of snow layers (K)
    REAL(KIND=r8) , INTENT(INOUT) :: hsno   (npoi,nsnolay) ! thickness of snow layers (m)

    !      include 'comveg.h'
    REAL(KIND=r8) , INTENT(INOUT) :: wliqu    (npoi)     ! intercepted liquid h2o on upper canopy leaf area (kg m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: wliqumax            ! maximum intercepted water on a unit upper canopy leaf area (kg m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: wsnou    (npoi)     ! intercepted frozen h2o (snow) on upper canopy leaf area (kg m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: wsnoumax            ! intercepted snow capacity for upper canopy leaves (kg m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: tu       (npoi)     ! temperature of upper canopy leaves (K)
    REAL(KIND=r8) , INTENT(INOUT) :: wliqs    (npoi)     ! intercepted liquid h2o on upper canopy stem area (kg m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: wliqsmax            ! maximum intercepted water on a unit upper canopy stem area (kg m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: wsnos    (npoi)     ! intercepted frozen h2o (snow) on upper canopy stem area (kg m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: wsnosmax            ! intercepted snow capacity for upper canopy stems (kg m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: ts       (npoi)     ! temperature of upper canopy stems (K)
    REAL(KIND=r8) , INTENT(INOUT) :: wliql    (npoi)     ! intercepted liquid h2o on lower canopy leaf and stem area (kg m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: wliqlmax            ! maximum intercepted water on a unit lower canopy stem & leaf area (kg m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: wsnol    (npoi)     ! intercepted frozen h2o (snow) on lower canopy leaf & stem area (kg m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: wsnolmax            ! intercepted snow capacity for lower canopy leaves & stems (kg m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: tl       (npoi)     ! temperature of lower canopy leaves & stems(K)
    REAL(KIND=r8) , INTENT(INOUT) :: topparu  (npoi)     ! total photosynthetically active raditaion absorbed by top leaves of upper canopy (W m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: topparl  (npoi)     ! total photosynthetically active raditaion absorbed by top leaves of lower canopy (W m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: fl       (npoi)     ! fraction of snow-free area covered by lower  canopy
    REAL(KIND=r8) , INTENT(IN   ) :: fu       (npoi)     ! fraction of overall area covered by upper canopy
    REAL(KIND=r8) , INTENT(IN   ) :: lai      (npoi,2)   ! canopy single-sided leaf area index (area leaf/area veg)
    REAL(KIND=r8) , INTENT(IN   ) :: sai      (npoi,2)   ! current single-sided stem area index
    REAL(KIND=r8) , INTENT(IN   ) :: rhoveg   (nband,2)  ! reflectance of an average leaf/stem
    REAL(KIND=r8) , INTENT(IN   ) :: tauveg   (nband,2)  ! transmittance of an average leaf/stem
    REAL(KIND=r8) , INTENT(IN   ) :: orieh    (2)        ! fraction of leaf/stems with horizontal orientation
    REAL(KIND=r8) , INTENT(IN   ) :: oriev    (2)        ! fraction of leaf/stems with vertical
    REAL(KIND=r8) , INTENT(INOUT) :: wliqmin             ! minimum intercepted water on unit vegetated area (kg m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: wsnomin             ! minimum intercepted snow on unit vegetated area (kg m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: t12      (npoi)     ! air temperature at z12 (K)
    REAL(KIND=r8) , INTENT(IN   ) :: tdripu              ! decay time for dripoff of liquid intercepted by upper canopy leaves (sec)
    REAL(KIND=r8) , INTENT(IN   ) :: tblowu              ! decay time for blowoff of snow intercepted by upper canopy leaves (sec)
    REAL(KIND=r8) , INTENT(IN   ) :: tdrips              ! decay time for dripoff of liquid intercepted by upper canopy stems (sec) 
    REAL(KIND=r8) , INTENT(IN   ) :: tblows              ! decay time for blowoff of snow intercepted by upper canopy stems (sec)
    REAL(KIND=r8) , INTENT(INOUT) :: t34      (npoi)     ! air temperature at z34 (K)
    REAL(KIND=r8) , INTENT(IN   ) :: tdripl              ! decay time for dripoff of liquid intercepted by lower canopy leaves & stem (sec)
    REAL(KIND=r8) , INTENT(IN   ) :: tblowl              ! decay time for blowoff of snow intercepted by lower canopy leaves & stems (sec)
    REAL(KIND=r8) , INTENT(IN   ) :: ztop     (npoi,2)   ! height of plant top above ground (m)
    REAL(KIND=r8) , INTENT(IN   ) :: alaiml              ! lower canopy leaf & stem maximum area (2 sided) for normalization of drag coefficient (m2 m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: zbot     (npoi,2)   ! height of lowest branches above ground (m)
    REAL(KIND=r8) , INTENT(IN   ) :: alaimu              ! upper canopy leaf & stem area (2 sided) for normalization of drag coefficient (m2 m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: froot    (npoi,nsoilay,2)! fraction of root in soil layer 
    REAL(KIND=r8) , INTENT(INOUT) :: q34      (npoi)     ! specific humidity of air at z34
    REAL(KIND=r8) , INTENT(INOUT) :: q12      (npoi)     ! specific humidity of air at z12
    REAL(KIND=r8) , INTENT(INOUT) :: su       (npoi)     ! air-vegetation transfer coefficients (*rhoa) for upper canopy leaves (m s-1 * kg m-3) (A39a Pollard & Thompson 1995)
    REAL(KIND=r8) , INTENT(IN   ) :: cleaf               ! empirical constant in upper canopy leaf-air aerodynamic transfer coefficient (m s-0.5) (A39a Pollard & Thompson 95)
    REAL(KIND=r8) , INTENT(IN   ) :: dleaf    (2)        ! typical linear leaf dimension in aerodynamic transfer coefficient (m)
    REAL(KIND=r8) , INTENT(INOUT) :: ss       (npoi)     ! air-vegetation transfer coefficients (*rhoa) for upper canopy stems (m s-1 * kg m-3) (A39a Pollard & Thompson 1995)
    REAL(KIND=r8) , INTENT(IN   ) :: cstem               ! empirical constant in upper canopy stem-air aerodynamic transfer coefficient (m s-0.5) (A39a Pollard & Thompson 95)
    REAL(KIND=r8) , INTENT(IN   ) :: dstem    (2)        ! typical linear stem dimension in aerodynamic transfer coefficient (m)
    REAL(KIND=r8) , INTENT(INOUT) :: sl       (npoi)     ! air-vegetation transfer coefficients (*rhoa) for lower canopy leaves & stems (m s-1*kg m-3) (A39a Pollard & Thompson 1995)
    REAL(KIND=r8) , INTENT(IN   ) :: cgrass              ! empirical constant in lower canopy-air aerodynamic transfer coefficient (m s-0.5) (A39a Pollard & Thompson 95)
    REAL(KIND=r8) , INTENT(INOUT) :: ciub     (npoi)     ! intercellular co2 concentration - broadleaf (mol_co2/mol_air)
    REAL(KIND=r8) , INTENT(INOUT) :: ciuc     (npoi)     ! intercellular co2 concentration - conifer   (mol_co2/mol_air)
    REAL(KIND=r8) , INTENT(IN   ) :: exist    (npoi,npft)! probability of existence of each plant functional type in a gridcell
    REAL(KIND=r8) , INTENT(INOUT) :: csub     (npoi)     ! leaf boundary layer co2 concentration - broadleaf (mol_co2/mol_air)
    REAL(KIND=r8) , INTENT(INOUT) :: gsub     (npoi)     ! upper canopy stomatal conductance - broadleaf  (mol_co2 m-2 s-1)
    REAL(KIND=r8) , INTENT(INOUT) :: csuc     (npoi)     ! leaf boundary layer co2 concentration - conifer   (mol_co2/mol_air)
    REAL(KIND=r8) , INTENT(INOUT) :: gsuc     (npoi)     ! upper canopy stomatal conductance - conifer    (mol_co2 m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: agcub    (npoi)     ! canopy average gross photosynthesis rate - broadleaf  (mol_co2 m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: agcuc    (npoi)     ! canopy average gross photosynthesis rate - conifer    (mol_co2 m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: ancub    (npoi)     ! canopy average net photosynthesis rate - broadleaf    (mol_co2 m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: ancuc    (npoi)     ! canopy average net photosynthesis rate - conifer      (mol_co2 m-2 s-1)
    REAL(KIND=r8) , INTENT(INOUT) :: totcondub(npoi)     ! 
    REAL(KIND=r8) , INTENT(INOUT) :: totconduc(npoi)     !
    REAL(KIND=r8) , INTENT(INOUT) :: cils     (npoi)     ! intercellular co2 concentration - shrubs    (mol_co2/mol_air)
    REAL(KIND=r8) , INTENT(INOUT) :: cil3     (npoi)     ! intercellular co2 concentration - c3 plants (mol_co2/mol_air)
    REAL(KIND=r8) , INTENT(INOUT) :: cil4     (npoi)     ! intercellular co2 concentration - c4 plants (mol_co2/mol_air)
    REAL(KIND=r8) , INTENT(INOUT) :: csls     (npoi)     ! leaf boundary layer co2 concentration - shrubs(mol_co2/mol_air)
    REAL(KIND=r8) , INTENT(INOUT) :: gsls     (npoi)     ! lower canopy stomatal conductance - shrubs     (mol_co2 m-2 s-1)
    REAL(KIND=r8) , INTENT(INOUT) :: csl3     (npoi)     ! leaf boundary layer co2 concentration - c3 plants (mol_co2/mol_air)
    REAL(KIND=r8) , INTENT(INOUT) :: gsl3     (npoi)     ! lower canopy stomatal conductance - c3 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8) , INTENT(INOUT) :: csl4     (npoi)     ! leaf boundary layer co2 concentration - c4 plants (mol_co2/mol_air)
    REAL(KIND=r8) , INTENT(INOUT) :: gsl4     (npoi)     ! lower canopy stomatal conductance - c4 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: agcls    (npoi)     ! canopy average gross photosynthesis rate - shrubs     (mol_co2 m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: agcl4    (npoi)     ! canopy average gross photosynthesis rate - c4 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: agcl3    (npoi)     ! canopy average gross photosynthesis rate - c3 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: ancls    (npoi)     ! canopy average net photosynthesis rate - shrubs       (mol_co2 m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: ancl4    (npoi)     ! canopy average net photosynthesis rate - c4 grasses   (mol_co2 m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: ancl3    (npoi)     ! canopy average net photosynthesis rate - c3 grasses   (mol_co2 m-2 s-1)
    REAL(KIND=r8) , INTENT(INOUT) :: totcondls(npoi)     ! 
    REAL(KIND=r8) , INTENT(INOUT) :: totcondl3(npoi)     !
    REAL(KIND=r8) , INTENT(INOUT) :: totcondl4(npoi)     !
    REAL(KIND=r8) , INTENT(IN   ) :: chu(1:nVegClass)    ! heat capacity of upper canopy leaves per unit leaf area (J kg-1 m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: chs(1:nVegClass)    ! heat capacity of upper canopy stems per unit stem area (J kg-1 m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: chl(1:nVegClass)    ! heat capacity of lower canopy leaves & stems per unit leaf/stem area (J kg-1 m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: frac     (npoi,npft)! fraction of canopy occupied by each plant functional type
    REAL(KIND=r8) , INTENT(INOUT) :: tlsub    (npoi)     ! temperature of lower canopy vegetation buried by snow (K)
    !      INCLUDE 'compft.h'
    REAL(KIND=r8) , INTENT(IN   ) :: vmax_pft(npft)  ! nominal vmax of top leaf at 15 C (mol-co2/m**2/s) [not used]
    REAL(KIND=r8) , INTENT(IN   ) :: tau15           ! co2/o2 specificity ratio at 15 degrees C (dimensionless)
    REAL(KIND=r8) , INTENT(IN   ) :: kc15            ! co2 kinetic parameter (mol/mol)
    REAL(KIND=r8) , INTENT(IN   ) :: ko15            ! o2 kinetic parameter (mol/mol) 
    REAL(KIND=r8) , INTENT(IN   ) :: cimax           ! maximum value for ci (needed for model stability)
    REAL(KIND=r8) , INTENT(IN   ) :: gammaub         ! leaf respiration coefficient
    REAL(KIND=r8) , INTENT(IN   ) :: alpha3          ! intrinsic quantum efficiency for C3 plants (dimensionless)
    REAL(KIND=r8) , INTENT(IN   ) :: theta3          ! photosynthesis coupling coefficient for C3 plants (dimensionless)
    REAL(KIND=r8) , INTENT(IN   ) :: beta3           ! photosynthesis coupling coefficient for C3 plants (dimensionless)
    REAL(KIND=r8) , INTENT(IN   ) :: coefmub         ! 'm' coefficient for stomatal conductance relationship
    REAL(KIND=r8) , INTENT(IN   ) :: coefbub         ! 'b' coefficient for stomatal conductance relationship
    REAL(KIND=r8) , INTENT(IN   ) :: gsubmin         ! absolute minimum stomatal conductance
    REAL(KIND=r8) , INTENT(IN   ) :: gammauc         ! leaf respiration coefficient
    REAL(KIND=r8) , INTENT(IN   ) :: coefmuc         ! 'm' coefficient for stomatal conductance relationship  
    REAL(KIND=r8) , INTENT(IN   ) :: coefbuc         ! 'b' coefficient for stomatal conductance relationship  
    REAL(KIND=r8) , INTENT(IN   ) :: gsucmin         ! absolute minimum stomatal conductance
    REAL(KIND=r8) , INTENT(IN   ) :: gammals         ! leaf respiration coefficient
    REAL(KIND=r8) , INTENT(IN   ) :: coefmls         ! 'm' coefficient for stomatal conductance relationship 
    REAL(KIND=r8) , INTENT(IN   ) :: coefbls         ! 'b' coefficient for stomatal conductance relationship 
    REAL(KIND=r8) , INTENT(IN   ) :: gslsmin         ! absolute minimum stomatal conductance
    REAL(KIND=r8) , INTENT(IN   ) :: gammal3         ! leaf respiration coefficient
    REAL(KIND=r8) , INTENT(IN   ) :: coefml3         ! 'm' coefficient for stomatal conductance relationship 
    REAL(KIND=r8) , INTENT(IN   ) :: coefbl3         ! 'b' coefficient for stomatal conductance relationship 
    REAL(KIND=r8) , INTENT(IN   ) :: gsl3min         ! absolute minimum stomatal conductance
    REAL(KIND=r8) , INTENT(IN   ) :: gammal4         ! leaf respiration coefficient
    REAL(KIND=r8) , INTENT(IN   ) :: alpha4          ! intrinsic quantum efficiency for C4 plants (dimensionless)
    REAL(KIND=r8) , INTENT(IN   ) :: theta4          ! photosynthesis coupling coefficient for C4 plants (dimensionless) 
    REAL(KIND=r8) , INTENT(IN   ) :: beta4           ! photosynthesis coupling coefficient for C4 plants (dimensionless)
    REAL(KIND=r8) , INTENT(IN   ) :: coefml4         ! 'm' coefficient for stomatal conductance relationship
    REAL(KIND=r8) , INTENT(IN   ) :: coefbl4         ! 'b' coefficient for stomatal conductance relationship
    REAL(KIND=r8) , INTENT(IN   ) :: gsl4min         ! absolute minimum stomatal conductance
    !      INCLUDE 'comsum.h'
    REAL(KIND=r8) , INTENT(INOUT) :: a10scalparamu(npoi)! 10-day average day-time scaling parameter - upper canopy (dimensionless)
    REAL(KIND=r8) , INTENT(INOUT) :: a10daylightu (npoi)! 10-day average day-time PAR - upper canopy (micro-Ein m-2 s-1)
    REAL(KIND=r8) , INTENT(INOUT) :: a10scalparaml(npoi)! 10-day average day-time scaling parameter - lower canopy (dimensionless)
    REAL(KIND=r8) , INTENT(INOUT) :: a10daylightl (npoi)! 10-day average day-time PAR - lower canopy (micro-Ein m-2 s-1)

    !      INCLUDE 'comhyd.h'
    REAL(KIND=r8) , INTENT(OUT  ) :: ginvap (npoi) ! total evaporation rate from all intercepted h2o (kg_h2o m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: gsuvap (npoi) ! total evaporation rate from surface (snow/soil) (kg_h2o m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: gtrans (npoi) ! total transpiration rate from all vegetation canopies (kg_h2o m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: gtransu(npoi) ! transpiration from upper canopy (kg_h2o m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: gtransl(npoi) ! transpiration from lower canopy (kg_h2o m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: grunof (npoi) ! surface runoff rate (kg_h2o m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: gdrain (npoi) ! drainage rate out of bottom of lowest soil layer (kg_h2o m-2 s-1)
    REAL(KIND=r8) , INTENT(INOUT) :: gadjust(npoi) ! h2o flux due to adjustments in subroutine wadjust (kg_h2o m-2 s-1)
    REAL(KIND=r8) , INTENT(IN   ) :: ux  (npoi)
    REAL(KIND=r8) , INTENT(IN   ) :: uy  (npoi)
    REAL(KIND=r8) , INTENT(OUT  ) :: taux(npoi)
    REAL(KIND=r8) , INTENT(OUT  ) :: tauy(npoi)
    REAL(KIND=r8) , INTENT(OUT  ) :: ts2 (npoi)
    REAL(KIND=r8) , INTENT(OUT  ) :: qs2 (npoi)
    REAL(KIND=r8) , INTENT(IN   ) :: vegtype0(npoi)      ! annual vegetation type - ibis classification
    REAL(KIND=r8) , INTENT(IN   ) :: stressfac(nVegClass)
    REAL(KIND=r8) , INTENT(IN   ) :: avmuir_factor(nVegClass,2)
    REAL(KIND=r8) , INTENT(OUT  ) :: bstar(npoi)
    REAL(KIND=r8)  :: qh(npoi)


    !
    ! Local variables
    !
    INTEGER :: ib     ! waveband number (1= visible, 2= near-IR)
    INTEGER :: i      ! loop indice
    !
    ! set physical soil quantities
    !
    CALL setsoi(npoi   , &! INTENT(IN   )
         nsoilay, &! INTENT(IN   )
         sand   , &! INTENT(IN   )
         clay   , &! INTENT(IN   )
         poros  , &! INTENT(IN   )
         wsoi   , &! INTENT(IN   )
         wisoi  , &! INTENT(IN   )
         consoi , &! INTENT(OUT  )
         zwpmax , &! INTENT(IN   )
         wpud   , &! INTENT(IN   )
         wipud  , &! INTENT(IN   )
         wpudmax, &! INTENT(IN   )
         qglif  , &! INTENT(OUT  )
         tsoi         , &! INTENT(IN   )
         hvasug , &! INTENT(OUt  )
         hvasui , &! INTENT(OUt  )
         tsno   , &! INTENT(IN   )
         ta     , &! INTENT(IN   )
         nsnolay, &! INTENT(IN   )
         hvap   , &! INTENT(IN   )
         cvap   , &! INTENT(IN   )
         ch2o   , &! INTENT(IN   )
         hsub   , &! INTENT(IN   )
         cice   , &! INTENT(IN   )
         tmelt  , &! INTENT(IN   )
         epsilon  )! INTENT(IN   )
    !
    ! calculate areal fractions wetted by intercepted h2o
    !
    CALL fwetcal(npoi    , &! INTENT(IN   )
         fwetu   , &! INTENT(OUT  )
         rliqu   , &! INTENT(OUT  )
         fwets   , &! INTENT(OUT  )
         rliqs   , &! INTENT(OUT  )
         fwetl   , &! INTENT(OUT  )
         rliql   , &! INTENT(OUT  )
         wliqu   , &! INTENT(IN   )
         wliqumax, &! INTENT(IN   ) ::
         wsnou   , &! INTENT(IN   ) ::
         wsnoumax, &! INTENT(IN   ) ::
         tu      , &! INTENT(IN   )
         wliqs   , &! INTENT(IN   )
         wliqsmax, &! INTENT(IN   )
         wsnos   , &! INTENT(IN   )
         wsnosmax, &! INTENT(IN   )
         ts      , &! INTENT(IN   )
         wliql   , &! INTENT(IN   )
         wliqlmax, &! INTENT(IN   )
         wsnol   , &! INTENT(IN   )
         wsnolmax, &! INTENT(IN   )
         tl      , &! INTENT(IN   )
         epsilon , &! INTENT(IN   )
         tmelt     )! INTENT(IN   )
    !
    ! set up for solar calculations
    !
    CALL solset(npoi     , &! INTENT(IN   )
         nsol     , &! INTENT(OUT  )
         nband    , &! INTENT(IN   )
         solu     , &! INTENT(OUT  )
         sols     , &! INTENT(OUT  )
         soll     , &! INTENT(OUT  )
         solg     , &! INTENT(OUT  )
         soli     , &! INTENT(OUT  )
         scalcoefl, &! INTENT(OUT  )
         scalcoefu, &! INTENT(OUT  )
         indsol   , &! INTENT(OUT  )
         topparu  , &! INTENT(OUT  )
         topparl  , &! INTENT(OUT  )
         asurd    , &! INTENT(OUT  )
         asuri    , &! INTENT(OUT  )
         coszen     )! INTENT(IN   )  
    !
    ! solar calculations for each waveband
    !
    DO  ib = 1, nband
       !
       ! solsur sets surface albedos for soil and snow
       ! solalb performs the albedo calculations
       ! solarf uses the unit-incident-flux results from solalb
       ! to obtain absorbed fluxes sol[u,s,l,g,i] and 
       ! incident pars sunp[u,l]
       !
       CALL solsur (ib       , &! INTENT(IN   )
            tmelt    , &! INTENT(IN   )
            nsol     , &! INTENT(IN   )
            albsod   , &! INTENT(OUt  )
            albsoi   , &! INTENT(OUt  )
            albsnd   , &! INTENT(OUt  )
            albsni   , &! INTENT(OUt  )
            indsol   , &! INTENT(IN   )
            wsoi     , &! INTENT(IN   )
            wisoi    , &! INTENT(IN   )
            albsav   , &! INTENT(IN   )
            albsan   , &! INTENT(IN   )
            tsno     , &! INTENT(IN   )
            coszen   , &! INTENT(IN   )
            npoi     , &! INTENT(IN   )
            nsoilay  , &! INTENT(IN   )
            nsnolay    )! INTENT(IN   )

       CALL solalb (ib       , &! INTENT(IN   )
            nVegClass , &! INTENT(IN   )
            vegtype0  , &! INTENT(IN   )! annual vegetation type - ibis classification
            avmuir_factor, &! INTENT(IN   )
            relod    , &! INTENT(OUT  )
            reloi    , &! INTENT(OUT  )
            indsol   , &! INTENT(IN   )
            reupd    , &! INTENT(OUT  )
            reupi    , &! INTENT(OUT  )
            albsnd   , &! INTENT(IN   )
            albsni   , &! INTENT(IN   )
            albsod   , &! INTENT(IN   )
            albsoi   , &! INTENT(IN   )
            fl       , &! INTENT(IN   )
            fu       , &! INTENT(IN   )
            fi       , &! INTENT(IN   )
            asurd    , &! INTENT(INOUT)! local
            asuri    , &! INTENT(INOUT)! local
            npoi     , &! INTENT(IN   )
            nband    , &! INTENT(IN   )
            nsol     , &! INTENT(IN   )
            ablod    , &! INTENT(OUT  )
            abloi    , &! INTENT(OUT  )
            flodd    , &! INTENT(OUT  )
            dummy    , &! INTENT(OUT  )
            flodi    , &! INTENT(OUT  )
            floii    , &! INTENT(OUT  )
            coszen   , &! INTENT(IN   )
            terml    , &! INTENT(OUT  )
            termu    , &! INTENT(OUT  )
            lai      , &! INTENT(IN   )
            sai      , &! INTENT(IN   )
            abupd    , &! INTENT(OUT  )
            abupi    , &! INTENT(OUT  )
            fupdd    , &! INTENT(OUT  )
            fupdi    , &! INTENT(OUT  )
            fupii    , &! INTENT(OUT  )
            fwetl    , &! INTENT(IN   )
            rliql    , &! INTENT(IN   )
            rliqu    , &! INTENT(IN   )
            rliqs    , &! INTENT(IN   )
            fwetu    , &! INTENT(IN   )
            fwets    , &! INTENT(IN   )
            rhoveg   , &! INTENT(IN   )
            tauveg   , &! INTENT(IN   )
            orieh    , &! INTENT(IN   )
            oriev    , &! INTENT(IN   )
            tl       , &! INTENT(IN   )
            ts       , &! INTENT(IN   )
            tu       , &! INTENT(IN   )
            pi       , &! INTENT(IN   )
            tmelt    , &! INTENT(IN   )
            epsilon    )! INTENT(IN   )

       CALL solarf (ib       , & ! INTENT(IN) 
            nsol     , & ! INTENT(IN) 
            solu     , & ! INTENT(INOUT) !global
            indsol   , & ! INTENT(IN) 
            abupd    , & ! INTENT(IN) 
            abupi    , & ! INTENT(IN) 
            sols     , & ! INTENT(INOUT) !global
            sol2d    , & ! INTENT(OUT  ) 
            fupdd    , & ! INTENT(IN) 
            sol2i    , & ! INTENT(OUT  ) 
            fupii    , & ! INTENT(IN) 
            fupdi    , & ! INTENT(IN) 
            sol3d    , & ! INTENT(OUT  ) 
            sol3i    , & ! INTENT(OUT  ) 
            soll     , & ! INTENT(INOUT) !global
            ablod    , & ! INTENT(IN) 
            abloi    , & ! INTENT(IN) 
            flodd    , & ! INTENT(IN) 
            flodi    , & ! INTENT(IN) 
            floii    , & ! INTENT(IN) 
            solg     , & ! INTENT(INOUT) !global
            albsod   , & ! INTENT(IN) 
            albsoi   , & ! INTENT(IN) 
            soli     , & ! INTENT(INOUT) !global
            albsnd   , & ! INTENT(IN) 
            albsni   , & ! INTENT(IN) 
            scalcoefu, & ! INTENT(OUT  ) 
            termu    , & ! INTENT(IN) 
            scalcoefl, & ! INTENT(OUT  ) 
            terml    , & ! INTENT(IN) 
            lai      , & ! INTENT(IN) 
            sai      , & ! INTENT(IN) 
            fu       , & ! INTENT(IN)
            fl       , & ! INTENT(IN)
            topparu  , & ! INTENT(OUT  ) 
            topparl  , & ! INTENT(OUT  ) 
            solad    , & ! INTENT(IN)
            solai    , & ! INTENT(IN)
            npoi     , & ! INTENT(IN) 
            nband    , & ! INTENT(IN) 
            epsilon    ) ! INTENT(IN) 
       !
    END DO
    !
    ! calculate ir fluxes
    !
    CALL irrad(npoi  , &! INTENT(IN) ::
         nsoilay, &! INTENT(IN   ) ::
         stef  , &! INTENT(IN) ::
         nVegClass , &! INTENT(IN   )
         vegtype0  , &! INTENT(IN   )! annual vegetation type - ibis classification
         avmuir_factor, &! INTENT(IN   )
         firb  , &! INTENT(OUT  ) ::
         firs  , &! INTENT(OUT  ) ::
         firu  , &! INTENT(OUT  ) ::
         firl  , &! INTENT(OUT  ) ::
         firg  , &! INTENT(OUT  ) ::
         firi  , &! INTENT(OUT  ) ::
         lai   , &! INTENT(IN   ) ::
         sai   , &! INTENT(IN   ) ::
         fu    , &! INTENT(IN   ) ::
         tu    , &! INTENT(IN   ) ::
         ts    , &! INTENT(IN   ) ::
         tl    , &! INTENT(IN   ) ::
         fl    , &! INTENT(IN   ) ::
         tg    , &! INTENT(IN   ) ::
         ti    , &! INTENT(IN   ) ::
         fi    , &! INTENT(IN   ) ::
         fira  , &! INTENT(IN   ) :: 
         poros , &! INTENT(IN   ) :: 
         wsoi    )! INTENT(IN   ) ::
    !
    ! step intercepted h2o
    !
    CALL cascade(npoi    , & ! INTENT(IN   )
         epsilon , & ! INTENT(IN   )
         dtime   , & ! INTENT(IN   )
         ch2o    , & ! INTENT(IN   )
         cice    , & ! INTENT(IN   )  
         tmelt   , & ! INTENT(IN   )
         hfus    , & ! INTENT(IN   )
         vzero   , & ! INTENT(IN   )
         snowg   , & ! INTENT(OUT  )
         tsnowg  , & ! INTENT(OUT  )
         tsnowl  , & ! INTENT(INOUT)
         pfluxl  , & ! INTENT(OUT  )
         raing   , & ! INTENT(OUT  )
         traing  , & ! INTENT(OUT  )
         trainl  , & ! INTENT(INOUT)
         snowl   , & ! INTENT(OUT  )
         tsnowu  , & ! INTENT(INOUT)
         pfluxu  , & ! INTENT(OUT  )
         rainu   , & ! INTENT(INOUT)
         trainu  , & ! INTENT(INOUT)
         snowu   , & ! INTENT(INOUT)
         pfluxs  , & ! INTENT(OUT  )
         rainl   , & ! INTENT(OUT  )
         wliqmin , & ! INTENT(INOUT)
         wliqumax, & ! INTENT(IN   )  
         wsnomin , & ! INTENT(INOUT)
         wsnoumax, & ! INTENT(IN   )  
         t12     , & ! INTENT(IN   )  
         lai     , & ! INTENT(IN   )  
         tu      , & ! INTENT(IN   ) 
         wliqu   , & ! INTENT(INOUT)
         wsnou   , & ! INTENT(INOUT)
         tdripu  , & ! INTENT(IN   ) 
         tblowu  , & ! INTENT(IN   ) 
         sai     , & ! INTENT(IN   ) 
         ts           , & ! INTENT(IN   ) 
         wliqs   , & ! INTENT(INOUT) :: 
         wsnos   , & ! INTENT(INOUT) :: 
         tdrips  , & ! INTENT(IN   ) 
         tblows  , & ! INTENT(IN   ) 
         wliqsmax, & ! INTENT(IN   ) 
         wsnosmax, & ! INTENT(IN   ) 
         fu      , & ! INTENT(IN   ) 
         t34     , & ! INTENT(IN   ) 
         tl      , & ! INTENT(IN   ) 
         wliql   , & ! INTENT(INOUT) :: 
         wsnol   , & ! INTENT(INOUT) :: 
         tdripl  , & ! INTENT(IN   ) 
         tblowl  , & ! INTENT(IN   ) 
         wliqlmax, & ! INTENT(IN   ) 
         wsnolmax, & ! INTENT(IN   ) 
         fl           , & ! INTENT(IN   ) 
         raina   , & ! INTENT(IN   ) 
         ta      , & ! INTENT(IN   ) 
         qa      , & ! INTENT(IN   ) 
         psurf   , & ! INTENT(IN   )
         snowa     ) ! INTENT(IN   )
    !
    ! re-calculate wetted fractions, changed by cascade
    !
    CALL fwetcal(npoi    , &! INTENT(IN   )
         fwetu   , &! INTENT(OUT  )
         rliqu   , &! INTENT(OUT  )
         fwets   , &! INTENT(OUT  )
         rliqs   , &! INTENT(OUT  )
         fwetl   , &! INTENT(OUT  )
         rliql   , &! INTENT(OUT  )
         wliqu   , &! INTENT(IN   )
         wliqumax, &! INTENT(IN   ) ::
         wsnou   , &! INTENT(IN   ) ::
         wsnoumax, &! INTENT(IN   ) ::
         tu      , &! INTENT(IN   )
         wliqs   , &! INTENT(IN   )
         wliqsmax, &! INTENT(IN   )
         wsnos   , &! INTENT(IN   )
         wsnosmax, &! INTENT(IN   )
         ts      , &! INTENT(IN   )
         wliql   , &! INTENT(IN   )
         wliqlmax, &! INTENT(IN   )
         wsnol   , &! INTENT(IN   )
         wsnolmax, &! INTENT(IN   )
         tl      , &! INTENT(IN   )
         epsilon , &! INTENT(IN   )
         tmelt     )! INTENT(IN   )
    !
    ! step vegetation canopy temperatures implicitly
    ! and calculate sensible heat and moisture fluxes
    !
    CALL canopy(bps   , &! INTENT(INOUT) !local
         rhoa         , &! INTENT(INOUT) !local
         cp           , &! INTENT(INOUT) !local
         za           , &! INTENT(INOUT) !local
         bdl          , &! INTENT(INOUT) !local
         dil          , &! INTENT(INOUT) !local
         z3               , &! INTENT(INOUT) !local
         z4           , &! INTENT(INOUT) !local
         z34          , &! INTENT(INOUT) !local
         exphl        , &! INTENT(INOUT) !local
         expl         , &! INTENT(INOUT) !local
         displ        , &! INTENT(OUT  ) !local
         bdu               , &! INTENT(INOUT) !local
         diu          , &! INTENT(INOUT) !local
         z1           , &! INTENT(INOUT) !local
         z2           , &! INTENT(INOUT) !local
         z12          , &! INTENT(INOUT) !local
         exphu        , &! INTENT(INOUT) !local
         expu         , &! INTENT(INOUT) !local
         dispu        , &! INTENT(OUT  ) !local
         alogg        , &! INTENT(OUT  ) !local
         alogi        , &! INTENT(OUT  ) !local
         alogav       , &! INTENT(INOUT) !local
         alog4        , &! INTENT(INOUT) !local
         alog3        , &! INTENT(INOUT) !local
         alog2        , &! INTENT(OUT  ) !local
         alog1        , &! INTENT(INOUT) !local
         aloga        , &! INTENT(INOUT) !local
         u2           , &! INTENT(INOUT) !local
         alogu        , &! INTENT(INOUT) !local
         alogl        , &! INTENT(INOUT) !local
         richl        , &! INTENT(OUT  ) !local
         straml       , &! INTENT(OUT  ) !local
         strahl       , &! INTENT(OUT  ) !local
         richu        , &! INTENT(OUT  ) !local
         stramu       , &! INTENT(OUT  ) !local
         strahu       , &! INTENT(OUT  ) !local
         u1           , &! INTENT(OUT  ) !local
         u12          , &! INTENT(OUT  ) !local
         u3           , &! INTENT(OUT  ) !local
         u34          , &! INTENT(OUT  ) !local
         u4           , &! INTENT(OUT  ) !local
         cu               , &! INTENT(INOUT) !local
         cl           , &! INTENT(INOUT) !local
         sg           , &! INTENT(INOUT) !local
         si           , &! INTENT(INOUT) !local
         fwetl        , &! INTENT(IN) !global
         scalcoefl    , &! INTENT(IN) !global
         scalcoefu    , &! INTENT(IN) !global
         fwetu        , &! INTENT(IN) !global
         termu        , &! INTENT(IN) !global
         fwetux       , &! INTENT(OUT  ) !local
         fwetsx       , &! INTENT(OUT  ) !local
         fwets        , &! INTENT(IN) !global
         fwetlx       , &! INTENT(OUT  ) !local
         solu               , &! INTENT(IN) !global
         firu         , &! INTENT(IN) !global
         sols         , &! INTENT(IN) !global
         firs         , &! INTENT(IN) !global
         soll         , &! INTENT(IN) !global
         firl         , &! INTENT(IN) !global
         rliqu        , &! INTENT(IN) !global
         rliqs        , &! INTENT(IN) !global
         rliql        , &! INTENT(IN) !global
         pfluxu       , &! INTENT(IN) !global
         pfluxs       , &! INTENT(IN) !global
         pfluxl       , &! INTENT(IN) !global
         solg         , &! INTENT(IN) !global
         firg         , &! INTENT(INOUT) !global
         soli         , &! INTENT(IN) !global
         firi         , &! INTENT(INOUT) !global
         fsena        , &! INTENT(OUT  ) !local
         fseng        , &! INTENT(OUT  ) !local
         fseni        , &! INTENT(OUT  ) !local
         fsenu        , &! INTENT(OUT  ) !local
         fsens        , &! INTENT(OUT  ) !local
         fsenl        , &! INTENT(OUT  ) !local
         fvapa        , &! INTENT(OUT  ) !local
         fvaput       , &! INTENT(OUT  ) !local
         fvaps        , &! INTENT(OUT  ) !local
         fvaplw       , &! INTENT(OUT  ) !local
         fvaplt       , &! INTENT(OUT  ) !local
         fvapg        , &! INTENT(OUT  ) !local
         fvapi        , &! INTENT(OUT  ) !local
         firb         , &! INTENT(INOUT) !global
         terml        , &! INTENT(IN) !global
         fvapuw       , &! INTENT(OUT  ) !local
         ztop         , &! INTENT(IN) !global
         fl           , &! INTENT(IN) !global
         lai          , &! INTENT(IN) !global
         sai               , &! INTENT(IN) !global
         alaiml       , &! INTENT(IN) !global
         zbot         , &! INTENT(IN) !global
         fu           , &! INTENT(IN) !global
         alaimu       , &! INTENT(IN) !global
         froot        , &! INTENT(IN) !global
         t34          , &! INTENT(INOUT) !global
         t12          , &! INTENT(INOUT) !global
         q34               , &! INTENT(INOUT) !global
         q12          , &! INTENT(INOUT) !global
         su           , &! INTENT(INOUT) !local
         cleaf        , &! INTENT(IN) !global
         dleaf        , &! INTENT(IN) !global
         ss           , &! INTENT(INOUT) !local
         cstem        , &! INTENT(IN) !global
         dstem        , &! INTENT(IN) !global
         sl           , &! INTENT(INOUT) !local
         cgrass       , &! INTENT(IN) !global
         tu           , &! INTENT(INOUT) !global
         ciub         , &! INTENT(INOUT) !global
         ciuc         , &! INTENT(INOUT) !global
         exist        , &! INTENT(IN) !global
         topparu      , &! INTENT(IN) !global
         csub         , &! INTENT(INOUT) !global
         gsub         , &! INTENT(INOUT) !global
         csuc         , &! INTENT(INOUT) !global
         gsuc         , &! INTENT(INOUT) !global
         agcub        , &! INTENT(OUT  ) !local
         agcuc        , &! INTENT(OUT  ) !local
         ancub        , &! INTENT(OUT  ) !local
         ancuc        , &! INTENT(OUT  ) !local
         totcondub    , &! INTENT(INOUT) !local
         totconduc    , &! INTENT(INOUT) !local
         tl           , &! INTENT(INOUT) !global
         cils         , &! INTENT(INOUT) !global
         cil3         , &! INTENT(INOUT) !global
         cil4         , &! INTENT(INOUT) !global
         topparl      , &! INTENT(IN) !global
         csls         , &! INTENT(INOUT) !global
         gsls         , &! INTENT(INOUT) !global
         csl3         , &! INTENT(INOUT) !global
         gsl3         , &! INTENT(INOUT) !global
         csl4         , &! INTENT(INOUT) !global
         gsl4         , &! INTENT(INOUT) !global
         agcls        , &! INTENT(OUT  ) !local
         agcl4        , &! INTENT(OUT  ) !local
         agcl3        , &! INTENT(OUT  ) !local
         ancls        , &! INTENT(OUT  ) !local
         ancl4        , &! INTENT(OUT  ) !local
         ancl3        , &! INTENT(OUT  ) !local
         totcondls    , &! INTENT(INOUT) !local
         totcondl3    , &! INTENT(INOUT) !local
         totcondl4    , &! INTENT(INOUT) !local
         chu(1:nVegClass), &! INTENT(IN    ) !global
         wliqu        , &! INTENT(IN) !global
         wsnou        , &! INTENT(IN) !global
         chs(1:nVegClass), &! INTENT(IN    ) !global
         wliqs        , &! INTENT(IN) !global
         wsnos        , &! INTENT(IN) !global
         chl(1:nVegClass), &! INTENT(IN    ) !global
         wliql        , &! INTENT(IN) !global
         wsnol        , &! INTENT(IN) !global
         ts           , &! INTENT(INOUT) !global
         frac         , &! INTENT(IN) !global
         z0soi        , &! INTENT(IN) !global
         wsoi         , &! INTENT(IN) !global
         wisoi        , &! INTENT(IN) !global
         swilt        , &! INTENT(IN) !global
         sfield       , &! INTENT(IN) !global
         stressl      , &! INTENT(INOUT) !local
         stressu      , &! INTENT(INOUT) !local
         stresstl     , &! INTENT(INOUT) !local
         stresstu     , &! INTENT(INOUT) !local
         poros        , &! INTENT(IN) !global
         wpud         , &! INTENT(IN) !global
         wipud        , &! INTENT(IN) !global
         csoi         , &! INTENT(IN) !global
         rhosoi       , &! INTENT(IN) !global
         hsoi         , &! INTENT(IN) !global
         consoi       , &! INTENT(IN) !global
         tg           , &! INTENT(INOUT) !global
         ti           , &! INTENT(INOUT) !global
         wpudmax      , &! INTENT(IN) !global
         suction      , &! INTENT(IN) !global
         bex          , &! INTENT(IN) !global
         hvasug       , &! INTENT(IN) !global
         tsoi         , &! INTENT(IN) !global
         hvasui       , &! INTENT(IN) !global
         upsoiu       , &! INTENT(OUT  ) !local
         upsoil       , &! INTENT(OUT  ) !local
         fi           , &! INTENT(IN) !global
         z0sno        , &! INTENT(IN) !global
         consno       , &! INTENT(IN) !global
         hsno         , &! INTENT(IN) !global
         hsnotop      , &! INTENT(IN) !global
         tsno         , &! INTENT(IN) !global
         psurf        , &! INTENT(IN) !global
         ta               , &! INTENT(IN) !global
         qa           , &! INTENT(IN) !global
         ua           , &! INTENT(IN) !global
         o2conc       , &! INTENT(IN) !global
         co2conc      , &! INTENT(IN) !global
         npoi         , &! INTENT(IN) !global
         nsoilay      , &! INTENT(IN) !global
         nsnolay      , &! INTENT(IN) !global
         npft         , &! INTENT(IN) !global
         vonk         , &! INTENT(IN) !global
         epsilon      , &! INTENT(IN) !global
         hvap         , &! INTENT(IN) !global
         ch2o         , &! INTENT(IN) !global
         hsub         , &! INTENT(IN) !global
         cice         , &! INTENT(IN) !global
         rhow         , &! INTENT(IN) !global
         stef         , &! INTENT(IN) !global
         tmelt        , &! INTENT(IN) !global
         hfus         , &! INTENT(IN) !global
         cappa        , &! INTENT(IN) !global
         rair         , &! INTENT(IN) !global
         rvap         , &! INTENT(IN) !global
         cair         , &! INTENT(IN) !global
         cvap         , &! INTENT(IN) !global
         grav         , &! INTENT(IN) !global
         dtime        , &! INTENT(IN) !global
         vmax_pft     , &! INTENT(IN) !global
         tau15        , &! INTENT(IN) !global
         kc15         , &! INTENT(IN) !global
         ko15         , &! INTENT(IN) !global
         cimax        , &! INTENT(IN) !global
         gammaub      , &! INTENT(IN) !global
         alpha3       , &! INTENT(IN) !global
         theta3       , &! INTENT(IN) !global
         beta3        , &! INTENT(IN) !global
         coefmub      , &! INTENT(IN) !global
         coefbub      , &! INTENT(IN) !global
         gsubmin      , &! INTENT(IN) !global
         gammauc      , &! INTENT(IN) !global
         coefmuc      , &! INTENT(IN) !global
         coefbuc      , &! INTENT(IN) !global
         gsucmin      , &! INTENT(IN) !global
         gammals      , &! INTENT(IN) !global
         coefmls      , &! INTENT(IN) !global
         coefbls      , &! INTENT(IN) !global
         gslsmin      , &! INTENT(IN) !global
         gammal3      , &! INTENT(IN) !global
         coefml3      , &! INTENT(IN) !global
         coefbl3      , &! INTENT(IN) !global
         gsl3min      , &! INTENT(IN) !global
         gammal4      , &! INTENT(IN) !global
         alpha4       , &! INTENT(IN) !global
         theta4       , &! INTENT(IN) !global
         beta4        , &! INTENT(IN) !global
         coefml4      , &! INTENT(IN) !global
         coefbl4      , &! INTENT(IN) !global
         gsl4min      , &! INTENT(IN) !global
         a10scalparamu, &! INTENT(INOUT) !global
         a10daylightu , &! INTENT(INOUT) !global
         a10scalparaml, &! INTENT(INOUT) !global
         a10daylightl , &! INTENT(INOUT) !global
         ginvap       , &! INTENT(OUT  ) !local
         gsuvap       , &! INTENT(OUT  ) !local
         gtrans       , &! INTENT(OUT  ) !local
         gtransu      , &! INTENT(OUT  ) !local
         gtransl      , &! INTENT(OUT  ) !local
         ux               , &! INTENT(IN) !global
         uy               , &! INTENT(IN) !global
         taux         , &! INTENT(OUT  ) !local
         tauy         , &! INTENT(OUT  ) !local
         ts2               , &! INTENT(OUT  ) !local
         qs2          ,& ! INTENT(OUT  ) !local
         vegtype0     , &! INTENT(in  ) !local
         stressfac    , &! INTENT(in  ) !local
         nVegClass  )

    !
    ! step intercepted h2o due to evaporation
    !
    CALL cascad2(rliqu , &! INTENT(IN   )
         fvapuw, &! INTENT(IN   )
         fvapa , &! INTENT(INOUT)
         fsena , &! INTENT(INOUT)
         rliqs , &! INTENT(IN   )
         fvaps , &! INTENT(IN   )
         rliql , &! INTENT(IN   )
         fvaplw, &! INTENT(IN   )
         ta    , &! INTENT(IN   )
         fu    , &! INTENT(IN   )
         lai   , &! INTENT(IN   )
         tu    , &! INTENT(INOUT)
         wliqu , &! INTENT(INOUT)
         wsnou , &! INTENT(INOUT)
         chu(1:nVegClass), &! INTENT(IN   )
         sai   , &! INTENT(IN   )
         ts    , &! INTENT(INOUT)
         wliqs , &! INTENT(INOUT)
         wsnos , &! INTENT(INOUT)
         chs(1:nVegClass), &! INTENT(IN   )
         fl    , &! INTENT(IN   )
         tl    , &! INTENT(INOUT)
         wliql , &! INTENT(INOUT)
         wsnol , &! INTENT(INOUT)
         chl(1:nVegClass), &! INTENT(IN   )
         fi    , &! INTENT(IN   )
         npoi  , &! INTENT(IN   )
         hvap  , &! INTENT(IN   )
         cvap  , &! INTENT(IN   )
         ch2o  , &! INTENT(IN   )
         hsub  , &! INTENT(IN   )
         cice  , &! INTENT(IN   )
         dtime , &! INTENT(IN   )
         hfus  , &! INTENT(IN   )
         vegtype0, &! INTENT(IN   )
         tmelt   ,&! INTENT(IN   )
         nVegClass)! INTENT(IN   )
    !
    ! arbitrarily set veg temps & intercepted h2o for no-veg locations
    !
    CALL noveg(lai  , &! INTENT(IN   )
         fu   , &! INTENT(IN   )
         tu   , &! INTENT(INOUT)
         wliqu, &! INTENT(INOUT)
         sai  , &! INTENT(IN   )
         ts   , &! INTENT(INOUT)
         wliqs, &! INTENT(INOUT)
         wsnos, &! INTENT(INOUT)
         fl   , &! INTENT(IN   )
         tl   , &! INTENT(INOUT)
         wliql, &! INTENT(INOUT)
         wsnol, &! INTENT(INOUT)
         wsnou, &! INTENT(INOUT)
         tg   , &! INTENT(IN   )
         ti   , &! INTENT(IN   )
         fi   , &! INTENT(IN   )
         npoi   )! INTENT(IN   )
    !
    ! set net surface heat fluxes for soil and snow models
    !
    DO i = 1, npoi
       !
       qh(i)=0.622e0_r8*EXP(21.65605e0_r8 -5418.0e0_r8 /ta(i))/(psurf(i)/100.0_r8)
       BSTAR(i) = (grav/(rhoa(i)*SQRT(cu(i)*MAX(ua(i),1.e-30_r8)/rhoa(i)))) * (sg(i)*(ts2(i)-ta(i))/ta(i)) 
       !(CT*(TH-TA-(MAPL_GRAV/MAPL_CP)*DZ)/TA + MAPL_VIREPS*CQ*(QH-QA))

       !bstar(i)=(sg(i))*grav*(sg(i)*(ts2(i)-ta(i))/ta(i))! &
       !+mapl_vireps*(cu(i)/(rhoa(i)*ua(i)))*(qh(i)-qa(i)))

       heatg(i) = solg(i) + firg(i) - fseng(i) - &
            hvasug(i)*fvapg(i)
       !
       heati(i) = soli(i) + firi(i) - fseni(i) - &
            hvasui(i)*fvapi(i)
       !
    END DO
    !
    ! step snow model
    !
    CALL snow(rainl   , &! INTENT(IN   )
         trainl  , &! INTENT(IN   )
         snowl   , &! INTENT(IN   )
         tsnowl  , &! INTENT(IN   )
         fvapi   , &! INTENT(IN   )
         snowg   , &! INTENT(IN   )
         tsnowg  , &! INTENT(IN   )
         solg    , &! INTENT(INOUT)
         fvapg   , &! INTENT(INOUT)
         raing   , &! INTENT(INOUT)
         traing  , &! INTENT(INOUT)
         fl      , &! INTENT(IN   )
         ztop    , &! INTENT(IN   )
         lai     , &! INTENT(IN   )
         sai     , &! INTENT(IN   )
         tlsub   , &! INTENT(INOUT)
         chl(1:nVegClass), &! INTENT(IN   )
         tl      , &! INTENT(IN   )
         wliql   , &! INTENT(INOUT)
         wsnol   , &! INTENT(INOUT)
         hsnomin , &! INTENT(IN   )
         hsnotop , &! INTENT(IN   )
         fi      , &! INTENT(INOUT)
         rhos    , &! INTENT(IN   )
         tsno    , &! INTENT(INOUT)
         hsno    , &! INTENT(INOUT)
         consno  , &! INTENT(IN   )
         fimax   , &! INTENT(IN   )
         fimin   , &! INTENT(IN   )
         heati   , &! INTENT(IN   )
         hsoi    , &! INTENT(IN   )
         consoi  , &! INTENT(IN   )
         tsoi    , &! INTENT(IN   )
         heatg   , &! INTENT(INOUT)
         npoi    , &! INTENT(IN   )
         nsoilay , &! INTENT(IN   )
         nsnolay , &! INTENT(IN   )
         dtime   , &! INTENT(IN   )
         cice    , &! INTENT(IN   )
         epsilon , &! INTENT(IN   )
         ch2o    , &! INTENT(IN   )
         tmelt   , &! INTENT(IN   )
         hfus    , &! INTENT(IN   )
         vegtype0, &! INTENT(IN   )
         vzero   , &
         nVegClass ) ! INTENT(IN   )
    !
    ! step soil model
    !
    CALL soilctl(raing   , &!  INTENT(IN   )
         fvapg   , &!  INTENT(IN   )
         traing  , &!  INTENT(IN   )
         tu      , &!  INTENT(IN   )
         tl      , &!  INTENT(IN   )
         wpud    , &!  INTENT(INOUT)
         wipud   , &!  INTENT(INOUT)
         wpudmax , &!  INTENT(IN   )
         tsoi    , &!  INTENT(INOUT)
         qglif   , &!  INTENT(IN   )
         wisoi   , &!  INTENT(INOUT)
         poros   , &!  INTENT(IN   )
         hsoi    , &!  INTENT(IN   )
         wsoi    , &!  INTENT(INOUT)
         hydraul , &!  INTENT(IN   )
         heatg   , &!  INTENT(IN   )
         porosflo, &!  INTENT(INOUT)
         upsoiu  , &!  INTENT(IN   )
         upsoil  , &!  INTENT(IN   )
         csoi    , &!  INTENT(IN   )
         rhosoi  , &!  INTENT(IN   )
         suction , &!  INTENT(IN   )
         bex     , &!  INTENT(IN   )
         ibex    , &!  INTENT(IN   )
         bperm   , &!  INTENT(IN   )
         consoi  , &!  INTENT(IN   )
         hflo    , &!  INTENT(INOUT)
         grunof  , &!  INTENT(OUT  )
         gadjust , &!  INTENT(INOUT)
         gdrain  , &!  INTENT(OUT  )
         npoi    , &!  INTENT(IN   )
         nsoilay , &!  INTENT(IN   )
         epsilon , &!  INTENT(IN   )
         dtime   , &!  INTENT(IN   )
         ch2o    , &!  INTENT(IN   )
         cice    , &!  INTENT(IN   )
         tmelt   , &!  INTENT(IN   )
         rhow    , &!  INTENT(IN   )
         hfus    , &!  INTENT(IN   )
         stressl , &! soil moisture stress factor for the lower canopy (dimensionless)
         stressu , &! soil moisture stress factor for the upper canopy (dimensionless)
         stresstl, &! sum of stressl over all 6 soil layers (dimensionless)
         stresstu  )! sum of stressu over all 6 soil layers (dimensionless)

    !
    ! return to main program
    !
    RETURN
  END  SUBROUTINE lsxmain
  !
  !  ####   #    #   ####   #    #
  ! #       ##   #  #    #  #    #
  !  ####   # #  #  #    #  #    #
  !      #  #  # #  #    #  # ## #
  ! #    #  #   ##  #    #  ##  ##
  !  ####   #    #   ####   #    #
  !


  SUBROUTINE snow (rainl  , & ! INTENT(IN   )
       trainl , & ! INTENT(IN   )
       snowl  , & ! INTENT(IN   )
       tsnowl , & ! INTENT(IN   )
       fvapi  , & ! INTENT(IN   )
       snowg  , & ! INTENT(IN   )
       tsnowg , & ! INTENT(IN   )
       solg   , & ! INTENT(INOUT)
       fvapg  , & ! INTENT(INOUT)
       raing  , & ! INTENT(INOUT)
       traing , & ! INTENT(INOUT)
       fl     , & ! INTENT(IN   )
       ztop   , & ! INTENT(IN   )
       lai    , & ! INTENT(IN   )
       sai    , & ! INTENT(IN   )
       tlsub  , & ! INTENT(INOUT)
       chl    , & ! INTENT(IN   )
       tl     , & ! INTENT(IN   )
       wliql  , & ! INTENT(INOUT)
       wsnol  , & ! INTENT(INOUT)
       hsnomin, & ! INTENT(IN   )
       hsnotop, & ! INTENT(IN   )
       fi     , & ! INTENT(INOUT)
       rhos   , & ! INTENT(IN   )
       tsno   , & ! INTENT(INOUT)
       hsno   , & ! INTENT(INOUT)
       consno , & ! INTENT(IN   )
       fimax  , & ! INTENT(IN   )
       fimin  , & ! INTENT(IN   )
       heati  , & ! INTENT(IN   )
       hsoi   , & ! INTENT(IN   )
       consoi , & ! INTENT(IN   )
       tsoi   , & ! INTENT(IN   )
       heatg  , & ! INTENT(INOUT)
       npoi   , & ! INTENT(IN   )
       nsoilay, & ! INTENT(IN   )
       nsnolay, & ! INTENT(IN   )
       dtime  , & ! INTENT(IN   )
       cice   , & ! INTENT(IN   )
       epsilon, & ! INTENT(IN   )
       ch2o   , & ! INTENT(IN   )
       tmelt  , & ! INTENT(IN   )
       hfus   , & ! INTENT(IN   )
       vegtype0, &! INTENT(IN   )
       vzero     ,&
       nVegClass ) ! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    ! steps snow model through one timestep
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN   ) :: nVegClass
    INTEGER, INTENT(IN   ) :: npoi    ! total number of land points
    INTEGER, INTENT(IN   ) :: nsoilay ! number of soil layers
    INTEGER, INTENT(IN   ) :: nsnolay ! number of snow layers
    REAL(KIND=r8)   , INTENT(IN   ) :: dtime   ! model timestep (seconds)
    REAL(KIND=r8)   , INTENT(IN   ) :: cice    ! specific heat of ice (J deg-1 kg-1)
    REAL(KIND=r8)   , INTENT(IN   ) :: epsilon ! small quantity to avoid zero-divides and other
    ! truncation or machine-limit troubles with small
    ! values. should be slightly greater than o(1)
    ! machine precision
    REAL(KIND=r8)   , INTENT(IN   ) :: ch2o    ! specific heat of liquid water (J deg-1 kg-1)
    REAL(KIND=r8)   , INTENT(IN   ) :: tmelt   ! freezing point of water (K)
    REAL(KIND=r8)   , INTENT(IN   ) :: hfus    ! latent heat of fusion of water (J kg-1)
    REAL(KIND=r8)   , INTENT(IN   ) :: vzero  (npoi)        ! a REAL(KIND=r8) array of zeros, of length npoi
    REAL(KIND=r8)    , INTENT(IN   ) :: heati  (npoi)              ! net heat flux into snow surface (W m-2)
    REAL(KIND=r8)    , INTENT(IN   ) :: hsoi   (npoi,nsoilay+1)   ! soil layer thickness (m)
    REAL(KIND=r8)    , INTENT(IN   ) :: consoi (npoi,nsoilay)! thermal conductivity of each soil layer (W m-1 K-1)
    REAL(KIND=r8)    , INTENT(IN   ):: tsoi   (npoi,nsoilay)! soil temperature for each layer (K)
    REAL(KIND=r8)    , INTENT(INOUT) :: heatg  (npoi)              ! net heat flux into soil surface (W m-2)
    REAL(KIND=r8)    , INTENT(IN   ) :: hsnomin              ! minimum total thickness of snow (m)
    REAL(KIND=r8)    , INTENT(IN   ) :: hsnotop              ! thickness of top snow layer (m)
    REAL(KIND=r8)    , INTENT(INOUT) :: fi     (npoi)        ! fractional snow cover
    REAL(KIND=r8)    , INTENT(IN   ) :: rhos                 ! density of snow (kg m-3)
    REAL(KIND=r8)    , INTENT(INOUT) :: tsno   (npoi,nsnolay)! temperature of snow layers (K)
    REAL(KIND=r8)    , INTENT(INOUT) :: hsno   (npoi,nsnolay)! thickness of snow layers (m)
    REAL(KIND=r8)    , INTENT(IN   ) :: consno               ! thermal conductivity of snow (W m-1 K-1)
    REAL(KIND=r8)    , INTENT(IN   ) :: fimax                ! maximum fractional snow cover
    REAL(KIND=r8)    , INTENT(IN   ) :: fimin                ! minimum fractional snow cover
    REAL(KIND=r8)    , INTENT(IN   ) :: fl    (npoi)             ! fraction of snow-free area covered by lower  canopy
    REAL(KIND=r8)    , INTENT(IN   ) :: ztop  (npoi,2)      ! height of plant top above ground (m)
    REAL(KIND=r8)    , INTENT(IN   ) :: lai   (npoi,2)      ! canopy single-sided leaf area index (area leaf/area veg)
    REAL(KIND=r8)    , INTENT(IN   ) :: sai   (npoi,2)      ! current single-sided stem area index
    REAL(KIND=r8)    , INTENT(INOUT) :: tlsub (npoi)             ! temperature of lower canopy vegetation buried by snow (K)
    REAL(KIND=r8)    , INTENT(IN   ) :: chl (1:nVegClass)          ! heat capacity of lower canopy leaves & stems per unit leaf/stem area (J kg-1 m-2)
    REAL(KIND=r8)    , INTENT(IN   ) :: tl    (npoi)             ! temperature of lower canopy leaves & stems(K)
    REAL(KIND=r8)    , INTENT(INOUT) :: wliql (npoi)             ! intercepted liquid h2o on lower canopy leaf and stem area (kg m-2)
    REAL(KIND=r8)    , INTENT(INOUT ):: wsnol (npoi)             ! intercepted frozen h2o (snow) on lower canopy leaf & stem area (kg m-2)
    REAL(KIND=r8)    , INTENT(IN   ) :: rainl (npoi)        ! rainfall rate below upper canopy (kg m-2 s-1)
    REAL(KIND=r8)    , INTENT(IN   ) :: trainl(npoi)        ! rainfall temperature below upper canopy (K)
    REAL(KIND=r8)    , INTENT(IN   ) :: snowl (npoi)        ! snowfall rate below upper canopy (kg h2o m-2 s-1)
    REAL(KIND=r8)    , INTENT(IN   ) :: tsnowl(npoi)        ! snowfall temperature below upper canopy (K)
    REAL(KIND=r8)    , INTENT(IN   ) :: fvapi (npoi)        ! h2o vapor flux (evaporation) between snow & air at z34 (kg m-2 s-1 / fi )
    REAL(KIND=r8)    , INTENT(IN   ) :: snowg (npoi)        ! snowfall rate at soil level (kg h2o m-2 s-1)
    REAL(KIND=r8)    , INTENT(IN   ) :: tsnowg(npoi)        ! snowfall temperature at soil level (K) 
    REAL(KIND=r8)    , INTENT(INOUT) :: solg  (npoi)        ! solar flux (direct + diffuse) absorbed by unit snow-free soil (W m-2)
    REAL(KIND=r8)    , INTENT(INOUT) :: fvapg (npoi)        ! h2o vapor flux (evaporation) between soil & air at z34
    ! (kg m-2 s-1/bare ground fraction)
    REAL(KIND=r8)    , INTENT(INOUT) :: raing (npoi)        ! rainfall rate at soil level (kg m-2 s-1)
    REAL(KIND=r8)    , INTENT(INOUT) :: traing(npoi)        ! rainfall temperature at soil level (K)
    REAL(KIND=r8)    , INTENT(IN   ) :: vegtype0(npoi)   

    REAL(KIND=r8)    :: traing2(npoi)        ! rainfall temperature at soil level (K)

    !
    ! local variables:
    !
    INTEGER :: i                   ! loop indices
    INTEGER :: k                   ! loop indices
    INTEGER :: npn                 ! index indsno, npcounter for pts with snow
    INTEGER :: iveg

    REAL(KIND=r8)    :: rwork               ! working vaiable
    REAL(KIND=r8)    :: rwork2              ! working vaiable
    REAL(KIND=r8)    :: finew               ! storing variable for fi
    REAL(KIND=r8)    :: zhh                 ! 0.5*hsnomin
    REAL(KIND=r8)    :: zdh                 ! max height of snow above hsnomin (?)     

    INTEGER :: indsno(npoi)        ! index of points with snow in current 1d strip
    !
    REAL(KIND=r8)    :: hinit(nsnolay)      ! initial layer thicknesses when snow first forms
    REAL(KIND=r8)    :: hsnoruf(npoi)       ! heigth of snow forced to cover lower canopy (?)
    REAL(KIND=r8)    :: fiold(npoi)         ! old fi at start of this timestep
    REAL(KIND=r8)    :: fhtop(npoi)         ! heat flux into upper snow surface
    REAL(KIND=r8)    :: sflo(npoi,nsnolay+2)! heat flux across snow and buried-lower-veg layer bdries
    REAL(KIND=r8)    :: zmelt(npoi)         ! liquid mass flux increments to soil, at temperature 
    ! tmelt, due to processes occuring during this step
    REAL(KIND=r8)    :: zheat(npoi)         ! heat flux to soil, due to processes occuring this step
    REAL(KIND=r8)    :: dfi(npoi)           ! change in fi
    REAL(KIND=r8)    :: xl(npoi)            ! lower veg density
    REAL(KIND=r8)    :: xh(npoi)            ! temporary arrays
    REAL(KIND=r8)    :: xm(npoi)            ! "
    REAL(KIND=r8)    :: ht(npoi)            ! "
    REAL(KIND=r8)    :: x1(npoi)            ! "
    REAL(KIND=r8)    :: x2(npoi)            ! "
    REAL(KIND=r8)    :: x3(npoi)            ! "
    !
    DO i = 1, npoi
       hsnoruf(i) =  MIN (0.700_r8, MAX (hsnomin+.050_r8, fl(i)*ztop(i,1)))
       xl(i) = fl(i) * 2.00_r8 * (lai(i,1) + sai(i,1))
       x1(i) = tlsub(i)
    END DO
    !
    hinit(1) = hsnotop
    !
    DO k = 2, nsnolay
       hinit(k) = (hsnomin - hsnotop) / (nsnolay-1)
    END DO
    !
    DO i = 1, npoi
       fiold(i) = fi(i)
    END DO
    !
    ! zero out arrays
    !

    DO k = 1, nsnolay+2
       DO i = 1, npoi
          sflo(i,k)  = 0.0_r8
       END DO
    END DO

    !      CALL const (sflo            , &  !INTENT(OUT  )
    !                  npoi*(nsnolay+2), &  !INTENT(IN   )
    !          0.0_r8               )  !INTENT(IN   )
    DO i = 1, npoi
       zmelt(i)  = 0.0_r8     
    END DO

    !      CALL const (zmelt           , &  !INTENT(OUT  )
    !                  npoi            , &  !INTENT(IN   )
    !          0.0_r8               )  !INTENT(IN   )

    DO i = 1, npoi
       zheat(i)  = 0.0_r8     
    END DO

    !      CALL const (zheat           , &  !INTENT(OUT  )
    !                  npoi            , &  !INTENT(IN   )
    !          0.0_r8               )  !INTENT(IN   )
    !
    ! set up index indsno, npn for pts with snow - indsno is used
    ! only by vadapt - elsewhere below, just test on npn > 0
    !
    npn = 0                                        
    !
    DO i = 1, npoi 
       IF (fi(i).GT.0.0_r8) THEN
          npn = npn + 1
          indsno(npn) = i
       END IF
    END DO
    !
    ! set surface heat flux fhtop and increment top layer thickness
    ! due to rainfall, snowfall and sublimation on existing snow
    !
    IF (npn.GT.0) THEN
       !
       rwork = dtime / rhos
       !
       DO i = 1, npoi
          !
          fhtop(i) = heati(i) +   &
               rainl(i) * (ch2o * (trainl(i) - tmelt)     + hfus + &
               cice * (tmelt     - tsno(i,1)))       + &
               snowl(i) *  cice * (tsnowl(i) - tsno(i,1))
          !
          IF (fi(i).GT.0.0_r8) hsno(i,1) = hsno(i,1) +   &
               (rainl(i) + snowl(i) - fvapi(i)) * rwork
          !
       END DO
       !
    END IF
    !
    ! step temperatures due to heat conduction, including buried
    ! lower-veg temperature tlsub
    !
    IF (npn.GT.0) THEN
       !
       DO i=1,npoi
          x1(i)=tlsub(i)
       END DO
       CALL snowheat (tlsub    , &! INTENT(INOUT)
            fhtop    , &! INTENT(IN   )
            sflo     , &! INTENT(out  )
            xl       , &! INTENT(IN   )
            chl(1:nVegClass), &! INTENT(IN   )
            hsoi     , &! INTENT(IN   )
            consoi   , &! INTENT(IN   )
            tsoi     , &! INTENT(IN   )
            tsno     , &! INTENT(INOUT)
            consno   , &! INTENT(IN   )
            hsno     , &! INTENT(IN   )
            rhos     , &! INTENT(IN   )
            fi       , &! INTENT(IN   )
            vegtype0 , &! INTENT(IN   )
            npoi     , &! INTENT(IN   )
            nsnolay  , &! INTENT(IN   )
            nsoilay  , &! INTENT(IN   )
            dtime    , &! INTENT(IN   )
            cice     , &! INTENT(IN   )
            nVegClass ) ! INTENT(IN   )
       !
    END IF
    !
    ! put snowfall from 1-fi snow-free area onto side of existing
    ! snow, or create new snow if current fi = 0. also reset index.
    ! (assumes total depth of newly created snow = hsnomin.)
    ! (fi will not become gt 1 here if one timestep's snowfall
    ! <= hsnomin, but protect against this anyway.)
    !
    ! if no adjacent snowfall or fi = 1, dfi = 0, so no effect
    !

    DO i=1,npoi
       ht(i) = 0.0_r8
    END DO

    !      CALL const (ht    , &!INTENT(OUT  )
    !                  npoi  , &!INTENT(IN)
    !          0.0_r8     )!INTENT(IN)
    DO k=1,nsnolay
       DO i=1,npoi
          ht(i) = ht(i) + hsno(i,k)
       END DO
    END DO
    !
    DO i=1,npoi
       IF (ht(i).EQ.0.0_r8) ht(i) = hsnomin
    END DO
    !
    rwork = dtime / rhos
    DO i=1,npoi
       dfi(i) = (1.0_r8-fi(i))*rwork*snowg(i) / ht(i)
       dfi(i) = MIN (dfi(i), 1.0_r8-fi(i))
    END DO
    !
    DO  k=1,nsnolay
       DO i=1,npoi
          IF (fi(i)+dfi(i).GT.0.0_r8)    &
               tsno(i,k) = (tsno(i,k)*fi(i) + tsnowg(i)*dfi(i))   &
               / (fi(i)+dfi(i))
          !
          ! set initial thicknesses for newly created snow
          !
          IF (fi(i).EQ.0.0_r8 .AND. dfi(i).GT.0.0_r8) hsno(i,k) = hinit(k)
       END DO
    END DO
    !
    npn = 0
    DO i=1,npoi
       fi(i) = fi(i) + dfi(i)
       IF (fi(i).GT.0.0_r8) THEN 
          npn = npn + 1
          indsno(npn) = i
       END IF
    END DO
    !
    ! melt from any layer (due to implicit heat conduction, any
    ! layer can exceed tmelt, not just the top layer), and reduce
    ! thicknesses (even to zero, and give extra heat to soil)
    !
    ! ok to do it for non-snow points, for which xh = xm = 0
    !
    IF (npn.GT.0) THEN
       !
       rwork = 1.0_r8 / rhos
       DO k=1,nsnolay
          DO i=1,npoi
             xh(i) = rhos*hsno(i,k)*cice * MAX(tsno(i,k)-tmelt, 0.0_r8)
             xm(i) = MIN (rhos*hsno(i,k), xh(i)/hfus)
             hsno(i,k) = hsno(i,k) - xm(i)*rwork
             tsno(i,k) = MIN (tsno(i,k),tmelt)
             zmelt(i) = zmelt(i) + fi(i)*xm(i)
             zheat(i) = zheat(i) + fi(i)*(xh(i)-hfus*xm(i))
          END DO
       END DO
       !
       ! adjust fi and thicknesses for coverage-vs-volume relation
       ! ie, total thickness = hsnomin for fi < fimax, and fi <= fimax.
       ! (ok to do it for no-snow points, for which ht=fi=finew=0.)
       !          
       DO i=1,npoi
          ht(i) = 0.0_r8 
       END DO

       !        CALL const (ht     , & !INTENT(OUT  )
       !                    npoi   , & !INTENT(IN   )
       !            0.0_r8      ) !INTENT(IN   )
       DO k=1,nsnolay
          DO i=1,npoi
             ht(i) = ht(i) + hsno(i,k)
          END DO
       END DO
       !
       ! linear variation  for 0 < fi < 1
       !
       zhh = 0.50_r8*hsnomin
       DO i=1,npoi
          zdh = hsnoruf(i)-hsnomin
          finew = ( -zhh + SQRT(zhh**2 + zdh*fi(i)*ht(i)) ) / zdh

          finew = MAX (0.0_r8, MIN (fimax, finew))
          x1(i) =  fi(i) / MAX (finew, epsilon)
          fi(i) =  finew
       END DO
       !
       DO k=1,nsnolay
          DO i=1,npoi
             hsno(i,k) = hsno(i,k) * x1(i)
          END DO
       END DO
       !
    END IF
    !
    ! re-adapt snow thickness profile, so top thickness = hsnotop
    ! and other thicknesses are equal
    !
    ! adjust temperature to conserve sensible heat
    !
    CALL vadapt (hsno   , &! INTENT(INOUT)
         tsno   , &! INTENT(INOUT)
         hsnotop, &! INTENT(IN   )
         indsno , &! INTENT(IN   )
         npn    , &! INTENT(IN   )
         nsnolay, &! INTENT(IN   )
         npoi   , &! INTENT(IN   )
         nsnolay, &! INTENT(IN   )
         epsilon  )! INTENT(IN   )
    !
    ! if fi is below fimin, melt all snow and adjust soil fluxes
    !
    IF (npn.GT.0) THEN

       DO  i = 1, npoi
          x1(i) = fi(i)
       END DO

       !        CALL scopy (npoi   , & ! INTENT(IN   )
       !                    fi     , & ! INTENT(IN   )
       !            x1       ) ! INTENT(OUT  )

       DO k=1,nsnolay
          DO i=1,npoi
             IF (x1(i).LT.fimin) THEN
                xm(i) = x1(i) * rhos * hsno(i,k)
                zmelt(i) = zmelt(i) + xm(i)
                zheat(i) = zheat(i) - xm(i)*(cice*(tmelt-tsno(i,k))+hfus)
                hsno(i,k) = 0.0_r8
                tsno(i,k) = tmelt
                fi(i) = 0.0_r8
             END IF
          END DO
       END DO
    END IF
    !
    ! adjust buried lower veg for fi changes. if fi has increased,
    ! incorporate newly buried intercepted h2o into bottom-layer 
    ! snow, giving associated heat increment to soil, and mix the
    ! specific heat of newly buried veg (at tl) into tlsub. if fi
    ! has decreased, change temp of newly exhumed veg to tl, giving
    ! assoc heat increment to soil, and smear out intercepted h2o
    !
    IF (npn.GT.0) THEN
       DO i=1,npoi
          iveg=vegtype0(i)   
          dfi(i) = fi(i) - fiold(i)
          !
          IF (dfi(i).GT.0.0_r8) THEN
             !
             ! factor of xl*chl has been divided out of next line
             !
             tlsub(i)= (tlsub(i)*fiold(i) + tl(i)*dfi(i)) / fi(i)
             zheat(i) = zheat(i) + dfi(i)*xl(i)              &
                  * ( wliql(i) * (ch2o*(tl(i)-tmelt) + hfus     &
                  +cice*(tmelt-tsno(i,nsnolay))) &
                  + wsnol(i) *  cice*(tl(i)-tsno(i,nsnolay)) )
             !
             hsno(i,nsnolay) = hsno(i,nsnolay)   &
                  + dfi(i)*xl(i)*(wliql(i)+wsnol(i)) &
                  / (rhos*fi(i))
          END IF
          !
          IF (dfi(i).LT.0.0_r8) THEN
             zheat(i) = zheat(i) - dfi(i)*xl(i)*chl(iveg)*(tlsub(i)-tl(i))
             rwork = (1.0_r8-fiold(i)) / (1.0_r8-fi(i))
             wliql(i) = wliql(i) * rwork
             wsnol(i) = wsnol(i) * rwork
          END IF
          !
       END DO
    END IF
    !
    ! areally average fluxes to be used by soil model. (don't use
    ! index due to mix call, but only need at all if npn > 0)
    !
    IF (npn.GT.0) THEN
       !
       rwork = 1.0_r8 / dtime
       DO i=1,npoi
          rwork2 = 1.0_r8 - fiold(i)
          heatg(i) = rwork2*heatg(i)  &
               + fiold(i)*sflo(i,nsnolay+2) &
               + zheat(i)*rwork
          solg(i)  = rwork2 * solg(i)
          fvapg(i) = rwork2 * fvapg(i)
          x1(i)    = rwork2 * raing(i)
          x2(i)    = zmelt(i)*rwork
          x3(i)    = tmelt
       END DO
       !
       traing2=traing
       CALL MixSnow (&
            raing   , &! INTENT(OUT  )
            traing  , &! INTENT(OUT  )
            x1      , &! INTENT(IN   )
            traing2 , &! INTENT(IN   )
            x2      , &! INTENT(IN   )
            x3      , &! INTENT(IN   )
            vzero   , &! INTENT(IN   )
            vzero   , &! INTENT(IN   )
            npoi    , &! INTENT(IN   )
            epsilon   )! INTENT(IN   )
       !
    END IF
    !
    RETURN
  END SUBROUTINE snow
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE snowheat (tlsub   , &! INTENT(INOUT)
       fhtop   , &! INTENT(IN   )
       sflo    , &! INTENT(out  )
       xl      , &! INTENT(IN   )
       chl     , &! INTENT(IN   )
       hsoi    , &! INTENT(IN   )
       consoi  , &! INTENT(IN   )
       tsoi    , &! INTENT(IN   )
       tsno    , &! INTENT(INOUT)
       consno  , &! INTENT(IN   )
       hsno    , &! INTENT(IN   )
       rhos    , &! INTENT(IN   )
       fi      , &! INTENT(IN   )
       vegtype0 , &! INTENT(IN   )
       npoi    , &! INTENT(IN   )
       nsnolay , &! INTENT(IN   )
       nsoilay , &! INTENT(IN   )
       dtime   , &! INTENT(IN   )
       cice    , &! INTENT(IN   )
       nVegClass ) ! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    ! sets up call to tridia to solve implicit snow heat conduction,
    ! using snow temperatures in tsno (in comsno). adds an extra
    ! buried-lower-veg layer to the bottom of the snow with 
    ! conduction coefficient conbur/xl and heat capacity chl*xl
    !

    IMPLICIT NONE
    !
    !include 'compar.h'
    INTEGER, INTENT(IN   ) ::       nVegClass
    INTEGER, INTENT(IN   ) :: npoi    ! total number of land points
    INTEGER, INTENT(IN   ) :: nsnolay ! number of snow layers
    INTEGER, INTENT(IN   ) :: nsoilay ! number of soil layers
    REAL(KIND=r8)   , INTENT(IN   ) :: dtime ! model timestep (seconds)
    REAL(KIND=r8)   , INTENT(IN   ) :: cice    ! specific heat of ice (J deg-1 kg-1)

    !      include 'comsno.h'
    REAL(KIND=r8)   , INTENT(INOUT) :: tsno  (npoi,nsnolay)   ! temperature of snow layers (K)
    REAL(KIND=r8)   , INTENT(IN   ) :: consno                        ! thermal conductivity of snow (W m-1 K-1)
    REAL(KIND=r8)   , INTENT(IN   ) :: hsno  (npoi,nsnolay)   ! thickness of snow layers (m)
    REAL(KIND=r8)   , INTENT(IN   ) :: rhos                          ! density of snow (kg m-3)
    REAL(KIND=r8)   , INTENT(IN   ) :: fi    (npoi)                ! fractional snow cover

    !      include 'comsoi.h'
    REAL(KIND=r8)    , INTENT(IN   ) :: hsoi  (npoi,nsoilay+1)          ! soil layer thickness (m)
    REAL(KIND=r8)    , INTENT(IN   ) :: consoi(npoi,nsoilay)     ! thermal conductivity of each soil layer (W m-1 K-1)
    REAL(KIND=r8)    , INTENT(IN   ) :: tsoi  (npoi,nsoilay)          ! soil temperature for each layer (K)

    !
    ! Arguments
    !
    REAL(KIND=r8)   , INTENT(IN   ) :: chl(1:nVegClass)             ! specific heat of lower veg per l/s area (supplied)
    REAL(KIND=r8)   , INTENT(INOUT) :: tlsub(npoi)           ! temperature of buried lower veg (supplied, returned)
    REAL(KIND=r8)   , INTENT(IN   ) :: fhtop(npoi)           ! heat flux into top snow layer from atmos (supplied)
    REAL(KIND=r8)   , INTENT(OUT  ) :: sflo (npoi,nsnolay+2) ! downward heat flow across layer boundaries (returned)
    REAL(KIND=r8)   , INTENT(IN   ) :: xl   (npoi)           ! (lai(i,1)+sai(i,1))*fl(i), lower-veg density(supplied)
    REAL(KIND=r8)   , INTENT(IN   ) :: vegtype0(npoi) 
    !
    ! Local variables
    !
    INTEGER :: iveg               ! loop indices 
    INTEGER :: k               ! loop indices 
    INTEGER :: i               ! loop indices
    INTEGER :: km1             ! used to avoid layer 0
    INTEGER :: kp1             ! used to avoid layer nsnolay+2
    !
    REAL(KIND=r8)   , PARAMETER :: rimp =1.00_r8           ! implicit fraction of the calculation (0 to 1)
    REAL(KIND=r8)   , PARAMETER :: conbur= 2.00_r8         ! conduction coeff of buried lower veg layer 
    ! for unit density xl=(lai+sai)*fl, in w m-2 k-1
    ! conbur (for xl=1) is chosen to be equiv to 10 cm of snow
    REAL(KIND=r8)   , PARAMETER :: hfake = 0.010_r8        ! arbitrary small thickness to allow processing 
    ! for zero snow. (doesn't use index since tridia
    ! not set up for index.)
    REAL(KIND=r8)    :: rwork                  ! to compute matrix diagonals and right-hand side
    REAL(KIND=r8)    :: dt         ! '
    REAL(KIND=r8)    :: dti         ! '
    !
    REAL(KIND=r8)    :: con (npoi,nsnolay+2)   ! conduction coefficents between layers
    REAL(KIND=r8)    :: temp(npoi,nsnolay+1)   ! combined snow and buried-veg temperatures
    REAL(KIND=r8)    :: d1  (npoi,nsnolay+1)   ! diagonals of tridiagonal systems of equations 
    REAL(KIND=r8)    :: d2  (npoi,nsnolay+1)   ! '
    REAL(KIND=r8)    :: d3  (npoi,nsnolay+1)   ! '
    REAL(KIND=r8)    :: rhs (npoi,nsnolay+1)   ! right-hand sides of systems of equations
    REAL(KIND=r8)    :: w1  (npoi,nsnolay+1)   ! work array needed by tridia
    REAL(KIND=r8)    :: w2  (npoi,nsnolay+1)   ! '
    !
    ! copy snow and buried-lower-veg temperatures into combined
    ! array temp
    !

    DO k=1,nsnolay
       DO i=1,npoi
          temp(i,k) =tsno(i,k) 
       END DO
    END DO
    !      CALL scopy (npoi*nsnolay     , & ! INTENT(IN   )
    !                  tsno             , & ! INTENT(IN   )
    !          temp               ) ! INTENT(OUT  )
    DO i=1,npoi     
       temp(i,nsnolay+1)=tlsub(i)
    END DO

    !      CALL scopy (npoi             , & ! INTENT(IN   )
    !                  tlsub            , & ! INTENT(IN   )
    !          temp(1,nsnolay+1)  ) ! INTENT(OUT  )
    !
    ! set conduction coefficients between layers
    !
    DO k=1,nsnolay+2
       IF (k.EQ.1) THEN
          DO i=1,npoi                  
             con(i,k) = 0.0_r8
          END DO
          !         CALL const (con(1,k)   , & !INTENT(OUT  )
          !                      npoi       , & !INTENT(IN   )
          !              0.0_r8          ) !INTENT(IN   )
          !
       ELSE IF (k.LE.nsnolay) THEN
          rwork = 0.50_r8 / consno
          DO i=1,npoi
             con(i,k) = 1.0_r8 / (   MAX(hsno(i,k-1),hfake)*rwork  &
                  + MAX(hsno(i,k)  ,hfake)*rwork )
          END DO
          !
       ELSE IF (k.EQ.nsnolay+1) THEN
          rwork = 0.50_r8 / consno
          DO  i=1,npoi
             con(i,k) = 1.0_r8 / (   MAX(hsno(i,k-1),hfake)*rwork  &
                  + 0.50_r8*xl(i)/conbur )
          END DO
          !
       ELSE IF (k.EQ.nsnolay+2) THEN
          rwork = 0.50_r8 / conbur
          DO i=1,npoi
             con(i,k) = 1.0_r8 / (   xl(i)*rwork                  & 
                  + 0.50_r8*hsoi(i,1) / consoi(i,1) )
          END DO
       END IF
    END DO
    !
    ! set matrix diagonals and right-hand side. for layer nsnolay+1
    ! (buried-lower-veg layer), use explicit contact with soil, and
    ! multiply eqn through by xl*chl/dtime to allow zero xl.
    !
    DO k=1,nsnolay+1
       km1 = MAX (k-1,1)
       kp1 = MIN (k+1,nsnolay+1)
       !
       IF (k.LE.nsnolay) THEN
          rwork = dtime /(rhos*cice)
          DO i=1,npoi
             dt = rwork / (MAX(hsno(i,k),hfake))
             d1(i,k) =    - dt*rimp* con(i,k)
             d2(i,k) = 1.0_r8 + dt*rimp*(con(i,k)+con(i,k+1))
             d3(i,k) =    - dt*rimp* con(i,k+1)
             !
             rhs(i,k) = temp(i,k) + dt                                  &
                  * ( (1.0_r8-rimp)*con(i,k  )*(temp(i,km1)-temp(i,k))  &
                  +   (1.0_r8-rimp)*con(i,k+1)*(temp(i,kp1)-temp(i,k)) )
          END DO
          !
          IF (k.EQ.1) THEN 
             rwork = dtime /(rhos*cice)
             DO i=1,npoi
                dt = rwork / (MAX(hsno(i,k),hfake))
                rhs(i,k) = rhs(i,k) + dt*fhtop(i)
             END DO
          END IF
          !
       ELSE IF (k.EQ.nsnolay+1) THEN
          !
          DO i=1,npoi
             iveg= vegtype0(i)
             rwork = chl(iveg) / dtime
             dti = xl(i)*rwork
             d1(i,k) =     -  rimp* con(i,k)
             d2(i,k) = dti +  rimp*(con(i,k)+con(i,k+1))
             d3(i,k) = 0.0_r8
             rhs(i,k) = dti*temp(i,k)                                 &
                  + ( (1.0_r8-rimp)*con(i,k)*(temp(i,km1)-temp(i,k))  &
                  + con(i,k+1)*(tsoi(i,1)-(1.0_r8-rimp)*temp(i,k)) )
          END DO
       END IF
    END DO
    !
    ! solve the tridiagonal systems
    !
    CALL tridia (npoi      , & ! INTENT(IN   )
         npoi      , & ! INTENT(IN   )
         nsnolay+1 , & ! INTENT(IN   )
         d1        , & ! INTENT(IN   )
         d2        , & ! INTENT(IN   )
         d3        , & ! INTENT(IN   )
         rhs       , & ! INTENT(IN   )
         temp      , & ! INTENT(INOUT)
         w1        , & ! INTENT(INOUT)
         w2          ) ! INTENT(INOUT)
    !
    ! deduce downward heat fluxes between layers
    !
    DO i=1, npoi
       sflo (i,1) = fhtop(i)
    END DO
    !      CALL scopy (npoi     , & ! INTENT(IN   )
    !                  fhtop    , & ! INTENT(IN   )
    !          sflo(1,1)  ) ! INTENT(OUT  )
    !
    DO k=1,nsnolay+1
       IF (k.LE.nsnolay) THEN
          rwork = rhos*cice/dtime
          DO i=1,npoi
             sflo(i,k+1) = sflo(i,k) - rwork*hsno(i,k)   &
                  &                          *(temp(i,k)-tsno(i,k))
          END DO
          !
       ELSE
          !
          DO i=1,npoi
             iveg= vegtype0(i)
             rwork = chl(iveg)/dtime

             sflo(i,k+1) = sflo(i,k)                                  &
                  - xl(i)*rwork*(temp(i,nsnolay+1)-tlsub(i))
          END DO
       END IF
    END DO
    !
    ! copy temperature solution to tsno and tlsub, but not for
    ! points with no snow
    !
    DO k=1,nsnolay
       DO i=1,npoi
          IF (fi(i).GT.0.0_r8) tsno(i,k) = temp(i,k) 
       END DO
    END DO
    !
    DO  i=1,npoi
       IF (fi(i).GT.0.0_r8) tlsub(i) = temp(i,nsnolay+1)
    END DO
    !
    RETURN
  END SUBROUTINE snowheat
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE vadapt (hcur   , &! INTENT(INOUT)
       tcur   , &! INTENT(INOUT)
       htop   , &! INTENT(IN   )
       indp   , &! INTENT(IN   )
       np     , &! INTENT(IN   )
       nlay   , &! INTENT(IN   )
       npoi   , &! INTENT(IN   )
       nsnolay, &! INTENT(IN   )
       epsilon  )! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    ! re-adapt snow layer thicknesses, so top thickness
    ! equals hsnotop and other thicknesses are equal
    !
    ! also adjusts profile of tracer field tcur so its vertical
    ! integral is conserved (eg, temperature)
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN   ) :: npoi               ! total number of land points
    INTEGER, INTENT(IN   ) :: nsnolay       ! number of snow layers
    REAL(KIND=r8)   , INTENT(IN   ) :: epsilon       ! small quantity to avoid zero-divides and other
    ! truncation or machine-limit troubles with small
    ! values. should be slightly greater than o(1)
    ! machine precision

    !
    ! Arguments
    !
    INTEGER, INTENT(IN   ) :: np                 ! number of snow pts in current strip (supplied)
    INTEGER, INTENT(IN   ) :: nlay               ! # of layer
    !
    INTEGER, INTENT(IN   ) :: indp (npoi)        ! index of snow pts in current strip (supplied)
    !
    REAL(KIND=r8)   , INTENT(IN   ) :: htop               ! prescribed top layer thickness (supplied)
    !
    REAL(KIND=r8)   , INTENT(INOUT) :: hcur (npoi,nlay)   ! layer thicknesses (supplied and returned)     
    REAL(KIND=r8)   , INTENT(INOUT) :: tcur (npoi,nlay)   ! tracer field (supplied and returned)
    !
    ! local variables
    !
    INTEGER :: i     ! loop indices
    INTEGER :: j     ! loop indices
    INTEGER :: k     ! loop indices
    INTEGER :: ko    ! loop indices
    !
    REAL(KIND=r8)    :: dz 
    REAL(KIND=r8)    :: rwork
    !
    REAL(KIND=r8)    :: ht   (npoi)         ! storing variable for zold        
    REAL(KIND=r8)    :: h1   (npoi)         ! to compute new layer thickness
    REAL(KIND=r8)    :: za   (npoi)         ! 
    REAL(KIND=r8)    :: zb   (npoi)         ! 
    REAL(KIND=r8)    :: zheat(npoi)         !
    !
    REAL(KIND=r8)    :: hnew (npoi,nsnolay)    ! new layer thickness
    REAL(KIND=r8)    :: tnew (npoi,nsnolay)    ! new temperatures of layers
    REAL(KIND=r8)    :: zold (npoi,nsnolay+1)  ! distances from surface to old layer boundaries
    !
    ! if no snow or seaice points in current 1d strip, return. note
    ! that the index is not used below (for cray vec and efficiency)
    ! except in the final loop setting the returned values
    !
    IF (np.EQ.0) RETURN
    !
    ! set distances zold from surface to old layer boundaries
    !
    DO i=1,npoi
       zold(i,1) = 0.0_r8  
    END DO

    !      CALL const (zold(1,1),& !INTENT(OUT  )
    !                  npoi     ,& !INTENT(IN   )
    !          0.0_r8       ) !INTENT(IN   )
    !
    DO k=1,nlay
       DO i=1,npoi
          zold(i,k+1) = zold(i,k) + hcur(i,k)
       END DO
    END DO
    !
    ! set new layer thicknesses hnew (tot thickness is unchanged).
    ! if total thickness is less than nlay*htop (which should be
    ! le hsnomin), make all new layers equal including
    ! top one, so other layers aren't so thin. use epsilon to 
    ! handle zero (snow) points
    !     
    DO i=1,npoi
       ht   (i) =zold(i,nlay+1)
    END DO
    !      CALL scopy (npoi          , & ! INTENT(IN   )
    !                  zold(1,nlay+1), & ! INTENT(IN   )
    !          ht              ) ! INTENT(OUT  )
    !
    rwork = nlay*htop
    DO i=1,npoi
       IF (ht(i).GE.rwork) THEN
          h1(i) = (ht(i)-htop)/(nlay-1)
       ELSE
          h1(i) = MAX (ht(i)/nlay, epsilon)
       END IF
    END DO
    !
    DO k=1,nlay
       DO i=1,npoi
          hnew(i,k) = h1(i)
       END DO
    END DO
    !
    rwork = nlay*htop
    DO i=1,npoi
       IF (ht(i).GE.rwork) hnew(i,1) = htop
    END DO
    !
    ! integrate old temperature profile (loop 410) over each
    ! new layer (loop 400), to get new field tnew
    !
    DO i=1,npoi
       zb   (i) =0.0_r8 
    END DO

    !      CALL const (zb     , & !INTENT(OUT  )
    !                  npoi   , & !INTENT(IN   )
    !          0.0_r8      ) !INTENT(IN   )
    !
    DO k=1,nlay
       !
       DO i=1,npoi
          za(i) = zb(i)
          zb(i) = za(i) + hnew(i,k)
       END DO

       DO i=1,npoi
          zheat   (i) =0.0_r8 
       END DO

       !        CALL const (zheat  , & !INTENT(OUT  )
       !                    npoi   , & !INTENT(IN   )
       !            0.0_r8      ) !INTENT(IN   )
       !
       DO ko=1,nlay
          DO i=1,npoi
             IF (za(i).LT.zold(i,ko+1) .AND. zb(i).GT.zold(i,ko)) THEN
                dz = MIN(zold(i,ko+1),zb(i)) - MAX(zold(i,ko),za(i))
                zheat(i) = zheat(i) + tcur(i,ko)*dz
             END IF
          END DO
       END DO
       !
       DO i=1,npoi
          tnew(i,k) = zheat(i) / hnew(i,k)
       END DO
       !
    END DO
    !
    ! use index for final copy to seaice or snow arrays, to avoid
    ! changing soil values (when called for seaice) and to avoid
    ! changing nominal snow values for no-snow points (when called
    ! for snow)
    !
    DO k=1,nlay
       DO j=1,np
          i = indp(j)
          hcur(i,k) = hnew(i,k)
          tcur(i,k) = tnew(i,k)
       END DO
    END DO
    !
    RETURN
  END SUBROUTINE vadapt
  ! ---------------------------------------------------------------------
  SUBROUTINE MixSnow (&
       xm    , &! INTENT(OUT  )
       tm    , &! INTENT(OUT  )
       x1    , &! INTENT(IN   )
       t1    , &! INTENT(IN   )
       x2    , &! INTENT(IN   )
       t2    , &! INTENT(IN   )
       x3    , &! INTENT(IN   )
       t3    , &! INTENT(IN   )
       npoi  , &! INTENT(IN   )
       epsilon )! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    ! calorimetrically mixes masses x1,x2,x3 with temperatures
    ! t1,t2,t3 into combined mass xm with temperature tm
    !
    ! xm,tm may be returned into same location as one of x1,t1,..,
    ! so hold result temporarily in xtmp,ttmp below
    !
    ! will work if some of x1,x2,x3 have opposite signs, but may 
    ! give unphysical tm's
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN   ) :: npoi    ! total number of land points
    REAL(KIND=r8)   , INTENT(IN   ) :: epsilon ! small quantity to avoid zero-divides and other
    ! truncation or machine-limit troubles with small
    ! values. should be slightly greater than o(1)
    ! machine precision
    !
    ! Arguments (input except for xm, tm)
    !
    REAL(KIND=r8)   , INTENT(OUT  ) :: xm(npoi)     ! resulting mass  
    REAL(KIND=r8)   , INTENT(OUT  ) :: tm(npoi)     ! resulting temp
    REAL(KIND=r8)   , INTENT(IN   ) :: x1(npoi)     ! mass 1
    REAL(KIND=r8)   , INTENT(IN   ) :: t1(npoi)     ! temp 1
    REAL(KIND=r8)   , INTENT(IN   ) :: x2(npoi)     ! mass 2
    REAL(KIND=r8)   , INTENT(IN   ) :: t2(npoi)     ! temp 2
    REAL(KIND=r8)   , INTENT(IN   ) :: x3(npoi)     ! mass 3
    REAL(KIND=r8)   , INTENT(IN   ) :: t3(npoi)     ! temp 3
    !
    ! local variables
    !
    INTEGER :: i            ! loop indice
    !
    REAL(KIND=r8) ::   xtmp         ! resulting mass (storing variable)
    REAL(KIND=r8) ::   ytmp         !  "
    REAL(KIND=r8) ::   ttmp         ! resulting temp
    !
    ! ---------------------------------------------------------------------
    !
    DO  i=1,npoi
       !
       xtmp = x1(i) + x2(i) + x3(i)
       !
       ytmp = SIGN (MAX (ABS(xtmp), epsilon), xtmp)
       !
       IF (ABS(xtmp).GE.epsilon) THEN
          ttmp = (t1(i)*x1(i) + t2(i)*x2(i) + t3(i)*x3(i)) / ytmp
       ELSE
          ttmp = 0.0_r8
          xtmp = 0.0_r8
       END IF
       !
       xm(i) = xtmp
       tm(i) = ttmp
       !
    END DO
    !
    RETURN
  END SUBROUTINE MixSnow
  ! ---------------------------------------------------------------------
  SUBROUTINE tridia (ns      , & ! INTENT(IN   ) :: ns ! number of systems to be solved.
       nd      , & ! INTENT(IN   ) :: nd ! first dimension of arrays (ge ns)
       ne      , & ! INTENT(IN   ) :: ne ! number of unknowns in each system. (>2)
       a       , & ! INTENT(IN   ) :: a(nd,ne)     ! subdiagonals of matrices stored in a(j,2)...a(j,ne).
       b       , & ! INTENT(IN   ) :: b(nd,ne)     ! main diagonals of matrices stored in b(j,1)...b(j,ne).
       c       , & ! INTENT(IN   ) :: c(nd,ne)     ! super-diagonals of matrices stored in c(j,1)...c(j,ne-1).
       y       , & ! INTENT(IN   ) :: y(nd,ne)     ! right hand side of equations stored in y(j,1)...y(j,ne).
       x       , & ! INTENT(INOUT) :: x(nd,ne)     ! solutions of the systems returned in x(j,1)...x(j,ne).
       alpha   , & ! INTENT(INOUT) :: alpha(nd,ne) ! work array 
       gamma     ) ! INTENT(INOUT) :: gamma(nd,ne) ! work array
    ! ---------------------------------------------------------------------
    !

    IMPLICIT NONE
    !
    !     purpose:
    !     to compute the solution of many tridiagonal linear systems.
    !
    !      arguments:
    !
    !      ns ..... the number of systems to be solved.
    !
    !      nd ..... first dimension of arrays (ge ns).
    !
    !      ne ..... the number of unknowns in each system.
    !               this must be > 2. second dimension of arrays.
    !
    !      a ...... the subdiagonals of the matrices are stored
    !               in locations a(j,2) through a(j,ne).
    !
    !      b ...... the main diagonals of the matrices are stored
    !               in locations b(j,1) through b(j,ne).
    !
    !      c ...... the super-diagonals of the matrices are stored in
    !               locations c(j,1) through c(j,ne-1).
    !
    !      y ...... the right hand side of the equations is stored in
    !               y(j,1) through y(j,ne).
    !
    !      x ...... the solutions of the systems are returned in
    !               locations x(j,1) through x(j,ne).
    !
    !      alpha .. work array dimensioned alpha(nd,ne)
    !
    !      gamma .. work array dimensioned gamma(nd,ne)
    !
    !       history:  based on a streamlined version of the old ncar
    !                 ulib subr trdi used in the phoenix climate
    !                 model of schneider and thompson (j.g.r., 1981).
    !                 revised by starley thompson to solve multiple
    !                 systems and vectorize well on the cray-1.
    !                 later revised to include a parameter statement
    !                 to define loop limits and thus enable cray short
    !                 vector loops.
    !
    !       algorithm:  lu decomposition followed by solution.
    !                   note: this subr executes satisfactorily
    !                   if the input matrix is diagonally dominant
    !                   and non-singular.  the diagonal elements are
    !                   used to pivot, and no tests are made to determine
    !                   singularity. if a singular or numerically singular
    !                   matrix is used as input a divide by zero or
    !                   floating point overflow will result.
    !
    !       last revision date:      4 february 1988
    !
    !
    ! Arguments
    !
    ! General information:

    !     Subroutine tridia solves a tridiagonal system of equations for the
    !     n+1 unknowns u(0), u(1), ..., u(n-1), u(n).  This system must
    !     consist of n+1 equations of the form:
    !
    !                   b(0) u(0) + c(0) u(1)   = d(0)   first equation
    !     a(i) u(i-1) + b(i) u(i) + c(i) u(i+1) = d(i)   for i = 1, 2, ..., n-1
    !     a(n) u(n-1) + b(n) u(n)               = d(n)   final equation

    INTEGER, INTENT(IN   ) :: ns ! number of systems to be solved.
    INTEGER, INTENT(IN   ) :: nd ! first dimension of arrays (ge ns)
    INTEGER, INTENT(IN   ) :: ne ! number of unknowns in each system. (>2)

    REAL(KIND=r8)   , INTENT(IN   ) :: a(nd,ne)     ! subdiagonals of matrices stored in a(j,2)...a(j,ne).
    REAL(KIND=r8)   , INTENT(IN   ) :: b(nd,ne)     ! main diagonals of matrices stored in b(j,1)...b(j,ne).
    REAL(KIND=r8)   , INTENT(IN   ) :: c(nd,ne)     ! super-diagonals of matrices stored in c(j,1)...c(j,ne-1).
    REAL(KIND=r8)   , INTENT(IN   ) :: y(nd,ne)     ! right hand side of equations stored in y(j,1)...y(j,ne).
    REAL(KIND=r8)   , INTENT(INOUT) :: x(nd,ne)     ! solutions of the systems returned in x(j,1)...x(j,ne).
    REAL(KIND=r8)   , INTENT(INOUT) :: alpha(nd,ne) ! work array 
    REAL(KIND=r8)   , INTENT(INOUT) :: gamma(nd,ne) ! work array
    !
    ! local variables
    !
    INTEGER :: nm1     ! loop indices
    INTEGER :: j       ! loop indices
    INTEGER :: i       ! loop indices
    INTEGER :: ib      ! loop indices

    !
    nm1 = ne-1
    !
    ! obtain the lu decompositions
    !
    DO j=1,ns
       alpha(j,1) = 1.0_r8/b(j,1)
       gamma(j,1) = c(j,1)*alpha(j,1)
    END DO
    DO i=2,nm1
       DO j=1,ns
          alpha(j,i) = 1.0_r8/(b(j,i)-a(j,i)*gamma(j,i-1))
          gamma(j,i) = c(j,i)*alpha(j,i)
       END DO
    END DO
    !
    ! solve
    !
    DO j=1,ns
       x(j,1) = y(j,1)*alpha(j,1)
    END DO

    DO i=2,nm1
       DO j=1,ns
          x(j,i) = (y(j,i)-a(j,i)*x(j,i-1))*alpha(j,i)
       END DO
    END DO

    DO j=1,ns
       x(j,ne) = (y(j,ne)-a(j,ne)*x(j,nm1))/     &
            (b(j,ne)-a(j,ne)*gamma(j,nm1))
    END DO
    DO i=1,nm1
       ib = ne-i
       DO j=1,ns
          x(j,ib) = x(j,ib)-gamma(j,ib)*x(j,ib+1)
       END DO
    END DO
    !
    RETURN
  END SUBROUTINE tridia



  !####### #     # ######     #####  #         # ####### #         #
  !#         ##    # #     #   #         # ##         # #         # #  #  #
  !#         # #   # #     #   #           # #   # #         # #  #  #
  !#####   #  #  # #     #    #####  #  #  # #         # #  #  #
  !#         #   # # #     #            # #   # # #         # #  #  #
  !#         #    ## #     #   #         # ### #         # #  #  #
  !####### #     # ######     #####  #         # #######  ## ##
  !






  !
  !  ####    ####      #    #
  ! #       #    #     #    #
  !  ####   #    #     #    #
  !      #  #    #     #    #
  ! #    #  #    #     #    #
  !  ####    ####      #    ######
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE setsoi(npoi   , &! INTENT(IN   )
       nsoilay, &! INTENT(IN   )
       sand   , &! INTENT(IN   )
       clay   , &! INTENT(IN   )
       poros  , &! INTENT(IN   )
       wsoi   , &! INTENT(IN   )
       wisoi  , &! INTENT(IN   )
       consoi , &! INTENT(OUT  )
       zwpmax , &! INTENT(IN   )
       wpud   , &! INTENT(IN   )
       wipud  , &! INTENT(IN   )
       wpudmax, &! INTENT(IN   )
       qglif  , &! INTENT(OUT  )
       tsoi   , &! INTENT(IN   )
       hvasug , &! INTENT(OUt  )
       hvasui , &! INTENT(OUt  )
       tsno   , &! INTENT(IN   )
       ta     , &! INTENT(IN   )
       nsnolay, &! INTENT(IN   )
       hvap   , &! INTENT(IN   )
       cvap   , &! INTENT(IN   )
       ch2o   , &! INTENT(IN   )
       hsub   , &! INTENT(IN   )
       cice   , &! INTENT(IN   )
       tmelt  , &! INTENT(IN   )
       epsilon  )! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    ! sets diagnostic soil quantities
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: npoi                 ! total number of land points
    INTEGER, INTENT(IN   ) :: nsoilay              ! number of soil layers
    INTEGER, INTENT(IN   ) :: nsnolay              ! number of snow layers
    REAL(KIND=r8)   , INTENT(IN   ) :: hvap                   ! latent heat of vaporization of water (J kg-1)
    REAL(KIND=r8)   , INTENT(IN   ) :: cvap                   ! specific heat of water vapor at constant pressure (J deg-1 kg-1)
    REAL(KIND=r8)   , INTENT(IN   ) :: ch2o                   ! specific heat of liquid water (J deg-1 kg-1)
    REAL(KIND=r8)   , INTENT(IN   ) :: hsub                   ! latent heat of sublimation of ice (J kg-1)
    REAL(KIND=r8)   , INTENT(IN   ) :: cice                   ! specific heat of ice (J deg-1 kg-1)
    REAL(KIND=r8)   , INTENT(IN   ) :: tmelt                   ! freezing point of water (K)
    REAL(KIND=r8)   , INTENT(IN   ) :: epsilon              ! small quantity to avoid zero-divides and other
    ! truncation or machine-limit troubles with small
    ! values. should be slightly greater than o(1)
    ! machine precision
    REAL(KIND=r8)   , INTENT(IN   ) :: sand   (npoi,nsoilay)! percent sand of soil
    REAL(KIND=r8)   , INTENT(IN   ) :: clay   (npoi,nsoilay)! percent clay of soil
    REAL(KIND=r8)   , INTENT(IN   ) :: poros  (npoi,nsoilay)! porosity (mass of h2o per unit vol at sat / rhow)
    REAL(KIND=r8)   , INTENT(IN   ) :: wsoi   (npoi,nsoilay)! fraction of soil pore space containing liquid water
    REAL(KIND=r8)   , INTENT(IN   ) :: wisoi  (npoi,nsoilay)! fraction of soil pore space containing ice
    REAL(KIND=r8)   , INTENT(OUT  ) :: consoi (npoi,nsoilay)! thermal conductivity of each soil layer (W m-1 K-1)
    REAL(KIND=r8)   , INTENT(IN   ) :: zwpmax               ! assumed maximum fraction of soil surface 
    REAL(KIND=r8)   , INTENT(IN   ) :: wpud   (npoi)        ! liquid content of puddles per soil area (kg m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: wipud  (npoi)        ! ice content of puddles per soil area (kg m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: wpudmax              ! normalization constant for puddles (kg m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: qglif  (npoi,4)      ! 1: fraction of soil evap (fvapg) from soil liquid
    ! 2: fraction of soil evap (fvapg) from soil ice
    ! 3: fraction of soil evap (fvapg) from puddle liquid
    ! 4: fraction of soil evap (fvapg) from puddle ice
    REAL(KIND=r8)   , INTENT(IN   ) :: tsoi   (npoi,nsoilay)! soil temperature for each layer (K)
    REAL(KIND=r8)   , INTENT(OUT  ) :: hvasug (npoi)        ! latent heat of vap/subl, for soil surface (J kg-1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: hvasui (npoi)        ! latent heat of vap/subl, for snow surface (J kg-1)
    REAL(KIND=r8)   , INTENT(IN   ) :: tsno   (npoi,nsnolay)! temperature of snow layers (K)
    REAL(KIND=r8)   , INTENT(IN   ) :: ta     (npoi)        ! air temperature (K)

    !
    ! Local variables
    !
    INTEGER i, k            ! loop indices
    INTEGER msand           ! % of sand in grid point
    INTEGER mclay           ! % of clay in grid point
    !
    REAL(KIND=r8)    fsand           ! fraction of sand in grid point
    REAL(KIND=r8)    fsilt           ! fraction of silt in grid point
    REAL(KIND=r8)    fclay           ! fraction of clay in grid point
    ! MEM: added forganic for organic soils.
    !      REAL(KIND=r8)    forganic        ! fraction of organic soil in grid point
    REAL(KIND=r8)    powliq              ! liquid water content in fraction of soil depth
    REAL(KIND=r8)    powice              ! ice water content in fraction of soil depth
    REAL(KIND=r8)    zcondry              ! dry-soil conductivity
    REAL(KIND=r8)    zvap              ! latent heat of vaporisation at soil temp
    REAL(KIND=r8)    zsub              ! latent heat of sublimation at soil temp
    REAL(KIND=r8)    zwpud              ! fraction of soil surface covered by puddle
    !             zwpmax              ! assumed maximum value of zwpud
    REAL(KIND=r8)    zwsoi              ! volumetric water content of top soil layer 
    REAL(KIND=r8)   ::  rwork
    REAL(KIND=r8)    rwork1
    REAL(KIND=r8)    rwork2
    REAL(KIND=r8)   zwtot
    rwork =zwpmax
    !
    !
    ! set soil layer quantities
    !
    DO  k = 1, nsoilay
       !
       DO  i = 1, npoi
          !
          ! Convert input sand and clay percents to fractions
          !
          msand = NINT(sand(i,k))
          mclay = NINT(clay(i,k))
          !
          fsand = 0.01_r8 * msand
          fclay = 0.01_r8 * mclay
          fsilt = MAX(0.01_r8 * (100 - msand - mclay),0.0007_r8)
          !
          ! update thermal conductivity (w m-1 k-1)
          !
          ! based on c = c1**v1 * c2**v2 * c3**v3 * c4**v4 where c1,c2..
          ! are conductivities of soil grains, air, liquid and ice
          ! respectively, and v1,v2... are their volume fractions 
          ! (so v1 = 1-p where p is the porosity, and v1+v2+v3+v4 = 1).
          ! then condry = c1**(1-p) * c2**p  is the dry-soil
          ! conductivity, and c = condry * (c3/c2)**v3 * (c4/c2)**v4, 
          ! where c2 = conductivity of air = .025 w m-1 k-1.
          ! however this formula agrees better with williams+smith
          ! table 4 for wet (unfrozen) sand and clay if c2 is decreased
          ! to ~.005. (for peat in next section, ok if c2 = .025).
          ! also see lachenbruch etal,1982,jgr,87,9301 and refs therein.
          !
          powliq = poros(i,k) * wsoi(i,k) * (1.0_r8 - wisoi(i,k))
          powice = poros(i,k) * wisoi(i,k)
          !
          zcondry = fsand * 0.300_r8 + fsilt * 0.265_r8 + fclay * 0.250_r8 ! +
          !zcondry = fsand * 0.500_r8 + fsilt * 0.265_r8 + fclay * 0.250_r8 ! +

          ! M. El Maayar added this to account for contribution of organic soils
          !     >              forganic * 0.026
          !
          consoi(i,k) = zcondry * ((0.56_r8*100.0_r8)**powliq)  &
               * ((2.24_r8*100.0_r8)**powice)
          !
       END DO
       !
    END DO


    !
    ! set qglif - the fraction of soil sfc evaporation from soil
    ! liquid (relative to total from liquid and ice)
    !
    ! 1: fraction of soil evap (fvapg) from soil liquid
    ! 2: fraction of soil evap (fvapg) from soil ice
    ! 3: fraction of soil evap (fvapg) from puddle liquid
    ! 4: fraction of soil evap (fvapg) from puddle ice
    !
    DO  i = 1, npoi        
    !              zwpmax = 0.5
       zwpud = max (0.0_r8, min (zwpmax,    &
             zwpmax*(wpud(i)+wipud(i))/wpudmax) )



       !
       zwtot = (1.00_r8 - wisoi(i,1)) * wsoi(i,1) + wisoi(i,1) &
            + (wpud(i) + wipud(i)) / wpudmax
       zwsoi = MAX((1.00_r8 - wisoi(i,1)) * wsoi(i,1) + wisoi(i,1),0.001_r8)

       !
       IF (zwtot.GE.epsilon) THEN
          !
          ! for a wet surface
          !
          rwork = 1.00_r8 / zwtot
          rwork1 = 1.0_r8/zwsoi

          !
          !qglif(i,1) = (1.0_r8 - wisoi(i,1)) * wsoi(i,1) * rwork1
          !qglif(i,2) =       wisoi(i,1)              * rwork1
          !qglif(i,3) = (wpud(i)  / wpudmax) * rwork1
          !qglif(i,4) = (wipud(i) / wpudmax) * rwork1
          IF (zwpud.ge.epsilon) THEN
            rwork2 = 1./(wpud(i) + wipud(i))
            qglif(i,1) = (1.0_r8 - zwpud) * (1.0_r8 - wisoi(i,1)) * wsoi(i,1) * rwork1
            qglif(i,2) = (1.0_r8 - zwpud) * wisoi(i,1) * rwork1
            qglif(i,3) = zwpud * wpud(i) * rwork2
            qglif(i,4) = zwpud * wipud(i) * rwork2
          ELSE
            qglif(i,1) = (1.0_r8 - wisoi(i,1)) * wsoi(i,1) * rwork1
            qglif(i,2) = wisoi(i,1) * rwork1
            qglif(i,3) = 0.00_r8
            qglif(i,4) = 0.00_r8
          END IF

          !
       ELSE
          !c
          !c for a 100% dry soil surface, assign all soil evap to the puddles.
          !c Note that for small puddle sizes, this could lead to negative
          !c puddle depths. However, for a 100% dry soil with small puddles,
          !c evaporation is likely to be very small or less than zero
          !c (condensation), so negative puddle depths are not likely to occur.
          !c
          IF (zwpud.ge.epsilon) THEN
            rwork2 = 1./(wpud(i) + wipud(i))
            qglif(i,1) = 0.0
            qglif(i,2) = 0.0
            qglif(i,3) = zwpud * wpud(i) * rwork2
            qglif(i,4) = zwpud * wipud(i) * rwork2
          ELSE

             IF (tsoi(i,1).GE.tmelt) THEN
                !
                ! above freezing
                !
                qglif(i,1) = 0.0_r8
                qglif(i,2) = 0.0_r8
                qglif(i,3) = 1.0_r8
                qglif(i,4) = 0.0_r8
                !
             ELSE
                !
                ! below freezing
                !
                qglif(i,1) = 0.0_r8
                qglif(i,2) = 0.0_r8
                qglif(i,3) = 0.0_r8
                qglif(i,4) = 1.0_r8
                !
             END IF
          !
          END IF
       END IF
       !
       ! set latent heat values
       !
       !        zvap = hvapf (tsoi(i,1), ta(i))
       zvap = hvap + cvap*(ta(i)-273.16_r8) - ch2o*(tsoi(i,1)-273.16_r8)

       !        zsub = hsubf (tsoi(i,1), ta(i))
       zsub = hsub + cvap*(ta(i)-273.16_r8) - cice*(tsoi(i,1)-273.16_r8)

       !
       !        hvasug(i) = (qglif(i,1) + qglif(i,3)) * zvap + &
       !                    (qglif(i,2) + qglif(i,4)) * zsub 
       !
       !        hvasui(i) = hsubf(tsno(i,1),ta(i))
       !       !zvap = hvapf (tsoi(i,1), ta(i),hvap,cvap,ch2o)
       !zsub = hsubf (tsoi(i,1), ta(i),hsub,cvap,cice)
       !
       hvasug(i) = (qglif(i,1) + qglif(i,3)) * zvap +  &
            (qglif(i,2) + qglif(i,4)) * zsub 
       !
       !       !hvasui(i) = hsubf(tsno(i,1),ta(i),hsub,cvap,cice)
       hvasui(i) = hsub + cvap*(ta(i)-273.16_r8) - cice*(tsno(i,1)-273.16_r8)

       !
    END DO

    !
    ! set qglif - the fraction of soil sfc evaporation from soil liquid,
    ! soil ice, puddle liquid, and puddle ice (relative to total sfc evap)
    !
    ! zwpud:   fraction of surface area covered by puddle (range: 0 - zwpmax)
    ! zwpmax:  maximum value of zwpud (currently assumed to be 0.5)
    ! 1-zwpud: fraction of surface area covered by soil (range: (1-zwpmax) - 1.0)
    ! zwsoi:   volumetric water content of top soil layer (range: 0 - 1.0)
    !
    ! qglif(i,1): fraction of soil evap (fvapg) from soil liquid
    ! qglif(i,2): fraction of soil evap (fvapg) from soil ice
    ! qglif(i,3): fraction of soil evap (fvapg) from puddle liquid
    ! qglif(i,4): fraction of soil evap (fvapg) from puddle ice
    !
    !     DO  i = 1, npoi
    !
    !        zwpmax = 0.5
    !       zwpud = max (0.0_r8, min (zwpmax,    &
    !             zwpmax*(wpud(i)+wipud(i))/wpudmax) )
    !       zwsoi = (1.0_r8 - wisoi(i,1)) * wsoi(i,1) + wisoi(i,1)
    !
    !       IF (zwsoi.ge.epsilon) THEN
    !
    !         rwork1 = 1.0_r8/zwsoi
    !
    !         IF (zwpud.ge.epsilon) THEN
    !           rwork2 = 1.0_r8/(wpud(i) + wipud(i))
    !           qglif(i,1) = (1.0_r8 - zwpud) * (1.0_r8 - wisoi(i,1)) *  &
    !               wsoi(i,1) * rwork1
    !           qglif(i,2) = (1.0_r8 - zwpud) * wisoi(i,1) * rwork1
    !           qglif(i,3) = zwpud * wpud(i) * rwork2
    !           qglif(i,4) = zwpud * wipud(i) * rwork2
    !         ELSE
    !           qglif(i,1) = (1.0_r8 - wisoi(i,1)) * wsoi(i,1) * rwork1
    !           qglif(i,2) = wisoi(i,1) * rwork1
    !           qglif(i,3) = 0.0_r8
    !           qglif(i,4) = 0.0_r8
    !         END IF
    !
    !       ELSE
    !
    ! for a 100% dry soil surface, assign all soil evap to the puddles.
    ! Note that for small puddle sizes, this could lead to negative
    ! puddle depths. However, for a 100% dry soil with small puddles,
    ! evaporation is likely to be very small or less than zero
    ! (condensation), so negative puddle depths are not likely to occur.
    !
    !         IF (zwpud.ge.epsilon) THEN
    !           rwork2 = 1.0_r8/(wpud(i) + wipud(i))
    !           qglif(i,1) = 0.0_r8
    !           qglif(i,2) = 0.0_r8
    !           qglif(i,3) = zwpud * wpud(i) * rwork2
    !           qglif(i,4) = zwpud * wipud(i) * rwork2
    !         ELSE
    !           IF (tsoi(i,1).ge.tmelt) THEN
    !
    ! above freezing
    !
    !             qglif(i,1) = 0.0_r8
    !             qglif(i,2) = 0.0_r8
    !             qglif(i,3) = 1.0_r8
    !             qglif(i,4) = 0.0_r8
    !
    !           ELSE
    !
    ! below freezing
    !
    !             qglif(i,1) = 0.0_r8
    !             qglif(i,2) = 0.0_r8
    !             qglif(i,3) = 0.0_r8
    !             qglif(i,4) = 1.0_r8
    !           END IF
    !         END IF
    !
    !       END IF
    !
    ! set latent heat values
    !
    !       !zvap = hvapf (tsoi(i,1), ta(i),hvap,cvap,ch2o)
    !        zvap = hvap + cvap*(ta(i)-273.16_r8) - ch2o*(tsoi(i,1)-273.16_r8)
    !zsub = hsubf (tsoi(i,1), ta(i),hsub,cvap,cice)
    !        zsub = hsub + cvap*(ta(i)-273.16_r8) - cice*(tsoi(i,1)-273.16_r8)
    !
    !       hvasug(i) = (qglif(i,1) + qglif(i,3)) * zvap +  &
    !                   (qglif(i,2) + qglif(i,4)) * zsub 
    !
    !       !hvasui(i) = hsubf(tsno(i,1),ta(i),hsub,cvap,cice)
    !hvasui(i) = hsub + cvap*(ta(i)-273.16_r8) - cice*(tsno(i,1)-273.16_r8)
    !
    !     END DO   
    !
    RETURN
  END   SUBROUTINE setsoi
  ! ---------------------------------------------------------------------
  SUBROUTINE soilctl(raing   ,&!  INTENT(IN   )
       fvapg   ,&!  INTENT(IN   )
       traing  ,&!  INTENT(IN   )
       tu      ,&!  INTENT(IN   )
       tl      ,&!  INTENT(IN   )
       wpud    ,&!  INTENT(INOUT)
       wipud   ,&!  INTENT(INOUT)
       wpudmax ,&!  INTENT(IN   )
       tsoi    ,&!  INTENT(INOUT)
       qglif   ,&!  INTENT(IN   )
       wisoi   ,&!  INTENT(INOUT)
       poros   ,&!  INTENT(IN   )
       hsoi    ,&!  INTENT(IN   )
       wsoi    ,&!  INTENT(INOUT)
       hydraul ,&!  INTENT(IN   )
       heatg   ,&!  INTENT(IN   )
       porosflo,&!  INTENT(INOUT)
       upsoiu  ,&!  INTENT(IN   )
       upsoil  ,&!  INTENT(IN   )
       csoi    ,&!  INTENT(IN   )
       rhosoi  ,&!  INTENT(IN   )
       suction ,&!  INTENT(IN   )
       bex     ,&!  INTENT(IN   )
       ibex    ,&!  INTENT(IN   )
       bperm   ,&!  INTENT(IN   )
       consoi  ,&!  INTENT(IN   )
       hflo    ,&!  INTENT(INOUT)
       grunof  ,&!  INTENT(OUT  )
       gadjust ,&!  INTENT(INOUT)
       gdrain  ,&!  INTENT(OUT  )
       npoi    ,&!  INTENT(IN   )
       nsoilay ,&!  INTENT(IN   )
       epsilon ,&!  INTENT(IN   )
       dtime   ,&!  INTENT(IN   )
       ch2o    ,&!  INTENT(IN   )
       cice    ,&!  INTENT(IN   )
       tmelt   ,&!  INTENT(IN   )
       rhow    ,&!  INTENT(IN   )
       hfus    ,&!  INTENT(IN   )
       stressl ,&! soil moisture stress factor for the lower canopy (dimensionless)
       stressu ,&! soil moisture stress factor for the upper canopy (dimensionless)
       stresstl,&! sum of stressl over all 6 soil layers (dimensionless)
       stresstu  )! sum of stressu over all 6 soil layers (dimensionless)

    ! ---------------------------------------------------------------------
    !
    ! steps soil/seaice model through one timestep
    !
    IMPLICIT NONE
    !
    INTEGER,  INTENT(IN   ) :: npoi          ! total number of land points
    INTEGER,  INTENT(IN   ) :: nsoilay       ! number of soil layers
    REAL(KIND=r8)   ,  INTENT(IN   ) :: epsilon       ! small quantity to avoid zero-divides and other
    ! truncation or machine-limit troubles with small
    ! values. should be slightly greater than o(1)
    ! machine precision
    REAL(KIND=r8)   ,  INTENT(IN   ) :: dtime         ! model timestep (seconds)
    REAL(KIND=r8)   ,  INTENT(IN   ) :: ch2o          ! specific heat of liquid water (J deg-1 kg-1)
    REAL(KIND=r8)   ,  INTENT(IN   ) :: cice          ! specific heat of ice (J deg-1 kg-1)
    REAL(KIND=r8)   ,  INTENT(IN   ) :: tmelt         ! freezing point of water (K)
    REAL(KIND=r8)   ,  INTENT(IN   ) :: rhow          ! density of liquid water (all types) (kg m-3)
    REAL(KIND=r8)   ,  INTENT(IN   ) :: hfus          ! latent heat of fusion of water (J kg-1)
    REAL(KIND=r8)   ,  INTENT(OUT  ) :: grunof (npoi)            ! surface runoff rate (kg_h2o m-2 s-1)
    REAL(KIND=r8)   ,  INTENT(INOUT) :: gadjust(npoi)            ! h2o flux due to adjustments in subroutine wadjust (kg_h2o m-2 s-1)
    REAL(KIND=r8)   ,  INTENT(OUT  ) :: gdrain (npoi)            ! drainage rate out of bottom of lowest soil layer (kg_h2o m-2 s-1)
    REAL(KIND=r8)   ,  INTENT(INOUT) :: wpud    (npoi)           ! liquid content of puddles per soil area (kg m-2)
    REAL(KIND=r8)   ,  INTENT(INOUT) :: wipud   (npoi)           ! ice content of puddles per soil area (kg m-2)
    REAL(KIND=r8)   ,  INTENT(IN   ) :: wpudmax                  ! normalization constant for puddles (kg m-2)
    REAL(KIND=r8)   ,  INTENT(INOUT) :: tsoi    (npoi,nsoilay)   ! soil temperature for each layer (K)
    REAL(KIND=r8)   ,  INTENT(IN   ) :: qglif   (npoi,4)         ! 1: fraction of soil evap (fvapg) from soil liquid
    ! 2: fraction of soil evap (fvapg) from soil ice
    ! 3: fraction of soil evap (fvapg) from puddle liquid
    ! 4: fraction of soil evap (fvapg) from puddle ice
    REAL(KIND=r8)   ,  INTENT(INOUT) :: wisoi   (npoi,nsoilay)   ! fraction of soil pore space containing ice
    REAL(KIND=r8)   ,  INTENT(IN   ) :: poros   (npoi,nsoilay)   ! porosity (mass of h2o per unit vol at sat / rhow)
    REAL(KIND=r8)   ,  INTENT(IN   ) :: hsoi    (npoi,nsoilay+1)      ! soil layer thickness (m)
    REAL(KIND=r8)   ,  INTENT(INOUT) :: wsoi    (npoi,nsoilay)   ! fraction of soil pore space containing liquid water
    REAL(KIND=r8)   ,  INTENT(IN   ) :: hydraul (npoi,nsoilay)   ! saturated hydraulic conductivity (m/s)
    REAL(KIND=r8)   ,  INTENT(IN   ) :: heatg   (npoi)           ! net heat flux into soil surface (W m-2)
    REAL(KIND=r8)   ,  INTENT(INOUT) :: porosflo(npoi,nsoilay)   ! porosity after reduction by ice content
    REAL(KIND=r8)   ,  INTENT(IN   ) :: upsoiu  (npoi,nsoilay)   ! soil water uptake from transpiration (kg_h2o m-2 s-1)
    REAL(KIND=r8)   ,  INTENT(IN   ) :: upsoil  (npoi,nsoilay)   ! soil water uptake from transpiration (kg_h2o m-2 s-1)
    REAL(KIND=r8)   ,  INTENT(IN   ) :: csoi    (npoi,nsoilay)   ! specific heat of soil, no pore spaces (J kg-1 deg-1)
    REAL(KIND=r8)   ,  INTENT(IN   ) :: rhosoi  (npoi,nsoilay)   ! soil density (without pores, not bulk) (kg m-3)
    REAL(KIND=r8)   ,  INTENT(IN   ) :: suction (npoi,nsoilay)   ! saturated matric potential (m-h2o)
    REAL(KIND=r8)   ,  INTENT(IN   ) :: bex     (npoi,nsoilay)   ! exponent "b" in soil water potential
    INTEGER,  INTENT(IN   ) :: ibex    (npoi,nsoilay)   ! nint(bex), used for cpu speed
    REAL(KIND=r8)   ,  INTENT(IN   ) :: bperm  (npoi)                   ! lower b.c. for soil profile drainage 
    ! (0.0 = impermeable; 1.0 = fully permeable)
    REAL(KIND=r8)   , INTENT(IN   ) :: consoi  (npoi,nsoilay)   ! thermal conductivity of each soil layer (W m-1 K-1)
    REAL(KIND=r8)   , INTENT(INOUT) :: hflo    (npoi,nsoilay+1) ! downward heat transport through soil layers (W m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: tu(npoi)! temperature of upper canopy leaves (K)
    REAL(KIND=r8)   , INTENT(IN   ) :: tl(npoi)                 ! temperature of lower canopy leaves & stems(K)
    REAL(KIND=r8)   , INTENT(IN   ) :: raing (npoi)! rainfall rate at soil level (kg m-2 s-1)
    REAL(KIND=r8)   , INTENT(IN   ) :: fvapg (npoi)! h2o vapor flux (evaporation) between soil & air at z34 
    !(kg m-2 s-1/bare ground fraction)
    REAL(KIND=r8)   , INTENT(IN   ) :: traing(npoi)! rainfall temperature at soil level (K)
    REAL(KIND=r8)   , INTENT(IN) :: stressl (npoi,nsoilay)     ! soil moisture stress factor for the lower canopy (dimensionless)
    REAL(KIND=r8)   , INTENT(IN) :: stressu (npoi,nsoilay)     ! soil moisture stress factor for the upper canopy (dimensionless)
    REAL(KIND=r8)   , INTENT(IN) :: stresstl(npoi)             ! sum of stressl over all 6 soil layers (dimensionless)
    REAL(KIND=r8)   , INTENT(IN) :: stresstu(npoi)             ! sum of stressu over all 6 soil layers (dimensionless)

    !
    INTEGER :: i                ! loop indices 
    INTEGER :: k                ! loop indices
    !
    REAL(KIND=r8)    :: zfrez            ! factor decreasing runoff fraction for tsoi < tmelt
    REAL(KIND=r8)    :: zrunf            ! fraction of rain that doesn't stay in puddle (runoff fraction)
    REAL(KIND=r8)    :: rwork            ! 
    REAL(KIND=r8)    :: wipre            ! storing variable
    REAL(KIND=r8)    :: zdpud            ! used to compute transfer from puddle to infiltration
    REAL(KIND=r8)    :: cx               ! average specific heat for soil, water and ice
    REAL(KIND=r8)    :: chav             ! average specific heat for water and ice
    REAL(KIND=r8)    :: zwsoi      
    !
    REAL(KIND=r8)    :: owsoi  (npoi,nsoilay)    ! old value of wsoi
    REAL(KIND=r8)    :: otsoi  (npoi,nsoilay)    ! old value of tsoi
    REAL(KIND=r8)    :: c0pud  (npoi,nsoilay)    ! layer heat capacity due to puddles (=0 except for top)
    REAL(KIND=r8)    :: c1pud  (npoi,nsoilay)    ! updated av. specifilayer heat capacity due to  puddle
    REAL(KIND=r8)    :: wflo   (npoi,nsoilay+1)  ! = drainage at the bottom, returned by soilh2o
    REAL(KIND=r8)    :: fwtop  (npoi)            ! evaporation rate from soil (for soilh2o)
    REAL(KIND=r8)    :: fhtop  (npoi)            ! heat flux through soil surface (for soilheat)
    REAL(KIND=r8)    :: fwpud  (npoi)            ! portion of puddle that infiltrates in soil (rate)
    REAL(KIND=r8)    :: fsqueez(npoi)            ! excess amount of water (soilh2o) 
    REAL(KIND=r8)    :: dh     (npoi)            ! correction if water at tsoi < tmelt or ice at temp > tmelt
    REAL(KIND=r8)    :: dw     (npoi)            ! '
    REAL(KIND=r8)    :: zporos (npoi)            
    REAL(KIND=r8)    :: weigth (npoi) 

    !
    DO k=1,nsoilay
       DO i=1,npoi
          c0pud  (i,k) =  0.0_r8
          c1pud  (i,k) =  0.0_r8
       END DO
    END DO
    !      CALL const (c0pud, npoi*nsoilay, 0.0_r8)
    !      CALL const (c1pud, npoi*nsoilay, 0.0_r8)
    !
    ! for soil, set soil infiltration rate fwtop (for 
    ! soilh2o) and upper heat flux fhtop (for soilheat)
    !
    ! also step puddle model wpud, wipud
    !
    ! procedure is:
    !
    !   (0) immediately transfer any excess puddle liq to runoff
    !
    !   (1) apportion raing btwn puddle liquid(wpud) or runoff(grunof)
    !
    !   (2) apportion evap/condens (fvapg) btwn infil rate(fwtop), soil
    !       ice (wisoi(i,1)), puddle liq (wpud), or puddle ice (wipud)
    !
    !   (3) transfer some puddle liquid to fwtop
    !
    !   (4) compute upper heat flx fhtop: includes fwtop*ch2o*tsoi(i,1)
    !       to be consistent with whflo in soilheat, and accounts for
    !       changing rain temp from traing to tsoi(i,1) and runoff temp
    !       from tsoi to max(tsoi(i,1),tmelt)
    !
    DO i = 1, npoi
       !
       ! (0) immediately transfer any excess puddle liq to runoff
       !
       ! the following runoff formulation could give rise to very
       ! small amounts of negative runoff
       !
       grunof(i) = MIN (wpud(i), MAX (0.0_r8, wpud(i) + wipud(i) - &
            wpudmax)) / dtime
       !
       wpud(i) = wpud(i) - grunof(i) * dtime
       !
       ! (1) apportion sfc-level rain between puddle liquid and runoff
       !
       ! linear dependence of runoff fraction on wpud+wipud assumes
       ! uniform statistical distribution of sub-grid puddle 
       ! capacities between 0 and wpudmax. runoff fraction is 
       ! reduced linearly for tsoi < tmelt (by zfrez) in which case
       ! any rain will increase wpud which will be frozen to wipud
       ! below
       !
       zfrez = MAX (0.0_r8, MIN (1.0_r8, (tsoi(i,1) - tmelt + .5_r8) * 2.0_r8))
       !
       ! always have some minimal amount of runoff (0.10) even if 
       ! puddles are dry or soil is cold
       !
       zrunf = zfrez * MAX (0.0_r8, MIN (1.0_r8, (wpud(i) + wipud(i)) / &
            wpudmax))
       !
       wpud(i) = wpud(i) + (1.0_r8 - zrunf) * raing(i) * dtime
       !
       grunof(i) = grunof(i) + zrunf * raing(i)
       !
       ! (2) apportion evaporation or condensation between 4 h2o stores:
       !
       rwork = fvapg(i) * dtime
       !
       IF (fvapg(i).GE.0.0_r8) THEN
          !
          ! evaporation: split according to qglif
          !
          fwtop(i)   =            - qglif(i,1)*fvapg(i)
          wpud(i)    = wpud(i)    - qglif(i,3)*rwork
          wipud(i)   = wipud(i)   - qglif(i,4)*rwork
          !
          wipre = wisoi(i,1)
          wisoi(i,1) = MAX (0.0_r8, wipre - qglif(i,2)*rwork / &
               (rhow*poros(i,1)*hsoi(i,1)))
          !
          IF (1.0_r8-wisoi(i,1).GT.epsilon) &
               wsoi(i,1) = wsoi(i,1)*(1.0_r8-wipre)/(1.0_r8-wisoi(i,1))
          !
       ELSE
          !
          ! condensation: give all to puddles (to avoid wsoi, wisoi > 1)
          !
          fwtop(i) = 0.0_r8
          wpud(i) = wpud(i)  - (qglif(i,1)+qglif(i,3))*rwork
          wipud(i)= wipud(i) - (qglif(i,2)+qglif(i,4))*rwork
          !
       END IF
       !
       ! (3) transfer some puddle liquid to infiltration; can lead
       !     to small amounts of negative wpud (in soilh2o) due to
       !     round-off error
       !
       weigth (i) = (1.0_r8 - ((stressu(i,1  ) / max (stresstu(i), epsilon)) + (stressl(i,1  ) / max (stresstl(i), epsilon))))

       zdpud = rhow * dtime * MAX (0.0_r8, 1.0_r8-wisoi(i,1))**2 *  &
            hydraul(i,1)*weigth (i)
       !
       fwpud(i) = MAX (0.0_r8, MIN (wpud(i), zdpud)) / dtime
       c0pud(i,1) = ch2o*wpud(i) + cice*wipud(i)
       !
       ! (4) compute upper soil heat flux
       !
       fhtop(i) = heatg(i)   &
            + raing(i)*ch2o*(traing(i)-tsoi(i,1))  &
            - grunof(i)*ch2o*MAX(tmelt-tsoi(i,1), 0.0_r8)
       !
       ! update diagnostic variables
       !
       gadjust(i) = 0.0_r8
       !
    END DO
    !
    ! reduce soil moisture due to transpiration (upsoi[u,l], from
    ! turvap).need to do that before other time stepping below since 
    ! specific heat of this transport is neglected
    !
    ! first set porosflo, reduced porosity due to ice content, used
    ! as the effective porosity for uptake here and liquid hydraulics
    ! later in soilh2o. to avoid divide-by-zeros, use small epsilon
    ! limit; this will always cancel with epsilon or 0 in numerators
    !
    ! also increment soil temperature to balance transpired water
    ! differential between temps of soil and leaf. physically
    ! should apply this to the tree, but would be awkward in turvap.
    ! 
    ! also, save old soil moisture owsoi and temperatures otsoi so
    ! implicit soilh2o and soilheat can aposteriori deduce fluxes.
    !
    DO k = 1, nsoilay
       !
       DO i = 1, npoi
          !
          porosflo(i,k) = poros(i,k) * MAX (epsilon, (1.0_r8-wisoi(i,k)))
          !
          ! next line just for ice whose poros(i,k) is 0.0
          !
          porosflo(i,k) = MAX (porosflo(i,k), epsilon)
          !
          wsoi(i,k) = wsoi(i,k) - dtime * (upsoiu(i,k) + upsoil(i,k)) / &
               (rhow * porosflo(i,k) * hsoi(i,k))
          !
          cx = c0pud(i,k) +     &
               (   (1.0_r8-poros(i,k))*csoi(i,k)*rhosoi(i,k)  &
               + poros(i,k)*(1.0_r8-wisoi(i,k))*wsoi(i,k)*ch2o*rhow  &
               + poros(i,k)*wisoi(i,k)*cice*rhow  &
               ) * hsoi(i,k)
          !          WRITE(*,*)'poros(i,k)=',poros(i,k)
          !
          tsoi(i,k) = tsoi(i,k) - dtime * ch2o *     &
               (  upsoiu(i,k)*(tu(i)-tsoi(i,k))   &
               + upsoil(i,k)*(tl(i)-tsoi(i,k)) ) / cx
          !
          owsoi(i,k)  = wsoi(i,k)
          otsoi(i,k)  = tsoi(i,k)
          !
       END DO
       !
    END DO
    !
    ! step soil moisture calculations
    !
    CALL soilh2o (owsoi   , &! INTENT(IN   )
         fwtop   , &! INTENT(IN   )
         fwpud   , &! INTENT(IN   )
         fsqueez , &! INTENT(OUT  )
         wflo    , &! INTENT(OUT  )
         wisoi   , &! INTENT(IN   )
         hsoi    , &! INTENT(IN   )
         wsoi    , &! INTENT(INOUT)
         hydraul , &! INTENT(IN   )
         suction , &! INTENT(IN   )
         bex     , &! INTENT(IN   )
         ibex    , &! INTENT(IN   )
         bperm   , &! INTENT(IN   )
         porosflo, &! INTENT(IN   )
         poros   , &! INTENT(IN   )
         wpud    , &! INTENT(INOUT)
         npoi    , &! INTENT(IN   )
         nsoilay , &! INTENT(IN   )
         dtime   , &! INTENT(IN   )
         rhow    , &! INTENT(IN   )
         stressl , &! INTENT(IN   ) 
         stressu , &! INTENT(IN   ) 
         stresstl, &! INTENT(IN   ) 
         stresstu, &! INTENT(IN   ) 
         epsilon   )! INTENT(IN   )
    !
    ! update drainage and puddle
    !
    DO i = 1, npoi
       !
       gdrain(i)  = wflo(i,nsoilay+1)
       c1pud(i,1) = ch2o*wpud(i) + cice*wipud(i)
       !
    END DO
    !
    ! step temperatures due to conductive heat transport
    !
    CALL soilheat (otsoi  , &! INTENT(IN   )
         owsoi  , &! INTENT(IN   )
         c0pud  , &! INTENT(IN   )
         fhtop  , &! INTENT(IN   )
         wflo   , &! INTENT(IN   )
         c1pud  , &! INTENT(IN   )
         tsoi   , &! INTENT(INOUT)
         hsoi   , &! INTENT(IN   )
         consoi , &! INTENT(IN   )
         poros  , &! INTENT(IN   )
         csoi   , &! INTENT(IN   )
         rhosoi , &! INTENT(IN   )
         wisoi  , &! INTENT(IN   )
         wsoi   , &! INTENT(IN   )
         hflo   , &! INTENT(INOUT)
         npoi   , &! INTENT(IN   )
         nsoilay, &! INTENT(IN   )
         ch2o   , &! INTENT(IN   )
         rhow   , &! INTENT(IN   )
         cice   , &! INTENT(IN   )
         dtime    )! INTENT(IN   )
    !
    ! set wsoi, wisoi to exactly 0 or 1 if differ by negligible 
    ! amount (needed to avoid epsilon errors in loop 400 below)
    !
    CALL wadjust (hsoi   , &! hsoi   , &! INTENT(IN   )
         poros  , &! poros  , &! INTENT(IN   )
         wisoi  , &! wisoi  , &! INTENT(INOUT)
         wsoi   , &! wsoi   , &! INTENT(INOUT)
         gadjust, &! gadjust, &! INTENT(INOUT)
         npoi   , &! npoi   , &! INTENT(IN   )
         nsoilay, &! nsoilay, &! INTENT(IN   )
         rhow   , &! rhow   , &! INTENT(IN   )
         epsilon, &! epsilon, &! INTENT(IN   )
         dtime    )! dtime         )! INTENT(IN   )
    !
    ! heat-conserving adjustment for liquid/ice below/above melt
    ! point. uses exactly the same logic as for intercepted veg h2o
    ! in steph2o2. we employ the fiction here that soil liquid and
    ! soil ice both have density rhow, to avoid "pot-hole"
    ! difficulties of expansion on freezing. this is done by 
    ! dividing all eqns through by rhow(*hsoi).
    !
    ! the factor (1-wsoi(old))/(1-wisoi(new)) in the wsoi increments
    ! results simply from conservation of h2o mass; recall wsoi is
    ! liquid content relative to ice-reduced pore space.
    !
    DO k = 1, nsoilay
       DO i = 1, npoi
          !
          ! next line is just to avoid divide-by-zero for ice with
          ! poros = 0
          !
          zporos(i) = MAX (poros(i,k), epsilon)
          rwork = c1pud(i,k)/rhow/hsoi(i,k)                 &
               + (1.0_r8-zporos(i))*csoi(i,k)*rhosoi(i,k)/rhow
          !
          chav = rwork                                  &
               + zporos(i)*(1.0_r8-wisoi(i,k))*wsoi(i,k)*ch2o  &
               + zporos(i)*wisoi(i,k)*cice
          !
          ! if liquid exists below melt point, freeze some to ice
          !
          ! (note that if tsoi>tmelt or wsoi=0, nothing changes.)
          ! (also note if resulting wisoi=1, either dw=0 and prev
          ! wisoi=1, or prev wsoi=1, so use of epsilon is ok.)
          !
          zwsoi = MIN (1.0_r8, wsoi(i,k))
          !
          dh(i) = chav * (tmelt-tsoi(i,k))
          dw(i) = MIN ( zporos(i)*(1.0_r8-wisoi(i,k))*zwsoi, &
               MAX (0.0_r8,dh(i)/hfus) )
          !
          wisoi(i,k) = wisoi(i,k) +  dw(i)/zporos(i)
          wsoi(i,k)  = wsoi(i,k)  - (dw(i)/zporos(i))*(1.0_r8-zwsoi)  &
               / MAX (epsilon,1.0_r8-wisoi(i,k))
          !
          chav = rwork    &
               + zporos(i)*(1.0_r8-wisoi(i,k))*wsoi(i,k)*ch2o  &
               + zporos(i)*wisoi(i,k)*cice
          !
          !IF (chav == 0) THEN 
          !   tsoi(i,k) = tmelt -5.0_r8 
          !ELSE
          tsoi(i,k) = tmelt - (dh(i)-hfus*dw(i)) / chav
          !END IF 

          !
          ! if ice exists above melt point, melt some to liquid
          !
          ! note that if tsoi<tmelt or wisoi=0, nothing changes
          !
          ! also note if resulting wisoi=1, dw=0 and prev wisoi=1,
          ! so use of epsilon is ok
          !
          dh(i) = chav * (tsoi(i,k) - tmelt)
          dw(i) = MIN ( zporos(i)*wisoi(i,k), MAX (0.0_r8, dh(i)/hfus) )
          !
          wisoi(i,k) = wisoi(i,k) -  dw(i)/zporos(i)
          wsoi(i,k)  = wsoi(i,k)  + (dw(i)/zporos(i))  &
               * (1.0_r8-wsoi(i,k)) / MAX(epsilon,1.0_r8-wisoi(i,k))
          !
          chav = rwork   &
               + zporos(i)*(1.0_r8-wisoi(i,k))*wsoi(i,k)*ch2o   &
               + zporos(i)*wisoi(i,k)*cice
          !
          !IF (chav == 0) THEN 
          !   tsoi(i,k) = tmelt  +5.0_r8
          !ELSE
          tsoi(i,k) = tmelt + (dh(i)-hfus*dw(i)) / chav
          !END IF 
          !
          ! reset porosflo (although not used after this)
          !
          porosflo(i,k) = zporos(i) * MAX (epsilon, 1.0_r8-wisoi(i,k))
          !
       END DO
    END DO
    ! 
    ! set wsoi, wisoi to exactly 0 or 1 if differ by negligible 
    ! amount (roundoff error in loop 400 above can produce very
    ! small negative amounts)
    !
    CALL wadjust(hsoi    , &! hsoi   , &! INTENT(IN   )
         poros   , &! poros  , &! INTENT(IN   )
         wisoi   , &! wisoi  , &! INTENT(INOUT)
         wsoi    , &! wsoi   , &! INTENT(INOUT)
         gadjust , &! gadjust, &! INTENT(INOUT)
         npoi    , &! npoi   , &! INTENT(IN   )
         nsoilay , &! nsoilay, &! INTENT(IN   )
         rhow    , &! rhow   , &! INTENT(IN   )
         epsilon , &! epsilon, &! INTENT(IN   )
         dtime     )! dtime)! INTENT(IN   )
    !
    ! repeat ice/liquid adjustment for upper-layer puddles (don't 
    ! divide through by rhow*hsoi). upper-layer moistures wsoi,wisoi
    ! are already consistent with tsoi(i,1) > or < tmelt, and will 
    ! remain consistent here since tsoi(i,1) will not cross tmelt
    !
    k = 1
    !
    DO i = 1, npoi
       !
       ! if any puddle liquid below tmelt, freeze some to puddle ice
       !
       rwork = ( (1.0_r8-poros(i,k))*csoi(i,k)*rhosoi(i,k)      &
            + poros(i,k)*(1.0_r8-wisoi(i,k))*wsoi(i,k)*ch2o*rhow  &
            + poros(i,k)*wisoi(i,k)*cice*rhow  &
            ) * hsoi(i,k)
       !
       chav = ch2o*wpud(i) + cice*wipud(i) + rwork
       !
       dh(i) = chav * (tmelt-tsoi(i,k))
       dw(i) = MIN (wpud(i), MAX (0.0_r8, dh(i)/hfus))
       wipud(i) = wipud(i) + dw(i)
       wpud(i)  = wpud(i)  - dw(i)
       chav = ch2o*wpud(i) + cice*wipud(i) + rwork
       !IF (chav == 0) THEN 
       !   tsoi(i,k) = tmelt -5.0_r8 
       !ELSE
       tsoi(i,k) = tmelt - (dh(i)-hfus*dw(i)) / chav
       !END IF 

       !
       ! if any puddle ice above tmelt, melt some to puddle liquid
       !
       dh(i) = chav * (tsoi(i,k)-tmelt)
       dw(i) = MIN (wipud(i), MAX (0.0_r8, dh(i)/hfus))
       wipud(i) = wipud(i) - dw(i)
       wpud(i)  = wpud(i)  + dw(i)
       chav = ch2o*wpud(i) + cice*wipud(i) + rwork

       !IF (chav == 0) THEN 
       !   tsoi(i,k) = tmelt  + 5.0
       !ELSE
       tsoi(i,k) = tmelt + (dh(i)-hfus*dw(i)) / chav
       !END IF 

       !  
    END DO
    !
    RETURN
  END SUBROUTINE soilctl
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE soilh2o (owsoi   , &! INTENT(IN   )
       fwtop   , &! INTENT(IN   )
       fwpud   , &! INTENT(IN   )
       fsqueez , &! INTENT(OUT  )
       wflo    , &! INTENT(OUT  )
       wisoi   , &! INTENT(IN   )
       hsoi    , &! INTENT(IN   )
       wsoi    , &! INTENT(INOUT)
       hydraul , &! INTENT(IN   )
       suction , &! INTENT(IN   )
       bex     , &! INTENT(IN   )
       ibex    , &! INTENT(IN   )
       bperm   , &! INTENT(IN   )
       porosflo, &! INTENT(IN   )
       poros   , &! INTENT(IN   )
       wpud    , &! INTENT(INOUT)
       npoi    , &! INTENT(IN   )
       nsoilay , &! INTENT(IN   )
       dtime   , &! INTENT(IN   )
       rhow    , &! INTENT(IN   )
       stressl , &! INTENT(IN   ) 
       stressu , &! INTENT(IN   ) 
       stresstl, &! INTENT(IN   ) 
       stresstu, &! INTENT(IN   ) 
       epsilon   )! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    ! sets up call to tridia to solve implicit soil moisture eqn,
    ! using soil temperatures in wsoi (in comsoi)
    !
    ! lower bc can be no h2o flow or free drainage, set by bperm below
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN   ) :: npoi     ! total number of land points
    INTEGER, INTENT(IN   ) :: nsoilay  ! number of soil layers
    REAL(KIND=r8)   , INTENT(IN   ) :: dtime    ! model timestep (seconds)
    REAL(KIND=r8)   , INTENT(IN   ) :: rhow     ! density of liquid water (all types) (kg m-3)
    REAL(KIND=r8)   , INTENT(IN   ) :: epsilon  ! small quantity to avoid zero-divides and other
    ! truncation or machine-limit troubles with small
    ! values. should be slightly greater than o(1)
    ! machine precision
    REAL(KIND=r8)   , INTENT(IN   ) :: wisoi   (npoi,nsoilay)   ! fraction of soil pore space containing ice
    REAL(KIND=r8)   , INTENT(IN   ) :: hsoi    (npoi,nsoilay+1)      ! soil layer thickness (m)
    REAL(KIND=r8)   , INTENT(INOUT) :: wsoi    (npoi,nsoilay)   ! fraction of soil pore space containing liquid water
    REAL(KIND=r8)   , INTENT(IN   ) :: hydraul (npoi,nsoilay)   ! saturated hydraulic conductivity (m/s)
    REAL(KIND=r8)   , INTENT(IN   ) :: suction (npoi,nsoilay)   ! saturated matric potential (m-h2o)
    REAL(KIND=r8)   , INTENT(IN   ) :: bex     (npoi,nsoilay)   ! exponent "b" in soil water potential
    INTEGER, INTENT(IN   ) :: ibex    (npoi,nsoilay)   ! nint(bex), used for cpu speed
    REAL(KIND=r8)   , INTENT(IN   ) :: bperm (npoi)                    ! lower b.c. for soil profile drainage 
    ! (0.0 = impermeable; 1.0 = fully permeable)      wisoi   
    REAL(KIND=r8)   , INTENT(IN   ) :: porosflo(npoi,nsoilay)   ! porosity after reduction by ice content
    REAL(KIND=r8)   , INTENT(IN   ) :: poros   (npoi,nsoilay)   ! porosity (mass of h2o per unit vol at sat / rhow)
    REAL(KIND=r8)   , INTENT(INOUT) :: wpud    (npoi)           ! liquid content of puddles per soil area (kg m-2)
    !
    ! Arguments : all are supplied except wflo (returned):
    !
    REAL(KIND=r8)   , INTENT(IN   ) :: owsoi(npoi,nsoilay)  ! soil moistures at start of timestep
    REAL(KIND=r8)   , INTENT(IN   ) :: fwtop(npoi)          ! evaporation from top soil layer
    REAL(KIND=r8)   , INTENT(IN   ) :: fwpud(npoi)          ! h2o flux into top layer (infiltration from puddle)
    REAL(KIND=r8)   , INTENT(OUT  ) :: fsqueez(npoi)        ! excess water at end of time step in soil column  
    REAL(KIND=r8)   , INTENT(OUT  ) :: wflo(npoi,nsoilay+1) ! downward h2o flow across layer boundaries
    REAL(KIND=r8)   , INTENT(IN   ) :: stressl (npoi,nsoilay)     ! soil moisture stress factor for the lower canopy (dimensionless)
    REAL(KIND=r8)   , INTENT(IN   ) :: stressu (npoi,nsoilay)     ! soil moisture stress factor for the upper canopy (dimensionless)
    REAL(KIND=r8)   , INTENT(IN   ) :: stresstl(npoi)             ! sum of stressl over all 6 soil layers (dimensionless)
    REAL(KIND=r8)   , INTENT(IN   ) :: stresstu(npoi)             ! sum of stressu over all 6 soil layers (dimensionless)

    !
    ! local variables
    !     
    INTEGER :: k    ! loop indices
    INTEGER :: i    ! loop indices
    INTEGER :: km1  ! loop indices
    INTEGER :: kka  ! loop indices
    INTEGER :: kkb  ! loop indices
    !
    INTEGER :: m(npoi)                ! exponents 
    INTEGER :: n(npoi)
    !
    REAL(KIND=r8)   , PARAMETER :: wmin=0.10_r8 ! minimum  soils moisture 
    REAL(KIND=r8)   , PARAMETER :: wmax=1.00_r8 ! maximum  soils moisture 

    REAL(KIND=r8)   , PARAMETER :: dmin= 1.e-9_r8 ! minimum diffusivity for dry soils (m**2 s-1) 
    REAL(KIND=r8)   , PARAMETER :: rimp= 1.0_r8   ! implicit fraction of the calculation (0 to 1)
    !      REAL(KIND=r8)    ::bperm,                 ! = 0 for impermeable (no drainage) lower bc,
    ! = 1 for free drainage lower bc
    REAL(KIND=r8)    :: zbex
    REAL(KIND=r8)    :: z
    REAL(KIND=r8)    :: dt
    REAL(KIND=r8)    :: zz
    ! 
    REAL(KIND=r8)    :: hsoim(npoi,nsoilay+1)  ! vertical distances between centers of layers
    !
    REAL(KIND=r8)    :: wsoim(npoi,nsoilay+1)    ! interpolated moisture values at layer boundaries
    REAL(KIND=r8)    :: wsoia(npoi,nsoilay+1)    ! 
    REAL(KIND=r8)    :: wsoib(npoi,nsoilay+1)    ! '
    REAL(KIND=r8)    :: weim (npoi,nsoilay+1)    ! '
    REAL(KIND=r8)    :: weip (npoi,nsoilay+1)    ! '
    REAL(KIND=r8)    :: a    (npoi)              ! intermediate terms (const for each pt)
    REAL(KIND=r8)    :: b    (npoi)              ! 
    REAL(KIND=r8)    :: bwn  (npoi)              ! 
    REAL(KIND=r8)    :: bwn1 (npoi)              ! 
    REAL(KIND=r8)    :: e    (npoi,nsoilay+1)    ! intermediate terms in algebraic devel 
    REAL(KIND=r8)    :: f    (npoi,nsoilay+1)    ! '
    REAL(KIND=r8)    :: g    (npoi,nsoilay+1)    ! '
    REAL(KIND=r8)    :: d1   (npoi,nsoilay)      ! diagonals of tridiagonal systems of equations 
    REAL(KIND=r8)    :: d2   (npoi,nsoilay)      !  '
    REAL(KIND=r8)    :: d3   (npoi,nsoilay)      !  '
    REAL(KIND=r8)    :: rhs  (npoi,nsoilay)      ! right-hand sides of systems of equations
    REAL(KIND=r8)    :: w1   (npoi,nsoilay)      ! work arrays needed by tridia
    REAL(KIND=r8)    :: w2   (npoi,nsoilay)      !  '
    REAL(KIND=r8)    :: weigth   (npoi)      !  '
    REAL(KIND=r8)    :: weigthp   (npoi)      !  '

    !
    ! set lower boundary condition for the soil
    ! (permeability of the base)
    !
    !     bperm = 0.00  ! e.g. fully impermeable base
    !     bperm = 1.00  ! e.g. fully permeable base
    !
    !      bperm = 0.10
    !
    ! set level vertical distances, interpolated moistures, and
    ! interpolation weights
    !
    ! top layer
    !
    k = 1
    !
    DO i = 1, npoi
       !
       hsoim(i,k) = 0.5_r8 * hsoi(i,k)
       !
       weim(i,k) = 0.0_r8
       weip(i,k) = 1.0_r8
       !
       wsoim(i,k) = MAX(MIN (wsoi(i,k) , wmax),wmin)
       wsoia(i,k) = MAX(MIN (wsoim(i,k), wmax),wmin)
       wsoib(i,k) = MAX(MIN (wsoim(i,k), wmax),wmin)
       !
    END DO
    !
    ! middle layers
    !
    DO k = 2, nsoilay
       !
       DO i = 1, npoi
          !
          hsoim(i,k) = 0.5_r8 * (hsoi(i,k-1) + hsoi(i,k))
          !
          !   --                                        --     --                                       --
          !  |  1        hsoi                             |   |     1       hsoi                          |
          ! =| ---- * ------------------------------------| + |1 - ---- * --------------------------------|
          !  |  2        0.5_r8 * (hsoi(k-1) + hsoi(k))   |   |     2      0.5_r8 * (hsoi(k-1) + hsoi(k)) |
          !   --                                        --     --                                       --
          !
          weim(i,k) = 0.5_r8 * hsoi(i,k) / hsoim(i,k)
          weip(i,k) = 1.0_r8 - weim(i,k)
          !
          wsoim(i,k) = weim(i,k) * wsoi(i,k-1) + weip(i,k) * wsoi(i,k)
          wsoia(i,k) = MAX(MIN (wsoim(i,k), wmax),wmin)
          wsoib(i,k) = MAX(MIN (wsoim(i,k), wmax),wmin)
          !
       END DO
       !
    END DO
    !
    ! bottom layer
    !
    k = nsoilay + 1
    !
    DO i = 1, npoi
       !
       hsoim(i,k) = 0.5_r8 * hsoi(i,k-1)
       !
       weim(i,k) = 1.0_r8
       weip(i,k) = 0.0_r8
       !
       wsoim(i,k) = MAX(MIN (wsoi(i,k-1), wmax),wmin)
       wsoia(i,k) = MAX(MIN (wsoim(i,k), wmax),wmin)
       wsoib(i,k) = MAX(MIN (wsoim(i,k), wmax),wmin)
       !
    END DO
    !
    ! set intermediate quantities e,f,g. these are terms in the
    ! expressions for the fluxes at boundaries between layers,
    ! so are zero for k=1. use bwn1 to account for minimum 
    ! diffusivity dmin. bperm is used for k=nsoilay+1 to set the
    ! type of the lower bc.
    !
    ! top layer
    !
    k = 1
    !
    DO i = 1, npoi
       e   (i,k) =0.0_r8
    END DO
    !      CALL const (e   (1,k),& !INTENT(OUT  )
    !                  npoi     ,& !INTENT(IN   )
    !          0.0_r8       ) !INTENT(IN   )
    DO i = 1, npoi
       f   (i,k) =0.0_r8
    END DO
    !      CALL const (f(1,k)   ,& !INTENT(OUT  )
    !                  npoi     ,& !INTENT(IN   )
    !          0.0_r8)              !INTENT(IN   )

    DO i = 1, npoi
       g   (i,k) =0.0_r8
    END DO

    !      CALL const (g(1,k)   ,& !INTENT(OUT  )
    !                  npoi     ,& !INTENT(IN   )
    !          0.0_r8       ) !INTENT(IN   )
    !
    !
    !                      -  - 2b +3               --            --
    !         dWl       d |    |                d  |     b+2    dWl |
    !       P---- = -K0---| Wl |    + K0*fhi0* ----| B Wl    * ---- |
    !         dt        dz|    |                dz |            dz  |
    !                      -  -                     --            --
    !
    !                      -  - 2b +3               --            --
    !         dWl       d |    |                d  |     b+2    dWl |
    !       P---- = -K0---| Wl |    + K0*fhi0* ----| B Wl    * ---- | + 
    !         dt        dz|    |                dz |            dz  |
    !                      -  -                     --            --
    !
    !               -    - 2b +3                     --            --
    !            d | dWl  |                   d  dWl|     b+2    dWl |
    !        -K0---| ---- |    + K0*fhi0*B * -------|   Wl    * ---- |
    !            dz| dWl  |                   dz dWl|            dz  |
    !               -    -                           --            --

    !
    !
    !                            -  - 2b+2                         --                                     --
    !         dWl               |    |     dWl                 d  |       b+2       dwl                 b+2 |
    !       P---- = -K0*(2b+3)* | Wl |   *------ + K0*fhi0*B* ----|(b+1) Wl      * ----  + K0*fhi0* B*Wl    | +
    !         dt                |    |     dz                  dz |                 dz                      |
    !                            -  -                              --                                     --
    !               -    - 2b +3                     --            --
    !            d | dWl  |                   d  dWl|     b+2    dWl |
    !        -K0---| ---- |    + K0*fhi0*B * -------|   Wl    * ---- |
    !            dz| dWl  |                   dz dWl|            dz  |
    !               -    -                           --            --

    !                            -  - 2b+2                           b+2                         b+2
    !         dWl               |    |     dWl                     dWl      dwl                dWl  
    !       P---- = -K0*(2b+3)* | Wl |   *------ + K0*fhi0*B* (b+1)----  * ----  + K0*fhi0* B* ----  +
    !         dt                |    |     dz                       dz      dz                  dz  
    !                            -  -                                                               
    !               -    - 2b +3                     --            --
    !            d | dWl  |                   d  dWl|     b+2    dWl |
    !        -K0---| ---- |    + K0*fhi0*B * -------|   Wl    * ---- |
    !            dz| dWl  |                   dz dWl|            dz  |
    !               -    -                           --            --

    !
    ! middle layers
    !
    DO  k = 2, nsoilay
       !
       DO  i = 1, npoi
          !
          ! now that hydraul, suction and ibex can vary with depth,
          ! use averages of surrounding mid-layer values
          !
          ! (see notes 8/27/93)
          !   stressu(i,k) / max (stresstu(i), epsilon)
          weigthp(i) = (1.0_r8 - ((stressu(i,k-1) / max (stresstu(i), epsilon)) + (stressl(i,k-1) / max (stresstl(i), epsilon))))
          weigth (i) = (1.0_r8 - ((stressu(i,k  ) / max (stresstu(i), epsilon)) + (stressl(i,k  ) / max (stresstl(i), epsilon))))
          !
          !  hsoi = DZ
          !
          !
          !        1       DZ                            1        DZ
          !term1= --- * -------------- * K(k-1) + ( 1 - ---- * -------------- ) *K(k)
          !        2      0.5_r8 * DZ                    2      0.5_r8 * DZ
          !  
          a(i) = weim(i,k) * hydraul(i,k-1)*weigthp(i) + &
                 weip(i,k) * hydraul(i,k  )*weigth(i) 
          !
          !FHYsat = suction  ! saturated matric potential (m-h2o)
          !          K0*Fhy0*B
          
          b(i) = weim(i,k) * hydraul(i,k-1)*weigthp(i) * &
                 suction(i,k-1) * bex(i,k-1) +  &
                 weip(i,k) * hydraul(i,k  )*weigth(i)  * &
                 suction(i,k  ) * bex(i,k  )
          !
          zbex = weim(i,k) * bex(i,k-1) + &
                 weip(i,k) * bex(i,k  ) 
          !
          !  http://www.ecmwf.int/research/ifsdocs/CY25r1/Physics/Physics-08-07.html
          !
          m(i) = 2 * NINT(zbex) + 3
          n(i) =     NINT(zbex) + 2
          !         O  
          ! wsoib =---
          !         Os 
          !                       B+1
          !       =  K0*Fhy0*B* Wl
          !
          bwn1(i) = b(i) * (wsoib(i,k)**(n(i)-1))
          !                       (B+1 )
          !       =  K0*Fhy0*B* Wl           * Wl

          bwn(i)  = bwn1(i) * wsoib(i,k)
          !
          IF (bwn(i).LT.dmin) bwn1(i) = 0.0_r8
          bwn(i) = MAX (bwn(i), dmin)
          !                                              --                                                            --
          !                                  -  - 2b +3 |                   -- -- b+1                                    |
          !         dWl                     |    |      |                  |     |                            b+1        |   DWl 
          !       P---- = (-1 + (2b+3))*K0* | Wl |    + |(1 - 1)*K0*fhi0*B*| Wl  |  * Wl  - (b+2)*K0*fhi0* B*Wl    * Wl  | *-----
          !         dt                      |    |      |                  |     |                                       |    Dz
          !                                  -  -       |                   -- --                                        |
          !                                              --                                                            --
          !
          !                                         2b+3  
          !          ( 2b+3  -1)        *   K0  *  Wl  

          e(i,k) =  (-1.0_r8 + rimp*m(i))*a(i)*(wsoia(i,k)**m(i))      &
                 + (( 1.0_r8-rimp)*bwn(i)         - rimp*n(i)*bwn1(i)*wsoib(i,k))  &
               * (wsoi(i,k)-wsoi(i,k-1)) / hsoim(i,k)
          !
          !
          !
          !
          !                            2B+3 - 1                              B+1 - 1
          !               (2B+3)*K0*Wl               + b+2*  +  K0*Fhy0*B* Wl
          f(i,k) = - rimp*m(i)*a(i)*(wsoia(i,k)**(m(i)-1)) &
                   + rimp*n(i)*bwn1(i) & 
                   * (wsoi(i,k)-wsoi(i,k-1)) / hsoim(i,k)
          !
          g(i,k) = rimp*bwn(i)
          !
       END DO
       !
    END DO
    !                      -  - 2b +3               --       --               
    !         dWl       d |    |                d  |     b+2   |   dwl                 b+2    dWl 
    !       P---- = -K0---| Wl |    + K0*fhi0* ----| B*Wl      |* ----  + K0*fhi0* B*Wl    * ------ 
    !         dt        dz|    |                dz |           |   dz                         dzdz  
    !                      -  -                     --       --         

    !
    ! bottom layer
    !
    k = nsoilay + 1
    !
    DO i = 1, npoi
       !
       weigthp(i) = (1.0_r8 - ((stressu(i,nsoilay) / max (stresstu(i), epsilon)) + (stressl(i,nsoilay) / max (stresstl(i), epsilon))))
       weigth (i) = (1.0_r8 - ((stressu(i,nsoilay) / max (stresstu(i), epsilon)) + (stressl(i,nsoilay) / max (stresstl(i), epsilon))))
       !
       a(i) = hydraul(i,nsoilay) * weigthp(i)
       b(i) = hydraul(i,nsoilay) * weigth (i)*suction(i,nsoilay)*ibex(i,nsoilay)
       !
       m(i) = 2*ibex(i,nsoilay) + 3
       n(i) = ibex(i,nsoilay)   + 2
       !
       e(i,k) =                -a(i)*(wsoia(i,k)**m(i))*bperm(i)
       f(i,k) = 0.0_r8
       g(i,k) = 0.0_r8
       !
    END DO
    !
    ! deduce all e,f,g in proportion to the minimum of the two 
    ! adjacent layers' (1-wisoi), to account for restriction of flow
    ! by soil ice. this will cancel in loop 300  with the factor 
    ! 1-wisoi in (one of) the layer's porosflo, even if wisoi=1 by 
    ! the use of epsilon limit. so a layer with wisoi=1 will form a 
    ! barrier to flow of liquid, but still have a predicted wsoi
    !
    DO k = 1, nsoilay+1
       !
       kka = MAX (k-1,1)
       kkb = MIN (k,nsoilay)
       !
       DO i=1,npoi
          !
          ! multiply by an additional factor of 1-wisoi for stability
          !
          z = MAX(0.0_r8,1.0_r8-MAX(wisoi(i,kka),wisoi(i,kkb)))**2
          !
          e(i,k) = z * e(i,k)
          f(i,k) = z * f(i,k)
          g(i,k) = z * g(i,k)
          !
       END DO
       !
    END DO
    !
    ! set matrix diagonals and right-hand sides
    !
    DO k = 1, nsoilay
       !
       DO  i = 1, npoi
          !
          dt = dtime / (porosflo(i,k)*hsoi(i,k))
          d1(i,k) = dt*(   f(i,k)*0.5_r8*hsoi(i,k)/hsoim(i,k) &
                         - g(i,k)/hsoim(i,k) )
          rhs(i,k) = wsoi(i,k) + dt*( e(i,k+1) - e(i,k) )
          !
       END DO
       !
       IF (k.EQ.1) THEN
          !
          DO i=1,npoi
             !
             dt = dtime / (porosflo(i,k)*hsoi(i,k))
             rhs(i,k) = rhs(i,k) + dt*(fwtop(i)+fwpud(i))/rhow
             !
          END DO
          !
       END IF
       !
       IF (k.LT.nsoilay) THEN
          !
          km1 = MAX (k-1,1)
          !
          DO i=1,npoi
             !
             dt = dtime / (porosflo(i,k)*hsoi(i,k))
             d2(i,k) = 1.0_r8 + dt*( - f(i,k+1)*0.5_r8*hsoi(i,k+1)/hsoim(i,k+1)  &
                  + f(i,k)  *0.5_r8*hsoi(i,km1)/hsoim(i,k)    &
                  + g(i,k+1)/hsoim(i,k+1)                & 
                  + g(i,k)  /hsoim(i,k) )
             d3(i,k) = dt*( - f(i,k+1)*0.5_r8*hsoi(i,k)/hsoim(i,k+1)         &
                  - g(i,k+1)            /hsoim(i,k+1) )
             !
          END DO
          !
       ELSE IF (k.EQ.nsoilay) THEN
          !
          DO i=1,npoi
             !
             dt = dtime / (porosflo(i,k)*hsoi(i,k))
             d2(i,k) = 1.0_r8 + dt*( - f(i,k+1)                         &
                  + f(i,k)  *0.5_r8*hsoi(i,k-1)/hsoim(i,k)  &
                  + g(i,k)  /hsoim(i,k) )
             d3(i,k) = 0.0_r8
             !
          END DO
          !
       END IF
       !
    END DO
    ! General information:
    
    !     Subroutine tridia solves a tridiagonal system of equations for the
    !     n+1 unknowns u(0), u(1), ..., u(n-1), u(n).  This system must
    !     consist of n+1 equations of the form:
    !
    !                   b(0) u(0) + c(0) u(1)   = d(0)   first equation
    !     a(i) u(i-1) + b(i) u(i) + c(i) u(i+1) = d(i)   for i = 1, 2, ..., n-1
    !     a(n) u(n-1) + b(n) u(n)               = d(n)   final equation

    !     Subroutine tridia uses the Thomas algorithm. This method does not
    !     use pivoting; it may crash or be inaccurate if the tridiagonal
    !     system is not strictly diagonally dominant or symmetric positive
    !     definite.

    !     Always check the error code ier after using tridia; if ier is 2 or
    !     more, there may be a problem.

    !     Copyright 1996 Leon van Dommelen
    !     Version 1.0 Leon van Dommelen 1/2/97

    ! Computational procedure:

    !     Subroutine tridia works by reducing the given system to a
    !     bidiagonal system of equations of the form
    !
    !                        u(0) + k(0) u(1)   = l(0)   first equation
    !                        u(i) + k(i) u(i+1) = l(i)   for i = 1, 2, ..., n-1
    !                        u(n)               = l(n)   final equation
    !
    !     and then solving this system in inverse order.
    !
    ! solve the systems of equations
    !
    CALL tridia (npoi    , &! INTENT(IN   ) :: ns ! number of systems to be solved.
         npoi    , &! INTENT(IN   ) :: nd ! first dimension of arrays (ge ns)
         nsoilay , &! INTENT(IN   ) :: ne ! number of unknowns in each system. (>2)
         d1      , &! INTENT(IN   ) :: a(nd,ne)     ! subdiagonals of matrices stored in a(j,2)...a(j,ne).
         d2      , &! INTENT(IN   ) :: b(nd,ne)     ! main diagonals of matrices stored in b(j,1)...b(j,ne).
         d3      , &! INTENT(IN   ) :: c(nd,ne)     ! super-diagonals of matrices stored in c(j,1)...c(j,ne-1).
         rhs     , &! INTENT(IN   ) :: y(nd,ne)     ! right hand side of equations stored in y(j,1)...y(j,ne).
         wsoi    , &! INTENT(INOUT) :: x(nd,ne)     ! solutions of the systems returned in x(j,1)...x(j,ne).
         w1      , &! INTENT(INOUT) :: alpha(nd,ne) ! work array 
         w2        )! INTENT(INOUT) :: gamma(nd,ne) ! work array
    !
    DO i = 1, npoi
       !
       fsqueez(i) = 0.0_r8
       wflo(i,nsoilay+1) = - rhow * e(i,nsoilay+1)
       !
    END DO
    !
    DO k = nsoilay, 1, -1
       !
       DO i = 1, npoi
          !
          zz = rhow * poros(i,k) *    &
               MAX(epsilon, (1.0_r8-wisoi(i,k))) * hsoi(i,k)
          !
          wsoi(i,k) = wsoi(i,k) + dtime * fsqueez(i) / zz 
          fsqueez(i) = MAX (wsoi(i,k)-1.0_r8,0.0_r8) * zz / dtime
          wsoi(i,k) = MIN (wsoi(i,k),wmax)
          wsoi(i,k) = MAX (wsoi(i,k),wmin)
          !
          wflo(i,k) = wflo(i,k+1) + (wsoi(i,k)-owsoi(i,k)) * zz / dtime
          !
       END DO
       !
    END DO
    !
    ! step puddle liquid due to fsqueez and fwpud
    !
    ! also subtract net puddle-to-top-layer flux from wflo(i,1),
    ! since puddle and top soil layer are lumped together in soilheat
    ! so upper wflo should be external flux only (evap/condens)

    DO i = 1, npoi
       !
       wpud(i)   = wpud(i)   + (fsqueez(i) - fwpud(i)) * dtime
       wflo(i,1) = wflo(i,1) + (fsqueez(i) - fwpud(i))
       !
    END DO
    !
    RETURN
  END SUBROUTINE soilh2o
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE soilheat (otsoi  , &! INTENT(IN   )
       owsoi  , &! INTENT(IN   )
       c0pud  , &! INTENT(IN   )
       fhtop  , &! INTENT(IN   )
       wflo   , &! INTENT(IN   )
       c1pud  , &! INTENT(IN   )
       tsoi   , &! INTENT(INOUT)
       hsoi   , &! INTENT(IN   )
       consoi , &! INTENT(IN   )
       poros  , &! INTENT(IN   )
       csoi   , &! INTENT(IN   )
       rhosoi , &! INTENT(IN   )
       wisoi  , &! INTENT(IN   )
       wsoi   , &! INTENT(IN   )
       hflo   , &! INTENT(INOUT)
       npoi   , &! INTENT(IN   )
       nsoilay, &! INTENT(IN   )
       ch2o   , &! INTENT(IN   )
       rhow   , &! INTENT(IN   )
       cice   , &! INTENT(IN   )
       dtime    )! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    !        sets up call to tridia to solve implicit soil/ice heat 
    !        conduction, using layer temperatures in tsoi (in comsoi).
    !        the heat flux due to liquid flow previously calculated
    !        in soilh2o is accounted for. lower bc is conductive flux = 0
    !        for soil (although the flux due to liquid drainage flow can
    !        be > 0)
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN   ) :: npoi                       ! total number of land points
    INTEGER, INTENT(IN   ) :: nsoilay                    ! number of soil layers
    REAL(KIND=r8)   , INTENT(IN   ) :: ch2o                       ! specific heat of liquid water (J deg-1 kg-1)
    REAL(KIND=r8)   , INTENT(IN   ) :: rhow                       ! density of liquid water (all types) (kg m-3)
    REAL(KIND=r8)   , INTENT(IN   ) :: cice                       ! specific heat of ice (J deg-1 kg-1)
    REAL(KIND=r8)   , INTENT(IN   ) :: dtime                      ! model timestep (seconds)
    REAL(KIND=r8)   , INTENT(INOUT) :: tsoi   (npoi,nsoilay)      ! soil temperature for each layer (K)
    REAL(KIND=r8)   , INTENT(IN   ) :: hsoi   (npoi,nsoilay+1)    ! soil layer thickness (m)
    REAL(KIND=r8)   , INTENT(IN   ) :: consoi (npoi,nsoilay)      ! thermal conductivity of each soil layer (W m-1 K-1)
    REAL(KIND=r8)   , INTENT(IN   ) :: poros  (npoi,nsoilay)      ! porosity (mass of h2o per unit vol at sat / rhow)
    REAL(KIND=r8)   , INTENT(IN   ) :: csoi   (npoi,nsoilay)      ! specific heat of soil, no pore spaces (J kg-1 deg-1)
    REAL(KIND=r8)   , INTENT(IN   ) :: rhosoi (npoi,nsoilay)      ! soil density (without pores, not bulk) (kg m-3)
    REAL(KIND=r8)   , INTENT(IN   ) :: wisoi  (npoi,nsoilay)      ! fraction of soil pore space containing ice
    REAL(KIND=r8)   , INTENT(IN   ) :: wsoi   (npoi,nsoilay)      ! fraction of soil pore space containing liquid water
    REAL(KIND=r8)   , INTENT(INOUT) :: hflo   (npoi,nsoilay+1)    ! downward heat transport through soil layers (W m-2)
    !
    ! Arguments
    !
    REAL(KIND=r8)   , INTENT(IN   ) :: otsoi(npoi,nsoilay)     ! soil/ice temperatures at start of timestep (redundant
    ! with tsoi, but passed to be consistent with soilh2o)
    REAL(KIND=r8)   , INTENT(IN   ) :: owsoi(npoi,nsoilay)     ! soil moistures at start of timestep (before soilh2o)
    REAL(KIND=r8)   , INTENT(IN   ) :: c0pud(npoi,nsoilay)     ! layer heat capacity due to puddles (=0 except for top)
    REAL(KIND=r8)   , INTENT(IN   ) :: c1pud(npoi,nsoilay)     ! updated c0pud
    REAL(KIND=r8)   , INTENT(IN   ) :: fhtop(npoi)             ! heat flux into top layer from atmos
    REAL(KIND=r8)   , INTENT(IN   ) :: wflo (npoi,nsoilay+1)   ! downward h2o  flow across layer boundaries
    !
    ! local variables
    !
    INTEGER :: k                       ! loop indices
    INTEGER :: i                       ! loop indices
    INTEGER :: km1                     ! loop indices
    INTEGER :: kp1                     ! loop indices
    !
    REAL(KIND=r8)   , PARAMETER :: rimp= 1.0_r8    ! implicit fraction of the calculation (0 to 1)
    ! = 0 for impermeable (no drainage) lower bc,
    ! = 1 for free drainage lower bc
    REAL(KIND=r8)    :: rwork                   ! work variables
    REAL(KIND=r8)    :: rwork1                  ! work variables
    REAL(KIND=r8)    :: rwork2                  ! work variables
    REAL(KIND=r8)    :: rwork3                  ! work variables
    REAL(KIND=r8)    :: t
    !
    REAL(KIND=r8)    :: whflo(npoi,nsoilay+1)   ! downward heat fluxes across layer bdries due to h2o
    ! movement calculated in soilh2o
    REAL(KIND=r8)    :: con  (npoi,nsoilay+1)   ! conduction coeffs between layers
    REAL(KIND=r8)    :: c0   (npoi,nsoilay)     ! specific heats at start of timestep
    REAL(KIND=r8)    :: c1   (npoi,nsoilay)     ! specific heats at end of timestep
    REAL(KIND=r8)    :: d1   (npoi,nsoilay)     ! diagonals of tridiagonal systems of equations
    REAL(KIND=r8)    :: d2   (npoi,nsoilay)     !  ''
    REAL(KIND=r8)    :: d3   (npoi,nsoilay)     !  ''
    REAL(KIND=r8)    :: rhs  (npoi,nsoilay)     ! right-hand sides of systems of equations
    REAL(KIND=r8)    :: w1   (npoi,nsoilay)     ! work arrays needed by tridia
    REAL(KIND=r8)    :: w2   (npoi,nsoilay)
    !
    ! set conduction coefficient between layers, and heat fluxes
    ! due to liquid transport
    !
    ! top layer
    !
    k = 1
    !
    DO i = 1, npoi
       con(i,k) = 0.0_r8
       whflo(i,k) = wflo(i,k) * ch2o * tsoi(i,k)
    END DO
    !
    ! middle layers
    !
    DO k = 2, nsoilay
       !
       DO i = 1, npoi
          !
          con(i,k) =  1.0_r8 / (0.5_r8 * (hsoi(i,k-1) / consoi(i,k-1) + &
               hsoi(i,k)   / consoi(i,k)))
          !
          t = (hsoi(i,k) * tsoi(i,k-1) + hsoi(i,k-1) * tsoi(i,k)) /  &
               (hsoi(i,k-1)             + hsoi(i,k))
          !
          whflo(i,k) = wflo(i,k) * ch2o * t
          !
       END DO
       !
    END DO
    !
    ! bottom layer
    !
    k = nsoilay + 1
    !
    DO i = 1, npoi
       con(i,k) = 0.0_r8
       whflo(i,k) = wflo(i,k) * ch2o * tsoi(i,k-1)
    END DO
    !
    ! set diagonals of matrix and right-hand side. use old and
    ! new heat capacities c0, c1 consistently with moisture fluxes
    ! whflo computed above, to conserve heat associated with 
    ! changing h2o amounts in each layer
    !
    DO k = 1, nsoilay
       !
       km1 = MAX (k-1,1)
       kp1 = MIN (k+1,nsoilay)
       !
       DO i=1,npoi
          !
          rwork1 = (1.0_r8-poros(i,k))*csoi(i,k)*rhosoi(i,k)
          rwork2 = poros(i,k)*(1.0_r8-wisoi(i,k))*ch2o*rhow
          rwork3 = poros(i,k)*wisoi(i,k)*cice*rhow
          !
          c0(i,k) = c0pud(i,k) +              &
               (   rwork1             &
               + rwork2 * owsoi(i,k)   &
               + rwork3             &
               ) * hsoi(i,k)
          !
          c1(i,k) = c1pud(i,k) +             &
               (   rwork1           &
               + rwork2 * wsoi(i,k)  &
               + rwork3           &
               ) * hsoi(i,k)
          !
          rwork = dtime/c1(i,k)
          !
          d1(i,k) =    - rwork * rimp * con(i,k)
          d2(i,k) = 1.0_r8 + rwork * rimp * (con(i,k)+con(i,k+1))
          d3(i,k) =    - rwork * rimp * con(i,k+1)
          !
          rhs(i,k) = (c0(i,k)/c1(i,k))*tsoi(i,k) + rwork               &
               * ( (1.0_r8-rimp)*con(i,k)  *(tsoi(i,km1)-tsoi(i,k))  &
               + (1.0_r8-rimp)*con(i,k+1)*(tsoi(i,kp1)-tsoi(i,k))  &
               + whflo(i,k) - whflo(i,k+1) )
!          rhs(i,k) = (c0(i,k)/c1(i,k))*tsoi(i,k) + DT/pcDZ               &
!                     * ( (1.0_r8-rimp)*kg/dz  *(tso i(i,km1)-tsoi(i,k)) + (1.0_r8-rimp)*kg/dz*(tsoi(i,kp1)-tsoi(i,k)) + whflo(i,k) - whflo(i,k+1) )

          !
       END DO
       !
       IF (k.EQ.1) THEN
          DO i=1,npoi
             rhs(i,k) = rhs(i,k) + (dtime/c1(i,k))*fhtop(i)
          END DO
       END IF
       !
    END DO
!                     -    -
!   d(p*c*T)    d kg |  dT  |    d (pl*cl*wl*Tl)
!   -------   = -----| ---- | + ------
!   dt          dz   |  dz  |    dz
!                     -    -

!                           -    -
!   (p*c*T)       dt*kg   |  d2T  |    dt        d (pl*cl*wl*Tl)
!   ------- =     --------| ----- | + -----    * ------
!                  p*c    |  d2z  |    p*c          dz
!                          -    -

!
! solve systems of equations
!
!  -         -     -  -      - -
! |b  c  0  0 |   |  y |    | x |
! |a  b  c  0 | * |  y |  = | x |
! |0  a  b  c |   |  y |    | x |
! |0  0  a  b |   |  Y |    | x |
!  -         -     -  -      -  -

!  x(1) = y(1)*b + y(2)*c + y(3)*0 + y(4)*0

    CALL tridia (npoi     , &! INTENT(IN   ) :: ns ! number of systems to be solved.
         npoi     , &! INTENT(IN   ) :: nd ! first dimension of arrays (ge ns)
         nsoilay  , &! INTENT(IN   ) :: ne ! number of unknowns in each system. (>2)
         d1       , &! INTENT(IN   ) :: a(nd,ne)     ! subdiagonals of matrices stored in a(j,2)...a(j,ne).
         d2       , &! INTENT(IN   ) :: b(nd,ne)     ! main diagonals of matrices stored in b(j,1)...b(j,ne).
         d3       , &! INTENT(IN   ) :: c(nd,ne)     ! super-diagonals of matrices stored in c(j,1)...c(j,ne-1).
         rhs      , &! INTENT(IN   ) :: y(nd,ne)     ! right hand side of equations stored in y(j,1)...y(j,ne).
         tsoi     , &! INTENT(INOUT) :: x(nd,ne)     ! solutions of the systems returned in x(j,1)...x(j,ne).
         w1       , &! INTENT(INOUT) :: alpha(nd,ne) ! work array 
         w2         )! INTENT(INOUT) :: gamma(nd,ne) ! work array
    !
    ! deduce downward heat fluxes between layers
    !
    DO i=1,npoi
       hflo(i,1) = fhtop(i)
    END DO

    !      CALL scopy (npoi     , & ! INTENT(IN   )
    !                  fhtop    , & ! INTENT(IN   )
    !          hflo(1,1)  ) ! INTENT(OUT  )
    !
    DO k=1,nsoilay
       DO i=1,npoi
          hflo(i,k+1) = hflo(i,k) - &
               (c1(i,k)*tsoi(i,k) - c0(i,k)*otsoi(i,k)) / dtime
       END DO
    END DO
    !
    RETURN
  END SUBROUTINE soilheat
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE wadjust(hsoi   , &! INTENT(IN   )
       poros  , &! INTENT(IN   )
       wisoi  , &! INTENT(INOUT)
       wsoi   , &! INTENT(INOUT)
       gadjust, &! INTENT(INOUT)
       npoi   , &! INTENT(IN   )
       nsoilay, &! INTENT(IN   )
       rhow   , &! INTENT(IN   )
       epsilon, &! INTENT(IN   )
       dtime    )! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    ! set wsoi, wisoi to exactly 0 if differ by negligible amount, 
    ! to protect epsilon logic in soilctl and soilh2o
    !
    ! ice-liquid transformations in soilctl loop 400 can produce very
    ! small -ve amounts due to roundoff error, and very small -ve or +ve
    ! amounts can cause (harmless) "underflow" fpes in soilh2o
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN   ) :: npoi     ! total number of land points
    INTEGER, INTENT(IN   ) :: nsoilay  ! number of soil layers
    REAL(KIND=r8)   , INTENT(IN   ) :: rhow     ! density of liquid water (all types) (kg m-3)
    REAL(KIND=r8)   , INTENT(IN   ) :: epsilon  ! small quantity to avoid zero-divides and other
    ! truncation or machine-limit troubles with small
    ! values. should be slightly greater than o(1)
    ! machine precision
    REAL(KIND=r8)   , INTENT(IN   ) :: dtime              ! model timestep (seconds)
    REAL(KIND=r8)   , INTENT(INOUT) :: gadjust(npoi)      ! h2o flux due to adjustments in 
    ! subroutine wadjust (kg_h2o m-2 s-1)
    REAL(KIND=r8)   , INTENT(IN   ) :: hsoi (npoi,nsoilay+1)! soil layer thickness (m)
    REAL(KIND=r8)   , INTENT(IN   ) :: poros(npoi,nsoilay)! porosity (mass of h2o per unit vol at sat / rhow)
    REAL(KIND=r8)   , INTENT(INOUT) :: wisoi(npoi,nsoilay)! fraction of soil pore space containing ice
    REAL(KIND=r8)   , INTENT(INOUT) :: wsoi (npoi,nsoilay)! fraction of soil pore space containing liquid water
    ! 
    ! local variables
    !
    INTEGER :: k
    INTEGER :: i
    REAL(KIND=r8)    :: ztot0 
    REAL(KIND=r8)    :: ztot1
    !
    DO k = 1, nsoilay
       DO i = 1, npoi
          !
          ! initial total soil water
          !
          ztot0 = hsoi(i,k) * poros(i,k) * rhow * &
               ((1.0_r8 - wisoi(i,k)) * wsoi(i,k) + wisoi(i,k))
          !
          ! set bounds on wsoi and wisoi
          !
          IF (wsoi(i,k).LT.epsilon)  wsoi(i,k)  = 0.0_r8
          IF (wisoi(i,k).LT.epsilon) wisoi(i,k) = 0.0_r8
          !
          wsoi(i,k)  = MIN (1.0_r8, wsoi(i,k))
          wisoi(i,k) = MIN (1.0_r8, wisoi(i,k))
          !
          IF (wisoi(i,k).GE.1-epsilon) wsoi(i,k) = 0.0_r8
          !
          ! for diagnosis of total adjustment
          !
          ztot1 = hsoi(i,k) * poros(i,k) * rhow * &
               ((1.0_r8 - wisoi(i,k)) * wsoi(i,k) + wisoi(i,k))
          !
          gadjust(i) = gadjust(i) + (ztot1 - ztot0) / dtime
          !
       END DO
    END DO
    !
    RETURN
  END SUBROUTINE wadjust


  !####### #     # ######     #####  #######   ###   #
  !#         ##    # #     #   #         # #         #    #    #
  !#         # #   # #     #   #           #         #    #    #
  !#####   #  #  # #     #    #####  #         #    #    #
  !#         #   # # #     #            # #         #    #    #
  !#         #    ## #     #   #         # #         #    #    #
  !####### #     # ######     #####  #######   ###   #######











  !
  !  ####     ##    #    #   ####   #####    #   #
  ! #    #   #  #   ##   #  #    #  #    #    # #
  ! #       #    #  # #  #  #    #  #    #     #
  ! #       ######  #  # #  #    #  #####      #
  ! #    #  #    #  #   ##  #    #  #          #
  !  ####   #    #  #    #   ####   #          #
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE canopy(bps         , &! INTENT(INOUT) !local
       rhoa         , &! INTENT(INOUT) !local
       cp           , &! INTENT(INOUT) !local
       za           , &! INTENT(INOUT) !local
       bdl          , &! INTENT(INOUT) !local
       dil          , &! INTENT(INOUT) !local
       z3           , &! INTENT(INOUT) !local
       z4           , &! INTENT(INOUT) !local
       z34          , &! INTENT(INOUT) !local
       exphl        , &! INTENT(INOUT) !local
       expl         , &! INTENT(INOUT) !local
       displ        , &! INTENT(OUT  ) !local
       bdu          , &! INTENT(INOUT) !local
       diu          , &! INTENT(INOUT) !local
       z1           , &! INTENT(INOUT) !local
       z2           , &! INTENT(INOUT) !local
       z12          , &! INTENT(INOUT) !local
       exphu        , &! INTENT(INOUT) !local
       expu         , &! INTENT(INOUT) !local
       dispu        , &! INTENT(OUT  ) !local
       alogg        , &! INTENT(OUT  ) !local
       alogi        , &! INTENT(OUT  ) !local
       alogav       , &! INTENT(INOUT) !local
       alog4        , &! INTENT(INOUT) !local
       alog3        , &! INTENT(INOUT) !local
       alog2        , &! INTENT(OUT  ) !local
       alog1        , &! INTENT(INOUT) !local
       aloga        , &! INTENT(INOUT) !local
       u2           , &! INTENT(INOUT) !local
       alogu        , &! INTENT(INOUT) !local
       alogl        , &! INTENT(INOUT) !local
       richl        , &! INTENT(OUT  ) !local
       straml       , &! INTENT(OUT  ) !local
       strahl       , &! INTENT(OUT  ) !local
       richu        , &! INTENT(OUT  ) !local
       stramu       , &! INTENT(OUT  ) !local
       strahu       , &! INTENT(OUT  ) !local
       u1           , &! INTENT(OUT  ) !local
       u12          , &! INTENT(OUT  ) !local
       u3           , &! INTENT(OUT  ) !local
       u34          , &! INTENT(OUT  ) !local
       u4           , &! INTENT(OUT  ) !local
       cu           , &! INTENT(INOUT) !local
       cl           , &! INTENT(INOUT) !local
       sg           , &! INTENT(INOUT) !local
       si           , &! INTENT(INOUT) !local
       fwetl        , &! INTENT(IN   ) !global
       scalcoefl    , &! INTENT(IN   ) !global
       scalcoefu    , &! INTENT(IN   ) !global
       fwetu        , &! INTENT(IN   ) !global
       termu        , &! INTENT(IN   ) !global
       fwetux       , &! INTENT(OUT  ) !local
       fwetsx       , &! INTENT(OUT  ) !local
       fwets        , &! INTENT(IN   ) !global
       fwetlx       , &! INTENT(OUT  ) !local
       solu             , &! INTENT(IN   ) !global
       firu             , &! INTENT(IN   ) !global
       sols             , &! INTENT(IN   ) !global
       firs         , &! INTENT(IN   ) !global
       soll         , &! INTENT(IN   ) !global
       firl         , &! INTENT(IN   ) !global
       rliqu        , &! INTENT(IN   ) !global
       rliqs        , &! INTENT(IN   ) !global
       rliql        , &! INTENT(IN   ) !global
       pfluxu       , &! INTENT(IN   ) !global
       pfluxs       , &! INTENT(IN   ) !global
       pfluxl       , &! INTENT(IN   ) !global
       solg         , &! INTENT(IN   ) !global
       firg         , &! INTENT(INOUT) !global
       soli             , &! INTENT(IN   ) !global
       firi         , &! INTENT(INOUT) !global
       fsena        , &! INTENT(OUT  ) !local
       fseng        , &! INTENT(OUT  ) !local
       fseni        , &! INTENT(OUT  ) !local
       fsenu        , &! INTENT(OUT  ) !local
       fsens        , &! INTENT(OUT  ) !local
       fsenl        , &! INTENT(OUT  ) !local
       fvapa        , &! INTENT(OUT  ) !local
       fvaput       , &! INTENT(OUT  ) !local
       fvaps        , &! INTENT(OUT  ) !local
       fvaplw       , &! INTENT(OUT  ) !local
       fvaplt       , &! INTENT(OUT  ) !local
       fvapg        , &! INTENT(OUT  ) !local
       fvapi        , &! INTENT(OUT  ) !local
       firb         , &! INTENT(INOUT) !global
       terml        , &! INTENT(IN   ) !global
       fvapuw       , &! INTENT(OUT  ) !local
       ztop         , &! INTENT(IN   ) !global
       fl           , &! INTENT(IN   ) !global
       lai          , &! INTENT(IN   ) !global
       sai              , &! INTENT(IN   ) !global
       alaiml       , &! INTENT(IN   ) !global
       zbot         , &! INTENT(IN   ) !global
       fu           , &! INTENT(IN   ) !global
       alaimu       , &! INTENT(IN   ) !global
       froot        , &! INTENT(IN   ) !global
       t34          , &! INTENT(INOUT) !global
       t12          , &! INTENT(INOUT) !global
       q34          , &! INTENT(INOUT) !global
       q12          , &! INTENT(INOUT) !global
       su           , &! INTENT(INOUT) !local
       cleaf        , &! INTENT(IN   ) !global
       dleaf        , &! INTENT(IN   ) !global
       ss           , &! INTENT(INOUT) !local
       cstem        , &! INTENT(IN   ) !global
       dstem        , &! INTENT(IN   ) !global
       sl           , &! INTENT(INOUT) !local
       cgrass       , &! INTENT(IN   ) !global
       tu             , &! INTENT(INOUT) !global
       ciub         , &! INTENT(INOUT) !global
       ciuc         , &! INTENT(INOUT) !global
       exist        , &! INTENT(IN   ) !global
       topparu      , &! INTENT(IN   ) !global
       csub         , &! INTENT(INOUT) !global
       gsub         , &! INTENT(INOUT) !global
       csuc         , &! INTENT(INOUT) !global
       gsuc         , &! INTENT(INOUT) !global
       agcub        , &! INTENT(OUT  ) !local
       agcuc        , &! INTENT(OUT  ) !local
       ancub        , &! INTENT(OUT  ) !local
       ancuc        , &! INTENT(OUT  ) !local
       totcondub    , &! INTENT(INOUT) !local
       totconduc    , &! INTENT(INOUT) !local
       tl           , &! INTENT(INOUT) !global
       cils         , &! INTENT(INOUT) !global
       cil3         , &! INTENT(INOUT) !global
       cil4         , &! INTENT(INOUT) !global
       topparl      , &! INTENT(IN   ) !global
       csls         , &! INTENT(INOUT) !global
       gsls         , &! INTENT(INOUT) !global
       csl3         , &! INTENT(INOUT) !global
       gsl3         , &! INTENT(INOUT) !global
       csl4         , &! INTENT(INOUT) !global
       gsl4         , &! INTENT(INOUT) !global
       agcls        , &! INTENT(OUT  ) !local
       agcl4        , &! INTENT(OUT  ) !local
       agcl3        , &! INTENT(OUT  ) !local
       ancls        , &! INTENT(OUT  ) !local
       ancl4        , &! INTENT(OUT  ) !local
       ancl3        , &! INTENT(OUT  ) !local
       totcondls    , &! INTENT(INOUT) !local
       totcondl3    , &! INTENT(INOUT) !local
       totcondl4    , &! INTENT(INOUT) !local
       chu          , &! INTENT(IN   ) !global
       wliqu        , &! INTENT(IN   ) !global
       wsnou        , &! INTENT(IN   ) !global
       chs          , &! INTENT(IN   ) !global
       wliqs        , &! INTENT(IN   ) !global
       wsnos        , &! INTENT(IN   ) !global
       chl          , &! INTENT(IN   ) !global
       wliql        , &! INTENT(IN   ) !global
       wsnol        , &! INTENT(IN   ) !global
       ts           , &! INTENT(INOUT) !global
       frac         , &! INTENT(IN   ) !global
       z0soi        , &! INTENT(IN   ) !global
       wsoi         , &! INTENT(IN   ) !global
       wisoi        , &! INTENT(IN   ) !global
       swilt        , &! INTENT(IN   ) !global
       sfield       , &! INTENT(IN   ) !global
       stressl      , &! INTENT(INOUT) !local
       stressu      , &! INTENT(INOUT) !local
       stresstl     , &! INTENT(INOUT) !local
       stresstu     , &! INTENT(INOUT) !local
       poros        , &! INTENT(IN   ) !global
       wpud         , &! INTENT(IN   ) !global
       wipud        , &! INTENT(IN   ) !global
       csoi         , &! INTENT(IN   ) !global
       rhosoi       , &! INTENT(IN   ) !global
       hsoi         , &! INTENT(IN   ) !global
       consoi       , &! INTENT(IN   ) !global
       tg           , &! INTENT(INOUT) !global
       ti           , &! INTENT(INOUT) !global
       wpudmax      , &! INTENT(IN   ) !global
       suction      , &! INTENT(IN   ) !global
       bex          , &! INTENT(IN   ) !global
       hvasug       , &! INTENT(IN   ) !global
       tsoi         , &! INTENT(IN   ) !global
       hvasui       , &! INTENT(IN   ) !global
       upsoiu       , &! INTENT(OUT  ) !local
       upsoil       , &! INTENT(OUT  ) !local
       fi           , &! INTENT(IN   ) !global
       z0sno        , &! INTENT(IN   ) !global
       consno       , &! INTENT(IN   ) !global
       hsno         , &! INTENT(IN   ) !global
       hsnotop      , &! INTENT(IN   ) !global
       tsno         , &! INTENT(IN   ) !global
       psurf        , &! INTENT(IN   ) !global
       ta           , &! INTENT(IN   ) !global
       qa           , &! INTENT(IN   ) !global
       ua           , &! INTENT(IN   ) !global
       o2conc       , &! INTENT(IN   ) !global
       co2conc      , &! INTENT(IN   ) !global
       npoi         , &! INTENT(IN   ) !global
       nsoilay      , &! INTENT(IN   ) !global
       nsnolay      , &! INTENT(IN   ) !global
       npft         , &! INTENT(IN   ) !global
       vonk         , &! INTENT(IN   ) !global
       epsilon      , &! INTENT(IN   ) !global
       hvap         , &! INTENT(IN   ) !global
       ch2o         , &! INTENT(IN   ) !global
       hsub         , &! INTENT(IN   ) !global
       cice         , &! INTENT(IN   ) !global
       rhow         , &! INTENT(IN   ) !global
       stef         , &! INTENT(IN   ) !global
       tmelt        , &! INTENT(IN   ) !global
       hfus         , &! INTENT(IN   ) !global
       cappa        , &! INTENT(IN   ) !global
       rair         , &! INTENT(IN   ) !global
       rvap         , &! INTENT(IN   ) !global
       cair         , &! INTENT(IN   ) !global
       cvap         , &! INTENT(IN   ) !global
       grav         , &! INTENT(IN   ) !global
       dtime        , &! INTENT(IN   ) !global
       vmax_pft     , &! INTENT(IN   ) !global
       tau15        , &! INTENT(IN   ) !global
       kc15         , &! INTENT(IN   ) !global
       ko15         , &! INTENT(IN   ) !global
       cimax        , &! INTENT(IN   ) !global
       gammaub      , &! INTENT(IN   ) !global
       alpha3       , &! INTENT(IN   ) !global
       theta3       , &! INTENT(IN   ) !global
       beta3        , &! INTENT(IN   ) !global
       coefmub      , &! INTENT(IN   ) !global
       coefbub      , &! INTENT(IN   ) !global
       gsubmin      , &! INTENT(IN   ) !global
       gammauc      , &! INTENT(IN   ) !global
       coefmuc      , &! INTENT(IN   ) !global
       coefbuc      , &! INTENT(IN   ) !global
       gsucmin      , &! INTENT(IN   ) !global
       gammals      , &! INTENT(IN   ) !global
       coefmls      , &! INTENT(IN   ) !global
       coefbls      , &! INTENT(IN   ) !global
       gslsmin      , &! INTENT(IN   ) !global
       gammal3      , &! INTENT(IN   ) !global
       coefml3      , &! INTENT(IN   ) !global
       coefbl3      , &! INTENT(IN   ) !global
       gsl3min      , &! INTENT(IN   ) !global
       gammal4      , &! INTENT(IN   ) !global
       alpha4       , &! INTENT(IN   ) !global
       theta4       , &! INTENT(IN   ) !global
       beta4        , &! INTENT(IN   ) !global
       coefml4      , &! INTENT(IN   ) !global
       coefbl4      , &! INTENT(IN   ) !global
       gsl4min      , &! INTENT(IN   ) !global
       a10scalparamu, &! INTENT(INOUT) !global
       a10daylightu , &! INTENT(INOUT) !global
       a10scalparaml, &! INTENT(INOUT) !global
       a10daylightl , &! INTENT(INOUT) !global
       ginvap       , &! INTENT(OUT  ) !local
       gsuvap       , &! INTENT(OUT  ) !local
       gtrans       , &! INTENT(OUT  ) !local
       gtransu      , &! INTENT(OUT  ) !local
       gtransl      , &! INTENT(OUT  ) !local
       ux           , &! INTENT(IN   ) !global
       uy           , &! INTENT(IN   ) !global
       taux         , &! INTENT(OUT  ) !local
       tauy         , &! INTENT(OUT  ) !local
       ts2          , &! INTENT(OUT  ) !local
       qs2          , &! INTENT(OUT  ) !local
       vegtype0     , &!
       stressfac    , &
       nVegClass      )




    ! ---------------------------------------------------------------------
    !
    ! calculates sensible heat and moisture flux coefficients,
    ! and steps canopy temperatures through one timestep
    !
    ! atmospheric conditions at za are supplied in comatm
    ! arrays ta, qa, psurf and scalar siga (p/ps)
    !
    ! downward sensible heat and moisture fluxes at za
    ! are returned in com1d arrays fsena, fvapa
    !
    ! sensible heat and moisture fluxes from solid objects to air
    ! are stored (for other models and budget) in com1d arrays
    ! fsen[u,s,l,g,i], fvap[u,s,l,g,i]
    !
    ! the procedure is first to compute wind speeds and aerodynamic
    ! transfer coefficients in turcof, then call turvap to solve an
    ! implicit linear system for temperatures and specific
    ! humidities and the corresponding fluxes - this is iterated
    ! niter times for non-linearities due to stratification,
    ! implicit/explicit (h2o phase), dew, vpd and max soil
    ! moisture uptake - t12 and q12 are changed each iteration,
    ! and tu, ts, tl, tg, ti can be adjusted too
    !
    ! initialize aerodynamic quantities
    !
    IMPLICIT NONE

    INTEGER, INTENT(IN   ) :: nVegClass


    INTEGER, INTENT(IN   ) :: npoi           ! total number of land points
    INTEGER, INTENT(IN   ) :: nsoilay        ! number of soil layers
    INTEGER, INTENT(IN   ) :: nsnolay        ! number of snow layers
    INTEGER, INTENT(IN   ) :: npft           ! number of plant functional types


    REAL(KIND=r8), INTENT(IN   ) :: vonk           ! von karman constant (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: epsilon        ! small quantity to avoid zero-divides and other
                                                   ! truncation or machine-limit troubles with small
                                                   ! values. should be slightly greater than o(1)
                                                   ! machine precision
    REAL(KIND=r8), INTENT(IN   ) :: hvap           ! latent heat of vaporization of water (J kg-1)
    REAL(KIND=r8), INTENT(IN   ) :: ch2o           ! specific heat of liquid water (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN   ) :: hsub           ! latent heat of sublimation of ice (J kg-1)
    REAL(KIND=r8), INTENT(IN   ) :: cice           ! specific heat of ice (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN   ) :: rhow           ! density of liquid water (all types) (kg m-3)
    REAL(KIND=r8), INTENT(IN   ) :: stef           ! stefan-boltzmann constant (W m-2 K-4)
    REAL(KIND=r8), INTENT(IN   ) :: tmelt          ! freezing point of water (K)
    REAL(KIND=r8), INTENT(IN   ) :: hfus           ! latent heat of fusion of water (J kg-1)
    REAL(KIND=r8), INTENT(IN   ) :: cappa               ! rair/cair
    REAL(KIND=r8), INTENT(IN   ) :: rair           ! gas constant for dry air (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN   ) :: rvap           ! gas constant for water vapor (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN   ) :: cair           ! specific heat of dry air at constant pressure (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN   ) :: cvap           ! specific heat of water vapor at constant pressure (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN   ) :: grav           ! gravitational acceleration (m s-2)
    REAL(KIND=r8), INTENT(IN   ) :: dtime          ! model timestep (seconds)
    REAL(KIND=r8), INTENT(IN   ) :: vmax_pft(npft) ! nominal vmax of top leaf at 15 C (mol-co2/m**2/s) [not used]
    REAL(KIND=r8), INTENT(IN   ) :: tau15          ! co2/o2 specificity ratio at 15 degrees C (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: kc15           ! co2 kinetic parameter (mol/mol)
    REAL(KIND=r8), INTENT(IN   ) :: ko15           ! o2 kinetic parameter (mol/mol) 
    REAL(KIND=r8), INTENT(IN   ) :: cimax          ! maximum value for ci (needed for model stability)
    REAL(KIND=r8), INTENT(IN   ) :: gammaub        ! leaf respiration coefficient
    REAL(KIND=r8), INTENT(IN   ) :: alpha3         ! intrinsic quantum efficiency for C3 plants (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: theta3         ! photosynthesis coupling coefficient for C3 plants (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: beta3          ! photosynthesis coupling coefficient for C3 plants (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: coefmub        ! 'm' coefficient for stomatal conductance relationship
    REAL(KIND=r8), INTENT(IN   ) :: coefbub        ! 'b' coefficient for stomatal conductance relationship
    REAL(KIND=r8), INTENT(IN   ) :: gsubmin        ! absolute minimum stomatal conductance
    REAL(KIND=r8), INTENT(IN   ) :: gammauc        ! leaf respiration coefficient
    REAL(KIND=r8), INTENT(IN   ) :: coefmuc        ! 'm' coefficient for stomatal conductance relationship  
    REAL(KIND=r8), INTENT(IN   ) :: coefbuc        ! 'b' coefficient for stomatal conductance relationship  
    REAL(KIND=r8), INTENT(IN   ) :: gsucmin        ! absolute minimum stomatal conductance
    REAL(KIND=r8), INTENT(IN   ) :: gammals        ! leaf respiration coefficient
    REAL(KIND=r8), INTENT(IN   ) :: coefmls        ! 'm' coefficient for stomatal conductance relationship 
    REAL(KIND=r8), INTENT(IN   ) :: coefbls        ! 'b' coefficient for stomatal conductance relationship 
    REAL(KIND=r8), INTENT(IN   ) :: gslsmin        ! absolute minimum stomatal conductance
    REAL(KIND=r8), INTENT(IN   ) :: gammal3        ! leaf respiration coefficient
    REAL(KIND=r8), INTENT(IN   ) :: coefml3        ! 'm' coefficient for stomatal conductance relationship 
    REAL(KIND=r8), INTENT(IN   ) :: coefbl3        ! 'b' coefficient for stomatal conductance relationship 
    REAL(KIND=r8), INTENT(IN   ) :: gsl3min        ! absolute minimum stomatal conductance
    REAL(KIND=r8), INTENT(IN   ) :: gammal4        ! leaf respiration coefficient
    REAL(KIND=r8), INTENT(IN   ) :: alpha4         ! intrinsic quantum efficiency for C4 plants (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: theta4         ! photosynthesis coupling coefficient for C4 plants (dimensionless) 
    REAL(KIND=r8), INTENT(IN   ) :: beta4          ! photosynthesis coupling coefficient for C4 plants (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: coefml4        ! 'm' coefficient for stomatal conductance relationship
    REAL(KIND=r8), INTENT(IN   ) :: coefbl4        ! 'b' coefficient for stomatal conductance relationship
    REAL(KIND=r8), INTENT(IN   ) :: gsl4min        ! absolute minimum stomatal conductance
    REAL(KIND=r8), INTENT(INOUT) :: a10scalparamu (npoi) ! 10-day average day-time scaling parameter - upper canopy (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: a10daylightu  (npoi) ! 10-day average day-time PAR - upper canopy (micro-Ein m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: a10scalparaml (npoi) ! 10-day average day-time scaling parameter - lower canopy (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: a10daylightl  (npoi) ! 10-day average day-time PAR - lower canopy (micro-Ein m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: ginvap (npoi)  ! total evaporation rate from all intercepted h2o (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  )  :: gsuvap (npoi)  ! total evaporation rate from surface (snow/soil) (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  )  :: gtrans (npoi)  ! total transpiration rate from all vegetation canopies (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  )  :: gtransu(npoi)  ! transpiration from upper canopy (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  )  :: gtransl(npoi)  ! transpiration from lower canopy (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(IN   ) :: psurf  (npoi)         ! surface pressure (Pa)  &
    REAL(KIND=r8), INTENT(IN   ) :: ta     (npoi)         ! air temperature (K)  &
    REAL(KIND=r8), INTENT(IN   ) :: qa     (npoi)         ! specific humidity (kg_h2o/kg_air)  &
    REAL(KIND=r8), INTENT(IN   ) :: ua     (npoi)         ! wind speed (m s-1)  &
    REAL(KIND=r8), INTENT(IN   ) :: o2conc           ! o2 concentration (mol/mol)
    REAL(KIND=r8), INTENT(IN   ) :: co2conc(npoi)          ! co2 concentration (mol/mol)  &
    REAL(KIND=r8), INTENT(IN   ) :: fi     (npoi)! fractional snow cover
    REAL(KIND=r8), INTENT(IN   ) :: z0sno! roughness length of snow surface (m)
    REAL(KIND=r8), INTENT(IN   ) :: consno  ! thermal conductivity of snow (W m-1 K-1)
    REAL(KIND=r8), INTENT(IN   ) :: hsno   (npoi,nsnolay)! thickness of snow layers (m)
    REAL(KIND=r8), INTENT(IN   ) :: hsnotop ! thickness of top snow layer (m)
    REAL(KIND=r8), INTENT(IN   ) :: tsno   (npoi,nsnolay)! temperature of snow layers (K)
    REAL(KIND=r8), INTENT(IN   ) :: z0soi(npoi)             ! roughness length of soil surface (m)
    REAL(KIND=r8), INTENT(IN   ) :: wsoi(npoi,nsoilay)     ! fraction of soil pore space containing liquid water
    REAL(KIND=r8), INTENT(IN   ) :: wisoi(npoi,nsoilay)     ! fraction of soil pore space containing ice
    REAL(KIND=r8), INTENT(IN   ) :: swilt(npoi,nsoilay)     ! wilting soil moisture value (fraction of pore space)
    REAL(KIND=r8), INTENT(IN   ) :: sfield  (npoi,nsoilay)     ! field capacity soil moisture value (fraction of pore space)
    REAL(KIND=r8), INTENT(INOUT) :: stressl (npoi,nsoilay)     ! soil moisture stress factor for the lower canopy (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: stressu (npoi,nsoilay)     ! soil moisture stress factor for the upper canopy (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: stresstl(npoi)             ! sum of stressl over all 6 soil layers (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: stresstu(npoi)             ! sum of stressu over all 6 soil layers (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: poros(npoi,nsoilay)     ! porosity (mass of h2o per unit vol at sat / rhow)
    REAL(KIND=r8), INTENT(IN   ) :: wpud(npoi)             ! liquid content of puddles per soil area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wipud(npoi)             ! ice content of puddles per soil area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: csoi(npoi,nsoilay)     ! specific heat of soil, no pore spaces (J kg-1 deg-1)
    REAL(KIND=r8), INTENT(IN   ) :: rhosoi  (npoi,nsoilay)     ! soil density (without pores, not bulk) (kg m-3)
    REAL(KIND=r8), INTENT(IN   ) :: hsoi    (npoi,nsoilay+1)   ! soil layer thickness (m)
    REAL(KIND=r8), INTENT(IN   ) :: consoi  (npoi,nsoilay)     ! thermal conductivity of each soil layer (W m-1 K-1)
    REAL(KIND=r8), INTENT(INOUT) :: tg(npoi)             ! soil skin temperature (K)
    REAL(KIND=r8), INTENT(INOUT) :: ti(npoi)             ! snow skin temperature (K)
    REAL(KIND=r8), INTENT(IN   ) :: wpudmax            ! normalization constant for puddles (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: suction (npoi,nsoilay)     ! saturated matric potential (m-h2o)
    REAL(KIND=r8), INTENT(IN   ) :: bex(npoi,nsoilay)     ! exponent "b" in soil water potential
    REAL(KIND=r8), INTENT(IN   ) :: hvasug  (npoi)             ! latent heat of vap/subl, for soil surface (J kg-1)
    REAL(KIND=r8), INTENT(IN   ) :: tsoi(npoi,nsoilay)     ! soil temperature for each layer (K)
    REAL(KIND=r8), INTENT(IN   ) :: hvasui  (npoi)             ! latent heat of vap/subl, for snow surface (J kg-1)
    REAL(KIND=r8), INTENT(OUT  ) :: upsoiu  (npoi,nsoilay)     ! soil water uptake from transpiration (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: upsoil  (npoi,nsoilay)     ! soil water uptake from transpiration (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(IN   ) :: ztop         (npoi,2)    ! height of plant top above ground (m)
    REAL(KIND=r8), INTENT(IN   ) :: fl         (npoi)      ! fraction of snow-free area covered by lower  canopy
    REAL(KIND=r8), INTENT(IN   ) :: lai         (npoi,2)    ! canopy single-sided leaf area index (area leaf/area veg)
    REAL(KIND=r8), INTENT(IN   ) :: sai         (npoi,2)    ! current single-sided stem area index
    REAL(KIND=r8), INTENT(IN   ) :: alaiml               ! lower canopy leaf & stem maximum area (2 sided) for 
    ! normalization of drag coefficient (m2 m-2)
    REAL(KIND=r8), INTENT(IN   ) :: zbot         (npoi,2)    ! height of lowest branches above ground (m)
    REAL(KIND=r8), INTENT(IN   ) :: fu         (npoi)      ! fraction of overall area covered by upper canopy
    REAL(KIND=r8), INTENT(IN   ) :: alaimu               ! upper canopy leaf & stem area (2 sided) for normalization of 
    ! drag coefficient (m2 m-2)
    REAL(KIND=r8), INTENT(IN   ) :: froot       (npoi,nsoilay,2) ! fraction of root in soil layer 
    REAL(KIND=r8), INTENT(INOUT) :: t34         (npoi)      ! air temperature at z34 (K)
    REAL(KIND=r8), INTENT(INOUT) :: t12         (npoi)      ! air temperature at z12 (K)
    REAL(KIND=r8), INTENT(INOUT) :: q34         (npoi)      ! specific humidity of air at z34
    REAL(KIND=r8), INTENT(INOUT) :: q12         (npoi)      ! specific humidity of air at z12
    REAL(KIND=r8), INTENT(INOUT) :: su         (npoi)      ! air-vegetation transfer coefficients (*rhoa) for upper 
    ! canopy leaves (m s-1 * kg m-3) (A39a Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(IN   ) :: cleaf             ! empirical constant in upper canopy leaf-air aerodynamic 
    ! transfer coefficient (m s-0.5) (A39a Pollard & Thompson 95)
    REAL(KIND=r8), INTENT(IN   ) :: dleaf         (2)             ! typical linear leaf dimension in aerodynamic transfer coefficient (m)
    REAL(KIND=r8), INTENT(INOUT) :: ss         (npoi)      ! air-vegetation transfer coefficients (*rhoa) for upper 
    ! canopy stems (m s-1 * kg m-3) (A39a Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(IN   ) :: cstem             ! empirical constant in upper canopy stem-air aerodynamic 
    ! transfer coefficient (m s-0.5) (A39a Pollard & Thompson 95)
    REAL(KIND=r8), INTENT(IN   ) :: dstem         (2)             ! typical linear stem dimension in aerodynamic transfer coefficient (m)
    REAL(KIND=r8), INTENT(INOUT) :: sl         (npoi)      ! air-vegetation transfer coefficients (*rhoa) for lower 
    ! canopy leaves & stems (m s-1*kg m-3) (A39a Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(IN   ) :: cgrass               ! empirical constant in lower canopy-air aerodynamic transfer
    ! coefficient (m s-0.5) (A39a Pollard & Thompson 95)
    REAL(KIND=r8), INTENT(INOUT) :: tu         (npoi)      ! temperature of upper canopy leaves (K)
    REAL(KIND=r8), INTENT(INOUT) :: ciub         (npoi)      ! intercellular co2 concentration - broadleaf (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: ciuc         (npoi)      ! intercellular co2 concentration - conifer   (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(IN   ) :: exist         (npoi,npft) ! probability of existence of each plant functional type in a gridcell
    REAL(KIND=r8), INTENT(IN   ) :: topparu  (npoi)      ! total photosynthetically active raditaion absorbed by 
    ! top leaves of upper canopy (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: csub         (npoi)      ! leaf boundary layer co2 concentration - broadleaf (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: gsub         (npoi)      ! upper canopy stomatal conductance - broadleaf  (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: csuc         (npoi)      ! leaf boundary layer co2 concentration - conifer   (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: gsuc         (npoi)      ! upper canopy stomatal conductance - conifer    (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: agcub         (npoi)      ! canopy average gross photosynthesis rate - broadleaf  (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: agcuc         (npoi)      ! canopy average gross photosynthesis rate - conifer    (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: ancub         (npoi)      ! canopy average net photosynthesis rate - broadleaf    (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: ancuc         (npoi)      ! canopy average net photosynthesis rate - conifer      (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: totcondub(npoi)      ! 
    REAL(KIND=r8), INTENT(INOUT) :: totconduc(npoi)      !
    REAL(KIND=r8), INTENT(INOUT) :: tl         (npoi)        ! temperature of lower canopy leaves & stems(K)
    REAL(KIND=r8), INTENT(INOUT) :: cils         (npoi)      ! intercellular co2 concentration - shrubs    (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: cil3         (npoi)      ! intercellular co2 concentration - c3 plants (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: cil4         (npoi)      ! intercellular co2 concentration - c4 plants (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(IN   ) :: topparl  (npoi)          ! total photosynthetically active raditaion absorbed by top 
                                                             ! leaves of lower canopy (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: csls         (npoi)      ! leaf boundary layer co2 concentration - shrubs         (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: gsls         (npoi)      ! lower canopy stomatal conductance - shrubs     (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: csl3         (npoi)      ! leaf boundary layer co2 concentration - c3 plants (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: gsl3         (npoi)      ! lower canopy stomatal conductance - c3 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: csl4         (npoi)      ! leaf boundary layer co2 concentration - c4 plants (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: gsl4         (npoi)      ! lower canopy stomatal conductance - c4 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: agcls         (npoi)      ! canopy average gross photosynthesis rate - shrubs     (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: agcl4         (npoi)      ! canopy average gross photosynthesis rate - c4 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: agcl3         (npoi)      ! canopy average gross photosynthesis rate - c3 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: ancls         (npoi)      ! canopy average net photosynthesis rate - shrubs       (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: ancl4         (npoi)      ! canopy average net photosynthesis rate - c4 grasses   (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: ancl3         (npoi)      ! canopy average net photosynthesis rate - c3 grasses   (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: totcondls(npoi)      ! 
    REAL(KIND=r8), INTENT(INOUT) :: totcondl3(npoi)      !
    REAL(KIND=r8), INTENT(INOUT) :: totcondl4(npoi)      !
    REAL(KIND=r8), INTENT(IN   ) :: chu(1:nVegClass)          ! heat capacity of upper canopy leaves per unit leaf area (J kg-1 m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wliqu         (npoi)      ! intercepted liquid h2o on upper canopy leaf area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wsnou         (npoi)      ! intercepted frozen h2o (snow) on upper canopy leaf area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: chs(1:nVegClass)          ! heat capacity of upper canopy stems per unit stem area (J kg-1 m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wliqs         (npoi)      ! intercepted liquid h2o on upper canopy stem area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wsnos         (npoi)      ! intercepted frozen h2o (snow) on upper canopy stem area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: chl(1:nVegClass)          ! heat capacity of lower canopy leaves & stems per unit leaf/stem 
                                                              ! area (J kg-1 m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wliql         (npoi)      ! intercepted liquid h2o on lower canopy leaf and stem area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wsnol         (npoi)      ! intercepted frozen h2o (snow) on lower canopy leaf & stem area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: ts            (npoi)      ! temperature of upper canopy stems (K)
    REAL(KIND=r8), INTENT(IN   ) :: frac          (npoi,npft) ! fraction of canopy occupied by each plant functional type
    REAL(KIND=r8), INTENT(IN   ) :: terml         (npoi,7)    ! term needed in lower canopy scaling
    REAL(KIND=r8), INTENT(OUT  ) :: fvapuw (npoi)              ! h2o vapor flux (evaporation from wet parts) between upper 
    ! canopy leaves and air at z12 (kg m-2 s-1/ LAI upper canopy/ fu)
    REAL(KIND=r8), INTENT(OUT  ) :: fvaput (npoi)          ! h2o vapor flux (transpiration from dry parts) between upper
    ! canopy leaves and air at z12 (kg m-2 s-1/ LAI upper canopy/ fu)
    REAL(KIND=r8), INTENT(IN   ) :: fwetl         (npoi)   ! fraction of lower canopy stem & leaf area wetted by
    ! intercepted liquid and/or snow
    REAL(KIND=r8), INTENT(IN   ) :: scalcoefl(npoi,4) ! term needed in lower canopy scaling
    REAL(KIND=r8), INTENT(IN   ) :: scalcoefu(npoi,4) ! term needed in upper canopy scalingterml  
    REAL(KIND=r8), INTENT(IN   ) :: fwetu         (npoi)   ! fraction of upper canopy leaf area wetted by intercepted liquid and/or snow
    REAL(KIND=r8), INTENT(IN   ) :: termu         (npoi,7) ! term needed in upper canopy scaling
    REAL(KIND=r8), INTENT(OUT  ) :: fwetux   (npoi)   ! fraction of upper canopy leaf area wetted if dew forms
    REAL(KIND=r8), INTENT(OUT  ) :: fwetsx   (npoi)   ! fraction of upper canopy stem area wetted if dew forms
    REAL(KIND=r8), INTENT(IN   ) :: fwets         (npoi)   ! fraction of upper canopy stem area wetted by intercepted liquid and/or snow
    REAL(KIND=r8), INTENT(OUT  ) :: fwetlx   (npoi)   ! fraction of lower canopy leaf and stem area wetted if dew forms
    REAL(KIND=r8), INTENT(IN   ) :: solu         (npoi)   ! solar flux (direct + diffuse) absorbed by upper
    ! canopy leaves per unit canopy area (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: firu         (npoi)   ! ir raditaion absorbed by upper canopy leaves (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: sols         (npoi)   ! solar flux (direct + diffuse) absorbed by upper canopy
    ! stems per unit canopy area (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: firs         (npoi)   ! ir radiation absorbed by upper canopy stems (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: soll         (npoi)   ! solar flux (direct + diffuse) absorbed by lower canopy
    ! leaves and stems per unit canopy area (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: firl         (npoi)   ! ir radiation absorbed by lower canopy leaves and stems (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: rliqu         (npoi)   ! proportion of fwetu due to liquid
    REAL(KIND=r8), INTENT(IN   ) :: rliqs         (npoi)   ! proportion of fwets due to liquid
    REAL(KIND=r8), INTENT(IN   ) :: rliql         (npoi)   ! proportion of fwetl due to liquid
    REAL(KIND=r8), INTENT(IN   ) :: pfluxu   (npoi)   ! heat flux on upper canopy leaves due to intercepted h2o (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: pfluxs   (npoi)   ! heat flux on upper canopy stems due to intercepted h2o (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: pfluxl   (npoi)   ! heat flux on lower canopy leaves & stems due to intercepted h2o (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: solg         (npoi)   ! solar flux (direct + diffuse) absorbed by unit snow-free soil (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: firg         (npoi)   ! ir radiation absorbed by soil/ice (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: soli         (npoi)   ! solar flux (direct + diffuse) absorbed by unit snow surface (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: firi         (npoi)   ! ir radiation absorbed by snow (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: fsena         (npoi)   ! downward sensible heat flux between za & z12 at za (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: fseng         (npoi)   ! upward sensible heat flux between soil surface & air at z34 (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: fseni         (npoi)   ! upward sensible heat flux between snow surface & air at z34 (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: fsenu         (npoi)   ! sensible heat flux from upper canopy leaves to air (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: fsens         (npoi)   ! sensible heat flux from upper canopy stems to air (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: fsenl         (npoi)   ! sensible heat flux from lower canopy to air (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: fvapa         (npoi)   ! downward h2o vapor flux between za & z12 at za (kg m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: fvaps         (npoi)   ! h2o vapor flux (evaporation from wet surface) between upper
    ! canopy stems and air at z12 (kg m-2 s-1 / SAI lower canopy / fu)
    REAL(KIND=r8), INTENT(OUT  ) :: fvaplw   (npoi)   ! h2o vapor flux (evaporation from wet surface) between lower
    ! canopy leaves & stems and air at z34 (kg m-2 s-1/ LAI lower canopy/ fl)
    REAL(KIND=r8), INTENT(OUT  ) :: fvaplt   (npoi)   ! h2o vapor flux (transpiration) between lower canopy &
    ! air at z34 (kg m-2 s-1 / LAI lower canopy / fl)
    REAL(KIND=r8), INTENT(OUT  ) :: fvapg         (npoi)   ! h2o vapor flux (evaporation) between soil & air at z34
    ! (kg m-2 s-1/bare ground fraction)
    REAL(KIND=r8), INTENT(OUT  ) :: fvapi         (npoi)   ! h2o vapor flux (evaporation) between snow & air at z34 (kg m-2 s-1 / fi )
    REAL(KIND=r8), INTENT(INOUT) :: firb         (npoi)   ! net upward ir radiation at reference atmospheric level za (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: bps   (npoi)! (ps/p) ** (rair/cair) for atmospheric level  (const)
    REAL(KIND=r8), INTENT(INOUT) :: rhoa  (npoi)! air density at za (allowing for h2o vapor) (kg m-3)
    REAL(KIND=r8), INTENT(INOUT) :: cp    (npoi)! specific heat of air at za (allowing for h2o vapor) (J kg-1 K-1)
    REAL(KIND=r8), INTENT(INOUT) :: za    (npoi)! height above the surface of atmospheric forcing (m)
    REAL(KIND=r8), INTENT(INOUT) :: bdl   (npoi)! aerodynamic coefficient ([(tau/rho)/u**2] for laower
    ! canopy (A31/A30 Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(INOUT) :: dil   (npoi)! inverse of momentum diffusion coefficient within lower canopy (m)
    REAL(KIND=r8), INTENT(INOUT) :: z3    (npoi)! effective top of the lower canopy (for momentum) (m)
    REAL(KIND=r8), INTENT(INOUT) :: z4    (npoi)! effective bottom of the lower canopy (for momentum) (m)
    REAL(KIND=r8), INTENT(INOUT) :: z34   (npoi)! effective middle of the lower canopy (for momentum) (m)
    REAL(KIND=r8), INTENT(INOUT) :: exphl (npoi)! exp(lamda/2*(z3-z4)) for lower canopy (A30 Pollard & Thompson)
    REAL(KIND=r8), INTENT(INOUT) :: expl  (npoi)! exphl**2
    REAL(KIND=r8), INTENT(OUT  ) :: displ (npoi)! zero-plane displacement height for lower canopy (m)
    REAL(KIND=r8), INTENT(INOUT) :: bdu   (npoi)! aerodynamic coefficient ([(tau/rho)/u**2] for upper canopy
    !  (A31/A30 Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(INOUT) :: diu   (npoi)! inverse of momentum diffusion coefficient within upper canopy (m)
    REAL(KIND=r8), INTENT(INOUT) :: z1    (npoi)! effective top of upper canopy (for momentum) (m)
    REAL(KIND=r8), INTENT(INOUT) :: z2    (npoi)! effective bottom of the upper canopy (for momentum) (m)
    REAL(KIND=r8), INTENT(INOUT) :: z12   (npoi)! effective middle of the upper canopy (for momentum) (m)
    REAL(KIND=r8), INTENT(INOUT) :: exphu (npoi)! exp(lamda/2*(z3-z4)) for upper canopy (A30 Pollard & Thompson)
    REAL(KIND=r8), INTENT(INOUT) :: expu  (npoi)! exphu**2
    REAL(KIND=r8), INTENT(OUT  ) :: dispu (npoi)! zero-plane displacement height for upper canopy (m)
    REAL(KIND=r8), INTENT(OUT  ) :: alogg (npoi)! log of soil roughness
    REAL(KIND=r8), INTENT(OUT  ) :: alogi (npoi)! log of snow roughness
    REAL(KIND=r8), INTENT(INOUT) :: alogav(npoi)! average of alogi and alogg 
    REAL(KIND=r8), INTENT(INOUT) :: alog4 (npoi)! log (max(z4, 1.1*z0sno, 1.1*z0soi)) 
    REAL(KIND=r8), INTENT(INOUT) :: alog3 (npoi)! log (z3 - displ)
    REAL(KIND=r8), INTENT(OUT  ) :: alog2 (npoi)! log (z2 - displ)
    REAL(KIND=r8), INTENT(INOUT) :: alog1 (npoi)! log (z1 - dispu) 
    REAL(KIND=r8), INTENT(INOUT) :: aloga (npoi)! log (za - dispu) 
    REAL(KIND=r8), INTENT(INOUT) :: u2    (npoi)! wind speed at level z2 (m s-1)
    REAL(KIND=r8), INTENT(INOUT) :: alogu (npoi)! log (roughness length of upper canopy)
    REAL(KIND=r8), INTENT(INOUT) :: alogl (npoi)! log (roughness length of lower canopy)
    REAL(KIND=r8), INTENT(OUT  ) :: richl (npoi)! richardson number for air above upper canopy (z3 to z2)
    REAL(KIND=r8), INTENT(OUT  ) :: straml(npoi)! momentum correction factor for stratif between upper
    ! & lower canopy (z3 to z2) (louis et al.)
    REAL(KIND=r8), INTENT(OUT  ) :: strahl(npoi)! heat/vap correction factor for stratif between upper
    ! & lower canopy (z3 to z2) (louis et al.)
    REAL(KIND=r8), INTENT(OUT  ) :: richu (npoi)! richardson number for air between upper & lower canopy (z1 to za)
    REAL(KIND=r8), INTENT(OUT  ) :: stramu(npoi)! momentum correction factor for stratif above upper
    ! canopy (z1 to za) (louis et al.)
    REAL(KIND=r8), INTENT(OUT  ) :: strahu(npoi)! heat/vap correction factor for stratif above upper
    ! canopy (z1 to za) (louis et al.)
    REAL(KIND=r8), INTENT(OUT  ) :: u1    (npoi)! wind speed at level z1 (m s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: u12   (npoi)! wind speed at level z12 (m s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: u3    (npoi)! wind speed at level z3 (m s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: u34   (npoi)! wind speed at level z34 (m s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: u4    (npoi)! wind speed at level z4 (m s-1)
    REAL(KIND=r8), INTENT(INOUT) :: cu    (npoi)! air transfer coefficient (*rhoa) (m s-1 kg m-3)
    ! for upper air region (z12 --> za) (A35 Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(INOUT) :: cl    (npoi)! air transfer coefficient (*rhoa) (m s-1 kg m-3)
    ! between the 2 canopies (z34 --> z12) (A36 Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(INOUT) :: sg    (npoi)! air-soil transfer coefficient
    REAL(KIND=r8), INTENT(INOUT) :: si    (npoi)! air-snow transfer coefficient
    REAL(KIND=r8), INTENT(IN   ) :: ux  (npoi)
    REAL(KIND=r8), INTENT(IN   ) :: uy  (npoi)
    REAL(KIND=r8), INTENT(OUT  ) :: taux(npoi)
    REAL(KIND=r8), INTENT(OUT  ) :: tauy(npoi)
    REAL(KIND=r8), INTENT(OUT  ) :: ts2 (npoi) 
    REAL(KIND=r8), INTENT(OUT  ) :: qs2 (npoi)
    REAL(KIND=r8), INTENT(IN   ) :: vegtype0(npoi)
    REAL(KIND=r8), INTENT(IN   ) :: stressfac(nVegClass)
    !
    ! Local variables
    !
    REAL(KIND=r8) :: xu    (npoi)  ! SAVE
    REAL(KIND=r8) :: xs    (npoi)  ! SAVE
    REAL(KIND=r8) :: xl    (npoi)  ! SAVE
    REAL(KIND=r8) :: chux  (npoi)  ! SAVE
    REAL(KIND=r8) :: chsx  (npoi)  ! SAVE
    REAL(KIND=r8) :: chlx  (npoi)  ! SAVE
    REAL(KIND=r8) :: chgx  (npoi)  ! SAVE
    REAL(KIND=r8) :: wlgx  (npoi)  ! SAVE
    REAL(KIND=r8) :: wigx  (npoi)  ! SAVE
    REAL(KIND=r8) :: cog   (npoi)  ! SAVE
    REAL(KIND=r8) :: coi   (npoi)  ! SAVE
    REAL(KIND=r8) :: zirg  (npoi)  ! SAVE
    REAL(KIND=r8) :: ziri  (npoi)  ! SAVE
    REAL(KIND=r8) :: wu    (npoi)  ! SAVE
    REAL(KIND=r8) :: ws    (npoi)  ! SAVE
    REAL(KIND=r8) :: wl    (npoi)  ! SAVE
    REAL(KIND=r8) :: wg    (npoi)  ! SAVE
    REAL(KIND=r8) :: wi    (npoi)  ! SAVE
    REAL(KIND=r8) :: tuold (npoi)  ! SAVE
    REAL(KIND=r8) :: tsold (npoi)  ! SAVE
    REAL(KIND=r8) :: tlold (npoi)  ! SAVE
    REAL(KIND=r8) :: tgold (npoi)  ! SAVE
    REAL(KIND=r8) :: tiold (npoi)  ! SAVE
    INTEGER :: niter         ! total number of ierations
    INTEGER :: iter          ! number of iteration
    REAL(KIND=r8) :: cdmaxa, cdmaxb, ctau,fmom,umom,vmom
    INTEGER :: i

    CALL canini(bps      , & ! INTENT(OUT  )
         rhoa      , & ! INTENT(OUT  )
         cp        , & ! INTENT(OUT  )
         za        , & ! INTENT(INOUT)
         bdl       , & ! INTENT(OUT  )
         dil       , & ! INTENT(OUT  )
         z3        , & ! INTENT(OUT  )
         z4        , & ! INTENT(OUT  )
         z34       , & ! INTENT(OUT  )
         exphl     , & ! INTENT(OUT  )
         expl      , & ! INTENT(OUT  )
         displ     , & ! INTENT(OUT  )
         bdu       , & ! INTENT(OUT  )
         diu       , & ! INTENT(OUT  )
         z1        , & ! INTENT(OUT  )
         z2        , & ! INTENT(OUT  )
         z12       , &! INTENT(OUT  )
         exphu     , & ! INTENT(OUT  )
         expu      , & ! INTENT(OUT  )
         dispu     , & ! INTENT(OUT  )
         alogg     , & ! INTENT(OUT  )
         alogi     , & ! INTENT(OUT  )
         alogav    , & ! INTENT(OUT  )
         alog4     , & ! INTENT(OUT  )
         alog3     , & ! INTENT(OUT  )
         alog2     , & ! INTENT(OUT  )
         alog1     , & ! INTENT(OUT  )
         aloga     , & ! INTENT(OUT  )
         u2        , & ! INTENT(OUT  )
         alogu     , & ! INTENT(OUT  )
         alogl     , & ! INTENT(OUT  )
         ztop      , & ! INTENT(IN   )
         fl        , & ! INTENT(IN   )
         lai       , & ! INTENT(IN   )
         sai       , & ! INTENT(IN   )
         alaiml    , & ! INTENT(IN   )
         zbot      , & ! INTENT(IN   )
         fu        , & ! INTENT(IN   )
         alaimu    , & ! INTENT(IN   )
         z0soi     , & ! INTENT(IN   )
         fi        , & ! INTENT(IN   )
         z0sno     , & ! INTENT(IN   )
         psurf     , & ! INTENT(IN   )
         ta        , & ! INTENT(IN   )
         qa        , & ! INTENT(IN   )
         ua        , & ! INTENT(IN   )
         npoi      , & ! INTENT(IN   ) 
         cappa     , & ! INTENT(IN   ) 
         rair      , & ! INTENT(IN   ) 
         rvap      , & ! INTENT(IN   ) 
         cair      , & ! INTENT(IN   ) 
         cvap      , & ! INTENT(IN   ) 
         grav        ) ! INTENT(IN   ) 
    !
    ! estimate soil moisture stress parameters
    !
    CALL drystress(froot   , &! INTENT(IN   )
         wsoi    , &! INTENT(IN   )
         wisoi   , &! INTENT(IN   )
         swilt   , &! INTENT(IN   )
         sfield  , &! INTENT(IN   ) 
         stressl , &! INTENT(OUt  ) 
         stressu , &! INTENT(OUt  ) 
         stresstl, &! INTENT(OUt  ) 
         stresstu, &! INTENT(OUt  ) 
         vegtype0, &! INTENT(IN   ) 
         stressfac, &! INTENT(IN   ) 
         nVegClass , &! INTENT(IN   ) 
         npoi    , &! INTENT(IN   )
         nsoilay   )! INTENT(IN   )
    !
    ! iterate the whole canopy physics solution niter times:niter = 3
    !
    niter = 3 !3 
    !
    DO iter = 1, niter
       !
       ! calculate wind speeds and aerodynamic transfer coeffs
       !
       CALL turcof (   iter      , &! INTENT(IN   )
            z3        , &! INTENT(IN   )
            z2        , &! INTENT(IN   )
            alogl     , &! INTENT(INOUT)
            u2        , &! INTENT(INOUT)
            richl     , &! INTENT(OUT  )
            straml    , &! INTENT(OUT  )
            strahl    , &! INTENT(OUT  )
            bps       , &! INTENT(IN   )
            z1        , &! INTENT(IN   )
            za        , &! INTENT(IN   )
            alogu     , &! INTENT(INOUT)
            aloga     , &! INTENT(IN   )
            richu     , &! INTENT(OUT  )
            stramu    , &! INTENT(OUT  )
            strahu    , &! INTENT(OUT  )
            alog4     , &! INTENT(OUT  )
            alogav    , &! INTENT(IN   )
            bdl       , &! INTENT(IN   )
            expl      , &! INTENT(IN   )
            alog3     , &! INTENT(IN   )
            bdu       , &! INTENT(IN   )
            expu      , &! INTENT(IN   )
            alog1     , &! INTENT(IN   )
            u1        , &! INTENT(OUT  )
            u12       , &! INTENT(OUT  )
            exphu     , &! INTENT(IN   )
            u3        , &! INTENT(OUT  )
            u34       , &! INTENT(OUT  )
            exphl     , &! INTENT(IN   )
            u4        , &! INTENT(OUT  )
            rhoa      , &! INTENT(IN   )
            diu       , &! INTENT(IN   )
            z12       , &! INTENT(IN   )
            dil       , &! INTENT(IN   )
            z34       , &! INTENT(IN   )
            z4        , &! INTENT(IN   )
            cu        , &! INTENT(OUT  )
            cl        , &! INTENT(OUT  )
            sg        , &! INTENT(OUT  )
            si        , &! INTENT(OUT  )
            alog2     , &! INTENT(OUT  )
            t34       , &! INTENT(IN   )
            t12       , &! INTENT(IN   )
            q34       , &! INTENT(IN   )
            q12       , &! INTENT(IN   )
            su        , &! INTENT(OUT  )
            cleaf     , &! INTENT(IN   )
            dleaf     , &! INTENT(IN   )
            ss        , &! INTENT(OUT  )
            cstem     , &! INTENT(IN   )
            dstem     , &! INTENT(IN   )
            sl        , &! INTENT(OUT  )
            cgrass    , &! INTENT(IN   )
            ua        , &! INTENT(IN   )
            ta        , &! INTENT(IN   )
            qa        , &! INTENT(IN   )
            npoi      , &! INTENT(IN   )
            dtime     , &! INTENT(IN   )
            vonk      , &! INTENT(IN   )
            grav        )! INTENT(IN   )
       !
       ! calculate canopy photosynthesis rates and conductance
       !
       CALL stomata(tau15          , &! INTENT(IN   )
            kc15         , &! INTENT(IN   )
            ko15         , &! INTENT(IN   )
            cimax        , &! INTENT(IN   )
            vmax_pft     , &! INTENT(IN   )
            gammaub      , &! INTENT(IN   )
            alpha3       , &! INTENT(IN   )
            theta3       , &! INTENT(IN   )
            beta3        , &! INTENT(IN   )
            coefmub      , &! INTENT(IN   )
            coefbub      , &! INTENT(IN   )
            gsubmin      , &! INTENT(IN   )
            gammauc      , &! INTENT(IN   )
            coefmuc      , &! INTENT(IN   )
            coefbuc      , &! INTENT(IN   )
            gsucmin      , &! INTENT(IN   )
            gammals      , &! INTENT(IN   )
            coefmls      , &! INTENT(IN   )
            coefbls      , &! INTENT(IN   )
            gslsmin      , &! INTENT(IN   )
            gammal3      , &! INTENT(IN   )
            coefml3      , &! INTENT(IN   )
            coefbl3      , &! INTENT(IN   )
            gsl3min      , &! INTENT(IN   )
            gammal4      , &! INTENT(IN   )
            alpha4       , &! INTENT(IN   )
            theta4       , &! INTENT(IN   )
            beta4        , &! INTENT(IN   )
            coefml4      , &! INTENT(IN   )
            coefbl4      , &! INTENT(IN   )
            gsl4min      , &! INTENT(IN   )
            a10scalparamu, &! INTENT(INOUT)
            a10daylightu , &! INTENT(INOUT)
            a10scalparaml, &! INTENT(INOUT)
            a10daylightl , &! INTENT(INOUT)
            fwetl        , &! INTENT(IN   )
            scalcoefl    , &! INTENT(IN   )
            terml        , &! INTENT(IN   )
            fwetu        , &! INTENT(IN   )
            scalcoefu    , &! INTENT(IN   )
            termu        , &! INTENT(IN   )
            tu           , &! INTENT(IN   )
            ciub         , &! INTENT(INOUT)
            ciuc         , &! INTENT(INOUT)
            su                  , &! INTENT(IN   )
            t12          , &! INTENT(IN   )
            q12          , &! INTENT(IN   )
            exist        , &! INTENT(IN   )
            topparu      , &! INTENT(IN   )
            csub         , &! INTENT(INOUT)
            gsub         , &! INTENT(INOUT)
            csuc         , &! INTENT(INOUT)
            gsuc          , &! INTENT(INOUT)
            lai          , &! INTENT(IN   )
            sai          , &! INTENT(IN   )
            agcub        , &! INTENT(OUT)
            agcuc        , &! INTENT(OUT)
            ancub        , &! INTENT(OUT)
            ancuc        , &! INTENT(OUT)
            totcondub    , &! INTENT(OUT)
            totconduc    , &! INTENT(OUT)
            tl           , &! INTENT(IN   )
            cils         , &! INTENT(INOUT)
            cil3         , &! INTENT(INOUT)
            cil4         , &! INTENT(INOUT)
            sl           , &! INTENT(IN   )
            t34          , &! INTENT(IN   )
            q34          , &! INTENT(IN   )
            topparl      , &! INTENT(IN   )
            csls         , &! INTENT(INOUT)
            gsls         , &! INTENT(INOUT)
            csl3         , &! INTENT(INOUT)
            gsl3         , &! INTENT(INOUT)
            csl4         , &! INTENT(INOUT)
            gsl4         , &! INTENT(INOUT)
            agcls        , &! INTENT(OUT)
            agcl4        , &! INTENT(OUT)
            agcl3        , &! INTENT(OUT)
            ancls        , &! INTENT(OUT)
            ancl4        , &! INTENT(OUT)
            ancl3        , &! INTENT(OUT)
            totcondls    , &! INTENT(OUT)
            totcondl3    , &! INTENT(OUT)
            totcondl4    , &! INTENT(OUT)
            stresstu     , &! INTENT(IN   )
            stresstl     , &! INTENT(IN   )
            o2conc       , &! INTENT(IN   )
            psurf        , &! INTENT(IN   )
            co2conc      , &! INTENT(IN   )
            npoi         , &! INTENT(IN   )
            npft         , &! INTENT(IN   )
            epsilon      , &! INTENT(IN   )
            dtime        , &! INTENT(IN   )
            niter	 , &
            iter	   )




       !
       ! solve implicit system of heat and water balance equations
       !
       CALL turvap (iter , &! INTENT(IN   )
            niter        , &! INTENT(IN   )
            cp           , &! INTENT(IN   )
            sg           , &! INTENT(IN   )
            fwetux       , &! INTENT(OUT  )
            fwetu        , &! INTENT(IN   )
            fwetsx       , &! INTENT(OUT  )
            fwets        , &! INTENT(IN   )
            fwetlx       , &! INTENT(OUT  )
            fwetl        , &! INTENT(IN   )
            solu         , &! INTENT(IN   )
            firu         , &! INTENT(IN   )
            sols         , &! INTENT(IN   )
            firs         , &! INTENT(IN   )
            soll         , &! INTENT(IN   )
            firl         , &! INTENT(IN   )
            rliqu        , &! INTENT(IN   )
            rliqs        , &! INTENT(IN   )
            rliql        , &! INTENT(IN   )
            pfluxu       , &! INTENT(IN   )
            pfluxs       , &! INTENT(IN   )
            pfluxl       , &! INTENT(IN   )
            cu           , &! INTENT(IN   )
            cl           , &! INTENT(IN   )
            bps          , &! INTENT(IN   )
            si           , &! INTENT(IN   )
            solg         , &! INTENT(IN   )
            firg         , &! INTENT(INOUT)
            soli         , &! INTENT(IN   )
            firi         , &! INTENT(INOUT)
            fsena        , &! INTENT(OUT  )
            fseng        , &! INTENT(OUT  )
            fseni        , &! INTENT(OUT  )
            fsenu        , &! INTENT(OUT  )
            fsens        , &! INTENT(OUT  )
            fsenl        , &! INTENT(OUT  )
            fvapa        , &! INTENT(OUT  )
            fvapuw       , &! INTENT(OUT  )
            fvaput       , &! INTENT(OUT  )
            fvaps        , &! INTENT(OUT  )
            fvaplw       , &! INTENT(OUT  )
            fvaplt       , &! INTENT(OUT  )
            fvapg        , &! INTENT(OUT  )
            fvapi        , &! INTENT(OUT  )
            firb         , &! INTENT(INOUT)
            lai          , &! INTENT(IN   )
            fu           , &! INTENT(IN   )
            sai          , &! INTENT(IN   )
            fl           , &! INTENT(IN   )
            chu(1:nVegClass)    , &! INTENT(IN   )
            wliqu        , &! INTENT(IN   )
            wsnou        , &! INTENT(IN   )
            chs(1:nVegClass)    , &! INTENT(IN   )
            wliqs        , &! INTENT(IN   )
            wsnos        , &! INTENT(IN   )
            chl(1:nVegClass)    , &! INTENT(IN   )
            wliql        , &! INTENT(IN   )
            wsnol        , &! INTENT(IN   )
            tu           , &! INTENT(INOUT)
            ts           , &! INTENT(INOUT)
            tl           , &! INTENT(INOUT)
            q34          , &! INTENT(INOUT)
            t34          , &! INTENT(INOUT)
            q12          , &! INTENT(INOUT)
            su           , &! INTENT(IN   )
            totcondub    , &! INTENT(IN   )
            frac         , &! INTENT(IN   )
            totconduc    , &! INTENT(IN   )
            ss           , &! INTENT(IN   )
            sl           , &! INTENT(IN   )
            totcondls    , &! INTENT(IN   )
            totcondl4    , &! INTENT(IN   )
            totcondl3    , &! INTENT(IN   )
            t12          , &! INTENT(INOUT)
            poros        , &! INTENT(IN   )
            wpud         , &! INTENT(IN   )
            wipud        , &! INTENT(IN   )
            csoi         , &! INTENT(IN   )
            rhosoi       , &! INTENT(IN   )
            wisoi        , &! INTENT(IN   )
            wsoi         , &! INTENT(IN   )
            hsoi         , &! INTENT(IN   )
            consoi       , &! INTENT(IN   )
            tg           , &! INTENT(INOUT)
            ti           , &! INTENT(INOUT)
            wpudmax      , &! INTENT(IN   )
            suction      , &! INTENT(IN   )
            bex          , &! INTENT(IN   )
            swilt        , &! INTENT(IN   )
            hvasug       , &! INTENT(IN   )
            tsoi         , &! INTENT(IN   )
            hvasui       , &! INTENT(IN   )
            upsoiu       , &! INTENT(OUT  )
            stressu      , &! INTENT(IN   )
            stresstu     , &! INTENT(IN   )
            upsoil       , &! INTENT(OUT  )
            stressl      , &! INTENT(IN   )
            stresstl     , &! INTENT(IN   )
            fi           , &! INTENT(IN   )
            consno       , &! INTENT(IN   )
            hsno         , &! INTENT(IN   )
            hsnotop      , &! INTENT(IN   )
            tsno         , &! INTENT(IN   )
            psurf        , &! INTENT(IN   )
            ta           , &! INTENT(IN   )
            qa           , &! INTENT(IN   )
            ginvap       , &! INTENT(OUT  )
            gsuvap       , &! INTENT(OUT  )
            gtrans       , &! INTENT(OUT  )
            gtransu      , &! INTENT(OUT  )
            gtransl      , &! INTENT(OUT  )
            xu           , &! INTENT(INOUT)
            xs           , &! INTENT(INOUT)
            xl           , &! INTENT(INOUT)
            chux         , &! INTENT(INOUT)
            chsx         , &! INTENT(INOUT)
            chlx         , &! INTENT(INOUT)
            chgx         , &! INTENT(INOUT)
            wlgx         , &! INTENT(INOUT)
            wigx         , &! INTENT(INOUT)
            cog          , &! INTENT(INOUT)
            coi          , &! INTENT(INOUT)
            zirg         , &! INTENT(INOUT)
            ziri         , &! INTENT(INOUT)
            wu           , &! INTENT(INOUT)
            ws           , &! INTENT(INOUT)
            wl           , &! INTENT(INOUT)
            wg           , &! INTENT(INOUT)
            wi           , &! INTENT(INOUT)
            tuold        , &! INTENT(INOUT)
            tsold        , &! INTENT(INOUT)
            tlold        , &! INTENT(INOUT)
            tgold        , &! INTENT(INOUT)
            tiold        , &! INTENT(INOUT)
            npoi         , &! INTENT(IN   )
            nsoilay      , &! INTENT(IN   )
            nsnolay      , &! INTENT(IN   )
            npft         , &! INTENT(IN   )
            hvap         , &! INTENT(IN   )
            cvap         , &! INTENT(IN   )
            ch2o         , &! INTENT(IN   )
            hsub         , &! INTENT(IN   )
            cice         , &! INTENT(IN   )
            rhow         , &! INTENT(IN   )
            stef         , &! INTENT(IN   )
            tmelt        , &! INTENT(IN   )
            hfus         , &! INTENT(IN   )
            epsilon      , &! INTENT(IN   )
            grav         , &! INTENT(IN   )
            rvap         , &! INTENT(IN   )
            vegtype0     , &! INTENT(IN   )
            dtime        , &! INTENT(IN   )
            nVegClass      )! INTENT(IN   )
       !
    END DO
    !
    cdmaxa = 300.0_r8/(2.0_r8*dtime)
    !cdmaxa = 0.0625_r8 
    !cdmaxa = 0.125_r8 
    cdmaxb = 1e20_r8

    DO i = 1, npoi
       ctau    = ua(i) * (vonk / (aloga(i) - alogu(i)))**2 * stramu(i)
       IF(ctau/cdmaxb /= -1.0_r8)THEN
          ctau    = MIN (cdmaxa, ctau / (1.0_r8 + ctau/cdmaxb))
       ELSE
          ctau    = cdmaxa
       END IF
       taux(i) = rhoa(i) * ctau * ux(i)
       tauy(i) = rhoa(i) * ctau * uy(i)
    END DO
    !
    ! Calculate 2-m surface air temperature (diagnostic, for history)
    ! Arguments are 1-dimensional --> can be passed only for the kpti-->kptj
    ! and tscreen doesn't use any array defined over all land points.
    !
    CALL tscreen  (ts2    , &! INTENT(OUT  )
         2.0_r8 , &! INTENT(IN   )
         za     , &! INTENT(IN   ) 
         z1     , &! INTENT(IN   )
         z12    , &! INTENT(IN   )
         z34    , &! INTENT(IN   )
         dispu  , &! INTENT(IN   )
         ta     , &! INTENT(IN   )
         t12    , &! INTENT(IN   )
         t34    , &! INTENT(IN   )
         npoi     )! INTENT(IN   )

    CALL tscreen  (qs2    , &! INTENT(OUT  )
         2.0_r8 , &! INTENT(IN   )
         za     , &! INTENT(IN   ) 
         z1     , &! INTENT(IN   )
         z12    , &! INTENT(IN   )
         z34    , &! INTENT(IN   )
         dispu  , &! INTENT(IN   )
         qa     , &! INTENT(IN   )
         q12    , &! INTENT(IN   )
         q34    , &! INTENT(IN   )
         npoi     )! INTENT(IN   )

    RETURN
  END SUBROUTINE canopy
  ! ---------------------------------------------------------------------
  SUBROUTINE tscreen(tscr    , &
       zscr    , &
       za      , &
       z1      , &
       z12     , &
       z34     , &
       dispu   , &
       ta      , &
       t12     , &
       t34     , &
       npoi       )
    ! ---------------------------------------------------------------------
    !        Interpolates diagnostic screen-height temperature tscr at
    !        height zscr.
    !
    !---------------------------------------------------------------

    IMPLICIT NONE
    !-----------------------------------------------------------------------
    !
    ! input variables
    !
    INTEGER, INTENT(IN   ) :: npoi          ! number of points in little vector 
    !
    REAL(KIND=r8), INTENT(IN   ) :: zscr            ! refernce height
    !
    REAL(KIND=r8), INTENT(OUT  ) :: tscr(npoi)
    REAL(KIND=r8), INTENT(IN   ) :: za   (npoi)
    REAL(KIND=r8), INTENT(IN   ) :: z1   (npoi)
    REAL(KIND=r8), INTENT(IN   ) :: z12  (npoi)
    REAL(KIND=r8), INTENT(IN   ) :: z34  (npoi)
    REAL(KIND=r8), INTENT(IN   ) :: dispu(npoi)
    REAL(KIND=r8), INTENT(IN   ) :: ta   (npoi) 
    REAL(KIND=r8), INTENT(IN   ) :: t12  (npoi)
    REAL(KIND=r8), INTENT(IN   ) :: t34  (npoi)
    !
    ! local variables
    !
    INTEGER :: i
    REAL(KIND=r8) :: w
    !
    DO i = 1, npoi

       IF (zscr.GT.z1(i)) THEN

          !         above upper canopy:
          w = LOG((zscr  -dispu(i)) / (z1(i)-dispu(i)))  &
               / LOG((za(i) -dispu(i)) / (z1(i)-dispu(i)))
          tscr(i) = w*ta(i) + (1.0_r8-w)*t12(i)

       ELSE IF (zscr.GT.z12(i)) THEN

          !         within top half of upper canopy:
          tscr(i) = t12(i)

       ELSE IF (zscr.GT.z34(i)) THEN

          !         between mid-points of canopies:
          tscr(i) = ( (zscr-z34(i))*t12(i) + (z12(i)-zscr)*t34(i) )  &
               / (z12(i)-z34(i))

       ELSE

          !         within or below lower canopy:
          tscr(i) = t34(i)

       END IF

    END DO

    RETURN
  END SUBROUTINE tscreen


  !
  ! #####   #    #   #   #   ####      #     ####   #        ####    ####    #   #
  ! #    #  #    #    # #   #          #    #    #  #       #    #  #    #    # #
  ! #    #  ######     #     ####      #    #    #  #       #    #  #          #
  ! #####   #    #     #         #     #    #    #  #       #    #  #  ###     #
  ! #       #    #     #    #    #     #    #    #  #       #    #  #    #     #
  ! #       #    #     #     ####      #     ####   ######   ####    ####      #
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE stomata(               &
       tau15        , &! INTENT(IN   )
       kc15              , &! INTENT(IN   )
       ko15              , &! INTENT(IN   )
       cimax        , &! INTENT(IN   )
       vmax_pft     , &! INTENT(IN   )
       gammaub      , &! INTENT(IN   )
       alpha3       , &! INTENT(IN   )
       theta3       , &! INTENT(IN   )
       beta3        , &! INTENT(IN   )
       coefmub      , &! INTENT(IN   )
       coefbub      , &! INTENT(IN   )
       gsubmin      , &! INTENT(IN   )
       gammauc      , &! INTENT(IN   )
       coefmuc      , &! INTENT(IN   )
       coefbuc      , &! INTENT(IN   )
       gsucmin      , &! INTENT(IN   )
       gammals      , &! INTENT(IN   )
       coefmls      , &! INTENT(IN   )
       coefbls      , &! INTENT(IN   )
       gslsmin      , &! INTENT(IN   )
       gammal3      , &! INTENT(IN   )
       coefml3      , &! INTENT(IN   )
       coefbl3      , &! INTENT(IN   )
       gsl3min      , &! INTENT(IN   )
       gammal4      , &! INTENT(IN   )
       alpha4       , &! INTENT(IN   )
       theta4       , &! INTENT(IN   )
       beta4        , &! INTENT(IN   )
       coefml4      , &! INTENT(IN   )
       coefbl4      , &! INTENT(IN   )
       gsl4min      , &! INTENT(IN   )
       a10scalparamu, &! INTENT(INOUT)
       a10daylightu , &! INTENT(INOUT)
       a10scalparaml, &! INTENT(INOUT)
       a10daylightl , &! INTENT(INOUT)
       fwetl        , &! INTENT(IN   )
       scalcoefl    , &! INTENT(IN   )
       terml              , &! INTENT(IN   )
       fwetu              , &! INTENT(IN   )
       scalcoefu    , &! INTENT(IN   )
       termu        , &! INTENT(IN   )
       tu           , &! INTENT(IN   )
       ciub         , &! INTENT(INOUT)
       ciuc         , &! INTENT(INOUT)
       su           , &! INTENT(IN   )
       t12          , &! INTENT(IN   )
       q12          , &! INTENT(IN   )
       exist              , &! INTENT(IN   )
       topparu      , &! INTENT(IN   )
       csub         , &! INTENT(INOUT)
       gsub         , &! INTENT(INOUT)
       csuc         , &! INTENT(INOUT)
       gsuc         , &! INTENT(INOUT)
       lai          , &! INTENT(IN   )
       sai          , &! INTENT(IN   )
       agcub        , &! INTENT(OUT)
       agcuc        , &! INTENT(OUT)
       ancub        , &! INTENT(OUT)
       ancuc        , &! INTENT(OUT)
       totcondub    , &! INTENT(OUT)
       totconduc    , &! INTENT(OUT)
       tl           , &! INTENT(IN   )
       cils         , &! INTENT(INOUT)
       cil3         , &! INTENT(INOUT)
       cil4         , &! INTENT(INOUT)
       sl           , &! INTENT(IN   )
       t34          , &! INTENT(IN   )
       q34          , &! INTENT(IN   )
       topparl      , &! INTENT(IN   )
       csls         , &! INTENT(INOUT)
       gsls         , &! INTENT(INOUT)
       csl3         , &! INTENT(INOUT)
       gsl3         , &! INTENT(INOUT)
       csl4         , &! INTENT(INOUT)
       gsl4         , &! INTENT(INOUT)
       agcls              , &! INTENT(OUT)
       agcl4        , &! INTENT(OUT)
       agcl3        , &! INTENT(OUT)
       ancls        , &! INTENT(OUT)
       ancl4        , &! INTENT(OUT)
       ancl3        , &! INTENT(OUT)
       totcondls    , &! INTENT(OUT)
       totcondl3    , &! INTENT(OUT)
       totcondl4    , &! INTENT(OUT)
       stresstu     , &! INTENT(IN   )
       stresstl     , &! INTENT(IN   )
       o2conc       , &! INTENT(IN   )
       psurf        , &! INTENT(IN   )
       co2conc      , &! INTENT(IN   )
       npoi         , &! INTENT(IN   )
       npft         , &! INTENT(IN   )
       epsilon      , &! INTENT(IN   )
       dtime	    , &! INTENT(IN   )
       niter	    , &
       iter	      )
    ! ---------------------------------------------------------------------
    !
    ! common blocks
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN   )  :: npoi          ! total number of land points
    INTEGER, INTENT(IN   )  :: npft          ! number of plant functional types
    INTEGER, INTENT(IN   )  :: niter
    INTEGER, INTENT(IN   )  :: iter
    REAL(KIND=r8), INTENT(IN   )  :: epsilon  ! small quantity to avoid zero-divides and other
    ! truncation or machine-limit troubles with small
    ! values. should be slightly greater than o(1)
    ! machine precision
    REAL(KIND=r8) , INTENT(IN   ) :: dtime    ! model timestep (seconds)
    REAL(KIND=r8) , INTENT(IN   ) :: o2conc             ! o2 concentration (mol/mol)
    REAL(KIND=r8) , INTENT(IN   ) :: psurf     (npoi)         ! surface pressure (Pa)  &
    REAL(KIND=r8) , INTENT(IN   ) :: co2conc  (npoi)          ! co2 concentration (mol/mol)  &
    REAL(KIND=r8) , INTENT(IN   ) :: stresstu(npoi)      ! sum of stressu over all 6 soil layers (dimensionless)
    REAL(KIND=r8) , INTENT(IN   ) :: stresstl(npoi)      ! sum of stressl over all 6 soil layers (dimensionless)
    REAL(KIND=r8) , INTENT(IN   ) :: tu       (npoi)     ! temperature of upper canopy leaves (K)
    REAL(KIND=r8) , INTENT(INOUT) :: ciub     (npoi)     ! intercellular co2 concentration - broadleaf (mol_co2/mol_air)
    REAL(KIND=r8) , INTENT(INOUT) :: ciuc     (npoi)     ! intercellular co2 concentration - conifer   (mol_co2/mol_air)
    REAL(KIND=r8) , INTENT(IN   ) :: su       (npoi)     ! air-vegetation transfer coefficients (*rhoa)
    ! for upper canopy leaves (m s-1 * kg m-3) (A39a Pollard & Thompson 1995)
    REAL(KIND=r8) , INTENT(IN   ) :: t12      (npoi)     ! air temperature at z12 (K)
    REAL(KIND=r8) , INTENT(IN   ) :: q12      (npoi)     ! specific humidity of air at z12
    REAL(KIND=r8) , INTENT(IN   ) :: exist    (npoi,npft)! probability of existence of each plant functional type in a gridcell
    REAL(KIND=r8) , INTENT(IN   ) :: topparu  (npoi)     ! total photosynthetically active raditaion
    ! absorbed by top leaves of upper canopy (W m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: csub     (npoi)     ! leaf boundary layer co2 concentration - broadleaf (mol_co2/mol_air)
    REAL(KIND=r8) , INTENT(INOUT) :: gsub     (npoi)     ! upper canopy stomatal conductance - broadleaf  (mol_co2 m-2 s-1)
    REAL(KIND=r8) , INTENT(INOUT) :: csuc     (npoi)     ! leaf boundary layer co2 concentration - conifer(mol_co2/mol_air)
    REAL(KIND=r8) , INTENT(INOUT) :: gsuc     (npoi)     ! upper canopy stomatal conductance - conifer    (mol_co2 m-2 s-1)
    REAL(KIND=r8) , INTENT(IN   ) :: lai      (npoi,2)   ! canopy single-sided leaf area index (area leaf/area veg)
    REAL(KIND=r8) , INTENT(IN   ) :: sai      (npoi,2)   ! current single-sided stem area index
    REAL(KIND=r8) , INTENT(OUT  ) :: agcub    (npoi)     ! canopy average gross photosynthesis rate - broadleaf  (mol_co2 m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: agcuc    (npoi)     ! canopy average gross photosynthesis rate - conifer    (mol_co2 m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: ancub    (npoi)     ! canopy average net photosynthesis rate - broadleaf    (mol_co2 m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: ancuc    (npoi)     ! canopy average net photosynthesis rate - conifer      (mol_co2 m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: totcondub(npoi)     ! 
    REAL(KIND=r8) , INTENT(OUT  ) :: totconduc(npoi)     !
    REAL(KIND=r8) , INTENT(IN   ) :: tl       (npoi)     ! temperature of lower canopy leaves & stems(K)
    REAL(KIND=r8) , INTENT(INOUT) :: cils     (npoi)     ! intercellular co2 concentration - shrubs    (mol_co2/mol_air)
    REAL(KIND=r8) , INTENT(INOUT) :: cil3     (npoi)     ! intercellular co2 concentration - c3 plants (mol_co2/mol_air)
    REAL(KIND=r8) , INTENT(INOUT) :: cil4     (npoi)     ! intercellular co2 concentration - c4 plants (mol_co2/mol_air)
    REAL(KIND=r8) , INTENT(IN   ) :: sl       (npoi)     ! air-vegetation transfer coefficients (*rhoa) for 
    ! lower canopy leaves & stems (m s-1*kg m-3) (A39a Pollard & Thompson 1995)
    REAL(KIND=r8) , INTENT(IN   ) :: t34      (npoi)     ! air temperature at z34 (K)
    REAL(KIND=r8) , INTENT(IN   ) :: q34      (npoi)     ! specific humidity of air at z34
    REAL(KIND=r8) , INTENT(IN   ) :: topparl  (npoi)     ! total photosynthetically active raditaion absorbed by 
    ! top leaves of lower canopy (W m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: csls     (npoi)     ! leaf boundary layer co2 concentration - shrubs    (mol_co2/mol_air)
    REAL(KIND=r8) , INTENT(INOUT) :: gsls     (npoi)     ! lower canopy stomatal conductance - shrubs     (mol_co2 m-2 s-1)
    REAL(KIND=r8) , INTENT(INOUT) :: csl3     (npoi)     ! leaf boundary layer co2 concentration - c3 plants (mol_co2/mol_air)
    REAL(KIND=r8) , INTENT(INOUT) :: gsl3     (npoi)     ! lower canopy stomatal conductance - c3 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8) , INTENT(INOUT) :: csl4     (npoi)     ! leaf boundary layer co2 concentration - c4 plants (mol_co2/mol_air)
    REAL(KIND=r8) , INTENT(INOUT) :: gsl4     (npoi)     ! lower canopy stomatal conductance - c4 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: agcls    (npoi)          ! canopy average gross photosynthesis rate - shrubs (mol_co2 m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: agcl4    (npoi)          ! canopy average gross photosynthesis rate - c4 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: agcl3    (npoi)          ! canopy average gross photosynthesis rate - c3 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: ancls    (npoi)          ! canopy average net photosynthesis rate - shrubs   (mol_co2 m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: ancl4    (npoi)          ! canopy average net photosynthesis rate - c4 grasses       (mol_co2 m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: ancl3    (npoi)          ! canopy average net photosynthesis rate - c3 grasses       (mol_co2 m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: totcondls(npoi)     ! 
    REAL(KIND=r8) , INTENT(OUT  ) :: totcondl3(npoi)     !
    REAL(KIND=r8) , INTENT(OUT  ) :: totcondl4(npoi)     !
    REAL(KIND=r8) , INTENT(IN   ) :: fwetl     (npoi)   ! fraction of lower canopy stem & leaf area wetted by 
    !intercepted liquid and/or snow
    REAL(KIND=r8) , INTENT(IN   ) :: scalcoefl (npoi,4) ! term needed in lower canopy scaling
    REAL(KIND=r8) , INTENT(IN   ) :: terml     (npoi,7) ! term needed in lower canopy scaling
    REAL(KIND=r8) , INTENT(IN   ) :: fwetu     (npoi)   ! fraction of upper canopy leaf area wetted by intercepted liquid and/or snow
    REAL(KIND=r8) , INTENT(IN   ) :: scalcoefu (npoi,4) ! term needed in upper canopy scaling
    REAL(KIND=r8) , INTENT(IN   ) :: termu     (npoi,7) ! term needed in upper canopy scaling
    REAL(KIND=r8) , INTENT(INOUT) :: a10scalparamu(npoi) ! 10-day average day-time scaling parameter - upper canopy (dimensionless)
    REAL(KIND=r8) , INTENT(INOUT) :: a10daylightu (npoi) ! 10-day average day-time PAR - upper canopy (micro-Ein m-2 s-1)
    REAL(KIND=r8) , INTENT(INOUT) :: a10scalparaml(npoi) ! 10-day average day-time scaling parameter - lower canopy (dimensionless)
    REAL(KIND=r8) , INTENT(INOUT) :: a10daylightl (npoi) ! 10-day average day-time PAR - lower canopy (micro-Ein m-2 s-1)
    REAL(KIND=r8) , INTENT(IN   ) :: tau15          ! co2/o2 specificity ratio at 15 degrees C (dimensionless)
    REAL(KIND=r8) , INTENT(IN   ) :: kc15           ! co2 kinetic parameter (mol/mol)
    REAL(KIND=r8) , INTENT(IN   ) :: ko15           ! o2 kinetic parameter (mol/mol) 
    REAL(KIND=r8) , INTENT(IN   ) :: cimax          ! maximum value for ci (needed for model stability)
    REAL(KIND=r8) , INTENT(IN   ) :: vmax_pft(npft) ! nominal vmax of top leaf at 15 C (mol-co2/m**2/s) [not used]
    REAL(KIND=r8) , INTENT(IN   ) :: gammaub        ! leaf respiration coefficient
    REAL(KIND=r8) , INTENT(IN   ) :: alpha3         ! intrinsic quantum efficiency for C3 plants (dimensionless)
    REAL(KIND=r8) , INTENT(IN   ) :: theta3         ! photosynthesis coupling coefficient for C3 plants (dimensionless)
    REAL(KIND=r8) , INTENT(IN   ) :: beta3          ! photosynthesis coupling coefficient for C3 plants (dimensionless)
    REAL(KIND=r8) , INTENT(IN   ) :: coefmub        ! 'm' coefficient for stomatal conductance relationship
    REAL(KIND=r8) , INTENT(IN   ) :: coefbub        ! 'b' coefficient for stomatal conductance relationship
    REAL(KIND=r8) , INTENT(IN   ) :: gsubmin        ! absolute minimum stomatal conductance
    REAL(KIND=r8) , INTENT(IN   ) :: gammauc        ! leaf respiration coefficient
    REAL(KIND=r8) , INTENT(IN   ) :: coefmuc        ! 'm' coefficient for stomatal conductance relationship  
    REAL(KIND=r8) , INTENT(IN   ) :: coefbuc        ! 'b' coefficient for stomatal conductance relationship  
    REAL(KIND=r8) , INTENT(IN   ) :: gsucmin        ! absolute minimum stomatal conductance
    REAL(KIND=r8) , INTENT(IN   ) :: gammals        ! leaf respiration coefficient
    REAL(KIND=r8) , INTENT(IN   ) :: coefmls        ! 'm' coefficient for stomatal conductance relationship 
    REAL(KIND=r8) , INTENT(IN   ) :: coefbls        ! 'b' coefficient for stomatal conductance relationship 
    REAL(KIND=r8) , INTENT(IN   ) :: gslsmin        ! absolute minimum stomatal conductance
    REAL(KIND=r8) , INTENT(IN   ) :: gammal3        ! leaf respiration coefficient
    REAL(KIND=r8) , INTENT(IN   ) :: coefml3        ! 'm' coefficient for stomatal conductance relationship 
    REAL(KIND=r8) , INTENT(IN   ) :: coefbl3        ! 'b' coefficient for stomatal conductance relationship 
    REAL(KIND=r8) , INTENT(IN   ) :: gsl3min        ! absolute minimum stomatal conductance
    REAL(KIND=r8) , INTENT(IN   ) :: gammal4        ! leaf respiration coefficient
    REAL(KIND=r8) , INTENT(IN   ) :: alpha4         ! intrinsic quantum efficiency for C4 plants (dimensionless)
    REAL(KIND=r8) , INTENT(IN   ) :: theta4         ! photosynthesis coupling coefficient for C4 plants (dimensionless) 
    REAL(KIND=r8) , INTENT(IN   ) :: beta4          ! photosynthesis coupling coefficient for C4 plants (dimensionless)
    REAL(KIND=r8) , INTENT(IN   ) :: coefml4        ! 'm' coefficient for stomatal conductance relationship
    REAL(KIND=r8) , INTENT(IN   ) :: coefbl4        ! 'b' coefficient for stomatal conductance relationship
    REAL(KIND=r8) , INTENT(IN   ) :: gsl4min        ! absolute minimum stomatal conductance      
    !
    ! local variables
    !
    INTEGER i
    REAL(KIND=r8) rwork  ! 3.47e-03 - 1. / tu(i)
    REAL(KIND=r8) tau    ! 
    REAL(KIND=r8) tleaf  ! leaf temp in celcius
    REAL(KIND=r8) tempvm(npoi) !
    REAL(KIND=r8) zweight!

    REAL(KIND=r8) esat12 ! vapor pressure in upper canopy air 
    REAL(KIND=r8) qsat12 ! specific humidity in upper canopy air
    REAL(KIND=r8) rh12   ! relative humidity in upper canopy air 
    REAL(KIND=r8) esat34 ! vapor pressure in lower canopy air
    REAL(KIND=r8) qsat34 ! specific humidity in lower canopy air 
    REAL(KIND=r8) rh34   ! relative humidity in lower canopy air 
    REAL(KIND=r8) gbco2u ! bound. lay. conductance for CO2 in upper canopy
    REAL(KIND=r8) gbco2l ! bound. lay. conductance for CO2 in lower canopy
    REAL(KIND=r8) gscub  ! 
    REAL(KIND=r8) gscuc  !
    REAL(KIND=r8) gscls  !
    REAL(KIND=r8) gscl3  !
    REAL(KIND=r8) gscl4  !
    REAL(KIND=r8) vmax (npoi)
    REAL(KIND=r8) vmaxub (npoi)
    REAL(KIND=r8) vmaxuc (npoi)
    REAL(KIND=r8) vmaxls (npoi)
    REAL(KIND=r8) vmaxl3 (npoi)
    REAL(KIND=r8) vmaxl4 (npoi)
    REAL(KIND=r8) rdarkub (npoi)
    REAL(KIND=r8) rdarkuc (npoi)
    REAL(KIND=r8) rdarkls (npoi)
    REAL(KIND=r8) rdarkl3 (npoi)
    REAL(KIND=r8) rdarkl4 (npoi)
    REAL(KIND=r8) agub (npoi)
    REAL(KIND=r8) aguc (npoi)
    REAL(KIND=r8) agls (npoi)
    REAL(KIND=r8) agl3 (npoi)
    REAL(KIND=r8) agl4 (npoi)
    REAL(KIND=r8) anub (npoi)
    REAL(KIND=r8) anuc (npoi)
    REAL(KIND=r8) anls (npoi)
    REAL(KIND=r8) anl3 (npoi)
    REAL(KIND=r8) anl4 (npoi)
    REAL(KIND=r8) duma 
    REAL(KIND=r8) dumb 
    REAL(KIND=r8) dumc 
    REAL(KIND=r8) dume 
    REAL(KIND=r8) dumq 
    REAL(KIND=r8) dump
    REAL(KIND=r8) pxaiu (npoi)
    REAL(KIND=r8) plaiu (npoi)
    REAL(KIND=r8) pxail (npoi)
    REAL(KIND=r8) plail (npoi)
    REAL(KIND=r8) cscub 
    REAL(KIND=r8) cscuc 
    REAL(KIND=r8) cscls 
    REAL(KIND=r8) cscl3 
    REAL(KIND=r8) cscl4
    REAL(KIND=r8) extpar 
    REAL(KIND=r8) scale
    !
    REAL(KIND=r8) kc     ! co2 kinetic parameter (mol/mol)
    REAL(KIND=r8) ko     ! o2  kinetic parameter (mol/mol)
    !*      REAL(KIND=r8) kc15
    !*      REAL(KIND=r8) ko15  ! o2  kinetic parameter (mol/mol) at 15 degrees C
    REAL(KIND=r8) kco2   ! initial c4 co2 efficiency (mol-co2/m**2/s)
    !REAL(KIND=r8) je     ! 'light limited' rate of photosynthesis (mol-co2/m**2/s)
    REAL(KIND=r8) je_ub (npoi)  ! 'light limited' rate of photosynthesis (mol-co2/m**2/s)
    REAL(KIND=r8) je_uc (npoi)  ! 'light limited' rate of photosynthesis (mol-co2/m**2/s)
    REAL(KIND=r8) je_ls (npoi)  ! 'light limited' rate of photosynthesis (mol-co2/m**2/s)
    REAL(KIND=r8) je_l3 (npoi)  ! 'light limited' rate of photosynthesis (mol-co2/m**2/s)
    REAL(KIND=r8) je_l4 (npoi)  ! 'light limited' rate of photosynthesis (mol-co2/m**2/s)

    !REAL(KIND=r8) jc     ! 'rubisco limited' rate of photosynthesis (mol-co2/m**2/s)
    REAL(KIND=r8) jc_ub (npoi)  ! 'rubisco limited' rate of photosynthesis (mol-co2/m**2/s)
    REAL(KIND=r8) jc_uc (npoi)  ! 'rubisco limited' rate of photosynthesis (mol-co2/m**2/s)
    REAL(KIND=r8) jc_ls (npoi)  ! 'rubisco limited' rate of photosynthesis (mol-co2/m**2/s)
    REAL(KIND=r8) jc_l3 (npoi)  ! 'rubisco limited' rate of photosynthesis (mol-co2/m**2/s)
    REAL(KIND=r8) jc_l4 (npoi)  ! 'rubisco limited' rate of photosynthesis (mol-co2/m**2/s)

    !REAL(KIND=r8) js     ! sucrose synthesis limitation
    REAL(KIND=r8) js_ub (npoi)  ! sucrose synthesis limitation
    REAL(KIND=r8) js_uc (npoi)  ! sucrose synthesis limitation
    REAL(KIND=r8) js_ls (npoi)  ! sucrose synthesis limitation
    REAL(KIND=r8) js_l3 (npoi)  ! sucrose synthesis limitation
    REAL(KIND=r8) js_l4 (npoi)  ! sucrose synthesis limitation

    !REAL(KIND=r8) ji     ! 'co2 limited' rate of photosynthesis (mol-co2/m**2/s)
    REAL(KIND=r8) ji_ub (npoi)  ! 'co2 limited' rate of photosynthesis (mol-co2/m**2/s)
    REAL(KIND=r8) ji_uc (npoi)  ! 'co2 limited' rate of photosynthesis (mol-co2/m**2/s)
    REAL(KIND=r8) ji_ls (npoi)  ! 'co2 limited' rate of photosynthesis (mol-co2/m**2/s)
    REAL(KIND=r8) ji_l3 (npoi)  ! 'co2 limited' rate of photosynthesis (mol-co2/m**2/s)
    REAL(KIND=r8) ji_l4 (npoi)  ! 'co2 limited' rate of photosynthesis (mol-co2/m**2/s)


    REAL(KIND=r8) jp     ! model-intermediate rate of photosynthesis (mol-co2/m**2/s)
    REAL(KIND=r8) gamstar! gamma*, the co2 compensation points for c3 plants
    !
    ! model parameters
    !
    ! intrinsic quantum efficiency for c3 and c4 plants (dimensionless)
    !
    !*      real alpha3, alpha4
    !
    !*      data alpha3 /0.060/
    !*      data alpha4 /0.050/
    !
    ! co2/o2 specificity ratio at 15 degrees C (dimensionless)
    !
    !*      real tau15
    !
    !*      data tau15 /4500.0/     
    !
    ! o2/co2 kinetic parameters (mol/mol)
    !
    !*      real kc15, ko12 
    !
    !*      data kc15 /1.5e-04/ 
    !*      data ko15 /2.5e-01/ 
    !
    ! leaf respiration coefficients
    !
    !*      real gammaub, gammauc, gammals, gammal3, gammal4
    !
    !*      data gammaub /0.0150/   ! broadleaf trees
    !*      data gammauc /0.0150/   ! conifer trees
    !*      data gammals /0.0150/   ! shrubs
    !*      data gammal3 /0.0150/   ! c3 grasses
    !*      data gammal4 /0.0300/   ! c4 grasses
    !
    ! 'm' coefficients for stomatal conductance relationship
    !
    !*      real coefmub, coefmuc, coefmls, coefml3, coefml4
    !
    !*      data coefmub /10.0/     ! broadleaf trees
    !*      data coefmuc / 6.0/     ! conifer trees
    !*      data coefmls / 9.0/     ! shrubs
    !*      data coefml3 / 9.0/     ! c3 grasses
    !*      data coefml4 / 4.0/     ! c4 grasses
    !
    ! 'b' coefficients for stomatal conductance relationship 
    ! (minimum conductance when net photosynthesis is zero)
    !
    !*      real coefbub, coefbuc, coefbls, coefbl3, coefbl4
    !
    !*      data coefbub /0.010/    ! broadleaf trees
    !*      data coefbuc /0.010/    ! conifer trees
    !*      data coefbls /0.010/    ! shrubs
    !*      data coefbl3 /0.010/    ! c3 grasses
    !*      data coefbl4 /0.040/    ! c4 grasses
    !
    ! absolute minimum stomatal conductances
    !
    !*      real gsubmin, gsucmin, gslsmin, gsl3min, gsl4min
    !
    !*      data gsubmin /0.00001/  ! broadleaf trees
    !*      data gsucmin /0.00001/  ! conifer trees
    !*      data gslsmin /0.00001/  ! shrubs
    !*      data gsl3min /0.00001/  ! c3 grasses
    !*      data gsl4min /0.00001/  ! c4 grasses
    !
    ! photosynthesis coupling coefficients (dimensionless)
    !
    !*      real theta3
    !
    !*      data theta3 /0.970/     ! c3 photosynthesis
    !
    !*      real theta4, beta4
    !
    !*      data theta4 /0.970/     ! c4 photosynthesis
    !*      data beta4  /0.800/     ! c4 photosynthesis
    !
    ! maximum values for ci (for model stability)
    !
    !*      real cimax
    !
    !*      data cimax /2000.e-06/  ! maximum values for ci
    !
    ! include water vapor functions
    !
    !
    !
    ! statement functions tsatl,tsati are used below so that lowe's
    ! polyomial for liquid is used if t gt 273.16, or for ice if 
    ! t lt 273.16. also impose range of validity for lowe's polys.
    !
    !      REAL(KIND=r8)    t        ! temperature argument of statement function 
    !      REAL(KIND=r8)    p1       ! pressure argument of function 
    !      REAL(KIND=r8)    e1       ! vapor pressure argument of function
    !      REAL(KIND=r8)    tsatl    ! statement function
    !      REAL(KIND=r8)    tsati    ! 
    !      REAL(KIND=r8)    esat     !
    !      REAL(KIND=r8)    qsat     ! 
    !      REAL(KIND=r8)    cvmgt    ! function
    !
    !      tsatl(t) = min (100., max (t-273.16, 0.))
    !      tsati(t) = max (-60., min (t-273.16, 0.))
    !
    ! statement function esat is svp in n/m**2, with t in deg k. 
    ! (100 * lowe's poly since 1 mb = 100 n/m**2.)
    !
    !      esat (t) =             &
    !       100.*(            &
    !              cvmgt (asat0, bsat0, t.ge.273.16)            &
    !              + tsatl(t)*(asat1 + tsatl(t)*(asat2 + tsatl(t)*(asat3            &
    !              + tsatl(t)*(asat4 + tsatl(t)*(asat5 + tsatl(t)* asat6))))) &
    !              + tsati(t)*(bsat1 + tsati(t)*(bsat2 + tsati(t)*(bsat3            &
    !              + tsati(t)*(bsat4 + tsati(t)*(bsat5 + tsati(t)* bsat6))))) &
    !       )
    !
    ! statement function qsat is saturation specific humidity,
    ! with svp e1 and ambient pressure p in n/m**2. impose an upper
    ! limit of 1 to avoid spurious values for very high svp
    ! and/or small p1
    !
    !       qsat (e1, p1) = 0.622 * e1 /  &
    !                     max ( p1 - (1.0 - 0.622) * e1, 0.622 * e1 )
    !
    !
    !
    !
    ! ---------------------------------------------------------------------
    ! * * * upper canopy physiology calculations * * *
    ! ---------------------------------------------------------------------
    !
    DO i = 1, npoi
       !
       ! calculate physiological parameter values which are a function of temperature
       !
       rwork = 3.47e-03_r8 - 1.0_r8 / MIN(MAX(tu(i),180.0_r8),360.0_r8)
       !
       tau = tau15 * EXP(-4500.0_r8 * rwork)
       kc  = kc15  * EXP( 6000.0_r8 * rwork)
       ko  = ko15  * EXP( 1500.0_r8 * rwork)
       !
       tleaf = tu(i) - 273.16_r8
       !
       tempvm(i) = EXP(3500.0_r8 * rwork ) /  &
            ((1.0_r8 + EXP(0.40_r8 * (  5.0_r8 - tleaf))) * &
             (1.0_r8 + EXP(0.40_r8 * (tleaf - 50.0_r8))))
       !
       ! upper canopy gamma-star values (mol/mol)
       ! is the compensation point for gross photosynthesis 
       !
       gamstar = o2conc / (2.0_r8 * tau)
       !
       ! constrain ci values to acceptable bounds -- to help ensure numerical stability
       !
       ciub(i) = MAX (1.05_r8 * gamstar, MIN (cimax, ciub(i)))
       ciuc(i) = MAX (1.05_r8 * gamstar, MIN (cimax, ciuc(i)))
       !
       ! calculate boundary layer parameters (mol/m**2/s) = su / 0.029 * 1.35
       !
       gbco2u = MIN (10.0_r8, MAX (0.1_r8, su(i) * 25.5_r8))
       ! 
       ! calculate the relative humidity in the canopy air space
       ! with a minimum value of 0.30_r8 to avoid errors in the 
       ! physiological calculations
       !
       esat12 = esat (t12(i))
       qsat12 = qsat (esat12, psurf(i))
       rh12   = MAX (0.30_r8, q12(i) / qsat12)
       !
       ! ---------------------------------------------------------------------
       ! broadleaf (evergreen & deciduous) tree physiology 
       ! ---------------------------------------------------------------------
       ! 
       ! nominal values for vmax of top leaf at 15 C (mol-co2/m**2/s)
       !
       ! tropical broadleaf trees          60.0 e-06 mol/m**2/sec
       ! warm-temperate broadleaf trees    40.0 e-06 mol/m**2/sec
       ! temperate broadleaf trees         25.0 e-06 mol/m**2/sec
       ! boreal broadleaf trees            25.0 e-06 mol/m**2/sec
       !_r8
       !*        if (exist(i,1).gt.0.5) then
       !*          vmaxub = 65.0e-06
       !*        else if (exist(i,3).gt.0.5) then
       !*          vmaxub = 40.0e-06
       !*        else 
       !*          vmaxub = 30.0e-06
       !*        endif
       !*
       !**** DTP 2001/06/06: Following code replaces above, making initialization
       !*                    dependent upon parameter values read in from external
       !*                    canopy parameter file "params.can".
       !*
       IF (exist(i,1).GT.0.5_r8) THEN
          vmaxub(i) = vmax_pft(1) ! 65.0e-06 ! Tropical broadleaf evergreen
       ELSE IF (exist(i,3).GT.0.5_r8) THEN
          vmaxub(i) = vmax_pft(3) ! 40.0e-06 ! Warm-temperate broadleaf evergreen
       ELSE 
          vmaxub(i) = vmax_pft(5) ! 30.0e-06 ! Temperate or boreal broadleaf cold deciduous
       END IF
       !
       ! vmax and dark respiration for current conditions
       !
       ! calculate (tempvm) physiological parameter values which are a function of temperature

       vmax(i)  = vmaxub(i) * tempvm(i) * stresstu(i)
       rdarkub(i) = gammaub * vmaxub(i) * tempvm(i)
       !
       ! 'light limited' rate of photosynthesis (mol/m**2/s)
       !  topparu is total photosynthetically active raditaion
       !   absorbed by top leaves of upper canopy (W m-2)
       !
       !   alpha3  instrinsic quantum efficiency of CO2 uptake in C3 plants 
       !
       je_ub(i) = topparu(i) * 4.59e-06_r8 * alpha3 * (ciub(i) - gamstar) / &
            (ciub(i) + 2.0_r8 * gamstar)
       !
       ! 'rubisco limited' rate of photosynthesis (mol/m**2/s)
       !
       jc_ub(i) = vmax(i) * (ciub(i) - gamstar) / &
            (ciub(i) + kc * (1.0_r8 + o2conc / ko))
       !
       ! solution to quadratic equation
       !
       !
       !       aX^2 + bX + c =0
       ! 
       !                    ________
       !                   /        \
       !            -b + \/     D
       !        X= --------------------
       !                  2a
       !
       !
       !      D = b^2 - 4ac 
       !
       duma = theta3
       dumb = je_ub(i) + jc_ub(i)
       dumc = je_ub(i) * jc_ub(i)
       !
       dume = MAX (dumb**2 - 4.0_r8 * duma * dumc, 0.0_r8)
       dumq = 0.5_r8 * (dumb + SQRT(dume)) + 1.e-15_r8
       !
       ! calculate the intermediate photosynthesis rate (mol/m**2/s)
       !
       jp = MIN (dumq/duma, dumc/dumq)
       !
       ! 'sucrose synthesis limited' rate of photosynthesis (mol/m**2/s)
       !
       js_ub(i) = vmax(i) / 2.2_r8
       !
       ! solution to quadratic equation
       !
       !
       !       aX^2 + bX + c =0
       ! 
       !                    ________
       !                   /        \
       !            -b + \/     D
       !        X= --------------------
       !                  2a
       !
       !
       !      D = b^2 - 4ac 
       !
       duma = beta3
       dumb = jp + js_ub(i)
       dumc = jp * js_ub(i)
       !
       dume = MAX (dumb**2 - 4.0_r8 * duma * dumc, 0.0_r8)
       !
       dumq = 0.5_r8 * (dumb + SQRT(dume)) + 1.e-15_r8
       !
       ! calculate the net photosynthesis rate (mol/m**2/s)
       !
       agub(i) = MIN (dumq/duma, dumc/dumq)
       anub(i) = agub(i) - rdarkub(i)
       !
       ! calculate co2 concentrations and stomatal condutance values
       ! using simple iterative procedure
       !
       ! weight results with the previous iteration's values -- this
       ! improves convergence by avoiding flip-flop between diffusion
       ! into and out of the stomatal cavities
       !
       ! calculate new value of cs using implicit scheme
       !
       csub(i) = 0.5_r8 * (csub(i) + co2conc(i) - anub(i) / gbco2u)
       csub(i) = MAX (1.05_r8 * gamstar, csub(i))
       !
       ! calculate new value of gs using implicit scheme
       !
       gsub(i) = 0.5_r8 * (gsub(i)  +  (coefmub * anub(i) * rh12 / csub(i) + coefbub * stresstu(i)))
       !
       gsub(i) = MAX (gsubmin, coefbub * stresstu(i), gsub(i))
       !
       ! calculate new value of ci using implicit scheme
       !
       ciub(i) = 0.5_r8 * (ciub(i) + csub(i) - 1.6_r8 * anub(i) / gsub(i))
       ciub(i) = MAX (1.05_r8 * gamstar, MIN (cimax, ciub(i)))
       !
       ! ---------------------------------------------------------------------
       ! conifer tree physiology 
       ! ---------------------------------------------------------------------
       ! 
       ! nominal values for vmax of top leaf at 15 C (mol-co2/m**2/s)
       !
       ! temperate conifer trees           30.0 e-06 mol/m**2/sec
       ! boreal conifer trees              20.0 e-06 mol/m**2/sec
       !
       !*        if (exist(i,4).gt.0.5) then
       !*          vmaxuc = 30.0e-06
       !*        else 
       !*          vmaxuc = 25.0e-06
       !*        endif
       !*
       !**** DTP 2001/06/06: Following code replaces above, making initialization
       !*                    dependent upon parameter values read in from external
       !*                    canopy parameter file "params.can".
       !*
       IF (exist(i,4).GT.0.5_r8) THEN
          vmaxuc(i) = vmax_pft(4) ! 30.0e-06 ! Temperate conifer
       ELSE 
          vmaxuc(i) = vmax_pft(6) ! 25.0e-06 ! Boreal conifer evergreen
       END IF

       !
       ! vmax and dark respiration for current conditions
       !
       vmax(i)  = vmaxuc(i) * tempvm(i) * stresstu(i)
       !gammauc         ! leaf respiration coefficient
       rdarkuc(i) = gammauc * vmaxuc(i) * tempvm(i)
       !
       ! 'light limited' rate of photosynthesis (mol/m**2/s)
       !
       je_uc(i) = topparu(i) * 4.59e-06_r8 * alpha3 * (ciuc(i) - gamstar) / &
            (ciuc(i) + 2.0_r8 * gamstar)
       !
       ! 'rubisco limited' rate of photosynthesis (mol/m**2/s)
       !
       jc_uc(i) = vmax(i) * (ciuc(i) - gamstar) /   &
            (ciuc(i) + kc * (1.0_r8 + o2conc / ko))
       !
       ! solution to quadratic equation
       !
       !
       !       aX^2 + bX + c =0
       ! 
       !                    ________
       !                   /        \
       !            -b + \/     D
       !        X= --------------------
       !                  2a
       !
       !
       !      D = b^2 - 4ac 
       !
       duma = theta3
       dumb = je_uc(i) + jc_uc(i)
       dumc = je_uc(i) * jc_uc(i)
       !
       dume = MAX (dumb**2 - 4.0_r8 * duma * dumc, 0.0_r8)
       dumq = 0.5_r8 * (dumb + SQRT(dume)) + 1.e-15_r8
       !
       ! calculate the intermediate photosynthesis rate (mol/m**2/s)
       !
       jp = MIN (dumq/duma, dumc/dumq)
       !
       ! 'sucrose synthesis limited' rate of photosynthesis (mol/m**2/s)
       !
       js_uc(i) = vmax(i) / 2.2_r8
       !
       ! solution to quadratic equation
       !
       !
       !       aX^2 + bX + c =0
       ! 
       !                    ________
       !                   /        \
       !            -b + \/     D
       !        X= --------------------
       !                  2a
       !
       !
       !      D = b^2 - 4ac 
       !
       duma = beta3
       dumb = jp + js_uc(i)
       dumc = jp * js_uc(i)
       !
       dume = MAX (dumb**2 - 4.0_r8 * duma * dumc, 0.0_r8)
       dumq = 0.5_r8 * (dumb + SQRT(dume)) + 1.e-15_r8
       !
       ! calculate the net photosynthesis rate (mol/m**2/s)
       !
       aguc(i) = MIN (dumq/duma, dumc/dumq) 
       anuc(i) = aguc(i) - rdarkuc(i)
       !
       ! calculate co2 concentrations and stomatal condutance values
       ! calculate co2 concentrations and stomatal condutance values
       ! using simple iterative procedure
       !
       ! weight results with the previous iteration's values -- this
       ! improves convergence by avoiding flip-flop between diffusion
       ! into and out of the stomatal cavities
       !
       ! calculate new value of cs using implicit scheme
       ! csuc  leaf boundary layer co2 concentration - conifer   (mol_co2/mol_air)
       ! gbco2u  ->  boundary layer parameters (mol/m**2/s) = su / 0.029 * 1.35
       !
       ! upper canopy gamma-star values (mol/mol)
       !
       ! gamstar = o2conc / (2.0_r8 * tau)

       csuc(i) = 0.5_r8 * (csuc(i) + co2conc(i) - anuc(i) / gbco2u)
       csuc(i) = MAX (1.05_r8 * gamstar, csuc(i))
       !
       ! calculate new value of gs using implicit scheme
       ! upper canopy stomatal conductance - conifer    (mol_co2 m-2 s-1)
       gsuc(i) = 0.5_r8 * (gsuc(i)  +  (coefmuc * anuc(i) * rh12 / csuc(i) +  &
            coefbuc * stresstu(i)))
       !
       gsuc(i) = MAX (gsucmin, coefbuc * stresstu(i), gsuc(i))
       !
       ! calculate new value of ci using implicit scheme
       !
       ciuc(i) = 0.5_r8 * (ciuc(i) + csuc(i) - 1.6_r8 * anuc(i) / gsuc(i))
       ciuc(i) = MAX (1.05_r8 * gamstar, MIN (cimax, ciuc(i)))
       !
       ! ---------------------------------------------------------------------
       ! upper canopy scaling
       ! ---------------------------------------------------------------------
       !
       ! the canopy scaling algorithm assumes that the net photosynthesis
       ! is proportional to absored par (apar) during the daytime. during night,
       ! the respiration is scaled using a 10-day running-average daytime canopy
       ! scaling parameter.
       !
       ! apar(x) = A exp(-k x) + B exp(-h x) + C exp(h x)
       ! an(x) is proportional to apar(x)
       !
       ! therefore, an(x) = an(0) * apar(x) / apar(0)
       ! an(x) = an(0) * (A exp(-k x) + B exp(-h x) + C exp(h x)) / 
       !                 (A + B + C)
       !
       ! this equation is further simplified to
       ! an(x) = an(0) * exp (-extpar * x)
       !
       ! an(0) is calculated for a sunlit leaf at the top of the canopy using
       ! the full-blown plant physiology model (Farquhar/Ball&Berry, Collatz).
       ! then the approximate par extinction coefficient (extpar) is calculated
       ! using parameters obtained from the two-stream radiation calculation.
       !
       ! an,canopy avg.= integral (an(x), from 0 to xai) / lai
       !               = an(0) * (1 - exp (-extpar * xai )) / (extpar * lai)
       !
       ! the term '(1 - exp (-extpar * xai )) / lai)' scales photosynthesis from leaf
       ! to canopy level (canopy average) at day time. A 10-day running mean of this
       ! scaling parameter (weighted by light) is then used to scale the respiration
       ! during night time.
       !
       ! once canopy average photosynthesis is calculated, then the canopy average
       ! stomatal conductance is calculated using the 'big leaf approach',i.e. 
       ! assuming that the canopy is a big leaf and applying the leaf-level stomatal
       ! conductance equations to the whole canopy.
       !
       ! calculate the approximate par extinction coefficient:
       !
       ! extpar = (k * A + h * B - h * C) / (A + B + C)
       !
       !WRITE(*,*)termu(i,6),termu(i,7),scalcoefu(i,3),scalcoefu(i,2),scalcoefu(i,1)
       extpar = (termu(i,6) * scalcoefu(i,1) +   &
            termu(i,7) * scalcoefu(i,2) -   &
            termu(i,7) * scalcoefu(i,3)) /  &
            MAX (scalcoefu(i,4), epsilon)
       !
       extpar = MAX (1.e-1_r8, MIN (1.e+1_r8, extpar))
       !
       ! calculate canopy average photosynthesis (per unit leaf area):
       !
       pxaiu(i) = extpar * (lai(i,2) + sai(i,2))
       plaiu(i) = extpar *  lai(i,2)
       !
       ! scale is the parameter that scales from leaf-level photosynthesis to
       ! canopy average photosynthesis
       ! CD : replaced 24 (hours) by 86400/dtime for use with other timestep
       !
       zweight = EXP(-1.0_r8 / (10.0_r8 * 86400.0_r8 / dtime))
       !  zweight = exp(-1.0_r8 / (10.0_r8 * 86400.0_r8 / 2400))
       !
       ! for non-zero lai
       !
       IF (plaiu(i).GT.0.0_r8) THEN
          !
          ! day-time conditions, use current scaling coefficient
          !
          IF (topparu(i).GT.10.0_r8) THEN
             !
             scale = (1.0_r8 - EXP(-pxaiu(i))) / plaiu(i)
             !
             ! update 10-day running mean of scale, weighted by light levels
             !
             a10scalparamu(i) = zweight * a10scalparamu(i) + &
                  (1.0_r8 - zweight) * scale * topparu(i)
             !
             a10daylightu(i)  = zweight * a10daylightu(i) +   &
                  (1.0_r8 - zweight) * topparu(i)
             !
             ! night-time conditions, use long-term day-time average scaling coefficient
             !
          ELSE
             !
             scale = a10scalparamu(i) / a10daylightu(i)
             !
          END IF
          !
          ! if no lai present
          !
       ELSE
          !
          scale = 0.0_r8
          !
       END IF
       !
       ! perform scaling on all carbon fluxes from upper canopy
       !
       agcub(i) = agub(i) * scale
       agcuc(i) = aguc(i) * scale
       !
       ancub(i) = anub(i) * scale
       ancuc(i) = anuc(i) * scale
       !
       ! calculate diagnostic canopy average surface co2 concentration 
       ! (big leaf approach)
       !
       cscub = MAX (1.05_r8 * gamstar, co2conc(i) - ancub(i) / gbco2u)
       cscuc = MAX (1.05_r8 * gamstar, co2conc(i) - ancuc(i) / gbco2u)
       !
       ! calculate diagnostic canopy average stomatal conductance (big leaf approach)
       !
       gscub = coefmub * ancub(i) * rh12 / cscub + coefbub * stresstu(i)
       !
       gscuc = coefmuc * ancuc(i) * rh12 / cscuc + coefbuc * stresstu(i)
       !
       gscub = MAX (gsubmin, coefbub * stresstu(i), gscub)
       gscuc = MAX (gsucmin, coefbuc * stresstu(i), gscuc)
       !
       ! calculate total canopy and boundary-layer total conductance for 
       ! water vapor diffusion
       !
       rwork = 1.0_r8 / su(i)
       dump  = 1.0_r8 / 0.029_r8
       !
       totcondub(i) = 1.0_r8 / (rwork + dump / gscub)
       totconduc(i) = 1.0_r8 / (rwork + dump / gscuc)
       !
       ! multiply canopy photosynthesis by wet fraction - this calculation is
       ! done here and not earlier to avoid using within canopy conductance
       !
       rwork = 1 - fwetu(i)
       !
       agcub(i) = rwork * agcub(i)
       agcuc(i) = rwork * agcuc(i)
       !
       ancub(i) = rwork * ancub(i)
       ancuc(i) = rwork * ancuc(i)
       !
    END DO
    !
    ! ---------------------------------------------------------------------
    ! * * * lower canopy physiology calculations * * *
    ! ---------------------------------------------------------------------
    !
    DO i = 1, npoi
       !
       ! calculate physiological parameter values which are a function of temperature
       !
       rwork = 3.47e-03_r8 - 1.0_r8 / MIN(MAX(tl(i),180.0_r8),360.0_r8)
       !rwork = 3.47e-03_r8 - 1.0_r8 / tl(i)
       !WRITE(*,*)tl(i),rwork
       !
       tau = tau15 * EXP(-4500.0_r8 * rwork)
       kc  = kc15  * EXP( 6000.0_r8 * rwork)
       ko  = ko15  * EXP( 1500.0_r8 * rwork)
       !
       tleaf = tl(i) - 273.16_r8
       !
       tempvm(i) = EXP(3500.0_r8 * rwork ) /  &
            ((1.0_r8 + EXP(0.40_r8 * (  5.0_r8 - tleaf))) *  &
            (1.0_r8 + EXP(0.40_r8 * (tleaf - 50.0_r8))))
       !
       ! lower canopy gamma-star values (mol/mol)
       !
       gamstar = o2conc / (2.0_r8 * tau)
       !
       ! constrain ci values to acceptable bounds -- to help ensure numerical stability
       !
       cils(i) = MAX (1.05_r8 * gamstar, MIN (cimax, cils(i)))
       cil3(i) = MAX (1.05_r8 * gamstar, MIN (cimax, cil3(i)))
       cil4(i) = MAX (0.0_r8           , MIN (cimax, cil4(i)))
       !
       ! calculate boundary layer parameters (mol/m**2/s) = su / 0.029 * 1.35
       !
       gbco2l = MIN (10.0_r8, MAX (0.1_r8, sl(i) * 25.5_r8))
       ! 
       ! calculate the relative humidity in the canopy air space
       ! with a minimum value of 0.30_r8 to avoid errors in the 
       ! physiological calculations
       !
       esat34 = esat (t34(i))
       qsat34 = qsat (esat34, psurf(i))
       rh34   = MAX (0.30_r8, q34(i) / qsat34)
       !
       ! ---------------------------------------------------------------------
       ! shrub physiology
       ! ---------------------------------------------------------------------
       ! 
       ! nominal values for vmax of top leaf at 15 C (mol-co2/m**2/s)
       !
       !*        vmaxls = 27.5e-06_r8
       !*
       !**** DTP 2001/06/06: Following code replaces above, making initialization
       !*                    dependent upon parameter values read in from external
       !*                    canopy parameter file "params.can".
       !*
       vmaxls(i) = vmax_pft(9) ! 27.5e-06 ! Shrubs (evergreen or cold deciduous) 
       ! 
       ! vmax and dark respiration for current conditions
       !
       vmax(i)  = vmaxls(i) * tempvm(i) * stresstl(i)
       !   gammals          leaf respiration coefficient
       rdarkls(i) = gammals * vmaxls(i) * tempvm(i)
       !
       ! 'light limited' rate of photosynthesis (mol/m**2/s)
       !
       je_ls(i) = topparl(i) * 4.59e-06_r8 * alpha3 * (cils(i) - gamstar) / &
            (cils(i) + 2.0_r8 * gamstar)
       !
       ! 'rubisco limited' rate of photosynthesis (mol/m**2/s)
       !
       jc_ls(i) = vmax(i) * (cils(i) - gamstar) /  &
            (cils(i) + kc * (1.0_r8 + o2conc / ko))
       !
       ! solution to quadratic equation
       !
       duma = theta3
       dumb = je_ls(i) + jc_ls(i)
       dumc = je_ls(i) * jc_ls(i)
       !
       dume = MAX (dumb**2 - 4.0_r8 * duma * dumc, 0.0_r8)
       dumq = 0.5_r8 * (dumb + SQRT(dume)) + 1.e-15_r8
       !
       ! calculate the intermediate photosynthesis rate (mol/m**2/s)
       !
       jp = MIN (dumq/duma, dumc/dumq)
       !
       ! 'sucrose synthesis limited' rate of photosynthesis (mol/m**2/s)
       !
       js_ls(i) = vmax(i) / 2.2_r8
       !
       ! solution to quadratic equation
       !
       !
       !       aX^2 + bX + c =0
       ! 
       !                    ________
       !                   /        \
       !            -b + \/     D
       !        X= --------------------
       !                  2a
       !
       !
       !      D = b^2 - 4ac 
       !
       duma = beta3
       dumb = jp + js_ls(i)
       dumc = jp * js_ls(i)
       !
       dume = MAX (dumb**2 - 4.0_r8 * duma * dumc, 0.0_r8)
       dumq = 0.5_r8 * (dumb + SQRT(dume)) + 1.e-15_r8
       !
       ! calculate the net photosynthesis rate (mol/m**2/s)
       !
       agls(i) = MIN (dumq/duma, dumc/dumq)
       anls(i) = agls(i) - rdarkls(i)
       !
       ! calculate co2 concentrations and stomatal condutance values
       ! using simple iterative procedure
       !
       ! weight results with the previous iteration's values -- this
       ! improves convergence by avoiding flip-flop between diffusion
       ! into and out of the stomatal cavities
       !
       ! calculate new value of cs using implicit scheme
       !
       csls(i) = 0.5_r8 * (csls(i) + co2conc(i) - anls(i) / gbco2l)
       csls(i) = MAX (1.05_r8 * gamstar, csls(i))
       !
       ! calculate new value of gs using implicit scheme
       !
       gsls(i) = 0.5_r8 * (gsls(i) + coefmls * anls(i) * rh34 / csls(i) + &
            coefbls * stresstl(i))
       !
       gsls(i) = MAX (gslsmin, coefbls * stresstl(i), gsls(i))
       !
       ! calculate new value of ci using implicit scheme
       !
       cils(i) = 0.5_r8 * (cils(i) + csls(i) - 1.6_r8 * anls(i) / gsls(i))
       cils(i) = MAX (1.05_r8 * gamstar, MIN (cimax, cils(i)))
       !
       ! ---------------------------------------------------------------------
       ! c3 grass physiology
       ! ---------------------------------------------------------------------
       ! 
       ! nominal values for vmax of top leaf at 15 C (mol-co2/m**2/s)
       !
       !*       vmaxl3 = 25.0e-06
       !*
       !**** DTP 2001/06/06: Following code replaces above, making initialization
       !*                    dependent upon parameter value read in from external
       !*                    canopy parameter file "params.can".
       !*
       vmaxl3(i) = vmax_pft(12) ! 25.0e-06 ! C3 grasses
       ! 
       ! vmax and dark respiration for current conditions
       !
       vmax(i)  = vmaxl3(i) * tempvm(i) * stresstl(i)
       ! gammal3   leaf respiration coefficient
       rdarkl3(i) = gammal3 * vmaxl3(i) * tempvm(i)
       !
       ! 'light limited' rate of photosynthesis (mol/m**2/s)
       !
       je_l3(i) = topparl(i) * 4.59e-06_r8 * alpha3 * (cil3(i) - gamstar) /   &
            (cil3(i) + 2.0_r8 * gamstar)
       !
       ! 'rubisco limited' rate of photosynthesis (mol/m**2/s)
       !
       jc_l3(i) = vmax(i) * (cil3(i) - gamstar) /      &
            (cil3(i) + kc * (1.0_r8 + o2conc / ko))
       !
       ! solution to quadratic equation
       !
       duma = theta3
       dumb = je_l3(i) + jc_l3(i)
       dumc = je_l3(i) * jc_l3(i)
       !
       dume = MAX (dumb**2 - 4.0_r8 * duma * dumc, 0.0_r8)
       dumq = 0.5_r8 * (dumb + SQRT(dume)) + 1.e-15_r8
       !
       ! calculate the intermediate photosynthesis rate (mol/m**2/s)
       !
       jp = MIN (dumq/duma, dumc/dumq)
       !
       ! 'sucrose synthesis limited' rate of photosynthesis (mol/m**2/s)
       !
       js_l3(i) = vmax(i) / 2.2_r8
       !
       ! solution to quadratic equation
       !
       !
       !       aX^2 + bX + c =0
       ! 
       !                    ________
       !                   /        \
       !            -b + \/     D
       !        X= --------------------
       !                  2a
       !
       !
       !      D = b^2 - 4ac 
       !
       duma = beta3
       dumb = jp + js_l3(i)
       dumc = jp * js_l3(i)
       !
       dume = MAX (dumb**2 - 4.0_r8 * duma * dumc, 0.0_r8)
       dumq = 0.5_r8 * (dumb + SQRT(dume)) + 1.e-15_r8
       !
       ! calculate the net photosynthesis rate (mol/m**2/s)
       !
       agl3(i) = MIN (dumq/duma, dumc/dumq)
       anl3(i) = agl3(i) - rdarkl3(i)
       !
       ! calculate co2 concentrations and stomatal condutance values
       ! using simple iterative procedure
       !
       ! weight results with the previous iteration's values -- this
       ! improves convergence by avoiding flip-flop between diffusion
       ! into and out of the stomatal cavities
       !
       ! calculate new value of cs using implicit scheme
       !
       csl3(i) = 0.5_r8 * (csl3(i) + co2conc(i) - anl3(i) / gbco2l)
       csl3(i) = MAX (1.05_r8 * gamstar, csl3(i))
       !
       ! calculate new value of gs using implicit scheme
       !
       gsl3(i) = 0.5_r8 * (gsl3(i) + coefml3 * anl3(i) * rh34 / csl3(i) +  &
            coefbl3 * stresstl(i))
       !
       gsl3(i) = MAX (gsl3min, coefbl3 * stresstl(i), gsl3(i))
       !
       ! calculate new value of ci using implicit scheme
       !
       cil3(i) = 0.5_r8 * (cil3(i) + csl3(i) - 1.6_r8 * anl3(i) / gsl3(i))
       cil3(i) = MAX (1.05_r8 * gamstar, MIN (cimax, cil3(i)))
       !
       ! ---------------------------------------------------------------------
       ! c4 grass physiology
       ! ---------------------------------------------------------------------
       !
       ! nominal values for vmax of top leaf at 15 C (mol-co2/m**2/s)
       !
       !*       vmaxl4 = 15.0e-06
       !*
       !**** DTP 2001/06/06: Following code replaces above, making initialization
       !*                    dependent upon parameter value read in from external
       !*                    canopy parameter file "params.can".
       !*
       vmaxl4(i) = vmax_pft(11) ! 15.0e-06 ! C4 grasses
       !
       ! calculate the parameter values which are a function of temperature
       !
       rwork = 3.47e-03_r8 - 1.0_r8 / MIN(MAX(tl(i),180.0_r8),360.0_r8)
       !rwork = 3.47e-03_r8 - 1.0_r8 / tl(i)
       !
       tleaf = tl(i) - 273.16_r8
       !
       tempvm(i) = EXP(3500.0_r8 * rwork ) /            &
            ((1.0_r8 + EXP(0.40_r8 * ( 10.0_r8 - tleaf))) * & 
            (1.0_r8 + EXP(0.40_r8 * (tleaf - 50.0_r8))))
       !
       ! vmax and dark respiration for current conditions
       !
       vmax(i)  = vmaxl4(i) * tempvm(i) * stresstl(i)
       ! gammal4   leaf respiration coefficient
       rdarkl4(i) = gammal4 * vmaxl4(i) * tempvm(i)
       !
       ! initial c4 co2 efficiency (mol/m**2/s)
       !
       kco2 = 18.0e+03_r8 * vmax(i)
       !
       ! 'light limited' rate of photosynthesis (mol/m**2/s)
       !
       je_l4(i) = topparl(i) * 4.59e-06_r8 * alpha4
       !
       ! 'rubisco limited' rate of photosynthesis
       !
       jc_l4(i) = vmax(i)
       !
       ! solve for intermediate photosynthesis rate
       !
       !
       !       aX^2 + bX + c =0
       ! 
       !                    ________
       !                   /        \
       !            -b + \/     D
       !        X= --------------------
       !                  2a
       !
       !
       !      D = b^2 - 4ac 
       !
       duma = theta4
       dumb = je_l4(i) + jc_l4(i)
       dumc = je_l4(i) * jc_l4(i)
       !
       dume = MAX (dumb**2 - 4.0_r8 * duma * dumc, 0.0_r8)
       dumq = 0.5_r8 * (dumb + SQRT(dume)) + 1.e-15_r8
       !
       jp = MIN (dumq/duma, dumc/dumq)
       !
       ! 'carbon dioxide limited' rate of photosynthesis (mol/m**2/s)
       !
       ji_l4(i) = kco2 * cil4(i)
       !
       ! solution to quadratic equation
       !
       !       aX^2 + bX + c =0
       ! 
       !                    ________
       !                   /        \
       !            -b + \/     D
       !        X= --------------------
       !                  2a
       !
       !
       !      D = b^2 - 4ac 
       !
       duma = beta4
       dumb = jp + ji_l4(i)
       dumc = jp * ji_l4(i)
       !
       dume = MAX (dumb**2 - 4.0_r8 * duma * dumc, 0.0_r8)
       dumq = 0.5_r8 * (dumb + SQRT(dume)) + 1.e-15_r8
       !
       ! calculate the net photosynthesis rate (mol/m**2/s)
       !
       agl4(i) = MIN (dumq/duma, dumc/dumq)
       anl4(i) = agl4(i) - rdarkl4(i)
       !
       ! calculate co2 concentrations and stomatal condutance values
       ! using simple iterative procedure
       !
       ! weight results with the previous iteration's values -- this
       ! improves convergence by avoiding flip-flop between diffusion
       ! into and out of the stomatal cavities
       !
       ! calculate new value of cs using implicit scheme
       ! CD: For numerical stability (to avoid division by zero in gsl4), 
       ! csl4 is limited to 1e-8 mol_co2/mol_air.
       !  
       csl4(i) = 0.5_r8 * (csl4(i) + co2conc(i) - anl4(i) / gbco2l)
       csl4(i) = MAX (1.e-8_r8, csl4(i))
       !
       ! calculate new value of gs using implicit scheme
       !
       gsl4(i) = 0.5_r8 * (gsl4(i) + coefml4 * anl4(i) * rh34 / csl4(i) +  &
            coefbl4 * stresstl(i))
       !
       gsl4(i) = MAX (gsl4min, coefbl4 * stresstl(i), gsl4(i))
       !
       ! calculate new value of ci using implicit scheme
       !
       cil4(i) = 0.5_r8 * (cil4(i) + csl4(i) - 1.6_r8 * anl4(i) / gsl4(i))
       cil4(i) = MAX (0.0_r8, MIN (cimax, cil4(i)))
       !
       ! ---------------------------------------------------------------------
       ! lower canopy scaling
       ! ---------------------------------------------------------------------
       !
       ! calculate the approximate extinction coefficient
       !
       extpar = (terml(i,6) * scalcoefl(i,1) +   &
            terml(i,7) * scalcoefl(i,2) -   &
            terml(i,7) * scalcoefl(i,3)) /  &
            MAX (scalcoefl(i,4), epsilon)
       !
       extpar = MAX (1.e-1_r8, MIN (1.e+1_r8, extpar))
       !
       ! calculate canopy average photosynthesis (per unit leaf area):
       !
       pxail(i) = extpar * (lai(i,1) + sai(i,1))
       plail(i) = extpar *  lai(i,1)
       !
       ! scale is the parameter that scales from leaf-level photosynthesis to
       ! canopy average photosynthesis
       ! CD : replaced 24 (hours) by 86400/dtime for use with other timestep
       !
       zweight = EXP(-1.0_r8 / (10.0_r8 * 86400.0_r8 / dtime))
       ! zweight = exp(-1.0_r8 / (10.0_r8 * 86400.0_r8 / 2400))
       !
       ! for non-zero lai
       !
       IF (plail(i).GT.0.0_r8) THEN
          !
          ! day-time conditions, use current scaling coefficient
          !
          IF (topparl(i).GT.10.0_r8) THEN
             !
             scale = (1.0_r8 - EXP(-pxail(i))) / plail(i)
             !
             ! update 10-day running mean of scale, weighted by light levels
             !
             a10scalparaml(i) = zweight * a10scalparaml(i) +    &
                  (1.0_r8 - zweight) * scale * topparl(i)
             !
             a10daylightl(i)  = zweight * a10daylightl(i) +  &
                  (1.0_r8 - zweight) * topparl(i)
             !
             ! night-time conditions, use long-term day-time average scaling coefficient
             !
          ELSE
             !
             scale = a10scalparaml(i) / a10daylightl(i)
             !
          END IF
          !
          ! if no lai present
          !
       ELSE
          !
          scale = 0.0_r8
          !
       END IF
       !       WRITE(*,*) scale,agls,agl4,agl3,anl4,anls

       !
       ! perform scaling on all carbon fluxes from upper canopy
       !
       agcls(i) = agls(i) * scale
       agcl4(i) = agl4(i) * scale
       agcl3(i) = agl3(i) * scale
       !
       ancls(i) = anls(i) * scale
       ancl4(i) = anl4(i) * scale
       ancl3(i) = anl3(i) * scale

       !
       ! calculate canopy average surface co2 concentration
       ! CD: For numerical stability (to avoid division by zero in gscl4),
       ! cscl4 is limited to 1e-8 mol_co2/mol_air.
       !
       cscls = MAX (1.05_r8 * gamstar, co2conc(i) - ancls(i) / gbco2l)
       cscl3 = MAX (1.05_r8 * gamstar, co2conc(i) - ancl3(i) / gbco2l)
       cscl4 = MAX (1.e-8_r8         , co2conc(i) - ancl4(i) / gbco2l)
       !
       ! calculate canopy average stomatal conductance
       !
       gscls = coefmls * ancls(i) * rh34 / cscls +  &
            coefbls * stresstl(i)
       !
       gscl3 = coefml3 * ancl3(i) * rh34 / cscl3 +   &
            coefbl3 * stresstl(i)
       !
       gscl4 = coefml4 * ancl4(i) * rh34 / cscl4 +   &
            coefbl4 * stresstl(i)
       !
       gscls = MAX (gslsmin, coefbls * stresstl(i), gscls)
       gscl3 = MAX (gsl3min, coefbl3 * stresstl(i), gscl3)
       gscl4 = MAX (gsl4min, coefbl4 * stresstl(i), gscl4)
       !
       ! calculate canopy and boundary-layer total conductance for water vapor diffusion
       !
       rwork = 1.0_r8 / sl(i)
       dump =  1.0_r8 / 0.029_r8
       !
       totcondls(i) = 1.0_r8 / (rwork + dump / gscls)
       totcondl3(i) = 1.0_r8 / (rwork + dump / gscl3)
       totcondl4(i) = 1.0_r8 / (rwork + dump / gscl4)
       !
       ! multiply canopy photosynthesis by wet fraction -- this calculation is
       ! done here and not earlier to avoid using within canopy conductance
       !
       rwork = 1.0_r8 - fwetl(i)
       !
       agcls(i) = rwork * agcls(i)
       agcl3(i) = rwork * agcl3(i)
       agcl4(i) = rwork * agcl4(i)
       !
       ancls(i) = rwork * ancls(i)
       ancl3(i) = rwork * ancl3(i)
       ancl4(i) = rwork * ancl4(i)
       !
    END DO
    !
    ! return to main program
    !
    RETURN
  END SUBROUTINE stomata


  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE co2 (npoi,co2initm, co2conc, iyear)
    ! 1.1 MIXING RATIO
    ! The mixing ratio CX of a gas X (equivalently called the mole fraction) is defined as the 
    ! number of moles of X per mole of air. It is given in units of mol/mol (abbreviation for 
    ! moles per mole), or equivalently in units of v/v (volume of gas per volume of air) since 
    ! the volume occupied by an ideal gas is proportional to the number of molecules. Pressures
    ! in the atmosphere are sufficiently low that the ideal gas law is always obeyed to within 1%.
    ! The mixing ratio of a gas has the virtue of remaining constant when the air density changes
    ! (as happens when the temperature or the pressure changes). Consider a balloon filled with
    ! room air and allowed to rise in the atmosphere. As the balloon rises it expands, so that 
    ! the number of molecules per unit volume inside the balloon decreases; however, the mixing
    ! ratios of the different gases in the balloon remain constant. The mixing ratio is therefore
    ! a robust measure of atmospheric composition.
  
    ! Table 1-1 lists the mixing ratios of some major atmospheric gases. The most abundant is 
    ! molecular nitrogen (N2) with a mixing ratio CN2 = 0.78 mol/mol; N2 accounts for 78% of 
    ! all molecules in the atmosphere. Next in abundance are molecular oxygen (O2) 
    ! with CO2 = 0.21 mol/mol, and argon (Ar) with CAr = 0.0093 mol/mol. The mixing ratios in 
    ! Table 1-1 are for dry air, excluding water vapor. Water vapor mixing ratios in the atmosphere
    ! are highly variable (10-6-10-2 mol/mol). This variability in water vapor is part of our 
    ! everyday experience as it affects the ability of sweat to evaporate and the drying rate of
    ! clothes on a line.
 
    ! Gases other than N2, O2, Ar, and H2O are present in the atmosphere at extremely low 
    ! concentrations and are called trace gases. Despite their low concentrations, 
    ! these trace gases can be of critical importance for the greenhouse effect, the ozone layer, 
    ! smog, and other environmental issues. Mixing ratios of trace gases are commonly given in 
    ! units of parts per million volume ( ppmv or simply ppm), parts per billion volume ( ppbv or ppb), 
    ! or parts per trillion volume ( pptv or ppt); 
    ! 1 ppmv = 1x10-6 mol/mol, 
    ! 1 ppbv = 1x10-9 mol/mol, and 
    ! 1 pptv = 1x10-12 mol/mol. 
    !
    ! For example, the present-day CO2 concentration is 365 ppmv (365x10-6 mol/mol).
    ! ---------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    ! Arguments 
    !
    INTEGER, INTENT(IN   ) :: npoi,iyear    ! current year
    !
    REAL(KIND=r8), INTENT(IN   ) :: co2initm(npoi)  ! input atmospheric co2 concentration
    REAL(KIND=r8), INTENT(OUT  ) :: co2conc(npoi)  ! output " for year iyear   
    INTEGER :: iyr
    iyr=iyear
    !
    ! calculate co2 concentration for this year
    !
    !     if (iyear.lt.1860) then
    !
    co2conc = co2initm
    !
    !     else
    !
    ! 1992 IPCC estimates
    !
    !       iyr = iyear - 1860 + 1
    !       co2conc = (297.12 - 0.26716 * iyr +
    !    >                      0.0015368 * iyr**2 +
    !    >                      3.451e-5 * iyr**3) * 1.e-6
    !
    !
    ! M. El Maayar: 1996 IPCC estimates
    !
    !       iyr = iyear - 1860 + 1
    !       co2conc = (303.514 - 0.57881 * iyr +
    !    >                      0.00622 * iyr**2 +
    !    >                      1.3e-5 * iyr**3) * 1.e-6
    !
    !
    !     end if
    !
    RETURN
  END SUBROUTINE co2

  ! ---------------------------------------------------------------------
  SUBROUTINE drystress(froot     , &! INTENT(IN   )
       wsoi      , &! INTENT(IN   )
       wisoi     , &! INTENT(IN   )
       swilt     , &! INTENT(IN   )
       sfield    , &! INTENT(IN   ) 
       stressl   , &! INTENT(OUt  ) 
       stressu   , &! INTENT(OUt  ) 
       stresstl  , &! INTENT(OUt  ) 
       stresstu  , &! INTENT(OUt  ) 
       vegtype0  , &! INTENT(IN   ) 
       stressfac , &! INTENT(IN   ) 
       nVegClass , &! INTENT(IN   ) 
       npoi      , &! INTENT(IN   )
       nsoilay     )! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    ! common blocks
    !
    IMPLICIT NONE
    !
    INTEGER , INTENT(IN   ) ::  nVegClass
    INTEGER , INTENT(IN   ) ::  npoi               ! total number of land points
    INTEGER , INTENT(IN   ) ::  nsoilay               ! number of soil layers
    REAL(KIND=r8) , INTENT(IN   ) :: wsoi    (npoi,nsoilay)! fraction of soil pore space containing liquid water
    REAL(KIND=r8) , INTENT(IN   ) :: wisoi   (npoi,nsoilay)! fraction of soil pore space containing ice
    REAL(KIND=r8) , INTENT(IN   ) :: swilt   (npoi,nsoilay)! wilting soil moisture value (fraction of pore space)
    REAL(KIND=r8) , INTENT(IN   ) :: sfield  (npoi,nsoilay)! field capacity soil moisture value (fraction of pore space)
    REAL(KIND=r8) , INTENT(OUt  ) :: stressl (npoi,nsoilay)! soil moisture stress factor for the lower canopy (dimensionless)
    REAL(KIND=r8) , INTENT(OUt  ) :: stressu (npoi,nsoilay)! soil moisture stress factor for the upper canopy (dimensionless)
    REAL(KIND=r8) , INTENT(OUt  ) :: stresstl(npoi)        ! sum of stressl over all 6 soil layers (dimensionless)
    REAL(KIND=r8) , INTENT(OUt  ) :: stresstu(npoi)        ! sum of stressu over all 6 soil layers (dimensionless)
    REAL(KIND=r8) , INTENT(IN   ) :: froot(npoi,nsoilay,2) ! fraction of root in soil layer 
    REAL(KIND=r8) , INTENT(IN   ) :: vegtype0(npoi) 
    REAL(KIND=r8) , INTENT(IN   ) :: stressfac(nVegClass)! to calculate moisture stress factor 
    !
    ! local variables
    !
    INTEGER i    ! loop indices
    INTEGER k    ! loop indices
    INTEGER iveg    ! loop indices
    !
    !REAL(KIND=r8) stressfac ! to calculate moisture stress factor 
    REAL(KIND=r8) awc! available water content (fraction)
    REAL(KIND=r8) znorm! normalizing factor
    REAL(KIND=r8) zwilt! function of awc, =1 if awc = 1 (no stress)
    !
    ! stressfac determines the 'strength' of the soil moisture
    ! stress on physiological processes
    !
    ! strictly speaking, stresst* is multiplied to the vmax
    ! parameters used in the photosynthesis calculations
    !
    ! stressfac determines the shape of the soil moisture response
    !
    !stressfac = -5.0_r8
    !
    !znorm = 1.0_r8 - EXP(stressfac)
    !
    DO i = 1, npoi

       iveg=int(vegtype0(i))

       znorm = 1.0_r8 - exp(stressfac(iveg))

       !
       ! initialize stress parameter
       !
       stresstl(i) = 0.0_r8
       stresstu(i) = 0.0_r8
       !
       ! fraction of soil water uptake in each layer
       !
       DO k = 1, nsoilay
          !
          ! plant available water content (fraction)
          !
          !wsoi    (npoi,nsoilay)    ! fraction of soil pore space containing liquid water
          !swilt   (npoi,nsoilay)    ! wilting soil moisture value        (fraction of pore space)
          !sfield  (npoi,nsoilay)    ! field capacity soil moisture value (fraction of pore space)

          IF((sfield(i,k) - swilt(i,k)) == 0.0_r8 ) THEN
             awc = 1

          ELSE

             awc = MIN (1.0_r8, MAX (0.0_r8,(wsoi(i,k)*(1 - wisoi(i,k)) - swilt(i,k)) / &
                                            (sfield(i,k) - swilt(i,k))))
          END IF
          !
          !                  1 - exp [ stressfac * awc]
          !         zwilt = ----------------------------
          !                    1 - exp [ stressfac ]
          !
          zwilt = (1.0_r8 - EXP(stressfac(iveg) * awc)) / znorm
          !
          ! update for each layer
          !
          stressl(i,k) = froot(i,k,1) * MAX (0.0_r8, MIN (1.0_r8, zwilt))
          stressu(i,k) = froot(i,k,2) * MAX (0.0_r8, MIN (1.0_r8, zwilt))

          !PRINT*,'pkubota',awc,sfield(i,k) , swilt(i,k),zwilt,stressl(i,k)
          !
          ! integral over rooting profile
          !
          stresstl(i) = stresstl(i) + stressl(i,k)
          stresstu(i) = stresstu(i) + stressu(i,k)
          !
       END DO
       !
    END DO
    !
    ! return to main program
    !              -
    !             |     Fc - Fi        Theta_sat - Tehta_ice
    ! wi=MAX(0,MIN| 1, ----------- * -------------------------)
    !             |_     Fc - Fo         Theta_sat  
    !
    !              Fc   
    RETURN
  END SUBROUTINE drystress
      !FHYsat = suction  ! saturated matric potential (m-h2o)
      !#poros   (npoi,nsoilay) ! porosity (mass of h2o per unit vol at sat / rhow)
      !
      !   Fi  is the soil waater matric potential (mm)
      !   Fc and Fo are the soil water potential (mm)
      !   Theta_sat and Tehta_ice  sre the saturated volumetric water e ice content (m3/m3)
! ---------------------------------------------------------------------
  SUBROUTINE canini (bps   , &! INTENT(OUT  )
       rhoa   , &! INTENT(OUT  )
       cp     , &! INTENT(OUT  )
       za     , &! INTENT(OUT  )
       bdl    , &! INTENT(OUT  )
       dil    , &! INTENT(OUT  )
       z3     , &! INTENT(OUT  )
       z4     , &! INTENT(OUT  )
       z34    , &! INTENT(OUT  )
       exphl  , &! INTENT(OUT  )
       expl   , &! INTENT(OUT  )
       displ  , &! INTENT(OUT  )
       bdu    , &! INTENT(OUT  )
       diu    , &! INTENT(OUT  )
       z1     , &! INTENT(OUT  )
       z2     , &! INTENT(OUT  )
       z12    , &! INTENT(OUT  )
       exphu  , &! INTENT(OUT  )
       expu   , &! INTENT(OUT  )
       dispu  , &! INTENT(OUT  )
       alogg  , &! INTENT(OUT  )
       alogi  , &! INTENT(OUT  )
       alogav , &! INTENT(OUT  )
       alog4  , &! INTENT(OUT  )
       alog3  , &! INTENT(OUT  )
       alog2  , &! INTENT(OUT  )
       alog1  , &! INTENT(OUT  )
       aloga  , &! INTENT(OUT  )
       u2     , &! INTENT(OUT  )
       alogu  , &! INTENT(OUT  )
       alogl  , &! INTENT(OUT  )
       ztop   , &! INTENT(IN   )
       fl     , &! INTENT(IN   )
       lai    , &! INTENT(IN   )
       sai    , &! INTENT(IN   )
       alaiml , &! INTENT(IN   )
       zbot   , &! INTENT(IN   )
       fu     , &! INTENT(IN   )
       alaimu , &! INTENT(IN   )
       z0soi  , &! INTENT(IN   )
       fi     , &! INTENT(IN   )
       z0sno  , &! INTENT(IN   )
       psurf  , &! INTENT(IN   )
       ta     , &! INTENT(IN   )
       qa     , &! INTENT(IN   )
       ua     , &! INTENT(IN   )
       npoi   , &! INTENT(IN   ) 
       cappa  , &! INTENT(IN   ) 
       rair   , &! INTENT(IN   ) 
       rvap   , &! INTENT(IN   ) 
       cair   , &! INTENT(IN   ) 
       cvap   , &! INTENT(IN   ) 
       grav     )! INTENT(IN   ) 
    ! ---------------------------------------------------------------------
    !
    ! initializes aerodynamic quantities that remain constant 
    ! through one timestep
    !
    ! note that some quantities actually are
    ! constant as long as the vegetation amounts and fractional
    ! coverage remain unchanged, so could re-arrange code for
    ! efficiency - currently all arrays initialized here are in
    ! com1d which can be overwritten elsewhere
    !
    ! rwork is used throughout as a scratch variable to reduce number of
    ! computations
    !
    !      include 'implicit.h'
    !
    IMPLICIT NONE 
    INTEGER , INTENT(IN   ) :: npoi            ! total number of land points
    REAL(KIND=r8) , INTENT(IN   ) :: cappa           ! rair/cair
    REAL(KIND=r8) , INTENT(IN   ) :: rair            ! gas constant for dry air (J deg-1 kg-1)
    REAL(KIND=r8) , INTENT(IN   ) :: rvap            ! gas constant for water vapor (J deg-1 kg-1)
    REAL(KIND=r8) , INTENT(IN   ) :: cair            ! specific heat of dry air at constant pressure (J deg-1 kg-1)
    REAL(KIND=r8) , INTENT(IN   ) :: cvap            ! specific heat of water vapor at constant pressure (J deg-1 kg-1)
    REAL(KIND=r8) , INTENT(IN   ) :: grav            ! gravitational acceleration (m s-2)
    REAL(KIND=r8) , INTENT(INOUT) :: bps(npoi) ! (ps/p) ** (rair/cair) for atmospheric level  (const)
    REAL(KIND=r8) , INTENT(IN   ) :: psurf (npoi)    ! surface pressure (Pa)
    REAL(KIND=r8) , INTENT(IN   ) :: ta    (npoi)    ! air temperature (K)
    REAL(KIND=r8) , INTENT(IN   ) ::  qa    (npoi)    ! specific humidity (kg_h2o/kg_air)
    REAL(KIND=r8) , INTENT(IN   ) ::  ua    (npoi)    ! wind speed (m s-1)
    REAL(KIND=r8) , INTENT(IN   ) :: fi    (npoi)    ! fractional snow cover
    REAL(KIND=r8) , INTENT(IN   ) :: z0sno               ! roughness length of snow surface (m)
    REAL(KIND=r8) , INTENT(IN   ) :: z0soi (npoi)    ! roughness length of soil surface (m)
    REAL(KIND=r8) , INTENT(IN   ) :: ztop  (npoi,2)  ! height of plant top above ground (m)
    REAL(KIND=r8) , INTENT(IN   ) :: fl    (npoi)    ! fraction of snow-free area covered by lower  canopy
    REAL(KIND=r8) , INTENT(IN   ) :: lai   (npoi,2)  ! canopy single-sided leaf area index (area leaf/area veg)
    REAL(KIND=r8) , INTENT(IN   ) :: sai   (npoi,2)  ! current single-sided stem area index
    REAL(KIND=r8) , INTENT(IN   ) :: alaiml           ! lower canopy leaf & stem maximum area (2 sided) 
    ! for normalization of drag coefficient (m2 m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: zbot  (npoi,2)  ! height of lowest branches above ground (m)
    REAL(KIND=r8) , INTENT(IN   ) :: fu    (npoi)    ! fraction of overall area covered by upper canopy
    REAL(KIND=r8) , INTENT(IN   ) :: alaimu          ! upper canopy leaf & stem area (2 sided) for 
    ! normalization of drag coefficient (m2 m-2)
    REAL(KIND=r8) , INTENT(OUT  ) :: rhoa  (npoi)    ! air density at za (allowing for h2o vapor) (kg m-3)
    REAL(KIND=r8) , INTENT(OUT  ) :: cp    (npoi)    ! specific heat of air at za (allowing for h2o vapor) (J kg-1 K-1)
    REAL(KIND=r8) , INTENT(INOUT) :: za    (npoi)    ! height above the surface of atmospheric forcing (m)
    REAL(KIND=r8) , INTENT(OUT  ) :: bdl   (npoi)    ! aerodynamic coefficient ([(tau/rho)/u**2] for
    ! laower canopy (A31/A30 Pollard & Thompson 1995)
    REAL(KIND=r8) , INTENT(OUT  ) :: dil   (npoi)    ! inverse of momentum diffusion coefficient within lower canopy (m)
    REAL(KIND=r8) , INTENT(OUT  ) :: z3    (npoi)    ! effective top of the lower canopy (for momentum) (m)
    REAL(KIND=r8) , INTENT(OUT  ) :: z4    (npoi)    ! effective bottom of the lower canopy (for momentum) (m)
    REAL(KIND=r8) , INTENT(OUT  ) :: z34   (npoi)    ! effective middle of the lower canopy (for momentum) (m)
    REAL(KIND=r8) , INTENT(OUT  ) :: exphl (npoi)    ! exp(lamda/2*(z3-z4)) for lower canopy (A30 Pollard & Thompson)
    REAL(KIND=r8) , INTENT(OUT  ) :: expl  (npoi)    ! exphl**2
    REAL(KIND=r8) , INTENT(OUT  ) :: displ (npoi)    ! zero-plane displacement height for lower canopy (m)
    REAL(KIND=r8) , INTENT(OUT  ) :: bdu   (npoi)    ! aerodynamic coefficient ([(tau/rho)/u**2] for upper
    ! canopy (A31/A30 Pollard & Thompson 1995)
    REAL(KIND=r8) , INTENT(OUT  ) :: diu   (npoi)    ! inverse of momentum diffusion coefficient within upper canopy (m)
    REAL(KIND=r8) , INTENT(OUT  ) :: z1    (npoi)    ! effective top of upper canopy (for momentum) (m)
    REAL(KIND=r8) , INTENT(OUT  ) :: z2    (npoi)    ! effective bottom of the upper canopy (for momentum) (m)
    REAL(KIND=r8) , INTENT(OUT  ) :: z12   (npoi)    ! effective middle of the upper canopy (for momentum) (m)
    REAL(KIND=r8) , INTENT(OUT  ) :: exphu (npoi)    ! exp(lamda/2*(z3-z4)) for upper canopy (A30 Pollard & Thompson)
    REAL(KIND=r8) , INTENT(OUT  ) :: expu  (npoi)    ! exphu**2
    REAL(KIND=r8) , INTENT(OUT  ) :: dispu (npoi)    ! zero-plane displacement height for upper canopy (m)
    REAL(KIND=r8) , INTENT(OUT  ) :: alogg (npoi)    ! log of soil roughness
    REAL(KIND=r8) , INTENT(OUT  ) :: alogi (npoi)    ! log of snow roughness
    REAL(KIND=r8) , INTENT(OUT  ) :: alogav(npoi)    ! average of alogi and alogg 
    REAL(KIND=r8) , INTENT(OUT  ) :: alog4 (npoi)    ! log (max(z4, 1.1*z0sno, 1.1*z0soi)) 
    REAL(KIND=r8) , INTENT(OUT  ) :: alog3 (npoi)    ! log (z3 - displ)
    REAL(KIND=r8) , INTENT(OUT  ) :: alog2 (npoi)    ! log (z2 - displ)
    REAL(KIND=r8) , INTENT(OUT  ) :: alog1 (npoi)    ! log (z1 - dispu) 
    REAL(KIND=r8) , INTENT(OUT  ) :: aloga (npoi)    ! log (za - dispu) 
    REAL(KIND=r8) , INTENT(OUT  ) :: u2    (npoi)    ! wind speed at level z2 (m s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: alogu (npoi)    ! log (roughness length of upper canopy)
    REAL(KIND=r8) , INTENT(OUT  ) :: alogl (npoi)    ! log (roughness length of lower canopy)
    !
    ! Local variables
    !
    REAL(KIND=r8) :: siga      ! sigma level of atmospheric data
    REAL(KIND=r8) ::   pa        ! pressure at level of atmospheric data
    REAL(KIND=r8) ::   x         ! density of vegetation (without distinction between
    ! lai,sai)
    REAL(KIND=r8) ::   x1        ! density of vegetation (different max)
    REAL(KIND=r8) ::   rwork     ! difference between top and bottom of canopy 
    REAL(KIND=r8) ::   cvegl     !
    REAL(KIND=r8) ::   dvegl     ! diffusion coefficient for lower canopy
    REAL(KIND=r8) ::   bvegl     ! e-folding depth in canopy for lower canopy
    REAL(KIND=r8) ::   cvegu     !
    REAL(KIND=r8) ::  dvegu     ! diffusion coefficient for upper canopy
    REAL(KIND=r8) ::  bvegu     ! e-folding depth in canopy for upper canopy

    INTEGER  :: i          ! loop indice

    !
    ! define sigma level of atmospheric data
    !
    ! currently, the value of siga is set to 0.999. This is roughly 10 meters
    ! above ground, which is the typical height for the CRU05 input wind speed data
    !
    siga = 0.999_r8
    !
    !tfac = 1.0_r8 / (siga**cappa)
    !
    ! atmospheric conditions at za
    ! za is variable, although siga = p/ps is constant
    !
    DO  i = 1, npoi
       ! PK bps(i) = 1.0_r8/(siga**cappa)!pkubota
       !
       pa = psurf(i) * siga
       !
       rhoa(i) = pa / ( rair * ta(i) *  &
            (1.0_r8 + (rvap / rair - 1.0_r8) * qa(i)) )
       !
       cp(i) = cair * (1.0_r8 + (cvap / cair - 1.0_r8) * qa(i))
       !
       !iPK za(i) = (psurf(i) - pa) / (rhoa(i) * grav)
       !
       ! make sure that atmospheric level is higher than canopy top
       !
       za(i) = MAX (za(i), ztop(i,2) + 1.0_r8)
       !
    END DO
    !
    ! aerodynamic coefficients for the lower story
    !
    ! cvegl (drag coeff for momentum) is proportional, and dvegl
    ! (diffusion coeff for momentum) inversely proportional,
    ! to x = density of vegetation (without distinction between
    ! lai,sai and fl*(1-fi)) - x is not allowed to be exactly
    ! zero to avoid divide-by-zeros, and for x>1 dvegl is 
    ! proportional to 1/x**2 so that roughness length tends to
    ! zero as x tends to infinity
    !
    ! also the top, bottom and displacement heights z3(i),z4(i),
    ! displ(i) tend to particular values as the density tends to
    ! zero, to give same results as equations for no veg at all.
    !
    DO i = 1, npoi
       !
       x = fl(i) * (1.0_r8 - fi(i)) * 2.0_r8 * (lai(i,1)  &
            + sai(i,1)) / alaiml
       !
       x  = MIN (x, 3.0_r8)
       x1 = MIN (x, 1.0_r8)
       !
       rwork = MAX(ztop(i,1)-zbot(i,1),0.01_r8)
       cvegl = (0.4_r8 / rwork) * MAX(1.e-5_r8, x)
       !
       dvegl = (0.1_r8 * rwork) / MAX(1.e-5_r8, x, x**2)
       !
       ! e-folding depth in canopy
       !
       bvegl = SQRT (2.0_r8 * cvegl / dvegl )
       !
       ! [(tau/rho)/u**2] for inf canopy
       !
       bdl(i) = 0.5_r8 * bvegl * dvegl
       !
       ! 1 / diffusion coefficient
       !
       dil(i) = 1.0_r8 / dvegl
       !
       rwork = (1.0_r8 - x1) * (MAX (z0soi(i),z0sno) + 0.01_r8) 
       !
       z3(i) = x1 * ztop(i,1) + rwork
       !
       z4(i) = x1 * zbot(i,1) + rwork
       !
       z34(i) = 0.5_r8 * (z3(i) + z4(i))
       !
       exphl(i) = EXP (0.5_r8 * bvegl * (z3(i)-z4(i)))
       expl(i)  = exphl(i)**2
       !
       displ(i) = x1 * 0.7_r8 * z3(i)
       !
    END DO
    !
    ! aerodynamic coefficients for the upper story
    ! same comments as for lower story
    !
    DO i = 1, npoi
       !
       x = fu(i) * 2.0_r8 * (lai(i,2)+sai(i,2)) / alaimu
       !
       x  = MIN (x, 3.0_r8)
       x1 = MIN (x, 1.0_r8)
       !
       rwork = MAX(ztop(i,2)-zbot(i,2),.01_r8)
       cvegu = (0.4_r8 / rwork) *  &
            MAX(1.e-5_r8,x)
       !
       dvegu = (0.1_r8 * rwork) /  &
            MAX(1.e-5_r8,x,x**2)
       !
       rwork = 1.0_r8 / dvegu
       bvegu  = SQRT (2.0_r8 * cvegu * rwork)
       bdu(i) = 0.5_r8 * bvegu * dvegu
       diu(i) = rwork
       !
       rwork = (1.0_r8 - x1) * (z3(i) + 0.01_r8)
       z1(i) = x1 * ztop(i,2) + rwork
       z2(i) = x1 * zbot(i,2) + rwork
       !
       z12(i) = 0.5_r8 * (z1(i) + z2(i))
       !
       exphu(i) = EXP (0.5_r8 * bvegu * (z1(i) - z2(i)))
       expu(i)  = exphu(i)**2
       !
       dispu(i) = x1 * 0.7_r8 * z1(i) + (1.0_r8 - x1) * displ(i)
       !
    END DO
    !
    ! mixing-length logarithms
    !
    DO i = 1, npoi
       !
       alogg(i)  = LOG (z0soi(i))
       alogi(i)  = LOG (z0sno)
       alogav(i) = (1.0_r8 - fi(i)) * alogg(i) + fi(i) * alogi(i)
       !
       ! alog4 must be > z0soi, z0sno to avoid possible problems later 
       !
       alog4(i) = LOG ( MAX (z4(i), 1.1_r8*z0soi(i), 1.1_r8*z0sno) )
       alog3(i) = LOG (z3(i)-displ(i))
       alog2(i) = LOG (z2(i)-displ(i))
       alog1(i) = LOG (z1(i)-dispu(i))
       aloga(i) = LOG (za(i)-dispu(i))
       !
       ! initialize u2, alogu, alogl for first iteration's fstrat
       !
       u2(i)    = ua(i)/exphu(i)
       alogu(i) = LOG (MAX(.01_r8, .1_r8*(z1(i)-z2(i))))
       alogl(i) = LOG (MAX(.01_r8, .1_r8*(z3(i)-z4(i))))
       !
    END DO
    !
    RETURN
  END SUBROUTINE canini
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE turcof (iter        , &! INTENT(IN   )
       z3          , &! INTENT(IN   )
       z2          , &! INTENT(IN   )
       alogl       , &! INTENT(INOUT)
       u2          , &! INTENT(INOUT)
       richl       , &! INTENT(OUT  )
       straml      , &! INTENT(OUT  )
       strahl      , &! INTENT(OUT  )
       bps         , &! INTENT(IN   )
       z1          , &! INTENT(IN   )
       za          , &! INTENT(IN   )
       alogu       , &! INTENT(INOUT)
       aloga       , &! INTENT(IN   )
       richu       , &! INTENT(OUT  )
       stramu      , &! INTENT(OUT  )
       strahu      , &! INTENT(OUT  )
       alog4       , &! INTENT(IN   )
       alogav      , &! INTENT(IN   )
       bdl         , &! INTENT(IN   )
       expl        , &! INTENT(IN   )
       alog3       , &! INTENT(IN   )
       bdu         , &! INTENT(IN   )
       expu        , &! INTENT(IN   )
       alog1       , &! INTENT(IN   )
       u1          , &! INTENT(OUT  )
       u12         , &! INTENT(OUT  )
       exphu       , &! INTENT(IN   )
       u3          , &! INTENT(OUT  )
       u34         , &! INTENT(OUT  )
       exphl       , &! INTENT(IN   )
       u4          , &! INTENT(OUT  )
       rhoa        , &! INTENT(IN   )
       diu         , &! INTENT(IN   )
       z12         , &! INTENT(IN   )
       dil         , &! INTENT(IN   )
       z34         , &! INTENT(IN   )
       z4          , &! INTENT(IN   )
       cu          , &! INTENT(OUT  )
       cl          , &! INTENT(OUT  )
       sg          , &! INTENT(OUT  )
       si          , &! INTENT(OUT  )
       alog2       , &! INTENT(OUT  )
       t34         , &! INTENT(IN   )
       t12         , &! INTENT(IN   )
       q34         , &! INTENT(IN   )
       q12         , &! INTENT(IN   )
       su          , &! INTENT(OUT  )
       cleaf       , &! INTENT(IN   )
       dleaf       , &! INTENT(IN   )
       ss          , &! INTENT(OUT  )
       cstem       , &! INTENT(IN   )
       dstem       , &! INTENT(IN   )
       sl          , &! INTENT(OUT  )
       cgrass      , &! INTENT(IN   )
       ua          , &! INTENT(IN   )
       ta          , &! INTENT(IN   )
       qa          , &! INTENT(IN   )
       npoi        , &! INTENT(IN   )
       dtime       , &! INTENT(IN   )
       vonk        , &! INTENT(IN   )
       grav          )! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    ! solves for wind speeds at various levels
    !
    ! also computes upper and lower-region air-air transfer coefficients
    ! and saves them in com1d arrays cu and cl for use by turvap,
    ! and similarly for the solid-air transfer coefficients
    ! su, ss, sl, sg and si
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN   ) :: npoi              ! total number of land points
    REAL(KIND=r8), INTENT(IN   ) :: dtime                ! model timestep (seconds)
    REAL(KIND=r8), INTENT(IN   ) :: vonk                ! von karman constant (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: grav              ! gravitational acceleration (m s-2)
    REAL(KIND=r8), INTENT(IN   ) :: ua(npoi)  ! wind speed (m s-1)
    REAL(KIND=r8), INTENT(IN   ) :: ta(npoi)  ! air temperature (K)
    REAL(KIND=r8), INTENT(IN   ) :: qa(npoi)  ! specific humidity (kg_h2o/kg_air)
    REAL(KIND=r8), INTENT(IN   ) :: t34   (npoi)      ! air temperature at z34 (K)
    REAL(KIND=r8), INTENT(IN   ) :: t12   (npoi)      ! air temperature at z12 (K)
    REAL(KIND=r8), INTENT(IN   ) :: q34   (npoi)      ! specific humidity of air at z34
    REAL(KIND=r8), INTENT(IN   ) :: q12   (npoi)      ! specific humidity of air at z12
    REAL(KIND=r8), INTENT(OUT  ) :: su    (npoi)      ! air-vegetation transfer coefficients (*rhoa) 
    !for upper canopy leaves (m s-1 * kg m-3) (A39a Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(IN   ) :: cleaf             ! empirical constant in upper canopy leaf-air aerodynamic
    ! transfer coefficient (m s-0.5) (A39a Pollard & Thompson 95)
    REAL(KIND=r8), INTENT(IN   ) :: dleaf (2)         ! typical linear leaf dimension in aerodynamic transfer coefficient (m)
    REAL(KIND=r8), INTENT(OUT  ) :: ss    (npoi)      ! air-vegetation transfer coefficients (*rhoa) for upper 
    ! canopy stems (m s-1 * kg m-3) (A39a Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(IN   ) :: cstem             ! empirical constant in upper canopy stem-air aerodynamic transfer
    ! coefficient (m s-0.5) (A39a Pollard & Thompson 95)
    REAL(KIND=r8), INTENT(IN   ) :: dstem (2)         ! typical linear stem dimension in aerodynamic transfer coefficient (m)
    REAL(KIND=r8), INTENT(OUT  ) :: sl    (npoi)      ! air-vegetation transfer coefficients (*rhoa) for lower canopy
    ! leaves & stems (m s-1*kg m-3) (A39a Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(IN   ) :: cgrass            ! empirical constant in lower canopy-air aerodynamic transfer
    ! coefficient (m s-0.5) (A39a Pollard & Thompson 95)
    REAL(KIND=r8), INTENT(IN   ) :: z3    (npoi)      ! effective top of the lower canopy (for momentum) (m)
    REAL(KIND=r8), INTENT(IN   ) :: z2    (npoi)      ! effective bottom of the upper canopy (for momentum) (m)
    REAL(KIND=r8), INTENT(INOUt) :: alogl (npoi)      ! log (roughness length of lower canopy)
    REAL(KIND=r8), INTENT(INOUt) :: u2    (npoi)      ! wind speed at level z2 (m s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: richl (npoi)      ! richardson number for air above upper canopy (z3 to z2)
    REAL(KIND=r8), INTENT(OUT  ) :: straml(npoi)      ! momentum correction factor for stratif between upper &
    ! lower canopy (z3 to z2) (louis et al.)
    REAL(KIND=r8), INTENT(OUT  ) :: strahl(npoi)      ! heat/vap correction factor for stratif between upper &
    ! lower canopy (z3 to z2) (louis et al.)
    REAL(KIND=r8), INTENT(INOUT) :: bps   (npoi)                ! (ps/p) ** (rair/cair) for atmospheric level  (const)
    REAL(KIND=r8), INTENT(IN   ) :: z1    (npoi)      ! effective top of upper canopy (for momentum) (m)
    REAL(KIND=r8), INTENT(IN   ) :: za    (npoi)      ! height above the surface of atmospheric forcing (m)
    REAL(KIND=r8), INTENT(INOUT) :: alogu (npoi)      ! log (roughness length of upper canopy)
    REAL(KIND=r8), INTENT(IN   ) :: aloga (npoi)      ! log (za - dispu) 
    REAL(KIND=r8), INTENT(OUT  ) :: richu (npoi)      ! richardson number for air between upper & lower canopy (z1 to za)
    REAL(KIND=r8), INTENT(OUT  ) :: stramu(npoi)      ! momentum correction factor for stratif above upper canopy
    ! (z1 to za) (louis et al.)
    REAL(KIND=r8), INTENT(OUT  ) :: strahu(npoi)      ! heat/vap correction factor for stratif above upper canopy 
    ! (z1 to za) (louis et al.)
    REAL(KIND=r8), INTENT(IN   ) :: alog4 (npoi)      ! log (max(z4, 1.1*z0sno, 1.1*z0soi)) 
    REAL(KIND=r8), INTENT(IN   ) :: alogav(npoi)      ! average of alogi and alogg 
    REAL(KIND=r8), INTENT(IN   ) :: bdl   (npoi)      ! aerodynamic coefficient ([(tau/rho)/u**2] for laower canopy
    ! (A31/A30 Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(IN   ) :: expl  (npoi)      ! exphl**2
    REAL(KIND=r8), INTENT(IN   ) :: alog3 (npoi)      ! log (z3 - displ)
    REAL(KIND=r8), INTENT(IN   ) :: bdu   (npoi)      ! aerodynamic coefficient ([(tau/rho)/u**2] for upper canopy
    ! (A31/A30 Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(IN   ) :: expu  (npoi)      ! exphu**2
    REAL(KIND=r8), INTENT(IN   ) :: alog1 (npoi)      ! log (z1 - dispu) 
    REAL(KIND=r8), INTENT(OUT  ) :: u1    (npoi)      ! wind speed at level z1 (m s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: u12   (npoi)      ! wind speed at level z12 (m s-1)
    REAL(KIND=r8), INTENT(IN   ) :: exphu (npoi)      ! exp(lamda/2*(z3-z4)) for upper canopy (A30 Pollard & Thompson)
    REAL(KIND=r8), INTENT(OUT  ) :: u3    (npoi)      ! wind speed at level z3 (m s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: u34   (npoi)      ! wind speed at level z34 (m s-1)
    REAL(KIND=r8), INTENT(IN   ) :: exphl (npoi)      ! exp(lamda/2*(z3-z4)) for lower canopy (A30 Pollard & Thompson)
    REAL(KIND=r8), INTENT(OUT  ) :: u4    (npoi)      ! wind speed at level z4 (m s-1)
    REAL(KIND=r8), INTENT(IN   ) :: rhoa  (npoi)      ! air density at za (allowing for h2o vapor) (kg m-3)
    REAL(KIND=r8), INTENT(IN   ) :: diu   (npoi)      ! inverse of momentum diffusion coefficient within upper canopy (m)
    REAL(KIND=r8), INTENT(IN   ) :: z12   (npoi)      ! effective middle of the upper canopy (for momentum) (m)
    REAL(KIND=r8), INTENT(IN   ) :: dil   (npoi)      ! inverse of momentum diffusion coefficient within lower canopy (m)
    REAL(KIND=r8), INTENT(IN   ) :: z34   (npoi)      ! effective middle of the lower canopy (for momentum) (m)
    REAL(KIND=r8), INTENT(IN   ) :: z4    (npoi)      ! effective bottom of the lower canopy (for momentum) (m)
    REAL(KIND=r8), INTENT(OUT  ) :: cu    (npoi)      ! air transfer coefficient (*rhoa) (m s-1 kg m-3) for upper
    ! air region (z12 --> za) (A35 Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(OUT  ) :: cl    (npoi)      ! air transfer coefficient (*rhoa) (m s-1 kg m-3) between
    ! the 2 canopies (z34 --> z12) (A36 Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(OUT  ) :: sg    (npoi)      ! air-soil transfer coefficient
    REAL(KIND=r8), INTENT(OUT  ) :: si    (npoi)      ! air-snow transfer coefficient
    REAL(KIND=r8), INTENT(IN   ) :: alog2 (npoi)      ! log (z2 - displ)
    !
    ! Arguments (input)
    !
    INTEGER, INTENT(IN   ) :: iter          !current iteration number
    !
    !
    ! Local variables
    !
    INTEGER i                ! loop indice
    !
    REAL(KIND=r8) xfac(npoi)             !
    REAL(KIND=r8) x                !
    REAL(KIND=r8) rwork            ! working variable
    REAL(KIND=r8) cdmax            ! max value for cd
    REAL(KIND=r8) tauu             !
    REAL(KIND=r8) a                !
    REAL(KIND=r8) b                !
    REAL(KIND=r8) c                !
    REAL(KIND=r8) d                !
    REAL(KIND=r8) taul             !
    REAL(KIND=r8) ca               ! to compute inverse air-air transfer coeffs
    REAL(KIND=r8) cai 
    REAL(KIND=r8) cbi 
    REAL(KIND=r8) cci              !
    REAL(KIND=r8) cdi 
    REAL(KIND=r8) cei 
    REAL(KIND=r8) cfi              !
    REAL(KIND=r8) sg0               ! to compute air-solid transfer coeff for soil
    REAL(KIND=r8) si0              ! to compute air-solid transfer coeff for ice
    !
    REAL(KIND=r8) yu(npoi) 
    REAL(KIND=r8) yl(npoi)
    !
    ! set stratification factors for lower and upper regions
    ! using values from the previous iteration
    !
    xfac = 1.0_r8
    !
    CALL fstrat (t34        , & ! INTENT(IN)
         t12        , & ! INTENT(IN)
         xfac       , & ! INTENT(IN)
         q34        , & ! INTENT(IN)
         q12        , & ! INTENT(IN)
         z3         , & ! INTENT(IN)
         z2         , & ! INTENT(IN)
         alogl      , & ! INTENT(IN)
         alogl      , & ! INTENT(IN)
         alog2      , & ! INTENT(IN)
         u2         , & ! INTENT(IN)
         richl      , & ! INTENT(OUT  )
         straml     , & ! INTENT(OUT  )
         strahl     , & ! INTENT(OUT  )
         iter       , & ! INTENT(IN)
         npoi       , & ! INTENT(IN)
         grav       , & ! INTENT(IN)
         vonk         ) ! INTENT(IN)
    !
    CALL fstrat (t12        , &! INTENT(IN   )
         ta         , &! INTENT(IN   )
         bps        , &! INTENT(IN   )
         q12        , &! INTENT(IN   )
         qa         , &! INTENT(IN   )
         z1         , &! INTENT(IN   )
         za         , &! INTENT(IN   )
         alogu      , &! INTENT(IN   )
         alogu      , &! INTENT(IN   )
         aloga      , &! INTENT(IN   )
         ua         , &! INTENT(IN   )
         richu      , &! INTENT(OUT  )
         stramu     , &! INTENT(OUT  )
         strahu     , &! INTENT(OUT  )
         iter       , &! INTENT(IN   )
         npoi       , &! INTENT(IN   )
         grav       , &! INTENT(IN   )
         vonk         )! INTENT(IN   )
    !
    ! eliminate c/d from eq (28), tau_l/rho from (26),(27), to get
    ! lower-story roughness alogl. yl/bdl is (tau_l/rho)/(c+d)
    !
    ! equation numbers correspond to lsx description section 4.e
    !
    DO i = 1, npoi
       !
       x = ((alog4(i)-alogav(i))/vonk)**2 * bdl(i)
       !
       rwork = 1.0_r8 / expl(i)
       yl(i) = ((x+1)*expl(i) + (x-1)*rwork)   &
            / ((x+1)*expl(i) - (x-1)*rwork)
       !
       alogl(i) = alog3(i) - vonk * SQRT(yl(i)/bdl(i))
       !
    END DO
    !
    ! eliminate tau_l/rho from (24),(25), tau_u/rho and a/b from
    ! (22),(23), to get upper-story roughness alogu
    ! 
    ! yu/bdu is (tau_u/rho)/(a+b)
    !
    DO i = 1, npoi
       !          
       x = ((alog2(i)-alogl(i))/vonk)**2 * bdu(i) / straml(i)
       !
       rwork = 1.0_r8 / expu(i)
       yu(i) = ((x+1)*expu(i) + (x-1)*rwork)    &
            / ((x+1)*expu(i) - (x-1)*rwork)
       !
       alogu(i) = alog1(i) - vonk * SQRT(yu(i)/bdu(i))
       !
    END DO
    !
    ! define the maximum value of cd
    !
    !cdmax = 300.0_r8 / (2.0_r8 * dtime)
    cdmax = 300.0_r8 / (1.0_r8 * dtime) ! 0.125 
    !cdmax = 0.125_r8 
    !
    ! get tauu (=tau_u/rho) from (21), a and b from (22),(23),
    ! taul (=tau_u/rho) from (25), c and d from (26),(27)
    !
    ! changed the following to eliminate small errors associated with
    ! moving this code to single precision - affected c and d,
    ! which made u_ become undefined, as well as affecting some
    ! other variables
    !
    DO i = 1, npoi
       !
       tauu = (ua(i) * vonk/(aloga(i)-alogu(i)))**2 * stramu(i)
       !
       a = 0.5_r8 * tauu * (yu(i)+1)/bdu(i)
       b = 0.5_r8 * tauu * (yu(i)-1)/bdu(i)
       !
       taul = bdu(i) * (a/expu(i) - b*expu(i))
       !
       c = 0.5_r8 * taul * (yl(i)+1)/bdl(i)
       d = 0.5_r8 * taul * (yl(i)-1)/bdl(i)
       !
       ! evaluate wind speeds at various levels, keeping a minimum 
       ! wind speed of 0.01_r8 m/s at all levels
       !   
       u1(i)  = MAX (0.01_r8, SQRT (MAX (0.0_r8, (a+b))))
       u12(i) = MAX (0.01_r8, SQRT (MAX (0.0_r8, (a/exphu(i)+b*exphu(i)))))
       u2(i)  = MAX (0.01_r8, SQRT (MAX (0.0_r8, (a/expu(i) +b*expu(i)))))
       u3(i)  = MAX (0.01_r8, SQRT (MAX (0.0_r8, (c+d))))
       u34(i) = MAX (0.01_r8, SQRT (MAX (0.0_r8, (c/exphl(i)+d*exphl(i)))))
       u4(i)  = MAX (0.01_r8, SQRT (MAX (0.0_r8, (c/expl(i) +d*expl(i)))))
       !
    END DO
    !
    ! compute inverse air-air transfer coeffs
    !
    ! use of inverse individual coeffs cai, cbi, cci, cdi, cei, cfi avoids
    ! divide-by-zero as vegetation vanishes - combine into
    ! upper-region coeff cu from za to z12, and lower-region coeff
    ! cl from z34 to z12, and also coeffs
    !
    DO i = 1, npoi
       !
       ca = ua(i)*strahu(i)*vonk**2  /   &
            ((aloga(i)-alogu(i)) * (aloga(i)-alog1(i)))
       !
       ca = MIN (cdmax, ca / (1.0_r8 + ca * 1.0e-20_r8))
       !
       cai = 1.0_r8 / (rhoa(i)*ca)
       !
       cbi = diu(i) * (z1(i)-z12(i)) / (rhoa(i) * 0.5_r8*(u1(i)+u12(i)))
       cci = diu(i) * (z12(i)-z2(i)) / (rhoa(i) * 0.5_r8*(u12(i)+u2(i)))
       !
       cdi = (alog2(i)-alogl(i)) * (alog2(i)-alog3(i)) /  &
            (rhoa(i)*u2(i)*strahl(i)*vonk**2)
       !
       cei = dil(i) * (z3(i)-z34(i)) / (rhoa(i) * 0.5_r8*(u3(i)+u34(i)))
       cfi = dil(i) * (z34(i)-z4(i)) / (rhoa(i) * 0.5_r8*(u34(i)+u4(i)))
       !
       cu(i) = 1.0_r8 / (cai + cbi)
       cl(i) = 1.0_r8 / (cci + cdi + cei)
       !
       ! compute air-solid transfer coeffs for upper leaves, upper
       ! stems, lower story (su,ss,sl)
       !
       su(i) = rhoa(i) * cleaf  * SQRT (u12(i) / dleaf(2))
       ss(i) = rhoa(i) * cstem  * SQRT (u12(i) / dstem(2))
       sl(i) = rhoa(i) * cgrass * SQRT (u34(i) / dleaf(1))
       !
       ! compute air-solid transfer coeffs for soil and snow (sg,si)
       !
       ! old technique
       !
       !       sg0 = rhoa(i) * u4(i) * (vonk/(alog4(i)-alogg(i)))**2
       !       si0 = rhoa(i) * u4(i) * (vonk/(alog4(i)-alogi(i)))**2
       !
       ! replace above formulations which depend on the log-wind profile
       ! (which may not work well below a canopy), with empirical formulation
       ! of Norman's. In the original LSX, turcof.f solves for the winds at
       ! the various levels from the momentum equations. This gives the transfer
       ! coefficients for heat and moisture. Heat and moisture eqns are then solved 
       ! in subroutine turvap. Using the empirical formulation of John Norman is 
       ! not consistent with the earlier solution for u4 (based on a logarithmic 
       ! profile just above the ground. However, this is used here because it 
       ! improved a lot simulations of the sensible heat flux over the 
       ! HAPEX-MOBILHY and FIFE sites
       !
       sg0 = rhoa(i) * (0.004_r8 + 0.012_r8 * u4(i))
       si0 = rhoa(i) * (0.003_r8 + 0.010_r8 * u4(i))
       !
       ! modify the cofficient to deal with cfi (see above)
       !
       sg(i) = 1.0_r8 / (cfi + 1.0_r8 / sg0)
       si(i) = 1.0_r8 / (cfi + 1.0_r8 / si0)
       !
    END DO
    !
    ! JAF:  not necessary 
    !
    ! if no veg, recalculate coefficients appropriately for a
    ! single logarithmic profile, and 2 fictitious levels just
    ! above soil/snow surface. these levels are arbitrary but are
    ! taken as z2 and z4, preset in vegdat to a few cm height
    ! for bare ground and ice. use strahu from above, which used
    ! t12 and alogu (ok after first iteration)
    !
    !     do 600 i = 1, npoi
    !
    !       if ((fu(i).eq.0.0).and.(fl(i).eq.0.0)) then
    !
    !         z = rhoa(i)*ua(i)*strahu(i)*vonk**2 / (aloga(i)-alogav(i))
    !
    !         ca    = z / (aloga(i)-alog2(i))
    !         cu(i) = rhoa(i)*min (cdmax,
    !    >                          ca / (1. + ca / 1.0e+20))
    !
    !         cl(i) = z / (alog2(i)-alog4(i))
    !
    !         sg(i) = z / (alog4(i)-alogg(i))
    !         si(i) = z / (alog4(i)-alogi(i))
    !
    !         alogu(i) = alogav(i)
    !
    !       endif
    !
    ! 600 continue
    !
    RETURN
  END SUBROUTINE turcof
  !
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE turvap(iter     , &! INTENT(IN   )
       niter    , &! INTENT(IN   )
       cp       , &! INTENT(IN   )
       sg       , &! INTENT(IN   )
       fwetux   , &! INTENT(OUT  )
       fwetu    , &! INTENT(IN   )
       fwetsx   , &! INTENT(OUT  )
       fwets    , &! INTENT(IN   )
       fwetlx   , &! INTENT(OUT  )
       fwetl    , &! INTENT(IN   )
       solu     , &! INTENT(IN   )
       firu     , &! INTENT(IN   )
       sols     , &! INTENT(IN   )
       firs     , &! INTENT(IN   )
       soll     , &! INTENT(IN   )
       firl     , &! INTENT(IN   )
       rliqu    , &! INTENT(IN   )
       rliqs    , &! INTENT(IN   )
       rliql    , &! INTENT(IN   )
       pfluxu   , &! INTENT(IN   )
       pfluxs   , &! INTENT(IN   )
       pfluxl   , &! INTENT(IN   )
       cu       , &! INTENT(IN   )
       cl       , &! INTENT(IN   )
       bps     , &! INTENT(IN   )
       si       , &! INTENT(IN   )
       solg     , &! INTENT(IN   )
       firg     , &! INTENT(INOUT)
       soli     , &! INTENT(IN   )
       firi     , &! INTENT(INOUT)
       fsena    , &! INTENT(OUT  )
       fseng    , &! INTENT(OUT  )
       fseni    , &! INTENT(OUT  )
       fsenu    , &! INTENT(OUT  )
       fsens    , &! INTENT(OUT  )
       fsenl    , &! INTENT(OUT  )
       fvapa    , &! INTENT(OUT  )
       fvapuw   , &! INTENT(OUT  )
       fvaput   , &! INTENT(OUT  )
       fvaps    , &! INTENT(OUT  )
       fvaplw   , &! INTENT(OUT  )
       fvaplt   , &! INTENT(OUT  )
       fvapg    , &! INTENT(OUT  )
       fvapi    , &! INTENT(OUT  )
       firb     , &! INTENT(INOUT)
       lai      , &! INTENT(IN   )
       fu       , &! INTENT(IN   )
       sai      , &! INTENT(IN   )
       fl       , &! INTENT(IN   )
       chu      , &! INTENT(IN   )
       wliqu    , &! INTENT(IN   )
       wsnou    , &! INTENT(IN   )
       chs      , &! INTENT(IN   )
       wliqs    , &! INTENT(IN   )
       wsnos    , &! INTENT(IN   )
       chl      , &! INTENT(IN   )
       wliql    , &! INTENT(IN   )
       wsnol    , &! INTENT(IN   )
       tu       , &! INTENT(INOUT)
       ts       , &! INTENT(INOUT)
       tl       , &! INTENT(INOUT)
       q34      , &! INTENT(INOUT)
       t34      , &! INTENT(INOUT)
       q12      , &! INTENT(INOUT)
       su       , &! INTENT(IN   )
       totcondub, &! INTENT(IN   )
       frac     , &! INTENT(IN   )
       totconduc, &! INTENT(IN   )
       ss       , &! INTENT(IN   )
       sl       , &! INTENT(IN   )
       totcondls, &! INTENT(IN   )
       totcondl4, &! INTENT(IN   )
       totcondl3, &! INTENT(IN   )
       t12      , &! INTENT(INOUT)
       poros    , &! INTENT(IN   )
       wpud     , &! INTENT(IN   )
       wipud    , &! INTENT(IN   )
       csoi     , &! INTENT(IN   )
       rhosoi   , &! INTENT(IN   )
       wisoi    , &! INTENT(IN   )
       wsoi     , &! INTENT(IN   )
       hsoi     , &! INTENT(IN   )
       consoi   , &! INTENT(IN   )
       tg       , &! INTENT(INOUT)
       ti       , &! INTENT(INOUT)
       wpudmax  , &! INTENT(IN   )
       suction  , &! INTENT(IN   )
       bex      , &! INTENT(IN   )
       swilt    , &! INTENT(IN   )
       hvasug   , &! INTENT(IN   )
       tsoi     , &! INTENT(IN   )
       hvasui   , &! INTENT(IN   )
       upsoiu   , &! INTENT(OUT  )
       stressu  , &! INTENT(IN   )
       stresstu , &! INTENT(IN   )
       upsoil   , &! INTENT(OUT  )
       stressl  , &! INTENT(IN   )
       stresstl , &! INTENT(IN   )
       fi       , &! INTENT(IN   )
       consno   , &! INTENT(IN   )
       hsno     , &! INTENT(IN   )
       hsnotop  , &! INTENT(IN   )
       tsno         , &! INTENT(IN   )
       psurf    , &! INTENT(IN   )
       ta       , &! INTENT(IN   )
       qa       , &! INTENT(IN   )
       ginvap   , &! INTENT(OUT  )
       gsuvap   , &! INTENT(OUT  )
       gtrans   , &! INTENT(OUT  )
       gtransu  , &! INTENT(OUT  )
       gtransl  , &! INTENT(OUT  )
       xu       , &! INTENT(INOUT)
       xs       , &! INTENT(INOUT)
       xl       , &! INTENT(INOUT)
       chux     , &! INTENT(INOUT)
       chsx     , &! INTENT(INOUT)
       chlx     , &! INTENT(INOUT)
       chgx     , &! INTENT(INOUT)
       wlgx     , &! INTENT(INOUT)
       wigx     , &! INTENT(INOUT)
       cog      , &! INTENT(INOUT)
       coi      , &! INTENT(INOUT)
       zirg     , &! INTENT(INOUT)
       ziri     , &! INTENT(INOUT)
       wu       , &! INTENT(INOUT)
       ws       , &! INTENT(INOUT)
       wl       , &! INTENT(INOUT)
       wg       , &! INTENT(INOUT)
       wi       , &! INTENT(INOUT)
       tuold    , &! INTENT(INOUT)
       tsold    , &! INTENT(INOUT)
       tlold    , &! INTENT(INOUT)
       tgold    , &! INTENT(INOUT)
       tiold    , &! INTENT(INOUT)
       npoi     , &! INTENT(IN   )
       nsoilay  , &! INTENT(IN   )
       nsnolay  , &! INTENT(IN   )
       npft     , &! INTENT(IN   )
       hvap     , &! INTENT(IN   )
       cvap     , &! INTENT(IN   )
       ch2o     , &! INTENT(IN   )
       hsub     , &! INTENT(IN   )
       cice     , &! INTENT(IN   )
       rhow     , &! INTENT(IN   )
       stef     , &! INTENT(IN   )
       tmelt    , &! INTENT(IN   )
       hfus     , &! INTENT(IN   )
       epsilon  , &! INTENT(IN   )
       grav     , &! INTENT(IN   )
       rvap     , &! INTENT(IN   )
       vegtype0 , &! INTENT(IN   )
       dtime     ,&! INTENT(IN   )
       nVegClass  )! INTENT(IN   )

    ! ---------------------------------------------------------------------
    !
    ! solves canopy system with linearized implicit sensible heat and
    ! moisture fluxes
    !
    ! first, assembles matrix arr of coeffs in linearized equations
    ! for tu,ts,tl,t12,t34,q12,q34,tg,ti and assembles the right hand
    ! sides in the rhs vector
    !
    ! then calls linsolve to solve this system, passing template mplate of
    ! zeros of arr 
    ! 
    ! finally calculates the implied fluxes and stores them 
    ! for the agcm, soil, snow models and budget calcs
    !
    ! common blocks
    !
    IMPLICIT NONE
    !
    INTEGER , INTENT(IN   ) :: nVegClass
    INTEGER , INTENT(IN   ) :: npoi    ! total number of land points
    INTEGER , INTENT(IN   ) :: nsoilay ! number of soil layers
    INTEGER , INTENT(IN   ) :: nsnolay ! number of snow layers
    INTEGER , INTENT(IN   ) :: npft    ! number of plant functional types
    REAL(KIND=r8) , INTENT(IN   ) :: hvap            ! latent heat of vaporization of water (J kg-1)
    REAL(KIND=r8) , INTENT(IN   ) :: cvap            ! specific heat of water vapor at constant pressure (J deg-1 kg-1)
    REAL(KIND=r8) , INTENT(IN   ) :: ch2o            ! specific heat of liquid water (J deg-1 kg-1)
    REAL(KIND=r8) , INTENT(IN   ) :: hsub            ! latent heat of sublimation of ice (J kg-1)
    REAL(KIND=r8) , INTENT(IN   ) :: cice            ! specific heat of ice (J deg-1 kg-1)
    REAL(KIND=r8) , INTENT(IN   ) :: rhow            ! density of liquid water (all types) (kg m-3)
    REAL(KIND=r8) , INTENT(IN   ) :: stef            ! stefan-boltzmann constant (W m-2 K-4)
    REAL(KIND=r8) , INTENT(IN   ) :: tmelt            ! freezing point of water (K)
    REAL(KIND=r8) , INTENT(IN   ) :: hfus            ! latent heat of fusion of water (J kg-1)
    REAL(KIND=r8) , INTENT(IN   ) :: epsilon            ! small quantity to avoid zero-divides and other
    ! truncation or machine-limit troubles with small
    ! values. should be slightly greater than o(1)
    ! machine precision      tmelthfus
    REAL(KIND=r8) , INTENT(IN   ) :: grav          ! gravitational acceleration (m s-2)
    REAL(KIND=r8) , INTENT(IN   ) :: rvap          ! gas constant for water vapor (J deg-1 kg-1)
    REAL(KIND=r8) , INTENT(IN   ) :: dtime         ! model timestep (seconds)
    REAL(KIND=r8) , INTENT(INOUT) :: xu    (npoi)
    REAL(KIND=r8) , INTENT(INOUT) :: xs    (npoi)
    REAL(KIND=r8) , INTENT(INOUT) :: xl    (npoi)
    REAL(KIND=r8) , INTENT(INOUT) :: chux  (npoi) 
    REAL(KIND=r8) , INTENT(INOUT) :: chsx  (npoi) 
    REAL(KIND=r8) , INTENT(INOUT) :: chlx  (npoi)
    REAL(KIND=r8) , INTENT(INOUT) :: chgx  (npoi)
    REAL(KIND=r8) , INTENT(INOUT) :: wlgx  (npoi)
    REAL(KIND=r8) , INTENT(INOUT) :: wigx  (npoi)
    REAL(KIND=r8) , INTENT(INOUT) :: cog   (npoi)
    REAL(KIND=r8) , INTENT(INOUT) :: coi   (npoi)
    REAL(KIND=r8) , INTENT(INOUT) :: zirg  (npoi)
    REAL(KIND=r8) , INTENT(INOUT) :: ziri  (npoi)
    REAL(KIND=r8) , INTENT(INOUT) :: wu    (npoi)
    REAL(KIND=r8) , INTENT(INOUT) :: ws    (npoi)
    REAL(KIND=r8) , INTENT(INOUT) :: wl    (npoi)              
    REAL(KIND=r8) , INTENT(INOUT) :: wg    (npoi)
    REAL(KIND=r8) , INTENT(INOUT) :: wi    (npoi)
    REAL(KIND=r8) , INTENT(INOUT) :: tuold (npoi)
    REAL(KIND=r8) , INTENT(INOUT) :: tsold (npoi)
    REAL(KIND=r8) , INTENT(INOUT) :: tlold (npoi)
    REAL(KIND=r8) , INTENT(INOUT) :: tgold (npoi)
    REAL(KIND=r8) , INTENT(INOUT) :: tiold (npoi)
    REAL(KIND=r8) , INTENT(OUT  ) :: ginvap (npoi)  ! total evaporation rate from all intercepted h2o (kg_h2o m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: gsuvap (npoi)  ! total evaporation rate from surface (snow/soil) (kg_h2o m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: gtrans (npoi)  ! total transpiration rate from all vegetation canopies (kg_h2o m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: gtransu(npoi)  ! transpiration from upper canopy (kg_h2o m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: gtransl(npoi)  ! transpiration from lower canopy (kg_h2o m-2 s-1)
    REAL(KIND=r8) , INTENT(IN   ) :: psurf(npoi)            ! surface pressure (Pa)  &
    REAL(KIND=r8) , INTENT(IN   ) :: ta   (npoi)            ! air temperature (K)  &
    REAL(KIND=r8) , INTENT(IN   ) :: qa   (npoi)            ! specific humidity (kg_h2o/kg_air)  &
    REAL(KIND=r8) , INTENT(IN   ) :: fi     (npoi)          ! fractional snow cover
    REAL(KIND=r8) , INTENT(IN   ) :: consno                 ! thermal conductivity of snow (W m-1 K-1)
    REAL(KIND=r8) , INTENT(IN   ) :: hsno   (npoi,nsnolay)  ! thickness of snow layers (m)
    REAL(KIND=r8) , INTENT(IN   ) :: hsnotop                ! thickness of top snow layer (m)
    REAL(KIND=r8) , INTENT(IN   ) :: tsno   (npoi,nsnolay)  ! temperature of snow layers (K)
    REAL(KIND=r8) , INTENT(IN   ) :: poros   (npoi,nsoilay) ! porosity (mass of h2o per unit vol at sat / rhow)
    REAL(KIND=r8) , INTENT(IN   ) :: wpud    (npoi)         ! liquid content of puddles per soil area (kg m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: wipud   (npoi)         ! ice content of puddles per soil area (kg m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: csoi    (npoi,nsoilay) ! specific heat of soil, no pore spaces (J kg-1 deg-1)
    REAL(KIND=r8) , INTENT(IN   ) :: rhosoi  (npoi,nsoilay) ! soil density (without pores, not bulk) (kg m-3)
    REAL(KIND=r8) , INTENT(IN   ) :: wisoi   (npoi,nsoilay) ! fraction of soil pore space containing ice
    REAL(KIND=r8) , INTENT(IN   ) :: wsoi    (npoi,nsoilay) ! fraction of soil pore space containing liquid water
    REAL(KIND=r8) , INTENT(IN   ) :: hsoi    (npoi,nsoilay+1)! soil layer thickness (m)
    REAL(KIND=r8) , INTENT(IN   ) :: consoi  (npoi,nsoilay) ! thermal conductivity of each soil layer (W m-1 K-1)
    REAL(KIND=r8) , INTENT(INOUT) :: tg      (npoi)         ! soil skin temperature (K)
    REAL(KIND=r8) , INTENT(INOUT) :: ti      (npoi)         ! snow skin temperature (K)
    REAL(KIND=r8) , INTENT(IN   ) :: wpudmax                ! normalization constant for puddles (kg m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: suction (npoi,nsoilay) ! saturated matric potential (m-h2o)
    REAL(KIND=r8) , INTENT(IN   ) :: bex     (npoi,nsoilay) ! exponent "b" in soil water potential
    REAL(KIND=r8) , INTENT(IN   ) :: swilt   (npoi,nsoilay) ! wilting soil moisture value (fraction of pore space)
    REAL(KIND=r8) , INTENT(IN   ) :: hvasug  (npoi)         ! latent heat of vap/subl, for soil surface (J kg-1)
    REAL(KIND=r8) , INTENT(IN   ) :: tsoi    (npoi,nsoilay) ! soil temperature for each layer (K)
    REAL(KIND=r8) , INTENT(IN   ) :: hvasui  (npoi)         ! latent heat of vap/subl, for snow surface (J kg-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: upsoiu  (npoi,nsoilay) ! soil water uptake from transpiration (kg_h2o m-2 s-1)
    REAL(KIND=r8) , INTENT(IN   ) :: stressu (npoi,nsoilay) ! soil moisture stress factor for the upper canopy (dimensionless)
    REAL(KIND=r8) , INTENT(IN   ) :: stresstu(npoi)         ! sum of stressu over all 6 soil layers (dimensionless)
    REAL(KIND=r8) , INTENT(OUT  ) :: upsoil  (npoi,nsoilay) ! soil water uptake from transpiration (kg_h2o m-2 s-1)
    REAL(KIND=r8) , INTENT(IN   ) :: stressl (npoi,nsoilay) ! soil moisture stress factor for the lower canopy (dimensionless)
    REAL(KIND=r8) , INTENT(IN   ) :: stresstl(npoi)         ! sum of stressl over all 6 soil layers (dimensionless)
    REAL(KIND=r8) , INTENT(IN   ) :: lai      (npoi,2) ! canopy single-sided leaf area index (area leaf/area veg)
    REAL(KIND=r8) , INTENT(IN   ) :: fu       (npoi)   ! fraction of overall area covered by upper canopy
    REAL(KIND=r8) , INTENT(IN   ) :: sai      (npoi,2) ! current single-sided stem area index
    REAL(KIND=r8) , INTENT(IN   ) :: fl       (npoi)   ! fraction of snow-free area covered by lower  canopy
    REAL(KIND=r8) , INTENT(IN   ) :: chu  (1:nVegClass)       ! heat capacity of upper canopy leaves per unit leaf area (J kg-1 m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: wliqu    (npoi)   ! intercepted liquid h2o on upper canopy leaf area (kg m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: wsnou    (npoi)   ! intercepted frozen h2o (snow) on upper canopy leaf area (kg m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: chs(1:nVegClass)         ! heat capacity of upper canopy stems per unit stem area (J kg-1 m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: wliqs    (npoi)   ! intercepted liquid h2o on upper canopy stem area (kg m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: wsnos    (npoi)   ! intercepted frozen h2o (snow) on upper canopy stem area (kg m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: chl (1:nVegClass)        ! heat capacity of lower canopy leaves & stems per unit
    ! leaf/stem area (J kg-1 m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: wliql    (npoi)   ! intercepted liquid h2o on lower canopy leaf and stem area (kg m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: wsnol    (npoi)   ! intercepted frozen h2o (snow) on lower canopy leaf & stem area (kg m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: tu       (npoi)   ! temperature of upper canopy leaves (K)
    REAL(KIND=r8) , INTENT(INOUT) :: ts       (npoi)   ! temperature of upper canopy stems (K)
    REAL(KIND=r8) , INTENT(INOUT) :: tl       (npoi)   ! temperature of lower canopy leaves & stems(K)
    REAL(KIND=r8) , INTENT(INOUT) :: q34      (npoi)   ! specific humidity of air at z34
    REAL(KIND=r8) , INTENT(INOUT) :: t34      (npoi)   ! air temperature at z34 (K)
    REAL(KIND=r8) , INTENT(INOUT) :: q12      (npoi)   ! specific humidity of air at z12
    REAL(KIND=r8) , INTENT(IN   ) :: su       (npoi)   ! air-vegetation transfer coefficients (*rhoa) for upper canopy
    ! leaves (m s-1 * kg m-3) (A39a Pollard & Thompson 1995)
    REAL(KIND=r8) , INTENT(IN   ) :: totcondub(npoi)   ! 
    REAL(KIND=r8) , INTENT(IN   ) :: frac     (npoi,npft)   ! fraction of canopy occupied by each plant functional type
    REAL(KIND=r8) , INTENT(IN   ) :: totconduc(npoi)   !
    REAL(KIND=r8) , INTENT(IN   ) :: ss       (npoi)   ! air-vegetation transfer coefficients (*rhoa) for upper canopy
    ! stems (m s-1 * kg m-3) (A39a Pollard & Thompson 1995)
    REAL(KIND=r8) , INTENT(IN   ) :: sl       (npoi)   ! air-vegetation transfer coefficients (*rhoa) for lower canopy
    ! leaves & stems (m s-1*kg m-3) (A39a Pollard & Thompson 1995)
    REAL(KIND=r8) , INTENT(IN   ) :: totcondls(npoi)   ! 
    REAL(KIND=r8) , INTENT(IN   ) :: totcondl4(npoi)   !
    REAL(KIND=r8) , INTENT(IN   ) :: totcondl3(npoi)   !
    REAL(KIND=r8) , INTENT(INOUT) :: t12      (npoi)   ! air temperature at z12 (K)
    REAL(KIND=r8) , INTENT(IN   ) :: cp     (npoi)     ! specific heat of air at za (allowing for h2o vapor) (J kg-1 K-1)
    REAL(KIND=r8) , INTENT(IN   ) :: sg     (npoi)     ! air-soil transfer coefficient
    REAL(KIND=r8) , INTENT(OUT  ) :: fwetux (npoi)     ! fraction of upper canopy leaf area wetted if dew forms
    REAL(KIND=r8) , INTENT(IN   ) :: fwetu  (npoi)     ! fraction of upper canopy leaf area wetted by intercepted liquid and/or snow
    REAL(KIND=r8) , INTENT(OUT  ) :: fwetsx (npoi)     ! fraction of upper canopy stem area wetted if dew forms
    REAL(KIND=r8) , INTENT(IN   ) :: fwets  (npoi)     ! fraction of upper canopy stem area wetted by intercepted liquid and/or snow
    REAL(KIND=r8) , INTENT(OUT  ) :: fwetlx (npoi)     ! fraction of lower canopy leaf and stem area wetted if dew forms
    REAL(KIND=r8) , INTENT(IN   ) :: fwetl  (npoi)     ! fraction of lower canopy stem & leaf area wetted by 
    ! intercepted liquid and/or snow
    REAL(KIND=r8) , INTENT(IN   ) :: solu   (npoi)     ! solar flux (direct + diffuse) absorbed by upper canopy 
    ! leaves per unit canopy area (W m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: firu   (npoi)     ! ir raditaion absorbed by upper canopy leaves (W m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: sols   (npoi)     ! solar flux (direct + diffuse) absorbed by upper canopy
    ! stems per unit canopy area (W m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: firs   (npoi)     ! ir radiation absorbed by upper canopy stems (W m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: soll   (npoi)     ! solar flux (direct + diffuse) absorbed by lower canopy
    ! leaves and stems per unit canopy area (W m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: firl   (npoi)     ! ir radiation absorbed by lower canopy leaves and stems (W m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: rliqu  (npoi)     ! proportion of fwetu due to liquid
    REAL(KIND=r8) , INTENT(IN   ) :: rliqs  (npoi)     ! proportion of fwets due to liquid
    REAL(KIND=r8) , INTENT(IN   ) :: rliql  (npoi)     ! proportion of fwetl due to liquid
    REAL(KIND=r8) , INTENT(IN   ) :: pfluxu (npoi)     ! heat flux on upper canopy leaves due to intercepted h2o (W m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: pfluxs (npoi)     ! heat flux on upper canopy stems due to intercepted h2o (W m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: pfluxl (npoi)     ! heat flux on lower canopy leaves & stems due to intercepted h2o (W m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: cu     (npoi)     ! air transfer coefficient (*rhoa) (m s-1 kg m-3) for upper air
    ! region (z12 --> za) (A35 Pollard & Thompson 1995)
    REAL(KIND=r8) , INTENT(IN   ) :: cl     (npoi)     ! air transfer coefficient (*rhoa) (m s-1 kg m-3) between the
    ! 2 canopies (z34 --> z12) (A36 Pollard & Thompson 1995)
    REAL(KIND=r8) , INTENT(INOUT) :: bps    (npoi)           ! (ps/p) ** (rair/cair) for atmospheric level  (const)
    REAL(KIND=r8) , INTENT(IN   ) :: si     (npoi)     ! air-snow transfer coefficient
    REAL(KIND=r8) , INTENT(IN   ) :: solg   (npoi)     ! solar flux (direct + diffuse) absorbed by unit snow-free soil (W m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: firg   (npoi)     ! ir radiation absorbed by soil/ice (W m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: soli   (npoi)     ! solar flux (direct + diffuse) absorbed by unit snow surface (W m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: firi   (npoi)     ! ir radiation absorbed by snow (W m-2)
    REAL(KIND=r8) , INTENT(OUT  ) :: fsena  (npoi)     ! downward sensible heat flux between za & z12 at za (W m-2)
    REAL(KIND=r8) , INTENT(OUT  ) :: fseng  (npoi)     ! upward sensible heat flux between soil surface & air at z34 (W m-2)
    REAL(KIND=r8) , INTENT(OUT  ) :: fseni  (npoi)     ! upward sensible heat flux between snow surface & air at z34 (W m-2)
    REAL(KIND=r8) , INTENT(OUT  ) :: fsenu  (npoi)     ! sensible heat flux from upper canopy leaves to air (W m-2)
    REAL(KIND=r8) , INTENT(OUT  ) :: fsens  (npoi)     ! sensible heat flux from upper canopy stems to air (W m-2)
    REAL(KIND=r8) , INTENT(OUT  ) :: fsenl  (npoi)     ! sensible heat flux from lower canopy to air (W m-2)
    REAL(KIND=r8) , INTENT(OUT  ) :: fvapa  (npoi)     ! downward h2o vapor flux between za & z12 at za (kg m-2 s-1)
    REAL(KIND=r8) , INTENT(OUT  ) :: fvapuw (npoi)     ! h2o vapor flux (evaporation from wet parts) between upper canopy
    ! leaves and air at z12 (kg m-2 s-1/ LAI upper canopy/ fu)
    REAL(KIND=r8) , INTENT(OUT  ) :: fvaput (npoi)     ! h2o vapor flux (transpiration from dry parts) between upper canopy 
    ! leaves and air at z12 (kg m-2 s-1/ LAI upper canopy/ fu)
    REAL(KIND=r8) , INTENT(OUT  ) :: fvaps  (npoi)     ! h2o vapor flux (evaporation from wet surface) between upper canopy 
    ! stems and air at z12 (kg m-2 s-1 / SAI lower canopy / fu)
    REAL(KIND=r8) , INTENT(OUT  ) :: fvaplw (npoi)     ! h2o vapor flux (evaporation from wet surface) between lower canopy
    ! leaves & stems and air at z34 (kg m-2 s-1/ LAI lower canopy/ fl)
    REAL(KIND=r8) , INTENT(OUT  ) :: fvaplt (npoi)     ! h2o vapor flux (transpiration) between lower canopy &
    ! air at z34 (kg m-2 s-1 / LAI lower canopy / fl)
    REAL(KIND=r8) , INTENT(OUT  ) :: fvapg  (npoi)     ! h2o vapor flux (evaporation) between soil & air
    ! at z34 (kg m-2 s-1/bare ground fraction)
    REAL(KIND=r8) , INTENT(OUT  ) :: fvapi  (npoi)     ! h2o vapor flux (evaporation) between snow & air at z34 (kg m-2 s-1 / fi )
    REAL(KIND=r8) , INTENT(INOUT) :: firb   (npoi)     ! net upward ir radiation at reference atmospheric level za (W m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: vegtype0(npoi)  
    !
    ! Arguments (input)
    !
    INTEGER, INTENT(IN   ) :: niter       ! total # of iteration
    INTEGER, INTENT(IN   ) :: iter        ! # of iteration
    !
    ! local variables
    !
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER :: iveg

    !
    REAL(KIND=r8) :: rwork
    REAL(KIND=r8) :: zwtot 
    REAL(KIND=r8) :: rwork2 
    REAL(KIND=r8) :: tgav 
    REAL(KIND=r8) :: tiav 
    REAL(KIND=r8) :: tuav 
    REAL(KIND=r8) :: tsav 
    REAL(KIND=r8) :: tlav 
    REAL(KIND=r8) :: quav 
    REAL(KIND=r8) :: qsav 
    REAL(KIND=r8) :: qlav 
    REAL(KIND=r8) :: qgav 
    REAL(KIND=r8) :: qiav 
    REAL(KIND=r8) :: zwpud 
    REAL(KIND=r8) :: zwsoi
    REAL(KIND=r8) :: psig 
    REAL(KIND=r8) :: hfac 
    REAL(KIND=r8) :: hfac2 
    REAL(KIND=r8) :: zwopt 
    REAL(KIND=r8) :: zwdry 
    REAL(KIND=r8) :: betaw 
    REAL(KIND=r8) :: emisoil 
    REAL(KIND=r8) :: e 
    REAL(KIND=r8) :: qs1
    REAL(KIND=r8) :: dqs1 
    REAL(KIND=r8) :: xnumer 
    REAL(KIND=r8) :: xdenom 
    REAL(KIND=r8) :: betafac 
    REAL(KIND=r8) :: betas
    !


    REAL(KIND=r8) :: fradu(npoi) 
    REAL(KIND=r8) :: frads(npoi) 
    REAL(KIND=r8) :: fradl(npoi)       
    REAL(KIND=r8) :: qu(npoi)
    REAL(KIND=r8) :: qs(npoi)  
    REAL(KIND=r8) :: ql(npoi)  
    REAL(KIND=r8) :: qg(npoi)  
    REAL(KIND=r8) :: qi(npoi)
    REAL(KIND=r8) :: dqu(npoi)
    REAL(KIND=r8) :: dqs(npoi) 
    REAL(KIND=r8) :: dql(npoi) 
    REAL(KIND=r8) :: dqg(npoi) 
    REAL(KIND=r8) :: dqi(npoi)

    REAL(KIND=r8) :: suw(npoi)
    REAL(KIND=r8) :: ssw(npoi)
    REAL(KIND=r8) :: slw(npoi)
    REAL(KIND=r8) :: sut(npoi)
    REAL(KIND=r8) :: slt(npoi)
    REAL(KIND=r8) :: slt0(npoi)
    REAL(KIND=r8) :: suh(npoi)
    REAL(KIND=r8) :: ssh(npoi) 
    REAL(KIND=r8) :: slh(npoi)
    REAL(KIND=r8) :: qgfac(npoi)
    REAL(KIND=r8) :: qgfac0(npoi)
    REAL(KIND=r8) :: tupre(npoi)
    REAL(KIND=r8) :: tspre(npoi)
    REAL(KIND=r8) :: tlpre(npoi)
    REAL(KIND=r8) :: tgpre(npoi)
    REAL(KIND=r8) :: tipre(npoi)
    !
    INTEGER, PARAMETER :: nqn=9
    !
    !
    REAL(KIND=r8) :: arr(npoi,nqn,nqn)      !    
    REAL(KIND=r8) :: rhs(npoi,nqn)          ! right hand side
    REAL(KIND=r8) :: vec(npoi,nqn)          ! 
    REAL(KIND=r8) :: levapr          ! 
    REAL(KIND=r8) :: hearcvap
    !
    INTEGER, PARAMETER :: mplate(1:nqn,1:nqn) = RESHAPE ((/&
         !                  tu  ts  tl t12 t34 q12 q34  tg  ti 
         !                  ----------------------------------
         1,  0,  0,  1,  0,  1,  0,  0,  0,&!tu
         0,  1,  0,  1,  0,  1,  0,  0,  0,&!ts
         0,  0,  1,  0,  1,  0,  1,  0,  0,&!tl
         1,  1,  0,  1,  1,  0,  0,  0,  0,&!t12
         0,  0,  1,  1,  1,  0,  0,  1,  1,&!t34
         1,  1,  0,  0,  0,  1,  1,  0,  0,&!q12
         0,  0,  1,  0,  0,  1,  1,  1,  1,&!q34
         0,  0,  0,  0,  1,  0,  1,  1,  0,&!tg
         0,  0,  0,  0,  1,  0,  1,  0,  1 &!ti
         /), (/nqn, nqn/) )
    !
    !include 'comsat.h'

    !
    ! ---------------------------------------------------------------------
    ! statement functions and associated parameters
    ! ---------------------------------------------------------------------
    !
    ! polynomials for svp(t), d(svp)/dt over water and ice are from
    ! lowe(1977),jam,16,101-103.
    !
    !
    REAL(KIND=r8) , PARAMETER :: asat0 =  6.1078000_r8
    REAL(KIND=r8) , PARAMETER :: asat1 =  4.4365185e-1_r8
    REAL(KIND=r8) , PARAMETER :: asat2 =  1.4289458e-2_r8
    REAL(KIND=r8) , PARAMETER :: asat3 =  2.6506485e-4_r8
    REAL(KIND=r8) , PARAMETER :: asat4 =  3.0312404e-6_r8
    REAL(KIND=r8) , PARAMETER :: asat5 =  2.0340809e-8_r8
    REAL(KIND=r8) , PARAMETER :: asat6 =  6.1368209e-11_r8 
    !
    !
    REAL(KIND=r8) , PARAMETER :: bsat0 =  6.1091780_r8
    REAL(KIND=r8) , PARAMETER :: bsat1 =  5.0346990e-1_r8
    REAL(KIND=r8) , PARAMETER :: bsat2 =  1.8860134e-2_r8
    REAL(KIND=r8) , PARAMETER :: bsat3 =  4.1762237e-4_r8
    REAL(KIND=r8) , PARAMETER :: bsat4 =  5.8247203e-6_r8
    REAL(KIND=r8) , PARAMETER :: bsat5 =  4.8388032e-8_r8
    REAL(KIND=r8) , PARAMETER :: bsat6 =  1.8388269e-10_r8
    !
    !
    REAL(KIND=r8) , PARAMETER :: csat0 =  4.4381000e-1_r8
    REAL(KIND=r8) , PARAMETER :: csat1 =  2.8570026e-2_r8
    REAL(KIND=r8) , PARAMETER :: csat2 =  7.9380540e-4_r8
    REAL(KIND=r8) , PARAMETER :: csat3 =  1.2152151e-5_r8
    REAL(KIND=r8) , PARAMETER :: csat4 =  1.0365614e-7_r8
    REAL(KIND=r8) , PARAMETER :: csat5 =  3.5324218e-10_r8
    REAL(KIND=r8) , PARAMETER :: csat6 = -7.0902448e-13_r8
    !
    !
    REAL(KIND=r8) , PARAMETER :: dsat0 =  5.0303052e-1_r8
    REAL(KIND=r8) , PARAMETER :: dsat1 =  3.7732550e-2_r8
    REAL(KIND=r8) , PARAMETER :: dsat2 =  1.2679954e-3_r8
    REAL(KIND=r8) , PARAMETER :: dsat3 =  2.4775631e-5_r8
    REAL(KIND=r8) , PARAMETER :: dsat4 =  3.0056931e-7_r8
    REAL(KIND=r8) , PARAMETER :: dsat5 =  2.1585425e-9_r8
    REAL(KIND=r8) , PARAMETER :: dsat6 =  7.1310977e-12_r8
    levapr=hsub
    hearcvap=cvap
    !
    ! statement functions tsatl,tsati are used below so that lowe's
    ! polyomial for liquid is used if t gt 273.16, or for ice if 
    ! t lt 273.16. also impose range of validity for lowe's polys.
    !
    !      REAL(KIND=r8) :: t        ! temperature argument of statement function 
    !      REAL(KIND=r8) :: tair      ! temperature argument of statement function 
    !      REAL(KIND=r8) :: p1       ! pressure argument of function 
    !      REAL(KIND=r8) :: e1       ! vapor pressure argument of function
    !      REAL(KIND=r8) :: q1       ! saturation specific humidity argument of function
    !REAL(KIND=r8) :: tsatl     ! statement function
    !REAL(KIND=r8) :: tsati     ! 
    !REAL(KIND=r8) :: esat     !
    !REAL(KIND=r8) :: desat    !
    !REAL(KIND=r8) :: qsat     ! 
    !REAL(KIND=r8) :: dqsat    ! 
    !REAL(KIND=r8) :: hvapf    ! 
    !REAL(KIND=r8) :: hsubf    !
    !REAL(KIND=r8) :: cvmgt    ! function
    !
    !tsatl(t) = min (100., max (t-273.16, 0.))
    !tsati(t) = max (-60., min (t-273.16, 0.))
    !
    ! statement function esat is svp in n/m**2, with t in deg k. 
    ! (100 * lowe's poly since 1 mb = 100 n/m**2.)
    !
    ! esat (t) =              &
    !  100.*(             &
    !         cvmgt (asat0, bsat0, t.ge.273.16)             &
    !         + tsatl(t)*(asat1 + tsatl(t)*(asat2 + tsatl(t)*(asat3             &
    !         + tsatl(t)*(asat4 + tsatl(t)*(asat5 + tsatl(t)* asat6)))))  &
    !         + tsati(t)*(bsat1 + tsati(t)*(bsat2 + tsati(t)*(bsat3             &
    !         + tsati(t)*(bsat4 + tsati(t)*(bsat5 + tsati(t)* bsat6)))))  &
    ! )
    !
    ! statement function desat is d(svp)/dt, with t in deg k.
    ! (100 * lowe's poly since 1 mb = 100 n/m**2.)
    !
    !desat (t) =             &
    ! 100.*(             &
    !         cvmgt (csat0, dsat0, t.ge.273.16)             &
    !         + tsatl(t)*(csat1 + tsatl(t)*(csat2 + tsatl(t)*(csat3             &
    !          + tsatl(t)*(csat4 + tsatl(t)*(csat5 + tsatl(t)* csat6)))))  &
    !         + tsati(t)*(dsat1 + tsati(t)*(dsat2 + tsati(t)*(dsat3             &
    !         + tsati(t)*(dsat4 + tsati(t)*(dsat5 + tsati(t)* dsat6)))))  &
    ! )
    !
    ! statement function qsat is saturation specific humidity,
    ! with svp e1 and ambient pressure p in n/m**2. impose an upper
    ! limit of 1 to avoid spurious values for very high svp
    ! and/or small p1
    !
    !       qsat (e1, p1) = 0.622 * e1 /  &
    !                     max ( p1 - (1.0 - 0.622) * e1, 0.622 * e1 )
    !
    ! statement function dqsat is d(qsat)/dt, with t in deg k and q1
    ! in kg/kg (q1 is *saturation* specific humidity)
    !
    !       dqsat (t, q1) = desat(t) * q1 * (1. + q1*(1./0.622 - 1.)) / &
    !                       esat(t)
    !
    ! statement functions hvapf, hsubf correct the latent heats of
    ! vaporization (liquid-vapor) and sublimation (ice-vapor) to
    ! allow for the concept that the phase change takes place at
    ! 273.16, and the various phases are cooled/heated to that 
    ! temperature before/after the change. this concept is not
    ! physical but is needed to balance the "black-box" energy 
    ! budget. similar correction is applied in convad in the agcm
    ! for precip. needs common comgrd for the physical constants.
    ! argument t is the temp of the liquid or ice, and tair is the
    ! temp of the delivered or received vapor.
    !
    !              hvapf(t,tair) = hvap + cvap*(tair-273.16) - ch2o*(t-273.16)
    !      hsubf(t,tair) = hsub + cvap*(tair-273.16) - cice*(t-273.16)
    !


    !
    ! if first iteration, save original canopy temps in t*old
    ! (can use tsoi,tsno for original soil/snow skin temps), for
    ! rhs heat capacity terms in matrix soln, and for adjustment
    ! of canopy temps after each iteration
    !
    ! also initialize soil/snow skin temps tg, ti to top-layer temps
    !
    ! the variables t12, t34, q12, q34, for the first iteration
    ! are saved via global arrays from the previous gcm timestep,
    ! this is worth doing only if the agcm forcing is
    ! smoothly varying from timestep to timestep
    !
    arr=0.0_r8
    rhs=0.0_r8
    vec=0.0_r8
    IF (iter.EQ.1) THEN
       !
       ! weights for canopy coverages
       !
       DO  i = 1, npoi
          xu(i) = 2.0_r8 * lai(i,2) * fu(i)
          xs(i) = 2.0_r8 * sai(i,2) * fu(i)
          xl(i) = 2.0_r8 * (lai(i,1) + sai(i,1)) * fl(i) * (1.0_r8 - fi(i))
       END DO
       !
       ! specific heats per leaf/stem area
       !
       DO i = 1, npoi
          iveg=vegtype0(i)
          chux(i) = chu(iveg) + ch2o * wliqu(i) + cice * wsnou(i)
          chsx(i) = chs(iveg) + ch2o * wliqs(i) + cice * wsnos(i)
          chlx(i) = chl(iveg) + ch2o * wliql(i) + cice * wsnol(i)
       END DO
       !
       DO i = 1, npoi 
          !
          rwork = poros(i,1) * rhow
          !
          chgx(i) = ch2o * wpud(i) + cice * wipud(i)              &
               + ((1.0_r8-poros(i,1))*csoi(i,1)*rhosoi(i,1)   &
               + rwork*(1.0_r8-wisoi(i,1))*wsoi(i,1)*ch2o     &
               + rwork*wisoi(i,1)*cice              &
               ) * hsoi(i,1)
          !
          wlgx(i) = wpud(i) +                 &
               rwork * (1.0_r8 - wisoi(i,1)) *         &
               wsoi(i,1) * hsoi(i,1)
          !
          wigx(i) = wipud(i) + rwork * wisoi(i,1) * hsoi(i,1)
          !
       END DO
       !
       ! conductivity coeffs between ground skin and first layer
       !
       DO i = 1, npoi
          cog(i) = consoi(i,1) / (0.5_r8 * hsoi(i,1))
          coi(i) = consno      / (0.5_r8 * MAX (hsno(i,1), hsnotop))
       END DO
       !
       ! d(ir emitted) / dt for soil
       !
       rwork = 4.0_r8 * 0.95_r8 * stef
       !
       DO i = 1, npoi
          zirg(i) = rwork * (tg(i)**3)
          ziri(i) = rwork * (ti(i)**3)
       END DO
       !
       ! updated temperature memory
       !
       DO i = 1, npoi
          tuold(i) = tu(i)
          tsold(i) = ts(i)
          tlold(i) = tl(i)
          tgold(i) = tg(i)
          tiold(i) = ti(i)
       END DO
       !
    END IF
    !
    ! set implicit/explicit factors w* (0 to 1) for this iteration
    ! w* is 1 for fully implicit, 0 for fully explicit
    ! for first iteration, impexp and impexp2 set w* to 1
    !
    CALL impexp (wu       , &! INTENT(INOUT)
         tu       , &! INTENT(IN   )
         chux     , &! INTENT(IN   )
         wliqu    , &! INTENT(IN   )
         wsnou    , &! INTENT(IN   )
         iter     , &! INTENT(IN   )
         npoi     , &! INTENT(IN   )
         tmelt    , &! INTENT(IN   )
         hfus     , &! INTENT(IN   )
         epsilon    )! INTENT(IN   )

    CALL impexp (ws       , &! INTENT(INOUT)
         ts       , &! INTENT(IN   )
         chsx     , &! INTENT(IN   )
         wliqs    , &! INTENT(IN   )
         wsnos    , &! INTENT(IN   )
         iter     , &! INTENT(IN   )
         npoi     , &! INTENT(IN   )
         tmelt    , &! INTENT(IN   )
         hfus     , &! INTENT(IN   )
         epsilon    )! INTENT(IN   )

    CALL impexp (wl       , &! INTENT(INOUT)
         tl       , &! INTENT(IN   )
         chlx     , &! INTENT(IN   )
         wliql    , &! INTENT(IN   )
         wsnol    , &! INTENT(IN   )
         iter     , &! INTENT(IN   )
         npoi     , &! INTENT(IN   )
         tmelt    , &! INTENT(IN   )
         hfus     , &! INTENT(IN   )
         epsilon    )! INTENT(IN   )

    CALL impexp (wg       , &! INTENT(INOUT)
         tg       , &! INTENT(IN   )
         chgx     , &! INTENT(IN   )
         wlgx     , &! INTENT(IN   )
         wigx     , &! INTENT(IN   )
         iter     , &! INTENT(IN   )
         npoi     , &! INTENT(IN   )
         tmelt    , &! INTENT(IN   )
         hfus     , &! INTENT(IN   )
         epsilon    )! INTENT(IN   )
    !
    ! call impexp2 for snow model
    !
    CALL impexp2 (wi      , &! INTENT(INOUT)
         ti      , &! INTENT(IN   )
         tiold   , &! INTENT(IN   )
         iter    , &! INTENT(IN   )
         npoi    , &! INTENT(IN   )
         tmelt   , &! INTENT(IN   )
                                !hfus    , &! INTENT(IN   )
         epsilon   )! INTENT(IN   )
    !
    ! adjust t* for this iteration 
    !
    ! in this routine we are free to choose them, 
    ! since they are just the central values about which the 
    ! equations are linearized - heat is conserved in the matrix
    ! solution because t*old are used for the rhs heat capacities
    !
    ! here, let t* represent the previous soln if it was fully
    ! implicit, but weight towards t*old depending on the amount
    ! (1-w*) the previous soln was explicit
    !
    ! this weighting is necessary for melting/freezing surfaces, for which t*
    ! is kept at t*old, presumably at or near tmelt
    !
    DO i = 1, npoi
       tu(i) = wu(i) * tu(i) + (1.0_r8 - wu(i)) * tuold(i)
       ts(i) = ws(i) * ts(i) + (1.0_r8 - ws(i)) * tsold(i)
       tl(i) = wl(i) * tl(i) + (1.0_r8 - wl(i)) * tlold(i)
       tg(i) = wg(i) * tg(i) + (1.0_r8 - wg(i)) * tgold(i)
       ti(i) = wi(i) * ti(i) + (1.0_r8 - wi(i)) * tiold(i)
    END DO
    !
    ! save current "central" values for final flux calculations
    !
    DO i = 1, npoi
       tupre(i) = tu(i)
       tspre(i) = ts(i)
       tlpre(i) = tl(i)
       tgpre(i) = tg(i)
       tipre(i) = ti(i)
    END DO
    !
    ! calculate various terms occurring in the linearized eqns,
    ! using values of t12, t34, q12, q34 from
    ! the previous iteration
    !
    ! specific humidities for canopy and ground, and derivs wrt t
    ! for canopy
    !
    ! limit derivs to avoid -ve implicit q's below,
    ! as long as d(temp)s in one iteration are le 10 deg k
    !
    DO i = 1, npoi
       !
       !       PRINT*,'esat= ', esat(tu(i)), 'es_Sat=  ', es_Sat   (tu(i), psurf(i))
       !       PRINT*,'desat=', desat(tu(i)),'esdT_Sat=', esdT_Sat (tu(i), psurf(i))
       !         e      = esat(tu(i))
       !        qu(i)  = qsat (e, psurf(i))
       !       PRINT*,'qsat= ', qsat(e,psurf(i)), 'qs_Sat=  ', qs_Sat   (tu(i), psurf(i))
       !       PRINT*,'dqsat=', dqsat(tu(i),qu(i)),'qsdT_Sat=', qsdT_Sat (tu(i), psurf(i))

       e      = esat(tu(i))
       qu(i)  = qsat (e, psurf(i))
       dqu(i) = dqsat (tu(i), qu(i))
       dqu(i) = MIN (dqu(i), qu(i) * 0.1_r8)
       !
       e      = esat(ts(i))
       qs(i)  = qsat (e, psurf(i))
       dqs(i) = dqsat (ts(i), qs(i))
       dqs(i) = MIN (dqs(i), qs(i) * 0.1_r8)
       !
       e      = esat(tl(i))
       ql(i)  = qsat (e, psurf(i))
       dql(i) = dqsat (tl(i), ql(i))
       dql(i) = MIN (dql(i), ql(i) * 0.1_r8)
       !
       e      = esat(tg(i))
       qg(i)  = qsat (e, psurf(i))
       dqg(i) = dqsat (tg(i), qg(i))
       dqg(i) = MIN (dqg(i), qg(i) * 0.1_r8)
       !
       e      = esat(ti(i))
       qi(i)  = qsat (e, psurf(i))
       dqi(i) = dqsat (ti(i), qi(i))
       dqi(i) = MIN (dqi(i), qi(i) * 0.1_r8)
       !
    END DO
    !
    ! set qgfac0, factor by which soil surface specific humidity
    ! is less than saturation
    !
    ! it is important to note that the qgfac expression should
    ! satisfy timestep cfl criterion for upper-layer soil moisture
    ! for small wsoi(i,1)
    !
    ! for each iteration, qgfac is set to qgfac0, or to 1 if
    ! condensation onto soil is anticipated (loop 110 in canopy.f)
    !
    ! Evaporation from bare soil is calculated using the "beta method"
    ! (e.g., eqns 5 & 7 of Mahfouf and Noilhan 1991, JAM 30 1354-1365),
    ! but converted to the "alpha method" (eqns 2 & 3 of M&N), to match
    ! the structure in IBIS. The conversion from the beta to alpha
    ! method is through the relationship:
    !   alpha * qgs - q34 = beta * (hfac * qgs - q34),
    ! from which one solves for alpha (which is equal to qgfac0):
    !   qgfac0 = alpha = (beta * hfac) + (1 - beta)*(q34/qgs)
    !
    DO i = 1, npoi
       !
       ! first calculate the total saturated fraction at the soil surface
       ! (including puddles ... see soil.f)
       !
       zwpud = MAX (0.0_r8, MIN (0.5_r8, 0.5_r8*(wpud(i)+wipud(i))/wpudmax) )
       zwsoi = (1.0_r8 - wisoi(i,1)) * wsoi(i,1) + wisoi(i,1)
       zwtot = zwpud + (1.0_r8 - zwpud) * zwsoi
       zwtot = MAX(zwtot,1.0e-12_r8)
       !
       ! next calculate the matric potential (from eqn 9.3 of Campbell and
       ! Norman), multiply by gravitational acceleration to get in units
       ! of J/kg, and calculate the relative humidity at the soil water
       ! surface (i.e., within the soil matrix), based on thermodynamic
       ! theory (eqn 4.13 of C&N)
       !
       psig = -grav * suction(i,1) * (zwtot ** (-bex(i,1)))
       hfac = EXP(psig/(rvap*tg(i)))
       !
       ! then calculate the relative humidity of the air (relative to
       ! saturation at the soil temperature). Note that if hfac2 > 1
       ! (which would imply condensation), then qgfac is set to 1
       ! later in the code (to allow condensation to proceed at the
       ! "potential rate")
       !
       hfac2 = q34(i)/qg(i)
       !
       ! set the "beta" factor and then calculate "alpha" (i.e., qgfac0)
       ! as the beta-weighted average of the soil water RH and the "air RH"
       ! First calculate beta_w:
       !
       zwopt = 1.0_r8
       zwdry = swilt(i,1)
       betaw = MAX(0.0_r8, MIN(1.0_r8, (zwtot - zwdry)/(zwopt - zwdry)) )
       !
       ! Next convert beta_w to beta_s (see Milly 1992, JClim 5 209-226):
       !
       emisoil = 0.95_r8
       e      = esat(t34(i))
       qs1    = qsat (e, psurf(i))
       dqs1   = dqsat (t34(i), qs1)
       xnumer = hvap * dqs1
       xdenom = cp(i) + (4.0_r8 * emisoil * stef * (t34(i))**3) / sg(i)
       betafac = xnumer / xdenom
       betas = betaw / (1.0_r8 + betafac * (1.0_r8 - betaw))
       !
       ! Combine hfac and hfac2 into qgfac0 ("alpha") using beta_s
       !
       qgfac0(i) = betas * hfac + (1.0_r8 - betas) * hfac2
    END DO
    !
    ! set fractions covered by intercepted h2o to 1 if dew forms
    !
    ! these fwet*x are used only in turvap, and are distinct from
    ! the real fractions fwet* that are set in fwetcal
    !
    ! they must be exactly 1 if q12 > qu or q34 > ql, to zero transpiration
    ! by the factor 1-fwet[u,l]x below, so preventing "-ve" transp
    !
    ! similarly, set qgfac, allowing for anticipated dew formation
    ! to avoid excessive dew formation (which then infiltrates) onto
    ! dry soils
    !
    DO i = 1, npoi
       !
       fwetux(i) = fwetu(i)
       IF (q12(i).GT.qu(i)) fwetux(i) = 1.0_r8
       !
       fwetsx(i) = fwets(i)
       IF (q12(i).GT.qs(i)) fwetsx(i) = 1.0_r8
       !
       fwetlx(i) = fwetl(i)
       IF (q34(i).GT.ql(i)) fwetlx(i) = 1.0_r8
       !
       qgfac(i) = qgfac0(i)
       IF (q34(i).GT.qg(i)) qgfac(i) = 1.0_r8
       !
       ! set net absorbed radiative fluxes for canopy components
       !
       fradu(i) = 0.0_r8
       !
       IF (lai(i,2).GT.epsilon) &
            fradu(i) = (solu(i) + firu(i)) / (2.0_r8 * lai(i,2))
       !
       frads(i) = 0.0_r8
       !
       IF (sai(i,2).GT.epsilon) &
            frads(i) = (sols(i) + firs(i)) / (2.0_r8 * sai(i,2))
       !
       fradl(i) = 0.0_r8
       !
       IF ((lai(i,1)+sai(i,1)).GT.epsilon)   &
            fradl(i) = (soll(i) + firl(i)) /   &
            (2.0_r8 * (lai(i,1) + sai(i,1)))
       !
    END DO
    !
    ! calculate canopy-air moisture transfer coeffs for wetted
    ! leaf/stem areas, and for dry (transpiring) leaf areas
    !
    ! the wetted-area coeffs suw,ssw,slw are constrained to be less
    ! than what would evaporate 0.8 * the intercepted h2o mass in 
    ! this timestep (using previous iteration's q* values)
    !
    ! this should virtually eliminate evaporation-overshoots and the need
    ! for the "negative intercepted h2o"  correction in steph2o2
    !        
    DO  i = 1, npoi
       !
       ! coefficient for evaporation from wet surfaces in the upper canopy:
       !
       suw(i) = MIN ( fwetux(i) * su(i),  &
            0.8_r8 * (wliqu(i) + wsnou(i)) /   &
            MAX (dtime * (qu(i) - q12(i)), epsilon))
       !
       ! coefficient for transpiration from average upper canopy leaves:
       !
       sut(i) = (1.0_r8 - fwetux(i)) * 0.5_r8 *     &
            ( totcondub(i) * frac(i,1) +   &
            totcondub(i) * frac(i,2) +   &
            totcondub(i) * frac(i,3) +   &
            totconduc(i) * frac(i,4) +   &
            totcondub(i) * frac(i,5) +   &
            totconduc(i) * frac(i,6) +   &
            totcondub(i) * frac(i,7) +   &
            totcondub(i) * frac(i,8) )
       !
       sut(i) = MAX (0.0_r8, sut(i))
       !
       ! coefficient for sensible heat flux from upper canopy:
       !
       suh(i) = suw(i) * (rliqu(i)  * hvapf(tu(i),ta(i))  +  &
            (1.0_r8-rliqu(i)) * hsubf(tu(i),ta(i))) +  &
            sut(i) *              hvapf(tu(i),ta(i))
       !
       ! coefficient for evaporation from wet surfaces on the stems:
       !
       ssw(i) = MIN (fwetsx(i) * ss(i),    &
            0.8_r8 * (wliqs(i) + wsnos(i))    &
            / MAX (dtime * (qs(i) - q12(i)), epsilon))
       !
       ! coefficient for sensible heat flux from stems:
       !
       ssh(i) = ssw(i) * (rliqs(i)  * hvapf(ts(i),ta(i)) +  &
            (1.0_r8-rliqs(i)) * hsubf(ts(i),ta(i)))
       !
       ! coefficient for evaporation from wet surfaces in the lower canopy:
       !
       slw(i) = MIN (fwetlx(i) * sl(i),         &
            0.8_r8 * (wliql(i) + wsnol(i))       &
            / MAX (dtime * (ql(i) - q34(i)), epsilon))
       !
       ! coefficient for transpiration from average lower canopy leaves:
       !
       !       WRITE(*,*) fwetlx(i),totcondls(i),totcondl4(i),totcondl3(i),frac(i,9),frac(i,10)

       slt0(i) = (1.0_r8 - fwetlx(i)) * 0.5_r8 *       &
            ( totcondls(i) * frac(i,9)  +  & 
            totcondls(i) * frac(i,10) +  &
            totcondl4(i) * frac(i,11) +  &
            totcondl3(i) * frac(i,12) )
       !
       slt0(i) = MAX (0.0_r8, slt0(i))
       !
       ! averaged over stems and lower canopy leaves:
       ! 
       slt(i) = slt0(i) * lai(i,1) / MAX (lai(i,1)+sai(i,1), epsilon)
       !
       ! coefficient for sensible heat flux from lower canopy:
       !
       slh(i) = slw(i) * (  rliql(i)  * hvapf(tl(i),ta(i))  +  &
            (1.0_r8-rliql(i)) * hsubf(tl(i),ta(i))) +  &
            slt(i) *                hvapf(tl(i),ta(i))
       !
    END DO
    !
    ! set the matrix of coefficients and the right-hand sides
    ! of the linearized equations
    !

    !
    DO k=1,nqn
       DO j=1,nqn
          DO i=1,npoi
             arr(i,j,k) = 0.0_r8
          END DO
       END DO
    END DO
    DO k=1,nqn     
       DO i=1,npoi
          rhs(i,k)  =   0.0_r8   
       END DO
    END DO

    !      CALL const(arr, npoi*nqn*nqn, 0.0_r8)
    !      CALL const(rhs, npoi*nqn, 0.0_r8)
    !
    rwork = 1.0_r8 / dtime
    !
    ! upper leaf temperature tu
    !
    DO i = 1, npoi
       !
       rwork2 = su(i)*cp(i)
       arr(i,1,1) = chux(i)*rwork &
            + wu(i)*rwork2   &
            + wu(i)*suh(i)*dqu(i)
       arr(i,1,4) = -rwork2
       arr(i,1,6) = -suh(i)
       rhs(i,1) = tuold(i)*chux(i)*rwork   &
            - (1.0_r8-wu(i))*rwork2*tu(i)   &
            - suh(i) * (qu(i)-wu(i)*dqu(i)*tu(i))  &
            + fradu(i) - pfluxu(i)
       ! 
    END DO
    !
    ! upper stem temperature ts
    !
    DO i = 1, npoi
       !
       rwork2 = ss(i)*cp(i)
       arr(i,2,2) = chsx(i)*rwork   &
            + ws(i)*rwork2     &
            + ws(i)*ssh(i)*dqs(i)
       arr(i,2,4) = -rwork2
       arr(i,2,6) = -ssh(i)
       rhs(i,2) = tsold(i)*chsx(i)*rwork       &
            - (1.0_r8-ws(i))*rwork2*ts(i)       &
            - ssh(i) * (qs(i)-ws(i)*dqs(i)*ts(i))   &
            + frads(i) - pfluxs(i)
       !
    END DO
    !
    ! lower veg temperature tl
    !
    DO i = 1, npoi
       !
       !       WRITE(*,*) si(i),cp(i),chlx(i),wl(i),slh(i),dql(i) 

       rwork2 = sl(i)*cp(i)
       arr(i,3,3) = chlx(i)*rwork            &
            + wl(i)*rwork2             &
            + wl(i)*slh(i)*dql(i)
       arr(i,3,5) = -rwork2
       arr(i,3,7) = -slh(i)
       rhs(i,3) = tlold(i)*chlx(i)*rwork     &
            - (1.0_r8-wl(i))*rwork2*tl(i)    &
            - slh(i) * (ql(i)-wl(i)*dql(i)*tl(i))  &
            + fradl(i) - pfluxl(i)

       !       WRITE(*,*) arr(i,3,3),arr(i,3,4),arr(i,3,7),arr(i,3,5),arr(i,3,8),arr(i,3,9) 

       !
    END DO
    !
    ! upper air temperature t12
    !
    DO i = 1, npoi
       !
       rwork = xu(i)*su(i)
       rwork2 = xs(i)*ss(i)
       arr(i,4,1) = -wu(i)*rwork
       arr(i,4,2) = -ws(i)*rwork2
       arr(i,4,4) = cu(i) + cl(i) + rwork + rwork2
       arr(i,4,5) = -cl(i)
       rhs(i,4) = cu(i)*ta(i)*bps(i)                &
            + (1.0_r8-wu(i))*rwork*tu(i)          &
            + (1.0_r8-ws(i))*rwork2*ts(i)

       !
    END DO
    !
    ! lower air temperature t34
    !
    DO i = 1, npoi
       !
       rwork = xl(i)*sl(i)
       rwork2 = fi(i)*si(i)
       arr(i,5,3) = -wl(i)*rwork
       arr(i,5,4) = -cl(i)
       arr(i,5,5) = cl(i) + rwork                     &
            + (1.0_r8-fi(i))*sg(i) + rwork2
       arr(i,5,8) = -wg(i)*(1.0_r8-fi(i))*sg(i)
       arr(i,5,9) = -wi(i)*rwork2
       rhs(i,5) = (1.0_r8-wl(i))*rwork           *tl(i)    &
            + (1.0_r8-wg(i))*(1.0_r8-fi(i))*sg(i)*tg(i)   &
            + (1.0_r8-wi(i))*rwork2          *ti(i)

       !
    END DO
    !
    ! upper air specific humidity q12
    !
    DO i = 1, npoi
       !
       rwork = xu(i)*(suw(i)+sut(i))
       rwork2 = xs(i)*ssw(i)
       arr(i,6,1) = -wu(i)*rwork *dqu(i)
       arr(i,6,2) = -ws(i)*rwork2*dqs(i)
       arr(i,6,6) = cu(i) + cl(i)   &
            + rwork + rwork2
       arr(i,6,7) = -cl(i)
       rhs(i,6) = cu(i)*qa(i)  &
            + rwork  * (qu(i)-wu(i)*dqu(i)*tu(i))  &
            + rwork2 * (qs(i)-ws(i)*dqs(i)*ts(i))
       !
    END DO
    !
    ! lower air specific humidity q34
    !
    DO i = 1, npoi
       !
       rwork  = xl(i)*(slw(i)+slt(i))
       rwork2 = (1.0_r8-fi(i))*sg(i)
       arr(i,7,3) = -wl(i)*rwork*dql(i)
       arr(i,7,6) = -cl(i)
       arr(i,7,7) = cl(i) + rwork    &
            + rwork2 +fi(i)*si(i)
       arr(i,7,8) = -wg(i)*rwork2*qgfac(i)*dqg(i)
       arr(i,7,9) = -wi(i)*fi(i)*si(i)*dqi(i)
       rhs(i,7)= rwork           *(ql(i)-wl(i)*dql(i)*tl(i))  &
            + rwork2*qgfac(i) *(qg(i)-wg(i)*dqg(i)*tg(i))  &
            + fi(i) *si(i)    *(qi(i)-wi(i)*dqi(i)*ti(i))
       !
    END DO
    !
    ! soil skin temperature
    !
    ! (there is no wg in this eqn since it solves for a fully
    ! implicit tg. wg can be thought of as the fractional soil
    ! area using a fully implicit soln, and 1-wg as that using a
    ! fully explicit soln. the combined soil temperature is felt
    ! by the lower air, so wg occurs in the t34,q34 eqns above.)
    !
    DO i = 1, npoi
       !
       rwork  = sg(i)*cp(i)
       rwork2 = sg(i)*hvasug(i)
       arr(i,8,5) = -rwork
       arr(i,8,7) = -rwork2
       arr(i,8,8) = rwork + rwork2*qgfac(i)*dqg(i)   &
            + cog(i) + zirg(i)
       rhs(i,8) = -rwork2*qgfac(i)*(qg(i)-dqg(i)*tg(i))   &
            + cog(i)*tsoi(i,1)   &
            + solg(i) + firg(i) + zirg(i) * tgold(i)
       !
    END DO
    !
    ! snow skin temperature
    !
    ! (there is no wi here, for the same reason as for wg above.)
    !
    DO i = 1, npoi
       !
       rwork  = si(i)*cp(i)
       rwork2 = si(i)*hvasui(i)
       arr(i,9,5) = -rwork
       arr(i,9,7) = -rwork2
       arr(i,9,9) = rwork + rwork2*dqi(i)  &
            + coi(i) + ziri(i)
       rhs(i,9) = -rwork2*(qi(i)-dqi(i)*ti(i))   &
            + coi(i)*tsno(i,1)              &
            + soli(i) + firi(i) + ziri(i) * tiold(i)
       !
    END DO
    !
    ! solve the systems of equations
    !
    DO  i = 1, npoi
       !
       ! WRITE(*,*)'pkubota1',rhs(i,1),rhs(i,2),rhs(i,3),rhs(i,4),rhs(i,5),rhs(i,8),rhs(i,9), rhs(i,6), rhs(i,7)
       !
    END DO

    CALL linsolve (arr   , &! INTENT(INOUT)
         rhs   , &! INTENT(INOUT)
         vec   , &! INTENT(INOUT)
         mplate, &! INTENT(IN   )
         nqn   , &! INTENT(IN   )
         npoi    )! INTENT(IN   )
    DO  i = 1, npoi
       !
       !WRITE(*,*)'pkubota2',vec(i,1),vec(i,2),vec(i,3),vec(i,4),vec(i,5),vec(i,8),vec(i,9), vec(i,6), vec(i,7)
       !
    END DO

    !
    ! copy this iteration's solution to t*, q12, q34
    !
    DO  i = 1, npoi
       !
       tu(i)  = MIN(MAX(vec(i,1),180.0_r8),340.0_r8)
       ts(i)  = MIN(MAX(vec(i,2),180.0_r8),340.0_r8)
       tl(i)  = MIN(MAX(vec(i,3),180.0_r8),340.0_r8)
       t12(i) = MIN(MAX(vec(i,4),180.0_r8),340.0_r8)
       t34(i) = MIN(MAX(vec(i,5),180.0_r8),340.0_r8)
       tg(i)  = MIN(MAX(vec(i,8),180.0_r8),340.0_r8)
       ti(i)  = MIN(MAX(vec(i,9),180.0_r8),340.0_r8)
       !
       q12(i) = MAX(vec(i,6),1.0e-12_r8)
       q34(i) = MAX(vec(i,7),1.0e-12_r8)
       !
    END DO
    !
    ! all done except for final flux calculations,
    ! so loop back for the next iteration (except the last)
    !
    IF (iter.LT.niter) RETURN
    !
    ! evaluate sensible heat and moisture fluxes (per unit
    ! leaf/stem/snow-free/snow-covered area as appropriate)
    !
    ! *******************************
    ! diagnostic sensible heat fluxes
    ! *******************************
    !
    DO i = 1, npoi
       !
       fsena(i) = cp(i) * cu(i) * (ta(i)*bps(i) - t12(i))
       !
       tgav = wg(i)*tg(i) + (1.0_r8-wg(i))*tgpre(i)
       fseng(i) = cp(i) * sg(i) * (tgav - t34(i))
       !
       tiav = wi(i)*ti(i) + (1.0_r8-wi(i))*tipre(i)
       fseni(i) = cp(i) * si(i) * (tiav - t34(i))

       tuav = wu(i)*tu(i) + (1.0_r8 - wu(i))*tupre(i)
       fsenu(i) = cp(i) * su(i) * (tuav - t12(i))
       !
       tsav = ws(i)*ts(i) + (1.0_r8 - ws(i))*tspre(i)
       fsens(i) = cp(i) * ss(i) * (tsav - t12(i))
       !
       tlav = wl(i)*tl(i) + (1.0_r8 - wl(i))*tlpre(i)
       fsenl(i) = cp(i) * sl(i) * (tlav - t12(i))
       !
    END DO
    !
    ! *************************
    ! calculate moisture fluxes
    ! *************************
    !
    DO i = 1, npoi
       !
       ! total evapotranspiration from the entire column
       !
       fvapa(i)  = cu(i) * (qa(i)-q12(i))
       !
       ! evaporation from wet surfaces in the upper canopy
       ! and transpiration per unit leaf area - upper canopy
       !
       quav = qu(i) + wu(i)*dqu(i)*(tu(i)-tupre(i))

       fvapuw(i) = suw(i) * (quav-q12(i))
       fvaput(i) = MAX (0.0_r8, sut(i) * (quav-q12(i)))
       !
       ! evaporation from wet surfaces on stems
       !
       qsav = qs(i) + ws(i)*dqs(i)*(ts(i)-tspre(i))
       fvaps(i) = ssw(i) * (qsav-q12(i))
       !
       ! evaporation from wet surfaces in the lower canopy
       ! and transpiration per unit leaf area - lower canopy
       !
       !        WRITE(*,*)qlav,q34(i),tu(i) 
       qlav = ql(i) + wl(i)*dql(i)*(tl(i)-tlpre(i))
       fvaplw(i) = slw(i)  * (qlav-q34(i))
       fvaplt(i) = MAX (0.0_r8, slt0(i) * (qlav-q34(i)))
       !
       ! evaporation from the ground
       !
       qgav = qg(i) + wg(i)*dqg(i)*(tg(i)-tgpre(i))
       fvapg(i) = sg(i) * (qgfac(i)*qgav - q34(i))
       !
       ! evaporation from the snow
       !
       qiav = qi(i) + wi(i)*dqi(i)*(ti(i)-tipre(i))
       fvapi(i) = si(i) * (qiav-q34(i))
       !
    END DO
    ! 
    ! adjust ir fluxes
    !
    DO i = 1, npoi
       !
       firg(i) = firg(i) - wg(i)*zirg(i)*(tg(i) - tgold(i))
       firi(i) = firi(i) - wi(i)*ziri(i)*(ti(i) - tiold(i))
       firb(i) = firb(i) + (1.0_r8-fi(i))*wg(i)*zirg(i)*(tg(i)-tgold(i)) &
            +     fi(i) *wi(i)*ziri(i)*(ti(i)-tiold(i))
       !
       ! impose constraint on skin temperature
       !
       ti(i) = MIN (ti(i), tmelt)
       !
    END DO
    !
    ! set upsoi[u,l], the actual soil water uptake rates from each
    ! soil layer due to transpiration in the upper and lower stories,
    ! for the soil model 
    !
    DO  k = 1, nsoilay
       DO  i = 1, npoi
          !
          !
          !          wsoi(i,k) = wsoi(i,k) - dtime * (upsoiu(i,k) + upsoil(i,k)) / &
          !                                  (rhow * porosflo(i,k) * hsoi(i,k))
          !
          upsoiu(i,k) = fvaput(i) * 2.0_r8 * lai(i,2) * fu(i) *  &
               stressu(i,k) / MAX (stresstu(i), epsilon)
          !
          upsoil(i,k) = fvaplt(i) * 2.0_r8 * lai(i,1) * fl(i) *  &
               (1.0_r8 - fi(i)) *                        &
               stressl(i,k) / MAX (stresstl(i), epsilon)
          !
       END DO
    END DO
    !
    ! set net evaporation from intercepted water, net evaporation
    ! from the surface, and net transpiration rates
    !
    DO i = 1, npoi
       !
       ! evaporation from intercepted water
       !
       ginvap(i) = fvapuw(i) * 2.0_r8 * lai(i,2) * fu(i) +      &
            fvaps (i) * 2.0_r8 * sai(i,2) * fu(i) +      &
            fvaplw(i) * 2.0_r8 * (lai(i,1) + sai(i,1)) * &
            fl(i) * (1.0_r8 - fi(i))
       !
       ! evaporation from soil and snow surfaces
       !
       gsuvap(i) = fvapg(i)  * (1.0_r8 - fi(i)) + fvapi(i)  * fi(i)
       !
       ! transpiration
       !
       gtrans(i) = fvaput(i) * 2.0_r8 * lai(i,2) * fu(i) +            &
            fvaplt(i) * 2.0_r8 * lai(i,1) * fl(i) * (1.0_r8-fi(i))
       !
       gtransu(i) = fvaput(i) * 2.0_r8 * lai(i,2) * fu(i)
       gtransl(i) = fvaplt(i) * 2.0_r8 * lai(i,1) * fl(i) * (1.0_r8-fi(i))
       !
    END DO
    !
    RETURN
  END SUBROUTINE turvap

  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE fstrat (tb      , &! INTENT(IN   )
       tt      , &! INTENT(IN   )
       bps   , &! INTENT(IN   )
       qb      , &! INTENT(IN   )
       qt      , &! INTENT(IN   )
       zb      , &! INTENT(IN   )
       zt      , &! INTENT(IN   )
       albm    , &! INTENT(IN   )
       albh    , &! INTENT(IN   )
       alt     , &! INTENT(IN   )
       u       , &! INTENT(IN   )
       rich    , &! INTENT(OUT  )
       stram   , &! INTENT(OUT  )
       strah   , &! INTENT(OUT  )
       iter    , &! INTENT(IN   )
       npoi    , &! INTENT(IN   )
       grav    , &! INTENT(IN   )
       vonk      )! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    ! computes mixing-length stratification correction factors
    ! for momentum and heat/vapor, for current 1d strip, using
    ! parameterizations in louis (1979),blm,17,187. first computes
    ! richardson numbers. sets an upper limit to richardson numbers
    ! so lower-veg winds don't become vanishingly small in very
    ! stable conditions (cf, carson and richards,1978,blm,14,68)
    !
    ! system (i) is as in louis(1979). system (vi) is improved as
    ! described in louis(1982), ecmwf workshop on planetary boundary
    ! layer parameterizations,november 1981,59-79 (qc880.4 b65w619)
    !
    ! common blocks
    !
    IMPLICIT NONE
    !
    !
    ! input variables
    !
    INTEGER, INTENT(IN   ) :: npoi         ! total number of land points
    REAL(KIND=r8), INTENT(IN   ) :: grav         ! gravitational acceleration (m s-2)
    REAL(KIND=r8), INTENT(IN   ) :: vonk         ! von karman constant (dimensionless)
    INTEGER, INTENT(IN   ) :: iter         ! current iteration number
    REAL(KIND=r8), INTENT(IN   ) :: bps (npoi)        ! pot. temp factor for ttop (relative to bottom,supplied)
    REAL(KIND=r8), INTENT(IN   ) :: tb    (npoi) ! bottom temperature (supplied)
    REAL(KIND=r8), INTENT(IN   ) :: tt    (npoi) ! top temperature (supplied)
    REAL(KIND=r8), INTENT(IN   ) :: qb    (npoi) ! bottom specific humidity (supplied)
    REAL(KIND=r8), INTENT(IN   ) :: qt    (npoi) ! top specific humidity (supplied)
    REAL(KIND=r8), INTENT(IN   ) :: zb    (npoi) ! height of bottom (supplied)
    REAL(KIND=r8), INTENT(IN   ) :: zt    (npoi) ! height of top (supplied)
    REAL(KIND=r8), INTENT(IN   ) :: albm  (npoi) ! log (bottom roughness length) for momentum (supplied)
    REAL(KIND=r8), INTENT(IN   ) :: albh  (npoi) ! log (bottom roughness length) for heat/h2o (supplied)
    REAL(KIND=r8), INTENT(IN   ) :: alt   (npoi) ! log (z at top) (supplied)
    REAL(KIND=r8), INTENT(IN   ) :: u     (npoi) ! wind speed at top (supplied)
    REAL(KIND=r8), INTENT(OUT  ) :: rich  (npoi) ! richardson number (returned)
    REAL(KIND=r8), INTENT(OUT  ) :: stram (npoi) ! stratification factor for momentum (returned)
    REAL(KIND=r8), INTENT(OUT  ) :: strah (npoi) ! stratification factor for heat/vap (returned)
    !
    ! local variables
    !
    INTEGER indp(npoi)   !
    INTEGER indq(npoi)   !
    REAL(KIND=r8) stramx(npoi) !
    REAL(KIND=r8) strahx(npoi) !

    !
    INTEGER i
    INTEGER j
    INTEGER np
    INTEGER nq
    !
    REAL(KIND=r8) zht
    REAL(KIND=r8) zhb
    REAL(KIND=r8) xm
    REAL(KIND=r8) xh
    REAL(KIND=r8) rwork
    REAL(KIND=r8) ym
    REAL(KIND=r8) yh
    REAL(KIND=r8) z
    REAL(KIND=r8) w
    ! ---------------------------------------------------------------------
    np = 0
    nq = 0
    stram=0.0_r8
    strah=0.0_r8
    !
    ! do for all points
    !
    DO i = 1, npoi
       !
       ! calculate richardson numbers
       !
       zht = tt(i)*bps(i)*(1._r8+.622_r8*qt(i))
       zhb = tb(i)*      (1._r8+.622_r8*qb(i))
       !
       rich(i) = grav * MAX (zt(i)-zb(i), 0.0_r8) &
            * (zht-zhb) / (0.5_r8*(zht+zhb) * u(i)**2)
       !
       ! bound richardson number between -2.0 (unstable) to 1.0 (stable)
       !
       rich(i) = MAX (-2.0_r8, MIN (rich(i), 1.0_r8))
       !
    END DO
    !
    ! set up indices for points with negative or positive ri
    !
    DO  i = 1, npoi
       !
       IF (rich(i).LE.0.0_r8) THEN
          np = np + 1
          indp(np) = i
       ELSE
          nq = nq + 1
          indq(nq) = i
       END IF
       !
    END DO
    !
    ! calculate momentum and heat/vapor factors for negative ri
    !
    IF (np.GT.0) THEN
       !
       DO j = 1, np
          !
          i = indp(j)
          !
          xm = MAX (alt(i)-albm(i), .5_r8)
          xh = MAX (alt(i)-albh(i), .5_r8)
          !
          rwork = SQRT(-rich(i))
          !
          ym = (vonk/xm)**2 * EXP (0.5_r8*xm) * rwork
          yh = (vonk/xh)**2 * EXP (0.5_r8*xh) * rwork
          !
          ! system (vi)
          !
          stramx(i) =   1.0_r8 - 2*5*rich(i) / (1.0_r8 + 75*ym)
          strahx(i) =   1.0_r8 - 3*5*rich(i) / (1.0_r8 + 75*yh)
          !
       END DO
       !
    END IF
    !
    ! calculate momentum and heat/vapor factors for positive ri
    !
    IF (nq.GT.0) THEN
       !
       DO j=1,nq
          !
          i = indq(j)
          !
          ! system (vi)
          !
          z = SQRT(1.0_r8 + 5 * rich(i))
          !
          stramx(i) = 1.0_r8 / (1.0_r8 + 2*5*rich(i) / z)
          strahx(i) = 1.0_r8 / (1.0_r8 + 3*5*rich(i) * z)
          !
       END DO
       !
    END IF
    !
    ! except for the first iteration, weight results with the
    ! previous iteration's values. this improves convergence by
    ! avoiding flip-flop between stable/unstable stratif, eg,
    ! with cold upper air and the lower surface being heated by
    ! solar radiation
    !
    IF (iter.EQ.1) THEN
       !
       DO i = 1, npoi
          !
          stram(i) = stramx(i)
          strah(i) = strahx(i)
          !
       END DO
       !
    ELSE
       !
       w = 0.5_r8
       !
       DO i = 1, npoi
          !
          stram(i) = w * stramx(i) + (1.0_r8 - w) * stram(i)
          strah(i) = w * strahx(i) + (1.0_r8 - w) * strah(i)
          !
       END DO
       !
    END IF
    !
    RETURN
  END SUBROUTINE fstrat
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE impexp (wimp     , &! INTENT(INOUT)
       tveg     , &! INTENT(IN   )
       ch       , &! INTENT(IN   )
       wliq     , &! INTENT(IN   )
       wsno     , &! INTENT(IN   )
       iter     , &! INTENT(IN   )
       npoi     , &! INTENT(IN   )
       tmelt    , &! INTENT(IN   )
       hfus     , &! INTENT(IN   )
       epsilon    )! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    ! sets the implicit vs explicit fraction in turvap calcs for
    ! upper leaves, upper stems or lower veg. this is to account for
    ! temperatures of freezing/melting intercepted h2o constrained
    ! at the melt point. if a purely implicit calc is used for such
    ! a surface, the predicted temperature would be nearly the atmos
    ! equil temp with little sensible heat input, so the amount of
    ! freezing or melting is underestimated. however, if a purely
    ! explicit calc is used with only a small amount of intercepted
    ! h2o, the heat exchange can melt/freeze all the h2o and cause
    ! an unrealistic huge change in the veg temp. the algorithm
    ! below attempts to avoid both pitfalls
    !
    ! common blocks
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: npoi
    REAL(KIND=r8), INTENT(IN   ) :: tmelt       ! freezing point of water (K)
    REAL(KIND=r8), INTENT(IN   ) :: hfus        ! latent heat of fusion of water (J kg-1)
    REAL(KIND=r8), INTENT(IN   ) :: epsilon     ! small quantity to avoid zero-divides and other
    ! truncation or machine-limit troubles with small
    ! values. should be slightly greater than o(1)
    ! machine precision

    !
    ! input/output variables
    !
    INTEGER, INTENT(IN   ) :: iter        ! current iteration number (supplied)
    !
    REAL(KIND=r8), INTENT(INOUT) ::  wimp(npoi) ! implicit/explicit fraction (0 to 1) (returned)
    REAL(KIND=r8), INTENT(IN   ) ::  tveg(npoi) ! temperature of veg (previous iteration's soln) (supp)
    REAL(KIND=r8), INTENT(IN   ) ::  ch  (npoi) ! heat capacity of veg (supplied)
    REAL(KIND=r8), INTENT(IN   ) ::  wliq(npoi) ! veg intercepted liquid (supplied)
    REAL(KIND=r8), INTENT(IN   ) ::  wsno(npoi) ! veg intercepted snow (supplied)
    !
    ! local variables
    !
    INTEGER ::  i
    !
    REAL(KIND=r8) ::  h
    REAL(KIND=r8) ::  z
    REAL(KIND=r8) ::  winew
    !
    ! for first iteration, set wimp to fully implicit, and return
    !
    IF (iter.EQ.1) THEN
       wimp=1.0_r8
       !CALL const(wimp, npoi, 1.0)
       RETURN
    END IF
    !
    ! for second and subsequent iterations, estimate wimp based on
    ! the previous iterations's wimp and its resulting tveg.
    !
    ! calculate h, the "overshoot" heat available to melt any snow
    ! or freeze any liquid. then the explicit fraction is taken to
    ! be the ratio of h to the existing h2o's latent heat (ie, 100%
    ! explicit calculation if not all of the h2o would be melted or
    ! frozen). so winew, the implicit amount, is 1 - that ratio.
    ! but since we are using the previous iteration's t* results
    ! for the next iteration, to ensure convergence we need to damp
    ! the returned estimate wimp by averaging winew with the 
    ! previous estimate. this works reasonably well even with a
    ! small number of iterations (3), since for instance with large
    ! amounts of h2o so that wimp should be 0., a good amount of 
    ! h2o is melted or frozen with wimp = .25
    !
    DO i = 1, npoi
       !
       h = ch(i) * (tveg(i) - tmelt)
       z = MAX (ABS(h), epsilon)
       !
       winew = 1.0_r8
       !
       IF (h.GT.epsilon)  winew = 1.0_r8 - MIN (1.0_r8, hfus * wsno(i) / z)
       IF (h.LT.-epsilon) winew = 1.0_r8 - MIN (1.0_r8, hfus * wliq(i) / z)
       !
       wimp(i) = 0.5_r8 * (wimp(i) + winew)
       !
    END DO
    !
    RETURN
  END SUBROUTINE impexp
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE impexp2 (wimp    , &! INTENT(INOUT)
       t       , &! INTENT(IN   )
       told    , &! INTENT(IN   )
       iter    , &! INTENT(IN   )
       npoi    , &! INTENT(IN   )
       tmelt   , &! INTENT(IN   )
       !  hfus    , &! INTENT(IN   )
       epsilon   )! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    ! sets the implicit vs explicit fraction in turvap calcs for
    ! seaice or snow skin temperatures, to account for temperatures
    ! of freezing/melting surfaces being constrained at the melt
    ! point
    !
    ! unlike impexp, don't have to allow for all h2o 
    ! vanishing within the timestep
    !
    ! wimp   = implicit fraction (0 to 1) (returned)
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: npoi
    REAL(KIND=r8), INTENT(IN   ) :: tmelt       ! freezing point of water (K)
    !      REAL(KIND=r8), INTENT(IN   ) :: hfus        ! latent heat of fusion of water (J kg-1)
    REAL(KIND=r8), INTENT(IN   ) :: epsilon     ! small quantity to avoid zero-divides and other
    ! truncation or machine-limit troubles with small
    ! values. should be slightly greater than o(1)
    ! machine precision

    !
    ! input variables
    !
    INTEGER , INTENT(IN   ) :: iter
    REAL(KIND=r8) , INTENT(INOUT) :: wimp(npoi) 
    REAL(KIND=r8) , INTENT(IN   ) :: t   (npoi) 
    REAL(KIND=r8) , INTENT(IN   ) :: told(npoi)
    !
    ! local variables
    !
    INTEGER :: i    ! loop indice
    !
    ! for first iteration, set wimp to fully implicit, and return
    !
    IF (iter.EQ.1) THEN
       wimp=1.0_r8
       !CALL const(wimp, npoi, 1.0)
       RETURN
    END IF
    !
    DO i = 1, npoi
       !
       IF ((t(i)-told(i)).GT.epsilon) wimp(i) = (tmelt - told(i)) /  &
            (t(i)  - told(i))
       wimp(i) = MAX (0.0_r8, MIN (1.0_r8, wimp(i)))
       !
    END DO
    !
    RETURN
  END SUBROUTINE impexp2
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE fwetcal(npoi     , &! INTENT(IN   )
       fwetu    , &! INTENT(OUT  )
       rliqu    , &! INTENT(OUT  )
       fwets    , &! INTENT(OUT  )
       rliqs    , &! INTENT(OUT  )
       fwetl    , &! INTENT(OUT  )
       rliql    , &! INTENT(OUT  )
       wliqu    , &! INTENT(IN   )
       wliqumax , &! INTENT(IN   ) ::
       wsnou    , &! INTENT(IN   ) ::
       wsnoumax , &! INTENT(IN   ) ::
       tu       , &! INTENT(IN   )
       wliqs    , &! INTENT(IN   )
       wliqsmax , &! INTENT(IN   )
       wsnos    , &! INTENT(IN   )
       wsnosmax , &! INTENT(IN   )
       ts       , &! INTENT(IN   )
       wliql    , &! INTENT(IN   )
       wliqlmax , &! INTENT(IN   )
       wsnol    , &! INTENT(IN   )
       wsnolmax , &! INTENT(IN   )
       tl       , &! INTENT(IN   )
       epsilon  , &! INTENT(IN   )
       tmelt      )! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    ! calculates fwet[u,s,l], the fractional areas wetted by 
    ! intercepted h2o (liquid and snow combined) -  the maximum value
    ! fmax (<1) allows some transpiration even in soaked conditions
    !
    ! use a linear relation between fwet* and wliq*,wsno* (at least
    ! for small values), so that the implied "thickness" is constant
    ! (equal to wliq*max, wsno*max as below) and the typical amount
    ! evaporated in one timestep in steph2o will not make wliq*,wsno*
    ! negative and thus cause a spurious unrecoverable h2o loss
    !
    ! (the max(w*max,.01) below numericaly allows w*max = 0 without
    ! blowup.) in fact evaporation in one timestep *does* sometimes
    ! exceed wliq*max (currently 1 kg/m2), so there is an additional
    ! safeguard in turvap that limits the wetted-area aerodynamic
    ! coefficients suw,ssw,slw -- if that too fails, there is an 
    ! ad-hoc adjustment in steph2o2 to reset negative wliq*,wsno*
    ! amounts to zero by taking some water vapor from the atmosphere.
    !
    ! also sets rliq[u,s,l], the proportion of fwet[u,s,l] due to
    ! liquid alone. fwet,rliq are used in turvap, rliq in steph2o. 
    ! (so rliq*fwet, (1-rliq)*fwet are the fractional areas wetted
    ! by liquid and snow individually.) if fwet is 0, choose rliq
    ! = 1 if t[u,s,l] ge tmelt or 0 otherwize, for use by turvap and
    ! steph2o in case of initial dew formation on dry surface.
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: npoi         ! total number of land points
    REAL(KIND=r8), INTENT(IN   ) :: epsilon      ! small quantity to avoid zero-divides and other
    ! truncation or machine-limit troubles with small
    ! values. should be slightly greater than o(1)
    ! machine precision
    REAL(KIND=r8), INTENT(IN   ) :: tmelt        ! freezing point of water (K)
    REAL(KIND=r8), INTENT(OUT  ) :: fwetu (npoi) ! fraction of upper canopy leaf area wetted by intercepted liquid and/or snow
    REAL(KIND=r8), INTENT(OUT  ) :: rliqu (npoi) ! proportion of fwetu due to liquid
    REAL(KIND=r8), INTENT(OUT  ) :: fwets (npoi) ! fraction of upper canopy stem area wetted by intercepted liquid and/or snow
    REAL(KIND=r8), INTENT(OUT  ) :: rliqs (npoi) ! proportion of fwets due to liquid
    REAL(KIND=r8), INTENT(OUT  ) :: fwetl (npoi) ! fraction of lower canopy stem & leaf area wetted by intercepted liquid and/or snow
    REAL(KIND=r8), INTENT(OUT  ) :: rliql (npoi) ! proportion of fwetl due to liquid
    REAL(KIND=r8), INTENT(IN   ) :: wliqu (npoi) ! intercepted liquid h2o on upper canopy leaf area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wliqumax     ! maximum intercepted water on a unit upper canopy leaf area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wsnou (npoi) ! intercepted frozen h2o (snow) on upper canopy leaf area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wsnoumax     ! intercepted snow capacity for upper canopy leaves (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: tu    (npoi) ! temperature of upper canopy leaves (K)
    REAL(KIND=r8), INTENT(IN   ) :: wliqs (npoi) ! intercepted liquid h2o on upper canopy stem area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wliqsmax     ! maximum intercepted water on a unit upper canopy stem area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wsnos (npoi) ! intercepted frozen h2o (snow) on upper canopy stem area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wsnosmax     ! intercepted snow capacity for upper canopy stems (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: ts    (npoi) ! temperature of upper canopy stems (K)
    REAL(KIND=r8), INTENT(IN   ) :: wliql (npoi) ! intercepted liquid h2o on lower canopy leaf and stem area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wliqlmax            ! maximum intercepted water on a unit lower canopy stem & leaf area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wsnol (npoi) ! intercepted frozen h2o (snow) on lower canopy leaf & stem area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wsnolmax     ! intercepted snow capacity for lower canopy leaves & stems (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: tl    (npoi) ! temperature of lower canopy leaves & stems(K)

    !INCLUDE 'com1d.h'
    !
    ! local variables
    !
    INTEGER :: i           ! loop indice
    !
    REAL(KIND=r8) ,PARAMETER :: fmax = 0.25_r8      ! maximum water cover on two-sided leaf
    REAL(KIND=r8) :: xliq        ! fraction of wetted leaf (liquid only)
    REAL(KIND=r8) :: xtot        ! fraction of wetted leaf (liquid and snow)
    !
    ! upper leaves
    !
    DO i = 1, npoi
       !
       xliq = wliqu(i) / MAX (wliqumax, 0.01_r8)
       xtot = xliq + wsnou(i) / MAX (wsnoumax, 0.01_r8)
       !
       fwetu(i) = MIN (fmax, xtot)
       rliqu(i) = xliq / MAX (xtot, epsilon)
       !
       IF (fwetu(i).EQ.0.0_r8) THEN
          rliqu(i) = 1.0_r8
          IF (tu(i).LT.tmelt) rliqu(i) = 0.0_r8
       END IF
       !
    END DO
    !
    ! upper stems
    !
    DO  i = 1, npoi
       !
       xliq = wliqs(i) / MAX (wliqsmax, 0.01_r8)
       xtot = xliq + wsnos(i) / MAX (wsnosmax, 0.01_r8)
       !
       fwets(i) = MIN (fmax, xtot)
       rliqs(i) = xliq / MAX (xtot, epsilon)
       !
       IF (fwets(i).EQ.0.0_r8) THEN
          rliqs(i) = 1.0_r8
          IF (ts(i).LT.tmelt) rliqs(i) = 0.0_r8
       END IF
       !
    END DO
    !
    ! lower veg
    !
    DO  i = 1, npoi
       !
       xliq = wliql(i) / MAX (wliqlmax, 0.01_r8)
       xtot = xliq + wsnol(i) / MAX (wsnolmax, 0.01_r8)
       !
       fwetl(i) = MIN (fmax, xtot)
       rliql(i) = xliq / MAX (xtot, epsilon)
       !
       IF (fwetl(i).EQ.0.0_r8) THEN
          rliql(i) = 1.0_r8
          IF (tl(i).LT.tmelt) rliql(i) = 0.0_r8
       END IF
       !
    END DO
    !
    RETURN
  END SUBROUTINE fwetcal
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE cascade(npoi    , &! INTENT(IN   )
       epsilon , &! INTENT(IN   )
       dtime   , &! INTENT(IN   )
       ch2o    , &! INTENT(IN   )
       cice    , &! INTENT(IN   )  
       tmelt   , &! INTENT(IN   )
       hfus    , &! INTENT(IN   )
       vzero   , &! INTENT(IN   )
       snowg   , &! INTENT(OUT  )
       tsnowg  , &! INTENT(OUT  )
       tsnowl  , &! INTENT(INOUT)
       pfluxl  , &! INTENT(OUT  )
       raing   , &! INTENT(OUT  )
       traing  , &! INTENT(OUT  )
       trainl  , &! INTENT(INOUT)
       snowl   , &! INTENT(OUT  )
       tsnowu  , &! INTENT(INOUT)
       pfluxu  , &! INTENT(OUT  )
       rainu   , &! INTENT(INOUT)
       trainu  , &! INTENT(INOUT)
       snowu   , &! INTENT(INOUT)
       pfluxs  , &! INTENT(OUT  )
       rainl   , &! INTENT(OUT  )
       wliqmin , &! INTENT(INOUT)
       wliqumax, &! INTENT(IN   )  
       wsnomin , &! INTENT(INOUT)
       wsnoumax, &! INTENT(IN   )  
       t12     , &! INTENT(IN   )  
       lai     , &! INTENT(IN   )  
       tu      , &! INTENT(IN   ) 
       wliqu   , &! INTENT(INOUT)
       wsnou   , &! INTENT(INOUT)
       tdripu  , &! INTENT(IN   ) 
       tblowu  , &! INTENT(IN   ) 
       sai     , &! INTENT(IN   ) 
       ts      , &! INTENT(IN   ) 
       wliqs   , &! INTENT(INOUT) :: 
       wsnos   , &! INTENT(INOUT) :: 
       tdrips  , &! INTENT(IN   ) 
       tblows  , &! INTENT(IN   ) 
       wliqsmax, &! INTENT(IN   ) 
       wsnosmax, &! INTENT(IN   ) 
       fu      , &! INTENT(IN   ) 
       t34     , &! INTENT(IN   ) 
       tl      , &! INTENT(IN   ) 
       wliql   , &! INTENT(INOUT) :: 
       wsnol   , &! INTENT(INOUT) :: 
       tdripl  , &! INTENT(IN   ) 
       tblowl  , &! INTENT(IN   ) 
       wliqlmax, &! INTENT(IN   ) 
       wsnolmax, &! INTENT(IN   ) 
       fl      , &! INTENT(IN   ) 
       raina   , &! INTENT(IN   ) 
       ta      , &! INTENT(IN   ) 
       qa      , &! INTENT(IN   ) 
       psurf   , &! INTENT(IN   )
       snowa   )  ! INTENT(IN   )




    ! ---------------------------------------------------------------------
    !
    ! steps intercepted h2o due to drip, precip, and min/max limits
    !
    ! calls steph2o for upper leaves, upper stems and lower veg in
    ! iurn, adjusting precips at each level
    !
    IMPLICIT NONE
    !


    INTEGER, INTENT(IN   ) :: npoi    ! total number of land points

    REAL(KIND=r8), INTENT(IN   ) :: epsilon ! small quantity to avoid zero-divides and other
    ! truncation or machine-limit troubles with small
    ! values. should be slightly greater than o(1)
    ! machine precision
    REAL(KIND=r8), INTENT(IN   ) :: dtime   ! model timestep (seconds)
    REAL(KIND=r8), INTENT(IN   ) :: ch2o    ! specific heat of liquid water (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN   ) :: cice    ! specific heat of ice (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN   ) :: tmelt   ! freezing point of water (K)
    REAL(KIND=r8), INTENT(IN   ) :: hfus    ! latent heat of fusion of water (J kg-1)
    REAL(KIND=r8), INTENT(IN   ) :: vzero (npoi)! a real array of zeros, of length npoi
    REAL(KIND=r8), INTENT(OUT  ) :: snowg (npoi)! snowfall rate at soil level (kg h2o m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: tsnowg(npoi)! snowfall temperature at soil level (K) 
    REAL(KIND=r8), INTENT(INOUT) :: tsnowl(npoi)! snowfall temperature below upper canopy (K)
    REAL(KIND=r8), INTENT(OUT  ) :: pfluxl(npoi)! heat flux on lower canopy leaves & stems due to intercepted h2o (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: raing (npoi)! rainfall rate at soil level (kg m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: traing(npoi)! rainfall temperature at soil level (K)
    REAL(KIND=r8), INTENT(INOUT) :: trainl(npoi)! rainfall temperature below upper canopy (K)
    REAL(KIND=r8), INTENT(OUT  ) :: snowl (npoi)! snowfall rate below upper canopy (kg h2o m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: tsnowu(npoi)! snowfall temperature above upper canopy (K)
    REAL(KIND=r8), INTENT(OUT  ) :: pfluxu(npoi)! heat flux on upper canopy leaves due to intercepted h2o (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: rainu (npoi)! rainfall rate above upper canopy (kg m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: trainu(npoi)! rainfall temperature above upper canopy (K)
    REAL(KIND=r8), INTENT(INOUT) :: snowu (npoi)! snowfall rate above upper canopy (kg h2o m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: pfluxs(npoi)! heat flux on upper canopy stems due to intercepted h2o (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: rainl (npoi)! rainfall rate below upper canopy (kg m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: wliqmin     ! minimum intercepted water on unit vegetated area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wliqumax    ! maximum intercepted water on a unit upper canopy leaf area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: wsnomin     ! minimum intercepted snow on unit vegetated area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wsnoumax    ! intercepted snow capacity for upper canopy leaves (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: t12(npoi)   ! air temperature at z12 (K)
    REAL(KIND=r8), INTENT(IN   ) :: lai(npoi,2) ! canopy single-sided leaf area index (area leaf/area veg)
    REAL(KIND=r8), INTENT(IN   ) :: tu(npoi)    ! temperature of upper canopy leaves (K)
    REAL(KIND=r8), INTENT(INOUT) :: wliqu(npoi) ! intercepted liquid h2o on upper canopy leaf area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: wsnou(npoi) ! intercepted frozen h2o (snow) on upper canopy leaf area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: tdripu      ! decay time for dripoff of liquid intercepted by upper canopy leaves (sec)
    REAL(KIND=r8), INTENT(IN   ) :: tblowu      ! decay time for blowoff of snow intercepted by upper canopy leaves (sec)
    REAL(KIND=r8), INTENT(IN   ) :: sai(npoi,2) ! current single-sided stem area index
    REAL(KIND=r8), INTENT(IN   ) :: ts(npoi)    ! temperature of upper canopy stems (K)
    REAL(KIND=r8), INTENT(INOUT) :: wliqs(npoi) ! intercepted liquid h2o on upper canopy stem area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: wsnos(npoi) ! intercepted frozen h2o (snow) on upper canopy stem area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: tdrips           ! decay time for dripoff of liquid intercepted by upper canopy stems (sec) 
    REAL(KIND=r8), INTENT(IN   ) :: tblows          ! decay time for blowoff of snow intercepted by upper canopy stems (sec)
    REAL(KIND=r8), INTENT(IN   ) :: wliqsmax           ! maximum intercepted water on a unit upper canopy stem area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wsnosmax    ! intercepted snow capacity for upper canopy stems (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: fu(npoi)           ! fraction of overall area covered by upper canopy
    REAL(KIND=r8), INTENT(IN   ) :: t34(npoi)          ! air temperature at z34 (K)
    REAL(KIND=r8), INTENT(IN   ) :: tl(npoi)           ! temperature of lower canopy leaves & stems(K)
    REAL(KIND=r8), INTENT(INOUT) :: wliql(npoi) ! intercepted liquid h2o on lower canopy leaf and stem area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: wsnol(npoi) ! intercepted frozen h2o (snow) on lower canopy leaf & stem area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: tdripl          ! decay time for dripoff of liquid intercepted by lower canopy leaves & stem (sec)
    REAL(KIND=r8), INTENT(IN   ) :: tblowl          ! decay time for blowoff of snow intercepted by lower canopy leaves & stems (sec)
    REAL(KIND=r8), INTENT(IN   ) :: wliqlmax           ! maximum intercepted water on a unit lower canopy stem & leaf area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wsnolmax           ! intercepted snow capacity for lower canopy leaves & stems (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: fl(npoi)           ! fraction of snow-free area covered by lower  canopy
    REAL(KIND=r8), INTENT(IN   ) :: raina(npoi) ! rainfall rate (mm/s or kg m-2 s-1)
    REAL(KIND=r8), INTENT(IN   ) :: ta(npoi)    ! air temperature (K)
    REAL(KIND=r8), INTENT(IN   ) :: qa(npoi)    ! specific humidity (kg_h2o/kg_air)
    REAL(KIND=r8), INTENT(IN   ) :: psurf(npoi) ! surface pressure (Pa)
    REAL(KIND=r8), INTENT(IN   ) :: snowa(npoi) ! snowfall rate (mm/s or kg m-2 s-1 of water)

    !
    ! local variables
    !
    INTEGER i           ! loop indice
    !
    !      REAL(KIND=r8) twet3       ! Function: wet bulb temperature (K)
    REAL(KIND=r8) twetbulb  ,es  ! wet bulb temperature (K)
    !    
    REAL(KIND=r8) xai(npoi)   ! lai and/or sai for veg component
    ! (allows steph2o to work on any veg component)
    REAL(KIND=r8) rain (npoi) ! rainfall at appropriate level (modified by steph2o)
    REAL(KIND=r8) train(npoi) ! temperature of rain (modified by steph2o)  
    REAL(KIND=r8) snow (npoi) ! snowfall at appropriate level (modified by steph2o)
    REAL(KIND=r8) tsnow(npoi) ! temperature of snow (modified by steph2o)
    REAL(KIND=r8) x1   (npoi) ! 
    REAL(KIND=r8) x2   (npoi) ! 
    REAL(KIND=r8) x3   (npoi) ! 
    REAL(KIND=r8) x4   (npoi) ! 
    REAL(KIND=r8) :: tsnowu2(npoi)! snowfall temperature above upper canopy (K)
    REAL(KIND=r8) :: snowu2 (npoi)! snowfall rate above upper canopy (kg h2o m-2 s-1)
    REAL(KIND=r8) :: rainu2 (npoi)! rainfall rate above upper canopy (kg m-2 s-1)
    REAL(KIND=r8) :: trainu2 (npoi)! rainfall temperature above upper canopy (K)
    !
    ! adjust rainfall and snowfall rates at above-tree level
    !
    ! set wliqmin, wsnomin -- unlike wliq*max, wsno*max, these are
    ! part of the lsx numerical method and not from the vegetation
    ! database, and they are the same for all veg components
    !
    ! the value 0.0010 should be small compared to typical precip rates
    ! times dtime to allow any intercepted h2o to be initiated, but
    ! not too small to allow evap rates to reduce wliq*, wsno* to
    ! that value in a reasonable number of time steps
    !
    wliqmin = 0.0010_r8 * (dtime/3600.0_r8) * (wliqumax / 0.2_r8)
    wsnomin = 0.0010_r8 * (dtime/3600.0_r8) * (wsnoumax / 2.0_r8)
    !
    DO i=1,npoi
       rainu(i) = raina(i)
       !
       ! set rain temperature to the wet bulb temperature
       !
       IF (ta(i) .GT. tmelt) THEN
          !twetbulb = twet3( ta(i), qa(i), psurf(i) )
          es = fpvs2es5(ta(i))
          twetbulb = ftdp(es)
       ELSE
          twetbulb = tmelt
       END IF
       trainu(i) = MAX (twetbulb, tmelt)
       x1(i) = 0.00_r8
       x2(i) = MAX (t12(i), tmelt)
    END DO
    !
    ! calorimetrically mixes masses x1,x2,x3 with temperatures
    ! t1,t2,t3 into combined mass xm with temperature tm
    !
    !
    rainu2=rainu
    trainu2=trainu
    CALL mix (&
         rainu    , & ! INTENT(OUT  )
         trainu   , & ! INTENT(OUT  )
         rainu2   , & ! INTENT(IN   )
         trainu2  , & ! INTENT(IN   )
         x1       , & ! INTENT(IN   )
         x2       , & ! INTENT(IN   )
         vzero    , & ! INTENT(IN   )
         vzero    , & ! INTENT(IN   )
         npoi     , & ! INTENT(IN   )
         epsilon    ) ! INTENT(IN   )
    !
    DO i=1,npoi
       snowu(i) = snowa(i)
       tsnowu(i) = MIN (ta(i), tmelt)
       x1(i) = 0.00_r8
       x2(i) = MIN (t12(i), tmelt)
    END DO
    !
    !
    ! calorimetrically mixes masses x1,x2,x3 with temperatures
    ! t1,t2,t3 into combined mass xm with temperature tm
    !
    tsnowu2=tsnowu
    snowu2=snowu
    CALL mix (&
         snowu   , & ! INTENT(OUT  )
         tsnowu  , & ! INTENT(OUT  )
         snowu2  , & ! INTENT(IN   )
         tsnowu2 , & ! INTENT(IN   )
         x1      , & ! INTENT(IN   )
         x2      , & ! INTENT(IN   )
         vzero   , & ! INTENT(IN   )
         vzero   , & ! INTENT(IN   )
         npoi    , & ! INTENT(IN   )
         epsilon   ) ! INTENT(IN   )

    !
    ! set up for upper leaves
    !
    DO i = 1, npoi
       xai(i)   = 2.00_r8 * lai(i,2)
       rain(i)  = rainu(i)
       train(i) = trainu(i)
       snow(i)  = snowu(i)
       tsnow(i) = tsnowu(i)
    END DO
    !
    ! step upper leaves
    !
    CALL steph2o   &
         (tu         , &! INTENT(IN   )
         wliqu      , &! INTENT(INOUT   )
         wsnou      , &! INTENT(INOUT   )
         xai        , &! INTENT(IN   )
         pfluxu     , &! INTENT(OUT  )
         rain       , &! INTENT(INOUT   )
         train      , &! INTENT(INOUT   )
         snow       , &! INTENT(INOUT   )
         tsnow      , &! INTENT(INOUT   )
         tdripu     , &! INTENT(IN   )
         tblowu     , &! INTENT(IN   )
         wliqumax   , &! INTENT(IN   )
         wsnoumax   , &! INTENT(IN   )
         wliqmin    , &! INTENT(IN   )
         wsnomin    , &! INTENT(IN   )
         npoi       , &! INTENT(IN   )
         epsilon    , &! INTENT(IN   )
         dtime      , &! INTENT(IN   )
         ch2o       , &! INTENT(IN   )
         cice       , &! INTENT(IN   )
         tmelt      , &! INTENT(IN   )
         hfus       , &! INTENT(IN   )
         vzero        )! INTENT(IN   )
    !
    ! set up for upper stems
    ! the upper stems get precip as modified by the upper leaves
    !
    DO  i=1,npoi
       xai(i) = 2.00_r8 * sai(i,2)
    END DO
    !
    ! step upper stems
    !
    !        WRITE(*,*)xai(npoi),wsnos(npoi),tblows,i

    CALL steph2o &
         (ts       , & ! INTENT(IN   )
         wliqs    , & ! INTENT(INOUT   )
         wsnos    , & ! INTENT(INOUT   )
         xai      , & ! INTENT(IN   )
         pfluxs   , & ! INTENT(OUT  )
         rain     , & ! INTENT(INOUT   )
         train    , & ! INTENT(INOUT   )
         snow     , & ! INTENT(INOUT   )
         tsnow    , & ! INTENT(INOUT   )
         tdrips   , & ! INTENT(IN   )
         tblows   , & ! INTENT(IN   )
         wliqsmax , & ! INTENT(IN   )
         wsnosmax , & ! INTENT(IN   )
         wliqmin  , & ! INTENT(IN   )
         wsnomin  , & ! INTENT(IN   )
         npoi     , & ! INTENT(IN   )
         epsilon  , & ! INTENT(IN   )
         dtime    , & ! INTENT(IN   )
         ch2o     , & ! INTENT(IN   )
         cice     , & ! INTENT(IN   )
         tmelt    , & ! INTENT(IN   )
         hfus     , & ! INTENT(IN   )
         vzero)             ! INTENT(IN   )
    !
    ! adjust rainfall and snowfall rates at below-tree level
    ! allowing for upper-veg interception/drip/belowoff
    !
    DO i=1,npoi
       x1(i) = fu(i)*rain(i)
       x2(i) = (1.0_r8-fu(i))*rainu(i)
       x3(i) = 0.00_r8
       x4(i) = MAX (t34(i), tmelt)
    END DO
    !
    CALL mix (&
         rainl   , &! INTENT(OUT  )
         trainl  , &! INTENT(OUT  )
         x1      , &! INTENT(IN   )
         train   , &! INTENT(IN   )
         x2      , &! INTENT(IN   )
         trainu  , &! INTENT(IN   )
         x3      , &! INTENT(IN   )
         x4      , &! INTENT(IN   )
         npoi    , &! INTENT(IN   )
         epsilon   )! INTENT(IN   )

    !
    DO  i=1,npoi
       x1(i) = fu(i)*snow(i)
       x2(i) = (1.0_r8-fu(i))*snowu(i)
       x3(i) = 0.00_r8
       x4(i) = MIN (t34(i), tmelt)
    END DO
    !
    CALL mix (&
         snowl   , &! INTENT(OUT  )
         tsnowl  , &! INTENT(OUT  )
         x1      , &! INTENT(IN   )
         tsnow   , &! INTENT(IN   )
         x2      , &! INTENT(IN   )
         tsnowu  , &! INTENT(IN   )
         x3      , &! INTENT(IN   )
         x4      , &! INTENT(IN   )
         npoi    , &! INTENT(IN   )
         epsilon   )! INTENT(IN   )

    !
    ! set up for lower veg
    !
    DO i = 1, npoi
       xai(i)   = 2.00_r8 * (lai(i,1) + sai(i,1))
       rain(i)  = rainl(i)
       train(i) = trainl(i)
       snow(i)  = snowl(i)
       tsnow(i) = tsnowl(i)
    END DO
    !
    ! step lower veg
    !
    !        WRITE(*,*)xai(npoi),wsnos(npoi),tblowl,i

    CALL steph2o  &
         (tl        , &! INTENT(IN   )
         wliql     , &! INTENT(INOUT   )
         wsnol     , &! INTENT(INOUT   )
         xai       , &! INTENT(IN   )
         pfluxl    , &! INTENT(OUT  )
         rain      , &! INTENT(INOUT   )
         train     , &! INTENT(INOUT   )
         snow      , &! INTENT(INOUT   )
         tsnow     , &! INTENT(INOUT   )
         tdripl    , &! INTENT(IN   )
         tblowl    , &! INTENT(IN   )
         wliqlmax  , &! INTENT(IN   )
         wsnolmax  , &! INTENT(IN   )
         wliqmin   , &! INTENT(IN   )
         wsnomin   , &! INTENT(IN   )
         npoi      , &! INTENT(IN   )
         epsilon   , &! INTENT(IN   )
         dtime     , &! INTENT(IN   )
         ch2o      , &! INTENT(IN   )
         cice      , &! INTENT(IN   )
         tmelt     , &! INTENT(IN   )
         hfus      , &! INTENT(IN   )
         vzero       )! INTENT(IN   )

    !
    ! adjust rainfall and  snowfall rates at soil level,
    ! allowing for lower-veg interception/drip/blowoff
    !
    DO i=1,npoi
       x1(i) = fl(i) * rain(i)
       x2(i) = (1.0_r8-fl(i)) * rainl(i)
    END DO
    !
    CALL mix (&
         raing   , &! INTENT(OUT  )
         traing  , &! INTENT(OUT  )
         x1      , &! INTENT(IN   )
         train   , &! INTENT(IN   )
         x2      , &! INTENT(IN   )
         trainl  , &! INTENT(IN   )
         vzero   , &! INTENT(IN   )
         vzero   , &! INTENT(IN   )
         npoi    , &! INTENT(IN   )
         epsilon   )! INTENT(IN   )

    !
    DO i=1,npoi
       x1(i) = fl(i) * snow(i)
       x2(i) = (1.0_r8-fl(i)) * snowl(i)
    END DO
    !
    CALL mix (&
         snowg    , &! INTENT(OUT  )
         tsnowg   , &! INTENT(OUT  )
         x1       , &! INTENT(IN   )
         tsnow    , &! INTENT(IN   )
         x2       , &! INTENT(IN   )
         tsnowl   , &! INTENT(IN   )
         vzero    , &! INTENT(IN   )
         vzero    , &! INTENT(IN   )
         npoi     , &! INTENT(IN   )
         epsilon    )! INTENT(IN   )

    !
    RETURN
  END SUBROUTINE cascade
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE steph2o( &
       tveg       , &! INTENT(IN   )
       wliq       , &! INTENT(INOUT   )
       wsno       , &! INTENT(INOUT   )
       xai        , &! INTENT(IN   )
       pflux      , &! INTENT(OUT  )
       rain       , &! INTENT(INOUT   )
       train      , &! INTENT(INOUT   )
       snow       , &! INTENT(INOUT   )
       tsnow      , &! INTENT(INOUT   )
       tdrip      , &! INTENT(IN   )
       tblow      , &! INTENT(IN   )
       wliqmax    , &! INTENT(IN   )
       wsnomax    , &! INTENT(IN   )
       wliqmin    , &! INTENT(IN   )
       wsnomin    , &! INTENT(IN   )
       npoi       , &! INTENT(IN   )
       epsilon    , &! INTENT(IN   )
       dtime      , &! INTENT(IN   )
       ch2o       , &! INTENT(IN   )
       cice       , &! INTENT(IN   )
       tmelt      , &! INTENT(IN   )
       hfus       , &! INTENT(IN   )
       vzero        )! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    ! steps intercepted h2o for one canopy component (upper leaves, 
    ! upper stems, or lower veg) through one lsx time step, adjusting
    ! for h2o sensible heat and phase changes. also modifies precip
    ! due to interception and drip,blowoff
    !
    ! 
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN   ) ::  npoi     ! total number of land points
    REAL(KIND=r8), INTENT(IN   ) ::  epsilon  ! small quantity to avoid zero-divides and other
    ! truncation or machine-limit troubles with small
    ! values. should be slightly greater than o(1)
    ! machine precision
    REAL(KIND=r8), INTENT(IN   ) ::  dtime    ! model timestep (seconds)
    REAL(KIND=r8), INTENT(IN   ) ::  ch2o     ! specific heat of liquid water (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN   ) ::  cice     ! specific heat of ice (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN   ) ::  tmelt    ! freezing point of water (K)
    REAL(KIND=r8), INTENT(IN   ) ::  hfus     ! latent heat of fusion of water (J kg-1)
    REAL(KIND=r8), INTENT(IN   ) ::  vzero(npoi) ! a real array of zeros, of length npoi

    !
    ! Arguments (all arguments are supplied (unchanged) unless otherwise noted
    !    
    REAL(KIND=r8) , INTENT(IN   ) :: tdrip       ! e-folding time of liquid drip  tdrip[u,s,l]
    REAL(KIND=r8) , INTENT(IN   ) :: tblow          ! e-folding time of snow blowoff tblow[u,s,l]
    REAL(KIND=r8) , INTENT(IN   ) :: wliqmax     ! max amount of intercepted liquid wliq[u,s,l]max
    REAL(KIND=r8) , INTENT(IN   ) :: wsnomax     ! max amount of intercepted snow   wsno[u,s,l]max
    REAL(KIND=r8) , INTENT(IN   ) :: wliqmin     ! min amount of intercepted liquid (same name for u,s,l)
    REAL(KIND=r8) , INTENT(IN   ) :: wsnomin     ! min amount of intercepted snow (same name for u,s,l)
    !
    REAL(KIND=r8) , INTENT(IN   ) :: tveg (npoi) ! temperature of veg component t[u,s,l]
    REAL(KIND=r8) , INTENT(INOUT) :: wliq (npoi) ! intercepted liquid amount wliq[u,s,l] (returned)
    REAL(KIND=r8) , INTENT(INOUT) :: wsno (npoi) ! intercepted snow amount wsno[u,s,l] (returned)
    REAL(KIND=r8) , INTENT(IN   ) :: xai  (npoi) ! lai, sai, lai+sai for upper leaves/stems,lower veg
    REAL(KIND=r8) , INTENT(OUT  ) :: pflux(npoi) ! ht flux due to adjust of intercep precip (returned)
    REAL(KIND=r8) , INTENT(INOUT) :: rain (npoi) ! rainfall rate. Input: above veg, Output: below veg
    REAL(KIND=r8) , INTENT(INOUT) :: train(npoi) ! temperature of rain. (returned)
    REAL(KIND=r8) , INTENT(INOUT) :: snow (npoi) ! snowfall rate. Input: above veg, output: below veg
    REAL(KIND=r8) , INTENT(INOUT) :: tsnow(npoi) ! temperature of snow (returned)
    !
    ! local variables:
    !
    INTEGER :: i       ! loop indice
    !
    REAL(KIND=r8) :: rwork   ! 1/dtime
    REAL(KIND=r8) :: x       ! work variable
    REAL(KIND=r8) :: rwork2  ! work variable: ch2o - cice
    REAL(KIND=r8) :: dw      ! correction: freezing liguid or melting snow
    !
    REAL(KIND=r8) :: fint(npoi)  ! precip fraction intercepted by unit leaf/stem area
    REAL(KIND=r8) :: drip(npoi)  ! rate of liquid drip
    REAL(KIND=r8) :: blow(npoi)  ! rate of snow blowoff
    REAL(KIND=r8) :: rain2 (npoi) ! rainfall rate. Input: above veg, Output: below veg
    REAL(KIND=r8) :: train2 (npoi)
    REAL(KIND=r8) :: snow2 (npoi) ! snowfall rate. Input: above veg, output: below veg
    REAL(KIND=r8) :: tsnow2(npoi) ! temperature of snow (returned)

    !
    ! ---------------------------------------------------------------------
    !
    ! calculate fint, the intercepted precip fraction per unit
    ! leaf/stem area -- note 0.5 * lai or sai (similar to irrad)
    ! 
    DO i = 1, npoi
       !
       IF (xai(i).GE.epsilon) THEN
          fint(i) = ( 1.0_r8-EXP(-0.50_r8*xai(i)) )/ xai(i)
       ELSE
          fint(i) = 0.50_r8
       END IF
       !
    END DO
    !
    ! step intercepted liquid and snow amounts due to drip/blow,
    ! intercepted rainfall/snowfall, and min/max limits. also 
    ! adjust temperature of intercepted precip to current veg temp,
    ! storing the heat needed to do this in pflux for use in turvap
    ! 
    ! without these pfluxes, the implicit turvap calcs could not
    ! account for the heat flux associated with precip adjustments,
    ! especially changes of phase (see below), and so could not
    ! handle equilibrium situations such as intercepted snowfall
    ! being continuously melted by warm atmos fluxes, with the veg 
    ! temp somewhat lower than the equil atmos temp to supply heat
    ! that melts the incoming snow; (turvap would just change veg 
    ! temp to atmos equil, with little sensible heat storage...then
    ! final phase adjustment would return veg temp to melt point)
    !
    ! the use of the current (ie, previous timestep's) veg temp 
    ! gives the best estimate of what this timestep's final temp
    ! will be, at least for steady conditions
    !
    rwork = 1.0_r8 / dtime
    !
    DO i=1,npoi
       !    
       ! liquid
       !
       drip(i) = xai(i)*wliq(i)/tdrip
      !wliq(i) = wliq(i) * (1.0_r8-dtime/tdrip)
       wliq(i) = wliq(i) * (1.0_r8-dtime/MAX(tdrip,dtime))
       !
       wliq(i) = wliq(i) + dtime*rain(i)*fint(i)
       pflux(i) = rain(i)*fint(i) * (tveg(i)-train(i))*ch2o
       rain(i) = rain(i)*(1.0_r8-xai(i)*fint(i))
       !
       x = wliq(i)
       wliq(i) = MIN (wliq(i), wliqmax)
       IF (wliq(i).LT.wliqmin) wliq(i) = 0.0_r8
       drip(i) = drip(i) + xai(i)*(x-wliq(i))*rwork
       !
       ! snow
       !
       blow(i) = xai(i)*wsno(i)/tblow
       wsno(i) = wsno(i) * (1.0_r8-dtime/tblow)
       !
       wsno(i) = wsno(i) + dtime*snow(i)*fint(i)
       pflux(i) = pflux(i) + snow(i)*fint(i) * (tveg(i)-tsnow(i))*cice
       snow(i) = snow(i)*(1.0_r8-xai(i)*fint(i))
       !
       x = wsno(i)
       wsno(i) = MIN (wsno(i), wsnomax)
       IF (wsno(i).LT.wsnomin) wsno(i) = 0.0_r8 
       blow(i) = blow(i) + xai(i)*(x-wsno(i))*rwork
       !
    END DO
    !
    ! change phase of liquid/snow below/above melt point, and add
    ! required heat to pflux (see comments above). this will only
    ! affect the precip intercepted in this timestep, since original
    ! wliq, wsno must have been ge/le melt point (ensured in later
    ! call to cascad2/steph2o2)
    !
    rwork2 = ch2o - cice
    !
    DO i=1,npoi
       !
       ! liquid below freezing
       !
       dw = 0.0_r8
       IF (tveg(i).LT.tmelt)  dw = wliq(i)
       !
       pflux(i) = pflux(i)  &
            + dw * (rwork2*(tmelt-tveg(i)) - hfus) * rwork
       wliq(i) = wliq(i) - dw
       wsno(i) = wsno(i) + dw
       !
       ! snow above freezing
       !
       dw = 0.0_r8
       IF (tveg(i).GT.tmelt)  dw = wsno(i)
       !
       pflux(i) = pflux(i)   &
            + dw * (rwork2*(tveg(i)-tmelt) + hfus) * rwork
       wsno(i) = wsno(i) - dw
       wliq(i) = wliq(i) + dw
       !
    END DO
    !
    ! adjust rainfall, snowfall below veg for interception 
    ! and drip, blowoff
    !
    rain2 = rain! rainfall rate. Input: above veg, Output: below veg
    train2=train
    CALL mix (&
         rain    , & ! INTENT(OUT  )
         train   , & ! INTENT(OUT  )
         rain2   , & ! INTENT(IN   )
         train2  , & ! INTENT(IN   )
         drip    , & ! INTENT(IN   )
         tveg    , & ! INTENT(IN   )
         vzero   , & ! INTENT(IN   )
         vzero   , & ! INTENT(IN   )
         npoi    , & ! INTENT(IN   )
         epsilon   ) ! INTENT(IN   )
 
    snow2  =snow
    tsnow2 =tsnow
 
    CALL mix (&
         snow    , & ! INTENT(OUT  )
         tsnow   , & ! INTENT(OUT  )
         snow2   , & ! INTENT(IN   )
         tsnow2  , & ! INTENT(IN   )
         blow    , & ! INTENT(IN   )
         tveg    , & ! INTENT(IN   )
         vzero   , & ! INTENT(IN   )
         vzero   , & ! INTENT(IN   )
         npoi    , & ! INTENT(IN   )
         epsilon   ) ! INTENT(IN   )

    !
    RETURN
  END SUBROUTINE steph2o
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE cascad2(rliqu , &! INTENT(IN   )
       fvapuw, &! INTENT(IN   )
       fvapa , &! INTENT(INOUT)
       fsena , &! INTENT(INOUT)
       rliqs , &! INTENT(IN   )
       fvaps , &! INTENT(IN   )
       rliql , &! INTENT(IN   )
       fvaplw, &! INTENT(IN   )
       ta    , &! INTENT(IN   )
       fu    , &! INTENT(IN   )
       lai   , &! INTENT(IN   )
       tu    , &! INTENT(INOUT)
       wliqu , &! INTENT(INOUT)
       wsnou , &! INTENT(INOUT)
       chu   , &! INTENT(IN   )
       sai   , &! INTENT(IN   )
       ts    , &! INTENT(INOUT)
       wliqs , &! INTENT(INOUT)
       wsnos , &! INTENT(INOUT)
       chs   , &! INTENT(IN   )
       fl    , &! INTENT(IN   )
       tl    , &! INTENT(INOUT)
       wliql , &! INTENT(INOUT)
       wsnol , &! INTENT(INOUT)
       chl   , &! INTENT(IN   )
       fi    , &! INTENT(IN   )
       npoi  , &! INTENT(IN   )
       hvap  , &! INTENT(IN   )
       cvap  , &! INTENT(IN   )
       ch2o  , &! INTENT(IN   )
       hsub  , &! INTENT(IN   )
       cice  , &! INTENT(IN   )
       dtime , &! INTENT(IN   )
       hfus  , &! INTENT(IN   )
       vegtype0 , &! INTENT(IN   )
       tmelt     ,&! INTENT(IN   )
       nVegClass)! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    ! at end of timestep, removes evaporation from intercepted h2o,
    ! and does final heat-conserving adjustment for any liquid/snow 
    ! below/above melt point. calls steph2o2 for upper leaves, 
    ! upper stems and lower veg in turn.
    !
    IMPLICIT NONE
    !
      INTEGER , INTENT(IN   ) :: nVegClass
    INTEGER , INTENT(IN   ) :: npoi           ! total number of land points
    REAL(KIND=r8) , INTENT(IN   ) :: hvap           ! latent heat of vaporization of water (J kg-1)
    REAL(KIND=r8) , INTENT(IN   ) :: cvap           ! specific heat of water vapor at constant pressure (J deg-1 
    REAL(KIND=r8) , INTENT(IN   ) :: ch2o           ! specific heat of liquid water (J deg-1 kg-1
    REAL(KIND=r8) , INTENT(IN   ) :: hsub           ! latent heat of sublimation of ice 
    REAL(KIND=r8) , INTENT(IN   ) :: cice           ! specific heat of ice (J deg-1 
    REAL(KIND=r8) , INTENT(IN   ) :: dtime          ! model timestep (seconds)
    REAL(KIND=r8) , INTENT(IN   ) :: hfus           ! latent heat of fusion of water (J 
    REAL(KIND=r8) , INTENT(IN   ) :: tmelt          ! freezing point of water (K)
    REAL(KIND=r8) , INTENT(IN   ) :: fi    (npoi)   ! fractional snow cover
    REAL(KIND=r8) , INTENT(IN   ) :: fu    (npoi)   ! fraction of overall area covered by upper canopy
    REAL(KIND=r8) , INTENT(IN   ) :: lai   (npoi,2) ! canopy single-sided leaf area index (area leaf/area veg)
    REAL(KIND=r8) , INTENT(INOUT) :: tu    (npoi)   ! temperature of upper canopy leaves (K)
    REAL(KIND=r8) , INTENT(INOUT) :: wliqu (npoi)   ! intercepted liquid h2o on upper canopy leaf area (kg m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: wsnou (npoi)   ! intercepted frozen h2o (snow) on upper canopy leaf area (kg m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: chu   (1:nVegClass)   ! heat capacity of upper canopy leaves per unit leaf area (J kg-1 m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: sai   (npoi,2) ! current single-sided stem area index
    REAL(KIND=r8) , INTENT(INOUT) :: ts    (npoi)   ! temperature of upper canopy stems (K)
    REAL(KIND=r8) , INTENT(INOUT) :: wliqs (npoi)   ! intercepted liquid h2o on upper canopy stem area (kg m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: wsnos (npoi)   ! intercepted frozen h2o (snow) on upper canopy stem area (kg m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: chs   (1:nVegClass)   ! heat capacity of upper canopy stems per unit stem area (J kg-1 m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: fl    (npoi)   ! fraction of snow-free area covered by lower  canopy
    REAL(KIND=r8) , INTENT(INOUT) :: tl    (npoi)   ! temperature of lower canopy leaves & stems(K)
    REAL(KIND=r8) , INTENT(INOUT) :: wliql (npoi)   ! intercepted liquid h2o on lower canopy leaf and stem area (kg m-2)
    REAL(KIND=r8) , INTENT(INOUT) :: wsnol (npoi)   ! intercepted frozen h2o (snow) on lower canopy leaf & stem area (kg m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: chl (1:nVegClass)     ! heat capacity of lower canopy leaves & stems per unit leaf/stem area (J kg-1 m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: rliqu (npoi)   ! proportion of fwetu due to liquid
    REAL(KIND=r8) , INTENT(IN   ) :: fvapuw(npoi)! h2o vapor flux (evaporation from wet parts) between upper canopy 
    ! leaves and air at z12 (kg m-2 s-1/ LAI upper canopy/ fu)
    REAL(KIND=r8) , INTENT(INOUT) :: fvapa (npoi)! downward h2o vapor flux between za & z12 at za (kg m-2 s-1)
    REAL(KIND=r8) , INTENT(INOUT) :: fsena (npoi)! downward sensible heat flux between za & z12 at za (W m-2)
    REAL(KIND=r8) , INTENT(IN   ) :: rliqs (npoi)! proportion of fwets due to liquid
    REAL(KIND=r8) , INTENT(IN   ) :: fvaps (npoi)! h2o vapor flux (evaporation from wet surface) between upper canopy
    !  stems and air at z12 (kg m-2 s-1 / SAI lower canopy / fu)
    REAL(KIND=r8) , INTENT(IN   ) :: rliql (npoi)! proportion of fwetl due to liquid
    REAL(KIND=r8) , INTENT(IN   ) :: fvaplw(npoi)   ! h2o vapor flux (evaporation from wet surface) between lower canopy
    !  leaves & stems and air at z34 (kg m-2 s-1/ LAI lower canopy/ fl)
    REAL(KIND=r8) , INTENT(IN   ) :: ta    (npoi)   ! air temperature (K)
    REAL(KIND=r8) , INTENT(IN   ) :: vegtype0(npoi)

    !
    ! local variables
    !
    INTEGER :: i               ! loop indice
    REAL(KIND=r8) :: fveg(npoi)    ! fractional areal coverage of veg component
    REAL(KIND=r8) :: xai (npoi)    ! lai and/or sai for veg component
    !
    ! ---------------------------------------------------------------------
    !
    ! set up for upper leaves
    !
    DO i=1,npoi
       fveg(i) = fu(i)
       xai(i) = 2.00_r8 * lai(i,2)
    END DO
    !
    ! step upper leaves
    !
    CALL steph2o2 (tu    ,& ! INTENT(INOUT)
         wliqu ,& ! INTENT(INOUT)
         wsnou ,& ! INTENT(INOUT)
         fveg  ,& ! INTENT(IN   )
         xai   ,& ! INTENT(IN   )
         rliqu ,& ! INTENT(IN   )
         fvapuw,& ! INTENT(IN   )
         chu   ,& ! INTENT(IN   )
         fvapa ,& ! INTENT(INOUT)
         fsena ,& ! INTENT(INOUT)
         ta    ,& ! INTENT(IN   )
         npoi  ,& ! INTENT(IN   )
         hvap  ,& ! INTENT(IN   )
         cvap  ,& ! INTENT(IN   )
         ch2o  ,& ! INTENT(IN   )
         hsub  ,& ! INTENT(IN   ) 
         cice  ,& ! INTENT(IN   )
         dtime ,& ! INTENT(IN   )
         vegtype0,& ! INTENT(IN   )
         hfus  ,& ! INTENT(IN   )
         tmelt , &! INTENT(IN   )
         nVegClass  )! INTENT(IN   )
    !
    ! set up for upper stems
    !
    DO  i=1,npoi
       fveg(i) = fu(i)
       xai(i) = 2.00_r8 * sai(i,2)
    END DO
    !
    ! step upper stems
    !
    CALL steph2o2 (ts    ,&! INTENT(INOUT)
         wliqs ,&! INTENT(INOUT)
         wsnos ,&! INTENT(INOUT)
         fveg  ,&! INTENT(IN   )
         xai   ,&! INTENT(IN   )
         rliqs ,&! INTENT(IN   )
         fvaps ,&! INTENT(IN   )
         chs(1:nVegClass),&! INTENT(IN   )
         fvapa ,&! INTENT(INOUT)
         fsena ,&! INTENT(INOUT)
         ta    ,&! INTENT(IN   )
         npoi  ,&! INTENT(IN   )
         hvap  ,&! INTENT(IN   )
         cvap  ,&! INTENT(IN   )
         ch2o  ,&! INTENT(IN   )
         hsub  ,&! INTENT(IN   ) 
         cice  ,&! INTENT(IN   )
         dtime ,&! INTENT(IN   )
         vegtype0,& ! INTENT(IN   )
         hfus  ,&! INTENT(IN   )
         tmelt , &! INTENT(IN   )
         nVegClass  )! INTENT(IN   )
    !
    ! set up for lower veg
    !
    DO i=1,npoi
       fveg(i) = (1.0_r8-fi(i))*fl(i)
       xai(i) = 2.00_r8 * (lai(i,1) + sai(i,1))
    END DO
    !
    ! step lower veg
    !
    CALL steph2o2 (tl    ,& ! INTENT(INOUT)
         wliql ,& ! INTENT(INOUT)
         wsnol ,& ! INTENT(INOUT)
         fveg  ,& ! INTENT(IN   )
         xai   ,& ! INTENT(IN   )
         rliql ,& ! INTENT(IN   )
         fvaplw,& ! INTENT(IN   )
         chl(1:nVegClass),& ! INTENT(IN   )
         fvapa ,& ! INTENT(INOUT)
         fsena ,& ! INTENT(INOUT)
         ta    ,& ! INTENT(IN   )
         npoi  ,& ! INTENT(IN   )
         hvap  ,& ! INTENT(IN   )
         cvap  ,& ! INTENT(IN   )
         ch2o  ,& ! INTENT(IN   )
         hsub  ,& ! INTENT(IN   ) 
         cice  ,& ! INTENT(IN   )
         dtime ,& ! INTENT(IN   )
         vegtype0,& ! INTENT(IN   )
         hfus  ,& ! INTENT(IN   )
         tmelt, &! INTENT(IN   )
         nVegClass  )! INTENT(IN   )
    !
    RETURN
  END SUBROUTINE cascad2
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE steph2o2 (tveg   , &! INTENT(INOUT)
       wliq   , &! INTENT(INOUT)
       wsno   , &! INTENT(INOUT)
       fveg   , &! INTENT(IN   )
       xai    , &! INTENT(IN   )
       rliq   , &! INTENT(IN   )
       fvapw  , &! INTENT(IN   )
       cveg   , &! INTENT(IN   )
       fvapa  , &! INTENT(INOUT)
       fsena  , &! INTENT(INOUT)
       ta     , &! INTENT(IN   )
       npoi   , &! INTENT(IN   )
       hvap   , &! INTENT(IN   )
       cvap   , &! INTENT(IN   )
       ch2o   , &! INTENT(IN   )
       hsub   , &! INTENT(IN   ) 
       cice   , &! INTENT(IN   )
       dtime  , &! INTENT(IN   )
       vegtype0,& ! INTENT(IN   )
       hfus   , &! INTENT(IN   )
       tmelt  , &! INTENT(IN   )
       nVegClass  )! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    ! removes evaporation from intercepted h2o, and does final
    ! heat-conserving adjustment for any liquid/snow below/above
    ! melt point, for one veg component
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN   ) :: nVegClass
    INTEGER, INTENT(IN   ) :: npoi           ! total number of land points
    REAL(KIND=r8), INTENT(IN   ) :: hvap           ! latent heat of vaporization of water (J kg-1)
    REAL(KIND=r8), INTENT(IN   ) :: cvap           ! specific heat of water vapor at constant pressure (J deg-1 
    REAL(KIND=r8), INTENT(IN   ) :: ch2o           ! specific heat of liquid water (J deg-1 kg-1
    REAL(KIND=r8), INTENT(IN   ) :: hsub           ! latent heat of sublimation of ice 
    REAL(KIND=r8), INTENT(IN   ) :: cice           ! specific heat of ice (J deg-1 
    REAL(KIND=r8), INTENT(IN   ) :: dtime          ! model timestep (seconds)
    REAL(KIND=r8), INTENT(IN   ) :: hfus           ! latent heat of fusion of water (J 
    REAL(KIND=r8), INTENT(IN   ) :: tmelt          ! freezing point of water (K)
    REAL(KIND=r8), INTENT(IN   ) :: ta    (npoi)   ! air temperature (K)
    REAL(KIND=r8), INTENT(INOUT) :: fvapa (npoi)   ! downward h2o vapor flux between za & z12 at za (kg m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: fsena (npoi)   ! downward sensible heat flux between za & z12 at za (W m-2)
    !
    ! Arguments (all arguments are supplied unless otherwise noted)
    !
    REAL(KIND=r8), INTENT(IN   ) :: cveg(nVegClass)           ! specific heat of veg component ch[u,s,l] 
    !
    REAL(KIND=r8), INTENT(INOUT) :: tveg (npoi)    ! temperature of veg component t[u,s,l] (returned)
    REAL(KIND=r8), INTENT(INOUT) :: wliq (npoi)    ! intercepted liquid amount wliq[u,s,l] (returned)
    REAL(KIND=r8), INTENT(INOUT) :: wsno (npoi)    ! intercepted snow amount wsno[u,s,l] (returned)
    REAL(KIND=r8), INTENT(IN   ) :: fveg (npoi)    ! fractional areal coverage, fu or (1-fi)*fl
    REAL(KIND=r8), INTENT(IN   ) :: xai  (npoi)    ! lai, sai, lai+sai for upper leaves/stems,lower veg
    REAL(KIND=r8), INTENT(IN   ) :: rliq (npoi)    ! ratio of area wetted by liquid to total wetted area
    REAL(KIND=r8), INTENT(IN   ) :: fvapw(npoi)    ! wetted evap h2o flx per leaf/stem area fvap[uw,s,lw]
    REAL(KIND=r8), INTENT(IN   ) :: vegtype0(npoi)    
    !
    ! local variables
    !
    INTEGER :: iveg              ! loopi indice
    INTEGER :: i              ! loopi indice
    !
    REAL(KIND=r8) :: zm                 ! to compute corrective fluxes
    REAL(KIND=r8) :: rwork          ! 1/specific heat of fusion 
    REAL(KIND=r8) :: chav                ! average specific heat for veg, liw and snow
    !
    REAL(KIND=r8) :: dh(npoi)       ! correct heat flux for liquid below melt point and opposite
    REAL(KIND=r8) :: dw(npoi)       ! correct water flux for liquid below melt point and opposite
    REAL(KIND=r8) :: hearhvap
    REAL(KIND=r8) :: heahsub
    REAL(KIND=r8) :: heacvap
    !
    !

    !
    ! statement functions tsatl,tsati are used below so that lowe's
    ! polyomial for liquid is used if t gt 273.16, or for ice if 
    ! t lt 273.16. also impose range of validity for lowe's polys.
    !
    !      REAL(KIND=r8) :: t             ! temperature argument of statement function 
    !      REAL(KIND=r8) :: tair          ! temperature argument of statement function 
    !      REAL(KIND=r8) :: hvapf         ! 
    !              REAL(KIND=r8) :: hsubf         !
    !
    ! statement functions hvapf, hsubf correct the latent heats of
    ! vaporization (liquid-vapor) and sublimation (ice-vapor) to
    ! allow for the concept that the phase change takes place at
    ! 273.16, and the various phases are cooled/heated to that 
    ! temperature before/after the change. this concept is not
    ! physical but is needed to balance the "black-box" energy 
    ! budget. similar correction is applied in convad in the agcm
    ! for precip. needs common comgrd for the physical constants.
    ! argument t is the temp of the liquid or ice, and tair is the
    ! temp of the delivered or received vapor.
    !
    !      hvapf(t,tair) = hvap + cvap*(tair-273.16) - ch2o*(t-273.16)
    !      hsubf(t,tair) = hsub + cvap*(tair-273.16) - cice*(t-273.16)
    !
    !
    ! ---------------------------------------------------------------------
    !
    ! step intercepted h2o due to evaporation/sublimation.
    ! (fvapw already has been multiplied by fwet factor in turvap,
    ! so it is per unit leaf/stem area.)
    !
    ! due to linear fwet factors (see comments in fwetcal) and
    ! the cap on suw,ssw,slw in turvap, evaporation in one timestep
    ! should hardly ever make wliq or wsno negative -- but if this
    ! happens, compensate by increasing vapor flux from atmosphere, 
    ! and decreasing sensib heat flux from atmos (the former is
    ! dangerous since it could suck moisture out of a dry atmos,
    ! and both are unphysical but do fix the budget) tveg in hvapf
    ! and hsubf should be pre-turvap-timestep values, but are not
    !
    hearhvap=hvap
    heahsub=hsub
    heacvap=cvap
    DO  i = 1, npoi
       !
       !WRITE(*,*)wliq(i), rliq(i),fvapw(i)

       wliq(i) = wliq(i) - dtime *     rliq(i)  * fvapw(i)
       wsno(i) = wsno(i) - dtime * (1.0_r8-rliq(i)) * fvapw(i)
       !
       ! check to see if predicted wliq or wsno are less than zero
       !
       IF ((wliq(i).LT.0.0_r8 .OR. wsno(i) .LT. 0.0_r8)  &
            .AND. fveg(i)*xai(i).GT.0.0_r8 )  THEN
          !
          !         write (*,9999) i, wliq(i), wsno(i)
          !9999     format(' ***warning: wliq<0 or wsno<0 -- steph2o2 9999',
          !    >           ' i, wliq, wsno:',i4, 2f12.6)
          !
          ! calculate corrective fluxes
          !
          zm = MAX (-wliq(i), 0.0_r8) * fveg(i) * xai(i) / dtime
          fvapa(i) = fvapa(i) + zm
          fsena(i) = fsena(i) - zm*hvapf(tveg(i),ta(i))
          wliq(i) = MAX (wliq(i), 0.0_r8)
          !
          zm = MAX (-wsno(i), 0.0_r8) * fveg(i) * xai(i) / dtime
          fvapa(i) = fvapa(i) + zm
          fsena(i) = fsena(i) - zm*hsubf(tveg(i),ta(i))
          wsno(i) = MAX (wsno(i), 0.0_r8)
          !
       END IF
       !
    END DO
    !
    ! final heat-conserving correction for liquid/snow below/above
    ! melting point
    !
    rwork = 1.0_r8 / hfus
    !
    DO i=1,npoi
        iveg=vegtype0(i)
       !
       chav = cveg(iveg) + ch2o*wliq(i) + cice*wsno(i)
       !
       ! correct for liquid below melt point
       !
       ! (nb: if tveg > tmelt or wliq = 0, nothing changes.)
       !
       IF (tveg(i).LT.tmelt .AND. wliq(i).GT.0.0_r8) THEN
          dh(i) = chav*(tmelt - tveg(i))
          dw(i) = MIN (wliq(i), MAX (0.0_r8, dh(i)*rwork))
          wliq(i) = wliq(i) - dw(i)
          wsno(i) = wsno(i) + dw(i) 
          chav = cveg(iveg) + ch2o*wliq(i) + cice*wsno(i)
          tveg(i) = tmelt - (dh(i)-hfus*dw(i))/chav
       END IF
       !
       ! correct for snow above melt point
       !
       ! (nb: if tveg < tmelt or wsno = 0, nothing changes.)
       !
       IF (tveg(i).GT.tmelt .AND. wsno(i).GT.0.0_r8) THEN
          dh(i) = chav*(tveg(i) - tmelt)
          dw(i) = MIN (wsno(i), MAX (0.0_r8, dh(i)*rwork))
          wsno(i) = wsno(i) - dw(i)
          wliq(i) = wliq(i) + dw(i)
          chav = cveg(iveg) + ch2o*wliq(i) + cice*wsno(i)
          tveg(i) = tmelt + (dh(i)-hfus*dw(i))/chav
       END IF
       !
    END DO
    !
    RETURN
  END SUBROUTINE steph2o2
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE mix (&
       xm       , &! INTENT(OUT  )
       tm       , &! INTENT(OUT  )
       x1       , &! INTENT(IN   )
       t1       , &! INTENT(IN   )
       x2       , &! INTENT(IN   )
       t2       , &! INTENT(IN   )
       x3       , &! INTENT(IN   )
       t3       , &! INTENT(IN   )
       npoi     , &! INTENT(IN   )
       epsilon    )! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    ! calorimetrically mixes masses x1,x2,x3 with temperatures
    ! t1,t2,t3 into combined mass xm with temperature tm
    !
    ! xm,tm may be returned into same location as one of x1,t1,..,
    ! so hold result temporarily in xtmp,ttmp below
    !
    ! will work if some of x1,x2,x3 have opposite signs, but may 
    ! give unphysical tm's
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN   ) :: npoi    ! total number of land points
    REAL(KIND=r8), INTENT(IN   ) :: epsilon ! small quantity to avoid zero-divides and other
    ! truncation or machine-limit troubles with small
    ! values. should be slightly greater than o(1)
    ! machine precision
    !
    ! Arguments (input except for xm, tm)
    !
    REAL(KIND=r8) , INTENT(OUT  ) :: xm(npoi)     ! resulting mass  
    REAL(KIND=r8) , INTENT(OUT  ) :: tm(npoi)     ! resulting temp
    REAL(KIND=r8) , INTENT(IN   ) :: x1(npoi)     ! mass 1
    REAL(KIND=r8) , INTENT(IN   ) :: t1(npoi)     ! temp 1
    REAL(KIND=r8) , INTENT(IN   ) :: x2(npoi)     ! mass 2
    REAL(KIND=r8) , INTENT(IN   ) :: t2(npoi)     ! temp 2
    REAL(KIND=r8) , INTENT(IN   ) :: x3(npoi)     ! mass 3
    REAL(KIND=r8) , INTENT(IN   ) :: t3(npoi)     ! temp 3
    !
    ! local variables
    !
    INTEGER :: i            ! loop indice
    !
    REAL(KIND=r8) :: xtmp         ! resulting mass (storing variable)
    REAL(KIND=r8) :: ytmp         !  "
    REAL(KIND=r8) :: ttmp         ! resulting temp
    !
    ! ---------------------------------------------------------------------
    !
    DO  i=1,npoi
       !
       xtmp = x1(i) + x2(i) + x3(i)
       !
       ytmp = SIGN (MAX (ABS(xtmp), epsilon), xtmp)
       !
       IF (ABS(xtmp).GE.epsilon) THEN
          ttmp = (t1(i)*x1(i) + t2(i)*x2(i) + t3(i)*x3(i)) / ytmp
       ELSE
          ttmp = 0.0_r8
          xtmp = 0.0_r8
       END IF
       !
       xm(i) = xtmp
       tm(i) = ttmp
       !
    END DO
    !
    RETURN
  END SUBROUTINE mix
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE noveg(lai  ,&! INTENT(IN   )
       fu   ,&! INTENT(IN   )
       tu   ,&! INTENT(INOUT)
       wliqu,&! INTENT(INOUT)
       sai  ,&! INTENT(IN   )
       ts   ,&! INTENT(INOUT)
       wliqs,&! INTENT(INOUT)
       wsnos,&! INTENT(INOUT)
       fl   ,&! INTENT(IN   )
       tl   ,&! INTENT(INOUT)
       wliql,&! INTENT(INOUT)
       wsnol,&! INTENT(INOUT)
       wsnou,&! INTENT(INOUT)
       tg   ,&! INTENT(IN   )
       ti   ,&! INTENT(IN   )
       fi   ,&! INTENT(IN   )
       npoi  )! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    ! if no veg surfaces exist, set prog vars to nominal values
    !
    ! (sensible fluxes fsen[u,s,l], latent fluxes fvap[u,s,l]*, 
    ! temperature t[u,s,l], and intercepted liquid, snow amounts 
    ! wliq[u,s,l], wsno[u,s,l] have been calculated for a unit 
    ! leaf/stem surface, whether or not one exists.)
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN   ) :: npoi          ! total number of land points
    REAL(KIND=r8), INTENT(IN   ) :: fi   (npoi)   ! fractional snow cover
    REAL(KIND=r8), INTENT(IN   ) :: tg   (npoi)   ! soil skin temperature (K)
    REAL(KIND=r8), INTENT(IN   ) :: ti   (npoi)   ! snow skin temperature (K)
    REAL(KIND=r8), INTENT(IN   ) :: lai  (npoi,2) ! canopy single-sided leaf area index (area leaf/area veg)
    REAL(KIND=r8), INTENT(IN   ) :: fu   (npoi)   ! fraction of overall area covered by upper canopy
    REAL(KIND=r8), INTENT(INOUT) :: tu   (npoi)   ! temperature of upper canopy leaves (K)
    REAL(KIND=r8), INTENT(INOUT) :: wliqu(npoi)   ! intercepted liquid h2o on upper canopy leaf area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: sai  (npoi,2) ! current single-sided stem area index
    REAL(KIND=r8), INTENT(INOUT) :: ts   (npoi)   ! temperature of upper canopy stems (K)
    REAL(KIND=r8), INTENT(INOUT) :: wliqs(npoi)   ! intercepted liquid h2o on upper canopy stem area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: wsnos(npoi)   ! intercepted frozen h2o (snow) on upper canopy stem area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: fl   (npoi)   ! fraction of snow-free area covered by lower  canopy
    REAL(KIND=r8), INTENT(INOUT) :: tl   (npoi)   ! temperature of lower canopy leaves & stems(K)
    REAL(KIND=r8), INTENT(INOUT) :: wliql(npoi)   ! intercepted liquid h2o on lower canopy leaf and stem area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: wsnol(npoi)   ! intercepted frozen h2o (snow) on lower canopy leaf & stem area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: wsnou(npoi)   ! intercepted frozen h2o (snow) on upper canopy leaf area (kg m-2)
    !
    ! local variables
    !
    INTEGER :: i   ! loop indice
    !
    REAL(KIND=r8) :: tav ! average temp for soil and snow 
    REAL(KIND=r8) :: x   ! total lai + sai
    REAL(KIND=r8) :: y   ! fraction of lower canopy not snow covered 
    !
    DO i = 1, npoi
       !
       tav = (1.0_r8-fi(i))*tg(i) + fi(i)*ti(i)
       !
       IF  (lai(i,2).EQ.0.0_r8 .OR. fu(i).EQ.0.0_r8) THEN
          tu(i) = tav
          wliqu(i) = 0.0_r8
          wsnou(i) = 0.0_r8
       END IF
       !
       IF (sai(i,2).EQ.0.0_r8 .OR. fu(i).EQ.0.0_r8) THEN
          ts(i) = tav
          wliqs(i) = 0.0_r8
          wsnos(i) = 0.0_r8
       END IF
       !
       x = 2.0_r8 * (lai(i,1) + sai(i,1))
       y = fl(i)*(1.0_r8-fi(i))
       !
       IF (x .EQ.0.0_r8 .OR. y.EQ.0.0_r8) THEN
          tl(i) = tav 
          wliql(i) = 0.0_r8
          wsnol(i) = 0.0_r8
       END IF
       !
    END DO
    !
    RETURN
  END SUBROUTINE noveg
  !
  ! ------------------------------------------------------------------------
  REAL(KIND=r8) FUNCTION twet3(tak, q, p)
    ! ------------------------------------------------------------------------
    !
    ! twet3.f last update 8/30/2000 C Molling
    !
    ! This function calculates the wet bulb temperature given
    ! air temp, specific humidity, and air pressure.  It needs the function esat
    ! in order to work (in comsat.h).  The function is an approximation to
    ! the actual wet bulb temperature relationship.  It agrees well with the
    ! formula in the Smithsonian Met. Tables for moderate humidities, but differs
    ! by as much as 1 K in extremely dry or moist environments.
    !
    ! INPUT
    !     tak - air temp in K
    !     q - specific humidity in kg/kg
    !     p - air pressure in Pa (Pa = 100 * mb)
    !
    ! OUTPUT
    !     twet3 - wet bulb temp in K, accuracy?
    !
    IMPLICIT NONE
    REAL(KIND=r8) , PARAMETER ::  hvap  = 2.5104e+6_r8 
    REAL(KIND=r8) , PARAMETER ::  hfus  = 0.3336e+6_r8   
    REAL(KIND=r8) , PARAMETER ::  cvap  = 1.81e+3_r8    
    REAL(KIND=r8) , PARAMETER ::  ch2o  = 4.218e+3_r8    
    REAL(KIND=r8) , PARAMETER ::  hsub  = hvap + hfus    
    REAL(KIND=r8) , PARAMETER ::  cice  = 2.106e+3_r8    
    REAL(KIND=r8) , PARAMETER ::  cair  = 1.00464e+3_r8 
    !
    INTEGER :: i
    !
    REAL(KIND=r8) :: tak 
    REAL(KIND=r8) :: q 
    REAL(KIND=r8) :: p 
    REAL(KIND=r8) :: ta 
    REAL(KIND=r8) :: twk 
    REAL(KIND=r8) :: twold 
    REAL(KIND=r8) :: diff
    !
    ! ------
    ! comsat
    ! ------
    !
    ! ---------------------------------------------------------------------
    ! statement functions and associated parameters
    ! ---------------------------------------------------------------------
    !
    ! polynomials for svp(t), d(svp)/dt over water and ice are from
    ! lowe(1977),jam,16,101-103.
    !
    !
    REAL(KIND=r8) , PARAMETER :: asat0 =  6.1078000_r8
    REAL(KIND=r8) , PARAMETER :: asat1 =  4.4365185e-1_r8
    REAL(KIND=r8) , PARAMETER :: asat2 =  1.4289458e-2_r8
    REAL(KIND=r8) , PARAMETER :: asat3 =  2.6506485e-4_r8
    REAL(KIND=r8) , PARAMETER :: asat4 =  3.0312404e-6_r8
    REAL(KIND=r8) , PARAMETER :: asat5 =  2.0340809e-8_r8
    REAL(KIND=r8) , PARAMETER :: asat6 =  6.1368209e-11_r8
    !
    REAL(KIND=r8) , PARAMETER :: bsat0 =  6.1091780_r8
    REAL(KIND=r8) , PARAMETER :: bsat1 =  5.0346990e-1_r8
    REAL(KIND=r8) , PARAMETER :: bsat2 =  1.8860134e-2_r8
    REAL(KIND=r8) , PARAMETER :: bsat3 =  4.1762237e-4_r8
    REAL(KIND=r8) , PARAMETER :: bsat4 =  5.8247203e-6_r8
    REAL(KIND=r8) , PARAMETER :: bsat5 =  4.8388032e-8_r8
    REAL(KIND=r8) , PARAMETER :: bsat6 =  1.8388269e-10_r8
    !
    REAL(KIND=r8) , PARAMETER :: csat0 =  4.4381000e-1_r8
    REAL(KIND=r8) , PARAMETER :: csat1 =  2.8570026e-2_r8
    REAL(KIND=r8) , PARAMETER :: csat2 =  7.9380540e-4_r8
    REAL(KIND=r8) , PARAMETER :: csat3 =  1.2152151e-5_r8
    REAL(KIND=r8) , PARAMETER :: csat4 =  1.0365614e-7_r8
    REAL(KIND=r8) , PARAMETER :: csat5 =  3.5324218e-10_r8
    REAL(KIND=r8) , PARAMETER :: csat6 = -7.0902448e-13_r8
    !
    REAL(KIND=r8) , PARAMETER :: dsat0 =  5.0303052e-1_r8
    REAL(KIND=r8) , PARAMETER :: dsat1 =  3.7732550e-2_r8
    REAL(KIND=r8) , PARAMETER :: dsat2 =  1.2679954e-3_r8
    REAL(KIND=r8) , PARAMETER :: dsat3 =  2.4775631e-5_r8
    REAL(KIND=r8) , PARAMETER :: dsat4 =  3.0056931e-7_r8
    REAL(KIND=r8) , PARAMETER :: dsat5 =  2.1585425e-9_r8
    REAL(KIND=r8) , PARAMETER :: dsat6 =  7.1310977e-12_r8
    !
    ! statement functions tsatl,tsati are used below so that lowe's
    ! polyomial for liquid is used if t gt 273.16, or for ice if 
    ! t lt 273.16. also impose range of validity for lowe's polys.
    !
    !      REAL(KIND=r8) :: t        ! temperature argument of statement function 
    !      REAL(KIND=r8) :: tair     ! temperature argument of statement function 
    !      REAL(KIND=r8) :: p1       ! pressure argument of function 
    !      REAL(KIND=r8) :: e1       ! vapor pressure argument of function
    !      REAL(KIND=r8) :: q1       ! saturation specific humidity argument of function
    ! REAL(KIND=r8) :: tsatl    ! statement function
    ! REAL(KIND=r8) :: tsati    ! 
    ! REAL(KIND=r8) :: esat     !
    ! REAL(KIND=r8) :: desat    !
    ! REAL(KIND=r8) :: qsat     ! 
    !  REAL(KIND=r8) :: dqsat    ! 
    !  REAL(KIND=r8) :: hvapf    ! 
    !  REAL(KIND=r8) :: hsubf    !
    !  REAL(KIND=r8) :: cvmgt    ! function
    !
    !tsatl(t) = min (100., max (t-273.16, 0.))
    !tsati(t) = max (-60., min (t-273.16, 0.))
    !
    ! statement function esat is svp in n/m**2, with t in deg k. 
    ! (100 * lowe's poly since 1 mb = 100 n/m**2.)
    !
    ! esat (t) =              &
    !  100.*(             &
    !         cvmgt (asat0, bsat0, t.ge.273.16)             &
    !         + tsatl(t)*(asat1 + tsatl(t)*(asat2 + tsatl(t)*(asat3             &
    !         + tsatl(t)*(asat4 + tsatl(t)*(asat5 + tsatl(t)* asat6)))))  &
    !         + tsati(t)*(bsat1 + tsati(t)*(bsat2 + tsati(t)*(bsat3             &
    !         + tsati(t)*(bsat4 + tsati(t)*(bsat5 + tsati(t)* bsat6)))))  &
    !  )
    !
    ! statement function desat is d(svp)/dt, with t in deg k.
    ! (100 * lowe's poly since 1 mb = 100 n/m**2.)
    !
    !desat (t) =              &
    ! 100.*(              &
    !         cvmgt (csat0, dsat0, t.ge.273.16)              &
    !          + tsatl(t)*(csat1 + tsatl(t)*(csat2 + tsatl(t)*(csat3              &
    !         + tsatl(t)*(csat4 + tsatl(t)*(csat5 + tsatl(t)* csat6)))))   &
    !         + tsati(t)*(dsat1 + tsati(t)*(dsat2 + tsati(t)*(dsat3              &
    !         + tsati(t)*(dsat4 + tsati(t)*(dsat5 + tsati(t)* dsat6)))))   &
    ! )
    !
    ! statement function qsat is saturation specific humidity,
    ! with svp e1 and ambient pressure p in n/m**2. impose an upper
    ! limit of 1 to avoid spurious values for very high svp
    ! and/or small p1
    !
    !       qsat (e1, p1) = 0.622 * e1 /  &
    !                    max ( p1 - (1.0 - 0.622) * e1, 0.622 * e1 )
    !
    ! statement function dqsat is d(qsat)/dt, with t in deg k and q1
    ! in kg/kg (q1 is *saturation* specific humidity)
    !
    !       dqsat (t, q1) = desat(t) * q1 * (1. + q1*(1./0.622 - 1.)) /  &
    !                       esat(t)
    !
    ! statement functions hvapf, hsubf correct the latent heats of
    ! vaporization (liquid-vapor) and sublimation (ice-vapor) to
    ! allow for the concept that the phase change takes place at
    ! 273.16, and the various phases are cooled/heated to that 
    ! temperature before/after the change. this concept is not
    ! physical but is needed to balance the "black-box" energy 
    ! budget. similar correction is applied in convad in the agcm
    ! for precip. needs common comgrd for the physical constants.
    ! argument t is the temp of the liquid or ice, and tair is the
    ! temp of the delivered or received vapor.
    !
    !      hvapf(t,tair) = hvap + cvap*(tair-273.16) - ch2o*(t-273.16)
    !      hsubf(t,tair) = hsub + cvap*(tair-273.16) - cice*(t-273.16)
    !

    !
    ! temperatures in twet3 equation must be in C
    ! pressure in qsat function must be in Pa
    ! temperatures in esat,hvapf functions must be in K
    !
    !     Air temp in C
    !     -------------
    ta = tak - 273.16_r8
    !
    !     First guess for wet bulb temp in C, K
    !     -------------------------------------
    twet3 = ta * q / qsat(esat(tak),p)
    twk = twet3 + 273.16_r8
    !
    !     Iterate to converge
    !     -------------------
    DO  i = 1, 50
       twold = twk - 273.16_r8
       twet3 = ta - (hvapf(twk,tak)/cair) * ( qsat( esat(twk),p )-q )
       diff = twet3 - twold
       !
       ! below, the 0.2 is the relaxation parameter that works up to 40C (at least)
       !
       twk = twold + 0.2_r8 * diff + 273.16_r8
       IF (ABS(twk-273.16_r8-twold) .LT. 0.02_r8) GO TO 999
    END DO
    !
    PRINT *, 'Warning, twet3 failed to converge after 20 iterations!'
    PRINT *, 'twet3, twetold: ', twk, twold+273.16_r8
    PRINT *, 'twetbulb is being set to the air temperature'
    !
    twet3 = tak
    !
    !     Return wet bulb temperature in K
    !     --------------------------------
999 twet3 = twk
    !
    RETURN
  END FUNCTION twet3

  !
  ! ---------------------------------------------------------------------
  SUBROUTINE linsolve (arr     , &! INTENT(INOUT)
       rhs     , &! INTENT(INOUT)
       vec     , &! INTENT(INOUT)
       mplate  , &! INTENT(IN   )
       nd      , &! INTENT(IN   )
       npoi      )! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    ! solves multiple linear systems of equations, vectorizing
    ! over the number of systems. basic gaussian elimination is 
    ! used, with no pivoting (relies on all diagonal elements
    ! being and staying significantly non-zero)
    !
    ! a template array mplate is used to detect when an operation 
    ! is not necessary (element already zero or would add zeros),
    ! assuming that every system has the same pattern of zero
    ! elements
    !
    ! this template is first copied to mplatex since it 
    ! must be updated during the procedure in case an original-zero
    ! pattern location becomes non-zero
    !
    ! the first subscript in arr, rhs, vec is over the multiple
    ! systems, and the others are the usual row, column subscripts
    !
    IMPLICIT NONE
    !
    ! Arguments (input-output)
    !
    INTEGER, INTENT(IN   ) :: npoi                ! total number of land points
    INTEGER, INTENT(IN   ) :: nd                  ! number of equations (supplied)
    !
    INTEGER, INTENT(IN   ) :: mplate(nd,nd)       ! pattern of zero elements of arr (supplied)
    !
    REAL(KIND=r8), INTENT(INOUT) :: arr(npoi,nd,nd)     ! equation coefficients (supplied, overwritten)
    REAL(KIND=r8), INTENT(INOUT) :: rhs(npoi,nd)        ! equation right-hand sides (supplied, overwritten) 
    REAL(KIND=r8), INTENT(INOUT) :: vec(npoi,nd)        ! solution (returned)
    ! 
    ! local variables
    !
    INTEGER :: ndx                 ! Max number of equations
    INTEGER :: j                   ! loop indices
    INTEGER :: i                   ! loop indices
    INTEGER :: id                  ! loop indices
    INTEGER :: m                   ! loop indices
    !
    PARAMETER (ndx=9)
    !
    INTEGER :: mplatex(ndx,ndx)
    !
    REAL(KIND=r8) :: f(npoi)
    !
    IF (nd.GT.ndx) THEN
       WRITE(*,900) nd, ndx
900    FORMAT(/' *** fatal error ***'/  &
            /' number of linsolve eqns',i4,' exceeds limit',i4)
       STOP 'have problem at linsolve'
    END IF
    f=0.0_r8
    !
    ! copy the zero template so it can be changed below
    !
    DO j=1,nd
       DO i=1,nd
          mplatex(i,j) = mplate(i,j)
       END DO
    END DO
    !
    ! zero all array elements below the diagonal, proceeding from
    ! the first row to the last. note that mplatex is set non-zero
    ! for changed (i,j) locations, in loop 20
    !
    DO id=1, nd-1
       DO i=id+1,nd
          !
          IF (mplatex(i,id).NE.0) THEN
             DO  m=1,npoi
                f(m) = arr(m,i,id) / arr(m,id,id)
             END DO
             !
             DO j=id,nd
                IF (mplatex(id,j).NE.0) THEN
                   DO  m=1,npoi
                      arr(m,i,j) = arr(m,i,j) - f(m)*arr(m,id,j)
                   END DO
                   mplatex(i,j) = 1
                END IF
             END DO
             !
             DO m=1,npoi
                rhs(m,i) = rhs(m,i) - f(m)*rhs(m,id)
             END DO
          END IF
          !
       END DO
    END DO
    !
    ! all array elements below the diagonal are zero, so can
    ! immediately solve the equations in reverse order
    !
    DO id=nd,1,-1
       !
       f =0.0_r8
       !call const (f, npoi, 0.0)

       IF (id.LT.nd) THEN
          DO  j=id+1,nd
             IF (mplatex(id,j).NE.0) THEN
                DO m=1,npoi
                   f(m) = f(m) + arr(m,id,j)*vec(m,j)
                END DO
             END IF
          END DO
       END IF
       !
       DO m=1,npoi
          vec(m,id) = (rhs(m,id) - f(m)) / arr(m,id,id)
       END DO
       !
    END DO
    !
    RETURN
  END SUBROUTINE linsolve
  ! ---------------------------------------------------------------------
  REAL(KIND=r8) FUNCTION qsat(e1, p1)
    ! ---------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN   ) :: e1 
    REAL(KIND=r8), INTENT(IN   ) :: p1 

    ! statement function qsat is saturation specific humidity,
    ! with svp e1 and ambient pressure p in n/m**2. impose an upper
    ! limit of 1 to avoid spurious values for very high svp
    ! and/or small p1
    !
    qsat = 0.622_r8 * e1 /  &
         MAX ( p1 - (1.0_r8 - 0.622_r8) * e1, 0.622_r8 * e1 )

  END  FUNCTION qsat

  ! ---------------------------------------------------------------------
  REAL(KIND=r8) FUNCTION desat(t)
  ! ---------------------------------------------------------------------
  !
  ! chooses between two things.  Used in canopy.f 
  !
    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN   ) :: t 

    desat  = 100.0_r8*( cvmgt (csat0, dsat0, t.GE.273.16_r8)             &
         + tsatl(t)*(csat1 + tsatl(t)*(csat2 + tsatl(t)*(csat3             &
         + tsatl(t)*(csat4 + tsatl(t)*(csat5 + tsatl(t)* csat6)))))  &
         + tsati(t)*(dsat1 + tsati(t)*(dsat2 + tsati(t)*(dsat3             &
         + tsati(t)*(dsat4 + tsati(t)*(dsat5 + tsati(t)* dsat6)))))  &
         )
  END  FUNCTION desat

  ! ---------------------------------------------------------------------
  REAL(KIND=r8) FUNCTION esat(t)
    ! ---------------------------------------------------------------------
    IMPLICIT NONE

    REAL(KIND=r8), INTENT(IN   ) :: t 

    esat = 100.0_r8*(cvmgt (asat0, bsat0, t.GE.273.16_r8)             &
         + tsatl(t)*(asat1 + tsatl(t)*(asat2 + tsatl(t)*(asat3             &
         + tsatl(t)*(asat4 + tsatl(t)*(asat5 + tsatl(t)* asat6)))))  &
         + tsati(t)*(bsat1 + tsati(t)*(bsat2 + tsati(t)*(bsat3             &
         + tsati(t)*(bsat4 + tsati(t)*(bsat5 + tsati(t)* bsat6)))))  &
         )
  END  FUNCTION esat
  ! 
  ! ---------------------------------------------------------------------
  REAL(KIND=r8) FUNCTION tsatl (t)
    ! ---------------------------------------------------------------------
    !
    ! chooses between two things.  Used in canopy.f
    !
    IMPLICIT NONE
    !
    REAL(KIND=r8), INTENT(IN   ) :: t
    !
    tsatl = MIN (100.0_r8, MAX (t-273.16_r8, 0.0_r8))
    !
    RETURN
  END  FUNCTION tsatl
  ! ---------------------------------------------------------------------
  REAL(KIND=r8) FUNCTION tsati (t)
    ! ---------------------------------------------------------------------
    !
    ! chooses between two things.  Used in canopy.f
    !
    IMPLICIT NONE
    !
    REAL(KIND=r8), INTENT(IN   ) :: t
    !
    tsati = MAX (-60.0_r8, MIN (t-273.16_r8, 0.0_r8))
    !
    RETURN
  END  FUNCTION tsati

  ! ---------------------------------------------------------------------
  REAL(KIND=r8) FUNCTION dqsat (t, q1)
    ! ---------------------------------------------------------------------
    !
    ! chooses between two things.  Used in canopy.f
    !
    IMPLICIT NONE      
    REAL(KIND=r8), INTENT(IN   ) :: t
    REAL(KIND=r8), INTENT(IN   ) :: q1

    !
    !
    ! statement function dqsat is d(qsat)/dt, with t in deg k and q1
    ! in kg/kg (q1 is *saturation* specific humidity)
    !
    dqsat = desat(t) * q1 * (1.0_r8 + q1*(1.0_r8/0.622_r8 - 1.0_r8)) / &
         esat(t)
    !
    RETURN
  END  FUNCTION dqsat
  ! ---------------------------------------------------------------------
  REAL(KIND=r8) FUNCTION hvapf (t,tair)
    ! ---------------------------------------------------------------------
    !
    ! chooses between two things.  Used in canopy.f
    !
    IMPLICIT NONE
    !
    REAL(KIND=r8), INTENT(IN   ) :: t
    REAL(KIND=r8), INTENT(IN   ) :: tair
    !
    !
    ! statement functions hvapf, hsubf correct the latent heats of
    ! vaporization (liquid-vapor) and sublimation (ice-vapor) to
    ! allow for the concept that the phase change takes place at
    ! 273.16, and the various phases are cooled/heated to that 
    ! temperature before/after the change. this concept is not
    ! physical but is needed to balance the "black-box" energy 
    ! budget. similar correction is applied in convad in the agcm
    ! for precip. needs common comgrd for the physical constants.
    ! argument t is the temp of the liquid or ice, and tair is the
    ! temp of the delivered or received vapor.
    !
    hvapf = hvap + cvap*(tair-273.16_r8) - ch2o*(t-273.16_r8)
    !
    RETURN
  END  FUNCTION hvapf

  ! ---------------------------------------------------------------------
  REAL(KIND=r8) FUNCTION hsubf (t,tair)
    ! ---------------------------------------------------------------------
    IMPLICIT NONE
    !
    REAL(KIND=r8), INTENT(IN   ) :: t
    REAL(KIND=r8), INTENT(IN   ) :: tair

    ! statement functions hvapf, hsubf correct the latent heats of
    ! vaporization (liquid-vapor) and sublimation (ice-vapor) to
    ! allow for the concept that the phase change takes place at
    ! 273.16, and the various phases are cooled/heated to that 
    ! temperature before/after the change. this concept is not
    ! physical but is needed to balance the "black-box" energy 
    ! budget. similar correction is applied in convad in the agcm
    ! for precip. needs common comgrd for the physical constants.
    ! argument t is the temp of the liquid or ice, and tair is the
    ! temp of the delivered or received vapor.

    hsubf = hsub + cvap*(tair-273.16_r8) - cice*(t-273.16_r8)

    RETURN
  END  FUNCTION hsubf

  ! ---------------------------------------------------------------------
  REAL(KIND=r8) FUNCTION cvmgt (x,y,l)
    ! ---------------------------------------------------------------------
    !
    ! chooses between two things.  Used in canopy.f
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN   ) :: l
    REAL(KIND=r8), INTENT(IN   ) :: x
    REAL(KIND=r8), INTENT(IN   ) :: y
    !
    IF (l) THEN
       cvmgt = x
    ELSE
       cvmgt = y
    END IF
    !
    RETURN
  END  FUNCTION cvmgt





      REAL(KIND=r8) FUNCTION es_Sat (T, p)
!
! !DESCRIPTION:
! Computes saturation mixing ratio and the change in saturation
! mixing ratio with respect to temperature.
! Reference:  Polynomial approximations from:
!             Piotr J. Flatau, et al.,1992:  Polynomial fits to saturation
!             vapor pressure.  Journal of Applied Meteorology, 31, 1507-1513.
!
! !USES:
!    use shr_kind_mod , only: r8 => shr_kind_r8
!    use shr_const_mod, only: SHR_CONST_TKFRZ
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: T        ! temperature (K)
    real(r8), intent(in)  :: p        ! surface atmospheric pressure (pa)
!    real(r8), intent(out) :: es       ! vapor pressure (pa)
!    real(r8), intent(out) :: esdT     ! d(es)/d(T)
!    real(r8), intent(out) :: qs       ! humidity (kg/kg)
!    real(r8), intent(out) :: qsdT     ! d(qs)/d(T)
!
! !CALLED FROM:
! subroutine Biogeophysics1 in module Biogeophysics1Mod
! subroutine BiogeophysicsLake in module BiogeophysicsLakeMod
! subroutine CanopyFluxesMod CanopyFluxesMod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
!
!
! !LOCAL VARIABLES:
!EOP
!
    real(r8) :: T_limit
    real(r8) :: td,vp,vp1,vp2,es
!
! For water vapor (temperature range 0C-100C)
!
    real(r8), parameter :: a0 =  6.11213476_r8
    real(r8), parameter :: a1 =  0.444007856_r8
    real(r8), parameter :: a2 =  0.143064234e-01_r8
    real(r8), parameter :: a3 =  0.264461437e-03_r8
    real(r8), parameter :: a4 =  0.305903558e-05_r8
    real(r8), parameter :: a5 =  0.196237241e-07_r8
    real(r8), parameter :: a6 =  0.892344772e-10_r8
    real(r8), parameter :: a7 = -0.373208410e-12_r8
    real(r8), parameter :: a8 =  0.209339997e-15_r8
!
! For derivative:water vapor
!
    real(r8), parameter :: b0 =  0.444017302_r8
    real(r8), parameter :: b1 =  0.286064092e-01_r8
    real(r8), parameter :: b2 =  0.794683137e-03_r8
    real(r8), parameter :: b3 =  0.121211669e-04_r8
    real(r8), parameter :: b4 =  0.103354611e-06_r8
    real(r8), parameter :: b5 =  0.404125005e-09_r8
    real(r8), parameter :: b6 = -0.788037859e-12_r8
    real(r8), parameter :: b7 = -0.114596802e-13_r8
    real(r8), parameter :: b8 =  0.381294516e-16_r8
!
! For ice (temperature range -75C-0C)
!
    real(r8), parameter :: c0 =  6.11123516_r8
    real(r8), parameter :: c1 =  0.503109514_r8
    real(r8), parameter :: c2 =  0.188369801e-01_r8
    real(r8), parameter :: c3 =  0.420547422e-03_r8
    real(r8), parameter :: c4 =  0.614396778e-05_r8
    real(r8), parameter :: c5 =  0.602780717e-07_r8
    real(r8), parameter :: c6 =  0.387940929e-09_r8
    real(r8), parameter :: c7 =  0.149436277e-11_r8
    real(r8), parameter :: c8 =  0.262655803e-14_r8
!
! For derivative:ice
!
    real(r8), parameter :: d0 =  0.503277922_r8
    real(r8), parameter :: d1 =  0.377289173e-01_r8
    real(r8), parameter :: d2 =  0.126801703e-02_r8
    real(r8), parameter :: d3 =  0.249468427e-04_r8
    real(r8), parameter :: d4 =  0.313703411e-06_r8
    real(r8), parameter :: d5 =  0.257180651e-08_r8
    real(r8), parameter :: d6 =  0.133268878e-10_r8
    real(r8), parameter :: d7 =  0.394116744e-13_r8
    real(r8), parameter :: d8 =  0.498070196e-16_r8
!-----------------------------------------------------------------------
    real(R8),parameter :: SHR_CONST_TKFRZ   = 273.15_r8       ! freezing T of fresh water          ~ K


    T_limit = T - SHR_CONST_TKFRZ
    if (T_limit > 100.0_r8) T_limit=100.0_r8
    if (T_limit < -75.0_r8) T_limit=-75.0_r8

    td       = T_limit
    if (td >= 0.0_r8) then
       es   = a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 &
            + td*(a5 + td*(a6 + td*(a7 + td*a8)))))))
!       esdT = b0 + td*(b1 + td*(b2 + td*(b3 + td*(b4 &
!            + td*(b5 + td*(b6 + td*(b7 + td*b8)))))))
    else
       es   = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 &
            + td*(c5 + td*(c6 + td*(c7 + td*c8)))))))
!       esdT = d0 + td*(d1 + td*(d2 + td*(d3 + td*(d4 &
!            + td*(d5 + td*(d6 + td*(d7 + td*d8)))))))
    endif

    es_Sat    = es    * 100._r8            ! pa
!    esdT  = esdT  * 100._r8            ! pa/K
!
!    vp    = 1.0_r8   / (p - 0.378_r8*es)
!    vp1   = 0.622_r8 * vp
!    vp2   = vp1   * vp
!
!    qs    = es    * vp1             ! kg/kg
!    qsdT  = esdT  * vp2 * p         ! 1 / K

  end FUNCTION es_Sat



      REAL(KIND=r8) FUNCTION esdT_Sat (T, p)
!
! !DESCRIPTION:
! Computes saturation mixing ratio and the change in saturation
! mixing ratio with respect to temperature.
! Reference:  Polynomial approximations from:
!             Piotr J. Flatau, et al.,1992:  Polynomial fits to saturation
!             vapor pressure.  Journal of Applied Meteorology, 31, 1507-1513.
!
! !USES:
!    use shr_kind_mod , only: r8 => shr_kind_r8
!    use shr_const_mod, only: SHR_CONST_TKFRZ
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: T        ! temperature (K)
    real(r8), intent(in)  :: p        ! surface atmospheric pressure (pa)
!    real(r8), intent(out) :: es       ! vapor pressure (pa)
!    real(r8), intent(out) :: esdT     ! d(es)/d(T)
!    real(r8), intent(out) :: qs       ! humidity (kg/kg)
!    real(r8), intent(out) :: qsdT     ! d(qs)/d(T)
!
! !CALLED FROM:
! subroutine Biogeophysics1 in module Biogeophysics1Mod
! subroutine BiogeophysicsLake in module BiogeophysicsLakeMod
! subroutine CanopyFluxesMod CanopyFluxesMod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
!
!
! !LOCAL VARIABLES:
!EOP
!
    real(r8) :: T_limit
    real(r8) :: td,vp,vp1,vp2,esdT
!
! For water vapor (temperature range 0C-100C)
!
    real(r8), parameter :: a0 =  6.11213476_r8
    real(r8), parameter :: a1 =  0.444007856_r8
    real(r8), parameter :: a2 =  0.143064234e-01_r8
    real(r8), parameter :: a3 =  0.264461437e-03_r8
    real(r8), parameter :: a4 =  0.305903558e-05_r8
    real(r8), parameter :: a5 =  0.196237241e-07_r8
    real(r8), parameter :: a6 =  0.892344772e-10_r8
    real(r8), parameter :: a7 = -0.373208410e-12_r8
    real(r8), parameter :: a8 =  0.209339997e-15_r8
!
! For derivative:water vapor
!
    real(r8), parameter :: b0 =  0.444017302_r8
    real(r8), parameter :: b1 =  0.286064092e-01_r8
    real(r8), parameter :: b2 =  0.794683137e-03_r8
    real(r8), parameter :: b3 =  0.121211669e-04_r8
    real(r8), parameter :: b4 =  0.103354611e-06_r8
    real(r8), parameter :: b5 =  0.404125005e-09_r8
    real(r8), parameter :: b6 = -0.788037859e-12_r8
    real(r8), parameter :: b7 = -0.114596802e-13_r8
    real(r8), parameter :: b8 =  0.381294516e-16_r8
!
! For ice (temperature range -75C-0C)
!
    real(r8), parameter :: c0 =  6.11123516_r8
    real(r8), parameter :: c1 =  0.503109514_r8
    real(r8), parameter :: c2 =  0.188369801e-01_r8
    real(r8), parameter :: c3 =  0.420547422e-03_r8
    real(r8), parameter :: c4 =  0.614396778e-05_r8
    real(r8), parameter :: c5 =  0.602780717e-07_r8
    real(r8), parameter :: c6 =  0.387940929e-09_r8
    real(r8), parameter :: c7 =  0.149436277e-11_r8
    real(r8), parameter :: c8 =  0.262655803e-14_r8
!
! For derivative:ice
!
    real(r8), parameter :: d0 =  0.503277922_r8
    real(r8), parameter :: d1 =  0.377289173e-01_r8
    real(r8), parameter :: d2 =  0.126801703e-02_r8
    real(r8), parameter :: d3 =  0.249468427e-04_r8
    real(r8), parameter :: d4 =  0.313703411e-06_r8
    real(r8), parameter :: d5 =  0.257180651e-08_r8
    real(r8), parameter :: d6 =  0.133268878e-10_r8
    real(r8), parameter :: d7 =  0.394116744e-13_r8
    real(r8), parameter :: d8 =  0.498070196e-16_r8
!-----------------------------------------------------------------------
    real(R8),parameter :: SHR_CONST_TKFRZ   = 273.15_r8       ! freezing T of fresh water          ~ K


    T_limit = T - SHR_CONST_TKFRZ
    if (T_limit > 100.0_r8) T_limit=100.0_r8
    if (T_limit < -75.0_r8) T_limit=-75.0_r8

    td       = T_limit
    if (td >= 0.0_r8) then
!       es   = a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 &
!            + td*(a5 + td*(a6 + td*(a7 + td*a8)))))))
       esdT = b0 + td*(b1 + td*(b2 + td*(b3 + td*(b4 &
            + td*(b5 + td*(b6 + td*(b7 + td*b8)))))))
    else
!       es   = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 &
!            + td*(c5 + td*(c6 + td*(c7 + td*c8)))))))
       esdT = d0 + td*(d1 + td*(d2 + td*(d3 + td*(d4 &
            + td*(d5 + td*(d6 + td*(d7 + td*d8)))))))
    endif

!    es    = es    * 100._r8            ! pa
    esdT_Sat  = esdT  * 100._r8            ! pa/K

!    vp    = 1.0_r8   / (p - 0.378_r8*es)
!    vp1   = 0.622_r8 * vp
!    vp2   = vp1   * vp
!
!    qs    = es    * vp1             ! kg/kg
!    qsdT  = esdT  * vp2 * p         ! 1 / K

  end FUNCTION esdT_Sat
  
  
      REAL(KIND=r8) FUNCTION qs_Sat (T, p)
!
! !DESCRIPTION:
! Computes saturation mixing ratio and the change in saturation
! mixing ratio with respect to temperature.
! Reference:  Polynomial approximations from:
!             Piotr J. Flatau, et al.,1992:  Polynomial fits to saturation
!             vapor pressure.  Journal of Applied Meteorology, 31, 1507-1513.
!
! !USES:
!    use shr_kind_mod , only: r8 => shr_kind_r8
!    use shr_const_mod, only: SHR_CONST_TKFRZ
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: T        ! temperature (K)
    real(r8), intent(in)  :: p        ! surface atmospheric pressure (pa)
!    real(r8), intent(out) :: es       ! vapor pressure (pa)
!    real(r8), intent(out) :: esdT     ! d(es)/d(T)
!    real(r8), intent(out) :: qs       ! humidity (kg/kg)
!    real(r8), intent(out) :: qsdT     ! d(qs)/d(T)
!
! !CALLED FROM:
! subroutine Biogeophysics1 in module Biogeophysics1Mod
! subroutine BiogeophysicsLake in module BiogeophysicsLakeMod
! subroutine CanopyFluxesMod CanopyFluxesMod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
!
!
! !LOCAL VARIABLES:
!EOP
!
    real(r8) :: T_limit
    real(r8) :: td,vp,vp1,vp2
    real(r8) :: es  ! vapor pressure (pa)
!    real(r8) :: esdT  ! d(es)/d(T)
!    real(r8) :: qs  ! humidity (kg/kg)
!    real(r8) :: qsdT  ! d(qs)/d(T)

!
! For water vapor (temperature range 0C-100C)
!
    real(r8), parameter :: a0 =  6.11213476_r8
    real(r8), parameter :: a1 =  0.444007856_r8
    real(r8), parameter :: a2 =  0.143064234e-01_r8
    real(r8), parameter :: a3 =  0.264461437e-03_r8
    real(r8), parameter :: a4 =  0.305903558e-05_r8
    real(r8), parameter :: a5 =  0.196237241e-07_r8
    real(r8), parameter :: a6 =  0.892344772e-10_r8
    real(r8), parameter :: a7 = -0.373208410e-12_r8
    real(r8), parameter :: a8 =  0.209339997e-15_r8
!
! For derivative:water vapor
!
    real(r8), parameter :: b0 =  0.444017302_r8
    real(r8), parameter :: b1 =  0.286064092e-01_r8
    real(r8), parameter :: b2 =  0.794683137e-03_r8
    real(r8), parameter :: b3 =  0.121211669e-04_r8
    real(r8), parameter :: b4 =  0.103354611e-06_r8
    real(r8), parameter :: b5 =  0.404125005e-09_r8
    real(r8), parameter :: b6 = -0.788037859e-12_r8
    real(r8), parameter :: b7 = -0.114596802e-13_r8
    real(r8), parameter :: b8 =  0.381294516e-16_r8
!
! For ice (temperature range -75C-0C)
!
    real(r8), parameter :: c0 =  6.11123516_r8
    real(r8), parameter :: c1 =  0.503109514_r8
    real(r8), parameter :: c2 =  0.188369801e-01_r8
    real(r8), parameter :: c3 =  0.420547422e-03_r8
    real(r8), parameter :: c4 =  0.614396778e-05_r8
    real(r8), parameter :: c5 =  0.602780717e-07_r8
    real(r8), parameter :: c6 =  0.387940929e-09_r8
    real(r8), parameter :: c7 =  0.149436277e-11_r8
    real(r8), parameter :: c8 =  0.262655803e-14_r8
!
! For derivative:ice
!
    real(r8), parameter :: d0 =  0.503277922_r8
    real(r8), parameter :: d1 =  0.377289173e-01_r8
    real(r8), parameter :: d2 =  0.126801703e-02_r8
    real(r8), parameter :: d3 =  0.249468427e-04_r8
    real(r8), parameter :: d4 =  0.313703411e-06_r8
    real(r8), parameter :: d5 =  0.257180651e-08_r8
    real(r8), parameter :: d6 =  0.133268878e-10_r8
    real(r8), parameter :: d7 =  0.394116744e-13_r8
    real(r8), parameter :: d8 =  0.498070196e-16_r8
!-----------------------------------------------------------------------
    real(R8),parameter :: SHR_CONST_TKFRZ   = 273.15_r8       ! freezing T of fresh water          ~ K


    T_limit = T - SHR_CONST_TKFRZ
    if (T_limit > 100.0_r8) T_limit=100.0_r8
    if (T_limit < -75.0_r8) T_limit=-75.0_r8

    td       = T_limit
    if (td >= 0.0_r8) then
       es   = a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 &
            + td*(a5 + td*(a6 + td*(a7 + td*a8)))))))
!       esdT = b0 + td*(b1 + td*(b2 + td*(b3 + td*(b4 &
!            + td*(b5 + td*(b6 + td*(b7 + td*b8)))))))
    else
       es   = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 &
            + td*(c5 + td*(c6 + td*(c7 + td*c8)))))))
!       esdT = d0 + td*(d1 + td*(d2 + td*(d3 + td*(d4 &
!            + td*(d5 + td*(d6 + td*(d7 + td*d8)))))))
    endif

    es    = es    * 100._r8            ! pa
!    esdT  = esdT  * 100._r8            ! pa/K

    vp    = 1.0_r8   / (p - 0.378_r8*es)
    vp1   = 0.622_r8 * vp
!    vp2   = vp1   * vp

    qs_Sat    = es    * vp1             ! kg/kg
!    qsdT  = esdT  * vp2 * p         ! 1 / K

  END FUNCTION qs_Sat

      REAL(KIND=r8) FUNCTION qsdT_Sat (T, p)
!
! !DESCRIPTION:
! Computes saturation mixing ratio and the change in saturation
! mixing ratio with respect to temperature.
! Reference:  Polynomial approximations from:
!             Piotr J. Flatau, et al.,1992:  Polynomial fits to saturation
!             vapor pressure.  Journal of Applied Meteorology, 31, 1507-1513.
!
! !USES:
!    use shr_kind_mod , only: r8 => shr_kind_r8
!    use shr_const_mod, only: SHR_CONST_TKFRZ
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: T        ! temperature (K)
    real(r8), intent(in)  :: p        ! surface atmospheric pressure (pa)
!    real(r8), intent(out) :: qsdT     ! d(qs)/d(T)
!
! !CALLED FROM:
! subroutine Biogeophysics1 in module Biogeophysics1Mod
! subroutine BiogeophysicsLake in module BiogeophysicsLakeMod
! subroutine CanopyFluxesMod CanopyFluxesMod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
!
!
! !LOCAL VARIABLES:
!EOP
!
    real(r8) :: T_limit
    real(r8) :: td,vp,vp1,vp2,qs,es,esdT
!
! For water vapor (temperature range 0C-100C)
!
    real(r8), parameter :: a0 =  6.11213476_r8
    real(r8), parameter :: a1 =  0.444007856_r8
    real(r8), parameter :: a2 =  0.143064234e-01_r8
    real(r8), parameter :: a3 =  0.264461437e-03_r8
    real(r8), parameter :: a4 =  0.305903558e-05_r8
    real(r8), parameter :: a5 =  0.196237241e-07_r8
    real(r8), parameter :: a6 =  0.892344772e-10_r8
    real(r8), parameter :: a7 = -0.373208410e-12_r8
    real(r8), parameter :: a8 =  0.209339997e-15_r8
!
! For derivative:water vapor
!
    real(r8), parameter :: b0 =  0.444017302_r8
    real(r8), parameter :: b1 =  0.286064092e-01_r8
    real(r8), parameter :: b2 =  0.794683137e-03_r8
    real(r8), parameter :: b3 =  0.121211669e-04_r8
    real(r8), parameter :: b4 =  0.103354611e-06_r8
    real(r8), parameter :: b5 =  0.404125005e-09_r8
    real(r8), parameter :: b6 = -0.788037859e-12_r8
    real(r8), parameter :: b7 = -0.114596802e-13_r8
    real(r8), parameter :: b8 =  0.381294516e-16_r8
!
! For ice (temperature range -75C-0C)
!
    real(r8), parameter :: c0 =  6.11123516_r8
    real(r8), parameter :: c1 =  0.503109514_r8
    real(r8), parameter :: c2 =  0.188369801e-01_r8
    real(r8), parameter :: c3 =  0.420547422e-03_r8
    real(r8), parameter :: c4 =  0.614396778e-05_r8
    real(r8), parameter :: c5 =  0.602780717e-07_r8
    real(r8), parameter :: c6 =  0.387940929e-09_r8
    real(r8), parameter :: c7 =  0.149436277e-11_r8
    real(r8), parameter :: c8 =  0.262655803e-14_r8
!
! For derivative:ice
!
    real(r8), parameter :: d0 =  0.503277922_r8
    real(r8), parameter :: d1 =  0.377289173e-01_r8
    real(r8), parameter :: d2 =  0.126801703e-02_r8
    real(r8), parameter :: d3 =  0.249468427e-04_r8
    real(r8), parameter :: d4 =  0.313703411e-06_r8
    real(r8), parameter :: d5 =  0.257180651e-08_r8
    real(r8), parameter :: d6 =  0.133268878e-10_r8
    real(r8), parameter :: d7 =  0.394116744e-13_r8
    real(r8), parameter :: d8 =  0.498070196e-16_r8
!-----------------------------------------------------------------------
    real(R8),parameter :: SHR_CONST_TKFRZ   = 273.15_r8       ! freezing T of fresh water          ~ K


    T_limit = T - SHR_CONST_TKFRZ
    if (T_limit > 100.0_r8) T_limit=100.0_r8
    if (T_limit < -75.0_r8) T_limit=-75.0_r8

    td       = T_limit
    if (td >= 0.0_r8) then
       es   = a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 &
            + td*(a5 + td*(a6 + td*(a7 + td*a8)))))))
       esdT = b0 + td*(b1 + td*(b2 + td*(b3 + td*(b4 &
            + td*(b5 + td*(b6 + td*(b7 + td*b8)))))))
    else
       es   = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 &
            + td*(c5 + td*(c6 + td*(c7 + td*c8)))))))
       esdT = d0 + td*(d1 + td*(d2 + td*(d3 + td*(d4 &
            + td*(d5 + td*(d6 + td*(d7 + td*d8)))))))
    endif

    es    = es    * 100._r8            ! pa
    esdT  = esdT  * 100._r8            ! pa/K

    vp    = 1.0_r8   / (p - 0.378_r8*es)
    vp1   = 0.622_r8 * vp
    vp2   = vp1   * vp

    qs    = es    * vp1             ! kg/kg
    qsdT_Sat  = esdT  * vp2 * p         ! 1 / K

  END FUNCTION qsdT_Sat

  ! ####### ## ######    #####     #    #         # ####### ######  #         #
  ! #          ### ##  #         #   # #   ##         # #         # #         #  ##
  ! #          # ## ##  #            ##  # #   # #         # #         #   # #
  ! #####   #  #  # ##  #           #         # #  #  # #         # ######     #
  ! #          #   # # ##  #           ####### #   # # #         # #              #
  ! #          #    ## ##  #         # #         # ### #         # #              #
  ! ####### ## ######    #####  #         # #         # ####### #              #

  !
  ! ---------------------------------------------------------------------

  !
  ! #####     ##    #####      #      ##     #####     #     ####   #    #
  ! #    #   #  #   #    #     #     #  #      #       #    #    #  ##   #
  ! #    #  #    #  #    #     #    #    #     #       #    #    #  # #  #
  ! #####   ######  #    #     #    ######     #       #    #    #  #  # #
  ! #   #   #    #  #    #     #    #    #     #       #    #    #  #   ##
  ! #    #  #    #  #####      #    #    #     #       #     ####   #    #
  !
  ! ---------------------------------------------------------------------





  !
  ! #####     ##    #####      #      ##     #####     #     ####   #    #
  ! #    #   #  #   #    #     #     #  #      #       #    #    #  ##   #
  ! #    #  #    #  #    #     #    #    #     #       #    #    #  # #  #
  ! #####   ######  #    #     #    ######     #       #    #    #  #  # #
  ! #   #   #    #  #    #     #    #    #     #       #    #    #  #   ##
  ! #    #  #    #  #####      #    #    #     #       #     ####   #    #
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE solset(npoi     , &! INTENT(IN   )
       nsol     , &! INTENT(OUT  )
       nband    , &! INTENT(IN   )
       solu     , &! INTENT(OUT  )
       sols         , &! INTENT(OUT  )
       soll         , &! INTENT(OUT  )
       solg     , &! INTENT(OUT  )
       soli     , &! INTENT(OUT  )
       scalcoefl, &! INTENT(OUT  )
       scalcoefu, &! INTENT(OUT  )
       indsol   , &! INTENT(OUT  )
       topparu  , &! INTENT(OUT  )
       topparl  , &! INTENT(OUT  )
       asurd    , &! INTENT(OUT  )
       asuri    , &! INTENT(OUT  )
       coszen     )! INTENT(IN  )  
    ! ---------------------------------------------------------------------
    !
    ! zeros albedos and internal absorbed solar fluxes, and sets
    ! index for other solar routines. the index indsol, with number
    ! of points nsol, points to current 1d strip arrays whose coszen 
    ! values are gt 0 (indsol, nsol are in com1d)
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN   ) :: npoi              ! total number of land points
    INTEGER, INTENT(OUT  ) :: nsol              ! number of points in indsol
    INTEGER, INTENT(IN   ) :: nband             ! number of solar radiation wavebands
    REAL(KIND=r8), INTENT(OUT  ) :: solu     (npoi)   ! solar flux (direct + diffuse) absorbed by upper
    ! canopy leaves per unit canopy area (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: sols     (npoi)   ! solar flux (direct + diffuse) absorbed by upper
    ! canopy stems per unit canopy area (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: soll     (npoi)   ! solar flux (direct + diffuse) absorbed by lower
    ! canopy leaves and stems per unit canopy area (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: solg     (npoi)   ! solar flux (direct + diffuse) absorbed by unit
    ! snow-free soil (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: soli     (npoi)   ! solar flux (direct + diffuse) absorbed by
    ! unit snow surface (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: scalcoefl(npoi,4) ! term needed in lower canopy scaling
    REAL(KIND=r8), INTENT(OUT  ) :: scalcoefu(npoi,4) ! term needed in upper canopy scaling
    INTEGER, INTENT(OUT  ) :: indsol   (npoi)   ! index of current strip for points with positive coszen
    REAL(KIND=r8), INTENT(OUT  ) :: topparu  (npoi)   ! total photosynthetically active raditaion
    ! absorbed by top leaves of upper canopy (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: topparl  (npoi)   ! total photosynthetically active raditaion absorbed
    ! by top leaves of lower canopy (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: asurd    (npoi,nband) ! direct albedo of surface system
    REAL(KIND=r8), INTENT(OUT  ) :: asuri    (npoi,nband) ! diffuse albedo of surface system 
    REAL(KIND=r8), INTENT(IN   ) :: coszen   (npoi)      ! cosine of solar zenith angle
    !
    INTEGER :: i,k
    !
    ! zero albedos returned just as a niceity
    !
    DO k=1,nband
       DO i=1, npoi
          asurd    (i,k) = 0.0_r8  
          asuri    (i,k) = 0.0_r8  
       END DO
    END DO

    !      CALL const2 (asurd     , & ! INTENT(OUT  ) ::  arr(nar)
    !                   npoi*nband, & ! INTENT(IN   ) ::  nar
    !           0.0_r8         ) ! INTENT(IN   ) ::  value
    !      CALL const2 (asuri     , & ! INTENT(OUT  ) ::  arr(nar)
    !                   npoi*nband, & ! INTENT(IN   ) ::  nar
    !           0.0_r8         ) ! INTENT(IN   ) ::  value
    !
    ! zeros absorbed solar fluxes sol[u,s,l,g,i]1 since only points
    ! with +ve coszen will be set in solarf, and since
    ! sol[u,l,s,g,i]1 are summed over wavebands in solarf
    !
    ! similarly zero par-related arrays set in solarf for turvap
    !
    DO i=1, npoi
       solu    (i) = 0.0_r8  
       sols    (i) = 0.0_r8  
       soll    (i) = 0.0_r8
       solg    (i) = 0.0_r8
       soli    (i) = 0.0_r8
       topparu (i) = 0.0_r8
       topparl (i) = 0.0_r8
    END DO

    !      CALL const2 (solu      , &! INTENT(OUT  ) :: arr(nar)
    !                   npoi      , &! INTENT(IN   ) :: nar
    !           0.0_r8         )! INTENT(IN   ) :: value
    !      CALL const2 (sols      , &! INTENT(OUT  ) :: arr(nar)
    !                   npoi      , &! INTENT(IN   ) :: nar
    !           0.0_r8         )! INTENT(IN   ) :: value
    !      CALL const2 (soll      , &! INTENT(OUT  ) :: arr(nar)
    !                   npoi      , &! INTENT(IN   ) :: nar
    !           0.0_r8         )! INTENT(IN   ) :: value
    !      CALL const2 (solg      , &! INTENT(OUT  ) :: arr(nar)
    !                   npoi      , &! INTENT(IN   ) :: nar
    !           0.0_r8         )! INTENT(IN   ) :: value
    !      CALL const2 (soli      , &! INTENT(OUT  ) :: arr(nar)
    !                   npoi      , &! INTENT(IN   ) :: nar
    !           0.0_r8         )! INTENT(IN   ) :: value
    !
    !      CALL const2 (topparu   , &! INTENT(OUT  ) :: arr(nar)
    !                   npoi      , &! INTENT(IN   ) :: nar
    !           0.0_r8         )! INTENT(IN   ) :: value
    !      CALL const2 (topparl   , &! INTENT(OUT  ) :: arr(nar)
    !                   npoi      , &! INTENT(IN   ) :: nar
    !           0.0_r8         )! INTENT(IN   ) :: value
    !
    ! set canopy scaling coefficients for night-time conditions
    !
    DO k=1,4
       DO i=1, npoi
          scalcoefl    (i,k) = 0.0_r8  
          scalcoefu    (i,k) = 0.0_r8  
       END DO
    END DO

    !      CALL const2 (scalcoefl , &! INTENT(OUT  ) :: arr(nar)
    !                   npoi*4    , &! INTENT(IN   ) :: nar
    !           0.0_r8         )! INTENT(IN   ) :: value
    !      CALL const2 (scalcoefu , &! INTENT(OUT  ) :: arr(nar)
    !                   npoi*4    , &! INTENT(IN   ) :: nar
    !           0.0_r8         )! INTENT(IN   ) :: value
    !
    ! set index of points with positive coszen
    !
    nsol = 0
    !
    DO i = 1, npoi
       IF (coszen(i).GT.0.0_r8) THEN
          nsol = nsol + 1
          indsol(nsol) = i
       END IF
    END DO
    !
    RETURN
  END SUBROUTINE solset
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE solsur (ib     , &! INTENT(IN   )
       tmelt  , &! INTENT(IN   )
       nsol   , &! INTENT(IN   )
       albsod , &! INTENT(OUt  )
       albsoi , &! INTENT(OUt  )
       albsnd , &! INTENT(OUt  )
       albsni , &! INTENT(OUt  )
       indsol , &! INTENT(IN   )
       wsoi   , &! INTENT(IN   )
       wisoi  , &! INTENT(IN   )
       albsav , &! INTENT(IN   )
       albsan , &! INTENT(IN   )
       tsno   , &! INTENT(IN   )
       coszen , &! INTENT(IN   )
       npoi   , &! INTENT(IN   )
       nsoilay, &! INTENT(IN   )
       nsnolay  )! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    ! sets surface albedos for soil and snow, prior to other
    ! solar calculations
    !
    ! ib = waveband number
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: npoi
    INTEGER, INTENT(IN   ) :: nsoilay
    INTEGER, INTENT(IN   ) :: nsnolay
    INTEGER, INTENT(IN   ) :: nsol 
    REAL(KIND=r8), INTENT(IN   ) :: tmelt
    REAL(KIND=r8), INTENT(OUt  ) :: albsod(npoi) ! direct  albedo for soil surface (visible or IR)
    REAL(KIND=r8), INTENT(OUt  ) :: albsoi(npoi) ! diffuse albedo for soil surface (visible or IR)
    REAL(KIND=r8), INTENT(OUt  ) :: albsnd(npoi) ! direct  albedo for snow surface (visible or IR)
    REAL(KIND=r8), INTENT(OUt  ) :: albsni(npoi) ! diffuse albedo for snow surface (visible or IR)
    INTEGER, INTENT(IN   ) :: indsol(npoi)        ! index of current strip for points with positive coszen
    REAL(KIND=r8), INTENT(IN   ) :: wsoi  (npoi,nsoilay)! fraction of soil pore space containing liquid water
    REAL(KIND=r8), INTENT(IN   ) :: wisoi (npoi,nsoilay)! fraction of soil pore space containing ice
    REAL(KIND=r8), INTENT(IN   ) :: albsav(npoi)          ! saturated soil surface albedo (visible waveband)
    REAL(KIND=r8), INTENT(IN   ) :: albsan(npoi)        ! saturated soil surface albedo (near-ir waveband)
    REAL(KIND=r8), INTENT(IN   ) :: tsno  (npoi,nsnolay)  ! temperature of snow layers (K)
    REAL(KIND=r8), INTENT(IN   ) :: coszen(npoi)        ! cosine of solar zenith angle

    !
    ! input variable
    !
    INTEGER, INTENT(IN   ) :: ib    ! waveband number. 1 = visible, 2 = near IR
    !
    ! local variables
    !     
    INTEGER j     ! loop indice on number of points with >0 coszen
    INTEGER i     ! indice of point in (1, npoi) array. 
    !
    REAL(KIND=r8), PARAMETER :: a7svlo=0.90_r8! snow albedo at low threshold temp., visible
    REAL(KIND=r8), PARAMETER :: a7svhi=0.70_r8! high              , visible
    REAL(KIND=r8), PARAMETER :: a7snlo=0.60_r8!                                   , near IR
    REAL(KIND=r8), PARAMETER :: a7snhi=0.40_r8!                                   , near-IR
    !    avisb(i)=icealv!icealv =  0.8e0_r8! constant icealv
    !    avisd(i)=icealv!icealv =  0.8e0_r8! constant icealv
    !    anirb(i)=icealn!icealn =  0.4e0_r8! constant icealn
    !    anird(i)=icealn!icealn =  0.4e0_r8! constant icealn
    REAL(KIND=r8) t7shi ! high threshold temperature for snow albed
    REAL(KIND=r8) t7slo ! low  threshold temperature for snow albedo
    REAL(KIND=r8) dinc  ! albedo correction du to soil moisture
    REAL(KIND=r8) zw    ! liquid moisture content

    REAL(KIND=r8) x   (npoi) 
    REAL(KIND=r8) zfac(npoi)
    !
    ! set the "standard" snow values:
    !
    !      DATA    a7svlo, a7svhi /0.90_r8, 0.70_r8/
    !      DATA    a7snlo, a7snhi /0.60_r8, 0.40_r8/
    !
    !     t7shi ... high threshold temperature for snow albedo
    !     t7slo ... low  threshold temperature for snow albedo
    !
    t7shi = tmelt
    t7slo = tmelt - 15.0_r8
    !
    ! do nothing if all points in current strip have coszen le 0
    !
    IF (nsol.EQ.0) THEN
       RETURN
    END IF
    !
    IF (ib.EQ.1) THEN
       !
       ! soil albedos (visible waveband)
       !
       DO  j = 1, nsol
          !
          i = indsol(j)
          !
          ! change the soil albedo as a function of soil moisture
          !
          zw = wsoi(i,1) * (1.0_r8-wisoi(i,1))
          !
          dinc = 1.0_r8 + 1.0_r8 * MIN (1.0_r8, MAX (0.0_r8, 1.0_r8 - (zw /.50_r8) ))
          !
          albsod(i) = MIN (albsav(i) * dinc, 0.80_r8)
          albsoi(i) = albsod(i)
          !
       END DO
       !
       ! snow albedos (visible waveband)
       !
       DO  j = 1, nsol
          !
          i = indsol(j)
          !
          x(i) = (a7svhi*(tsno(i,1)-t7slo) + a7svlo*(t7shi-tsno(i,1)))   &
               / (t7shi-t7slo)
          !
          x(i) = MIN (a7svlo, MAX (a7svhi, x(i)))
          !
          zfac(i)   = MAX ( 0.0_r8, 1.5_r8 / (1.0_r8 + 4.0_r8*coszen(i)) - 0.5_r8 )
          albsnd(i) = MIN (0.99_r8, x(i) + (1.0_r8-x(i))*zfac(i))
          albsni(i) = MIN (1.0_r8, x(i))
          !
       END DO
       !
    ELSE
       !
       ! soil albedos (near-ir waveband)
       !
       DO  j = 1, nsol
          i = indsol(j)
          !
          ! lsx.2 formulation (different from lsx.1)
          !
          zw = wsoi(i,1) * (1.0_r8 - wisoi(i,1))
          !
          dinc = 1.0_r8 + 1.0_r8 * MIN (1.0_r8, MAX (0.0_r8, 1.0_r8 - (zw / .50_r8)  ))
          !
          albsod(i) = MIN (albsan(i) * dinc, .80_r8)
          albsoi(i) = albsod(i)
          !
       END DO
       !
       ! snow albedos (near-ir waveband)
       !
       DO  j = 1, nsol
          !
          i = indsol(j)
          !
          x(i) = (a7snhi*(tsno(i,1)-t7slo) + a7snlo*(t7shi-tsno(i,1)))  &
               / (t7shi-t7slo)
          x(i) = MIN (a7snlo, MAX (a7snhi, x(i)))
          !
          zfac(i) = MAX ( 0.0_r8, 1.5_r8/(1.0_r8+4.0_r8*coszen(i)) - 0.5_r8 )
          !
          albsnd(i) = MIN (0.99_r8, x(i) + (1.0_r8-x(i))*zfac(i))
          albsni(i) = MIN (1.0_r8, x(i))
          !
       END DO
       !
    END IF
    !
    RETURN
  END SUBROUTINE solsur
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE solalb (ib     , &! INTENT(IN   )
       nVegClass , &! INTENT(IN   )
       vegtype0    , &! INTENT(IN   )   ! annual vegetation type - ibis classification
       avmuir_factor , &! INTENT(IN   )
       relod  , &! INTENT(OUT  )
       reloi  , &! INTENT(OUT  )
       indsol , &! INTENT(IN   )
       reupd  , &! INTENT(OUT  )
       reupi  , &! INTENT(OUT  )
       albsnd , &! INTENT(IN   )
       albsni , &! INTENT(IN   )
       albsod , &! INTENT(IN   )
       albsoi , &! INTENT(IN   )
       fl     , &! INTENT(IN   )
       fu     , &! INTENT(IN   )
       fi     , &! INTENT(IN   )
       asurd  , &! INTENT(INOUT)! local
       asuri  , &! INTENT(INOUT)! local
       npoi   , &! INTENT(IN   )
       nband  , &! INTENT(IN   )
       nsol   , &! INTENT(IN   )
       ablod  , &! INTENT(OUT  )
       abloi  , &! INTENT(OUT  )
       flodd  , &! INTENT(OUT  )
       dummy  , &! INTENT(OUT  )
       flodi  , &! INTENT(OUT  )
       floii  , &! INTENT(OUT  )
       coszen , &! INTENT(IN   )
       terml  , &! INTENT(OUT  )
       termu  , &! INTENT(OUT  )
       lai    , &! INTENT(IN   )
       sai    , &! INTENT(IN   )
       abupd  , &! INTENT(OUT  )
       abupi  , &! INTENT(OUT  )
       fupdd  , &! INTENT(OUT  )
       fupdi  , &! INTENT(OUT  )
       fupii  , &! INTENT(OUT  )
       fwetl  , &! INTENT(IN   )
       rliql  , &! INTENT(IN   )
       rliqu  , &! INTENT(IN   )
       rliqs  , &! INTENT(IN   )
       fwetu  , &! INTENT(IN   )
       fwets  , &! INTENT(IN   )
       rhoveg , &! INTENT(IN   )
       tauveg , &! INTENT(IN   )
       orieh  , &! INTENT(IN   )
       oriev  , &! INTENT(IN   )
       tl     , &! INTENT(IN   )
       ts     , &! INTENT(IN   )
       tu     , &! INTENT(IN   )
       pi     , &! INTENT(IN   )
       tmelt  , &! INTENT(IN   )
       epsilon  )! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    ! calculates effective albedos of the surface system,
    ! separately for unit incoming direct and diffuse flux -- the 
    ! incoming direct zenith angles are supplied in comatm array 
    ! coszen, and the effective albedos are returned in comatm
    ! arrays asurd, asuri -- also detailed absorbed and reflected flux
    ! info is stored in com1d arrays, for later use by solarf
    !
    ! the procedure is first to calculate the grass+soil albedos,
    ! then the tree + (grass+soil+snow) albedos. the labels
    ! (a) to (d) correspond to those in the description doc
    !
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: npoi
    INTEGER, INTENT(IN   ) :: nsol
    INTEGER, INTENT(IN   ) :: nband
    INTEGER, INTENT(IN   ) :: nVegClass ! INTENT(IN   )
    REAL(KIND=r8), INTENT(IN   ) :: vegtype0 (npoi)   ! INTENT(IN   )   ! annual vegetation type - ibis classification
    REAL(KIND=r8), INTENT(IN   ) :: avmuir_factor(nVegClass,2 )! INTENT(IN   )

    REAL(KIND=r8), INTENT(IN   ) :: fwetl (npoi)     ! fraction of lower canopy stem & leaf area wetted
    ! by intercepted liquid and/or snow
    REAL(KIND=r8), INTENT(IN   ) :: rliql (npoi)     ! proportion of fwetl due to liquid
    REAL(KIND=r8), INTENT(IN   ) :: rliqu (npoi)     ! proportion of fwetu due to liquid 
    REAL(KIND=r8), INTENT(IN   ) :: rliqs (npoi)     ! proportion of fwets due to liquid
    REAL(KIND=r8), INTENT(IN   ) :: fwetu (npoi)     ! fraction of upper canopy leaf area wetted by
    ! intercepted liquid and/or snow
    REAL(KIND=r8), INTENT(IN   ) :: fwets (npoi)     ! fraction of upper canopy stem area wetted by
    ! intercepted liquid and/or snow
    REAL(KIND=r8), INTENT(IN   ) :: rhoveg(nband,2)  ! reflectance of an average leaf/stem
    REAL(KIND=r8), INTENT(IN   ) :: tauveg(nband,2)  ! transmittance of an average leaf/stem
    REAL(KIND=r8), INTENT(IN   ) :: orieh (2)        ! fraction of leaf/stems with horizontal orientation
    REAL(KIND=r8), INTENT(IN   ) :: oriev (2)        ! fraction of leaf/stems with vertical
    REAL(KIND=r8), INTENT(IN   ) :: tl    (npoi)     ! temperature of lower canopy leaves & stems(K)
    REAL(KIND=r8), INTENT(IN   ) :: ts    (npoi)     ! temperature of upper canopy stems (K)
    REAL(KIND=r8), INTENT(IN   ) :: tu    (npoi)     ! temperature of upper canopy leaves (K)
    REAL(KIND=r8), INTENT(IN   ) :: pi        
    REAL(KIND=r8), INTENT(IN   ) :: tmelt            ! freezing point of water (K) 
    REAL(KIND=r8), INTENT(IN   ) :: epsilon          ! small quantity to avoid zero-divides and other
    ! truncation or machine-limit troubles with small
    ! values. should be slightly greater than o(1)
    ! machine precision


    REAL(KIND=r8), INTENT(INOUT) :: asurd(npoi,nband)! direct albedo of surface system
    REAL(KIND=r8), INTENT(INOUT) :: asuri(npoi,nband)! diffuse albedo of surface system 
    REAL(KIND=r8), INTENT(OUT  ) :: abupd(npoi)      ! fraction of direct  radiation absorbed by upper canopy
    REAL(KIND=r8), INTENT(OUT  ) :: abupi(npoi)      ! fraction of diffuse radiation absorbed by upper canopy
    REAL(KIND=r8), INTENT(OUT  ) :: fupdd(npoi)      ! downward direct radiation per unit incident direct
    ! beam on upper canopy (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: fupdi(npoi)      ! downward diffuse radiation per unit icident direct
    ! radiation on upper canopy (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: fupii(npoi)      ! downward diffuse radiation per unit incident diffuse
    ! radiation on upper canopy (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: fi(npoi)         ! fractional snow cover
    REAL(KIND=r8), INTENT(IN   ) :: fl(npoi)         ! fraction of snow-free area covered by lower  canopy
    REAL(KIND=r8), INTENT(IN   ) :: fu(npoi)         ! fraction of overall area covered by upper canopy
    REAL(KIND=r8), INTENT(IN   ) :: lai(npoi,2)      ! canopy single-sided leaf area index (area leaf/area veg)
    REAL(KIND=r8), INTENT(IN   ) :: sai(npoi,2)      ! current single-sided stem area index
    REAL(KIND=r8), INTENT(OUT  ) :: relod (npoi)     ! upward direct radiation per unit icident direct beam
    ! on lower canopy (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: reloi (npoi)     ! upward diffuse radiation per unit incident diffuse
    ! radiation on lower canopy (W m-2)
    INTEGER, INTENT(IN   ) :: indsol(npoi)     ! index of current strip for points with positive coszen
    REAL(KIND=r8), INTENT(OUT  ) :: reupd (npoi)     ! upward direct radiation per unit incident direct
    ! radiation on upper canopy (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: reupi (npoi)     ! upward diffuse radiation per unit incident diffuse
    ! radiation on upper canopy (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: albsnd(npoi)     ! direct  albedo for snow surface (visible or IR)
    REAL(KIND=r8), INTENT(IN   ) :: albsni(npoi)     ! diffuse albedo for snow surface (visible or IR)
    REAL(KIND=r8), INTENT(IN   ) :: albsod(npoi)     ! direct  albedo for soil surface (visible or IR)
    REAL(KIND=r8), INTENT(IN   ) :: albsoi(npoi)     ! diffuse albedo for soil surface (visible or IR)
    REAL(KIND=r8), INTENT(OUT  ) :: ablod (npoi)     ! fraction of direct  radiation absorbed by lower canopy
    REAL(KIND=r8), INTENT(OUT  ) :: abloi (npoi)     ! fraction of diffuse radiation absorbed by lower canopy
    REAL(KIND=r8), INTENT(OUT  ) :: flodd (npoi)     ! downward direct radiation per unit incident direct
    ! radiation on lower canopy (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: dummy (npoi)     ! placeholder, always = 0: no direct flux produced for diffuse incident
    REAL(KIND=r8), INTENT(OUT  ) :: flodi (npoi)     ! downward diffuse radiation per unit incident direct
    ! radiation on lower canopy (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: floii (npoi)     ! downward diffuse radiation per unit incident diffuse
    ! radiation on lower canopy
    REAL(KIND=r8), INTENT(IN   ) :: coszen(npoi)     ! cosine of solar zenith angle
    REAL(KIND=r8), INTENT(OUT  ) :: terml (npoi,7)   ! term needed in lower canopy scaling
    REAL(KIND=r8), INTENT(OUT  ) :: termu (npoi,7)   ! term needed in upper canopy scaling
    ! 
    ! Arguments
    ! 
    INTEGER, INTENT(IN   ) :: ib     ! waveband number (1= visible, 2= near-IR)
    !
    ! local variables
    !     
    INTEGER :: j      ! loop indice on number of points with >0 coszen
    INTEGER :: i      ! indice of point in (1, npoi) array. 
    !
    ! do nothing if all points in current strip have coszen le 0
    !
    IF (nsol.EQ.0) RETURN
    !
    ! (a) obtain albedos, etc, for two-stream lower veg + soil
    !     system, for direct and diffuse incoming unit flux
    !
    DO  j = 1, nsol
       !
       i = indsol(j)
       !
       asurd(i,ib) = albsod(i)
       asuri(i,ib) = albsoi(i)
       !
    END DO
    !

    CALL twostr (       &! INTENT(IN   )
         nVegClass    , &! INTENT(IN   )
         vegtype0     , &! INTENT(IN   )   ! annual vegetation type - ibis classification
         avmuir_factor,&! INTENT(IN   )
         ablod , &! INTENT(OUT  )
         abloi , &! INTENT(OUT  )
         relod , &! INTENT(OUT  )
         reloi , &! INTENT(OUT  )
         flodd , &! INTENT(OUT  )
         dummy , &! INTENT(OUT  )
         flodi , &! INTENT(OUT  )
         floii , &! INTENT(OUT  )
         asurd , &! INTENT(IN   )
         asuri , &! INTENT(IN   )
         1     , &! INTENT(IN   )
         coszen, &! INTENT(IN   )
         ib    , &! INTENT(IN   )
         indsol, &! INTENT(IN   )
         terml , &! INTENT(OUT  )
         termu , &! INTENT(OUT  )
         lai   , &! INTENT(IN   )
         sai   , &! INTENT(IN   )
         npoi  , &! INTENT(IN   )
         nband , &! INTENT(IN   )
         nsol  , &! INTENT(IN   )
         fwetl , &! INTENT(IN   )
         rliql , &! INTENT(IN   )
         rliqu , &! INTENT(IN   )
         rliqs , &! INTENT(IN   )
         fwetu , &! INTENT(IN   )
         fwets , &! INTENT(IN   )
         rhoveg, &! INTENT(IN   )
         tauveg, &! INTENT(IN   )
         orieh , &! INTENT(IN   )
         oriev , &! INTENT(IN   )
         tl    , &! INTENT(IN   )
         ts    , &! INTENT(IN   )
         tu    , &! INTENT(IN   )
         pi    , &! INTENT(IN   )
         tmelt , &! INTENT(IN   )
         epsilon )! INTENT(IN   )
    !
    ! (b) areally average surface albedos (lower veg, soil, snow)
    !
    DO  j = 1, nsol
       !
       i = indsol(j)
       !
       asurd(i,ib) = fl(i)*(1.0_r8-fi(i))*relod(i)       &
            + (1.0_r8-fl(i))*(1.0_r8-fi(i))*albsod(i)  &
            + fi(i)*albsnd(i)
       !
       asuri(i,ib) = fl(i)*(1.0_r8-fi(i))*reloi(i)        &
            + (1.0_r8-fl(i))*(1.0_r8-fi(i))*albsoi(i)   &
            + fi(i)*albsni(i)
       !
    END DO
    !
    ! (c) obtain albedos, etc, for two-stream upper veg + surface
    !     system, for direct and diffuse incoming unit flux
    !
    CALL twostr (       &! INTENT(IN   )
         nVegClass    , &! INTENT(IN   )
         vegtype0     , &! INTENT(IN   )   ! annual vegetation type - ibis classification
         avmuir_factor,&! INTENT(IN   )
         abupd  , &! INTENT(OUT  )
         abupi  , &! INTENT(OUT  )
         reupd  , &! INTENT(OUT  )
         reupi  , &! INTENT(OUT  )
         fupdd  , &! INTENT(OUT  )
         dummy  , &! INTENT(OUT  )
         fupdi  , &! INTENT(OUT  )
         fupii  , &! INTENT(OUT  )
         asurd  , &! INTENT(IN   )
         asuri  , &! INTENT(IN   )
         2      , &! INTENT(IN   )
         coszen , &! INTENT(IN   )
         ib     , &! INTENT(IN   )
         indsol , &! INTENT(IN   )
         terml  , &! INTENT(OUT  )
         termu  , &! INTENT(OUT  )
         lai    , &! INTENT(IN   )
         sai    , &! INTENT(IN   )
         npoi   , &! INTENT(IN   )
         nband  , &! INTENT(IN   )
         nsol   , &! INTENT(IN   )
         fwetl  , &! INTENT(IN   )
         rliql  , &! INTENT(IN   )
         rliqu  , &! INTENT(IN   )
         rliqs  , &! INTENT(IN   )
         fwetu  , &! INTENT(IN   )
         fwets  , &! INTENT(IN   )
         rhoveg , &! INTENT(IN   )
         tauveg , &! INTENT(IN   )
         orieh  , &! INTENT(IN   )
         oriev  , &! INTENT(IN   )
         tl     , &! INTENT(IN   )
         ts     , &! INTENT(IN   )
         tu     , &! INTENT(IN   )
         pi     , &! INTENT(IN   )
         tmelt  , &! INTENT(IN   )
         epsilon  )! INTENT(IN   )
    !
    ! (d) calculate average overall albedos 
    !
    DO  j = 1, nsol
       !
       i = indsol(j)
       !  number of solar radiation wavebands : vis, nir
       !  REAL (KIND=r8), PARAMETER   :: icealv =         0.8e0_r8! constant icealv
       !  REAL (KIND=r8), PARAMETER   :: icealn =         0.4e0_r8! constant icealn
       !         avisb(i)=asuri(ncount,1,jb)           !asurd  (npoi,nband)   ! local  ! direct albedo of surface system
       !         avisd(i)=asurd(ncount,1,jb)           !asuri  (npoi,nband)   ! local  ! diffuse albedo of surface system 
       !         anirb(i)=asuri(ncount,2,jb)
       !         anird(i)=asurd(ncount,2,jb)
       asurd(i,ib) = fu(i)*reupd(i) + (1.0_r8-fu(i))*asurd(i,ib)
       asuri(i,ib) = fu(i)*reupi(i) + (1.0_r8-fu(i))*asuri(i,ib)

       IF(ib == 1) THEN
          !vis
          asurd(i,ib) = MIN(0.8e0_r8,MAX(0.00_r8,asurd(i,ib)))
          asuri(i,ib) = MIN(0.8e0_r8,MAX(0.00_r8,asuri(i,ib)))
       ELSE
          !nir
          asurd(i,ib) = MIN(0.4e0_r8,MAX(0.00_r8,asurd(i,ib)))
          asuri(i,ib) = MIN(0.4e0_r8,MAX(0.00_r8,asuri(i,ib)))
       END IF
       !
    END DO
    !
    RETURN
  END SUBROUTINE solalb
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE solarf (ib       , &! INTENT(IN   ) 
       nsol     , &! INTENT(IN   ) 
       solu     , &! INTENT(INOUT) !global
       indsol   , &! INTENT(IN   ) 
       abupd    , &! INTENT(IN   ) 
       abupi    , &! INTENT(IN   ) 
       sols          , &! INTENT(INOUT) !global
       sol2d    , &! INTENT(OUT  ) 
       fupdd    , &! INTENT(IN   ) 
       sol2i    , &! INTENT(OUT  ) 
       fupii    , &! INTENT(IN   ) 
       fupdi    , &! INTENT(IN   ) 
       sol3d    , &! INTENT(OUT  ) 
       sol3i    , &! INTENT(OUT  ) 
       soll     , &! INTENT(INOUT) !global
       ablod    , &! INTENT(IN   ) 
       abloi          , &! INTENT(IN   ) 
       flodd    , &! INTENT(IN   ) 
       flodi    , &! INTENT(IN   ) 
       floii    , &! INTENT(IN   ) 
       solg     , &! INTENT(INOUT) !global
       albsod   , &! INTENT(IN   ) 
       albsoi   , &! INTENT(IN   ) 
       soli     , &! INTENT(INOUT) !global
       albsnd   , &! INTENT(IN   ) 
       albsni   , &! INTENT(IN   ) 
       scalcoefu, &! INTENT(OUT  ) 
       termu    , &! INTENT(IN   ) 
       scalcoefl, &! INTENT(OUT  ) 
       terml    , &! INTENT(IN   ) 
       lai      , &! INTENT(IN   ) 
       sai          , &! INTENT(IN   ) 
       fu          , &! INTENT(IN   )
       fl          , &! INTENT(IN   )
       topparu  , &! INTENT(OUT  ) 
       topparl  , &! INTENT(OUT  ) 
       solad    , &! INTENT(IN   )
       solai          , &! INTENT(IN   )
       npoi     , &! INTENT(IN   ) 
       nband    , &! INTENT(IN   ) 
       epsilon    )! INTENT(IN   ) 
    ! ---------------------------------------------------------------------
    !
    ! calculates solar fluxes absorbed by upper and lower stories,
    ! soil and snow
    !
    ! zenith angles are in comatm array coszen, and must be the same
    ! as supplied earlier to solalb
    !
    ! solarf uses the results obtained earlier by solalb and 
    ! stored in com1d arrays. the absorbed fluxes are returned in
    ! com1d arrays sol[u,s,l,g,i]
    !
    ! the procedure is first to calculate the upper-story absorbed
    ! fluxes and fluxes below the upper story, then the lower-story
    ! absorbed fluxes and fluxes below the lower story, then fluxes
    ! absorbed by the soil and snow
    !
    ! ib = waveband number
    !
    IMPLICIT NONE 
    ! 
    ! Arguments
    ! 
    INTEGER, INTENT(IN  ) ::  ib     ! waveband number (1= visible, 2= near-IR)
    INTEGER, INTENT(IN  ) ::  nsol   ! number of points in indsol
    INTEGER, INTENT(IN  ) ::  npoi   ! total number of land points
    INTEGER, INTENT(IN  ) ::  nband  ! number of solar radiation wavebands
    REAL(KIND=r8), INTENT(IN  ) ::  epsilon! small quantity to avoid zero-divides and other
    ! truncation or machine-limit troubles with small
    ! values. should be slightly greater than o(1)
    ! machine precision      
    REAL(KIND=r8), INTENT(INOUT) :: solu  (npoi) ! solar flux (direct + diffuse) absorbed by upper canopy
    ! leaves per unit canopy area (W m-2)
    INTEGER, INTENT(IN   ) :: indsol(npoi) ! index of current strip for points with positive coszen
    REAL(KIND=r8), INTENT(IN   ) :: abupd (npoi) ! fraction of direct  radiation absorbed by upper canopy
    REAL(KIND=r8), INTENT(IN   ) :: abupi (npoi) ! fraction of diffuse radiation absorbed by upper canopy
    REAL(KIND=r8), INTENT(INOUT) :: sols  (npoi) ! solar flux (direct + diffuse) absorbed by upper canopy
    ! stems per unit canopy area (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: sol2d (npoi) ! direct downward radiation  out of upper canopy per unit
    ! vegetated (upper) area (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: fupdd (npoi) ! downward direct radiation per unit incident direct beam
    ! on upper canopy (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: sol2i (npoi) ! diffuse downward radiation out of upper canopy per unit
    ! vegetated (upper) area(W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: fupii (npoi) ! downward diffuse radiation per unit incident diffuse
    ! radiation on upper canopy (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: fupdi (npoi) ! downward diffuse radiation per unit icident direct radiation
    ! on upper canopy (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: sol3d (npoi) ! direct downward radiation  out of upper canopy + gaps per
    ! unit grid cell area (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: sol3i (npoi) ! diffuse downward radiation out of upper canopy + gaps per
    ! unit grid cell area (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: soll  (npoi) ! solar flux (direct + diffuse) absorbed by lower canopy
    ! leaves and stems per unit canopy area (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: ablod (npoi) ! fraction of direct  radiation absorbed by lower canopy
    REAL(KIND=r8), INTENT(IN   ) :: abloi (npoi) ! fraction of diffuse radiation absorbed by lower canopy
    REAL(KIND=r8), INTENT(IN   ) :: flodd (npoi) ! downward direct radiation per unit incident direct radiation
    ! on lower canopy (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: flodi (npoi) ! downward diffuse radiation per unit incident direct radiation
    ! on lower canopy (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: floii (npoi) ! downward diffuse radiation per unit incident diffuse radiation
    ! on lower canopy
    REAL(KIND=r8), INTENT(INOUT) :: solg  (npoi) ! solar flux (direct + diffuse) absorbed by unit snow-free soil (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: albsod(npoi) ! direct  albedo for soil surface (visible or IR)
    REAL(KIND=r8), INTENT(IN   ) :: albsoi(npoi) ! diffuse albedo for soil surface (visible or IR)
    REAL(KIND=r8), INTENT(INOUT) :: soli  (npoi) ! solar flux (direct + diffuse) absorbed by unit snow surface (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: albsnd(npoi) ! direct  albedo for snow surface (visible or IR)
    REAL(KIND=r8), INTENT(IN   ) :: albsni(npoi) ! diffuse albedo for snow surface (visible or IR)
    REAL(KIND=r8), INTENT(OUT  ) :: scalcoefu (npoi,4)! term needed in upper canopy scaling
    REAL(KIND=r8), INTENT(IN   ) :: termu     (npoi,7)! term needed in upper canopy scaling
    REAL(KIND=r8), INTENT(OUT  ) :: scalcoefl (npoi,4)! term needed in lower canopy scaling
    REAL(KIND=r8), INTENT(IN   ) :: terml     (npoi,7)! term needed in lower canopy scaling
    REAL(KIND=r8), INTENT(IN   ) :: lai       (npoi,2)! canopy single-sided leaf area index (area leaf/area veg)
    REAL(KIND=r8), INTENT(IN   ) :: sai       (npoi,2)! current single-sided stem area index
    REAL(KIND=r8), INTENT(IN   ) :: fu        (npoi)  ! fraction of overall area covered by upper canopy
    REAL(KIND=r8), INTENT(IN   ) :: fl     (npoi) ! fraction of snow-free area covered by lower  canopy
    REAL(KIND=r8), INTENT(OUT  ) ::  topparu(npoi) ! total photosynthetically active raditaion absorbed
    ! by top leaves of upper canopy (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) ::  topparl(npoi) ! total photosynthetically active raditaion absorbed by
    ! top leaves of lower canopy (W m-2)
    REAL(KIND=r8), INTENT(IN   ) ::  solad(npoi,nband) ! direct downward solar flux (W m-2)
    REAL(KIND=r8), INTENT(IN   ) ::  solai(npoi,nband) ! diffuse downward solar flux (W m-2)
    !
    ! local variables
    !     
    INTEGER :: j     ! loop indice on number of points with >0 coszen
    INTEGER :: i     ! indice of point in (1, npoi) array. 
    !
    REAL(KIND=r8) :: x 
    REAL(KIND=r8) :: y 
    REAL(KIND=r8) :: xd 
    REAL(KIND=r8) :: xi 
    REAL(KIND=r8) :: xaiu! total single-sided lai+sai, upper
    REAL(KIND=r8) :: xail! total single-sided lai+sai, lower
    !
    ! do nothing if all points in current strip have coszen le 0
    !
    IF (nsol.EQ.0) RETURN
    !
    ! (f) calculate fluxes absorbed by upper leaves and stems,
    !     and downward fluxes below upper veg, using unit-flux
    !     results of solalb(c) (apportion absorbed flux between
    !     leaves and stems in proportion to their lai and sai)
    !
    DO j=1,nsol
       !
       i = indsol(j)
       x = solad(i,ib)*abupd(i) + solai(i,ib)*abupi(i)
       y = lai(i,2) / MAX (lai(i,2)+sai(i,2), epsilon)
       solu(i)  = solu (i) + x * y
       sols(i)  = sols (i) + x * (1.0_r8-y)
       sol2d(i) = solad(i,ib)*fupdd(i)
       sol2i(i) = solad(i,ib)*fupdi(i) + solai(i,ib)*fupii(i)
       !
    END DO
    !
    ! (g) areally average fluxes to lower veg, soil, snow
    !
    DO j=1,nsol
       !
       i = indsol(j)
       sol3d(i) = fu(i)*sol2d(i) + (1.0_r8-fu(i))*solad(i,ib)
       sol3i(i) = fu(i)*sol2i(i) + (1.0_r8-fu(i))*solai(i,ib)
       !
    END DO
    !
    ! (h,i) calculate fluxes absorbed by lower veg, snow-free soil
    !       and snow, using results of (g) and unit-flux results
    !       of solalb(a)
    !
    DO  j=1,nsol
       !
       i = indsol(j)
       soll(i) = soll(i) + sol3d(i)*ablod(i) + sol3i(i)*abloi(i)
       !
       xd = (fl(i)*flodd(i) + 1.0_r8-fl(i)) * sol3d(i)
       !
       xi = fl(i)*(sol3d(i)*flodi(i) + sol3i(i)*floii(i))   &
            + (1.0_r8-fl(i)) * sol3i(i)
       !
       solg(i) = solg(i)             &
            + (1.0_r8-albsod(i))*xd + (1.0_r8-albsoi(i))*xi
       !
       soli(i) = soli(i)                    &
            + (1.0_r8-albsnd(i))*sol3d(i)  &
            + (1.0_r8-albsni(i))*sol3i(i)
       !
    END DO
    !
    ! estimate absorbed pars at top of canopy, toppar[u,l] and
    ! some canopy scaling parameters
    !
    ! this neglects complications due to differing values of dead vs 
    ! live elements, averaged into rhoveg, tauveg in vegdat, and 
    ! modifications of omega due to intercepted snow in twoset
    !
    ! do only for visible band (ib=1)
    !
    IF (ib.EQ.1) THEN
       !
       DO j = 1, nsol
          !
          i = indsol(j)
          !
          ! the canopy scaling algorithm assumes that the net photosynthesis
          ! is proportional to absored par (apar) during the daytime. during night,
          ! the respiration is scaled using a 10-day running-average daytime canopy
          ! scaling parameter.
          !
          ! apar(x) = A exp(-k x) + B exp(-h x) + C exp(h x)
          !
          ! some of the required terms (i.e. term[u,l] are calculated in the subroutine 'twostr'.
          ! in the equations below, 
          !
          !   A = scalcoefu(i,1) = term[u,l](i,1) * ipardir(0)
          !   B = scalcoefu(i,2) = term[u,l](i,2) * ipardir(0) + term[u,l](i,3) * ipardif(0)
          !   C = scalcoefu(i,3) = term[u,l](i,4) * ipardir(0) + term[u,l](i,5) * ipardif(0)
          !   A + B + C = scalcoefu(i,4) = also absorbed par at canopy of canopy by leaves & stems
          !
          ! upper canopy:
          !
          ! total single-sided lai+sai
          !
          xaiu = MAX (lai(i,2)+sai(i,2), epsilon)
          !
          ! some terms required for use in canopy scaling:
          !
          scalcoefu(i,1) = termu(i,1) * solad(i,ib)
          !
          scalcoefu(i,2) = termu(i,2) * solad(i,ib) +    &
               termu(i,3) * solai(i,ib)
          !
          scalcoefu(i,3) = termu(i,4) * solad(i,ib) +    &
               termu(i,5) * solai(i,ib)
          !
          scalcoefu(i,4) = scalcoefu(i,1) +   &
               scalcoefu(i,2) +   &
               scalcoefu(i,3)
          !
          ! apar of the "top" leaves of the canopy
          !
          topparu(i) = scalcoefu(i,4) * lai(i,2) / xaiu
          !
          ! lower canopy:
          !
          ! total single-sided lai+sai
          !
          xail = MAX (lai(i,1)+sai(i,1), epsilon)
          !
          ! some terms required for use in canopy scaling:
          !
          scalcoefl(i,1) = terml(i,1) * sol3d(i)
          !
          scalcoefl(i,2) = terml(i,2) * sol3d(i) +  &
               terml(i,3) * sol3i(i)
          !
          scalcoefl(i,3) = terml(i,4) * sol3d(i) +  &
               terml(i,5) * sol3i(i)
          !
          scalcoefl(i,4) = scalcoefl(i,1) +   &
               scalcoefl(i,2) +   &
               scalcoefl(i,3)
          !
          ! apar of the "top" leaves of the canopy
          !
          topparl(i) = scalcoefl(i,4) * lai(i,1) / xail
          !
       END DO
       !
    END IF
    !
    RETURN
  END SUBROUTINE solarf
  !
  !
  ! ------------------------------------------------------------------------
  SUBROUTINE twostr (       &! INTENT(IN   )
       nVegClass    , &! INTENT(IN   )
       vegtype0     , &! INTENT(IN   )   ! annual vegetation type - ibis classification
       avmuir_factor,&! INTENT(IN   )
       abvegd, &! INTENT(OUT  )
       abvegi, &! INTENT(OUT  )
       refld , &! INTENT(OUT  )
       refli , &! INTENT(OUT  )
       fbeldd, &! INTENT(OUT  )
       fbeldi, &! INTENT(OUT  )
       fbelid, &! INTENT(OUT  )
       fbelii, &! INTENT(OUT  )
       asurd , &! INTENT(IN   )
       asuri , &! INTENT(IN   )
       iv    , &! INTENT(IN   )
       coszen, &! INTENT(IN   )
       ib    , &! INTENT(IN   )
       indsol, &! INTENT(IN   )
       terml , &! INTENT(OUT  )
       termu , &! INTENT(OUT  )
       lai   , &! INTENT(IN   )
       sai   , &! INTENT(IN   )
       npoi  , &! INTENT(IN   )
       nband , &! INTENT(IN   )
       nsol  , &! INTENT(IN   )
       fwetl , &! INTENT(IN   )
       rliql , &! INTENT(IN   )
       rliqu , &! INTENT(IN   )
       rliqs , &! INTENT(IN   )
       fwetu , &! INTENT(IN   )
       fwets , &! INTENT(IN   )
       rhoveg, &! INTENT(IN   )
       tauveg, &! INTENT(IN   )
       orieh , &! INTENT(IN   )
       oriev , &! INTENT(IN   )
       tl    , &! INTENT(IN   )
       ts    , &! INTENT(IN   )
       tu    , &! INTENT(IN   )
       pi    , &! INTENT(IN   )    
       tmelt , &! INTENT(IN   )
       epsilon )! INTENT(IN   )
    ! ------------------------------------------------------------------------
    !
    ! solves canonical radiative transfer problem of two-stream veg
    ! layer + underlying surface of known albedo, for unit incoming
    ! direct or diffuse flux. returns flux absorbed within layer,
    ! reflected flux, and downward fluxes below layer. note that all
    ! direct fluxes are per unit horizontal zrea, ie, already 
    ! including a factor cos (zenith angle)
    !
    ! the solutions for the twostream approximation follow Sellers (1985),
    ! and Bonan (1996) (the latter being the LSM documentation)
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: npoi
    INTEGER, INTENT(IN   ) :: nband
    INTEGER, INTENT(IN   ) :: nsol      
    INTEGER, INTENT(IN   ) :: nVegClass   ! INTENT(IN   )
    REAL(KIND=r8), INTENT(IN   ) :: vegtype0  (npoi)   ! INTENT(IN   )   ! annual vegetation type - ibis classification
    REAL(KIND=r8), INTENT(IN   ) :: avmuir_factor(nVegClass,2)! INTENT(IN   )
    REAL(KIND=r8), INTENT(IN   ) :: fwetl (npoi)          ! fraction of lower canopy stem & leaf area 
    ! wetted by intercepted liquid and/or snow
    REAL(KIND=r8), INTENT(IN   ) :: rliql (npoi)     ! proportion of fwetl due to liquid
    REAL(KIND=r8), INTENT(IN   ) :: rliqu (npoi)     ! proportion of fwetu due to liquid 
    REAL(KIND=r8), INTENT(IN   ) :: rliqs (npoi)     ! proportion of fwets due to liquid
    REAL(KIND=r8), INTENT(IN   ) :: fwetu (npoi)     ! fraction of upper canopy leaf area wetted by
    ! intercepted liquid and/or snow
    REAL(KIND=r8), INTENT(IN   ) :: fwets (npoi)     ! fraction of upper canopy stem area wetted by
    ! intercepted liquid and/or snow

    REAL(KIND=r8), INTENT(IN   ) :: rhoveg(nband,2)  ! reflectance of an average leaf/stem
    REAL(KIND=r8), INTENT(IN   ) :: tauveg(nband,2)  ! transmittance of an average leaf/stem
    REAL(KIND=r8), INTENT(IN   ) :: orieh (2)        ! fraction of leaf/stems with horizontal orientation
    REAL(KIND=r8), INTENT(IN   ) :: oriev (2)        ! fraction of leaf/stems with vertical
    REAL(KIND=r8), INTENT(IN   ) :: tl    (npoi)     ! temperature of lower canopy leaves & stems(K)
    REAL(KIND=r8), INTENT(IN   ) :: ts    (npoi)     ! temperature of upper canopy stems (K)
    REAL(KIND=r8), INTENT(IN   ) :: tu    (npoi)     ! temperature of upper canopy leaves (K)

    REAL(KIND=r8), INTENT(IN   ) :: pi        
    REAL(KIND=r8), INTENT(IN   ) :: tmelt            ! freezing point of water (K) 
    REAL(KIND=r8), INTENT(IN   ) :: epsilon          ! small quantity to avoid zero-divides and other
    ! truncation or machine-limit troubles with small
    ! values. should be slightly greater than o(1)
    ! machine precision

    !
    ! 1: lower canopy
    ! 2: upper canopy
    !
    REAL(KIND=r8), INTENT(IN   ) :: lai   (npoi,2)   ! canopy single-sided leaf area index (area leaf/area veg)
    REAL(KIND=r8), INTENT(IN   ) :: sai   (npoi,2)   ! current single-sided stem area index

    INTEGER, INTENT(IN   ) :: indsol(npoi)     ! index of current strip for points with positive coszen
    REAL(KIND=r8), INTENT(OUT  ) :: terml (npoi,7)   ! term needed in lower canopy scaling
    REAL(KIND=r8), INTENT(OUT  ) :: termu (npoi,7)   ! term needed in upper canopy scaling
    !
    ! Arguments
    !
    INTEGER, INTENT(IN   ) :: ib               ! waveband number (1= visible, 2= near-IR)
    INTEGER, INTENT(IN   ) :: iv               ! 1 for lower, 2 for upper story params (supplied)
    !
    REAL(KIND=r8), INTENT(OUT  ) :: abvegd(npoi)     ! direct flux absorbed by two-stream layer (returned)
    REAL(KIND=r8), INTENT(OUT  ) :: abvegi(npoi)     ! diffuse flux absorbed by two-stream layer (returned)
    REAL(KIND=r8), INTENT(OUT  ) :: refld (npoi)      ! direct flux reflected above two-stream layer (returned)
    REAL(KIND=r8), INTENT(OUT  ) :: refli(npoi)      ! diffuse flux reflected above two-stream layer (returned)
    REAL(KIND=r8), INTENT(OUT  ) :: fbeldd(npoi)     ! downward direct  flux below two-stream layer(returned)
    REAL(KIND=r8), INTENT(OUT  ) :: fbeldi(npoi)     ! downward direct  flux below two-stream layer(returned)
    REAL(KIND=r8), INTENT(OUT  ) :: fbelid(npoi)     ! downward diffuse flux below two-stream layer(returned)
    REAL(KIND=r8), INTENT(OUT  ) :: fbelii(npoi)     ! downward diffuse flux below two-stream layer(returned)
    REAL(KIND=r8), INTENT(IN   ) :: asurd(npoi,nband)! direct  albedo of underlying surface (supplied)
    REAL(KIND=r8), INTENT(IN   ) :: asuri(npoi,nband)! diffuse albedo of underlying surface (supplied)
    REAL(KIND=r8), INTENT(IN   ) :: coszen(npoi)     ! cosine of direct zenith angle (supplied, must be gt 0)
    !
    ! local variables
    !
    INTEGER :: j                ! loop indice on number of points with >0 coszen
    INTEGER :: i                ! indice of point in (1, npoi) array. 
    !
    REAL(KIND=r8) b, c, c0, d, f, h, k, q, p, sigma
    !
    REAL(KIND=r8) ud1, ui1, ud2, ui2, ud3, xai, s1
    REAL(KIND=r8) s2, p1, p2, p3, p4 
    REAL(KIND=r8) rwork, dd1, di1, dd2, di2, h1, h2
    REAL(KIND=r8) h3, h4, h5, h6, h7, h8
    REAL(KIND=r8) h9, h10, absurd, absuri
    !
    ! [d,i] => per unit incoming direct, diffuse (indirect) flux
    !
    REAL(KIND=r8) omega(npoi)       !
    REAL(KIND=r8) betad(npoi)       !
    REAL(KIND=r8) betai(npoi)       !
    REAL(KIND=r8) avmu(npoi)        !
    REAL(KIND=r8) gdir(npoi)        !
    REAL(KIND=r8) tmp0(npoi)        !

    !
    ! do nothing if all points in current strip have coszen le 0
    !
    IF (nsol.EQ.0) RETURN
    !
    ! calculate two-stream parameters omega, betad, betai, avmu, gdir
    !
    CALL twoset (       &! INTENT(IN   )
         nVegClass    , &! INTENT(IN   )
         vegtype0     , &! INTENT(IN   )   ! annual vegetation type - ibis classification
         avmuir_factor, &! INTENT(IN   )
         omega  , &! INTENT(OUT   )
         betad  , &! INTENT(OUT   )
         betai  , &! INTENT(OUT   )
         avmu   , &! INTENT(OUT   )
         gdir   , &! INTENT(OUT   )
         coszen , &! INTENT(IN    )
         iv     , &! INTENT(IN    )
         ib     , &! INTENT(IN    )
         indsol , &! INTENT(IN    )
         fwetl  , &! INTENT(IN    )
         rliql  , &! INTENT(IN    )
         rliqu  , &! INTENT(IN    )
         rliqs  , &! INTENT(IN    )
         fwetu  , &! INTENT(IN    )
         fwets  , &! INTENT(IN    )
         rhoveg , &! INTENT(IN    )
         tauveg , &! INTENT(IN    )
         oriev  , &! INTENT(IN    )
         orieh  , &! INTENT(IN    )
         tl          , &! INTENT(IN    )
         lai    , &! INTENT(IN    )
         sai    , &! INTENT(IN    )
         tu          , &! INTENT(IN    )
         ts          , &! INTENT(IN    )
         nband  , &! INTENT(IN   )
         npoi   , &! INTENT(IN   )
         nsol   , &! INTENT(IN   )
         pi          , &! INTENT(IN   )
         tmelt  , &! INTENT(IN   )
         epsilon  )! INTENT(IN   ) 
    !
    DO  j=1,nsol
       !
       i = indsol(j)
       !
       ! the notations used here are taken from page 21 of Bonan's LSM documentation:
       ! Bonan, 1996: A Land Surface Model (LSM version 1.0) for ecological, hydrological,
       ! and atmospheric studies: Technical description and user's guide. NCAR Technical
       ! Note. NCAR/TN-417+STR, January 1996.
       !
       ! some temporary variables are also introduced, which are from the original
       ! lsx model.
       !
       b = 1.0_r8 - omega(i) * (1.0_r8-betai(i))
       c = omega(i) * betai(i)
       !
       tmp0(i) = b*b-c*c!pkubota deve ser alterado
       !
       q = SQRT ( MAX(0.000000000001_r8, tmp0(i)) )
       k = gdir(i) / MAX(coszen(i), 0.01_r8)
       p = avmu(i) * k
!!!WRITE(*,*)coszen(i),k,p,avmu(i)
       !
       ! next line perturbs p if p = q
       !
       IF ( ABS(p-q) .LT. .001_r8*p ) p = (1.0_r8+SIGN(.001_r8,p-q)) * p
       !
       c0 = omega(i) * p
       d = c0 * betad(i)
       f = c0 * (1.0_r8-betad(i))
       h = q / avmu(i)
       !
       sigma = p*p - tmp0(i)
       !
       ! direct & diffuse parameters are separately calculated
       !
       ud1 = b - c/asurd(i,ib)
       ui1 = b - c/asuri(i,ib)
       ud2 = b - c*asurd(i,ib)
       ui2 = b - c*asuri(i,ib)
       ud3 = f + c*asurd(i,ib)
       !
       xai = MAX (lai(i,iv) + sai(i,iv), epsilon)
       !
       s1 = EXP(-1.0_r8*h*xai)
       s2 = EXP(-1.0_r8*k*xai)
       !
       p1 = b + q
       p2 = b - q
       p3 = b + p
       p4 = b - p
       rwork = 1.0_r8/s1
       !
       ! direct & diffuse parameters are separately calculated
       !

       dd1 = p1*(ud1-q)*rwork - p2*(ud1+q)*s1
       di1 = p1*(ui1-q)*rwork - p2*(ui1+q)*s1
       dd2 = (ud2+q)*rwork - (ud2-q)*s1
       di2 = (ui2+q)*rwork - (ui2-q)*s1
       h1 = -1.0_r8*d*p4 - c*f
       rwork = s2*(d-c-h1*(ud1+p)/sigma)
       h2 = 1.0_r8/dd1*( (d-h1*p3/sigma)*(ud1-q)/s1 -   &
            p2*rwork )
       h3 = -1.0_r8/dd1*( (d-h1*p3/sigma)*(ud1+q)*s1 -   &
            p1*rwork )
       h4 = -1.0_r8*f*p3 - c*d
       rwork = s2*(ud3-h4*(ud2-p)/sigma)
       h5 = -1.0_r8/dd2*( h4*(ud2+q)/(sigma*s1) +  &
            rwork )
       h6 = 1.0_r8/dd2*( h4*s1*(ud2-q)/sigma +    &
            rwork )
       h7 = c*(ui1-q)/(di1*s1)
       h8 = -1.0_r8*c*s1*(ui1+q)/di1
       h9 = (ui2+q)/(di2*s1)
       h10= -1.0_r8*s1*(ui2-q)/di2
       !
       ! save downward direct, diffuse fluxes below two-stream layer
       !
       fbeldd(i) = s2
       fbeldi(i) = 0.0_r8
       fbelid(i) = h4/sigma*s2 + h5*s1 + h6/s1
       fbelii(i) = h9*s1 + h10/s1
       !
       ! save reflected flux, and flux absorbed by two-stream layer
       !
       refld(i) = h1/sigma + h2 + h3
       refli(i) = h7 + h8
       absurd = (1.0_r8-asurd(i,ib)) * fbeldd(i)  &
            + (1.0_r8-asuri(i,ib)) * fbelid(i)
       absuri = (1.0_r8-asuri(i,ib)) * fbelii(i)
       !
       abvegd(i) = MAX (0.0_r8, 1.0_r8 - refld(i) - absurd)
       abvegi(i) = MAX (0.0_r8, 1.0_r8 - refli(i) - absuri)
       !
       ! if no veg, make sure abveg (flux absorbed by veg) is exactly zero
       ! if this is not done, roundoff error causes small (+/-)
       ! sols, soll values in solarf and subsequent problems in turvap
       ! via stomata
       !
       IF (xai.LT.epsilon) abvegd(i) = 0.0_r8
       IF (xai.LT.epsilon) abvegi(i) = 0.0_r8
       !
       ! some terms needed in canopy scaling
       ! the canopy scaling algorithm assumes that the net photosynthesis
       ! is proportional to absored par (apar) during the daytime. during night,
       ! the respiration is scaled using a 10-day running-average daytime canopy
       ! scaling parameter.
       !
       ! apar(x) = A exp(-k x) + B exp(-h x) + C exp(h x)
       !
       ! in the equations below, 
       !
       !   k = term[u,l](i,6)
       !   h = term[u,l](i,7)
       !
       !   A = term[u,l](i,1) * ipardir(0)
       !   B = term[u,l](i,2) * ipardir(0) + term[u,l](i,3) * ipardif(0)
       !   C = term[u,l](i,4) * ipardir(0) + term[u,l](i,5) * ipardif(0)
       !
       ! calculations performed only for visible (ib=1)
       !
       IF (ib.EQ.1) THEN
          !
          IF (iv.EQ.1) THEN
             terml(i,1) = k * (1.0_r8 + (h4-h1) / sigma)
             terml(i,2) = h * (h5 - h2)
             terml(i,3) = h * (h9 - h7)
             terml(i,4) = h * (h3 - h6)
             terml(i,5) = h * (h8 - h10)
             terml(i,6) = k
             terml(i,7) = h
          ELSE
             termu(i,1) = k * (1.0_r8 + (h4-h1) / sigma)
             termu(i,2) = h * (h5 - h2)
             termu(i,3) = h * (h9 - h7)
             termu(i,4) = h * (h3 - h6)
             termu(i,5) = h * (h8 - h10)
             termu(i,6) = k
             termu(i,7) = h
          END IF
          !
       END IF
       !
    END DO
    !
    RETURN
  END SUBROUTINE twostr
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE twoset (       &! INTENT(IN   )
       nVegClass    , &! INTENT(IN   )
       vegtype0     , &! INTENT(IN   )   ! annual vegetation type - ibis classification
       avmuir_factor,&! INTENT(IN   )
       omega  , &! INTENT(OUT   )
       betad  , &! INTENT(OUT   )
       betai  , &! INTENT(OUT   )
       avmu   , &! INTENT(OUT   )
       gdir   , &! INTENT(OUT   )
       coszen , &! INTENT(IN    )
       iv     , &! INTENT(IN    )
       ib     , &! INTENT(IN    )
       indsol , &! INTENT(IN    )
       fwetl  , &! INTENT(IN    )
       rliql  , &! INTENT(IN    )
       rliqu  , &! INTENT(IN    )
       rliqs  , &! INTENT(IN    )
       fwetu  , &! INTENT(IN    )
       fwets  , &! INTENT(IN    )
       rhoveg , &! INTENT(IN    )
       tauveg , &! INTENT(IN    )
       oriev  , &! INTENT(IN    )
       orieh  , &! INTENT(IN    )
       tl     , &! INTENT(IN    )
       lai    , &! INTENT(IN    )
       sai    , &! INTENT(IN    )
       tu     , &! INTENT(IN    )
       ts     , &! INTENT(IN    )
       nband  , &! INTENT(IN   )
       npoi   , &! INTENT(IN   )
       nsol   , &! INTENT(IN   )
       pi     , &! INTENT(IN   )
       tmelt  , &! INTENT(IN   )
       epsilon  )! INTENT(IN   ) 
    ! ---------------------------------------------------------------------
    !
    ! sets two-stream parameters, given single-element transmittance
    ! and reflectance, leaf orientation weights, and cosine of the
    ! zenith angle, then adjusts for amounts of intercepted snow
    !
    ! the two-stream parameters omega,betad,betai are weighted 
    ! combinations of the "exact" values for the 3 orientations:
    ! all vertical, all horizontal, or all random (ie, spherical)
    !
    ! the vertical, horizontal weights are in oriev,orieh (comveg)
    !
    ! the "exact" expressions are as derived in my notes(8/6/91,p.6).
    ! note that values for omega*betad and omega*betai are calculated
    ! and then divided by the new omega, since those products are 
    ! actually used in twostr. also those depend *linearly* on the
    ! single-element transmittances and reflectances tauveg, rhoveg,
    ! which are themselves linear weights of leaf and stem values 
    !
    ! for random orientation, omega*betad depends on coszen according
    ! to the function in array tablemu
    !
    ! the procedure is approximate since omega*beta[d,i] and gdir
    ! should depend non-linearly on the complete leaf-angle
    ! distribution. then we should also treat leaf and stem angle
    ! distributions separately, and allow for the cylindrical
    ! shape of stems (norman and jarvis, app.b; the expressions 
    ! below are appropriate for flat leaves)
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: nband
    INTEGER, INTENT(IN   ) :: npoi
    INTEGER, INTENT(IN   ) :: nsol
    REAL(KIND=r8), INTENT(IN   ) :: pi
    REAL(KIND=r8), INTENT(IN   ) :: tmelt
    REAL(KIND=r8), INTENT(IN   ) :: epsilon        ! small quantity to avoid zero-divides and other
    ! truncation or machine-limit troubles with small
    ! values. should be slightly greater than o(1)
    ! machine precision
    INTEGER      , INTENT(IN   ) :: nVegClass               ! INTENT(IN   )
    REAL(KIND=r8), INTENT(IN   ) :: vegtype0 (npoi)         ! INTENT(IN   )   ! annual vegetation type - ibis classification
    REAL(KIND=r8), INTENT(IN   ) :: avmuir_factor(nVegClass,2)! INTENT(IN   )
    REAL(KIND=r8), INTENT(IN   ) :: rhoveg(nband,2)! reflectance of an average leaf/stem
    REAL(KIND=r8), INTENT(IN   ) :: tauveg(nband,2)! transmittance of an average leaf/stem
    REAL(KIND=r8), INTENT(IN   ) :: oriev (2)             ! fraction of leaf/stems with vertical
    REAL(KIND=r8), INTENT(IN   ) :: orieh (2)      ! fraction of leaf/stems with horizontal orientation
    REAL(KIND=r8), INTENT(IN   ) :: tl    (npoi)   ! temperature of lower canopy leaves & stems(K)
    REAL(KIND=r8), INTENT(IN   ) :: lai   (npoi,2) ! canopy single-sided leaf area index (area leaf/area veg)
    REAL(KIND=r8), INTENT(IN   ) :: sai   (npoi,2) ! current single-sided stem area index
    REAL(KIND=r8), INTENT(IN   ) :: tu    (npoi)   ! temperature of upper canopy leaves (K)
    REAL(KIND=r8), INTENT(IN   ) :: ts    (npoi)   ! temperature of upper canopy stems (K)

    INTEGER, INTENT(IN   ) :: indsol(npoi) ! index of current strip for points with positive coszen
    REAL(KIND=r8), INTENT(IN   ) :: fwetl (npoi) ! fraction of lower canopy stem & leaf area wetted by intercepted liquid and/or snow
    REAL(KIND=r8), INTENT(IN   ) :: rliql (npoi) ! proportion of fwetl due to liquid
    REAL(KIND=r8), INTENT(IN   ) :: rliqu (npoi) ! proportion of fwetu due to liquid
    REAL(KIND=r8), INTENT(IN   ) :: rliqs (npoi) ! proportion of fwets due to liquid
    REAL(KIND=r8), INTENT(IN   ) :: fwetu (npoi) ! fraction of upper canopy leaf area wetted by intercepted liquid and/or snow
    REAL(KIND=r8), INTENT(IN   ) :: fwets (npoi) ! fraction of upper canopy stem area wetted by intercepted liquid and/or snow
    !
    ! Arguments (all quantities are returned unless otherwise note)
    !
    INTEGER, INTENT(IN    ) :: ib              ! waveband number (1= visible, 2= near-IR)
    INTEGER, INTENT(IN    ) :: iv              ! 1 for lower, 2 for upper story params (supplied)
    !
    REAL(KIND=r8), INTENT(OUT   ) :: omega(npoi)     ! fraction of intercepted radiation that is scattered
    REAL(KIND=r8), INTENT(OUT   ) :: betad(npoi)     ! fraction of scattered *direct* radiation that is
    !  scattered into upwards hemisphere
    REAL(KIND=r8), INTENT(OUT   ) :: betai(npoi)     ! fraction of scattered downward *diffuse* radiation
    ! that is scattered into upwards hemisphere (or fraction
    ! of scattered upward diffuse rad. into downwards hemis)
    REAL(KIND=r8), INTENT(OUT   ) :: avmu(npoi)      ! average diffuse optical depth
    REAL(KIND=r8), INTENT(OUT   ) :: gdir(npoi)      ! average projected leaf area into solar direction
    REAL(KIND=r8), INTENT(IN    ) :: coszen(npoi)    ! cosine of solar zenith angle (supplied)
    !
    ! local variables
    !
    INTEGER j     ! loop indice on number of points with >0 coszen
    INTEGER i     ! indice of point in (1, npoi) array. 
    INTEGER, PARAMETER :: ntmu=100  !
    INTEGER, PARAMETER :: nbandas=2
    INTEGER itab,iveg
    !
    REAL(KIND=r8) zrho, ztau, orand, ztab, rwork, y
    REAL(KIND=r8) o, x
    !
    REAL(KIND=r8) otmp(npoi)
    !
    !REAL(KIND=r8) tablemu(ntmu+1)
    !REAL(KIND=r8) omegasno(nbandas)

    !SAVE    tablemu, betadsno, betaisno,omegasno
    !
    REAL(KIND=r8), PARAMETER :: tablemu(1:ntmu+1) = (/&
         0.5000_r8, 0.4967_r8, 0.4933_r8, 0.4900_r8, 0.4867_r8, 0.4833_r8, 0.4800_r8, 0.4767_r8,&
         0.4733_r8, 0.4700_r8, 0.4667_r8, 0.4633_r8, 0.4600_r8, 0.4567_r8, 0.4533_r8, 0.4500_r8,&
         0.4467_r8, 0.4433_r8, 0.4400_r8, 0.4367_r8, 0.4333_r8, 0.4300_r8, 0.4267_r8, 0.4233_r8,&
         0.4200_r8, 0.4167_r8, 0.4133_r8, 0.4100_r8, 0.4067_r8, 0.4033_r8, 0.4000_r8, 0.3967_r8,&
         0.3933_r8, 0.3900_r8, 0.3867_r8, 0.3833_r8, 0.3800_r8, 0.3767_r8, 0.3733_r8, 0.3700_r8,&
         0.3667_r8, 0.3633_r8, 0.3600_r8, 0.3567_r8, 0.3533_r8, 0.3500_r8, 0.3467_r8, 0.3433_r8,&
         0.3400_r8, 0.3367_r8, 0.3333_r8, 0.3300_r8, 0.3267_r8, 0.3233_r8, 0.3200_r8, 0.3167_r8,&
         0.3133_r8, 0.3100_r8, 0.3067_r8, 0.3033_r8, 0.3000_r8, 0.2967_r8, 0.2933_r8, 0.2900_r8,&
         0.2867_r8, 0.2833_r8, 0.2800_r8, 0.2767_r8, 0.2733_r8, 0.2700_r8, 0.2667_r8, 0.2633_r8,&
         0.2600_r8, 0.2567_r8, 0.2533_r8, 0.2500_r8, 0.2467_r8, 0.2433_r8, 0.2400_r8, 0.2367_r8,&
         0.2333_r8, 0.2300_r8, 0.2267_r8, 0.2233_r8, 0.2200_r8, 0.2167_r8, 0.2133_r8, 0.2100_r8,&
         0.2067_r8, 0.2033_r8, 0.2000_r8, 0.1967_r8, 0.1933_r8, 0.1900_r8, 0.1867_r8, 0.1833_r8,&
         0.1800_r8, 0.1767_r8, 0.1733_r8, 0.1700_r8, 0.1667_r8 /)
    !
    REAL(KIND=r8), PARAMETER :: omegasno(1:nbandas)=(/0.9_r8, 0.7_r8/)
    REAL(KIND=r8), PARAMETER :: betadsno=0.5_r8
    REAL(KIND=r8), PARAMETER :: betaisno=0.5_r8
    !DATA betadsno, betaisno /0.5, 0.5/
    !
    ! set two-stream parameters omega, betad, betai, gdir and avmu
    ! as weights of those for 100% vert,horiz,random orientations
    !
    DO  j=1,nsol
       i = indsol(j)
       !
       zrho = rhoveg(ib,iv)
       ztau = tauveg(ib,iv)
       !
       ! weight for random orientation is 1 - those for vert and horiz
       !
       orand = 1.0_r8 - oriev(iv) - orieh(iv)
       !
       omega(i) = zrho + ztau
       !
       ! ztab is transmittance coeff - for random-orientation omega*betad,
       ! given by tablemu as a function of coszen
       !
       itab = NINT (coszen(i)*ntmu + 1)
       ztab = tablemu(itab)
       rwork = 1.0_r8/omega(i)
       !
       betad(i) = (  oriev(iv) * 0.5_r8*(zrho + ztau)   &
            + orieh(iv) * zrho                &
            + orand       * ((1.0_r8-ztab)*zrho + ztab*ztau) )   &
            * rwork
       !
       betai(i) = (  oriev(iv) * 0.5_r8*(zrho + ztau)   &
            + orieh(iv) * zrho                &
            + orand       * ((2.0_r8/3.0_r8)*zrho + (1.0_r8/3.0_r8)*ztau) )  &
            * rwork
       !
       gdir(i)  = oriev(iv) * (2.0_r8/pi) *              &
            SQRT ( MAX (0.0_r8, 1.0_r8-coszen(i)*coszen(i)) )   &
            + orieh(iv) * coszen(i)   &
            + orand       * 0.5_r8
    !
    !
    !
    !      emu(i) - 1.0_r8 = - exp ( -lai(i,2) / avmuir_local )
    !
    !      exp ( -lai(i,2) / avmuir_local ) = 1.0_r8 - emu(i)
    !
    !      -lai(i,2) / avmuir_local   =log (max(1.0_r8 - emu(i),0.0000001_r8))

    !    avmu(i)   = - max(lai(i,2)+lai(i,1),0.01_r8)/log (max(1.0_r8 - 0.987_r8,0.0000001_r8))

    !PK        avmu(i) = 1.0_r8
            iveg=vegtype0 (i) 
            avmu(i) = avmuir_factor(iveg,1)
    END DO
    !
    ! adjust omega, betad and betai for amounts of intercepted snow
    ! (omegasno decreases to .6 of cold values within 1 deg of tmelt)
    !
    IF (iv.EQ.1) THEN
       !
       ! lower story 
       !
       DO  j=1,nsol
          i = indsol(j)
          y = fwetl(i)*(1.0_r8-rliql(i))
          o = omegasno(ib)*(0.6_r8 + 0.4_r8*MAX(0.0_r8,MIN(1.0_r8,(tmelt-tl(i))/1.0_r8)))
          otmp(i)  = omega(i)
          rwork = y * o
          omega(i) =  (1-y)*otmp(i)          + rwork
          betad(i) = ((1-y)*otmp(i)*betad(i) + rwork*betadsno) /  &
               omega(i)  
          betai(i) = ((1-y)*otmp(i)*betai(i) + rwork*betaisno) /  &
               omega(i)  
       END DO
       !
    ELSE
       !
       ! upper story
       !
       DO  j=1,nsol
          i = indsol(j)
          x = lai(i,iv) / MAX (lai(i,iv)+sai(i,iv), epsilon)
          y = x * fwetu(i)*(1.0_r8-rliqu(i)) + (1-x) *fwets(i)*(1.0_r8-rliqs(i))
          o = (     x  * MIN (1.0_r8, MAX (0.6_r8, (tmelt-tu(i))/0.1_r8))   &
               + (1-x) * MIN (1.0_r8, MAX (0.6_r8, (tmelt-ts(i))/0.1_r8)) )  &
               *  omegasno(ib) 
          !
          otmp(i)  = omega(i)
          rwork = y * o
          omega(i) =  (1-y)*otmp(i)          + rwork
          betad(i) = ((1-y)*otmp(i)*betad(i) + rwork*betadsno) / &
               omega(i)
          betai(i) = ((1-y)*otmp(i)*betai(i) + rwork*betaisno) /  &
               omega(i)
          !
       END DO
       !
    END IF
    !
    RETURN
  END SUBROUTINE twoset
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE irrad(npoi  , &! INTENT(IN   ) ::
       nsoilay, &! INTENT(IN   ) ::
       stef  , &! INTENT(IN   ) ::
       nVegClass , &! INTENT(IN   )
       vegtype0  , &! INTENT(IN   )! annual vegetation type - ibis classification
       avmuir_factor, &! INTENT(IN   )
       firb  , &! INTENT(OUT  ) ::
       firs  , &! INTENT(OUT  ) ::
       firu  , &! INTENT(OUT  ) ::
       firl  , &! INTENT(OUT  ) ::
       firg  , &! INTENT(OUT  ) ::
       firi  , &! INTENT(OUT  ) ::
       lai   , &! INTENT(IN   ) ::
       sai   , &! INTENT(IN   ) ::
       fu    , &! INTENT(IN   ) ::
       tu    , &! INTENT(IN   ) ::
       ts    , &! INTENT(IN   ) ::
       tl    , &! INTENT(IN   ) ::
       fl    , &! INTENT(IN   ) ::
       tg    , &! INTENT(IN   ) ::
       ti    , &! INTENT(IN   ) ::
       fi    , &! INTENT(IN   ) ::
       fira  , &! INTENT(IN   ) ::
       poros , &! INTENT(IN   ) ::
       wsoi    )! INTENT(IN   ) ::
    ! ---------------------------------------------------------------------
    !
    ! calculates overall emitted ir flux, and net absorbed minus
    ! emitted ir fluxes for upper leaves, upper stems, lower story,
    ! soil and snow. assumes upper leaves, upper stems and lower
    ! story each form a semi-transparent plane, with the upper-leaf
    ! plane just above the upper-stem plane. the soil and snow 
    ! surfaces have emissivities of 0.95.
    !
    ! the incoming flux is supplied in comatm array fira
    !
    ! the emitted ir flux by overall surface system is returned in
    ! com1d array firb - the ir fluxes absorbed by upper leaves,
    ! upper stems, lower veg, soil and snow are returned in com1d 
    ! arrays firu, firs, firl, firg and firi
    ! 
    ! other com1d arrays used are:
    !
    ! emu, ems, eml  = emissivities of the vegetation planes
    ! fup, fdown     = upward and downward fluxes below tree level
    !
    IMPLICIT NONE
    !
    !include 'compar.h'
    !include 'comatm.h'
    !include 'comsno.h'
    !include 'comsoi.h'
    !include 'comveg.h'
    INTEGER, INTENT(IN   ) :: npoi         ! total number of land points      
    INTEGER, INTENT(IN   ) :: nsoilay   
    REAL(KIND=r8), INTENT(IN   ) :: stef         ! stefan-boltzmann constant (W m-2 K-4)
    INTEGER, INTENT(IN   ) :: nVegClass ! INTENT(IN   )
    REAL(KIND=r8), INTENT(IN   ) :: vegtype0(npoi)  ! INTENT(IN   )! annual vegetation type - ibis classification
    REAL(KIND=r8), INTENT(IN   ) :: avmuir_factor(nVegClass,2)! INTENT(IN   )

    REAL(KIND=r8), INTENT(OUT  ) :: firb (npoi)  ! net upward ir radiation at reference atmospheric level za (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: firs (npoi)  ! ir radiation absorbed by upper canopy stems (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: firu (npoi)  ! ir raditaion absorbed by upper canopy leaves (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: firl (npoi)  ! ir radiation absorbed by lower canopy leaves and stems (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: firg (npoi)  ! ir radiation absorbed by soil/ice (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: firi (npoi)  ! ir radiation absorbed by snow (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: lai  (npoi,2)! canopy single-sided leaf area index (area leaf/area veg)
    REAL(KIND=r8), INTENT(IN   ) :: sai  (npoi,2)! current single-sided stem area index
    REAL(KIND=r8), INTENT(IN   ) :: fu   (npoi)  ! fraction of overall area covered by upper canopy
    REAL(KIND=r8), INTENT(IN   ) :: tu   (npoi)  ! temperature of upper canopy leaves (K)
    REAL(KIND=r8), INTENT(IN   ) :: ts   (npoi)  ! temperature of upper canopy stems (K)
    REAL(KIND=r8), INTENT(IN   ) :: tl   (npoi)  ! temperature of lower canopy leaves & stems(K)
    REAL(KIND=r8), INTENT(IN   ) :: fl   (npoi)  ! fraction of snow-free area covered by lower  canopy
    REAL(KIND=r8), INTENT(IN   ) :: tg   (npoi)  ! soil skin temperature (K)
    REAL(KIND=r8), INTENT(IN   ) :: ti   (npoi)  ! snow skin temperature (K)
    REAL(KIND=r8), INTENT(IN   ) :: fi   (npoi)  ! fractional snow cover
    REAL(KIND=r8), INTENT(IN   ) :: fira (npoi)  ! incoming ir flux (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: poros (npoi,nsoilay) ! INTENT(IN   ) ::
    REAL(KIND=r8), INTENT(IN   ) :: wsoi (npoi,nsoilay) ! INTENT(IN   ) :: !    ! fraction of soil pore space containing liquid water

    !      include 'com1d.h'
    !
    ! Local arrays:
    !
    INTEGER i,iveg            ! loop indice
    !
    !
    REAL(KIND=r8) emu   (npoi) ! ir emissivity of upper-leaves veg plane
    REAL(KIND=r8) ems   (npoi) ! ir emissivity of upper-stems veg plane
    REAL(KIND=r8) eml   (npoi) ! ir emissivity of lower-story veg plane
    REAL(KIND=r8) emg   (npoi) ! ir emissivity (gray) of soil surface
    REAL(KIND=r8) emi   (npoi) ! ir emissivity (gray) of snow surface
    REAL(KIND=r8) fdown (npoi) ! downward ir flux below tree level per overall area
    REAL(KIND=r8) fdowng(npoi) ! upward   ir flux below tree level per overall area
    REAL(KIND=r8) fup   (npoi) ! downward ir flux below lower-story veg
    REAL(KIND=r8) fupg  (npoi) ! upward   ir flux below lower-story veg
    REAL(KIND=r8) fupgb (npoi) ! upward   ir flux above bare soil surface
    REAL(KIND=r8) fupi  (npoi) ! upward   ir flux above snow surface
    REAL(KIND=r8) avmu (npoi)
    !
    ! set emissivities of soil and snow
    !
    REAL(KIND=r8), PARAMETER :: emisoil  =  0.95_r8  ! soil emissivity
    REAL(KIND=r8), PARAMETER :: emisnow  =  0.95_r8  ! snow emissivity

    !      DATA emisoil, emisnow
    !     &       /0.95, 0.95/
    !
    ! use uniform value 1.0 for average diffuse optical depth
    ! (although an array for solar, all values are set to 1 in twoset).
    ! The typical values of emissivity are 0.80-0.95 for bare soil,
    ! 0.95-0.97 for vegetated areas and 0.99 for snow 
    ! (Wilber et al. 1999; Jin 2004; Jin and Liang 2004, manuscript submitted to J. Climate)
    !
    !
    REAL(KIND=r8), PARAMETER :: avmuir_local = 1.0_r8      ! average diffuse optical depth
      
!M. Mira1, E. Valor1, R. Boluda2, V. Caselles1 and C. Coll1 Influence of the soil moisture effect on the thermal infraredemissivity
!Tethys,4,3-9,2007
      REAL(KIND=r8), PARAMETER :: par_A = -0.010_r8     ! coeff A*10-3
      REAL(KIND=r8), PARAMETER :: Par_B =  0.08_r8      ! coeff B*10-2
      REAL(KIND=r8), PARAMETER :: Par_C =  emisoil      ! coeff C

    !SAVE avmuir
    !DATA avmuir /1./
    !
    DO i=1,npoi
       iveg=vegtype0(i)  ! INTENT(IN   )! annual vegetation type - ibis classification
       !
       !
       !      emu(i) - 1.0_r8 = - exp ( -lai(i,2) / avmuir_local )
       !
       !      exp ( -lai(i,2) / avmuir_local ) = 1.0_r8 - emu(i)
       !
       !      -lai(i,2) / avmuir_local   =log (max(1.0_r8 - emu(i),0.0000001_r8))

       !       avmuir_local   = - lai(i,2)/log (max(1.0_r8 - emu(i),0.0000001_r8))
       ! avmu(i)   = - max(lai(i,2)+lai(i,1),0.01_r8)/log (max(1.0_r8 - 0.987_r8,0.0000001_r8))
       !
       emu(i) = 1.0_r8 - EXP ( -lai(i,2) / avmuir_factor(iveg,2) )
       ems(i) = 1.0_r8 - EXP ( -sai(i,2) / avmuir_factor(iveg,2) )
       eml(i) = 1.0_r8 - EXP ( -(lai(i,1)+sai(i,1)) / avmuir_factor(iveg,2))
       !
       emg(i) = emisoil
       !emg(i) = par_A*(((wsoi (i,1)*poros (i,1) + wsoi (i,2)*poros (i,2) + wsoi (i,3)*poros (i,3))/3.0_r8)**2) + &
       !         Par_B*(((wsoi (i,1)*poros (i,1) + wsoi (i,2)*poros (i,2) + wsoi (i,3)*poros (i,3))/3.0_r8)   ) + Par_C
       !emg(i) = MIN(MAX(emg(i),emisoil), 1.0_r8)
       emi(i) = emisnow
       !
       !        fu   ! fraction of overall area covered by upper canopy
       fdown(i) =  (1.0_r8-fu(i)) * fira(i)&
            + fu(i) * ( (1.0_r8-emu(i))*(1.0_r8-ems(i))*fira(i)&
            +         emu(i)* (1.0_r8-ems(i))*stef*(tu(i)**4)&
            +         ems(i)*stef*(ts(i)**4) )
       !
       fdowng(i) = (1.0_r8-eml(i))*fdown(i)  + eml(i)*stef*(tl(i)**4)
       !
       fupg(i)   = (1.0_r8-emg(i))*fdowng(i) + emg(i)*stef*(tg(i)**4)
       !
       fupgb(i)  = (1.0_r8-emg(i))*fdown(i)  + emg(i)*stef*(tg(i)**4)
       !
       fupi(i)   = (1.0_r8-emi(i))*fdown(i)  + emi(i)*stef*(ti(i)**4)
       !
       fup(i) = (1.0_r8-fi(i))*(      fl(i)*(       eml(i) *stef*(tl(i)**4) &
            + (1.0_r8-eml(i))*fupg(i) )         &
            +(1.0_r8-fl(i))*fupgb(i)         &
            )         &
            +     fi(i) * fupi(i)
       !
       firb(i) =   (1.0_r8-fu(i)) * fup(i) &
            + fu(i)  * ( (1.0_r8-emu(i))*(1.0_r8-ems(i))*fup(i)&
            +    emu(i)*stef*(tu(i)**4)&
            +    ems(i)*(1.0_r8-emu(i))*stef*(ts(i)**4) )
       !
       firu(i) =   emu(i)*ems(i)*stef*(ts(i)**4)  &
            + emu(i)*(1.0_r8-ems(i))*fup(i)           &
            + emu(i)*fira(i)           &
            - 2*emu(i)*stef*(tu(i)**4)
       !
       firs(i) =   ems(i)*emu(i)*stef*(tu(i)**4)   &
            + ems(i)*fup(i)                      &
            + ems(i)*(1.0_r8-emu(i))*fira(i)          &
            - 2*ems(i)*stef*(ts(i)**4)
       !
       firl(i) =   eml(i)*fdown(i)        &
            + eml(i)*fupg(i)         &
            - 2*eml(i)*stef*(tl(i)**4)
       !
       firg(i) =       fl(i)  * (fdowng(i) - fupg(i))   &
            + (1.0_r8-fl(i)) * (fdown(i)  - fupgb(i))
       !
       firi(i) =   fdown(i) - fupi(i)
       !
    END DO
    !
    RETURN
  END  SUBROUTINE irrad

  !
  ! ######  #    #  #####  #####     ##         #####      #           ##          #####     #          ####   #    #
  ! #          ##   #  #    # #    #   #  #   #    #     #          #  #      #            #         #    #  ##   #
  ! #####   # #  #  #    # #    #  #    #  #    #     #         #    #     #            #         #    #  # #  #
  ! #          #  # #  #    # #####   ######  #    #     #         ######     #            #         #    #  #  # #
  ! #          #   ##  #    # #   #   #    #  #    #     #         #    #     #            #         #    #  #   ##
  ! ######  #    #  #####  #    #  #    #  #####      #         #    #     #            #          ####   #    #
  !
  ! ---------------------------------------------------------------------


  ! ---------------------------------------------------------------------
  SUBROUTINE const (arr   , & ! INTENT(OUT  )
       nar   , & ! INTENT(IN)
       value   ) ! INTENT(IN)
    ! ---------------------------------------------------------------------
    !
    ! sets all elements of REAL(KIND=r8) vector arr to value
    !
    IMPLICIT NONE
    !
    ! Arguments
    !
    INTEGER, INTENT(IN   ) :: nar
    REAL(KIND=r8)   , INTENT(IN   ) :: value
    REAL(KIND=r8)   , INTENT(OUT  ) :: arr(nar)
    !
    ! Local variables
    !
    INTEGER :: j
    !
    DO j = 1, nar
       arr(j) = value
    END DO
    !
    RETURN
  END SUBROUTINE const

  ! ---------------------------------------------------------------------
  SUBROUTINE scopy (nt       , & ! INTENT(IN   )
       arr      , & ! INTENT(IN   )
       brr        ) ! INTENT(OUT  )
    ! ---------------------------------------------------------------------
    !
    ! copies array arr to brr,for 1st nt words of arr
    !
    IMPLICIT NONE
    !
    ! Arguments
    !
    INTEGER, INTENT(IN   ) ::  nt     
    REAL(KIND=r8)   , INTENT(IN   ) ::  arr(nt)     ! input
    REAL(KIND=r8)   , INTENT(OUT  ) ::  brr(nt)     ! output
    !
    ! Local variables
    !
    INTEGER  ia
    !
    DO  ia = 1, nt
       brr(ia) = arr(ia)
    END DO
    !
    RETURN
  END SUBROUTINE scopy
END MODULE Sfc_Ibis_LsxMain
