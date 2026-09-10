!+---+-----------------------------------------------------------------+
!+---+-----------------------------------------------------------------+
!..This set of routines facilitates computing radar reflectivity.
!.. This module is more library code whereas the individual microphysics
!.. schemes contains specific details needed for the final computation,
!.. so refer to location within each schemes calling the routine named
!.. rayleigh_soak_wetgraupel.
!.. The bulk of this code originated from Ulrich Blahak (Germany) and
!.. was adapted to WRF by G. Thompson.  This version of code is only
!.. intended for use when Rayleigh scattering principles dominate and
!.. is not intended for wavelengths in which Mie scattering is a
!.. significant portion.  Therefore, it is well-suited to use with
!.. 5 or 10 cm wavelength like USA NEXRAD radars.
!.. This code makes some rather simple assumptions about water
!.. coating on outside of frozen species (snow/graupel).  Fraction of
!.. meltwater is simply the ratio of mixing ratio below melting level
!.. divided by mixing ratio at level just above highest T>0C.  Also,
!.. immediately 90% of the melted water exists on the ice's surface
!.. and 10% is embedded within ice.  No water is "shed" at all in these
!.. assumptions. The code is quite slow because it does the reflectivity
!.. calculations based on 50 individual size bins of the distributions.
!+---+-----------------------------------------------------------------+

MODULE module_mp_radar
  USE Parallelism,Only: MsgOne

  IMPLICIT NONE
SAVE

  !      USE module_wrf_error

  PUBLIC :: rayleigh_soak_wetgraupel
  PUBLIC :: radar_init
  PRIVATE :: m_complex_water_ray
  PRIVATE :: m_complex_ice_maetzler
  PRIVATE :: m_complex_maxwellgarnett
  PRIVATE :: get_m_mix_nested
  PRIVATE :: get_m_mix
  PRIVATE :: WGAMMA
  PRIVATE :: GAMMLN

 ! INTEGER      , PARAMETER :: r8  = SELECTED_REAL_KIND(P=13,R=300)
  INTEGER,PRIVATE      , PARAMETER :: r8  = SELECTED_REAL_KIND(15)

  INTEGER, PARAMETER, PUBLIC:: nrbins = 50
  REAL(KIND=r8), DIMENSION(nrbins+1), PUBLIC:: xxDx
  REAL(KIND=r8), DIMENSION(nrbins), PUBLIC:: xxDs,xdts,xxDg,xdtg
  REAL(KIND=r8), PARAMETER, PUBLIC:: lamda_radar = 0.10_r8           ! in meters
  REAL(KIND=r8), PUBLIC:: K_w, PI5, lamda4
  COMPLEX*16, PUBLIC:: m_w_0, m_i_0
  REAL(KIND=r8), DIMENSION(nrbins+1), PUBLIC:: simpson
  REAL(KIND=r8), DIMENSION(3), PARAMETER, PUBLIC:: basis =       &
       (/1.0_r8/3.0_r8, 4.0_r8/3.0_r8, 1.0_r8/3.0_r8/)
  REAL(KIND=r8), DIMENSION(4), PUBLIC:: xcre, xcse, xcge, xcrg, xcsg, xcgg
  REAL(KIND=r8), PUBLIC:: xam_r, xbm_r, xmu_r, xobmr
  REAL(KIND=r8), PUBLIC:: xam_s, xbm_s, xmu_s, xoams, xobms, xocms
  REAL(KIND=r8), PUBLIC:: xam_g, xbm_g, xmu_g, xoamg, xobmg, xocmg
  REAL(KIND=r8), PUBLIC:: xorg2, xosg2, xogg2

  INTEGER, PARAMETER, PUBLIC:: slen = 20
  CHARACTER(len=slen), PUBLIC::                                     &
       mixingrulestring_s, matrixstring_s, inclusionstring_s,    &
       hoststring_s, hostmatrixstring_s, hostinclusionstring_s,  &
       mixingrulestring_g, matrixstring_g, inclusionstring_g,    &
       hoststring_g, hostmatrixstring_g, hostinclusionstring_g

  !..Single melting snow/graupel particle 90% meltwater on external sfc
  REAL(KIND=r8), PARAMETER:: melt_outside_s = 0.90_r8
  REAL(KIND=r8), PARAMETER:: melt_outside_g = 0.90_r8

  CHARACTER*256:: radar_debug

CONTAINS

  !+---+-----------------------------------------------------------------+
  !+---+-----------------------------------------------------------------+
  !+---+-----------------------------------------------------------------+

  SUBROUTINE radar_init()

    IMPLICIT NONE
    INTEGER:: n
    PI5 = 3.14159_r8*3.14159_r8*3.14159_r8*3.14159_r8*3.14159_r8
    lamda4 = lamda_radar*lamda_radar*lamda_radar*lamda_radar
    m_w_0 = m_complex_water_ray (lamda_radar, 0.00_r8)
    m_i_0 = m_complex_ice_maetzler (lamda_radar, 0.00_r8)
    K_w = (ABS( (m_w_0*m_w_0 - 1.0_r8) /(m_w_0*m_w_0 + 2.0_r8) ))**2

    DO n = 1, nrbins+1
       simpson(n) = 0.00_r8
    ENDDO
    DO n = 1, nrbins-1, 2
       simpson(n) = simpson(n) + basis(1)
       simpson(n+1) = simpson(n+1) + basis(2)
       simpson(n+2) = simpson(n+2) + basis(3)
    ENDDO

    DO n = 1, slen
       mixingrulestring_s(n:n) = CHAR(0)
       matrixstring_s(n:n) = CHAR(0)
       inclusionstring_s(n:n) = CHAR(0)
       hoststring_s(n:n) = CHAR(0)
       hostmatrixstring_s(n:n) = CHAR(0)
       hostinclusionstring_s(n:n) = CHAR(0)
       mixingrulestring_g(n:n) = CHAR(0)
       matrixstring_g(n:n) = CHAR(0)
       inclusionstring_g(n:n) = CHAR(0)
       hoststring_g(n:n) = CHAR(0)
       hostmatrixstring_g(n:n) = CHAR(0)
       hostinclusionstring_g(n:n) = CHAR(0)
    ENDDO

    mixingrulestring_s = 'maxwellgarnett'
    hoststring_s = 'air'
    matrixstring_s = 'water'
    inclusionstring_s = 'spheroidal'
    hostmatrixstring_s = 'icewater'
    hostinclusionstring_s = 'spheroidal'

    mixingrulestring_g = 'maxwellgarnett'
    hoststring_g = 'air'
    matrixstring_g = 'water'
    inclusionstring_g = 'spheroidal'
    hostmatrixstring_g = 'icewater'
    hostinclusionstring_g = 'spheroidal'

    !..Create bins of snow (from 100 microns up to 2 cm).
    xxDx(1) = 100.e-6_r8
    xxDx(nrbins+1) = 0.020_r8
    DO n = 2, nrbins
       xxDx(n) = DEXP(DFLOAT(n-1)/DFLOAT(nrbins) &
            *DLOG(xxDx(nrbins+1)/xxDx(1)) +DLOG(xxDx(1)))
    ENDDO
    DO n = 1, nrbins
       xxDs(n) = DSQRT(xxDx(n)*xxDx(n+1))
       xdts(n) = xxDx(n+1) - xxDx(n)
    ENDDO

    !..Create bins of graupel (from 100 microns up to 5 cm).
    xxDx(1) = 100.e-6
    xxDx(nrbins+1) = 0.050_r8
    DO n = 2, nrbins
       xxDx(n) = DEXP(DFLOAT(n-1)/DFLOAT(nrbins) &
            *DLOG(xxDx(nrbins+1)/xxDx(1)) +DLOG(xxDx(1)))
    ENDDO
    DO n = 1, nrbins
       xxDg(n) = DSQRT(xxDx(n)*xxDx(n+1))
       xdtg(n) = xxDx(n+1) - xxDx(n)
    ENDDO


    !..The calling program must set the m(D) relations and gamma shape
    !.. parameter mu for rain, snow, and graupel.  Easily add other types
    !.. based on the template here.  For majority of schemes with simpler
    !.. exponential number distribution, mu=0.

    xcre(1) = 1.0_r8 + xbm_r
    xcre(2) = 1.0_r8 + xmu_r
    xcre(3) = 1.0_r8 + xbm_r + xmu_r
    xcre(4) = 1.0_r8 + 2.0_r8*xbm_r + xmu_r
    DO n = 1, 4
       xcrg(n) = WGAMMA(xcre(n))
    ENDDO
    xorg2 = 1.0_r8/xcrg(2)

    xcse(1) = 1.0_r8 + xbm_s
    xcse(2) = 1.0_r8 + xmu_s
    xcse(3) = 1.0_r8 + xbm_s + xmu_s
    xcse(4) = 1.0_r8 + 2.0_r8*xbm_s + xmu_s
    DO n = 1, 4
       xcsg(n) = WGAMMA(xcse(n))
    ENDDO
    xosg2 = 1.0_r8/xcsg(2)

    xcge(1) = 1.0_r8 + xbm_g
    xcge(2) = 1.0_r8 + xmu_g
    xcge(3) = 1.0_r8 + xbm_g + xmu_g
    xcge(4) = 1.0_r8 + 2.0_r8*xbm_g + xmu_g
    DO n = 1, 4
       xcgg(n) = WGAMMA(xcge(n))
    ENDDO
    xogg2 = 1.0_r8/xcgg(2)

    xobmr = 1.0_r8/xbm_r
    xoams = 1.0_r8/xam_s
    xobms = 1.0_r8/xbm_s
    xocms = xoams**xobms
    xoamg = 1.0_r8/xam_g
    xobmg = 1.0_r8/xbm_g
    xocmg = xoamg**xobmg


  END SUBROUTINE radar_init

  !+---+-----------------------------------------------------------------+
  !+---+-----------------------------------------------------------------+

  COMPLEX*16 FUNCTION m_complex_water_ray(lambda,T)

    !      Complex refractive Index of Water as function of Temperature T
    !      [deg C] and radar wavelength lambda [m]; valid for
    !      lambda in [0.001,1.0] m; T in [-10.0,30.0] deg C
    !      after Ray (1972)

    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN):: T,lambda
    REAL(KIND=r8):: epsinf,epss,epsr,epsi
    REAL(KIND=r8):: alpha,lambdas,nenner
!    COMPLEX*16, PARAMETER:: i = (0d0,1d0)
    REAL(KIND=r8), PARAMETER:: PIx=3.1415926535897932384626434_r8

    epsinf  = 5.27137_r8 + 0.02164740_r8 * T - 0.00131198_r8 * T*T
    epss    = 78.54e+0_r8 * (1.0_r8 - 4.579e-3_r8 * (T - 25.0_r8)                 &
         + 1.190e-5_r8 * (T - 25.0_r8)*(T - 25.0_r8)                        &
         - 2.800e-8_r8 * (T - 25.0_r8)*(T - 25.0_r8)*(T - 25.0_r8))
    alpha   = -16.8129_r8/(T+273.16_r8) + 0.0609265_r8
    lambdas = 0.00033836_r8 * EXP(2513.98_r8/(T+273.16_r8)) * 1e-2

    nenner = 1.e0_r8+2.e0_r8*(lambdas/lambda)**(1e0_r8-alpha)*SIN(alpha*PIx*0.5_r8) &
         + (lambdas/lambda)**(2e0_r8-2e0_r8*alpha)
    epsr = epsinf + ((epss-epsinf) * ((lambdas/lambda)**(1e0_r8-alpha)   &
         * SIN(alpha*PIx*0.5_r8)+1e0_r8)) / nenner
    epsi = ((epss-epsinf) * ((lambdas/lambda)**(1e0_r8-alpha)            &
         * COS(alpha*PIx*0.5_r8)+0e0_r8)) / nenner                           &
         + lambda*1.25664_r8/1.88496_r8

    m_complex_water_ray = SQRT(CMPLX(epsr,-epsi))

  END FUNCTION m_complex_water_ray

  !+---+-----------------------------------------------------------------+

  COMPLEX*16 FUNCTION m_complex_ice_maetzler(lambda,T)

    !      complex refractive index of ice as function of Temperature T
    !      [deg C] and radar wavelength lambda [m]; valid for
    !      lambda in [0.0001,30] m; T in [-250.0,0.0] C
    !      Original comment from the Matlab-routine of Prof. Maetzler:
    !      Function for calculating the relative permittivity of pure ice in
    !      the microwave region, according to C. Maetzler, "Microwave
    !      properties of ice and snow", in B. Schmitt et al. (eds.) Solar
    !      System Ices, Astrophys. and Space Sci. Library, Vol. 227, Kluwer
    !      Academic Publishers, Dordrecht, pp. 241-257 (1998). Input:
    !      TK = temperature (K), range 20 to 273.15
    !      f = frequency in GHz, range 0.01 to 3000

    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN):: T,lambda
    REAL(KIND=r8):: f,c,TK,B1,B2,b,deltabeta,betam,beta,theta,alfa

    c = 2.99e8_r8
    TK = T + 273.16_r8
    f = c / lambda * 1e-9_r8

    B1 = 0.0207_r8
    B2 = 1.16e-11_r8
    b = 335.0e0_r8
    deltabeta = EXP(-10.02_r8 + 0.0364_r8*(TK-273.16_r8))
    betam = (B1/TK) * ( EXP(b/TK) / ((EXP(b/TK)-1)**2) ) + B2*f*f
    beta = betam + deltabeta
    theta = 300.0_r8 / TK - 1.0_r8
    alfa = (0.00504_r8 + 0.0062_r8*theta) * EXP(-22.1_r8*theta)
    m_complex_ice_maetzler = 3.1884_r8 + 9.1e-4_r8*(TK-273.16_r8)
    m_complex_ice_maetzler = m_complex_ice_maetzler                   &
         + CMPLX(0.0_r8, (alfa/f + beta*f)) 
    m_complex_ice_maetzler = SQRT(CONJG(m_complex_ice_maetzler))

  END FUNCTION m_complex_ice_maetzler

  !+---+-----------------------------------------------------------------+

  SUBROUTINE rayleigh_soak_wetgraupel (x_g, a_geo, b_geo, fmelt,    &
       meltratio_outside, m_w, m_i, C_back,       &
       mixingrule,matrix,inclusion,                       &
       host,hostmatrix,hostinclusion)

    IMPLICIT NONE

    REAL(KIND=r8), INTENT(in):: x_g, a_geo, b_geo, fmelt,  &
         meltratio_outside
    REAL(KIND=r8), INTENT(out):: C_back
    COMPLEX*16, INTENT(in):: m_w, m_i
    CHARACTER(len=*), INTENT(in):: mixingrule, matrix, inclusion,     &
         host, hostmatrix, hostinclusion

    COMPLEX*16:: m_core, m_air
    REAL(KIND=r8):: D_large, D_g, rhog, x_w, fm, fmgrenz,    &
         volg, vg, volair, volice, volwater,            &
         meltratio_outside_grenz, mra
    INTEGER:: error
    REAL(KIND=r8), PARAMETER:: PIx=3.1415926535897932384626434_r8

    !     refractive index of air:
    m_air = (1.0e0_r8,0.0e0_r8)

    !     Limiting the degree of melting --- for safety: 
    fm = DMAX1(DMIN1(fmelt, 1.0e0_r8), 0.0e0_r8)
    !     Limiting the ratio of (melting on outside)/(melting on inside):
    mra = DMAX1(DMIN1(meltratio_outside, 1.0e0_r8), 0.0e0_r8)

    !    ! The relative portion of meltwater melting at outside should increase
    !    ! from the given input value (between 0 and 1)
    !    ! to 1 as the degree of melting approaches 1,
    !    ! so that the melting particle "converges" to a water drop.
    !    ! Simplest assumption is linear:
    mra = mra + (1.0e0_r8-mra)*fm

    x_w = x_g * fm

    D_g = a_geo * x_g**b_geo

    IF (D_g .GE. 1e-12_r8) THEN

       vg = PIx/6.0_r8 * D_g**3
       rhog = DMAX1(DMIN1(x_g / vg, 900.0e0_r8), 10.0e0_r8)
       vg = x_g / rhog

       meltratio_outside_grenz = 1.0e0_r8 - rhog / 1000.0_r8

       IF (mra .LE. meltratio_outside_grenz) THEN
          !..In this case, it cannot happen that, during melting, all the
          !.. air inclusions within the ice particle get filled with
          !.. meltwater. This only happens at the end of all melting.
          volg = vg * (1.0e0_r8 - mra * fm)

       ELSE
          !..In this case, at some melting degree fm, all the air
          !.. inclusions get filled with meltwater.
          fmgrenz=(900.0_r8-rhog)/(mra*900.0_r8-rhog+900.0_r8*rhog/1000.0_r8)

          IF (fm .LE. fmgrenz) THEN
             !.. not all air pockets are filled:
             volg = (1.0_r8 - mra * fm) * vg
          ELSE
             !..all air pockets are filled with meltwater, now the
             !.. entire ice sceleton melts homogeneously:
             volg = (x_g - x_w) / 900.0_r8 + x_w / 1000.0_r8
          ENDIF

       ENDIF

       D_large  = (6.0_r8 / PIx * volg) ** (1.0_r8/3.0_r8)
       volice = (x_g - x_w) / (volg * 900.0_r8)
       volwater = x_w / (1000.0_r8 * volg)
       volair = 1.0_r8 - volice - volwater

       !..complex index of refraction for the ice-air-water mixture
       !.. of the particle:
       m_core = get_m_mix_nested (m_air, m_i, m_w, volair, volice,      &
            volwater, mixingrule, host, matrix, inclusion, &
            hostmatrix, hostinclusion, error)
       IF (error .NE. 0) THEN
          C_back = 0.0e0_r8
          RETURN
       ENDIF

       !..Rayleigh-backscattering coefficient of melting particle: 
       C_back = (ABS((m_core**2-1.0e0_r8)/(m_core**2+2.0e0_r8)))**2           &
            * PI5 * D_large**6 / lamda4

    ELSE
       C_back = 0.0e0_r8
    ENDIF

  END SUBROUTINE rayleigh_soak_wetgraupel

  !+---+-----------------------------------------------------------------+

  COMPLEX*16 FUNCTION get_m_mix_nested (m_a, m_i, m_w, volair,      &
       volice, volwater, mixingrule, host, matrix,        &
       inclusion, hostmatrix, hostinclusion, cumulerror)

    IMPLICIT NONE

    REAL(KIND=r8), INTENT(in):: volice, volair, volwater
    COMPLEX*16, INTENT(in):: m_a, m_i, m_w
    CHARACTER(len=*), INTENT(in):: mixingrule, host, matrix,          &
         inclusion, hostmatrix, hostinclusion
    INTEGER, INTENT(out):: cumulerror

    REAL(KIND=r8):: vol1, vol2
    COMPLEX*16:: mtmp
    INTEGER:: error

    !..Folded: ( (m1 + m2) + m3), where m1,m2,m3 could each be
    !.. air, ice, or water

    cumulerror = 0
    get_m_mix_nested = CMPLX(1.0e0_r8,0.0e0_r8)

    IF (host .EQ. 'air') THEN

       IF (matrix .EQ. 'air') THEN
          WRITE(radar_debug,*) 'GET_M_MIX_NESTED: bad matrix: ', matrix
          CALL wrf_debug_radar('...get_m_mix_nested...', radar_debug)
          cumulerror = cumulerror + 1
       ELSE
          vol1 = volice / MAX(volice+volwater,1e-10_r8)
          vol2 = 1.0e0_r8 - vol1
          mtmp = get_m_mix (m_a, m_i, m_w, 0.0e0_r8, vol1, vol2,             &
               mixingrule, matrix, inclusion, error)
          cumulerror = cumulerror + error

          IF (hostmatrix .EQ. 'air') THEN
             get_m_mix_nested = get_m_mix (m_a, mtmp, 2.0_r8*m_a,              &
                  volair, (1.0e0_r8-volair), 0.0e0_r8, mixingrule,     &
                  hostmatrix, hostinclusion, error)
             cumulerror = cumulerror + error
          ELSEIF (hostmatrix .EQ. 'icewater') THEN
             get_m_mix_nested = get_m_mix (m_a, mtmp, 2.0_r8*m_a,              &
                  volair, (1.0e0_r8-volair), 0.0e0_r8, mixingrule,     &
                  'ice', hostinclusion, error)
             cumulerror = cumulerror + error
          ELSE
             WRITE(radar_debug,*) 'GET_M_MIX_NESTED: bad hostmatrix: ',        &
                  hostmatrix
             CALL wrf_debug_radar('...get_m_mix_nested...', radar_debug)
             cumulerror = cumulerror + 1
          ENDIF
       ENDIF

    ELSEIF (host .EQ. 'ice') THEN

       IF (matrix .EQ. 'ice') THEN
          WRITE(radar_debug,*) 'GET_M_MIX_NESTED: bad matrix: ', matrix
          CALL wrf_debug_radar('...get_m_mix_nested...', radar_debug)
          cumulerror = cumulerror + 1
       ELSE
          vol1 = volair / MAX(volair+volwater,1e-10_r8)
          vol2 = 1.0e0_r8 - vol1
          mtmp = get_m_mix (m_a, m_i, m_w, vol1, 0.0e0_r8, vol2,             &
               mixingrule, matrix, inclusion, error)
          cumulerror = cumulerror + error

          IF (hostmatrix .EQ. 'ice') THEN
             get_m_mix_nested = get_m_mix (mtmp, m_i, 2.0_r8*m_a,              &
                  (1.0e0_r8-volice), volice, 0.0e0_r8, mixingrule,     &
                  hostmatrix, hostinclusion, error)
             cumulerror = cumulerror + error
          ELSEIF (hostmatrix .EQ. 'airwater') THEN
             get_m_mix_nested = get_m_mix (mtmp, m_i, 2.0_r8*m_a,              &
                  (1.0e0_r8-volice), volice, 0.0e0_r8, mixingrule,     &
                  'air', hostinclusion, error)
             cumulerror = cumulerror + error          
          ELSE
             WRITE(radar_debug,*) 'GET_M_MIX_NESTED: bad hostmatrix: ',        &
                  hostmatrix
             CALL wrf_debug_radar('...get_m_mix_nested...', radar_debug)
             cumulerror = cumulerror + 1
          ENDIF
       ENDIF

    ELSEIF (host .EQ. 'water') THEN

       IF (matrix .EQ. 'water') THEN
          WRITE(radar_debug,*) 'GET_M_MIX_NESTED: bad matrix: ', matrix
          CALL wrf_debug_radar('...get_m_mix_nested...', radar_debug)
          cumulerror = cumulerror + 1
       ELSE
          vol1 = volair / MAX(volice+volair,1e-10_r8)
          vol2 = 1.0e0_r8 - vol1
          mtmp = get_m_mix (m_a, m_i, m_w, vol1, vol2, 0.0e0_r8,             &
               mixingrule, matrix, inclusion, error)
          cumulerror = cumulerror + error

          IF (hostmatrix .EQ. 'water') THEN
             get_m_mix_nested = get_m_mix (2*m_a, mtmp, m_w,                &
                  0.0e0_r8, (1.0e0_r8-volwater), volwater, mixingrule, &
                  hostmatrix, hostinclusion, error)
             cumulerror = cumulerror + error
          ELSEIF (hostmatrix .EQ. 'airice') THEN
             get_m_mix_nested = get_m_mix (2*m_a, mtmp, m_w,                &
                  0.0e0_r8, (1.0e0_r8-volwater), volwater, mixingrule, &
                  'ice', hostinclusion, error)
             cumulerror = cumulerror + error          
          ELSE
             WRITE(radar_debug,*) 'GET_M_MIX_NESTED: bad hostmatrix: ',         &
                  hostmatrix
             CALL wrf_debug_radar('...get_m_mix_nested...', radar_debug)
             cumulerror = cumulerror + 1
          ENDIF
       ENDIF

    ELSEIF (host .EQ. 'none') THEN

       get_m_mix_nested = get_m_mix (m_a, m_i, m_w,                     &
            volair, volice, volwater, mixingrule,            &
            matrix, inclusion, error)
       cumulerror = cumulerror + error

    ELSE
       WRITE(radar_debug,*) 'GET_M_MIX_NESTED: unknown matrix: ', host
       CALL wrf_debug_radar('...get_m_mix_nested...', radar_debug)
       cumulerror = cumulerror + 1
    ENDIF

    IF (cumulerror .NE. 0) THEN
       WRITE(radar_debug,*) 'GET_M_MIX_NESTED: error encountered'
       CALL wrf_debug_radar('...get_m_mix_nested...', radar_debug)
       get_m_mix_nested = CMPLX(1.0e0_r8,0.0e0_r8)    
    ENDIF

  END FUNCTION get_m_mix_nested

  !+---+-----------------------------------------------------------------+

  COMPLEX*16 FUNCTION get_m_mix (m_a, m_i, m_w, volair, volice,     &
       volwater, mixingrule, matrix, inclusion, error)

    IMPLICIT NONE

    REAL(KIND=r8), INTENT(in):: volice, volair, volwater
    COMPLEX*16, INTENT(in):: m_a, m_i, m_w
    CHARACTER(len=*), INTENT(in):: mixingrule, matrix, inclusion
    INTEGER, INTENT(out):: error

    error = 0
    get_m_mix = CMPLX(1.0e0_r8,0.0e0_r8)

    IF (mixingrule .EQ. 'maxwellgarnett') THEN
       IF (matrix .EQ. 'ice') THEN
          get_m_mix = m_complex_maxwellgarnett(volice, volair, volwater,  &
               m_i, m_a, m_w, inclusion, error)
       ELSEIF (matrix .EQ. 'water') THEN
          get_m_mix = m_complex_maxwellgarnett(volwater, volair, volice,  &
               m_w, m_a, m_i, inclusion, error)
       ELSEIF (matrix .EQ. 'air') THEN
          get_m_mix = m_complex_maxwellgarnett(volair, volwater, volice,  &
               m_a, m_w, m_i, inclusion, error)
       ELSE
          WRITE(radar_debug,*) 'GET_M_MIX: unknown matrix: ', matrix
          CALL wrf_debug_radar('...get_m_mix...', radar_debug)
          error = 1
       ENDIF

    ELSE
       WRITE(radar_debug,*) 'GET_M_MIX: unknown mixingrule: ', mixingrule
       CALL wrf_debug_radar('...get_m_mix...', radar_debug)
       error = 2
    ENDIF

    IF (error .NE. 0) THEN
       WRITE(radar_debug,*) 'GET_M_MIX: error encountered'
       CALL wrf_debug_radar('...get_m_mix...', radar_debug)
    ENDIF

  END FUNCTION get_m_mix

  !+---+-----------------------------------------------------------------+

  COMPLEX*16 FUNCTION m_complex_maxwellgarnett(vol1, vol2, vol3,    &
       m1, m2, m3, inclusion, error)

    IMPLICIT NONE

    COMPLEX*16 :: m1, m2, m3
    REAL(KIND=r8) :: vol1, vol2, vol3
    CHARACTER(len=*) :: inclusion

    COMPLEX*16 :: beta2, beta3, m1t, m2t, m3t
    INTEGER, INTENT(out) :: error

    error = 0

    IF (DABS(vol1+vol2+vol3-1.0e0_r8) .GT. 1e-6_r8) THEN
       WRITE(radar_debug,*) 'M_COMPLEX_MAXWELLGARNETT: sum of the ',       &
            'partial volume fractions is not 1...ERROR'
       CALL wrf_debug_radar('...m_complex_maxwellgarnett...', radar_debug)
       m_complex_maxwellgarnett=CMPLX(-999.99e0_r8,-999.99e0_r8)
       error = 1
       RETURN
    ENDIF

    m1t = m1**2
    m2t = m2**2
    m3t = m3**2

    IF (inclusion .EQ. 'spherical') THEN
       beta2 = 3.0e0_r8*m1t/(m2t+2.0e0_r8*m1t)
       beta3 = 3.0e0_r8*m1t/(m3t+2.0e0_r8*m1t)
    ELSEIF (inclusion .EQ. 'spheroidal') THEN
       beta2 = 2.0e0_r8*m1t/(m2t-m1t) * (m2t/(m2t-m1t)*LOG(m2t/m1t)-1.0e0_r8)
       beta3 = 2.0e0_r8*m1t/(m3t-m1t) * (m3t/(m3t-m1t)*LOG(m3t/m1t)-1.0e0_r8)
    ELSE
       WRITE(radar_debug,*) 'M_COMPLEX_MAXWELLGARNETT: ',                  &
            'unknown inclusion: ', inclusion
       CALL wrf_debug_radar('...m_complex_maxwellgarnett...', radar_debug)
       m_complex_maxwellgarnett=DCMPLX(-999.99e0_r8,-999.99e0_r8)
       error = 1
       RETURN
    ENDIF

    m_complex_maxwellgarnett = &
         SQRT(((1.0e0_r8-vol2-vol3)*m1t + vol2*beta2*m2t + vol3*beta3*m3t) / &
         (1.0e0_r8-vol2-vol3+vol2*beta2+vol3*beta3))

  END FUNCTION m_complex_maxwellgarnett

  !+---+-----------------------------------------------------------------+
  REAL(KIND=r8) FUNCTION GAMMLN(XX)
    !     --- RETURNS THE VALUE LN(GAMMA(XX)) FOR XX > 0.
    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN):: XX
    REAL(KIND=r8), PARAMETER:: STP = 2.5066282746310005_r8
    REAL(KIND=r8), DIMENSION(6), PARAMETER:: &
         COF = (/76.18009172947146_r8, -86.50532032941677_r8, &
         24.01409824083091_r8, -1.231739572450155_r8, &
         0.1208650973866179e-2_r8, -.5395239384953e-5_r8/)
    REAL(KIND=r8):: SER,TMP,X,Y
    INTEGER:: J

    X=XX
    Y=X
    TMP=X+5.5e0_r8
    TMP=(X+0.5e0_r8)*LOG(TMP)-TMP
    SER=1.000000000190015e0_r8
    DO  J=1,6
       Y=Y+1.e0_r8
       SER=SER+COF(J)/Y
       !11     CONTINUE
    END DO
    GAMMLN=TMP+LOG(STP*SER/X)
  END FUNCTION GAMMLN
  !  (C) Copr. 1986-92 Numerical Recipes Software 2.02
  !+---+-----------------------------------------------------------------+
  REAL(KIND=r8)  FUNCTION WGAMMA(y)

    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN):: y

    WGAMMA = EXP(GAMMLN(y))

  END FUNCTION WGAMMA

  !+---+-----------------------------------------------------------------+
  !wrf_debug_radar
  !+---+-----------------------------------------------------------------+


  SUBROUTINE wrf_debug_radar(str,str2)
    IMPLICIT NONE
    CHARACTER(LEN=*) , INTENT(IN   ) :: str
    CHARACTER(LEN=*) , INTENT(IN   ) :: str2

    CALL MsgOne(str,str2)
  END  SUBROUTINE wrf_debug_radar



  !+---+-----------------------------------------------------------------+
END MODULE module_mp_radar
!+---+-----------------------------------------------------------------+


!+---+-----------------------------------------------------------------+
!.. This subroutine computes the moisture tendencies of water vapor,
!.. cloud droplets, rain, cloud ice (pristine), snow, and graupel.
!.. Prior to WRFv2.2 this code was based on Reisner et al (1998), but
!.. few of those pieces remain.  A complete description is now found in
!.. Thompson, G., P. R. Field, R. M. Rasmussen, and W. D. Hall, 2008:
!.. Explicit Forecasts of winter precipitation using an improved bulk
!.. microphysics scheme. Part II: Implementation of a new snow
!.. parameterization.  Mon. Wea. Rev., 136, 5095-5115.
!.. Prior to WRFv3.1, this code was single-moment rain prediction as
!.. described in the reference above, but in v3.1 and higher, the
!.. scheme is two-moment rain (predicted rain number concentration).
!..
!.. Beginning with WRFv3.6, this is also the "aerosol-aware" scheme as
!.. described in Thompson, G. and T. Eidhammer, 2014:  A study of
!.. aerosol impacts on clouds and precipitation development in a large
!.. winter cyclone.  J. Atmos. Sci., 71, 3636-3658.  Setting WRF
!.. namelist option mp_physics=8 utilizes the older one-moment cloud
!.. water with constant droplet concentration set as Nt_c (found below)
!.. while mp_physics=28 uses double-moment cloud droplet number
!.. concentration, which is not permitted to exceed Nt_c_max below.
!..
!.. Most importantly, users may wish to modify the prescribed number of
!.. cloud droplets (Nt_c; see guidelines mentioned below).  Otherwise,
!.. users may alter the rain and graupel size distribution parameters
!.. to use exponential (Marshal-Palmer) or generalized gamma shape.
!.. The snow field assumes a combination of two gamma functions (from
!.. Field et al. 2005) and would require significant modifications
!.. throughout the entire code to alter its shape as well as accretion
!.. rates.  Users may also alter the constants used for density of rain,
!.. graupel, ice, and snow, but the latter is not constant when using
!.. Paul Field's snow distribution and moments methods.  Other values
!.. users can modify include the constants for mass and/or velocity
!.. power law relations and assumed capacitances used in deposition/
!.. sublimation/evaporation/melting.
!.. Remaining values should probably be left alone.
!..
!..Author: Greg Thompson, NCAR-RAL, gthompsn@ucar.edu, 303-497-2805
!..Last modified: 11 Feb 2015   Aerosol additions to v3.5.1 code 9/2013
!..                 Cloud fraction additions 11/2014 part of pre-v3.7
!+---+-----------------------------------------------------------------+
!wrft:model_layer:physics
!+---+-----------------------------------------------------------------+
!
MODULE Micro_GTHOMPSON
  USE Parallelism,Only: MsgOne,FatalError
  !      USE module_wrf_error
  USE module_mp_radar
  !#if ( defined( DM_PARALLEL ) && ( ! defined( STUBMPI ) ) )
  !      USE module_dm, ONLY : wrf_dm_max_real
  !#endif

  IMPLICIT NONE
  PRIVATE

  INTEGER      , PARAMETER :: r8  = SELECTED_REAL_KIND(15)

  LOGICAL, PARAMETER, PRIVATE:: iiwarm = .FALSE.
  LOGICAL, PRIVATE           :: is_aerosol_aware = .TRUE.
  LOGICAL, PARAMETER, PRIVATE:: dustyIce = .TRUE.
  LOGICAL, PARAMETER, PRIVATE:: homogIce = .TRUE.

  INTEGER, PARAMETER, PRIVATE:: IFDRY = 0
  REAL(KIND=r8), PARAMETER, PRIVATE:: T_0 = 273.15_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: PI = 3.1415926536_r8

  REAL(R8),PARAMETER :: SHR_CONST_MWDAIR = 28.966_R8       ! molecular weight dry air ~ kg/kmole  
  REAL(r8),PARAMETER :: SHR_CONST_MWWV   = 18.016_r8       ! molecular weight water vapor
  REAL(R8),PARAMETER :: SHR_CONST_AVOGAD = 6.02214e26_R8   ! Avogadro's number ~ molecules/kmole  
  REAL(R8),PARAMETER :: SHR_CONST_BOLTZ  = 1.38065e-23_R8  ! Boltzmann's constant ~ J/K/molecule
  REAL(R8),PARAMETER :: SHR_CONST_G      = 9.80616_R8      ! acceleration of gravity ~ m/s^2

  REAL(R8),PARAMETER :: SHR_CONST_RGAS   = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ ! Universal gas constant ~ J/K/kmole

  REAL(R8),PARAMETER :: SHR_CONST_RDAIR  = SHR_CONST_RGAS/SHR_CONST_MWDAIR  ! Dry air gas constant ~ J/K/kg
  REAL(r8),PARAMETER :: SHR_CONST_RWV    = SHR_CONST_RGAS/SHR_CONST_MWWV    ! Water vapor gas constant ~ J/K/kg

  REAL(r8),PARAMETER :: rair   = SHR_CONST_RDAIR    ! Gas constant for dry air (J/K/kg)
  REAL(r8),PARAMETER :: gravit = SHR_CONST_G      ! gravitational acceleration
  REAL(r8),PARAMETER :: zvir   = SHR_CONST_RWV/SHR_CONST_RDAIR - 1          ! rh2o/rair - 1
  !LOGICAL, PARAMETER :: f_qndrop=.FALSE.

!  REAL(KIND=r8), PARAMETER, PRIVATE :: rair=
  !..Densities of rain, snow, graupel, and cloud ice.
  REAL(KIND=r8), PARAMETER, PRIVATE:: rho_w = 1000.0_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: rho_s = 100.0_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: rho_g = 500.0_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: rho_i = 890.0_r8


  !..Prescribed number of cloud droplets.  Set according to known data or
  !.. roughly 100 per cc (100.E6 m^-3) for Maritime cases and
  !.. 300 per cc (300.E6 m^-3) for Continental.  Gamma shape parameter,
  !.. mu_c, calculated based on Nt_c is important in autoconversion
  !.. scheme.  In 2-moment cloud water, Nt_c represents a maximum of
  !.. droplet concentration and nu_c is also variable depending on local
  !.. droplet number concentration.
  REAL(KIND=r8), PARAMETER, PRIVATE:: Nt_c = 100.E6_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: Nt_c_max = 1999.E6_r8

  !..Declaration of constants for assumed CCN/IN aerosols when none in
  !.. the input data.  Look inside the init routine for modifications
  !.. due to surface land-sea points or vegetation characteristics.
  REAL(KIND=r8), PARAMETER, PRIVATE:: naIN0 = 1.5E6_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: naIN1 = 0.5E6_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: naCCN0 = 300.0E6_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: naCCN1 = 50.0E6_r8

  !..Generalized gamma distributions for rain, graupel and cloud ice.
  !.. N(D) = N_0 * D**mu * exp(-lamda*D);  mu=0 is exponential.
  REAL(KIND=r8), PARAMETER, PRIVATE:: mu_r = 0.0_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: mu_g = 0.0_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: mu_i = 0.0_r8
  REAL(KIND=r8), PRIVATE:: mu_c

  !..Sum of two gamma distrib for snow (Field et al. 2005).
  !.. N(D) = M2**4/M3**3 * [Kap0*exp(-M2*Lam0*D/M3)
  !..    + Kap1*(M2/M3)**mu_s * D**mu_s * exp(-M2*Lam1*D/M3)]
  !.. M2 and M3 are the (bm_s)th and (bm_s+1)th moments respectively
  !.. calculated as function of ice water content and temperature.
  REAL(KIND=r8), PARAMETER, PRIVATE:: mu_s = 0.6357_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: Kap0 = 490.6_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: Kap1 = 17.46_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: Lam0 = 20.78_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: Lam1 = 3.29_r8

  !..Y-intercept parameter for graupel is not constant and depends on
  !.. mixing ratio.  Also, when mu_g is non-zero, these become equiv
  !.. y-intercept for an exponential distrib and proper values are
  !.. computed based on same mixing ratio and total number concentration.
  REAL(KIND=r8), PARAMETER, PRIVATE:: gonv_min = 1.E4_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: gonv_max = 3.E6_r8

  !..Mass power law relations:  mass = am*D**bm
  !.. Snow from Field et al. (2005), others assume spherical form.
  REAL(KIND=r8), PARAMETER, PRIVATE:: am_r = PI*rho_w/6.0_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: bm_r = 3.0_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: am_s = 0.069_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: bm_s = 2.0_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: am_g = PI*rho_g/6.0_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: bm_g = 3.0_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: am_i = PI*rho_i/6.0_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: bm_i = 3.0_r8

  !..Fallspeed power laws relations:  v = (av*D**bv)*exp(-fv*D)
  !.. Rain from Ferrier (1994), ice, snow, and graupel from
  !.. Thompson et al (2008). Coefficient fv is zero for graupel/ice.
  REAL(KIND=r8), PARAMETER, PRIVATE:: av_r = 4854.0_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: bv_r = 1.0_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: fv_r = 195.0_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: av_s = 40.0_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: bv_s = 0.55_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: fv_s = 100.0_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: av_g = 442.0_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: bv_g = 0.89_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: av_i = 1847.5_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: bv_i = 1.0_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: av_c = 0.316946E8_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: bv_c = 2.0_r8

  !..Capacitance of sphere and plates/aggregates: D**3, D**2
  REAL(KIND=r8), PARAMETER, PRIVATE:: C_cube = 0.5_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: C_sqrd = 0.3_r8

  !..Collection efficiencies.  Rain/snow/graupel collection of cloud
  !.. droplets use variables (Ef_rw, Ef_sw, Ef_gw respectively) and
  !.. get computed elsewhere because they are dependent on stokes
  !.. number.
  REAL(KIND=r8), PARAMETER, PRIVATE:: Ef_si = 0.05_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: Ef_rs = 0.95_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: Ef_rg = 0.75_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: Ef_ri = 0.95_r8

  !..Minimum microphys values
  !.. R1 value, 1.E-12, cannot be set lower because of numerical
  !.. problems with Paul Field's moments and should not be set larger
  !.. because of truncation problems in snow/ice growth.
  REAL(KIND=r8), PARAMETER, PRIVATE:: R1 = 1.E-12_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: R2 = 1.E-6_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: eps = 1.E-15_r8

  !..Constants in Cooper curve relation for cloud ice number.
  REAL(KIND=r8), PARAMETER, PRIVATE:: TNO = 5.0_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: ATO = 0.304_r8

  !..Rho_not used in fallspeed relations (rho_not/rho)**.5 adjustment.
  REAL(KIND=r8), PARAMETER, PRIVATE:: rho_not = 101325.0_r8/(287.05_r8*298.0_r8)

  !..Schmidt number
  REAL(KIND=r8), PARAMETER, PRIVATE:: Sc = 0.632_r8
  REAL(KIND=r8), PRIVATE:: Sc3

  !..Homogeneous freezing temperature
  REAL(KIND=r8), PARAMETER, PRIVATE:: HGFR = 235.16_r8

  !..Water vapor and air gas constants at constant pressure
  REAL(KIND=r8), PARAMETER, PRIVATE:: Rv = 461.5_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: oRv = 1.0_r8/Rv
  REAL(KIND=r8), PARAMETER, PRIVATE:: R = 287.04_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: Cp = 1004.0_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: R_uni = 8.314_r8                           ! J (mol K)-1

  REAL(KIND=r8), PARAMETER, PRIVATE:: k_b = 1.38065E-23_r8           ! Boltzmann constant [J/K]
  REAL(KIND=r8), PARAMETER, PRIVATE:: M_w = 18.01528E-3_r8           ! molecular mass of water [kg/mol]
  REAL(KIND=r8), PARAMETER, PRIVATE:: M_a = 28.96E-3_r8              ! molecular mass of air [kg/mol]
  REAL(KIND=r8), PARAMETER, PRIVATE:: N_avo = 6.022E23_r8            ! Avogadro number [1/mol]
  REAL(KIND=r8), PARAMETER, PRIVATE:: ma_w = M_w / N_avo          ! mass of water molecule [kg]
  REAL(KIND=r8), PARAMETER, PRIVATE:: ar_volume = 4.0_r8/3.0_r8*PI*(2.5e-6_r8)**3        ! assume radius of 0.025 micrometer, 2.5e-6 cm

  !..Enthalpy of sublimation, vaporization, and fusion at 0C.
  REAL(KIND=r8), PARAMETER, PRIVATE:: lsub = 2.834E6_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: lvap0 = 2.5E6_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: lfus = lsub - lvap0
  REAL(KIND=r8), PARAMETER, PRIVATE:: olfus = 1.0_r8/lfus

  !..Ice initiates with this mass (kg), corresponding diameter calc.
  !..Min diameters and mass of cloud, rain, snow, and graupel (m, kg).
  REAL(KIND=r8), PARAMETER, PRIVATE:: xm0i = 1.E-12_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: D0c = 1.E-6_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: D0r = 50.E-6_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: D0s = 200.E-6_r8
  REAL(KIND=r8), PARAMETER, PRIVATE:: D0g = 250.E-6_r8
  REAL(KIND=r8), PRIVATE:: D0i, xm0s, xm0g

  !..Lookup table dimensions
  INTEGER, PARAMETER, PRIVATE:: nbins = 100
  INTEGER, PARAMETER, PRIVATE:: nbc = nbins
  INTEGER, PARAMETER, PRIVATE:: nbi = nbins
  INTEGER, PARAMETER, PRIVATE:: nbr = nbins
  INTEGER, PARAMETER, PRIVATE:: nbs = nbins
  INTEGER, PARAMETER, PRIVATE:: nbg = nbins
  INTEGER, PARAMETER, PRIVATE:: ntb_c = 37
  INTEGER, PARAMETER, PRIVATE:: ntb_i = 64
  INTEGER, PARAMETER, PRIVATE:: ntb_r = 37
  INTEGER, PARAMETER, PRIVATE:: ntb_s = 28
  INTEGER, PARAMETER, PRIVATE:: ntb_g = 28
  INTEGER, PARAMETER, PRIVATE:: ntb_g1 = 28
  INTEGER, PARAMETER, PRIVATE:: ntb_r1 = 37
  INTEGER, PARAMETER, PRIVATE:: ntb_i1 = 55
  INTEGER, PARAMETER, PRIVATE:: ntb_t = 9
  INTEGER, PRIVATE:: nic1, nic2, nii2, nii3, nir2, nir3, nis2, nig2, nig3
  INTEGER, PARAMETER, PRIVATE:: ntb_arc = 7
  INTEGER, PARAMETER, PRIVATE:: ntb_arw = 9
  INTEGER, PARAMETER, PRIVATE:: ntb_art = 7
  INTEGER, PARAMETER, PRIVATE:: ntb_arr = 5
  INTEGER, PARAMETER, PRIVATE:: ntb_ark = 4
  INTEGER, PARAMETER, PRIVATE:: ntb_IN = 55
  INTEGER, PRIVATE:: niIN2

  REAL(KIND=r8), DIMENSION(nbins+1):: xDx
  REAL(KIND=r8), DIMENSION(nbc):: Dc, dtc
  REAL(KIND=r8), DIMENSION(nbi):: Di, dti
  REAL(KIND=r8), DIMENSION(nbr):: Dr, dtr
  REAL(KIND=r8), DIMENSION(nbs):: Ds, dts
  REAL(KIND=r8), DIMENSION(nbg):: Dg, dtg
  REAL(KIND=r8), DIMENSION(nbc):: t_Nc

  !..Lookup tables for cloud water content (kg/m**3).
  REAL(KIND=r8), DIMENSION(ntb_c), PARAMETER, PRIVATE:: &
       r_c = (/1.e-6_r8,2.e-6_r8,3.e-6_r8,4.e-6_r8,5.e-6_r8,6.e-6_r8,7.e-6_r8,8.e-6_r8,9.e-6_r8, &
               1.e-5_r8,2.e-5_r8,3.e-5_r8,4.e-5_r8,5.e-5_r8,6.e-5_r8,7.e-5_r8,8.e-5_r8,9.e-5_r8, &
               1.e-4_r8,2.e-4_r8,3.e-4_r8,4.e-4_r8,5.e-4_r8,6.e-4_r8,7.e-4_r8,8.e-4_r8,9.e-4_r8, &
               1.e-3_r8,2.e-3_r8,3.e-3_r8,4.e-3_r8,5.e-3_r8,6.e-3_r8,7.e-3_r8,8.e-3_r8,9.e-3_r8, &
               1.e-2_r8/)

  !..Lookup tables for cloud ice content (kg/m**3).
  REAL(KIND=r8), DIMENSION(ntb_i), PARAMETER, PRIVATE:: &
       r_i = (/1.e-10_r8,2.e-10_r8,3.e-10_r8,4.e-10_r8, &
       5.e-10_r8,6.e-10_r8,7.e-10_r8,8.e-10_r8,9.e-10_r8, &
       1.e-9_r8,2.e-9_r8,3.e-9_r8,4.e-9_r8,5.e-9_r8,6.e-9_r8,7.e-9_r8,8.e-9_r8,9.e-9_r8, &
       1.e-8_r8,2.e-8_r8,3.e-8_r8,4.e-8_r8,5.e-8_r8,6.e-8_r8,7.e-8_r8,8.e-8_r8,9.e-8_r8, &
       1.e-7_r8,2.e-7_r8,3.e-7_r8,4.e-7_r8,5.e-7_r8,6.e-7_r8,7.e-7_r8,8.e-7_r8,9.e-7_r8, &
       1.e-6_r8,2.e-6_r8,3.e-6_r8,4.e-6_r8,5.e-6_r8,6.e-6_r8,7.e-6_r8,8.e-6_r8,9.e-6_r8, &
       1.e-5_r8,2.e-5_r8,3.e-5_r8,4.e-5_r8,5.e-5_r8,6.e-5_r8,7.e-5_r8,8.e-5_r8,9.e-5_r8, &
       1.e-4_r8,2.e-4_r8,3.e-4_r8,4.e-4_r8,5.e-4_r8,6.e-4_r8,7.e-4_r8,8.e-4_r8,9.e-4_r8, &
       1.e-3_r8/)

  !..Lookup tables for rain content (kg/m**3).
  REAL(KIND=r8), DIMENSION(ntb_r), PARAMETER, PRIVATE:: &
       r_r = (/1.e-6_r8,2.e-6_r8,3.e-6_r8,4.e-6_r8,5.e-6_r8,6.e-6_r8,7.e-6_r8,8.e-6_r8,9.e-6_r8, &
       1.e-5_r8,2.e-5_r8,3.e-5_r8,4.e-5_r8,5.e-5_r8,6.e-5_r8,7.e-5_r8,8.e-5_r8,9.e-5_r8, &
       1.e-4_r8,2.e-4_r8,3.e-4_r8,4.e-4_r8,5.e-4_r8,6.e-4_r8,7.e-4_r8,8.e-4_r8,9.e-4_r8, &
       1.e-3_r8,2.e-3_r8,3.e-3_r8,4.e-3_r8,5.e-3_r8,6.e-3_r8,7.e-3_r8,8.e-3_r8,9.e-3_r8, &
       1.e-2_r8/)

  !..Lookup tables for graupel content (kg/m**3).
  REAL(KIND=r8), DIMENSION(ntb_g), PARAMETER, PRIVATE:: &
       r_g = (/1.e-5_r8,2.e-5_r8,3.e-5_r8,4.e-5_r8,5.e-5_r8,6.e-5_r8,7.e-5_r8,8.e-5_r8,9.e-5_r8, &
       1.e-4_r8,2.e-4_r8,3.e-4_r8,4.e-4_r8,5.e-4_r8,6.e-4_r8,7.e-4_r8,8.e-4_r8,9.e-4_r8, &
       1.e-3_r8,2.e-3_r8,3.e-3_r8,4.e-3_r8,5.e-3_r8,6.e-3_r8,7.e-3_r8,8.e-3_r8,9.e-3_r8, &
       1.e-2_r8/)

  !..Lookup tables for snow content (kg/m**3).
  REAL(KIND=r8), DIMENSION(ntb_s), PARAMETER, PRIVATE:: &
       r_s = (/1.e-5_r8,2.e-5_r8,3.e-5_r8,4.e-5_r8,5.e-5_r8,6.e-5_r8,7.e-5_r8,8.e-5_r8,9.e-5_r8, &
       1.e-4_r8,2.e-4_r8,3.e-4_r8,4.e-4_r8,5.e-4_r8,6.e-4_r8,7.e-4_r8,8.e-4_r8,9.e-4_r8, &
       1.e-3_r8,2.e-3_r8,3.e-3_r8,4.e-3_r8,5.e-3_r8,6.e-3_r8,7.e-3_r8,8.e-3_r8,9.e-3_r8, &
       1.e-2_r8/)

  !..Lookup tables for rain y-intercept parameter (/m**4).
  REAL(KIND=r8), DIMENSION(ntb_r1), PARAMETER, PRIVATE:: &
       N0r_exp = (/1.e6_r8,2.e6_r8,3.e6_r8,4.e6_r8,5.e6_r8,6.e6_r8,7.e6_r8,8.e6_r8,9.e6_r8, &
       1.e7_r8,2.e7_r8,3.e7_r8,4.e7_r8,5.e7_r8,6.e7_r8,7.e7_r8,8.e7_r8,9.e7_r8, &
       1.e8_r8,2.e8_r8,3.e8_r8,4.e8_r8,5.e8_r8,6.e8_r8,7.e8_r8,8.e8_r8,9.e8_r8, &
       1.e9_r8,2.e9_r8,3.e9_r8,4.e9_r8,5.e9_r8,6.e9_r8,7.e9_r8,8.e9_r8,9.e9_r8, &
       1.e10_r8/)

  !..Lookup tables for graupel y-intercept parameter (/m**4).
  REAL(KIND=r8), DIMENSION(ntb_g1), PARAMETER, PRIVATE:: &
       N0g_exp = (/1.e4_r8,2.e4_r8,3.e4_r8,4.e4_r8,5.e4_r8,6.e4_r8,7.e4_r8,8.e4_r8,9.e4_r8, &
       1.e5_r8,2.e5_r8,3.e5_r8,4.e5_r8,5.e5_r8,6.e5_r8,7.e5_r8,8.e5_r8,9.e5_r8, &
       1.e6_r8,2.e6_r8,3.e6_r8,4.e6_r8,5.e6_r8,6.e6_r8,7.e6_r8,8.e6_r8,9.e6_r8, &
       1.e7_r8/)

  !..Lookup tables for ice number concentration (/m**3).
  REAL(KIND=r8), DIMENSION(ntb_i1), PARAMETER, PRIVATE:: &
       Nt_i = (/1.0_r8,2.0_r8,3.0_r8,4.0_r8,5.0_r8,6.0_r8,7.0_r8,8.0_r8,9.0_r8, &
       1.e1_r8,2.e1_r8,3.e1_r8,4.e1_r8,5.e1_r8,6.e1_r8,7.e1_r8,8.e1_r8,9.e1_r8, &
       1.e2_r8,2.e2_r8,3.e2_r8,4.e2_r8,5.e2_r8,6.e2_r8,7.e2_r8,8.e2_r8,9.e2_r8, &
       1.e3_r8,2.e3_r8,3.e3_r8,4.e3_r8,5.e3_r8,6.e3_r8,7.e3_r8,8.e3_r8,9.e3_r8, &
       1.e4_r8,2.e4_r8,3.e4_r8,4.e4_r8,5.e4_r8,6.e4_r8,7.e4_r8,8.e4_r8,9.e4_r8, &
       1.e5_r8,2.e5_r8,3.e5_r8,4.e5_r8,5.e5_r8,6.e5_r8,7.e5_r8,8.e5_r8,9.e5_r8, &
       1.e6_r8/)

  !..Aerosol table parameter: Number of available aerosols, vertical
  !.. velocity, temperature, aerosol mean radius, and hygroscopicity.
  REAL(KIND=r8), DIMENSION(ntb_arc), PARAMETER, PRIVATE:: &
       ta_Na = (/10.0_r8, 31.6_r8, 100.0_r8, 316.0_r8, 1000.0_r8, 3160.0_r8, 10000.0_r8/)
  REAL(KIND=r8), DIMENSION(ntb_arw), PARAMETER, PRIVATE:: &
       ta_Ww = (/0.01_r8, 0.0316_r8, 0.1_r8, 0.316_r8, 1.0_r8, 3.16_r8, 10.0_r8, 31.6_r8, 100.0_r8/)
  REAL(KIND=r8), DIMENSION(ntb_art), PARAMETER, PRIVATE:: &
       ta_Tk = (/243.15_r8, 253.15_r8, 263.15_r8, 273.15_r8, 283.15_r8, 293.15_r8, 303.15_r8/)
  REAL(KIND=r8), DIMENSION(ntb_arr), PARAMETER, PRIVATE:: &
       ta_Ra = (/0.01_r8, 0.02_r8, 0.04_r8, 0.08_r8, 0.16_r8/)
  REAL(KIND=r8), DIMENSION(ntb_ark), PARAMETER, PRIVATE:: &
       ta_Ka = (/0.2_r8, 0.4_r8, 0.6_r8, 0.8_r8/)

  !..Lookup tables for IN concentration (/m**3) from 0.001 to 1000/Liter.
  REAL(KIND=r8), DIMENSION(ntb_IN), PARAMETER, PRIVATE:: &
       Nt_IN = (/1.0_r8,2.0_r8,3.0_r8,4.0_r8,5.0_r8,6.0_r8,7.0_r8,8.0_r8,9.0_r8, &
       1.e1_r8,2.e1_r8,3.e1_r8,4.e1_r8,5.e1_r8,6.e1_r8,7.e1_r8,8.e1_r8,9.e1_r8, &
       1.e2_r8,2.e2_r8,3.e2_r8,4.e2_r8,5.e2_r8,6.e2_r8,7.e2_r8,8.e2_r8,9.e2_r8, &
       1.e3_r8,2.e3_r8,3.e3_r8,4.e3_r8,5.e3_r8,6.e3_r8,7.e3_r8,8.e3_r8,9.e3_r8, &
       1.e4_r8,2.e4_r8,3.e4_r8,4.e4_r8,5.e4_r8,6.e4_r8,7.e4_r8,8.e4_r8,9.e4_r8, &
       1.e5_r8,2.e5_r8,3.e5_r8,4.e5_r8,5.e5_r8,6.e5_r8,7.e5_r8,8.e5_r8,9.e5_r8, &
       1.e6_r8/)

  !..For snow moments conversions (from Field et al. 2005)
  REAL(KIND=r8), DIMENSION(10), PARAMETER, PRIVATE:: &
       sa = (/ 5.065339_r8, -0.062659_r8, -3.032362_r8, 0.029469_r8, -0.000285_r8, &
       0.31255_r8,   0.000204_r8,  0.003199_r8, 0.0_r8,      -0.015952_r8/)
  REAL(KIND=r8), DIMENSION(10), PARAMETER, PRIVATE:: &
       sb = (/ 0.476221_r8, -0.015896_r8,  0.165977_r8, 0.007468_r8, -0.000141_r8, &
       0.060366_r8,  0.000079_r8,  0.000594_r8, 0.0_r8,      -0.003577_r8/)

  !..Temperatures (5 C interval 0 to -40) used in lookup tables.
  REAL(KIND=r8), DIMENSION(ntb_t), PARAMETER, PRIVATE:: &
       Tc = (/-0.01_r8, -5.0_r8, -10.0_r8, -15.0_r8, -20.0_r8, -25.0_r8, -30.0_r8, -35.0_r8, -40.0_r8/)

  !..Lookup tables for various accretion/collection terms.
  !.. ntb_x refers to the number of elements for rain, snow, graupel,
  !.. and temperature array indices.  Variables beginning with t-p/c/m/n
  !.. represent lookup tables.  Save compile-time memory by making
  !.. allocatable (2009Jun12, J. Michalakes).
  INTEGER, PARAMETER, PRIVATE:: R8SIZE = 8
  INTEGER, PARAMETER, PRIVATE:: R4SIZE = 4
  REAL (KIND=R8SIZE), ALLOCATABLE, DIMENSION(:,:,:,:)::             &
       tcg_racg, tmr_racg, tcr_gacr, tmg_gacr,                 &
       tnr_racg, tnr_gacr
  REAL (KIND=R8SIZE), ALLOCATABLE, DIMENSION(:,:,:,:)::             &
       tcs_racs1, tmr_racs1, tcs_racs2, tmr_racs2,             &
       tcr_sacr1, tms_sacr1, tcr_sacr2, tms_sacr2,             &
       tnr_racs1, tnr_racs2, tnr_sacr1, tnr_sacr2
  REAL (KIND=R8SIZE), ALLOCATABLE, DIMENSION(:,:,:,:)::             &
       tpi_qcfz, tni_qcfz
  REAL (KIND=R8SIZE), ALLOCATABLE, DIMENSION(:,:,:,:)::             &
       tpi_qrfz, tpg_qrfz, tni_qrfz, tnr_qrfz
  REAL (KIND=R8SIZE), ALLOCATABLE, DIMENSION(:,:)::                 &
       tps_iaus, tni_iaus, tpi_ide
  REAL (KIND=R8SIZE), ALLOCATABLE, DIMENSION(:,:):: t_Efrw
  REAL (KIND=R8SIZE), ALLOCATABLE, DIMENSION(:,:):: t_Efsw
  REAL (KIND=R8SIZE), ALLOCATABLE, DIMENSION(:,:,:):: tnr_rev
  REAL (KIND=R8SIZE), ALLOCATABLE, DIMENSION(:,:,:)::               &
       tpc_wev, tnc_wev
  REAL (KIND=R4SIZE), ALLOCATABLE, DIMENSION(:,:,:,:,:):: tnccn_act

  !..Variables holding a bunch of exponents and gamma values (cloud water,
  !.. cloud ice, rain, snow, then graupel).
  REAL(KIND=r8), DIMENSION(5,15), PRIVATE:: cce, ccg
  REAL(KIND=r8), DIMENSION(15), PRIVATE::  ocg1, ocg2
  REAL(KIND=r8), DIMENSION(7), PRIVATE:: cie, cig
  REAL(KIND=r8), PRIVATE:: oig1, oig2, obmi
  REAL(KIND=r8), DIMENSION(13), PRIVATE:: cre, crg
  REAL(KIND=r8), PRIVATE:: ore1, org1, org2, org3, obmr
  REAL(KIND=r8), DIMENSION(18), PRIVATE:: cse, csg
  REAL(KIND=r8), PRIVATE:: oams, obms, ocms
  REAL(KIND=r8), DIMENSION(12), PRIVATE:: cge, cgg
  REAL(KIND=r8), PRIVATE:: oge1, ogg1, ogg2, ogg3, oamg, obmg, ocmg

  !..Declaration of precomputed constants in various rate eqns.
  REAL(KIND=r8) :: t1_qr_qc, t1_qr_qi, t2_qr_qi, t1_qg_qc, t1_qs_qc, t1_qs_qi
  REAL(KIND=r8) :: t1_qr_ev, t2_qr_ev
  REAL(KIND=r8) :: t1_qs_sd, t2_qs_sd, t1_qg_sd, t2_qg_sd
  REAL(KIND=r8) :: t1_qs_me, t2_qs_me, t1_qg_me, t2_qg_me
  LOGICAL       :: restart_mic
  LOGICAL       :: First_mic=.TRUE.
  INTEGER       :: nLev_mic
  !+---+
  !+---+-----------------------------------------------------------------+
  !..END DECLARATIONS
  !+---+-----------------------------------------------------------------+
  !+---+
  !ctrlL
  PUBLIC :: Init_Micro_thompson
  PUBLIC :: RunMicro_thompson
CONTAINS

  SUBROUTINE Init_Micro_thompson(kMax,path_in,restart,a_hybr,b_hybr)

    IMPLICIT NONE


    !..OPTIONAL variables that control application of aerosol-aware scheme
    CHARACTER(LEN=*), INTENT(IN   ) :: path_in
    LOGICAL         , INTENT(IN   ) :: restart
    INTEGER         , INTENT(IN   ) :: kMax
    REAL(KIND=r8)   , INTENT(IN   ) :: a_hybr   (kMax+1)
    REAL(KIND=r8)   , INTENT(IN   ) :: b_hybr   (kMax+1)
    REAL(KIND=r8) :: hypm (kMax+1)

    !
    !  LOCAL VARIABLE
    !
!PK -      CHARACTER(LEN=256) :: mp_debug
    CHARACTER(LEN=600) :: wrf_err_message
    REAL(KIND=r8) :: ps0
    INTEGER       :: i
    INTEGER       :: j
    INTEGER       :: k
    INTEGER       :: l
    INTEGER       :: m
    INTEGER       :: n
!PK -      REAL(KIND=r8) :: h_01
!PK -      REAL(KIND=r8) :: niIN3
!PK -      REAL(KIND=r8) :: niCCN3
!PK -      REAL(KIND=r8) :: max_test
    LOGICAL       :: micro_init
!PK -    LOGICAL       :: has_CCN
!PK -    LOGICAL       :: has_IN
!PK -    nLev_mic=0
!PK -    DO k=1,Size(si,1)
!PK -      nLev_mic=nLev_mic+1
!PK -      IF(si(k)*100000.0_r8 < 5000.0_r8)exit ! 50 mb
!PK -    END DO



    ! hypm     reference state midpoint pressures
    ps0    = 1.0e5_r8            ! Base state surface pressure (pascals)
    DO k=Size(a_hybr,1),1,-1
!  SB
!      hypm(k) =  ps0*sig(k)
       l = kMax+1-k  ! inversion from top to bottom to bottom to top
       hypm(k) =  0.5_r8 * (ps0*(b_hybr(l)+b_hybr(l+1)) + a_hybr(l) + &
                            a_hybr(l+1) )
    END DO

    nLev_mic = 0
    DO k=kMax,1,-1
       IF (hypm(k) >= 5.e4_r8) THEN!40000
          nLev_mic = nLev_mic + 1
       END IF
    END DO
    nLev_mic = MAX(nLev_mic,1)
    
        
!PK -    is_aerosol_aware = .FALSE.
    micro_init = .FALSE.
    restart_mic=restart
    First_mic=.true.
!PK -    has_CCN    = .FALSE.
!PK -    has_IN     = .FALSE.

!PK -    WRITE(mp_debug,*) ' DEBUG  checking column of hgt ', 1+1,1+1
!PK -    CALL wrf_debug(0, mp_debug)
!PK -    DO k = 1, kMAx
!PK -       WRITE(mp_debug,*) ' DEBUGT  k, hgt = ', k, hgt(1+1,k,1+1)
!PK -       CALL wrf_debug(0, mp_debug)
!PK -    ENDDO
!PK -
!PK -    IF (PRESENT(nwfa2d) .AND. PRESENT(nwfa) .AND. PRESENT(nifa)) is_aerosol_aware = .TRUE.
!PK -
!PK -    IF (is_aerosol_aware) THEN
!PK -
!PK -      !..Check for existing aerosol data, both CCN and IN aerosols.  If missing
!PK -      !.. fill in just a basic vertical profile, somewhat boundary-layer following.
!PK -
!PK -      !#if ( defined( DM_PARALLEL ) && ( ! defined( STUBMPI ) ) )
!PK -      !      max_test = wrf_dm_max_real ( MAXVAL(nwfa(1:ibMax-1,:,1:jbMax-1)) )
!PK -      !#else
!PK -      max_test = MAXVAL ( nwfa(1:ibMax-1,:,1:jbMax-1) )
!PK -      !#endif
!PK -
!PK -      IF (max_test .LT. eps) THEN
!PK -         WRITE(mp_debug,*) ' Apparently there are no initial CCN aerosols.'
!PK -         CALL wrf_debug(0, mp_debug)
!PK -         WRITE(mp_debug,*) '   checked column at point (i,j) = ', 1,1
!PK -         CALL wrf_debug(0, mp_debug)
!PK -         DO j = 1, jbMax
!PK -            DO i =1, ibMax
!PK -               IF (hgt(i,1,j).LE.1000.0_r8) THEN
!PK -                  h_01 = 0.8_r8
!PK -               ELSEIF (hgt(i,1,j).GE.2500.0_r8) THEN
!PK -                  h_01 = 0.01_r8
!PK -               ELSE
!PK -                  h_01 = 0.8_r8*COS(hgt(i,1,j)*0.001_r8 - 1.0_r8)
!PK -               ENDIF
!PK -               niCCN3 = -1.0_r8*log(naCCN1/naCCN0)/h_01
!PK -               nwfa(i,1,j) = naCCN1+naCCN0*EXP(-((hgt(i,2,j)-hgt(i,1,j))/1000.0_r8)*niCCN3)
!PK -               DO k = 2, kMAx
!PK -                  nwfa(i,k,j) = naCCN1+naCCN0*EXP(-((hgt(i,k,j)-hgt(i,1,j))/1000.0_r8)*niCCN3)
!PK -               ENDDO
!PK -            ENDDO
!PK -         ENDDO
!PK -      ELSE
!PK -         has_CCN    = .TRUE.
!PK -         WRITE(mp_debug,*) ' Apparently initial CCN aerosols are present.'
!PK -         CALL wrf_debug(0, mp_debug)
!PK -         WRITE(mp_debug,*) '   column sum at point (i,j) = ', 1,1, SUM(nwfa(1,:,1))
!PK -         CALL wrf_debug(0, mp_debug)
!PK -      ENDIF
!PK -
!PK -
!PK -      !#if ( defined( DM_PARALLEL ) && ( ! defined( STUBMPI ) ) )
!PK -      !      max_test = wrf_dm_max_real ( MAXVAL(nifa(1:ibMax-1,:,1:jbMax-1)) )
!PK -      !#else
!PK -      max_test = MAXVAL ( nifa(1:ibMax-1,:,1:jbMax-1) )
!PK -      !#endif
!PK -
!PK -      IF (max_test .LT. eps) THEN
!PK -         WRITE(mp_debug,*) ' Apparently there are no initial IN aerosols.'
!PK -         CALL wrf_debug(0, mp_debug)
!PK -         WRITE(mp_debug,*) '   checked column at point (i,j) = ', 1,1
!PK -         CALL wrf_debug(0, mp_debug)
!PK -         DO j = 1, jbMax
!PK -            DO i =1,  ibMax
!PK -               IF (hgt(i,1,j).LE.1000.0_r8) THEN
!PK -                  h_01 = 0.8_r8
!PK -               ELSEIF (hgt(i,1,j).GE.2500.0_r8) THEN
!PK -                  h_01 = 0.01_r8
!PK -               ELSE
!PK -                  h_01 = 0.8_r8*COS(hgt(i,1,j)*0.001_r8 - 1.0_r8)
!PK -               ENDIF
!PK -               niIN3 = -1.0_r8*log(naIN1/naIN0)/h_01
!PK -               nifa(i,1,j) = naIN1+naIN0*EXP(-((hgt(i,2,j)-hgt(i,1,j))/1000.0_r8)*niIN3)
!PK -               DO k = 2, kMAx
!PK -                  nifa(i,k,j) = naIN1+naIN0*EXP(-((hgt(i,k,j)-hgt(i,1,j))/1000.0_r8)*niIN3)
!PK -               ENDDO
!PK -            ENDDO
!PK -         ENDDO
!PK -      ELSE
!PK -         has_IN     = .TRUE.
!PK -         WRITE(mp_debug,*) ' Apparently initial IN aerosols are present.'
!PK -         CALL wrf_debug(0, mp_debug)
!PK -         WRITE(mp_debug,*) '   column sum at point (i,j) = ', 1,1, SUM(nifa(1,:,1))
!PK -         CALL wrf_debug(0, mp_debug)
!PK -      ENDIF
!PK -
!PK -      !..Capture initial state lowest level CCN aerosol data in 2D array.
!PK -
!PK -      !     do j = 1, jbMax
!PK -      !     do i = 1, ibMax
!PK -      !        nwfa2d(i,j) = nwfa(i,1,j)
!PK -      !     enddo
!PK -      !     enddo
!PK -
!PK -      !..Scale the lowest level aerosol data into an emissions rate.  This is
!PK -      !.. very far from ideal, but need higher emissions where larger amount
!PK -      !.. of existing and lesser emissions where not already lots of aerosols
!PK -      !.. for first-order simplistic approach.  Later, proper connection to
!PK -      !.. emission inventory would be better, but, for now, scale like this:
!PK -      !.. where: Nwfa=50 per cc, emit 0.875E4 aerosols per kg per second
!PK -      !..        Nwfa=500 per cc, emit 0.875E5 aerosols per kg per second
!PK -      !..        Nwfa=5000 per cc, emit 0.875E6 aerosols per kg per second
!PK -      !.. for a grid with 20km spacing and scale accordingly for other spacings.
!PK -
!PK -      IF (is_start) THEN
!PK -         IF (SQRT(DX*DY)/20000.0_r8 .GE. 1.0_r8) THEN
!PK -            h_01 = 0.875_r8
!PK -         ELSE
!PK -            h_01 = (0.875_r8 + 0.125_r8*((20000.0_r8-SQRT(DX*DY))/16000.0_r8)) * SQRT(DX*DY)/20000.0_r8
!PK -         ENDIF
!PK -         WRITE(mp_debug,*) '   aerosol surface flux emission scale factor is: ', h_01
!PK -         CALL wrf_debug(0, mp_debug)
!PK -         DO j = 1, jbMax
!PK -            DO i = 1, ibMax
!PK -               nwfa2d(i,j) = 10.0_r8**(LOG10(nwfa(i,1,j)*1.E-6_r8)-3.69897_r8)
!PK -               nwfa2d(i,j) = nwfa2d(i,j)*h_01 * 1.E6_r8
!PK -            ENDDO
!PK -         ENDDO
!PK -         !     else
!PK -         !        write(mp_debug,*) '   sample (lower-left) aerosol surface flux emission rate: ', nwfa2d(1,1)
!PK -         !        CALL wrf_debug(0, mp_debug)
!PK -      ENDIF
!PK -
!PK -   ENDIF


    !..Allocate space for lookup tables (J. Michalakes 2009Jun08).

    IF (.NOT. ALLOCATED(tcg_racg) ) THEN
       ALLOCATE(tcg_racg(ntb_g1,ntb_g,ntb_r1,ntb_r))
       micro_init = .TRUE.
    ENDIF

    IF (.NOT. ALLOCATED(tmr_racg)) ALLOCATE(tmr_racg(ntb_g1,ntb_g,ntb_r1,ntb_r))
    IF (.NOT. ALLOCATED(tcr_gacr)) ALLOCATE(tcr_gacr(ntb_g1,ntb_g,ntb_r1,ntb_r))
    IF (.NOT. ALLOCATED(tmg_gacr)) ALLOCATE(tmg_gacr(ntb_g1,ntb_g,ntb_r1,ntb_r))
    IF (.NOT. ALLOCATED(tnr_racg)) ALLOCATE(tnr_racg(ntb_g1,ntb_g,ntb_r1,ntb_r))
    IF (.NOT. ALLOCATED(tnr_gacr)) ALLOCATE(tnr_gacr(ntb_g1,ntb_g,ntb_r1,ntb_r))

    IF (.NOT. ALLOCATED(tcs_racs1)) ALLOCATE(tcs_racs1(ntb_s,ntb_t,ntb_r1,ntb_r))
    IF (.NOT. ALLOCATED(tmr_racs1)) ALLOCATE(tmr_racs1(ntb_s,ntb_t,ntb_r1,ntb_r))
    IF (.NOT. ALLOCATED(tcs_racs2)) ALLOCATE(tcs_racs2(ntb_s,ntb_t,ntb_r1,ntb_r))
    IF (.NOT. ALLOCATED(tmr_racs2)) ALLOCATE(tmr_racs2(ntb_s,ntb_t,ntb_r1,ntb_r))
    IF (.NOT. ALLOCATED(tcr_sacr1)) ALLOCATE(tcr_sacr1(ntb_s,ntb_t,ntb_r1,ntb_r))
    IF (.NOT. ALLOCATED(tms_sacr1)) ALLOCATE(tms_sacr1(ntb_s,ntb_t,ntb_r1,ntb_r))
    IF (.NOT. ALLOCATED(tcr_sacr2)) ALLOCATE(tcr_sacr2(ntb_s,ntb_t,ntb_r1,ntb_r))
    IF (.NOT. ALLOCATED(tms_sacr2)) ALLOCATE(tms_sacr2(ntb_s,ntb_t,ntb_r1,ntb_r))
    IF (.NOT. ALLOCATED(tnr_racs1)) ALLOCATE(tnr_racs1(ntb_s,ntb_t,ntb_r1,ntb_r))
    IF (.NOT. ALLOCATED(tnr_racs2)) ALLOCATE(tnr_racs2(ntb_s,ntb_t,ntb_r1,ntb_r))
    IF (.NOT. ALLOCATED(tnr_sacr1)) ALLOCATE(tnr_sacr1(ntb_s,ntb_t,ntb_r1,ntb_r))
    IF (.NOT. ALLOCATED(tnr_sacr2)) ALLOCATE(tnr_sacr2(ntb_s,ntb_t,ntb_r1,ntb_r))

    IF (.NOT. ALLOCATED(tpi_qcfz)) ALLOCATE(tpi_qcfz(ntb_c,nbc,45,ntb_IN))
    IF (.NOT. ALLOCATED(tni_qcfz)) ALLOCATE(tni_qcfz(ntb_c,nbc,45,ntb_IN))

    IF (.NOT. ALLOCATED(tpi_qrfz)) ALLOCATE(tpi_qrfz(ntb_r,ntb_r1,45,ntb_IN))
    IF (.NOT. ALLOCATED(tpg_qrfz)) ALLOCATE(tpg_qrfz(ntb_r,ntb_r1,45,ntb_IN))
    IF (.NOT. ALLOCATED(tni_qrfz)) ALLOCATE(tni_qrfz(ntb_r,ntb_r1,45,ntb_IN))
    IF (.NOT. ALLOCATED(tnr_qrfz)) ALLOCATE(tnr_qrfz(ntb_r,ntb_r1,45,ntb_IN))

    IF (.NOT. ALLOCATED(tps_iaus)) ALLOCATE(tps_iaus(ntb_i,ntb_i1))
    IF (.NOT. ALLOCATED(tni_iaus)) ALLOCATE(tni_iaus(ntb_i,ntb_i1))
    IF (.NOT. ALLOCATED(tpi_ide)) ALLOCATE(tpi_ide(ntb_i,ntb_i1))

    IF (.NOT. ALLOCATED(t_Efrw)) ALLOCATE(t_Efrw(nbr,nbc))
    IF (.NOT. ALLOCATED(t_Efsw)) ALLOCATE(t_Efsw(nbs,nbc))

    IF (.NOT. ALLOCATED(tnr_rev)) ALLOCATE(tnr_rev(nbr, ntb_r1, ntb_r))
    IF (.NOT. ALLOCATED(tpc_wev)) ALLOCATE(tpc_wev(nbc,ntb_c,nbc))
    IF (.NOT. ALLOCATED(tnc_wev)) ALLOCATE(tnc_wev(nbc,ntb_c,nbc))

    IF (.NOT. ALLOCATED(tnccn_act))                                   &
         ALLOCATE(tnccn_act(ntb_arc,ntb_arw,ntb_art,ntb_arr,ntb_ark))
    CALL MsgOne('******Init_Micro_thompson*******','ALLOCATE VAR')
    IF (micro_init) THEN

       !..From Martin et al. (1994), assign gamma shape parameter mu for cloud
       !.. drops according to general dispersion characteristics (disp=~0.25
       !.. for Maritime and 0.45 for Continental).
       !.. disp=SQRT((mu+2)/(mu+1) - 1) so mu varies from 15 for Maritime
       !.. to 2 for really dirty air.  This not used in 2-moment cloud water
       !.. scheme and nu_c used instead and varies from 2 to 15 (integer-only).
       mu_c = MIN(15.0_r8, (1000.E6_r8/Nt_c + 2.0_r8))

       !..Schmidt number to one-third used numerous times.
       Sc3 = Sc**(1.0_r8/3.0_r8)

       !..Compute min ice diam from mass, min snow/graupel mass from diam.
       D0i = (xm0i/am_i)**(1.0_r8/bm_i)
       xm0s = am_s * D0s**bm_s
       xm0g = am_g * D0g**bm_g

       !..These constants various exponents and gamma() assoc with cloud,
       !.. rain, snow, and graupel.
       DO n = 1, 15
          cce(1,n) = n + 1.0_r8
          cce(2,n) = bm_r + n + 1.0_r8
          cce(3,n) = bm_r + n + 4.0_r8
          cce(4,n) = n + bv_c + 1.0_r8
          cce(5,n) = bm_r + n + bv_c + 1.0_r8
          ccg(1,n) = WGAMMA(cce(1,n))
          ccg(2,n) = WGAMMA(cce(2,n))
          ccg(3,n) = WGAMMA(cce(3,n))
          ccg(4,n) = WGAMMA(cce(4,n))
          ccg(5,n) = WGAMMA(cce(5,n))
          ocg1(n) = 1.0_r8/ccg(1,n)
          ocg2(n) = 1.0_r8/ccg(2,n)
       ENDDO

       cie(1) = mu_i + 1.0_r8
       cie(2) = bm_i + mu_i + 1.0_r8
       cie(3) = bm_i + mu_i + bv_i + 1.0_r8
       cie(4) = mu_i + bv_i + 1.0_r8
       cie(5) = mu_i + 2.0_r8
       cie(6) = bm_i*0.5_r8 + mu_i + bv_i + 1.0_r8
       cie(7) = bm_i*0.5_r8 + mu_i + 1.0_r8
       cig(1) = WGAMMA(cie(1))
       cig(2) = WGAMMA(cie(2))
       cig(3) = WGAMMA(cie(3))
       cig(4) = WGAMMA(cie(4))
       cig(5) = WGAMMA(cie(5))
       cig(6) = WGAMMA(cie(6))
       cig(7) = WGAMMA(cie(7))
       oig1 = 1.0_r8/cig(1)
       oig2 = 1.0_r8/cig(2)
       obmi = 1.0_r8/bm_i

       cre(1) = bm_r + 1.0_r8
       cre(2) = mu_r + 1.0_r8
       cre(3) = bm_r + mu_r + 1.0_r8
       cre(4) = bm_r*2.0_r8 + mu_r + 1.0_r8
       cre(5) = mu_r + bv_r + 1.0_r8
       cre(6) = bm_r + mu_r + bv_r + 1.0_r8
       cre(7) = bm_r*0.5_r8 + mu_r + bv_r + 1.0_r8
       cre(8) = bm_r + mu_r + bv_r + 3.0_r8
       cre(9) = mu_r + bv_r + 3.0_r8
       cre(10) = mu_r + 2.0_r8
       cre(11) = 0.5_r8*(bv_r + 5.0_r8 + 2.0_r8*mu_r)
       cre(12) = bm_r*0.5_r8 + mu_r + 1.0_r8
       cre(13) = bm_r*2.0_r8 + mu_r + bv_r + 1.0_r8
       DO n = 1, 13
          crg(n) = WGAMMA(cre(n))
       ENDDO
       obmr = 1.0_r8/bm_r
       ore1 = 1.0_r8/cre(1)
       org1 = 1.0_r8/crg(1)
       org2 = 1.0_r8/crg(2)
       org3 = 1.0_r8/crg(3)

       cse(1) = bm_s + 1.0_r8
       cse(2) = bm_s + 2.0_r8
       cse(3) = bm_s*2.0_r8
       cse(4) = bm_s + bv_s + 1.0_r8
       cse(5) = bm_s*2.0_r8 + bv_s + 1.0_r8
       cse(6) = bm_s*2.0_r8 + 1.0_r8
       cse(7) = bm_s + mu_s + 1.0_r8
       cse(8) = bm_s + mu_s + 2.0_r8
       cse(9) = bm_s + mu_s + 3.0_r8
       cse(10) = bm_s + mu_s + bv_s + 1.0_r8
       cse(11) = bm_s*2.0_r8 + mu_s + bv_s + 1.0_r8
       cse(12) = bm_s*2.0_r8 + mu_s + 1.0_r8
       cse(13) = bv_s + 2.0_r8
       cse(14) = bm_s + bv_s
       cse(15) = mu_s + 1.0_r8
       cse(16) = 1.0_r8 + (1.0_r8 + bv_s)/2.0_r8
       cse(17) = cse(16) + mu_s + 1.0_r8
       cse(18) = bv_s + mu_s + 3.0_r8
       DO n = 1, 18
          csg(n) = WGAMMA(cse(n))
       ENDDO
       oams = 1.0_r8/am_s
       obms = 1.0_r8/bm_s
       ocms = oams**obms

       cge(1) = bm_g + 1.0_r8
       cge(2) = mu_g + 1.0_r8
       cge(3) = bm_g + mu_g + 1.0_r8
       cge(4) = bm_g*2.0_r8 + mu_g + 1.0_r8
       cge(5) = bm_g*2.0_r8 + mu_g + bv_g + 1.0_r8
       cge(6) = bm_g + mu_g + bv_g + 1.0_r8
       cge(7) = bm_g + mu_g + bv_g + 2.0_r8
       cge(8) = bm_g + mu_g + bv_g + 3.0_r8
       cge(9) = mu_g + bv_g + 3.0_r8
       cge(10) = mu_g + 2.0_r8
       cge(11) = 0.5_r8*(bv_g + 5.0_r8 + 2.0_r8*mu_g)
       cge(12) = 0.5_r8*(bv_g + 5.0_r8) + mu_g
       DO n = 1, 12
          cgg(n) = WGAMMA(cge(n))
       ENDDO
       oamg = 1.0_r8/am_g
       obmg = 1.0_r8/bm_g
       ocmg = oamg**obmg
       oge1 = 1.0_r8/cge(1)
       ogg1 = 1.0_r8/cgg(1)
       ogg2 = 1.0_r8/cgg(2)
       ogg3 = 1.0_r8/cgg(3)

       !+---+-----------------------------------------------------------------+
       !..Simplify various rate eqns the best we can now.
       !+---+-----------------------------------------------------------------+

       !..Rain collecting cloud water and cloud ice
       t1_qr_qc = PI*0.25_r8*av_r * crg(9)
       t1_qr_qi = PI*0.25_r8*av_r * crg(9)
       t2_qr_qi = PI*0.25_r8*am_r*av_r * crg(8)

       !..Graupel collecting cloud water
       t1_qg_qc = PI*0.25_r8*av_g * cgg(9)

       !..Snow collecting cloud water
       t1_qs_qc = PI*0.25_r8*av_s

       !..Snow collecting cloud ice
       t1_qs_qi = PI*0.25_r8*av_s

       !..Evaporation of rain; ignore depositional growth of rain.
       t1_qr_ev = 0.78_r8 * crg(10)
       t2_qr_ev = 0.308_r8*Sc3*SQRT(av_r) * crg(11)

       !..Sublimation/depositional growth of snow
       t1_qs_sd = 0.86_r8
       t2_qs_sd = 0.28_r8*Sc3*SQRT(av_s)

       !..Melting of snow
       t1_qs_me = PI*4.0_r8*C_sqrd*olfus * 0.86_r8
       t2_qs_me = PI*4.0_r8*C_sqrd*olfus * 0.28_r8*Sc3*SQRT(av_s)

       !..Sublimation/depositional growth of graupel
       t1_qg_sd = 0.86_r8 * cgg(10)
       t2_qg_sd = 0.28_r8*Sc3*SQRT(av_g) * cgg(11)

       !..Melting of graupel
       t1_qg_me = PI*4.0_r8*C_cube*olfus * 0.86_r8 * cgg(10)
       t2_qg_me = PI*4.0_r8*C_cube*olfus * 0.28_r8*Sc3*SQRT(av_g) * cgg(11)

       !..Constants for helping find lookup table indexes.
       nic2 = NINT(log10(r_c(1)))
       nii2 = NINT(log10(r_i(1)))
       nii3 = NINT(log10(Nt_i(1)))
       nir2 = NINT(log10(r_r(1)))
       nir3 = NINT(log10(N0r_exp(1)))
       nis2 = NINT(log10(r_s(1)))
       nig2 = NINT(log10(r_g(1)))
       nig3 = NINT(log10(N0g_exp(1)))
       niIN2 = NINT(log10(Nt_IN(1)))

       !..Create bins of cloud water (from min diameter up to 100 microns).
       Dc(1) = D0c*1.0e0_r8
       dtc(1) = D0c*1.0e0_r8
       DO n = 2, nbc
          Dc(n) = Dc(n-1) + 1.0e-6_r8
          dtc(n) = (Dc(n) - Dc(n-1))
       ENDDO

       !..Create bins of cloud ice (from min diameter up to 5x min snow size).
       xDx(1) = D0i*1.0e0_r8
       xDx(nbi+1) = 5.0e0_r8*D0s
       DO n = 2, nbi
          xDx(n) = DEXP(DFLOAT(n-1)/DFLOAT(nbi) &
               *DLOG(xDx(nbi+1)/xDx(1)) +DLOG(xDx(1)))
       ENDDO
       DO n = 1, nbi
          Di(n) = DSQRT(xDx(n)*xDx(n+1))
          dti(n) = xDx(n+1) - xDx(n)
       ENDDO

       !..Create bins of rain (from min diameter up to 5 mm).
       xDx(1) = D0r*1.0e0_r8
       xDx(nbr+1) = 0.005e0_r8
       DO n = 2, nbr
          xDx(n) = DEXP(DFLOAT(n-1)/DFLOAT(nbr) &
               *DLOG(xDx(nbr+1)/xDx(1)) +DLOG(xDx(1)))
       ENDDO
       DO n = 1, nbr
          Dr(n) = DSQRT(xDx(n)*xDx(n+1))
          dtr(n) = xDx(n+1) - xDx(n)
       ENDDO

       !..Create bins of snow (from min diameter up to 2 cm).
       xDx(1) = D0s*1.0e0_r8
       xDx(nbs+1) = 0.02e0_r8
       DO n = 2, nbs
          xDx(n) = DEXP(DFLOAT(n-1)/DFLOAT(nbs) &
               *DLOG(xDx(nbs+1)/xDx(1)) +DLOG(xDx(1)))
       ENDDO
       DO n = 1, nbs
          Ds(n) = DSQRT(xDx(n)*xDx(n+1))
          dts(n) = xDx(n+1) - xDx(n)
       ENDDO

       !..Create bins of graupel (from min diameter up to 5 cm).
       xDx(1) = D0g*1.0e0_r8
       xDx(nbg+1) = 0.05e0_r8
       DO n = 2, nbg
          xDx(n) = DEXP(DFLOAT(n-1)/DFLOAT(nbg) &
               *DLOG(xDx(nbg+1)/xDx(1)) +DLOG(xDx(1)))
       ENDDO
       DO n = 1, nbg
          Dg(n) = DSQRT(xDx(n)*xDx(n+1))
          dtg(n) = xDx(n+1) - xDx(n)
       ENDDO

       !..Create bins of cloud droplet number concentration (1 to 3000 per cc).
       xDx(1) = 1.0e0_r8
       xDx(nbc+1) = 3000.0e0_r8
       DO n = 2, nbc
          xDx(n) = DEXP(DFLOAT(n-1)/DFLOAT(nbc)                          &
               *DLOG(xDx(nbc+1)/xDx(1)) +DLOG(xDx(1)))
       ENDDO
       DO n = 1, nbc
          t_Nc(n) = DSQRT(xDx(n)*xDx(n+1)) * 1.e6_r8
       ENDDO
       nic1 = INT( LOG(t_Nc(nbc)/t_Nc(1)))

       !+---+-----------------------------------------------------------------+
       !..Create lookup tables for most costly calculations.
       !+---+-----------------------------------------------------------------+

       DO m = 1, ntb_r
          DO k = 1, ntb_r1
             DO j = 1, ntb_g
                DO i = 1, ntb_g1
                   tcg_racg(i,j,k,m) = 0.0e0_r8
                   tmr_racg(i,j,k,m) = 0.0e0_r8
                   tcr_gacr(i,j,k,m) = 0.0e0_r8
                   tmg_gacr(i,j,k,m) = 0.0e0_r8
                   tnr_racg(i,j,k,m) = 0.0e0_r8
                   tnr_gacr(i,j,k,m) = 0.0e0_r8
                ENDDO
             ENDDO
          ENDDO
       ENDDO

       DO m = 1, ntb_r
          DO k = 1, ntb_r1
             DO j = 1, ntb_t
                DO i = 1, ntb_s
                   tcs_racs1(i,j,k,m) = 0.0e0_r8
                   tmr_racs1(i,j,k,m) = 0.0e0_r8
                   tcs_racs2(i,j,k,m) = 0.0e0_r8
                   tmr_racs2(i,j,k,m) = 0.0e0_r8
                   tcr_sacr1(i,j,k,m) = 0.0e0_r8
                   tms_sacr1(i,j,k,m) = 0.0e0_r8
                   tcr_sacr2(i,j,k,m) = 0.0e0_r8
                   tms_sacr2(i,j,k,m) = 0.0e0_r8
                   tnr_racs1(i,j,k,m) = 0.0e0_r8
                   tnr_racs2(i,j,k,m) = 0.0e0_r8
                   tnr_sacr1(i,j,k,m) = 0.0e0_r8
                   tnr_sacr2(i,j,k,m) = 0.0e0_r8
                ENDDO
             ENDDO
          ENDDO
       ENDDO

       DO m = 1, ntb_IN
          DO k = 1, 45
             DO j = 1, ntb_r1
                DO i = 1, ntb_r
                   tpi_qrfz(i,j,k,m) = 0.0e0_r8
                   tni_qrfz(i,j,k,m) = 0.0e0_r8
                   tpg_qrfz(i,j,k,m) = 0.0e0_r8
                   tnr_qrfz(i,j,k,m) = 0.0e0_r8
                ENDDO
             ENDDO
             DO j = 1, nbc
                DO i = 1, ntb_c
                   tpi_qcfz(i,j,k,m) = 0.0e0_r8
                   tni_qcfz(i,j,k,m) = 0.0e0_r8
                ENDDO
             ENDDO
          ENDDO
       ENDDO

       DO j = 1, ntb_i1
          DO i = 1, ntb_i
             tps_iaus(i,j) = 0.0e0_r8
             tni_iaus(i,j) = 0.0e0_r8
             tpi_ide(i,j) = 0.0e0_r8
          ENDDO
       ENDDO

       DO j = 1, nbc
          DO i = 1, nbr
             t_Efrw(i,j) = 0.0_r8
          ENDDO
          DO i = 1, nbs
             t_Efsw(i,j) = 0.0_r8
          ENDDO
       ENDDO

       DO k = 1, ntb_r
          DO j = 1, ntb_r1
             DO i = 1, nbr
                tnr_rev(i,j,k) = 0.0e0_r8
             ENDDO
          ENDDO
       ENDDO

       DO k = 1, nbc
          DO j = 1, ntb_c
             DO i = 1, nbc
                tpc_wev(i,j,k) = 0.0e0_r8
                tnc_wev(i,j,k) = 0.0e0_r8
             ENDDO
          ENDDO
       ENDDO

       DO m = 1, ntb_ark
          DO l = 1, ntb_arr
             DO k = 1, ntb_art
                DO j = 1, ntb_arw
                   DO i = 1, ntb_arc
                      tnccn_act(i,j,k,l,m) = 1.0_r8
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       CALL wrf_debug('******Init_Micro_thompson*******', 'CREATING MICROPHYSICS LOOKUP TABLES ... ')
       WRITE (wrf_err_message, '(a, f5.2, a, f5.2, a, f5.2, a, f5.2)') &
            ' using: mu_c=',mu_c,' mu_i=',mu_i,' mu_r=',mu_r,' mu_g=',mu_g
       CALL wrf_debug('******Init_Micro_thompson*******', wrf_err_message)

       !..Read a static file containing CCN activation of aerosols. The
       !.. data were created from a parcel model by Feingold & Heymsfield with
       !.. further changes by Eidhammer and Kriedenweis.
       IF (is_aerosol_aware) THEN
          CALL wrf_debug('******Init_Micro_thompson*******', '  calling table_ccnAct routine')
          CALL table_ccnAct(path_in)
       ENDIF

       !..Collision efficiency between rain/snow and cloud water.
       CALL wrf_debug('******Init_Micro_thompson*******', '  creating qc collision eff tables')
       CALL table_Efrw()
       CALL table_Efsw()

       !..Drop evaporation.
       CALL wrf_debug('******Init_Micro_thompson*******', '  creating rain evap table')
       CALL table_dropEvap()

       !..Initialize various constants for computing radar reflectivity.
       xam_r = am_r
       xbm_r = bm_r
       xmu_r = mu_r
       xam_s = am_s
       xbm_s = bm_s
       xmu_s = mu_s
       xam_g = am_g
       xbm_g = bm_g
       xmu_g = mu_g

       CALL radar_init()


       IF (.NOT. iiwarm) THEN

          !..Rain collecting graupel & graupel collecting rain.
          CALL wrf_debug('******Init_Micro_thompson*******', '  creating rain collecting graupel table')
          CALL qr_acr_qg(path_in)

          !..Rain collecting snow & snow collecting rain.
          CALL wrf_debug('******Init_Micro_thompson*******', '  creating rain collecting snow table')
          CALL qr_acr_qs(path_in)

          !..Cloud water and rain freezing (Bigg, 1953).
          CALL wrf_debug('******Init_Micro_thompson*******', '  creating freezing of water drops table')
          CALL freezeH2O(path_in)

          !..Conversion of some ice mass into snow category.
          CALL wrf_debug('******Init_Micro_thompson*******', '  creating ice converting to snow table')
          CALL qi_aut_qs()

       ENDIF

       CALL wrf_debug('******Init_Micro_thompson*******', ' ... DONE microphysical lookup tables')

    ENDIF

  END SUBROUTINE Init_Micro_thompson

  


  SUBROUTINE RunMicro_thompson( &
       nCols       , &!INTEGER      , INTENT(IN   ) :: nCols
       kMax        , &!INTEGER      , INTENT(IN   ) :: kMax 
       prsi        , &
       prsl        , &
       tc          , &!REAL(KIND=r8), INTENT(INOUT) :: Tc (1:nCols, 1:kMax)
       QV          , &!REAL(KIND=r8), INTENT(INOUT) :: qv (1:nCols, 1:kMax)
       QC          , &!REAL(KIND=r8), INTENT(INOUT) :: qc (1:nCols, 1:kMax)
       QR          , &!REAL(KIND=r8), INTENT(INOUT) :: qr (1:nCols, 1:kMax)
       QI          , &!REAL(KIND=r8), INTENT(INOUT) :: qi (1:nCols, 1:kMax)
       QS          , &!REAL(KIND=r8), INTENT(INOUT) :: qs (1:nCols, 1:kMax)
       QG          , &!REAL(KIND=r8), INTENT(INOUT) :: qg (1:nCols, 1:kMax)
       NI          , &!REAL(KIND=r8), INTENT(INOUT) :: ni (1:nCols, 1:kMax)
       NR          , &!REAL(KIND=r8), INTENT(INOUT) :: nr (1:nCols, 1:kMax)
       nifa        , &!REAL(KIND=r8), INTENT(INOUT) :: ns (1:nCols, 1:kMax)
       nwfa        , &!REAL(KIND=r8), INTENT(INOUT) :: NG (1:nCols, 1:kMax)   
       nc          , &!REAL(KIND=r8), INTENT(INOUT) :: NC (1:nCols, 1:kMax)   
       dTcdt       , &!
       dqvdt       , &!
       dqcdt       , &!
       dqrdt       , &!
       dqidt       , &!
       dqsdt       , &!
       dqgdt       , &!
       dnidt       , &!
       dnrdt       , &!
       dnifadt     , &!
       dnwfadt     , &!
       dncdt       , &!
       TKE         , &!REAL(KIND=r8), INTENT(IN   ) :: TKE (1:nCols, 1:kMax)   
       KZH         , &!REAL(KIND=r8), INTENT(IN   ) :: KZH (1:nCols, 1:kMax)   
       DT_IN       , &!REAL(KIND=r8), INTENT(IN   ) :: dt_in
       omega       , &!REAL(KIND=r8), INTENT(IN   ) :: omega  ! omega (Pa/s)
       EFFCS       , &!REAL(KIND=r8), INTENT(OUT  ) :: EFFCS (1:nCols, 1:kMax)   ! EFFCS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
       EFFIS       , &!REAL(KIND=r8), INTENT(OUT  ) :: EFFIS (1:nCols, 1:kMax)   ! EFFIS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
       LSRAIN      , &!REAL(KIND=r8), INTENT(OUT) :: LSRAIN(1:nCols)
       LSSNOW        )!REAL(KIND=r8), INTENT(OUT) :: LSSNOW(1:nCols)

    IMPLICIT NONE
    INTEGER      , INTENT(IN   ) :: nCols
    INTEGER      , INTENT(IN   ) :: kMax 
    REAL(KIND=r8), INTENT(IN   ) :: prsi       (1:nCols,1:kMax+1)  
    REAL(KIND=r8), INTENT(IN   ) :: prsl       (1:nCols,1:kMax)    
    ! Temporary changed from INOUT to IN
    REAL(KIND=r8), INTENT(INOUT) :: Tc   (1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qv   (1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qc   (1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qr   (1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qi   (1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qs   (1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qg   (1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: ni   (1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: nr   (1:nCols, 1:kMax)    
    REAL(KIND=r8), INTENT(INOUT) :: nifa (1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: nwfa (1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: nc   (1:nCols, 1:kMax)

    REAL(KIND=r8), INTENT(INOUT) :: dTcdt   (1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: dqvdt   (1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: dqcdt   (1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: dqrdt   (1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: dqidt   (1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: dqsdt   (1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: dqgdt   (1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: dnidt   (1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: dnrdt   (1:nCols, 1:kMax)    
    REAL(KIND=r8), INTENT(INOUT) :: dnifadt (1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: dnwfadt (1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: dncdt   (1:nCols, 1:kMax)


    REAL(KIND=r8), INTENT(IN   ) :: tke(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(IN   ) :: kzh(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(IN   ) :: dt_in
    REAL(KIND=r8), INTENT(IN   ) :: omega (1:nCols, 1:kMax) ! omega (Pa/s)
    REAL(KIND=r8), INTENT(OUT  ) :: LSRAIN(1:nCols)
    REAL(KIND=r8), INTENT(OUT  ) :: LSSNOW(1:nCols)
    REAL(KIND=r8), INTENT(OUT  ) :: EFFCS (1:nCols, 1:kMax)   ! EFFCS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
    ! note: effis not currently passed out of microphysics (no coupling of ice eff rad with radiation)
    REAL(KIND=r8), INTENT(OUT  ) :: EFFIS (1:nCols, 1:kMax)   ! EFFIS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)


!nc          (1:nCols, 1:kMax)  , & ! REAL(KIND=r8)  , OPTIONAL, INTENT(INOUT) :: nc  (1:nCols, 1:kMax)
!nwfa          (1:nCols, 1:kMax)  , & ! REAL(KIND=r8)  , OPTIONAL, INTENT(INOUT) :: nwfa(1:nCols, 1:kMax)
!nifa          (1:nCols, 1:kMax)  , & ! REAL(KIND=r8)  , OPTIONAL, INTENT(INOUT) :: nifa(1:nCols, 1:kMax)
    REAL(KIND=r8) :: nwfa2d(1:nCols)   !, OPTIONAL, INTENT(IN   )
    REAL(KIND=r8) :: SR        (1:nCols)
    REAL(KIND=r8) :: refl_10cm (1:nCols, 1:kMax)
    REAL(KIND=r8) :: re_snow   (1:nCols, 1:kMax)
    ! add cumulus tendencies

    REAL(KIND=r8) :: qrcuten (1:nCols, 1:kMax)
    REAL(KIND=r8) :: qscuten (1:nCols, 1:kMax)
    REAL(KIND=r8) :: qicuten (1:nCols, 1:kMax)
    REAL(KIND=r8) :: mu      (1:nCols)

    LOGICAL       :: diagflag,is_start
    INTEGER       :: do_radar_ref  ! GT added for reflectivity calcs
!    LOGICAL       :: F_QNDROP      ! wrf-chem
    REAL(KIND=r8) :: qndrop(1:nCols, 1:kMax) ! hm added, wrf-chem 

    REAL(KIND=r8) :: flip_pint  (nCols,kMax+1)   ! Interface pressures  
    REAL(KIND=r8) :: flip_pmid  (nCols,kMax)! Midpoint pressures 
    REAL(KIND=r8) :: flip_t     (nCols,kMax)! temperature
    REAL(KIND=r8) :: flip_q     (nCols,kMax)! specific humidity
    REAL(KIND=r8) :: flip_pdel  (nCols,kMax)! layer thickness
    REAL(KIND=r8) :: flip_rpdel (nCols,kMax)! inverse of layer thickness
    REAL(KIND=r8) :: flip_lnpmid(nCols,kMax)! Log Midpoint pressures    
    REAL(KIND=r8) :: flip_lnpint(nCols,kMax+1)   ! Log interface pressures
    REAL(KIND=r8) :: flip_zi    (nCols,kMax+1)! Height above surface at interfaces 
    REAL(KIND=r8) :: flip_zm    (nCols,kMax)  ! Geopotential height at mid level

    REAL(KIND=r8) :: zi         (nCols,kMax+1)     ! Height above surface at interfaces
    REAL(KIND=r8) :: zm         (nCols,kMax)        ! Geopotential height at mid level
    REAL(KIND=r8) :: p          (nCols,kMax) ! pressure at all points, on u,v levels (Pa).
    REAL(KIND=r8) :: dz         (nCols,kMax)
    REAL(KIND=r8) :: RAINNC (1:nCols)
    REAL(KIND=r8) :: RAINNCV(1:nCols)
    REAL(KIND=r8) :: SNOWNC    (1:nCols)
    REAL(KIND=r8) :: SNOWNCV   (1:nCols)
    REAL(KIND=r8) :: GRAUPELNC (1:nCols)
    REAL(KIND=r8) :: GRAUPELNCV(1:nCols)
    REAL(KIND=r8) :: rho(nCols,kMax)     
    REAL(KIND=r8) :: w  (1:nCols, 1:kMax) !, tke, nctend, nnColsnd,kzh
    REAL(KIND=r8) :: rainprod  (1:nCols, 1:kMax)
    REAL(KIND=r8) :: evapprod  (1:nCols, 1:kMax)


    REAL(KIND=r8) :: Tc_mic   (1:nCols, 1:kMax)
    REAL(KIND=r8) :: qv_mic   (1:nCols, 1:kMax)
    REAL(KIND=r8) :: qc_mic   (1:nCols, 1:kMax)
    REAL(KIND=r8) :: qr_mic   (1:nCols, 1:kMax)
    REAL(KIND=r8) :: qi_mic   (1:nCols, 1:kMax)
    REAL(KIND=r8) :: qs_mic   (1:nCols, 1:kMax)
    REAL(KIND=r8) :: qg_mic   (1:nCols, 1:kMax)
    REAL(KIND=r8) :: ni_mic   (1:nCols, 1:kMax)
    REAL(KIND=r8) :: nr_mic   (1:nCols, 1:kMax)    
    REAL(KIND=r8) :: nifa_mic (1:nCols, 1:kMax)
    REAL(KIND=r8) :: nwfa_mic (1:nCols, 1:kMax)
    REAL(KIND=r8) :: nc_mic   (1:nCols, 1:kMax)

    INTEGER       :: has_reqc   =1                 ! , & ! INTEGER             , INTENT(IN   ) :: has_reqc
    INTEGER       :: has_reqi   =1                 ! , & ! INTEGER             , INTENT(IN   ) :: has_reqi
    INTEGER       :: has_reqs   =1                 ! , & ! INTEGER             , INTENT(IN   ) :: has_reqs
    CHARACTER(LEN=256) :: mp_debug
    REAL(KIND=r8) :: h_01,max_test,niCCN3,niIN3
    INTEGER       :: I,K,kflip
    nwfa2d=0.0_r8
    refl_10cm=0.0_r8
    RAINNC =0.0_r8
    rainprod =0.0_r8
    evapprod =0.0_r8
    RAINNCV=0.0_r8
    SNOWNC    =0.0_r8
    SNOWNCV   =0.0_r8
    GRAUPELNC =0.0_r8
    GRAUPELNCV=0.0_r8
    qrcuten=0.0_r8
    qscuten=0.0_r8
    qicuten=0.0_r8
    mu     =1.0_r8
!    F_QNDROP=.FALSE.
    qndrop=0.0_r8
    EFFCS =0.0_r8
    EFFIS =0.0_r8
    diagflag=.TRUE.
    is_start=.TRUE.
    do_radar_ref=1
    DO i=1,nCols
       !flip_pint       (i,kMax+1) = gps(i)*si(1) ! gps --> Pa
       flip_pint       (i,kMax+1) = prsi(i,1)
    END DO
    DO k=kMax,1,-1
       kflip=kMax+2-k
       DO i=1,nCols
         ! flip_pint    (i,k)      = MAX(si(kflip)*gps(i) ,1.0e-12_r8)
          flip_pint    (i,k)      = MAX(prsi(i,kflip)    ,1.0e-12_r8)
       END DO
    END DO
    DO k=1,kMax
       kflip=kMax+1-k
       DO i=1,nCols

          dTcdt   (i, k)=0.0_r8
          dqvdt   (i, k)=0.0_r8
          dqcdt   (i, k)=0.0_r8
          dqrdt   (i, k)=0.0_r8
          dqidt   (i, k)=0.0_r8
          dqsdt   (i, k)=0.0_r8
          dqgdt   (i, k)=0.0_r8
          dnidt   (i, k)=0.0_r8
          dnrdt   (i, k)=0.0_r8
          dnifadt (i, k)=0.0_r8
          dnwfadt (i, k)=0.0_r8
          dncdt   (i, k)=0.0_r8

          Tc_mic   (i, k) = Tc   (i, k)
          qv_mic   (i, k) = qv   (i, k)
          qc_mic   (i, k) = qc   (i, k)
          qr_mic   (i, k) = qr   (i, k)
          qi_mic   (i, k) = qi   (i, k)
          qs_mic   (i, k) = qs   (i, k)
          qg_mic   (i, k) = qg   (i, k)
          ni_mic   (i, k) = ni   (i, k)
          nr_mic   (i, k) = nr   (i, k)
          nifa_mic (i, k) = nifa (i, k)
          nwfa_mic (i, k) = nwfa (i, k)
          nc_mic   (i, k) = nc   (i, k)


          flip_t   (i,kflip) =  TC_mic (i,k)
          flip_q   (i,kflip) =  qv_mic (i,k)
          !flip_pmid(i,kflip) =  sl(  k)*gps (i)
          flip_pmid(i,kflip) =  prsl(i,k)
       END DO
    END DO
    DO k=1,kMax
       DO i=1,nCols    
          flip_pdel    (i,k) = MAX(flip_pint(i,k+1) - flip_pint(i,k),1.0e-12_r8)
          flip_rpdel   (i,k) = 1.0_r8/MAX((flip_pint(i,k+1) - flip_pint(i,k)),1.0e-12_r8)
          flip_lnpmid  (i,k) = LOG(flip_pmid(i,k))
       END DO
    END DO
    DO k=1,kMax+1
       DO i=1,nCols
          flip_lnpint(i,k) =  LOG(flip_pint  (i,k))
       END DO
    END DO

    !
    !..delsig     k=2  ****gu,gv,gt,gq,gyu,gyv,gtd,gqd,sl*** } delsig(2)
    !             k=3/2----si,ric,rf,km,kh,b,l -----------
    !             k=1  ****gu,gv,gt,gq,gyu,gyv,gtd,gqd,sl*** } delsig(1)
    !             k=1/2----si ----------------------------

    ! Derive new temperature and geopotential fields

    CALL geopotential_t(                                 &
         flip_lnpint(1:nCols,1:kMax+1)   , flip_pint  (1:nCols,1:kMax+1)  , &
         flip_pmid  (1:nCols,1:kMax)     , flip_pdel  (1:nCols,1:kMax)   , flip_rpdel(1:nCols,1:kMax)   , &
         flip_t     (1:nCols,1:kMax)     , flip_q     (1:nCols,1:kMax)   , rair   , gravit , zvir   ,&
         flip_zi    (1:nCols,1:kMax+1)   , flip_zm    (1:nCols,1:kMax)   , nCols   ,nCols, kMax)
    DO i=1,nCols
       zi(i,1) = flip_zi    (i,kMax+1)
    END DO
    DO k=1,kMax
       kflip=kMax+1-k
       DO i=1,nCols
          zi (i,k+1) = flip_zi    (i,kflip)
          zm (i,k  ) = flip_zm    (i,kflip)
          p  (i,k  ) = flip_pmid  (i,kflip)
       END DO
    END DO
    DO k=1,kMax
       DO i=1,nCols
          dz (i,k  ) = MAX(zi(i,k+1)-zi(i,k),1.0e-12)
          !j/kg/kelvin
          !
          ! P = rho * R * T
          !
          !            P
          ! rho  = -------
          !          R * T
          !
          rho   (i,k) =  (p(i,k)/(R*tc_mic(i,k)))       ! density
          w     (i,k) = -omega(i,k)/(rho(i,k)*gravit) ! (Pa/s)  - (m/s)
       END DO
    END DO
    IF(First_mic)THEN
       IF(.not.restart_mic)THEN
       IF (is_aerosol_aware) THEN

          !..Check for existing aerosol data, both CCN and IN aerosols.  If missing
          !.. fill in just a basic vertical profile, somewhat boundary-layer following.

          max_test = MAXVAL ( nwfa_mic(1:nCols-1,:))
          IF (max_test .LT. eps) THEN
             DO i =1, nCols
                IF (zm(i,1).LE.1000.0_r8) THEN
                   h_01 = 0.8_r8
                ELSEIF (zm(i,1).GE.2500.0_r8) THEN
                   h_01 = 0.01_r8
                ELSE
                   h_01 = 0.8_r8*COS(zm(i,1)*0.001_r8 - 1.0_r8)
                ENDIF
                niCCN3 = -1.0_r8*log(naCCN1/naCCN0)/h_01
                nwfa_mic(i,1) = naCCN1+naCCN0*EXP(-((zm(i,2)-zm(i,1))/1000.0_r8)*niCCN3)
                DO k = 2, kMAx
                   nwfa_mic(i,k) = naCCN1+naCCN0*EXP(-((zm(i,k)-zm(i,1))/1000.0_r8)*niCCN3)
                ENDDO
             ENDDO
          ENDIF

          max_test = MAXVAL ( nifa_mic(1:nCols-1,:) )

         IF (max_test .LT. eps) THEN
            DO i =1,  nCols
               IF (zm(i,1).LE.1000.0_r8) THEN
                  h_01 = 0.8_r8
               ELSEIF (zm(i,1).GE.2500.0_r8) THEN
                  h_01 = 0.01_r8
               ELSE
                  h_01 = 0.8_r8*COS(zm(i,1)*0.001_r8 - 1.0_r8)
               ENDIF
               niIN3 = -1.0_r8*log(naIN1/naIN0)/h_01
               nifa_mic(i,1) = naIN1+naIN0*EXP(-((zm(i,2)-zm(i,1))/1000.0_r8)*niIN3)
               DO k = 2, kMAx
                  nifa_mic(i,k) = naIN1+naIN0*EXP(-((zm(i,k)-zm(i,1))/1000.0_r8)*niIN3)
               ENDDO
            ENDDO
         ENDIF

         !..Capture initial state lowest level CCN aerosol data in 2D array.

        DO i = 1, nCols
           nwfa2d(i) = nwfa_mic(i,1)
        END DO

        !..Scale the lowest level aerosol data into an emissions rate.  This is
        !.. very far from ideal, but need higher emissions where larger amount
        !.. of existing and lesser emissions where not already lots of aerosols
        !.. for first-order simplistic approach.  Later, proper connection to
        !.. emission inventory would be better, but, for now, scale like this:
        !.. where: Nwfa=50 per cc, emit 0.875E4 aerosols per kg per second
        !..        Nwfa=500 per cc, emit 0.875E5 aerosols per kg per second
        !..        Nwfa=5000 per cc, emit 0.875E6 aerosols per kg per second
        !.. for a grid with 20km spacing and scale accordingly for other spacings.

        IF (is_start) THEN
           !IF (SQRT(DX*DY)/20000.0_r8 .GE. 1.0_r8) THEN
               h_01 = 0.875_r8
           !ELSE
           !    h_01 = (0.875_r8 + 0.125_r8*((20000.0_r8-SQRT(DX*DY))/16000.0_r8)) * SQRT(DX*DY)/20000.0_r8
           !ENDIF
           ! WRITE(mp_debug,*) '   aerosol surface flux emission scale factor is: ', h_01
           ! CALL wrf_debug(0, mp_debug)
           DO i = 1, nCols
              nwfa2d(i) = 10.0_r8**(LOG10(nwfa_mic(i,1)*1.E-6_r8)-3.69897_r8)
              nwfa2d(i) = nwfa2d(i)*h_01 * 1.E6_r8
           END DO
        ENDIF
      END IF
      ENDIF
      First_mic=.FALSE.
   END IF
   flip_t=tc_mic
   flip_q=qv_mic
   CALL mp_gt_driver(                                  &                        ! memory dims
                    nCols                            , & ! INTEGER , INTENT(IN   ) :: nCols
                    nLev_mic                         , & ! INTEGER , INTENT(IN   ) :: kMax
                    tc_mic    (1:nCols, 1:nLev_mic)  , & ! REAL(KIND=r8)  , INTENT(INOUT) :: tc_mic(1:nCols, 1:kMax)
                    qv_mic    (1:nCols, 1:nLev_mic)  , & ! REAL(KIND=r8)  , INTENT(INOUT) :: qv_mic(1:nCols, 1:kMax)
                    qc_mic    (1:nCols, 1:nLev_mic)  , & ! REAL(KIND=r8)  , INTENT(INOUT) :: qc_mic(1:nCols, 1:kMax)
                    qr_mic    (1:nCols, 1:nLev_mic)  , & ! REAL(KIND=r8)  , INTENT(INOUT) :: qr_mic(1:nCols, 1:kMax)
                    qi_mic    (1:nCols, 1:nLev_mic)  , & ! REAL(KIND=r8)  , INTENT(INOUT) :: qi_mic(1:nCols, 1:kMax)
                    qs_mic    (1:nCols, 1:nLev_mic)  , & ! REAL(KIND=r8)  , INTENT(INOUT) :: qs_mic(1:nCols, 1:kMax)
                    qg_mic    (1:nCols, 1:nLev_mic)  , & ! REAL(KIND=r8)  , INTENT(INOUT) :: qg_mic(1:nCols, 1:kMax)
                    ni_mic    (1:nCols, 1:nLev_mic)  , & ! REAL(KIND=r8)  , INTENT(INOUT) :: ni_mic(1:nCols, 1:kMax)
                    nr_mic    (1:nCols, 1:nLev_mic)  , & ! REAL(KIND=r8)  , INTENT(INOUT) :: nr_mic(1:nCols, 1:kMax)
                    p         (1:nCols, 1:nLev_mic)  , & ! REAL(KIND=r8)  , INTENT(IN   ) :: p  (1:nCols, 1:kMax)
                    w         (1:nCols, 1:nLev_mic)  , & ! REAL(KIND=r8)  , INTENT(IN   ) :: w  (1:nCols, 1:kMax)
                    dz        (1:nCols, 1:nLev_mic)  , & ! REAL(KIND=r8)  , INTENT(IN   ) :: dz (1:nCols, 1:kMax)
                    dt_in                            , & ! REAL(KIND=r8)  , INTENT(IN   ):: dt_in
                    RAINNC    (1:nCols)              , & ! REAL(KIND=r8)  , INTENT(INOUT) :: RAINNC  (1:nCols)
                    RAINNCV   (1:nCols)              , & ! REAL(KIND=r8)  , INTENT(INOUT) :: RAINNCV (1:nCols)
                    SNOWNC    (1:nCols)              , & ! REAL(KIND=r8)  , OPTIONAL, INTENT(INOUT) :: SNOWNC    (1:nCols)
                    SNOWNCV   (1:nCols)              , & ! REAL(KIND=r8)  , OPTIONAL, INTENT(INOUT) :: SNOWNCV   (1:nCols)
                    GRAUPELNC (1:nCols)              , & ! REAL(KIND=r8)  , OPTIONAL, INTENT(INOUT) :: GRAUPELNC (1:nCols)
                    GRAUPELNCV(1:nCols)              , & ! REAL(KIND=r8)  , OPTIONAL, INTENT(INOUT) :: GRAUPELNCV(1:nCols)
                    SR        (1:nCols)              , & ! REAL(KIND=r8)  , INTENT(INOUT) :: SR   (1:nCols)
                                                         !#if ( WRF_CHEM == 1 )
                    rainprod  (1:nCols, 1:nLev_mic)  , & ! REAL(KIND=r8)  , INTENT(INOUT) :: rainprod(1:nCols, 1:kMax)
                    evapprod  (1:nCols, 1:nLev_mic)  , & ! REAL(KIND=r8)  , INTENT(INOUT) :: evapprod(1:nCols, 1:kMax)
                                                         !#endif
                    refl_10cm (1:nCols, 1:nLev_mic)  , & ! REAL(KIND=r8)  , INTENT(INOUT):: refl_10cm(1:nCols, 1:kMax)
                    diagflag                         , & ! LOGICAL             , OPTIONAL, INTENT(IN) :: diagflag
                    do_radar_ref                     , & ! INTEGER             , OPTIONAL, INTENT(IN) :: do_radar_ref
                    EFFCS     (1:nCols, 1:nLev_mic)  , & ! REAL(KIND=r8)  , INTENT(INOUT) :: re_cloud(1:nCols, 1:kMax)
                    EFFIS     (1:nCols, 1:nLev_mic)  , & ! REAL(KIND=r8)  , INTENT(INOUT) :: re_ice  (1:nCols, 1:kMax)
                    re_snow   (1:nCols, 1:nLev_mic)  , & ! REAL(KIND=r8)  , INTENT(INOUT) :: re_snow (1:nCols, 1:kMax)
                    has_reqc                         , & ! INTEGER             , INTENT(IN   ) :: has_reqc
                    has_reqi                         , & ! INTEGER             , INTENT(IN   ) :: has_reqi
                    has_reqs                         , & ! INTEGER             , INTENT(IN   ) :: has_reqs
                    nc_mic        (1:nCols, 1:nLev_mic)  , & ! REAL(KIND=r8)  , OPTIONAL, INTENT(INOUT) :: nc_mic  (1:nCols, 1:kMax)
                    nwfa_mic      (1:nCols, 1:nLev_mic)  , & ! REAL(KIND=r8)  , OPTIONAL, INTENT(INOUT) :: nwfa_mic(1:nCols, 1:kMax)
                    nifa_mic      (1:nCols, 1:nLev_mic)  , & ! REAL(KIND=r8)  , OPTIONAL, INTENT(INOUT) :: nifa_mic(1:nCols, 1:kMax)
                    nwfa2d    (1:nCols)           )  ! REAL(KIND=r8)  , OPTIONAL, INTENT(IN   ) :: nwfa2d(1:nCols)


    DO k=1,kMax
       DO i=1,nCols

          dTcdt   (i, k)=(Tc_mic   (i, k) -  Tc   (i, k))/dt_in
          dqvdt   (i, k)=(qv_mic   (i, k) -  qv   (i, k))/dt_in
          dqcdt   (i, k)=(qc_mic   (i, k) -  qc   (i, k))/dt_in
          dqrdt   (i, k)=(qr_mic   (i, k) -  qr   (i, k))/dt_in
          dqidt   (i, k)=(qi_mic   (i, k) -  qi   (i, k))/dt_in
          dqsdt   (i, k)=(qs_mic   (i, k) -  qs   (i, k))/dt_in
          dqgdt   (i, k)=(qg_mic   (i, k) -  qg   (i, k))/dt_in
          dnidt   (i, k)=(ni_mic   (i, k) -  ni   (i, k))/dt_in
          dnrdt   (i, k)=(nr_mic   (i, k) -  nr   (i, k))/dt_in
          dnifadt (i, k)=(nifa_mic (i, k) -  nifa (i, k))/dt_in
          dnwfadt (i, k)=(nwfa_mic (i, k) -  nwfa (i, k))/dt_in
          dncdt   (i, k)=(nc_mic   (i, k) -  nc   (i, k))/dt_in


          Tc   (i, k) = Tc_mic   (i, k)
          qv   (i, k) = qv_mic   (i, k)
          qc   (i, k) = qc_mic   (i, k)
          qr   (i, k) = qr_mic   (i, k)
          qi   (i, k) = qi_mic   (i, k)
          qs   (i, k) = qs_mic   (i, k)
          qg   (i, k) = qg_mic   (i, k)
          ni   (i, k) = ni_mic   (i, k)
          nr   (i, k) = nr_mic   (i, k)
          nifa (i, k) = nifa_mic (i, k)
          nwfa (i, k) = nwfa_mic (i, k)
          nc   (i, k) = nc_mic   (i, k)
       END DO
    END DO   

    !DO k=1,kMax
    !   Do i=1,nCols
    !      PRINT*,qc(i,k),qr(i,k),(qs(i,k)),(qv(i,k)-flip_q(i,k))/dt_in,p(i,k),rainnc(i)
    !   END DO
    !END DO   
    LSRAIN(1:nCols)=0.5_r8*RAINNC(1:nCols)/1000.0_r8  !(mm)->m
    LSSNOW(1:nCols)=0.5_r8*SNOWNC(1:nCols)/1000.0_r8  !(mm)->m
  END SUBROUTINE RunMicro_thompson

  !+---+-----------------------------------------------------------------+
  !ctrlL
  !+---+-----------------------------------------------------------------+
  !..This is a wrapper routine designed to transfer values from 3D to 1D.
  !+---+-----------------------------------------------------------------+
  SUBROUTINE mp_gt_driver(              &                              ! memory dims
                          nCols       , & ! INTEGER , INTENT(IN   ) :: nCols
                          kMax        , & ! INTEGER , INTENT(IN   ) :: kMax
                          tc          , & ! REAL(KIND=r8)  , INTENT(INOUT) :: tc(1:nCols, 1:kMax)
                          qv          , & ! REAL(KIND=r8)  , INTENT(INOUT) :: qv(1:nCols, 1:kMax)
                          qc          , & ! REAL(KIND=r8)  , INTENT(INOUT) :: qc(1:nCols, 1:kMax)
                          qr          , & ! REAL(KIND=r8)  , INTENT(INOUT) :: qr(1:nCols, 1:kMax)
                          qi          , & ! REAL(KIND=r8)  , INTENT(INOUT) :: qi(1:nCols, 1:kMax)
                          qs          , & ! REAL(KIND=r8)  , INTENT(INOUT) :: qs(1:nCols, 1:kMax)
                          qg          , & ! REAL(KIND=r8)  , INTENT(INOUT) :: qg(1:nCols, 1:kMax)
                          ni          , & ! REAL(KIND=r8)  , INTENT(INOUT) :: ni(1:nCols, 1:kMax)
                          nr          , & ! REAL(KIND=r8)  , INTENT(INOUT) :: nr(1:nCols, 1:kMax)
                          p           , & ! REAL(KIND=r8)  , INTENT(IN   ) :: p  (1:nCols, 1:kMax)
                          w           , & ! REAL(KIND=r8)  , INTENT(IN   ) :: w  (1:nCols, 1:kMax)
                          dz          , & ! REAL(KIND=r8)  , INTENT(IN   ) :: dz (1:nCols, 1:kMax)
                          dt_in       , & ! REAL(KIND=r8)  , INTENT(IN   ):: dt_in
                          RAINNC      , & ! REAL(KIND=r8)  , INTENT(INOUT) :: RAINNC  (1:nCols)
                          RAINNCV     , & ! REAL(KIND=r8)  , INTENT(INOUT) :: RAINNCV (1:nCols)
                          SNOWNC      , & ! REAL(KIND=r8)  , OPTIONAL, INTENT(INOUT) :: SNOWNC    (1:nCols)
                          SNOWNCV     , & ! REAL(KIND=r8)  , OPTIONAL, INTENT(INOUT) :: SNOWNCV   (1:nCols)
                          GRAUPELNC   , & ! REAL(KIND=r8)  , OPTIONAL, INTENT(INOUT) :: GRAUPELNC (1:nCols)
                          GRAUPELNCV  , & ! REAL(KIND=r8)  , OPTIONAL, INTENT(INOUT) :: GRAUPELNCV(1:nCols)
                          SR          , & ! REAL(KIND=r8)  , INTENT(INOUT) :: SR   (1:nCols)
                                          !#if ( WRF_CHEM == 1 )
                          rainprod    , & ! REAL(KIND=r8)  , INTENT(INOUT) :: rainprod(1:nCols, 1:kMax)
                          evapprod    , & ! REAL(KIND=r8)  , INTENT(INOUT) :: evapprod(1:nCols, 1:kMax)
                                          !#endif
                          refl_10cm   , & ! REAL(KIND=r8)  , INTENT(INOUT):: refl_10cm(1:nCols, 1:kMax)
                          diagflag    , & ! LOGICAL        , OPTIONAL, INTENT(IN) :: diagflag
                          do_radar_ref, & ! INTEGER        , OPTIONAL, INTENT(IN) :: do_radar_ref
                          re_cloud    , & ! REAL(KIND=r8)  , INTENT(INOUT) :: re_cloud(1:nCols, 1:kMax)
                          re_ice      , & ! REAL(KIND=r8)  , INTENT(INOUT) :: re_ice  (1:nCols, 1:kMax)
                          re_snow     , & ! REAL(KIND=r8)  , INTENT(INOUT) :: re_snow (1:nCols, 1:kMax)
                          has_reqc    , & ! INTEGER        , INTENT(IN   ) :: has_reqc
                          has_reqi    , & ! INTEGER        , INTENT(IN   ) :: has_reqi
                          has_reqs    , & ! INTEGER        , INTENT(IN   ) :: has_reqs
                          nc          , & ! REAL(KIND=r8)  , OPTIONAL, INTENT(INOUT) :: nc  (1:nCols, 1:kMax)
                          nwfa        , & ! REAL(KIND=r8)  , OPTIONAL, INTENT(INOUT) :: nwfa(1:nCols, 1:kMax)
                          nifa        , & ! REAL(KIND=r8)  , OPTIONAL, INTENT(INOUT) :: nifa(1:nCols, 1:kMax)
                          nwfa2d         ) ! REAL(KIND=r8)  , OPTIONAL, INTENT(IN   ) :: nwfa2d(1:nCols)

    IMPLICIT NONE

    !..Subroutine arguments

    INTEGER        , INTENT(IN   ) :: nCols
    INTEGER        , INTENT(IN   ) :: kMax
    REAL(KIND=r8)  , INTENT(INOUT) :: qv(1:nCols, 1:kMax)
    REAL(KIND=r8)  , INTENT(INOUT) :: qc(1:nCols, 1:kMax)
    REAL(KIND=r8)  , INTENT(INOUT) :: qr(1:nCols, 1:kMax)
    REAL(KIND=r8)  , INTENT(INOUT) :: qi(1:nCols, 1:kMax)
    REAL(KIND=r8)  , INTENT(INOUT) :: qs(1:nCols, 1:kMax)
    REAL(KIND=r8)  , INTENT(INOUT) :: qg(1:nCols, 1:kMax)
    REAL(KIND=r8)  , INTENT(INOUT) :: ni(1:nCols, 1:kMax)
    REAL(KIND=r8)  , INTENT(INOUT) :: nr(1:nCols, 1:kMax)
    REAL(KIND=r8)  , INTENT(INOUT) :: tc(1:nCols, 1:kMax)
    REAL(KIND=r8)  , INTENT(IN   ) :: p  (1:nCols, 1:kMax)
    REAL(KIND=r8)  , INTENT(IN   ) :: w  (1:nCols, 1:kMax)
    REAL(KIND=r8)  , INTENT(IN   ) :: dz (1:nCols, 1:kMax)
    REAL(KIND=r8)  , INTENT(IN   ) :: dt_in
    REAL(KIND=r8)  , INTENT(INOUT) :: RAINNC  (1:nCols)
    REAL(KIND=r8)  , INTENT(INOUT) :: RAINNCV (1:nCols)
    REAL(KIND=r8)  , INTENT(INOUT) :: SNOWNC  (1:nCols)
    REAL(KIND=r8)  , INTENT(INOUT) :: SNOWNCV   (1:nCols)
    REAL(KIND=r8)  , INTENT(INOUT) :: GRAUPELNC (1:nCols)
    REAL(KIND=r8)  , INTENT(INOUT) :: GRAUPELNCV(1:nCols)
    REAL(KIND=r8)  , INTENT(INOUT) :: SR      (1:nCols)
    !#if ( WRF_CHEM == 1 )
    REAL(KIND=r8)  , INTENT(INOUT) :: rainprod(1:nCols, 1:kMax)
    REAL(KIND=r8)  , INTENT(INOUT) :: evapprod(1:nCols, 1:kMax)
    !#endif
    REAL(KIND=r8)  , INTENT(INOUT):: refl_10cm(1:nCols, 1:kMax)
    LOGICAL        , INTENT(IN) :: diagflag
    INTEGER        , INTENT(IN) :: do_radar_ref
    REAL(KIND=r8)  , INTENT(INOUT) :: re_cloud(1:nCols, 1:kMax)
    REAL(KIND=r8)  , INTENT(INOUT) :: re_ice  (1:nCols, 1:kMax)
    REAL(KIND=r8)  , INTENT(INOUT) :: re_snow (1:nCols, 1:kMax)
    INTEGER        , INTENT(IN   ) :: has_reqc
    INTEGER        , INTENT(IN   ) :: has_reqi
    INTEGER        , INTENT(IN   ) :: has_reqs
    REAL(KIND=r8)  , OPTIONAL, INTENT(INOUT) :: nc  (1:nCols, 1:kMax)
    REAL(KIND=r8)  , OPTIONAL, INTENT(INOUT) :: nwfa(1:nCols, 1:kMax)
    REAL(KIND=r8)  , OPTIONAL, INTENT(INOUT) :: nifa(1:nCols, 1:kMax)
    REAL(KIND=r8)  , OPTIONAL, INTENT(IN   ) :: nwfa2d(1:nCols)

    !..Local variables
    REAL(KIND=r8) :: qv1d(1:kMax)
    REAL(KIND=r8) :: qc1d(1:kMax)
    REAL(KIND=r8) :: qi1d(1:kMax)
    REAL(KIND=r8) :: qr1d(1:kMax)
    REAL(KIND=r8) :: qs1d(1:kMax)
    REAL(KIND=r8) :: qg1d(1:kMax)
    REAL(KIND=r8) :: ni1d(1:kMax)
    REAL(KIND=r8) :: nr1d(1:kMax)
    REAL(KIND=r8) :: nc1d(1:kMax)
    REAL(KIND=r8) :: nwfa1d(1:kMax)
    REAL(KIND=r8) :: nifa1d(1:kMax)
    REAL(KIND=r8) :: t1d(1:kMax)
    REAL(KIND=r8) :: p1d(1:kMax)
    REAL(KIND=r8) :: w1d(1:kMax)
    REAL(KIND=r8) :: dz1d(1:kMax)
    REAL(KIND=r8) :: dBZ(1:kMax)
    REAL(KIND=r8) :: re_qc1d(1:kMax)
    REAL(KIND=r8) :: re_qi1d(1:kMax)
    REAL(KIND=r8) :: re_qs1d(1:kMax)
    !#if ( WRF_CHEM == 1 )
    REAL(KIND=r8):: rainprod1d(1:kMax)
    REAL(KIND=r8):: evapprod1d(1:kMax)
    !#endif
    REAL(KIND=r8):: pcp_ra(1:nCols)
    REAL(KIND=r8):: pcp_sn(1:nCols)
    REAL(KIND=r8):: pcp_gr(1:nCols)
    REAL(KIND=r8):: pcp_ic(1:nCols)
    REAL(KIND=r8):: dt
    REAL(KIND=r8):: pptrain
    REAL(KIND=r8):: pptsnow
    REAL(KIND=r8):: pptgraul
    REAL(KIND=r8):: pptice
    REAL(KIND=r8):: qc_max
    REAL(KIND=r8):: qr_max
    REAL(KIND=r8):: qs_max
    REAL(KIND=r8):: qi_max
    REAL(KIND=r8):: qg_max
    REAL(KIND=r8):: ni_max
    REAL(KIND=r8):: nr_max
    REAL(KIND=r8):: nwfa1
    INTEGER :: i
!    INTEGER :: j
    INTEGER :: k
    INTEGER :: imax_qc
    INTEGER :: imax_qr
    INTEGER :: imax_qi
    INTEGER :: imax_qs
    INTEGER :: imax_qg
    INTEGER :: imax_ni
    INTEGER :: imax_nr
    INTEGER :: jmax_qc
    INTEGER :: jmax_qr
    INTEGER :: jmax_qi
    INTEGER :: jmax_qs
    INTEGER :: jmax_qg
    INTEGER :: jmax_ni
    INTEGER :: jmax_nr
    INTEGER :: kmax_qc
    INTEGER :: kmax_qr
    INTEGER :: kmax_qi
    INTEGER :: kmax_qs
    INTEGER :: kmax_qg
    INTEGER :: kmax_ni
    INTEGER :: kmax_nr
    INTEGER :: i_start
!    INTEGER :: j_start
    INTEGER :: i_end
!    INTEGER :: j_end
    CHARACTER*256:: mp_debug

    !+---+

    i_start = 1
    !j_start = jts
    i_end   = nCols
    !j_end   = MIN(jte, jde-1)

    !..For idealized testing by developer.
    !     if ( (ide-ids+1).gt.4 .and. (jde+1).lt.4 .and.                &
    !          ids.eq.its.and.ide.eq.nCols.and.eq.and.jde.eq.jte) then
    !        i_start = its + 2
    !        i_end   = ite
    !        j_start 
    !        j_end   = jte
    !     endif

    dt = dt_in

    qc_max = 0.0_r8
    qr_max = 0.0_r8
    qs_max = 0.0_r8
    qi_max = 0.0_r8
    qg_max = 0
    ni_max = 0.0_r8
    nr_max = 0.0_r8
    imax_qc = 0
    imax_qr = 0
    imax_qi = 0
    imax_qs = 0
    imax_qg = 0
    imax_ni = 0
    imax_nr = 0
    jmax_qc = 0
    jmax_qr = 0
    jmax_qi = 0
    jmax_qs = 0
    jmax_qg = 0
    jmax_ni = 0
    jmax_nr = 0
    kmax_qc = 0
    kmax_qr = 0
    kmax_qi = 0
    kmax_qs = 0
    kmax_qg = 0
    kmax_ni = 0
    kmax_nr = 0
    DO i = 1, 256
       mp_debug(i:i) = CHAR(0)
    ENDDO

    IF (.NOT. is_aerosol_aware .AND. PRESENT(nc) .AND. PRESENT(nwfa)  &
         .AND. PRESENT(nifa) .AND. PRESENT(nwfa2d)) THEN
       WRITE(mp_debug,*) 'WARNING, nc-nwfa-nifa-nwfa2d present but is_aerosol_aware is FALSE'
       CALL wrf_debug('..mp_gt_driver..', mp_debug)
    ENDIF

    !      j_loop:  do j = j_start, j_end
    i_loop:  DO i = 1, nCols

       pptrain = 0.0_r8
       pptsnow = 0.0_r8
       pptgraul = 0.0_r8
       pptice = 0.0_r8
       RAINNCV(i) = 0.0_r8
      ! IF ( PRESENT (snowncv) ) THEN
          SNOWNCV(i) = 0.0_r8
      ! ENDIF
      ! IF ( PRESENT (graupelncv) ) THEN
          GRAUPELNCV(i) = 0.0_r8
      ! ENDIF
       SR(i) = 0.0_r8

       DO k = 1, kMax
          IF(tc(i,k) < 0.0_r8)THEN
              PRINT*,i,k,tc(i,k)
           ELSe 
          t1d(k) = tc(i,k)!*pii(i,k)
          END IF
          p1d(k) = p(i,k)
          w1d(k) = w(i,k)
          dz1d(k) = dz(i,k)
          qv1d(k) = qv(i,k)
          qc1d(k) = qc(i,k)
          qi1d(k) = qi(i,k)
          qr1d(k) = qr(i,k)
          qs1d(k) = qs(i,k)
          qg1d(k) = qg(i,k)
          ni1d(k) = ni(i,k)
          nr1d(k) = nr(i,k)
       ENDDO
       IF (is_aerosol_aware) THEN
          DO k = 1, kMax
             nc1d(k)   = nc(i,k)
             nwfa1d(k) = nwfa(i,k)
             nifa1d(k) = nifa(i,k)
          ENDDO
          nwfa1 = nwfa2d(i)
       ELSE
          DO k = 1, kMax
             nc1d(k) = Nt_c
             nwfa1d(k) = 11.1E6_r8
             nifa1d(k) = naIN1*0.01_r8
          ENDDO
          nwfa1 = 11.1E6_r8
       ENDIF

       CALL mp_thompson(&
            qv1d      (1:kMax), &  !  REAL(KIND=r8), INTENT(INOUT) :: qv1d(kts:kte)
            qc1d      (1:kMax), &  !  REAL(KIND=r8), INTENT(INOUT) :: qc1d(kts:kte)
            qi1d      (1:kMax), &  !  REAL(KIND=r8), INTENT(INOUT) :: qi1d(kts:kte)
            qr1d      (1:kMax), &  !  REAL(KIND=r8), INTENT(INOUT) :: qr1d(kts:kte)
            qs1d      (1:kMax), &  !  REAL(KIND=r8), INTENT(INOUT) :: qs1d(kts:kte)
            qg1d      (1:kMax), &  !  REAL(KIND=r8), INTENT(INOUT) :: qg1d(kts:kte)
            ni1d      (1:kMax), &  !  REAL(KIND=r8), INTENT(INOUT) :: ni1d(kts:kte)
            nr1d      (1:kMax), &  !  REAL(KIND=r8), INTENT(INOUT) :: nr1d(kts:kte)
            nc1d      (1:kMax), &  !  REAL(KIND=r8), INTENT(INOUT) :: nc1d(kts:kte)
            nwfa1d    (1:kMax), &  !  REAL(KIND=r8), INTENT(INOUT) :: nwfa1d(kts:kte)
            nifa1d    (1:kMax), &  !  REAL(KIND=r8), INTENT(INOUT) :: nifa1d(kts:kte)
            t1d       (1:kMax), &  !  REAL(KIND=r8), INTENT(INOUT) :: t1d(kts:kte)
            p1d       (1:kMax), &  !REAL(KIND=r8), INTENT(IN  ) :: p1d(kts:kte)
            w1d       (1:kMax), &  !REAL(KIND=r8), INTENT(IN  ) :: w1d(kts:kte)
            dz1d      (1:kMax), &  !REAL(KIND=r8), INTENT(IN  ) :: dzq(kts:kte)
            pptrain           , &  !REAL(KIND=r8), INTENT(INOUT) :: pptrain
            pptsnow           , &  !REAL(KIND=r8), INTENT(INOUT) :: pptsnow
            pptgraul          , &  !REAL(KIND=r8), INTENT(INOUT) :: pptgraul
            pptice            , &  !REAL(KIND=r8), INTENT(INOUT) :: pptice
                                !#if ( WRF_CHEM ==     !#if ( WRF_CHEM == 1 )
            rainprod1d(1:kMax), &  ! REAL(KIND=r8), INTENT(INOUT) :: rainprod(kts:kte)
            evapprod1d(1:kMax), &  !REAL(KIND=r8), INTENT(INOUT) :: evapprod(kts:kte)
                                !#endif        !#endif
            1                 , &   !  INTEGER, INTENT(IN) :: kts
            kMax              , &  ! INTEGER, INTENT(IN) :: kte
            dt                , &  !REAL(KIND=r8), INTENT(IN   ) :: dt
            i                   )  !INTEGER, INTENT(IN) :: ii

       pcp_ra(i) = pptrain
       pcp_sn(i) = pptsnow
       pcp_gr(i) = pptgraul
       pcp_ic(i) = pptice
       RAINNCV(i) = pptrain + pptsnow + pptgraul + pptice
       RAINNC(i) = RAINNC(i) + pptrain + pptsnow + pptgraul + pptice
       !IF ( PRESENT(snowncv) .AND. PRESENT(snownc) ) THEN
          SNOWNCV(i) = pptsnow + pptice
          SNOWNC(i) = SNOWNC(i) + pptsnow + pptice
       !ENDIF
       !IF ( PRESENT(graupelncv) .AND. PRESENT(graupelnc) ) THEN
          GRAUPELNCV(i) = pptgraul
          GRAUPELNC(i) = GRAUPELNC(i) + pptgraul
       !ENDIF
       SR(i) = (pptsnow + pptgraul + pptice)/(RAINNCV(i)+1.e-12_r8)



       !..Reset lowest model level to initial state aerosols (fake sfc source).
       !.. Changed 13 May 2013 to fake emissions in which nwfa2d is aerosol
       !.. number tendency (number per kg per second).
       IF (is_aerosol_aware) THEN
          !-GT        nwfa1d(1) = nwfa1
          nwfa1d(1) = nwfa1d(1) + nwfa2d(i)*dt_in

          DO k = 1, kMax
             nc(i,k) = nc1d(k)
             nwfa(i,k) = nwfa1d(k)
             nifa(i,k) = nifa1d(k)
          ENDDO
       ENDIF

       DO k = 1, kMax
          qv(i,k) = qv1d(k)
          qc(i,k) = qc1d(k)
          qi(i,k) = qi1d(k)
          qr(i,k) = qr1d(k)
          qs(i,k) = qs1d(k)
          qg(i,k) = qg1d(k)
          ni(i,k) = ni1d(k)
          nr(i,k) = nr1d(k)
          tc(i,k) = t1d(k)!/pii(i,k)
          !if ( WRF_CHEM == 1 )
          rainprod(i,k) = rainprod1d(k)
          evapprod(i,k) = evapprod1d(k)
          !endif
          IF (qc1d(k) .GT. qc_max) THEN
             imax_qc = i
             kmax_qc = k
             qc_max = qc1d(k)
          ELSEIF (qc1d(k) .LT. 0.0_r8) THEN
             WRITE(mp_debug,*) 'WARNING, negative qc ', qc1d(k),        &
                  ' at i,k=', i,k
             CALL wrf_debug('..mp_gt_driver..', mp_debug)
          ENDIF
          IF (qr1d(k) .GT. qr_max) THEN
             imax_qr = i
             kmax_qr = k
             qr_max = qr1d(k)
          ELSEIF (qr1d(k) .LT. 0.0_r8) THEN
             WRITE(mp_debug,*) 'WARNING, negative qr ', qr1d(k),        &
                  ' at i,k=', i,k
             CALL wrf_debug('..mp_gt_driver..', mp_debug)
          ENDIF
          IF (nr1d(k) .GT. nr_max) THEN
             imax_nr = i
             kmax_nr = k
             nr_max = nr1d(k)
          ELSEIF (nr1d(k) .LT. 0.0_r8) THEN
             WRITE(mp_debug,*) 'WARNING, negative nr ', nr1d(k),        &
                  ' at i,k=', i,k
             CALL wrf_debug('..mp_gt_driver..', mp_debug)
          ENDIF
          IF (qs1d(k) .GT. qs_max) THEN
             imax_qs = i
             kmax_qs = k
             qs_max = qs1d(k)
          ELSEIF (qs1d(k) .LT. 0.0_r8) THEN
             WRITE(mp_debug,*) 'WARNING, negative qs ', qs1d(k),        &
                  ' at i,k=', i,k
             CALL wrf_debug('..mp_gt_driver..', mp_debug)
          ENDIF
          IF (qi1d(k) .GT. qi_max) THEN
             imax_qi = i
             kmax_qi = k
             qi_max = qi1d(k)
          ELSEIF (qi1d(k) .LT. 0.0_r8) THEN
             WRITE(mp_debug,*) 'WARNING, negative qi ', qi1d(k),        &
                  ' at i,k=', i,k
             CALL wrf_debug('..mp_gt_driver..', mp_debug)
          ENDIF
          IF (qg1d(k) .GT. qg_max) THEN
             imax_qg = i
             kmax_qg = k
             qg_max = qg1d(k)
          ELSEIF (qg1d(k) .LT. 0.0_r8) THEN
             WRITE(mp_debug,*) 'WARNING, negative qg ', qg1d(k),        &
                  ' at i,k=', i,k
             CALL wrf_debug('..mp_gt_driver..', mp_debug)
          ENDIF
          IF (ni1d(k) .GT. ni_max) THEN
             imax_ni = i
             kmax_ni = k
             ni_max = ni1d(k)
          ELSEIF (ni1d(k) .LT. 0.0_r8) THEN
             WRITE(mp_debug,*) 'WARNING, negative ni ', ni1d(k),        &
                  ' at i,k=', i,k
             CALL wrf_debug('..mp_gt_driver..', mp_debug)
          ENDIF
          IF (qv1d(k) .LT. 0.0_r8) THEN
             WRITE(mp_debug,*) 'WARNING, negative qv ', qv1d(k),        &
                  ' at i,k=', i,k
             CALL wrf_debug('..mp_gt_driver..', mp_debug)
             IF (k.LT.kMax-2 .AND. k.GT.1+1) THEN
                WRITE(mp_debug,*) '   below and above are: ', qv(i,k-1), qv(i,k+1)
                CALL wrf_debug('..mp_gt_driver..', mp_debug)
                qv(i,k) = MAX(1.E-7_r8, 0.5_r8*(qv(i,k-1) + qv(i,k+1)))
             ELSE
                qv(i,k) = 1.E-7_r8
             ENDIF
          ENDIF
       ENDDO

 !      IF ( PRESENT (diagflag) ) THEN
          IF (diagflag .AND. do_radar_ref == 1) THEN
             CALL calc_refl10cm (qv1d, qc1d, qr1d, nr1d, qs1d, qg1d,       &
                  t1d, p1d, dBZ, 1, kMax)
             DO k = 1, kMax
                refl_10cm(i,k) = MAX(-35.0_r8, dBZ(k))
             ENDDO
          ENDIF
!       ENDIF

       IF (has_reqc.NE.0 .AND. has_reqi.NE.0 .AND. has_reqs.NE.0) THEN
          DO k = 1, kMax
             re_qc1d(k) = 2.51E-6_r8
             re_qi1d(k) = 10.01E-6_r8
             re_qs1d(k) = 10.01E-6_r8
          ENDDO
          CALL calc_effectRad (t1d, p1d, qv1d, qc1d, nc1d, qi1d, ni1d, qs1d,  &
               re_qc1d, re_qi1d, re_qs1d, 1, kMax)
          DO k = 1, kMax
             re_cloud(i,k) = MAX( 2.51E-6_r8, MIN(re_qc1d(k), 50.E-6_r8))
             re_ice  (i,k) = MAX(10.01E-6_r8, MIN(re_qi1d(k), 125.E-6_r8))
             re_snow (i,k) = MAX(10.01E-6_r8, MIN(re_qs1d(k), 999.E-6_r8))
          ENDDO
       ENDIF

    ENDDO i_loop
    !      enddo j_loop

    ! DEBUG - GT
!    WRITE(mp_debug,'(a,7(a,e13.6,1x,a,i3,a,i3,a,i3,a,1x))') 'MP-GT:', &
!         'qc: ', qc_max, '(', imax_qc, ',', jmax_qc, ',', kmax_qc, ')', &
!         'qr: ', qr_max, '(', imax_qr, ',', jmax_qr, ',', kmax_qr, ')', &
!         'qi: ', qi_max, '(', imax_qi, ',', jmax_qi, ',', kmax_qi, ')', &
!         'qs: ', qs_max, '(', imax_qs, ',', jmax_qs, ',', kmax_qs, ')', &
!         'qg: ', qg_max, '(', imax_qg, ',', jmax_qg, ',', kmax_qg, ')', &
!         'ni: ', ni_max, '(', imax_ni, ',', jmax_ni, ',', kmax_ni, ')', &
!         'nr: ', nr_max, '(', imax_nr, ',', jmax_nr, ',', kmax_nr, ')'
!    CALL wrf_debug('..mp_gt_driver..', mp_debug)
    ! END DEBUG - GT

!    DO i = 1, 256
!       mp_debug(i:i) = CHAR(0)
!    ENDDO

  END SUBROUTINE mp_gt_driver

  !+---+-----------------------------------------------------------------+
  !ctrlL
  !+---+-----------------------------------------------------------------+
  !+---+-----------------------------------------------------------------+
  !.. This subroutine computes the moisture tendencies of water vapor,
  !.. cloud droplets, rain, cloud ice (pristine), snow, and graupel.
  !.. Previously this code was based on Reisner et al (1998), but few of
  !.. those pieces remain.  A complete description is now found in
  !.. Thompson et al. (2004, 2008).
  !+---+-----------------------------------------------------------------+
  !
  SUBROUTINE mp_thompson (&
       qv1d       , &!  REAL(KIND=r8), INTENT(INOUT) :: qv1d(kts:kte)
       qc1d       , &!  REAL(KIND=r8), INTENT(INOUT) :: qc1d(kts:kte)
       qi1d       , &!  REAL(KIND=r8), INTENT(INOUT) :: qi1d(kts:kte)
       qr1d       , &!  REAL(KIND=r8), INTENT(INOUT) :: qr1d(kts:kte)
       qs1d       , &!  REAL(KIND=r8), INTENT(INOUT) :: qs1d(kts:kte)
       qg1d       , &!  REAL(KIND=r8), INTENT(INOUT) :: qg1d(kts:kte)
       ni1d       , &!  REAL(KIND=r8), INTENT(INOUT) :: ni1d(kts:kte)
       nr1d       , &!  REAL(KIND=r8), INTENT(INOUT) :: nr1d(kts:kte)
       nc1d       , &!  REAL(KIND=r8), INTENT(INOUT) :: nc1d(kts:kte)
       nwfa1d     , &!  REAL(KIND=r8), INTENT(INOUT) :: nwfa1d(kts:kte)
       nifa1d     , &!  REAL(KIND=r8), INTENT(INOUT) :: nifa1d(kts:kte)
       t1d        , &!REAL(KIND=r8), INTENT(INOUT) :: t1d(kts:kte)
       p1d        , &!REAL(KIND=r8), INTENT(IN) :: p1d(kts:kte)
       w1d        , &!REAL(KIND=r8), INTENT(IN) :: w1d(kts:kte)
       dzq        , &!REAL(KIND=r8), INTENT(IN) :: dzq(kts:kte)
       pptrain    , &!REAL(KIND=r8), INTENT(INOUT) :: pptrain
       pptsnow    , &!REAL(KIND=r8), INTENT(INOUT) :: pptsnow
       pptgraul   , &!REAL(KIND=r8), INTENT(INOUT) :: pptgraul
       pptice     , &!REAL(KIND=r8), INTENT(INOUT) :: pptice
                                !#if ( WRF_CHEM == 1 )
       rainprod   , &! REAL(KIND=r8), INTENT(INOUT) :: rainprod(kts:kte)
       evapprod   , &!REAL(KIND=r8), INTENT(INOUT) :: evapprod(kts:kte)
                                !#endif
       kts        , &!  INTEGER, INTENT(IN) :: kts
       kte        , &! INTEGER, INTENT(IN) :: kte
       dt         , &!REAL(KIND=r8), INTENT(IN   ) :: dt
       ii           )!INTEGER, INTENT(IN) :: ii

    IMPLICIT NONE

    !..Sub arguments
    INTEGER, INTENT(IN) :: kts
    INTEGER, INTENT(IN) :: kte
    INTEGER, INTENT(IN) :: ii
    REAL(KIND=r8), INTENT(INOUT) :: qv1d(kts:kte)
    REAL(KIND=r8), INTENT(INOUT) :: qc1d(kts:kte)
    REAL(KIND=r8), INTENT(INOUT) :: qi1d(kts:kte)
    REAL(KIND=r8), INTENT(INOUT) :: qr1d(kts:kte)
    REAL(KIND=r8), INTENT(INOUT) :: qs1d(kts:kte)
    REAL(KIND=r8), INTENT(INOUT) :: qg1d(kts:kte)
    REAL(KIND=r8), INTENT(INOUT) :: ni1d(kts:kte)
    REAL(KIND=r8), INTENT(INOUT) :: nr1d(kts:kte)
    REAL(KIND=r8), INTENT(INOUT) :: nc1d(kts:kte)
    REAL(KIND=r8), INTENT(INOUT) :: nwfa1d(kts:kte)
    REAL(KIND=r8), INTENT(INOUT) :: nifa1d(kts:kte)
    REAL(KIND=r8), INTENT(INOUT) :: t1d(kts:kte)

    REAL(KIND=r8), INTENT(IN   ) :: p1d(kts:kte)
    REAL(KIND=r8), INTENT(IN   ) :: w1d(kts:kte)
    REAL(KIND=r8), INTENT(IN   ) :: dzq(kts:kte)
    REAL(KIND=r8), INTENT(INOUT) :: pptrain
    REAL(KIND=r8), INTENT(INOUT) :: pptsnow
    REAL(KIND=r8), INTENT(INOUT) :: pptgraul
    REAL(KIND=r8), INTENT(INOUT) :: pptice
    REAL(KIND=r8), INTENT(IN   ) :: dt
    !#if ( WRF_CHEM == 1 )
    REAL(KIND=r8), INTENT(INOUT) :: rainprod(kts:kte)
    REAL(KIND=r8), INTENT(INOUT) :: evapprod(kts:kte)
    !#endif

    !..Local variables
    REAL(KIND=r8):: tten(kts:kte)
    REAL(KIND=r8):: qvten(kts:kte)
    REAL(KIND=r8):: qcten(kts:kte)
    REAL(KIND=r8):: qiten(kts:kte)
    REAL(KIND=r8):: qrten(kts:kte)
    REAL(KIND=r8):: qsten(kts:kte)
    REAL(KIND=r8):: qgten(kts:kte)
    REAL(KIND=r8):: niten(kts:kte)
    REAL(KIND=r8):: nrten(kts:kte)
    REAL(KIND=r8):: ncten(kts:kte)
    REAL(KIND=r8):: nwfaten(kts:kte)
    REAL(KIND=r8):: nifaten(kts:kte)

    REAL(KIND=r8)  :: prw_vcd(kts:kte)

    REAL(KIND=r8)  :: pnc_wcd(kts:kte)
    REAL(KIND=r8)  :: pnc_wau(kts:kte)
    REAL(KIND=r8)  :: pnc_rcw(kts:kte)
    REAL(KIND=r8)  :: pnc_scw(kts:kte)
    REAL(KIND=r8)  :: pnc_gcw(kts:kte)

    REAL(KIND=r8) :: pna_rca(kts:kte)
    REAL(KIND=r8) :: pna_sca(kts:kte)
    REAL(KIND=r8) :: pna_gca(kts:kte)
    REAL(KIND=r8) :: pnd_rcd(kts:kte)
    REAL(KIND=r8) :: pnd_scd(kts:kte)
    REAL(KIND=r8) :: pnd_gcd(kts:kte)

    REAL(KIND=r8) :: prr_wau(kts:kte)
    REAL(KIND=r8) :: prr_rcw(kts:kte)
    REAL(KIND=r8) :: prr_rcs(kts:kte)
    REAL(KIND=r8) :: prr_rcg(kts:kte)
    REAL(KIND=r8) :: prr_sml(kts:kte)
    REAL(KIND=r8) :: prr_gml(kts:kte)
    REAL(KIND=r8) :: prr_rci(kts:kte)
    REAL(KIND=r8) :: prv_rev(kts:kte)
    REAL(KIND=r8) :: pnr_wau(kts:kte)
    REAL(KIND=r8) :: pnr_rcs(kts:kte)
    REAL(KIND=r8) :: pnr_rcg(kts:kte)
    REAL(KIND=r8) :: pnr_rci(kts:kte)
    REAL(KIND=r8) :: pnr_sml(kts:kte)
    REAL(KIND=r8) :: pnr_gml(kts:kte)
    REAL(KIND=r8) :: pnr_rev(kts:kte)
    REAL(KIND=r8) :: pnr_rcr(kts:kte)
    REAL(KIND=r8) :: pnr_rfz(kts:kte)

    REAL(KIND=r8) :: pri_inu(kts:kte)
    REAL(KIND=r8) :: pni_inu(kts:kte)
    REAL(KIND=r8) :: pri_ihm(kts:kte)
    REAL(KIND=r8) :: pni_ihm(kts:kte)
    REAL(KIND=r8) :: pri_wfz(kts:kte)
    REAL(KIND=r8) :: pni_wfz(kts:kte)
    REAL(KIND=r8) :: pri_rfz(kts:kte)
    REAL(KIND=r8) :: pni_rfz(kts:kte)
    REAL(KIND=r8) :: pri_ide(kts:kte)
    REAL(KIND=r8) :: pni_ide(kts:kte)
    REAL(KIND=r8) :: pri_rci(kts:kte)
    REAL(KIND=r8) :: pni_rci(kts:kte)
    REAL(KIND=r8) :: pni_sci(kts:kte)
    REAL(KIND=r8) :: pni_iau(kts:kte)
    REAL(KIND=r8) :: pri_iha(kts:kte)
    REAL(KIND=r8) :: pni_iha(kts:kte)

    REAL(KIND=r8) :: prs_iau(kts:kte)
    REAL(KIND=r8) :: prs_sci(kts:kte)
    REAL(KIND=r8) :: prs_rcs(kts:kte)
    REAL(KIND=r8) :: prs_scw(kts:kte)
    REAL(KIND=r8) :: prs_sde(kts:kte)
    REAL(KIND=r8) :: prs_ihm(kts:kte)
    REAL(KIND=r8) :: prs_ide(kts:kte)

    REAL(KIND=r8) :: prg_scw(kts:kte)
    REAL(KIND=r8) :: prg_rfz(kts:kte)
    REAL(KIND=r8) :: prg_gde(kts:kte)
    REAL(KIND=r8) :: prg_gcw(kts:kte)
    REAL(KIND=r8) :: prg_rci(kts:kte)
    REAL(KIND=r8) :: prg_rcs(kts:kte)
    REAL(KIND=r8) :: prg_rcg(kts:kte)
    REAL(KIND=r8) :: prg_ihm(kts:kte)

    REAL(KIND=r8), PARAMETER:: zeroD0 = 0.0e0_r8

    REAL(KIND=r8):: temp(kts:kte)
    REAL(KIND=r8):: pres(kts:kte)
    REAL(KIND=r8):: qv  (kts:kte)
    REAL(KIND=r8):: rc  (kts:kte)
    REAL(KIND=r8):: ri  (kts:kte)
    REAL(KIND=r8):: rr  (kts:kte)
    REAL(KIND=r8):: rs  (kts:kte)
    REAL(KIND=r8):: rg  (kts:kte)
    REAL(KIND=r8):: ni  (kts:kte)
    REAL(KIND=r8):: nr  (kts:kte)
    REAL(KIND=r8):: nc  (kts:kte)
    REAL(KIND=r8):: nwfa(kts:kte)
    REAL(KIND=r8):: nifa(kts:kte)
    REAL(KIND=r8):: rho (kts:kte)
    REAL(KIND=r8):: rhof(kts:kte)
    REAL(KIND=r8):: rhof2 (kts:kte)
    REAL(KIND=r8):: qvs   (kts:kte)
    REAL(KIND=r8):: qvsi  (kts:kte)
    REAL(KIND=r8):: delQvs(kts:kte)
    REAL(KIND=r8):: satw  (kts:kte)
    REAL(KIND=r8):: sati  (kts:kte)
    REAL(KIND=r8):: ssatw (kts:kte)
    REAL(KIND=r8):: ssati (kts:kte)
    REAL(KIND=r8):: diffu (kts:kte)
    REAL(KIND=r8):: visco (kts:kte)
    REAL(KIND=r8):: vsc2  (kts:kte)
    REAL(KIND=r8):: tcond (kts:kte)
    REAL(KIND=r8):: lvap  (kts:kte)
    REAL(KIND=r8):: ocp   (kts:kte)
    REAL(KIND=r8):: lvt2  (kts:kte)

    REAL(KIND=r8) :: ilamr(kts:kte)
    REAL(KIND=r8) :: ilamg(kts:kte)
    REAL(KIND=r8) :: N0_r (kts:kte)
    REAL(KIND=r8) :: N0_g (kts:kte)
    REAL(KIND=r8):: mvd_r(kts:kte)
    REAL(KIND=r8):: mvd_c(kts:kte)
    REAL(KIND=r8):: smob (kts:kte)
    REAL(KIND=r8):: smo2 (kts:kte)
    REAL(KIND=r8):: smo1 (kts:kte)
    REAL(KIND=r8):: smo0 (kts:kte)
    REAL(KIND=r8):: smoc (kts:kte)
    REAL(KIND=r8):: smod (kts:kte)
    REAL(KIND=r8):: smoe (kts:kte)
    REAL(KIND=r8):: smof (kts:kte)

    REAL(KIND=r8):: sed_r(kts:kte)
    REAL(KIND=r8):: sed_s(kts:kte)
    REAL(KIND=r8):: sed_g(kts:kte)
    REAL(KIND=r8):: sed_i(kts:kte)
    REAL(KIND=r8):: sed_n(kts:kte)
    REAL(KIND=r8):: sed_c(kts:kte)

!    REAL(KIND=r8):: rgvm
    REAL(KIND=r8):: delta_tp
    REAL(KIND=r8):: orho
    REAL(KIND=r8):: lfus2
    REAL(KIND=r8) :: onstep(5)
    REAL(KIND=r8) :: N0_exp
    REAL(KIND=r8) :: N0_min
    REAL(KIND=r8) :: lam_exp
    REAL(KIND=r8) :: lamc
    REAL(KIND=r8) :: lamr
    REAL(KIND=r8) :: lamg
    REAL(KIND=r8) :: lami
    REAL(KIND=r8) :: ilami
    REAL(KIND=r8) :: ilamc
    REAL(KIND=r8):: xDc
    REAL(KIND=r8):: Dc_b
    REAL(KIND=r8):: Dc_g
    REAL(KIND=r8):: xDi
!    REAL(KIND=r8):: xDr
    REAL(KIND=r8):: xDs
    REAL(KIND=r8):: xDg
!    REAL(KIND=r8):: Ds_m
!    REAL(KIND=r8):: Dg_m
!    REAL(KIND=r8) :: Dr_star
    REAL(KIND=r8) ::  Dc_star
    REAL(KIND=r8):: zeta1
    REAL(KIND=r8):: zeta
    REAL(KIND=r8):: taud
    REAL(KIND=r8):: tau
!    REAL(KIND=r8):: stoke_r
!    REAL(KIND=r8):: stoke_s
    REAL(KIND=r8):: stoke_g
!    REAL(KIND=r8):: stoke_i
    REAL(KIND=r8):: vti
    REAL(KIND=r8):: vtr
    REAL(KIND=r8):: vts
    REAL(KIND=r8):: vtg
    REAL(KIND=r8):: vtc
    REAL(KIND=r8):: vtik (kts:kte+1)
    REAL(KIND=r8):: vtnik(kts:kte+1)
    REAL(KIND=r8):: vtrk (kts:kte+1)
    REAL(KIND=r8):: vtnrk(kts:kte+1)
    REAL(KIND=r8):: vtsk (kts:kte+1)
    REAL(KIND=r8):: vtgk (kts:kte+1)
    REAL(KIND=r8):: vtck (kts:kte+1)
    REAL(KIND=r8):: vtnck(kts:kte+1)

    REAL(KIND=r8):: vts_boost(kts:kte)
    REAL(KIND=r8):: Mrat
    REAL(KIND=r8):: ils1
    REAL(KIND=r8):: ils2
    REAL(KIND=r8):: t1_vts
    REAL(KIND=r8):: t2_vts
    REAL(KIND=r8):: t3_vts
    REAL(KIND=r8):: t4_vts
    REAL(KIND=r8):: C_snow
    REAL(KIND=r8):: a_
    REAL(KIND=r8):: b_
    REAL(KIND=r8):: loga_
!    REAL(KIND=r8):: A1
!    REAL(KIND=r8):: A2
    REAL(KIND=r8):: tf
    REAL(KIND=r8):: tempc
    REAL(KIND=r8):: tc0
!    REAL(KIND=r8):: r_mvd1
!    REAL(KIND=r8):: r_mvd2
!    REAL(KIND=r8):: xkrat
    REAL(KIND=r8):: xnc
    REAL(KIND=r8):: xri
    REAL(KIND=r8):: xni
    REAL(KIND=r8):: xmi
    REAL(KIND=r8):: oxmi
    REAL(KIND=r8):: xrc
    REAL(KIND=r8):: xrr
    REAL(KIND=r8):: xnr
    REAL(KIND=r8):: xsat
    REAL(KIND=r8):: rate_max
    REAL(KIND=r8):: sump
    REAL(KIND=r8):: ratio
    REAL(KIND=r8):: clap
    REAL(KIND=r8):: fcd
    REAL(KIND=r8):: dfcd
    REAL(KIND=r8):: otemp
    REAL(KIND=r8):: rvs
    REAL(KIND=r8):: rvs_p
    REAL(KIND=r8):: rvs_pp
    REAL(KIND=r8):: gamsc
    REAL(KIND=r8):: alphsc
    REAL(KIND=r8):: t1_evap
    REAL(KIND=r8):: t1_subl
    REAL(KIND=r8):: r_frac
    REAL(KIND=r8):: g_frac
    REAL(KIND=r8):: Ef_rw
    REAL(KIND=r8):: Ef_sw
    REAL(KIND=r8):: Ef_gw
    REAL(KIND=r8):: Ef_rr
    REAL(KIND=r8):: Ef_ra
    REAL(KIND=r8):: Ef_sa
    REAL(KIND=r8):: Ef_ga
    REAL(KIND=r8):: dtsave
    REAL(KIND=r8):: odts
    REAL(KIND=r8):: odt
    REAL(KIND=r8):: odzq
    REAL(KIND=r8):: hgt_agl
    REAL(KIND=r8):: xslw1
    REAL(KIND=r8):: ygra1
    REAL(KIND=r8):: zans1
    REAL(KIND=r8):: eva_factor
!    INTEGER :: i
    INTEGER :: k
!    INTEGER :: k2
    INTEGER :: n
    INTEGER :: nn
    INTEGER :: nstep
!    INTEGER :: k_0
!    INTEGER :: kbot
    INTEGER :: IT
    INTEGER :: iexfrq
    INTEGER :: ksed1(5)
    INTEGER :: nir
    INTEGER :: nis
    INTEGER :: nig
    INTEGER :: nii
    INTEGER :: nic
    INTEGER :: niin
    INTEGER :: idx_tc
    INTEGER :: idx_t
    INTEGER :: idx_s
    INTEGER :: idx_g1
    INTEGER :: idx_g
    INTEGER :: idx_r1
    INTEGER :: idx_r
    INTEGER :: idx_i1
    INTEGER :: idx_i
    INTEGER :: idx_c
    INTEGER :: idx
    INTEGER :: idx_d
    INTEGER :: idx_n
    INTEGER :: idx_in

!    LOGICAL :: melti
    LOGICAL :: no_micro
    LOGICAL :: L_qc(kts:kte)
    LOGICAL :: L_qi(kts:kte)
    LOGICAL :: L_qr(kts:kte)
    LOGICAL :: L_qs(kts:kte)
    LOGICAL :: L_qg(kts:kte)
    LOGICAL :: debug_flag
    CHARACTER*256:: mp_debug
    INTEGER:: nu_c

    !+---+

    debug_flag = .FALSE.
    !     if (ii.eq.901 .and.eq.379) debug_flag = .true.
    IF(debug_flag) THEN
       WRITE(mp_debug, *) 'DEBUG INFO, mp_thompson at (i,j) ', ii, ', '
       CALL wrf_debug('mp_thompson..', mp_debug)
    ENDIF

    no_micro = .TRUE.
    dtsave = dt
    odt = 1.0_r8/dt
    odts = 1.0_r8/dtsave
    iexfrq = 1

    !+---+-----------------------------------------------------------------+
    !.. Source/sink terms.  First 2 chars: "pr" represents source/sink of
    !.. mass while "pn" represents source/sink of number.  Next char is one
    !.. of "v" for water vapor, "r" for rain, "i" for cloud ice, "w" for
    !.. cloud water, "s" for snow, and "g" for graupel.  Next chars
    !.. represent processes: "de" for sublimation/deposition, "ev" for
    !.. evaporation, "fz" for freezing, "ml" for melting, "au" for
    !.. autoconversion, "nu" for ice nucleation, "hm" for Hallet/Mossop
    !.. secondary ice production, and "c" for collection followed by the
    !.. character for the species being collected.  ALL of these terms are
    !.. positive (except for deposition/sublimation terms which can switch
    !.. signs based on super/subsaturation) and are treated as negatives
    !.. where necessary in the tendency equations.
    !+---+-----------------------------------------------------------------+

    DO k = kts, kte
       tten(k) = 0.0_r8
       qvten(k) = 0.0_r8
       qcten(k) = 0.0_r8
       qiten(k) = 0.0_r8
       qrten(k) = 0.0_r8
       qsten(k) = 0.0_r8
       qgten(k) = 0.0_r8
       niten(k) = 0.0_r8
       nrten(k) = 0.0_r8
       ncten(k) = 0.0_r8
       nwfaten(k) = 0.0_r8
       nifaten(k) = 0.0_r8

       prw_vcd(k) = 0.0_r8

       pnc_wcd(k) = 0.0_r8
       pnc_wau(k) = 0.0_r8
       pnc_rcw(k) = 0.0_r8
       pnc_scw(k) = 0.0_r8
       pnc_gcw(k) = 0.0_r8

       prv_rev(k) = 0.0_r8
       prr_wau(k) = 0.0_r8
       prr_rcw(k) = 0.0_r8
       prr_rcs(k) = 0.0_r8
       prr_rcg(k) = 0.0_r8
       prr_sml(k) = 0.0_r8
       prr_gml(k) = 0.0_r8
       prr_rci(k) = 0.0_r8
       pnr_wau(k) = 0.0_r8
       pnr_rcs(k) = 0.0_r8
       pnr_rcg(k) = 0.0_r8
       pnr_rci(k) = 0.0_r8
       pnr_sml(k) = 0.0_r8
       pnr_gml(k) = 0.0_r8
       pnr_rev(k) = 0.0_r8
       pnr_rcr(k) = 0.0_r8
       pnr_rfz(k) = 0.0_r8

       pri_inu(k) = 0.0_r8
       pni_inu(k) = 0.0_r8
       pri_ihm(k) = 0.0_r8
       pni_ihm(k) = 0.0_r8
       pri_wfz(k) = 0.0_r8
       pni_wfz(k) = 0.0_r8
       pri_rfz(k) = 0.0_r8
       pni_rfz(k) = 0.0_r8
       pri_ide(k) = 0.0_r8
       pni_ide(k) = 0.0_r8
       pri_rci(k) = 0.0_r8
       pni_rci(k) = 0.0_r8
       pni_sci(k) = 0.0_r8
       pni_iau(k) = 0.0_r8
       pri_iha(k) = 0.0_r8
       pni_iha(k) = 0.0_r8

       prs_iau(k) = 0.0_r8
       prs_sci(k) = 0.0_r8
       prs_rcs(k) = 0.0_r8
       prs_scw(k) = 0.0_r8
       prs_sde(k) = 0.0_r8
       prs_ihm(k) = 0.0_r8
       prs_ide(k) = 0.0_r8

       prg_scw(k) = 0.0_r8
       prg_rfz(k) = 0.0_r8
       prg_gde(k) = 0.0_r8
       prg_gcw(k) = 0.0_r8
       prg_rci(k) = 0.0_r8
       prg_rcs(k) = 0.0_r8
       prg_rcg(k) = 0.0_r8
       prg_ihm(k) = 0.0_r8

       pna_rca(k) = 0.0_r8
       pna_sca(k) = 0.0_r8
       pna_gca(k) = 0.0_r8

       pnd_rcd(k) = 0.0_r8
       pnd_scd(k) = 0.0_r8
       pnd_gcd(k) = 0.0_r8
    ENDDO
    !#if ( WRF_CHEM == 1 )
    DO k = kts, kte
       rainprod(k) = 0.0_r8
       evapprod(k) = 0.0_r8
    ENDDO
    !#endif

    !+---+-----------------------------------------------------------------+
    !..Put column of data into local arrays.
    !+---+-----------------------------------------------------------------+
    DO k = kts, kte
       temp(k) = t1d(k)
       qv(k) = MAX(1.E-10_r8, qv1d(k))
       pres(k) = p1d(k)

       rho(k) = 0.622_r8*pres(k)/(R*temp(k)*(qv(k)+0.622_r8))
       nwfa(k) = MAX(11.1E6_r8, MIN(9999.E6_r8, nwfa1d(k)*rho(k)))
       nifa(k) = MAX(naIN1*0.01_r8, MIN(9999.E6_r8, nifa1d(k)*rho(k)))

       IF (qc1d(k) .GT. R1) THEN
          no_micro = .FALSE.
          rc(k) = qc1d(k)*rho(k)
          nc(k) = MAX(2.0_r8, nc1d(k)*rho(k))
          L_qc(k) = .TRUE.
          nu_c = MIN(15, NINT(1000.E6_r8/nc(k)) + 2)
          IF((nc(k)*am_r*ccg(2,nu_c)*ocg1(nu_c)/rc(k)) <0.0_r8)THEN
             PRINT*,k,pres(k),R,temp(k),qv(k),rc(k)

            ! PRINT*,nc(k),am_r,ccg(2,nu_c),ocg1(nu_c),rc(k)
             lamc=1.e-12_r8
          ELSE
             lamc = (nc(k)*am_r*ccg(2,nu_c)*ocg1(nu_c)/rc(k))**obmr
          END IF

          xDc = (bm_r + nu_c + 1.0_r8) / lamc
          IF (xDc.LT. D0c) THEN
             lamc = cce(2,nu_c)/D0c
          ELSEIF (xDc.GT. D0r*2.0_r8) THEN
             lamc = cce(2,nu_c)/(D0r*2.0_r8)
          ENDIF
          nc(k) = MIN( DBLE(Nt_c_max), ccg(1,nu_c)*ocg2(nu_c)*rc(k)   &
               / am_r*lamc**bm_r)
          IF (.NOT. is_aerosol_aware) nc(k) = Nt_c
       ELSE
          qc1d(k) = 0.0_r8
          nc1d(k) = 0.0_r8
          rc(k) = R1
          nc(k) = 2.0_r8
          L_qc(k) = .FALSE.
       ENDIF

       IF (qi1d(k) .GT. R1) THEN
          no_micro = .FALSE.
          ri(k) = qi1d(k)*rho(k)
          ni(k) = MAX(R2, ni1d(k)*rho(k))
          L_qi(k) = .TRUE.
          lami = (am_i*cig(2)*oig1*ni(k)/ri(k))**obmi
          ilami = 1.0_r8/lami
          xDi = (bm_i + mu_i + 1.0_r8) * ilami
          IF (xDi.LT. 20.E-6_r8) THEN
             lami = cie(2)/20.E-6_r8
             ni(k) = MIN(499.e3_r8, cig(1)*oig2*ri(k)/am_i*lami**bm_i)
          ELSEIF (xDi.GT. 300.E-6_r8) THEN
             lami = cie(2)/300.E-6_r8
             ni(k) = cig(1)*oig2*ri(k)/am_i*lami**bm_i
          ENDIF
       ELSE
          qi1d(k) = 0.0_r8
          ni1d(k) = 0.0_r8
          ri(k) = R1
          ni(k) = R2
          L_qi(k) = .FALSE.
       ENDIF

       IF (qr1d(k) .GT. R1) THEN
          no_micro = .FALSE.
          rr(k) = qr1d(k)*rho(k)
          nr(k) = MAX(R2, nr1d(k)*rho(k))
          L_qr(k) = .TRUE.
          lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
          mvd_r(k) = (3.0_r8 + mu_r + 0.672_r8) / lamr
          IF (mvd_r(k) .GT. 2.5E-3_r8) THEN
             mvd_r(k) = 2.5E-3_r8
             lamr = (3.0_r8 + mu_r + 0.672_r8) / mvd_r(k)
             nr(k) = crg(2)*org3*rr(k)*lamr**bm_r / am_r
          ELSEIF (mvd_r(k) .LT. D0r*0.75_r8) THEN
             mvd_r(k) = D0r*0.75_r8
             lamr = (3.0_r8 + mu_r + 0.672_r8) / mvd_r(k)
             nr(k) = crg(2)*org3*rr(k)*lamr**bm_r / am_r
          ENDIF
       ELSE
          qr1d(k) = 0.0_r8
          nr1d(k) = 0.0_r8
          rr(k) = R1
          nr(k) = R2
          L_qr(k) = .FALSE.
       ENDIF
       IF (qs1d(k) .GT. R1) THEN
          no_micro = .FALSE.
          rs(k) = qs1d(k)*rho(k)
          L_qs(k) = .TRUE.
       ELSE
          qs1d(k) = 0.0_r8
          rs(k) = R1
          L_qs(k) = .FALSE.
       ENDIF
       IF (qg1d(k) .GT. R1) THEN
          no_micro = .FALSE.
          rg(k) = qg1d(k)*rho(k)
          L_qg(k) = .TRUE.
       ELSE
          qg1d(k) = 0.0_r8
          rg(k) = R1
          L_qg(k) = .FALSE.
       ENDIF
    ENDDO

    !+---+-----------------------------------------------------------------+
    !     if (debug_flag) then
    !      write(mp_debug,*) 'DEBUG-VERBOSE at (i,j) ', ii, ', '
    !      CALL wrf_debug(550, mp_debug)
    !      do k = kts, kte
    !        write(mp_debug, '(a,i3,f8.2,1x,f7.2,1x, 11(1x,e13.6))')        &
    !    &              'VERBOSE: ', k, pres(k)*0.01, temp(k)-273.15, qv(k), rc(k), rr(k), ri(k), rs(k), rg(k), nc(k), nr(k), ni(k), nwfa(k), nifa(k)
    !        CALL wrf_debug(550, mp_debug)
    !      enddo
    !     endif
    !+---+-----------------------------------------------------------------+

    !+---+-----------------------------------------------------------------+
    !..Derive various thermodynamic variables frequently used.
    !.. Saturation vapor pressure (mixing ratio) over liquid/ice comes from
    !.. Flatau et al. 1992; enthalpy (latent heat) of vaporization from
    !.. Bohren & Albrecht 1998; others from Pruppacher & Klett 1978.
    !+---+-----------------------------------------------------------------+
    DO k = kts, kte
       tempc = temp(k) - 273.15_r8
       rhof(k) = SQRT(RHO_NOT/rho(k))
       rhof2(k) = SQRT(rhof(k))
       qvs(k) = rslf(pres(k), temp(k))
       delQvs(k) = MAX(0.0_r8, rslf(pres(k), 273.15_r8)-qv(k))
       IF (tempc .LE. 0.0_r8) THEN
          qvsi(k) = rsif(pres(k), temp(k))
       ELSE
          qvsi(k) = qvs(k)
       ENDIF
       satw(k) = qv(k)/qvs(k)
       sati(k) = qv(k)/qvsi(k)
       ssatw(k) = satw(k) - 1.0_r8
       ssati(k) = sati(k) - 1.0_r8
       IF (ABS(ssatw(k)).LT. eps) ssatw(k) = 0.0_r8
       IF (ABS(ssati(k)).LT. eps) ssati(k) = 0.0_r8
       IF (no_micro .AND. ssati(k).GT. 0.0_r8) no_micro = .FALSE.
       diffu(k) = 2.11E-5_r8*(temp(k)/273.15_r8)**1.94_r8 * (101325.0_r8/pres(k))
       IF (tempc .GE. 0.0_r8) THEN
          visco(k) = (1.718_r8+0.0049_r8*tempc)*1.0E-5_r8
       ELSE
          visco(k) = (1.718_r8+0.0049_r8*tempc-1.2E-5_r8*tempc*tempc)*1.0E-5_r8
       ENDIF
       ocp(k) = 1.0_r8/(Cp*(1.0_r8+0.887_r8*qv(k)))
       vsc2(k) = SQRT(rho(k)/visco(k))
       lvap(k) = lvap0 + (2106.0_r8 - 4218.0_r8)*tempc
       tcond(k) = (5.69_r8 + 0.0168_r8*tempc)*1.0E-5_r8 * 418.936_r8
    ENDDO

    !+---+-----------------------------------------------------------------+
    !..If no existing hydrometeor species and no chance to initiate ice or
    !.. condense cloud water, just exit quickly!
    !+---+-----------------------------------------------------------------+

    IF (no_micro) RETURN

    !+---+-----------------------------------------------------------------+
    !..Calculate y-intercept, slope, and useful moments for snow.
    !+---+-----------------------------------------------------------------+
    IF (.NOT. iiwarm) THEN
       DO k = kts, kte
          IF (.NOT. L_qs(k)) CYCLE
          tc0 = MIN(-0.1_r8, temp(k)-273.15_r8)
          smob(k) = rs(k)*oams

          !..All other moments based on reference, 2nd moment.  If bm_s.ne.2,
          !.. then we must compute actual 2nd moment and use as reference.
          IF (bm_s.GT.(2.0_r8-1.e-3_r8) .AND. bm_s.LT.(2.0_r8+1.e-3_r8)) THEN
             smo2(k) = smob(k)
          ELSE
             loga_ = sa(1) + sa(2)*tc0 + sa(3)*bm_s &
                  + sa(4)*tc0*bm_s + sa(5)*tc0*tc0 &
                  + sa(6)*bm_s*bm_s + sa(7)*tc0*tc0*bm_s &
                  + sa(8)*tc0*bm_s*bm_s + sa(9)*tc0*tc0*tc0 &
                  + sa(10)*bm_s*bm_s*bm_s
             a_ = 10.0_r8**loga_
             b_ = sb(1) + sb(2)*tc0 + sb(3)*bm_s &
                  + sb(4)*tc0*bm_s + sb(5)*tc0*tc0 &
                  + sb(6)*bm_s*bm_s + sb(7)*tc0*tc0*bm_s &
                  + sb(8)*tc0*bm_s*bm_s + sb(9)*tc0*tc0*tc0 &
                  + sb(10)*bm_s*bm_s*bm_s
             smo2(k) = (smob(k)/a_)**(1.0_r8/b_)
          ENDIF

          !..Calculate 0th moment.  Represents snow number concentration.
          loga_ = sa(1) + sa(2)*tc0 + sa(5)*tc0*tc0 + sa(9)*tc0*tc0*tc0
          a_ = 10.0_r8**loga_
          b_ = sb(1) + sb(2)*tc0 + sb(5)*tc0*tc0 + sb(9)*tc0*tc0*tc0
          smo0(k) = a_ * smo2(k)**b_

          !..Calculate 1st moment.  Useful for depositional growth and melting.
          loga_ = sa(1) + sa(2)*tc0 + sa(3) &
               + sa(4)*tc0 + sa(5)*tc0*tc0 &
               + sa(6) + sa(7)*tc0*tc0 &
               + sa(8)*tc0 + sa(9)*tc0*tc0*tc0 &
               + sa(10)
          a_ = 10.0_r8**loga_
          b_ = sb(1)+ sb(2)*tc0 + sb(3) + sb(4)*tc0 &
               + sb(5)*tc0*tc0 + sb(6) &
               + sb(7)*tc0*tc0 + sb(8)*tc0 &
               + sb(9)*tc0*tc0*tc0 + sb(10)
          smo1(k) = a_ * smo2(k)**b_

          !..Calculate bm_s+1 (th) moment.  Useful for diameter calcs.
          loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(1) &
               + sa(4)*tc0*cse(1) + sa(5)*tc0*tc0 &
               + sa(6)*cse(1)*cse(1) + sa(7)*tc0*tc0*cse(1) &
               + sa(8)*tc0*cse(1)*cse(1) + sa(9)*tc0*tc0*tc0 &
               + sa(10)*cse(1)*cse(1)*cse(1)
          a_ = 10.0_r8**loga_
          b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(1) + sb(4)*tc0*cse(1) &
               + sb(5)*tc0*tc0 + sb(6)*cse(1)*cse(1) &
               + sb(7)*tc0*tc0*cse(1) + sb(8)*tc0*cse(1)*cse(1) &
               + sb(9)*tc0*tc0*tc0 + sb(10)*cse(1)*cse(1)*cse(1)
          smoc(k) = a_ * smo2(k)**b_

          !..Calculate bv_s+2 (th) moment.  Useful for riming.
          loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(13) &
               + sa(4)*tc0*cse(13) + sa(5)*tc0*tc0 &
               + sa(6)*cse(13)*cse(13) + sa(7)*tc0*tc0*cse(13) &
               + sa(8)*tc0*cse(13)*cse(13) + sa(9)*tc0*tc0*tc0 &
               + sa(10)*cse(13)*cse(13)*cse(13)
          a_ = 10.0_r8**loga_
          b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(13) + sb(4)*tc0*cse(13) &
               + sb(5)*tc0*tc0 + sb(6)*cse(13)*cse(13) &
               + sb(7)*tc0*tc0*cse(13) + sb(8)*tc0*cse(13)*cse(13) &
               + sb(9)*tc0*tc0*tc0 + sb(10)*cse(13)*cse(13)*cse(13)
          smoe(k) = a_ * smo2(k)**b_

          !..Calculate 1+(bv_s+1)/2 (th) moment.  Useful for depositional growth.
          loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(16) &
               + sa(4)*tc0*cse(16) + sa(5)*tc0*tc0 &
               + sa(6)*cse(16)*cse(16) + sa(7)*tc0*tc0*cse(16) &
               + sa(8)*tc0*cse(16)*cse(16) + sa(9)*tc0*tc0*tc0 &
               + sa(10)*cse(16)*cse(16)*cse(16)
          a_ = 10.0_r8**loga_
          b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(16) + sb(4)*tc0*cse(16) &
               + sb(5)*tc0*tc0 + sb(6)*cse(16)*cse(16) &
               + sb(7)*tc0*tc0*cse(16) + sb(8)*tc0*cse(16)*cse(16) &
               + sb(9)*tc0*tc0*tc0 + sb(10)*cse(16)*cse(16)*cse(16)
          smof(k) = a_ * smo2(k)**b_

       ENDDO

       !+---+-----------------------------------------------------------------+
       !..Calculate y-intercept, slope values for graupel.
       !+---+-----------------------------------------------------------------+
       N0_min = gonv_max
       DO k = kte, kts, -1
          IF (temp(k).LT.270.65_r8 .AND. L_qr(k) .AND. mvd_r(k).GT.100.E-6_r8) THEN
             xslw1 = 4.01_r8 + log10(mvd_r(k))
          ELSE
             xslw1 = 0.01_r8
          ENDIF
          ygra1 = 4.31_r8 + log10(MAX(5.E-5_r8, rg(k)))
          zans1 = 3.1_r8 + (100.0_r8/(300.0_r8*xslw1*ygra1/(10.0_r8/xslw1+1.0_r8+0.25_r8*ygra1)+30.0_r8+10.0_r8*ygra1))
          N0_exp = 10.0_r8**(zans1)
          N0_exp = MAX(DBLE(gonv_min), MIN(N0_exp, DBLE(gonv_max)))
          N0_min = MIN(N0_exp, N0_min)
          N0_exp = N0_min
          lam_exp = (N0_exp*am_g*cgg(1)/rg(k))**oge1
          lamg = lam_exp * (cgg(3)*ogg2*ogg1)**obmg
          ilamg(k) = 1.0_r8/lamg
          N0_g(k) = N0_exp/(cgg(2)*lam_exp) * lamg**cge(2)
       ENDDO

    ENDIF

    !+---+-----------------------------------------------------------------+
    !..Calculate y-intercept, slope values for rain.
    !+---+-----------------------------------------------------------------+
    DO k = kte, kts, -1
       lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
       ilamr(k) = 1.0_r8/lamr
       mvd_r(k) = (3.0_r8 + mu_r + 0.672_r8) / lamr
       N0_r(k) = nr(k)*org2*lamr**cre(2)
    ENDDO

    !+---+-----------------------------------------------------------------+
    !..Compute warm-rain process terms (except evap done later).
    !+---+-----------------------------------------------------------------+

    DO k = kts, kte

       !..Rain self-collection follows Seifert, 1994 and drop break-up
       !.. follows Verlinde and Cotton, 1993.                                        RAIN2M
       IF (L_qr(k) .AND. mvd_r(k).GT. D0r) THEN
          !-GT      Ef_rr = 1.0
          !-GT      if (mvd_r(k) .gt. 1500.0E-6) then
          Ef_rr = 2.0_r8 - EXP(2300.0_r8*(mvd_r(k)-1600.0E-6_r8))
          !-GT      endif
          pnr_rcr(k) = Ef_rr * 0.5_r8*nr(k)*rr(k)
       ENDIF

       mvd_c(k) = D0c
       IF (L_qc(k)) THEN
          nu_c = MIN(15, NINT(1000.E6_r8/nc(k)) + 2)
          xDc = MAX(D0c*1.E6_r8, ((rc(k)/(am_r*nc(k)))**obmr) * 1.E6_r8)
          lamc = (nc(k)*am_r* ccg(2,nu_c) * ocg1(nu_c) / rc(k))**obmr
          mvd_c(k) = (3.0_r8+nu_c+0.672_r8) / lamc
       ENDIF

       !..Autoconversion follows Berry & Reinhardt (1974) with characteristic
       !.. diameters correctly computed from gamma distrib of cloud droplets.
       IF (rc(k).GT. 0.01e-3_r8) THEN
          Dc_g = ((ccg(3,nu_c)*ocg2(nu_c))**obmr / lamc) * 1.E6_r8
          Dc_b = (xDc*xDc*xDc*Dc_g*Dc_g*Dc_g - xDc*xDc*xDc*xDc*xDc*xDc) &
               **(1.0_r8/6.0_r8)
          zeta1 = 0.5_r8*((6.25E-6_r8*xDc*Dc_b*Dc_b*Dc_b - 0.4_r8) &
               + ABS(6.25E-6_r8*xDc*Dc_b*Dc_b*Dc_b - 0.4_r8))
          zeta = 0.027_r8*rc(k)*zeta1
          taud = 0.5_r8*((0.5_r8*Dc_b - 7.5_r8) + ABS(0.5_r8*Dc_b - 7.5_r8)) + R1
          tau  = 3.72_r8/(rc(k)*taud)
          prr_wau(k) = zeta/tau
          prr_wau(k) = MIN(DBLE(rc(k)*odts), prr_wau(k))
          pnr_wau(k) = prr_wau(k) / (am_r*nu_c*D0r*D0r*D0r)              ! RAIN2M
          pnc_wau(k) = MIN(DBLE(nc(k)*odts), prr_wau(k)                 &
               / (am_r*mvd_c(k)*mvd_c(k)*mvd_c(k)))                   ! Qc2M
       ENDIF

       !..Rain collecting cloud water.  In CE, assume Dc<<Dr and vtc=~0.
       IF (L_qr(k) .AND. mvd_r(k).GT. D0r .AND. mvd_c(k).GT. D0c) THEN
          lamr = 1.0_r8/ilamr(k)
          idx = 1 + INT(nbr*DLOG(mvd_r(k)/Dr(1))/DLOG(Dr(nbr)/Dr(1)))
          idx = MIN(idx, nbr)
          Ef_rw = t_Efrw(idx, INT(mvd_c(k)*1.E6_r8))
          prr_rcw(k) = rhof(k)*t1_qr_qc*Ef_rw*rc(k)*N0_r(k) &
               *((lamr+fv_r)**(-cre(9)))
          prr_rcw(k) = MIN(DBLE(rc(k)*odts), prr_rcw(k))
          pnc_rcw(k) = rhof(k)*t1_qr_qc*Ef_rw*nc(k)*N0_r(k)             &
               *((lamr+fv_r)**(-cre(9)))                          ! Qc2M
          pnc_rcw(k) = MIN(DBLE(nc(k)*odts), pnc_rcw(k))
       ENDIF

       !..Rain collecting aerosols, wet scavenging.
       IF (L_qr(k) .AND. mvd_r(k).GT. D0r) THEN
          Ef_ra = Eff_aero(mvd_r(k),0.04E-6_r8,visco(k),rho(k),temp(k),'r')
          lamr = 1.0_r8/ilamr(k)
          pna_rca(k) = rhof(k)*t1_qr_qc*Ef_ra*nwfa(k)*N0_r(k)           &
               *((lamr+fv_r)**(-cre(9)))
          pna_rca(k) = MIN(DBLE(nwfa(k)*odts), pna_rca(k))

          Ef_ra = Eff_aero(mvd_r(k),0.8E-6_r8,visco(k),rho(k),temp(k),'r')
          pnd_rcd(k) = rhof(k)*t1_qr_qc*Ef_ra*nifa(k)*N0_r(k)           &
               *((lamr+fv_r)**(-cre(9)))
          pnd_rcd(k) = MIN(DBLE(nifa(k)*odts), pnd_rcd(k))
       ENDIF

    ENDDO

    !+---+-----------------------------------------------------------------+
    !..Compute all frozen hydrometeor species' process terms.
    !+---+-----------------------------------------------------------------+
    IF (.NOT. iiwarm) THEN
       DO k = kts, kte
          vts_boost(k) = 1.5_r8

          !..Temperature lookup table indexes.
          tempc = temp(k) - 273.15_r8
          idx_tc = MAX(1, MIN(NINT(-tempc), 45) )
          idx_t = INT( (tempc-2.5_r8)/5.0_r8 ) - 1
          idx_t = MAX(1, -idx_t)
          idx_t = MIN(idx_t, ntb_t)
          IT = MAX(1, MIN(NINT(-tempc), 31) )

          !..Cloud water lookup table index.
          IF (rc(k).GT. r_c(1)) THEN
             nic = NINT(log10(rc(k)))
             DO nn = nic-1, nic+1
                n = nn
                IF ( (rc(k)/10.0_r8**nn).GE.1.0_r8 .AND. &
                     (rc(k)/10.0_r8**nn).LT.10.0_r8) GOTO 141
             ENDDO
141          CONTINUE
             idx_c = INT(rc(k)/10.0_r8**n) + 10*(n-nic2) - (n-nic2)
             idx_c = MAX(1, MIN(idx_c, ntb_c))
          ELSE
             idx_c = 1
          ENDIF

          !..Cloud droplet number lookup table index.
          idx_n = NINT(1.0_r8 + FLOAT(nbc) * DLOG(nc(k)/t_Nc(1)) / nic1)
          idx_n = MAX(1, MIN(idx_n, nbc))

          !..Cloud ice lookup table indexes.
          IF (ri(k).GT. r_i(1)) THEN
             nii = NINT(log10(ri(k)))
             DO nn = nii-1, nii+1
                n = nn
                IF ( (ri(k)/10.0_r8**nn).GE.1.0_r8 .AND. &
                     (ri(k)/10.0_r8**nn).LT.10.0_r8) GOTO 142
             ENDDO
142          CONTINUE
             idx_i = INT(ri(k)/10.0_r8**n) + 10*(n-nii2) - (n-nii2)
             idx_i = MAX(1, MIN(idx_i, ntb_i))
          ELSE
             idx_i = 1
          ENDIF

          IF (ni(k).GT. Nt_i(1)) THEN
             nii = NINT(log10(ni(k)))
             DO nn = nii-1, nii+1
                n = nn
                IF ( (ni(k)/10.0_r8**nn).GE.1.0_r8 .AND. &
                     (ni(k)/10.0_r8**nn).LT.10.0_r8) GOTO 143
             ENDDO
143          CONTINUE
             idx_i1 = INT(ni(k)/10.0_r8**n) + 10*(n-nii3) - (n-nii3)
             idx_i1 = MAX(1, MIN(idx_i1, ntb_i1))
          ELSE
             idx_i1 = 1
          ENDIF

          !..Rain lookup table indexes.
          IF (rr(k).GT. r_r(1)) THEN
             nir = NINT(log10(rr(k)))
             DO nn = nir-1, nir+1
                n = nn
                IF ( (rr(k)/10.0_r8**nn).GE.1.0_r8 .AND. &
                     (rr(k)/10.0_r8**nn).LT.10.0_r8) GOTO 144
             ENDDO
144          CONTINUE
             idx_r = INT(rr(k)/10.0_r8**n) + 10*(n-nir2) - (n-nir2)
             idx_r = MAX(1, MIN(idx_r, ntb_r))

             lamr = 1.0_r8/ilamr(k)
             lam_exp = lamr * (crg(3)*org2*org1)**bm_r
             N0_exp = org1*rr(k)/am_r * lam_exp**cre(1)
             nir = NINT(DLOG10(N0_exp))
             DO nn = nir-1, nir+1
                n = nn
                IF ( (N0_exp/10.0_r8**nn).GE.1.0_r8 .AND. &
                     (N0_exp/10.0_r8**nn).LT.10.0_r8) GOTO 145
             ENDDO
145          CONTINUE
             idx_r1 = INT(N0_exp/10.0_r8**n) + 10*(n-nir3) - (n-nir3)
             idx_r1 = MAX(1, MIN(idx_r1, ntb_r1))
          ELSE
             idx_r = 1
             idx_r1 = ntb_r1
          ENDIF

          !..Snow lookup table index.
          IF (rs(k).GT. r_s(1)) THEN
             nis = NINT(log10(rs(k)))
             DO nn = nis-1, nis+1
                n = nn
                IF ( (rs(k)/10.0_r8**nn).GE.1.0_r8 .AND. &
                     (rs(k)/10.0_r8**nn).LT.10.0_r8) GOTO 146
             ENDDO
146          CONTINUE
             idx_s = INT(rs(k)/10.0_r8**n) + 10*(n-nis2) - (n-nis2)
             idx_s = MAX(1, MIN(idx_s, ntb_s))
          ELSE
             idx_s = 1
          ENDIF

          !..Graupel lookup table index.
          IF (rg(k).GT. r_g(1)) THEN
             nig = NINT(log10(rg(k)))
             DO nn = nig-1, nig+1
                n = nn
                IF ( (rg(k)/10.0_r8**nn).GE.1.0_r8 .AND. &
                     (rg(k)/10.0_r8**nn).LT.10.0_r8) GOTO 147
             ENDDO
147          CONTINUE
             idx_g = INT(rg(k)/10.0_r8**n) + 10*(n-nig2) - (n-nig2)
             idx_g = MAX(1, MIN(idx_g, ntb_g))

             lamg = 1.0_r8/ilamg(k)
             lam_exp = lamg * (cgg(3)*ogg2*ogg1)**bm_g
             N0_exp = ogg1*rg(k)/am_g * lam_exp**cge(1)
             nig = NINT(DLOG10(N0_exp))
             DO nn = nig-1, nig+1
                n = nn
                IF ( (N0_exp/10.0_r8**nn).GE.1.0_r8 .AND. &
                     (N0_exp/10.0_r8**nn).LT.10.0_r8) GOTO 148
             ENDDO
148          CONTINUE
             idx_g1 = INT(N0_exp/10.0_r8**n) + 10*(n-nig3) - (n-nig3)
             idx_g1 = MAX(1, MIN(idx_g1, ntb_g1))
          ELSE
             idx_g = 1
             idx_g1 = ntb_g1
          ENDIF

          !..Deposition/sublimation prefactor (from Srivastava & Coen 1992).
          otemp = 1.0_r8/temp(k)
          rvs = rho(k)*qvsi(k)
          rvs_p = rvs*otemp*(lsub*otemp*oRv - 1.0_r8)
          rvs_pp = rvs * ( otemp*(lsub*otemp*oRv - 1.0_r8) &
               *otemp*(lsub*otemp*oRv - 1.0_r8) &
               + (-2.0_r8*lsub*otemp*otemp*otemp*oRv) &
               + otemp*otemp)
          gamsc = lsub*diffu(k)/tcond(k) * rvs_p
          alphsc = 0.5_r8*(gamsc/(1.0_r8+gamsc))*(gamsc/(1.0_r8+gamsc)) &
               * rvs_pp/rvs_p * rvs/rvs_p
          alphsc = MAX(1.E-9_r8, alphsc)
          xsat = ssati(k)
          IF (ABS(xsat).LT. 1.E-9_r8) xsat=0.0_r8
          t1_subl = 4.0_r8*PI*( 1.0_r8 - alphsc*xsat &
               + 2.0_r8*alphsc*alphsc*xsat*xsat &
               - 5.0_r8*alphsc*alphsc*alphsc*xsat*xsat*xsat ) &
               / (1.0_r8+gamsc)

          !..Snow collecting cloud water.  In CE, assume Dc<<Ds and vtc=~0.
          IF (L_qc(k) .AND. mvd_c(k).GT. D0c) THEN
             xDs = 0.0_r8
             IF (L_qs(k)) xDs = smoc(k) / smob(k)
             IF (xDs .GT. D0s) THEN
                idx = 1 + INT(nbs*DLOG(xDs/Ds(1))/DLOG(Ds(nbs)/Ds(1)))
                idx = MIN(idx, nbs)
                Ef_sw = t_Efsw(idx, INT(mvd_c(k)*1.E6_r8))
                prs_scw(k) = rhof(k)*t1_qs_qc*Ef_sw*rc(k)*smoe(k)
                pnc_scw(k) = rhof(k)*t1_qs_qc*Ef_sw*nc(k)*smoe(k)                ! Qc2M
                pnc_scw(k) = MIN(DBLE(nc(k)*odts), pnc_scw(k))
             ENDIF

             !..Graupel collecting cloud water.  In CE, assume Dc<<Dg and vtc=~0.
             IF (rg(k).GE. r_g(1) .AND. mvd_c(k).GT. D0c) THEN
                xDg = (bm_g + mu_g + 1.0_r8) * ilamg(k)
                vtg = rhof(k)*av_g*cgg(6)*ogg3 * ilamg(k)**bv_g
                stoke_g = mvd_c(k)*mvd_c(k)*vtg*rho_w/(9.0_r8*visco(k)*xDg)
                IF (xDg.GT. D0g) THEN
                   IF (stoke_g.GE.0.4_r8 .AND. stoke_g.LE.10.0_r8) THEN
                      Ef_gw = 0.55_r8*log10(2.51_r8*stoke_g)
                   ELSEIF (stoke_g.LT.0.4_r8) THEN
                      Ef_gw = 0.0_r8
                   ELSEIF (stoke_g.GT.10) THEN
                      Ef_gw = 0.77_r8
                   ENDIF
                   prg_gcw(k) = rhof(k)*t1_qg_qc*Ef_gw*rc(k)*N0_g(k) &
                        *ilamg(k)**cge(9)
                   pnc_gcw(k) = rhof(k)*t1_qg_qc*Ef_gw*nc(k)*N0_g(k)           &
                        *ilamg(k)**cge(9)                                 ! Qc2M
                   pnc_gcw(k) = MIN(DBLE(nc(k)*odts), pnc_gcw(k))
                ENDIF
             ENDIF
          ENDIF

          !..Snow and graupel collecting aerosols, wet scavenging.
          IF (rs(k) .GT. r_s(1)) THEN
             xDs = smoc(k) / smob(k)
             Ef_sa = Eff_aero(xDs,0.04E-6_r8,visco(k),rho(k),temp(k),'s')
             pna_sca(k) = rhof(k)*t1_qs_qc*Ef_sa*nwfa(k)*smoe(k)
             pna_sca(k) = MIN(DBLE(nwfa(k)*odts), pna_sca(k))

             Ef_sa = Eff_aero(xDs,0.8E-6_r8,visco(k),rho(k),temp(k),'s')
             pnd_scd(k) = rhof(k)*t1_qs_qc*Ef_sa*nifa(k)*smoe(k)
             pnd_scd(k) = MIN(DBLE(nifa(k)*odts), pnd_scd(k))
          ENDIF
          IF (rg(k) .GT. r_g(1)) THEN
             xDg = (bm_g + mu_g + 1.0_r8) * ilamg(k)
             Ef_ga = Eff_aero(xDg,0.04E-6_r8,visco(k),rho(k),temp(k),'g')
             pna_gca(k) = rhof(k)*t1_qg_qc*Ef_ga*nwfa(k)*N0_g(k)           &
                  *ilamg(k)**cge(9)
             pna_gca(k) = MIN(DBLE(nwfa(k)*odts), pna_gca(k))

             Ef_ga = Eff_aero(xDg,0.8E-6_r8,visco(k),rho(k),temp(k),'g')
             pnd_gcd(k) = rhof(k)*t1_qg_qc*Ef_ga*nifa(k)*N0_g(k)           &
                  *ilamg(k)**cge(9)
             pnd_gcd(k) = MIN(DBLE(nifa(k)*odts), pnd_gcd(k))
          ENDIF

          !..Rain collecting snow.  Cannot assume Wisner (1972) approximation
          !.. or Mizuno (1990) approach so we solve the CE explicitly and store
          !.. results in lookup table.
          IF (rr(k).GE. r_r(1)) THEN
             IF (rs(k).GE. r_s(1)) THEN
                IF (temp(k).LT.T_0) THEN
                   prr_rcs(k) = -(tmr_racs2(idx_s,idx_t,idx_r1,idx_r) &
                        + tcr_sacr2(idx_s,idx_t,idx_r1,idx_r) &
                        + tmr_racs1(idx_s,idx_t,idx_r1,idx_r) &
                        + tcr_sacr1(idx_s,idx_t,idx_r1,idx_r))
                   prs_rcs(k) = tmr_racs2(idx_s,idx_t,idx_r1,idx_r) &
                        + tcr_sacr2(idx_s,idx_t,idx_r1,idx_r) &
                        - tcs_racs1(idx_s,idx_t,idx_r1,idx_r) &
                        - tms_sacr1(idx_s,idx_t,idx_r1,idx_r)
                   prg_rcs(k) = tmr_racs1(idx_s,idx_t,idx_r1,idx_r) &
                        + tcr_sacr1(idx_s,idx_t,idx_r1,idx_r) &
                        + tcs_racs1(idx_s,idx_t,idx_r1,idx_r) &
                        + tms_sacr1(idx_s,idx_t,idx_r1,idx_r)
                   prr_rcs(k) = MAX(DBLE(-rr(k)*odts), prr_rcs(k))
                   prs_rcs(k) = MAX(DBLE(-rs(k)*odts), prs_rcs(k))
                   prg_rcs(k) = MIN(DBLE((rr(k)+rs(k))*odts), prg_rcs(k))
                   pnr_rcs(k) = tnr_racs1(idx_s,idx_t,idx_r1,idx_r)            &   ! RAIN2M
                        + tnr_racs2(idx_s,idx_t,idx_r1,idx_r)          &
                        + tnr_sacr1(idx_s,idx_t,idx_r1,idx_r)          &
                        + tnr_sacr2(idx_s,idx_t,idx_r1,idx_r)
                ELSE
                   prs_rcs(k) = -tcs_racs1(idx_s,idx_t,idx_r1,idx_r)           &
                        - tms_sacr1(idx_s,idx_t,idx_r1,idx_r)          &
                        + tmr_racs2(idx_s,idx_t,idx_r1,idx_r)          &
                        + tcr_sacr2(idx_s,idx_t,idx_r1,idx_r)
                   prs_rcs(k) = MAX(DBLE(-rs(k)*odts), prs_rcs(k))
                   prr_rcs(k) = -prs_rcs(k)
                   pnr_rcs(k) = tnr_racs2(idx_s,idx_t,idx_r1,idx_r)            &   ! RAIN2M
                        + tnr_sacr2(idx_s,idx_t,idx_r1,idx_r)
                ENDIF
                pnr_rcs(k) = MIN(DBLE(nr(k)*odts), pnr_rcs(k))
             ENDIF

             !..Rain collecting graupel.  Cannot assume Wisner (1972) approximation
             !.. or Mizuno (1990) approach so we solve the CE explicitly and store
             !.. results in lookup table.
             IF (rg(k).GE. r_g(1)) THEN
                IF (temp(k).LT.T_0) THEN
                   prg_rcg(k) = tmr_racg(idx_g1,idx_g,idx_r1,idx_r) &
                        + tcr_gacr(idx_g1,idx_g,idx_r1,idx_r)
                   prg_rcg(k) = MIN(DBLE(rr(k)*odts), prg_rcg(k))
                   prr_rcg(k) = -prg_rcg(k)
                   pnr_rcg(k) = tnr_racg(idx_g1,idx_g,idx_r1,idx_r)            &   ! RAIN2M
                        + tnr_gacr(idx_g1,idx_g,idx_r1,idx_r)
                   pnr_rcg(k) = MIN(DBLE(nr(k)*odts), pnr_rcg(k))
                ELSE
                   prr_rcg(k) = tcg_racg(idx_g1,idx_g,idx_r1,idx_r)
                   prr_rcg(k) = MIN(DBLE(rg(k)*odts), prr_rcg(k))
                   prg_rcg(k) = -prr_rcg(k)
                   !..Put in explicit drop break-up due to collisions.
                   pnr_rcg(k) = -5.0_r8*tnr_gacr(idx_g1,idx_g,idx_r1,idx_r)         ! RAIN2M
                ENDIF
             ENDIF
          ENDIF

          !+---+-----------------------------------------------------------------+
          !..Next IF block handles only those processes below 0C.
          !+---+-----------------------------------------------------------------+

          IF (temp(k).LT.T_0) THEN

             vts_boost(k) = 1.0_r8
             rate_max = (qv(k)-qvsi(k))*rho(k)*odts*0.999_r8

             !+---+---------------- BEGIN NEW ICE NUCLEATION -----------------------+
             !..Freezing of supercooled water (rain or cloud) is influenced by dust
             !.. but still using Bigg 1953 with a temperature adjustment of a few
             !.. degrees depending on dust concentration.  A default value by way
             !.. of idx_IN is 1.0 per Liter of air is used when dustyIce flag is
             !.. false.  Next, a combination of deposition/condensation freezing
             !.. using DeMott et al (2010) dust nucleation when water saturated or
             !.. Phillips et al (2008) when below water saturation; else, without
             !.. dustyIce flag, use the previous Cooper (1986) temperature-dependent
             !.. value.  Lastly, allow homogeneous freezing of deliquesced aerosols
             !.. following Koop et al. (2001, Nature).
             !.. Implemented by T. Eidhammer and G. Thompson 2012Dec18
             !+---+-----------------------------------------------------------------+

             IF (dustyIce .AND. is_aerosol_aware) THEN
                xni = iceDeMott(tempc,qvs(k),qvs(k),qvsi(k),rho(k),nifa(k))
             ELSE
                xni = 1.0_r8 *1000.0_r8                                               ! Default is 1.0_r8 per Liter
             ENDIF

             !..Ice nuclei lookup table index.
             IF (xni.GT. Nt_IN(1)) THEN
                niin = NINT(log10(xni))
                DO nn = niin-1, niin+1
                   n = nn
                   IF ( (xni/10.0_r8**nn).GE.1.0_r8 .AND. &
                        (xni/10.0_r8**nn).LT.10.0_r8) GOTO 149
                ENDDO
149             CONTINUE
                idx_IN = INT(xni/10.0_r8**n) + 10*(n-niin2) - (n-niin2)
                idx_IN = MAX(1, MIN(idx_IN, ntb_IN))
             ELSE
                idx_IN = 1
             ENDIF

             !..Freezing of water drops into graupel/cloud ice (Bigg 1953).
             IF (rr(k).GT. r_r(1)) THEN
                prg_rfz(k) = tpg_qrfz(idx_r,idx_r1,idx_tc,idx_IN)*odts
                pri_rfz(k) = tpi_qrfz(idx_r,idx_r1,idx_tc,idx_IN)*odts
                pni_rfz(k) = tni_qrfz(idx_r,idx_r1,idx_tc,idx_IN)*odts
                pnr_rfz(k) = tnr_qrfz(idx_r,idx_r1,idx_tc,idx_IN)*odts          ! RAIN2M
                pnr_rfz(k) = MIN(DBLE(nr(k)*odts), pnr_rfz(k))
             ELSEIF (rr(k).GT. R1 .AND. temp(k).LT.HGFR) THEN
                pri_rfz(k) = rr(k)*odts
                pnr_rfz(k) = nr(k)*odts                                         ! RAIN2M
                pni_rfz(k) = pnr_rfz(k)
             ENDIF

             IF (rc(k).GT. r_c(1)) THEN
                pri_wfz(k) = tpi_qcfz(idx_c,idx_n,idx_tc,idx_IN)*odts
                pri_wfz(k) = MIN(DBLE(rc(k)*odts), pri_wfz(k))
                pni_wfz(k) = tni_qcfz(idx_c,idx_n,idx_tc,idx_IN)*odts
                pni_wfz(k) = MIN(DBLE(nc(k)*odts), pri_wfz(k)/(2.0_r8*xm0i),     &
                     pni_wfz(k))
             ELSEIF (rc(k).GT. R1 .AND. temp(k).LT.HGFR) THEN
                pri_wfz(k) = rc(k)*odts
                pni_wfz(k) = nc(k)*odts
             ENDIF

             !..Deposition nucleation of dust/mineral from DeMott et al (2010)
             !.. we may need to relax the temperature and ssati constraints.
             IF ( (ssati(k).GE. 0.25_r8) .OR. (ssatw(k).GT. eps &
                  .AND. temp(k).LT.261.15_r8) ) THEN
                IF (dustyIce .AND. is_aerosol_aware) THEN
                   xnc = iceDeMott(tempc,qv(k),qvs(k),qvsi(k),rho(k),nifa(k))
                ELSE
                   xnc = MIN(250.E3_r8, TNO*EXP(ATO*(T_0-temp(k))))
                ENDIF
                xni = ni(k) + (pni_rfz(k)+pni_wfz(k))*dtsave
                pni_inu(k) = 0.5_r8*(xnc-xni + ABS(xnc-xni))*odts
                pri_inu(k) = MIN(DBLE(rate_max), xm0i*pni_inu(k))
                pni_inu(k) = pri_inu(k)/xm0i
             ENDIF

             !..Freezing of aqueous aerosols based on Koop et al (2001, Nature)
             xni = smo0(k)+ni(k) + (pni_rfz(k)+pni_wfz(k)+pni_inu(k))*dtsave
             IF (is_aerosol_aware .AND. homogIce .AND. (xni.LE.500.E3)     &
                  &                .AND.(temp(k).LT.238).AND.(ssati(k).GE.0.4_r8) ) THEN
                xnc = iceKoop(temp(k),qv(k),qvs(k),nwfa(k), dtsave)
                pni_iha(k) = xnc*odts
                pri_iha(k) = MIN(DBLE(rate_max), xm0i*0.1_r8*pni_iha(k))
                pni_iha(k) = pri_iha(k)/(xm0i*0.1_r8)
             ENDIF
             !+---+------------------ END NEW ICE NUCLEATION -----------------------+


             !..Deposition/sublimation of cloud ice (Srivastava & Coen 1992).
             IF (L_qi(k)) THEN
                lami = (am_i*cig(2)*oig1*ni(k)/ri(k))**obmi
                ilami = 1.0_r8/lami
                xDi = MAX(DBLE(D0i), (bm_i + mu_i + 1.0_r8) * ilami)
                xmi = am_i*xDi**bm_i
                oxmi = 1.0_r8/xmi
                pri_ide(k) = C_cube*t1_subl*diffu(k)*ssati(k)*rvs &
                     *oig1*cig(5)*ni(k)*ilami

                IF (pri_ide(k) .LT. 0.0_r8) THEN
                   pri_ide(k) = MAX(DBLE(-ri(k)*odts), pri_ide(k), DBLE(rate_max))
                   pni_ide(k) = pri_ide(k)*oxmi
                   pni_ide(k) = MAX(DBLE(-ni(k)*odts), pni_ide(k))
                ELSE
                   pri_ide(k) = MIN(pri_ide(k), DBLE(rate_max))
                   prs_ide(k) = (1.0e0_r8-tpi_ide(idx_i,idx_i1))*pri_ide(k)
                   pri_ide(k) = tpi_ide(idx_i,idx_i1)*pri_ide(k)
                ENDIF

                !..Some cloud ice needs to move into the snow category.  Use lookup
                !.. table that resulted from explicit bin representation of distrib.
                IF ( (idx_i.EQ. ntb_i) .OR. (xDi.GT. 5.0_r8*D0s) ) THEN
                   prs_iau(k) = ri(k)*0.99_r8*odts
                   pni_iau(k) = ni(k)*0.95_r8*odts
                ELSEIF (xDi.LT. 0.1_r8*D0s) THEN
                   prs_iau(k) = 0.0_r8
                   pni_iau(k) = 0.0_r8
                ELSE
                   prs_iau(k) = tps_iaus(idx_i,idx_i1)*odts
                   prs_iau(k) = MIN(DBLE(ri(k)*0.99_r8*odts), prs_iau(k))
                   pni_iau(k) = tni_iaus(idx_i,idx_i1)*odts
                   pni_iau(k) = MIN(DBLE(ni(k)*0.95_r8*odts), pni_iau(k))
                ENDIF
             ENDIF

             !..Deposition/sublimation of snow/graupel follows Srivastava & Coen
             !.. (1992).
             IF (L_qs(k)) THEN
                C_snow = C_sqrd + (tempc+15.0_r8)*(C_cube-C_sqrd)/(-30.0_r8+15.0_r8)
                C_snow = MAX(C_sqrd, MIN(C_snow, C_cube))
                prs_sde(k) = C_snow*t1_subl*diffu(k)*ssati(k)*rvs &
                     * (t1_qs_sd*smo1(k) &
                     + t2_qs_sd*rhof2(k)*vsc2(k)*smof(k))
                IF (prs_sde(k).LT. 0.0_r8) THEN
                   prs_sde(k) = MAX(DBLE(-rs(k)*odts), prs_sde(k), DBLE(rate_max))
                ELSE
                   prs_sde(k) = MIN(prs_sde(k), DBLE(rate_max))
                ENDIF
             ENDIF

             IF (L_qg(k) .AND. ssati(k).LT. -eps) THEN
                prg_gde(k) = C_cube*t1_subl*diffu(k)*ssati(k)*rvs &
                     * N0_g(k) * (t1_qg_sd*ilamg(k)**cge(10) &
                     + t2_qg_sd*vsc2(k)*rhof2(k)*ilamg(k)**cge(11))
                IF (prg_gde(k).LT. 0.0_r8) THEN
                   prg_gde(k) = MAX(DBLE(-rg(k)*odts), prg_gde(k), DBLE(rate_max))
                ELSE
                   prg_gde(k) = MIN(prg_gde(k), DBLE(rate_max))
                ENDIF
             ENDIF

             !..Snow collecting cloud ice.  In CE, assume Di<<Ds and vti=~0.
             IF (L_qi(k)) THEN
                lami = (am_i*cig(2)*oig1*ni(k)/ri(k))**obmi
                ilami = 1.0_r8/lami
                xDi = MAX(DBLE(D0i), (bm_i + mu_i + 1.0_r8) * ilami)
                xmi = am_i*xDi**bm_i
                oxmi = 1.0_r8/xmi
                IF (rs(k).GE. r_s(1)) THEN
                   prs_sci(k) = t1_qs_qi*rhof(k)*Ef_si*ri(k)*smoe(k)
                   pni_sci(k) = prs_sci(k) * oxmi
                ENDIF

                !..Rain collecting cloud ice.  In CE, assume Di<<Dr and vti=~0.
                IF (rr(k).GE. r_r(1) .AND. mvd_r(k).GT. 4.0_r8*xDi) THEN
                   lamr = 1.0_r8/ilamr(k)
                   pri_rci(k) = rhof(k)*t1_qr_qi*Ef_ri*ri(k)*N0_r(k) &
                        *((lamr+fv_r)**(-cre(9)))
                   pnr_rci(k) = rhof(k)*t1_qr_qi*Ef_ri*ni(k)*N0_r(k)           &   ! RAIN2M
                        *((lamr+fv_r)**(-cre(9)))
                   pni_rci(k) = pri_rci(k) * oxmi
                   prr_rci(k) = rhof(k)*t2_qr_qi*Ef_ri*ni(k)*N0_r(k) &
                        *((lamr+fv_r)**(-cre(8)))
                   prr_rci(k) = MIN(DBLE(rr(k)*odts), prr_rci(k))
                   prg_rci(k) = pri_rci(k) + prr_rci(k)
                ENDIF
             ENDIF

             !..Ice multiplication from rime-splinters (Hallet & Mossop 1974).
             IF (prg_gcw(k).GT. eps .AND. tempc.GT.-8.0_r8) THEN
                tf = 0.0_r8
                IF (tempc.GE.-5.0_r8 .AND. tempc.LT.-3.0_r8) THEN
                   tf = 0.5_r8*(-3.0_r8 - tempc)
                ELSEIF (tempc.GT.-8.0_r8 .AND. tempc.LT.-5.0_r8) THEN
                   tf = 0.33333333_r8*(8.0_r8 + tempc)
                ENDIF
                pni_ihm(k) = 3.5E8_r8*tf*prg_gcw(k)
                pri_ihm(k) = xm0i*pni_ihm(k)
                prs_ihm(k) = prs_scw(k)/(prs_scw(k)+prg_gcw(k)) &
                     * pri_ihm(k)
                prg_ihm(k) = prg_gcw(k)/(prs_scw(k)+prg_gcw(k)) &
                     * pri_ihm(k)
             ENDIF

             !..A portion of rimed snow converts to graupel but some remains snow.
             !.. Interp from 5 to 75% as riming factor increases from 5.0 to 30.0
             !.. 0.028 came from (.75-.05)/(30.-5.).  This remains ad-hoc and should
             !.. be revisited.
             IF (prs_scw(k).GT.5.0_r8*prs_sde(k) .AND. &
                  prs_sde(k).GT.eps) THEN
                r_frac = MIN(30.0e0_r8, prs_scw(k)/prs_sde(k))
                g_frac = MIN(0.75_r8, 0.05_r8 + (r_frac-5.0_r8)*0.028_r8)
                vts_boost(k) = MIN(1.5_r8, 1.1_r8 + (r_frac-5.0_r8)*0.016_r8)
                prg_scw(k) = g_frac*prs_scw(k)
                prs_scw(k) = (1.0_r8 - g_frac)*prs_scw(k)
             ENDIF

          ELSE

             !..Melt snow and graupel and enhance from collisions with liquid.
             !.. We also need to sublimate snow and graupel if subsaturated.
             IF (L_qs(k)) THEN
                prr_sml(k) = (tempc*tcond(k)-lvap0*diffu(k)*delQvs(k))       &
                     * (t1_qs_me*smo1(k) + t2_qs_me*rhof2(k)*vsc2(k)*smof(k))
                prr_sml(k) = prr_sml(k) + 4218.0_r8*olfus*tempc &
                     * (prr_rcs(k)+prs_scw(k))
                prr_sml(k) = MIN(DBLE(rs(k)*odts), MAX(0.e0_r8, prr_sml(k)))
                pnr_sml(k) = smo0(k)/rs(k)*prr_sml(k) * 10.0_r8**(-0.75_r8*tempc)      ! RAIN2M
                pnr_sml(k) = MIN(DBLE(smo0(k)*odts), pnr_sml(k))
                IF (tempc.GT.3.5_r8 .OR. rs(k).LT.0.005E-3_r8) pnr_sml(k)=0.0_r8

                IF (ssati(k).LT. 0.0_r8) THEN
                   prs_sde(k) = C_cube*t1_subl*diffu(k)*ssati(k)*rvs &
                        * (t1_qs_sd*smo1(k) &
                        + t2_qs_sd*rhof2(k)*vsc2(k)*smof(k))
                   prs_sde(k) = MAX(DBLE(-rs(k)*odts), prs_sde(k))
                ENDIF
             ENDIF

             IF (L_qg(k)) THEN
                prr_gml(k) = (tempc*tcond(k)-lvap0*diffu(k)*delQvs(k))       &
                     * N0_g(k)*(t1_qg_me*ilamg(k)**cge(10)             &
                     + t2_qg_me*rhof2(k)*vsc2(k)*ilamg(k)**cge(11))
                !-GT       prr_gml(k) = prr_gml(k) + 4218.*olfus*tempc &
                !-GT                               * (prr_rcg(k)+prg_gcw(k))
                prr_gml(k) = MIN(DBLE(rg(k)*odts), MAX(0.e0_r8, prr_gml(k)))
                pnr_gml(k) = N0_g(k)*cgg(2)*ilamg(k)**cge(2) / rg(k)         &   ! RAIN2M
                     * prr_gml(k) * 10.0_r8**(-0.25_r8*tempc)
                IF (tempc.GT.7.5_r8 .OR. rg(k).LT.0.005E-3_r8) pnr_gml(k)=0.0_r8

                IF (ssati(k).LT. 0.0_r8) THEN
                   prg_gde(k) = C_cube*t1_subl*diffu(k)*ssati(k)*rvs &
                        * N0_g(k) * (t1_qg_sd*ilamg(k)**cge(10) &
                        + t2_qg_sd*vsc2(k)*rhof2(k)*ilamg(k)**cge(11))
                   prg_gde(k) = MAX(DBLE(-rg(k)*odts), prg_gde(k))
                ENDIF
             ENDIF

             !.. This change will be required if users run adaptive time step that
             !.. results in delta-t that is generally too long to allow cloud water
             !.. collection by snow/graupel above melting temperature.
             !.. Credit to Bjorn-Egil Nygaard for discovering.
             IF (dt .GT. 120.0_r8) THEN
                prr_rcw(k)=prr_rcw(k)+prs_scw(k)+prg_gcw(k)
                prs_scw(k)=0.0_r8
                prg_gcw(k)=0.0_r8
             ENDIF

          ENDIF

       ENDDO
    ENDIF

    !+---+-----------------------------------------------------------------+
    !..Ensure we do not deplete more hydrometeor species than exists.
    !+---+-----------------------------------------------------------------+
    DO k = kts, kte

       !..If ice supersaturated, ensure sum of depos growth terms does not
       !.. deplete more vapor than possibly exists.  If subsaturated, limit
       !.. sum of sublimation terms such that vapor does not reproduce ice
       !.. supersat again.
       sump = pri_inu(k) + pri_ide(k) + prs_ide(k) &
            + prs_sde(k) + prg_gde(k) + pri_iha(k)
       rate_max = (qv(k)-qvsi(k))*odts*0.999_r8
       IF ( (sump.GT. eps .AND. sump.GT. rate_max) .OR. &
            (sump.LT. -eps .AND. sump.LT. rate_max) ) THEN
          ratio = rate_max/sump
          pri_inu(k) = pri_inu(k) * ratio
          pri_ide(k) = pri_ide(k) * ratio
          pni_ide(k) = pni_ide(k) * ratio
          prs_ide(k) = prs_ide(k) * ratio
          prs_sde(k) = prs_sde(k) * ratio
          prg_gde(k) = prg_gde(k) * ratio
          pri_iha(k) = pri_iha(k) * ratio
       ENDIF

       !..Cloud water conservation.
       sump = -prr_wau(k) - pri_wfz(k) - prr_rcw(k) &
            - prs_scw(k) - prg_scw(k) - prg_gcw(k)
       rate_max = -rc(k)*odts
       IF (sump.LT. rate_max .AND. L_qc(k)) THEN
          ratio = rate_max/sump
          prr_wau(k) = prr_wau(k) * ratio
          pri_wfz(k) = pri_wfz(k) * ratio
          prr_rcw(k) = prr_rcw(k) * ratio
          prs_scw(k) = prs_scw(k) * ratio
          prg_scw(k) = prg_scw(k) * ratio
          prg_gcw(k) = prg_gcw(k) * ratio
       ENDIF

       !..Cloud ice conservation.
       sump = pri_ide(k) - prs_iau(k) - prs_sci(k) &
            - pri_rci(k)
       rate_max = -ri(k)*odts
       IF (sump.LT. rate_max .AND. L_qi(k)) THEN
          ratio = rate_max/sump
          pri_ide(k) = pri_ide(k) * ratio
          prs_iau(k) = prs_iau(k) * ratio
          prs_sci(k) = prs_sci(k) * ratio
          pri_rci(k) = pri_rci(k) * ratio
       ENDIF

       !..Rain conservation.
       sump = -prg_rfz(k) - pri_rfz(k) - prr_rci(k) &
            + prr_rcs(k) + prr_rcg(k)
       rate_max = -rr(k)*odts
       IF (sump.LT. rate_max .AND. L_qr(k)) THEN
          ratio = rate_max/sump
          prg_rfz(k) = prg_rfz(k) * ratio
          pri_rfz(k) = pri_rfz(k) * ratio
          prr_rci(k) = prr_rci(k) * ratio
          prr_rcs(k) = prr_rcs(k) * ratio
          prr_rcg(k) = prr_rcg(k) * ratio
       ENDIF

       !..Snow conservation.
       sump = prs_sde(k) - prs_ihm(k) - prr_sml(k) &
            + prs_rcs(k)
       rate_max = -rs(k)*odts
       IF (sump.LT. rate_max .AND. L_qs(k)) THEN
          ratio = rate_max/sump
          prs_sde(k) = prs_sde(k) * ratio
          prs_ihm(k) = prs_ihm(k) * ratio
          prr_sml(k) = prr_sml(k) * ratio
          prs_rcs(k) = prs_rcs(k) * ratio
       ENDIF

       !..Graupel conservation.
       sump = prg_gde(k) - prg_ihm(k) - prr_gml(k) &
            + prg_rcg(k)
       rate_max = -rg(k)*odts
       IF (sump.LT. rate_max .AND. L_qg(k)) THEN
          ratio = rate_max/sump
          prg_gde(k) = prg_gde(k) * ratio
          prg_ihm(k) = prg_ihm(k) * ratio
          prr_gml(k) = prr_gml(k) * ratio
          prg_rcg(k) = prg_rcg(k) * ratio
       ENDIF

       !..Re-enforce proper mass conservation for subsequent elements in case
       !.. any of the above terms were altered.  Thanks P. Blossey. 2009Sep28
       pri_ihm(k) = prs_ihm(k) + prg_ihm(k)
       ratio = MIN( ABS(prr_rcg(k)), ABS(prg_rcg(k)) )
       prr_rcg(k) = ratio * SIGN(1.0_r8, REAL(SNGL(prr_rcg(k)),kind=r8))
       prg_rcg(k) = -prr_rcg(k)
       IF (temp(k).GT.T_0) THEN
          ratio = MIN( ABS(prr_rcs(k)), ABS(prs_rcs(k)) )
          prr_rcs(k) = ratio * SIGN(1.0_r8, REAL(SNGL(prr_rcs(k)),kind=r8))
          prs_rcs(k) = -prr_rcs(k)
       ENDIF

    ENDDO

    !+---+-----------------------------------------------------------------+
    !..Calculate tendencies of all species but constrain the number of ice
    !.. to reasonable values.
    !+---+-----------------------------------------------------------------+
    DO k = kts, kte
       orho = 1.0_r8/rho(k)
       lfus2 = lsub - lvap(k)

       !..Aerosol number tendency
       IF (is_aerosol_aware) THEN
          nwfaten(k) = nwfaten(k) - (pna_rca(k) + pna_sca(k)          &
               + pna_gca(k) + pni_iha(k)) * orho
          nifaten(k) = nifaten(k) - (pnd_rcd(k) + pnd_scd(k)          &
               + pnd_gcd(k)) * orho
          IF (dustyIce) THEN
             nifaten(k) = nifaten(k) - pni_inu(k)*orho
          ELSE
             nifaten(k) = 0.0_r8
          ENDIF
       ENDIF

       !..Water vapor tendency
       qvten(k) = qvten(k) + (-pri_inu(k) - pri_iha(k) - pri_ide(k)   &
            - prs_ide(k) - prs_sde(k) - prg_gde(k)) &
            * orho

       !..Cloud water tendency
       qcten(k) = qcten(k) + (-prr_wau(k) - pri_wfz(k) &
            - prr_rcw(k) - prs_scw(k) - prg_scw(k) &
            - prg_gcw(k)) &
            * orho

       !..Cloud water number tendency
       ncten(k) = ncten(k) + (-pnc_wau(k) - pnc_rcw(k) &
            - pni_wfz(k) - pnc_scw(k) - pnc_gcw(k)) &
            * orho

       !..Cloud water mass/number balance; keep mass-wt mean size between
       !.. 1 and 50 microns.  Also no more than Nt_c_max drops total.
       xrc=MAX(R1, (qc1d(k) + qcten(k)*dtsave)*rho(k))
       xnc=MAX(2.0_r8, (nc1d(k) + ncten(k)*dtsave)*rho(k))
       IF (xrc .GT. R1) THEN
          nu_c = MIN(15, NINT(1000.E6_r8/xnc) + 2)
          lamc = (xnc*am_r*ccg(2,nu_c)*ocg1(nu_c)/rc(k))**obmr
          xDc = (bm_r + nu_c + 1.0_r8) / lamc
          IF (xDc.LT. D0c) THEN
             lamc = cce(2,nu_c)/D0c
             xnc = ccg(1,nu_c)*ocg2(nu_c)*xrc/am_r*lamc**bm_r
             ncten(k) = (xnc-nc1d(k)*rho(k))*odts*orho
          ELSEIF (xDc.GT. D0r*2.0_r8) THEN
             lamc = cce(2,nu_c)/(D0r*2.0_r8)
             xnc = ccg(1,nu_c)*ocg2(nu_c)*xrc/am_r*lamc**bm_r
             ncten(k) = (xnc-nc1d(k)*rho(k))*odts*orho
          ENDIF
       ELSE
          ncten(k) = -nc1d(k)*odts
       ENDIF
       xnc=MAX(0.0_r8,(nc1d(k) + ncten(k)*dtsave)*rho(k))
       IF (xnc.GT.Nt_c_max) &
            ncten(k) = (Nt_c_max-nc1d(k)*rho(k))*odts*orho

       !..Cloud ice mixing ratio tendency
       qiten(k) = qiten(k) + (pri_inu(k) + pri_iha(k) + pri_ihm(k)    &
            + pri_wfz(k) + pri_rfz(k) + pri_ide(k) &
            - prs_iau(k) - prs_sci(k) - pri_rci(k)) &
            * orho

       !..Cloud ice number tendency.
       niten(k) = niten(k) + (pni_inu(k) + pni_iha(k) + pni_ihm(k)    &
            + pni_wfz(k) + pni_rfz(k) + pni_ide(k) &
            - pni_iau(k) - pni_sci(k) - pni_rci(k)) &
            * orho

       !..Cloud ice mass/number balance; keep mass-wt mean size between
       !.. 20 and 300 microns.  Also no more than 250 xtals per liter.
       xri=MAX(R1,(qi1d(k) + qiten(k)*dtsave)*rho(k))
       xni=MAX(R2,(ni1d(k) + niten(k)*dtsave)*rho(k))
       IF (xri.GT. R1) THEN
          lami = (am_i*cig(2)*oig1*xni/xri)**obmi
          ilami = 1.0_r8/lami
          xDi = (bm_i + mu_i + 1.0_r8) * ilami
          IF (xDi.LT. 20.E-6_r8) THEN
             lami = cie(2)/20.E-6_r8
             xni = MIN(499.e3_r8, cig(1)*oig2*xri/am_i*lami**bm_i)
             niten(k) = (xni-ni1d(k)*rho(k))*odts*orho
          ELSEIF (xDi.GT. 300.E-6_r8) THEN
             lami = cie(2)/300.E-6_r8
             xni = cig(1)*oig2*xri/am_i*lami**bm_i
             niten(k) = (xni-ni1d(k)*rho(k))*odts*orho
          ENDIF
       ELSE
          niten(k) = -ni1d(k)*odts
       ENDIF
       xni=MAX(0.0_r8,(ni1d(k) + niten(k)*dtsave)*rho(k))
       IF (xni.GT.499.E3_r8) &
            niten(k) = (499.E3_r8-ni1d(k)*rho(k))*odts*orho

       !..Rain tendency
       qrten(k) = qrten(k) + (prr_wau(k) + prr_rcw(k) &
            + prr_sml(k) + prr_gml(k) + prr_rcs(k) &
            + prr_rcg(k) - prg_rfz(k) &
            - pri_rfz(k) - prr_rci(k)) &
            * orho

       !..Rain number tendency
       nrten(k) = nrten(k) + (pnr_wau(k) + pnr_sml(k) + pnr_gml(k)    &
            - (pnr_rfz(k) + pnr_rcr(k) + pnr_rcg(k)           &
            + pnr_rcs(k) + pnr_rci(k)) )                      &
            * orho

       !..Rain mass/number balance; keep median volume diameter between
       !.. 37 microns (D0r*0.75) and 2.5 mm.
       xrr=MAX(R1,(qr1d(k) + qrten(k)*dtsave)*rho(k))
       xnr=MAX(R2,(nr1d(k) + nrten(k)*dtsave)*rho(k))
       IF (xrr.GT. R1) THEN
          lamr = (am_r*crg(3)*org2*xnr/xrr)**obmr
          mvd_r(k) = (3.0_r8 + mu_r + 0.672_r8) / lamr
          IF (mvd_r(k) .GT. 2.5E-3_r8) THEN
             mvd_r(k) = 2.5E-3_r8
             lamr = (3.0_r8 + mu_r + 0.672_r8) / mvd_r(k)
             xnr = crg(2)*org3*xrr*lamr**bm_r / am_r
             nrten(k) = (xnr-nr1d(k)*rho(k))*odts*orho
          ELSEIF (mvd_r(k) .LT. D0r*0.75_r8) THEN
             mvd_r(k) = D0r*0.75_r8
             lamr = (3.0_r8 + mu_r + 0.672_r8) / mvd_r(k)
             xnr = crg(2)*org3*xrr*lamr**bm_r / am_r
             nrten(k) = (xnr-nr1d(k)*rho(k))*odts*orho
          ENDIF
       ELSE
          qrten(k) = -qr1d(k)*odts
          nrten(k) = -nr1d(k)*odts
       ENDIF

       !..Snow tendency
       qsten(k) = qsten(k) + (prs_iau(k) + prs_sde(k) &
            + prs_sci(k) + prs_scw(k) + prs_rcs(k) &
            + prs_ide(k) - prs_ihm(k) - prr_sml(k)) &
            * orho

       !..Graupel tendency
       qgten(k) = qgten(k) + (prg_scw(k) + prg_rfz(k) &
            + prg_gde(k) + prg_rcg(k) + prg_gcw(k) &
            + prg_rci(k) + prg_rcs(k) - prg_ihm(k) &
            - prr_gml(k)) &
            * orho

       !..Temperature tendency
       IF (temp(k).LT.T_0) THEN
          tten(k) = tten(k) &
               + ( lsub*ocp(k)*(pri_inu(k) + pri_ide(k) &
               + prs_ide(k) + prs_sde(k) &
               + prg_gde(k) + pri_iha(k)) &
               + lfus2*ocp(k)*(pri_wfz(k) + pri_rfz(k) &
               + prg_rfz(k) + prs_scw(k) &
               + prg_scw(k) + prg_gcw(k) &
               + prg_rcs(k) + prs_rcs(k) &
               + prr_rci(k) + prg_rcg(k)) &
               )*orho * (1-IFDRY)
       ELSE
          tten(k) = tten(k) &
               + ( lfus*ocp(k)*(-prr_sml(k) - prr_gml(k) &
               - prr_rcg(k) - prr_rcs(k)) &
               + lsub*ocp(k)*(prs_sde(k) + prg_gde(k)) &
               )*orho * (1-IFDRY)
       ENDIF

    ENDDO

    !+---+-----------------------------------------------------------------+
    !..Update variables for TAU+1 before condensation & sedimention.
    !+---+-----------------------------------------------------------------+
    DO k = kts, kte
       temp(k) = t1d(k) + DT*tten(k)
       otemp = 1.0_r8/temp(k)
       tempc = temp(k) - 273.15_r8
       qv(k) = MAX(1.E-10_r8, qv1d(k) + DT*qvten(k))
       rho(k) = 0.622_r8*pres(k)/(R*temp(k)*(qv(k)+0.622_r8))
       rhof(k) = SQRT(RHO_NOT/rho(k))
       rhof2(k) = SQRT(rhof(k))
       qvs(k) = rslf(pres(k), temp(k))
       ssatw(k) = qv(k)/qvs(k) - 1.0_r8
       IF (ABS(ssatw(k)).LT. eps) ssatw(k) = 0.0_r8
       diffu(k) = 2.11E-5_r8*(temp(k)/273.15_r8)**1.94_r8 * (101325.0_r8/pres(k))
       IF (tempc .GE. 0.0_r8) THEN
          visco(k) = (1.718_r8+0.0049_r8*tempc)*1.0E-5_r8
       ELSE
          visco(k) = (1.718_r8+0.0049_r8*tempc-1.2E-5_r8*tempc*tempc)*1.0E-5_r8
       ENDIF
       vsc2(k) = SQRT(rho(k)/visco(k))
       lvap(k) = lvap0 + (2106.0_r8 - 4218.0_r8)*tempc
       tcond(k) = (5.69_r8 + 0.0168_r8*tempc)*1.0E-5_r8 * 418.936_r8
       ocp(k) = 1.0_r8/(Cp*(1.0_r8+0.887_r8*qv(k)))
       lvt2(k)=lvap(k)*lvap(k)*ocp(k)*oRv*otemp*otemp

       nwfa(k) = MAX(11.1E6_r8, (nwfa1d(k) + nwfaten(k)*DT)*rho(k))

       IF ((qc1d(k) + qcten(k)*DT) .GT. R1) THEN
          rc(k) = (qc1d(k) + qcten(k)*DT)*rho(k)
          nc(k) = MAX(2.0_r8, (nc1d(k) + ncten(k)*DT)*rho(k))
          IF (.NOT. is_aerosol_aware) nc(k) = Nt_c
          L_qc(k) = .TRUE.
       ELSE
          rc(k) = R1
          nc(k) = 2.0_r8
          L_qc(k) = .FALSE.
       ENDIF

       IF ((qi1d(k) + qiten(k)*DT) .GT. R1) THEN
          ri(k) = (qi1d(k) + qiten(k)*DT)*rho(k)
          ni(k) = MAX(R2, (ni1d(k) + niten(k)*DT)*rho(k))
          L_qi(k) = .TRUE. 
       ELSE
          ri(k) = R1
          ni(k) = R2
          L_qi(k) = .FALSE.
       ENDIF

       IF ((qr1d(k) + qrten(k)*DT) .GT. R1) THEN
          rr(k) = (qr1d(k) + qrten(k)*DT)*rho(k)
          nr(k) = MAX(R2, (nr1d(k) + nrten(k)*DT)*rho(k))
          L_qr(k) = .TRUE.
          lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
          mvd_r(k) = (3.0_r8 + mu_r + 0.672_r8) / lamr
          IF (mvd_r(k) .GT. 2.5E-3_r8) THEN
             mvd_r(k) = 2.5E-3_r8
             lamr = (3.0_r8 + mu_r + 0.672_r8) / mvd_r(k)
             nr(k) = crg(2)*org3*rr(k)*lamr**bm_r / am_r
          ELSEIF (mvd_r(k) .LT. D0r*0.75_r8) THEN
             mvd_r(k) = D0r*0.75_r8
             lamr = (3.0_r8 + mu_r + 0.672_r8) / mvd_r(k)
             nr(k) = crg(2)*org3*rr(k)*lamr**bm_r / am_r
          ENDIF
       ELSE
          rr(k) = R1
          nr(k) = R2
          L_qr(k) = .FALSE.
       ENDIF

       IF ((qs1d(k) + qsten(k)*DT) .GT. R1) THEN
          rs(k) = (qs1d(k) + qsten(k)*DT)*rho(k)
          L_qs(k) = .TRUE.
       ELSE
          rs(k) = R1
          L_qs(k) = .FALSE.
       ENDIF

       IF ((qg1d(k) + qgten(k)*DT) .GT. R1) THEN
          rg(k) = (qg1d(k) + qgten(k)*DT)*rho(k)
          L_qg(k) = .TRUE.
       ELSE
          rg(k) = R1
          L_qg(k) = .FALSE.
       ENDIF
    ENDDO

    !+---+-----------------------------------------------------------------+
    !..With tendency-updated mixing ratios, recalculate snow moments and
    !.. intercepts/slopes of graupel and rain.
    !+---+-----------------------------------------------------------------+
    IF (.NOT. iiwarm) THEN
       DO k = kts, kte
          IF (.NOT. L_qs(k)) CYCLE
          tc0 = MIN(-0.1_r8, temp(k)-273.15_r8)
          smob(k) = rs(k)*oams

          !..All other moments based on reference, 2nd moment.  If bm_s.ne.2,
          !.. then we must compute actual 2nd moment and use as reference.
          IF (bm_s.GT.(2.0_r8-1.e-3_r8) .AND. bm_s.LT.(2.0_r8+1.e-3_r8)) THEN
             smo2(k) = smob(k)
          ELSE
             loga_ = sa(1) + sa(2)*tc0 + sa(3)*bm_s &
                  + sa(4)*tc0*bm_s + sa(5)*tc0*tc0 &
                  + sa(6)*bm_s*bm_s + sa(7)*tc0*tc0*bm_s &
                  + sa(8)*tc0*bm_s*bm_s + sa(9)*tc0*tc0*tc0 &
                  + sa(10)*bm_s*bm_s*bm_s
             a_ = 10.0_r8**loga_
             b_ = sb(1) + sb(2)*tc0 + sb(3)*bm_s &
                  + sb(4)*tc0*bm_s + sb(5)*tc0*tc0 &
                  + sb(6)*bm_s*bm_s + sb(7)*tc0*tc0*bm_s &
                  + sb(8)*tc0*bm_s*bm_s + sb(9)*tc0*tc0*tc0 &
                  + sb(10)*bm_s*bm_s*bm_s
             smo2(k) = (smob(k)/a_)**(1.0_r8/b_)
          ENDIF

          !..Calculate bm_s+1 (th) moment.  Useful for diameter calcs.
          loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(1) &
               + sa(4)*tc0*cse(1) + sa(5)*tc0*tc0 &
               + sa(6)*cse(1)*cse(1) + sa(7)*tc0*tc0*cse(1) &
               + sa(8)*tc0*cse(1)*cse(1) + sa(9)*tc0*tc0*tc0 &
               + sa(10)*cse(1)*cse(1)*cse(1)
          a_ = 10.0_r8**loga_
          b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(1) + sb(4)*tc0*cse(1) &
               + sb(5)*tc0*tc0 + sb(6)*cse(1)*cse(1) &
               + sb(7)*tc0*tc0*cse(1) + sb(8)*tc0*cse(1)*cse(1) &
               + sb(9)*tc0*tc0*tc0 + sb(10)*cse(1)*cse(1)*cse(1)
          smoc(k) = a_ * smo2(k)**b_

          !..Calculate bm_s+bv_s (th) moment.  Useful for sedimentation.
          loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(14) &
               + sa(4)*tc0*cse(14) + sa(5)*tc0*tc0 &
               + sa(6)*cse(14)*cse(14) + sa(7)*tc0*tc0*cse(14) &
               + sa(8)*tc0*cse(14)*cse(14) + sa(9)*tc0*tc0*tc0 &
               + sa(10)*cse(14)*cse(14)*cse(14)
          a_ = 10.0_r8**loga_
          b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(14) + sb(4)*tc0*cse(14) &
               + sb(5)*tc0*tc0 + sb(6)*cse(14)*cse(14) &
               + sb(7)*tc0*tc0*cse(14) + sb(8)*tc0*cse(14)*cse(14) &
               + sb(9)*tc0*tc0*tc0 + sb(10)*cse(14)*cse(14)*cse(14)
          smod(k) = a_ * smo2(k)**b_
       ENDDO

       !+---+-----------------------------------------------------------------+
       !..Calculate y-intercept, slope values for graupel.
       !+---+-----------------------------------------------------------------+
       N0_min = gonv_max
       DO k = kte, kts, -1
          IF (temp(k).LT.270.65_r8 .AND. L_qr(k) .AND. mvd_r(k).GT.100.E-6_r8) THEN
             xslw1 = 4.01_r8 + log10(mvd_r(k))
          ELSE
             xslw1 = 0.01_r8
          ENDIF
          ygra1 = 4.31_r8 + log10(MAX(5.E-5_r8, rg(k)))
          zans1 = 3.1_r8 + (100.0_r8/(300.0_r8*xslw1*ygra1/(10.0_r8/xslw1+1.0_r8+0.25_r8*ygra1)+30.0_r8+10.0_r8*ygra1))
          N0_exp = 10.0_r8**(zans1)
          N0_exp = MAX(DBLE(gonv_min), MIN(N0_exp, DBLE(gonv_max)))
          N0_min = MIN(N0_exp, N0_min)
          N0_exp = N0_min
          lam_exp = (N0_exp*am_g*cgg(1)/rg(k))**oge1
          lamg = lam_exp * (cgg(3)*ogg2*ogg1)**obmg
          ilamg(k) = 1.0_r8/lamg
          N0_g(k) = N0_exp/(cgg(2)*lam_exp) * lamg**cge(2)
       ENDDO

    ENDIF

    !+---+-----------------------------------------------------------------+
    !..Calculate y-intercept, slope values for rain.
    !+---+-----------------------------------------------------------------+
    DO k = kte, kts, -1
       lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
       ilamr(k) = 1.0_r8/lamr
       mvd_r(k) = (3.0_r8 + mu_r + 0.672_r8) / lamr
       N0_r(k) = nr(k)*org2*lamr**cre(2)
    ENDDO

    !+---+-----------------------------------------------------------------+
    !..Cloud water condensation and evaporation.  Nucleate cloud droplets
    !.. using explicit CCN aerosols with hygroscopicity like sulfates using
    !.. parcel model lookup table results provided by T. Eidhammer.  Evap
    !.. drops using calculation of max drop size capable of evaporating in
    !.. single timestep and explicit number of drops smaller than Dc_star
    !.. from lookup table.
    !+---+-----------------------------------------------------------------+
    DO k = kts, kte
       orho = 1.0_r8/rho(k)
       IF ( (ssatw(k).GT. eps) .OR. (ssatw(k).LT. -eps .AND. &
            L_qc(k)) ) THEN
          clap = (qv(k)-qvs(k))/(1.0_r8 + lvt2(k)*qvs(k))
          DO n = 1, 3
             fcd = qvs(k)* EXP(lvt2(k)*clap) - qv(k) + clap
             dfcd = qvs(k)*lvt2(k)* EXP(lvt2(k)*clap) + 1.0_r8
             clap = clap - fcd/dfcd
          ENDDO
          xrc = rc(k) + clap*rho(k)
          xnc = 0.0_r8
          IF (xrc.GT. R1) THEN
             prw_vcd(k) = clap*odt
             !+---+-----------------------------------------------------------------+ !  DROPLET NUCLEATION
             IF (clap .GT. eps) THEN
                IF (is_aerosol_aware) THEN
                   xnc = MAX(2.0_r8, activ_ncloud(temp(k), w1d(k), nwfa(k)))
                ELSE
                   xnc = Nt_c
                ENDIF
                pnc_wcd(k) = 0.5_r8*(xnc-nc(k) + ABS(xnc-nc(k)))*odts*orho

                !+---+-----------------------------------------------------------------+ !  EVAPORATION
             ELSEIF (clap .LT. -eps .AND. ssatw(k).LT.-1.E-6_r8 .AND. is_aerosol_aware) THEN
                tempc = temp(k) - 273.15_r8
                otemp = 1.0_r8/temp(k)
                rvs = rho(k)*qvs(k)
                rvs_p = rvs*otemp*(lvap(k)*otemp*oRv - 1.0_r8)
                rvs_pp = rvs * ( otemp*(lvap(k)*otemp*oRv - 1.0_r8) &
                     *otemp*(lvap(k)*otemp*oRv - 1.0_r8) &
                     + (-2.0_r8*lvap(k)*otemp*otemp*otemp*oRv) &
                     + otemp*otemp)
                gamsc = lvap(k)*diffu(k)/tcond(k) * rvs_p
                alphsc = 0.5_r8*(gamsc/(1.0_r8+gamsc))*(gamsc/(1.0_r8+gamsc)) &
                     * rvs_pp/rvs_p * rvs/rvs_p
                alphsc = MAX(1.E-9_r8, alphsc)
                xsat = ssatw(k)
                IF (ABS(xsat).LT. 1.E-9_r8) xsat=0.0_r8
                t1_evap = 2.0_r8*PI*( 1.0_r8 - alphsc*xsat  &
                     + 2.0_r8*alphsc*alphsc*xsat*xsat  &
                     - 5.0_r8*alphsc*alphsc*alphsc*xsat*xsat*xsat ) &
                     / (1.0_r8+gamsc)

                Dc_star = DSQRT(-2.e0_r8*DT * t1_evap/(2.0_r8*PI) &
                     * 4.0_r8*diffu(k)*ssatw(k)*rvs/rho_w)
                idx_d = MAX(1, MIN(INT(1.E6_r8*Dc_star), nbc))

                idx_n = NINT(1.0_r8 + FLOAT(nbc) * DLOG(nc(k)/t_Nc(1)) / nic1)
                idx_n = MAX(1, MIN(idx_n, nbc))

                !..Cloud water lookup table index.
                IF (rc(k).GT. r_c(1)) THEN
                   nic = NINT(log10(rc(k)))
                   DO nn = nic-1, nic+1
                      n = nn
                      IF ( (rc(k)/10.0_r8**nn).GE.1.0_r8 .AND. &
                           (rc(k)/10.0_r8**nn).LT.10.0_r8) GOTO 159
                   ENDDO
159                CONTINUE
                   idx_c = INT(rc(k)/10.0_r8**n) + 10*(n-nic2) - (n-nic2)
                   idx_c = MAX(1, MIN(idx_c, ntb_c))
                ELSE
                   idx_c = 1
                ENDIF

                !prw_vcd(k) = MAX(DBLE(-rc(k)*orho*odt),                     &
                !           -tpc_wev(idx_d, idx_c, idx_n)*orho*odt)
                prw_vcd(k) = MAX(DBLE(-rc(k)*0.99_r8*orho*odt), prw_vcd(k))
                pnc_wcd(k) = MAX(DBLE(-nc(k)*0.99_r8*orho*odt),                &
                     -tnc_wev(idx_d, idx_c, idx_n)*orho*odt)

             ENDIF
          ELSE
             prw_vcd(k) = -rc(k)*orho*odt
             pnc_wcd(k) = -nc(k)*orho*odt
          ENDIF

          !+---+-----------------------------------------------------------------+

          qvten(k) = qvten(k) - prw_vcd(k)
          qcten(k) = qcten(k) + prw_vcd(k)
          ncten(k) = ncten(k) + pnc_wcd(k)
          nwfaten(k) = nwfaten(k) - pnc_wcd(k)
          tten(k) = tten(k) + lvap(k)*ocp(k)*prw_vcd(k)*(1-IFDRY)
          rc(k) = MAX(R1, (qc1d(k) + DT*qcten(k))*rho(k))
          nc(k) = MAX(2.0_r8, (nc1d(k) + DT*ncten(k))*rho(k))
          IF (.NOT. is_aerosol_aware) nc(k) = Nt_c
          qv(k) = MAX(1.E-10_r8, qv1d(k) + DT*qvten(k))
          temp(k) = t1d(k) + DT*tten(k)
          rho(k) = 0.622_r8*pres(k)/(R*temp(k)*(qv(k)+0.622_r8))
          qvs(k) = rslf(pres(k), temp(k))
          ssatw(k) = qv(k)/qvs(k) - 1.0_r8
       ENDIF
    ENDDO

    !+---+-----------------------------------------------------------------+
    !.. If still subsaturated, allow rain to evaporate, following
    !.. Srivastava & Coen (1992).
    !+---+-----------------------------------------------------------------+
    DO k = kts, kte
       IF ( (ssatw(k).LT. -eps) .AND. L_qr(k) &
            .AND. (.NOT.(prw_vcd(k).GT. 0.0_r8)) ) THEN
          tempc = temp(k) - 273.15_r8
          otemp = 1.0_r8/temp(k)
          orho = 1.0_r8/rho(k)
          rhof(k) = SQRT(RHO_NOT*orho)
          rhof2(k) = SQRT(rhof(k))
          diffu(k) = 2.11E-5_r8*(temp(k)/273.15_r8)**1.94_r8 * (101325.0_r8/pres(k))
          IF (tempc .GE. 0.0_r8) THEN
             visco(k) = (1.718_r8+0.0049_r8*tempc)*1.0E-5_r8
          ELSE
             visco(k) = (1.718_r8+0.0049_r8*tempc-1.2E-5_r8*tempc*tempc)*1.0E-5_r8
          ENDIF
          vsc2(k) = SQRT(rho(k)/visco(k))
          lvap(k) = lvap0 + (2106.0_r8 - 4218.0_r8)*tempc
          tcond(k) = (5.69_r8 + 0.0168_r8*tempc)*1.0E-5_r8 * 418.936_r8
          ocp(k) = 1.0_r8/(Cp*(1.0_r8+0.887_r8*qv(k)))

          rvs = rho(k)*qvs(k)
          rvs_p = rvs*otemp*(lvap(k)*otemp*oRv - 1.0_r8)
          rvs_pp = rvs * ( otemp*(lvap(k)*otemp*oRv - 1.0_r8) &
               *otemp*(lvap(k)*otemp*oRv - 1.0_r8) &
               + (-2.0_r8*lvap(k)*otemp*otemp*otemp*oRv) &
               + otemp*otemp)
          gamsc = lvap(k)*diffu(k)/tcond(k) * rvs_p
          alphsc = 0.5_r8*(gamsc/(1.0_r8+gamsc))*(gamsc/(1.0_r8+gamsc)) &
               * rvs_pp/rvs_p * rvs/rvs_p
          alphsc = MAX(1.E-9_r8, alphsc)
          xsat   = MIN(-1.E-9_r8, ssatw(k))
          t1_evap = 2.0_r8*PI*( 1.0_r8 - alphsc*xsat  &
               + 2.0_r8*alphsc*alphsc*xsat*xsat  &
               - 5.0_r8*alphsc*alphsc*alphsc*xsat*xsat*xsat ) &
               / (1.0_r8+gamsc)

          lamr = 1.0_r8/ilamr(k)
          !..Rapidly eliminate near zero values when low humidity (<95%)
          IF (qv(k)/qvs(k) .LT. 0.95_r8 .AND. rr(k)*orho.LE.1.E-8_r8) THEN
             prv_rev(k) = rr(k)*orho*odts
          ELSE
             prv_rev(k) = t1_evap*diffu(k)*(-ssatw(k))*N0_r(k)*rvs &
                  * (t1_qr_ev*ilamr(k)**cre(10) &
                  + t2_qr_ev*vsc2(k)*rhof2(k)*((lamr+0.5_r8*fv_r)**(-cre(11))))
             rate_max = MIN((rr(k)*orho*odts), (qvs(k)-qv(k))*odts)
             prv_rev(k) = MIN(DBLE(rate_max), prv_rev(k)*orho)

             !..TEST: G. Thompson  10 May 2013
             !..Reduce the rain evaporation in same places as melting graupel occurs. 
             !..Rationale: falling and simultaneous melting graupel in subsaturated
             !..regions will not melt as fast because particle temperature stays
             !..at 0C.  Also not much shedding of the water from the graupel so
             !..likely that the water-coated graupel evaporating much slower than
             !..if the water was immediately shed off.
             IF (prr_gml(k).GT.0.0_r8) THEN
                eva_factor = MIN(1.0_r8, 0.01_r8+(0.99_r8-0.01_r8)*(tempc/20.0_r8))
                prv_rev(k) = prv_rev(k)*eva_factor
             ENDIF
          ENDIF

          pnr_rev(k) = MIN(DBLE(nr(k)*0.99_r8*orho*odts),                  &   ! RAIN2M
               prv_rev(k) * nr(k)/rr(k))

          qrten(k) = qrten(k) - prv_rev(k)
          qvten(k) = qvten(k) + prv_rev(k)
          nrten(k) = nrten(k) - pnr_rev(k)
          nwfaten(k) = nwfaten(k) + pnr_rev(k)
          tten(k) = tten(k) - lvap(k)*ocp(k)*prv_rev(k)*(1-IFDRY)

          rr(k) = MAX(R1, (qr1d(k) + DT*qrten(k))*rho(k))
          qv(k) = MAX(1.E-10_r8, qv1d(k) + DT*qvten(k))
          nr(k) = MAX(R2, (nr1d(k) + DT*nrten(k))*rho(k))
          temp(k) = t1d(k) + DT*tten(k)
          rho(k) = 0.622_r8*pres(k)/(R*temp(k)*(qv(k)+0.622_r8))
       ENDIF
    ENDDO
    !#if ( WRF_CHEM == 1 )
    DO k = kts, kte
       evapprod(k) = prv_rev(k) - (MIN(zeroD0,prs_sde(k)) + &
            MIN(zeroD0,prg_gde(k)))
       rainprod(k) = prr_wau(k) + prr_rcw(k) + prs_scw(k) + &
            prg_scw(k) + prs_iau(k) + &
            prg_gcw(k) + prs_sci(k) + &
            pri_rci(k)
    ENDDO
    !#endif

    !+---+-----------------------------------------------------------------+
    !..Find max terminal fallspeed (distribution mass-weighted mean
    !.. velocity) and use it to determine if we need to split the timestep
    !.. (var nstep>1).  Either way, only bother to do sedimentation below
    !.. 1st level that contains any sedimenting particles (k=ksed1 on down).
    !.. New in v3.0+ is computing separate for rain, ice, snow, and
    !.. graupel species thus making code faster with credit to J. Schmidt.
    !+---+-----------------------------------------------------------------+
    nstep = 0
    onstep(:) = 1.0_r8
    ksed1(:) = 1
    DO k = kte+1, kts, -1
       vtrk(k) = 0._r8
       vtnrk(k) = 0._r8
       vtik(k) = 0._r8
       vtnik(k) = 0._r8
       vtsk(k) = 0._r8
       vtgk(k) = 0._r8
       vtck(k) = 0._r8
       vtnck(k) = 0._r8
    ENDDO
    DO k = kte, kts, -1
       vtr = 0.0_r8
       rhof(k) = SQRT(RHO_NOT/rho(k))

       IF (rr(k).GT. R1) THEN
          lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
          vtr = rhof(k)*av_r*crg(6)*org3 * lamr**cre(3)                 &
               *((lamr+fv_r)**(-cre(6)))
          vtrk(k) = vtr
          ! First below is technically correct:
          !         vtr = rhof(k)*av_r*crg(5)*org2 * lamr**cre(2)                 &
          !                     *((lamr+fv_r)**(-cre(5)))
          ! Test: make number fall faster (but still slower than mass)
          ! Goal: less prominent size sorting
          vtr = rhof(k)*av_r*crg(7)/crg(12) * lamr**cre(12)             &
               *((lamr+fv_r)**(-cre(7)))
          vtnrk(k) = vtr
       ELSE
          vtrk(k) = vtrk(k+1)
          vtnrk(k) = vtnrk(k+1)
       ENDIF

       IF (MAX(vtrk(k),vtnrk(k)) .GT. 1.E-3_r8) THEN
          ksed1(1) = MAX(ksed1(1), k)
          delta_tp = dzq(k)/(MAX(vtrk(k),vtnrk(k)))
          nstep = MAX(nstep, INT(DT/delta_tp + 1.0_r8))
       ENDIF
    ENDDO
    IF (ksed1(1) .EQ. kte) ksed1(1) = kte-1
    IF (nstep .GT. 0) onstep(1) = 1.0_r8/REAL(nstep)

    !+---+-----------------------------------------------------------------+

    hgt_agl = 0.0_r8
    DO k = kts, kte-1
       IF (rc(k) .GT. R2) ksed1(5) = k
       hgt_agl = hgt_agl + dzq(k)
       IF (hgt_agl .GT. 500.0_r8) GOTO 151
    ENDDO
151 CONTINUE

    DO k = ksed1(5), kts, -1
       vtc = 0.0_r8
       IF (rc(k) .GT. R1 .AND. w1d(k) .LT. 1.E-1_r8) THEN
          nu_c = MIN(15, NINT(1000.E6_r8/nc(k)) + 2)
          lamc = (nc(k)*am_r*ccg(2,nu_c)*ocg1(nu_c)/rc(k))**obmr
          ilamc = 1.0_r8/lamc
          vtc = rhof(k)*av_c*ccg(5,nu_c)*ocg2(nu_c) * ilamc**bv_c
          vtck(k) = vtc
          vtc = rhof(k)*av_c*ccg(4,nu_c)*ocg1(nu_c) * ilamc**bv_c
          vtnck(k) = vtc
       ENDIF
    ENDDO

    !+---+-----------------------------------------------------------------+

    IF (.NOT. iiwarm) THEN

       nstep = 0
       DO k = kte, kts, -1
          vti = 0.0_r8

          IF (ri(k).GT. R1) THEN
             lami = (am_i*cig(2)*oig1*ni(k)/ri(k))**obmi
             ilami = 1.0_r8/lami
             vti = rhof(k)*av_i*cig(3)*oig2 * ilami**bv_i
             vtik(k) = vti
             ! First below is technically correct:
             !          vti = rhof(k)*av_i*cig(4)*oig1 * ilami**bv_i
             ! Goal: less prominent size sorting
             vti = rhof(k)*av_i*cig(6)/cig(7) * ilami**bv_i
             vtnik(k) = vti
          ELSE
             vtik(k) = vtik(k+1)
             vtnik(k) = vtnik(k+1)
          ENDIF

          IF (vtik(k) .GT. 1.E-3_r8) THEN
             ksed1(2) = MAX(ksed1(2), k)
             delta_tp = dzq(k)/vtik(k)
             nstep = MAX(nstep, INT(DT/delta_tp + 1.0_r8))
          ENDIF
       ENDDO
       IF (ksed1(2) .EQ. kte) ksed1(2) = kte-1
       IF (nstep .GT. 0) onstep(2) = 1.0_r8/REAL(nstep)

       !+---+-----------------------------------------------------------------+

       nstep = 0
       DO k = kte, kts, -1
          vts = 0.0_r8

          IF (rs(k).GT. R1) THEN
             xDs = smoc(k) / smob(k)
             Mrat = 1.0_r8/xDs
             ils1 = 1.0_r8/(Mrat*Lam0 + fv_s)
             ils2 = 1.0_r8/(Mrat*Lam1 + fv_s)
             t1_vts = Kap0*csg(4)*ils1**cse(4)
             t2_vts = Kap1*Mrat**mu_s*csg(10)*ils2**cse(10)
             ils1 = 1.0_r8/(Mrat*Lam0)
             ils2 = 1.0_r8/(Mrat*Lam1)
             t3_vts = Kap0*csg(1)*ils1**cse(1)
             t4_vts = Kap1*Mrat**mu_s*csg(7)*ils2**cse(7)
             vts = rhof(k)*av_s * (t1_vts+t2_vts)/(t3_vts+t4_vts)
             IF (temp(k).GT. T_0) THEN
                vtsk(k) = MAX(vts*vts_boost(k), vtrk(k))
             ELSE
                vtsk(k) = vts*vts_boost(k)
             ENDIF
          ELSE
             vtsk(k) = vtsk(k+1)
          ENDIF

          IF (vtsk(k) .GT. 1.E-3_r8) THEN
             ksed1(3) = MAX(ksed1(3), k)
             delta_tp = dzq(k)/vtsk(k)
             nstep = MAX(nstep, INT(DT/delta_tp + 1.0_r8))
          ENDIF
       ENDDO
       IF (ksed1(3) .EQ. kte) ksed1(3) = kte-1
       IF (nstep .GT. 0) onstep(3) = 1.0_r8/REAL(nstep)

       !+---+-----------------------------------------------------------------+

       nstep = 0
       DO k = kte, kts, -1
          vtg = 0.0_r8

          IF (rg(k).GT. R1) THEN
             vtg = rhof(k)*av_g*cgg(6)*ogg3 * ilamg(k)**bv_g
             IF (temp(k).GT. T_0) THEN
                vtgk(k) = MAX(vtg, vtrk(k))
             ELSE
                vtgk(k) = vtg
             ENDIF
          ELSE
             vtgk(k) = vtgk(k+1)
          ENDIF

          IF (vtgk(k) .GT. 1.E-3_r8) THEN
             ksed1(4) = MAX(ksed1(4), k)
             delta_tp = dzq(k)/vtgk(k)
             nstep = MAX(nstep, INT(DT/delta_tp + 1.0_r8))
          ENDIF
       ENDDO
       IF (ksed1(4) .EQ. kte) ksed1(4) = kte-1
       IF (nstep .GT. 0) onstep(4) = 1.0_r8/REAL(nstep)
    ENDIF

    !+---+-----------------------------------------------------------------+
    !..Sedimentation of mixing ratio is the integral of v(D)*m(D)*N(D)*dD,
    !.. whereas neglect m(D) term for number concentration.  Therefore,
    !.. cloud ice has proper differential sedimentation.
    !.. New in v3.0+ is computing separate for rain, ice, snow, and
    !.. graupel species thus making code faster with credit to J. Schmidt.
    !.. Bug fix, 2013Nov01 to tendencies using rho(k+1) correction thanks to
    !.. Eric Skyllingstad.
    !+---+-----------------------------------------------------------------+

    nstep = NINT(1.0_r8/onstep(1))
    DO n = 1, nstep
       DO k = kte, kts, -1
          sed_r(k) = vtrk(k)*rr(k)
          sed_n(k) = vtnrk(k)*nr(k)
       ENDDO
       k = kte
       odzq = 1.0_r8/dzq(k)
       orho = 1.0_r8/rho(k)
       qrten(k) = qrten(k) - sed_r(k)*odzq*onstep(1)*orho
       nrten(k) = nrten(k) - sed_n(k)*odzq*onstep(1)*orho
       rr(k) = MAX(R1, rr(k) - sed_r(k)*odzq*DT*onstep(1))
       nr(k) = MAX(R2, nr(k) - sed_n(k)*odzq*DT*onstep(1))
       DO k = ksed1(1), kts, -1
          odzq = 1.0_r8/dzq(k)
          orho = 1.0_r8/rho(k)
          qrten(k) = qrten(k) + (sed_r(k+1)-sed_r(k))                 &
               *odzq*onstep(1)*orho
          nrten(k) = nrten(k) + (sed_n(k+1)-sed_n(k))                 &
               *odzq*onstep(1)*orho
          rr(k) = MAX(R1, rr(k) + (sed_r(k+1)-sed_r(k)) &
               *odzq*DT*onstep(1))
          nr(k) = MAX(R2, nr(k) + (sed_n(k+1)-sed_n(k)) &
               *odzq*DT*onstep(1))
       ENDDO

       IF (rr(kts).GT.R1*10.0_r8) &
            pptrain = pptrain + sed_r(kts)*DT*onstep(1)
    ENDDO

    !+---+-----------------------------------------------------------------+

    DO k = kte, kts, -1
       sed_c(k) = vtck(k)*rc(k)
       sed_n(k) = vtnck(k)*nc(k)
    ENDDO
    DO k = ksed1(5), kts, -1
       odzq = 1.0_r8/dzq(k)
       orho = 1.0_r8/rho(k)
       qcten(k) = qcten(k) + (sed_c(k+1)-sed_c(k)) *odzq*orho
       ncten(k) = ncten(k) + (sed_n(k+1)-sed_n(k)) *odzq*orho
       rc(k) = MAX(R1, rc(k) + (sed_c(k+1)-sed_c(k)) *odzq*DT)
       nc(k) = MAX(10.0_r8, nc(k) + (sed_n(k+1)-sed_n(k)) *odzq*DT)
    ENDDO

    !+---+-----------------------------------------------------------------+

    nstep = NINT(1.0_r8/onstep(2))
    DO n = 1, nstep
       DO k = kte, kts, -1
          sed_i(k) = vtik(k)*ri(k)
          sed_n(k) = vtnik(k)*ni(k)
       ENDDO
       k = kte
       odzq = 1.0_r8/dzq(k)
       orho = 1.0_r8/rho(k)
       qiten(k) = qiten(k) - sed_i(k)*odzq*onstep(2)*orho
       niten(k) = niten(k) - sed_n(k)*odzq*onstep(2)*orho
       ri(k) = MAX(R1, ri(k) - sed_i(k)*odzq*DT*onstep(2))
       ni(k) = MAX(R2, ni(k) - sed_n(k)*odzq*DT*onstep(2))
       DO k = ksed1(2), kts, -1
          odzq = 1.0_r8/dzq(k)
          orho = 1.0_r8/rho(k)
          qiten(k) = qiten(k) + (sed_i(k+1)-sed_i(k))                 &
               *odzq*onstep(2)*orho
          niten(k) = niten(k) + (sed_n(k+1)-sed_n(k))                 &
               *odzq*onstep(2)*orho
          ri(k) = MAX(R1, ri(k) + (sed_i(k+1)-sed_i(k)) &
               *odzq*DT*onstep(2))
          ni(k) = MAX(R2, ni(k) + (sed_n(k+1)-sed_n(k)) &
               *odzq*DT*onstep(2))
       ENDDO

       IF (ri(kts).GT.R1*10.0_r8) &
            pptice = pptice + sed_i(kts)*DT*onstep(2)
    ENDDO

    !+---+-----------------------------------------------------------------+

    nstep = NINT(1.0_r8/onstep(3))
    DO n = 1, nstep
       DO k = kte, kts, -1
          sed_s(k) = vtsk(k)*rs(k)
       ENDDO
       k = kte
       odzq = 1.0_r8/dzq(k)
       orho = 1.0_r8/rho(k)
       qsten(k) = qsten(k) - sed_s(k)*odzq*onstep(3)*orho
       rs(k) = MAX(R1, rs(k) - sed_s(k)*odzq*DT*onstep(3))
       DO k = ksed1(3), kts, -1
          odzq = 1.0_r8/dzq(k)
          orho = 1.0_r8/rho(k)
          qsten(k) = qsten(k) + (sed_s(k+1)-sed_s(k))                 &
               *odzq*onstep(3)*orho
          rs(k) = MAX(R1, rs(k) + (sed_s(k+1)-sed_s(k)) &
               *odzq*DT*onstep(3))
       ENDDO

       IF (rs(kts).GT.R1*10.0_r8) &
            pptsnow = pptsnow + sed_s(kts)*DT*onstep(3)
    ENDDO

    !+---+-----------------------------------------------------------------+

    nstep = NINT(1.0_r8/onstep(4))
    DO n = 1, nstep
       DO k = kte, kts, -1
          sed_g(k) = vtgk(k)*rg(k)
       ENDDO
       k = kte
       odzq = 1.0_r8/dzq(k)
       orho = 1.0_r8/rho(k)
       qgten(k) = qgten(k) - sed_g(k)*odzq*onstep(4)*orho
       rg(k) = MAX(R1, rg(k) - sed_g(k)*odzq*DT*onstep(4))
       DO k = ksed1(4), kts, -1
          odzq = 1.0_r8/dzq(k)
          orho = 1.0_r8/rho(k)
          qgten(k) = qgten(k) + (sed_g(k+1)-sed_g(k))                 &
               *odzq*onstep(4)*orho
          rg(k) = MAX(R1, rg(k) + (sed_g(k+1)-sed_g(k)) &
               *odzq*DT*onstep(4))
       ENDDO

       IF (rg(kts).GT.R1*10.0_r8) &
            pptgraul = pptgraul + sed_g(kts)*DT*onstep(4)
    ENDDO

    !+---+-----------------------------------------------------------------+
    !.. Instantly melt any cloud ice into cloud water if above 0C and
    !.. instantly freeze any cloud water found below HGFR.
    !+---+-----------------------------------------------------------------+
    IF (.NOT. iiwarm) THEN
       DO k = kts, kte
          xri = MAX(0.0_r8, qi1d(k) + qiten(k)*DT)
          IF ( (temp(k).GT. T_0) .AND. (xri.GT. 0.0_r8) ) THEN
             qcten(k) = qcten(k) + xri*odt
             ncten(k) = ncten(k) + ni1d(k)*odt
             qiten(k) = qiten(k) - xri*odt
             niten(k) = -ni1d(k)*odt
             tten(k) = tten(k) - lfus*ocp(k)*xri*odt*(1-IFDRY)
          ENDIF

          xrc = MAX(0.0_r8, qc1d(k) + qcten(k)*DT)
          IF ( (temp(k).LT. HGFR) .AND. (xrc.GT. 0.0_r8) ) THEN
             lfus2 = lsub - lvap(k)
             xnc = nc1d(k) + ncten(k)*DT
             qiten(k) = qiten(k) + xrc*odt
             niten(k) = niten(k) + xnc*odt
             qcten(k) = qcten(k) - xrc*odt
             ncten(k) = ncten(k) - xnc*odt
             tten(k) = tten(k) + lfus2*ocp(k)*xrc*odt*(1-IFDRY)
          ENDIF
       ENDDO
    ENDIF

    !+---+-----------------------------------------------------------------+
    !.. All tendencies computed, apply and pass back final values to parent.
    !+---+-----------------------------------------------------------------+
    DO k = kts, kte
       t1d(k)  = t1d(k) + tten(k)*DT
       qv1d(k) = MAX(1.E-10_r8, qv1d(k) + qvten(k)*DT)
       qc1d(k) = qc1d(k) + qcten(k)*DT
       nc1d(k) = MAX(2.0_r8/rho(k), nc1d(k) + ncten(k)*DT)
       nwfa1d(k) = MAX(11.1E6_r8/rho(k), MIN(9999.E6_r8/rho(k),             &
            (nwfa1d(k)+nwfaten(k)*DT)))
       nifa1d(k) = MAX(naIN1*0.01_r8, MIN(9999.E6_r8/rho(k),                &
            (nifa1d(k)+nifaten(k)*DT)))

       IF (qc1d(k) .LE. R1) THEN
          qc1d(k) = 0.0_r8
          nc1d(k) = 0.0_r8
       ELSE
          nu_c = MIN(15, NINT(1000.E6_r8/(nc1d(k)*rho(k))) + 2)
          lamc = (am_r*ccg(2,nu_c)*ocg1(nu_c)*nc1d(k)/qc1d(k))**obmr
          xDc = (bm_r + nu_c + 1.0_r8) / lamc
          IF (xDc.LT. D0c) THEN
             lamc = cce(2,nu_c)/D0c
          ELSEIF (xDc.GT. D0r*2.) THEN
             lamc = cce(2,nu_c)/(D0r*2.0_r8)
          ENDIF
          nc1d(k) = MIN(ccg(1,nu_c)*ocg2(nu_c)*qc1d(k)/am_r*lamc**bm_r,&
               DBLE(Nt_c_max)/rho(k))
       ENDIF

       qi1d(k) = qi1d(k) + qiten(k)*DT
       ni1d(k) = MAX(R2/rho(k), ni1d(k) + niten(k)*DT)
       IF (qi1d(k) .LE. R1) THEN
          qi1d(k) = 0.0_r8
          ni1d(k) = 0.0_r8
       ELSE
          lami = (am_i*cig(2)*oig1*ni1d(k)/qi1d(k))**obmi
          ilami = 1.0_r8/lami
          xDi = (bm_i + mu_i + 1.0_r8) * ilami
          IF (xDi.LT. 20.E-6_r8) THEN
             lami = cie(2)/20.E-6_r8
          ELSEIF (xDi.GT. 300.E-6_r8) THEN
             lami = cie(2)/300.E-6_r8
          ENDIF
          ni1d(k) = MIN(cig(1)*oig2*qi1d(k)/am_i*lami**bm_i,           &
               499.e3_r8/rho(k))
       ENDIF
       qr1d(k) = qr1d(k) + qrten(k)*DT
       nr1d(k) = MAX(R2/rho(k), nr1d(k) + nrten(k)*DT)
       IF (qr1d(k) .LE. R1) THEN
          qr1d(k) = 0.0_r8
          nr1d(k) = 0.0_r8
       ELSE
          lamr = (am_r*crg(3)*org2*nr1d(k)/qr1d(k))**obmr
          mvd_r(k) = (3.0_r8 + mu_r + 0.672_r8) / lamr
          IF (mvd_r(k) .GT. 2.5E-3_r8) THEN
             mvd_r(k) = 2.5E-3_r8
          ELSEIF (mvd_r(k) .LT. D0r*0.75_r8) THEN
             mvd_r(k) = D0r*0.75_r8
          ENDIF
          lamr = (3.0_r8 + mu_r + 0.672_r8) / mvd_r(k)
          nr1d(k) = crg(2)*org3*qr1d(k)*lamr**bm_r / am_r
       ENDIF
       qs1d(k) = qs1d(k) + qsten(k)*DT
       IF (qs1d(k) .LE. R1) qs1d(k) = 0.0_r8
       qg1d(k) = qg1d(k) + qgten(k)*DT
       IF (qg1d(k) .LE. R1) qg1d(k) = 0.0_r8
    ENDDO

  END SUBROUTINE mp_thompson

  !+---+-----------------------------------------------------------------+
  !+---+-----------------------------------------------------------------+
  !ctrlL
  !+---+-----------------------------------------------------------------+
  !..Creation of the lookup tables and support functions found below here.
  !+---+-----------------------------------------------------------------+
  !..Rain collecting graupel (and inverse).  Explicit CE integration.
  !+---+-----------------------------------------------------------------+

  SUBROUTINE qr_acr_qg(path_in)

    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN   ) ::path_in
    !..Local variables
    INTEGER:: i, j, k, m, n, n2
    INTEGER:: km, km_s, km_e
    REAL(KIND=r8), DIMENSION(nbg):: vg, N_g
    REAL(KIND=r8), DIMENSION(nbr):: vr, N_r
    REAL(KIND=r8):: N0_r, N0_g, lam_exp, lamg, lamr
    REAL(KIND=r8):: massg, massr, dvg, dvr, t1, t2, z1, z2, y1, y2
    LOGICAL :: force_read_thompson=.FALSE.
    LOGICAL :: write_thompson_tables=.TRUE.
    LOGICAL :: wrf_dm_on_monitor=.TRUE.
    LOGICAL lexist,lopen,opened
    INTEGER good,iunit_mp_th1
    !      LOGICAL, EXTERNAL :: wrf_dm_on_monitor
    iunit_mp_th1 = -1
    IF ( wrf_dm_on_monitor ) THEN
       DO i = 7,99
          INQUIRE ( i , OPENED = opened )
          IF ( .NOT. opened ) THEN
             iunit_mp_th1 = i
             GOTO 2010
          ENDIF
       ENDDO
2010   CONTINUE
    ENDIF
    !+---+

    ! CALL nl_get_force_read_thompson(1,force_read_thompson)
    ! CALL nl_get_write_thompson_tables(1,write_thompson_tables)

    good = 0
    IF ( wrf_dm_on_monitor ) THEN
       INQUIRE(FILE=TRIM(path_in)//'/'//"qr_acr_qg.dat",EXIST=lexist)
       IF ( lexist ) THEN
          CALL wrf_message("ThompMP: read qr_acr_qg.dat stead of computing")
          OPEN(iunit_mp_th1,file=TRIM(path_in)//'/'//"qr_acr_qg.dat",form="unformatted",err=1234)
          READ(iunit_mp_th1,err=1234) tcg_racg
          READ(iunit_mp_th1,err=1234) tmr_racg
          READ(iunit_mp_th1,err=1234) tcr_gacr
          READ(iunit_mp_th1,err=1234) tmg_gacr
          READ(iunit_mp_th1,err=1234) tnr_racg
          READ(iunit_mp_th1,err=1234) tnr_gacr
          good = 1
          CLOSE(iunit_mp_th1)
1234      CONTINUE
          IF ( good .NE. 1 ) THEN
             INQUIRE(iunit_mp_th1,opened=lopen)
             IF (lopen) THEN
                IF( force_read_thompson ) THEN
                   CALL wrf_error_fatal("Error reading qr_acr_qg.dat. Aborting because force_read_thompson is .true.")
                ENDIF
                CLOSE(iunit_mp_th1)
             ELSE
                IF( force_read_thompson ) THEN
                   CALL wrf_error_fatal("Error opening qr_acr_qg.dat. Aborting because force_read_thompson is .true.")
                ENDIF
             ENDIF
          ENDIF
       ELSE
          IF( force_read_thompson ) THEN
             CALL wrf_error_fatal("Non-existent qr_acr_qg.dat. Aborting because force_read_thompson is .true.")
          ENDIF
       ENDIF
    ENDIF
    !#if defined(DM_PARALLEL) && !defined(STUBMPI)
    !      CALL wrf_dm_bcast_integer(good,1)
    !#endif

    IF ( good .EQ. 1 ) THEN
       !#if defined(DM_PARALLEL) && !defined(STUBMPI)
       !        CALL wrf_dm_bcast_double(tcg_racg,SIZE(tcg_racg))
       !        CALL wrf_dm_bcast_double(tmr_racg,SIZE(tmr_racg))
       !        CALL wrf_dm_bcast_double(tcr_gacr,SIZE(tcr_gacr))
       !        CALL wrf_dm_bcast_double(tmg_gacr,SIZE(tmg_gacr))
       !        CALL wrf_dm_bcast_double(tnr_racg,SIZE(tnr_racg))
       !        CALL wrf_dm_bcast_double(tnr_gacr,SIZE(tnr_gacr))
       !#endif
    ELSE
       CALL wrf_message("ThompMP: computing qr_acr_qg")
       DO n2 = 1, nbr
          !        vr(n2) = av_r*Dr(n2)**bv_r * DEXP(-fv_r*Dr(n2))
          vr(n2) = -0.1021_r8 + 4.932E3_r8*Dr(n2) - 0.9551E6_r8*Dr(n2)*Dr(n2)     &
               + 0.07934E9_r8*Dr(n2)*Dr(n2)*Dr(n2)                          &
               - 0.002362E12_r8*Dr(n2)*Dr(n2)*Dr(n2)*Dr(n2)
       ENDDO
       DO n = 1, nbg
          vg(n) = av_g*Dg(n)**bv_g
       ENDDO

       !..Note values returned from wrf_dm_decomp1d are zero-based, add 1 for
       !.. fortran indices.  J. Michalakes, 2009Oct30.

       !#if ( defined( DM_PARALLEL ) && ( ! defined( STUBMPI ) ) )
       !        CALL wrf_dm_decomp1d ( ntb_r*ntb_r1, km_s, km_e )
       !#else
       km_s = 0
       km_e = ntb_r*ntb_r1 - 1
       !#endif

       DO km = km_s, km_e
          m = km / ntb_r1 + 1
          k = MOD( km , ntb_r1 ) + 1

          lam_exp = (N0r_exp(k)*am_r*crg(1)/r_r(m))**ore1
          lamr = lam_exp * (crg(3)*org2*org1)**obmr
          N0_r = N0r_exp(k)/(crg(2)*lam_exp) * lamr**cre(2)
          DO n2 = 1, nbr
             N_r(n2) = N0_r*Dr(n2)**mu_r *EXP(-lamr*Dr(n2))*dtr(n2)
          ENDDO

          DO j = 1, ntb_g
             DO i = 1, ntb_g1
                lam_exp = (N0g_exp(i)*am_g*cgg(1)/r_g(j))**oge1
                lamg = lam_exp * (cgg(3)*ogg2*ogg1)**obmg
                N0_g = N0g_exp(i)/(cgg(2)*lam_exp) * lamg**cge(2)
                DO n = 1, nbg
                   N_g(n) = N0_g*Dg(n)**mu_g * EXP(-lamg*Dg(n))*dtg(n)
                ENDDO

                t1 = 0.0e0_r8
                t2 = 0.0e0_r8
                z1 = 0.0e0_r8
                z2 = 0.0e0_r8
                y1 = 0.0e0_r8
                y2 = 0.0e0_r8
                DO n2 = 1, nbr
                   massr = am_r * Dr(n2)**bm_r
                   DO n = 1, nbg
                      massg = am_g * Dg(n)**bm_g

                      dvg = 0.5e0_r8*((vr(n2) - vg(n)) + ABS(vr(n2)-vg(n)))
                      dvr = 0.5e0_r8*((vg(n) - vr(n2)) + ABS(vg(n)-vr(n2)))

                      t1 = t1+ PI*0.25_r8*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                           *dvg*massg * N_g(n)* N_r(n2)
                      z1 = z1+ PI*0.25_r8*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                           *dvg*massr * N_g(n)* N_r(n2)
                      y1 = y1+ PI*0.25_r8*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                           *dvg       * N_g(n)* N_r(n2)

                      t2 = t2+ PI*0.25_r8*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                           *dvr*massr * N_g(n)* N_r(n2)
                      y2 = y2+ PI*0.25_r8*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                           *dvr       * N_g(n)* N_r(n2)
                      z2 = z2+ PI*0.25_r8*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                           *dvr*massg * N_g(n)* N_r(n2)
                   ENDDO
!97                 CONTINUE
                ENDDO
                tcg_racg(i,j,k,m) = t1
                tmr_racg(i,j,k,m) = MIN(z1, r_r(m)*1.0e0_r8)
                tcr_gacr(i,j,k,m) = t2
                tmg_gacr(i,j,k,m) = z2
                tnr_racg(i,j,k,m) = y1
                tnr_gacr(i,j,k,m) = y2
             ENDDO
          ENDDO
       ENDDO

       !..Note wrf_dm_gatherv expects zero-based km_s, km_e (J. Michalakes, 2009Oct30).

       !#if ( defined( DM_PARALLEL ) && ( ! defined( STUBMPI ) ) )
       !        CALL wrf_dm_gatherv(tcg_racg, ntb_g*ntb_g1, km_s, km_e, R8SIZE)
       !        CALL wrf_dm_gatherv(tmr_racg, ntb_g*ntb_g1, km_s, km_e, R8SIZE)
       !        CALL wrf_dm_gatherv(tcr_gacr, ntb_g*ntb_g1, km_s, km_e, R8SIZE)
       !        CALL wrf_dm_gatherv(tmg_gacr, ntb_g*ntb_g1, km_s, km_e, R8SIZE)
       !        CALL wrf_dm_gatherv(tnr_racg, ntb_g*ntb_g1, km_s, km_e, R8SIZE)
       !        CALL wrf_dm_gatherv(tnr_gacr, ntb_g*ntb_g1, km_s, km_e, R8SIZE)
       !#endif
       IF ( write_thompson_tables .AND. wrf_dm_on_monitor ) THEN
          CALL wrf_message("Writing qr_acr_qg.dat in Thompson MP init")
          OPEN(iunit_mp_th1,file=TRIM(path_in)//'/'//"qr_acr_qg.dat",form="unformatted",err=9234)
          WRITE(iunit_mp_th1,err=9234) tcg_racg
          WRITE(iunit_mp_th1,err=9234) tmr_racg
          WRITE(iunit_mp_th1,err=9234) tcr_gacr
          WRITE(iunit_mp_th1,err=9234) tmg_gacr
          WRITE(iunit_mp_th1,err=9234) tnr_racg
          WRITE(iunit_mp_th1,err=9234) tnr_gacr
          CLOSE(iunit_mp_th1)
          RETURN    ! ----- RETURN
9234      CONTINUE!
          CALL wrf_error_fatal("Error writing qr_acr_qg.dat")!
       ENDIF

    ENDIF

  END SUBROUTINE qr_acr_qg

  !+---+-----------------------------------------------------------------+
  !ctrlL
  !+---+-----------------------------------------------------------------+
  !..Rain collecting snow (and inverse).  Explicit CE integration.
  !+---+-----------------------------------------------------------------+

  SUBROUTINE qr_acr_qs(path_in)

    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN   ) :: path_in
    !..Local variables
    INTEGER:: i, j, k, m, n, n2
    INTEGER:: km, km_s, km_e
    REAL(KIND=r8), DIMENSION(nbr):: vr, D1, N_r
    REAL(KIND=r8), DIMENSION(nbs):: vs, N_s
    REAL(KIND=r8):: loga_, a_, b_, second2, M0, M2, M3, Mrat, oM3
    REAL(KIND=r8):: N0_r, lam_exp, lamr, slam1, slam2
    REAL(KIND=r8):: dvs, dvr, masss, massr
    REAL(KIND=r8):: t1, t2, t3, t4, z1, z2, z3, z4
    REAL(KIND=r8):: y1, y2, y3, y4
    LOGICAL :: force_read_thompson=.FALSE.
    LOGICAL :: write_thompson_tables=.TRUE.
    LOGICAL :: wrf_dm_on_monitor=.TRUE.
    LOGICAL lexist,lopen,opened
    INTEGER good,iunit_mp_th1

    iunit_mp_th1 = -1
    IF ( wrf_dm_on_monitor ) THEN
       DO i = 7,99
          INQUIRE ( i , OPENED = opened )
          IF ( .NOT. opened ) THEN
             iunit_mp_th1 = i
             GOTO 2010
          ENDIF
       ENDDO
2010   CONTINUE
    ENDIF
    !+---+

    ! CALL nl_get_force_read_thompson(1,force_read_thompson)
    ! CALL nl_get_write_thompson_tables(1,write_thompson_tables)

    good = 0
    IF ( wrf_dm_on_monitor ) THEN
       INQUIRE(FILE=TRIM(path_in)//'/'//"qr_acr_qs.dat",EXIST=lexist)
       IF ( lexist ) THEN
          CALL wrf_message("ThompMP: read qr_acr_qs.dat instead of computing")
          OPEN(iunit_mp_th1,file=TRIM(path_in)//'/'//"qr_acr_qs.dat",form="unformatted",err=1234)
          READ(iunit_mp_th1,err=1234)tcs_racs1
          READ(iunit_mp_th1,err=1234)tmr_racs1
          READ(iunit_mp_th1,err=1234)tcs_racs2
          READ(iunit_mp_th1,err=1234)tmr_racs2
          READ(iunit_mp_th1,err=1234)tcr_sacr1
          READ(iunit_mp_th1,err=1234)tms_sacr1
          READ(iunit_mp_th1,err=1234)tcr_sacr2
          READ(iunit_mp_th1,err=1234)tms_sacr2
          READ(iunit_mp_th1,err=1234)tnr_racs1
          READ(iunit_mp_th1,err=1234)tnr_racs2
          READ(iunit_mp_th1,err=1234)tnr_sacr1
          READ(iunit_mp_th1,err=1234)tnr_sacr2
          good = 1
          CLOSE(iunit_mp_th1)
1234      CONTINUE
          IF ( good .NE. 1 ) THEN
             INQUIRE(iunit_mp_th1,opened=lopen)
             IF (lopen) THEN
                IF( force_read_thompson ) THEN
                   CALL wrf_error_fatal("Error reading qr_acr_qs.dat. Aborting because force_read_thompson is .true.")
                ENDIF
                CLOSE(iunit_mp_th1)
             ELSE
                IF( force_read_thompson ) THEN
                   CALL wrf_error_fatal("Error opening qr_acr_qs.dat. Aborting because force_read_thompson is .true.")
                ENDIF
             ENDIF
          ENDIF
       ELSE
          IF( force_read_thompson ) THEN
             CALL wrf_error_fatal("Non-existent qr_acr_qs.dat. Aborting because force_read_thompson is .true.")
          ENDIF
       ENDIF
    ENDIF
    !#if defined(DM_PARALLEL) && !defined(STUBMPI)
    !      CALL wrf_dm_bcast_integer(good,1)
    !#endif

    IF ( good .EQ. 1 ) THEN
       !#if defined(DM_PARALLEL) && !defined(STUBMPI)
       !        CALL wrf_dm_bcast_double(tcs_racs1,SIZE(tcs_racs1))
       !        CALL wrf_dm_bcast_double(tmr_racs1,SIZE(tmr_racs1))
       !        CALL wrf_dm_bcast_double(tcs_racs2,SIZE(tcs_racs2))
       !        CALL wrf_dm_bcast_double(tmr_racs2,SIZE(tmr_racs2))
       !        CALL wrf_dm_bcast_double(tcr_sacr1,SIZE(tcr_sacr1))
       !        CALL wrf_dm_bcast_double(tms_sacr1,SIZE(tms_sacr1))
       !        CALL wrf_dm_bcast_double(tcr_sacr2,SIZE(tcr_sacr2))
       !        CALL wrf_dm_bcast_double(tms_sacr2,SIZE(tms_sacr2))
       !        CALL wrf_dm_bcast_double(tnr_racs1,SIZE(tnr_racs1))
       !        CALL wrf_dm_bcast_double(tnr_racs2,SIZE(tnr_racs2))
       !        CALL wrf_dm_bcast_double(tnr_sacr1,SIZE(tnr_sacr1))
       !        CALL wrf_dm_bcast_double(tnr_sacr2,SIZE(tnr_sacr2))
       !#endif
    ELSE
       CALL wrf_message("ThompMP: computing qr_acr_qs")
       DO n2 = 1, nbr
          !        vr(n2) = av_r*Dr(n2)**bv_r * DEXP(-fv_r*Dr(n2))
          vr(n2) = -0.1021_r8 + 4.932E3_r8*Dr(n2) - 0.9551E6_r8*Dr(n2)*Dr(n2)     &
               + 0.07934E9_r8*Dr(n2)*Dr(n2)*Dr(n2)                          &
               - 0.002362E12_r8*Dr(n2)*Dr(n2)*Dr(n2)*Dr(n2)
          D1(n2) = (vr(n2)/av_s)**(1.0_r8/bv_s)
       ENDDO
       DO n = 1, nbs
          vs(n) = 1.5_r8*av_s*Ds(n)**bv_s * DEXP(-fv_s*Ds(n))
       ENDDO

       !..Note values returned from wrf_dm_decomp1d are zero-based, add 1 for
       !.. fortran indices.  J. Michalakes, 2009Oct30.

       !#if ( defined( DM_PARALLEL ) && ( ! defined( STUBMPI ) ) )
       !        CALL wrf_dm_decomp1d ( ntb_r*ntb_r1, km_s, km_e )
       !#else
       km_s = 0
       km_e = ntb_r*ntb_r1 - 1
       !#endif

       DO km = km_s, km_e
          m = km / ntb_r1 + 1
          k = MOD( km , ntb_r1 ) + 1

          lam_exp = (N0r_exp(k)*am_r*crg(1)/r_r(m))**ore1
          lamr = lam_exp * (crg(3)*org2*org1)**obmr
          N0_r = N0r_exp(k)/(crg(2)*lam_exp) * lamr**cre(2)
          DO n2 = 1, nbr
             N_r(n2) = N0_r*Dr(n2)**mu_r * DEXP(-lamr*Dr(n2))*dtr(n2)
          ENDDO

          DO j = 1, ntb_t
             DO i = 1, ntb_s

                !..From the bm_s moment, compute plus one moment.  If we are not
                !.. using bm_s=2, then we must transform to the pure 2nd moment
                !.. (variable called "second2") and then to the bm_s+1 moment.

                M2 = r_s(i)*oams *1.0e0_r8
                IF (bm_s.GT.2.0_r8-1.E-3_r8 .AND. bm_s.LT.2.0_r8+1.E-3_r8) THEN
                   loga_ = sa(1) + sa(2)*Tc(j) + sa(3)*bm_s &
                        + sa(4)*Tc(j)*bm_s + sa(5)*Tc(j)*Tc(j) &
                        + sa(6)*bm_s*bm_s + sa(7)*Tc(j)*Tc(j)*bm_s &
                        + sa(8)*Tc(j)*bm_s*bm_s + sa(9)*Tc(j)*Tc(j)*Tc(j) &
                        + sa(10)*bm_s*bm_s*bm_s
                   a_ = 10.0_r8**loga_
                   b_ = sb(1) + sb(2)*Tc(j) + sb(3)*bm_s &
                        + sb(4)*Tc(j)*bm_s + sb(5)*Tc(j)*Tc(j) &
                        + sb(6)*bm_s*bm_s + sb(7)*Tc(j)*Tc(j)*bm_s &
                        + sb(8)*Tc(j)*bm_s*bm_s + sb(9)*Tc(j)*Tc(j)*Tc(j) &
                        + sb(10)*bm_s*bm_s*bm_s
                   second2 = (M2/a_)**(1.0_r8/b_)
                ELSE
                   second2 = M2
                ENDIF

                loga_ = sa(1) + sa(2)*Tc(j) + sa(3)*cse(1) &
                     + sa(4)*Tc(j)*cse(1) + sa(5)*Tc(j)*Tc(j) &
                     + sa(6)*cse(1)*cse(1) + sa(7)*Tc(j)*Tc(j)*cse(1) &
                     + sa(8)*Tc(j)*cse(1)*cse(1) + sa(9)*Tc(j)*Tc(j)*Tc(j) &
                     + sa(10)*cse(1)*cse(1)*cse(1)
                a_ = 10.0_r8**loga_
                b_ = sb(1)+sb(2)*Tc(j)+sb(3)*cse(1) + sb(4)*Tc(j)*cse(1) &
                     + sb(5)*Tc(j)*Tc(j) + sb(6)*cse(1)*cse(1) &
                     + sb(7)*Tc(j)*Tc(j)*cse(1) + sb(8)*Tc(j)*cse(1)*cse(1) &
                     + sb(9)*Tc(j)*Tc(j)*Tc(j)+sb(10)*cse(1)*cse(1)*cse(1)
                M3 = a_ * second2**b_

                oM3 = 1.0_r8/M3
                Mrat = M2*(M2*oM3)*(M2*oM3)*(M2*oM3)
                M0   = (M2*oM3)**mu_s
                slam1 = M2 * oM3 * Lam0
                slam2 = M2 * oM3 * Lam1

                DO n = 1, nbs
                   N_s(n) = Mrat*(Kap0*DEXP(-slam1*Ds(n)) &
                        + Kap1*M0*Ds(n)**mu_s * DEXP(-slam2*Ds(n)))*dts(n)
                ENDDO

                t1 = 0.0e0_r8
                t2 = 0.0e0_r8
                t3 = 0.0e0_r8
                t4 = 0.0e0_r8
                z1 = 0.0e0_r8
                z2 = 0.0e0_r8
                z3 = 0.0e0_r8
                z4 = 0.0e0_r8
                y1 = 0.0e0_r8
                y2 = 0.0e0_r8
                y3 = 0.0e0_r8
                y4 = 0.0e0_r8
                DO n2 = 1, nbr
                   massr = am_r * Dr(n2)**bm_r
                   DO n = 1, nbs
                      masss = am_s * Ds(n)**bm_s

                      dvs = 0.5e0_r8*((vr(n2) - vs(n)) + DABS(vr(n2)-vs(n)))
                      dvr = 0.5e0_r8*((vs(n) - vr(n2)) + DABS(vs(n)-vr(n2)))

                      IF (massr .GT. 1.5_r8*masss) THEN
                         t1 = t1+ PI*0.25_r8*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                              *dvs*masss * N_s(n)* N_r(n2)
                         z1 = z1+ PI*0.25_r8*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                              *dvs*massr * N_s(n)* N_r(n2)
                         y1 = y1+ PI*0.25_r8*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                              *dvs       * N_s(n)* N_r(n2)
                      ELSE
                         t3 = t3+ PI*0.25_r8*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                              *dvs*masss * N_s(n)* N_r(n2)
                         z3 = z3+ PI*0.25_r8*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                              *dvs*massr * N_s(n)* N_r(n2)
                         y3 = y3+ PI*0.25_r8*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                              *dvs       * N_s(n)* N_r(n2)
                      ENDIF

                      IF (massr .GT. 1.5_r8*masss) THEN
                         t2 = t2+ PI*0.25_r8*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                              *dvr*massr * N_s(n)* N_r(n2)
                         y2 = y2+ PI*0.25_r8*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                              *dvr       * N_s(n)* N_r(n2)
                         z2 = z2+ PI*0.25_r8*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                              *dvr*masss * N_s(n)* N_r(n2)
                      ELSE
                         t4 = t4+ PI*0.25_r8*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                              *dvr*massr * N_s(n)* N_r(n2)
                         y4 = y4+ PI*0.25_r8*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                              *dvr       * N_s(n)* N_r(n2)
                         z4 = z4+ PI*0.25_r8*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                              *dvr*masss * N_s(n)* N_r(n2)
                      ENDIF

                   ENDDO
                ENDDO
                tcs_racs1(i,j,k,m) = t1
                tmr_racs1(i,j,k,m) = DMIN1(z1, r_r(m)*1.0e0_r8)
                tcs_racs2(i,j,k,m) = t3
                tmr_racs2(i,j,k,m) = z3
                tcr_sacr1(i,j,k,m) = t2
                tms_sacr1(i,j,k,m) = z2
                tcr_sacr2(i,j,k,m) = t4
                tms_sacr2(i,j,k,m) = z4
                tnr_racs1(i,j,k,m) = y1
                tnr_racs2(i,j,k,m) = y3
                tnr_sacr1(i,j,k,m) = y2
                tnr_sacr2(i,j,k,m) = y4
             ENDDO
          ENDDO
       ENDDO

       !..Note wrf_dm_gatherv expects zero-based km_s, km_e (J. Michalakes, 2009Oct30).

       !#if ( defined( DM_PARALLEL ) && ( ! defined( STUBMPI ) ) )
       !        CALL wrf_dm_gatherv(tcs_racs1, ntb_s*ntb_t, km_s, km_e, R8SIZE)
       !        CALL wrf_dm_gatherv(tmr_racs1, ntb_s*ntb_t, km_s, km_e, R8SIZE)
       !        CALL wrf_dm_gatherv(tcs_racs2, ntb_s*ntb_t, km_s, km_e, R8SIZE)
       !        CALL wrf_dm_gatherv(tmr_racs2, ntb_s*ntb_t, km_s, km_e, R8SIZE)
       !        CALL wrf_dm_gatherv(tcr_sacr1, ntb_s*ntb_t, km_s, km_e, R8SIZE)
       !        CALL wrf_dm_gatherv(tms_sacr1, ntb_s*ntb_t, km_s, km_e, R8SIZE)
       !        CALL wrf_dm_gatherv(tcr_sacr2, ntb_s*ntb_t, km_s, km_e, R8SIZE)
       !        CALL wrf_dm_gatherv(tms_sacr2, ntb_s*ntb_t, km_s, km_e, R8SIZE)
       !        CALL wrf_dm_gatherv(tnr_racs1, ntb_s*ntb_t, km_s, km_e, R8SIZE)
       !        CALL wrf_dm_gatherv(tnr_racs2, ntb_s*ntb_t, km_s, km_e, R8SIZE)
       !        CALL wrf_dm_gatherv(tnr_sacr1, ntb_s*ntb_t, km_s, km_e, R8SIZE)
       !        CALL wrf_dm_gatherv(tnr_sacr2, ntb_s*ntb_t, km_s, km_e, R8SIZE)
       !#endif

       IF ( write_thompson_tables .AND. wrf_dm_on_monitor ) THEN
          CALL wrf_message("Writing qr_acr_qs.dat in Thompson MP init")
          OPEN(iunit_mp_th1,file=TRIM(path_in)//'/'//"qr_acr_qs.dat",form="unformatted",err=9234)
          WRITE(iunit_mp_th1,err=9234)tcs_racs1
          WRITE(iunit_mp_th1,err=9234)tmr_racs1
          WRITE(iunit_mp_th1,err=9234)tcs_racs2
          WRITE(iunit_mp_th1,err=9234)tmr_racs2
          WRITE(iunit_mp_th1,err=9234)tcr_sacr1
          WRITE(iunit_mp_th1,err=9234)tms_sacr1
          WRITE(iunit_mp_th1,err=9234)tcr_sacr2
          WRITE(iunit_mp_th1,err=9234)tms_sacr2
          WRITE(iunit_mp_th1,err=9234)tnr_racs1
          WRITE(iunit_mp_th1,err=9234)tnr_racs2
          WRITE(iunit_mp_th1,err=9234)tnr_sacr1
          WRITE(iunit_mp_th1,err=9234)tnr_sacr2
          CLOSE(iunit_mp_th1)
          RETURN    ! ----- RETURN
9234      CONTINUE
          CALL wrf_error_fatal("Error writing qr_acr_qs.dat")
       ENDIF
    ENDIF

  END SUBROUTINE qr_acr_qs
  !+---+-----------------------------------------------------------------+
  !ctrlL
  !+---+-----------------------------------------------------------------+
  !..This is a literal adaptation of Bigg (1954) probability of drops of
  !..a particular volume freezing.  Given this probability, simply freeze
  !..the proportion of drops summing their masses.
  !+---+-----------------------------------------------------------------+

  SUBROUTINE freezeH2O(path_in)

    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN   ) :: path_in
    !..Local variables
    INTEGER:: i, j, k, m, n, n2
    REAL(KIND=r8), DIMENSION(nbr):: N_r, massr
    REAL(KIND=r8), DIMENSION(nbc):: N_c, massc
    REAL(KIND=r8):: sum1, sum2, sumn1, sumn2, &
         prob, vol, Texp, orho_w, &
         lam_exp, lamr, N0_r, lamc, N0_c!, y
    INTEGER:: nu_c
    REAL(KIND=r8) :: T_adjust
    LOGICAL :: force_read_thompson=.FALSE.
    LOGICAL :: write_thompson_tables=.TRUE.
    LOGICAL :: wrf_dm_on_monitor=.TRUE.

    LOGICAL lexist,lopen,opened
    INTEGER good,iunit_mp_th1

    iunit_mp_th1 = -1
    IF ( wrf_dm_on_monitor ) THEN
       DO i = 7,99
          INQUIRE ( i , OPENED = opened )
          IF ( .NOT. opened ) THEN
             iunit_mp_th1 = i
             GOTO 2010
          ENDIF
       ENDDO
2010   CONTINUE
    ENDIF
    !+---+
    !CALL nl_get_force_read_thompson(1,force_read_thompson)
    !CALL nl_get_write_thompson_tables(1,write_thompson_tables)

    good = 0
    IF ( wrf_dm_on_monitor) THEN
       INQUIRE(FILE=TRIM(path_in)//'/'//"freezeH2O.dat",EXIST=lexist)
       IF ( lexist ) THEN
          CALL wrf_message("ThompMP: read freezeH2O.dat stead of computing")
          OPEN(iunit_mp_th1,file=TRIM(path_in)//'/'//"freezeH2O.dat",form="unformatted",err=1234)
          READ(iunit_mp_th1,err=1234)tpi_qrfz
          READ(iunit_mp_th1,err=1234)tni_qrfz
          READ(iunit_mp_th1,err=1234)tpg_qrfz
          READ(iunit_mp_th1,err=1234)tnr_qrfz
          READ(iunit_mp_th1,err=1234)tpi_qcfz
          READ(iunit_mp_th1,err=1234)tni_qcfz
          good = 1
          CLOSE(iunit_mp_th1)
1234      CONTINUE
          IF ( good .NE. 1 ) THEN
             INQUIRE(iunit_mp_th1,opened=lopen)
             IF (lopen) THEN
                IF( force_read_thompson ) THEN
                   CALL wrf_error_fatal("Error reading freezeH2O.dat. Aborting because force_read_thompson is .true.")
                ENDIF
                CLOSE(iunit_mp_th1)
             ELSE
                IF( force_read_thompson ) THEN
                   CALL wrf_error_fatal("Error opening freezeH2O.dat. Aborting because force_read_thompson is .true.")
                ENDIF
             ENDIF
          ENDIF
       ELSE
          IF( force_read_thompson ) THEN
             CALL wrf_error_fatal("Non-existent freezeH2O.dat. Aborting because force_read_thompson is .true.")
          ENDIF
       ENDIF
    ENDIF

    !#if defined(DM_PARALLEL) && !defined(STUBMPI)
    !      CALL wrf_dm_bcast_integer(good,1)
    !#endif

    IF ( good .EQ. 1 ) THEN
       !#if defined(DM_PARALLEL) && !defined(STUBMPI)
       !        CALL wrf_dm_bcast_double(tpi_qrfz,SIZE(tpi_qrfz))
       !        CALL wrf_dm_bcast_double(tni_qrfz,SIZE(tni_qrfz))
       !        CALL wrf_dm_bcast_double(tpg_qrfz,SIZE(tpg_qrfz))
       !        CALL wrf_dm_bcast_double(tnr_qrfz,SIZE(tnr_qrfz))
       !        CALL wrf_dm_bcast_double(tpi_qcfz,SIZE(tpi_qcfz))
       !        CALL wrf_dm_bcast_double(tni_qcfz,SIZE(tni_qcfz))
       !#endif
    ELSE
       CALL wrf_message("ThompMP: computing freezeH2O")

       orho_w = 1.0_r8/rho_w

       DO n2 = 1, nbr
          massr(n2) = am_r*Dr(n2)**bm_r
       ENDDO
       DO n = 1, nbc
          massc(n) = am_r*Dc(n)**bm_r
       ENDDO

       !..Freeze water (smallest drops become cloud ice, otherwise graupel).
       DO m = 1, ntb_IN
          T_adjust = MAX(-3.0_r8, MIN(3.0_r8 - log10(Nt_IN(m)), 3.0_r8))
          DO k = 1, 45
             !         print*, ' Freezing water for temp = ', -k
             Texp = DEXP( DFLOAT(k) - T_adjust*1.0e0_r8 ) - 1.0e0_r8
             DO j = 1, ntb_r1
                DO i = 1, ntb_r
                   lam_exp = (N0r_exp(j)*am_r*crg(1)/r_r(i))**ore1
                   lamr = lam_exp * (crg(3)*org2*org1)**obmr
                   N0_r = N0r_exp(j)/(crg(2)*lam_exp) * lamr**cre(2)
                   sum1 = 0.0e0_r8
                   sum2 = 0.0e0_r8
                   sumn1 = 0.0e0_r8
                   sumn2 = 0.0e0_r8
                   DO n2 = nbr, 1, -1
                      N_r(n2) = N0_r*Dr(n2)**mu_r*DEXP(-lamr*Dr(n2))*dtr(n2)
                      vol = massr(n2)*orho_w
                      prob = 1.0e0_r8 - DEXP(-120.0e0_r8*vol*5.2e-4_r8 * Texp)
                      IF (massr(n2) .LT. xm0g) THEN
                         sumn1 = sumn1 + prob*N_r(n2)
                         sum1 = sum1 + prob*N_r(n2)*massr(n2)
                      ELSE
                         sumn2 = sumn2 + prob*N_r(n2)
                         sum2 = sum2 + prob*N_r(n2)*massr(n2)
                      ENDIF
                      IF ((sum1+sum2).GE.r_r(i)) EXIT
                   ENDDO
                   tpi_qrfz(i,j,k,m) = sum1
                   tni_qrfz(i,j,k,m) = sumn1
                   tpg_qrfz(i,j,k,m) = sum2
                   tnr_qrfz(i,j,k,m) = sumn2
                ENDDO
             ENDDO

             DO j = 1, nbc
                nu_c = MIN(15, NINT(1000.E6_r8/t_Nc(j)) + 2)
                DO i = 1, ntb_c
                   lamc = (t_Nc(j)*am_r* ccg(2,nu_c) * ocg1(nu_c) / r_c(i))**obmr
                   N0_c = t_Nc(j)*ocg1(nu_c) * lamc**cce(1,nu_c)
                   sum1 = 0.0e0_r8
                   sumn2 = 0.0e0_r8
                   DO n = nbc, 1, -1
                      vol = massc(n)*orho_w
                      prob = 1.0e0_r8 - DEXP(-120.0e0_r8*vol*5.2e-4_r8 * Texp)
                      N_c(n) = N0_c*Dc(n)**nu_c*EXP(-lamc*Dc(n))*dtc(n)
                      sumn2 = MIN(t_Nc(j), sumn2 + prob*N_c(n))
                      sum1 = sum1 + prob*N_c(n)*massc(n)
                      IF (sum1 .GE. r_c(i)) EXIT
                   ENDDO
                   tpi_qcfz(i,j,k,m) = sum1
                   tni_qcfz(i,j,k,m) = sumn2
                ENDDO
             ENDDO
          ENDDO
       ENDDO

       IF ( write_thompson_tables .AND. wrf_dm_on_monitor ) THEN
          CALL wrf_message("Writing freezeH2O.dat in Thompson MP init")
          OPEN(iunit_mp_th1,file=TRIM(path_in)//'/'//"freezeH2O.dat",form="unformatted",err=9234)
          WRITE(iunit_mp_th1,err=9234)tpi_qrfz
          WRITE(iunit_mp_th1,err=9234)tni_qrfz
          WRITE(iunit_mp_th1,err=9234)tpg_qrfz
          WRITE(iunit_mp_th1,err=9234)tnr_qrfz
          WRITE(iunit_mp_th1,err=9234)tpi_qcfz
          WRITE(iunit_mp_th1,err=9234)tni_qcfz
          CLOSE(iunit_mp_th1)
          RETURN    ! ----- RETURN
9234      CONTINUE
          CALL wrf_error_fatal("Error writing freezeH2O.dat")
       ENDIF
    ENDIF

  END SUBROUTINE freezeH2O
  !+---+-----------------------------------------------------------------+
  !ctrlL
  !+---+-----------------------------------------------------------------+
  !..Cloud ice converting to snow since portion greater than min snow
  !.. size.  Given cloud ice content (kg/m**3), number concentration
  !.. (#/m**3) and gamma shape parameter, mu_i, break the distrib into
  !.. bins and figure out the mass/number of ice with sizes larger than
  !.. D0s.  Also, compute incomplete gamma function for the integration
  !.. of ice depositional growth from diameter=0 to D0s.  Amount of
  !.. ice depositional growth is this portion of distrib while larger
  !.. diameters contribute to snow growth (as in Harrington et al. 1995).
  !+---+-----------------------------------------------------------------+

  SUBROUTINE qi_aut_qs

    IMPLICIT NONE

    !..Local variables
    INTEGER                          :: i, j, n2
    REAL(KIND=r8), DIMENSION(nbi) :: N_i
    REAL(KIND=r8)                 :: N0_i, lami, Di_mean, t1, t2
    REAL(KIND=r8) :: xlimit_intg

    !+---+

    DO j = 1, ntb_i1
       DO i = 1, ntb_i
          lami = (am_i*cig(2)*oig1*Nt_i(j)/r_i(i))**obmi
          Di_mean = (bm_i + mu_i + 1.0_r8) / lami
          N0_i = Nt_i(j)*oig1 * lami**cie(1)
          t1 = 0.0e0_r8
          t2 = 0.0e0_r8
          IF (SNGL(Di_mean) .GT. 5.0_r8*D0s) THEN
             t1 = r_i(i)
             t2 = Nt_i(j)
             tpi_ide(i,j) = 0.0e0_r8
          ELSEIF (SNGL(Di_mean) .LT. D0i) THEN
             t1 = 0.0e0_r8
             t2 = 0.0e0_r8
             tpi_ide(i,j) = 1.0e0_r8
          ELSE
             xlimit_intg = lami*D0s
             tpi_ide(i,j) = GAMMP(mu_i+2.0_r8, xlimit_intg) * 1.0e0_r8
             DO n2 = 1, nbi
                N_i(n2) = N0_i*Di(n2)**mu_i * DEXP(-lami*Di(n2))*dti(n2)
                IF (Di(n2).GE.D0s) THEN
                   t1 = t1 + N_i(n2) * am_i*Di(n2)**bm_i
                   t2 = t2 + N_i(n2)
                ENDIF
             ENDDO
          ENDIF
          tps_iaus(i,j) = t1
          tni_iaus(i,j) = t2
       ENDDO
    ENDDO

  END SUBROUTINE qi_aut_qs


  !ctrlL
  !+---+-----------------------------------------------------------------+
  !..Variable collision efficiency for rain collecting cloud water using
  !.. method of Beard and Grover, 1974 if a/A less than 0.25; otherwise
  !.. uses polynomials to get close match of Pruppacher & Klett Fig 14-9.
  !+---+-----------------------------------------------------------------+

  SUBROUTINE table_Efrw()

    IMPLICIT NONE

    !..Local variables
    REAL(KIND=r8):: vtr, stokes, reynolds, Ef_rw
    REAL(KIND=r8):: p, yc0, F, G, H, z, K0, X
    INTEGER:: i, j

    DO j = 1, nbc
       DO i = 1, nbr
          Ef_rw = 0.0_r8
          p = Dc(j)/Dr(i)
          IF (Dr(i).LT.50.E-6_r8 .OR. Dc(j).LT.3.E-6_r8) THEN
             t_Efrw(i,j) = 0.0_r8
          ELSEIF (p.GT.0.25_r8) THEN
             X = Dc(j)*1.e6_r8
             IF (Dr(i) .LT. 75.e-6_r8) THEN
                Ef_rw = 0.026794_r8*X - 0.20604_r8
             ELSEIF (Dr(i) .LT. 125.e-6_r8) THEN
                Ef_rw = -0.00066842_r8*X*X + 0.061542_r8*X - 0.37089_r8
             ELSEIF (Dr(i) .LT. 175.e-6_r8) THEN
                Ef_rw = 4.091e-06_r8*X*X*X*X - 0.00030908_r8*X*X*X               &
                     + 0.0066237_r8*X*X - 0.0013687_r8*X - 0.073022_r8
             ELSEIF (Dr(i) .LT. 250.e-6_r8) THEN
                Ef_rw = 9.6719e-5_r8*X*X*X - 0.0068901_r8*X*X + 0.17305_r8*X        &
                     - 0.65988_r8
             ELSEIF (Dr(i) .LT. 350.e-6_r8) THEN
                Ef_rw = 9.0488e-5_r8*X*X*X - 0.006585_r8*X*X + 0.16606_r8*X         &
                     - 0.56125_r8
             ELSE
                Ef_rw = 0.00010721_r8*X*X*X - 0.0072962_r8*X*X + 0.1704_r8*X        &
                     - 0.46929_r8
             ENDIF
          ELSE
             vtr = -0.1021_r8 + 4.932E3_r8*Dr(i) - 0.9551E6_r8*Dr(i)*Dr(i) &
                  + 0.07934E9_r8*Dr(i)*Dr(i)*Dr(i) &
                  - 0.002362E12_r8*Dr(i)*Dr(i)*Dr(i)*Dr(i)
             stokes = Dc(j)*Dc(j)*vtr*rho_w/(9.0_r8*1.718E-5_r8*Dr(i))
             reynolds = 9.0_r8*stokes/(p*p*rho_w)

             F = DLOG(reynolds)
             G = -0.1007e0_r8 - 0.358e0_r8*F + 0.0261e0_r8*F*F
             K0 = DEXP(G)
             z = DLOG(stokes/(K0+1.e-15_r8))
             H = 0.1465e0_r8 + 1.302e0_r8*z - 0.607e0_r8*z*z + 0.293e0_r8*z*z*z
             yc0 = 2.0e0_r8/PI * ATAN(H)
             Ef_rw = (yc0+p)*(yc0+p) / ((1.0_r8+p)*(1.0_r8+p))

          ENDIF

          t_Efrw(i,j) = MAX(0.0_r8, MIN(SNGL(Ef_rw), 0.95_r8))

       ENDDO
    ENDDO

  END SUBROUTINE table_Efrw
  !ctrlL
  !+---+-----------------------------------------------------------------+
  !..Variable collision efficiency for snow collecting cloud water using
  !.. method of Wang and Ji, 2000 except equate melted snow diameter to
  !.. their "effective collision cross-section."
  !+---+-----------------------------------------------------------------+

  SUBROUTINE table_Efsw()

    IMPLICIT NONE

    !..Local variables
    REAL(KIND=r8):: Ds_m, vts, vtc, stokes, reynolds, Ef_sw
    REAL(KIND=r8):: p, yc0, F, G, H, z, K0
    INTEGER:: i, j

    DO j = 1, nbc
       vtc = 1.19e4_r8 * (1.0e4_r8*Dc(j)*Dc(j)*0.25e0_r8)
       DO i = 1, nbs
          vts = av_s*Ds(i)**bv_s * DEXP(-fv_s*Ds(i)) - vtc
          Ds_m = (am_s*Ds(i)**bm_s / am_r)**obmr
          p = Dc(j)/Ds_m
          IF (p.GT.0.25_r8 .OR. Ds(i).LT.D0s .OR. Dc(j).LT.6.E-6_r8 &
               .OR. vts.LT.1.E-3_r8) THEN
             t_Efsw(i,j) = 0.0_r8
          ELSE
             stokes = Dc(j)*Dc(j)*vts*rho_w/(9.0_r8*1.718E-5_r8*Ds_m)
             reynolds = 9.0_r8*stokes/(p*p*rho_w)

             F = DLOG(reynolds)
             G = -0.1007e0_r8 - 0.358e0_r8*F + 0.0261e0_r8*F*F
             K0 = DEXP(G)
             z = DLOG(stokes/(K0+1.e-15_r8))
             H = 0.1465e0_r8 + 1.302e0_r8*z - 0.607e0_r8*z*z + 0.293e0_r8*z*z*z
             yc0 = 2.0e0_r8/PI * ATAN(H)
             Ef_sw = (yc0+p)*(yc0+p) / ((1.0_r8+p)*(1.0_r8+p))

             t_Efsw(i,j) = MAX(0.0_r8, MIN(SNGL(Ef_sw), 0.95_r8))
          ENDIF

       ENDDO
    ENDDO

  END SUBROUTINE table_Efsw


  !ctrlL
  !+---+-----------------------------------------------------------------+
  !..Function to compute collision efficiency of collector species (rain,
  !.. snow, graupel) of aerosols.  Follows Wang et al, 2010, ACP, which
  !.. follows Slinn (1983).
  !+---+-----------------------------------------------------------------+

  REAL(KIND=r8)FUNCTION Eff_aero(D, Da, visc,rhoa,Temp,species)

    IMPLICIT NONE
    REAL(KIND=r8) :: D, Da, visc, rhoa, Temp
    CHARACTER(LEN=1):: species
    REAL(KIND=r8) :: aval, Cc, diff, Re, Sc, St, St2, vt, Eff
    REAL(KIND=r8), PARAMETER:: boltzman = 1.3806503E-23_r8
    REAL(KIND=r8), PARAMETER:: meanPath = 0.0256E-6_r8

    vt = 1.0_r8
    IF (species .EQ. 'r') THEN
       vt = -0.1021_r8 + 4.932E3_r8*D - 0.9551E6_r8*D*D                        &
            + 0.07934E9_r8*D*D*D - 0.002362E12_r8*D*D*D*D
    ELSEIF (species .EQ. 's') THEN
       vt = av_s*D**bv_s
    ELSEIF (species .EQ. 'g') THEN
       vt = av_g*D**bv_g
    ENDIF

    Cc    = 1.0_r8 + 2.0_r8*meanPath/Da *(1.257_r8+0.4_r8*EXP(-0.55_r8*Da/meanPath))
    diff  = boltzman*Temp*Cc/(3.0_r8*PI*visc*Da)

    Re    = 0.5_r8*rhoa*D*vt/visc
    Sc    = visc/(rhoa*diff)

    St    = Da*Da*vt*1000.0_r8/(9.0_r8*visc*D)
    aval  = 1.0_r8+LOG(1.0_r8+Re)
    St2   = (1.2_r8 + 1.0_r8/12.0_r8*aval)/(1.0_r8+aval)

    Eff = 4.0_r8/(Re*Sc) * (1.0_r8 + 0.4_r8*SQRT(Re)*Sc**0.3333_r8                  &
         + 0.16_r8*SQRT(Re)*SQRT(Sc))                  &
         + 4.0_r8*Da/D * (0.02_r8 + Da/D*(1.0_r8+2.0_r8*SQRT(Re)))

    IF (St.GT.St2) Eff = Eff  + ( (St-St2)/(St-St2+0.666667_r8))**1.5_r8
    Eff_aero = MAX(1.E-5_r8, MIN(Eff, 1.0_r8))

  END FUNCTION Eff_aero

  !ctrlL
  !+---+-----------------------------------------------------------------+
  !..Integrate rain size distribution from zero to D-star to compute the
  !.. number of drops smaller than D-star that evaporate in a single
  !.. timestep.  Drops larger than D-star dont evaporate entirely so do
  !.. not affect number concentration.
  !+---+-----------------------------------------------------------------+

  SUBROUTINE table_dropEvap()

    IMPLICIT NONE

    !..Local variables
    INTEGER:: i, j, k, n
    REAL(KIND=r8), DIMENSION(nbc):: N_c, massc
    REAL(KIND=r8):: summ, summ2, lamc, N0_c
    INTEGER:: nu_c
    !      REAL(KIND=r8):: Nt_r, N0, lam_exp, lam
    !      REAL(KIND=r8) :: xlimit_intg

    DO n = 1, nbc
       massc(n) = am_r*Dc(n)**bm_r
    ENDDO

    DO k = 1, nbc
       nu_c = MIN(15, NINT(1000.E6_r8/t_Nc(k)) + 2)
       DO j = 1, ntb_c
          lamc = (t_Nc(k)*am_r* ccg(2,nu_c)*ocg1(nu_c) / r_c(j))**obmr
          N0_c = t_Nc(k)*ocg1(nu_c) * lamc**cce(1,nu_c)
          DO i = 1, nbc
             !-GT           tnc_wev(i,j,k) = GAMMP(nu_c+1., SNGL(Dc(i)*lamc))*t_Nc(k)
             N_c(i) = N0_c* Dc(i)**nu_c*EXP(-lamc*Dc(i))*dtc(i)
             !     if(j.eq.18 .and. k.eq.50) print*, ' N_c = ', N_c(i)
             summ = 0._r8
             summ2 = 0._r8
             DO n = 1, i
                summ = summ + massc(n)*N_c(n)
                summ2 = summ2 + N_c(n)
             ENDDO
             !      if(j.eq.18 .and. k.eq.50) print*, '  DEBUG-TABLE: ', r_c(j), t_Nc(k), summ2, summ
             tpc_wev(i,j,k) = summ
             tnc_wev(i,j,k) = summ2
          ENDDO
       ENDDO
    ENDDO

    !
    !..To do the same thing for rain.
    !
    !     do k = 1, ntb_r
    !     do j = 1, ntb_r1
    !        lam_exp = (N0r_exp(j)*am_r*crg(1)/r_r(k))**ore1
    !        lam = lam_exp * (crg(3)*org2*org1)**obmr
    !        N0 = N0r_exp(j)/(crg(2)*lam_exp) * lam**cre(2)
    !        Nt_r = N0 * crg(2) / lam**cre(2)
    !        do i = 1, nbr
    !           xlimit_intg = lam*Dr(i)
    !           tnr_rev(i,j,k) = GAMMP(mu_r+1.0, xlimit_intg) * Nt_r
    !        enddo
    !     enddo
    !     enddo

    ! TO APPLY TABLE ABOVE
    !..Rain lookup table indexes.
    !         Dr_star = DSQRT(-2.D0*DT * t1_evap/(2.*PI) &
    !                 * 0.78*4.*diffu(k)*xsat*rvs/rho_w)
    !         idx_d = NINT(1.0 + FLOAT(nbr) * DLOG(Dr_star/D0r)             &
    !               / DLOG(Dr(nbr)/D0r))
    !         idx_d = MAX(1, MIN(idx_d, nbr))
    !
    !         nir = NINT(log10(rr(k)))
    !         do nn = nir-1, nir+1
    !            n = nn
    !            if ( (rr(k)/10.**nn).ge.1.0 .and. &
    !                 (rr(k)/10.**nn).lt.10.0) goto 154
    !         enddo
    !154      continue
    !         idx_r = INT(rr(k)/10.**n) + 10*(n-nir2) - (n-nir2)
    !         idx_r = MAX(1, MIN(idx_r, ntb_r))
    !
    !         lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
    !         lam_exp = lamr * (crg(3)*org2*org1)**bm_r
    !         N0_exp = org1*rr(k)/am_r * lam_exp**cre(1)
    !         nir = NINT(DLOG10(N0_exp))
    !         do nn = nir-1, nir+1
    !            n = nn
    !            if ( (N0_exp/10.**nn).ge.1.0 .and. &
    !                 (N0_exp/10.**nn).lt.10.0) goto 155
    !         enddo
    !155      continue
    !         idx_r1 = INT(N0_exp/10.**n) + 10*(n-nir3) - (n-nir3)
    !         idx_r1 = MAX(1, MIN(idx_r1, ntb_r1))
    !
    !         pnr_rev(k) = MIN(nr(k)*odts, SNGL(tnr_rev(idx_d,idx_r1,idx_r) &   ! RAIN2M
    !                    * odts))

  END SUBROUTINE table_dropEvap


  !
  !ctrlL
  !+---+-----------------------------------------------------------------+
  !..Fill the table of CCN activation data created from parcel model run
  !.. by Trude Eidhammer with inputs of aerosol number concentration,
  !.. vertical velocity, temperature, lognormal mean aerosol radius, and
  !.. hygroscopicity, kappa.  The data are read from external file and
  !.. contain activated fraction of CCN for given conditions.
  !+---+-----------------------------------------------------------------+

  SUBROUTINE table_ccnAct(path_in)

    !USE module_domain
    !USE module_dm
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN   ) :: path_in
    LOGICAL :: wrf_dm_on_monitor=.TRUE.

    !..Local variables
    INTEGER:: iunit_mp_th1, i
    LOGICAL:: opened
    CHARACTER*64 errmess

    iunit_mp_th1 = -1
    IF ( wrf_dm_on_monitor ) THEN
       DO i = 7,99
          INQUIRE ( i , OPENED = opened )
          IF ( .NOT. opened ) THEN
             iunit_mp_th1 = i
             GOTO 2010
          ENDIF
       ENDDO
2010   CONTINUE
    ENDIF
    !#if defined(DM_PARALLEL) && !defined(STUBMPI)
    !      CALL wrf_dm_bcast_bytes ( iunit_mp_th1 , IWORDSIZE )
    !#endif
    IF ( iunit_mp_th1 < 0 ) THEN
       CALL wrf_error_fatal ( 'module_mp_thompson: table_ccnAct: '//   &
            'Can not find unused fortran unit to read in lookup table.')
    ENDIF

    IF ( wrf_dm_on_monitor ) THEN
       WRITE(errmess, '(A,I2)') 'module_mp_thompson: opening CCN_ACTIVATE.BIN on unit ',iunit_mp_th1
       CALL wrf_debug('table_ccnAct...', errmess)
       OPEN(iunit_mp_th1,FILE=TRIM(path_in)//'/'//'CCN_ACTIVATE.BIN',                      &
            FORM='UNFORMATTED',STATUS='OLD',ERR=9009)
    ENDIF

    !#define DM_BCAST_MACRO(A) CALL wrf_dm_bcast_bytes(A, size(A)*R4SIZE)

    IF ( wrf_dm_on_monitor ) THEN
       READ(iunit_mp_th1,ERR=9010) tnccn_act
       CLOSE(iunit_mp_th1,STATUS='KEEP')
    END IF
    !#if defined(DM_PARALLEL) && !defined(STUBMPI)
    !      DM_BCAST_MACRO(tnccn_act)
    !#endif


    RETURN
9009 CONTINUE
    WRITE( errmess , '(A,I2)' ) 'module_mp_thompson: error opening CCN_ACTIVATE.BIN on unit ',iunit_mp_th1
    CALL wrf_error_fatal(errmess)
    RETURN
9010 CONTINUE
    WRITE( errmess , '(A,I2)' ) 'module_mp_thompson: error reading CCN_ACTIVATE.BIN on unit ',iunit_mp_th1
    CALL wrf_error_fatal(errmess)

  END SUBROUTINE table_ccnAct



  !^L
  !+---+-----------------------------------------------------------------+
  !..Retrieve fraction of CCN that gets activated given the model temp,
  !.. vertical velocity, and available CCN concentration.  The lookup
  !.. table (read from external file) has CCN concentration varying the
  !.. quickest, then updraft, then temperature, then mean aerosol radius,
  !.. and finally hygroscopicity, kappa.
  !.. TO_DO ITEM:  For radiation cooling producing fog, in which case the
  !.. updraft velocity could easily be negative, we could use the temp
  !.. and its tendency to diagnose a pretend postive updraft velocity.
  !+---+-----------------------------------------------------------------+
  REAL(KIND=r8)FUNCTION activ_ncloud(Tt, Ww, NCCN)

    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN):: Tt, Ww, NCCN
    REAL(KIND=r8) :: n_local, w_local
    INTEGER:: i, j, k, l, m, n
    REAL(KIND=r8) :: A, B, C, D, t, u, x1, x2, y1, y2, nx, wy, fraction


    !     ta_Na = (/10.0, 31.6, 100.0, 316.0, 1000.0, 3160.0, 10000.0/)  ntb_arc
    !     ta_Ww = (/0.01, 0.0316, 0.1, 0.316, 1.0, 3.16, 10.0, 31.6, 100.0/)  ntb_arw
    !     ta_Tk = (/243.15, 253.15, 263.15, 273.15, 283.15, 293.15, 303.15/)  ntb_art
    !     ta_Ra = (/0.01, 0.02, 0.04, 0.08, 0.16/)  ntb_arr
    !     ta_Ka = (/0.2, 0.4, 0.6, 0.8/)  ntb_ark

    n_local = NCCN * 1.E-6_r8
    w_local = Ww

    IF (n_local .GE. ta_Na(ntb_arc)) THEN
       n_local = ta_Na(ntb_arc) - 1.0_r8
    ELSEIF (n_local .LE. ta_Na(1)) THEN
       n_local = ta_Na(1) + 1.0_r8
    ENDIF
    DO n = 2, ntb_arc
       IF (n_local.GE.ta_Na(n-1) .AND. n_local.LT.ta_Na(n)) GOTO 8003
    ENDDO
8003 CONTINUE
    i = n
    x1 = LOG(ta_Na(i-1))
    x2 = LOG(ta_Na(i))

    IF (w_local .GE. ta_Ww(ntb_arw)) THEN
       w_local = ta_Ww(ntb_arw) - 1.0_r8
    ELSEIF (w_local .LE. ta_Ww(1)) THEN
       w_local = ta_Ww(1) + 0.001_r8
    ENDIF
    DO n = 2, ntb_arw
       IF (w_local.GE.ta_Ww(n-1) .AND. w_local.LT.ta_Ww(n)) GOTO 8005
    ENDDO
8005 CONTINUE
    j = n
    y1 = LOG(ta_Ww(j-1))
    y2 = LOG(ta_Ww(j))

    k = MAX(1, MIN( NINT( (Tt - ta_Tk(1))*0.1_r8) + 1, ntb_art))

    !..The next two values are indexes of mean aerosol radius and
    !.. hygroscopicity.  Currently these are constant but a future version
    !.. should implement other variables to allow more freedom such as
    !.. at least simple separation of tiny size sulfates from larger
    !.. sea salts.
    l = 3
    m = 2

    A = tnccn_act(i-1,j-1,k,l,m)
    B = tnccn_act(i,j-1,k,l,m)
    C = tnccn_act(i,j,k,l,m)
    D = tnccn_act(i-1,j,k,l,m)
    nx = LOG(n_local)
    wy = LOG(w_local)

    t = (nx-x1)/(x2-x1)
    u = (wy-y1)/(y2-y1)

    !     t = (n_local-ta(Na(i-1))/(ta_Na(i)-ta_Na(i-1))
    !     u = (w_local-ta_Ww(j-1))/(ta_Ww(j)-ta_Ww(j-1))

    fraction = (1.0_r8-t)*(1.0_r8-u)*A + t*(1.0_r8-u)*B + t*u*C + (1.0_r8-t)*u*D

    !     if (NCCN*fraction .gt. 0.75*Nt_c_max) then
    !        write(*,*) ' DEBUG-GT ', n_local, w_local, Tt, i, j, k
    !     endif

    activ_ncloud = NCCN*fraction

  END FUNCTION activ_ncloud
  !+---+-----------------------------------------------------------------+
  !+---+-----------------------------------------------------------------+
  SUBROUTINE GCF(GAMMCF,A,X,GLN)
    !     --- RETURNS THE INCOMPLETE GAMMA FUNCTION Q(A,X) EVALUATED BY ITS
    !     --- CONTINUED FRACTION REPRESENTATION AS GAMMCF.  ALSO RETURNS
    !     --- LN(GAMMA(A)) AS GLN.  THE CONTINUED FRACTION IS EVALUATED BY
    !     --- A MODIFIED LENTZ METHOD.
    !     --- USES GAMMLN
    IMPLICIT NONE
    INTEGER, PARAMETER:: ITMAX=100
    REAL(KIND=r8), PARAMETER:: gEPS=3.E-7_r8
    REAL(KIND=r8), PARAMETER:: FPMIN=1.E-30_r8
    REAL(KIND=r8), INTENT(IN):: A, X
    REAL(KIND=r8) :: GAMMCF,GLN
    INTEGER:: I
    REAL(KIND=r8) :: AN,B,C,D,DEL,H
    GLN=GAMMLN(A)
    B=X+1.0_r8-A
    C=1.0_r8/FPMIN
    D=1.0_r8/B
    H=D
    DO  I=1,ITMAX
       AN=-I*(I-A)
       B=B+2.0_r8
       D=AN*D+B
       IF(ABS(D).LT.FPMIN)D=FPMIN
       C=B+AN/C
       IF(ABS(C).LT.FPMIN)C=FPMIN
       D=1.0_r8/D
       DEL=D*C
       H=H*DEL
       IF(ABS(DEL-1.0_r8).LT.gEPS)GOTO 1
       !11        CONTINUE
    END DO
    PRINT *, 'A TOO LARGE, ITMAX TOO SMALL IN GCF'
1   GAMMCF=EXP(-X+A*LOG(X)-GLN)*H
  END SUBROUTINE GCF
  !  (C) Copr. 1986-92 Numerical Recipes Software 2.02
  !+---+-----------------------------------------------------------------+
  SUBROUTINE GSER(GAMSER,A,X,GLN)
    !     --- RETURNS THE INCOMPLETE GAMMA FUNCTION P(A,X) EVALUATED BY ITS
    !     --- ITS SERIES REPRESENTATION AS GAMSER.  ALSO RETURNS LN(GAMMA(A)) 
    !     --- AS GLN.
    !     --- USES GAMMLN
    IMPLICIT NONE
    INTEGER, PARAMETER:: ITMAX=100
    REAL(KIND=r8), PARAMETER:: gEPS=3.E-7_r8
    REAL(KIND=r8), INTENT(IN):: A, X
    REAL(KIND=r8) :: GAMSER,GLN
    INTEGER:: N
    REAL(KIND=r8) :: AP,DEL,SUM
    GLN=GAMMLN(A)
    IF(X.LE.0.0_r8)THEN
       IF(X.LT.0.0_r8) PRINT *, 'X < 0 IN GSER'
       GAMSER=0.0_r8
       RETURN
    ENDIF
    AP=A
    SUM=1.0_r8/A
    DEL=SUM
    DO  N=1,ITMAX
       AP=AP+1.0_r8
       DEL=DEL*X/AP
       SUM=SUM+DEL
       IF(ABS(DEL).LT.ABS(SUM)*gEPS)GOTO 1
       !11           CONTINUE
    END DO
    PRINT *,'A TOO LARGE, ITMAX TOO SMALL IN GSER'
1   GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)
  END SUBROUTINE GSER
  !  (C) Copr. 1986-92 Numerical Recipes Software 2.02
  !+---+-----------------------------------------------------------------+
  REAL(KIND=r8)FUNCTION GAMMLN(XX)
    !     --- RETURNS THE VALUE LN(GAMMA(XX)) FOR XX > 0.
    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN):: XX
    REAL(KIND=r8), PARAMETER:: STP = 2.5066282746310005e0_r8
    REAL(KIND=r8), DIMENSION(6), PARAMETER:: &
         COF = (/76.18009172947146e0_r8, -86.50532032941677e0_r8, &
         24.01409824083091e0_r8, -1.231739572450155e0_r8, &
         0.1208650973866179e-2_r8, -0.5395239384953e-5_r8/)
    REAL(KIND=r8):: SER,TMP,X,Y
    INTEGER:: J

    X=XX
    Y=X
    TMP=X+5.5e0_r8
    TMP=(X+0.5e0_r8)*LOG(TMP)-TMP
    SER=1.000000000190015e0_r8
    DO  J=1,6
       Y=Y+1.e0_r8
       SER=SER+COF(J)/Y
       !11              CONTINUE
    END DO
    GAMMLN=TMP+LOG(STP*SER/X)
  END FUNCTION GAMMLN
  !  (C) Copr. 1986-92 Numerical Recipes Software 2.02
  !+---+-----------------------------------------------------------------+
  REAL(KIND=r8)FUNCTION GAMMP(A,X)
    !     --- COMPUTES THE INCOMPLETE GAMMA FUNCTION P(A,X)
    !     --- SEE ABRAMOWITZ AND STEGUN 6.5.1
    !     --- USES GCF,GSER
    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN):: A,X
    REAL(KIND=r8) :: GAMMCF,GAMSER,GLN
    GAMMP = 0.0_r8
    IF((X.LT.0.0_r8) .OR. (A.LE.0.0_r8)) THEN
       PRINT *, 'BAD ARGUMENTS IN GAMMP'
       RETURN
    ELSEIF(X.LT.A+1.0_r8)THEN
       CALL GSER(GAMSER,A,X,GLN)
       GAMMP=GAMSER
    ELSE
       CALL GCF(GAMMCF,A,X,GLN)
       GAMMP=1.0_r8-GAMMCF
    ENDIF
  END FUNCTION GAMMP
  !  (C) Copr. 1986-92 Numerical Recipes Software 2.02
  !+---+-----------------------------------------------------------------+
  REAL(KIND=r8)FUNCTION WGAMMA(y)

    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN):: y

    WGAMMA = EXP(GAMMLN(y))

  END FUNCTION WGAMMA
  !+---+-----------------------------------------------------------------+
  ! THIS FUNCTION CALCULATES THE LIQUID SATURATION VAPOR MIXING RATIO AS
  ! A FUNCTION OF TEMPERATURE AND PRESSURE
  !
  REAL(KIND=r8)FUNCTION RSLF(P,T)

    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN):: P, T
    REAL(KIND=r8) :: ESL,X
    REAL(KIND=r8), PARAMETER:: C0= .611583699E03_r8
    REAL(KIND=r8), PARAMETER:: C1= .444606896E02_r8
    REAL(KIND=r8), PARAMETER:: C2= .143177157E01_r8
    REAL(KIND=r8), PARAMETER:: C3= .264224321E-1_r8
    REAL(KIND=r8), PARAMETER:: C4= .299291081E-3_r8
    REAL(KIND=r8), PARAMETER:: C5= .203154182E-5_r8
    REAL(KIND=r8), PARAMETER:: C6= .702620698E-8_r8
    REAL(KIND=r8), PARAMETER:: C7= .379534310E-11_r8
    REAL(KIND=r8), PARAMETER:: C8=-.321582393E-13_r8

    X=MAX(-80.0_r8,T-273.16_r8)

    !      ESL=612.2*EXP(17.67*X/(T-29.65))
    ESL=C0+X*(C1+X*(C2+X*(C3+X*(C4+X*(C5+X*(C6+X*(C7+X*C8)))))))
    RSLF=0.622_r8*ESL/(P-ESL)

    !    ALTERNATIVE
    !  ; Source: Murphy and Koop, Review of the vapour pressure of ice and
    !             supercooled water for atmospheric applications, Q. J. R.
    !             Meteorol. Soc (2005), 131, pp. 1539-1565.
    !    ESL = EXP(54.842763 - 6763.22 / T - 4.210 * log(T) + 0.000367 * T
    !        + TANH(0.0415 * (T - 218.8)) * (53.878 - 1331.22
    !        / T - 9.44523 * log(T) + 0.014025 * T))

  END FUNCTION RSLF
  !+---+-----------------------------------------------------------------+
  ! THIS FUNCTION CALCULATES THE ICE SATURATION VAPOR MIXING RATIO AS A
  ! FUNCTION OF TEMPERATURE AND PRESSURE
  !
  REAL(KIND=r8)FUNCTION RSIF(P,T)

    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN):: P, T
    REAL(KIND=r8) :: ESI,X
    REAL(KIND=r8), PARAMETER:: C0= .609868993E03_r8
    REAL(KIND=r8), PARAMETER:: C1= .499320233E02_r8
    REAL(KIND=r8), PARAMETER:: C2= .184672631E01_r8
    REAL(KIND=r8), PARAMETER:: C3= .402737184E-1_r8
    REAL(KIND=r8), PARAMETER:: C4= .565392987E-3_r8
    REAL(KIND=r8), PARAMETER:: C5= .521693933E-5_r8
    REAL(KIND=r8), PARAMETER:: C6= .307839583E-7_r8
    REAL(KIND=r8), PARAMETER:: C7= .105785160E-9_r8
    REAL(KIND=r8), PARAMETER:: C8= .161444444E-12_r8

    X=MAX(-80.0_r8,T-273.16_r8)
    ESI=C0+X*(C1+X*(C2+X*(C3+X*(C4+X*(C5+X*(C6+X*(C7+X*C8)))))))
    RSIF=0.622_r8*ESI/(P-ESI)

    !    ALTERNATIVE
    !  ; Source: Murphy and Koop, Review of the vapour pressure of ice and
    !             supercooled water for atmospheric applications, Q. J. R.
    !             Meteorol. Soc (2005), 131, pp. 1539-1565.
    !     ESI = EXP(9.550426 - 5723.265/T + 3.53068*log(T) - 0.00728332*T)

  END FUNCTION RSIF

  !+---+-----------------------------------------------------------------+
  REAL(KIND=r8)FUNCTION iceDeMott(tempc, qv, qvs, qvsi, rho, nifa)
    IMPLICIT NONE

    REAL(KIND=r8), INTENT(IN):: tempc, qv, qvs, qvsi, rho, nifa

    !..Local vars
    REAL(KIND=r8) :: satw, sati, siw, p_x, si0x, dtt, dsi, dsw, dab, fc, hx
    REAL(KIND=r8) :: ntilde, n_in, nmax, nhat, mux, xni, nifa_cc
    REAL(KIND=r8), PARAMETER:: p_c1    = 1000.0_r8
    REAL(KIND=r8), PARAMETER:: p_rho_c = 0.76_r8
    REAL(KIND=r8), PARAMETER:: p_alpha = 1.0_r8
    REAL(KIND=r8), PARAMETER:: p_gam   = 2.0_r8
    REAL(KIND=r8), PARAMETER:: delT    = 5.0_r8
    REAL(KIND=r8), PARAMETER:: T0x     = -40.0_r8
    REAL(KIND=r8), PARAMETER:: Sw0x    = 0.97_r8
    REAL(KIND=r8), PARAMETER:: delSi   = 0.1_r8
    REAL(KIND=r8), PARAMETER:: hdm     = 0.15_r8
    REAL(KIND=r8), PARAMETER:: p_psi   = 0.058707_r8*p_gam/p_rho_c
    REAL(KIND=r8), PARAMETER:: aap     = 1.0_r8
    REAL(KIND=r8), PARAMETER:: bbp     = 0.0_r8
    REAL(KIND=r8), PARAMETER:: y1p     = -35.0_r8
    REAL(KIND=r8), PARAMETER:: y2p     = -25.0_r8
    REAL(KIND=r8), PARAMETER:: rho_not0 = 101325.0_r8/(287.05_r8*273.15_r8)

    !+---+

    xni = 0.0_r8
    satw = qv/qvs
    sati = qv/qvsi
    siw = qvs/qvsi
    p_x = -1.0261_r8+(3.1656e-3_r8*tempc)+(5.3938e-4_r8*(tempc*tempc))         &
         +  (8.2584e-6_r8*(tempc*tempc*tempc))
    si0x = 1.0_r8+(10.0_r8**p_x)
    IF (sati.GE.si0x .AND. satw.LT.0.985_r8) THEN
       dtt = delta_p (tempc, T0x, T0x+delT, 1.0_r8, hdm)
       dsi = delta_p (sati, Si0x, Si0x+delSi, 0.0_r8, 1.0_r8)
       dsw = delta_p (satw, Sw0x, 1.0_r8, 0.0_r8, 1.0_r8)
       fc = dtt*dsi*0.5_r8
       hx = MIN(fc+((1.0_r8-fc)*dsw), 1.0_r8)
       ntilde = p_c1*p_gam*((EXP(12.96_r8*(sati-1.1_r8)))**0.3_r8) / p_rho_c
       IF (tempc .LE. y1p) THEN
          n_in = ntilde
       ELSEIF (tempc .GE. y2p) THEN
          n_in = p_psi*p_c1*EXP(12.96_r8*(sati-1.0_r8)-0.639_r8)
       ELSE
          IF (tempc .LE. -30.0_r8) THEN
             nmax = p_c1*p_gam*(EXP(12.96_r8*(siw-1.1_r8)))**0.3_r8/p_rho_c
          ELSE
             nmax = p_psi*p_c1*EXP(12.96_r8*(siw-1.0_r8)-0.639_r8)
          ENDIF
          ntilde = MIN(ntilde, nmax)
          nhat = MIN(p_psi*p_c1*EXP(12.96_r8*(sati-1.)-0.639_r8), nmax)
          dab = delta_p (tempc, y1p, y2p, aap, bbp)
          n_in = MIN(nhat*(ntilde/nhat)**dab, nmax)
       ENDIF
       mux = hx*p_alpha*n_in*rho
       xni = mux*((6700.0_r8*nifa)-200.0_r8)/((6700.0_r8*5.E5_r8)-200.0_r8)
    ELSEIF (satw.GE.0.985_r8 .AND. tempc.GT.HGFR-273.15_r8) THEN
       nifa_cc = nifa*RHO_NOT0*1.E-6_r8/rho
       xni  = 3.0_r8*nifa_cc**(1.25_r8)*EXP((0.46_r8*(-tempc))-11.6_r8)             !  [DeMott, 2015]
       !        xni = (5.94e-5*(-tempc)**3.33)                                 &
       !                   * (nifa_cc**((-0.0264*(tempc))+0.0033))
       xni = xni*rho/RHO_NOT0 * 1000.0_r8
    ENDIF

    iceDeMott = MAX(0.0_r8, xni)

  END FUNCTION iceDeMott

  !+---+-----------------------------------------------------------------+
  !..Newer research since Koop et al (2001) suggests that the freezing
  !.. rate should be lower than original paper, so J_rate is reduced
  !.. by two orders of magnitude.

  REAL(KIND=r8)FUNCTION iceKoop(temp, qv, qvs, naero, dt)
    IMPLICIT NONE

    REAL(KIND=r8), INTENT(IN):: temp, qv, qvs, naero, DT
    REAL(KIND=r8) :: mu_diff, a_w_i, delta_aw, log_J_rate, J_rate, prob_h, satw
    REAL(KIND=r8) :: xni

    xni = 0.0_r8
    satw = qv/qvs
    mu_diff    = 210368.0_r8 + (131.438_r8*temp) - (3.32373E6_r8/temp)         &
         &           - (41729.1_r8*log(temp))
    a_w_i      = EXP(mu_diff/(R_uni*temp))
    delta_aw   = satw - a_w_i
    log_J_rate = -906.7_r8 + (8502.0_r8*delta_aw)                           &
         &           - (26924.0_r8*delta_aw*delta_aw)                          &
         &           + (29180.0_r8*delta_aw*delta_aw*delta_aw)
    log_J_rate = MIN(20.0_r8, log_J_rate)
    J_rate     = 0.01_r8*(10.0_r8**log_J_rate)                                ! cm-3 s-1
    prob_h     = MIN(1.0_r8-EXP(-J_rate*ar_volume*DT), 1.0_r8)
    IF (prob_h .GT. 0.0_r8) THEN
       xni     = MIN(prob_h*naero, 1000.E3_r8)
    ENDIF

    iceKoop = MAX(0.0_r8, xni)

  END FUNCTION iceKoop

  !+---+-----------------------------------------------------------------+
  !.. Helper routine for Phillips et al (2008) ice nucleation.  Trude

  REAL(KIND=r8)FUNCTION delta_p (yy, y1, y2, aa, bb)
    IMPLICIT NONE

    REAL(KIND=r8), INTENT(IN):: yy, y1, y2, aa, bb
    REAL(KIND=r8) :: dab, A, B, a0, a1, a2, a3

    A   = 6.0_r8*(aa-bb)/((y2-y1)*(y2-y1)*(y2-y1))
    B   = aa+(A*y1*y1*y1/6.0_r8)-(A*y1*y1*y2*0.5_r8)
    a0  = B
    a1  = A*y1*y2
    a2  = -A*(y1+y2)*0.5_r8
    a3  = A/3.0_r8

    IF (yy.LE.y1) THEN 
       dab = aa
    ELSE IF (yy.GE.y2) THEN
       dab = bb
    ELSE
       dab = a0+(a1*yy)+(a2*yy*yy)+(a3*yy*yy*yy)
    ENDIF

    IF (dab.LT.aa) THEN 
       dab = aa
    ENDIF
    IF (dab.GT.bb) THEN 
       dab = bb
    ENDIF
    delta_p = dab

  END FUNCTION delta_p

  !+---+-----------------------------------------------------------------+
  !ctrlL

  !+---+-----------------------------------------------------------------+
  !..Compute _radiation_ effective radii of cloud water, ice, and snow.
  !.. These are entirely consistent with microphysics assumptions, not
  !.. constant or otherwise ad hoc as is internal to most radiation
  !.. schemes.  Since only the smallest snowflakes should impact
  !.. radiation, compute from first portion of complicated Field number
  !.. distribution, not the second part, which is the larger sizes.
  !+---+-----------------------------------------------------------------+

  SUBROUTINE calc_effectRad (t1d, p1d, qv1d, qc1d, nc1d, qi1d, ni1d, qs1d,   &
       &                re_qc1d, re_qi1d, re_qs1d, kts, kte)

    IMPLICIT NONE

    !..Sub arguments
    INTEGER, INTENT(IN):: kts, kte
    REAL(KIND=r8), DIMENSION(kts:kte), INTENT(IN)::                            &
         &                    t1d, p1d, qv1d, qc1d, nc1d, qi1d, ni1d, qs1d
    REAL(KIND=r8), DIMENSION(kts:kte), INTENT(INOUT):: re_qc1d, re_qi1d, re_qs1d
    !..Local variables
    INTEGER:: k
    REAL(KIND=r8), DIMENSION(kts:kte):: rho, rc, nc, ri, ni, rs
    REAL(KIND=r8) :: smo2, smob, smoc
    REAL(KIND=r8) :: tc0, loga_, a_, b_
    REAL(KIND=r8):: lamc, lami
    LOGICAL:: has_qc, has_qi, has_qs
    INTEGER:: inu_c
    REAL(KIND=r8), DIMENSION(15), PARAMETER:: g_ratio = (/24.0_r8,60.0_r8,120.0_r8,210.0_r8,336.0_r8,   &
         &                504.0_r8,720.0_r8,990.0_r8,1320.0_r8,1716.0_r8,2184.0_r8,2730.0_r8,3360.0_r8,4080.0_r8,4896.0_r8/)

    has_qc = .FALSE.
    has_qi = .FALSE.
    has_qs = .FALSE.

    DO k = kts, kte
       rho(k) = 0.622_r8*p1d(k)/(R*t1d(k)*(qv1d(k)+0.622_r8))
       rc(k) = MAX(R1, qc1d(k)*rho(k))
       nc(k) = MAX(R2, nc1d(k)*rho(k))
       IF (.NOT. is_aerosol_aware) nc(k) = Nt_c
       IF (rc(k).GT.R1 .AND. nc(k).GT.R2) has_qc = .TRUE.
       ri(k) = MAX(R1, qi1d(k)*rho(k))
       ni(k) = MAX(R2, ni1d(k)*rho(k))
       IF (ri(k).GT.R1 .AND. ni(k).GT.R2) has_qi = .TRUE.
       rs(k) = MAX(R1, qs1d(k)*rho(k))
       IF (rs(k).GT.R1) has_qs = .TRUE.
    ENDDO

    IF (has_qc) THEN
       DO k = kts, kte
          IF (rc(k).LE.R1 .OR. nc(k).LE.R2) CYCLE
          IF (nc(k).LT.100) THEN
             inu_c = 15
          ELSEIF (nc(k).GT.1.E10_r8) THEN
             inu_c = 2
          ELSE
             inu_c = MIN(15, NINT(1000.E6_r8/nc(k)) + 2)
          ENDIF
          lamc = (nc(k)*am_r*g_ratio(inu_c)/rc(k))**obmr
          re_qc1d(k) = MAX(2.51E-6_r8, MIN(SNGL(0.5e0_r8 * DBLE(3.0_r8+inu_c)/lamc), 50.E-6_r8))
       ENDDO
    ENDIF

    IF (has_qi) THEN
       DO k = kts, kte
          IF (ri(k).LE.R1 .OR. ni(k).LE.R2) CYCLE
          lami = (am_i*cig(2)*oig1*ni(k)/ri(k))**obmi
          re_qi1d(k) = MAX(10.01E-6_r8, MIN(SNGL(0.5e0_r8 * DBLE(3.0_r8+mu_i)/lami), 125.E-6_r8))
       ENDDO
    ENDIF

    IF (has_qs) THEN
       DO k = kts, kte
          IF (rs(k).LE.R1) CYCLE
          tc0 = MIN(-0.1_r8, t1d(k)-273.15_r8)
          smob = rs(k)*oams

          !..All other moments based on reference, 2nd moment.  If bm_s.ne.2,
          !.. then we must compute actual 2nd moment and use as reference.
          IF (bm_s.GT.(2.0_r8-1.e-3_r8) .AND. bm_s.LT.(2.0_r8+1.e-3_r8)) THEN
             smo2 = smob
          ELSE
             loga_ = sa(1) + sa(2)*tc0 + sa(3)*bm_s &
                  &         + sa(4)*tc0*bm_s + sa(5)*tc0*tc0 &
                  &         + sa(6)*bm_s*bm_s + sa(7)*tc0*tc0*bm_s &
                  &         + sa(8)*tc0*bm_s*bm_s + sa(9)*tc0*tc0*tc0 &
                  &         + sa(10)*bm_s*bm_s*bm_s
             a_ = 10.0_r8**loga_
             b_ = sb(1) + sb(2)*tc0 + sb(3)*bm_s &
                  &         + sb(4)*tc0*bm_s + sb(5)*tc0*tc0 &
                  &         + sb(6)*bm_s*bm_s + sb(7)*tc0*tc0*bm_s &
                  &         + sb(8)*tc0*bm_s*bm_s + sb(9)*tc0*tc0*tc0 &
                  &         + sb(10)*bm_s*bm_s*bm_s
             smo2 = (smob/a_)**(1.0_r8/b_)
          ENDIF
          !..Calculate bm_s+1 (th) moment.  Useful for diameter calcs.
          loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(1) &
               &         + sa(4)*tc0*cse(1) + sa(5)*tc0*tc0 &
               &         + sa(6)*cse(1)*cse(1) + sa(7)*tc0*tc0*cse(1) &
               &         + sa(8)*tc0*cse(1)*cse(1) + sa(9)*tc0*tc0*tc0 &
               &         + sa(10)*cse(1)*cse(1)*cse(1)
          a_ = 10.0_r8**loga_
          b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(1) + sb(4)*tc0*cse(1) &
               &        + sb(5)*tc0*tc0 + sb(6)*cse(1)*cse(1) &
               &        + sb(7)*tc0*tc0*cse(1) + sb(8)*tc0*cse(1)*cse(1) &
               &        + sb(9)*tc0*tc0*tc0 + sb(10)*cse(1)*cse(1)*cse(1)
          smoc = a_ * smo2**b_
          re_qs1d(k) = MAX(10.E-6_r8, MIN(0.5_r8*(smoc/smob), 999.E-6_r8))
       ENDDO
    ENDIF

  END SUBROUTINE calc_effectRad

  !+---+-----------------------------------------------------------------+
  !..Compute radar reflectivity assuming 10 cm wavelength radar and using
  !.. Rayleigh approximation.  Only complication is melted snow/graupel
  !.. which we treat as water-coated ice spheres and use Uli Blahak's
  !.. library of routines.  The meltwater fraction is simply the amount
  !.. of frozen species remaining from what initially existed at the
  !.. melting level interface.
  !+---+-----------------------------------------------------------------+

  SUBROUTINE calc_refl10cm (qv1d, qc1d, qr1d, nr1d, qs1d, qg1d,     &
       t1d, p1d, dBZ, kts, kte)

    IMPLICIT NONE

    !..Sub arguments
    INTEGER, INTENT(IN):: kts, kte
    REAL(KIND=r8), DIMENSION(kts:kte), INTENT(IN)::                            &
         qv1d, qc1d, qr1d, nr1d, qs1d, qg1d, t1d, p1d
    REAL(KIND=r8), DIMENSION(kts:kte), INTENT(INOUT):: dBZ
    !     REAL(KIND=r8), DIMENSION(kts:kte), INTENT(INOUT):: vt_dBZ

    !..Local variables
    REAL(KIND=r8), DIMENSION(kts:kte):: temp, pres, qv, rho, rhof
    REAL(KIND=r8), DIMENSION(kts:kte):: rc, rr, nr, rs, rg

    REAL(KIND=r8), DIMENSION(kts:kte):: ilamr, ilamg, N0_r, N0_g
    REAL(KIND=r8), DIMENSION(kts:kte):: mvd_r
    REAL(KIND=r8), DIMENSION(kts:kte):: smob, smo2, smoc, smoz
    REAL(KIND=r8) :: oM3, M0, Mrat, slam1, slam2!, xDs
 !   REAL(KIND=r8) :: ils1, ils2!, t1_vts, t2_vts, t3_vts, t4_vts
 !   REAL(KIND=r8) :: vtr_dbz_wt, vts_dbz_wt, vtg_dbz_wt

    REAL(KIND=r8), DIMENSION(kts:kte):: ze_rain, ze_snow, ze_graupel

    REAL(KIND=r8):: N0_exp, N0_min, lam_exp, lamr, lamg
    REAL(KIND=r8) :: a_, b_, loga_, tc0
    REAL(KIND=r8):: fmelt_s, fmelt_g

    INTEGER:: k, k_0, n
    LOGICAL:: melti
    LOGICAL, DIMENSION(kts:kte):: L_qr, L_qs, L_qg

    REAL(KIND=r8):: cback, x, eta, f_d
    REAL(KIND=r8) :: xslw1, ygra1, zans1

    !+---+

    DO k = kts, kte
       dBZ(k) = -35.0_r8
    ENDDO

    !+---+-----------------------------------------------------------------+
    !..Put column of data into local arrays.
    !+---+-----------------------------------------------------------------+
    DO k = kts, kte
       temp(k) = t1d(k)
       qv(k) = MAX(1.E-10_r8, qv1d(k))
       pres(k) = p1d(k)
       rho(k) = 0.622_r8*pres(k)/(R*temp(k)*(qv(k)+0.622_r8))
       rhof(k) = SQRT(RHO_NOT/rho(k))
       rc(k) = MAX(R1, qc1d(k)*rho(k))
       IF (qr1d(k) .GT. R1) THEN
          rr(k) = qr1d(k)*rho(k)
          nr(k) = MAX(R2, nr1d(k)*rho(k))
          lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
          ilamr(k) = 1.0_r8/lamr
          N0_r(k) = nr(k)*org2*lamr**cre(2)
          mvd_r(k) = (3.0_r8 + mu_r + 0.672_r8) * ilamr(k)
          L_qr(k) = .TRUE.
       ELSE
          rr(k) = R1
          nr(k) = R1
          mvd_r(k) = 50.E-6_r8
          L_qr(k) = .FALSE.
       ENDIF
       IF (qs1d(k) .GT. R2) THEN
          rs(k) = qs1d(k)*rho(k)
          L_qs(k) = .TRUE.
       ELSE
          rs(k) = R1
          L_qs(k) = .FALSE.
       ENDIF
       IF (qg1d(k) .GT. R2) THEN
          rg(k) = qg1d(k)*rho(k)
          L_qg(k) = .TRUE.
       ELSE
          rg(k) = R1
          L_qg(k) = .FALSE.
       ENDIF
    ENDDO

    !+---+-----------------------------------------------------------------+
    !..Calculate y-intercept, slope, and useful moments for snow.
    !+---+-----------------------------------------------------------------+
    DO k = kts, kte
       tc0 = MIN(-0.1_r8, temp(k)-273.15_r8)
       smob(k) = rs(k)*oams

       !..All other moments based on reference, 2nd moment.  If bm_s.ne.2,
       !.. then we must compute actual 2nd moment and use as reference.
       IF (bm_s.GT.(2.0_r8-1.e-3_r8) .AND. bm_s.LT.(2.0_r8+1.e-3_r8)) THEN
          smo2(k) = smob(k)
       ELSE
          loga_ = sa(1) + sa(2)*tc0 + sa(3)*bm_s &
               &         + sa(4)*tc0*bm_s + sa(5)*tc0*tc0 &
               &         + sa(6)*bm_s*bm_s + sa(7)*tc0*tc0*bm_s &
               &         + sa(8)*tc0*bm_s*bm_s + sa(9)*tc0*tc0*tc0 &
               &         + sa(10)*bm_s*bm_s*bm_s
          a_ = 10.0_r8**loga_
          b_ = sb(1) + sb(2)*tc0 + sb(3)*bm_s &
               &         + sb(4)*tc0*bm_s + sb(5)*tc0*tc0 &
               &         + sb(6)*bm_s*bm_s + sb(7)*tc0*tc0*bm_s &
               &         + sb(8)*tc0*bm_s*bm_s + sb(9)*tc0*tc0*tc0 &
               &         + sb(10)*bm_s*bm_s*bm_s
          smo2(k) = (smob(k)/a_)**(1.0_r8/b_)
       ENDIF

       !..Calculate bm_s+1 (th) moment.  Useful for diameter calcs.
       loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(1) &
            &         + sa(4)*tc0*cse(1) + sa(5)*tc0*tc0 &
            &         + sa(6)*cse(1)*cse(1) + sa(7)*tc0*tc0*cse(1) &
            &         + sa(8)*tc0*cse(1)*cse(1) + sa(9)*tc0*tc0*tc0 &
            &         + sa(10)*cse(1)*cse(1)*cse(1)
       a_ = 10.0_r8**loga_
       b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(1) + sb(4)*tc0*cse(1) &
            &        + sb(5)*tc0*tc0 + sb(6)*cse(1)*cse(1) &
            &        + sb(7)*tc0*tc0*cse(1) + sb(8)*tc0*cse(1)*cse(1) &
            &        + sb(9)*tc0*tc0*tc0 + sb(10)*cse(1)*cse(1)*cse(1)
       smoc(k) = a_ * smo2(k)**b_

       !..Calculate bm_s*2 (th) moment.  Useful for reflectivity.
       loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(3) &
            &         + sa(4)*tc0*cse(3) + sa(5)*tc0*tc0 &
            &         + sa(6)*cse(3)*cse(3) + sa(7)*tc0*tc0*cse(3) &
            &         + sa(8)*tc0*cse(3)*cse(3) + sa(9)*tc0*tc0*tc0 &
            &         + sa(10)*cse(3)*cse(3)*cse(3)
       a_ = 10.0_r8**loga_
       b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(3) + sb(4)*tc0*cse(3) &
            &        + sb(5)*tc0*tc0 + sb(6)*cse(3)*cse(3) &
            &        + sb(7)*tc0*tc0*cse(3) + sb(8)*tc0*cse(3)*cse(3) &
            &        + sb(9)*tc0*tc0*tc0 + sb(10)*cse(3)*cse(3)*cse(3)
       smoz(k) = a_ * smo2(k)**b_
    ENDDO

    !+---+-----------------------------------------------------------------+
    !..Calculate y-intercept, slope values for graupel.
    !+---+-----------------------------------------------------------------+

    N0_min = gonv_max
    DO k = kte, kts, -1
       IF (temp(k).LT.270.65_r8 .AND. L_qr(k) .AND. mvd_r(k).GT.100.E-6_r8) THEN
          xslw1 = 4.01_r8 + log10(mvd_r(k))
       ELSE
          xslw1 = 0.01_r8
       ENDIF
       ygra1 = 4.31_r8 + log10(MAX(5.E-5_r8, rg(k)))
       zans1 = 3.1_r8 + (100.0_r8/(300.0_r8*xslw1*ygra1/(10.0_r8/xslw1+1.0_r8+0.25_r8*ygra1)+30.0_r8+10.0_r8*ygra1))
       N0_exp = 10.0_r8**(zans1)
       N0_exp = MAX(DBLE(gonv_min), MIN(N0_exp, DBLE(gonv_max)))
       N0_min = MIN(N0_exp, N0_min)
       N0_exp = N0_min
       lam_exp = (N0_exp*am_g*cgg(1)/rg(k))**oge1
       lamg = lam_exp * (cgg(3)*ogg2*ogg1)**obmg
       ilamg(k) = 1.0_r8/lamg
       N0_g(k) = N0_exp/(cgg(2)*lam_exp) * lamg**cge(2)
    ENDDO

    !+---+-----------------------------------------------------------------+
    !..Locate K-level of start of melting (k_0 is level above).
    !+---+-----------------------------------------------------------------+
    melti = .FALSE.
    k_0 = kts
    DO k = kte-1, kts, -1
       IF ( (temp(k).GT.273.15_r8) .AND. L_qr(k)                         &
            &                            .AND. (L_qs(k+1).OR.L_qg(k+1)) ) THEN
          k_0 = MAX(k+1, k_0)
          melti=.TRUE.
          GOTO 195
       ENDIF
    ENDDO
195 CONTINUE

    !+---+-----------------------------------------------------------------+
    !..Assume Rayleigh approximation at 10 cm wavelength. Rain (all temps)
    !.. and non-water-coated snow and graupel when below freezing are
    !.. simple. Integrations of m(D)*m(D)*N(D)*dD.
    !+---+-----------------------------------------------------------------+

    DO k = kts, kte
       ze_rain(k) = 1.e-22_r8
       ze_snow(k) = 1.e-22_r8
       ze_graupel(k) = 1.e-22_r8
       IF (L_qr(k)) ze_rain(k) = N0_r(k)*crg(4)*ilamr(k)**cre(4)
       IF (L_qs(k)) ze_snow(k) = (0.176_r8/0.93_r8) * (6.0_r8/PI)*(6.0_r8/PI)     &
            &                           * (am_s/900.0_r8)*(am_s/900.0_r8)*smoz(k)
       IF (L_qg(k)) ze_graupel(k) = (0.176_r8/0.93_r8) * (6.0_r8/PI)*(6.0_r8/PI)  &
            &                              * (am_g/900.0_r8)*(am_g/900.0_r8)         &
            &                              * N0_g(k)*cgg(4)*ilamg(k)**cge(4)
    ENDDO

    !+---+-----------------------------------------------------------------+
    !..Special case of melting ice (snow/graupel) particles.  Assume the
    !.. ice is surrounded by the liquid water.  Fraction of meltwater is
    !.. extremely simple based on amount found above the melting level.
    !.. Uses code from Uli Blahak (rayleigh_soak_wetgraupel and supporting
    !.. routines).
    !+---+-----------------------------------------------------------------+

    IF (.NOT. iiwarm .AND. melti .AND. k_0.GE.2) THEN
       DO k = k_0-1, kts, -1

          !..Reflectivity contributed by melting snow
          IF (L_qs(k) .AND. L_qs(k_0) ) THEN
             fmelt_s = MAX(0.05e0_r8, MIN(1.0e0_r8-rs(k)/rs(k_0), 0.99e0_r8))
             eta = 0.e0_r8
             oM3 = 1.0_r8/smoc(k)
             M0 = (smob(k)*oM3)
             Mrat = smob(k)*M0*M0*M0
             slam1 = M0 * Lam0
             slam2 = M0 * Lam1
             DO n = 1, nrbins
                x = am_s * xxDs(n)**bm_s
                CALL rayleigh_soak_wetgraupel (x, DBLE(ocms), DBLE(obms), &
                     &              fmelt_s, melt_outside_s, m_w_0, m_i_0, &
                     &              CBACK, mixingrulestring_s, matrixstring_s,          &
                     &              inclusionstring_s, hoststring_s,                    &
                     &              hostmatrixstring_s, hostinclusionstring_s)
                f_d = Mrat*(Kap0*DEXP(-slam1*xxDs(n))                     &
                     &              + Kap1*(M0*xxDs(n))**mu_s * DEXP(-slam2*xxDs(n)))
                eta = eta + f_d * CBACK * simpson(n) * xdts(n)
             ENDDO
             ze_snow(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
          ENDIF

          !..Reflectivity contributed by melting graupel

          IF (L_qg(k) .AND. L_qg(k_0) ) THEN
             fmelt_g = MAX(0.05e0_r8, MIN(1.0e0_r8-rg(k)/rg(k_0), 0.99e0_r8))
             eta = 0.e0_r8
             lamg = 1.0_r8/ilamg(k)
             DO n = 1, nrbins
                x = am_g * xxDg(n)**bm_g
                CALL rayleigh_soak_wetgraupel (x, DBLE(ocmg), DBLE(obmg), &
                     &              fmelt_g, melt_outside_g, m_w_0, m_i_0, &
                     &              CBACK, mixingrulestring_g, matrixstring_g,          &
                     &              inclusionstring_g, hoststring_g,                    &
                     &              hostmatrixstring_g, hostinclusionstring_g)
                f_d = N0_g(k)*xxDg(n)**mu_g * DEXP(-lamg*xxDg(n))
                eta = eta + f_d * CBACK * simpson(n) * xdtg(n)
             ENDDO
             ze_graupel(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
          ENDIF

       ENDDO
    ENDIF

    DO k = kte, kts, -1
       dBZ(k) = 10.0_r8*LOG10((ze_rain(k)+ze_snow(k)+ze_graupel(k))*1.e18_r8)
    ENDDO


    !..Reflectivity-weighted terminal velocity (snow, rain, graupel, mix).
    !     do k = kte, kts, -1
    !        vt_dBZ(k) = 1.E-3
    !        if (rs(k).gt.R2) then
    !         Mrat = smob(k) / smoc(k)
    !         ils1 = 1./(Mrat*Lam0 + fv_s)
    !         ils2 = 1./(Mrat*Lam1 + fv_s)
    !         t1_vts = Kap0*csg(5)*ils1**cse(5)
    !         t2_vts = Kap1*Mrat**mu_s*csg(11)*ils2**cse(11)
    !         ils1 = 1./(Mrat*Lam0)
    !         ils2 = 1./(Mrat*Lam1)
    !         t3_vts = Kap0*csg(6)*ils1**cse(6)
    !         t4_vts = Kap1*Mrat**mu_s*csg(12)*ils2**cse(12)
    !         vts_dbz_wt = rhof(k)*av_s * (t1_vts+t2_vts)/(t3_vts+t4_vts)
    !         if (temp(k).ge.273.15 .and. temp(k).lt.275.15) then
    !            vts_dbz_wt = vts_dbz_wt*1.5
    !         elseif (temp(k).ge.275.15) then
    !            vts_dbz_wt = vts_dbz_wt*2.0
    !         endif
    !        else
    !         vts_dbz_wt = 1.E-3
    !        endif

    !        if (rr(k).gt.R1) then
    !         lamr = 1./ilamr(k)
    !         vtr_dbz_wt = rhof(k)*av_r*crg(13)*(lamr+fv_r)**(-cre(13))      &
    !    &               / (crg(4)*lamr**(-cre(4)))
    !        else
    !         vtr_dbz_wt = 1.E-3
    !        endif

    !        if (rg(k).gt.R2) then
    !         lamg = 1./ilamg(k)
    !         vtg_dbz_wt = rhof(k)*av_g*cgg(5)*lamg**(-cge(5))               &
    !    &               / (cgg(4)*lamg**(-cge(4)))
    !        else
    !         vtg_dbz_wt = 1.E-3
    !        endif

    !        vt_dBZ(k) = (vts_dbz_wt*ze_snow(k) + vtr_dbz_wt*ze_rain(k)      &
    !    &                + vtg_dbz_wt*ze_graupel(k))                        &
    !    &                / (ze_rain(k)+ze_snow(k)+ze_graupel(k))
    !     enddo

  END SUBROUTINE calc_refl10cm
  !
  !+---+-----------------------------------------------------------------+

  !+---+-----------------------------------------------------------------+
  !wrf_debug
  !+---+-----------------------------------------------------------------+


  SUBROUTINE wrf_debug(str,str2)
    IMPLICIT NONE
    CHARACTER(LEN=*) , INTENT(IN   ) :: str
    CHARACTER(LEN=*) , INTENT(IN   ) :: str2

    CALL MsgOne(str,str2)
  END  SUBROUTINE wrf_debug

  !+---+-----------------------------------------------------------------+
  !wrf_message
  !+---+-----------------------------------------------------------------+
  SUBROUTINE wrf_message(str2)
    IMPLICIT NONE
    CHARACTER(LEN=*) , INTENT(IN   ) :: str2
    CHARACTER(LEN=21)  ::     str='.....wrf_message.....'
    CALL MsgOne(str,str2)
  END  SUBROUTINE wrf_message
  !+---+-----------------------------------------------------------------+
  !wrf_error_fatal
  !+---+-----------------------------------------------------------------+
  SUBROUTINE wrf_error_fatal(str)
    IMPLICIT NONE
    CHARACTER(LEN=*) , INTENT(IN   ) :: str
    CALL  FatalError(str)
  END  SUBROUTINE wrf_error_fatal


  !+---+-----------------------------------------------------------------+

  !===============================================================================
  SUBROUTINE geopotential_t(                                 &
       piln   ,  pint   , pmid   , pdel   , rpdel  , &
       t      , q      , rair   , gravit , zvir   ,          &
       zi     , zm     , ncol   ,nCols, kMax)

    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Compute the geopotential height (above the surface) at the midpoints and 
    ! interfaces using the input temperatures and pressures.
    !
    !-----------------------------------------------------------------------

    IMPLICIT NONE

    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    INTEGER, INTENT(in) :: ncol                  ! Number of longitudes
    INTEGER, INTENT(in) :: nCols
    INTEGER, INTENT(in) :: kMax
    REAL(r8), INTENT(in) :: piln (nCols,kMax+1)   ! Log interface pressures
    REAL(r8), INTENT(in) :: pint (nCols,kMax+1)   ! Interface pressures
    REAL(r8), INTENT(in) :: pmid (nCols,kMax)    ! Midpoint pressures
    REAL(r8), INTENT(in) :: pdel (nCols,kMax)    ! layer thickness
    REAL(r8), INTENT(in) :: rpdel(nCols,kMax)    ! inverse of layer thickness
    REAL(r8), INTENT(in) :: t    (nCols,kMax)    ! temperature
    REAL(r8), INTENT(in) :: q    (nCols,kMax)    ! specific humidity
    REAL(r8), INTENT(in) :: rair                 ! Gas constant for dry air
    REAL(r8), INTENT(in) :: gravit               ! Acceleration of gravity
    REAL(r8), INTENT(in) :: zvir                 ! rh2o/rair - 1

    ! Output arguments

    REAL(r8), INTENT(out) :: zi(nCols,kMax+1)     ! Height above surface at interfaces
    REAL(r8), INTENT(out) :: zm(nCols,kMax)      ! Geopotential height at mid level
    !
    !---------------------------Local variables-----------------------------
    !
    LOGICAL  :: fvdyn              ! finite volume dynamics
    INTEGER  :: i,k                ! Lon, level indices
    REAL(r8) :: hkk(nCols)         ! diagonal element of hydrostatic matrix
    REAL(r8) :: hkl(nCols)         ! off-diagonal element
    REAL(r8) :: rog                ! Rair / gravit
    REAL(r8) :: tv                 ! virtual temperature
    REAL(r8) :: tvfac              ! Tv/T
    zi= 0.0_r8;    zm= 0.0_r8;    hkk= 0.0_r8;hkl= 0.0_r8
    rog= 0.0_r8;tv = 0.0_r8;tvfac= 0.0_r8
    !
    !-----------------------------------------------------------------------
    !
    rog = rair/gravit

    ! Set dynamics flag

    fvdyn = .FALSE.!dycore_is ('LR')

    ! The surface height is zero by definition.

    DO i = 1,ncol
       zi(i,kMax+1) = 0.0_r8
    END DO

    ! Compute zi, zm from bottom up. 
    ! Note, zi(i,k) is the interface above zm(i,k)

    DO k = kMax, 1, -1

       ! First set hydrostatic elements consistent with dynamics

       IF (fvdyn) THEN
          DO i = 1,ncol
             hkl(i) = piln(i,k+1) - piln(i,k)
             hkk(i) = 1.0_r8 - pint(i,k) * hkl(i) * rpdel(i,k)
          END DO
       ELSE
          DO i = 1,ncol
             hkl(i) = pdel(i,k) / pmid(i,k)
             hkk(i) = 0.5_r8 * hkl(i)
          END DO
       END IF

       ! Now compute tv, zm, zi

       DO i = 1,ncol
          tvfac   = 1.0_r8 + zvir * q(i,k)
          tv      = t(i,k) * tvfac

          zm(i,k) = zi(i,k+1) + rog * tv * hkk(i)
          zi(i,k) = zi(i,k+1) + rog * tv * hkl(i)
       END DO
    END DO

    RETURN
  END SUBROUTINE geopotential_t

END MODULE Micro_GTHOMPSON
!PROGRAM Main
!  USE module_mp_thompson
!END PROGRAM Main
