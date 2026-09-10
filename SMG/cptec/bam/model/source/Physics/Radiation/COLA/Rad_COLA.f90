!
!  $Author: pkubota $
!  $Date: 2011/04/26 17:52:17 $
!  $Revision: 1.4 $
!
MODULE Rad_COLA

  ! InitRadiation
  !
  ! swrad   ------| setsw ------| clear
  !                             |
  !                             | cloudy
  !
  ! lwrad   ------| lwflux -----| crunch
  !               !
  !               ! cldslw

  USE Constants, ONLY :     &
       r8, i8

    IMPLICIT NONE
  SAVE


  PRIVATE

  PUBLIC :: InitRadCOLA
  PUBLIC :: swrad
  PUBLIC :: lwrad

  ! Tabulated planck functions for LW
  REAL(KIND=r8) :: b2502 (32) ! water vapor
  REAL(KIND=r8) :: b2501 (32) ! water vapor
  REAL(KIND=r8) :: blkwin(32) ! ozone
  REAL(KIND=r8) :: blkco2(32) ! co2

  ! Water vapor absorption for SW
  INTEGER, PARAMETER :: nwaterbd = 10 ! number of sw water vapor bands
  REAL(KIND=r8)      :: xk(nwaterbd)  ! coefficients
  REAL(KIND=r8)      :: fk(nwaterbd)  ! weights

CONTAINS

  SUBROUTINE InitRadCOLA()

    ! planck function table for water vapor bands (center)
    b2501(:) = (/ &
         16.280e0_r8, 17.471e0_r8, 18.701e0_r8, 19.974e0_r8, 21.292e0_r8, &
         22.661e0_r8, 24.086e0_r8, 25.575e0_r8, 27.135e0_r8, 28.775e0_r8, &
         30.506e0_r8, 32.339e0_r8, 34.286e0_r8, 36.361e0_r8, 38.578e0_r8, &
         40.954e0_r8, 43.505e0_r8, 46.248e0_r8, 49.203e0_r8, 52.388e0_r8, &
         55.824e0_r8, 59.532e0_r8, 63.533e0_r8, 67.849e0_r8, 72.502e0_r8, &
         77.516e0_r8, 82.913e0_r8, 88.717e0_r8, 94.952e0_r8, 101.64e0_r8, &
         108.806e0_r8,116.472e0_r8 /)

    ! planck function table for water vapor bands (wing)
    b2502(:) = (/ &
         16.379e0_r8, 18.744e0_r8, 21.345e0_r8, 24.195e0_r8,  27.311e0_r8, &
         30.708e0_r8, 34.405e0_r8, 38.417e0_r8, 42.763e0_r8,  47.461e0_r8, &
         52.529e0_r8, 57.985e0_r8, 63.850e0_r8, 70.141e0_r8,  76.880e0_r8, &
         84.088e0_r8, 91.784e0_r8,  99.99e0_r8,108.726e0_r8, 118.016e0_r8, &
         127.881e0_r8,138.344e0_r8,149.429e0_r8,161.160e0_r8, 173.561e0_r8, &
         186.659e0_r8,200.478e0_r8,215.046e0_r8,230.390e0_r8, 246.539e0_r8, &
         263.523e0_r8,281.368e0_r8 /)

    ! planck function table for ozone band.
    blkwin(:) = (/ &
         0.593e0_r8, 0.774e0_r8, 0.993e0_r8, 1.258e0_r8, 1.573e0_r8, &
         1.944e0_r8, 2.377e0_r8, 2.877e0_r8, 3.450e0_r8, 4.102e0_r8, &
         4.838e0_r8, 5.664e0_r8, 6.585e0_r8, 7.606e0_r8, 8.733e0_r8, &
         9.969e0_r8,11.320e0_r8,12.788e0_r8,14.380e0_r8,16.097e0_r8, &
         17.944e0_r8,19.923e0_r8,22.038e0_r8,24.292e0_r8,26.685e0_r8, &
         29.221e0_r8,31.902e0_r8,34.729e0_r8,37.703e0_r8,40.825e0_r8, &
         44.097e0_r8,47.520e0_r8   /)

    ! planck function table for co2 bands.
    blkco2(:) = (/ &
         8.789e0_r8, 10.385e0_r8, 12.159e0_r8, 14.117e0_r8, 16.264e0_r8, &
         18.606e0_r8, 21.145e0_r8, 23.884e0_r8, 26.826e0_r8, 29.973e0_r8, &
         33.325e0_r8, 36.883e0_r8, 40.647e0_r8, 44.617e0_r8, 48.792e0_r8, &
         53.170e0_r8, 57.750e0_r8, 62.530e0_r8, 67.509e0_r8, 72.683e0_r8, &
         78.050e0_r8, 83.609e0_r8, 89.354e0_r8, 95.285e0_r8,101.397e0_r8, &
         107.688e0_r8,114.155e0_r8,120.794e0_r8,127.601e0_r8,134.574e0_r8, &
         141.710e0_r8,149.004e0_r8 /)


    ! Include global version 2.2 - Increase number of spectral radiation bands
    
    ! Ramaswamy &  Friedenreich data

    ! Water vapor absorption coefficients
    xk(:) = (/0.0002e1_r8,0.0035e1_r8,0.0377e1_r8,0.195e1_r8,0.940e1_r8, &
         4.46e1_r8,19.0e1_r8,98.9e1_r8,270.60e1_r8,3901.1e1_r8/)

    ! Water vapor absorption function weights
    fk(:) = (/0.0698e0_r8,0.1558e0_r8,0.0631e0_r8,0.0362e0_r8,0.0243e0_r8, &
         0.0158e0_r8,0.0087e0_r8,0.001467e0_r8,0.002342e0_r8,0.001075e0_r8/)

    ! Original COLA data
    !fk(:) = (/0.107e0_r8, 0.104e0_r8, 0.073e0_r8, 0.044e0_r8,  0.025e0_r8/)
    !xk(:) = (/0.005e0_r8, 0.041e0_r8, 0.416e0_r8, 4.752e0_r8, 72.459e0_r8/)

  END SUBROUTINE InitRadCOLA

  !------------------------------------------------------------------------
  ! LONG WAVE FLUXES CALCULATION
  !    Original Paper:
  !      Harshvardan et al, 1987: "A fast radiation parameterization for 
  !      atmospheric circulation models", J. Geophys. Res., v92, 1009-1016.
  !------------------------------------------------------------------------

  ! crunch: Computation of the gaseous transmission functions
  SUBROUTINE crunch(indx1 ,indx2 ,ncols ,kmax  ,h0p   ,h1p   ,ozone ,txuf  , &
       tv1   ,tv2   ,tui   ,tui2  ,x1    ,x2    ,cc    ,rawi  , &
       x3    ,x4    ,ch    ,css   ,ccu   ,shi   ,shu   ,wdel  , &
       fw    ,pai   ,tai   ,ozai  ,ubar  ,vbar  ,wbar  ,ubarm , &
       vbarm ,wbarm ,fluxu ,fluxd )
    !
    !
    !  Input formal parameters :: index1, index2
    !  parameters used to calculate trnasmission functions
    !  From one level to another
    !==========================================================================
    ! ncols......Number of grid points on a gaussian latitude circle
    ! kmax......Number of grid points at vertical
    ! indx1.....parameters used to calculate trnasmission functions
    ! indx2.....parameters used to calculate trnasmission functions
    ! h0p.......constant h0p = 0.0e0 fac converts to degrees / time step
    ! h1p.......constant h1p = 1.0e0 fac converts to degrees / time step
    ! ozone.....set ozone logical variable  ozone = (.NOT. noz)
    ! txuf......1.used as matrix of g-functions for paths from each level
    !             to all other layers.
    !           2.used for transmission in co2 band.
    !           3.used for transmission in ozone band.
    !           4.in cldslw used for probability of clear line-of-sight
    !             from each level to all other layers for max overlap.
    ! tv1.......Working dimension
    ! tv2 ......Working dimension
    ! tui.......Working dimension
    ! tui2......Working dimension
    ! x1........path water vapor(e-type) and working dimension
    ! x2........path water vapor(band-center) and working dimension
    ! cc........planck function at level temperature for co2 bands.
    ! rawi......water vapor amount in layer.
    ! x3........path water vapor (band-wings) and working dimension
    ! x4........Working dimension
    ! ch........Probability of clear line-of-sight from level to top of
    !           the atmosphere.
    ! css.......Large scale cloud amount and working dimension
    ! ccu.......Cumulus cloud amount and working dimension
    ! shi.......Total transmission function (water vapor + CO2 + ozone)
    !           g-function for a path from level to top of atmosphere.
    ! shu.......Total transmission function (water vapor + CO2 + ozone)
    !           g-function for a path from level  of atmosphere to surface
    ! wdel......Ozone path, water vapor (e-TYPE) transmission function in
    !           9.6 mcm band
    ! fw........Ozone path multiplied bye pressure
    ! pai.......Pressure at middle of layer
    ! tai.......Temperature at middle of layer
    ! ozai......ozone amount in layer.
    ! ubar......scaled water vapor path length in window.
    ! vbar......scaled water vapor path length in center.
    ! wbar......scaled water vapor path length in wing.
    ! ubarm.....ubarm(i,2) = (ubar(i,2) + ubar(i,1)) * hp5
    ! vbarm.... planck function at level temperature for ozone band.
    ! wbarm.... ubarm(i,2) = (ubar(i,2) + ubar(i,1)) * hp5
    ! fluxu.....Ozone path
    ! fluxd.....Ozone path mutiplicated by pressure
    !
    !==========================================================================
    INTEGER, INTENT(IN   ) :: ncols
    INTEGER, INTENT(IN   ) :: kmax
    INTEGER, INTENT(IN   ) :: indx1
    INTEGER, INTENT(IN   ) :: indx2
    REAL(KIND=r8),    INTENT(IN   ) :: h0p
    REAL(KIND=r8),    INTENT(IN   ) :: h1p
    LOGICAL, INTENT(IN   ) :: ozone
    REAL(KIND=r8),    INTENT(INOUT  ) :: txuf  (ncols,kmax+2,kmax+2)
    REAL(KIND=r8),    INTENT(IN   ) :: tv1   (ncols,kmax+3)
    REAL(KIND=r8),    INTENT(IN   ) :: tv2   (ncols,kmax+3)
    REAL(KIND=r8),    INTENT(IN   ) :: tui   (ncols,kmax+3)
    REAL(KIND=r8),    INTENT(IN   ) :: tui2  (ncols,kmax+3)
    REAL(KIND=r8),    INTENT(INOUT  ) :: x1    (ncols,kmax+3)
    REAL(KIND=r8),    INTENT(INOUT  ) :: x2    (ncols,kmax+3)
    REAL(KIND=r8),    INTENT(IN   ) :: cc    (ncols,kmax+3)
    REAL(KIND=r8),    INTENT(IN   ) :: rawi  (ncols,kmax+3)
    REAL(KIND=r8),    INTENT(INOUT  ) :: x3    (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(INOUT  ) :: x4    (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(IN   ) :: ch    (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(IN   ) :: css   (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(IN   ) :: ccu   (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(INOUT  ) :: shi   (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(INOUT  ) :: shu   (ncols,kmax+1)
    REAL(KIND=r8),    INTENT(INOUT  ) :: wdel  (ncols,kmax+1)
    REAL(KIND=r8),    INTENT(INOUT  ) :: fw    (ncols,kmax+1)

    ! Local Variables --->> Global Variables

    REAL(KIND=r8),    INTENT(IN   ) :: pai  (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(IN   ) :: tai  (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(IN   ) :: ozai (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(IN   ) :: ubar (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(IN   ) :: vbar (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(IN   ) :: wbar (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(IN   ) :: ubarm(ncols,kmax+2)
    REAL(KIND=r8),    INTENT(IN   ) :: vbarm(ncols,kmax+2)
    REAL(KIND=r8),    INTENT(IN   ) :: wbarm(ncols,kmax+2)
    REAL(KIND=r8),    INTENT(IN   ) :: fluxu(ncols,kmax+2)
    REAL(KIND=r8),    INTENT(IN   ) :: fluxd(ncols,kmax+2)

    REAL(KIND=r8)                   :: adel  (ncols,kmax+1)
    REAL(KIND=r8)                   :: bdel  (ncols,kmax+1)
    REAL(KIND=r8)                   :: yv    (ncols,kmax+1)
    REAL(KIND=r8)                   :: zv    (ncols,kmax+1)
    REAL(KIND=r8)                   :: wv    (ncols,kmax+1)
    REAL(KIND=r8)                   :: fu    (ncols,kmax+1)
    REAL(KIND=r8)                   :: yw    (ncols,kmax+1)
    REAL(KIND=r8)                   :: zw    (ncols,kmax+1)
    REAL(KIND=r8)                   :: ww    (ncols,kmax+1)

    !  Local Parameter

    REAL(KIND=r8), PARAMETER        :: h9p8  = -9.79e0_r8
    REAL(KIND=r8), PARAMETER        :: h27p  = -27.0e0_r8
    REAL(KIND=r8), PARAMETER        :: hp61  = 6.15384615e-1_r8
    REAL(KIND=r8), PARAMETER        :: h15p1 = 15.1e0_r8
    REAL(KIND=r8), PARAMETER        :: h3p1  = -3.1e0_r8
    REAL(KIND=r8), PARAMETER        :: hp9   = 0.9e0_r8
    REAL(KIND=r8), PARAMETER        :: hp04  = -0.04e0_r8
    REAL(KIND=r8), PARAMETER        :: h16p  = 16.0e0_r8
    REAL(KIND=r8), PARAMETER        :: h6p7  = -6.7e0_r8
    REAL(KIND=r8), PARAMETER        :: h1013 = 1013.25e0_r8
    REAL(KIND=r8), PARAMETER        :: h1381 = 1381.12e0_r8
    REAL(KIND=r8), PARAMETER        :: hp88  = 0.8796e0_r8
    REAL(KIND=r8), PARAMETER        :: hp677 = 6.7675e-1_r8
    REAL(KIND=r8), PARAMETER        :: h4p4  = -4.398e0_r8
    REAL(KIND=r8), PARAMETER        :: hp38  = 3.84615384e-1_r8

    INTEGER                :: ip
    !INTEGER                :: ix
    INTEGER                :: i  ! loop indices
    INTEGER                :: k  ! loop indices
    !INTEGER                :: isub
    !INTEGER                :: ksub
    !INTEGER                :: imtn
    INTEGER                :: ip1
    INTEGER                ::  n
    INTEGER                :: i0
    INTEGER                :: i1
    INTEGER                :: i2

    REAL(KIND=r8)                   :: x1_s
    REAL(KIND=r8)                   :: x2_s
    REAL(KIND=r8)                   :: x3_s
    REAL(KIND=r8)                   :: adel_s
    REAL(KIND=r8)                   :: bdel_s
    REAL(KIND=r8)                   :: wdel_s
    REAL(KIND=r8)                   :: fw_s
    REAL(KIND=r8)                   :: fu_s
    REAL(KIND=r8)                   :: yw_s
    REAL(KIND=r8)                   :: ww_s
    REAL(KIND=r8)                   :: yv_s
    REAL(KIND=r8)                   :: wv_s
    REAL(KIND=r8)                   :: zv_s
    REAL(KIND=r8)                   :: zw_s

    IF (ozone) THEN
       IF (indx1 == indx2) THEN
          IF (indx2 == 1) THEN
             DO k = 1, kMax+1
                DO i = 1, ncols
                   x1(i,k) = ABS(ubar(i,k+1))
                   x2(i,k) = ABS(vbar(i,k+1))
                   x3(i,k) = ABS(wbar(i,k+1))
                   adel(i,k) = ABS(tai(i,k+1))
                   bdel(i,k) = ABS(ch(i,k+1))
                END DO
             END DO
          ELSE
             DO k = 1, kmax+1
                DO i = 1, ncols
                   x1(i,k)   = ubar(i,kmax+2) - ubar(i,k)
                   x2(i,k)   = vbar(i,kmax+2) - vbar(i,k)
                   x3(i,k)   = wbar(i,kmax+2) - wbar(i,k)
                   adel(i,k) = tai(i,kmax+2)  - tai(i,k)
                   bdel(i,k) = ch(i,kmax+2)   - ch(i,k)
                   wdel(i,k) = fluxu(i,kmax+2) - fluxu(i,k)
                   fw(i,k) = fluxd(i,kmax+2) - fluxd(i,k)
                   x1(i,k) = ABS(x1(i,k))
                   x2(i,k) = ABS(x2(i,k))
                   x3(i,k) = ABS(x3(i,k))
                   adel(i,k) = ABS(adel(i,k))
                   bdel(i,k) = ABS(bdel(i,k))
                   wdel(i,k) = ABS(wdel(i,k))
                   fw  (i,k) = ABS(fw  (i,k))
                END DO
             END DO
          END IF

          DO k=1,kmax+1
             DO i = 1, ncols
                fu_s = h9p8 * x1(i,k)
                yw_s = EXP(LOG(x1(i,k) + 1.0e-100_r8)*.83_r8)
                ww_s = EXP(LOG(x3(i,k) + 1.0e-100_r8)*.6_r8)
                ww_s = h1p + h16p  * ww_s
                ww_s = h6p7 * x3(i,k) / ww_s
                ww_s = h27p * yw_s + ww_s
                yw_s = EXP(LOG(adel(i,k) + 1.0e-100_r8)*.56_r8)
                yw_s = h1p  + h15p1     * yw_s
                yv_s = h3p1 * adel(i,k) / yw_s
                yv_s = yv_s       + ww_s
                yw_s = EXP(LOG(bdel(i,k) + 1.0e-100_r8)*.57_r8)
                yw_s = h1p  + hp9       * yw_s
                wv_s = hp04 * bdel(i,k) / yw_s
                wv_s = wv_s + ww_s
                fw_s   = fw(i,k)/ (wdel(i,k) * h1013)
                zv_s   = h1381 * wdel(i,k)  / (hp88 * fw_s)
                adel_s = h1p + zv_s
                adel_s = SQRT(adel_s)
                adel_s = h4p4 * fw_s * (adel_s - h1p)
                zv_s = EXP(adel_s)
                fw_s = h1p - hp677 * (h1p - zv_s)
                adel_s = EXP(yv_s)
                bdel_s = EXP(wv_s)
                wdel_s = EXP(fu_s)
                yw_s = SQRT(x2(i,k))
                zw_s = SQRT(x3(i,k))
                ww_s = ((h1p   + 32.2095e0_r8  * x1(i,k)) &
                     /  (h1p        + 52.85e0_r8      * x1(i,k)) &
                     +  (0.534874e0_r8 + 199.0e0_r8      * x1(i,k) &
                     -  1990.63e0_r8  * x1(i,k)      * x1(i,k)) &
                     *  zw_s &
                     / (h1p         + 333.244e0_r8    * x1(i,k))) &
                     / ((h1p        + 74.144e0_r8     * x1(i,k)) &
                     / (0.43368e0_r8   + 24.7442e0_r8    * x1(i,k)) &
                     * zw_s      + h1p)
                wv_s = (h1p   + 9.22411e0_r8  * yw_s &
                     + 33.1236e0_r8    * x2(i,k) &
                     + 176.396e0_r8    * x2(i,k)      * x2(i,k))
                wv_s = h1p   / wv_s
                ww_s = MAX(ww_s, h0p)
                wv_s = MAX(wv_s, h0p)
                x1_s = MIN(x1(i,k), 0.06e0_r8)
                x2_s = MIN(x2(i,k), 2.0e0_r8)
                x3_s = MIN(x3(i,k), 8.0e0_r8)
                yw_s = SQRT(x2_s)
                zw_s = SQRT(x3_s)
                fu_s = x1_s * x1_s
                yv_s = (0.0851069e0_r8 * yw_s        &
                     - 0.187096e0_r8   * x2_s  * yw_s &
                     + 0.323105e0_r8   * x2_s) * 0.1e0_r8
                zv_s = 0.239186e0_r8   * x2_s        &
                     - 0.0922289e0_r8  * x2_s  * yw_s &
                     - 0.0167413e0_r8  * x2_s  * x2_s
                zv_s = zv_s * 1.0e-3_r8
                yw_s = (5.6383e-4_r8    + 1.05173e0_r8  * x1_s &
                     - 39.0722e0_r8 * fu_s) &
                     / (h1p   + 202.357e0_r8  * x1_s) &
                     + (0.0779555e0_r8  + 4.40720e0_r8  * x1_s &
                     + 3.15851e0_r8 * fu_s)   * zw_s &
                     / (h1p   + 40.2298e0_r8  * x1_s) &
                     + (-0.0381305e0_r8 - 3.63684e0_r8  * x1_s &
                     + 7.98951e0_r8 * fu_s)   * x3_s &
                     / (h1p   + 62.5692e0_r8  * x1_s) &
                     + (6.21039e-3_r8 + 0.710061e0_r8 * x1_s &
                     - 2.85241e0_r8 * fu_s)   * x3_s &
                     / (h1p   + 70.2912e0_r8  * x1_s) &
                     * zw_s
                yw_s = 0.1e0_r8    * yw_s
                zw_s = (-2.99542e-4_r8 + 0.238219e0_r8 * x1_s &
                     + 0.519264e0_r8   * fu_s) &
                     / (h1p         + 10.7775e0_r8  * x1_s) &
                     + (-2.91325e-2_r8 - 2.30007e0_r8  * x1_s &
                     + 10.946e0_r8     * fu_s)   * zw_s &
                     / (h1p         + 63.519e0_r8   * x1_s) &
                     + (1.43812e-2_r8  + 1.80265e0_r8  * x1_s &
                     - 10.1311e0_r8    * fu_s)   * x3_s &
                     / (h1p         + 98.4758e0_r8  * x1_s) &
                     + (-2.39016e-3_r8 - 3.71427e-1_r8 * x1_s &
                     + 2.35443e0_r8    * fu_s)   * x3_s &
                     / (h1p         + 120.228e0_r8  * x1_s) &
                     * zw_s
                zw_s = 1.0e-3_r8   * zw_s
                adel_s = hp38 * adel_s + hp61 * bdel_s
                fw_s   = fw_s * wdel_s
                !
                yw(i,k) = yw_s
                ww(i,k) = ww_s
                yv(i,k) = yv_s
                wv(i,k) = wv_s
                fw(i,k) = fw_s
                zv(i,k) = zv_s
                adel(i,k)= adel_s
                zw(i,k) = zw_s
                x1(i,k) = x1_s
                x2(i,k) = x2_s
                x3(i,k) = x3_s
             END DO
          END DO

          IF (indx2 == 1) THEN
             DO k = 2, kmax+2
                DO i = 1, ncols
                   x1(i,k)   = wv(i,k-1) * tv1(i,1)
                   x2(i,k)   = yv(i,k-1) * tui(i,1) + h1p &
                        + zv(i,k-1)   * tui2(i,1)
                   x3(i,k)   = ww(i,k-1) * tv2(i,1)
                   x4(i,k)   = yw(i,k-1) * tui(i,1) + h1p &
                        + zw(i,k-1)   * tui2(i,1)
                   fw(i,k-1) = adel(i,k-1) * cc(i,1) &
                        + fw(i,k-1)   * rawi(i,1)
                END DO
             END DO

             DO k = 2, kmax+2
                DO i = 1, ncols
                   shi(i,k) = x1(i,k)*x2(i,k) + x3(i,k)*x4(i,k) + fw(i,k-1)
                END DO
             END DO

          ELSE

             DO k = 1, kmax+1
                DO i = 1, ncols
                   x1(i,k) =  wv(i,k) * tv1(i,(kmax+3)) &
                        * (yv(i,k) * tui(i,(kmax+3))  + h1p &
                        + zv(i,k)  * tui2(i,(kmax+3)))
                   x2(i,k) = ww(i,k)  * tv2(i,(kmax+3)) &
                        * (yw(i,k) * tui(i,(kmax+3))  + h1p &
                        +  zw(i,k) * tui2(i,(kmax+3)))
                   x3(i,k) =  wv(i,k) * tv1(i,(kmax+2)) &
                        * (yv(i,k) * tui(i,(kmax+2))    + h1p &
                        + zv(i,k)  * tui2(i,(kmax+2)))
                   x4(i,k) = ww(i,k)  * tv2(i,(kmax+2)) &
                        * (yw(i,k) * tui(i,(kmax+2))    + h1p &
                        +  zw(i,k) * tui2(i,(kmax+2)))
                   shu(i,k) = (cc(i,(kmax+3))-cc(i,(kmax+2)))*adel(i,k) &
                        + (rawi(i,(kmax+3)) - rawi(i,(kmax+2))) * fw(i,k)
                   shu(i,k) = x1(i,k) + x2(i,k) - x3(i,k) - x4(i,k) &
                        + shu(i,k)
                END DO
             END DO

          END IF
       ELSE
          IF (indx2 == (kmax+2)) THEN
             DO ip = indx1, indx2
                DO k = indx1, ip
                   DO i = 1, (ncols)
                      x1(i,k-1)   = ubar(i,ip) - ubarm(i,k)
                      x2(i,k-1)   = vbar(i,ip) - vbarm(i,k)
                      x3(i,k-1)   = wbar(i,ip) - wbarm(i,k)
                      adel(i,k-1) = tai(i,ip)  - css(i,k)
                      bdel(i,k-1) = ch(i,ip)   - ccu(i,k)
                      wdel(i,k-1) = fluxu(i,ip) - ozai(i,k)
                      fw(i,k-1)   = fluxd(i,ip) - pai(i,k)
                   END DO
                END DO
                DO k = 1, ip-1
                   DO i  = 1, ncols
                      x1_s = ABS(x1(i,k))
                      x2_s = ABS(x2(i,k))
                      x3_s = ABS(x3(i,k))
                      adel_s = ABS(adel(i,k))
                      bdel_s = ABS(bdel(i,k))
                      wdel_s = ABS(wdel(i,k))
                      fw_s = ABS(fw(i,k))
                      fu_s = h9p8 * x1_s
                      yw_s = EXP(LOG(x1_s + 1.0e-100_r8)*.83_r8)
                      ww_s = EXP(LOG(x3_s + 1.0e-100_r8)*.6_r8)
                      ww_s = h1p + h16p  * ww_s
                      ww_s = h6p7 * x3_s / ww_s
                      ww_s = h27p * yw_s + ww_s
                      yw_s = EXP(LOG(adel_s + 1.0e-100_r8)*.56_r8)
                      yw_s = h1p  + h15p1     * yw_s
                      yv_s = h3p1 * adel_s / yw_s
                      yv_s = yv_s       + ww_s
                      yw_s = EXP(LOG(bdel_s + 1.0e-100_r8)*.57_r8)
                      yw_s = h1p  + hp9       * yw_s
                      wv_s = hp04 * bdel_s / yw_s
                      wv_s = wv_s + ww_s
                      fw_s = fw_s / (wdel_s * h1013)
                      zv_s = h1381   * wdel_s  / (hp88 * fw_s)
                      adel_s = h1p   + zv_s
                      adel_s = SQRT(adel_s)
                      adel_s = h4p4 * fw_s * (adel_s - h1p)
                      zv_s = EXP(adel_s)
                      fw_s = h1p - hp677 * (h1p - zv_s)
                      adel_s = EXP(yv_s)
                      bdel_s = EXP(wv_s)
                      wdel_s = EXP(fu_s)
                      yw_s = SQRT(x2_s)
                      zw_s = SQRT(x3_s)
                      ww_s = ((h1p   + 32.2095e0_r8  * x1_s) &
                           /  (h1p        + 52.85e0_r8      * x1_s) &
                           +  (0.534874e0_r8 + 199.0e0_r8      * x1_s  &
                           -  1990.63e0_r8   * x1_s      * x1_s) &
                           *  zw_s     &
                           / (h1p         + 333.244e0_r8    * x1_s))&
                           / ((h1p        + 74.144e0_r8     * x1_s) &
                           / (0.43368e0_r8   + 24.7442e0_r8    * x1_s) &
                           * zw_s      + h1p)
                      wv_s = (h1p   + 9.22411e0_r8  * yw_s  &
                           + 33.1236e0_r8    * x2_s   &
                           + 176.396e0_r8    * x2_s      * x2_s)
                      wv_s = h1p   / wv_s
                      ww_s = MAX(ww_s, h0p)
                      wv_s = MAX(wv_s, h0p)
                      x1_s = MIN(x1_s, 0.06e0_r8)
                      x2_s = MIN(x2_s, 2.0e0_r8)
                      x3_s = MIN(x3_s, 8.0e0_r8)
                      yw_s = SQRT(x2_s)
                      zw_s = SQRT(x3_s)
                      fu_s = x1_s * x1_s
                      yv_s = (0.0851069e0_r8  * yw_s &
                           -  0.187096e0_r8 * x2_s  * yw_s &
                           +  0.323105e0_r8 * x2_s) * 0.1e0_r8
                      zv_s =  0.239186e0_r8   * x2_s &
                           -  0.0922289e0_r8  * x2_s  * yw_s &
                           -  0.0167413e0_r8  * x2_s  * x2_s
                      zv_s =  zv_s * 1.0e-3_r8
                      yw_s = (5.6383e-4_r8    + 1.05173e0_r8  * x1_s &
                           - 39.0722e0_r8 * fu_s) &
                           / (h1p   + 202.357e0_r8  * x1_s) &
                           + (0.0779555e0_r8  + 4.40720e0_r8  * x1_s  &
                           + 3.15851e0_r8 * fu_s)   * zw_s  &
                           / (h1p   + 40.2298e0_r8  * x1_s) &
                           + (-0.0381305e0_r8 - 3.63684e0_r8  * x1_s  &
                           + 7.98951e0_r8 * fu_s)   * x3_s  &
                           / (h1p   + 62.5692e0_r8  * x1_s) &
                           + (6.21039e-3_r8 + 0.710061e0_r8 * x1_s  &
                           - 2.85241e0_r8 * fu_s)   * x3_s  &
                           / (h1p   + 70.2912e0_r8  * x1_s) &
                           * zw_s
                      yw_s = 0.1e0_r8    * yw_s
                      zw_s = (-2.99542e-4_r8 + 0.238219e0_r8 * x1_s &
                           + 0.519264e0_r8   * fu_s) &
                           / (h1p         + 10.7775e0_r8  * x1_s) &
                           + (-2.91325e-2_r8 - 2.30007e0_r8  * x1_s &
                           + 10.946e0_r8     * fu_s)   * zw_s &
                           / (h1p         + 63.519e0_r8   * x1_s) &
                           + (1.43812e-2_r8  + 1.80265e0_r8  * x1_s &
                           - 10.1311e0_r8    * fu_s)   * x3_s &
                           / (h1p         + 98.4758e0_r8  * x1_s) &
                           + (-2.39016e-3_r8 - 3.71427e-1_r8 * x1_s &
                           + 2.35443e0_r8    * fu_s)   * x3_s &
                           / (h1p         + 120.228e0_r8  * x1_s) &
                           * zw_s
                      zw_s = 1.0e-3_r8   * zw_s
                      !
                      adel(i,k) = adel_s
                      bdel(i,k) = bdel_s
                      wdel(i,k) = wdel_s
                      fw(i,k)   = fw_s
                      !
                      ! fu(i,k) nao usa depois
                      yw(i,k)   = yw_s
                      ww(i,k)   = ww_s
                      yv(i,k)   = yv_s
                      wv(i,k)   = wv_s
                      zv(i,k)   = zv_s
                      zw(i,k)   = zw_s
                   END DO
                END DO
                DO k = 1, ip-1
                   DO i = 1, ncols
                      x1_s =  wv(i,k)  * tv1 (i,k) &
                           * (yv(i,k)  * tui (i,k)    + h1p &
                           +  zv(i,k)  * tui2(i,k))
                      x2_s =  ww(i,k)  * tv2 (i,k) &
                           * (yw(i,k)  * tui (i,k)    + h1p &
                           +  zw(i,k)  * tui2(i,k))
                      x3_s =  wv(i,k)  * tv1 (i,2+k-1) &
                           * (yv(i,k)  * tui (i,2+k-1) + h1p &
                           +  zv(i,k)  * tui2(i,2+k-1))
                      x4(i,k) =  ww(i,k)  * tv2 (i,2+k-1) &
                           * (yw(i,k)  * tui (i,2+k-1)   + h1p &
                           +  zw(i,k)  * tui2(i,2+k-1))
                      txuf(i,k,ip) =  x1_s  + x2_s &
                           -  x3_s      - x4(i,k) &
                           + (cc(i,1+k-1)    - cc(i,2+k-1)) &
                           * (hp38       * adel(i,k) &
                           +  hp61       * bdel(i,k)) &
                           + (rawi(i,1+k-1)  - rawi(i,2+k-1)) &
                           *  fw(i,k)      * wdel(i,k)
                      !
                      x1(i,k) = x1_s
                      x2(i,k) = x2_s
                      x3(i,k) = x3_s
                   END DO
                END DO
             END DO
          END IF


          IF (indx2 /= (kmax+2)) THEN
             DO ip = indx1, indx2
                DO k = indx1, (kmax+2)-ip
                   DO i = 1, (ncols)
                      x1(i,k-0)   = ubar(i,ip) - ubarm(i,k+ip)
                      x2(i,k-0)   = vbar(i,ip) - vbarm(i,k+ip)
                      x3(i,k-0)   = wbar(i,ip) - wbarm(i,k+ip)
                      adel(i,k-0) = tai(i,ip)  - css(i,k+ip)
                      bdel(i,k-0) = ch(i,ip) - ccu(i,k+ip)
                      wdel(i,k-0) = fluxu(i,ip) - ozai(i,k+ip)
                      fw(i,k-0)   = fluxd(i,ip) - pai(i,k+ip)
                   END DO
                END DO

                DO k = 1, (kmax+2)-ip
                   DO i  = 1, ncols
                      x1(i,k) = ABS(x1(i,k))
                      x2(i,k) = ABS(x2(i,k))
                      x3(i,k) = ABS(x3(i,k))
                      adel(i,k) = ABS(adel(i,k))
                      bdel(i,k) = ABS(bdel(i,k))
                      wdel(i,k) = ABS(wdel(i,k))
                      fw(i,k) = ABS(fw(i,k))
                      fu(i,k) = h9p8 * x1(i,k)
                      yw(i,k) = EXP(LOG(x1(i,k) + 1.0e-100_r8)*.83_r8)
                      ww(i,k) = EXP(LOG(x3(i,k) + 1.0e-100_r8)*.6_r8)
                      ww(i,k) = h1p + h16p  * ww(i,k)
                      ww(i,k) = h6p7 * x3(i,k) / ww(i,k)
                      ww(i,k) = h27p * yw(i,k) + ww(i,k)
                      yw(i,k) = EXP(LOG(adel(i,k) + 1.0e-100_r8)*.56_r8)
                      yw(i,k) = h1p  + h15p1  * yw(i,k)
                      yv(i,k) = h3p1 * adel(i,k) / yw(i,k)
                      yv(i,k) = yv(i,k)   + ww(i,k)
                      yw(i,k) = EXP(LOG(bdel(i,k) + 1.0e-100_r8)*.57_r8)
                      yw(i,k) = h1p  + hp9  * yw(i,k)
                      wv(i,k) = hp04 * bdel(i,k) / yw(i,k)
                      wv(i,k) = wv(i,k) + ww(i,k)
                      fw(i,k)    = fw(i,k) / (wdel(i,k) * h1013)
                      zv(i,k)    = h1381   * wdel(i,k)  / (hp88 * fw(i,k))
                      adel(i,k) = h1p      + zv(i,k)
                      adel(i,k) = SQRT(adel(i,k))
                      adel(i,k) = h4p4 * fw(i,k) * (adel(i,k) - h1p)
                      zv(i,k) = EXP(adel(i,k))
                      fw(i,k) = h1p - hp677 * (h1p - zv(i,k))
                      adel(i,k) = EXP(yv(i,k))
                      bdel(i,k) = EXP(wv(i,k))
                      wdel(i,k) = EXP(fu(i,k))
                      yw(i,k) = SQRT(x2(i,k))
                      zw(i,k) = SQRT(x3(i,k))
                      ww(i,k) = ((h1p      + 32.2095e0_r8    * x1(i,k)) &
                           /  (h1p   + 52.85e0_r8  * x1(i,k)) &
                           +  (0.534874e0_r8 + 199.0e0_r8  * x1(i,k)  &
                           -  1990.63e0_r8   * x1(i,k)  * x1(i,k)) &
                           *  zw(i,k)        &
                           / (h1p   + 333.244e0_r8  * x1(i,k)))&
                           / ((h1p   + 74.144e0_r8  * x1(i,k)) &
                           / (0.43368e0_r8   + 24.7442e0_r8  * x1(i,k)) &
                           * zw(i,k)   + h1p)
                      wv(i,k) = (h1p      + 9.22411e0_r8    * yw(i,k)  &
                           + 33.1236e0_r8    * x2(i,k)      &
                           + 176.396e0_r8    * x2(i,k)  * x2(i,k))
                      wv(i,k) = h1p      / wv(i,k)
                      ww(i,k) = MAX(ww(i,k), h0p)
                      wv(i,k) = MAX(wv(i,k), h0p)
                      x1(i,k) = MIN(x1(i,k), 0.06e0_r8)
                      x2(i,k) = MIN(x2(i,k), 2.0e0_r8)
                      x3(i,k) = MIN(x3(i,k), 8.0e0_r8)
                      yw(i,k) = SQRT(x2(i,k))
                      zw(i,k) = SQRT(x3(i,k))
                      fu(i,k) = x1(i,k) * x1(i,k)
                      yv(i,k) = (0.0851069e0_r8  * yw(i,k) &
                           -  0.187096e0_r8   * x2(i,k)  * yw(i,k) &
                           +  0.323105e0_r8   * x2(i,k)) * 0.1e0_r8
                      zv(i,k) =  0.239186e0_r8   * x2(i,k) &
                           -  0.0922289e0_r8  * x2(i,k)  * yw(i,k) &
                           -  0.0167413e0_r8  * x2(i,k)  * x2(i,k)
                      zv(i,k) =  zv(i,k) * 1.0e-3_r8

                      yw(i,k) = (5.6383e-4_r8    + 1.05173e0_r8  * x1(i,k) &
                           - 39.0722e0_r8     * fu(i,k)) &
                           / (h1p    + 202.357e0_r8  * x1(i,k)) &
                           + (0.0779555e0_r8  + 4.40720e0_r8  * x1(i,k)  &
                           + 3.15851e0_r8     * fu(i,k)) * zw(i,k)  &
                           / (h1p    + 40.2298e0_r8  * x1(i,k)) &
                           + (-0.0381305e0_r8 - 3.63684e0_r8  * x1(i,k)  &
                           + 7.98951e0_r8     * fu(i,k)) * x3(i,k)  &
                           / (h1p    + 62.5692e0_r8  * x1(i,k)) &
                           + (6.21039e-3_r8   + 0.710061e0_r8 * x1(i,k)  &
                           - 2.85241e0_r8     * fu(i,k)) * x3(i,k)  &
                           / (h1p    + 70.2912e0_r8  * x1(i,k)) &
                           * zw(i,k)
                      yw(i,k) = 0.1e0_r8       * yw(i,k)
                      zw(i,k) = (-2.99542e-4_r8 + 0.238219e0_r8 * x1(i,k) &
                           + 0.519264e0_r8   * fu(i,k)) &
                           / (h1p   + 10.7775e0_r8  * x1(i,k)) &
                           + (-2.91325e-2 - 2.30007e0_r8  * x1(i,k) &
                           + 10.946e0_r8   * fu(i,k))   * zw(i,k) &
                           / (h1p   + 63.519e0_r8   * x1(i,k)) &
                           + (1.43812e-2_r8  + 1.80265e0_r8  * x1(i,k) &
                           - 10.1311e0_r8    * fu(i,k))   * x3(i,k) &
                           / (h1p   + 98.4758e0_r8  * x1(i,k)) &
                           + (-2.39016e-3_r8 - 3.71427e-1_r8 * x1(i,k) &
                           + 2.35443e0_r8    * fu(i,k))   * x3(i,k) &
                           / (h1p   + 120.228e0_r8  * x1(i,k)) &
                           * zw(i,k)
                      zw(i,k) = 1.0e-3_r8       * zw(i,k)
                   END DO
                END DO

                DO k = 1, (kmax+2)-ip
                   DO i = 1, ncols
                      x1(i,k) =  wv(i,k)  * tv1 (i,ip+k) &
                           * (yv(i,k)  * tui (i,ip+k)  + h1p &
                           +  zv(i,k)  * tui2(i,ip+k))
                      x2(i,k) =  ww(i,k)  * tv2 (i,ip+k) &
                           * (yw(i,k)  * tui (i,ip+k)  + h1p &
                           +  zw(i,k)  * tui2(i,ip+k))
                      x3(i,k) =  wv(i,k)  * tv1 (i,ip+0+k-1) &
                           * (yv(i,k)  * tui (i,ip+0+k-1) + h1p &
                           +  zv(i,k)  * tui2(i,ip+0+k-1))
                      x4(i,k) =  ww(i,k)  * tv2 (i,ip+0+k-1) &
                           * (yw(i,k)  * tui (i,ip+0+k-1)   + h1p &
                           +  zw(i,k)  * tui2(i,ip+0+k-1))
                      txuf(i,ip+k,ip) =  x1(i,k)    + x2(i,k) &
                           -  x3(i,k)  - x4(i,k) &
                           + (cc(i,ip+1+k-1) - cc(i,ip+0+k-1)) &
                           * (hp38  * adel(i,k) &
                           +  hp61  * bdel(i,k)) &
                           + (rawi(i,ip+1+k-1)  - rawi(i,ip+0+k-1)) &
                           *  fw(i,k)  * wdel(i,k)
                   END DO
                END DO
             END DO
          END IF

       END IF
    ENDIF

    IF (.not.ozone) THEN
       IF (indx1 == indx2) THEN
          IF (indx2 == 1) THEN
             DO k = 1, kMax+1
                DO i = 1, ncols
                   x1(i,k) = ABS(ubar(i,k+1))
                   x2(i,k) = ABS(vbar(i,k+1))
                   x3(i,k) = ABS(wbar(i,k+1))
                   adel(i,k) = ABS(tai(i,k+1))
                   bdel(i,k) = ABS(ch(i,k+1))
                END DO
             END DO
          ELSE
             DO k = 1, kmax+1
                DO i = 1, ncols
                   x1(i,k)   = ubar(i,kmax+2) - ubar(i,k)
                   x2(i,k)   = vbar(i,kmax+2) - vbar(i,k)
                   x3(i,k)   = wbar(i,kmax+2) - wbar(i,k)
                   adel(i,k) = tai(i,kmax+2)  - tai(i,k)
                   bdel(i,k) = ch(i,kmax+2)   - ch(i,k)
                END DO
             END DO

             DO k = 1, kMax+1
                DO i = 1, ncols
                   x1(i,k) = ABS(x1(i,k))
                   x2(i,k) = ABS(x2(i,k))
                   x3(i,k) = ABS(x3(i,k))
                   adel(i,k) = ABS(adel(i,k))
                   bdel(i,k) = ABS(bdel(i,k))
                END DO
             END DO


          END IF
          DO k=1,(kmax+1)
             DO i = 1, (ncols)
                fu(i,k) = h9p8 * x1(i,k)
                yw(i,k) = EXP(LOG(x1(i,k) + 1.0e-100_r8)*.83_r8)
                ww(i,k) = EXP(LOG(x3(i,k) + 1.0e-100_r8)*.6_r8)
                ww(i,k) = h1p + h16p  * ww(i,k)
                ww(i,k) = h6p7 * x3(i,k) / ww(i,k)
                ww(i,k) = h27p * yw(i,k) + ww(i,k)
                yw(i,k) = EXP(LOG(adel(i,k) + 1.0e-100_r8)*.56_r8)
                yw(i,k) = h1p  + h15p1     * yw(i,k)
                yv(i,k) = h3p1 * adel(i,k) / yw(i,k)
                yv(i,k) = yv(i,k)       + ww(i,k)
                yw(i,k) = EXP(LOG(bdel(i,k) + 1.0e-100_r8)*.57_r8)
                yw(i,k) = h1p  + hp9       * yw(i,k)
                wv(i,k) = hp04 * bdel(i,k) / yw(i,k)
                wv(i,k) = wv(i,k) + ww(i,k)
             END DO
          END DO

          DO k=1,(kmax+1)
             DO i = 1, (ncols)
                fw(i,k) = h1p
             END DO
          END DO


          DO k = 1, kMax+1
             DO i = 1, ncols
                adel(i,k) = EXP(yv(i,k))
                bdel(i,k) = EXP(wv(i,k))
                wdel(i,k) = EXP(fu(i,k))
             END DO
          END DO

          DO k = 1, kMax+1
             DO i = 1, ncols
                yw(i,k) = SQRT(x2(i,k))
                zw(i,k) = SQRT(x3(i,k))
             END DO
          END DO

          DO k=1,(kmax+1)
             DO i = 1, (ncols)
                ww(i,k) = ((h1p   + 32.2095e0_r8  * x1(i,k)) &
                     /  (h1p        + 52.85e0_r8      * x1(i,k)) &
                     +  (0.534874e0_r8 + 199.0e0_r8      * x1(i,k) &
                     -  1990.63e0_r8   * x1(i,k)      * x1(i,k)) &
                     *  zw(i,k) &
                     / (h1p         + 333.244e0_r8    * x1(i,k))) &
                     / ((h1p        + 74.144e0_r8     * x1(i,k)) &
                     / (0.43368e0_r8   + 24.7442e0_r8    * x1(i,k)) &
                     * zw(i,k)      + h1p)
                wv(i,k) = (h1p   + 9.22411e0_r8  * yw(i,k) &
                     + 33.1236e0_r8    * x2(i,k) &
                     + 176.396e0_r8    * x2(i,k)      * x2(i,k))
                wv(i,k) = h1p   / wv(i,k)
             END DO
          END DO

          DO k = 1, kMax+1
             DO i = 1, ncols
                ww(i,k) = MAX(ww(i,k), h0p)
                wv(i,k) = MAX(wv(i,k), h0p)
                x1(i,k) = MIN(x1(i,k), 0.06e0_r8)
                x2(i,k) = MIN(x2(i,k), 2.0e0_r8)
                x3(i,k) = MIN(x3(i,k), 8.0e0_r8)
                yw(i,k) = SQRT(x2(i,k))
                zw(i,k) = SQRT(x3(i,k))
             END DO
          END DO
          DO k=1,(kmax+1)
             DO i = 1, (ncols)
                fu(i,k) = x1(i,k) * x1(i,k)
                yv(i,k) = (0.0851069e0_r8 * yw(i,k)        &
                     - 0.187096e0_r8   * x2(i,k)  * yw(i,k) &
                     + 0.323105e0_r8   * x2(i,k)) * 0.1e0_r8
                zv(i,k) = 0.239186e0_r8   * x2(i,k)        &
                     - 0.0922289e0_r8  * x2(i,k)  * yw(i,k) &
                     - 0.0167413e0_r8  * x2(i,k)  * x2(i,k)
                zv(i,k) = zv(i,k) * 1.0e-3_r8
                yw(i,k) = (5.6383e-4_r8    + 1.05173e0_r8  * x1(i,k) &
                     - 39.0722e0_r8 * fu(i,k)) &
                     / (h1p   + 202.357e0_r8  * x1(i,k)) &
                     + (0.0779555e0_r8  + 4.40720e0_r8  * x1(i,k) &
                     + 3.15851e0_r8 * fu(i,k))   * zw(i,k) &
                     / (h1p   + 40.2298e0_r8  * x1(i,k)) &
                     + (-0.0381305e0_r8 - 3.63684e0_r8  * x1(i,k) &
                     + 7.98951e0_r8 * fu(i,k))   * x3(i,k) &
                     / (h1p   + 62.5692e0_r8  * x1(i,k)) &
                     + (6.21039e-3_r8 + 0.710061e0_r8 * x1(i,k) &
                     - 2.85241e0_r8 * fu(i,k))   * x3(i,k) &
                     / (h1p   + 70.2912e0_r8  * x1(i,k)) &
                     * zw(i,k)
                yw(i,k) = 0.1e0_r8    * yw(i,k)
                zw(i,k) = (-2.99542e-4_r8 + 0.238219e0_r8 * x1(i,k) &
                     + 0.519264e0_r8   * fu(i,k)) &
                     / (h1p         + 10.7775e0_r8  * x1(i,k)) &
                     + (-2.91325e-2_r8 - 2.30007e0_r8  * x1(i,k) &
                     + 10.946e0_r8     * fu(i,k))   * zw(i,k) &
                     / (h1p         + 63.519e0_r8   * x1(i,k)) &
                     + (1.43812e-2_r8  + 1.80265e0_r8  * x1(i,k) &
                     - 10.1311e0_r8    * fu(i,k))   * x3(i,k) &
                     / (h1p         + 98.4758e0_r8  * x1(i,k)) &
                     + (-2.39016e-3_r8 - 3.71427e-1_r8 * x1(i,k) &
                     + 2.35443e0_r8    * fu(i,k))   * x3(i,k) &
                     / (h1p         + 120.228e0_r8  * x1(i,k)) &
                     * zw(i,k)
                zw(i,k) = 1.0e-3_r8   * zw(i,k)
             END DO
          END DO
          DO k=1,(kmax+1)
             DO i = 1, (ncols)
                adel(i,k) = hp38 * adel(i,k) + hp61 * bdel(i,k)
                fw(i,k)   = fw(i,k) * wdel(i,k)
             END DO
          END DO

          IF (indx2 == 1) THEN

             DO k = 2, (kmax+2)
                DO i = 1, (ncols)
                   x1(i,k)   = wv(i,k-1) * tv1(i,1)
                   x2(i,k)   = yv(i,k-1) * tui(i,1) + h1p &
                        + zv(i,k-1)   * tui2(i,1)
                   x3(i,k)   = ww(i,k-1) * tv2(i,1)
                   x4(i,k)   = yw(i,k-1) * tui(i,1) + h1p &
                        + zw(i,k-1)   * tui2(i,1)
                   fw(i,k-1) = adel(i,k-1) * cc(i,1) &
                        + fw(i,k-1)   * rawi(i,1)
                END DO
             END DO

             DO k = 2, (kmax+2)
                DO i = 1, (ncols)
                   shi(i,k) = x1(i,k)*x2(i,k) + x3(i,k)*x4(i,k) + fw(i,k-1)
                END DO
             END DO

          ELSE

             DO k = 1, (kmax+1)
                DO i = 1, (ncols)
                   x1(i,k) =  wv(i,k) * tv1(i,(kmax+3)) &
                        * (yv(i,k) * tui(i,(kmax+3))  + h1p &
                        + zv(i,k)  * tui2(i,(kmax+3)))
                   x2(i,k) = ww(i,k)  * tv2(i,(kmax+3)) &
                        * (yw(i,k) * tui(i,(kmax+3))  + h1p &
                        +  zw(i,k) * tui2(i,(kmax+3)))
                   x3(i,k) =  wv(i,k) * tv1(i,(kmax+2)) &
                        * (yv(i,k) * tui(i,(kmax+2))    + h1p &
                        + zv(i,k)  * tui2(i,(kmax+2)))
                   x4(i,k) = ww(i,k)  * tv2(i,(kmax+2)) &
                        * (yw(i,k) * tui(i,(kmax+2))    + h1p &
                        +  zw(i,k) * tui2(i,(kmax+2)))
                   shu(i,k) = (cc(i,(kmax+3))-cc(i,(kmax+2)))*adel(i,k) &
                        + (rawi(i,(kmax+3)) - rawi(i,(kmax+2))) * fw(i,k)
                END DO
             END DO
             DO k=1,(kmax+1)
                DO i = 1, (ncols)
                   shu(i,k) = x1(i,k) + x2(i,k) - x3(i,k) - x4(i,k) &
                        + shu(i,k)
                END DO
             END DO

          END IF
       ELSE
          DO ip = indx1, indx2

             IF (indx2 == (kmax+2)) THEN
                ip1 = ip
                n   = ip - 1
                i0  = 0
                i1  = 1
                i2  = 2
             ELSE
                ip1 = (kmax+2) - ip
                n   = ip1
                i0  = ip
                i1  = 0
                i2  = 0
             END IF

             DO k = indx1, ip1
                DO i = 1, (ncols)
                   x1(i,k-i1)   = ubar(i,ip) - ubarm(i,k+i0)
                   x2(i,k-i1)   = vbar(i,ip) - vbarm(i,k+i0)
                   x3(i,k-i1)   = wbar(i,ip) - wbarm(i,k+i0)
                   adel(i,k-i1) = tai(i,ip)  - css(i,k+i0)
                   bdel(i,k-i1) = ch(i,ip)   - ccu(i,k+i0)
                END DO
             END DO


             DO k = 1, n
                DO i  = 1, ncols
                   x1(i,k) = ABS(x1(i,k))
                   x2(i,k) = ABS(x2(i,k))
                   x3(i,k) = ABS(x3(i,k))
                   adel(i,k) = ABS(adel(i,k))
                   bdel(i,k) = ABS(bdel(i,k))
                END DO
             END DO

             DO k= 1, n
                DO i= 1, ncols
                   fu(i,k) = h9p8 * x1(i,k)
                   yw(i,k) = EXP(LOG(x1(i,k) + 1.0e-100_r8)*.83_r8)
                   ww(i,k) = EXP(LOG(x3(i,k) + 1.0e-100_r8)*.6_r8)
                   ww(i,k) = h1p + h16p  * ww(i,k)
                   ww(i,k) = h6p7 * x3(i,k) / ww(i,k)
                   ww(i,k) = h27p * yw(i,k) + ww(i,k)
                   yw(i,k) = EXP(LOG(adel(i,k) + 1.0e-100_r8)*.56_r8)
                   yw(i,k) = h1p  + h15p1     * yw(i,k)
                   yv(i,k) = h3p1 * adel(i,k) / yw(i,k)
                   yv(i,k) = yv(i,k)       + ww(i,k)
                   yw(i,k) = EXP(LOG(bdel(i,k) + 1.0e-100_r8)*.57_r8)
                   yw(i,k) = h1p  + hp9       * yw(i,k)
                   wv(i,k) = hp04 * bdel(i,k) / yw(i,k)
                   wv(i,k) = wv(i,k) + ww(i,k)
                END DO
             END DO


             DO k=1, n
                DO i = 1, ncols
                   fw(i,k) = h1p
                END DO
             END DO

             DO k = 1, n
                DO i = 1, ncols
                   adel(i,k) = EXP(yv(i,k))
                   bdel(i,k) = EXP(wv(i,k))
                   wdel(i,k) = EXP(fu(i,k))
                END DO
             END DO
             DO k = 1, n
                DO i = 1, ncols
                   yw(i,k) = SQRT(x2(i,k))
                   zw(i,k) = SQRT(x3(i,k))
                END DO
             END DO

             DO k = 1, n
                DO i = 1, ncols
                   ww(i,k) = ((h1p   + 32.2095e0_r8  * x1(i,k)) &
                        /  (h1p        + 52.85e0_r8      * x1(i,k)) &
                        +  (0.534874e0_r8 + 199.0e0_r8      * x1(i,k)  &
                        -  1990.63e0_r8   * x1(i,k)      * x1(i,k)) &
                        *  zw(i,k)     &
                        / (h1p         + 333.244e0_r8    * x1(i,k)))&
                        / ((h1p        + 74.144e0_r8     * x1(i,k)) &
                        / (0.43368e0_r8   + 24.7442e0_r8    * x1(i,k)) &
                        * zw(i,k)      + h1p)
                   wv(i,k) = (h1p   + 9.22411e0_r8  * yw(i,k)  &
                        + 33.1236e0_r8    * x2(i,k)   &
                        + 176.396e0_r8    * x2(i,k)      * x2(i,k))
                   wv(i,k) = h1p   / wv(i,k)
                END DO
             END DO

             DO k = 1, n
                DO i = 1, ncols
                   ww(i,k) = MAX(ww(i,k), h0p)
                   wv(i,k) = MAX(wv(i,k), h0p)
                   x1(i,k) = MIN(x1(i,k), 0.06e0_r8)
                   x2(i,k) = MIN(x2(i,k), 2.0e0_r8)
                   x3(i,k) = MIN(x3(i,k), 8.0e0_r8)
                   yw(i,k) = SQRT(x2(i,k))
                   zw(i,k) = SQRT(x3(i,k))
                END DO
             END DO

             DO k =1, n
                DO i = 1, ncols
                   fu(i,k) = x1(i,k) * x1(i,k)
                   yv(i,k) = (0.0851069e0_r8  * yw(i,k) &
                        -  0.187096e0_r8 * x2(i,k)  * yw(i,k) &
                        +  0.323105e0_r8 * x2(i,k)) * 0.1e0_r8
                   zv(i,k) =  0.239186e0_r8   * x2(i,k) &
                        -  0.0922289e0_r8  * x2(i,k)  * yw(i,k) &
                        -  0.0167413e0_r8  * x2(i,k)  * x2(i,k)
                   zv(i,k) =  zv(i,k) * 1.0e-3_r8

                   yw(i,k) = (5.6383e-4_r8    + 1.05173e0_r8  * x1(i,k) &
                        - 39.0722e0_r8 * fu(i,k)) &
                        / (h1p   + 202.357e0_r8  * x1(i,k)) &
                        + (0.0779555e0_r8  + 4.40720e0_r8  * x1(i,k)  &
                        + 3.15851e0_r8 * fu(i,k))   * zw(i,k)  &
                        / (h1p   + 40.2298e0_r8  * x1(i,k)) &
                        + (-0.0381305e0_r8 - 3.63684e0_r8  * x1(i,k)  &
                        + 7.98951e0_r8 * fu(i,k))   * x3(i,k)  &
                        / (h1p   + 62.5692e0_r8  * x1(i,k)) &
                        + (6.21039e-3_r8 + 0.710061e0_r8 * x1(i,k)  &
                        - 2.85241e0_r8 * fu(i,k))   * x3(i,k)  &
                        / (h1p   + 70.2912e0_r8  * x1(i,k)) &
                        * zw(i,k)
                   yw(i,k) = 0.1e0_r8    * yw(i,k)
                   zw(i,k) = (-2.99542e-4_r8 + 0.238219e0_r8 * x1(i,k) &
                        + 0.519264e0_r8   * fu(i,k)) &
                        / (h1p         + 10.7775e0_r8  * x1(i,k)) &
                        + (-2.91325e-2_r8 - 2.30007e0_r8  * x1(i,k) &
                        + 10.946e0_r8     * fu(i,k))   * zw(i,k) &
                        / (h1p         + 63.519e0_r8   * x1(i,k)) &
                        + (1.43812e-2_r8  + 1.80265e0_r8  * x1(i,k) &
                        - 10.1311e0_r8    * fu(i,k))   * x3(i,k) &
                        / (h1p         + 98.4758e0_r8  * x1(i,k)) &
                        + (-2.39016e-3_r8 - 3.71427e-1_r8 * x1(i,k) &
                        + 2.35443e0_r8    * fu(i,k))   * x3(i,k) &
                        / (h1p         + 120.228e0_r8  * x1(i,k)) &
                        * zw(i,k)
                   zw(i,k) = 1.0e-3_r8   * zw(i,k)
                END DO
             END DO

             DO k = 1, n
                DO i = 1, ncols
                   x1(i,k) =  wv(i,k)  * tv1 (i,i0+k) &
                        * (yv(i,k)  * tui (i,i0+k)    + h1p &
                        +  zv(i,k)  * tui2(i,i0+k))
                   x2(i,k) =  ww(i,k)  * tv2 (i,i0+k) &
                        * (yw(i,k)  * tui (i,i0+k)    + h1p &
                        +  zw(i,k)  * tui2(i,i0+k))
                   x3(i,k) =  wv(i,k)  * tv1 (i,i0+i2+k-1) &
                        * (yv(i,k)  * tui (i,i0+i2+k-1) + h1p &
                        +  zv(i,k)  * tui2(i,i0+i2+k-1))
                   x4(i,k) =  ww(i,k)  * tv2 (i,i0+i2+k-1) &
                        * (yw(i,k)  * tui (i,i0+i2+k-1)   + h1p &
                        +  zw(i,k)  * tui2(i,i0+i2+k-1))
                   txuf(i,i0+k,ip) =  x1(i,k)  + x2(i,k) &
                        -  x3(i,k)      - x4(i,k) &
                        + (cc(i,i0+1+k-1)    - cc(i,i0+i2+k-1)) &
                        * (hp38       * adel(i,k) &
                        +  hp61       * bdel(i,k)) &
                        + (rawi(i,i0+1+k-1)  - rawi(i,i0+i2+k-1)) &
                        *  fw(i,k)      * wdel(i,k)
                END DO
             END DO
          END DO
       END IF
    END IF
  END SUBROUTINE crunch


  ! lwflux :computes carbon dioxide and ozone fluxes and upward and downward
  !         fluxes txuf for water vapor; the transmittances are also
  !         calculated.
  SUBROUTINE lwflux(pai   ,tai   ,ozai  ,ubar  ,vbar  ,wbar  ,ubarm ,vbarm , &
       wbarm ,fluxu ,fluxd ,txuf  ,tv1   ,tv2   ,tui   ,x1    , &
       x2    ,cc    ,rawi  ,x3    ,x4    ,ch    ,dp    ,css   , &
       ccu   ,shi   ,shh   ,shu   ,sumsav,h0p   ,h1p   ,h1p5  , &
       hp5   ,dtb   ,dtbinv,pr    ,ntm1  ,ozone ,co2m, &
       ncols ,kmax  )
    !
    !==========================================================================
    !
    !   ncols......Number of grid points on a gaussian latitude circle
    !   kmax......Number of grid points at vertical
    !   co2m....co2val is wgne standard value in ppm "co2val = /345.0/
    !   h0p.......constant h0p = 0.0e0 fac converts to degrees / time step
    !   h1p.......constant h1p = 1.0e0 fac converts to degrees / time step
    !   h1p5......Fact converts absorption to rate in degrees/ time step
    !             constant h1p5   = 1.5e0
    !   hp5.......constant hp5    = 0.5e0
    !   dtb.......temperature increment in b250.  Constant dtb  = 5.0e0
    !   dtbinv....constant dtbinv = h1p / dtb
    !   pr(1).....constant pr(1)  = h1p / 3.0e01
    !   pr(1).....constant pr(2)  = h1p / 3.0e02
    !   ntm1......number of rows in b250 - 1. constant  ntm1   = 31
    !   ozone.....set ozone logical variable  ozone = (.NOT. noz)
    !             true if there is ozone  absorption computation
    !   txuf......1.used as matrix of g-functions for paths from each level
    !               to all other layers.
    !             2.used for transmission in co2 band.
    !             3.used for transmission in ozone band.
    !             4.in cldslw used for probability of clear line-of-sight
    !               from each level to all other layers for max overlap.
    !   tv1.......Working dimension
    !   tv2 ......Working dimension
    !   tui.......Working dimension
    !   x1........path water vapor(e-type) and working dimension
    !   x2........path water vapor(band-center) and working dimension
    !   cc........planck function at level temperature for co2 bands.
    !   rawi......water vapor amount in layer
    !   x3........path water vapor (band-wings) and working dimension
    !   x4........Working dimension
    !   ch........Probability of clear line-of-sight from level to top of
    !             the atmosphere.
    !   dp........Pressure difference between levels
    !   css.......Large scale cloud amount and working dimension
    !   ccu.......Cumulus cloud amount and working dimension
    !   shi.......Total transmission function (water vapor + CO2 + ozone)
    !             g-function for a path from level to top of atmosphere.
    !   shh.......planck function at level temperature for water vapor
    !              bands.
    !   shu.......Total transmission function (water vapor + CO2 + ozone)
    !             g-function for a path from level  of atmosphere to surface
    !   sumsav....
    !   pai.......Pressure at middle of layer
    !   tai.......Temperature at middle of layer
    !   ozai......ozone amount in layer.
    !   ubar......scaled water vapor path length in window.
    !   vbar......scaled water vapor path length in center.
    !   wbar......scaled water vapor path length in wing.
    !   ubarm.....ubarm(i,2) = (ubar(i,2) + ubar(i,1)) * hp5
    !   vbarm.... planck function at level temperature for ozone band.
    !   wbarm.... ubarm(i,2) = (ubar(i,2) + ubar(i,1)) * hp5
    !   fluxu.....Ozone path
    !   fluxd.....Ozone path mutiplicated by pressure
    !
    !==========================================================================
    INTEGER, INTENT(in   ) :: ncols
    INTEGER, INTENT(in   ) :: kmax
    REAL(KIND=r8),    INTENT(in   ) :: h0p
    REAL(KIND=r8),    INTENT(in   ) :: h1p
    REAL(KIND=r8),    INTENT(in   ) :: h1p5
    REAL(KIND=r8),    INTENT(in   ) :: hp5
    REAL(KIND=r8),    INTENT(in   ) :: dtb
    REAL(KIND=r8),    INTENT(in   ) :: dtbinv
    REAL(KIND=r8),    INTENT(in   ) :: pr(2)
    INTEGER, INTENT(in   ) :: ntm1
    LOGICAL, INTENT(in   ) :: ozone
    REAL(KIND=r8),    INTENT(inout  ) :: txuf  (ncols,kmax+2,kmax+2)
    REAL(KIND=r8),    INTENT(inout  ) :: tv1   (ncols,kmax+3)
    REAL(KIND=r8),    INTENT(inout  ) :: tv2   (ncols,kmax+3)
    REAL(KIND=r8),    INTENT(inout) :: tui   (ncols,kmax+3)
    REAL(KIND=r8),    INTENT(inout  ) :: x1    (ncols,kmax+3)
    REAL(KIND=r8),    INTENT(inout  ) :: x2    (ncols,kmax+3)
    REAL(KIND=r8),    INTENT(inout  ) :: cc    (ncols,kmax+3)
    REAL(KIND=r8),    INTENT(inout  ) :: rawi  (ncols,kmax+3)
    REAL(KIND=r8),    INTENT(inout  ) :: x3    (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(inout  ) :: x4    (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(inout  ) :: ch    (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(in   ) :: dp    (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(inout  ) :: css   (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(inout  ) :: ccu   (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(inout  ) :: shi   (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(inout  ) :: shh   (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(inout  ) :: shu   (ncols,kmax+1)
    REAL(KIND=r8),    INTENT(inout  ) :: sumsav(ncols)

    REAL(KIND=r8),    INTENT(INOUT) :: pai   (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(INOUT) :: tai   (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(INOUT) :: ozai  (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(IN   ) :: ubar  (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(IN   ) :: vbar  (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(IN   ) :: wbar  (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(IN   ) :: ubarm (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(IN   ) :: vbarm (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(IN   ) :: wbarm (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(INOUT) :: fluxu (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(INOUT) :: fluxd (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(in   ) :: co2m  (nCols,kMax) !mol/mol

    REAL(KIND=r8)                   :: auxco2  (ncols,kmax+2)
    REAL(KIND=r8)                   :: tui2  (ncols,kmax+3)
    REAL(KIND=r8)                   :: wdel  (ncols,kmax+1)
    REAL(KIND=r8)                   :: fw   (ncols,kmax+1)
    INTEGER                :: it    (ncols,kmax+3)


    !
    !     REAL(KIND=r8) constants
    !
    REAL(KIND=r8), PARAMETER :: hp25=0.26_r8
    REAL(KIND=r8), PARAMETER :: temp1=165.0_r8
    REAL(KIND=r8), PARAMETER :: h250p=250.0_r8
    REAL(KIND=r8), PARAMETER :: pp=1.0_r8

    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER :: ix
    INTEGER :: ip
    INTEGER :: indx1
    INTEGER :: indx2

    DO i = 1, ncols
       shi(i,1) = h0p
    END DO
    DO k=1,kmax
       DO i = 1, ncols
          auxco2(i,k) = co2m(i,k)*1e6_r8 !convet mol/mol to ppm
       END DO
    END DO
    DO i = 1, ncols
       auxco2(i,kmax+1) = co2m(i,kmax)*1e6_r8!convet mol/mol to ppm
       auxco2(i,kmax+2) = co2m(i,kmax)*1e6_r8!convet mol/mol to ppm
   END DO

    DO k=1,(kmax+2)
       DO j=1,(kmax+2)
          DO i = 1, ncols
             txuf(i,j,k) = h0p
          END DO
       END DO
    END DO

    DO k=1,kmax+3
       DO i = 1, ncols
          rawi(i,k) = (tui(i,k) - temp1) * dtbinv + h1p5
          it(i,k)   = rawi(i,k)
       END DO
    END DO

    DO k = 1, kmax+3
       DO i = 1, ncols
          it(i,k) = MAX(1,MIN(it(i,k), ntm1))
       END DO
    END DO

    DO k=1,(kmax+3)
       DO i = 1, (ncols)
          rawi(i,k) = it  (i,k)  -  h1p
          rawi(i,k) = tui (i,k)  - (temp1 + rawi(i,k) * dtb)
          rawi(i,k) = rawi(i,k)  *  dtbinv
       END DO
    END DO

    DO k = 1, kmax+3
       DO i = 1, ncols
          IF(it(i,k) .LE. ntm1 )THEN
             x1(i,k)=b2501(it(i,k))
             x2(i,k)=b2501(it(i,k)+1)
          ENDIF
       END DO
    END DO

    DO k=1,(kmax+3)
       DO i = 1, (ncols)
          tv1(i,k) = x1(i,k) + (x2(i,k) - x1(i,k)) * rawi(i,k)
       END DO
    END DO

    DO k = 1, kmax+3
       DO i = 1, ncols
          IF(it(i,k) .LE. ntm1 )THEN
             x2(i,k)=b2502(it(i,k)+1)
             x1(i,k)=b2502(it(i,k))
          ENDIF
       END DO
    END DO

    DO k=1,(kmax+3)
       DO i = 1, (ncols)
          tv2(i,k) = x1(i,k) + (x2(i,k) - x1(i,k)) * rawi(i,k)
       END DO
    END DO

    DO k=1,(kmax+2)
       DO i = 1, (ncols)
          shh(i,k) = tv1(i,k) + tv2(i,k)
       END DO
    END DO

    DO i = 1, ncols
       sumsav(i) = tv1(i,(kmax+3)) + tv2(i,(kmax+3))
    END DO

    DO k=1,(kmax+3)
       DO i = 1, (ncols)
          tui (i,k)  = tui(i,k) - h250p
          tui2(i,k)  = tui(i,k) * tui(i,k)
       END DO
    END DO
    !
    !     carbon dioxide and ozone fluxes are calculated here.
    !
    DO k = 1, kmax+3
       DO i = 1, ncols
          IF(it(i,k) .LE. ntm1 )THEN
             x1(i,k)=blkco2(it(i,k)+1)
             x2(i,k)=blkco2(it(i,k))
          ENDIF
       END DO
    END DO

    DO k=1,(kmax+3)
       DO i = 1, (ncols)
          cc(i,k) = x2(i,k) + (x1(i,k) - x2(i,k)) * rawi(i,k)
       ENDDO
    END DO

    DO k = 1, kmax+3
       DO i = 1, ncols
          IF(it(i,k) .LE. ntm1 )THEN
             x2(i,k)=blkwin(it(i,k))
             x1(i,k)=blkwin(it(i,k)+1)
          ENDIF
       END DO
    END DO

    DO k=1,(kmax+3)
       DO i = 1, (ncols)
          rawi(i,k) = x2(i,k) + (x1(i,k) - x2(i,k)) * rawi(i,k)
       END DO
    END DO
    !
    !     compute transmittances in the 15 micron and 9.6 micron bands.
    !
    DO k=1,(kmax+1)
       DO i = 1, (ncols)
          x1(i,k) = 2.5e-2_r8 * (tai(i,k+1) - 240.0e0_r8)
       END DO
    END DO

    ch(1:ncols,2:(kmax+2))=exp(x1(1:ncols,1:(kmax+1)))

    DO k=1,(kmax+1)
       DO i = 1, (ncols)
          x1(i,k) = 8.9e-3_r8 * (tai(i,k+1) - 240.0e0_r8)
       END DO
    END DO

    tai(1:ncols,2:(kmax+2))=exp(x1(1:ncols,1:(kmax+1)))

    DO i = 1, ncols
       tai(i,1) = h0p
       ch(i,1)  = h0p
    END DO

    DO k=1,(kmax+1)
       DO i = 1, (ncols)
          x1(i,k) = pai(i,k)
       END DO
    END DO

    DO k=1,kmax+1
       DO i = 1, ncols
          x1(i,k)=MAX(pp,x1(i,k))
       END DO
    END DO

    DO k=1,(kmax+1)
       DO i = 1, (ncols)
          fw(i,k) = x1(i,k) * pr(1)
       END DO
    END DO

    DO k=1,(kmax+1)
       DO i = 1, (ncols)
          x2(i,k) = EXP(0.85_r8* LOG(fw(i,k)))
       ENDDO
    END DO

    DO k=1,(kmax+1)
       DO i = 1, (ncols)
          tai(i,k+1) = x2(i,k) * tai(i,k+1)
          x1(i,k)  = x1(i,k) * pr(2)
       END DO
    END DO

    x1(1:ncols,1:kmax+1) = sqrt(x1(1:ncols,1:kmax+1))

    DO k=1,(kmax+1)
       DO i = 1, (ncols)
          ch(i,k+1) = x1(i,k) * ch(i,k+1)
          !
          !     hp25 = 0.26 for co2 = 330 ppmv
          !
          tai(i,k+1) = dp(i,k+1) * tai(i,k+1) * hp25 * auxco2(i,k) / 330.0_r8
          ch(i,k+1)  = dp(i,k+1) * ch (i,k+1) * hp25 * auxco2(i,k) / 330.0_r8
       END DO
    END DO

    DO ip = 2, (kmax+2)
       DO i = 1, ncols
          tai(i,ip) = tai(i,ip-1) + tai(i,ip)
          ch(i,ip)  = ch(i,ip-1)  + ch(i,ip)
       END DO
    END DO

    DO k=1,(kmax+1)
       DO i = 1, (ncols)
          css(i,k+1) = (tai(i,k+1)  + tai(i,k))  * hp5
          ccu(i,k+1) = ( ch(i,k+1)  +  ch(i,k))  * hp5
       END DO
    END DO

    IF ( ozone ) THEN

       DO i = 1, ncols
          fluxd(i,1) = h0p
       END DO

       DO k=1,(kmax+2)
          DO i = 1, (ncols)
             fluxu(i,k) = h0p
          END DO
       END DO

       DO k=1,(kmax+1)
          DO i = 1, (ncols)
             fluxd(i,k+1) = pai(i,k) * ozai(i,k+1)
          END DO
       END DO

       DO ix = 2, (kmax+2)
          DO i = 1, ncols
             fluxd(i,ix) = fluxd(i,ix-1) + fluxd(i,ix)
             fluxu(i,ix) = fluxu(i,ix-1) + ozai(i,ix)
          END DO
       END DO

       DO k = 1, (kmax+1)
          DO i = 1, (ncols)
             pai (i,k+1)  = (fluxd(i,k) + fluxd(i,k+1)) * hp5
             ozai(i,k+1)  = (fluxu(i,k) + fluxu(i,k+1)) * hp5
          ENDDO
       END DO

       DO k=1,kmax+1
          DO i=1, ncols
             wdel(i,k)=ABS(fluxu(i,k+1))
             fw(i,k)= ABS(fluxd(i,k+1))
          END DO
       END DO
    END IF
    indx1 = 1
    indx2 = 1
    CALL crunch(indx1 ,indx2 ,ncols  ,kmax  ,h0p   ,h1p   ,ozone ,txuf  , &
         tv1   ,tv2   ,tui   ,tui2  ,x1    ,x2    ,cc    ,rawi  , &
         x3    ,x4    ,ch    ,css   ,ccu   ,shi   ,shu   ,wdel  , &
         fw    ,pai   ,tai   ,ozai  ,ubar  ,vbar  ,wbar  ,ubarm , &
         vbarm ,wbarm ,fluxu ,fluxd )
    indx1 = 2
    indx2 = 2
    CALL crunch(indx1 ,indx2 ,ncols ,kmax  ,h0p   ,h1p   ,ozone ,txuf  , &
         tv1   ,tv2   ,tui   ,tui2  ,x1    ,x2    ,cc    ,rawi  , &
         x3    ,x4    ,ch    ,css   ,ccu   ,shi   ,shu   ,wdel  , &
         fw    ,pai   ,tai   ,ozai  ,ubar  ,vbar  ,wbar  ,ubarm , &
         vbarm ,wbarm ,fluxu ,fluxd )

    !
    !     downward flux txuf for water vapor
    !

    indx1 = 2
    indx2 = (kmax+2)
    CALL crunch(indx1 ,indx2 ,ncols ,kmax  ,h0p   ,h1p   ,ozone ,txuf  , &
         tv1   ,tv2   ,tui   ,tui2  ,x1    ,x2    ,cc    ,rawi  , &
         x3    ,x4    ,ch    ,css   ,ccu   ,shi   ,shu   ,wdel  , &
         fw    ,pai   ,tai   ,ozai  ,ubar  ,vbar  ,wbar  ,ubarm , &
         vbarm ,wbarm ,fluxu ,fluxd )
    !
    !     upward flux txuf for water vapor
    !
    indx1 = 1
    indx2 = (kmax+1)

    CALL crunch(indx1 ,indx2 ,ncols  ,kmax  ,h0p   ,h1p   ,ozone ,txuf  , &
         tv1   ,tv2   ,tui   ,tui2  ,x1    ,x2    ,cc    ,rawi  , &
         x3    ,x4    ,ch    ,css   ,ccu   ,shi   ,shu   ,wdel  , &
         fw    ,pai   ,tai   ,ozai  ,ubar  ,vbar  ,wbar  ,ubarm , &
         vbarm ,wbarm ,fluxu ,fluxd )

    DO k = 1,(kmax+2)
       DO i = 1, (ncols)
          shh(i,k) = shh(i,k) + rawi(i,k) + cc(i,k)
       END DO
    END DO

    DO i = 1, ncols
       sumsav(i) = sumsav(i) + rawi(i,(kmax+3)) + cc(i,(kmax+3))
    END DO
  END SUBROUTINE lwflux



  !cldslw:  Estimate a contribution from cloudiness
  !    Using cloud amount of large scale(css) and cumulus(ccu) clouds
  SUBROUTINE cldslw(ncols ,kmax  ,nlcs  ,h1p   ,cs    ,x1    , &
       x2    ,cc    ,x3    ,x4    ,ch    ,css   ,ccu)
    !
    !==========================================================================
    !  imax.......Number of grid points on a gaussian latitude circle
    !  kmax.......Number of grid points at vertical
    !  nlcs.......nlcs=30
    !  h1p........h1p    = 1.0e0      fac converts to degrees / time step
    !  cs.........probability of clear line-of-sight from each level to
    !             all other layers.
    !  x1.........Water vapor path (e-type ) and working dimension
    !  x2.........Water vapor path (band-center) and working dimension
    !  cc.........planck function at level temperature for co2 bands.
    !  x3.........water vapor path (band-wings), working dimension
    !  x4.........Working dimension
    !  ch.........probability of clear line-of-sight from level to top of
    !             the atmosphere.
    !  css........Large scale cloud amount
    !             css=css*(1-exp(-0.01*dp)) for ice cloud t < 253.0
    !  ccu........Cumulus cloud amount
    !=========================================================================!
    ! >>> icld=1     : old cloud emisivity setting                           !
    !       ccu = ccu*(1-exp(-0.05*dp))                                      !
    !       css = css*(1-exp(-0.01*dp))          for ice cloud t<253.0       !
    !       css = css*(1-exp(-0.05*dp))          for     cloud t>253.0       !
    ! >>> icld=2     : new cloud emisivity setting                           !
    !       ccu = 1.0-exp(-0.12*ccu*dp)                                      !
    !       css = 0.0                                    for  t<-82.5c        !
    !       css = 1-exp(-1.5e-6*(t-tcrit)**2*css*dp)     for -82.5<t<-10.0    !
    !       css = 1-exp(-5.2e-3*(t-273.)-0.06)*css*dp)   for -10.0<t< 0.0     !
    !       css = 1-exp(-0.06*css*dp)                    for t> 0.0c          !
    ! >>> icld = 3   : ccm3 based cloud emisivity                             !
    !=========================================================================!
    !==========================================================================

    INTEGER, INTENT(in   ) :: ncols
    INTEGER, INTENT(in   ) :: kmax
    INTEGER, INTENT(in   ) :: nlcs
    REAL(KIND=r8),    INTENT(in   ) :: h1p
    REAL(KIND=r8),    INTENT(inout  ) :: cs  (ncols,kmax+2,nlcs)
    REAL(KIND=r8),    INTENT(inout  ) :: x1  (ncols,kmax+3)
    REAL(KIND=r8),    INTENT(inout  ) :: x2  (ncols,kmax+3)
    REAL(KIND=r8),    INTENT(inout  ) :: cc  (ncols,kmax+3)
    REAL(KIND=r8),    INTENT(inout  ) :: x3  (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(inout  ) :: x4  (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(inout  ) :: ch  (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(in   ) :: css (ncols,kmax+2)
    REAL(KIND=r8),    INTENT(in   ) :: ccu (ncols,kmax+2)

    INTEGER :: i
    INTEGER :: k
    INTEGER :: j
    INTEGER :: ip
    INTEGER :: il
    INTEGER :: ix
    INTEGER :: ipm1
    INTEGER :: ipp1

    DO k=1,(kmax+2)
       DO i=1,(ncols)
          ch(i,k) = h1p
          cc(i,k) = h1p
          x1(i,k) = h1p
          x2(i,k) = h1p
          x3(i,k) = h1p
       END DO
    END DO

    DO k=1,(kmax+2)
       DO j=1,(kmax+2)
          DO i=1,(ncols)
             cs(i,j,k) = h1p
          END DO
       END DO
    END DO

    DO ip=2,(kmax+2)
       DO i=1,ncols
          x1(i,ip) = h1p - ccu(i,ip-1)
       END DO

       DO i=1,ncols
          x1(i,ip)=MIN(x1(i,ip),x1(i,ip-1))
       END DO
       DO i=1,ncols
          ch(i,ip) = ch(i,ip-1) * (h1p - css(i,ip-1))
       END DO
    END DO

    DO k=1,(kmax+1)
       DO i=1,(ncols)
          ch(i,k+1) = ch(i,k+1) * x1(i,k+1)
       END DO
    END DO

    DO il=1,(kmax+1)
       ip=(kmax+2)-il
       DO i=1,ncols
          x2(i,ip) = h1p - ccu(i,ip)
       END DO

       DO i=1,ncols
          x2(i,ip)=MIN(x2(i,ip),x1(i,ip+1))
       END DO

       DO i=1,ncols
          cc(i,ip) = cc(i,ip+1) * (h1p - css(i,ip))
       END DO
    END DO

    DO k=1,(kmax+1)
       DO i=1,(ncols)
          cc(i,k) = cc(i,k) * x2(i,k)
       END DO
    END DO

    DO ip=2,(kmax+2)
       ipm1   = ip - 1

       DO ix=1,ipm1
          DO i=1,ncols
             x4(i,ix)    = h1p - ccu(i,ip-1)
             cs(i,ix,ip) = cs(i,ix,ip-1) * (h1p - css(i,ip-1))
          END DO
       END DO

       IF (ip .GT. 2) THEN
          DO j=1,ipm1
             DO i=1,ncols
                cs(i,j,ip-1) = cs(i,j,ip-1) * x3(i,j)
             END DO
          END DO
       END IF

       DO k=1,ipm1
          DO i=1,ncols
             x4(i,k)=MIN(x4(i,k),x3(i,k))
             x3(i,k) = x4(i,k)
          END DO
       END DO

    END DO

    DO j=1,ipm1
       DO i=1,ncols
          cs(i,j,(kmax+2)) = cs(i,j,(kmax+2)) * x4(i,j)
       END DO
    END DO

    DO k=1,kmax+2
       DO i=1,ncols
          x3(i,k) = h1p
       END DO
    END DO

    DO il=1,(kmax+1)
       ip   = (kmax+2) - il
       ipp1 = ip + 1

       DO ix=ipp1,(kmax+2)
          DO i=1,ncols
             x4(i,ix)    = h1p - ccu(i,ip)
             cs(i,ix,ip) = cs(i,ix,ip+1) * (h1p - css(i,ip))
          END DO
       END DO

       IF (il .GT. 1) THEN
          DO k=ip,((kmax+1)-ip)! DO k=ip,((kmax+2)-ip)
             DO  i=1, ncols
                cs(i,k+2,ip+1) = cs(i,k+2,ip+1) * x3(i,k+2)
             END DO
          END DO
       END IF

       DO k=ipp1,((kmax+2)-ip)
          DO i=1,ncols
             x4(i,j) = MIN(x4(i,j),x3(i,j))
             x3(i,k) = x4(i,k)
          END DO
       END DO

    END DO

    DO k=ipp1,((kmax+2)-ip)
       DO i=1, ncols
          cs(i,k,ip) = cs(i,k,ip) * x4(i,k)
       END DO
    END DO

  END SUBROUTINE cldslw


  ! lwrad  :compute upward and downward fluxes.
  SUBROUTINE lwrad( &
       ! Model Info and flags
       ncols ,kmax  ,nls   ,nlcs  , noz   ,icld  ,&
       ! Atmospheric fields
       pl20  ,pl    ,tl    ,ql    , o3l   ,tg    ,&
       co2m,                                    &
       ! LW Radiation fields 
       ulwclr,ulwtop,atlclr,atl   ,rsclr , rs    ,&
       dlwclr,dlwbot,                             &
       ! Cloud field and Microphysics
       cld   ,clu   ,clwp  ,fice  ,rei   ,emisd     )
    IMPLICIT NONE
    ! input variables
    !
    !     noz,  pl,  pl20, tl,  tg,  ql,  o3l,  cld,  clu
    !
    !     output variables
    !
    !     atl, atlclr,  rs,  rsclr, ulwtop, ulwclr, dlwbot, dlwclr
    !
    !     parameter list variables
    !
    !     nim......number of grid points around a latitude circle.
    !     nlm......number of model layers.
    !     nlmp1....nlm plus one.
    !     nlmp2....nlm plus two.
    !     nlmp3....nlm plus three.
    !     nls......number of layers in the stratosphere.
    !
    ! local variables
    !
    !     b250.....planck function table for water vapor bands; center and
    !              wing.
    !     ntm1.....number of rows in b250 - 1.
    !     temp1....lowest temperature for which b250 is tabulated.
    !     dtb......temperature increment in b250.
    !     nup1.....number of columns in gl and coeff.
    !     ubar.....scaled water vapor path length in window.
    !     vbar.....scaled water vapor path length in center.
    !     wbar.....scaled water vapor path length in wing.
    !     rawi.....water vapor amount in layer.
    !     ozai.....ozone amount in layer.
    !     shh......planck function at level temperature for water vapor
    !              bands.
    !     shi......g-function for a path from level to top of atmosphere.
    !     txuf.....1.used as matrix of g-functions for paths from each level
    !                to all other layers.
    !              2.used for transmission in co2 band.
    !              3.used for transmission in ozone band.
    !              4.in cldslw used for probability of clear line-of-sight
    !                from each level to all other layers for max overlap.
    !     wv,ww....interpolated value of gl for band center and wing
    !              respectively.
    !     yv,yw....linear term in temperature correction for band center
    !              and wing respectively.
    !     zv,zw....quadratic term in temperature correction for band center
    !              and wing respectively.
    !     blkco2...planck function table for co2 bands.
    !     blkwin...planck function table for ozone band.
    !     cc.......planck function at level temperature for co2 bands.
    !     vbarm....planck function at level temperature for ozone band.
    !     pp.......doppler broadening cut-off.
    !     pscalv...scaled co2 amount in band center.
    !     pscalw...scaled co2 amount in band wing.
    !     dx.......parameterized optical depth of water vapor line in co2
    !              band wing.
    !     dy.......parameterized optical depth of water vapor continuum in
    !              co2 band wing, also used for ozone band wing.
    !     sv,sw....minimum in table of log water vapor amount in band center
    !              and wing respectively.
    !     cs.......probability of clear line-of-sight from each level to
    !              all other layers.
    !     ch.......probability of clear line-of-sight from level to top of
    !              the atmosphere.
    !     cc.......probability of clear line-of-sight from level to surface.
    !     ct.......probability of clear line-of-sight from level to top of
    !              atmosphere for maximun overlap.
    !     cu.......probability of clear line-of-sight from level to surface
    !              for maximum overlap.
    !
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! >>> icld=1     : old cloud emisivity setting                      c
    !       ccu = ccu*(1-exp(-0.05*dp))                                 c
    !       css = css*(1-exp(-0.01*dp))          for ice cloud t<253.0  c
    !       css = css*(1-exp(-0.05*dp))          for     cloud t>253.0  c
    ! >>> icld=2     : new cloud emisivity setting                      c
    !       ccu = 1.0-exp(-0.12*ccu*dp)                                 c
    !       css = 0.0                                for  t<-82.5c      c
    !       css = 1-exp(-1.5e-6*(t-tcrit)**2*css*dp) for -82.5<t<-10.0c c
    !       css = 1-exp(-5.2e-3*(t-273.)-0.06)*css*dp)for -10.0<t< 0.0c c
    !       css = 1-exp(-0.06*css*dp)                 for t> 0.0c       c
    ! >>> icld = 3   : ccm3 based cloud emisivity                       c
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    ! ncols......Number of grid points on a gaussian latitude circle
    ! kmax......Number of grid points at vertical
    ! nls.......number of layers in the stratosphere.
    ! nlcs......nlcs=30
    ! noz.......constant logical noz = .FALSE.
    ! tg........Surface Temperature (K)
    ! pl20......pl20(i,k)=gps(i)*sigml(kflip) where
    !                         gps   =  surface pressure   (mb)
    !                         sigml =  sigma coordinate at bottom of layer
    ! pl........Pressure at Middle of Layer(mb)
    ! tl........Temperature at middle of Layer (K)
    ! ql........Specific Humidity at middle of layer (g/g)
    ! o3l.......Ozone Mixing ratio at middle of layer (g/g)
    ! cld.......Large scale cloud amount in layers
    ! clu.......Cumulus cloud amount in layers
    ! ulwclr....Upward flux at top in clear case (W/m2)
    ! ulwtop....Upward flux at top (W/m2)
    ! atlclr....Heating rate in clear case (K/s)
    ! atl.......Heating rate (K/s)
    ! rsclr.....Net surface flux in clear case (W/m2 )
    ! rs........Net surface flux
    ! dlwclr....Downward flux at surface in clear case (W/m2 )
    ! dlwbot....Downward flux at surface (W/m2 )
    ! clwp
    ! fice......controle of change of fase of water
    ! rei.......determine rei as function of normalized pressure
    ! emisd.....emis(i,kflip) = 1.- EXP(-1.66*rkabs(i,k)*clwp(i,k))
    ! co2m....co2val is wgne standard value in ppm "co2val = /345.0/
    !

    ! Model Info and flags
    INTEGER, INTENT(in   ) :: ncols
    INTEGER, INTENT(in   ) :: kmax
    INTEGER, INTENT(in   ) :: nls
    INTEGER, INTENT(in   ) :: nlcs
    LOGICAL, INTENT(in   ) :: noz
    INTEGER, INTENT(in   ) :: icld

    ! Atmospheric fields
    REAL(KIND=r8),    INTENT(in   ) :: pl20  (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: pl    (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: tl    (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: ql    (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: o3l   (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: tg    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: co2m(ncols,kMax)  !mol/mol 
    ! LW Radiation fields 
    REAL(KIND=r8),    INTENT(inout  ) :: ulwclr(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: ulwtop(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: atlclr(ncols,kmax)
    REAL(KIND=r8),    INTENT(inout  ) :: atl   (ncols,kmax)
    REAL(KIND=r8),    INTENT(inout  ) :: rsclr (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: rs    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: dlwclr(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: dlwbot(ncols)

    ! Cloud field and Microphysics
    REAL(KIND=r8),    INTENT(in   ) :: cld   (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: clu   (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: clwp  (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: fice  (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: rei   (ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: emisd (ncols,kmax)

    REAL(KIND=r8)    :: txuf  (ncols,kmax+2,kmax+2)
    REAL(KIND=r8)    :: cs    (ncols,kmax+2,nlcs)
    REAL(KIND=r8)    :: tv1   (ncols,kmax+3)
    REAL(KIND=r8)    :: tv2   (ncols,kmax+3)
    REAL(KIND=r8)    :: tui   (ncols,kmax+3)
    REAL(KIND=r8)    :: x1    (ncols,kmax+3)
    REAL(KIND=r8)    :: x2    (ncols,kmax+3)
    REAL(KIND=r8)    :: cc    (ncols,kmax+3)
    REAL(KIND=r8)    :: rawi  (ncols,kmax+3)
    REAL(KIND=r8)    :: x3    (ncols,kmax+2)
    REAL(KIND=r8)    :: x4    (ncols,kmax+2)
    REAL(KIND=r8)    :: ch    (ncols,kmax+2)
    REAL(KIND=r8)    :: dp    (ncols,kmax+2)
    REAL(KIND=r8)    :: css   (ncols,kmax+2)
    REAL(KIND=r8)    :: ccu   (ncols,kmax+2)
    REAL(KIND=r8)    :: shi   (ncols,kmax+2)
    REAL(KIND=r8)    :: shh   (ncols,kmax+2)
    REAL(KIND=r8)    :: shu   (ncols,kmax+1)
    REAL(KIND=r8)    :: suma  (ncols)
    REAL(KIND=r8)    :: sumsav(ncols)
    LOGICAL :: bitx  (ncols*(kmax+3))


    REAL(KIND=r8) :: pai   (ncols,kmax+2)
    REAL(KIND=r8) :: tai   (ncols,kmax+2)
    REAL(KIND=r8) :: ozai  (ncols,kmax+2)
    REAL(KIND=r8) :: ubar  (ncols,kmax+2)
    REAL(KIND=r8) :: vbar  (ncols,kmax+2)
    REAL(KIND=r8) :: wbar  (ncols,kmax+2)
    REAL(KIND=r8) :: ubarm (ncols,kmax+2)
    REAL(KIND=r8) :: vbarm (ncols,kmax+2)
    REAL(KIND=r8) :: wbarm (ncols,kmax+2)
    REAL(KIND=r8) :: fluxu (ncols,kmax+2)
    REAL(KIND=r8) :: fluxd (ncols,kmax+2)


    REAL(KIND=r8)      :: h0p
    REAL(KIND=r8)      :: h1p
    REAL(KIND=r8)      :: h1p5
    REAL(KIND=r8)      :: hp5
    REAL(KIND=r8)      :: dtb
    REAL(KIND=r8)      :: dtbinv
    REAL(KIND=r8)      :: pr(2)
    INTEGER   :: ntm1
    INTEGER   :: imnpnp
    LOGICAL   :: ozone

    REAL(KIND=r8)      :: emis (ncols,kmax)
    REAL(KIND=r8)      :: emis1(ncols,kmax+2)
    REAL(KIND=r8)      :: rkabs(ncols,kmax)

    REAL(KIND=r8)      :: h1p02
    REAL(KIND=r8)      :: h6p08
    REAL(KIND=r8)      :: tice

    INTEGER  :: ls1
    INTEGER  :: ls2
    INTEGER  :: imls1
    INTEGER  :: imlsm1
    INTEGER  :: imlm1
    INTEGER  :: imt2
    INTEGER  :: npmls1
    INTEGER  :: npmls2
    INTEGER  :: i
    INTEGER  :: ip
    INTEGER  :: ipm1
    INTEGER  :: ix
    INTEGER  :: ipp1
    INTEGER  :: l
    INTEGER  :: k
    INTEGER  :: kflip
    REAL(KIND=r8)     :: fac
    REAL(KIND=r8)     :: h3ppm
    REAL(KIND=r8)     :: tcrit
    REAL(KIND=r8)     :: ecrit
    REAL(KIND=r8)     :: d642
    REAL(KIND=r8)     :: pre   (2)

    h1p02 =   1.02e0_r8
    h6p08 =   6.0811e0_r8
    tice  = 273.16_r8
    !
    !     copy parameters into local variables
    !
    ls1    = nls+1
    ls2    = nls+2
    imls1  = ncols*ls1
    imlsm1 = ncols*(nls-1)
    imnpnp = (ncols*(kmax+2))*(kmax+2)
    imlm1  = ncols*(kmax-1)
    imt2   = ncols*2
    npmls1 = ncols* ((kmax+2) - ls1)
    npmls2 = ncols* ((kmax+2) - ls2)
    !
    !     fac converts to degrees / time step
    !
    fac    = 9.8e-1_r8 / 1.0030e04_r8
    h1p    = 1.0e0_r8
    h0p    = 0.0e0_r8
    h1p5   = 1.5e0_r8
    hp5    = 0.5e0_r8
    h3ppm  = 3.0e-6_r8
    tcrit  = tice - 82.5_r8
    ecrit  = 0.007884375_r8
    ntm1   = 31
    dtb    = 5.0e0_r8
    dtbinv = h1p / dtb
    pre(1) = h1p / 2.75e02_r8
    pre(2) = h1p / 5.5e02_r8
    d642   = h1p / 6.426e02_r8
    pr(1)  = h1p / 3.0e01_r8
    pr(2)  = h1p / 3.0e02_r8
    emis   = 0.0e0_r8
    emis1  = 0.0e0_r8
    emisd  = 0.0e0_r8
    !
    !     set ozone logical variable
    !
    ozone = (.NOT. noz)
    !
    !     ptop and dp at top don't change
    !
    DO i = 1, ncols
       dp (i,1)     = h0p
       x1 (i,1)     = h1p
       dp (i,2)     = h1p
       pai(i,1)     = hp5
    END DO

    DO k=1,((kmax+2) - (nls+1))
       DO i = 1, ncols
          rawi(i,nls+1+k) = ql(i,nls+k-1)
       END DO
    END DO

    IF (nls > 1) THEN
       DO k = 1,(nls-1)
          DO i = 1, ncols
             rawi(i,2+k) = h3ppm
          END DO
       END DO
    END IF

    DO k=1,kmax
       DO i = 1, ncols
          pai(i,k+1) = pl(i,k)
          tai(i,k+2) = tl(i,k)
          x1 (i,k+1) = pl20(i,k)
       END DO
    END DO

    DO i = 1, ncols
       pai(i,(kmax+2)) = pl20(i,kmax)
    END DO

    IF ( ozone ) THEN
       DO k = 1,kmax
          DO i = 1, (ncols)
             ozai(i,k+2) = o3l(i,k)
          END DO
       END DO
    END IF

    DO k = 1,kmax
       DO i = 1, ncols
          dp(i,k+2) = x1(i,k+1) - x1(i,k)
       END DO
    END DO
    !
    !     temperature and humidity interpolations
    !
    DO k = 1,kmax
       DO i = 1, ncols
          rawi(i,k+2)=MAX(0.1e-22_r8,rawi(i,k+2))
       END DO
    END DO

    IF ( ozone ) THEN

       DO k = 1,kmax
          DO i = 1, ncols
             ozai(i,k+2)=MAX(0.1e-9_r8,ozai(i,k+2))
          END DO
       END DO

       DO i = 1, ncols
          ozai(i,1) = h0p
          ozai(i,2) = ozai(i,3) * dp(i,2) * h1p02
       END DO

       DO k = 1,kmax
          DO  i = 1, ncols
             ozai(i,k+2) = ozai(i,k+2) * dp(i,k+2)  * h1p02
          END DO
       END DO

    END IF
    !
    !     do temperature interpolation
    !
    DO i = 1, ncols
       rawi(i,1) = h0p
       rawi(i,2) = h3ppm * dp(i,2) * h1p02
    END DO

    DO k=1,kmax
       DO i = 1, ncols
          rawi(i,k+2) = rawi(i,k+2) * dp(i,k+2) * h1p02
          x1(i,k)   = x1(i,k+1) / pai(i,k+1)
       END DO
    END DO

    x2(1:ncols,1:kmax)=LOG(x1(1:ncols,1:kmax))

    DO k = 1,kmax
       DO i = 1, ncols
          x1(i,k)   = pai(i,k+2) / pai(i,k+1)
       END DO
    END DO

    x4(1:ncols,1:kmax)=LOG(x1(1:ncols,1:kmax))

    DO k=1,(kmax-1)
       DO i = 1, ncols
          tui(i,k+2) = tai(i,k+2)  + x2(i,k) / x4(i,k) &
               * (tai(i,k+3) - tai(i,k+2))
       END DO
    END DO
    !
    !     set surface air temperature as mean of lowest layer t and tg
    !
    DO i = 1, ncols
       tui(i,(kmax+2)) = hp5 * ( tai(i,(kmax+2)) + tg(i) )
       tui(i,(kmax+3)) = tg(i)
       tai(i,1)        = tai(i,3)
       tai(i,2)        = tai(i,3)
    END DO

    DO k=1,2
       DO i = 1, ncols
          tui(i,k)     = tai(i,k)
       END DO
    END DO

    !
    !     compute scaled water vapor amounts
    !

    DO i = 1, ncols
       ubar(i,1) = h0p
       vbar(i,1) = h0p
       wbar(i,1) = h0p
    END DO

    DO k=1,(kmax+1)
       DO i = 1, (ncols)
          x3(i,k) = 0.5e-2_r8 * (tai(i,k+1) - 225.0e0_r8)
       END DO
    END DO

    x2(1:ncols,1:kmax+1)=EXP(x3(1:ncols,1:kmax+1))

    DO k=1,(kmax+1)
       DO i = 1, (ncols)
          vbar(i,k+1) = rawi(i,k+1) * x2(i,k) &
               * (pai(i,k) * pre(1))
          x3(i,k)     = 1.6e-2_r8 * (tai(i,k+1) - 256.0e0_r8)
       END DO
    END DO

    x2(1:ncols,1:kmax+1)=EXP(x3(1:ncols,1:kmax+1))

    DO k = 1,(kmax+1)
       DO i = 1, (ncols)
          wbar(i,k+1) = rawi(i,k+1) * x2(i,k) &
               * (pai(i,k) * pre(2))
          x1(i,k)     = 1.8e03_r8    / tai(i,k+1) - h6p08
       END DO
    END DO

    x2(1:ncols,1:kmax+1)=EXP(x1(1:ncols,1:kmax+1))

    DO k = 1, (kmax+1)
       DO i = 1, (ncols)
          x2(i,k) = x2(i,k) * pai(i,k)  * d642
          x2(i,k) = x2(i,k) * rawi(i,k+1) &
               * rawi(i,k+1) / dp(i,k+1)
       END DO
    END DO

    DO ip = 2, (kmax+2)
       DO i = 1, ncols
          ubar(i,ip) = ubar(i,ip-1) + x2(i,ip-1)
          vbar(i,ip) = vbar(i,ip-1) + vbar(i,ip)
          wbar(i,ip) = wbar(i,ip-1) + wbar(i,ip)
       END DO
    END DO

    DO k=1,(kmax+1)
       DO i = 1, (ncols)
          ubarm(i,k+1) = (ubar(i,k+1) + ubar(i,k)) * hp5
          vbarm(i,k+1) = (vbar(i,k+1) + vbar(i,k)) * hp5
          wbarm(i,k+1) = (wbar(i,k+1) + wbar(i,k)) * hp5
       END DO
    END DO

    CALL lwflux(pai   ,tai   ,ozai  ,ubar  ,vbar  ,wbar  ,ubarm ,vbarm , &
         wbarm ,fluxu ,fluxd ,txuf  ,tv1   ,tv2   ,tui   ,x1    , &
         x2    ,cc    ,rawi  ,x3    ,x4    ,ch    ,dp    ,css   , &
         ccu   ,shi   ,shh   ,shu   ,sumsav,h0p   ,h1p   ,h1p5  , &
         hp5   ,dtb   ,dtbinv,pr    ,ntm1  ,ozone ,co2m, &
         ncols ,kmax  )
    !
    !     compute clear sky fluxes only for cloud radiative forcing
    !     zero out downward flux accumulators at l=1
    !
    DO i = 1, ncols
       fluxd(i,1) = h0p
       !
       !     upward flux at the surface
       !
       fluxu(i,(kmax+2)) = sumsav(i)
    END DO
    !
    !     compute downward fluxes
    !
    DO ip = 2, (kmax+2)
       DO i = 1, ncols
          suma(i) = h0p
       END DO
       ipm1 = ip - 1
       DO ix = 1, ipm1
          DO i = 1, ncols
             suma(i) = suma(i) + txuf(i,ix,ip)
          END DO
       END DO
       DO i = 1, ncols
          fluxd(i,ip) = suma(i) - shi(i,ip) + shh(i,ip)
       END DO
    END DO
    !
    !     compute upward fluxes
    !
    DO ip = 1, (kmax+1)
       ipp1      = ip + 1
       DO i = 1, ncols
          suma(i) = h0p
       END DO
       DO ix = ipp1, (kmax+2)
          DO i = 1, ncols
             suma(i) = suma(i) + txuf(i,ix,ip)
          END DO
       END DO
       DO i = 1, ncols
          fluxu(i,ip) = suma(i)  + shu(i,ip) + shh(i,ip)
       END DO
    END DO

    DO k=1,(kmax+2)
       DO  i = 1, (ncols)
          x1(i,k) = fluxd(i,k) - fluxu(i,k)
       END DO
    END DO

    DO i = 1, ncols
       ulwclr(i)  = -x1(i,2)
       rsclr (i)  = -x1(i,(kmax+2))
       dlwclr(i)  = fluxd(i,(kmax+2))
    END DO

    DO k = 1,kmax
       DO i = 1, (ncols)
          atlclr(i,k) = (x1(i,k+1) - x1(i,k+2)) * fac / dp(i,k+2)
       END DO
    END DO

    !
    !     we don't allow clouds between gcm levels in the stratosphere
    !

    DO k=1,kmax
       DO i = 1, (ncols)
          css(i,k+1)   = cld(i,k)
          ccu(i,k+1)   = clu(i,k)
       END DO
    END DO

    DO k=1,(nls+1)
       DO i = 1, ncols
          css(i,k) = h0p
          ccu(i,k) = h0p
       END DO
    END DO
    !
    !     icld = 1     : old cloud emissivity setting
    !     icld = 2     : new cloud emissivity setting
    !
    IF (icld == 1) THEN

       DO k=1,((kmax+2) - (nls+2))
          DO i = 1, ncols
             tv1(i,nls+k+1) = -0.05e0_r8 *  dp(i,nls+k+2)
          END DO
       END DO

       DO k=1,((kmax+2) - (nls+2))
          DO i = 1,ncols
             tv2(i,nls+k+1)=EXP(tv1(i,nls+k+1))
             ccu(i,nls+k+1) = ccu(i,nls+k+1) * ( h1p - tv2(i,nls+k+1) )
          END DO
       END DO

       DO l = ls2, (kmax+1)
          DO i = 1, ncols
             !
             !     bash down cloud emissivities
             !
             IF (tl(i,l-1) < 253.0_r8) THEN
                tv1(i,l) = -0.01e0_r8 * dp(i,l+1)
             END IF
          END DO
       END DO

       DO k=1,((kmax+2) - (nls+2))
          DO i = 1,  ncols
             tv2(i,nls+k+1)=EXP(tv1(i,nls+k+1))
             css(i,nls+k+1) = css(i,nls+k+1) * ( h1p - tv2(i,nls+k+1) )
          END DO
       END DO

    ELSE IF (icld == 2) THEN

       DO k=1,((kmax+2) - (nls+2))
          DO i = 1,  ncols
             tv1(i,nls+k+1) = -0.12e0_r8 * ccu(i,nls+k+1) * dp(i,nls+k+2)
          END DO
       END DO

       DO k=1,((kmax+2) - (nls+2))
          DO i = 1,  ncols
             tv2(i,nls+k+1)=EXP(tv1(i,nls+k+1))
             ccu(i,nls+k+1) =  h1p - tv2(i,nls+k+1)
          END DO
       END DO

       DO l = ls2, (kmax+1)
          DO i = 1, ncols
             tv2(i,l) =  tl(i,l-1) - tcrit
          END DO
       END DO

       DO k=1,((kmax+2) - (nls+2))
          DO i = 1,  ncols
             tv2(i,nls+k+1) = MAX(tv2(i,nls+k+1),h0p)
             tv1(i,nls+k+1) = MIN(1.5e-6_r8 * tv2(i,nls+k+1) * tv2(i,nls+k+1),ecrit)
          END DO
       END DO

       DO l = ls2, (kmax+1)
          DO i = 1, ncols
             tv2(i,l) =  tl(i,l-1) - tice
          END DO
       END DO

       DO k=1,((kmax+2) - (nls+2))
          DO i = 1,  ncols
             tv2(i,nls+k+1) = MIN(MAX(5.2115625e-3_r8 * tv2(i,nls+k+1) + 0.06e0_r8 ,ecrit),0.06_r8)
             bitx(k) = (tv1(i,k).eq.ecrit)
             IF (bitx(k)) tv1(i,k)=tv2(i,k)
          END DO
       END DO

       DO l = ls2, (kmax+1)
          DO i = 1, ncols
             tv1(i,l) =  -tv1(i,l) * css(i,l) * dp(i,l+1)
          END DO
       END DO

       DO k=1,((kmax+2) - (nls+2))
          DO i = 1,  ncols
             tv2(i,nls+k+1)=EXP(tv1(i,nls+k+1))
             css(i,nls+k+1) =  h1p - tv2(i,nls+k+1)
          END DO
       END DO

    ELSE IF ((icld == 3).OR.(icld == 4).OR.(icld == 5).OR.(icld == 6) .OR.(icld == 7)) THEN

       DO k = 1, kmax
          kflip=kmax+1-k
          DO i = 1, ncols
             rkabs(i,k) = 0.090361_r8*(1.0_r8-fice(i,k)) + &
                  (0.005_r8 + 1.0_r8/rei(i,k))*fice(i,k)
             emis(i,k) =MIN(MAX(1.0_r8- EXP(-1.66_r8*rkabs(i,k)*clwp(i,k)),0.0_r8),1.0_r8)
             emisd(i,k)=MIN(MAX(emis(i,kflip),0.0_r8),1.0_r8)
          END DO
       END DO

       DO k=1,kmax
          DO i = 1, (ncols)
             emis1(i,k+1) = emis(i,k)
          END DO
       END DO

       DO k=1,nls+1
          DO i = 1,  ncols
             emis1(i,k) = h0p
          END DO
       END DO

       DO k=nls+2,kmax+1
          DO i = 1, ncols
             css(i,k) = css(i,k) * emis1(i,k)
             ccu(i,k) = ccu(i,k) * emis1(i,k)
          END DO
       END DO

    END IF
    !
    !     get the contribution from cloudiness
    !
    CALL cldslw(ncols ,kmax  ,nlcs  ,h1p   ,cs    ,x1    ,x2    , &
         cc    ,x3    ,x4    ,ch    ,css   ,ccu   )

    DO k=1,(kmax+2)
       DO l=1,(kmax+2)
          DO i = 1, (ncols)
             txuf(i,l,k) = txuf(i,l,k) * cs(i,l,k)
          END DO
       END DO
    END DO

    DO k=1,(kmax+2)
       DO i = 1, (ncols)
          shi(i,k) = shi(i,k) * ch(i,k)
       END DO
    END DO

    DO k=1,(kmax+1)
       DO i = 1, (ncols)
          shu(i,k) = shu(i,k) * cc(i,k)
          shu(i,k) = shu(i,k) + shh(i,k)
       END DO
    END DO
    !
    !     zero out downward flux accumulators at l=1
    !
    DO i = 1, ncols
       fluxd(i,1) = h0p
    END DO
    !
    !     upward flux at the surface
    !
    DO i = 1, ncols
       fluxu(i,(kmax+2)) = sumsav(i)
    END DO
    !
    !     upward flux at the surface was computed before call to cldslw
    !     compute downward fluxes
    !
    DO ip = 2, (kmax+2)

       DO i = 1, ncols
          suma(i) = h0p
       END DO

       ipm1 = ip - 1

       DO ix = 1, ipm1
          DO i = 1, ncols
             suma(i) = suma(i) + txuf(i,ix,ip)
          END DO
       END DO

       DO i = 1, ncols
          fluxd(i,ip) = suma(i) - shi(i,ip) + shh(i,ip)
       END DO

    END DO
    !
    !     compute upward fluxes
    !
    DO ip = 1, (kmax+1)
       ipp1 = ip + 1

       DO i = 1, ncols
          suma(i) = h0p
       END DO

       DO ix = ipp1, (kmax+2)
          DO i = 1, ncols
             suma(i) = suma(i) + txuf(i,ix,ip)
          END DO
       END DO

       DO i = 1, ncols
          fluxu(i,ip) = suma(i) + shu(i,ip)
       END DO

    END DO

    DO k=1,(kmax+2)
       DO i = 1, (ncols)
          x1(i,k) = fluxd(i,k) - fluxu(i,k)
       END DO
    END DO

    DO i = 1, ncols
       ulwtop(i) = -x1(i,2)
       rs(i)     = -x1(i,(kmax+2))
       dlwbot(i) = fluxd(i,(kmax+2))
    END DO

    DO k=1,kmax
       DO i = 1, (ncols)
          atl(i,k) = (x1(i,k+1) - x1(i,k+2)) * fac / dp(i,k+2)
       END DO
    END DO

  END SUBROUTINE lwrad


  !------------------------------------------------------------------------
  ! SHORT WAVE FLUXES CALCULATION
  !    Original Paper:
  !       Lacis and Hansen, 1974: "A Parameterization for the Absorption of 
  !       Solar Radiation in the Earth's Atmosphere", J. Atmos. Sci., v31, 
  !       118-133
  !------------------------------------------------------------------------


  !----------------------------------------------------------------------C
  ! Subroutine: CLOUDY
  !
  ! CALCULATES:
  !    computes  ozone heating, downflux and ground absorption, water
  !    vapor heating and computes for cycles over five bands total
  !    optical depth, reflection and transmission for diffuse
  !    incidence,direct transmission, diffuse reflection and
  !    transmission for direct beam, upward and dwonward fluxes at
  !    layer boundaries, absorption in the column.
  !
  !----------------------------------------------------------------------C

  SUBROUTINE cloudy (ncols ,kmax  ,scosc ,cmuc  ,csmcld,dscld ,rvbc  ,rvdc  , &
       rnbc  ,rndc  ,agvcd ,agncd ,rsvcd ,rsncd ,sc    ,rco   , &
       rcg   ,taut  ,rc2   ,tr1   ,tr2   ,tr3   ,ta    ,wa    , &
       oa    ,e0    ,pu    ,ozcd  ,swale ,swil  ,css   ,acld  , &
       dpc   ,swilc ,ccu   ,tauc  ,litd  ,sqrt3 ,gg    ,ggp   , &
       ggsq  ,np    ,lmp1  ,nsol  ,ncld  ,nclmp1,ncldnp,dooz    )
    !
    !==========================================================================
    !     imax......Number of grid points on a gaussian latitude circle
    !     kmax......Number of grid points at vertical
    !     scosc.....scosz(i)   = s0     * cmu(i)
    !                         where s0  is constant solar
    !                               cmu is cosine of solar zenith angle
    !     cmuc......cmuc is cosine of solar zenith angle
    !     csmcld....csmcld = cosmag(i)  = 1224.0 * cmu(i) * cmu(i) + 1.0
    !     dscld.....Total absorption in atmosphere and ground
    !     rvbc......Visible beam cloudy flux
    !     rvdc......Visible diffuse cloudy flux
    !     rnbc......Nir beam cloudy flux
    !     rndc......Nir diffuse cloudy flux
    !     agvcd.....Ground Visible Diffuse Albedo
    !     agncd.....Ground near Infrared Diffuse albedo
    !     rsvcd.....Ground Viseible beam albedo
    !     rsncd.....Ground near infrared beam albedo
    !     sc........Ground absorption
    !     rco.......Cloudy an ground reflection for ozone absorption comp.
    !     rcg.......Cloudy and ground reflection for ground absoption
    !               computation in visible region of spectrum
    !     taut......Downward visible beam
    !     rc2.......cloudy and ground reflection for ground absoption
    !               computation in region of spectrum from 0.7 to 0.9 mcm
    !     tr1.......Extinction of beam radiation (clear)
    !     tr2.......Extinction of beam radiation (cloud)
    !     tr3.......Extinction of beam radiation (Cloud, Raley)
    !     ta........Layer Temperature in DLGP
    !     wa........Layer specific humidity in DLGP
    !     oa........Layer ozone mixing ratio in DLGP
    !     e0........cloud optical depth
    !     pu........pressure at botton of layer in DLGP
    !     ozcd......ozone amount in column in CDLGP
    !     swale.....water vapor amount in column in DLGP
    !     swil......water vapor amount in layer in DLGP
    !     css.......cloud amount
    !     acld......Heating rate (cloudy)
    !     dpc.......pressure defference
    !     swilc.....water vapor amount in layer in CDLGP
    !     ccu.......cumulus cloud amount in DLGP
    !     tauc......Cloud optical depth in cloudy DLGP
    !     litd......numbers de CDLGP in all layers
    !     sqrt3.....Magification factor for diffuse reflected radiation sqrt3  = SQRT(3.0)
    !     gg........Asymmetry Factor      = 0.85
    !     ggp.......ggp    = gg  / (1.0 + gg) * 0.75
    !     ggsq......ggsq   = gg  * gg
    !     np........np = (kmax+2)
    !     lmp1......lmp1   = (kmax+1)
    !     nsol......Number of solar latitude grid points (cosz>0.01)
    !     ncld......Number of cloudy DLGP
    !     nclmp1....NCLD*(KMAX+1), where NCLD is number of cloudy DLGP
    !     ncldnp....NCLD*(KMAX+2), where NCLD is number of cloudy DLGP
    !     dooz......dooz   = (.NOT. noz) where : noz = .FALSE.
    !               (do ozone computation )
    !==========================================================================
    !
    INTEGER, INTENT(in   ) :: ncols
    INTEGER, INTENT(in   ) :: kmax
    REAL(KIND=r8),    INTENT(in   ) :: scosc (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: cmuc  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: csmcld(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: dscld (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: rvbc  (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: rvdc  (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: rnbc  (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: rndc  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: agvcd (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: agncd (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: rsvcd (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: rsncd (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: sc    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: rco   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: rcg   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: taut  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: rc2   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: tr1   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: tr2   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: tr3   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: ta    (ncols*(kmax+1))
    REAL(KIND=r8),    INTENT(inout  ) :: wa    (ncols*(kmax+1))
    REAL(KIND=r8),    INTENT(inout  ) :: oa    (ncols*(kmax+1))
    REAL(KIND=r8),    INTENT(inout  ) :: e0    (ncols*(kmax+1))
    REAL(KIND=r8),    INTENT(inout  ) :: pu    (ncols*(kmax+2))
    REAL(KIND=r8),    INTENT(in   ) :: ozcd  (ncols*(kmax+2))
    REAL(KIND=r8),    INTENT(inout  ) :: swale (ncols*(kmax+2))
    REAL(KIND=r8),    INTENT(inout  ) :: swil  (ncols*(kmax+1))
    REAL(KIND=r8),    INTENT(inout  ) :: css   (ncols*(kmax+1))
    REAL(KIND=r8),    INTENT(inout  ) :: acld  (ncols*(kmax+1))
    REAL(KIND=r8),    INTENT(inout  ) :: dpc   (ncols*(kmax+1))
    REAL(KIND=r8),    INTENT(in   ) :: swilc (ncols*(kmax+1))
    REAL(KIND=r8),    INTENT(inout  ) :: ccu   (ncols*(kmax+1))
    REAL(KIND=r8),    INTENT(in   ) :: tauc  (ncols*(kmax+1))
    INTEGER, INTENT(in   ) :: litd  (ncols*(kmax+2))
    REAL(KIND=r8),    INTENT(in   ) :: sqrt3
    REAL(KIND=r8),    INTENT(in   ) :: gg
    REAL(KIND=r8),    INTENT(in   ) :: ggp
    REAL(KIND=r8),    INTENT(in   ) :: ggsq
    INTEGER, INTENT(in   ) :: np
    INTEGER, INTENT(in   ) :: lmp1
    INTEGER, INTENT(in   ) :: nsol
    INTEGER, INTENT(in   ) :: ncld
    INTEGER, INTENT(in   ) :: nclmp1
    INTEGER, INTENT(in   ) :: ncldnp
    LOGICAL, INTENT(in   ) :: dooz



    REAL(KIND=r8)                 :: w     (ncols*(kmax+2))
    REAL(KIND=r8)                 :: dn    (ncols*(kmax+2))
    REAL(KIND=r8)                 :: up    (ncols*(kmax+2))
    REAL(KIND=r8)                 :: vm    (ncols*(kmax+2))
    REAL(KIND=r8)                 :: vp    (ncols*(kmax+2))
    REAL(KIND=r8)                 :: cm    (ncols*(kmax+2))
    REAL(KIND=r8)                 :: dl    (ncols*(kmax+2))
    REAL(KIND=r8)                 :: den   (ncols*(kmax+1))

    REAL(KIND=r8)                 :: f2   (ncols*(kmax+1))
    REAL(KIND=r8)                 :: cr   (ncols*(kmax+2))
    REAL(KIND=r8)                 :: e2   (ncols*(kmax+1))
    REAL(KIND=r8)                 :: ul   (ncols*(kmax+2))
    REAL(KIND=r8)                 :: u1   (ncols*(kmax+1))
    REAL(KIND=r8)                 :: ak   (ncols*(kmax+1))
    REAL(KIND=r8)                 :: da   (ncols*(kmax+1))
    REAL(KIND=r8)                 :: db   (ncols*(kmax+1))
    REAL(KIND=r8)                 :: alf1 (ncols*(kmax+1))
    REAL(KIND=r8)                 :: alf2 (ncols*(kmax+1))
    REAL(KIND=r8)                 :: sol  (ncols*(kmax+2))

    REAL(KIND=r8) :: avbc (ncols)
    REAL(KIND=r8) :: anbc (ncols)
    REAL(KIND=r8) :: rvnbc(ncols)

    REAL(KIND=r8) :: expcut
    REAL(KIND=r8) :: ggpp1
    REAL(KIND=r8) :: ggpm1
    REAL(KIND=r8) :: gsqm1

    INTEGER :: lp
    INTEGER :: i
    INTEGER :: l
    INTEGER, DIMENSION(lmp1) :: l0a
    INTEGER, DIMENSION(np)   :: l0b
    INTEGER, DIMENSION(lmp1) :: l0c
    INTEGER, DIMENSION(lmp1) :: l1a
    INTEGER, DIMENSION(np)   :: l1b
    INTEGER, DIMENSION(lmp1) :: l1c
    INTEGER :: ik
    INTEGER :: k


    !----------------------------------------------------------------------C
    !-----PARAMETERS USING ASYMMETRY FACTOR OF CLOUD PARTICLE--------------C
    !-----------------------------------PHASE FUNCTION g=0.85--------------C
    !----------------------------------------------------------------------C
    ggpp1 = 1.0_r8 + ggp         !  1 + (g/(1+g))*0.75_r8
    ggpm1 = ggp - 1.0_r8         ! -1 + (g/(1+g))*0.75_r8
    gsqm1 = 1.0_r8 - ggsq        !  1 - g*g
    expcut=- LOG(1.0e53_r8)
    lp = lmp1 * ncld
    DO i= 1, ncld
       dn(i)    = 0.0_r8
       tr1(i)   = 0.0_r8
       rvnbc(i) = 0.0_r8
    END DO
    !----------------------------------------------------------------------C
    !-------DO OZONE ABSORPTION COMPUTATIONS ON THE WAY DOWN---------------C
    !----------------------------------------------------------------------C
    IF (dooz) THEN

       DO l = 1, lmp1
          l0a(l) = (l-1) * ncld
          l1a(l) =  l    * ncld
       END DO
       DO l = 2, np
          l0b(l) = (l-2) * ncld
          l1b(l) = (l-1) * ncld
       END DO
       DO k = 1, lmp1
          l  = np-k
          l0c(k) = (l-1) * ncld
          l1c(k) =  l    * ncld
       END DO
       DO l = 2, np
          DO i = 1, ncld
             w(l0b(l)+i) = ozcd(l1b(l)+i) * csmcld(i)  !Ozone path on the way down
          END DO
       END DO
       !----------------------------------------------------------------------C
       !--COMPUTE OZONE  ABSORPTION (DN) TILL THE LAYER ON THE WAY DOWN-------C
       !--USING OZONE ABSORPTION FUNCTION DUE TO Lacis and Hansen (1974)------C
       !--OZONE AMOUNT BELOW TOP CLOUD LAYER IS EQUAL ZERO--------------------C
       !----------------------------------------------------------------------C
       DO i = 1, nclmp1
          e0(i) = 103.63_r8 * w(i)
          e0(i) = e0(i)  * e0(i) * e0(i)
          vp(i) = 1.0_r8 + 138.57_r8   * w(i)
          cm(i) = vp(i) ** 0.805_r8
          dn(ncld+i) = 0.02118_r8 * w(i) &
               / (1.0_r8 + 0.042_r8 * w(i) &
               + 0.000323_r8 * w(i) * w(i)) &
               + 1.08173_r8  * w(i) &
               / cm(i) + 0.0658_r8  * w(i) &
               / (1.0_r8   + e0(i))
       END DO
       !
       !     compute ozone heating on the way down
       !
       !----------------------------------------------------------------------C
       !--COMPUTE OZONE ABSORPTION (ACLD) ON THE WAY DOWN IN EACH LAYER-------C
       !--BELOW TOP CLOUD LAYER OZONE ABSORPTION IS EQUAL ZERO----------------C
       !----------------------------------------------------------------------C
       DO i = 1, nclmp1
          acld(i) = dn(ncld+i) - dn(i)
       END DO

       DO l = 1, lmp1
          DO i = 1, ncld
             acld(l0a(l)+i) = acld(l0a(l)+i) * scosc(i)
          END DO
       END DO
    ELSE
       DO i = 1, ncld
          dn(lp+i) = 0.0_r8
       END DO
    END IF
    !
    !     downflux and ground absorption for lambda under 0.7 microns
    !     as above but for lamda between 0.7 and 0.9 microns
    !
    !----------------------------------------------------------------------C
    !-----------DOWNWARD FLUXES AT GROUND AND GROUND ABSORPTION------------C
    !--FOR WAVELENGTH LESS THAN 0.7 MICRONS--------------------------------C
    !--AS ABOVE AND FOR WAVELENGTH BETWEEN 0.7 AND 0.9 MICRONS-------------C
    !--TAKING INTO ACCOUNT OZONE ABSORPTION, MOLECULAR SCATTERING AND------C
    !--GROUND REFLECTION.--------------------------------------------------C
    !--0.5 is the part of Solar Constant in visible region of spectrum-----C
    !--0.147 is the part of Solar Constant in spectral region--------------C
    !---------from 0.7 to 0.9 mcm------------------------------------------C
    !--TR2(I)=EXP(-(1-g*g)TAUcld*M), TAUcld is CLOUD OPTICAL DEPTH,--------C
    !-M IS MAGNIFICATION FACTOR, g IS ASYMMETRY FACTOR OF CLOUD PARTICLES--C
    !--TR3(I)=EXP(-(1-g*g)(TAUcld+TAURay)*M), TAURay=0.15746,--------------C
    !--RSVCD...GROUND VISIBLE BEAM ALBEDO----------------------------------C
    !--AGVCD...GROUND VISIBLE DIFFUSE ALBEDO-------------------------------C
    !--RSNCD...GROUND NEAR INFRARED BEAM ALBEDO----------------------------C
    !--AGNCD...GROUND NEAR INFRARED DIFFUSE ALBEDO-------------------------C
    !--RCG.....CLOUDY AND GROUND REFLECTION FOR GROUND ABSORPTION----------C
    !----------COMPUTATION IN VISIBLE REGION OF SPECTRUM-------------------C
    !--RC2.....CLOUDY AND GROUND REFLECTION FOR GROUND ABSORPTION----------C
    !----------COMPUTATION IN REGION OF SPECTRUM FROM 0.7 TO 0.9 mcm-------C
    !----------------------------------------------------------------------C
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !----------------------------------------------------------------------C
    !-----------SPECTRAL REGION FROM 0.2 TO 0.7 mcm------------------------C
    !----------------------------------------------------------------------C
    DO i = 1, ncld
       sc(i)   = 0.5_r8     - dn(lp+i) ! (DN)..Total ozone absorption=0?
       taut(i) = sc(i)   * tr3(i) * scosc(i)  !Downward visible beam
       avbc(i) = taut(i) * (1.0_r8 - rsvcd(i))  !Net visible beam
       rvbc(litd(i)) = taut(i)
    END DO
    DO i = 1, ncld
       sc(i)   = sc(i)   * (1.0_r8 - rcg(i)) * scosc(i)
       taut(i) = (sc(i)  - avbc(i)) &
            / (1.0_r8    - agvcd(i))     !Net visible diffuse
       rvdc(litd(i)) = taut(i)
    END DO
    DO i = 1, ncld
       !rnbc(i) = 0.147_r8   * tr2(i)         * scosc(i)!Downward NearIR beam
       !anbc(i) = (1.0_r8    - rsncd(i))      * rnbc(i) !Net NearIR beam
       !rndc(i) = (0.147_r8  * (1.0_r8 - rc2(i)) * scosc(i) &!Downward NearIR dif.
       !     - anbc(i))/ (1.0_r8 - agncd(i))
       !sc(i)   = sc(i)   + (1.0_r8 - rc2(i)) * scosc(i) &!Ground absorption
       !     * 0.147_r8

       ! Include global version 2.2
       rnbc(i) = 0.1214_r8  * tr2(i)         * scosc(i)   !Downward NearIR beam
       anbc(i) = (1.0_r8    - rsncd(i))      * rnbc(i)    !Net NearIR beam
       rndc(i) = (0.1214_r8  * (1.0_r8 - rc2(i)) * scosc(i) &!Downward NearIR dif.
            - anbc(i))/ (1.0_r8 - agncd(i))
       sc(i)   = sc(i)   + (1.0_r8 - rc2(i)) * scosc(i) & !Ground absorption
            * 0.1214_r8
    END DO
    !----------------------------------------------------------------------C
    !------DO OZONE ABSORPTION COMPUTATIONS ON THE WAY DOWN AND UP---------C
    !----------------------------------------------------------------------C
    IF (dooz) THEN
       !----------------------------------------------------------------------C
       !---COMPUTE OZONE PATH (W) TILL THE LAYER  ON THE WAY DOWN AND UP------C
       !----------------------------------------------------------------------C
       DO i = 1, ncld
          up(lp+i) = dn(lp+i)
          w(lp+i)  = ozcd(lp+i) * (1.9_r8 + csmcld(i))
       END DO
       DO i = 1, nclmp1
          w(i) = 1.9_r8 * ozcd(i)
       END DO

       DO l = 1, lmp1
          !cdir nodep
          DO i = 1, ncld
             w(l0a(l)+i) = w(lp+i) - w(l0a(l)+i)
          END DO
       END DO
       !----------------------------------------------------------------------C
       !--OZONE  ABSORPTION (UP) TILL THE LAYER ON THE WAY DOWN AND UP--------C
       !--USING OZONE ABSORPTION FUNCTION DUE TO Lacis and Hansen (1974)------C
       !----------------------------------------------------------------------C
       DO i = 1, nclmp1
          e0(i) = 103.63_r8       * w(i)
          e0(i) = e0(i) * e0(i) * e0(i)
          vp(i) = 1.0_r8 + 138.57_r8 * w(i)
          cm(i) = vp(i) ** 0.805_r8
          up(i) = 0.02118_r8 * w(i) &
               / (1.0_r8 + 0.042_r8 * w(i) &
               + 0.000323_r8 * w(i) * w(i)) &
               + 1.08173_r8  * w(i) &
               / cm(i) + 0.0658_r8  * w(i) &
               / (1.0_r8 + e0(i))
       END DO
       !
       !     compute ozone heating on the way up
       !
       !----------------------------------------------------------------------C
       !-----OZONE ABSORPTION (W) IN EACH LAYER ON WAY DOWN AND UP,-----------C
       !-----TOTAL OZONE ABSORPTION (ACLD) IN EACH LAYER----------------------C
       !--RCO....CLOUDY AND GROUND REFLECTION FOR OZONE ABSORPTION COMP.------C
       !----------------------------------------------------------------------C
       DO i = 1, nclmp1
          w(i) = up(i) - up(ncld+i)
       END DO
       DO l = 1, lmp1
          DO i=1,ncld
             acld(l0a(l)+i) = acld(l0a(l)+i) + w(l0a(l)+i) * rco(i) &
                  * scosc(i)
          END DO
       END DO
       acld(1:nclmp1) = MAX(0.0_r8,acld(1:nclmp1))
    END IF
    DO l = 1, lmp1
       DO i=1,ncld
          oa(l0a(l)+i)  = cmuc(i)          !Cosine of solar zenith angle
       END DO
    END DO
    DO i = 1, nclmp1                     !NCLD*(KMAX+1)
       ta(i) = 0.5_r8  - ggp * oa(i)        !0.5_r8-gD*0.75_r8*mu0
       wa(i) = 1.0_r8  - ta(i)              !gD=g/(1+g)
    END DO
    !----------------------------------------------------------------------C
    !-------CYCLE OVER TEN WATER VAPOR ABSORPTION COEFFICIENTS (XK)--------C
    !-----------(TAKING FROM EXPONENTIAL EXPANSION OF WATER VAPOR----------C
    !---------- ABSORPTION FUNCTION WITH THE WEIGHTS (FK))-----------------C
    !----------------------------------------------------------------------C
    DO ik=1,nwaterbd

       !
       !     compute total optical depth and single scattering albedo
       !     pizero of clouds is 0.99
       !
       DO i = 1, nclmp1
          ccu(i) = tauc(i) + swilc(i) * xk(ik)
          cr(i)  = tauc(i) / ccu(i)   * 0.99_r8
          !
          !     compute reflection and transmission for diffuse incidence
          !-----IN TWO STREAM APPROACH DUE TO Sagan and Pollack (1967)-----------C
          !-----GG....g=0.85  ASYMMETRY FACTOR-----------------------------------C
          !
          vp(i)  = 1.0_r8 - cr(i) * gg
          vm(i)  = 1.0_r8 - cr(i)
       END DO
       DO i = 1, nclmp1
          vp(i)  = SQRT(vp(i))
          vm(i)  = SQRT(vm(i))
          u1(i)  = vp(i) / vm(i)
          w(i)   = u1(i)  + 1.0_r8
          w(i)   = w(i)   * w(i)
          cm(i)  = u1(i)  - 1.0_r8
          cm(i)  = cm(i)  * cm(i)
          dpc(i) = sqrt3 * vp(i) * vm(i) * ccu(i)
          e2(i)  = -dpc(i)
          IF (ik >= 4) e2(i) = MAX(e2(i),expcut)
          e2(i)  = EXP(e2(i))
          cm(i)  = w(i)   - cm(i) * e2(i) * e2(i)
          pu(i)  = (u1(i) * u1(i) - 1.0_r8) &           ! Reflection
               * (1.0_r8 - e2(i) * e2(i)) / cm(i)
          css(i) =  4.0_r8 * u1(i) * e2(i)  / cm(i)     ! Transmission
          !----------------------------------------------------------------------C
          !-----DELTA-SCALING OF OPTICAL DEPTH AND SINGLE SCATTERING ALBEDO------C
          !----------------------------------------------------------------------C
          vp(i)  = 1.0_r8 - ggsq    * cr(i)
          ccu(i) = ccu(i) * vp(i)                    ! Optical depth
          cr(i)  = cr(i)  * gsqm1 / vp(i)            ! Single scattering albedo
       END DO
       !----------------------------------------------------------------------C
       !----COMPUTE DIRECT TRANSMISSION (EO) FOR DELTA-SCALING OPTICAL DEPTH--C
       !----------------------------------------------------------------------C
       DO i = 1, ncld
          sol(i)  = scosc(i) * fk(ik)
       END DO
       DO l = 1, lmp1
          DO i=1,ncld
             dpc(l0a(l)+i) = ccu(l0a(l)+i) * csmcld(i) ! Delta Tau*Magnif. factor
          END DO
       END DO
       DO i = 1, nclmp1
          e0(i)  = -dpc(i)
          IF (ik >= 4) e0(i) = MAX(e0(i), expcut)
          e0(i)  = EXP(e0(i))
       END DO
       !
       !     compute diffuse reflection and transmission for direct beam
       !     set up delta-eddington parameters
       !
       DO l = 1, lmp1
          !cdir nodep
          DO i = 1, ncld
             sol(l1a(l)+i) = sol(l0a(l)+i)  * e0(l0a(l)+i) ! Beam rad. at all levels
          END DO
       END DO
       DO i = 1, ncld
          rvnbc(i) = rvnbc(i) + sol(lp+i)  ! Beam rad. at the ground
       END DO

       DO i = 1, nclmp1
          swil(i)  =  1.75_r8 - cr(i) * ggpp1  ! 1.75-omegaD*(1+gD*0.75)
          swale(i) = -0.25_r8 - cr(i) * ggpm1  ! -0.25- omegaD*(gD*0.75-1)
          ak(i)    =  swil(i)  * swil(i) &
               -  swale(i) * swale(i)
       END DO
       DO i = 1, nclmp1
          ak(i) = SQRT(ak(i))
          alf1(i)  = swil(i)  * wa(i) &
               + swale(i) * ta(i)
          alf2(i)  = swil(i)  * ta(i) &
               + swale(i) * wa(i)
          vp(i) = 1.0_r8 - ak(i) * oa(i)
          vm(i) = 1.0_r8 + ak(i) * oa(i)
          ul(i) = ta(i)   - alf2(i) * oa(i)
          u1(i) = wa(i)   + alf1(i) * oa(i)
          cm(i) = alf2(i) + ak(i)   * ta(i)
          dl(i) = alf2(i) - ak(i)   * ta(i)
          e2(i) = alf1(i) - ak(i)   * wa(i)
          f2(i) = alf1(i) + ak(i)   * wa(i)
          den(i)  = 1.0_r8 - ak(i)  * ak(i) &
               *       oa(i)  * oa(i)
       END DO
       DO i = 1, nclmp1
          den(i) = MAX(den(i),.000001e0_r8)
          dpc(i) = ak(i) * ccu(i)
          da(i)  = -dpc(i)
          IF (ik >= 4) da(i) = MAX(da(i), expcut)
          da(i) = EXP(da(i))
          dpc(i) = 2.0_r8 * dpc(i)
          db(i)  = -dpc(i)
          IF (ik >= 4) db(i) = MAX(db(i), expcut)
          db(i) = EXP(db(i))
          den(i) = den(i) * (ak(i)   + swil(i) &
               + (ak(i) - swil(i)) * db(i))
          !
          !     compute upward and downward fluxes at layer boundaries
          !
          up(i) = cr(i)  * (vp(i) * cm(i) &
               - vm(i)  * dl(i)  * db(i) &
               - 2.0_r8           * ak(i)  * ul(i) &
               * e0(i)  * da(i)) / den(i)
          up(i) = up(i)  * sol(i)
          dn(ncld+i) = -cr(i)  * (vm(i) * f2(i) &
               * e0(i)   - vp(i)  * e2(i) &
               * e0(i)   * db(i)  - 2.0_r8 &
               * ak(i)   * u1(i)  * da(i)) &
               / den(i)
          dn(ncld+i) = sol(i)  * dn(ncld+i)
       END DO
       !
       !     fill up boundary terms
       !
       DO i = 1, ncld
          up(lp+i) = rsncd(i) * sol(lp+i)  !RSNCD...ground NIR beam albedo
          pu(lp+i) = agncd(i)              !AGNCD...ground vis. beam albedo
          dn(i) = 0.0_r8
          cr(i) = 0.0_r8
          cm(i) = 1.0_r8
          vp(i) = 0.0_r8
          vm(i) = up(i)
       END DO
       !
       !     safety net
       !
       DO i = 1, ncldnp
          up (i) = MAX(up (i), 0.0_r8)
          dn (i) = MAX(dn (i), 0.0_r8)
          sol(i) = MAX(sol(i), 0.0_r8)
       END DO
       DO i = 1, nclmp1
          pu (i) = MAX(pu (i), 0.0_r8)
          css(i) = MAX(css(i), 0.0_r8)
       END DO
       !----------------------------------------------------------------------C
       !-------COMPUTE FLUXES AT THE BOUNDARIES OF ADDING LAYERS--------------C
       !----------------------------------------------------------------------C
       !---------COMPUTE MAGNIFICATION FACTORS....CM=1/(1-R1*R2)--------------C
       !----------------------------------------------------------------------C
       DO i = 1, nclmp1
          e0(i) = css(i) * css(i)
       END DO
       DO l = 1, lmp1
          !cdir nodep
          DO i = 1, ncld
             cr(l1a(l)+i) = pu(l0a(l)+i) + e0(l0a(l)+i) * cr(l0a(l)+i) &
                  * cm(l0a(l)+i)
             cm(l1a(l)+i) = 1.0_r8      - cr(l1a(l)+i) * pu(l1a(l)+i)
             cm(l1a(l)+i) = 1.0_r8      / cm(l1a(l)+i)
          END DO
       END DO
       !
       !     compute fluxes on the way down
       !
       DO i = 1, nclmp1
          db(i) = cr(i)  * up(i)
          e0(i) = css(i) * cm(i)
       END DO
       DO l = 1, lmp1
          !cdir nodep
          DO i = 1, ncld
             vp(l1a(l)+i) = e0(l0a(l)+i) * (vp(l0a(l)+i) + db(l0a(l)+i)) &
                  + dn(l1a(l)+i)
          END DO
       END DO
       DO i = 1, nclmp1
          vm(ncld+i) = cm(ncld+i) *  &
               (vp(ncld+i) * pu(ncld+i) + up(ncld+i))
       END DO
       !----------------------------------------------------------------------C
       !-------COMPUTATION OF THE FLUXES AT THE GROUND------------------------C
       !----------------------------------------------------------------------C
       DO i = 1, ncld
          ul(lp+i) = vm(lp+i)
          dl(lp+i) = vp(lp+i)  + ul(lp+i) * cr(lp+i) + sol(lp+i)
          tr1(i)   = tr1(i)    + dl(lp+i) - ul(lp+i)
       END DO
       !
       !     compute fluxes on the way up
       !
       DO i = 1, nclmp1
          e0(i) = css(i) * cm(i)
       END DO

       DO k = 1, lmp1
          !cdir nodep
          DO i = 1, ncld
             ul(l0c(k)+i) = vm(l0c(k)+i) + ul(l1c(k)+i) * e0(l0c(k)+i)
          END DO
       END DO
       DO i = 1, nclmp1
          dl(i)   = vp(i)   + ul(i) * cr(i) + sol(i)
       END DO
       DO i = 1, nclmp1
          acld(i) = acld(i) + dl(i) - ul(i) &
               - (dl(ncld+i)     -   ul(ncld+i))
       END DO
    END DO
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !CC-END OF CYCLE OVER FIVE WATER VAPOR ABSORPTION COEFFICIENTS (XK)-CCCC
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    DO i = 1, ncld
       taut(i) = rnbc(i) + rvnbc(i)
    END DO
    !----------------------------------------------------------------------C
    !----RNBC...DOWNWARD NIR BEAM CLOUDY GROUND FLUX-----------------------C
    !----------------------------------------------------------------------C
    DO i = 1, ncld
       rnbc(litd(i)) = taut(i)
    END DO
    DO i = 1, ncld
       taut(i) = rndc(i) + (tr1(i) - rvnbc(i) &
            * (1.0_r8 - rsncd(i))) / (1.0_r8 - agncd(i))
    END DO
    !----------------------------------------------------------------------C
    !----RNDC...DOWNWARD NIR DIFFUSE CLOUDY GROUND FLUX--------------------C
    !----------------------------------------------------------------------C
    DO i = 1, ncld
       rndc(litd(i)) = taut(i)
    END DO
    DO i = 1, ncld
       taut(i) = sc(i) + tr1(i)
    END DO
    !----------------------------------------------------------------------C
    !-----SC...GROUND ABSORPTION-------------------------------------------C
    !----------------------------------------------------------------------C
    DO i = 1, ncld
       sc(litd(i)) = taut(i)
    END DO
    !
    !     add up absorption in column
    !
    DO i = 1, nclmp1
       e0(i) = acld(i)
    END DO
    DO i = 1, nclmp1
       acld(litd(i)) = e0(i)
    END DO
    !----------------------------------------------------------------------C
    !-----DSCLD...TOTAL ABSORPTION IN ATMOSPHERE AND GROUND----------------C
    !----------------------------------------------------------------------C
    DO i = 1, nsol
       dscld(i) = sc(i)
    END DO
    DO l = 1, lmp1
       DO i = 1, nsol
          dscld(i) = dscld(i) + acld(l0a(l)+i)
       END DO
    END DO
  END SUBROUTINE cloudy

  !-----------------------------------------------------------------------
  ! Subroutine: CLEAR
  !
  ! CALCULATES:
  !     OZONE PATH TILL THE LAYER ILLUMINATED FROM ABOVE AND FROM BELOW
  !     OZONE ABSORPTION IN EACH LAYER AND THE GROUND
  !     WATER VAPOR  PATH TILL THE LAYER ILLUMINATED FROM ABOVE AND FROM BELOW
  !     WATER VAPOR ABSORPTION IN EACH LAYER AND THE GROUND
  !     DOWNWARD AND UPWARD FLUXES AT THE GROUND (VISIBLE AND NEAR INFRARED)
  !
  ! MAIN OUTPUT:
  !    DOWNWARD GROUND FLUXES AND ABSORPTIONS, IN CLEAR CASE AND
  !    IN ALL DAYTIME LATITUDE GRID POINTS (NSOL)
  !
  ! COMPUTATION IS DONE 10 BANDS
  !---------------------------------------------------------------------C
  SUBROUTINE clear(ncols ,kmax  ,sqrt3 ,np    ,lmp1  ,nsol  ,nslmp1,nsolnp, &
       dooz  ,scosz ,cosmag,dsclr ,rvbl  ,rvdl  ,rnbl  ,rndl  , &
       agv   ,rsurfv,rsurfn,sl    ,rlo   ,rlg   ,tr1   ,e0    , &
       ozale ,swale ,aclr                                       )
    !==========================================================================
    !    imax......Number of grid points on a gaussian latitude circle
    !    kmax......Number of grid points at vertical
    !    sqrt3.....Magification factor for diffuse reflected radiation
    !              sqrt3  = SQRT(3.0)
    !    np........np = (kmax+2)
    !    lmp1......lmp1   = (kmax+1)
    !    nsol......Number of solar latitude grid points (cosz>0.01)
    !    nslmp1....nsol*(kmax +1 ), where nsol is number of solar
    !              latitude grid points (cosz>0.01)
    !    nsolnp....nsol*(kmax +2 ), where nsol is number of solar
    !              latitude grid points (cosz>0.01)
    !    dooz......dooz   = (.NOT. noz) where : noz = .FALSE.
    !               (do ozone computation )
    !    scosz.....scosz(i)   = s0     * cmu(i)
    !                         where s0  is constant solar
    !                               cmu is cosine of solar zenith angle
    !    cosmag....Magnification factor in DLGP
    !              csmcld = cosmag(i)  = 1224.0 * cmu(i) * cmu(i) + 1.0
    !    dsclr.....Absorption of clear atmosphere and ground clear
    !    rvbl......Visible beam clear  "Downward ground fluxes in DLGP"
    !    rvdl......Visible diffuse clear "Downward ground fluxes in DLGP"
    !    rnbl......NearIR beam clear "Downward ground fluxes in DLGP"
    !    rndl......NearIR diffuse clear  "Downward ground fluxes in DLGP"
    !    agv.......Ground visible diffuse albedo in DLGP
    !    rsurfv....Ground visible beam albedo in DLGP
    !    rsurfn....Ground near IR beam albedo in DLGP
    !    sl........Net ground clear flux
    !    rlo.......Clear sky and ground reflection in DLGP for
    !              ozone computation
    !    rlg.......Clear sky and ground reflection in DLGP for
    !              ground absorption computation
    !    tr1.......Extinction of beam radiation (clear)
    !    e0........Cloud optical depth in DLGP
    !    ozale.....Ozone amount in column in DLGP
    !    swale.....water vapor amount in column in DLGP
    !    aclr......Absorption in clear ATM
    !==========================================================================
    INTEGER, INTENT(in   ) :: ncols
    INTEGER, INTENT(in   ) :: kmax
    REAL(KIND=r8),    INTENT(in   ) :: sqrt3
    INTEGER, INTENT(in   ) :: np
    INTEGER, INTENT(in   ) :: lmp1
    INTEGER, INTENT(in   ) :: nsol
    INTEGER, INTENT(in   ) :: nslmp1
    INTEGER, INTENT(in   ) :: nsolnp
    LOGICAL, INTENT(in   ) :: dooz
    REAL(KIND=r8),    INTENT(in   ) :: scosz (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: cosmag(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: dsclr (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: rvbl  (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: rvdl  (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: rnbl  (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: rndl  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: agv   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: rsurfv(ncols)
    REAL(KIND=r8),    INTENT(in   ) :: rsurfn(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: sl    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: rlo   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: rlg   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: tr1   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: e0    (ncols*(kmax+1))
    REAL(KIND=r8),    INTENT(in   ) :: ozale (ncols*(kmax+2))
    REAL(KIND=r8),    INTENT(in   ) :: swale (ncols*(kmax+2))
    REAL(KIND=r8),    INTENT(inout  ) :: aclr  (ncols*(kmax+1))

    REAL(KIND=r8)                  :: den   (ncols*(kmax+1))
    REAL(KIND=r8)                  :: w     (ncols*(kmax+2))
    REAL(KIND=r8)                  :: dn    (ncols*(kmax+2))
    REAL(KIND=r8)                  :: up    (ncols*(kmax+2))
    REAL(KIND=r8)                  :: vp    (ncols*(kmax+2))
    REAL(KIND=r8)                  :: cm    (ncols*(kmax+2))

    REAL(KIND=r8) :: avbl(ncols)
    REAL(KIND=r8) :: dum(1)
    REAL(KIND=r8) :: expcut

    INTEGER :: i
    INTEGER :: l
    INTEGER :: l0
    INTEGER :: l1
    INTEGER :: lp
    INTEGER :: ik

    expcut=- LOG(1.0e53_r8)
    lp = lmp1 * nsol           ! (KMAX+1)*NSOL
    DO i = 1, nsol
       sl(i) = 0.0_r8
       dn(i) = 0.0_r8
    END DO
    DO i = 1, nslmp1           ! NSOL*(KMAX+1)
       aclr(i) = 0.0_r8
    END DO
    DO l = 1, lmp1             ! (KMAX+1)
       l0 = (l-1) * nsol
       DO i = 1, nsol
          !----------------------------------------------------------------------C
          !---SOLAR DOWNWARD FLUX (DEN) WITHOUT ANY EXTINCTION AT ALL LEVELS-----C
          !----------------------------------------------------------------------C
          den(l0+i) = scosz(i)
       END DO
    END DO
    !----------------------------------------------------------------------C
    !-------DO OZONE ABSORPTION COMPUTATIONS ON THE WAY DOWN---------------C
    !----------------------------------------------------------------------C
    IF (dooz) THEN
       DO l = 2, np            ! (KMAX+2)
          l0 = (l-2) * nsol
          l1 = (l-1) * nsol
          DO i = 1, nsol
             !----------------------------------------------------------------------C
             !---COMPUTE OZONE PATH (W) TILL THE LAYER ON THE WAY DOWN--------------C
             !----------------------------------------------------------------------C
             w(l0+i) = ozale(l1+i) * cosmag(i)   !Magnification factor
          END DO
       END DO
       !----------------------------------------------------------------------C
       !--COMPUTE OZONE  ABSORPTION (DN) TILL THE LAYER ON THE WAY DOWN-------C
       !--USING OZONE ABSORPTION FUNCTION DUE TO Lacis and Hansen (1974)------C
       !----------------------------------------------------------------------C
       DO i = 1, nslmp1        ! NSOL*(KMAX+1)
          e0(i) = 103.63_r8       * w(i)
          e0(i) = e0(i) * e0(i) * e0(i)
          vp(i) = 1.0_r8 + 138.57_r8 * w(i)
          cm(i) = EXP(LOG(vp(i) + 1.0e-100_r8)* 0.805_r8)
       END DO
       DO i = 1, nslmp1        ! NSOL*(KMAX+1)
          dn(nsol+i) = 0.02118_r8 * w(i) &
               / (1.0_r8 + 0.042_r8 * w(i) &
               + 0.000323_r8 * w(i) * w(i)) &
               + 1.08173_r8  * w(i) &
               / cm(i) + 0.0658_r8  * w(i) &
               / (1.0_r8   + e0(i))
       END DO
       !----------------------------------------------------------------------C
       !--COMPUTE OZONE ABSORPTION (ACLR) ON THE WAY DOWN IN EACH LAYER-------C
       !----------------------------------------------------------------------C
       DO i = 1, nslmp1                  !NSOL*(KMAX+1)
          aclr(i) = (dn(nsol+i) - dn(i)) * den(i)
       END DO
    ELSE
       DO i = 1, nsol
          dn(lp+i) = 0.0_r8
       END DO
    END IF
    !
    !     downflux and ground absorption for lamda less than 0.7 microns
    !     as above but for lambda between 0.7 and 0.9 microns
    !
    !----------------------------------------------------------------------C
    !-----DOWNWARD AND UPWARD FLUXES AT GROUND AND GROUND ABSORPTION-------C
    !--FOR WAVELENGTH LESS THAN 0.7 MICRONS--------------------------------C
    !--AS ABOVE AND FOR WAVELENGTH BETWEEN 0.7 AND 0.9 MICRONS-------------C
    !--TAKING INTO ACCOUNT OZONE ABSORPTION, MOLECULAR SCATTERING AND------C
    !--GROUND REFLECTION.--------------------------------------------------C
    !--0.5 is the part of Solar Constant in visible region of spectrum-----C
    !--0.147 is the part of Solar Constant in spectral region--------------C
    !---------from 0.7 to 0.9 mcm------------------------------------------C
    !--TR1(I)=EXP(-TAURAY*M), TAURAY=0.15746, M is Magnification Factor----C
    !--RSURFV..GROUND VISIBLE BEAM ALBEDO----------------------------------C
    !--RLG.....CLEAR SKY AND GROUND REFLECTION FOR GROUND ABSORPTION COMP.-C
    !--AGV.....GROUND VISIBLE DIFFUSE ALBEDO-------------------------------C
    !--RSURFN..GROUND NEAR INFRARED BEAM ALBEDO----------------------------C
    !----------------------------------------------------------------------C
    DO i = 1, nsol
       sl(i)   = 0.5_r8        - dn(lp+i)          !LP=NSOL*(KMAX+1)
       rvbl(i) = sl(i) * tr1(i) * scosz(i)      !Downward visible beam
       avbl(i) = (1.0_r8 - rsurfv(i)) * rvbl(i)    !Net visible beam
       sl(i)   = (1.0_r8 - rlg(i))    * scosz(i) * sl(i)
       rvdl(i) = (sl(i) - avbl(i)) / (1.0_r8 - agv(i))   !Downward visible dif.
       !rnbl(i) = 0.147_r8 * scosz(i)                    !Downward NearIR beam
       !rndl(i) = 0.0_r8*sl(i)                           !Downward NearIR diffuse
       ! Include global version 2.2
       rnbl(i) = 0.1214_r8 * scosz(i)
       rndl(i) = 0.0_r8
       sl(i)   = sl(i) + rnbl(i)  * (1.0_r8 - rsurfn(i)) !Ground absorption
    END DO
    !----------------------------------------------------------------------C
    !------DO OZONE ABSORPTION COMPUTATIONS ON THE WAY DOWN AND UP---------C
    !----------------------------------------------------------------------C
    IF (dooz) THEN
       DO i = 1, nsol
          up(lp+i) = dn(lp+i)
          w(lp+i)  = ozale(lp+i) * (1.90_r8 + cosmag(i))
       END DO
       DO i = 1, nslmp1                !NSOL*(KMAX+1)
          w(i) = 1.9_r8 * ozale(i)
       END DO
       DO l = 1, lmp1                  !(KMAX+1)
          l0 = (l-1) * nsol
          !cdir nodep
          DO i = 1, nsol
             w(l0+i) = w(lp+i) - w(l0+i)
          END DO
       END DO
       DO i = 1, nslmp1                !NSOL*(KMAX+1)
          e0(i) = 103.63_r8 * w(i)
          e0(i) = e0(i)  * e0(i) * e0(i)
          vp(i) = 1.0_r8 + 138.57_r8   * w(i)
          cm(i) = EXP(LOG(vp(i) + 1.0e-100_r8)* 0.805_r8)
       END DO
       DO i = 1, nslmp1                !NSOL*(KMAX+1)
          up(i) = 0.02118_r8 * w(i) &
               / (1.0_r8 + 0.042_r8 * w(i) &
               + 0.000323_r8 * w(i) * w(i)) &
               + 1.08173_r8  * w(i) &
               / cm(i) + 0.0658_r8  * w(i) &
               / (1.0_r8   + e0(i))
       END DO
       !----------------------------------------------------------------------C
       !-----OZONE ABSORPTION (W) IN EACH LAYER ON WAY UP---------------------C
       !----------------------------------------------------------------------C
       DO i = 1, nslmp1
          w(i) = up(i) - up(nsol+i)
       END DO
       !----------------------------------------------------------------------C
       !-------COMPUTE TOTAL OZONE ABSORPTION (ACLR) IN EACH LAYER------------C
       !--RLO....CLEAR SKY AND GROUND REFLECTION FOR OZONE ABSORPTION COMP.---C
       !----------------------------------------------------------------------C
       DO l = 1, lmp1
          l0 = (l-1) * nsol
          DO i = 1, nsol
             aclr(l0+i) = aclr(l0+i) + w(l0+i) * rlo(i) &
                  * scosz(i)
          END DO
       END DO
       !-----------------------SET ACLR>0-------------------------------------C
       dum(1)=0.0e0_r8
       aclr(1:nslmp1) = MAX(aclr(1:nslmp1), 0.0_r8)
    END IF
    !----------------------------------------------------------------------C
    !-------CYCLE OVER TEN WATER VAPOR ABSORPTION COEFFICIENTS (XK)--------C
    !-----------(TAKING FROM EXPONENTIAL EXPANSION OF WATER VAPOR----------C
    !---------- ABSORPTION FUNCTION WITH THE WEIGHTS (FK))-----------------C
    !----------------------------------------------------------------------C
    DO ik=1,nwaterbd
       DO l = 1, np                        !KMAX+2
          l0 = (l-1) * nsol
          DO i = 1, nsol
             !----------------------------------------------------------------------C
             !--COMPUTE WATER VAPOR PATH TILL THE LAYER (W) ON THE WAY DOWN---------C
             !----------------------------------------------------------------------C
             w(l0+i) = swale(l0+i) * cosmag(i)
          END DO
       END DO
       IF (ik >= 4) THEN
          DO i = 1, nsolnp                 !NSOL*(KMAX+2)
             dn(i) = -w(i) * xk(ik)
             dn(i) = MAX(dn(i), expcut)
             dn(i) = EXP(dn(i))
          END DO
       ELSE
          DO i = 1, nsolnp                 !NSOL*(KMAX+2)
             !----------------------------------------------------------------------C
             !---WATER VAPOR ABSORPTION TILL THE LAYER (DN) ON THE WAY DOWN---------C
             !----------------------------------------------------------------------C
             dn(i) = EXP(-w(i) * xk(ik))
          END DO
       END IF
       !----------------------------------------------------------------------C
       !--WATER VAPOR ABSOPTION (W) ON THE WAY DOWN IN EACH LAYER-------------C
       !---ACLR..TOTAL OZONE AND WATER VAPOR  ABSORPTIONON THE WAY DOWN IN EACH LAYER
       !---DEN...DOWNWARD FLUX AT THE TOP AT ALL LEVELS
       !----------------------------------------------------------------------C
       DO i = 1, nslmp1             !NSOL*(KMAX+1)
          w(i)    = (dn(i)  - dn(nsol+i)) * fk(ik)
          w(i)    = w(i)    * den(i)
          aclr(i) = aclr(i) + w(i)
       END DO
       !----------------------------------------------------------------------C
       !-----GROUND ABSORPTION (SL) FOR WAVELENGTHS OVER 0.9 MICRONS.---------C
       !------RNBL.....DOWNWARD GROUND NEAR INFRARED DIRECT BEAM FLUX---------C
       !----------------------------------------------------------------------C
       DO i = 1, nsol
          sl(i)   = sl(i)   + dn(lp+i) * (1.0_r8 - rsurfn(i)) &
               * scosz(i) * fk(ik)
          rnbl(i) = rnbl(i) + dn(lp+i) * scosz(i) * fk(ik)
          !----------------------------------------------------------------------C
          !----COMPUTE WATER VAPOR ABSORPTION  ON THE WAY DOWN AND UP,-----------C
          !---UP....WATER VAPOR ABSOPTION ON THE WAY DOWN AND UP TILL THE LAYER--C
          !----------------------------------------------------------------------C
          up(lp+i) = dn(lp+i)
          !----------------------------------------------------------------------C
          !-----W.....WATER VAPOR PATH ON THE WAY DOWN AND UP--------------------C
          !-----MULTIPLIED BY  WATER VAPOR ABSOPTION COEFFICIENT,----------------C
          !-----SQRT3....MAGNIFICATION FACTOR FOR DIFFUSE REFLECTED RADIATION----C
          !----------------------------------------------------------------------C
          w(lp+i)  = swale(lp+i) * (sqrt3 + cosmag(i)) * xk(ik)
       END DO
       DO i = 1, nslmp1
          w(i) = sqrt3 * swale(i) * xk(ik)
       END DO
       DO l = 1, lmp1
          l0 = (l-1) * nsol
          !cdir nodep
          DO i = 1, nsol
             w(l0+i) = w(l0+i) - w(lp+i)
          END DO
       END DO
       up(1:nslmp1) = EXP(w(1:nslmp1))
       !----------------------------------------------------------------------C
       !--WATER VAPOR ABSOPTION (W) ON THE WAY UP IN EACH LAYER---------------C
       !----------------------------------------------------------------------C
       DO i = 1, nslmp1
          w(i) = (up(nsol+i) - up(i)) * fk(ik)
       END DO
       DO l = 1, lmp1
          l0 = (l-1) * nsol
          !cdir nodep
          DO i = 1, nsol
             !----------------------------------------------------------------------C
             !---ACLR......OZONE AND WATER VAPOR ABSORPTION IN EACH LAYER,----------C
             !---RSURFN....GROUND NEAR INFRARED DIRECT BEAM ALBEDO.-----------------C
             !----------------------------------------------------------------------C
             aclr(l0+i) = aclr(l0+i) + w(l0+i) * rsurfn(i) * scosz(i)
          END DO
       END DO
    END DO
    !----------------------------------------------------------------------C
    !----COMPUTE TOTAL GROUND AND CLEAR SKY ABSORPTION (DSCLR)-------------C
    !----------------------------------------------------------------------C
    DO i = 1, nsol
       dsclr(i) = sl(i)
    END DO
    DO l = 1, lmp1
       l0 = (l-1) * nsol
       DO i = 1, nsol
          dsclr(i) = dsclr(i) + aclr(l0+i)
       END DO
    END DO
  END SUBROUTINE clear




  !-----------------------------------------------------------------------
  ! Subroutine: SETSW
  !
  ! CALL TWO MAIN SUBROUTINES: CLEAR AND CLOUDY
  !   Calculates
  !     ozone and water vapor amounts, clear sky
  !     reflectivities, clear sky ozone reflection, clear sky ground
  !     reflection; finds cloud top; computes cloud optical depth,
  !     cloudy sky reflectivities, cloudy sky ozone reflection, cloudy
  !     sky ground reflection.
  !     NUMBER OF CLOUDY DLGP,
  !     NUMBER OF TOP CLOUD LAYER IN EACH DLGP,
  !     COMPRESS OUT CLOUDY DAYTIME LATITUDE GRID POINT (DLGP) VALUES
  !         reflectivities, clear sky ozone reflection, clear sky ground
  !         reflection; finds cloud top; computes cloud optical depth,
  !         cloudy sky reflectivities, cloudy sky ozone reflection, cloudy
  !         sky ground reflection.
  !---------------------------------------------------------------------C
  ! ACRONYMS:
  !   CDLGP...CLOUDY DAYTIME LATITUDE GRID POINTS
  !   DLGP....DAYTIME LATITUDE GRID POINTS
  !   LGP.....LATITUDE GRID POINTS
  !   NSOL....NUMBER OF DAYTIME LATITUDE GRID POINTS
  !   NCLD....NUMBER OF CLOUDY DLGP
  !-----------------------------------------------------------------------

  SUBROUTINE setsw(ncols ,kmax  ,tice  ,icld  ,tauc  , &
       scosz ,cmu   ,cosmag,dsclr ,rvbl  ,scosc ,cmuc  , &
       csmcld,dscld ,rvbc  ,rvdl  ,rnbl  ,rndl  ,agv   ,agn   , &
       rvdc  ,rnbc  ,rndc  ,agncd ,rsurfv,rsurfn,sl    ,sc    , &
       ta    ,wa    ,oa    ,pu    ,aclr  ,dp    ,css   ,acld  , &
       dpc   ,ccu   ,listim,bitd  ,sqrt3 ,gg    ,ggp   , &
       ggsq  ,athrd ,tthrd ,rcn1  ,rcn2  ,tcrit ,ecrit ,np    , &
       lmp1  ,nsol  ,nslmp1,nsolnp,ncld  ,ncldp1,nclmp1, &
       ncldnp,dooz  )  !hmjb

    !==========================================================================
    ! INPUT AND OUTPUT VARIABLES
    !
    !    ncols.....Number of grid points on a gaussian latitude circle
    !    kmax......Number of grid points at vertical
    !    tice......tice=273.16 zero grau absoluto
    !    icld......>>> icld=1     : old cloud emisivity setting
    !                               ccu = ccu*(1-exp(-0.05*dp))
    !                               css = css*(1-exp(-0.01*dp))
    !                                     for ice cloud t<253.0
    !                               css = css*(1-exp(-0.05*dp))
    !                                     for     cloud t>253.0
    !              >>> icld=2     : new cloud emisivity setting
    !                               ccu = 1.0-exp(-0.12*ccu*dp)
    !                               css = 0.0 for      t<-82.5c
    !                               css = 1-exp(-1.5e-6*(t-tcrit)**2*css*dp)
    !                                     for -82.5<t<-10.0
    !                               css = 1-exp(-5.2e-3*(t-273.)-0.06)*css*dp)!
    !                                     for -10.0<t< 0.0
    !                               css = 1-exp(-0.06*css*dp)
    !                                     for t> 0.0c
    !             >>> icld = 3   : ccm3 based cloud emisivity
    !    clwp
    !    fice......controle of change of fase of water
    !    rei.......Ice particle Effective Radius
    !    rel.......Liquid particle Effective Radius
    !    tauc......Shortwave cloud optical depth
    !              only works if usind icld=3 or 4 (i.e. arakawa or clirad with ccm3)
    !    tsea
    !    scosz.....scosz(i)   = s0     * cmu(i)
    !                         where s0  is constant solar
    !                               cmu is cosine of solar zenith angle
    !    cmu.......is cosine of solar zenith angle
    !    cosmag....Magnification factor in DLGP
    !    dsclr.....Absorption of clear atmosphere and ground clear
    !    rvbl......Visible beam clear  "Downward ground fluxes in DLGP"
    !    scosc.....scosz(i)   = s0     * cmu(i)
    !                         where s0  is constant solar
    !                               cmu is cosine of solar zenith angle
    !    cmuc......cmuc is cosine of solar zenith angle
    !    csmcld....csmcld = cosmag(i)  = 1224.0 * cmu(i) * cmu(i) + 1.0
    !    dscld.....Total absorption in atmosphere and ground
    !    rvbc......Visible beam cloudy flux
    !    rvdl......Visible diffuse clear
    !    rnbl......NearIR beam clear
    !    rndl......NearIR diffuse clear
    !    agv.......Ground visible diffuse albedo in DLGP
    !    agn.......Ground near IR diffuse albedo in DLGP
    !    rvdc......Visible diffuse cloudy
    !    rnbc......NearIR beam cloudy
    !    rndc......NearIR diffuse cloudy
    !    agncd.....Ground near IR diffuse albedo in CDLGP
    !    rsurfv....Ground visible beam albedo in DLGP
    !    rsurfn....Ground near IR beam albedo in DLGP
    !    sl........Net ground clear flux
    !    sc........Ground absorption
    !    ta........Layer Temperature in DLGP
    !    wa........Layer specific humidity in DLGP
    !    oa........Layer ozone mixing ratio in DLGP
    !    pu........pressure at botton of layer in DLGP
    !    aclr......Absorption in clear ATM
    !    dp........Pressure difference in DLGP
    !    css.......Large scale cloud anount in DLGP
    !    acld......heatinf rate (cloudy)
    !    dpc.......Pressure difference in DLGP
    !    ccu.......cumulus cloud amount in DLGP
    !    litx......Numbers of DLGP in all layers
    !    listim....1,2,3...imax*(kmax+1)
    !    bitd......(.true.) in CDLGP
    !    sqrt3.....Magification factor for diffuse reflected radiation
    !              sqrt3  = SQRT(3.0)
    !    gg........Asymmetry Factor      = 0.85
    !    ggp.......ggp    = gg  / (1.0 + gg) * 0.75
    !    ggsq......ggsq   = gg  * gg
    !    athrd.....constant athrd  = 1.0 / 3.0
    !    tthrd.....constant tthrd  = 2.0 / 3.0
    !    rcn1......constant  rcn1   = 1.0 / 6.55
    !    rcn2......constant  rcn2   = 1.0 / 4.47238
    !    tcrit.....constant tcrit  = tice - 82.5
    !    ecrit.....ecrit  = 0.0105125
    !    np........np     = (kmax+2)
    !    lmp1......lmp1   = (kmax+1)
    !    nsol......Number of solar latitude grid points (cosz>0.01)
    !    nsolp1....nsolp1=nsol+1
    !    nslmp1....nslmp1=nsol*(kmax +1 ), where nsol is number of solar
    !              latitude grid points (cosz>0.01)
    !    nsolnp....nsolnp=nsol*(kmax +2 ), where nsol is number of solar
    !              latitude grid points (cosz>0.01)
    !    ncld......Number of cloudy DLGP
    !    ncldp1....ncldp1=NCLD+1
    !    nclmp1....nclmp1=NCLD*(kmax+1)
    !    ncldnp....ncldnp=NCLD*(kmax+2)
    !    dooz......dooz   = (.NOT. noz) where : noz = .FALSE.
    !               (do ozone computation )
    !-----------------------------------------------------------------------

    INTEGER, INTENT(in   ) :: ncols
    INTEGER, INTENT(in   ) :: kmax
    REAL(KIND=r8),    INTENT(in   ) :: tice
    REAL(KIND=r8),    INTENT(inout) :: tauc((ncols*(kmax+1)))
    INTEGER, INTENT(in   ) :: icld


    REAL(KIND=r8),    INTENT(in   ) :: sqrt3
    REAL(KIND=r8),    INTENT(in   ) :: gg
    REAL(KIND=r8),    INTENT(in   ) :: ggp
    REAL(KIND=r8),    INTENT(in   ) :: ggsq
    REAL(KIND=r8),    INTENT(in   ) :: athrd
    REAL(KIND=r8),    INTENT(in   ) :: tthrd
    REAL(KIND=r8),    INTENT(in   ) :: rcn1
    REAL(KIND=r8),    INTENT(in   ) :: rcn2
    REAL(KIND=r8),    INTENT(in   ) :: tcrit
    REAL(KIND=r8),    INTENT(in   ) :: ecrit
    INTEGER, INTENT(in   ) :: np
    INTEGER, INTENT(in   ) :: lmp1
    INTEGER, INTENT(in   ) :: nsol
    INTEGER, INTENT(in   ) :: nslmp1
    INTEGER, INTENT(in   ) :: nsolnp
    INTEGER, INTENT(inout  ) :: ncld
    INTEGER, INTENT(inout  ) :: ncldp1
    INTEGER, INTENT(inout  ) :: nclmp1
    INTEGER, INTENT(inout  ) :: ncldnp
    LOGICAL, INTENT(in   ) :: dooz


    REAL(KIND=r8),    INTENT(in   ) :: scosz (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: cmu   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: cosmag(ncols)
    REAL(KIND=r8),    INTENT(inout) :: dsclr (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rvbl  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: scosc (ncols)
    REAL(KIND=r8),    INTENT(inout) :: cmuc  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: csmcld(ncols)
    REAL(KIND=r8),    INTENT(inout) :: dscld (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rvbc  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rvdl  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rnbl  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rndl  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: agv   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: agn   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rvdc  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rnbc  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rndc  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: agncd (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: rsurfv(ncols)
    REAL(KIND=r8),    INTENT(in   ) :: rsurfn(ncols)
    REAL(KIND=r8),    INTENT(inout) :: sl    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: sc    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ta    ((ncols*(kmax+1)))
    REAL(KIND=r8),    INTENT(inout) :: wa    ((ncols*(kmax+1)))
    REAL(KIND=r8),    INTENT(inout) :: oa    ((ncols*(kmax+1)))
    REAL(KIND=r8),    INTENT(inout) :: pu    ((ncols*(kmax+2)))
    REAL(KIND=r8),    INTENT(inout) :: aclr  ((ncols*(kmax+1)))
    REAL(KIND=r8),    INTENT(in   ) :: dp    ((ncols*(kmax+1)))
    REAL(KIND=r8),    INTENT(inout) :: css   ((ncols*(kmax+1)))
    REAL(KIND=r8),    INTENT(inout) :: acld  ((ncols*(kmax+1)))
    REAL(KIND=r8),    INTENT(inout) :: dpc   ((ncols*(kmax+1)))
    REAL(KIND=r8),    INTENT(inout) :: ccu   ((ncols*(kmax+1)))
    INTEGER, INTENT(in   ) :: listim((ncols*(kmax+2)))
    LOGICAL, INTENT(inout) :: bitd  ((ncols*(kmax+2)))

    !---- LOCAL VARIABLES --------------------------------------------------

    REAL(KIND=r8)    :: agvcd (ncols) ! AGVCD (NCLD)...GROUND VISIBLE DIFFUSE
    ! ALBEDO IN CDLGP (SETSW)
    REAL(KIND=r8)    :: rlo   (ncols) ! Clear sky and ground reflection in
    ! DLGP for ozone computation
    REAL(KIND=r8)    :: rlg   (ncols) ! Clear sky and ground reflection in
    ! DLGP for ground absorption computation
    REAL(KIND=r8)    :: rsvcd (ncols) ! Ground Viseible beam albedo
    REAL(KIND=r8)    :: rsncd (ncols) ! Ground near infrared beam albedo
    REAL(KIND=r8)    :: rco   (ncols) ! Cloudy an ground reflection for ozone
    ! absorption comp.
    REAL(KIND=r8)    :: rcg   (ncols) ! Cloudy and ground reflection for ground
    ! absoption computation in visible region of spectrum
    REAL(KIND=r8)    :: taut  (ncols) ! Downward visible beam
    REAL(KIND=r8)    :: rc2   (ncols) ! Cloudy and ground reflection for ground absoption
    ! computation in region of spectrum from 0.7 to 0.9 mcm
    REAL(KIND=r8)    :: tr1   (ncols) ! Extinction of beam radiation (clear)
    REAL(KIND=r8)    :: tr2   (ncols) ! Extinction of beam radiation (cloud)
    REAL(KIND=r8)    :: tr3   (ncols) ! Extinction of beam radiation (Cloud, Raley)
    REAL(KIND=r8)    :: e0    ((ncols*(kmax+1))) ! cloud optical depth
    REAL(KIND=r8)    :: ozale ((ncols*(kmax+2))) ! Ozone amount in column in DLGP
    REAL(KIND=r8)    :: ozcd  ((ncols*(kmax+2))) ! ozone amount in column in CDLGP
    REAL(KIND=r8)    :: swale ((ncols*(kmax+2))) ! water vapor amount in column in DLGP
    REAL(KIND=r8)    :: swil  ((ncols*(kmax+1))) ! water vapor amount in layer in DLGP
    REAL(KIND=r8)    :: swilc ((ncols*(kmax+1))) ! water vapor amount in layer in CDLGP
    INTEGER :: litd  ((ncols*(kmax+2))) ! numbers de CDLGP in all layers
    LOGICAL :: bitc  ((ncols*(kmax+2))) ! Working logical dimension


    INTEGER :: icl (ncols)
    INTEGER :: icc (ncols)
    REAL(KIND=r8)    :: dum (1)
    REAL(KIND=r8)    :: arg1(ncols)
    REAL(KIND=r8)    :: a1  (ncols)
    REAL(KIND=r8)    :: e1  (ncols)
    REAL(KIND=r8)    :: upim(ncols)
    REAL(KIND=r8)    :: dnim(ncols)
    REAL(KIND=r8)    :: rc1 (ncols)
    REAL(KIND=r8)    :: tlim(ncols)
    REAL(KIND=r8)    :: g   (ncols)
    REAL(KIND=r8)    :: taui(ncols)
    REAL(KIND=r8)    :: b1  (ncols)
    REAL(KIND=r8)    :: c1  (ncols)

    REAL(KIND=r8)                   :: expcut
    INTEGER                :: i
    INTEGER                :: l
    INTEGER                :: l0
    INTEGER                :: l1
    INTEGER                :: k

    !---- SUBROUTINE STARTS HERE -------------------------------------------

    expcut=- LOG(1.0e53_r8)

    !-----------------------------------------------------------------------
    !   ZERO OUT OZONE BELOW 500 MB
    !-----------------------------------------------------------------------
    IF (dooz) THEN
       !---------------------------------------------------------------------C
       !------SET BITC=.TRUE. IF LAYER PRESSURE PU>500 mb--------------------C
       !---------------------------------------------------------------------C
       DO i = 1, nslmp1
          bitc(i) = pu(nsol+i).gt.500.0e0_r8
       END DO
       !---------------------------------------------------------------------C
       !------SET OZONE MIXING RATIO OA=0 IF LAYER PRESSURE PU>500 mb--------C
       !---------------------------------------------------------------------C
       WHERE(bitc(1:nslmp1)) oa(1:nslmp1)=0.0e0_r8
       !-----------------------------------------------------------------------
       !      compute ozone and water vapor amounts
       !-----------------------------------------------------------------------
       DO i=1,nsol
          ozale(i) = 0.0_r8
       END DO
       DO i=1,nslmp1                     ! NSOL*(KMAX+1)
          oa(i) = 476.0_r8 * oa(i) * dp(i)  ! Ozone amount in layer
       END DO                            ! DP...Pressure difference
       DO l=2,np
          l0 = (l-2) * nsol
          l1 = (l-1) * nsol
          !cdir nodep
          DO i=1,nsol
             ozale(l1+i) = ozale(l0+i) + oa(l0+i)
          END DO
       END DO
    END IF
    !---------------------------------------------------------------------C
    !--COMPUTE WATER VAPOR AMOUNT IN cm FROM SPECIFIC HUMIDITY IN g/g-----C
    !--SWALE (NSOL(KMAX+2)------------------------------------------------C
    !------WATER VAPOR AMOUNT IN ATMOSPHERE COLUMN ABOVE THE LAYER--------C
    !---------------------------------------------------------------------C
    DO i=1,nsol
       swale(i) = 0.0_r8
    END DO
    DO i=1,nsolnp     ! NSOL*(KMAX+2)
       pu(i) = pu(i) * pu(i)
    END DO
    swil(1:nslmp1)=sqrt(ta(1:nslmp1))
    DO i=1,nslmp1
       swil(i) = 120.1612_r8    * swil(i)
       swil(i) = (pu(nsol+i) - pu(i)) &
            * wa(i)       / swil(i) ! WA is specific humidity
    END DO
    DO l=2,np  !KMAX+2
       l0 = (l-2) * nsol
       l1 = (l-1) * nsol
       !cdir nodep
       DO i=1,nsol
          swale(l1+i) = swale(l0+i) + swil(l0+i)
       END DO
    END DO
    !-------------------------------------------------------------------------------C
    !     compute clear sky reflectivities
    !     clear sky ozone reflection, tauray=0.85
    !-------------------------------------------------------------------------------C
    !  RLO (NSOL)..CLEAR SKY AND GROUND REFLECTION FOR OZONE ABSORPTION COMPUTATION
    !  TAURAY=0.85 OPTICAL THICKNESS OF MOLECULAR SCATTERING IN UV SPECTRUM REGION
    !-------------------------------------------------------------------------------C
    DO i=1,nsol
       e1(i) = -0.85_r8 * cosmag(i)  !Magnification factor
       e1(i) = EXP(e1(i))
    END DO
    DO i=1,nsol
       b1(i)   = 3.0_r8 * cmu(i)    !Cosine solar zenith angle
       c1(i)   = 2.0_r8 - b1(i)
       b1(i)   = 2.0_r8 + b1(i)
       ! Eq. (3.6) - Rad. Doc. by J.Chagas & T.Tarasova
       a1(i)   = (b1(i) + c1(i) * e1(i)) * rcn1  ! 1.0_r8/6.55_r8
       upim(i) = 1.0_r8 - a1(i)
       ! Eq. (3.7) - Rad. Doc. by J.Chagas & T.Tarasova
       dnim(i) = a1(i) - e1(i)
       ! Eq. (3.5) - Rad. Doc. by J.Chagas & T.Tarasova
       rlo(i)  = upim(i) + (e1(i) * rsurfv(i) & !Ground,visible,beam alb.
            + dnim(i) * agv(i)) &               !Ground,visible,diffuse albedo
            * 0.576004_r8 / (1.0_r8 - 0.423996_r8 * agv(i))
       ! Eq. (3.4) - Rad. Doc. by J.Chagas & T.Tarasova
       rlo(i)  = tthrd * rsurfv(i) + athrd  * rlo(i)   !ATHRD=1.0/3.0
       !TTHRD=2.0/3.0
       !-------------------------------------------------------------------------------C
       !     clear sky ground reflection, tauray=0.15746
       !-------------------------------------------------------------------------------C
       !  RLG(NSOL)..CLEAR SKY AND GROUND REFLECTION FOR GROUND ABSORPTION COMPUTATION
       !  TAURAY=0.15746..OPTICAL THICKNESS OF MOLECULAR SCATTERING IN VISIBLE SPECTRUM REGION
       !-------------------------------------------------------------------------------C
       e1(i)   = -0.15746_r8 * cosmag(i)
       e1(i)   = EXP(e1(i))
    END DO
    DO i=1,nsol
       tr1(i)  = e1(i)
       a1(i)   = (b1(i) + c1(i) * e1(i)) * rcn2  ! 1.0/4.47238
       upim(i) = 1.0_r8 - a1(i)
       dnim(i) = a1(i) - e1(i)
       rlg(i)  = upim(i) + (e1(i)  * rsurfv(i) &
            + dnim(i) * agv(i)) * 0.88_r8 &
            / (1.0_r8 - 0.12_r8  * agv(i))
    END DO
    !-----------------------------------------------------------------------
    ! CALL SUBROUTINE CLEAR TO CALCULATE CLEAR SKY FLUXES
    !   computes for cycles over five bands ozone heating, radiational
    !   downflux and ground absorption, and water vapor heating.
    !-----------------------------------------------------------------------
    CALL clear (ncols ,kmax  ,sqrt3 ,np    ,lmp1  ,nsol  ,nslmp1,nsolnp, &
         dooz  ,scosz ,cosmag,dsclr ,rvbl  ,rvdl  ,rnbl  ,rndl  , &
         agv   ,rsurfv,rsurfn,sl    ,rlo   ,rlg   ,tr1   ,e0    , &
         ozale ,swale ,aclr                                       )

    !-----------------------------------------------------------------------
    ! FIND CLOUD TOP
    !-----------------------------------------------------------------------
    ! SET CLOUD AMOUNT as Maximum cloudiness
    DO i = 1, nslmp1
       css(i)=MAX(ccu(i),css(i))
    END DO

    !-----------------------------------------------------------------------
    ! FIND NUMBER OF TOP CLOUD LAYER ICC(NSOL) IN EACH LGP
    ! After this, icc stores, for each DLGP the layer id of the
    ! highest layer with clouds. If no clouds, icc stores np=kmax+2
    !-----------------------------------------------------------------------
    DO i=1,nsol
       icc(i) = np     !KMAX+2
    END DO
    DO k=1,lmp1              ! loop over all layers
       l=np-k                ! layer id goes from bottom (kmax+1) to TOA (1)
       l0 = (l-1) * nsol + 1 ! first position (I=1) of this layer in css matrix
       !-----------------------------------------------------------------------
       ! SET LOGICAL BITC(NSOL*(KMAX+2)=.TRUE. IF CLOUD AMOUNT CSS>0
       ! in the l-th layer, i.e., for all DLGP in this layer
       !-----------------------------------------------------------------------
       dum(1)=0.0e0_r8
       DO i = 1, nsol
          bitc(i) = css(l0+i-1).gt.0.0e0_r8
       END DO

       ! if there is a cloudy in this layer, over some of the DLGP,
       ! then we save the layer id in vector icc
       WHERE (bitc(1:nsol)) icc(1:nsol) = l
    END DO

    !-----------------------------------------------------------------------
    ! SET BITD(NSOL)=.TRUE. IF CLOUDS ARE AT ANY LAYER IN LGP
    !   NCLD........NUMBER OF CLOUDY LGP
    !   ICC(NSOL)...NUMBER OF TOP CLOUD LAYER IN EACH LGP
    ! Since icc has the index of the top cloud layer for each DLGP, it
    ! is only necessary to check which of the indexes are different
    ! from kmax+2=np to know if there is any cloud at all over this DLGP
    !-----------------------------------------------------------------------
    bitd(1:nsol) = icc(1:nsol) /= np
    ncld=COUNT(bitd(1:nsol))

    ! If there are no clouds, then copy clear fluxes over cloudy ones
    ! and exit from subroutine
    IF (ncld == 0) THEN    ! Cloudy fluxes are equal clear sky fluxes
       DO i=1,nslmp1       ! NSOL*(KMAX+1)
          acld(i) = aclr(i)
       END DO
       DO i=1,nsol
          sc(i)    = sl(i)
          dscld(i) = dsclr(i)
          rndc(i)  = rndl(i)
          rvdc(i)  = rvdl(i)
          rnbc(i)  = rnbl(i)
          rvbc(i)  = rvbl(i)
       END DO
       RETURN
    END IF
    !-----------------------------------------------------------------------
    ! SET NEW PARAMETERS FROM NCLD
    !-----------------------------------------------------------------------
    ncldp1 = ncld + 1
    nclmp1 = ncld * lmp1   !NCLD *(KMAX+1)
    ncldnp = ncld * np     !NCLD *(KMAX+2)
    !-----------------------------------------------------------------------
    ! COMPLETE BITD(NSOL(KMAX+1))
    ! Repeat the values we assigned for the first row (1:nsol) over
    ! all the other rows (2:kmax+2)
    !-----------------------------------------------------------------------
    DO l=1,lmp1       ! KMAX+1
       l1 = l * nsol + 1
       bitd(l1:l1+nsol-1) = icc(1:nsol) /= np
    END DO
    !-----------------------------------------------------------------------
    ! compress out cloudy grid point values for ozone amt., water amt.scosz, etc.
    !   NCLDNP....NCLD(KMAX+2)
    !   NSOLNP....NSOL(KMAX+2)
    !   LITD(NCLD*(KMAX+2))....NUMBERS OF CLOUDY DLGP
    !-----------------------------------------------------------------------
    litd(1:COUNT(bitd(1:nsolnp))) = PACK(listim(1:nsolnp), bitd(1:nsolnp))
    DO i = 1, ncld
       IF (litd(i)<=nsol) THEN
          icl(i)=icc(litd(i))
          cmuc(i)=cmu(litd(i))
          agvcd(i)=agv(litd(i))
          agncd(i)=agn(litd(i))
          rsvcd(i)=rsurfv(litd(i))
          rsncd(i)=rsurfn(litd(i))
          scosc(i)=scosz(litd(i))
          csmcld(i)=cosmag(litd(i))
       END IF
    END DO
    DO i = 1, nclmp1
       IF (litd(i)<=nslmp1) &
            swilc(i)=swil(litd(i))
    END DO
    IF (dooz) THEN
       DO i = 1, ncldnp
          IF(litd(i).le.nsolnp) ozcd(i)=ozale(litd(i))
       END DO
       !-----------------------------------------------------------------------
       !      make ozone amount constant below cloud top
       !-----------------------------------------------------------------------
       DO l=1,lmp1
          bitc(1:ncld) = icl(1:ncld) <= l
          l0 = (l-1) * ncld
          l1 =  l    * ncld
          DO i = 1, ncld
             IF(bitc(i)) ozcd(l1+i)=ozcd(l0+i)
          END DO
       END DO
    END IF
    !-----------------------------------------------------------------------
    ! VALUES IN CLOUDY DAYTIME LATITUDE GRID POINTS (DLGP)
    !-----------------------------------------------------------------------
    ! ICL......NUMBER OF TOP CLOUD LAYER
    ! CMUC.....COSINE SOLAR ZENITH ANGLE
    ! AGVCD....GROUND VISIBLE DIFFUSE ALBEDO
    ! AGNCD....GROUND NEAR IR DIFFUSE ALBEDO
    ! RSVCD....GROUND VISIBLE BEAM ALBEDO
    ! RSNCD....GROUND NEAR IR BEAM ALBEDO
    ! SCOSC....SOLAR FLUX AT ATMOSPHERE TOP
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !     compute cloud optical depth
    !     icld = 1   : old cloud emisivity setting
    !     icld = 2   : new cloud emisivity setting
    !-----------------------------------------------------------------------
    IF (icld == 1) THEN
       DO i=1,nslmp1
          e0(i)   = 0.05_r8
          bitc(i) = (ta(i)  .LT. 253.0_r8).AND. (ccu(i) .EQ. 0.0_r8)
       END DO
       WHERE (bitc(1:nslmp1)) e0(1:nslmp1)=0.025_r8

    ELSE IF (icld == 2) THEN
       DO i=1,nslmp1
          e0(i) = (ta(i) - tcrit)
       END DO
       DO i = 1, nslmp1
          e0(i)=MAX(1.0e0_r8,e0(i))
       END DO

       DO i=1,nslmp1
          tauc(i) = 2.0e-6_r8 * e0(i) * e0(i)
       END DO
       DO i = 1, nslmp1
          tauc(i)=MIN(ecrit,tauc(i))
       END DO

       DO i=1,nslmp1
          e0(i) = 6.94875e-3_r8 * (ta(i) - tice) + 0.08_r8
       END DO

       DO i = 1, nslmp1
          e0(i)=min(0.08_r8,max(ecrit,e0(i)))
          bitc(i)=e0(i).eq.ecrit
       END DO

       WHERE (bitc(1:nslmp1)) e0(1:nslmp1)=tauc(1:nslmp1)
       bitc(1:nslmp1)=ccu(1:nslmp1)>0.0e0_r8
       WHERE(bitc(1:nslmp1)) e0(1:nslmp1)=0.16e0_r8

    ELSE IF ((icld == 3).OR.(icld == 4).OR.(icld == 5).OR.(icld == 6).OR.(icld == 7)) THEN
       DO i=1,nslmp1
          ! the extra cloud fraction idea from ncar
          e0(i) = SQRT(css(i))*css(i)*tauc(i)
       END DO
    END IF
    !-----------------------------------------------------------------------
    ! DPC...PRESSURE DIFFERENSE
    ! CSS...CLOUD AMOUNT
    ! E0....CLOUD OPTICAL DEPTH
    !-----------------------------------------------------------------------
    IF (icld /= 3) THEN
       DO i=1,nslmp1
          e0(i) = e0(i) * dpc(i) * css(i)
       END DO
    END IF
    !-----------------------------------------------------------------------
    ! TRANSFORM E0 TO TAUC IN CLOUDY DLGP
    !-----------------------------------------------------------------------

    DO i = 1, nclmp1
       IF (litd(i).le.nslmp1) tauc(i)=e0(litd(i))
    END DO

    DO i=1,ncld
       taut(i) = 0.0_r8
    END DO
    DO l=1,lmp1
       l0 = (l-1) * ncld
       DO i=1,ncld
          taut(i) = taut(i) + tauc(l0+i)
       END DO
    END DO
    !
    !     compute cloudy sky reflectivities
    !     cloudy sky ozone reflection, tauray=0.0_r8, gcld=0.85_r8
    !
    !-----CLOUDY SKY AND GROUND REFLECTION FOR GROUND ABSORPTION ---------C
    !-----------COMPUTATION IN SPECTRUM REGION FROM 0.7 TO 0.9 mcm--------C
    !-------------------------RC2(NCLD)-----------------------------------C
    !-----TAURAY=0.0....MOLECULAR OPTICAL DEPTH IN NEAR INFRARED REGION---C
    !-----GCLD=0.85...ASYMMETRY FACTOR OF CLOUD PARTICLES PHASE FUNCTION--C
    DO i=1,ncld
       tlim(i) = 1.0_r8 / (1.0_r8 + 0.1299_r8 * taut(i))
       arg1(i) = -0.2775_r8 * taut(i) * csmcld(i)  ! Magnification factor in CDLGP
       e1(i)   = EXP(arg1(i))
    END DO
    DO i=1,ncld
       tr2(i)  = e1(i)
       b1(i)   = 3.0_r8 * cmuc(i)
       c1(i)   = 2.0_r8 - b1(i)
       b1(i)   = 2.0_r8 + b1(i)
       a1(i)   = b1(i)   + c1(i)  * e1(i)
       a1(i)   = a1(i)   / (4.0_r8 + 0.45_r8 * taut(i))
       upim(i) = 1.0_r8          - a1(i)
       dnim(i) = a1(i)   - e1(i)
       rc1(i)  = upim(i) + (e1(i)    * rsvcd(i) &
            + dnim(i) * agvcd(i)) * tlim(i) &
            / (1.0_r8 - (1.0_r8  - tlim(i))  * agvcd(i))
       rc2(i)  = upim(i) + (e1(i)    * rsncd(i) &
            + dnim(i) * agncd(i)) * tlim(i) &
            / (1.0_r8 - (1.0_r8  - tlim(i))  * agncd(i))
       !
       !     cloudy sky ozone reflection, tauray=0.85, gcld=0.85
       !
       !--CLOUDY SKY AND GROUND REFLECTION FOR OZONE ABSORPTION COMPUTATION--C
       !-------------------RCO(NCLD)-----------------------------------------C
       !-----TAURAY=0.85..MOLECULAR OPTICAL DEPTH IN UV REGION OF SPECTRUM---C
       !-----GCLD=0.85....ASYMMETRY FACTOR OF CLOUD PARTICLE PHASE FUNCTION--C
       taui(i) = taut(i) + 0.85_r8
       g(i)    = 0.85_r8 * taut(i) / taui(i)
       tlim(i) = 1.0_r8  / (1.0_r8 + 0.866_r8 * (1.0_r8 - g(i)) &
            * taui(i))
       arg1(i) = - (1.0_r8 - g(i) * g(i)) &
            * taui(i)     * csmcld(i)
       e1(i)   = EXP(arg1(i))
    END DO
    DO i=1,ncld
       a1(i)   = b1(i) + c1(i) * e1(i)
       a1(i)   = a1(i) / (4.0_r8 + 3.0_r8 * (1.0_r8 - g(i)) &
            * taui(i))
       upim(i) = 1.0_r8          - a1(i)
       dnim(i) = a1(i)   - e1(i)
       rco(i)  = upim(i) + (e1(i)    * rsvcd(i) &
            + dnim(i) * agvcd(i)) * tlim(i) &
            / (1.0_r8 - (1.0_r8  - tlim(i))  * agvcd(i))
       rco(i)  = tthrd * rc1(i) + athrd * rco(i)
       !
       !     cloudy sky ground reflection, tauray=0.15746, gcld=0.85
       !
       !-----------CLOUDY SKY AND GROUND REFLECTION--------------------------C
       !-----------FOR GROUND ABSORPTION COMPUTATION IN VISIBLE SPECTRUM-----C
       !------------------------RCG(NCLD)------------------------------------C
       !--TAURAY=0.15746..MOLECULAR OPTICAL DEPTH IN VISIBLE SPECTRUM--------C
       !--GCLD=0.85..ASYMMETRY FACTOR OF CLOUD PARTICLE PHASE FUNCTION-------C
       taui(i) = taut(i)        + 0.15746_r8
       g(i)    = 0.85_r8 * taut(i) / taui(i)
       tlim(i) = 1.0_r8  / (1.0_r8 + 0.866_r8 * (1.0_r8 - g(i)) &
            * taui(i))
       arg1(i) = -(1.0_r8 - g(i) * g(i)) &
            * taui(i)     * csmcld(i)
       e1(i)   = EXP(arg1(i))
    END DO
    DO i=1,ncld
       tr3(i)  = e1(i)
       a1(i)   = b1(i) + c1(i) * e1(i)
       a1(i)   = a1(i) / (4.0_r8 + 3.0_r8 * (1.0_r8 - g(i)) &
            * taui(i))
       upim(i) = 1.0_r8 - a1(i)
       dnim(i) = a1(i)   - e1(i)
       rcg(i)  = upim(i) + (e1(i)    * rsvcd(i) &
            + dnim(i) * agvcd(i)) * tlim(i) &
            / (1.0_r8 - (1.0_r8  - tlim(i))  * agvcd(i))
    END DO

    !-----------------------------------------------------------------------
    ! CALL SUBROUTINE CLOUDY TO CALCULATE CLOUDY SKY FLUXES
    !   computes  ozone heating, downflux and ground absorption, water
    !   vapor heating and computes for cycles over five bands total
    !   optical depth, reflection and transmission for diffuse
    !   incidence,direct transmission, diffuse reflection and
    !   transmission for direct beam, upward and dwonward fluxes at
    !   layer boundaries, absorption in the column.
    !-----------------------------------------------------------------------

    CALL cloudy (ncols ,kmax  ,scosc ,cmuc  ,csmcld,dscld ,rvbc  ,rvdc  , &
         rnbc  ,rndc  ,agvcd ,agncd ,rsvcd ,rsncd ,sc    ,rco   , &
         rcg   ,taut  ,rc2   ,tr1   ,tr2   ,tr3   ,ta    ,wa    , &
         oa    ,e0    ,pu    ,ozcd  ,swale ,swil  ,css   ,acld  , &
         dpc   ,swilc ,ccu   ,tauc  ,litd  ,sqrt3 ,gg    ,ggp   , &
         ggsq  ,np    ,lmp1  ,nsol  ,ncld  ,nclmp1,ncldnp,dooz    )

  END SUBROUTINE setsw


  !-----------------------------------------------------------------------
  ! Subroutine: SWRAD
  !
  ! MAIN SUBROUTINE FOR SOLAR RADIATION COMPUTATIONS
  ! CALL ONE MAIN SUBROUTINE:  SETSW
  !   COMPRESSES TWO-SIZE INPUT ARRAYES
  !   TO ONE-SIZE ARRAYES IN DAYTIME LATITUDE GRID POINTS (DLGP)
  !   CALCULATE DIRECT SURFACE ALBEDO FROM DIFFUSE ONES,
  !   CALCULATE DIFFUSE MAGNIFICATION FACTOR (Rodgers, 1967)
  !   SET OUTPUT SHORTWAVE RADIATIVE FLUXES IN ALL LGP.
  !---------------------------------------------------------------------C
  ! ACRONYMS:
  !   CDLGP...CLOUDY DAYTIME LATITUDE GRID POINTS
  !   DLGP....DAYTIME LATITUDE GRID POINTS
  !   LGP.....LATITUDE GRID POINTS
  !   NSOL....NUMBER OF DAYTIME LATITUDE GRID POINTS
  !   NCLD....NUMBER OF CLOUDY DLGP
  !-----------------------------------------------------------------------
  SUBROUTINE swrad( &
       ! Model Info and flags
       ncols , kmax  , nls   , noz   , &
       icld  , inalb , s0    , cosz  , &
       ! Atmospheric fields
       pl20  , dpl   , tl    , ql    , &
       o3l   ,                         &
       ! SURFACE:  albedo
       alvdf , alndf , alvdr , alndr , &
       ! SW Radiation fields 
       swinc ,                         &
       radvbc, radvdc, radnbc, radndc, &
       radvbl, radvdl, radnbl, radndl, &
       dswclr, dswtop, ssclr , ss    , &
       aslclr, asl   ,                 &
       ! Cloud field
       cld   , clu   , taud               )
    !
    ! >>> inalb= 1    : input two  types surfc albedo (2 diffused)
    !                   direct beam albedos are calculated by the subr.
    ! >>> inalb= 2    : input four types surfc albedo (2 diff,2 direct)
    ! >>> icld = 1    : old cloud emisivity (optical depth) setting
    !             ccu :     0.05 *dp
    !             css :     0.025*dp             for ice cloud t<253.0
    !                       0.05 *dp             for ice cloud t>253.0
    ! >>> icld = 2    : new cloud emisivity (optical depth) setting
    !             ccu :     (0.16)*dp
    !             css :      0.0                         t<-82.5c
    !                       (2.0e-6*(t-tcrit)**2)*dp    -82.5<t<-10.0c
    !                       (6.949e-3*(t-273)+.08)*dp   -10.0<t< 0.0c
    !                       (0.08)*dp                   -10.0<t< 0.0c
    ! >>> icld = 3    : ccm3 based cloud emisivity
    !
    !==========================================================================
    !   imax......Number of grid points on a gaussian latitude circle
    !   kmax......Number of grid points at vertical
    !   nls.......number of layers in the stratosphere.
    !   noz.......Logical (true when no ozone computation)
    !   icld......Input two types of cloud emissivity
    !             =1 : old cloud emissivity setting
    !             =2 : new cloud emissivity setting
    !   s0........Solar constant  at proper sun-earth distance
    !   inalb.....Input two types of surface albedo
    !             =1 : two diffuse (beam albedos are calculated)
    !             =2 : two diffuse , and two beam albedos
    !                  avisd,anird, avisb,anirb
    !   alvdf.....visible diffuse surface albedo
    !   alndf.....near-ir diffuse surface albedo
    !   alvdr.....visible beam surface albedo
    !   alndr.....near-ir beam surface albedo
    !   cosz......Cosine of zenith angle
    !   pl20......Flip array of pressure at bottom of layers (mb)
    !             flip arrays (k=1 means top of atmosphere)
    !             pl20(i,k)=gps(i)*sigml(kflip) where
    !                         gps   =  surface pressure   (mb)
    !                         sigml =  sigma coordinate at bottom of layer
    !   dpl.......Flip array  of pressure difference bettween levels
    !             flip arrays (k=1 means top of atmosphere)
    !   tl........Flip array of temperature in kelvin
    !   ql........Flip array of specific humidity in g/g
    !   o3l.......Ozone mixing ratio (g/g) in 18 layers and in all latitude grids
    !   cld.......Large scale cloud amount
    !   clu.......cumulus cloud amount
    !   swinc.....Solar flux at top of atmosphere
    !   dswclr....Absorption in the clear atmosphere and at the ground
    !   dswtop....Absorption in the cloudy atmosphere and at the ground
    !   ssclr.....Absorption  at the ground in clear case
    !   ss........Absorption at the ground in cloudy case
    !   aslclr....Heating rate (clear case) (K/s)
    !   asl.......Heating rate (cloudy case) (K/s)
    !   radvbl....Downward Surface shortwave fluxe visible beam (clear)
    !   radvdl....Downward Surface shortwave fluxe visible diffuse (clear)
    !   radnbl....Downward Surface shortwave fluxe Near-IR beam (clear)
    !   radndl....Downward Surface shortwave fluxe Near-IR diffuse (clear)
    !   radvbc....Downward Surface shortwave fluxe visible beam (cloudy)
    !   radvdc....Downward Surface shortwave fluxe visible diffuse (cloudy)
    !   radnbc....Downward Surface shortwave fluxe Near-IR beam (cloudy)
    !   radndc....Downward Surface shortwave fluxe Near-IR diffuse (cloudy)
    !   taud......Shortwave cloud optical depth
    !             only works if usind icld=3 or 4 (i.e. arakawa or clirad with ccm3)
    !==========================================================================

    ! Model Info and flags
    INTEGER      ,    INTENT(in   ) :: ncols
    INTEGER      ,    INTENT(in   ) :: kmax
    INTEGER      ,    INTENT(in   ) :: nls
    LOGICAL      ,    INTENT(in   ) :: noz
    INTEGER      ,    INTENT(in   ) :: icld
    INTEGER      ,    INTENT(in   ) :: inalb
    REAL(KIND=r8),    INTENT(in   ) :: s0
    REAL(KIND=r8),    INTENT(in   ) :: cosz   (ncols)

    ! Atmospheric fields
    REAL(KIND=r8),    INTENT(in   ) :: pl20   (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: dpl    (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: tl     (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: ql     (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: o3l    (ncols,kmax)

    ! SURFACE:  albedo
    REAL(KIND=r8),    INTENT(in   ) :: alvdf  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: alndf  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: alvdr  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: alndr  (ncols)

    ! SW Radiation fields 
    REAL(KIND=r8),    INTENT(out) :: swinc  (ncols)
    REAL(KIND=r8),    INTENT(out) :: radvbc (ncols)
    REAL(KIND=r8),    INTENT(out) :: radvdc (ncols)
    REAL(KIND=r8),    INTENT(out) :: radnbc (ncols)
    REAL(KIND=r8),    INTENT(out) :: radndc (ncols)
    REAL(KIND=r8),    INTENT(out) :: radvbl (ncols)
    REAL(KIND=r8),    INTENT(out) :: radvdl (ncols)
    REAL(KIND=r8),    INTENT(out) :: radnbl (ncols)
    REAL(KIND=r8),    INTENT(out) :: radndl (ncols)
    REAL(KIND=r8),    INTENT(out) :: dswclr (ncols)
    REAL(KIND=r8),    INTENT(out) :: dswtop (ncols)
    REAL(KIND=r8),    INTENT(out) :: ssclr  (ncols)
    REAL(KIND=r8),    INTENT(out) :: ss     (ncols)
    REAL(KIND=r8),    INTENT(out) :: aslclr (ncols,kmax)
    REAL(KIND=r8),    INTENT(out) :: asl    (ncols,kmax)

    ! Cloud field and Microphysics
    REAL(KIND=r8),    INTENT(in   ) :: cld    (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: clu    (ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: taud   (ncols,kmax)

    !---- LOCAL VARIABLES --------------------------------------------------

    INTEGER :: im      ! =IMAX
    INTEGER :: np      ! =KMAX+2
    INTEGER :: npp1    ! =KMAX+3
    INTEGER :: lmp1    ! =KMAX+1
    INTEGER :: ls1
    INTEGER :: ls
    INTEGER :: nsol    ! NUMBER OF SOLAR LATITUDE GRID POINTS (COSZ>0.01)
    INTEGER :: nsolp1  ! NSOL+1
    INTEGER :: nslmp1  ! NSOL*(KMAX+1)
    INTEGER :: nsolnp  ! NSOL*(KMAX+2)
    INTEGER :: ncld    ! NUMBER OF CLOUDY DLGP
    INTEGER :: ncldp1  ! NCLD+1
    INTEGER :: nclmp1  ! NCLD*(KMAX+1)
    INTEGER :: ncldnp  ! NCLD*(KMAX+2)
    LOGICAL :: dooz    ! DO OZONE CALCULATION (.NOT.NOZ)

    REAL(KIND=r8)    :: scosz (ncols) ! DOWNWARD FLUX AT TOP IN DLGP
    REAL(KIND=r8)    :: cmu   (ncols) ! COSINE OF SOLAR ZENITH ANGLE IN DLGP
    REAL(KIND=r8)    :: cosmag(ncols) ! MAGNIFICATION FACTOR IN DLGP
    REAL(KIND=r8)    :: dsclr (ncols) ! ABSORPTION OF CLEAR ATMOSPHERE AND GROUND-CLEAR
    REAL(KIND=r8)    :: scosc (ncols) ! DOWNWARD FLUX AT TOP IN CLOUDY DLGP
    REAL(KIND=r8)    :: cmuc  (ncols) ! COSINE OF SOLAR ZENITH ANGLE IN CDLGP
    REAL(KIND=r8)    :: csmcld(ncols) ! MAGNIFICATION FACTOR IN CDLGP
    REAL(KIND=r8)    :: dscld (ncols)
    REAL(KIND=r8)    :: rvbl  (ncols) ! VISIBLE, BEAM, CLEAR                      (CLEAR )
    REAL(KIND=r8)    :: rvbc  (ncols) ! VISIBLE, BEAM, CLOUDY                     (CLOUDY)
    REAL(KIND=r8)    :: rvdl  (ncols) ! VISIBLE, DIFFUSE, CLEAR                   (CLEAR )
    REAL(KIND=r8)    :: rnbl  (ncols) ! NearIR, BEAM, CLEAR                       (CLEAR )
    REAL(KIND=r8)    :: rndl  (ncols) ! NearIR, DIFFUSE, CLEAR                    (CLEAR )
    REAL(KIND=r8)    :: agv   (ncols) ! GROUND VISIBLE DIFFUSE ALBEDO IN DLGP     (SWRAD )
    REAL(KIND=r8)    :: agn   (ncols) ! GROUND NEAR IR DIFFUSE ALBEDO IN DLGP     (SWRAD )
    REAL(KIND=r8)    :: rvdc  (ncols) ! VISIBLE, DIFFUSE, CLOUDY                  (CLOUDY)
    REAL(KIND=r8)    :: rnbc  (ncols) ! NearIR, BEAM, CLOUDY                      (CLOUDY)
    REAL(KIND=r8)    :: rndc  (ncols) ! NearIR, DIFFUSE, CLOUDY                   (CLOUDY)
    REAL(KIND=r8)    :: agncd (ncols) ! GROUND NEAR IR DIFFUSE ALBEDO IN CDLGP    (SETSW )
    REAL(KIND=r8)    :: rsurfv(ncols) ! GROUND VISIBLE BEAM ALBEDO IN DLGP        (SWRAD )
    REAL(KIND=r8)    :: rsurfn(ncols) ! GROUND NEAR IR BEAM ALBEDO IN DLGP        (SWRAD )
    REAL(KIND=r8)    :: sl    (ncols) ! NET GROUND CLEAR FLUX                     (CLEAR )
    REAL(KIND=r8)    :: sc    (ncols)
    REAL(KIND=r8)    :: ta    ((ncols*(kmax+1))) ! LAYER TEMPERATURE IN DLGP           (SWRAD)
    REAL(KIND=r8)    :: wa    ((ncols*(kmax+1))) ! LAYER SPECIFIC HUMIDITY IN DLGP     (SWRAD)
    REAL(KIND=r8)    :: oa    ((ncols*(kmax+1))) ! LAYER OZONE MIXING RATIO IN DLGP    (SWRAD)
    REAL(KIND=r8)    :: pu    ((ncols*(kmax+2))) ! PRESSURE AT BOTTOM OF LAYER IN DLGP (SWRAD)
    REAL(KIND=r8)    :: tauc  ((ncols*(kmax+1))) ! CLOUD OPTICAL DEPTH IN CLOUDY DLGP
    REAL(KIND=r8)    :: aclr  ((ncols*(kmax+1)))! ABSORPTION IN CLEAR ATM.             (CLEAR)
    REAL(KIND=r8)    :: dp    ((ncols*(kmax+1)))! PRESSURE DIFFERENCE IN DLGP          (SWRAD)
    REAL(KIND=r8)    :: css   ((ncols*(kmax+1)))! LARGE SCALE CLOUD AMOUNT IN DLGP     (SWRAD)
    REAL(KIND=r8)    :: acld  ((ncols*(kmax+1)))! HEATING RATE                         (CLOUDY)
    REAL(KIND=r8)    :: dpc   ((ncols*(kmax+1)))! PRESSURE DIFFERENCE IN DLGP          (SWRAD)
    REAL(KIND=r8)    :: ccu   ((ncols*(kmax+1)))! CUMULUS CLOUD AMOUNT IN DLGP         (SWRAD)
    INTEGER :: litx  ((ncols*(kmax+1)))! NUMBERS OF DLGP IN ALL LAYERS        (SWRAD)
    INTEGER :: listim((ncols*(kmax+2)))! =1,2,3...IMAX*(KMAX+1)               (SWRAD)
    LOGICAL :: bitx  ((ncols*(kmax+1)))! (.TRUE.) IN SOLAR LATITUDE GRID POINTS
    LOGICAL :: bitd  ((ncols*(kmax+2)))! (.TRUE.) IN CLOUDY DLGP
    LOGICAL :: bitn((ncols*(kmax+2)))

    !-------------------------------------------------------------------C
    !-----FAC1...COEFFICIENT FOR COMPUTING HEATING RATE-----------------C
    !-----DAY....NUMBER OF SECONDS IN DAY-------------------------------C
    !-------------------------------------------------------------------C
    REAL(KIND=r8), PARAMETER    :: fac1  = 8.441874377_r8
    REAL(KIND=r8), PARAMETER    :: day   = 86400.0_r8
    REAL(KIND=r8), PARAMETER    :: fac   = fac1 / day
    REAL(KIND=r8), PARAMETER    :: pai   = 3.141592653589793_r8
    REAL(KIND=r8), PARAMETER    :: tice  = 273.16_r8
    REAL(KIND=r8), PARAMETER    :: gg    = 0.85_r8
    REAL(KIND=r8), PARAMETER    :: ggp   = GG/(1.0_r8+GG)*0.75_r8
    REAL(KIND=r8), PARAMETER    :: ggsq  = GG*GG
    REAL(KIND=r8), PARAMETER    :: athrd = 1.0_r8/3.0_r8
    REAL(KIND=r8), PARAMETER    :: tthrd = 2.0_r8/3.0_r8
    REAL(KIND=r8), PARAMETER    :: rcn1  = 1.0_r8/6.55_r8
    REAL(KIND=r8), PARAMETER    :: rcn2  = 1.0_r8/4.47238_r8
    REAL(KIND=r8), PARAMETER    :: tcrit = TICE-82.5_r8
    REAL(KIND=r8), PARAMETER    :: ecrit = 0.0105125_r8
    REAL(KIND=r8), PARAMETER    :: p1em9 = 0.1e-9_r8       
    REAL(KIND=r8), PARAMETER    :: p1em22= 0.1e-22_r8      

    REAL(KIND=r8) :: sqrt3 
    REAL(KIND=r8) :: expcut
    INTEGER                :: ik
    INTEGER                :: lm
    INTEGER                :: i
    INTEGER                :: k
    INTEGER                :: il
    INTEGER                :: nsol2
    INTEGER                :: nsl2p1
    INTEGER                :: nsollm
    INTEGER                :: nzercd
    INTEGER                :: nlimwa
    INTEGER                :: nrstwa
    INTEGER                :: nwa1
    INTEGER                :: l
    INTEGER                :: l1


    !---- SUBROUTINE STARTS HERE -------------------------------------------
    dooz   = (.NOT. noz)
    lm     = kmax
    im     = ncols
    np     = (kmax+2)
    npp1   = np+1
    lmp1   = (kmax+1)
    ls1    = nls+1
    ls     = nls
    !-----------------------------------------------------------------------
    ! FAC...Coefficient for heating rate calculation in K/s
    !-----------------------------------------------------------------------
    sqrt3 = SQRT(3.0_r8)   
    expcut= -LOG(1.0e53_r8)
    acld  = 0.0_r8
    sc    = 0.0_r8
    !-----------------------------------------------------------------------
    ! SET ARRAY LISTIM = I, WHEN I=1,ncols*(kmax+2)
    !-----------------------------------------------------------------------
    DO i = 1, ncols*(kmax+2)
       listim(i)=i
    END DO

    !-----------------------------------------------------------------------
    ! set bits for daytime grid points
    !-----------------------------------------------------------------------
    bitx(1:ncols)=cosz(1:ncols).ge.0.01e0_r8
    !-----------------------------------------------------------------------
    ! CALCULATE NSOL = NUMBER OF DAYTIME LATITUDE GRID POINTS
    !-----------------------------------------------------------------------
    nsol=COUNT(bitx(1:ncols))
    !-----------------------------------------------------------------------
    ! SET ZERO TO ALL LATITUDE GRIDS SURFACE FLUXES
    !-----------------------------------------------------------------------
    DO i = 1, im
       swinc(i)  = 0.0_r8
       ss(i)     = 0.0_r8
       ssclr(i)  = 0.0_r8
       dswtop(i) = 0.0_r8
       dswclr(i) = 0.0_r8
       radvbl(i) = 0.0_r8
       radvdl(i) = 0.0_r8
       radnbl(i) = 0.0_r8
       radndl(i) = 0.0_r8
       radvbc(i) = 0.0_r8
       radvdc(i) = 0.0_r8
       radnbc(i) = 0.0_r8
       radndc(i) = 0.0_r8
    END DO
    DO k=1,kmax
       DO il = 1,ncols
          asl(il,k)    = 0.0_r8
          aslclr(il,k) = 0.0_r8
       END DO
    END DO

    !-----------------------------------------------------------------------
    ! IF THERE ARE NO DAYTIME POINTS THEN RETURN
    !-----------------------------------------------------------------------
    IF (nsol == 0) RETURN

    nsolp1     = nsol   + 1
    nsol2      = nsol   + nsol
    nsl2p1     = nsol2  + 1
    nslmp1     = nsol   * lmp1
    nsollm     = nsol   * lm
    nsolnp     = nsol   * np
    nzercd     = nsol   * ls1
    nlimwa     = nsol   * ls
    nrstwa     = nslmp1 - nlimwa
    nwa1       = nlimwa + 1
    !-----------------------------------------------------------------------
    ! SET BITX IN ALL GRID POINTS AT ALL LEVELS AS AT FIRST LEVEL
    !   BITX IS (.TRUE.) IN DAYTIME GRID POINTS
    !   SIZE OF BITX IS IMLMP1=NCOLS*(KMAX+1)
    !-----------------------------------------------------------------------
    DO l = 1, lm
       l1 = l * im
       DO i = 1, im
          bitx(l1+i)=bitx(i)
       END DO
    END DO
    !-----------------------------------------------------------------------
    ! SET PRESSURE AT FIRST AND SECOND LEVELS FOR RADIATION ONLY
    !-----------------------------------------------------------------------
    DO i = 1, nsol
       pu(i) = 0.0_r8
    END DO
    DO i = nsolp1, nsol2    !NSOL+1....2*NSOL
       pu(i) = 0.5_r8
    END DO
    !-----------------------------------------------------------------------
    ! SET INTEGER ARRAY LITX (NSOL*(KMAX+1))
    ! NUMBERS OF LATITUDE DAYTIME GRID POINTS AT ALL LEVELS, USING
    ! INTEGER ARRAY LISTIM (NCOLS*(KMAX+1)) = 1,2,....NCOLS*(KMAX+1)
    !-----------------------------------------------------------------------
    litx(1:COUNT(bitx(1:(ncols*(kmax+1))))) = &
         PACK(listim(1:(ncols*(kmax+1))), bitx(1:(ncols*(kmax+1))))
    !-----------------------------------------------------------------------
    !  TRANSFORM  TWO-SIZE  INPUT ARRAYS:
    !        PL20(NCOLS,KMAX),DPL,TL,QL,CLD,CLU
    !  IN ONE-SIZE ARRAYS:
    !        PU(NSOL*(KMAX+2)),DP(NSOL*(KMAX+1)),TA,WA,CSS,CCU
    !  IN DAYTIME LATITUDE GRID POINTS AT ALL LEVELS -> LITX(NSOL*KMAX)
    !-----------------------------------------------------------------------
!    DO i = 1, nsollm ! nsol   * lm
!       IF (litx(i).le.(ncols*kmax)) THEN
!          pu  (nsol2+i) =pl20 (litx(i),1)
!          dp  (nsol+i ) =dpl  (litx(i),1)
!          ta  (nsol+i ) =tl   (litx(i),1)
!          wa  (nsol+i ) =ql   (litx(i),1)
!          css (nsol+i ) =cld  (litx(i),1)
!          ccu (nsol+i ) =clu  (litx(i),1)
!          tauc(nsol+i ) =taud (litx(i),1)
!       END IF
!    END DO
    ik=0
    DO k=1,lm
       DO i = 1, im
          IF (bitx(i)) THEN
             ik=ik+1
             pu  (nsol2+ik ) =pl20 (i,k)
             dp  (nsol +ik ) =dpl  (i,k)
             ta  (nsol +ik ) =tl   (i,k)
             wa  (nsol +ik ) =ql   (i,k)
             css (nsol +ik ) =cld  (i,k)
             ccu (nsol +ik ) =clu  (i,k)
             tauc(nsol +ik ) =taud (i,k)
          END IF
       END DO
    END DO

    !-----------------------------------------------------------------------
    ! IF OZONE IS INCLUDED
    ! TRANSFORM  TWO-SIZE  INPUT ARRAY O3L(NCOLS,KMAX)
    ! IN ONE-SIZE ARRAY OA(NSOL*KMAX) IN DAYTIME LATITUDE GRID POINTS
    !-----------------------------------------------------------------------
    IF (dooz) THEN
!       DO i = 1, nsollm
!          IF (litx(i).le.(ncols*kmax)) oa(nsol+i)=o3l(litx(i),1)
!       END DO
       ik=0
       DO k=1,lm
          DO i = 1, im
             IF (bitx(i)) THEN
                ik=ik+1
	        oa(nsol+ik)=o3l(i,k)
             END IF
          END DO
       END DO

       DO i = 1, nsol
          oa(i) = oa(nsol+i)
       END DO
       DO i = 1, nslmp1
          oa(i)=MAX(p1em9,oa(i))
       END DO
    END IF
    !-----------------------------------------------------------------------
    ! THE SAME TRANSFORMATION AS MENTIONED ABOVE FOR:
    !   VISIBLE SURFACE ALBEDO....ALVDF to AGV
    !   NearIR SURFACE ALBEDO.....ALNDF to AGN
    !   COSINE OF SOLAR ZENITH ANGLE..COSZ to CMU
    !-----------------------------------------------------------------------
    DO i = 1, nsol
       IF (litx(i).le.im) THEN
          agv(i)=alvdf(litx(i))
          agn(i)=alndf(litx(i))
          cmu(i)=cosz(litx(i))
       END IF
    END DO
    !-----------------------------------------------------------------------
    ! IF DIRECT BEAM ALBEDOS ARE GIVEN THEN
    ! ALVDR TRANSFORM TO RSURFV(NSOL) AND ALNDR TO RSURFN(NSOL)
    ! IN DAYTIME GRID POINTS
    !-----------------------------------------------------------------------
    IF (inalb == 2) THEN
       DO i = 1, nsol
          IF (litx(i).le.im) THEN
             rsurfv(i)=alvdr(litx(i))
             rsurfn(i)=alndr(litx(i))
          END IF
       END DO
    ELSE
       !-----------------------------------------------------------------------
       ! IF DIRECT BEAM ALBEDOS ARE NOT GIVEN THEN DO THE REVERSE
       ! CALCULATE DIRECT BEAM SURFACE ALBEDO
       !-----------------------------------------------------------------------
       rvbl(1:nsol)=acos(cmu(1:nsol)) ! RVBL... solar zenith angle
       DO i = 1, nsol
          rvdc(i)  =  -18.0_r8 * (0.5_r8 * pai - rvbl(i)) / pai
          rvbc(i)  =  EXP(rvdc(i))
       END DO
       DO i = 1, nsol
          rvdc(i)  = (agv(i) - 0.054313_r8) / 0.945687_r8
          rndc(i)  = (agn(i) - 0.054313_r8) / 0.945687_r8
          rsurfv(i) = rvdc(i)+(1.0_r8-rvdc(i))*rvbc(i)
          rsurfn(i) = rndc(i)+(1.0_r8-rndc(i))*rvbc(i)
       END DO
       DO i = 1, im
          alvdr(i) = 0.0_r8
          alndr(i) = 0.0_r8
       END DO
       DO i = 1, nsol
          alvdr(litx(i))=rsurfv(i)
          alndr(litx(i))=rsurfn(i)
       END DO
    END IF
    !-----------------------------------------------------------------------
    ! SET SOME PARAMETERS AT FIRST, SECOND AND THIRD LEVELS
    !-----------------------------------------------------------------------
    DO i = 1, nsol
       dp(i) = pu(nsol+i)    ! pressure differense
       ta(i) = ta(nsol+i)    ! temperature
       wa(i) = wa(nsol+i)    ! specific humidity
       tauc(i) = 0.0_r8      ! cloud optical depth
    END DO
    DO i = 1, nslmp1         ! NSOL*(KMAX+1)
       dpc(i) = dp(i)        ! DPC...pressure differense
    END DO
    DO i = 1, nzercd         ! NSOL*3
       css(i) = 0.0_r8
       ccu(i) = 0.0_r8
    END DO
    DO i = 1, nlimwa         ! NSOL*2
       wa(i) = 3.0e-6_r8
    END DO
    DO i = nwa1,(nwa1+nrstwa-1)
       wa(i)=MAX(p1em22,wa(i)) ! nrstwa=nsol*(kmax+1)-nsol*2
    END DO

    !-----------------------------------------------------------------------
    ! CALCULATION OF MAGNIFICATION FACTOR (Rodgers, 1967)
    ! CMU.......COSINE OF SOLAR ZENITH ANGLE AT DLGP
    !-----------------------------------------------------------------------
    DO i = 1, nsol
       cosmag(i)  = 1224.0_r8 * cmu(i) * cmu(i) + 1.0_r8
    END DO
    cosmag(1:nsol) = sqrt(cosmag(1:nsol))
    DO i = 1, nsol
       cosmag(i)  = 35.0_r8   / cosmag(i)
       scosz(i)   = s0     * cmu(i)  ! DOWNWARD SOLAR FLUX AT TOP
    END DO
    !-----------------------------------------------------------------------
    ! TRANSFORM SCOSZ(NSOL) TO SWINC(NCOLS) AT ALL LGP
    !-----------------------------------------------------------------------
    DO i = 1,nsol
       swinc(litx(i))=scosz(i)
    END DO
    !-----------------------------------------------------------------------
    ! CALL SUBROUTINE SETSW TO CALCULATE SOLAR FLUXES
    ! SETSW CALLS CLEAR AND CLOUDY TO CALCULATE:
    !   CLEAR:  dsclr,sl,aclr,rvbl,rvdl,rnbl,rndl
    !   CLOUDY: dscld,sc,acld,rvbc,rvdc,rnbc,rndc
    ! The values are packed at the begining of the arrays.
    ! Instead of occupying 1..ncols, they cover only the range 1..nsol
    !-----------------------------------------------------------------------
    CALL setsw(ncols ,kmax  ,tice  ,icld  ,tauc  , &
         scosz , cmu  ,cosmag,dsclr ,rvbl  ,scosc ,cmuc  , &
         csmcld,dscld ,rvbc  ,rvdl  ,rnbl  ,rndl  ,agv   ,agn   , &
         rvdc  ,rnbc  ,rndc  ,agncd ,rsurfv,rsurfn,sl    ,sc    , &
         ta    ,wa    ,oa    ,pu    ,aclr  ,dp    ,css   ,acld  , &
         dpc   ,ccu   ,listim,bitd  ,sqrt3 ,gg    ,ggp   , &
         ggsq  ,athrd ,tthrd ,rcn1  ,rcn2  ,tcrit ,ecrit ,np    , &
         lmp1  ,nsol  ,nslmp1,nsolnp,ncld  ,ncldp1,nclmp1, &
         ncldnp,dooz  )
    ik=0
    DO k=1,lm
       DO i = 1, im
          IF (bitx(i)) THEN
             ik=ik+1
              taud (i,k)=tauc(nsol +ik )
          END IF
       END DO
    END DO

    !-----------------------------------------------------------------------
    ! SET SOLAR FLUXES AT ALL LATITUDE GRID POINTS
    ! For CLEAR it is already set.
    ! If there are any clouds, then copy the clear values over the
    ! zero's in the cloudy vectors, to complete the cloudy fields.
    !   BITD(nslmp1)=.TRUE. IF CLOUD IN LAYER OVER A DLGP
    !   BITN(nslmp1)=.TRUE. IF CLEAR SKY IN LAYER OVER A DLGP
    !-----------------------------------------------------------------------

    !hmjb??? What happen if ncld=0? I think that cloudy-arrays will
    ! have null values. Maybe clear->cloudy should be done always.

    IF (ncld /= 0) THEN
       bitn(1:nslmp1)=.not.bitd(1:nslmp1)
       DO i = 1,nslmp1
          IF (bitn(i)) acld(i)=aclr(i)
       END DO
       DO i = 1,nsol
          IF (bitn(i)) THEN
             rvbc(i)=rvbl(i)
             rvdc(i)=rvdl(i)
             rnbc(i)=rnbl(i)
             rndc(i)=rndl(i)
             sc(i)=sl(i)
             dscld(i)=dsclr(i)
          END IF
       END DO
    END IF
    !-----------------------------------------------------------------------
    ! SET SOLAR FLUXES IN ALL GRID POINTS
    ! All values are nsol-packed and need to be unpacked
    ! This is done by copying values from positions (1:nsol) to
    ! positions litx(1:nsol).
    !-----------------------------------------------------------------------
    DO i = 1,nsol
       ! clear
       ssclr(litx(i))=sl(i)
       dswclr(litx(i))=dsclr(i)
       radvbl(litx(i))=rvbl(i)
       radvdl(litx(i))=rvdl(i)
       radnbl(litx(i))=rnbl(i)
       radndl(litx(i))=rndl(i)

       ! cloudy
       ss(litx(i))=sc(i)
       dswtop(litx(i))=dscld(i)
       radvbc(litx(i))=rvbc(i)
       radvdc(litx(i))=rvdc(i)
       radnbc(litx(i))=rnbc(i)
       radndc(litx(i))=rndc(i)
    END DO
    ik=0
    DO k=1,lm
       DO i = 1, im
          IF (bitx(i)) THEN
             ik=ik+1
             aslclr(i,k) =aclr(nsol+ik)
             asl   (i,k) =acld(nsol+ik)
          END IF
       END DO
    END DO


    !-----------------------------------------------------------------------
    ! CALCULATION OF SOLAR HEATING RATE IN K/s
    !-----------------------------------------------------------------------
    DO k=1,kmax
       DO i = 1, ncols
          IF (aslclr(i,k) < 1.e-22_r8) aslclr(i,k)=0.0_r8
          aslclr(i,k) = aslclr(i,k) * fac / dpl(i,k)
          IF (asl(i,k) < 1.e-22_r8) asl(i,k) = 0.0_r8
          asl(i,k)    = asl(i,k)    * fac / dpl(i,k)
       END DO
    END DO

  END SUBROUTINE swrad



END MODULE Rad_COLA
