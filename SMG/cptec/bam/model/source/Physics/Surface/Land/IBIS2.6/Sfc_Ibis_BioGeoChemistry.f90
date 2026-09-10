MODULE Sfc_Ibis_Vegetation

  IMPLICIT NONE
SAVE

  PRIVATE
  INTEGER, PARAMETER :: r4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers 
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers

  INTEGER, PUBLIC, PARAMETER :: nVegClass=15

  PUBLIC :: sumnow
  PUBLIC :: sumday
  PUBLIC :: summonth
  PUBLIC :: sumyear
  PUBLIC :: gdiag
  PUBLIC :: vdiag
  PUBLIC :: climanl2
  PUBLIC :: soilbgc

  PUBLIC :: pheno
  PUBLIC :: dynaveg1 
  PUBLIC :: dynaveg2
  PUBLIC :: DailyDynaVeg
CONTAINS
  !
  ! #    #  ######   ####   ######   #####    ##     #####     #     ####   #    #
  ! #    #  #       #    #  #          #     #  #      #       #    #    #  ##   #
  ! #    #  #####   #       #####      #    #    #     #       #    #    #  # #  #
  ! #    #  #       #  ###  #          #    ######     #       #    #    #  #  # #
  !  #  #   #       #    #  #          #    #    #     #       #    #    #  #   ##
  !   ##    ######   ####   ######     #    #    #     #       #     ####   #    #
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE pheno(tc      , &! INTENT(IN   )
       agddu   , &! INTENT(INOUT) global
       tempu   , &! INTENT(INOUT) global
       agddl   , &! INTENT(INOUT) global
       templ   , &! INTENT(INOUT) global
       dropu   , &! INTENT(INOUT) global
       dropls  , &! INTENT(INOUT) global
       dropl4  , &! INTENT(INOUT) global
       dropl3  , &! INTENT(INOUT) global
       vegtype0, &! INTENT(IN   ) global
       froot   , &! INTENT(INOUT) global
       hsoi    , &! INTENT(IN   )
       beta1   , &! INTENT(IN   )
       beta2   , &! INTENT(IN   )
       plai    , &! INTENT(IN   )
       adplai  , &! INTENT(IN   )
       frac    , &! INTENT(OUT  )
       lai     , &! INTENT(OUT  )
       fl      , &! INTENT(IN   )
       fu      , &! INTENT(IN   )
       zbot    , &! INTENT(INOUT  )
       ztop    , &! INTENT(INOUT  )
       a10td   , &! INTENT(IN   )
       a10ancub, &! INTENT(IN   )
       a10ancls, &! INTENT(IN   )
       a10ancl4, &! INTENT(IN   )
       a10ancl3, &! INTENT(IN   )
       td      , &! INTENT(IN   )
       tthreshold , &! INTENT(IN   )
       gthreshold , &! INTENT(IN   )
       avglaiu    , &! INTENT(IN   )
       avglail    , &! INTENT(IN   )
       adnpp     , &! INTENT(IN   )
       adtsoi    , &! INTENT(IN   )
       adwsoi    , &! INTENT(IN   )
       adwisoi   , &! INTENT(IN   )
       poros     , &! INTENT(IN   )
       rhow      , &! INTENT(IN   )
       npoi    , &! INTENT(IN   )
       npft    , &! INTENT(IN   )
       nsoilay    , &! INTENT(IN   )
       nVegClass, &! INTENT(IN   )
       rootmode, &! INTENT(IN   )
       epsilon   )! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    ! common blocks
    !
    IMPLICIT NONE
    !
    INTEGER      , INTENT(IN   ) :: nVegClass! INTENT(IN   )
    CHARACTER(len=*) , INTENT(IN   ) ::  rootmode
    INTEGER      , INTENT(IN   ) :: npoi                ! total number of land points
    INTEGER      , INTENT(IN   ) :: npft                ! number of plant functional types
    INTEGER      , INTENT(IN   ) :: nsoilay!global   ! number of soil layers
    REAL(KIND=r8), INTENT(IN   ) :: epsilon             ! small quantity to avoid zero-divides and other
    ! truncation or machine-limit troubles with small
    ! values. should be slightly greater than o(1)
    ! machine precision
    REAL(KIND=r8), INTENT(IN   ) :: td       (npoi)  ! daily average temperature (K)
    REAL(KIND=r8), INTENT(IN   ) :: a10td    (npoi)  ! 10-day average daily air temperature (K)
    REAL(KIND=r8), INTENT(IN   ) :: a10ancub (npoi)  ! 10-day average canopy photosynthesis rate - broadleaf (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(IN   ) :: a10ancls (npoi)  ! 10-day average canopy photosynthesis rate - shrubs (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(IN   ) :: a10ancl4 (npoi)  ! 10-day average canopy photosynthesis rate - c4 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(IN   ) :: a10ancl3 (npoi)  ! 10-day average canopy photosynthesis rate - c3 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(IN   ) :: vegtype0 (npoi)      ! annual vegetation type - ibis classification
    REAL(KIND=r8), INTENT(INOUT) :: froot    (npoi,nsoilay,2)! global! fraction of root in soil layer 
    REAL(KIND=r8), INTENT(IN   ) :: hsoi     (npoi,nsoilay+1)   ! global ! soil layer thickness (m)
    REAL(KIND=r8), INTENT(IN   ) :: beta1(nVegClass)
    REAL(KIND=r8), INTENT(IN   ) :: beta2(nVegClass)
    REAL(KIND=r8), INTENT(IN   ) :: tc       (npoi)  ! coldest monthly temperature (C)
    REAL(KIND=r8), INTENT(INOUT) :: agddu    (npoi)  ! annual accumulated growing degree days for bud
    ! burst, upper canopy (day-degrees)
    REAL(KIND=r8), INTENT(INOUT) :: tempu    (npoi)  ! cold-phenology trigger for trees (non-dimensional)
    REAL(KIND=r8), INTENT(INOUT) :: agddl    (npoi)  ! annual accumulated growing degree days for bud burst,
    ! lower canopy (day-degrees)
    REAL(KIND=r8), INTENT(INOUT) :: templ    (npoi)  ! cold-phenology trigger for grasses/shrubs (non-dimensional)
    REAL(KIND=r8), INTENT(INOUT) :: dropu    (npoi)  ! drought-phenology trigger for trees (non-dimensional)
    REAL(KIND=r8), INTENT(INOUT) :: dropls   (npoi)  ! drought-phenology trigger for shrubs (non-dimensional)
    REAL(KIND=r8), INTENT(INOUT) :: dropl4   (npoi)  ! drought-phenology trigger for c4 grasses (non-dimensional)
    REAL(KIND=r8), INTENT(INOUT) :: dropl3   (npoi)  ! drought-phenology trigger for c3 grasses (non-dimensional)
    REAL(KIND=r8), INTENT(IN   ) :: plai     (npoi,npft)! total leaf area index of each plant functional type
    REAL(KIND=r8), INTENT(IN   ) :: adplai   (npoi,npft)! global  ! total leaf area index of each plant functional type (non-dimensional)
    REAL(KIND=r8), INTENT(OUT  ) :: frac     (npoi,npft)! fraction of canopy occupied by each plant functional type
    REAL(KIND=r8), INTENT(OUT  ) :: lai      (npoi,2)  ! canopy single-sided leaf area index (area leaf/area veg)
    REAL(KIND=r8), INTENT(IN   ) :: fl       (npoi)  ! fraction of snow-free area covered by lower  canopy
    REAL(KIND=r8), INTENT(IN   ) :: fu       (npoi)  ! fraction of overall area covered by upper canopy
    REAL(KIND=r8), INTENT(INOUT  ) :: zbot     (npoi,2)  ! height of lowest branches above ground (m)
    REAL(KIND=r8), INTENT(INOUT  ) :: ztop     (npoi,2)  ! height of plant top above ground (m)   
    REAL(KIND=r8), INTENT(OUT  ) :: tthreshold (npoi)   ! temperature threshold for budburst and senescence
    REAL(KIND=r8), INTENT(OUT  ) :: gthreshold (npoi)   ! temperature threshold for budburst and senescence
    REAL(KIND=r8), INTENT(OUT  ) :: avglaiu (npoi)  ! average lai of upper canopy 
    REAL(KIND=r8), INTENT(OUT  ) :: avglail (npoi)! average lai of lower canopy 
    REAL(KIND=r8), INTENT(IN   ) :: adnpp     (npoi,npft)! local ! instantaneous NPP for each pft (mol-CO2 / m-2 / second)
    REAL(KIND=r8), INTENT(IN   ) :: adtsoi    (npoi)! global -daily average ! soil temperature for each layer (C)    ! averages for first 2 layers (0.3 meters) of soil  	 
    REAL(KIND=r8), INTENT(IN   ) :: adwsoi    (npoi)! global -daily average ! fraction of soil pore space containing liquid water ! averages for first 2 layers (0.3 meters) of soil(fraction) 
    REAL(KIND=r8), INTENT(IN   ) :: adwisoi   (npoi)! global -daily average ! fraction of soil pore space containing ice    ! averages for first 2 layers (0.3 meters) of soil  	 
    REAL(KIND=r8), INTENT(IN   ) :: rhow           !global  ! density of liquid water (all types) (kg m-3)
    REAL(KIND=r8), INTENT(IN   ) :: poros   (npoi,nsoilay)! global ! porosity (mass of h2o per unit vol at sat / rhow)

 
 
 
    ! REAL(KIND=r8)    :: totdepth(npoi)
    !REAL(KIND=r8)    :: frootnorm1(npoi)
    !REAL(KIND=r8)    :: frootnorm2(npoi)
    !REAL(KIND=r8)    :: depth(npoi,nsoilay)   ! soil layer depth (cm)
    !
    ! local variables
    !
    INTEGER :: i      
    !
    ! INTEGER :: imonth     !
    ! INTEGER :: iday       !
    !
    REAL(KIND=r8) , PARAMETER   :: uplimit=1.0_r8      !
    REAL(KIND=r8) , PARAMETER   :: dwlimit=1.0_r8      !
    REAL(KIND=r8)    :: ddays      !
    REAL(KIND=r8)    :: ddfac      !
    !**** DTP 2000/06/28 Modified this following discussion with Navin. We now
    !*    retain fu(i) derived from dynaveg and constrain it to a local value
    !*    "fu_phys" in the range 0.25 to 0.975 used in the canopy physics calcs.
    REAL(KIND=r8)    :: fu_phys    ! Local value of fu(i) constrained to range 0.25 to 0.975
    ! to keep physics calculations stable.
    !
    ! define 'drop days' -- number of days to affect phenology change
    !
    ddays = 15.00_r8
    ddfac = 1.00_r8 / ddays
    !
    ! begin global grid
    !
    DO i = 1, npoi
       !
       ! ---------------------------------------------------------------------
       ! * * * upper canopy winter phenology * * *
       ! ---------------------------------------------------------------------
       !
       ! temperature threshold for budburst and senescence
       !
       ! temperature threshold is assumed to be 0 degrees C 
       ! or 5 degrees warmer than the coldest monthly temperature
       !
        tthreshold(i) = max (0.00_r8         + 273.160_r8,  &
                             tc(i) + 5.00_r8 + 273.160_r8)
       !
       ! gdd threshold temperature for leaf budburst
       ! with a growing degree threshold of 100 units
       !
        gthreshold(i) = 0.00_r8 + 273.160_r8
       ! 
       ! determine if growing degree days are initiated
       !
        IF (a10td(i).lt.gthreshold(i)) THEN
          agddu(i)  = 0.00_r8
        ELSE
          agddu(i) = agddu(i) + td(i) - gthreshold(i)
        END IF
       !
       ! determine leaf display
       !
        IF (a10td(i).lt.tthreshold(i)) THEN
          tempu(i)  = max (0.00_r8, tempu(i) - ddfac)
        ELSE
          tempu(i) = min (1.0_r8, max (0.00_r8, agddu(i) - 100.00_r8) / 50.00_r8)
        END IF
       !
       ! ---------------------------------------------------------------------
       ! * * * lower canopy winter phenology * * *
       ! ---------------------------------------------------------------------
       !
       ! temperature threshold for budburst and senescence
       !
       ! temperature threshold is assumed to be 0 degrees C 
       !
        tthreshold(i) = 0.00_r8 + 273.160_r8
       !
       ! gdd threshold temperature for leaf budburst
       ! with a growing degree threshold of 150 units
       !
        gthreshold(i) = -5.00_r8 + 273.160_r8
       ! 
       ! determine if growing degree days are initiated
       !
        IF (a10td(i).lt.gthreshold(i)) THEN
          agddl(i)  = 0.00_r8
       ELSE
          agddl(i) = agddl(i) + td(i) - gthreshold(i)
       END IF
       !
       ! determine leaf display
       !
        IF (a10td(i).lt.tthreshold(i)) THEN
          templ(i)  = max (0.00_r8, templ(i) - ddfac)
        ELSE
          templ(i) = min (1.0_r8, max (0.00_r8, agddl(i) - 150.00_r8) / 50.00_r8)
        END IF
       !
       ! ---------------------------------------------------------------------
       ! * * * drought canopy winter phenology * * *
       ! ---------------------------------------------------------------------
       !
       IF (a10ancub(i).lt.0.00_r8) dropu(i) = max (0.10_r8, dropu(i) - ddfac)
       IF (a10ancub(i).ge.0.00_r8) dropu(i) = min (1.00_r8, dropu(i) + ddfac)
       !
       IF (a10ancls(i).lt.0.00_r8) dropls(i) = max (0.10_r8, dropls(i) - ddfac)
       IF (a10ancls(i).ge.0.00_r8) dropls(i) = min (1.00_r8, dropls(i) + ddfac)
       !
       IF (a10ancl4(i).lt.0.00_r8) dropl4(i) = max (0.10_r8, dropl4(i) - ddfac)
       IF (a10ancl4(i).ge.0.00_r8) dropl4(i) = min (1.00_r8, dropl4(i) + ddfac)
       !
       IF (a10ancl3(i).lt.0.00_r8) dropl3(i) = max (0.10_r8, dropl3(i) - ddfac)
       IF (a10ancl3(i).ge.0.00_r8) dropl3(i) = min (1.00_r8, dropl3(i) + ddfac)
       !
       ! ---------------------------------------------------------------------
       ! * * * update lai and canopy fractions * * *
       ! ---------------------------------------------------------------------
       !
       ! upper canopy single sided leaf area index (area-weighted)
       !
       ! ---------------------------------------------------
       !
       ! these classes consist of some combination of 
       ! plant functional types:
       !
       ! ---------------------------------------------------
       !  1: tropical broadleaf evergreen trees
       !  2: tropical broadleaf drought-deciduous trees
       !  3: warm-temperate broadleaf evergreen trees
       !  4: temperate conifer evergreen trees
       !  5: temperate broadleaf cold-deciduous trees
       !  6: boreal conifer evergreen trees
       !  7: boreal broadleaf cold-deciduous trees
       !  8: boreal conifer cold-deciduous trees
       !  9: evergreen shrubs
       ! 10: cold-deciduous shrubs
       ! 11: warm (c4) grasses
       ! 12: cool (c3) grasses
       ! ---------------------------------------------------
       !        avglaiu = plai(i,1)             +  &
       !                  plai(i,2) * dropu(i)  +  &
       !                  plai(i,3)             +  &
       !                  plai(i,4)             +  &
       !                  plai(i,5) * tempu(i)  +  &
       !                  plai(i,6)             +  &
       !                  plai(i,7) * tempu(i)  +  &
       !                  plai(i,8) * tempu(i)

        avglaiu(i) = MIN(MAX(adplai(i,1),dwlimit*plai(i,1)),uplimit*plai(i,1))	          +  &
                     MIN(MAX(adplai(i,2),dwlimit*plai(i,2)),uplimit*plai(i,2)) * dropu(i)  +  &
                     MIN(MAX(adplai(i,3),dwlimit*plai(i,3)),uplimit*plai(i,3))	          +  &
                     MIN(MAX(adplai(i,4),dwlimit*plai(i,4)),uplimit*plai(i,4))	          +  &
                     MIN(MAX(adplai(i,5),dwlimit*plai(i,5)),uplimit*plai(i,5)) * tempu(i)  +  &
                     MIN(MAX(adplai(i,6),dwlimit*plai(i,6)),uplimit*plai(i,6))	          +  &
                     MIN(MAX(adplai(i,7),dwlimit*plai(i,7)),uplimit*plai(i,7)) * tempu(i)  +  &
                     MIN(MAX(adplai(i,8),dwlimit*plai(i,8)),uplimit*plai(i,8)) * tempu(i)
!
! upper canopy fractions
!
!        frac(i,1) = plai(i,1)            / max (avglaiu, epsilon)
!        frac(i,2) = plai(i,2) * dropu(i) / max (avglaiu, epsilon)
!        frac(i,3) = plai(i,3)            / max (avglaiu, epsilon)
!        frac(i,4) = plai(i,4)            / max (avglaiu, epsilon)
!        frac(i,5) = plai(i,5) * tempu(i) / max (avglaiu, epsilon)
!        frac(i,6) = plai(i,6)            / max (avglaiu, epsilon)
!        frac(i,7) = plai(i,7) * tempu(i) / max (avglaiu, epsilon)
!        frac(i,8) = plai(i,8) * tempu(i) / max (avglaiu, epsilon)
       !
       ! upper canopy fractions
       !
        frac(i,1) = MIN(MAX(adplai(i,1), dwlimit*plai(i,1)),uplimit*plai(i,1))	         / max (avglaiu(i), epsilon)
        frac(i,2) = MIN(MAX(adplai(i,2), dwlimit*plai(i,2)),uplimit*plai(i,2)) * dropu(i) / max (avglaiu(i), epsilon)
        frac(i,3) = MIN(MAX(adplai(i,3), dwlimit*plai(i,3)),uplimit*plai(i,3))	         / max (avglaiu(i), epsilon)
        frac(i,4) = MIN(MAX(adplai(i,4), dwlimit*plai(i,4)),uplimit*plai(i,4))	         / max (avglaiu(i), epsilon)
        frac(i,5) = MIN(MAX(adplai(i,5), dwlimit*plai(i,5)),uplimit*plai(i,5)) * tempu(i) / max (avglaiu(i), epsilon)
        frac(i,6) = MIN(MAX(adplai(i,6), dwlimit*plai(i,6)),uplimit*plai(i,6))	         / max (avglaiu(i), epsilon)
        frac(i,7) = MIN(MAX(adplai(i,7), dwlimit*plai(i,7)),uplimit*plai(i,7)) * tempu(i) / max (avglaiu(i), epsilon)
        frac(i,8) = MIN(MAX(adplai(i,8), dwlimit*plai(i,8)),uplimit*plai(i,8)) * tempu(i) / max (avglaiu(i), epsilon)
       !
       ! lower canopy single sided leaf area index (area-weighted)
       !
!        avglail = plai(i,9)                              + &
!                  plai(i,10) * min (templ(i), dropls(i)) + &
!                  plai(i,11) * min (templ(i), dropl4(i)) + &
!                  plai(i,12) * min (templ(i), dropl3(i))

        avglail(i) = MIN(MAX(adplai(i, 9),dwlimit*plai(i, 9)),uplimit*plai(i, 9))		              + &
                     MIN(MAX(adplai(i,10),dwlimit*plai(i,10)),uplimit*plai(i,10)) * min (templ(i), dropls(i)) + &
                     MIN(MAX(adplai(i,11),dwlimit*plai(i,11)),uplimit*plai(i,11)) * min (templ(i), dropl4(i)) + &
                     MIN(MAX(adplai(i,12),dwlimit*plai(i,12)),uplimit*plai(i,12)) * min (templ(i), dropl3(i))

       !
       ! lower canopy fractions
       !
!        frac(i,9)  = plai(i,9)                              /  &
!                     max (avglail(i), epsilon)
!!
!        frac(i,10) = plai(i,10) * min (templ(i), dropls(i)) /  &
!                     max (avglail(i), epsilon)
!!
!        frac(i,11) = plai(i,11) * min (templ(i), dropl4(i)) /  &
!                     max (avglail(i), epsilon)
!!
!        frac(i,12) = plai(i,12) * min (templ(i), dropl3(i)) /  &
!                     max (avglail(i), epsilon)

        frac(i,9)  = MIN(MAX(adplai(i, 9),dwlimit*plai(i, 9)),uplimit*plai(i, 9))                             /  max (avglail(i), epsilon)
        frac(i,10) = MIN(MAX(adplai(i,10),dwlimit*plai(i,10)),uplimit*plai(i,10)) * min (templ(i), dropls(i)) /  max (avglail(i), epsilon)
        frac(i,11) = MIN(MAX(adplai(i,11),dwlimit*plai(i,11)),uplimit*plai(i,11)) * min (templ(i), dropl4(i)) /  max (avglail(i), epsilon)
        frac(i,12) = MIN(MAX(adplai(i,12),dwlimit*plai(i,12)),uplimit*plai(i,12)) * min (templ(i), dropl3(i)) /  max (avglail(i), epsilon)

       !
       ! calculate the canopy leaf area index using the fractional vegetation cover
       !
        lai(i,1) = avglail(i) / fl(i)

       !**** DTP 2000/06/28 Modified this following discussion with Navin. We now
       !*    retain fu(i) derived from dynaveg and constrain it to a local value
       !*    "fu_phys" in the range 0.25 to 0.975 used in the canopy physics calcs.

        fu_phys = max (0.250_r8, min (0.9750_r8, fu(i)))
        lai(i,2) = avglaiu(i) / fu_phys
        lai(i,2) = avglaiu(i) / fu(i)
       !
       ! put a fix on canopy lais to avoid problems in physics
       !
        lai(i,1) = min (lai(i,1), 12.00_r8)
        lai(i,2) = min (lai(i,2), 12.00_r8)
       !
       ! ---------------------------------------------------------------------
       ! * * * update canopy height parameters * * *
       ! ---------------------------------------------------------------------
       !
       ! update lower canopy height parameters
       !
       ! note that they are based on vegetation fraction and not
       ! averaged over the entire gridcell
       !
        zbot(i,1)   =  0.050_r8
        ztop(i,1)   =  max (0.250_r8, lai(i,1) * 0.250_r8)
       !        
       ! constrain ztop to be at least 0.5 meter lower than 
       ! zbot for upper canopy
       !
        ztop(i,1) = min (ztop(i,1), zbot(i,2) - 0.50_r8)
       !
       ! end of loop
       !
    END DO !i = 1, npoi

    !totdepth = 0.0_r8
    !DO k = 1, nsoilay
    !      DO i=1,npoi
    !         totdepth(i) = totdepth(i) + hsoi(i,k) * 100.0_r8
    !      END DO  
    !END DO
    !
    ! normalization factors
    !
    !DO i=1,npoi
    !   inveg = NINT (vegtype0(i))
    !   frootnorm1(i) = 1.0_r8 - beta1(inveg) ** totdepth(i)
    !   frootnorm2(i) = 1.0_r8 - beta2(inveg) ** totdepth(i)
    !END DO 


!    !
!    ! calculate rooting profiles
!    !
!    DO k = 1, nsoilay
!       !
!          DO i=1,npoi
!             inveg = NINT (vegtype0(i))
!             IF (k.EQ.1) THEN
!                !
!                depth(i,k) = hsoi(i,k) * 100.0_r8
!                !
!                froot(i,k,1) = 1.0_r8 - beta1(inveg) ** depth(i,k)
!                froot(i,k,2) = 1.0_r8 - beta2(inveg) ** depth(i,k)
!                !
!             ELSE
!                !
!                depth(i,k) = depth(i,k-1) + hsoi(i,k) * 100.0_r8
!                !
!                froot(i,k,1) = (1.0_r8 - beta1(inveg) ** depth(i,k)) -  &
!                               (1.0_r8 - beta1(inveg) ** depth(i,k-1)) 
!                !
!                froot(i,k,2) = (1.0_r8 - beta2(inveg) ** depth(i,k)) -   & 
!                               (1.0_r8 - beta2(inveg) ** depth(i,k-1)) 
!                !
!             END IF
!             !
!             froot(i,k,1) = froot(i,k,1) / frootnorm1(i)
!             froot(i,k,2) = froot(i,k,2) / frootnorm2(i)
!             !
!          END DO 
!    END DO
      IF(TRIM(rootmode) == 'JACKSON')THEN
         CALL RootingProfilesJackson(nVegClass,nsoilay,npoi,hsoi,vegtype0,beta1,beta2,froot)
      ELSE  IF(TRIM(rootmode) == 'MILENA')THEN
         CALL RootingProfilesMilena(nsoilay,npoi,npft,hsoi,adnpp ,adtsoi,adwsoi,adwisoi,vegtype0,poros,rhow,froot)
      ELSE
         PRINT*, 'ERROR at rootmode parameter',TRIM(rootmode)
         STOP
      END IF 
      
    !
    ! return to main program
    ! 
    RETURN
  END SUBROUTINE pheno

!
      SUBROUTINE RootingProfilesMilena(nsoilay,npoi,npft,hsoi,adnpp ,adtsoi,adwsoi,adwisoi,vegtype0,poros,rhow,froot)
        IMPLICIT  NONE
        INTEGER, INTENT(IN   ) :: npoi
        INTEGER, INTENT(IN   ) :: npft
        INTEGER, INTENT(IN   ) :: nsoilay
        REAL(KIND=r8), INTENT(IN   ) :: hsoi         (npoi,nsoilay+1)! soil layer thickness (m)
        REAL(KIND=r8), INTENT(IN   ) :: vegtype0     (npoi) ! fixed vegetation map
        REAL(KIND=r8), INTENT(IN   ) :: adnpp     (npoi,npft)! global! daily  total npp for each plant type (kg-C/m**2/day)
        REAL(KIND=r8), INTENT(IN   ) :: adtsoi    (npoi)! global -daily average ! soil temperature for each layer (C)    ! averages for first 2 layers (0.3 meters) of soil  	 
        REAL(KIND=r8), INTENT(IN   ) :: adwsoi    (npoi)! global -daily average ! fraction of soil pore space containing liquid water ! averages for first 2 layers (0.3 meters) of soil(fraction) 
        REAL(KIND=r8), INTENT(IN   ) :: adwisoi   (npoi)! global -daily average ! fraction of soil pore space containing ice    ! averages for first 2 layers (0.3 meters) of soil  	 
        REAL(KIND=r8), INTENT(IN   ) :: poros     (npoi,nsoilay)! global ! porosity (mass of h2o per unit vol at sat / rhow)
        REAL(KIND=r8), INTENT(IN   ) :: rhow !global  ! density of liquid water (all types) (kg m-3)
        REAL(KIND=r8), INTENT(INOUT) :: froot    (npoi,nsoilay,2)! global! fraction of root in soil layer 
        !
        ! LOCAL VARIABEL
        !
        REAL(KIND=r8)    :: totdepth  (npoi)  ! total soil depth
        REAL(KIND=r8)    :: frootnorm1(npoi)  ! normalization factor for Jackson rooting profile,low
        REAL(KIND=r8)    :: frootnorm2 (npoi) ! normalization factor for Jackson rooting profile, up
        REAL(KIND=r8)    :: depth(nsoilay)   ! soil layer depth (cm)
        REAL(KIND=r8)    :: depth_aux(npoi,nsoilay)   ! soil layer depth (cm)

        REAL(KIND=r8), PARAMETER    :: beta_lw_canopy(nVegClass)     = RESHAPE ( (/ &
        !  beta_lw_canopy        !                                                  beta_lw_canopy
        2.100_r8, &    !  1: tropical evergreen forest / woodland                 0.962
        1.800_r8, &    !  2: tropical deciduous forest / woodland                 0.961
        2.000_r8, &    !  3: temperate evergreen broadleaf forest / woodland      0.966
        2.700_r8, &    !  4: temperate evergreen conifer forest / woodland        0.966
        2.100_r8, &    !  5: temperate deciduous forest / woodland                0.965
        1.500_r8, &    !  6: boreal evergreen forest / woodland                   0.960
        1.400_r8, &    !  7: boreal deciduous forest / woodland                   0.950
        1.800_r8, &    !  8: mixed forest / woodland                              0.960
        1.300_r8, &    !  9: savanna                                              0.962
        1.000_r8, &    ! 10: grassland / steppe                                   0.952
        1.200_r8, &    ! 11: dense shrubland                                      0.970
        1.100_r8, &    ! 12: open shrubland                                       0.950
        1.000_r8, &    ! 13: tundra                                               0.914
        2.000_r8, &    ! 14: desert                                               0.970
        1.900_r8  &    ! 15: polar desert / rock / ice                            0.970
        /), (/nVegClass/) )!---->  grassland / shrub systems  

        REAL(KIND=r8), PARAMETER    :: beta_up_canopy(nVegClass)     = RESHAPE ( (/ &
        !  beta_up_canopy        !                                                 beta_up_canopy
        2.600_r8, &    !  1: tropical evergreen forest / woodland                 0.962
        2.300_r8, &    !  2: tropical deciduous forest / woodland                 0.961
        2.300_r8, &    !  3: temperate evergreen broadleaf forest / woodland      0.966
        2.400_r8, &    !  4: temperate evergreen conifer forest / woodland        0.966
        2.200_r8, &    !  5: temperate deciduous forest / woodland                0.965
        1.300_r8, &    !  6: boreal evergreen forest / woodland                   0.960
        1.300_r8, &    !  7: boreal deciduous forest / woodland                   0.950
        2.500_r8, &    !  8: mixed forest / woodland                              0.960
        1.500_r8, &    !  9: savanna                                              0.962
        1.200_r8, &    ! 10: grassland / steppe                                   0.952
        1.200_r8, &    ! 11: dense shrubland                                      0.970
        1.100_r8, &    ! 12: open shrubland                                       0.950
        1.000_r8, &    ! 13: tundra                                               0.914
        2.100_r8, &    ! 14: desert                                               0.970
        2.000_r8  &    ! 15: polar desert / rock / ice                            0.970
        /), (/nVegClass/) )!---->  grassland / shrub systems  

        REAL(KIND=r8), PARAMETER    :: theta_lw_canopy(nVegClass)     = RESHAPE ( (/ &
        !  theta_lw_canopy        !                                                theta_lw_canopy
        50.0_r8, &    !  1: tropical evergreen forest / woodland                 0.962
        40.0_r8, &    !  2: tropical deciduous forest / woodland                 0.961
        45.0_r8, &    !  3: temperate evergreen broadleaf forest / woodland      0.966
        33.0_r8, &    !  4: temperate evergreen conifer forest / woodland        0.966
        28.0_r8, &    !  5: temperate deciduous forest / woodland                0.965
        29.0_r8, &    !  6: boreal evergreen forest / woodland                   0.960
        27.0_r8, &    !  7: boreal deciduous forest / woodland                   0.950
        28.0_r8, &    !  8: mixed forest / woodland                              0.960
        50.0_r8, &    !  9: savanna                                              0.962
        18.0_r8, &    ! 10: grassland / steppe                                   0.952
        40.0_r8, &    ! 11: dense shrubland                                      0.970
        40.0_r8, &    ! 12: open shrubland                                       0.950
        15.0_r8, &    ! 13: tundra                                               0.914
        13.0_r8, &    ! 14: desert                                               0.970
        13.0_r8  &    ! 15: polar desert / rock / ice                            0.970
        /), (/nVegClass/) )!---->  grassland / shrub systems  

        REAL(KIND=r8), PARAMETER    :: theta_up_canopy(nVegClass)     = RESHAPE ( (/ &
        !  theta_up_canopy        !                                                theta_up_canopy
        200.0_r8, &    !  1: tropical evergreen forest / woodland                 0.962
        150.0_r8, &    !  2: tropical deciduous forest / woodland                 0.961
        150.0_r8, &    !  3: temperate evergreen broadleaf forest / woodland      0.966
        120.0_r8, &    !  4: temperate evergreen conifer forest / woodland        0.966
         60.0_r8, &    !  5: temperate deciduous forest / woodland                0.965
         35.0_r8, &    !  6: boreal evergreen forest / woodland                   0.960
         40.0_r8, &    !  7: boreal deciduous forest / woodland                   0.950
         50.0_r8, &    !  8: mixed forest / woodland                              0.960
         70.0_r8, &    !  9: savanna                                              0.962
         23.0_r8, &    ! 10: grassland / steppe                                   0.952
         60.0_r8, &    ! 11: dense shrubland                                      0.970
         50.0_r8, &    ! 12: open shrubland                                       0.950
         16.0_r8, &    ! 13: tundra                                               0.914
         13.0_r8, &    ! 14: desert                                               0.970
         13.0_r8  &    ! 15: polar desert / rock / ice                            0.970
         /), (/nVegClass/) )!---->  grassland / shrub systems  

         REAL(KIND=r8), PARAMETER    :: delta_lw_canopy(nVegClass)     = RESHAPE ( (/ &
         !  delta_lw_canopy        !                                                delta_lw_canopy
         1.000_r8, &    !  1: tropical evergreen forest / woodland                 0.962
         1.000_r8, &    !  2: tropical deciduous forest / woodland                 0.961
         1.000_r8, &    !  3: temperate evergreen broadleaf forest / woodland      0.966
         1.000_r8, &    !  4: temperate evergreen conifer forest / woodland        0.966
         1.000_r8, &    !  5: temperate deciduous forest / woodland                0.965
         1.000_r8, &    !  6: boreal evergreen forest / woodland                   0.960
         1.000_r8, &    !  7: boreal deciduous forest / woodland                   0.950
         1.000_r8, &    !  8: mixed forest / woodland                              0.960
         1.000_r8, &    !  9: savanna                                              0.962
         1.000_r8, &    ! 10: grassland / steppe                                   0.952
         1.000_r8, &    ! 11: dense shrubland                                      0.970
         1.000_r8, &    ! 12: open shrubland                                       0.950
         1.000_r8, &    ! 13: tundra                                               0.914
         1.000_r8, &    ! 14: desert                                               0.970
         1.000_r8  &    ! 15: polar desert / rock / ice                            0.970
         /), (/nVegClass/) )!---->  grassland / shrub systems  

         REAL(KIND=r8), PARAMETER    :: delta_up_canopy(nVegClass)     = RESHAPE ( (/ &
         !  delta_up_canopy        !                                                delta_up_canopy
         1.000_r8, &    !  1: tropical evergreen forest / woodland                 0.962
         1.000_r8, &    !  2: tropical deciduous forest / woodland                 0.961
         1.000_r8, &    !  3: temperate evergreen broadleaf forest / woodland      0.966
         1.000_r8, &    !  4: temperate evergreen conifer forest / woodland        0.966
         1.000_r8, &    !  5: temperate deciduous forest / woodland                0.965
         1.000_r8, &    !  6: boreal evergreen forest / woodland                   0.960
         1.000_r8, &    !  7: boreal deciduous forest / woodland                   0.950
         1.000_r8, &    !  8: mixed forest / woodland                              0.960
         1.000_r8, &    !  9: savanna                                              0.962
         1.000_r8, &    ! 10: grassland / steppe                                   0.952
         1.000_r8, &    ! 11: dense shrubland                                      0.970
         1.000_r8, &    ! 12: open shrubland                                       0.950
         1.000_r8, &    ! 13: tundra                                               0.914
         1.000_r8, &    ! 14: desert                                               0.970
         1.000_r8  &    ! 15: polar desert / rock / ice                            0.970
         /), (/nVegClass/) )!---->  grassland / shrub systems  
         REAL(KIND=r8)    :: x,maxdepth,nlayers

         REAL(KIND=r8)    :: hsoi_limit =0.3_r8      !       [m]    ! (Sistema Internacional at S.I.)
         !         
	 REAL(KIND=r8)    :: zdepth  ! total soil depth
         REAL(KIND=r8)    :: zdepth30cm  (npoi)  ! total soil depth
         REAL(KIND=r8)    :: poros30cm   (npoi)  ! total soil poros
         REAL(KIND=r8)    :: Theta     (npoi)  !  !Water Content of Soil Layer   [ kg m-2 ]
         REAL(KIND=r8)    :: adnpptot  (npoi)! ! daily  total npp (kg-C/m**2/day)
         INTEGER :: i,j,k,inveg,lrec
         !
         !Water Content of Soil Layer   =>  Theta  [ kg m-2 ]

         !                  Theta                            [ kg m-2 ]
         !wsoi=---------------------------- =    ---------------------------
         !              poros*hsoi*rhow            [ %]    *   [m] *  [kg m-3]

         !Water Content of Soil Layer    Theta  [ kg m-2 ]  =   poros * wsoi * hsoi * rhow

         !                            rho_h20_soil
         !poros = [%]   =----------------------------------                 ! Porosity : porosity (volume fraction -> porosity (mass of h2o per unit vol at sat / rhow))
         !                          rho_h20_soil_sat.

         !wsoi  =  [%]

         !rhow  = 1.0e+3_r8   [kg m-3]  !(Sistema Internacional at S.I.)

         !hsoi =0.3_r8             [m]    ! (Sistema Internacional at S.I.)
         !
         !
         ! averages for first 2 layers (0.3 meters) of soil
         !
         nlayers=0
         zdepth=0.0_r8
         zdepth30cm=0.0_r8
         poros30cm =0.0_r8
         DO k = 1, nsoilay
            DO i=1,npoi
               zdepth = zdepth + hsoi(i,k)
               IF(zdepth30cm(i)<=0.30_r8)THEN!m            
                  nlayers       =  nlayers + 1.0_r8
                  zdepth30cm(i) =  zdepth30cm(i) + hsoi(i,k)
                  poros30cm(i)  =  poros30cm(i) + poros(i,k)
               END IF
            END DO
         END DO            
         DO i=1,npoi
            poros30cm(i)=poros30cm(i)/nlayers
            
            !Water Content of Soil Layer    Theta  [ kg m-2 ]  =   poros * wsoi * hsoi * rhow
         
            Theta(i)  =   poros30cm(i) * adwsoi(i) * zdepth30cm(i) * rhow
         END DO 
         !
         ! daily  total npp for each plant type (kg-C/m**2/day)
         !
         ! determine total ecosystem positive npp (changed by exist at begin
         ! of subroutine). Different from sum of monthly and daily npp)
         !
         DO i=1,npoi

            adnpptot(i) = max(0.0_r8,adnpp(i,1))  + max(0.0_r8,adnpp(i,2)) + &
                          max(0.0_r8,adnpp(i,3))  + max(0.0_r8,adnpp(i,4)) + &
                          max(0.0_r8,adnpp(i,5))  + max(0.0_r8,adnpp(i,6)) + &
                          max(0.0_r8,adnpp(i,7))  + max(0.0_r8,adnpp(i,8)) + &
                          max(0.0_r8,adnpp(i,9))  + max(0.0_r8,adnpp(i,10)) + &
                          max(0.0_r8,adnpp(i,11)) + max(0.0_r8,adnpp(i,12))
         END DO 
         ! ************************************************************************
         ! define rooting profiles
         ! ************************************************************************
         !
         ! define rooting profiles based upon data published in:
         !
         ! Milena Dantas et al., 2020: 
         !
         ! and
         !
         ! Jackson et al., 1997:  A global budget for fine root biomass, 
         ! surface area, and nutrient contents, Proceedings of the National
         ! Academy of Sciences, 94, 7362-7366.
         !
         DO k = 1, nsoilay
            !
               DO i=1,npoi
                  inveg = NINT (vegtype0(i))
                  IF (k.eq.1) THEN
                     depth(k) = hsoi(i,k) * 100.0_r8
                  ELSE
                     depth(k) = depth(k-1) + hsoi(i,k) * 100.0_r8
                  END IF
                  ! beta_lw_canopy  é o parametro de forma
                  ! theta_lw_canopy é o parametro de escala 
                  ! delta_lw_canopy é o parametro de localização
                  maxdepth= depth(k)
               END DO
         END DO

         frootnorm1=0.0
         frootnorm2=0.0
         DO k = 1, nsoilay
            !
               DO i=1,npoi
                  inveg = NINT (vegtype0(i))
                  IF (k.eq.1) THEN
                     depth(k) = hsoi(i,k) * 100.0_r8
                  ELSE
                     depth(k) = depth(k-1) + hsoi(i,k) * 100.0_r8
                  END IF
                  ! beta_lw_canopy  é o parametro de forma
                  ! theta_lw_canopy é o parametro de escala 
                  ! delta_lw_canopy é o parametro de localização
                  x=depth(k)
                  IF( x < delta_lw_canopy (inveg))STOP 'ERROR depth < delta_lw_canopy (inveg)'
                  frootnorm1(i) =frootnorm1(i) +  (beta_lw_canopy(inveg)/theta_lw_canopy(inveg))* &
                   (((x-delta_lw_canopy(inveg))/theta_lw_canopy(inveg))**(beta_lw_canopy(inveg)-1.0))* &
                   (exp(-(((x-delta_lw_canopy(inveg))/theta_lw_canopy(inveg))**beta_lw_canopy(inveg))))
                  frootnorm2(i) =frootnorm2(i) +  (beta_up_canopy(inveg)/theta_up_canopy(inveg))* &
                   (((x-delta_up_canopy(inveg))/theta_up_canopy(inveg))**(beta_up_canopy(inveg)-1.0))* &
                  (exp(-(((x-delta_up_canopy(inveg))/theta_up_canopy(inveg))**beta_up_canopy(inveg))))
               END DO
         END DO

         DO k = 1, nsoilay
            !
               DO i=1,npoi
                  inveg = NINT (vegtype0(i))
                  IF (k.eq.1) THEN
                     depth(k) = hsoi(i,k) * 100.0_r8
                  ELSE
                     depth(k) = depth(k-1) + hsoi(i,k) * 100.0_r8
                  END IF
                  ! beta_lw_canopy  é o parametro de forma
                  ! theta_lw_canopy é o parametro de escala 
                  ! delta_lw_canopy é o parametro de localização
                  x=depth(k)
                  IF( x < delta_lw_canopy (inveg))STOP 'ERROR depth < delta_lw_canopy (inveg)'
                  froot(i,k,1) = (beta_lw_canopy(inveg)/theta_lw_canopy(inveg))* &
                 (((x-delta_lw_canopy(inveg))/theta_lw_canopy(inveg))**(beta_lw_canopy(inveg)-1.0))* &
                 (exp(-(((x-delta_lw_canopy(inveg))/theta_lw_canopy(inveg))**beta_lw_canopy(inveg))))
                  froot(i,k,2) = (beta_up_canopy(inveg)/theta_up_canopy(inveg))* &
                  (((x-delta_up_canopy(inveg))/theta_up_canopy(inveg))**(beta_up_canopy(inveg)-1.0))* &
                   (exp(-(((x-delta_up_canopy(inveg))/theta_up_canopy(inveg))**beta_up_canopy(inveg))))

                  froot(i,k,1) = froot(i,k,1) / frootnorm1(i)
                  froot(i,k,2) = froot(i,k,2) / frootnorm2(i)
                  !
               END DO
         END DO
      !
      ! return to main program
      !
      END SUBROUTINE RootingProfilesMilena
      !


      SUBROUTINE RootingProfilesJackson(nVegClass,nsoilay,npoi,hsoi,vegtype0,beta1,beta2,froot)
       IMPLICIT  NONE
       INTEGER      , INTENT(IN   ) :: nVegClass
       INTEGER      , INTENT(IN   ) :: nsoilay
       INTEGER      , INTENT(IN   ) :: npoi
       REAL(KIND=r8), INTENT(IN   ) :: hsoi     (npoi,nsoilay+1)   ! global ! soil layer thickness (m)
       REAL(KIND=r8), INTENT(IN   ) :: vegtype0 (npoi)      ! annual vegetation type - ibis classification
       REAL(KIND=r8), INTENT(IN   ) :: beta1    (nVegClass)
       REAL(KIND=r8), INTENT(IN   ) :: beta2    (nVegClass)
       REAL(KIND=r8), INTENT(INOUT) :: froot    (npoi,nsoilay,2)! global! fraction of root in soil layer 
       !
       ! LOCAL VARIABEL
       !
      REAL(KIND=r8)    :: totdepth(npoi)        ! total soil depth
      REAL(KIND=r8)    :: frootnorm1(npoi)      ! normalization factor for Jackson rooting profile,low
      REAL(KIND=r8)    :: frootnorm2(npoi)      ! normalization factor for Jackson rooting profile, up
      REAL(KIND=r8)    :: depth(npoi,nsoilay)   ! soil layer depth (cm)

      INTEGER :: i,j,k,inveg

      ! ************************************************************************
      ! define rooting profiles
      ! ************************************************************************
      !
      ! define rooting profiles based upon data published in:
      !
      ! Jackson et al., 1996:  A global analysis of root distributions
      ! for terrestrial biomes, Oecologia, 108, 389-411.
      !
      ! and
      !
      ! Jackson et al., 1997:  A global budget for fine root biomass, 
      ! surface area, and nutrient contents, Proceedings of the National
      ! Academy of Sciences, 94, 7362-7366.
      !
      ! rooting profiles are defined by the "beta" parameter
      !
      ! beta1 is assigned to the lower vegetation layer (grasses and shrubs)
      ! beta2 is assigned to the upper vegetation layer (trees)
      !
      ! according to Jackson et al. (1996, 1997), the values of beta
      ! typically fall in the following range
      !
      ! note that the 1997 paper specifically discusses the distribution
      ! of *fine roots* (instead of total root biomass), which may be more
      ! important for water and nutrient uptake
      !
      ! --------------                 ------------   ------------
      ! forest systems                 beta2 (1996)   beta2 (1997)
      ! --------------                 ------------   ------------
      ! tropical evergreen forest:        0.962          0.972
      ! tropical deciduous forest:        0.961          0.982
      ! temperate conifer forest:         0.976          0.980
      ! temperate broadleaf forest:       0.966          0.967
      ! all tropical/temperate forest:    0.970  
      ! boreal forest:                    0.943          0.943
      ! all trees:                                       0.976
      !
      ! -------------------------      ------------   ------------
      ! grassland / shrub systems      beta1 (1996)   beta1 (1997)
      ! -------------------------      ------------   ------------
      ! tropical grassland / savanna:     0.972          0.972
      ! temperate grassland:              0.943          0.943
      ! all grasses:                      0.952          0.952
      ! schlerophyllous shrubs:           0.964          0.950
      ! all shrubs:                       0.978          0.975
      ! crops:                            0.961
      ! desert:                           0.975          0.970
      ! tundra:                           0.914
      !
      ! --------------                 ------------
      ! all ecosystems                 beta  (1996)
      ! --------------                 ------------
      ! all ecosystems:                   0.966
      !
      ! for global simulations, we typically assign the following
      ! values to the beta parameters
      !
      ! beta1 = 0.950, which is typical for tropical/temperate grasslands
      ! beta2 = 0.970, which is typical for tropical/temperate forests
      !
      ! however, these values could be (and should be) further refined
      ! when using the model for specific regions
      ! 
      !      beta1 = 0.950  ! for lower layer herbaceous plants
      !      beta2 = 0.975  ! for upper layer trees
      !
      ! calculate total depth in centimeters
      !
      !
      totdepth = 0.0_r8
      DO k = 1, nsoilay
         DO i=1,npoi
            totdepth(i) = totdepth(i) + hsoi(i,k) * 100.0_r8
         END DO  
      END DO
    !
    ! normalization factors
    !
    DO i=1,npoi
       inveg = NINT (vegtype0(i))
       frootnorm1(i) = 1.0_r8 - beta1(inveg) ** totdepth(i)
       frootnorm2(i) = 1.0_r8 - beta2(inveg) ** totdepth(i)
    END DO 
    !
    ! calculate rooting profiles
    !
    DO k = 1, nsoilay
       !
       DO i=1,npoi

         inveg = NINT (vegtype0(i))

         IF (k.eq.1) THEN
!
            depth(i,k) = hsoi(i,k) * 100.0_r8
!
            froot(i,k,1) = 1.0_r8 - beta1(inveg) ** depth(i,k)
            froot(i,k,2) = 1.0_r8 - beta2(inveg) ** depth(i,k)
!
         ELSE
!
            depth(i,k) = depth(i,k-1) + hsoi(i,k) * 100.0_r8

!
            froot(i,k,1) = (1.0_r8 - beta1(inveg) ** depth(i,k)) -  &
                           (1.0_r8 - beta1(inveg) ** depth(i,k-1)) 
!
            froot(i,k,2) = (1.0_r8 - beta2(inveg) ** depth(i,k)) -   & 
                           (1.0_r8 - beta2(inveg) ** depth(i,k-1)) 
!
         END IF
!
         froot(i,k,1) = froot(i,k,1) / frootnorm1(i)
         froot(i,k,2) = froot(i,k,2) / frootnorm2(i)
!
       END DO
    END DO
!
! return to main program
!
END SUBROUTINE RootingProfilesJackson


!
!
! ---------------------------------------------------------------------
  SUBROUTINE dynaveg1 (isimfire , &! INTENT(IN   )
       tauwood0 , &! INTENT(IN   )
       tauwood  , &! INTENT(OUT  )
       tauleaf  , &! INTENT(IN   )
       tauroot  , &! INTENT(IN   )
       xminlai  , &! INTENT(IN   )
       falll    , &! INTENT(OUT  )
       fallr    , &! INTENT(OUT  )
       fallw    , &! INTENT(OUT  )
       cdisturb , &! INTENT(OUT  )
       exist    , &! INTENT(IN   )
       aleaf    , &! INTENT(IN   )
       awood    , &! INTENT(IN   )
       cbiol    , &! INTENT(INOUT) global
       cbior    , &! INTENT(INOUT) global
       cbiow    , &! INTENT(INOUT) global
       aroot    , &! INTENT(IN   )
       disturbf , &! INTENT(OUT  )
       disturbo , &! INTENT(OUT  )
       firefac  , &! INTENT(IN   )
       totlit   , &! INTENT(IN   )
       specla   , &! INTENT(IN   )
       plai     , &! INTENT(INOUT) local
       biomass  , &! INTENT(OUT  )
       totlaiu  , &! INTENT(INOUT) local
       totlail  , &! INTENT(INOUT) local
       totbiou  , &! INTENT(INOUT) local
       totbiol  , &! INTENT(OUT  )
       fu   , &! INTENT(OUT  )
       woodnorm , &! INTENT(IN   )
       fl   , &! INTENT(OUT  )
       zbot     , &! INTENT(OUT  )
       ztop     , &! INTENT(OUT  )
       sai      , &! INTENT(OUT  )
       sapfrac  , &! INTENT(OUT  )
       vegtype0 , &! INTENT(OUT  )
       gdd5     , &! INTENT(IN   )
       gdd0     , &! INTENT(IN   )
       aynpp    , &! INTENT(INOUT) global
       ayanpp   , &! INTENT(OUT  )
       ayneetot , &! INTENT(INOUT) global
       ayanpptot, &! INTENT(OUT  )
       npoi     , &!
       npft       )! , isim_ac, year)
    ! ---------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN   ) :: npoi                   ! total number of land points
    INTEGER, INTENT(IN   ) :: npft                   ! number of plant functional types
    REAL(KIND=r8)   , INTENT(INOUT) :: aynpp    (npoi,npft)   ! annual total npp for each plant type(kg-c/m**2/yr)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ayanpp   (npoi,npft)   ! annual above-ground npp for each plant type(kg-c/m**2/yr)
    REAL(KIND=r8)   , INTENT(INOUT) :: ayneetot (npoi)        ! annual total NEE for ecosystem (kg-C/m**2/yr)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ayanpptot(npoi)        ! annual above-ground npp for ecosystem (kg-c/m**2/yr)
    REAL(KIND=r8)   , INTENT(OUT  ) :: falll   (npoi)         ! annual leaf litter fall                      (kg_C m-2/year)
    REAL(KIND=r8)   , INTENT(OUT  ) :: fallr   (npoi)         ! annual root litter input                     (kg_C m-2/year)
    REAL(KIND=r8)   , INTENT(OUT  ) :: fallw   (npoi)         ! annual wood litter fall                      (kg_C m-2/year)
    REAL(KIND=r8)   , INTENT(OUT  ) :: cdisturb(npoi)         ! annual amount of vegetation carbon lost 
    ! to atmosphere due to fire  (biomass burning) (kg_C m-2/year)
    REAL(KIND=r8)   , INTENT(IN   ) :: exist   (npoi,npft)    ! probability of existence of each plant functional type in a gridcell
    REAL(KIND=r8)   , INTENT(IN   ) :: aleaf   (npft)         ! carbon allocation fraction to leaves
    REAL(KIND=r8)   , INTENT(IN   ) :: awood   (npft)         ! carbon allocation fraction to wood 
    REAL(KIND=r8)   , INTENT(INOUT) :: cbiol   (npoi,npft)    ! carbon in leaf biomass pool (kg_C m-2)
    REAL(KIND=r8)   , INTENT(INOUT) :: cbior   (npoi,npft)    ! carbon in fine root biomass pool (kg_C m-2)
    REAL(KIND=r8)   , INTENT(INOUT) :: cbiow   (npoi,npft)    ! carbon in woody biomass pool (kg_C m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: aroot   (npft)         ! carbon allocation fraction to fine roots
    REAL(KIND=r8)   , INTENT(OUT  ) :: disturbf(npoi)         ! annual fire disturbance regime (m2/m2/yr)
    REAL(KIND=r8)   , INTENT(OUT  ) :: disturbo(npoi)         ! fraction of biomass pool lost every year to disturbances other than fire
    REAL(KIND=r8)   , INTENT(IN   ) :: firefac (npoi)         ! factor that respresents the annual average fuel
    ! dryness of a grid cell, and hence characterizes the readiness to burn
    REAL(KIND=r8)   , INTENT(IN   ) :: totlit  (npoi)         ! total carbon in all litter pools (kg_C m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: specla  (npft)         ! specific leaf area (m**2/kg) 
    REAL(KIND=r8)   , INTENT(INOUT) :: plai    (npoi,npft)    ! total leaf area index of each plant functional type
    REAL(KIND=r8)   , INTENT(OUT  ) :: biomass (npoi,npft)    ! total biomass of each plant functional type  (kg_C m-2)
    REAL(KIND=r8)   , INTENT(INOUT) :: totlaiu (npoi)         ! total leaf area index for the upper canopy
    REAL(KIND=r8)   , INTENT(INOUT) :: totlail (npoi)         ! total leaf area index for the lower canopy
    REAL(KIND=r8)   , INTENT(INOUT) :: totbiou (npoi)         ! total biomass in the upper canopy (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: totbiol (npoi)         ! total biomass in the lower canopy (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: fu      (npoi)         ! fraction of overall area covered by upper canopy
    REAL(KIND=r8)   , INTENT(IN   ) :: woodnorm       ! value of woody biomass for upper canopy closure
    ! (ie when wood = woodnorm fu = 1.0) (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: fl      (npoi)         ! fraction of snow-free area covered by lower  canopy
    REAL(KIND=r8)   , INTENT(INOUT  ) :: zbot    (npoi,2)       ! height of lowest branches above ground (m)
    REAL(KIND=r8)   , INTENT(INOUT  ) :: ztop    (npoi,2)       ! height of plant top above ground (m)
    REAL(KIND=r8)   , INTENT(OUT  ) :: sai     (npoi,2)       ! current single-sided stem area index
    REAL(KIND=r8)   , INTENT(OUT  ) :: sapfrac (npoi)         ! fraction of woody biomass that is in sapwood
    REAL(KIND=r8)   , INTENT(INOUT) :: vegtype0(npoi)         ! annual vegetation type - ibis classification
    REAL(KIND=r8)   , INTENT(IN   ) :: gdd5    (npoi)         ! growing degree days > 5C
    REAL(KIND=r8)    :: gdd0    (npoi)         ! growing degree days > 0C 
    REAL(KIND=r8)   , INTENT(IN   ) :: tauwood0(npft)   ! normal (unstressed) turnover time for wood biomass (years)
    REAL(KIND=r8)   , INTENT(OUT  ) :: tauwood (npoi,npft)   ! wood biomass turnover time constant (years)
    REAL(KIND=r8)   , INTENT(IN   ) :: tauleaf (npft)   ! foliar biomass turnover time constant (years)
    REAL(KIND=r8)   , INTENT(IN   ) :: tauroot (npft)   ! fine root biomass turnover time constant (years)
    REAL(KIND=r8)   , INTENT(IN   ) :: xminlai          ! Minimum LAI for each existing PFT
    !
    ! Arguments
    !
    INTEGER, INTENT(IN   ) :: isimfire  ! fire switch
    ! isim_ac   ! age-class dynamics switch
    ! year      ! year of simulation

    !REAL(KIND=r8)    :: pdist      ! probability of other disturbance types....

    REAL(KIND=r8)    , PARAMETER :: pfire = 1.00_r8 ! probability of fire -- should be determined externally.
    ! for now we just assume it occurs all the time


    !
    ! local variables
    !
    INTEGER :: i     ! gridcell counter      
    INTEGER :: j     ! gridcell counter
    !
    REAL(KIND=r8)    :: sapspeed      ! in mm/day
    REAL(KIND=r8)    :: trans        ! (2.5 mm/day) 
    REAL(KIND=r8)    :: saparea       ! in m**2
    REAL(KIND=r8)    :: sapvolume     ! in m**3
    REAL(KIND=r8)    :: denswood      ! kg/m**3
    REAL(KIND=r8)    :: wood          ! total amount of woody biomass in gridcell
    REAL(KIND=r8)    :: taufin        !
    !*      REAL(KIND=r8)    :: xminlai               !
    !
    !*      REAL(KIND=r8)
    !*      REAL(KIND=r8)    :: aleaf(npft),   ! allocation fraction to leaves
    !*      REAL(KIND=r8)    :: aroot(npft),   ! allocation fraction to fine roots
    !*      REAL(KIND=r8)    :: awood(npft),   ! allocation fraction to wood
    !*      REAL(KIND=r8)    :: tauleaf(npft), ! turnover time of carbon in leaves (years)
    !*      REAL(KIND=r8)    :: tauroot(npft), ! turnover time of carbon in fine roots (years)
    !*      REAL(KIND=r8)    :: tauwood(npft)   ! turnover time of carbon in wood (years)
    !*      REAL(KIND=r8)    :: tauwood0(npft) ! normal (unstressed) turnover time
    !
    ! ibis uses a small number of plant functional types:
    !
    !  1: tropical broadleaf evergreen tree
    !  2: tropical broadleaf drought-deciduous trees
    !  3: warm-temperate broadleaf evergreen tree
    !  4: temperate conifer evergreen tree
    !  5: temperate broadleaf cold-deciduous tree
    !  6: boREAL(KIND=r8) conifer evergreen tree
    !  7: boREAL(KIND=r8) broadleaf cold-deciduous tree
    !  8: boREAL(KIND=r8) conifer cold-deciduous tree
    !  9: evergreen shrub
    ! 10: deciduous shrub
    ! 11: warm (c4) grass
    ! 12: cool (c3) grass
    !
    ! ---------------------------------------------------------------------
    ! * * * specify biomass turnover parameters (years) * * *
    ! ---------------------------------------------------------------------
    !
    !      data tauleaf / 1.01,   ! tropical broadleaf evergreen trees
    !     >               1.00,   ! tropical broadleaf drought-deciduous trees
    !     >               1.00,   ! warm-temperate broadleaf evergreen trees
    !     >               2.00,   ! temperate conifer evergreen trees
    !     >               1.00,   ! temperate broadleaf cold-deciduous trees
    !     >               2.50,   ! boREAL(KIND=r8) conifer evergreen trees
    !     >               1.00,   ! boREAL(KIND=r8) broadleaf cold-deciduous trees
    !     >               1.00,   ! boREAL(KIND=r8) conifer cold-deciduous trees
    !     >               1.50,   ! evergreen shrubs
    !     >               1.00,   ! deciduous shrubs
    !     >               1.25,   ! warm (c4) grasses
    !     >               1.50 /  ! cool (c3) grasses
    !
    !      data tauwood0 / 25.0,  ! tropical broadleaf evergreen trees
    !     >                25.0,  ! tropical broadleaf drought-deciduous trees
    !     >                25.0,  ! warm-temperate broadleaf evergreen trees
    !     >                50.0,  ! temperate conifer evergreen trees
    !     >                50.0,  ! temperate broadleaf cold-deciduous trees
    !     >               100.0,  ! boREAL(KIND=r8) conifer evergreen trees
    !     >               100.0,  ! boREAL(KIND=r8) broadleaf cold-deciduous trees
    !     >               100.0,  ! boREAL(KIND=r8) conifer cold-deciduous trees
    !     >                 5.0,  ! evergreen shrubs
    !     >                 5.0,  ! deciduous shrubs
    !     >               999.0,  ! warm (c4) grasses
    !     >               999.0 / ! cool (c3) grasses
    !
    ! ---------------------------------------------------------------------
    ! * * * apply disturbances * * *
    ! ---------------------------------------------------------------------
    !
    ! set fixed disturbance regime
    !    
    DO i = 1, npoi
       disturbf(i) = 0.0050_r8
       disturbo(i) = 0.0050_r8
    END DO
    !**** DTP 2000/08/10. One can do a decent test of ACME by setting
    !*    these disturbance rates to zero. With these values, the area
    !*    disturbed each year will be zero so the distribution of biomass
    !*    and PFTs across the domain should be identical to those 
    !*    resulting from a run of standard IBIS with zero disturbance

    !*        disturbf(i) = 0.0  ! Test with zero disturbance rate 
    !*        disturbo(i) = 0.0  ! (This should equal reference sim).

    !
    ! call fire disturbance routine
    !**** DTP 2001/03/06: In general isimfire should be set to zero if isim_ac
    !*    is set to 1 (but what does isimfire REAL(KIND=r8)ly do?)
    !
    IF (isimfire.eq.1) THEN
       CALL fire(npoi      , & ! INTENT(IN)
            firefac   , & ! INTENT(IN)
            totlit    , & ! INTENT(IN)
            disturbf  ,&
            vegtype0   ) ! INTENT(OUT  )
    END IF


    !
    ! begin global grid
    !
    DO i = 1, npoi
       !
       ! ---------------------------------------------------------------------
       ! * * * initialize vegetation dynamics pools * * *
       ! ---------------------------------------------------------------------
       !
       ! zero out litter fall fields
       !
       falll(i) = 0.00_r8
       fallr(i) = 0.00_r8
       fallw(i) = 0.00_r8
       !
       ! zero out carbon lost due to disturbance
       ! 
       cdisturb(i) = 0.00_r8
       !
       wood = 0.0010_r8
       !
       ! ---------------------------------------------------------------------
       ! * * * update npp, and pool losses  * * *
       ! ---------------------------------------------------------------------
       !
       ! go through all the pfts
       !
       DO j = 1, npft
          !
          ! apply this year's existence arrays to npp
          !
          aynpp(i,j)  = exist(i,j) * aynpp(i,j)
          !
          ! determine above-ground npp for each plant type
          !
          ayanpp(i,j) = (aleaf(j) + awood(j)) * aynpp(i,j)
          !
          ! determine turnover rates for woody biomass:
          !
          ! if pft can exist,    then tauwood = tauwood0 (normal turnover),
          ! if pft cannot exist, then tauwood = taufin years (to kill off trees)
          !
          !          taufin     = 5.00_r8
          taufin     = tauwood0(j)/2.00_r8
          !
          tauwood(i,j) = tauwood0(j) - (tauwood0(j) - taufin) *  &
               (1.00_r8 - exist(i,j))
          !
          ! assume a constant fine root turnover time
          !
          !          tauroot(j) = 1.00_r8
          !
          ! determine litter fall rates
          !
          falll(i) = falll(i) + cbiol(i,j) / tauleaf(j)
          fallr(i) = fallr(i) + cbior(i,j) / tauroot(j)
          fallw(i) = fallw(i) + cbiow(i,j) / tauwood(i,j)
          !
          ! ---------------------------------------------------------------------
          ! * * * update biomass pools  * * *
          ! ---------------------------------------------------------------------
          !
          ! update carbon reservoirs using an analytical solution
          ! to the original carbon balance differential equation
          !
          cbiol(i,j) = cbiol(i,j) * exp( -1.0_r8/tauleaf(j) ) + &
               aleaf(j) * tauleaf(j)   *  max (0.0_r8, aynpp(i,j)) *  (1.0_r8 - exp(-1.0_r8/tauleaf(j)))
          !
          cbiow(i,j) = cbiow(i,j) * exp(-1.0_r8/tauwood(i,j)) + awood(j) * tauwood(i,j) * max (0.0_r8, aynpp(i,j)) *  &
               (1.0_r8 - exp(-1.0_r8/tauwood(i,j)))
          !
          cbior(i,j) = cbior(i,j) * exp( -1.0_r8/tauroot(j) ) + aroot(j) * tauroot(j)   * max (0.0_r8, aynpp(i,j)) *   &
               (1.0_r8 - exp(-1.0_r8/tauroot(j)))
          !
          IF (j.le.8) wood = wood + max (0.00_r8, cbiow(i,j))
          !
       END DO
       !
       ! ---------------------------------------------------------------------
       ! * * * apply disturbances * * *
       ! ---------------------------------------------------------------------
       !
       ! set fixed disturbance regime
       !
       !        disturbf(i) = 0.0050_r8
       !        disturbo(i) = 0.0050_r8
       !        
       !**** DTP 2000/08/10. One can do a decent test of ACME by setting
       !*    these disturbance rates to zero. With these values, the area
       !*    disturbed each year will be zero so the distribution of biomass
       !*    and PFTs across the domain should be identical to those 
       !*    resulting from a run of standard IBIS with zero disturbance

       !*        disturbf(i) = 0.0  ! Test with zero disturbance rate 
       !*        disturbo(i) = 0.0  ! (This should equal reference sim).
       !
       !
       ! call fire disturbance routine
       !**** DTP 2001/03/06: In general isimfire should be set to zero if isim_ac
       !*    is set to 1 (but what does isimfire REAL(KIND=r8)ly do?)
       !
       !        IF (isimfire.eq.1) THEN
       !  CALL fire(npoi      , & ! INTENT(IN)
       !            firefac   , & ! INTENT(IN)
       !    totlit    , & ! INTENT(IN)
       !    disturbf    ) ! INTENT(OUT  )
       !        END IF

       DO j = 1, npft 
          !
          ! calculate biomass (vegetations) carbon lost to atmosphere   
          ! used to balance net ecosystem exchange  
          !
          ! ---------------------------------------------------------------------
          !**** DTP 2000/04/22 QUESTION: 
          ! ---------------------------------------------------------------------
          !* Shouldn't a portion of the destroyed material be added to litter fall?

          cdisturb(i) = cdisturb(i) +  &
               cbiol(i,j) * (disturbf(i) + disturbo(i)) + &
               cbiow(i,j) * (disturbf(i) + disturbo(i)) + &
               cbior(i,j) * (disturbf(i) + disturbo(i))                  
          !          
          ! adjust biomass pools due to disturbances
          !
          cbiol(i,j) = cbiol(i,j) * (1.0_r8 - disturbf(i) - disturbo(i))
          cbiow(i,j) = cbiow(i,j) * (1.0_r8 - disturbf(i) - disturbo(i))
          cbior(i,j) = cbior(i,j) * (1.0_r8 - disturbf(i) - disturbo(i))
          !
          ! constrain biomass fields to be positive
          !
          cbiol(i,j) = max (0.00_r8, cbiol(i,j))
          cbiow(i,j) = max (0.00_r8, cbiow(i,j))
          cbior(i,j) = max (0.00_r8, cbior(i,j))

       END DO



       ! ---------------------------------------------------------------------
       ! * * * check and update biomass pools following disturbance * * *
       ! ---------------------------------------------------------------------
       !
       DO  j = 1, npft
          !
          ! maintain minimum value of leaf carbon in areas that plants exist
          !
          !          xminlai = 0.010
          !
          !
          ! initialize specific leaf area values specific leaf area (m**2/kg) 
          !
          !      data specla  / 25.0,  ! tropical broadleaf evergreen trees
          !     >               25.0,  ! tropical broadleaf drought-deciduous trees
          !     >               25.0,  ! warm-temperate broadleaf evergreen trees
          !     >               12.5,  ! temperate conifer evergreen trees
          !     >               25.0,  ! temperate broadleaf cold-deciduous trees
          !     >               12.5,  ! boreal conifer evergreen trees
          !     >               25.0,  ! boreal broadleaf cold-deciduous trees  
          !     >               25.0,  ! boreal conifer cold-deciduous trees
          !     >               12.5,  ! evergreen shrubs 
          !     >               25.0,  ! deciduous shrubs 
          !     >               20.0,  ! warm (c4) grasses
          !     >               20.0 / ! cool (c3) grasses
          !

          cbiol(i,j) = max (exist(i,j) * xminlai / specla(j), cbiol(i,j))
          !
          ! update vegetation's physical characteristics
          !
          plai(i,j)    = cbiol(i,j) * specla(j)
          biomass(i,j) = cbiol(i,j) + cbiow(i,j) + cbior(i,j)
          !
       END DO
       !
       ! ---------------------------------------------------------------------
       ! * * * update annual npp, lai, and biomass * * *
       ! ---------------------------------------------------------------------
       !
       ! adjust annual net ecosystem exchange (calculated in stats.f) 
       ! by loss of carbon to atmosphere due to biomass burning (fire)
       !
       ayneetot(i) = ayneetot(i) - cdisturb(i)
       !
       ! determine total ecosystem above-ground npp
       !     
       ayanpptot(i) = ayanpp(i,1)  + ayanpp(i,2) +  &
            ayanpp(i,3)  + ayanpp(i,4) +  &
            ayanpp(i,5)  + ayanpp(i,6) +  &
            ayanpp(i,7)  + ayanpp(i,8) +  &
            ayanpp(i,9)  + ayanpp(i,10) + &
            ayanpp(i,11) + ayanpp(i,12)
       !
       ! update total canopy leaf area
       !
       totlaiu(i) = plai(i,1)  + plai(i,2) +   &
            plai(i,3)  + plai(i,4) +   &
            plai(i,5)  + plai(i,6) +   &
            plai(i,7)  + plai(i,8)
       !
       totlail(i) = plai(i,9)  + plai(i,10) +  &
            plai(i,11) + plai(i,12)
       !
       ! update total biomass
       !
       totbiou(i) = biomass(i,1) +  &
            biomass(i,2) +  &
            biomass(i,3) +  &
            biomass(i,4) +  &
            biomass(i,5) +  &
            biomass(i,6) +  &
            biomass(i,7) +  &
            biomass(i,8)
       !
       totbiol(i) = biomass(i,9)  +  &
            biomass(i,10) +  &
            biomass(i,11) +  &
            biomass(i,12)
       !
       ! ---------------------------------------------------------------------
       ! * * * update fractional cover and vegetation height parameters * * *
       ! ---------------------------------------------------------------------
       !
       !

       !**** Added these in temporarily for comparison with original code.
       !**** Delete these from production version....
       !
              fu(i) = (1.00_r8 - exp(-wood)) / (1.00_r8 - exp(-woodnorm))
              fu(i) = fu(i) * (1.0_r8 - disturbf(i) - disturbo(i))
       !IF((totlaiu(i,j)) <= 0.0_r8)THEN
       !    fu(i,j) = 0.0_r8
       !ELSE
       !    fu(i,j) = (1.0_r8 - exp(-totlaiu(i,j))) / (1.0_r8 - exp(-(totlaiu(i,j)+totlail(i,j))))
       !END IF
       !IF((totlail(i,j)) <= 0.0_r8)THEN
       !    fl(i,j) = 0.0_r8
       !ELSE
       !    fl(i,j) = (1.0_r8 - exp(-totlail(i,j))) / (1.0_r8 - exp(-(totlaiu(i,j)+totlail(i,j))))
       !END IF
       !
       ! constrain the fractional cover (upper canopy)
       !
       fu(i) = max (0.250_r8, min (0.9750_r8, fu(i)))
       !
       ! update fractional cover of herbaceous (lower) canopy:
       ! 
       fl(i) = totlail(i) / 1.00_r8
       !
       ! apply disturbances to fractional cover (lower canopy)
       !
       fl(i) = fl(i) * (1.0_r8 - disturbf(i) - disturbo(i))
       !
       ! constrain the fractional cover (lower canopy)
       !
       fl(i) = max (0.250_r8, min (0.9750_r8, fl(i)))
       !
       !
       ! annual update upper canopy height parameters
       ! should be calculated based on vegetative fraction and not the
       ! average over the entire grid cell
       !
       zbot(i,2) = 3.00_r8
       ztop(i,2) = max(zbot(i,2) + 1.000_r8, 2.500_r8 * totbiou(i) / fu(i) * 0.750_r8)
       !
       ! ---------------------------------------------------------------------
       ! * * * update stem area index and sapwood fraction * * *
       ! ---------------------------------------------------------------------
       !
       ! estimate stem area index (sai) as a fraction of the lai
       !
       sai(i,1) = 0.0500_r8 * totlail(i)
       sai(i,2) = 0.2500_r8 * totlaiu(i)
       !
       ! estimate sapwood fraction of woody biomass
       !
       sapspeed  = 25.00_r8                        ! (m/day)
       trans     = 0.00250_r8                      ! (2.5 mm/day) 
       saparea   = (trans / sapspeed)          ! m**2
       !
       sapvolume = saparea * ztop(i,2) * 0.750_r8  ! m**3
       !
       denswood  = 400.00_r8                       ! kg/m**3
       !
       sapfrac(i) = min (0.500_r8, max (0.050_r8, sapvolume * denswood / wood))
       !
    END DO ! DO 100 i = 1, npoi
    !
    ! ---------------------------------------------------------------------
    ! * * * map out vegetation classes for this year * * *
    ! ---------------------------------------------------------------------
    !
    CALL vegmap(totlaiu , &! INTENT(IN   )
         plai    , &! INTENT(IN   )
         totlail , &! INTENT(IN   )
         vegtype0, &! INTENT(OUT  )
         gdd5    , &! INTENT(IN   )
         gdd0    , &! INTENT(IN   )
         npoi    , &! INTENT(IN   )
         npft      )! INTENT(IN   )
    !
    !
    ! return to the main program
    !
    RETURN
  END SUBROUTINE dynaveg1 ! DYNAVEG


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      SUBROUTINE DailyDynaVeg(isimfire ,&
                              npoi     ,&!  INTEGER      , INTENT(IN   ) :: npoi                 ! total number of land points
                              npft     ,&!  INTEGER      , INTENT(IN   ) :: npft                 ! number of plant functional types
                              woodnorm ,&!  REAL(KIND=r8), INTENT(IN   ) :: woodnorm                ! value of woody biomass for upper canopy closure! (ie when wood = woodnorm fu = 1.0) (kg_C m-2)
                              xminlai  ,&!  REAL(KIND=r8), INTENT(IN   ) :: xminlai                ! Minimum LAI for each existing PFT
                              specla   ,&!  REAL(KIND=r8), INTENT(IN   ) :: specla    (npft)         ! specific leaf area (m**2/kg) 
                              aleaf    ,&!  REAL(KIND=r8), INTENT(IN   ) :: aleaf     (npft)         ! carbon allocation fraction to leaves
                              awood    ,&!  REAL(KIND=r8), INTENT(IN   ) :: awood     (npft)         ! carbon allocation fraction to wood 
                              tauwood0 ,&!  REAL(KIND=r8), INTENT(IN   ) :: tauwood0  (npft)         ! normal (unstressed) turnover time for wood biomass (years)
                              tauleaf  ,&!  REAL(KIND=r8), INTENT(IN   ) :: tauleaf   (npft)         ! foliar biomass turnover time constant (years)
                              tauroot  ,&!  REAL(KIND=r8), INTENT(IN   ) :: tauroot   (npft)         ! fine root biomass turnover time constant (years)
                              aroot    ,&!  REAL(KIND=r8), INTENT(IN   ) :: aroot     (npft)         ! carbon allocation fraction to fine roots
                              exist    ,&!  REAL(KIND=r8), INTENT(IN   ) :: exist     (npoi,npft)  ! probability of existence of each plant functional type in a gridcell
                              adco2mic ,&!  REAL(KIND=r8), INTENT(IN   ) :: adco2mic  (npoi)         ! global! daily accumulated co2 respiration from microbes (kg_C m-2 /day)
                              adnpp    ,&!  REAL(KIND=r8), INTENT(INOUT) :: adnpp     (npoi,npft)  ! annual total npp for each plant type(kg-c/m**2/yr)
                              adcbiol  ,&!  REAL(KIND=r8), INTENT(INOUT) :: adcbiol   (npoi,npft)  ! carbon in leaf biomass pool (kg_C m-2)
                              adcbior  ,&!  REAL(KIND=r8), INTENT(INOUT) :: adcbior   (npoi,npft)  ! carbon in fine root biomass pool (kg_C m-2)
                              adcbiow  ,&!  REAL(KIND=r8), INTENT(INOUT) :: adcbiow   (npoi,npft)  ! carbon in woody biomass pool (kg_C m-2)
                              adplai   ,&!  REAL(KIND=r8), INTENT(INOUT) :: adplai    (npoi,npft)  ! global  ! total leaf area index of each plant functional type
                              adfalll  ,&!  REAL(KIND=r8), INTENT(OUT  ) :: adfalll   (npoi)         ! global ! annual leaf litter fall (kg_C m-2/year)
                              adfallr  ,&!  REAL(KIND=r8), INTENT(OUT  ) :: adfallr   (npoi)         ! global ! annual root litter input(kg_C m-2/year)
                              adfallw  ,&!  REAL(KIND=r8), INTENT(OUT  ) :: adfallw   (npoi)         ! global ! annual wood litter fall (kg_C m-2/year)
                              fu       ,&!  REAL(KIND=r8), INTENT(OUT  ) :: fu        (npoi) ! fraction of overall area covered by upper canopy
                              fl       ,&!  REAL(KIND=r8), INTENT(OUT  ) :: fl        (npoi) ! fraction of snow-free area covered by lower  canopy
                              zbot     ,&!  REAL(KIND=r8), INTENT(OUT  ) :: zbot      (npoi,2) ! height of lowest branches above ground (m)
                              ztop     ,&!  REAL(KIND=r8), INTENT(OUT  ) :: ztop      (npoi,2) ! height of plant top above ground (m)
                              sai      ,&!  REAL(KIND=r8), INTENT(OUT  ) :: sai       (npoi,2) ! current single-sided stem area index
                              sapfrac  ,&!  REAL(KIND=r8), INTENT(OUT  ) :: sapfrac   (npoi) ! fraction of woody biomass that is in sapwood
                              vegtype0 ,&!REAL(KIND=r8), INTENT(INOUT) :: vegtype0 (npoi)      ! annual vegetation type - ibis classification
                              totlit   ,&
                              firefac  ,&
                              totlaiu  ,&
                              totlail  ,&
                              totbiou  ,&
                              totbiol   )
                               !  REAL(KIND=r8), INTENT(OUT  ) :: cdisturb  (npoi)         ! annual amount of vegetation carbon lost to atmosphere due to fire  (biomass burning) (kg_C m-2/year)

      IMPLICIT NONE
      INTEGER, INTENT(IN   ) :: isimfire  ! fire switch
                                          ! isim_ac   ! age-class dynamics switch
                                          ! year      ! year of simulation
      INTEGER      , INTENT(IN   ) :: npoi                   ! total number of land points
      INTEGER      , INTENT(IN   ) :: npft                   ! number of plant functional types
      REAL(KIND=r8), INTENT(IN   ) :: woodnorm               ! value of woody biomass for upper canopy closure! (ie when wood = woodnorm fu = 1.0) (kg_C m-2)
      REAL(KIND=r8), INTENT(IN   ) :: xminlai                ! Minimum LAI for each existing PFT
      REAL(KIND=r8), INTENT(IN   ) :: specla    (npft)       ! specific leaf area (m**2/kg) 
      REAL(KIND=r8), INTENT(IN   ) :: aleaf     (npft)       ! carbon allocation fraction to leaves
      REAL(KIND=r8), INTENT(IN   ) :: awood     (npft)       ! carbon allocation fraction to wood 
      REAL(KIND=r8), INTENT(IN   ) :: tauwood0  (npft)       ! normal (unstressed) turnover time for wood biomass (years)
      REAL(KIND=r8), INTENT(IN   ) :: tauleaf   (npft)       ! foliar biomass turnover time constant (years)
      REAL(KIND=r8), INTENT(IN   ) :: tauroot   (npft)       ! fine root biomass turnover time constant (years)
      REAL(KIND=r8), INTENT(IN   ) :: aroot     (npft)       ! carbon allocation fraction to fine roots
      REAL(KIND=r8), INTENT(IN   ) :: exist     (npoi,npft)  ! probability of existence of each plant functional type in a gridcell
      REAL(KIND=r8), INTENT(IN   ) :: adco2mic  (npoi)       ! global! daily accumulated co2 respiration from microbes (kg_C m-2 /day)
      REAL(KIND=r8), INTENT(INOUT) :: adnpp     (npoi,npft)  ! annual total npp for each plant type(kg-c/m**2/yr)
      REAL(KIND=r8), INTENT(INOUT) :: adcbiol   (npoi,npft)  ! carbon in leaf biomass pool (kg_C m-2)
      REAL(KIND=r8), INTENT(INOUT) :: adcbior   (npoi,npft)  ! carbon in fine root biomass pool (kg_C m-2)
      REAL(KIND=r8), INTENT(INOUT) :: adcbiow   (npoi,npft)  ! carbon in woody biomass pool (kg_C m-2)
      REAL(KIND=r8), INTENT(INOUT) :: adplai    (npoi,npft)  ! global  ! total leaf area index of each plant functional type
      REAL(KIND=r8), INTENT(OUT  ) :: adfalll   (npoi)       ! global ! annual leaf litter fall (kg_C m-2/year)
      REAL(KIND=r8), INTENT(OUT  ) :: adfallr   (npoi)       ! global ! annual root litter input(kg_C m-2/year)
      REAL(KIND=r8), INTENT(OUT  ) :: adfallw   (npoi)       ! global ! annual wood litter fall (kg_C m-2/year)
      REAL(KIND=r8), INTENT(OUT  ) :: fu        (npoi)       ! fraction of overall area covered by upper canopy
      REAL(KIND=r8), INTENT(OUT  ) :: fl        (npoi)       ! fraction of snow-free area covered by lower  canopy
      REAL(KIND=r8), INTENT(INOUT) :: zbot      (npoi,2)     ! height of lowest branches above ground (m)
      REAL(KIND=r8), INTENT(INOUT) :: ztop      (npoi,2)     ! height of plant top above ground (m)
      REAL(KIND=r8), INTENT(OUT  ) :: sai       (npoi,2)     ! current single-sided stem area index
      REAL(KIND=r8), INTENT(OUT  ) :: sapfrac   (npoi)       ! fraction of woody biomass that is in sapwood
      REAL(KIND=r8), INTENT(INOUT) :: vegtype0 (npoi)      ! annual vegetation type - ibis classification      
      REAL(KIND=r8), INTENT(IN   ) :: totlit   (npoi)         ! total carbon in all litter pools (kg_C m-2)
      REAL(KIND=r8), INTENT(IN   ) :: firefac  (npoi)         ! factor that respresents the annual average fuel
                                                              ! dryness of a grid cell, and hence characterizes the readiness to burn
      REAL(KIND=r8), INTENT(IN   ) :: totlaiu	(npoi)       ! total leaf area index for the upper canopy
      REAL(KIND=r8), INTENT(IN   ) :: totlail	(npoi)       ! total leaf area index for the lower canopy
      REAL(KIND=r8), INTENT(IN   ) :: totbiou	(npoi)       ! total biomass in the upper canopy (kg_C m-2)
      REAL(KIND=r8), INTENT(IN   ) :: totbiol	(npoi)       ! total biomass in the lower canopy (kg_C m-2)
 
      REAL(KIND=r8) :: cdisturb  (npoi)       ! annual amount of vegetation carbon lost to atmosphere due to fire  (biomass burning) (kg_C m-2/year)

!
! local variables
!
      REAL(KIND=r8)    :: adneetot  (npoi)        ! annual total NEE for ecosystem (kg-C/m**2/yr)
      REAL(KIND=r8)    :: adnpptot  (npoi)
      REAL(KIND=r8)    :: totlaiu_local   (npoi)       ! total leaf area index for the upper canopy
      REAL(KIND=r8)    :: totlail_local   (npoi)       ! total leaf area index for the lower canopy
      REAL(KIND=r8)    :: totbiou_local   (npoi)       ! total biomass in the upper canopy (kg_C m-2)
      REAL(KIND=r8)    :: totbiol_local   (npoi)       ! total biomass in the lower canopy (kg_C m-2)
       
      REAL(KIND=r8)    :: biomass (npoi,npft)    ! total biomass of each plant functional type  (kg_C m-2)
      REAL(KIND=r8)    :: tauwood (npoi,npft)    ! wood biomass turnover time constant (years)
      REAL(KIND=r8)    :: adanpp  (npoi,npft)    ! annual above-ground npp for each plant type(kg-c/m**2/yr)
      REAL(KIND=r8)    :: disturbf(npoi)         ! annual fire disturbance regime (m2/m2/yr)
      REAL(KIND=r8)    :: disturbo(npoi)         ! fraction of biomass pool lost every year to disturbances other than fire
      REAL(KIND=r8)    :: caccount(npoi)     
      REAL(KIND=r8)    :: cbiolmin(npoi,npft) ! minimum leaf biomass used as seed.
      REAL(KIND=r8)    :: wood           ! total amount of woody biomass in gridcell
      REAL(KIND=r8)    :: taufin        !
      REAL(KIND=r8)    :: seedbio
      REAL(KIND=r8)    :: sapspeed      ! in mm/day
      REAL(KIND=r8)    :: trans         ! (2.5 mm/day) 
      REAL(KIND=r8)    :: saparea       ! in m**2
      REAL(KIND=r8)    :: sapvolume     ! in m**3
      REAL(KIND=r8)    :: denswood      ! kg/m**3

      REAL(KIND=r8)    :: rwork
      INTEGER          :: niter,nit
      INTEGER :: i     ! gridcell counter      
      INTEGER :: j     ! gridcell counter

!
! iteration
!
      niter = 10
      rwork = 1.0_r8 / float(niter)

!
! set fixed disturbance regime
!      
      DO i = 1, npoi
         disturbf(i) = 0.0100_r8 / 365.0_r8
         disturbo(i) = 0.0100_r8 / 365.0_r8
      END DO

!
! call fire disturbance routine
!**** DTP 2001/03/06: In general isimfire should be set to zero if isim_ac
!*    is set to 1 (but what does isimfire REAL(KIND=r8)ly do?)
!
!      IF (isimfire.eq.1) THEN
!           CALL fire(npoi      , & ! INTENT(IN   )
!                     firefac   , & ! INTENT(IN   )
!                     totlit    , & ! INTENT(IN   )
!                     disturbf    ) ! INTENT(OUT  )
!         DO i = 1, npoi
!            disturbf(i) = disturbf(i) / 365.0_r8
!         END DO
!      END IF
!
! begin global grid
!
      DO i = 1, npoi
! 
! initialize wood for gridcell
!
         wood = 0.001_r8
!
! ---------------------------------------------------------------------
! * * * initialize vegetation dynamics pools * * *
! ---------------------------------------------------------------------
!
! zero out litter fall fields
!
        adfalll(i) = 0.00_r8
        adfallr(i) = 0.00_r8
        adfallw(i) = 0.00_r8
!
	caccount(i) = 0.00_r8    
!
! zero out carbon lost due to disturbance
! 
        cdisturb(i) = 0.00_r8        
!	cdistinit = 0.0_r8
!
! ---------------------------------------------------------------------
! * * * apply disturbances * * *
! ---------------------------------------------------------------------
!
! set fixed disturbance regime
!
!        disturbf(i) = 0.0050_r8
!        disturbo(i) = 0.0050_r8
!        
!**** DTP 2000/08/10. One can do a decent test of ACME by setting
!*    these disturbance rates to zero. With these values, the area
!*    disturbed each year will be zero so the distribution of biomass
!*    and PFTs across the domain should be identical to those 
!*    resulting from a run of standard IBIS with zero disturbance

!*        disturbf(i) = 0.0  ! Test with zero disturbance rate 
!*        disturbo(i) = 0.0  ! (This should equal reference sim).

!
!
! ---------------------------------------------------------------------
! * * * update npp, and pool losses  * * *
! ---------------------------------------------------------------------
!
! initialize specific leaf area values specific leaf area (m**2/kg) 
!
!      data specla  / 25.0,  ! tropical broadleaf evergreen trees
!     >               25.0,  ! tropical broadleaf drought-deciduous trees
!     >               25.0,  ! warm-temperate broadleaf evergreen trees
!     >               12.5,  ! temperate conifer evergreen trees
!     >               25.0,  ! temperate broadleaf cold-deciduous trees
!     >               12.5,  ! boREAL(KIND=r8) conifer evergreen trees
!     >               25.0,  ! boreal broadleaf cold-deciduous trees  
!     >               25.0,  ! boreal conifer cold-deciduous trees
!     >               12.5,  ! evergreen shrubs 
!     >               25.0,  ! deciduous shrubs 
!     >               20.0,  ! warm (c4) grasses
!     >               20.0 / ! cool (c3) grasses
!
!      woodnorm = 7.5
!

!
! go through all the pfts
!
        DO j = 1, npft
!
! maintain minimum value of leaf carbon in areas where plants exist
!
          cbiolmin(i,j) = exist(i,j)*xminlai/specla(j)
!
! apply this year's existence arrays to npp
!
          adnpp(i,j)  = exist(i,j) * adnpp(i,j)
!
! determine above-ground npp for each plant type
!
          adanpp(i,j) = (aleaf(j) + awood(j)) * adnpp(i,j)

! determine turnover rates for woody biomass:
!
! if pft can exist,    then tauwood = tauwood0 (normal turnover),
! if pft cannot exist, then tauwood = taufin years (to kill off trees)
!
!          taufin     = 5.00_r8
           taufin     = (tauwood0(j))/2.00_r8
!
          tauwood(i,j) = (tauwood0(j)) - ((tauwood0(j)) - taufin) * (1.00_r8 - exist(i,j))
!
! assume a constant fine root turnover time
!
!          tauroot(j) = 1.00_r8
! calculate carbon lost to atmosphere by disturbance (non iterated) :
! corresponds to value calculated by sumnow and used in the instantaneous nee
! used to balance carbon
!
!          cdistinit = cdistinit +  &
!                        ((cbiol(i,j) - cbiolmin(i,j)) * &
!                                      (disturbf(i) + disturbo(i)) + &
!                         cbiow(i,j) * (disturbf(i) + disturbo(i)) + &
!                         cbior(i,j) * (disturbf(i) + disturbo(i))) 
!
! iteration loop
!
          DO nit = 1, niter
          
!
! determine litter fall rates
!
!          falll(i) = falll(i) + cbiol(i,j) / tauleaf(j)
!          fallr(i) = fallr(i) + cbior(i,j) / tauroot(j)
!          fallw(i) = fallw(i) + cbiow(i,j) / tauwood(j)
           adfalll (i) = adfalll(i) + (adcbiol(i,j) - cbiolmin(i,j)) / (tauleaf(j  )*365.0) * rwork
           adfallr (i) = adfallr(i) + (adcbior(i,j)                ) / (tauroot(j  )*365.0) * rwork
           adfallw (i) = adfallw(i) + (adcbiow(i,j)                ) / (tauwood(i,j)*365.0) * rwork
!
! ---------------------------------------------------------------------
! * * * apply disturbances * * *
! ---------------------------------------------------------------------
!
! calculate biomass (vegetations) carbon lost to atmosphere   
! used to balance net ecosystem exchange  
!
          cdisturb(i) = cdisturb(i) + ((adcbiol(i,j) - cbiolmin(i,j)) * (disturbf(i) + disturbo(i)) +   &
                                                         adcbiow(i,j) * (disturbf(i) + disturbo(i)) +   &
                                                         adcbior(i,j) * (disturbf(i) + disturbo(i))) * rwork

!
!
! ---------------------------------------------------------------------
! * * * update biomass pools  * * *
! ---------------------------------------------------------------------
!
! update carbon reservoirs using an analytical solution
! to the original carbon balance differential equation
!
          adcbiol(i,j) = adcbiol(i,j) + (aleaf(j) * max (0.0_r8, adnpp(i,j)) - &
                        (adcbiol(i,j) - cbiolmin(i,j)) / (tauleaf(j)*365.0) -       &
                        (disturbf(i ) + disturbo(i)) *                    &
                        (adcbiol(i,j) - cbiolmin(i,j))) * rwork 
!
          adcbiow(i,j) = adcbiow(i,j) + (awood(j) * max (0.0_r8, adnpp(i,j)) -  &
                         adcbiow(i,j) / (tauwood(i,j)*365.0) -                        &
                        (disturbf(i) + disturbo(i)) * adcbiow(i,j)) * rwork
!
          adcbior(i,j) = adcbior(i,j) + (aroot(j) * max (0.0_r8, adnpp(i,j)) - &
                         adcbior(i,j) / (tauroot(j)*365.0) -                       &
                       (disturbf(i) + disturbo(i)) * adcbior(i,j)) * rwork
!
          END DO
!
! end of iteration loop
!
!
          IF (j.le.8) wood = wood + max (0.00_r8, adcbiow(i,j))

!
          seedbio = max(0.0_r8,(cbiolmin(i,j) - adcbiol(i,j)))
! account for negative biomass in nee (caccount > 0: carbon that has been accounted for 
! as absorbed in the fluxes but that is not accounted for in the calculation of biomass 
! pools ==> has to be released to atmosphere)
!
          caccount(i) = caccount(i) + seedbio -          & 
                     min (0.0_r8, adcbiow(i,j)) - &
                     min (0.0_r8, adcbior(i,j)) 

!
! constrain biomass fields to be positive
!
          adcbiol(i,j) = max (cbiolmin(i,j), adcbiol(i,j))
          adcbiow(i,j) = max (0.0_r8, adcbiow(i,j))
          adcbior(i,j) = max (0.0_r8, adcbior(i,j))
!
! update vegetation's physical characteristics
!
!
! initialize specific leaf area values specific leaf area (m**2/kg) 
!
!      data specla  / 25.0,  ! tropical broadleaf evergreen trees
!     >               25.0,  ! tropical broadleaf drought-deciduous trees
!     >               25.0,  ! warm-temperate broadleaf evergreen trees
!     >               12.5,  ! temperate conifer evergreen trees
!     >               25.0,  ! temperate broadleaf cold-deciduous trees
!     >               12.5,  ! boreal conifer evergreen trees
!     >               25.0,  ! boreal broadleaf cold-deciduous trees  
!     >               25.0,  ! boreal conifer cold-deciduous trees
!     >               12.5,  ! evergreen shrubs 
!     >               25.0,  ! deciduous shrubs 
!     >               20.0,  ! warm (c4) grasses
!     >               20.0 / ! cool (c3) grasses
!
!   cbiol==> carbon in leaf biomass pool (kg_C m-2)

          adplai (i,j)    = adcbiol(i,j) * specla(j)  !(kg_C m-2) * (m**2/kg) 

          biomass(i,j)    = adcbiol(i,j) + adcbiow(i,j) + adcbior(i,j)
!

        END DO

!
! ---------------------------------------------------------------------
! * * * update annual npp, lai, and biomass * * *
! ---------------------------------------------------------------------
!
! Disturbance can't result in negative biomass. caccount account for the 
! carbon not to be removed by the disturbance.
!
       cdisturb(i) = cdisturb(i) - caccount(i)
!
! determine total ecosystem positive npp (changed by exist at begin
! of subroutine). Different from sum of monthly and daily npp)
!
        adnpptot(i) = max(0.0_r8,adnpp(i,1))  + max(0.0_r8,adnpp(i,2)) + &
                      max(0.0_r8,adnpp(i,3))  + max(0.0_r8,adnpp(i,4)) + &
                      max(0.0_r8,adnpp(i,5))  + max(0.0_r8,adnpp(i,6)) + &
                      max(0.0_r8,adnpp(i,7))  + max(0.0_r8,adnpp(i,8)) + &
                      max(0.0_r8,adnpp(i,9))  + max(0.0_r8,adnpp(i,10)) + &
                      max(0.0_r8,adnpp(i,11)) + max(0.0_r8,adnpp(i,12))
!
! adjust annual net ecosystem exchange (calculated in stats.f) 
! by new value of npp (depending on exist), andloss of carbon to
! atmosphere due to biomass burning (fire)
!
        adneetot(i) = adnpptot(i) - adco2mic(i) - cdisturb(i)
!
! determine total ecosystem above-ground npp
!
!        adanpptot(i) = adanpp(i,1)  + adanpp(i,2) +  &
!                       adanpp(i,3)  + adanpp(i,4) +  &
!                       adanpp(i,5)  + adanpp(i,6) +  &
!                       adanpp(i,7)  + adanpp(i,8) +  &
!                       adanpp(i,9)  + adanpp(i,10) + &
!                       adanpp(i,11) + adanpp(i,12)
!
!
! update total canopy leaf area
!
        totlaiu_local(i) = adplai(i,1)  + adplai(i,2) + &
                     adplai(i,3)  + adplai(i,4) + & 
                     adplai(i,5)  + adplai(i,6) + &
                     adplai(i,7)  + adplai(i,8)
!
        totlail_local(i) = adplai(i,9)  + adplai(i,10) + &
                     adplai(i,11) + adplai(i,12)

!
! update total biomass
!
        totbiou_local(i) = biomass(i,1) +  & 
                     biomass(i,2) +  & 
                     biomass(i,3) +  &
                     biomass(i,4) +  &
                     biomass(i,5) +  &
                     biomass(i,6) +  &
                     biomass(i,7) +  &
                     biomass(i,8)
!
        totbiol_local(i) = biomass(i,9)  + &
                     biomass(i,10) + &
                     biomass(i,11) + &
                     biomass(i,12)
!
! ---------------------------------------------------------------------
! * * * update fractional cover and vegetation height parameters * * *
! ---------------------------------------------------------------------
!
! update fractional cover of forest and herbaceous canopies:
! 
!PK        IF(vegtype0(i) == 15)THEN
           fu(i) = (1.0_r8 - exp(-wood)) / (1.0_r8 - exp(-woodnorm))
           fl(i) = totlail_local(i) / 1.0_r8
!PK           ELSE
!PK              IF((totlaiu_local(i) ) <= 0.0_r8)THEN
!PK                 fu(i) = 0.0_r8
!PK              ELSE
!PK                 fu(i) = (1.0_r8 - exp(-totlaiu_local(i))) / (1.0_r8 - exp(-(totlaiu_local(i)+totlail_local(i))))
!PK               END IF
!PK              IF((totlail_local(i)) <= 0.0_r8)THEN
!PK                 fl(i) = 0.0_r8
!PK              ELSE
!PK                 fl(i) = (1.0_r8 - exp(-totlail_local(i))) / (1.0_r8 - exp(-(totlaiu_local(i)+totlail_local(i))))
!PK               END IF
!PK            END IF
!
        fu(i) = max (0.25_r8, min (0.975_r8, fu(i)))
        fl(i) = max (0.25_r8, min (0.975_r8, fl(i)))

!
! apply disturbances to fractional cover
!
        fu(i) = fu(i) * (1.0_r8 - disturbf(i) - disturbo(i))
        fl(i) = fl(i) * (1.0_r8 - disturbf(i) - disturbo(i))
!
! constrain the fractional cover
!
        fu(i) = max (0.25_r8, min (0.975_r8, fu(i)))
        fl(i) = max (0.25_r8, min (0.975_r8, fl(i)))

!
! annual update upper canopy height parameters
! should be calculated based on vegetative fraction and not the
! average over the entire grid cell
!
        zbot(i,2) = 3.0_r8
!PK   check Agro
 !       ztop(i,2) = max(zbot(i,2) + 1.00_r8, 2.50_r8 * totbiou_local(i) / fu(i) * 0.75_r8)

        ztop(i,2) = max(zbot(i,2) + 1.00_r8, 2.50_r8 * totbiou(i) / fu(i) * 0.75_r8)

!PK
!
! ---------------------------------------------------------------------
! * * * update stem area index and sapwood fraction * * *
! ---------------------------------------------------------------------
!
! estimate stem area index (sai) as a fraction of the lai
!
        sai(i,1) = 0.050_r8 * totlail_local(i)
        sai(i,2) = 0.250_r8 * totlaiu_local(i)
!
! estimate sapwood fraction of woody biomass
!
        sapspeed  = 25.0_r8                        ! (m/day)
        trans     = 0.0025_r8                      ! (2.5 mm/day) 
        saparea   = (trans / sapspeed)          ! m**2
!
        sapvolume = saparea * ztop(i,2) * 0.75_r8  ! m**3
!
        denswood  = 400.0_r8                       ! kg/m**3
!
        sapfrac(i) = min (0.50_r8, max (0.05_r8, sapvolume * denswood / wood))
!
!
! ---------------------------------------------------------------------
! * * * update annual npp, lai, and biomass * * *
! ---------------------------------------------------------------------
      END DO ! DO 100 i = 1, npoi
!
! ---------------------------------------------------------------------
! * * * map out vegetation classes for this year * * *
! ---------------------------------------------------------------------
!

      RETURN
      END SUBROUTINE DailyDynaVeg ! DYNAVEG

  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE dynaveg2(isimfire , &! INTENT(IN   )
       tauwood0 , &! INTENT(IN   )
       tauwood  , &! INTENT(OUT  )
       tauleaf  , &! INTENT(IN   )
       tauroot  , &! INTENT(IN   )
       xminlai  , &! INTENT(IN   )
       falll    , &! INTENT(OUT  )
       fallr    , &! INTENT(OUT  )
       fallw    , &! INTENT(OUT  )
       cdisturb , &! INTENT(OUT  )
       exist    , &! INTENT(IN   )
       aleaf    , &! INTENT(IN   )
       awood    , &! INTENT(IN   )
       cbiol    , &! INTENT(INOUT) global
       cbior    , &! INTENT(INOUT) global
       cbiow    , &! INTENT(INOUT) global
       aroot    , &! INTENT(IN   )
       disturbf , &! INTENT(OUT  )
       disturbo , &! INTENT(OUT  )
       firefac  , &! INTENT(IN   )
       totlit   , &! INTENT(IN   )
       specla   , &! INTENT(IN   )
       plai     , &! INTENT(INOUT) local
       biomass  , &! INTENT(OUT  )
       totlaiu  , &! INTENT(INOUT) local
       totlail  , &! INTENT(INOUT) local
       totbiou  , &! INTENT(INOUT) local
       totbiol  , &! INTENT(OUT  )
       fu   , &! INTENT(OUT  )
       woodnorm , &! INTENT(IN   )
       fl   , &! INTENT(OUT  )
       zbot     , &! INTENT(OUT  )
       ztop     , &! INTENT(OUT  )
       sai      , &! INTENT(OUT  )
       sapfrac  , &! INTENT(OUT  )
       vegtype0 , &! INTENT(OUT  )
       gdd5     , &! INTENT(IN   )
       gdd0     , &! INTENT(IN   )
       aynpp    , &! INTENT(INOUT) global
       ayanpp   , &! INTENT(OUT  )
       ayneetot , &! INTENT(INOUT) global
       ayanpptot, &! INTENT(OUT  )
       aynpptot, &! INTENT(OUT  )
       ayco2mic, &! INTENT(IN  )
       npoi     , &!
       npft       )! , isim_ac, year)
    ! ---------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN   ) :: npoi                   ! total number of land points
    INTEGER, INTENT(IN   ) :: npft                   ! number of plant functional types
    REAL(KIND=r8)   , INTENT(INOUT) :: aynpp    (npoi,npft)   ! annual total npp for each plant type(kg-c/m**2/yr)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ayanpp   (npoi,npft)   ! annual above-ground npp for each plant type(kg-c/m**2/yr)
    REAL(KIND=r8)   , INTENT(INOUT) :: ayneetot (npoi)        ! annual total NEE for ecosystem (kg-C/m**2/yr)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ayanpptot(npoi)        ! annual above-ground npp for ecosystem (kg-c/m**2/yr)
    REAL(KIND=r8)   , INTENT(OUT  ) :: falll   (npoi)         ! annual leaf litter fall                      (kg_C m-2/year)
    REAL(KIND=r8)   , INTENT(OUT  ) :: fallr   (npoi)         ! annual root litter input                     (kg_C m-2/year)
    REAL(KIND=r8)   , INTENT(OUT  ) :: fallw   (npoi)         ! annual wood litter fall                      (kg_C m-2/year)
    REAL(KIND=r8)   , INTENT(OUT  ) :: cdisturb(npoi)         ! annual amount of vegetation carbon lost 
    ! to atmosphere due to fire  (biomass burning) (kg_C m-2/year)
    REAL(KIND=r8)   , INTENT(IN   ) :: exist   (npoi,npft)    ! probability of existence of each plant functional type in a gridcell
    REAL(KIND=r8)   , INTENT(IN   ) :: aleaf   (npft)         ! carbon allocation fraction to leaves
    REAL(KIND=r8)   , INTENT(IN   ) :: awood   (npft)         ! carbon allocation fraction to wood 
    REAL(KIND=r8)   , INTENT(INOUT) :: cbiol   (npoi,npft)    ! carbon in leaf biomass pool (kg_C m-2)
    REAL(KIND=r8)   , INTENT(INOUT) :: cbior   (npoi,npft)    ! carbon in fine root biomass pool (kg_C m-2)
    REAL(KIND=r8)   , INTENT(INOUT) :: cbiow   (npoi,npft)    ! carbon in woody biomass pool (kg_C m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: aroot   (npft)         ! carbon allocation fraction to fine roots
    REAL(KIND=r8)   , INTENT(OUT  ) :: disturbf(npoi)         ! annual fire disturbance regime (m2/m2/yr)
    REAL(KIND=r8)   , INTENT(OUT  ) :: disturbo(npoi)         ! fraction of biomass pool lost every year to disturbances other than fire
    REAL(KIND=r8)   , INTENT(IN   ) :: firefac (npoi)         ! factor that respresents the annual average fuel
    ! dryness of a grid cell, and hence characterizes the readiness to burn
    REAL(KIND=r8)   , INTENT(IN   ) :: totlit  (npoi)         ! total carbon in all litter pools (kg_C m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: specla  (npft)         ! specific leaf area (m**2/kg) 
    REAL(KIND=r8)   , INTENT(INOUT) :: plai    (npoi,npft)    ! total leaf area index of each plant functional type
    REAL(KIND=r8)   , INTENT(OUT  ) :: biomass (npoi,npft)    ! total biomass of each plant functional type  (kg_C m-2)
    REAL(KIND=r8)   , INTENT(INOUT) :: totlaiu (npoi)         ! total leaf area index for the upper canopy
    REAL(KIND=r8)   , INTENT(INOUT) :: totlail (npoi)         ! total leaf area index for the lower canopy
    REAL(KIND=r8)   , INTENT(INOUT) :: totbiou (npoi)         ! total biomass in the upper canopy (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: totbiol (npoi)         ! total biomass in the lower canopy (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: fu      (npoi)         ! fraction of overall area covered by upper canopy
    REAL(KIND=r8)   , INTENT(IN   ) :: woodnorm       ! value of woody biomass for upper canopy closure
    ! (ie when wood = woodnorm fu = 1.0) (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: fl      (npoi)         ! fraction of snow-free area covered by lower  canopy
    REAL(KIND=r8)   , INTENT(INOUT) :: zbot    (npoi,2)       ! height of lowest branches above ground (m)
    REAL(KIND=r8)   , INTENT(INOUT) :: ztop    (npoi,2)       ! height of plant top above ground (m)
    REAL(KIND=r8)   , INTENT(OUT  ) :: sai     (npoi,2)       ! current single-sided stem area index
    REAL(KIND=r8)   , INTENT(OUT  ) :: sapfrac (npoi)         ! fraction of woody biomass that is in sapwood
    REAL(KIND=r8)   , INTENT(INOUT) :: vegtype0(npoi)         ! annual vegetation type - ibis classification
    REAL(KIND=r8)   , INTENT(IN   ) :: gdd5    (npoi)         ! growing degree days > 5C
    REAL(KIND=r8)   , INTENT(IN   ) :: gdd0    (npoi)         ! growing degree days > 0C 
    REAL(KIND=r8)   , INTENT(IN   ) :: tauwood0(npft)   ! normal (unstressed) turnover time for wood biomass (years)
    REAL(KIND=r8)   , INTENT(OUT  ) :: tauwood (npoi,npft)   ! wood biomass turnover time constant (years)
    REAL(KIND=r8)   , INTENT(IN   ) :: tauleaf (npft)   ! foliar biomass turnover time constant (years)
    REAL(KIND=r8)   , INTENT(IN   ) :: tauroot (npft)   ! fine root biomass turnover time constant (years)
    REAL(KIND=r8)   , INTENT(IN   ) :: xminlai          ! Minimum LAI for each existing PFT
    REAL(KIND=r8)   , INTENT(OUT  ) :: aynpptot (npoi)
    REAL(KIND=r8)   , INTENT(IN   ) :: ayco2mic (npoi)        ! global! annual total CO2 flux from microbial respiration (kg-C/m**2/yr)
    REAL(KIND=r8)  :: caccount(npoi)     
    REAL(KIND=r8)  :: cbiolmin(npoi,npft) ! minimum leaf biomass used as seed.

    !
    ! Arguments
    !
    INTEGER, INTENT(IN   ) :: isimfire  ! fire switch
    ! isim_ac   ! age-class dynamics switch
    ! year      ! year of simulation

    !REAL(KIND=r8)    :: pdist      ! probability of other disturbance types....

    REAL(KIND=r8)    , PARAMETER :: pfire = 1.00_r8 ! probability of fire -- should be determined externally.
    ! for now we just assume it occurs all the time

    !
    ! local variables
    !
    INTEGER :: i     ! gridcell counter      
    INTEGER :: j     ! gridcell counter
    !
    !      REAL(KIND=r8)    :: cdistinit
    REAL(KIND=r8)    :: seedbio
    REAL(KIND=r8)    :: sapspeed      ! in mm/day
    REAL(KIND=r8)    :: trans        ! (2.5 mm/day) 
    REAL(KIND=r8)    :: saparea       ! in m**2
    REAL(KIND=r8)    :: sapvolume     ! in m**3
    REAL(KIND=r8)    :: denswood      ! kg/m**3
    REAL(KIND=r8)    :: wood          ! total amount of woody biomass in gridcell
    REAL(KIND=r8)    :: taufin        !
    INTEGER          :: niter,nit
    REAL(KIND=r8)    :: rwork
    !*      REAL(KIND=r8)    :: xminlai               !
    !
    !*      REAL(KIND=r8)
    !*      REAL(KIND=r8)    :: aleaf(npft),   ! allocation fraction to leaves
    !*      REAL(KIND=r8)    :: aroot(npft),   ! allocation fraction to fine roots
    !*      REAL(KIND=r8)    :: awood(npft),   ! allocation fraction to wood
    !*      REAL(KIND=r8)    :: tauleaf(npft), ! turnover time of carbon in leaves (years)
    !*      REAL(KIND=r8)    :: tauroot(npft), ! turnover time of carbon in fine roots (years)
    !*      REAL(KIND=r8)    :: tauwood(npft)   ! turnover time of carbon in wood (years)
    !*      REAL(KIND=r8)    :: tauwood0(npft) ! normal (unstressed) turnover time
    !
    ! ibis uses a small number of plant functional types:
    !
    !  1: tropical broadleaf evergreen tree
    !  2: tropical broadleaf drought-deciduous trees
    !  3: warm-temperate broadleaf evergreen tree
    !  4: temperate conifer evergreen tree
    !  5: temperate broadleaf cold-deciduous tree
    !  6: boREAL(KIND=r8) conifer evergreen tree
    !  7: boREAL(KIND=r8) broadleaf cold-deciduous tree
    !  8: boREAL(KIND=r8) conifer cold-deciduous tree
    !  9: evergreen shrub
    ! 10: deciduous shrub
    ! 11: warm (c4) grass
    ! 12: cool (c3) grass
    !
    ! ---------------------------------------------------------------------
    ! * * * specify biomass turnover parameters (years) * * *
    ! ---------------------------------------------------------------------
    !
    !      data tauleaf / 1.01,   ! tropical broadleaf evergreen trees
    !     >               1.00,   ! tropical broadleaf drought-deciduous trees
    !     >               1.00,   ! warm-temperate broadleaf evergreen trees
    !     >               2.00,   ! temperate conifer evergreen trees
    !     >               1.00,   ! temperate broadleaf cold-deciduous trees
    !     >               2.50,   ! boREAL(KIND=r8) conifer evergreen trees
    !     >               1.00,   ! boREAL(KIND=r8) broadleaf cold-deciduous trees
    !     >               1.00,   ! boREAL(KIND=r8) conifer cold-deciduous trees
    !     >               1.50,   ! evergreen shrubs
    !     >               1.00,   ! deciduous shrubs
    !     >               1.25,   ! warm (c4) grasses
    !     >               1.50 /  ! cool (c3) grasses
    !
    !      data tauwood0 / 25.0,  ! tropical broadleaf evergreen trees
    !     >                25.0,  ! tropical broadleaf drought-deciduous trees
    !     >                25.0,  ! warm-temperate broadleaf evergreen trees
    !     >                50.0,  ! temperate conifer evergreen trees
    !     >                50.0,  ! temperate broadleaf cold-deciduous trees
    !     >               100.0,  ! boREAL(KIND=r8) conifer evergreen trees
    !     >               100.0,  ! boREAL(KIND=r8) broadleaf cold-deciduous trees
    !     >               100.0,  ! boREAL(KIND=r8) conifer cold-deciduous trees
    !     >                 5.0,  ! evergreen shrubs
    !     >                 5.0,  ! deciduous shrubs
    !     >               999.0,  ! warm (c4) grasses
    !     >               999.0 / ! cool (c3) grasses
    !
    ! iteration
    !
    niter = 10
    rwork = 1.0_r8 / REAL(niter,kind=r8)
    !
    ! set fixed disturbance regime
    !      
    DO i = 1, npoi
       IF(vegtype0(i) /= 15.0_r8)THEN
          disturbf(i) = 0.0050_r8
          disturbo(i) = 0.0050_r8
       END IF
    END DO
    !
    ! call fire disturbance routine
    !**** DTP 2001/03/06: In general isimfire should be set to zero if isim_ac
    !*    is set to 1 (but what does isimfire REAL(KIND=r8)ly do?)
    !
    IF (isimfire.eq.1) THEN
       CALL fire(npoi      , & ! INTENT(IN)
            firefac   , & ! INTENT(IN)
            totlit    , & ! INTENT(IN)
            disturbf  ,&
            vegtype0  ) ! INTENT(OUT  )
    END IF

    !
    ! begin global grid
    !
    DO i = 1, npoi
       IF(vegtype0(i) /= 15.0_r8)THEN
       ! 
       ! initialize wood for gridcell
       !
       wood = 0.001_r8
       !
       ! ---------------------------------------------------------------------
       ! * * * initialize vegetation dynamics pools * * *
       ! ---------------------------------------------------------------------
       !
       ! zero out litter fall fields
       !
       falll(i) = 0.00_r8
       fallr(i) = 0.00_r8
       fallw(i) = 0.00_r8
       !
       caccount(i) = 0.00_r8    
       !
       ! zero out carbon lost due to disturbance
       ! 
       cdisturb(i) = 0.00_r8        
       !cdistinit = 0.0_r8
       !
       ! ---------------------------------------------------------------------
       ! * * * apply disturbances * * *
       ! ---------------------------------------------------------------------
       !
       ! set fixed disturbance regime
       !
       !        disturbf(i) = 0.0050_r8
       !        disturbo(i) = 0.0050_r8
       !        
       !**** DTP 2000/08/10. One can do a decent test of ACME by setting
       !*    these disturbance rates to zero. With these values, the area
       !*    disturbed each year will be zero so the distribution of biomass
       !*    and PFTs across the domain should be identical to those 
       !*    resulting from a run of standard IBIS with zero disturbance

       !*        disturbf(i) = 0.0  ! Test with zero disturbance rate 
       !*        disturbo(i) = 0.0  ! (This should equal reference sim).

       !
       !
       ! ---------------------------------------------------------------------
       ! * * * update npp, and pool losses  * * *
       ! ---------------------------------------------------------------------
       !
       ! initialize specific leaf area values specific leaf area (m**2/kg) 
       !
       !      data specla  / 25.0,  ! tropical broadleaf evergreen trees
       !     >               25.0,  ! tropical broadleaf drought-deciduous trees
       !     >               25.0,  ! warm-temperate broadleaf evergreen trees
       !     >               12.5,  ! temperate conifer evergreen trees
       !     >               25.0,  ! temperate broadleaf cold-deciduous trees
       !     >               12.5,  ! boREAL(KIND=r8) conifer evergreen trees
       !     >               25.0,  ! boreal broadleaf cold-deciduous trees  
       !     >               25.0,  ! boreal conifer cold-deciduous trees
       !     >               12.5,  ! evergreen shrubs 
       !     >               25.0,  ! deciduous shrubs 
       !     >               20.0,  ! warm (c4) grasses
       !     >               20.0 / ! cool (c3) grasses
       !
       !      woodnorm = 7.5
       !
       
       !
       ! go through all the pfts
       !
       DO j = 1, npft
          !
          ! maintain minimum value of leaf carbon in areas where plants exist
          !
          cbiolmin(i,j) = exist(i,j)*xminlai/specla(j)
          !
          ! apply this year's existence arrays to npp
          !
          aynpp(i,j)  = exist(i,j) * aynpp(i,j)
          !
          ! determine above-ground npp for each plant type
          !
          ayanpp(i,j) = (aleaf(j) + awood(j)) * aynpp(i,j)
          !
          ! determine turnover rates for woody biomass:
          !
          ! if pft can exist,    then tauwood = tauwood0 (normal turnover),
          ! if pft cannot exist, then tauwood = taufin years (to kill off trees)
          !
          !          taufin     = 5.00_r8
          taufin     = tauwood0(j)/2.00_r8
          !
          tauwood(i,j) = tauwood0(j) - (tauwood0(j) - taufin) * (1.00_r8 - exist(i,j))
          !
          ! assume a constant fine root turnover time
          !
          !          tauroot(j) = 1.00_r8
          ! calculate carbon lost to atmosphere by disturbance (non iterated) :
          ! corresponds to value calculated by sumnow and used in the instantaneous nee
          ! used to balance carbon
          !
          !          cdistinit = cdistinit +  &
          !                        ((cbiol(i,j) - cbiolmin(i,j)) * &
          !                                      (disturbf(i) + disturbo(i)) + &
          !                         cbiow(i,j) * (disturbf(i) + disturbo(i)) + &
          !                         cbior(i,j) * (disturbf(i) + disturbo(i))) 
          !
          ! iteration loop
          !
          DO nit = 1, niter

             !
             ! determine litter fall rates
             !
             !          falll(i) = falll(i) + cbiol(i,j) / tauleaf(j)
             !          fallr(i) = fallr(i) + cbior(i,j) / tauroot(j)
             !          fallw(i) = fallw(i) + cbiow(i,j) / tauwood(j)
             falll (i) = falll(i) + (cbiol(i,j) - cbiolmin(i,j)) / tauleaf(j) * rwork
             fallr (i) = fallr(i) + cbior(i,j) / tauroot(j) * rwork
             fallw (i) = fallw(i) + cbiow(i,j) / tauwood(i,j) * rwork
             !
             ! ---------------------------------------------------------------------
             ! * * * apply disturbances * * *
             ! ---------------------------------------------------------------------
             !
             ! calculate biomass (vegetations) carbon lost to atmosphere   
             ! used to balance net ecosystem exchange  
             !
             cdisturb(i) = cdisturb(i) + ((cbiol(i,j) - cbiolmin(i,j)) * &
                  (disturbf(i) + disturbo(i)) +   &
                  cbiow(i,j) * (disturbf(i) + disturbo(i)) +   &
                  cbior(i,j) * (disturbf(i) + disturbo(i))) * rwork                  

             !
             ! ---------------------------------------------------------------------
             ! * * * update biomass pools  * * *
             ! ---------------------------------------------------------------------
             !
             ! update carbon reservoirs using an analytical solution
             ! to the original carbon balance differential equation
             !
             cbiol(i,j) = cbiol(i,j) + (aleaf(j) * max (0.0_r8, aynpp(i,j)) - &
                  (cbiol(i,j) - cbiolmin(i,j)) / tauleaf(j) -       &
                  (disturbf(i) + disturbo(i)) *                    &
                  (cbiol(i,j) - cbiolmin(i,j))) * rwork 
             !
             cbiow(i,j) = cbiow(i,j) + (awood(j) * max (0.0_r8, aynpp(i,j)) -  &
                  cbiow(i,j) / tauwood(i,j) -                        &
                  (disturbf(i) + disturbo(i)) * cbiow(i,j)) * rwork
             !
             cbior(i,j) = cbior(i,j) + (aroot(j) * max (0.0_r8, aynpp(i,j)) - &
                  cbior(i,j) / tauroot(j) -                       &
                  (disturbf(i) + disturbo(i)) * cbior(i,j)) * rwork
             !
          END DO
          !
          ! end of iteration loop
          !
          !
          IF (j.le.8) wood = wood + max (0.00_r8, cbiow(i,j))
          !
          seedbio = max(0.0_r8,(cbiolmin(i,j) - cbiol(i,j)))
          ! account for negative biomass in nee (caccount > 0: carbon that has been accounted for 
          ! as absorbed in the fluxes but that is not accounted for in the calculation of biomass 
          ! pools ==> has to be released to atmosphere)
          !
          caccount(i) = caccount(i) + seedbio -          & 
               min (0.0_r8, cbiow(i,j)) - &
               min (0.0_r8, cbior(i,j)) 
          !
          ! constrain biomass fields to be positive
          !
          cbiol(i,j) = max (cbiolmin(i,j), cbiol(i,j))
          cbiow(i,j) = max (0.0_r8, cbiow(i,j))
          cbior(i,j) = max (0.0_r8, cbior(i,j))
          !
          ! update vegetation's physical characteristics
          !
          !
          ! initialize specific leaf area values specific leaf area (m**2/kg) 
          !
          !      data specla  / 25.0,  ! tropical broadleaf evergreen trees
          !     >               25.0,  ! tropical broadleaf drought-deciduous trees
          !     >               25.0,  ! warm-temperate broadleaf evergreen trees
          !     >               12.5,  ! temperate conifer evergreen trees
          !     >               25.0,  ! temperate broadleaf cold-deciduous trees
          !     >               12.5,  ! boreal conifer evergreen trees
          !     >               25.0,  ! boreal broadleaf cold-deciduous trees  
          !     >               25.0,  ! boreal conifer cold-deciduous trees
          !     >               12.5,  ! evergreen shrubs 
          !     >               25.0,  ! deciduous shrubs 
          !     >               20.0,  ! warm (c4) grasses
          !     >               20.0 / ! cool (c3) grasses
          !
          !   cbiol==> carbon in leaf biomass pool (kg_C m-2)

          plai(i,j)    = cbiol(i,j) * specla(j)
          biomass(i,j) = cbiol(i,j) + cbiow(i,j) + cbior(i,j)
          !

       END DO
       !
       ! ---------------------------------------------------------------------
       ! * * * update annual npp, lai, and biomass * * *
       ! ---------------------------------------------------------------------
       !
       ! Disturbance can't result in negative biomass. caccount account for the 
       ! carbon not to be removed by the disturbance.
       !
       cdisturb(i) = cdisturb(i) - caccount(i)
       !
       ! determine total ecosystem positive npp (changed by exist at begin
       ! of subroutine). Different from sum of monthly and daily npp)
       !
       aynpptot(i) = max(0.0_r8,aynpp(i,1))  + max(0.0_r8,aynpp(i,2)) + &
            max(0.0_r8,aynpp(i,3))  + max(0.0_r8,aynpp(i,4)) + &
            max(0.0_r8,aynpp(i,5))  + max(0.0_r8,aynpp(i,6)) + &
            max(0.0_r8,aynpp(i,7))  + max(0.0_r8,aynpp(i,8)) + &
            max(0.0_r8,aynpp(i,9))  + max(0.0_r8,aynpp(i,10)) + &
            max(0.0_r8,aynpp(i,11)) + max(0.0_r8,aynpp(i,12))
       !
       ! adjust annual net ecosystem exchange (calculated in stats.f) 
       ! by new value of npp (depending on exist), andloss of carbon to
       ! atmosphere due to biomass burning (fire)
       !
       ayneetot(i) = aynpptot(i) - ayco2mic(i) - cdisturb(i)
       !
       ! determine total ecosystem above-ground npp
       !
       ayanpptot(i) = ayanpp(i,1)  + ayanpp(i,2) +  &
            ayanpp(i,3)  + ayanpp(i,4) +  &
            ayanpp(i,5)  + ayanpp(i,6) +  &
            ayanpp(i,7)  + ayanpp(i,8) +  &
            ayanpp(i,9)  + ayanpp(i,10) + &
            ayanpp(i,11) + ayanpp(i,12)
       !
       ! update total canopy leaf area
       !
       totlaiu(i) = plai(i,1)  + plai(i,2) + &
            plai(i,3)  + plai(i,4) + & 
            plai(i,5)  + plai(i,6) + &
            plai(i,7)  + plai(i,8)
       !
       totlail(i) = plai(i,9)  + plai(i,10) + &
            plai(i,11) + plai(i,12)
       !
       ! update total biomass
       !
       totbiou(i) = biomass(i,1) +  & 
            biomass(i,2) +  & 
            biomass(i,3) +  &
            biomass(i,4) +  &
            biomass(i,5) +  &
            biomass(i,6) +  &
            biomass(i,7) +  &
            biomass(i,8)
       !
       totbiol(i) = biomass(i,9)  + &
            biomass(i,10) + &
            biomass(i,11) + &
            biomass(i,12)

       !
       ! ---------------------------------------------------------------------
       ! * * * update fractional cover and vegetation height parameters * * *
       ! ---------------------------------------------------------------------
       !
       ! update fractional cover of forest and herbaceous canopies:
       ! 
       !PKfu(i) = (1.0_r8 - exp(-wood)) / (1.0_r8 - exp(-woodnorm))
       !
       !PKfl(i) = totlail(i) / 1.0_r8
       !
       ! apply disturbances to fractional cover
       !
       !PKfu(i) = fu(i) * (1.0_r8 - disturbf(i) - disturbo(i))
       !PKfl(i) = fl(i) * (1.0_r8 - disturbf(i) - disturbo(i))
       !
       ! constrain the fractional cover
       !
       !PKfu(i) = max (0.25_r8, min (0.975_r8, fu(i)))
       !PKfl(i) = max (0.25_r8, min (0.975_r8, fl(i)))
       !
       ! annual update upper canopy height parameters
       ! should be calculated based on vegetative fraction and not the
       ! average over the entire grid cell
       !
       !PKzbot(i,2) = 3.0_r8
       !PKztop(i,2) = max(zbot(i,2) + 1.00_r8, 2.50_r8 * totbiou(i) / fu(i) * 0.75_r8)
       !
       ! ---------------------------------------------------------------------
       ! * * * update stem area index and sapwood fraction * * *
       ! ---------------------------------------------------------------------
       !
       ! estimate stem area index (sai) as a fraction of the lai
       !
       !PKsai(i,1) = 0.050_r8 * totlail(i)
       !PKsai(i,2) = 0.250_r8 * totlaiu(i)
       !
       ! estimate sapwood fraction of woody biomass
       !
       sapspeed  = 25.0_r8                        ! (m/day)
       trans     = 0.0025_r8                      ! (2.5 mm/day) 
       saparea   = (trans / sapspeed)          ! m**2
       !
       sapvolume = saparea * ztop(i,2) * 0.75_r8  ! m**3
       !
       denswood  = 400.0_r8                       ! kg/m**3
       !
       !PKsapfrac(i) = min (0.50_r8, max (0.05_r8, sapvolume * denswood / wood))
       !
       END IF 
    END DO ! DO 100 i = 1, npoi
    !
    ! ---------------------------------------------------------------------
    ! * * * map out vegetation classes for this year * * *
    ! ---------------------------------------------------------------------
    !
    CALL vegmap(totlaiu , &! INTENT(IN   )
         plai    , &! INTENT(IN   )
         totlail , &! INTENT(IN   )
         vegtype0, &! INTENT(OUT  )
         gdd5    , &! INTENT(IN   )
         gdd0    , &! INTENT(IN   )
         npoi    , &! INTENT(IN   )
         npft      )! INTENT(IN   )
    !
    !
    ! return to the main program
    !
    RETURN
  END SUBROUTINE dynaveg2 ! DYNAVEG
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE fire(npoi, &! INTENT(IN   )
       firefac   , &! INTENT(IN   )
       totlit    , &! INTENT(IN   )
       disturbf  , &
       vegtype0  )! INTENT(OUT  )
    ! ---------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN   ) :: npoi                ! total number of land points
    REAL(KIND=r8)   , INTENT(IN   ) :: firefac  (npoi)     ! factor that respresents the annual average fuel
    ! dryness of a grid cell, and hence characterizes the readiness to burn
    REAL(KIND=r8)   , INTENT(IN   ) :: totlit   (npoi)     ! total carbon in all litter pools (kg_C m-2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: disturbf (npoi)     ! annual fire disturbance regime (m2/m2/yr)
    REAL(KIND=r8)   , INTENT(IN   ) :: vegtype0(npoi)  
    !
    ! local variables
    !
    INTEGER :: i
    !
    REAL(KIND=r8)    :: burn
    !
    ! begin global grid
    !
    DO i = 1, npoi
       IF(vegtype0(i) /= 15.0_r8)THEN
          !
          burn = firefac(i) * min (1.00_r8, totlit(i) / 0.2000_r8)
          !
          disturbf(i) = 1.00_r8 - exp(-0.50_r8 * burn)
          !
          disturbf(i) = max (0.00_r8, min (1.00_r8, disturbf(i)))
          !
       END IF
    END DO
    !
    RETURN
  END SUBROUTINE fire
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE vegmap(totlaiu , &! INTENT(IN   )
       plai    , &! INTENT(IN   )
       totlail , &! INTENT(IN   )
       vegtype0, &! INTENT(OUT  )
       gdd5, &! INTENT(IN   )
       gdd0    , &! INTENT(IN   )
       npoi    , &! INTENT(IN   )
       npft      )! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN   ) :: npoi                ! total number of land points
    INTEGER, INTENT(IN   ) :: npft                ! number of plant functional types
    REAL(KIND=r8)   , INTENT(IN   ) :: totlaiu (npoi)      ! total leaf area index for the upper canopy
    REAL(KIND=r8)   , INTENT(IN   ) :: plai    (npoi,npft) ! total leaf area index of each plant functional type
    REAL(KIND=r8)   , INTENT(IN   ) :: totlail (npoi)      ! total leaf area index for the lower canopy
    REAL(KIND=r8)   , INTENT(INOUT) :: vegtype0(npoi)      ! annual vegetation type - ibis classification
    REAL(KIND=r8)   , INTENT(IN   ) :: gdd5    (npoi)     ! growing degree days > 5C
    REAL(KIND=r8)   , INTENT(IN   ) :: gdd0    (npoi)     ! growing degree days > 0C 
    !
    ! local variables
    !
    INTEGER :: i 
    INTEGER :: j          ! loop indice
    INTEGER :: domtree    ! dominant tree
    !
    REAL(KIND=r8)    :: maxlai     ! maximum lai
    REAL(KIND=r8)    :: totlai    ! total ecosystem lai
    !      REAL(KIND=r8)    :: grassfrac  ! fraction of total lai in grasses
    !      REAL(KIND=r8)    :: treefrac   ! fraction of total lai in trees
    REAL(KIND=r8)    :: treelai    ! lai of trees
    REAL(KIND=r8)    :: shrublai   ! lai of shrubs
    REAL(KIND=r8)    :: grasslai   ! lai of grass
    REAL(KIND=r8)    :: ratio
    !
    ! classify vegetation cover into standard ibis vegetation classes 
    !
    ! ---------------------------------------------------
    !  1: tropical evergreen forest / woodland
    !  2: tropical deciduous forest / woodland
    !  3: temperate evergreen broadleaf forest / woodland
    !  4: temperate evergreen conifer forest / woodland
    !  5: temperate deciduous forest / woodland
    !  6: boREAL(KIND=r8) evergreen forest / woodland
    !  7: boREAL(KIND=r8) deciduous forest / woodland
    !  8: mixed forest / woodland
    !  9: savanna
    ! 10: grassland / steppe 
    ! 11: dense shrubland
    ! 12: open shrubland
    ! 13: tundra
    ! 14: desert 
    ! 15: polar desert / rock / ice
    ! ---------------------------------------------------
    !
    ! begin global grid
    !
    DO i = 1, npoi
       !
       ! determine total lai and tree, shrub, and grass fractions
       !
       treelai   = totlaiu(i) 
       shrublai  = plai(i,9)  + plai(i,10)
       grasslai  = plai(i,11) + plai(i,12)
       !
       totlai    = max (0.010_r8, totlail(i) + totlaiu(i))
       !
       ! determine dominant tree type by lai dominance
       !
       domtree = 0
       maxlai = 0.00_r8
       !
       DO j = 1, 8
          IF (plai(i,j).gt.maxlai) THEN
             domtree = j
             maxlai = plai(i,j)
          END IF
       END DO
       !
       ! assign initial vegetation type
       !
       !vegtype0(i) = -999.990_r8
       !
       ! dominant type:  tropical broadleaf evergreen tree
       !
       IF (domtree.eq.1) THEN
          IF (treelai.gt.2.50_r8)         vegtype0(i) =  1.00_r8  ! tropical evergreen forest / woodland
          IF (treelai.le.2.50_r8)         vegtype0(i) =  9.00_r8  ! savanna
          IF (treelai.le.0.50_r8) THEN
             IF (grasslai.ge.shrublai) vegtype0(i) = 10.00_r8  ! grassland
             IF (shrublai.ge.grasslai) vegtype0(i) = 11.00_r8  ! closed shrubland
          END IF
       END IF
       !
       ! dominant type:  tropical broadleaf drought-deciduous tree
       !
       IF (domtree.eq.2) THEN
          IF (treelai.gt.2.50_r8)         vegtype0(i) =  2.00_r8  ! tropical deciduous forest / woodland
          IF (treelai.le.2.50_r8)         vegtype0(i) =  9.00_r8  ! savanna
          IF (treelai.le.0.50_r8) THEN
             IF (grasslai.ge.shrublai) vegtype0(i) = 10.00_r8  ! grassland
             IF (shrublai.ge.grasslai) vegtype0(i) = 11.00_r8  ! closed shrubland
          END IF
       END IF
       !
       ! dominant type:  warm-temperate broadleaf evergreen tree
       !
       IF (domtree.eq.3) THEN
          IF (treelai.gt.2.50_r8)         vegtype0(i) =  3.00_r8  ! temperate evergreen broadleaf forest / woodland
          IF (treelai.le.2.50_r8)         vegtype0(i) =  9.00_r8  ! savanna
          IF (treelai.le.0.50_r8) THEN
             IF (grasslai.ge.shrublai) vegtype0(i) = 10.00_r8  ! grassland
             IF (shrublai.ge.grasslai) vegtype0(i) = 11.00_r8  ! closed shrubland
          END IF
       END IF
       !
       ! dominant type:  temperate conifer evergreen tree
       !
       IF (domtree.eq.4) THEN
          IF (treelai.gt.1.50_r8)         vegtype0(i) =  4.00_r8  ! temperate evergreen conifer forest / woodland
          IF (treelai.le.1.50_r8)         vegtype0(i) =  9.00_r8  ! savanna
          IF (treelai.le.0.50_r8) THEN
             IF (grasslai.ge.shrublai) vegtype0(i) = 10.00_r8  ! grassland
             IF (shrublai.ge.grasslai) vegtype0(i) = 11.00_r8  ! closed shrubland
          END IF
       END IF
       !
       ! dominant type:  temperate broadleaf deciduous tree
       !
       IF (domtree.eq.5) THEN
          IF (treelai.gt.1.50_r8)         vegtype0(i) =  5.00_r8  ! temperate deciduous forest / woodland
          IF (treelai.le.1.50_r8)         vegtype0(i) =  9.00_r8  ! savanna
          IF (treelai.le.0.50_r8) THEN
             IF (grasslai.ge.shrublai) vegtype0(i) = 10.00_r8  ! grassland
             IF (shrublai.ge.grasslai) vegtype0(i) = 11.00_r8  ! closed shrubland
          END IF
       END IF
       !
       ! dominant type:  boreal conifer evergreen tree
       !
       IF (domtree.eq.6)             vegtype0(i) =  6.00_r8  ! boreal evergreen forest / woodland
       !
       !       if (domtree.eq.6) then
       !         if (treelai.gt.1.0)         vegtype0(i) =  6.0  ! boreal evergreen forest / woodland
       !         if (treelai.le.1.0) then
       !           if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
       !           if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
       !         endif
       !       endif
       !
       ! dominant type:  boreal broadleaf cold-deciduous tree
       !
       IF (domtree.eq.7)             vegtype0(i) =  7.00_r8  ! boreal deciduous forest / woodland
       !
       !       if (domtree.eq.7) then
       !         if (treelai.gt.1.0)         vegtype0(i) =  7.0  ! boreal deciduous forest / woodland
       !         if (treelai.le.1.0) then
       !           if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
       !           if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
       !         endif
       !       endif
       !
       ! dominant type:  boreal conifer cold-deciduous tree
       !
       IF (domtree.eq.8)             vegtype0(i) =  7.00_r8  ! boreal deciduous forest / woodland
       !
       !       if (domtree.eq.8) then
       !         if (treelai.gt.1.0)         vegtype0(i) =  7.0  ! boreal deciduous forest / woodland
       !         if (treelai.le.1.0) then
       !           if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
       !           if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
       !         endif
       !       endif
       !
       ! temperate/boreal forest mixtures
       !
       IF ((domtree.ge.4).and.(domtree.le.8)) THEN
          ratio = (plai(i,5) + plai(i,7) + plai(i,8)) /  &
               (plai(i,4) + plai(i,5) + plai(i,6) +    &
               plai(i,7) + plai(i,8))
          IF (treelai.gt.1.00_r8) THEN
             IF ((ratio.gt.0.450_r8).and.(ratio.lt.0.550_r8)) vegtype0(i) = 8.0_r8
          END IF
          IF ((domtree.le.5).and.(treelai.le.1.00_r8)) THEN
             IF (grasslai.ge.shrublai) vegtype0(i) = 10.00_r8  ! grassland
             IF (shrublai.ge.grasslai) vegtype0(i) = 11.00_r8  ! closed shrubland
          END IF
       END IF
       !
       ! no tree is dominant
       !
       IF (domtree.eq.0) THEN
          IF (treelai.gt.1.00_r8)         vegtype0(i) =  9.00_r8  ! savanna
          IF (treelai.le.1.00_r8) THEN
             IF (grasslai.ge.shrublai) vegtype0(i) = 10.00_r8  ! grassland
             IF (shrublai.ge.grasslai) vegtype0(i) = 11.00_r8  ! closed shrubland
          END IF
       END IF
       !
       ! overriding vegtation classifications
       !
       IF (totlai.lt.1.00_r8)            vegtype0(i) = 12.00_r8  ! open shrubland
       IF (totlai.le.0.40_r8)            vegtype0(i) = 14.00_r8  ! desert
       !
       ! overriding climatic rules
       !
       IF (gdd5(i).lt.350.00_r8) THEN
          IF (totlai.ge.0.40_r8)          vegtype0(i) = 13.00_r8  ! tundra
          IF (totlai.lt.0.40_r8)          vegtype0(i) = 15.00_r8  ! polar desert
       END IF
       !
       IF (gdd0(i).lt.100.00_r8)         vegtype0(i) = 15.00_r8  ! polar desert
       !
    END DO! END DO i = 1, npoi
    !
    ! return to the main program
    !
    RETURN
  END SUBROUTINE vegmap

  !
  !  ####    #####    ##     #####   ####
  ! #          #     #  #      #    #
  !  ####      #    #    #     #     ####
  !      #     #    ######     #         #
  ! #    #     #    #    #     #    #    #
  !  ####      #    #    #     #     ####
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE sumnow(a10td, &! INTENT(INOUT) !global
       a10ancub, &! INTENT(INOUT) !global
       a10ancuc, &! INTENT(INOUT) !global
       a10ancls, &! INTENT(INOUT) !global
       a10ancl3, &! INTENT(INOUT) !global
       a10ancl4, &! INTENT(INOUT) !global
       nppdummy, &! INTENT(OUT  ) !local
       frac    , &! INTENT(IN   ) !global
       ancub, &! INTENT(IN   ) !global
       lai     , &! INTENT(IN   ) !global
       fu      , &! INTENT(IN   ) !global
       ancuc   , &! INTENT(IN   ) !global
       ancls   , &! INTENT(IN   ) !global
       fl      , &! INTENT(IN   ) !global
       ancl4   , &! INTENT(IN   ) !global
       ancl3   , &! INTENT(IN   ) !global
       tgpp    , &! INTENT(OUT  ) !local
       agcub   , &! INTENT(IN   ) !global
       agcuc   , &! INTENT(IN   ) !global
       agcls   , &! INTENT(IN   ) !global
       agcl4   , &! INTENT(IN   ) !global
       agcl3   , &! INTENT(IN   ) !global
       tgpptot , &! INTENT(OUT  ) !local
       ts      , &! INTENT(IN   ) !global
       froot   , &! INTENT(IN   ) !global
       tnpp    , &! INTENT(OUT  ) !local
       cbiow, &! INTENT(IN   ) !global
       sapfrac , &! INTENT(IN   ) !global
       cbior   , &! INTENT(IN   ) !global
       tnpptot , &! INTENT(OUT  ) !local
       tco2root, &! INTENT(OUT  ) !local
       tneetot , &! INTENT(OUT  ) !local
       tco2mic , &! INTENT(IN   ) !global
       tsoi    , &! INTENT(IN   ) !global
       fi      , &! INTENT(IN   ) !global
       td      , &! INTENT(IN   ) !global
       npoi    , &! INTENT(IN   ) !global
       nsoilay , &! INTENT(IN   ) !global
       npft, &! INTENT(IN   ) !global
       ndaypy  , &! INTENT(IN   ) !global
       dtime     )! INTENT(IN   ) !global
    ! ---------------------------------------------------------------------
    !
    ! common blocks
    !
    IMPLICIT NONE
    !
    INTEGER , INTENT(IN   ) :: npoi    ! total number of land points
    INTEGER , INTENT(IN   ) :: nsoilay ! number of soil layers
    INTEGER , INTENT(IN   ) :: npft ! number of plant functional types
    INTEGER , INTENT(IN   ) :: ndaypy ! number of days per year
    REAL(KIND=r8)    , INTENT(IN   ) :: dtime   ! model timestep (seconds)
    REAL(KIND=r8)    , INTENT(IN   ) :: td      (npoi)           ! daily average temperature (K)
    REAL(KIND=r8)    , INTENT(IN   ) :: fi      (npoi)         ! fractional snow cover
    REAL(KIND=r8)    , INTENT(IN   ) :: tsoi    (npoi,nsoilay)  ! soil temperature for each layer (K)
    REAL(KIND=r8)    , INTENT(OUT  ) :: nppdummy(npoi,npft) ! canopy NPP before accounting for stem and root respiration
    REAL(KIND=r8)    , INTENT(IN   ) :: frac    (npoi,npft) ! fraction of canopy occupied by each plant functional type
    REAL(KIND=r8)    , INTENT(IN   ) :: ancub   (npoi)      ! canopy average net photosynthesis rate - broadleaf    (mol_co2 m-2 s-1)
    REAL(KIND=r8)    , INTENT(IN   ) :: lai     (npoi,2)    ! canopy single-sided leaf area index (area leaf/area veg)
    REAL(KIND=r8)    , INTENT(IN   ) :: fu      (npoi)      ! fraction of overall area covered by upper canopy
    REAL(KIND=r8)    , INTENT(IN   ) :: ancuc   (npoi)      ! canopy average net photosynthesis rate - conifer      (mol_co2 m-2 s-1)
    REAL(KIND=r8)    , INTENT(IN   ) :: ancls   (npoi)      ! canopy average net photosynthesis rate - shrubs       (mol_co2 m-2 s-1)
    REAL(KIND=r8)    , INTENT(IN   ) :: fl      (npoi)      ! fraction of snow-free area covered by lower  canopy
    REAL(KIND=r8)    , INTENT(IN   ) :: ancl4   (npoi)      ! canopy average net photosynthesis rate - c4 grasses   (mol_co2 m-2 s-1)
    REAL(KIND=r8)    , INTENT(IN   ) :: ancl3   (npoi)      ! canopy average net photosynthesis rate - c3 grasses   (mol_co2 m-2 s-1)
    REAL(KIND=r8)    , INTENT(OUT  ) :: tgpp    (npoi,npft) ! instantaneous GPP for each pft (mol-CO2 / m-2 / second)
    REAL(KIND=r8)    , INTENT(IN   ) :: agcub   (npoi)      ! canopy average gross photosynthesis rate - broadleaf  (mol_co2 m-2 s-1)
    REAL(KIND=r8)    , INTENT(IN   ) :: agcuc   (npoi)      ! canopy average gross photosynthesis rate - conifer    (mol_co2 m-2 s-1)
    REAL(KIND=r8)    , INTENT(IN   ) :: agcls   (npoi)      ! canopy average gross photosynthesis rate - shrubs     (mol_co2 m-2 s-1)
    REAL(KIND=r8)    , INTENT(IN   ) :: agcl4   (npoi)      ! canopy average gross photosynthesis rate - c4 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8)    , INTENT(IN   ) :: agcl3   (npoi)      ! canopy average gross photosynthesis rate - c3 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8)    , INTENT(OUT  ) :: tgpptot (npoi)      ! instantaneous gpp (mol-CO2 / m-2 / second)
    REAL(KIND=r8)    , INTENT(IN   ) :: ts      (npoi)      ! temperature of upper canopy stems (K)
    REAL(KIND=r8)    , INTENT(IN   ) :: froot   (npoi,nsoilay,2) ! fraction of root in soil layer 
    REAL(KIND=r8)    , INTENT(OUT  ) :: tnpp    (npoi,npft) ! instantaneous NPP for each pft (mol-CO2 / m-2 / second)
    REAL(KIND=r8)    , INTENT(IN   ) :: cbiow   (npoi,npft) ! carbon in woody biomass pool (kg_C m-2)
    REAL(KIND=r8)    , INTENT(IN   ) :: sapfrac (npoi)      ! fraction of woody biomass that is in sapwood
    REAL(KIND=r8)    , INTENT(IN   ) :: cbior   (npoi,npft) ! carbon in fine root biomass pool (kg_C m-2)
    REAL(KIND=r8)    , INTENT(OUT  ) :: tnpptot (npoi)      ! instantaneous npp (mol-CO2 / m-2 / second)
    REAL(KIND=r8)    , INTENT(OUT  ) :: tco2root(npoi)      ! instantaneous fine co2 flux from soil (mol-CO2 / m-2 / second)
    REAL(KIND=r8)    , INTENT(OUT  ) :: tneetot (npoi)      ! instantaneous net ecosystem exchange of co2 per timestep (mol-CO2 / m-2 / second)
    REAL(KIND=r8)    , INTENT(IN   ) :: tco2mic (npoi)      ! instantaneous microbial co2 flux from soil (mol-CO2 / m-2 / second)
    REAL(KIND=r8)    , INTENT(INOUT) :: a10td    (npoi)     ! 10-day average daily air temperature (K)
    REAL(KIND=r8)    , INTENT(INOUT) :: a10ancub (npoi)     ! 10-day average canopy photosynthesis rate - broadleaf (mol_co2 m-2 s-1)
    REAL(KIND=r8)    , INTENT(INOUT) :: a10ancuc (npoi)     ! 10-day average canopy photosynthesis rate - conifer (mol_co2 m-2 s-1)
    REAL(KIND=r8)    , INTENT(INOUT) :: a10ancls (npoi)     ! 10-day average canopy photosynthesis rate - shrubs (mol_co2 m-2 s-1)
    REAL(KIND=r8)    , INTENT(INOUT) :: a10ancl3 (npoi)     ! 10-day average canopy photosynthesis rate - c3 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8)    , INTENT(INOUT) :: a10ancl4 (npoi)     ! 10-day average canopy photosynthesis rate - c4 grasses (mol_co2 m-2 s-1)
    !
    ! local variables
    !
    INTEGER :: i         ! loop indices
    INTEGER :: k         ! loop indices
    !
    REAL(KIND=r8)    :: tgpptot2 (npoi)      ! instantaneous gpp (mol-CO2 / m-2 / second)
    REAL(KIND=r8)    :: tnpp2    (npoi,npft) ! instantaneous NPP for each pft (mol-CO2 / m-2 / second)
    REAL(KIND=r8)    :: tnpptot2 (npoi)      ! instantaneous npp (mol-CO2 / m-2 / second)
    REAL(KIND=r8)    :: rwood    ! maintenance respiration coefficient for wood (/s)
    REAL(KIND=r8)    :: rroot    ! maintenance respiration coefficient for root (/s)
    REAL(KIND=r8)    :: rgrowth  ! growth respiration coefficient (fraction)
    REAL(KIND=r8)    :: stemtemp ! stem temperature
    REAL(KIND=r8)    :: roottemp ! average root temperature for all roots
    REAL(KIND=r8)    :: funca    ! temperature function for aboveground biomass (stems)
    REAL(KIND=r8)    :: funcb    ! temperature function for belowground biomass (roots)
    REAL(KIND=r8)    :: zweight  ! 10-day time averaging factor
    REAL(KIND=r8)    :: smask    ! 1 - fi
    !
    !
    ! ---------------------------------------------------------------------
    ! * * * define working variables * * *
    ! ---------------------------------------------------------------------
    !
    ! maintenance respiration coefficients (per second)
    !
    ! initially, we pick values for respiration coefficients that
    ! defined in units of  / year
    !
    !   rwood ~ 0.0125 
    !   rroot ~ 1.2500
    !
    ! however, we convert the unitsconvert to have resulting respiration
    ! fluxes in units of mol-C / m**2 / second
    !
    ! this requires we convert the time unit to seconds and add an additional
    ! factor to convert biomass units from kilograms to moles
    !
    rwood   = 0.0125_r8 / (ndaypy * 86400.0_r8) * (1000.0_r8 / 12.0_r8)
    rroot   = 1.2500_r8 / (ndaypy * 86400.0_r8) * (1000.0_r8 / 12.0_r8)
    !
    ! growth respiration coefficient (fraction)
    !
    rgrowth = 0.30_r8
    !
    ! 10-day time averaging factor
    !
    zweight = exp(-1.0_r8 / (10.0_r8 * 86400.0_r8 / dtime)) 
    !
    ! begin global grid
    !
    DO i = 1, npoi
       !
       ! calculate instantaneous carbon flux parameters, including
       ! npp (net primary production) and nee (net ecosystem exchange)
       !
       ! in this routine, all of the fluxes are calculated in the units
       ! of mol-C / m**2 / sec
       !
       ! ---------------------------------------------------------------------
       ! * * * calculate instantaneous GPP * * *
       ! ---------------------------------------------------------------------
       !
       ! snow masking for lower canopy vegetation
       !
       smask = 1.0_r8 - fi(i)
       !
       ! note that the following plants types follow different physiological paths
       !
       !   - broadleaf trees   :  types 1, 2, 3, 5, 7, 8 
       !   - conifer   trees   :  types 4, 6
       !   - shrubs            :  types 9, 10
       !   - c4 grasses        :  type 11
       !   - c3 grasses        :  type 12
       !
       ! note that plant type 8 is actually a deciduous conifer (e.g., Larix), but
       ! we are assuming that it's physiological behavior is like a broadleaf tree
       !
       ! nppdummy is canopy npp before accounting for stem & root respirtation
       ! Navin Sept 02
       !
       nppdummy(i,1)  = frac(i,1)  * ancub(i) * lai(i,2) * fu(i)
       nppdummy(i,2)  = frac(i,2)  * ancub(i) * lai(i,2) * fu(i)
       nppdummy(i,3)  = frac(i,3)  * ancub(i) * lai(i,2) * fu(i)
       nppdummy(i,4)  = frac(i,4)  * ancuc(i) * lai(i,2) * fu(i)
       nppdummy(i,5)  = frac(i,5)  * ancub(i) * lai(i,2) * fu(i)
       nppdummy(i,6)  = frac(i,6)  * ancuc(i) * lai(i,2) * fu(i)
       nppdummy(i,7)  = frac(i,7)  * ancub(i) * lai(i,2) * fu(i)
       nppdummy(i,8)  = frac(i,8)  * ancub(i) * lai(i,2) * fu(i)
       nppdummy(i,9)  = frac(i,9)  * ancls(i) * lai(i,1) * fl(i) * smask 
       nppdummy(i,10) = frac(i,10) * ancls(i) * lai(i,1) * fl(i) * smask
       nppdummy(i,11) = frac(i,11) * ancl4(i) * lai(i,1) * fl(i) * smask
       nppdummy(i,12) = frac(i,12) * ancl3(i) * lai(i,1) * fl(i) * smask
       !
       ! Navin's correction to compute npp using tgpp via agXXX
       ! agXXX should be used 
       !
       tgpp(i,1)  = frac(i,1)  * agcub(i) * lai(i,2) * fu(i)
       tgpp(i,2)  = frac(i,2)  * agcub(i) * lai(i,2) * fu(i)
       tgpp(i,3)  = frac(i,3)  * agcub(i) * lai(i,2) * fu(i)
       tgpp(i,4)  = frac(i,4)  * agcuc(i) * lai(i,2) * fu(i)
       tgpp(i,5)  = frac(i,5)  * agcub(i) * lai(i,2) * fu(i)
       tgpp(i,6)  = frac(i,6)  * agcuc(i) * lai(i,2) * fu(i)
       tgpp(i,7)  = frac(i,7)  * agcub(i) * lai(i,2) * fu(i)
       tgpp(i,8)  = frac(i,8)  * agcub(i) * lai(i,2) * fu(i)
       tgpp(i,9)  = frac(i,9)  * agcls(i) * lai(i,1) * fl(i) * smask 
       tgpp(i,10) = frac(i,10) * agcls(i) * lai(i,1) * fl(i) * smask
       tgpp(i,11) = frac(i,11) * agcl4(i) * lai(i,1) * fl(i) * smask
       tgpp(i,12) = frac(i,12) * agcl3(i) * lai(i,1) * fl(i) * smask
       !
       ! calculate total gridcell gpp
       !
       tgpptot(i) = 0.0_r8
       !
       DO k = 1, npft
          tgpptot(i) = tgpptot(i) + tgpp(i,k)
       END DO
       !
       ! calculate total gridcell gpp
       !
       tgpptot2(i) = 0.0_r8
       !
       DO k = 1, npft
          tgpptot2(i) = tgpptot2(i) + nppdummy(i,k)
       END DO

       !
       ! ---------------------------------------------------------------------
       ! * * * calculate temperature functions for respiration * * *
       ! ---------------------------------------------------------------------
       !
       ! calculate the stem temperature
       !
       stemtemp = MIN(MAX(ts(i),180.0_r8),350.0_r8)
       !
       ! calculate average root temperature (average of all roots)
       !
       roottemp = 0.0_r8
       !
       DO  k = 1, nsoilay
          roottemp = roottemp + tsoi(i,k) * 0.5_r8 *  &
               (froot(i,k,1) + froot(i,k,2))
       END DO
       roottemp = MIN(MAX(roottemp,180.0_r8),350.0_r8)
       !
       ! calculate respiration terms on a 15 degree base
       ! following respiration parameterization of Lloyd and Taylor
       !
       !        WRITE(*,*)ts(i)
       funca = exp(3500.0_r8 * (1.0_r8 / 288.16_r8 - 1.0_r8 / stemtemp))
       funcb = exp(3500.0_r8 * (1.0_r8 / 288.16_r8 - 1.0_r8 / roottemp))
       !
       ! ---------------------------------------------------------------------
       ! * * * calculate instantaneous NPP * * *
       ! ---------------------------------------------------------------------
       !
       ! the basic equation for npp is
       !
       !   npp = (1 - growth respiration term) * (gpp - maintenance respiration terms)
       !
       ! here the respiration terms are simulated as
       !
       !   growth respiration = rgrowth * (gpp - maintenance respiration terms)
       !
       ! where
       !
       !   rgrowth is the construction cost of new tissues
       !
       ! and
       !
       !   root respiration = rroot * cbior(i,k) * funcb
       !   wood respiration = rwood * cbiow(i,k) * funca * sapwood fraction
       !
       ! where
       ! 
       !   funca = temperature function for aboveground biomass (stems)
       !   funcb = temperature function for belowground biomass (roots)
       !
       ! note that we assume the sapwood fraction for shrubs is 1.0
       !
       ! also note that we apply growth respiration, (1 - rgrowth), 
       ! throughout the year; this may cause problems when comparing
       ! these npp values with flux tower measurements
       !
       ! also note that we need to convert the mass units of wood and
       ! root biomass from kilograms of carbon to moles of carbon
       ! to maintain consistent units (done in rwood, rroot)
       !
       ! finally, note that growth respiration is only applied to 
       ! positive carbon gains (i.e., when gpp-rmaint is positive)
       !
       ! Navin fix Sept 02 using nppdummy
       tnpp(i,1)  = tgpp(i,1)                           -   &
            rwood * cbiow(i,1) * sapfrac(i) * funca -   &
            rroot * cbior(i,1)              * funcb
       !
       tnpp(i,2)  = tgpp(i,2)                           -  &
            rwood * cbiow(i,2) * sapfrac(i) * funca -  &
            rroot * cbior(i,2)              * funcb
       !
       tnpp(i,3)  = tgpp(i,3)                           - &
            rwood * cbiow(i,3) * sapfrac(i) * funca - &
            rroot * cbior(i,3)              * funcb
       !
       tnpp(i,4)  = tgpp(i,4)                           - &
            rwood * cbiow(i,4) * sapfrac(i) * funca - &
            rroot * cbior(i,4)              * funcb
       !
       tnpp(i,5)  = tgpp(i,5)                           - &
            rwood * cbiow(i,5) * sapfrac(i) * funca - &
            rroot * cbior(i,5)              * funcb
       !
       tnpp(i,6)  = tgpp(i,6)                           -  &
            rwood * cbiow(i,6) * sapfrac(i) * funca -  &
            rroot * cbior(i,6)              * funcb
       !
       tnpp(i,7)  = tgpp(i,7)                           -   &
            rwood * cbiow(i,7) * sapfrac(i) * funca -   &
            rroot * cbior(i,7)              * funcb
       !
       tnpp(i,8)  = tgpp(i,8)                           -  &
            rwood * cbiow(i,8) * sapfrac(i) * funca -  &
            rroot * cbior(i,8)              * funcb
       !
       tnpp(i,9)  = tgpp(i,9)                           -  &
            rwood * cbiow(i,9)              * funca -  &
            rroot * cbior(i,9)              * funcb
       !
       tnpp(i,10) = tgpp(i,10) -                            &
            rwood * cbiow(i,10)             * funca -   & 
            rroot * cbior(i,10)             * funcb
       !
       tnpp(i,11) = tgpp(i,11) -   &
            rroot * cbior(i,11)            * funcb
       !
       tnpp(i,12) = tgpp(i,12) - &
            rroot * cbior(i,12)            * funcb
 
 !!!!!
        tnpp2(i,1)  = nppdummy(i,1)                           -   &
            rwood * cbiow(i,1) * sapfrac(i) * funca -   &
            rroot * cbior(i,1)              * funcb
       !
       tnpp2(i,2)  = nppdummy(i,2)                           -  &
            rwood * cbiow(i,2) * sapfrac(i) * funca -  &
            rroot * cbior(i,2)              * funcb
       !
       tnpp2(i,3)  = nppdummy(i,3)                           - &
            rwood * cbiow(i,3) * sapfrac(i) * funca - &
            rroot * cbior(i,3)              * funcb
       !
       tnpp2(i,4)  = nppdummy(i,4)                           - &
            rwood * cbiow(i,4) * sapfrac(i) * funca - &
            rroot * cbior(i,4)              * funcb
       !
       tnpp2(i,5)  = nppdummy(i,5)                           - &
            rwood * cbiow(i,5) * sapfrac(i) * funca - &
            rroot * cbior(i,5)              * funcb
       !
       tnpp2(i,6)  = nppdummy(i,6)                           -  &
            rwood * cbiow(i,6) * sapfrac(i) * funca -  &
            rroot * cbior(i,6)              * funcb
       !
       tnpp2(i,7)  = nppdummy(i,7)                           -   &
            rwood * cbiow(i,7) * sapfrac(i) * funca -   &
            rroot * cbior(i,7)              * funcb
       !
       tnpp2(i,8)  = nppdummy(i,8)                           -  &
            rwood * cbiow(i,8) * sapfrac(i) * funca -  &
            rroot * cbior(i,8)              * funcb
       !
       tnpp2(i,9)  = nppdummy(i,9)                           -  &
            rwood * cbiow(i,9)              * funca -  &
            rroot * cbior(i,9)              * funcb
       !
       tnpp2(i,10) = nppdummy(i,10) -                            &
            rwood * cbiow(i,10)             * funca -   & 
            rroot * cbior(i,10)             * funcb
       !
       tnpp2(i,11) = nppdummy(i,11) -   &
            rroot * cbior(i,11)            * funcb
       !
       tnpp2(i,12) = nppdummy(i,12) - &
            rroot * cbior(i,12)            * funcb
       !
       ! apply growth respiration and calculate total gridcell npp
       !
       tnpptot(i) = 0.0_r8
       !
       DO k = 1, npft
          IF (tnpp(i,k).gt.0.0_r8) THEN
             tnpp(i,k) = tnpp(i,k)  * (1.0_r8 - rgrowth)
          END IF
          tnpptot(i) = tnpptot(i) + tnpp(i,k)!! instantaneous npp (mol-CO2 / m-2 / second)
       END DO
!!
       !
       ! apply growth respiration and calculate total gridcell npp
       !
       tnpptot2(i) = 0.0_r8
       !
       DO k = 1, npft
          IF (tnpp2(i,k).gt.0.0_r8) THEN
             tnpp2(i,k) = tnpp2(i,k)  * (1.0_r8 - rgrowth)
          END IF
          tnpptot2(i) = tnpptot2(i) + tnpp2(i,k)!! instantaneous npp (mol-CO2 / m-2 / second)
       END DO

       !
       ! ---------------------------------------------------------------------
       ! * * * calculate total fine root respiration * * *
       ! ---------------------------------------------------------------------
       !
       tco2root(i) = 0.0_r8
       !
       DO k = 1, npft
                  !cbior -> carbon in fine root biomass pool (kg_C m-2)
                  !funcb -> temperature function for belowground biomass (roots)
                  !funcb -> exp(3500.0_r8 * (1.0_r8 / 288.16_r8 - 1.0_r8 / roottemp))
                  !roottemp -> average root temperature for all roots
                  
                  ! rroot   = 1.2500_r8 / (ndaypy * 86400.0_r8) * (1000.0_r8 / 12.0_r8) 
                  !         =  1/s                                  0.012kg  --> 

          tco2root(i) = tco2root(i) + rroot * cbior(i,k) * funcb !(mol-CO2 / m-2 / second)
       END DO
       !
       ! ---------------------------------------------------------------------
       ! * * * calculate instantaneous NEE * * *
       ! ---------------------------------------------------------------------
       !
       ! microbial respiration is calculated in biogeochem.f
       ! tco2mic(i)  -> instantaneous microbial co2 flux from soil (mol-CO2 / m-2 / second)
       !        WRITE(*,*)tnpptot(i) , tco2mic(i)
       tneetot(i) = tnpptot2(i) - tco2mic(i)!! instantaneous net ecosystem exchange of co2 per timestep  (mol-CO2 / m-2 / second)
       !
       ! ---------------------------------------------------------------------
       ! * * * update 10-day running-mean parameters * * *
       ! ---------------------------------------------------------------------
       !
       ! 10-day daily air temperature
       !
       a10td(i)    = zweight * a10td(i)    + (1.0_r8 - zweight) * td(i)
       !
       ! 10-day canopy photosynthesis rates
       !
       a10ancub(i) = zweight * a10ancub(i) + (1.0_r8 - zweight) * ancub(i)
       a10ancuc(i) = zweight * a10ancuc(i) + (1.0_r8 - zweight) * ancuc(i)
       a10ancls(i) = zweight * a10ancls(i) + (1.0_r8 - zweight) * ancls(i)
       a10ancl3(i) = zweight * a10ancl3(i) + (1.0_r8 - zweight) * ancl3(i)
       a10ancl4(i) = zweight * a10ancl4(i) + (1.0_r8 - zweight) * ancl4(i)
       !
    END DO  !DO  100 i = 1, npoi
    !
    ! return to main program
    !
    RETURN
  END SUBROUTINE sumnow
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE sumday (adnpp     , &! INTENT(INOUT) global
       tnpp      , &! INTENT(IN  ) !local
       raina     , &! INTENT(IN   )
       snowa     , &! INTENT(IN   )
       fvapa     , &! INTENT(IN   )
       grunof    , &! INTENT(IN   )
       gdrain    , &! INTENT(IN   )
       hsno      , &! INTENT(IN   )
       fi        , &! INTENT(IN   )
       hsoi      , &! INTENT(IN   )
       tsoi      , &! INTENT(IN   )
       wsoi      , &! INTENT(IN   )
       wisoi     , &! INTENT(IN   )
       ndtimes   , &! INTENT(INOUT) global
       adrain    , &! INTENT(INOUT) global
       adsnow    , &! INTENT(INOUT) global
       adaet     , &! INTENT(INOUT) global
       adtrunoff , &! INTENT(INOUT) global
       adsrunoff , &! INTENT(INOUT) global
       addrainage, &! INTENT(INOUT) global
       adrh      , &! INTENT(INOUT) global
       adsnod    , &! INTENT(INOUT) global
       adsnof    , &! INTENT(INOUT) global
       adwsoi    , &! INTENT(INOUT) global
       adtsoi    , &! INTENT(INOUT) global
       adwisoi   , &! INTENT(INOUT) global
       adtlaysoi , &! INTENT(INOUT) global
       adwlaysoi , &! INTENT(INOUT) global
       adwsoic   , &! INTENT(INOUT) global
       adtsoic   , &! INTENT(INOUT) global
       adco2mic  , &! INTENT(INOUT) global
       adco2root , &! INTENT(INOUT) global
       adco2soi  , &! INTENT(INOUT) global
       adco2ratio, &! INTENT(INOUT) global
       adnmintot , &! INTENT(INOUT) global
       froot     , &! INTENT(IN   )
       tco2mic   , &! INTENT(IN   )
       tco2root  , &! INTENT(IN   )
       decompl   , &! INTENT(INOUT) global
       decomps   , &! INTENT(INOUT) global
       tnmin     , &! INTENT(IN   )
       npoi      , &! INTENT(IN   )
       npft      , &! INTENT(IN   )
       nsoilay   , &! INTENT(IN   )
       nsnolay   , &! INTENT(IN   )
       dtime     , &! INTENT(IN   )
       td        , &! INTENT(INOUT)
       gdd0this  , &! INTENT(INOUT) global
       gdd5this  , &! INTENT(INOUT) global
       ts2       , &! INTENT(IN   )
       mcsec      ) ! INTENT(INOUT) global
    ! ---------------------------------------------------------------------
    !
    ! common blocks
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN   )  :: npoi ! total number of land points
    INTEGER, INTENT(IN   )  :: npft
    INTEGER, INTENT(IN   )  :: nsoilay ! number of soil layers
    INTEGER, INTENT(IN   )  :: nsnolay ! number of snow layers
    REAL(KIND=r8)   , INTENT(IN   )  :: dtime   ! model timestep (seconds)
    REAL(KIND=r8)   , INTENT(IN   ) :: tnpp  (npoi,npft)     ! instantaneous NPP for each pft (mol-CO2 / m-2 / second)
    REAL(KIND=r8)   , INTENT(IN   ) :: raina (npoi)       ! rainfall rate (mm/s or kg m-2 s-1)
    REAL(KIND=r8)   , INTENT(IN   ) :: snowa (npoi)    ! snowfall rate (mm/s or kg m-2 s-1 of water)
    REAL(KIND=r8)   , INTENT(IN   ) :: fvapa (npoi)       ! downward h2o vapor flux between za & z12 at za (kg m-2 s-1)
    REAL(KIND=r8)   , INTENT(IN   ) :: grunof(npoi)  ! surface runoff rate (kg_h2o m-2 s-1)
    REAL(KIND=r8)   , INTENT(IN   ) :: gdrain(npoi)  ! drainage rate out of bottom of lowest soil layer (kg_h2o m-2 s-1)
    REAL(KIND=r8)   , INTENT(IN   ) :: hsno  (npoi,nsnolay)   ! thickness of snow layers (m)
    REAL(KIND=r8)   , INTENT(IN   ) :: fi    (npoi)      ! fractional snow cover
    REAL(KIND=r8)   , INTENT(IN   ) :: hsoi  (npoi,nsoilay+1)  ! soil layer thickness (m)
    REAL(KIND=r8)   , INTENT(IN   ) :: tsoi  (npoi,nsoilay)  ! soil temperature for each layer (K)
    REAL(KIND=r8)   , INTENT(IN   ) :: wsoi  (npoi,nsoilay)  ! fraction of soil pore space containing liquid water
    REAL(KIND=r8)   , INTENT(IN   ) :: wisoi (npoi,nsoilay)  ! fraction of soil pore space containing ice
    INTEGER          , INTENT(INOUT) :: ndtimes    (npoi)          ! counter for daily average calculations

    REAL(KIND=r8)   , INTENT(INOUT) :: adrain    (npoi)    ! daily average rainfall rate (mm/day)
    REAL(KIND=r8)   , INTENT(INOUT) :: adsnow    (npoi)    ! daily average snowfall rate (mm/day)
    REAL(KIND=r8)   , INTENT(INOUT) :: adaet     (npoi)    ! daily average aet (mm/day)
    REAL(KIND=r8)   , INTENT(INOUT) :: adtrunoff (npoi)    ! daily average total runoff (mm/day)
    REAL(KIND=r8)   , INTENT(INOUT) :: adsrunoff (npoi)    ! daily average surface runoff (mm/day)
    REAL(KIND=r8)   , INTENT(INOUT) :: addrainage(npoi)    ! daily average drainage (mm/day)
    REAL(KIND=r8)   , INTENT(INOUT) :: adrh      (npoi)    ! daily average rh (percent)
    REAL(KIND=r8)   , INTENT(INOUT) :: adsnod    (npoi)    ! daily average snow depth (m)
    REAL(KIND=r8)   , INTENT(INOUT) :: adsnof    (npoi)    ! daily average snow fraction (fraction)
    REAL(KIND=r8)   , INTENT(INOUT) :: adwsoi    (npoi)    ! daily average soil moisture (fraction)
    REAL(KIND=r8)   , INTENT(INOUT) :: adtsoi    (npoi)    ! daily average soil temperature (c)
    REAL(KIND=r8)   , INTENT(INOUT) :: adwisoi   (npoi)    ! daily average soil ice (fraction)
    REAL(KIND=r8)   , INTENT(INOUT) :: adtlaysoi (npoi)    ! daily average soil temperature (c) of top layer
    REAL(KIND=r8)   , INTENT(INOUT) :: adwlaysoi (npoi)    ! daily average soil moisture of top layer(fraction)
    REAL(KIND=r8)   , INTENT(INOUT) :: adwsoic   (npoi)    ! daily average soil moisture using root profile weighting (fraction)
    REAL(KIND=r8)   , INTENT(INOUT) :: adtsoic   (npoi)    ! daily average soil temperature (c) using profile weighting
    REAL(KIND=r8)   , INTENT(INOUT) :: adco2mic  (npoi)    ! daily accumulated co2 respiration from microbes (kg_C m-2 /day)
    REAL(KIND=r8)   , INTENT(INOUT) :: adco2root (npoi)    ! daily accumulated co2 respiration from roots (kg_C m-2 /day)
    REAL(KIND=r8)   , INTENT(INOUT) :: adco2soi  (npoi)    ! daily accumulated co2 respiration from soil(total) (kg_C m-2 /day)
    REAL(KIND=r8)   , INTENT(INOUT) :: adco2ratio(npoi)    ! ratio of root to total co2 respiration
    REAL(KIND=r8)   , INTENT(INOUT) :: adnmintot (npoi)    ! daily accumulated net nitrogen mineralization (kg_N m-2 /day)
    REAL(KIND=r8)   , INTENT(INOUT) :: adnpp     (npoi,npft)! global! monthly total npp for each plant type (kg-C/m**2/day)

    REAL(KIND=r8)   , INTENT(IN   ) :: froot   (npoi,nsoilay,2) ! fraction of root in soil layer 
    REAL(KIND=r8)   , INTENT(IN   ) :: tco2mic (npoi)      ! instantaneous microbial co2 flux from soil (mol-CO2 / m-2 / second)
    REAL(KIND=r8)   , INTENT(IN   ) :: tco2root(npoi)      ! instantaneous fine co2 flux from soil (mol-CO2 / m-2 / second)
    REAL(KIND=r8)   , INTENT(INOUT) :: decompl (npoi)      ! litter decomposition factor    (dimensionless)
    REAL(KIND=r8)   , INTENT(INOUT) :: decomps (npoi)      ! soil organic matter decomposition factor     (dimensionless)
    REAL(KIND=r8)   , INTENT(IN   ) :: tnmin   (npoi)     ! instantaneous nitrogen mineralization (kg_N m-2/timestep)
    REAL(KIND=r8)   , INTENT(INOUT) :: td      (npoi)
    REAL(KIND=r8)   , INTENT(INOUT) :: gdd0this(npoi)       
    REAL(KIND=r8)   , INTENT(INOUT) :: gdd5this(npoi) 
    REAL(KIND=r8)   , INTENT(IN   ) :: ts2     (npoi) 
    REAL(KIND=r8)   , INTENT(IN   ) :: mcsec ! current seconds in day (0 - (86400 - dtime))     
    !
    ! Arguments
    !
    !      INTEGER, INTENT(IN   ) :: istep      ! daily timestep number (passed in)
    !
    ! local variables
    !
    INTEGER :: i          ! loop indices
    INTEGER :: k          ! loop indices
    !
    REAL(KIND=r8)    :: rwork      !working time variable
    REAL(KIND=r8)    :: rwork2     ! "
    REAL(KIND=r8)    :: rwork3     ! " 
    REAL(KIND=r8)    :: rwork4     ! "
    REAL(KIND=r8)    :: tconst     ! constant for Lloyd and Taylor (1994) function
    REAL(KIND=r8)    :: bconst     ! base temperature used for carbon decomposition
    REAL(KIND=r8)    :: btemp      ! maximum value of decomposition factor
    REAL(KIND=r8)    :: depth  (npoi)       ! total depth of the 4 1st soil layers
    REAL(KIND=r8)    :: depth2 (npoi)      ! total depth of the 2 1st soil layers
    REAL(KIND=r8)    :: zdepth       ! total depth of the 2 1st soil layers

    REAL(KIND=r8)    :: rdepth     ! total depth of the 4 1st soil layers
    REAL(KIND=r8)    :: rdepth2    ! total depth of the 2 1st soil layers
    REAL(KIND=r8)    :: snodpth    ! total snow depth
    REAL(KIND=r8)    :: soiltemp   ! average soil temp for 2 1st layers
    REAL(KIND=r8)    :: soilmois   ! average soil moisture (fraction of porosity) for 2 1st layers
    REAL(KIND=r8)    :: soilice    ! average soil ice for 2 1st layers
    REAL(KIND=r8)    :: soitempc   ! average soil temp over 6 layers
    REAL(KIND=r8)    :: soimoisc   ! average soil moisture over 6 layers
    REAL(KIND=r8)    :: factor     ! temperature decomposition factor for ltter/soil carbon
    REAL(KIND=r8)    :: wfps       ! water filled pore space
    REAL(KIND=r8)    :: moist      ! moisture effect on decomposition
    !      REAL(KIND=r8)    :: precipfac
    !
    !      INTEGER :: niter
    !
    ! ---------------------------------------------------------------------
    ! * * * update counters and working variables * * *
    ! ---------------------------------------------------------------------
    !
    ! reset sumday if the first timestep of the day 
    ! different from off-line IBIS where istep=1 :  1st timestep in the day
    !
    !      IF (istep .eq. 1) ndtimes = 0
    !
    rwork=dtime
    adrh=0.0_r8
    IF (mcsec .eq. 0.0_r8) THEN
       DO i = 1, npoi 
          ndtimes (i) = 0
          gdd0this(i) = gdd0this(i) + max(0.0_r8, (td(i) - 273.16_r8))
          gdd5this(i) = gdd5this(i) + max(0.0_r8, (td(i) - 278.16_r8))
       END DO
    END IF
    depth=0.0_r8
    depth2=0.0_r8
    DO  k = 1, nsoilay
        DO i = 1, npoi
           IF(depth(i)<=1.0_r8)THEN!m
             depth (i)=depth (i)+hsoi(i,k)
           END IF
           IF(depth2(i)<=0.30_r8)THEN!m
             depth2(i)=depth2(i)+hsoi(i,k)
           END IF
        END DO
    END DO
      
    DO i = 1, npoi
       !
       ! accumulate daily output (at this point for soil decomposition)
       !
       ndtimes(i) = ndtimes(i) + 1




       !
       ! working variables
       !
       rwork  = 1.0_r8 / real(ndtimes(i),kind=r8)
       rwork2 = 86400.0_r8
       rwork3 = 86400.0_r8 * 12.e-3_r8
       rwork4 = 86400.0_r8 * 14.e-3_r8
       !
       ! constants used in temperature function for c decomposition
       ! (arrhenius function constant) 
       !
       tconst  = 344.00_r8  ! constant for Lloyd and Taylor (1994) function
       btemp   = 288.16_r8  ! base temperature used for carbon decomposition
       !
       bconst  = 10.0_r8    ! maximum value of decomposition factor
       !
       ! soil weighting factors
       !
       rdepth  = 1.0_r8 / (depth(i))
       rdepth2 = 1.0_r8 / (depth2(i))
       !PK    rdepth  = 1.0_r8 / (hsoi(i,1) + hsoi(i,2) + hsoi(i,3) + hsoi(i,4))
       !PK    rdepth2 = 1.0_r8 / (hsoi(i,1) + hsoi(i,2))
       !
       ! begin global grid
       !
       !      DO i = 1, npoi
       !
       ! ---------------------------------------------------------------------
       ! * * * daily water budget terms * * *
       ! ---------------------------------------------------------------------
       !
       adrain(i)     = ((ndtimes(i)-1) * adrain(i) + raina(i) * 86400.0_r8) * rwork
       adsnow(i)     = ((ndtimes(i)-1) * adsnow(i) + snowa(i) * 86400.0_r8) * rwork
       adaet(i)      = ((ndtimes(i)-1) * adaet(i)  - fvapa(i) * 86400.0_r8) * rwork
       adtrunoff(i)  = ((ndtimes(i)-1) * adtrunoff(i)  +  (grunof(i) + gdrain(i)) * 86400.0_r8) * rwork
       adsrunoff(i)  = ((ndtimes(i)-1) * adsrunoff(i)  +   grunof(i)              * 86400.0_r8) * rwork
       addrainage(i) = ((ndtimes(i)-1) * addrainage(i) +   gdrain(i)              * 86400.0_r8) * rwork
       !
       ! ---------------------------------------------------------------------
       ! * * * daily atmospheric terms * * *
       ! ---------------------------------------------------------------------
       ! Different from off-line IBIS where td comes from climatology
       ! Daily mean temperature used for phenology based on 2-m screen
       ! temperature instead of 1st atmospheric level (~ 70 m)
       !
       !       td(i)      = ((ndtimes(i)-1) * td(i) + ta(i)) * rwork
       !
       td(i)   = ((ndtimes(i)-1) * td(i) + ts2(i)) * rwork
       !adrh(i) = ((ndtimes(i)-1) * adrh(i) + rh(i)) * rwork
       !
       ! ---------------------------------------------------------------------
       ! * * * daily snow parameters * * *
       ! ---------------------------------------------------------------------
       !
       snodpth = hsno(i,1) + hsno(i,2) + hsno(i,3)
       !
       adsnod(i) = ((ndtimes(i)-1) * adsnod(i) + snodpth) * rwork
       adsnof(i) = ((ndtimes(i)-1) * adsnof(i) + fi(i))   * rwork
       !
       ! ---------------------------------------------------------------------
       ! * * * soil parameters * * *
       ! ---------------------------------------------------------------------
       !
       ! initialize average soil parameters
       !
       soiltemp = 0.0_r8
       soilmois = 0.0_r8
       soilice  = 0.0_r8
       !
       soitempc = 0.0_r8
       soimoisc = 0.0_r8
       !
       ! averages for first 2 layers of soil
       !
        zdepth=0
        DO k = 1, nsoilay
           zdepth=zdepth+ hsoi(i,k)
            IF(zdepth<=0.30_r8)THEN!m
               soiltemp =  soiltemp + tsoi(i,k)  * hsoi(i,k)
               soilmois =  soilmois + wsoi(i,k)  * hsoi(i,k)
               soilice  =  soilice  + wisoi(i,k) * hsoi(i,k)
           END IF
        END DO
       !
       ! weighting on just thickness of each layer
       !
       soilmois = soilmois * rdepth2
       soilice  = soilice  * rdepth2
       soiltemp = soiltemp * rdepth2
       !
       ! calculate average root temperature, soil temperature and moisture and 
       ! ice content based on rooting profiles (weighted) from jackson et al
       ! 1996
       !
       ! these soil moisture and temperatures are used in biogeochem.f 
       ! we assume that the rooting profiles approximate
       ! where carbon resides in the soil
       !
       DO  k = 1, nsoilay

          soitempc = soitempc + tsoi(i,k)  * 0.5_r8 *  &
               (froot(i,k,1) + froot(i,k,2)) 
          soimoisc = soimoisc + wsoi(i,k)  * 0.5_r8 *  &
               (froot(i,k,1) + froot(i,k,2))

       END DO
       !
       ! calculate daily average soil moisture and soil ice
       ! using thickness of each layer as weighting function
       !
       adwsoi(i)  = ((ndtimes(i)-1) * adwsoi(i)  + soilmois) * rwork
       adtsoi(i)  = ((ndtimes(i)-1) * adtsoi(i)  + soiltemp) * rwork
       adwisoi(i) = ((ndtimes(i)-1) * adwisoi(i) + soilice)  * rwork
       !
       ! calculate daily average for soil temp/moisture of top layer
       !
       adtlaysoi(i) = ((ndtimes(i)-1) * adtlaysoi(i) + tsoi(i,1)) * rwork

       adwlaysoi(i) = ((ndtimes(i)-1) * adwlaysoi(i) + wsoi(i,1)) * rwork

       !
       ! calculate separate variables to keep track of weighting using 
       ! rooting profile information
       !
       ! note that these variables are only used for diagnostic purposes
       ! and that they are not needed in the biogeochemistry code
       !
       adwsoic(i)  = ((ndtimes(i)-1) * adwsoic(i) + soimoisc) * rwork
       adtsoic(i)  = ((ndtimes(i)-1) * adtsoic(i) + soitempc) * rwork
       !
       ! ---------------------------------------------------------------------
       ! * * * calculate daily soil co2 fluxes * * *
       ! ---------------------------------------------------------------------
       !
       ! increment daily total co2 respiration from microbes
       ! tco2mic is instantaneous value of co2 flux calculated in biogeochem.f
       !
       adco2mic(i) = ((ndtimes(i)-1) * adco2mic(i) + tco2mic(i) * rwork3) * rwork
       !
       ! increment daily total co2 respiration from fine roots
       ! tco2root is instantaneous value of co2 flux calculated in stats.f
       !
       adco2root(i) = ((ndtimes(i)-1) * adco2root(i) +    &
            tco2root(i) * rwork3) * rwork

       ! 
       ! calculate daily total co2 respiration from soil
       !
       adco2soi(i)  = adco2root(i) + adco2mic(i)
       !
       ! calculate daily ratio of total root to total co2 respiration
       !
       IF (adco2soi(i).gt.0.0_r8) THEN
          adco2ratio(i) = adco2root(i) / adco2soi(i)
       ELSE
          adco2ratio(i) = -999.99_r8
       END IF
       !
       ! ---------------------------------------------------------------------
       ! * * * calculate daily litter decomposition parameters * * *
       ! ---------------------------------------------------------------------
       !
       ! calculate litter carbon decomposition factors
       ! using soil temp, moisture and ice for top soil layer
       !
       ! calculation of soil biogeochemistry decomposition factors 
       ! based on moisture and temperature affects on microbial
       ! biomass dynamics
       !
       ! moisture function based on water-filled pore space (wfps)  
       ! williams et al., 1992 and friend et al., 1997 used in the
       ! hybrid 4.0 model; this is based on linn and doran, 1984
       !
       ! temperature functions are derived from arrhenius function
       ! found in lloyd and taylor, 1994 with a 15 c base 
       !
       ! calculate temperature decomposition factor
       ! CD impose lower limit to avoid division by zero at tsoi=227.13
       !
       IF (tsoi(i,1) .gt. 237.13_r8) THEN
          factor = min (exp(tconst * ((1.0_r8 / (btemp - 227.13_r8)) - (1.0_r8 /  &
               (tsoi(i,1)-227.13_r8)))), bconst)
       ELSE
          factor = exp(tconst * ((1.0_r8 / (btemp - 227.13_r8)) - (1.0_r8 /   &
               (237.13_r8-227.13_r8))))
       END IF
       !
       ! calculate water-filled pore space (in percent)
       !
       ! wsoi is relative to pore space not occupied by ice and water
       ! thus must include the ice fraction in the calculation
       !
       wfps = (1.0_r8 - wisoi(i,1)) * wsoi(i,1) * 100.0_r8
       !
       ! calculate moisture decomposition factor
       !
       IF (wfps .ge. 60.0_r8) THEN

          moist = 0.000371_r8 * (wfps**2) - (0.0748_r8 * wfps) + 4.13_r8

       ELSE

          moist = exp((wfps - 60.0_r8)**2 / (-800.0_r8))

       END IF
       !
       ! calculate combined temperature / moisture decomposition factor
       !
       factor = max (0.001_r8, min (bconst, factor * moist))
       !
       ! calculate daily average litter decomposition factor
       !
       decompl(i) = ((ndtimes(i)-1) * decompl(i) + factor) * rwork
       ! ---------------------------------------------------------------------
       ! * * * calculate daily soil carbon decomposition parameters * * *
       ! ---------------------------------------------------------------------
       !
       ! calculate soil carbon decomposition factors
       ! using soil temp, moisture and ice weighted by rooting profile scheme 
       !
       ! calculation of soil biogeochemistry decomposition factors 
       ! based on moisture and temperature affects on microbial
       ! biomass dynamics
       !
       ! moisture function based on water-filled pore space (wfps)  
       ! williams et al., 1992 and friend et al., 1997 used in the
       ! hybrid 4.0 model; this is based on linn and doran, 1984
       !
       ! temperature functions are derived from arrhenius function
       ! found in lloyd and taylor, 1994 with a 15 c base 
       !
       ! calculate temperature decomposition factor
       !
       IF (soiltemp .gt. 237.13_r8) THEN
          factor = min (exp(tconst * ((1.0_r8 / (btemp - 227.13_r8)) - (1.0_r8 /   &
               (soiltemp - 227.13_r8)))), bconst)
       ELSE
          factor = exp(tconst * ((1.0_r8 / (btemp - 227.13_r8)) - (1.0_r8 /   &
               (237.13_r8-227.13_r8))))
       END IF
       !
       ! calculate water-filled pore space (in percent)
       !
       ! wsoi is relative to pore space not occupied by ice and water
       ! thus must include the ice fraction in the calculation
       !
       wfps = (1.0_r8 - soilice) * soilmois * 100.0_r8
       !
       ! calculate moisture decomposition factor
       !
       IF (wfps .ge. 60.0_r8) THEN

          moist = 0.000371_r8 * (wfps**2) - (0.0748_r8 * wfps) + 4.13_r8

       ELSE

          moist = exp((wfps - 60.0_r8)**2 / (-800.0_r8))

       END IF
       !
       ! calculate combined temperature / moisture decomposition factor
       !
       factor = max (0.001_r8, min (bconst, factor * moist))
       !
       ! calculate daily average soil decomposition factor
       !
       decomps(i) = ((ndtimes(i)-1) * decomps(i) + factor) * rwork
       !
       ! ---------------------------------------------------------------------
       ! * * * calculate other daily biogeochemical parameters * * *
       ! ---------------------------------------------------------------------
       !
       ! increment daily total of net nitrogen mineralization
       ! value for tnmin is calculated in biogeochem.f
       !
       adnmintot(i) = ((ndtimes(i)-1) * adnmintot(i) +   &
            tnmin(i) * rwork4) * rwork



       ! ---------------------------------------------------------------------
       ! * * * determine daily npp * * *
       ! ---------------------------------------------------------------------
       !         adwsoi(i)  = ((ndtimes(i)-1) * adwsoi(i)  + soilmois) * rwork   ! average

       adnpp(i,1)  = ((ndtimes(i)-1) * adnpp(i,1)  + tnpp(i,1)  * rwork3) * rwork
       adnpp(i,2)  = ((ndtimes(i)-1) * adnpp(i,2)  + tnpp(i,2)  * rwork3) * rwork
       adnpp(i,3)  = ((ndtimes(i)-1) * adnpp(i,3)  + tnpp(i,3)  * rwork3) * rwork
       adnpp(i,4)  = ((ndtimes(i)-1) * adnpp(i,4)  + tnpp(i,4)  * rwork3) * rwork
       adnpp(i,5)  = ((ndtimes(i)-1) * adnpp(i,5)  + tnpp(i,5)  * rwork3) * rwork
       adnpp(i,6)  = ((ndtimes(i)-1) * adnpp(i,6)  + tnpp(i,6)  * rwork3) * rwork
       adnpp(i,7)  = ((ndtimes(i)-1) * adnpp(i,7)  + tnpp(i,7)  * rwork3) * rwork
       adnpp(i,8)  = ((ndtimes(i)-1) * adnpp(i,8)  + tnpp(i,8)  * rwork3) * rwork
       adnpp(i,9)  = ((ndtimes(i)-1) * adnpp(i,9)  + tnpp(i,9)  * rwork3) * rwork
       adnpp(i,10) = ((ndtimes(i)-1) * adnpp(i,10) + tnpp(i,10) * rwork3) * rwork
       adnpp(i,11) = ((ndtimes(i)-1) * adnpp(i,11) + tnpp(i,11) * rwork3) * rwork
       adnpp(i,12) = ((ndtimes(i)-1) * adnpp(i,12) + tnpp(i,12) * rwork3) * rwork

       !       adnpptot(i) = adnpp(i,1)  + adnpp(i,2)  + adnpp(i,3)  +  &
       !                     adnpp(i,4)  + adnpp(i,5)  + adnpp(i,6)  +  &
       !                     adnpp(i,7)  + adnpp(i,8)  + adnpp(i,9)  +  &
       !                     adnpp(i,10) + adnpp(i,11) + adnpp(i,12)
       

       !
    END DO !DO i = 1, npoi
    !
    ! return to main program
    !
    RETURN
  END SUBROUTINE sumday
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE summonth(&
       dtime     , &! INTENT(IN   )!global
       mcsec     , &! INTENT(IN   )!global
       iday      , &! INTENT(IN   )!global
       imonth    , &! INTENT(IN   )!global
       nmtimes   , &! INTENT(INOUT)!global
       amrain    , &! INTENT(INOUT)!global
       amsnow    , &! INTENT(INOUT)!global
       amaet     , &! INTENT(INOUT)!global
       amtrunoff , &! INTENT(INOUT)!global
       amsrunoff , &! INTENT(INOUT)!global
       amdrainage, &! INTENT(INOUT)!global
       amtemp    , &! INTENT(INOUT)!global
       amqa      , &! INTENT(INOUT)!global
       amsolar   , &! INTENT(INOUT)!global
       amirup    , &! INTENT(INOUT)!global
       amirdown  , &! INTENT(INOUT)!global
       amsens    , &! INTENT(INOUT)!global
       amlatent  , &! INTENT(INOUT)!global
       amlaiu    , &! INTENT(INOUT)!global
       amlail    , &! INTENT(INOUT)!global
       amtsoi    , &! INTENT(INOUT)!global
       amwsoi    , &! INTENT(INOUT)!global
       amwisoi   , &! INTENT(INOUT)!global
       amvwc    , &! INTENT(INOUT)!global
       amawc     , &! INTENT(INOUT)!global
       amsnod    , &! INTENT(INOUT)!global
       amsnof    , &! INTENT(INOUT)!global
       amnpp    , &! INTENT(INOUT)!global
       amnpptot  , &! INTENT(OUT  )!local
       amco2mic  , &! INTENT(INOUT)!global
       amco2root , &! INTENT(INOUT)!global
       amco2soi  , &! INTENT(OUT  )!local
       amco2ratio, &! INTENT(OUT  )!local
       amneetot  , &! INTENT(OUT  )!local
       amnmintot , &! INTENT(INOUT)!global
       amts2     , &! INTENT(INOUT)!global
       amtransu  , &! INTENT(INOUT)!global
       amtransl  , &! INTENT(INOUT)!global
       amsuvap   , &! INTENT(INOUT)!global
       aminvap   , &! INTENT(INOUT)!global
       amalbedo  , &! INTENT(INOUT)!global
       amtsoil   , &! INTENT(INOUT)!global
       amwsoil   , &! INTENT(INOUT)!global
       amwisoil  , &! INTENT(INOUT)!global
       ts2       , &! INTENT(INOUT)!global
       fu    , &! INTENT(IN   )!global
       lai       , &! INTENT(IN   )!global
       fl    , &! INTENT(IN   )!global
       tnpp      , &! INTENT(IN   )!global
       tco2mic   , &! INTENT(IN   )!global
       tco2root  , &! INTENT(IN   )!global
       tnmin     , &! INTENT(IN   )!global
       hsoi      , &! INTENT(IN   )!global
       tsoi    , &! INTENT(IN   )!global
       wsoi      , &! INTENT(IN   )!global
       wisoi     , &! INTENT(IN   )!global
       poros     , &! INTENT(IN   )!global
       swilt    , &! INTENT(IN   )!global
       hsno      , &! INTENT(IN   )!global
       fi    , &! INTENT(IN   )!global
       grunof    , &! INTENT(IN   )!global
       gdrain    , &! INTENT(IN   )!global
       gtransu   , &! INTENT(IN   )!global
       gtransl   , &! INTENT(IN   )!global
       gsuvap    , &! INTENT(IN   )!global
       ginvap    , &! INTENT(IN   )!global
       asurd     , &! INTENT(IN   )!global
       asuri     , &! INTENT(IN   )!global
       fvapa     , &! INTENT(IN   )!global
       firb      , &! INTENT(IN   )!global
       fsena     , &! INTENT(IN   )!global
       raina     , &! INTENT(IN   )!global
       snowa     , &! INTENT(IN   )!global
       ta        , &! INTENT(IN   )!global
       qa        , &! INTENT(IN   )!global
       solad     , &! INTENT(IN   )!global
       solai     , &! INTENT(IN   )!global
       fira      , &! INTENT(IN   )!global
       npoi      , &! INTENT(IN   )!global
       nband     , &! INTENT(IN   )!global
       nsoilay   , &! INTENT(IN   )!global
       nsnolay   , &! INTENT(IN   )!global
       npft      , &! INTENT(IN   )!global
       ndaypm    , &! INTENT(IN   )!global
       hvap        )! INTENT(IN   )!global
    ! ---------------------------------------------------------------------
    !
    ! first convert to units that make sense for output
    !
    !   - convert all temperatures to deg c
    !   - convert all liquid or vapor fluxes to mm/day
    !   - redefine upwd directed heat fluxes as positive
    !
    ! common blocks
    ! 
    IMPLICIT NONE
    !
    REAL(KIND=r8)   , INTENT(IN   ) :: dtime
    REAL(KIND=r8)   , INTENT(IN   ) :: mcsec     ! current seconds in day (0 - (86400 - dtime))
    INTEGER          , INTENT(IN   ) :: npoi      ! total number of land points
    INTEGER          , INTENT(IN   ) :: nband     ! number of solar radiation wavebands
    INTEGER          , INTENT(IN   ) :: nsoilay   ! number of soil layers
    INTEGER          , INTENT(IN   ) :: nsnolay   ! number of snow layers
    INTEGER          , INTENT(IN   ) :: npft   ! number of plant functional types
    INTEGER          , INTENT(IN   ) :: ndaypm(12)! number of days per month
    REAL(KIND=r8)   , INTENT(IN   ) :: hvap      ! latent heat of vaporization of water (J kg-1)
    REAL(KIND=r8)   , INTENT(IN   ) :: raina (npoi)       ! rainfall rate (mm/s or kg m-2 s-1)
    REAL(KIND=r8)   , INTENT(IN   ) :: snowa (npoi)    ! snowfall rate (mm/s or kg m-2 s-1 of water)
    REAL(KIND=r8)   , INTENT(IN   ) :: ta    (npoi)    ! air temperature (K)
    REAL(KIND=r8)   , INTENT(IN   ) :: qa    (npoi)    ! specific humidity (kg_h2o/kg_air)
    REAL(KIND=r8)   , INTENT(IN   ) :: solad (npoi,nband) ! direct downward solar flux (W m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: solai (npoi,nband) ! diffuse downward solar flux (W m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: fira  (npoi)    ! incoming ir flux (W m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: fvapa (npoi)      ! downward h2o vapor flux between za & z12 at za (kg m-2 s-1)
    REAL(KIND=r8)   , INTENT(IN   ) :: firb  (npoi)      ! net upward ir radiation at reference atmospheric level za (W m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: fsena (npoi)      ! downward sensible heat flux between za & z12 at za (W m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: grunof(npoi)  ! surface runoff rate (kg_h2o m-2 s-1)
    REAL(KIND=r8)   , INTENT(IN   ) :: gdrain(npoi)  ! drainage rate out of bottom of lowest soil layer (kg_h2o m-2 s-1)
    REAL(KIND=r8)   , INTENT(IN   ) :: gtransu (npoi)
    REAL(KIND=r8)   , INTENT(IN   ) :: gtransl (npoi)
    REAL(KIND=r8)   , INTENT(IN   ) :: gsuvap  (npoi)
    REAL(KIND=r8)   , INTENT(IN   ) :: ginvap  (npoi)
    REAL(KIND=r8)   , INTENT(IN   ) :: asurd   (npoi,nband) 
    REAL(KIND=r8)   , INTENT(IN   ) :: asuri   (npoi,nband) 
    REAL(KIND=r8)   , INTENT(IN   ) :: hsno(npoi,nsnolay)   ! thickness of snow layers (m)
    REAL(KIND=r8)   , INTENT(IN   ) :: fi  (npoi)      ! fractional snow cover
    REAL(KIND=r8)   , INTENT(IN   ) :: hsoi (npoi,nsoilay+1)        ! soil layer thickness (m)
    REAL(KIND=r8)   , INTENT(IN   ) :: tsoi (npoi,nsoilay)     ! soil temperature for each layer (K)
    REAL(KIND=r8)   , INTENT(IN   ) :: wsoi (npoi,nsoilay)     ! fraction of soil pore space containing liquid water
    REAL(KIND=r8)   , INTENT(IN   ) :: wisoi(npoi,nsoilay)     ! fraction of soil pore space containing ice
    REAL(KIND=r8)   , INTENT(IN   ) :: poros(npoi,nsoilay)     ! porosity (mass of h2o per unit vol at sat / rhow)
    REAL(KIND=r8)   , INTENT(IN   ) :: swilt(npoi,nsoilay)     ! wilting soil moisture value (fraction of pore space)
    REAL(KIND=r8)   , INTENT(IN   ) :: fu      (npoi)          ! fraction of overall area covered by upper canopy
    REAL(KIND=r8)   , INTENT(IN   ) :: lai     (npoi,2)        ! canopy single-sided leaf area index (area leaf/area veg)
    REAL(KIND=r8)   , INTENT(IN   ) :: fl      (npoi)          ! fraction of snow-free area covered by lower  canopy
    REAL(KIND=r8)   , INTENT(IN   ) :: tnpp    (npoi,npft)     ! instantaneous NPP for each pft (mol-CO2 / m-2 / second)
    REAL(KIND=r8)   , INTENT(IN   ) :: tco2mic (npoi)          ! instantaneous microbial co2 flux from soil (mol-CO2 / m-2 / second)
    REAL(KIND=r8)   , INTENT(IN   ) :: tco2root(npoi)          ! instantaneous fine co2 flux from soil (mol-CO2 / m-2 / second)
    REAL(KIND=r8)   , INTENT(IN   ) :: tnmin   (npoi)          ! instantaneous nitrogen mineralization (kg_N m-2/timestep)

    INTEGER          , INTENT(INOUT) :: nmtimes   (npoi)     ! counter for monthly average calculations
    REAL(KIND=r8)   , INTENT(INOUT) :: amrain    (npoi)     ! monthly average rainfall rate (mm/day)
    REAL(KIND=r8)   , INTENT(INOUT) :: amsnow    (npoi)     ! monthly average snowfall rate (mm/day)
    REAL(KIND=r8)   , INTENT(INOUT) :: amaet     (npoi)     ! monthly average aet (mm/day)
    REAL(KIND=r8)   , INTENT(INOUT) :: amtrunoff (npoi)     ! monthly average total runoff (mm/day)
    REAL(KIND=r8)   , INTENT(INOUT) :: amsrunoff (npoi)     ! monthly average surface runoff (mm/day)
    REAL(KIND=r8)   , INTENT(INOUT) :: amdrainage(npoi)     ! monthly average drainage (mm/day)
    REAL(KIND=r8)   , INTENT(INOUT) :: amtemp    (npoi)     ! monthly average air temperature (C)
    REAL(KIND=r8)   , INTENT(INOUT) :: amqa      (npoi)     ! monthly average specific humidity (kg-h2o/kg-air)
    REAL(KIND=r8)   , INTENT(INOUT) :: amsolar   (npoi)     ! monthly average incident solar radiation (W/m**2)
    REAL(KIND=r8)   , INTENT(INOUT) :: amirup    (npoi)     ! monthly average upward ir radiation (W/m**2)
    REAL(KIND=r8)   , INTENT(INOUT) :: amirdown  (npoi)     ! monthly average downward ir radiation (W/m**2)
    REAL(KIND=r8)   , INTENT(INOUT) :: amsens    (npoi)     ! monthly average sensible heat flux (W/m**2)
    REAL(KIND=r8)   , INTENT(INOUT) :: amlatent  (npoi)     ! monthly average latent heat flux (W/m**2)
    REAL(KIND=r8)   , INTENT(INOUT) :: amlaiu    (npoi)     ! monthly average lai for upper canopy (m**2/m**2)
    REAL(KIND=r8)   , INTENT(INOUT) :: amlail    (npoi)     ! monthly average lai for lower canopy (m**2/m**2)
    REAL(KIND=r8)   , INTENT(INOUT) :: amtsoi    (npoi)     ! monthly average 1m soil temperature (C)
    REAL(KIND=r8)   , INTENT(INOUT) :: amwsoi    (npoi)     ! monthly average 1m soil moisture (fraction)
    REAL(KIND=r8)   , INTENT(INOUT) :: amwisoi   (npoi)     ! monthly average 1m soil ice (fraction)
    REAL(KIND=r8)   , INTENT(INOUT) :: amvwc     (npoi)     ! monthly average 1m volumetric water content (fraction)
    REAL(KIND=r8)   , INTENT(INOUT) :: amawc     (npoi)     ! monthly average 1m plant-available water content (fraction)
    REAL(KIND=r8)   , INTENT(INOUT) :: amsnod    (npoi)     ! monthly average snow depth (m)
    REAL(KIND=r8)   , INTENT(INOUT) :: amsnof    (npoi)     ! monthly average snow fraction (fraction)
    REAL(KIND=r8)   , INTENT(INOUT) :: amnpp     (npoi,npft)! monthly total npp for each plant type (kg-C/m**2/month)
    REAL(KIND=r8)   , INTENT(OUT  ) :: amnpptot  (npoi)     ! monthly total npp for ecosystem (kg-C/m**2/month)
    REAL(KIND=r8)   , INTENT(INOUT) :: amco2mic  (npoi)     ! monthly total CO2 flux from microbial respiration (kg-C/m**2/month)
    REAL(KIND=r8)   , INTENT(INOUT) :: amco2root (npoi)     ! monthly total CO2 flux from soil due to root respiration (kg-C/m**2/month)
    REAL(KIND=r8)   , INTENT(OUT  ) :: amco2soi  (npoi)     ! monthly total soil CO2 flux from microbial
    ! and root respiration (kg-C/m**2/month)
    REAL(KIND=r8)   , INTENT(OUT  ) :: amco2ratio(npoi)     ! monthly ratio of root to total co2 flux
    REAL(KIND=r8)   , INTENT(OUT  ) :: amneetot  (npoi)     ! monthly total net ecosystem exchange of CO2 (kg-C/m**2/month)
    REAL(KIND=r8)   , INTENT(INOUT) :: amnmintot (npoi)     ! monthly total N mineralization from microbes (kg-N/m**2/month)
    REAL(KIND=r8)   , INTENT(INOUT) :: amts2     (npoi)     ! monthly average 2-m surface-air temperature 
    REAL(KIND=r8)   , INTENT(INOUT) :: amtransu  (npoi)     !
    REAL(KIND=r8)   , INTENT(INOUT) :: amtransl  (npoi)     !
    REAL(KIND=r8)   , INTENT(INOUT) :: amsuvap   (npoi)     !
    REAL(KIND=r8)   , INTENT(INOUT) :: aminvap   (npoi)     !
    REAL(KIND=r8)   , INTENT(INOUT) :: amalbedo  (npoi)     
    REAL(KIND=r8)   , INTENT(INOUT) :: amtsoil   (npoi, nsoilay) 
    REAL(KIND=r8)   , INTENT(INOUT) :: amwsoil   (npoi, nsoilay) 
    REAL(KIND=r8)   , INTENT(INOUT) :: amwisoil  (npoi, nsoilay)
    REAL(KIND=r8)   , INTENT(INOUT) :: ts2       (npoi)     ! monthly average 2-m surface-air temperature 


    !
    ! Arguments (input)
    !
    !      INTEGER, INTENT(IN   ) :: istep      ! daily timestep number (passed in)
    INTEGER, INTENT(IN   ) :: iday       ! day number  (passed in)
    INTEGER, INTENT(IN   ) :: imonth     ! month number (passed in)
    !
    ! local variables
    !
    INTEGER :: i 
    INTEGER :: k          ! loop indices
    !
    REAL(KIND=r8)    :: rwork     ! time work variable
    REAL(KIND=r8)    :: rwork2    !
    REAL(KIND=r8)    :: rwork3    !
    REAL(KIND=r8)    :: rwork4    !
    REAL(KIND=r8)    :: rdepth    ! 1/total soil depth over 4 1st layers
    REAL(KIND=r8)    :: solartot  ! total incoming radiation (direct + diffuse, visible + nearIR)
    REAL(KIND=r8)    :: soiltemp  ! average soil temp for 4 1st layers
    REAL(KIND=r8)    :: soilmois  ! average soil moisture for 4 1st layers 
    REAL(KIND=r8)    :: soilice   ! average soil ice for 4 1st layers 
    REAL(KIND=r8)    :: vwc       ! total liquid + ice content of 4 1st layers
    REAL(KIND=r8)    :: awc       ! total available water (+ ice) content of 4 1st layer
    REAL(KIND=r8)    :: snodpth   ! total snow depth
    REAL(KIND=r8)    :: depth (npoi)  
    REAL(KIND=r8)    :: depth2(npoi)  
    REAL(KIND=r8)    :: zdepth    ! total soil depth over 4 1st layers

    !
    REAL(KIND=r8)    :: albedotot
    !
    ! ---------------------------------------------------------------------
    ! * * * update counters and working variables * * *
    ! ---------------------------------------------------------------------
    ! 
    ! if the first timestep of the month then reset averages
    !
      depth =0.0_r8
      depth2=0.0_r8
      DO  k = 1, nsoilay
         DO i = 1, npoi
            IF(depth(i)<=1.0_r8)THEN!m
              depth (i)=depth (i)+hsoi(i,k)
            END IF
            IF(depth2(i)<=0.30_r8)THEN!m
              depth2(i)=depth2(i)+hsoi(i,k)
            END IF
         END DO
      END DO

    !IF ((istep.eq.1).and.(iday.eq.1)) nmtimes = 0
    albedotot=ta(1)*0.0_r8
    rwork=1.0_r8/dtime
    amtemp=amtemp
    DO i=1, npoi
       IF ((mcsec .eq. 0.0_r8) .and. (iday .eq. 1)) nmtimes(i) = 0
       !
       ! accumulate terms
       !

       !
       ! working variables
       !
       nmtimes(i) = nmtimes(i) + 1
       !
       ! rwork4 for conversion of nitrogen mineralization (moles)
       !
       rwork  = 1.0_r8 / real(nmtimes(i),kind=r8)
       rwork2 = real(ndaypm(imonth),kind=r8) * 86400.0_r8
       rwork3 = real(ndaypm(imonth),kind=r8) * 86400.0_r8 * 12.e-3_r8
       rwork4 = real(ndaypm(imonth),kind=r8) * 86400.0_r8 * 14.e-3_r8
       !
       !PK      rdepth = 1.0_r8 / (hsoi(i,1) + hsoi(i,2) + hsoi(i,3) + hsoi(i,4))
       rdepth = 1.0_r8 / (depth (i))
       !
       ! begin global grid
       !
       !do i = 1, npoi
       !
       ! monthly average temperature
       ! Different from offline IBIS where average T is from climatology
       !
       amts2(i) = ((nmtimes(i)-1) * amts2(i) + ts2(i)) * rwork
       !
       !      end do

       !     DO i = 1, npoi
       !
       ! ---------------------------------------------------------------------
       ! * * * monthly water budget terms * * *
       ! ---------------------------------------------------------------------
       ! 
       amrain(i)    = ((nmtimes(i)-1) * amrain(i) + raina(i) * 86400.0_r8) * rwork

       amsnow(i)    = ((nmtimes(i)-1) * amsnow(i) + snowa(i) * 86400.0_r8) * rwork

       amaet(i)     = ((nmtimes(i)-1) * amaet(i)  - fvapa(i) * 86400.0_r8) * rwork

       amtransu(i)     = ((nmtimes(i)-1) * amtransu(i) + gtransu(i) * 86400.0_r8) *rwork

       amtransl(i)     = ((nmtimes(i)-1) * amtransl(i) + gtransl(i) * 86400.0_r8) *rwork

       amsuvap(i)     = ((nmtimes(i)-1) * amsuvap(i) + gsuvap(i) * 86400.0_r8) *rwork

       aminvap(i)     = ((nmtimes(i)-1) * aminvap(i) + ginvap(i) * 86400.0_r8) *rwork

       amtrunoff(i)  = ((nmtimes(i)-1) * amtrunoff(i)  +        &
            (grunof(i) + gdrain(i)) * 86400.0_r8) * rwork

       amsrunoff(i)  = ((nmtimes(i)-1) * amsrunoff(i)  +        &
            grunof(i)              * 86400.0_r8) * rwork

       amdrainage(i) = ((nmtimes(i)-1) * amdrainage(i) +        &
            gdrain(i)              * 86400.0_r8) * rwork
       !
       ! ---------------------------------------------------------------------
       ! * * * monthly atmospheric terms * * *
       ! ---------------------------------------------------------------------
       !
       amqa(i)    = ((nmtimes(i)-1) * amqa(i)    + qa(i)) * rwork
       !
       ! ---------------------------------------------------------------------
       ! * * * energy budget terms * * *
       ! ---------------------------------------------------------------------
       !
       solartot = solad(i,1) + solad(i,2) + solai(i,1) + solai(i,2)
       !
       amsolar(i)  = ((nmtimes(i)-1) * amsolar(i)  +   &
            solartot)         * rwork
       amirup(i)   = ((nmtimes(i)-1) * amirup(i)   +   &
            firb(i))         * rwork
       amirdown(i) = ((nmtimes(i)-1) * amirdown(i) +   &
            fira(i))         * rwork
       amsens(i)   = ((nmtimes(i)-1) * amsens(i)   -   &
            fsena(i))        * rwork
       amlatent(i) = ((nmtimes(i)-1) * amlatent(i) -   &
            fvapa(i) * hvap) * rwork
       ! ---------------------------------------------------------------------
       ! ******* albedo calculations
       ! ---------------------------------------------------------------------
       albedotot = asurd(i,1) * solad(i,1) + &
            asurd(i,2) * solad(i,2) + &
            asuri(i,1) * solai(i,1) + &
            asuri(i,2) * solai(i,2)

       amalbedo(i) = ((nmtimes(i)-1) * amalbedo(i) + albedotot) * rwork
       !
       ! ---------------------------------------------------------------------
       ! * * * monthly vegetation parameters * * *
       ! ---------------------------------------------------------------------
       !
       amlaiu(i) = ((nmtimes(i)-1) * amlaiu(i) + fu(i) * lai(i,2)) * rwork
       amlail(i) = ((nmtimes(i)-1) * amlail(i) + fl(i) * lai(i,1)) * rwork
       !
       ! ---------------------------------------------------------------------
       ! * * * monthly soil parameters * * *
       ! ---------------------------------------------------------------------
       !
       soiltemp = 0.00_r8
       soilmois = 0.00_r8
       soilice  = 0.00_r8
       !
       vwc = 0.00_r8
       awc = 0.00_r8
       !
       ! averages for first 4 layers of soil (assumed to add to 1 meter depth)
       !
       zdepth=0.0_r8
       DO k = 1, nsoilay
          zdepth=zdepth+hsoi(i,k)
          IF(zdepth<=1.0_r8)THEN
             soiltemp =  soiltemp + tsoi(i,k)  * hsoi(i,k)
             soilmois =  soilmois + wsoi(i,k)  * hsoi(i,k)
             soilice  =  soilice  + wisoi(i,k) * hsoi(i,k)
       !
             vwc = vwc + (wisoi(i,k) + (1.0_r8 - wisoi(i,k)) * wsoi(i,k)) *  &
                         hsoi(i,k) * poros(i,k)
       !
             awc = awc + max (0.00_r8, (wisoi(i,k) +   &
                         (1.00_r8 - wisoi(i,k)) * wsoi(i,k)) - swilt(i,k)) *  &
                         hsoi(i,k) * poros(i,k) * 100.00_r8
       !
          END IF
       END DO
       !
       soiltemp = soiltemp * rdepth - 273.160_r8
       soilmois = soilmois * rdepth
       soilice  = soilice  * rdepth
       !
       vwc = vwc * rdepth
       awc = awc * rdepth
       !---------------------------------------------------------------------
       ! monthly average soil parameters:
       !---------------------------------------------------------------------
       amtsoi(i)  = ((nmtimes(i)-1) * amtsoi(i)  + soiltemp) * rwork
       amwsoi(i)  = ((nmtimes(i)-1) * amwsoi(i)  + soilmois) * rwork
       amwisoi(i) = ((nmtimes(i)-1) * amwisoi(i) + soilice)  * rwork
       amvwc(i)   = ((nmtimes(i)-1) * amvwc(i)   + vwc)      * rwork
       amawc(i)   = ((nmtimes(i)-1) * amawc(i)   + awc)      * rwork
       !
       ! Monthly averages per layer
       ! amalbedo amtsoil(i,k) amwsoil(i,k) amwisoil(i,k)
       do k = 1, nsoilay
          !
          amtsoil(i,k) = ((nmtimes(i)-1)*amtsoil(i,k)  &
               + tsoi(i,k)) *rwork
          amwsoil(i,k) = ((nmtimes(i)-1)*amwsoil(i,k)   &
               + wsoi(i,k)) * rwork
          amwisoil(i,k) = ((nmtimes(i)-1)*amwisoil(i,k)  &
               + wisoi(i,k)) * rwork
       end do
       !
       ! ---------------------------------------------------------------------
       ! * * * snow parameters * * *
       ! ---------------------------------------------------------------------
       !
       snodpth = hsno(i,1) + hsno(i,2) + hsno(i,3)
       !
       amsnod(i) = ((nmtimes(i)-1) * amsnod(i) + snodpth) * rwork
       amsnof(i) = ((nmtimes(i)-1) * amsnof(i) + fi(i))   * rwork
       !
       ! ---------------------------------------------------------------------
       ! * * * determine monthly npp * * *
       ! ---------------------------------------------------------------------
       !
       amnpp(i,1)  = ((nmtimes(i)-1) * amnpp(i,1)  + tnpp(i,1)  * rwork3) * rwork
       amnpp(i,2)  = ((nmtimes(i)-1) * amnpp(i,2)  + tnpp(i,2)  * rwork3) * rwork
       amnpp(i,3)  = ((nmtimes(i)-1) * amnpp(i,3)  + tnpp(i,3)  * rwork3) * rwork
       amnpp(i,4)  = ((nmtimes(i)-1) * amnpp(i,4)  +   &
            tnpp(i,4)  * rwork3) * rwork
       amnpp(i,5)  = ((nmtimes(i)-1) * amnpp(i,5)  +   &
            tnpp(i,5)  * rwork3) * rwork
       amnpp(i,6)  = ((nmtimes(i)-1) * amnpp(i,6)  +   &
            tnpp(i,6)  * rwork3) * rwork
       amnpp(i,7)  = ((nmtimes(i)-1) * amnpp(i,7)  +   &
            tnpp(i,7)  * rwork3) * rwork
       amnpp(i,8)  = ((nmtimes(i)-1) * amnpp(i,8)  +   &
            tnpp(i,8)  * rwork3) * rwork
       amnpp(i,9)  = ((nmtimes(i)-1) * amnpp(i,9)  +   &
            tnpp(i,9)  * rwork3) * rwork
       amnpp(i,10) = ((nmtimes(i)-1) * amnpp(i,10) +   &
            tnpp(i,10) * rwork3) * rwork
       amnpp(i,11) = ((nmtimes(i)-1) * amnpp(i,11) +   &
            tnpp(i,11) * rwork3) * rwork
       amnpp(i,12) = ((nmtimes(i)-1) * amnpp(i,12) +   &
            tnpp(i,12) * rwork3) * rwork
       !
       amnpptot(i) = amnpp(i,1)  + amnpp(i,2)  + amnpp(i,3)  +  &
            amnpp(i,4)  + amnpp(i,5)  + amnpp(i,6)  +  &
            amnpp(i,7)  + amnpp(i,8)  + amnpp(i,9)  +  &
            amnpp(i,10) + amnpp(i,11) + amnpp(i,12)
       !
       ! ---------------------------------------------------------------------
       ! * * * monthly biogeochemistry parameters * * *
       ! ---------------------------------------------------------------------
       !
       ! increment monthly total co2 respiration from microbes
       ! tco2mic is instantaneous value of co2 flux calculated in biogeochem.f
       !
       amco2mic(i) = ((nmtimes(i)-1) * amco2mic(i) +   &
            tco2mic(i) * rwork3) * rwork
       !
       ! increment monthly total co2 respiration from roots
       ! tco2root is instantaneous value of co2 flux calculated in stats.f
       !
       amco2root(i) = ((nmtimes(i)-1) * amco2root(i) +    &
            tco2root(i) * rwork3) * rwork
       !
       ! calculate average total co2 respiration from soil
       !
       amco2soi(i)  = amco2root(i) + amco2mic(i)
       !  
       !  calculate ratio of root to total co2 respiration
       !
       IF (amco2soi(i).gt.0.00_r8) THEN
          amco2ratio(i) = amco2root(i) / amco2soi(i)
       ELSE
          amco2ratio(i) = -999.990_r8
       END IF
       ! 
       !  monthly net ecosystem co2 flux -- npp total minus microbial respiration 
       !  the npp total includes losses from root respiration
       !
       amneetot(i)  = amnpptot(i) - amco2mic(i) 
       !
       ! increment monthly total of net nitrogen mineralization
       ! value for tnmin is calculated in biogeochem.f
       !
       amnmintot(i) = ((nmtimes(i)-1) * amnmintot(i) + tnmin(i) *  &
            rwork4) * rwork
       !
    END DO !DO 100 i = 1, npoi
    !
    ! return to main program
    !
    RETURN
  END SUBROUTINE summonth
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE sumyear(&
       dtime     , &! INTENT(IN   )
       mcsec     , &! INTENT(IN   )
       iday      , &! INTENT(IN   )
       imonth    , &! INTENT(IN   )
       wliqu     , &! INTENT(IN   )
       wsnou     , &! INTENT(IN   )
       fu   , &! INTENT(IN   )
       lai       , &! INTENT(IN   )
       wliqs     , &! INTENT(IN   )
       wsnos     , &! INTENT(IN   )
       sai       , &! INTENT(IN   )
       wliql     , &! INTENT(IN   )
       wsnol     , &! INTENT(IN   )
       fl        , &! INTENT(IN   )
       tgpp      , &! INTENT(IN   )
       tnpp      , &! INTENT(IN   )
       firefac   , &! INTENT(INOUT) global
       tco2mic   , &! INTENT(IN   )
       tco2root  , &! INTENT(IN   )
       cbior     , &! INTENT(IN   )
       tnmin     , &! INTENT(IN   )
       totalit   , &! INTENT(IN   )
       totrlit   , &! INTENT(IN   )
       totcsoi   , &! INTENT(IN   )
       totcmic   , &! INTENT(IN   )
       totanlit  , &! INTENT(IN   )
       totrnlit  , &! INTENT(IN   )
       totnsoi   , &! INTENT(IN   )
       nytimes   , &! INTENT(INOUT) global
       aysolar   , &! INTENT(INOUT) global
       ayirup    , &! INTENT(INOUT) global
       ayirdown  , &! INTENT(INOUT) global
       aysens    , &! INTENT(INOUT) global
       aylatent  , &! INTENT(INOUT) global
       ayprcp    , &! INTENT(INOUT) global
       ayaet     , &! INTENT(INOUT) global
       aytrans   , &! INTENT(INOUT) global
       aytrunoff , &! INTENT(INOUT) global
       aysrunoff , &! INTENT(INOUT) global
       aydrainage, &! INTENT(INOUT) global
       aydwtot   , &! INTENT(INOUT) global
       aywsoi    , &! INTENT(INOUT) global
       aywisoi   , &! INTENT(INOUT) global
       aytsoi    , &! INTENT(INOUT) global
       ayvwc     , &! INTENT(INOUT) global
       ayawc     , &! INTENT(INOUT) global
       aystresstu, &! INTENT(INOUT) global
       aystresstl, &! INTENT(INOUT) global
       aygpp     , &! INTENT(INOUT) global
       aygpptot  , &! INTENT(OUT  ) local
       aynpp     , &! INTENT(INOUT) global
       aynpptot  , &! INTENT(OUT  ) local
       ayco2mic  , &! INTENT(INOUT) global
       ayco2root , &! INTENT(INOUT) global
       ayco2soi  , &! INTENT(OUT  ) global
       ayneetot  , &! INTENT(OUT  ) global
       ayrootbio , &! INTENT(INOUT) global
       aynmintot , &! INTENT(INOUT) global
       ayalit    , &! INTENT(INOUT) global
       ayblit    , &! INTENT(INOUT) global
       aycsoi    , &! INTENT(INOUT) global
       aycmic    , &! INTENT(INOUT) global
       ayanlit   , &! INTENT(INOUT) global
       aybnlit   , &! INTENT(INOUT) global
       aynsoi    , &! INTENT(INOUT) global
       ayalbedo  , &! INTENT(INOUT) global
       hsoi   , &! INTENT(IN   ) global
       wpud      , &! INTENT(IN   ) global
       wipud     , &! INTENT(IN   ) global
       poros     , &! INTENT(IN   ) global
       wsoi   , &! INTENT(IN   ) global
       wisoi     , &! INTENT(IN   ) global
       tsoi      , &! INTENT(IN   ) global
       swilt     , &! INTENT(IN   ) global
       stresstu  , &! INTENT(IN   ) global
       stresstl  , &! INTENT(IN   ) global
       fi        , &! INTENT(IN   ) global
       rhos      , &! INTENT(IN   ) global
       hsno   , &! INTENT(IN   ) global
       gtrans    , &! INTENT(IN   ) global
       grunof    , &! INTENT(IN   ) global
       gdrain    , &! INTENT(IN   ) global
       wtot   , &! INTENT(INOUT) global
       firb      , &! INTENT(IN   ) global
       fsena     , &! INTENT(IN   ) global
       fvapa     , &! INTENT(IN   ) global
       solad     , &! INTENT(IN   ) global
       solai     , &! INTENT(IN   ) global
       fira      , &! INTENT(IN   ) global
       raina     , &! INTENT(IN   ) global
       snowa     , &! INTENT(IN   ) global
       asurd     , &! INTENT(IN   ) global
       asuri     , &! INTENT(IN   ) global
       npoi      , &! INTENT(IN   ) global
       nband     , &! INTENT(IN   ) global
       nsoilay   , &! INTENT(IN   ) global
       nsnolay   , &! INTENT(IN   ) global
       npft      , &! INTENT(IN   ) global
       ndaypy    , &! INTENT(IN   ) global
       hvap      , &! INTENT(IN   ) global
       rhow     )! INTENT(IN   ) global
    ! ---------------------------------------------------------------------
    !
    ! common blocks
    !
    IMPLICIT NONE
    !
    !      include 'compar.h'
    REAL(KIND=r8)   , INTENT(IN   ) :: dtime
    REAL(KIND=r8)   , INTENT(IN   ) :: mcsec     ! current seconds in day (0 - (86400 - dtime))
    INTEGER          , INTENT(IN   ) :: npoi ! total number of land points
    INTEGER          , INTENT(IN   ) :: nband ! number of solar radiation wavebands
    INTEGER          , INTENT(IN   ) :: nsoilay  ! number of soil layers
    INTEGER          , INTENT(IN   ) :: nsnolay  ! number of snow layers
    INTEGER          , INTENT(IN   ) :: npft ! number of plant functional types
    INTEGER          , INTENT(IN   ) :: ndaypy   ! number of days per year
    REAL(KIND=r8)   , INTENT(IN   ) :: hvap ! latent heat of vaporization of water (J kg-1)
    REAL(KIND=r8)   , INTENT(IN   ) :: rhow ! density of liquid water (all types) (kg m-3)
    !      include 'comatm.h'
    REAL(KIND=r8)   , INTENT(IN   ) :: solad(npoi,nband)  ! direct downward solar flux (W m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: solai(npoi,nband)  ! diffuse downward solar flux (W m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: fira (npoi)    ! incoming ir flux (W m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: raina(npoi)    ! rainfall rate (mm/s or kg m-2 s-1)
    REAL(KIND=r8)   , INTENT(IN   ) :: snowa(npoi)    ! snowfall rate (mm/s or kg m-2 s-1 of water)
    REAL(KIND=r8)   , INTENT(IN   ) :: asurd(npoi,nband)
    REAL(KIND=r8)   , INTENT(IN   ) :: asuri(npoi,nband)
    !      include 'com1d.h'
    REAL(KIND=r8)   , INTENT(IN   ) :: firb  (npoi)    ! net upward ir radiation at reference atmospheric level za (W m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: fsena (npoi)      ! downward sensible heat flux between za & z12 at za (W m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: fvapa (npoi)      ! downward h2o vapor flux between za & z12 at za (kg m-2 s-1)
    !      include 'comhyd.h'
    REAL(KIND=r8)   , INTENT(IN   ) :: gtrans (npoi)  ! total transpiration rate from all vegetation canopies (kg_h2o m-2 s-1)
    REAL(KIND=r8)   , INTENT(IN   ) :: grunof (npoi)  ! surface runoff rate (kg_h2o m-2 s-1)
    REAL(KIND=r8)   , INTENT(IN   ) :: gdrain (npoi)  ! drainage rate out of bottom of lowest soil layer (kg_h2o m-2 s-1)
    REAL(KIND=r8)   , INTENT(INOUT) :: wtot   (npoi)  ! total amount of water stored in snow, soil, puddels, and on vegetation (kg_h2o)
    !      include 'comsno.h'
    REAL(KIND=r8)   , INTENT(IN   ) :: fi    (npoi)        ! fractional snow cover
    REAL(KIND=r8)   , INTENT(IN   ) :: rhos                ! density of snow (kg m-3)
    REAL(KIND=r8)   , INTENT(IN   ) :: hsno  (npoi,nsnolay)! thickness of snow layers (m)
    !      include 'comsoi.h'
    REAL(KIND=r8)   , INTENT(IN   ) :: hsoi    (npoi,nsoilay+1)     ! soil layer thickness (m)
    REAL(KIND=r8)   , INTENT(IN   ) :: wpud    (npoi)          ! liquid content of puddles per soil area (kg m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: wipud   (npoi)          ! ice content of puddles per soil area (kg m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: poros   (npoi,nsoilay)  ! porosity (mass of h2o per unit vol at sat / rhow)
    REAL(KIND=r8)   , INTENT(IN   ) :: wsoi    (npoi,nsoilay)  ! fraction of soil pore space containing liquid water
    REAL(KIND=r8)   , INTENT(IN   ) :: wisoi   (npoi,nsoilay)  ! fraction of soil pore space containing ice
    REAL(KIND=r8)   , INTENT(IN   ) :: tsoi    (npoi,nsoilay)  ! soil temperature for each layer (K)
    REAL(KIND=r8)   , INTENT(IN   ) :: swilt   (npoi,nsoilay)  ! wilting soil moisture value (fraction of pore space)
    REAL(KIND=r8)   , INTENT(IN   ) :: stresstu(npoi)         ! sum of stressu over all 6 soil layers (dimensionless)
    REAL(KIND=r8)   , INTENT(IN   ) :: stresstl(npoi)         ! sum of stressl over all 6 soil layers (dimensionless)
    !      include 'comsum.h'
    INTEGER          , INTENT(INOUT) :: nytimes   (npoi)             ! counter for yearly average calculations
    REAL(KIND=r8)   , INTENT(INOUT) :: aysolar   (npoi)     ! annual average incident solar radiation (w/m**2)
    REAL(KIND=r8)   , INTENT(INOUT) :: ayirup    (npoi)     ! annual average upward ir radiation (w/m**2)
    REAL(KIND=r8)   , INTENT(INOUT) :: ayirdown  (npoi)     ! annual average downward ir radiation (w/m**2)
    REAL(KIND=r8)   , INTENT(INOUT) :: aysens    (npoi)     ! annual average sensible heat flux (w/m**2)
    REAL(KIND=r8)   , INTENT(INOUT) :: aylatent  (npoi)     ! annual average latent heat flux (w/m**2)
    REAL(KIND=r8)   , INTENT(INOUT) :: ayprcp    (npoi)     ! annual average precipitation (mm/yr)
    REAL(KIND=r8)   , INTENT(INOUT) :: ayaet     (npoi)     ! annual average aet (mm/yr)
    REAL(KIND=r8)   , INTENT(INOUT) :: aytrans   (npoi)     ! annual average transpiration (mm/yr)
    REAL(KIND=r8)   , INTENT(INOUT) :: aytrunoff (npoi)     ! annual average total runoff (mm/yr)
    REAL(KIND=r8)   , INTENT(INOUT) :: aysrunoff (npoi)     ! annual average surface runoff (mm/yr)
    REAL(KIND=r8)   , INTENT(INOUT) :: aydrainage(npoi)     ! annual average drainage (mm/yr)
    REAL(KIND=r8)   , INTENT(INOUT) :: aydwtot   (npoi)     ! annual average soil+vegetation+snow water recharge (mm/yr or kg_h2o/m**2/yr)
    REAL(KIND=r8)   , INTENT(INOUT) :: aywsoi    (npoi)     ! annual average 1m soil moisture (fraction)
    REAL(KIND=r8)   , INTENT(INOUT) :: aywisoi   (npoi)     ! annual average 1m soil ice (fraction)
    REAL(KIND=r8)   , INTENT(INOUT) :: aytsoi    (npoi)     ! annual average 1m soil temperature (C)
    REAL(KIND=r8)   , INTENT(INOUT) :: ayvwc     (npoi)     ! annual average 1m volumetric water content (fraction)
    REAL(KIND=r8)   , INTENT(INOUT) :: ayawc     (npoi)     ! annual average 1m plant-available water content (fraction)
    REAL(KIND=r8)   , INTENT(INOUT) :: aystresstu(npoi)     ! annual average soil moisture stress 
    ! parameter for upper canopy (dimensionless)
    REAL(KIND=r8)   , INTENT(INOUT) :: aystresstl(npoi)     ! annual average soil moisture stress 
    ! parameter for lower canopy (dimensionless)
    REAL(KIND=r8)   , INTENT(INOUT) :: aygpp     (npoi,npft)! annual gross npp for each plant type(kg-c/m**2/yr)
    REAL(KIND=r8)   , INTENT(OUT  ) :: aygpptot  (npoi)     ! annual total gpp for ecosystem (kg-c/m**2/yr)
    REAL(KIND=r8)   , INTENT(INOUT) :: aynpp     (npoi,npft)! annual total npp for each plant type(kg-c/m**2/yr)
    REAL(KIND=r8)   , INTENT(OUT  ) :: aynpptot  (npoi)     ! annual total npp for ecosystem (kg-c/m**2/yr)
    REAL(KIND=r8)   , INTENT(INOUT) :: ayco2mic  (npoi)     ! annual total CO2 flux from microbial respiration (kg-C/m**2/yr)
    REAL(KIND=r8)   , INTENT(INOUT) :: ayco2root (npoi)     ! annual total CO2 flux from soil due to root respiration (kg-C/m**2/yr)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ayco2soi  (npoi)     ! annual total soil CO2 flux from microbial and root respiration (kg-C/m**2/yr)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ayneetot  (npoi)     ! annual total NEE for ecosystem (kg-C/m**2/yr)
    REAL(KIND=r8)   , INTENT(INOUT) :: ayrootbio (npoi)     ! annual average live root biomass (kg-C / m**2)
    REAL(KIND=r8)   , INTENT(INOUT) :: aynmintot (npoi)     ! annual total nitrogen mineralization (kg-N/m**2/yr)
    REAL(KIND=r8)   , INTENT(INOUT) :: ayalit    (npoi)     ! aboveground litter (kg-c/m**2)
    REAL(KIND=r8)   , INTENT(INOUT) :: ayblit    (npoi)     ! belowground litter (kg-c/m**2)
    REAL(KIND=r8)   , INTENT(INOUT) :: aycsoi    (npoi)     ! total soil carbon (kg-c/m**2)
    REAL(KIND=r8)   , INTENT(INOUT) :: aycmic    (npoi)     ! total soil carbon in microbial biomass (kg-c/m**2)
    REAL(KIND=r8)   , INTENT(INOUT) :: ayanlit   (npoi)     ! aboveground litter nitrogen (kg-N/m**2)
    REAL(KIND=r8)   , INTENT(INOUT) :: aybnlit   (npoi)     ! belowground litter nitrogen (kg-N/m**2)
    REAL(KIND=r8)   , INTENT(INOUT) :: aynsoi    (npoi)     ! total soil nitrogen (kg-N/m**2)
    REAL(KIND=r8)   , INTENT(INOUT) :: ayalbedo  (npoi)  
    !      include 'comveg.h'
    REAL(KIND=r8)   , INTENT(IN   ) :: wliqu   (npoi)        ! intercepted liquid h2o on upper canopy leaf area (kg m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: wsnou   (npoi)        ! intercepted frozen h2o (snow) on upper canopy leaf area (kg m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: fu      (npoi)        ! fraction of overall area covered by upper canopy
    REAL(KIND=r8)   , INTENT(IN   ) :: lai     (npoi,2)      ! canopy single-sided leaf area index (area leaf/area veg)
    REAL(KIND=r8)   , INTENT(IN   ) :: wliqs   (npoi)        ! intercepted liquid h2o on upper canopy stem area (kg m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: wsnos   (npoi)        ! intercepted frozen h2o (snow) on upper canopy stem area (kg m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: sai     (npoi,2)      ! current single-sided stem area index
    REAL(KIND=r8)   , INTENT(IN   ) :: wliql   (npoi)        ! intercepted liquid h2o on lower canopy leaf and stem area (kg m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: wsnol   (npoi)        ! intercepted frozen h2o (snow) on lower canopy leaf & stem area (kg m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: fl      (npoi)        ! fraction of snow-free area covered by lower  canopy
    REAL(KIND=r8)   , INTENT(IN   ) :: tgpp    (npoi,npft)   ! instantaneous GPP for each pft (mol-CO2 / m-2 / second)
    REAL(KIND=r8)   , INTENT(IN   ) :: tnpp    (npoi,npft)   ! instantaneous NPP for each pft (mol-CO2 / m-2 / second)
    REAL(KIND=r8)   , INTENT(INOUT) :: firefac (npoi)        ! factor that respresents the annual average
    ! fuel dryness of a grid cell, and hence characterizes the readiness to burn
    REAL(KIND=r8)   , INTENT(IN   ) :: tco2mic (npoi)        ! instantaneous microbial co2 flux from soil (mol-CO2 / m-2 / second)
    REAL(KIND=r8)   , INTENT(IN   ) :: tco2root(npoi)        ! instantaneous fine co2 flux from soil (mol-CO2 / m-2 / second)
    REAL(KIND=r8)   , INTENT(IN   ) :: cbior   (npoi,npft)   ! carbon in fine root biomass pool (kg_C m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: tnmin   (npoi)        ! instantaneous nitrogen mineralization (kg_N m-2/timestep)
    REAL(KIND=r8)   , INTENT(IN   ) :: totalit (npoi)        ! total standing aboveground litter (kg_C m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: totrlit (npoi)        ! total root litter carbon belowground (kg_C m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: totcsoi (npoi)        ! total carbon in all soil pools (kg_C m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: totcmic (npoi)        ! total carbon residing in microbial pools (kg_C m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: totanlit(npoi)        ! total standing aboveground nitrogen in litter (kg_N m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: totrnlit(npoi)        ! total root litter nitrogen belowground (kg_N m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: totnsoi (npoi)        ! total nitrogen in soil (kg_N m-2)
    !
    ! Arguments (input)
    !
    !      INTEGER, INTENT(IN   ) :: istep     ! daily timestep number  (passed in)
    INTEGER, INTENT(IN   ) :: iday      ! day number  (passed in)
    INTEGER, INTENT(IN   ) :: imonth    ! month number (passed in)
    !
    ! local variables
    !
    INTEGER :: i        ! loop indices
    INTEGER :: k        ! loop indices
    !
    REAL(KIND=r8)    :: rwork     !
    REAL(KIND=r8)    :: rwork2    !
    REAL(KIND=r8)    :: rwork3    !
    REAL(KIND=r8)    :: rwork4    !
    REAL(KIND=r8)    :: rdepth    ! 1/total soil depth over 4 1st layers
    REAL(KIND=r8)    :: solartot  ! total incoming radiation (direct + diffuse, visible + nearIR)
    REAL(KIND=r8)    :: soiltemp  ! average soil temp for 4 1st layers
    REAL(KIND=r8)    :: soilmois  ! average soil moisture for 4 1st layers 
    REAL(KIND=r8)    :: soilice   ! average soil ice for 4 1st layers 
    REAL(KIND=r8)    :: vwc       ! total liquid + ice content of 4 1st layers
    REAL(KIND=r8)    :: awc       ! total available water (+ ice) content of 4 1st layer
    REAL(KIND=r8)    :: water     ! fire factor: total water content of 1st layer (liquid+ice)
    REAL(KIND=r8)    :: waterfrac ! fire factor: available water content of 1st layer
    REAL(KIND=r8)    :: fueldry   ! fire factor
    REAL(KIND=r8)    :: allroots  ! annual average root biomass
    REAL(KIND=r8)    :: wtotp     ! total water stored in soil+vegetation+snow at previous timestep
    REAL(KIND=r8)    :: albedotot
    REAL(KIND=r8)    :: depth (npoi) 
    REAL(KIND=r8)    :: depth2(npoi) 
    REAL(KIND=r8)    :: zdepth    ! total soil depth over 4 1st layers
    !
    ! ---------------------------------------------------------------------
    ! * * * update counters and working variables * * *
    ! ---------------------------------------------------------------------
    !
    ! reset sumyear if the first timestep of the year
    !
    !      IF ((istep.eq.1).and.(iday.eq.1).and.(imonth.eq.1)) nytimes = 0
    rwork=1.0_r8/dtime

    depth =0.0_r8
    depth2=0.0_r8
    DO  k = 1, nsoilay
       DO i = 1, npoi
          IF(depth(i)<=1.0_r8)THEN!m
            depth (i)=depth (i)+hsoi(i,k)
          END IF
          IF(depth2(i)<=0.30_r8)THEN!m
            depth2(i)=depth2(i)+hsoi(i,k)
          END IF
       END DO
    END DO
    
    DO i=1,npoi
       IF ((mcsec.eq.0.0_r8).and.(iday.eq.1).and.(imonth.eq.1)) nytimes(i) = 0

       !
       ! accumulate yearly output
       !
       nytimes(i) = nytimes(i) + 1
       !
       ! working variables
       !
       ! rwork4 is for nitrogen mineralization conversion
       !
       rwork  = 1.0_r8 / real(nytimes(i),kind=r8)
       rwork2 = real(ndaypy,kind=r8) * 86400.0_r8
       rwork3 = real(ndaypy,kind=r8) * 86400.0_r8 * 12.e-3_r8
       rwork4 = real(ndaypy,kind=r8) * 86400.0_r8 * 14.e-3_r8
       !
       !PK       rdepth = 1.0_r8 / (hsoi(i,1) + hsoi(i,2) + hsoi(i,3) + hsoi(i,4))
       rdepth = 1.0_r8 / (depth (i))
       !
       ! begin global grid
       !
       !DO i = 1, npoi
       !
       ! ---------------------------------------------------------------------
       ! * * * annual energy budget terms * * *
       ! ---------------------------------------------------------------------
       !
       solartot = solad(i,1) + solad(i,2) + solai(i,1) + solai(i,2)
       !
       albedotot = asurd(i,1) * solad(i,1) + &
            asurd(i,2) * solad(i,2) + &
            asuri(i,1) * solai(i,1) + &
            asuri(i,2) * solai(i,2)

       aysolar(i)  = ((nytimes(i)-1) * aysolar(i)  + solartot) * rwork  

       ayalbedo(i) = ((nytimes(i)-1) * ayalbedo(i) + albedotot) * rwork

       ayirup(i)   = ((nytimes(i)-1) * ayirup(i)   + firb(i))  * rwork

       ayirdown(i) = ((nytimes(i)-1) * ayirdown(i) + fira(i))  * rwork
       aysens(i)   = ((nytimes(i)-1) * aysens(i)   - fsena(i)) * rwork
       aylatent(i) = ((nytimes(i)-1) * aylatent(i) - fvapa(i)  * hvap)*  &
            rwork
       !
       ! ---------------------------------------------------------------------
       ! * * * annual water budget terms * * *
       ! ---------------------------------------------------------------------
       !
       ayprcp(i)     = ((nytimes(i)-1) * ayprcp(i)  +                  &
            (raina(i) + snowa(i)) * rwork2) * rwork

       ayaet(i)      = ((nytimes(i)-1) * ayaet(i)   -                  &
            fvapa(i)             * rwork2) * rwork

       aytrans(i)    = ((nytimes(i)-1) * aytrans(i) +                  &
            gtrans(i)            * rwork2) * rwork
       !
       aytrunoff(i)  = ((nytimes(i)-1) * aytrunoff(i)  +               &
            (grunof(i) + gdrain(i)) * rwork2) * rwork

       aysrunoff(i)  = ((nytimes(i)-1) * aysrunoff(i)  +               &   
            grunof(i)              * rwork2) * rwork

       aydrainage(i) = ((nytimes(i)-1) * aydrainage(i) +               &
            gdrain(i)  * rwork2)   * rwork
       !
       !---------------------------------------------------------------------
       ! CD
       ! estimate the change in soil-vegetation water content. Used to check 
       ! mass conservation
       !---------------------------------------------------------------------
       !
       wtotp = wtot(i)
       !
       wtot(i) = (wliqu(i)+wsnou(i)) * fu(i) * 2.00_r8 * lai(i,2) +   &
            (wliqs(i)+wsnos(i)) * fu(i) * 2.00_r8 * sai(i,2) +   &
            (wliql(i)+wsnol(i)) * fl(i) * 2.00_r8 *              &
            (lai(i,1) + sai(i,1)) * (1.0_r8 - fi(i))
       !
       wtot(i) = wtot(i) + wpud(i) + wipud(i)
       !
       DO  k = 1, nsoilay
          wtot(i) = wtot(i) +    &
               poros(i,k)*wsoi(i,k)*(1.0_r8-wisoi(i,k))*hsoi(i,k)*rhow+   &
               poros(i,k)*wisoi(i,k)*hsoi(i,k)*rhow
       END DO
       !
       DO k = 1, nsnolay
          wtot(i) = wtot(i) + fi(i)*rhos*hsno(i,k)
       END DO
       !
       aydwtot(i) = ((nytimes(i)-1) * aydwtot(i) +   &
            wtot(i) - wtotp) * rwork
       !
       ! ---------------------------------------------------------------------
       ! * * * annual soil parameters * * *
       ! ---------------------------------------------------------------------
       !
       soiltemp = 0.00_r8
       soilmois = 0.00_r8
       soilice  = 0.00_r8
       !
       vwc = 0.00_r8
       awc = 0.00_r8
       !
       ! averages for first 4 layers of soil
       !
        zdepth=0.0_r8
        DO k = 1, nsoilay
            zdepth=zdepth+hsoi(i,k)
            IF(zdepth<=1.0_r8)THEN!m
               soiltemp =  soiltemp + tsoi(i,k)  * hsoi(i,k)
               soilmois =  soilmois + wsoi(i,k)  * hsoi(i,k)
               soilice  =  soilice  + wisoi(i,k) * hsoi(i,k)
        !
               vwc = vwc + (wisoi(i,k) + (1.00_r8 - wisoi(i,k)) * wsoi(i,k)) *  &
                         hsoi(i,k) * poros(i,k)
        !
               awc = awc + max (0.00_r8, (wisoi(i,k) +                          &
                      (1.0_r8 - wisoi(i,k)) * wsoi(i,k)) - swilt(i,k)) *   &
                      hsoi(i,k) * poros(i,k) * 100.00_r8
        !
           END IF
        END DO
       !
       ! average soil and air temperatures
       !
       soiltemp = soiltemp * rdepth - 273.160_r8
       soilmois = soilmois * rdepth
       soilice  = soilice  * rdepth
       !
       vwc = vwc * rdepth
       awc = awc * rdepth
       !
       ! annual average soil moisture and soil ice
       !
       aywsoi(i)  = ((nytimes(i)-1) * aywsoi(i)  + soilmois) * rwork
       aywisoi(i) = ((nytimes(i)-1) * aywisoi(i) + soilice)  * rwork
       aytsoi(i)  = ((nytimes(i)-1) * aytsoi(i)  + soiltemp) * rwork
       ayvwc(i)   = ((nytimes(i)-1) * ayvwc(i)   + vwc)      * rwork
       ayawc(i)   = ((nytimes(i)-1) * ayawc(i)   + awc)      * rwork
       !
       ! soil moisture stress
       !
       aystresstu(i) = rwork * ((nytimes(i)-1) * aystresstu(i) + stresstu(i))
       !
       aystresstl(i) = rwork * ((nytimes(i)-1) * aystresstl(i) + stresstl(i))
       !
       ! ---------------------------------------------------------------------
       ! * * * determine annual gpp * * *
       ! ---------------------------------------------------------------------
       !
       ! gross primary production of each plant type
       !
       aygpp(i,1)  = ((nytimes(i)-1) * aygpp(i,1)  + tgpp(i,1)  * rwork3) * rwork
       aygpp(i,2)  = ((nytimes(i)-1) * aygpp(i,2)  + tgpp(i,2)  * rwork3) * rwork
       aygpp(i,3)  = ((nytimes(i)-1) * aygpp(i,3)  + tgpp(i,3)  * rwork3) * rwork
       aygpp(i,4)  = ((nytimes(i)-1) * aygpp(i,4)  + tgpp(i,4)  * rwork3) * rwork
       aygpp(i,5)  = ((nytimes(i)-1) * aygpp(i,5)  + tgpp(i,5)  * rwork3) * rwork
       aygpp(i,6)  = ((nytimes(i)-1) * aygpp(i,6)  + tgpp(i,6)  * rwork3) * rwork
       aygpp(i,7)  = ((nytimes(i)-1) * aygpp(i,7)  + tgpp(i,7)  * rwork3) * rwork
       aygpp(i,8)  = ((nytimes(i)-1) * aygpp(i,8)  + tgpp(i,8)  * rwork3) * rwork
       aygpp(i,9)  = ((nytimes(i)-1) * aygpp(i,9)  + tgpp(i,9)  * rwork3) * rwork
       aygpp(i,10) = ((nytimes(i)-1) * aygpp(i,10) + tgpp(i,10) * rwork3) * rwork
       aygpp(i,11) = ((nytimes(i)-1) * aygpp(i,11) + tgpp(i,11) * rwork3) * rwork
       aygpp(i,12) = ((nytimes(i)-1) * aygpp(i,12) + tgpp(i,12) * rwork3) * rwork
       !
       ! gross primary production of the entire gridcell
       !
       aygpptot(i) = aygpp(i,1)  + aygpp(i,2)  + aygpp(i,3)  +  &
            aygpp(i,4)  + aygpp(i,5)  + aygpp(i,6)  +  &
            aygpp(i,7)  + aygpp(i,8)  + aygpp(i,9)  +  &
            aygpp(i,10) + aygpp(i,11) + aygpp(i,12)
       !
       ! ---------------------------------------------------------------------
       ! * * * determine annual npp * * *
       ! ---------------------------------------------------------------------
       !
       ! net primary production of each plant type
       !
       aynpp(i,1)  = ((nytimes(i)-1) * aynpp(i,1) + tnpp(i,1)  * rwork3) * rwork
       aynpp(i,2)  = ((nytimes(i)-1) * aynpp(i,2) + tnpp(i,2)  * rwork3) * rwork
       aynpp(i,3)  = ((nytimes(i)-1) * aynpp(i,3) + tnpp(i,3)  * rwork3) * rwork
       aynpp(i,4)  = ((nytimes(i)-1) * aynpp(i,4) + tnpp(i,4)  * rwork3) * rwork
       aynpp(i,5)  = ((nytimes(i)-1) * aynpp(i,5) + tnpp(i,5)  * rwork3) * rwork
       aynpp(i,6)  = ((nytimes(i)-1) * aynpp(i,6) + tnpp(i,6)  * rwork3) * rwork
       aynpp(i,7)  = ((nytimes(i)-1) * aynpp(i,7) + tnpp(i,7)  * rwork3) * rwork
       aynpp(i,8)  = ((nytimes(i)-1) * aynpp(i,8) + tnpp(i,8)  * rwork3) * rwork
       aynpp(i,9)  = ((nytimes(i)-1) * aynpp(i,9) + tnpp(i,9)  * rwork3) * rwork
       aynpp(i,10) = ((nytimes(i)-1) * aynpp(i,10) + tnpp(i,10) * rwork3) * rwork
       aynpp(i,11) = ((nytimes(i)-1) * aynpp(i,11) + tnpp(i,11) * rwork3) * rwork
       aynpp(i,12) = ((nytimes(i)-1) * aynpp(i,12) + tnpp(i,12) * rwork3) * rwork
       !
       ! net primary production of the entire gridcell
       !
       aynpptot(i) = aynpp(i,1)  + aynpp(i,2)  + aynpp(i,3)  +  &
            aynpp(i,4)  + aynpp(i,5)  + aynpp(i,6)  +  &
            aynpp(i,7)  + aynpp(i,8)  + aynpp(i,9)  +  &
            aynpp(i,10) + aynpp(i,11) + aynpp(i,12)
       !
       ! ---------------------------------------------------------------------
       ! * * * annual carbon budget terms * * *
       ! ---------------------------------------------------------------------
       !
       ! fire factor used in vegetation dynamics calculations
       !
       water     = wisoi(i,1) + (1.0_r8 - wisoi(i,1)) * wsoi(i,1)
       waterfrac = (water - swilt(i,1)) / (1.0_r8 - swilt(i,1))
       !
       fueldry = max (0.00_r8, min (1.00_r8, -2.00_r8 * (waterfrac - 0.50_r8)))
       !
       firefac(i) = ((nytimes(i)-1) * firefac(i) + fueldry) * rwork
       !
       ! increment annual total co2 respiration from microbes
       ! tco2mic is instantaneous value of co2 flux calculated in biogeochem.f
       !
       ayco2mic(i) = ((nytimes(i)-1) * ayco2mic(i) +   &
            tco2mic(i) * rwork3) * rwork
       !
       ! increment annual total co2 respiration from roots
       !
       ayco2root(i) = ((nytimes(i)-1) * ayco2root(i) +  &
            tco2root(i) * rwork3) * rwork
       !
       ! calculate annual total co2 respiration from soil
       !
       ayco2soi(i)  = ayco2root(i) + ayco2mic(i)
       !  
       ! annual net ecosystem co2 flux -- npp total minus microbial respiration 
       ! the npp total includes losses from root respiration
       !
       ayneetot(i)  = aynpptot(i) - ayco2mic(i) 
       !
       ! annual average root biomass
       !
       allroots = cbior(i,1)  + cbior(i,2)  + cbior(i,3)  +  &
            cbior(i,4)  + cbior(i,5)  + cbior(i,6)  +  &
            cbior(i,7)  + cbior(i,8)  + cbior(i,9)  +  &
            cbior(i,10) + cbior(i,11) + cbior(i,12)
       !
       ayrootbio(i) =((nytimes(i)-1) * ayrootbio(i) + allroots) * rwork
       !
       ! ---------------------------------------------------------------------
       ! * * * annual biogeochemistry terms * * *
       ! ---------------------------------------------------------------------
       !
       ! increment annual total of net nitrogen mineralization
       ! value for tnmin is calculated in biogeochem.f
       !
       aynmintot(i) = ((nytimes(i)-1) * aynmintot(i) +  &
            tnmin(i) * rwork4) * rwork
       !
       ! other biogeochemistry variables
       !
       ayalit(i)  = ((nytimes(i)-1) * ayalit(i)  + totalit(i))  * rwork
       ayblit(i)  = ((nytimes(i)-1) * ayblit(i)  + totrlit(i))  * rwork
       aycsoi(i)  = ((nytimes(i)-1) * aycsoi(i)  + totcsoi(i))  * rwork
       aycmic(i)  = ((nytimes(i)-1) * aycmic(i)  + totcmic(i))  * rwork
       ayanlit(i) = ((nytimes(i)-1) * ayanlit(i) + totanlit(i)) * rwork
       aybnlit(i) = ((nytimes(i)-1) * aybnlit(i) + totrnlit(i)) * rwork
       aynsoi(i)  = ((nytimes(i)-1) * aynsoi(i)  + totnsoi(i))  * rwork
       !
    END DO! DO 100 i = 1, npoi
    !
    RETURN
  END SUBROUTINE sumyear
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE gdiag (iyear     , &! INTENT(IN   )
       iyear0    , &! INTENT(IN   )
       totbiou   , &! INTENT(IN   )
       totbiol   , &! INTENT(IN   )
       ayrratio  , &! INTENT(OUT  )
       aytrunoff , &! INTENT(IN   )
       ayprcp    , &! INTENT(IN   )
       aytratio  , &! INTENT(OUT  )
       aytrans   , &! INTENT(IN   )  
       ayaet  , &! INTENT(IN   )
       ayneetot  , &! INTENT(IN   )
       aynpptot  , &! INTENT(IN   )
       aygpptot  , &! INTENT(IN   )  
       ayalit  , &! INTENT(IN   )
       ayblit    , &! INTENT(IN   )
       aycsoi    , &! INTENT(IN   )
       ayco2soi  , &! INTENT(IN   ) 
       ayanlit   , &! INTENT(IN   )
       aybnlit   , &! INTENT(IN   )
       aynsoi    , &! INTENT(IN   )
       aysrunoff , &! INTENT(IN   )
       aydrainage, &! INTENT(IN   )
       aydwtot   , &! INTENT(IN   )
       nytimes   , &! INTENT(IN   )
       garea     , &! INTENT(IN   )
       npoi        )! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    ! common blocks
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN   ) :: npoi               ! total number of land points
    REAL(KIND=r8)   , INTENT(IN   ) :: garea     (npoi)   ! area of each gridcell (m**2)

    REAL(KIND=r8)   , INTENT(OUT  ) :: ayrratio  (npoi)   ! annual average runoff ratio (fraction)
    REAL(KIND=r8)   , INTENT(IN   ) :: aytrunoff (npoi)   ! annual average total runoff (mm/yr)
    REAL(KIND=r8)   , INTENT(IN   ) :: ayprcp    (npoi)   ! annual average precipitation (mm/yr)
    REAL(KIND=r8)   , INTENT(OUT  ) :: aytratio  (npoi)   ! annual average transpiration ratio (fraction)
    REAL(KIND=r8)   , INTENT(IN   ) :: aytrans   (npoi)   ! annual average transpiration (mm/yr)
    REAL(KIND=r8)   , INTENT(IN   ) :: ayaet     (npoi)   ! annual average aet (mm/yr)
    REAL(KIND=r8)   , INTENT(IN   ) :: ayneetot  (npoi)   ! annual total NEE for ecosystem (kg-C/m**2/yr)
    REAL(KIND=r8)   , INTENT(IN   ) :: aynpptot  (npoi)   ! annual total npp for ecosystem (kg-c/m**2/yr)
    REAL(KIND=r8)   , INTENT(IN   ) :: aygpptot  (npoi)   ! annual total gpp for ecosystem (kg-c/m**2/yr)
    REAL(KIND=r8)   , INTENT(IN   ) :: ayalit    (npoi)   ! aboveground litter (kg-c/m**2)
    REAL(KIND=r8)   , INTENT(IN   ) :: ayblit    (npoi)   ! belowground litter (kg-c/m**2)
    REAL(KIND=r8)   , INTENT(IN   ) :: aycsoi    (npoi)   ! total soil carbon (kg-c/m**2)
    REAL(KIND=r8)   , INTENT(IN   ) :: ayco2soi  (npoi)   ! annual total soil CO2 flux from microbial and root respiration (kg-C/m**2/yr)
    REAL(KIND=r8)   , INTENT(IN   ) :: ayanlit   (npoi)   ! aboveground litter nitrogen (kg-N/m**2)
    REAL(KIND=r8)   , INTENT(IN   ) :: aybnlit   (npoi)   ! belowground litter nitrogen (kg-N/m**2)
    REAL(KIND=r8)   , INTENT(IN   ) :: aynsoi    (npoi)   ! total soil nitrogen (kg-N/m**2)
    REAL(KIND=r8)   , INTENT(IN   ) :: aysrunoff (npoi)   ! annual average surface runoff (mm/yr)
    REAL(KIND=r8)   , INTENT(IN   ) :: aydrainage(npoi)   ! annual average drainage (mm/yr)
    REAL(KIND=r8)   , INTENT(IN   ) :: aydwtot   (npoi)   ! annual average soil+vegetation+snow water recharge (mm/yr or kg_h2o/m**2/yr)
    INTEGER, INTENT(IN   ) :: nytimes     (npoi)         ! counter for yearly average calculations

    REAL(KIND=r8)   , INTENT(IN   ) :: totbiou(npoi)      ! total biomass in the upper canopy (kg_C m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: totbiol(npoi)      ! total biomass in the lower canopy (kg_C m-2)

    !
    ! Arguments (input)
    !
    INTEGER, INTENT(IN   ) :: iyear    ! year counter
    INTEGER, INTENT(IN   ) :: iyear0   ! first year of simulation
    !
    ! local variables
    !
    INTEGER :: i        ! loop indice
    !
    REAL(KIND=r8)    :: gnee       ! domain total nee (gt-c/yr)
    REAL(KIND=r8)    :: gnpp       ! domain total npp (gt-c/yr)
    REAL(KIND=r8)    :: ggpp       ! domain total gpp (gt-c/yr)
    REAL(KIND=r8)    :: gbiomass   ! domain total biomass (gt-c)
    REAL(KIND=r8)    :: galitc     ! domain total aboveground litter carbon (gt-c)
    REAL(KIND=r8)    :: gblitc     ! domain total belowground litter carbon (gt-c)
    REAL(KIND=r8)    :: gsoic      ! domain total soil carbon (gt-c)
    REAL(KIND=r8)    :: gco2soi    ! domain total soil surface co2 flux (gt-c)
    REAL(KIND=r8)    :: galitn     ! domain total aboveground litter nitrogen (gt-c)
    REAL(KIND=r8)    :: gblitn     ! domain total belowground litter nitrogen (gt-c)
    REAL(KIND=r8)    :: gsoin      ! domain total soil nitrogen (gt-c)
    REAL(KIND=r8)    :: gprcp      ! domain average annual precipitation (mm/yr)
    REAL(KIND=r8)    :: gaet       ! domain average annual evapotranspiration (mm/yr)
    REAL(KIND=r8)    :: gt         ! domain average annual transpiration (mm/yr)
    REAL(KIND=r8)    :: gtrunoff   ! domain average total runoff (mm/yr)
    REAL(KIND=r8)    :: gsrunoff   ! domain average surface runoff (mm/yr)
    REAL(KIND=r8)    :: gdrainage  ! domain average drainage (mm/yr)
    REAL(KIND=r8)    :: gdwtot     !   "      "     water recharge (mm/yr)
    REAL(KIND=r8)    :: gtarea     ! total land area of the domain (m**2)
    REAL(KIND=r8)    :: aratio     ! aet / prcp ratio
    REAL(KIND=r8)    :: rratio     ! runoff / prcp ratio
    REAL(KIND=r8)    :: sratio     ! surface runoff / drainage ratio
    REAL(KIND=r8)    :: tratio     ! transpiration / aet ratio
    !
    ! initialize variables
    !
    gtarea    = 0.00_r8
    gnee      = 0.00_r8
    gnpp      = 0.00_r8
    ggpp      = 0.00_r8
    gbiomass  = 0.00_r8
    galitc    = 0.00_r8
    gblitc    = 0.00_r8
    gsoic     = 0.00_r8
    gco2soi   = 0.00_r8
    galitn    = 0.00_r8
    gblitn    = 0.00_r8
    gsoin     = 0.00_r8
    gprcp     = 0.00_r8
    gaet      = 0.00_r8
    gt        = 0.00_r8
    gtrunoff  = 0.00_r8
    gsrunoff  = 0.00_r8
    gdrainage = 0.00_r8
    gdwtot    = 0.00_r8
    !
    DO i = 1, npoi
       !
       ayrratio(i) = min (1.00_r8, max (0.00_r8, aytrunoff(i)) /  &
            max (0.10_r8, ayprcp(i)))
       !
       aytratio(i) = min (1.00_r8, max (0.00_r8, aytrans(i))   /  &
            max (0.10_r8, ayaet(i)))
       !
       gtarea    = gtarea    + garea(i)
       !
       gnee      = gnee      + garea(i) * ayneetot(i) * 1.e-12_r8
       gnpp      = gnpp      + garea(i) * aynpptot(i) * 1.e-12_r8
       ggpp      = ggpp      + garea(i) * aygpptot(i) * 1.e-12_r8
       gbiomass  = gbiomass  + garea(i) * totbiou(i)  * 1.e-12_r8 &
            + garea(i) * totbiol(i)  * 1.e-12_r8
       galitc    = galitc    + garea(i) * ayalit(i)   * 1.e-12_r8
       gblitc    = gblitc    + garea(i) * ayblit(i)   * 1.e-12_r8
       gsoic     = gsoic     + garea(i) * aycsoi(i)   * 1.e-12_r8
       gco2soi   = gco2soi   + garea(i) * ayco2soi(i) * 1.e-12_r8
       galitn    = galitn    + garea(i) * ayanlit(i)  * 1.e-12_r8
       gblitn    = gblitn    + garea(i) * aybnlit(i)  * 1.e-12_r8
       gsoin     = gsoin     + garea(i) * aynsoi(i)   * 1.e-12_r8
       !
       gprcp     = gprcp     + garea(i) * ayprcp(i)
       gaet      = gaet      + garea(i) * ayaet(i)
       gt        = gt        + garea(i) * aytrans(i)
       gtrunoff  = gtrunoff  + garea(i) * aytrunoff(i)
       gsrunoff  = gsrunoff  + garea(i) * aysrunoff(i)
       gdrainage = gdrainage + garea(i) * aydrainage(i)
       gdwtot    = gdwtot    + garea(i) * aydwtot(i)*nytimes(i)
       !
    END DO !DO i = 1, npoi
    !
    gprcp     = gprcp     / gtarea
    gaet      = gaet      / gtarea
    gt        = gt        / gtarea
    gtrunoff  = gtrunoff  / gtarea
    gsrunoff  = gsrunoff  / gtarea
    gdrainage = gdrainage / gtarea
    gdwtot    = gdwtot    / gtarea
    !
    aratio   = gaet     / gprcp
    rratio   = gtrunoff / gprcp
    sratio   = gsrunoff / gtrunoff
    tratio   = gt       / gaet
    !
    WRITE (*,*) ' '
    WRITE (*,*) '* * * annual diagnostic fields * * *'
    WRITE (*,*) ' '
    WRITE (*,9001) gnee
    WRITE (*,9000) gnpp
    WRITE (*,9002) ggpp
    WRITE (*,9010) gbiomass
    WRITE (*,9020) galitc
    WRITE (*,9021) gblitc
    WRITE (*,9030) gsoic
    WRITE (*,9032) gco2soi
    WRITE (*,9034) galitn
    WRITE (*,9036) gblitn
    WRITE (*,9038) gsoin
    WRITE (*,*) ' '
    WRITE (*,9040) gprcp
    WRITE (*,9050) gaet
    WRITE (*,9060) gt
    WRITE (*,9070) gtrunoff
    WRITE (*,9080) gsrunoff
    WRITE (*,9090) gdrainage
    WRITE (*,9095) gdwtot
    WRITE (*,*) ' '
    WRITE (*,9100) aratio
    WRITE (*,9110) rratio
    WRITE (*,*) ' '
    WRITE (*,9120) tratio
    WRITE (*,9130) sratio
    WRITE (*,*) ' '
    !
    ! WRITE some diagnostic output to history file
    !
    IF (iyear.eq.iyear0) THEN
       !
       OPEN (20,file='ibis.out.global',status='unknown')
       !
       WRITE (20,*) ' '
       WRITE (20,*) '* * * annual diagnostic fields * * *'
       WRITE (20,*) ' '
       WRITE (20,*) &
            'year       nee       npp       gpp   biomass   scarbon '// &
            'snitrogen   alitter   blitter    co2soi    '  //&
            'aratio    rratio    tratio'
       !
    END IF
    !
    WRITE (20,9500) iyear, gnee, gnpp, ggpp, gbiomass, gsoic, &
         gsoin, galitc,&
         gblitc, gco2soi, 100.00_r8 * aratio, 100.00_r8 * rratio,&
         100.00_r8 * tratio
    !
    CALL flush (20)
    !
    !     close (20)
    !
9000 FORMAT (1x,'total npp             of the domain (gt-c/yr) : ', &
         f12.3)
9001 FORMAT (1x,'total nee             of the domain (gt-c/yr) : ', &
         f12.5)
9002 FORMAT (1x,'total gpp             of the domain (gt-c/yr) : ', &
         f12.3)
9010 FORMAT (1x,'total biomass         of the domain (gt-c)    : ', &
         f12.3)
9020 FORMAT (1x,'aboveground litter    of the domain (gt-c)    : ', &
         f12.3)
9021 FORMAT (1x,'belowground litter    of the domain (gt-c)    : ', &
         f12.3)
9030 FORMAT (1x,'total soil carbon     of the domain (gt-c)    : ', &
         f12.3)
9032 FORMAT (1x,'total soil co2 flux   of the domain (gt-c)    : ', &
         f12.3)
9034 FORMAT (1x,'aboveground litter n  of the domain (gt-c)    : ', &
         f12.3)
9036 FORMAT (1x,'belowground litter n  of the domain (gt-c)    : ', &
         f12.3)
9038 FORMAT (1x,'total soil nitrogen   of the domain (gt-c)    : ', &
         f12.3)
9040 FORMAT (1x,'average precipitation of the domain (mm/yr)   : ', &
         f12.3)
9050 FORMAT (1x,'average aet           of the domain (mm/yr)   : ', &
         f12.3)
9060 FORMAT (1x,'average transpiration of the domain (mm/yr)   : ', &
         f12.3)
9070 FORMAT (1x,'average runoff        of the domain (mm/yr)   : ', &
         f12.3)
9080 FORMAT (1x,'average surf runoff   of the domain (mm/yr)   : ', &
         f12.3)
9090 FORMAT (1x,'average drainage      of the domain (mm/yr)   : ', &
         f12.3)
9095 FORMAT (1x,'average moisture recharge of the domain (mm/yr) : ', &
         f12.3)
9100 FORMAT (1x,'total aet      / precipitation                : ', &
         f12.3)
9110 FORMAT (1x,'total runoff   / precipitation                : ', &
         f12.3)
9120 FORMAT (1x,'transpiration  / total aet                    : ', &
         f12.3)
9130 FORMAT (1x,'surface runoff / total runoff                 : ', &
         f12.3)
9500 FORMAT (1x,i4,12f10.2)
    !
    ! return to main program
    !
    RETURN
  END SUBROUTINE gdiag
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE vdiag (iyear    , &! INTENT(IN   )
       iyear0   , &! INTENT(IN   )
       vegtype0 , &! INTENT(IN   )
       totbiou  , &! INTENT(IN   )
       totbiol  , &! INTENT(IN   )
       totlaiu  , &! INTENT(IN   )
       totlail  , &! INTENT(IN   )
       ayneetot , &! INTENT(IN   )
       aynpptot , &! INTENT(IN   )
       aygpptot , &! INTENT(IN   )
       aycsoi   , &! INTENT(IN   )
       aytrunoff, &! INTENT(IN   )
       garea    , &! INTENT(IN   )
       npoi     , &! INTENT(IN   )
       nVegClass  )! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    ! common blocks
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN   ) :: nVegClass
    INTEGER, INTENT(IN   ) :: npoi                ! total number of land points
    REAL(KIND=r8)   , INTENT(IN   ) :: garea    (npoi)     ! area of each gridcell (m**2)

    REAL(KIND=r8)   , INTENT(IN   ) :: ayneetot (npoi)     ! annual total NEE for ecosystem (kg-C/m**2/yr)
    REAL(KIND=r8)   , INTENT(IN   ) :: aynpptot (npoi)     ! annual total npp for ecosystem (kg-c/m**2/yr)
    REAL(KIND=r8)   , INTENT(IN   ) :: aygpptot (npoi)     ! annual total gpp for ecosystem (kg-c/m**2/yr)
    REAL(KIND=r8)   , INTENT(IN   ) :: aycsoi   (npoi)     ! total soil carbon (kg-c/m**2)
    REAL(KIND=r8)   , INTENT(IN   ) :: aytrunoff(npoi)     ! annual average total runoff (mm/yr)

    REAL(KIND=r8)   , INTENT(IN   ) :: vegtype0(npoi)      ! annual vegetation type - ibis classification
    REAL(KIND=r8)   , INTENT(IN   ) :: totbiou (npoi)      ! total biomass in the upper canopy (kg_C m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: totbiol (npoi)      ! total biomass in the lower canopy (kg_C m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: totlaiu (npoi)      ! total leaf area index for the upper canopy
    REAL(KIND=r8)   , INTENT(IN   ) :: totlail (npoi)      ! total leaf area index for the lower canopy

    !
    ! Arguments (input)
    !
    INTEGER, INTENT(IN   ) :: iyear       ! year counter
    INTEGER, INTENT(IN   ) :: iyear0      ! first year of simulation
    !
    ! local variables
    !
    INTEGER :: i        ! loop indices
    INTEGER :: k        ! loop indices
    !
    REAL(KIND=r8)    :: vtarea(nVegClass)    ! total area of the vegetation type (m**2)
    REAL(KIND=r8)    :: vnee(nVegClass)      ! vegetation type average nee (kg-c/m**2/yr)
    REAL(KIND=r8)    :: vnpp(nVegClass)      ! vegetation type average npp (kg-c/m**2/yr)
    REAL(KIND=r8)    :: vgpp(nVegClass)      ! vegetation type average gpp (kg-c/m**2/yr)
    REAL(KIND=r8)    :: vbiomass(nVegClass)  ! vegetation type average biomass (kg-c/m**2)
    REAL(KIND=r8)    :: vlai(nVegClass)      ! vegetation type average lai (m**2/m**2)
    REAL(KIND=r8)    :: vsoic(nVegClass)     ! vegetation type average soil carbon (kg-c/m**2)
    REAL(KIND=r8)    :: vrunoff(nVegClass)   ! vegetation type average runoff (mm/yr)
    !
    !
    ! initialize variables
    !
    DO k = 1, nVegClass
       !
       vtarea(k)   = 0.00_r8
       vnee(k)     = 0.00_r8
       vnpp(k)     = 0.00_r8
       vgpp(k)     = 0.00_r8
       vbiomass(k) = 0.00_r8
       vlai(k)     = 0.00_r8
       vsoic(k)    = 0.00_r8
       vrunoff(k)  = 0.00_r8
       !
    END DO
    !
    ! sum ecosystem properties over each vegetation type
    !
    DO i = 1, npoi
       !
       k = int (max (1.00_r8, min (15.00_r8, vegtype0(i))))
       !
       vtarea(k)   = vtarea(k)   + garea(i)
       !
       vnee(k)     = vnee(k)     + garea(i) * ayneetot(i)
       vnpp(k)     = vnpp(k)     + garea(i) * aynpptot(i)
       vgpp(k)     = vgpp(k)     + garea(i) * aygpptot(i)
       !
       vbiomass(k) = vbiomass(k) + garea(i) * totbiou(i)  &
            + garea(i) * totbiol(i)
       !
       vlai(k)     = vlai(k)     + garea(i) * totlaiu(i)   &
            + garea(i) * totlail(i)
       !
       vsoic(k)    = vsoic(k)    + garea(i) * aycsoi(i) 
       !
       vrunoff(k)  = vrunoff(k)  + garea(i) * aytrunoff(i)
       !
    END DO
    !
    ! calculate area averages
    !
    DO k = 1, nVegClass
       !
       vnee(k)     = vnee(k)     / max (1.00_r8, vtarea(k))
       vnpp(k)     = vnpp(k)     / max (1.00_r8, vtarea(k))
       vgpp(k)     = vgpp(k)     / max (1.00_r8, vtarea(k))
       vbiomass(k) = vbiomass(k) / max (1.00_r8, vtarea(k))
       vlai(k)     = vlai(k)     / max (1.00_r8, vtarea(k))
       vsoic(k)    = vsoic(k)    / max (1.00_r8, vtarea(k))
       vrunoff(k)  = vrunoff(k)  / max (1.00_r8, vtarea(k))
       !
    END DO
    !
    ! write some diagnostic output to history file
    !
    IF (iyear.eq.iyear0) THEN
       OPEN (30,file='ibis.out.vegtype',status='unknown')
    END IF
    !
    WRITE (30,*) ' '
    WRITE (30,*) '* * annual diagnostic fields by vegetation type * *'
    WRITE (30,*) ' '
    WRITE (30,*) &
         'year    veg           area       nee       npp       gpp   '// &
         'biomass       lai   scarbon    runoff '
    !
    DO k = 1, nVegClass
       !
       WRITE (30,9000) &
            iyear, k, vtarea(k) / 1.0e+06_r8, vnee(k), vnpp(k), vgpp(k), &
            vbiomass(k), vlai(k), vsoic(k), vrunoff(k)
       !
    END DO
    !
    CALL flush (30)
    !
    ! FORMAT statements
    !
9000 FORMAT (1x,i4,5x,i2,5x,1e10.3,7f10.3)
    !
    ! return to main program
    !
    RETURN
  END SUBROUTINE vdiag

  !
  !  ####   #          #    #    #    ##     #####  ######
  ! #    #  #          #    ##  ##   #  #      #    #
  ! #       #          #    # ## #  #    #     #    #####
  ! #       #          #    #    #  ######     #    #
  ! #    #  #          #    #    #  #    #     #    #
  !  ####   ######     #    #    #  #    #     #    ######
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE climanl2(TminL     , &! INTENT(IN   )
       TminU     , &! INTENT(IN   )
       Twarm     , &! INTENT(IN   )
       GDD       , &! INTENT(IN   )
       gdd0      , &! INTENT(INOUT)
       gdd0this  , &! INTENT(IN   )
       tc        , &! INTENT(INOUT)
       tw        , &! INTENT(INOUT)
       tcthis    , &! INTENT(IN   )
       twthis    , &! INTENT(IN   )
       tcmin     , &! INTENT(INOUT) local
       gdd5      , &! INTENT(INOUT) local
       gdd5this  , &! INTENT(IN   )
       exist     , &! INTENT(OUT  )
       deltat    , &! INTENT(IN   )
       npoi      , &! INTENT(IN   )
       npft        )! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    ! this subroutine updates the growing degree days, coldest temp, and
    ! warmest temp if monthly anomalies or daily values are used
    !
    ! common blocks
    !
    IMPLICIT NONE 
    !
    INTEGER, INTENT(IN   ) :: npoi                ! total number of land points
    INTEGER, INTENT(IN   ) :: npft                ! number of plant functional types
    REAL(KIND=r8), INTENT(IN   ) :: deltat  (npoi)      ! absolute minimum temperature -
    ! temp on average of coldest month (C)
    REAL(KIND=r8), INTENT(INOUT) :: gdd0    (npoi)      ! growing degree days > 0C 
    REAL(KIND=r8), INTENT(INOUT) :: gdd0this(npoi)      ! annual total growing degree days for current year
    REAL(KIND=r8), INTENT(INOUT) :: tc      (npoi)      ! coldest monthly temperature (C)
    REAL(KIND=r8), INTENT(INOUT) :: tw      (npoi)      ! warmest monthly temperature (C)
    REAL(KIND=r8), INTENT(INOUT) :: tcthis  (npoi)      ! coldest monthly temperature of current year (C)
    REAL(KIND=r8), INTENT(INOUT) :: twthis  (npoi)      ! warmest monthly temperature of current year (C)
    REAL(KIND=r8), INTENT(INOUT) :: tcmin   (npoi)      ! coldest daily temperature of current year (C)
    REAL(KIND=r8), INTENT(INOUT) :: gdd5    (npoi)      ! growing degree days > 5C
    REAL(KIND=r8), INTENT(INOUT) :: gdd5this(npoi)      ! annual total growing degree days for current year
    REAL(KIND=r8), INTENT(INOUT) :: exist   (npoi,npft) ! probability of existence of each plant functional type in a gridcell
    REAL(KIND=r8), INTENT(IN   ) :: TminL(npft)      ! Absolute minimum temperature -- lower limit (upper canopy PFTs)
    REAL(KIND=r8), INTENT(IN   ) :: TminU(npft)      ! Absolute minimum temperature -- upper limit (upper canopy PFTs)
    REAL(KIND=r8), INTENT(IN   ) :: Twarm(npft)      ! Temperature of warmest month (lower canopy PFTs)
    REAL(KIND=r8), INTENT(IN   ) :: GDD(npft)      ! minimum GDD needed (base 5 C for upper canopy PFTs, 
    ! base 0 C for lower canopy PFTs)

    ! 
    ! local variables
    !
    INTEGER :: i             ! loop indice
    !
    REAL(KIND=r8):: zweigc        ! 30-year e-folding time-avarage
    REAL(KIND=r8):: zweigw        ! 30-year e-folding time-avarage
    REAL(KIND=r8):: rworkc        ! 30-year e-folding time-avarage
    REAL(KIND=r8):: rworkw 
    !
    ! calculate a 30-year e-folding time-avarage
    !
    !      zweigc = exp(-1.0_r8/30.0_r8)
    !      zweigw = exp(-1.0_r8/30.0_r8)
    !
    !
    !     The filtering of the growing degree days and the climatic limits
    !     for existence of pft is done over 5 years instead of 30 in off
    !     -line IBIS.
    !
    !      zweigc = exp(-1.0_r8/30.0_r8)
    !      zweigw = exp(-1.0_r8/30.0_r8)
    zweigc = exp(-1.0_r8/5.0_r8)
    zweigw = exp(-1.0_r8/5.0_r8)

    rworkc = 1.0_r8 - zweigc
    rworkw = 1.0_r8 - zweigw
    !
    ! update critical climatic parameters with running average
    !
    DO  i = 1, npoi
       !
       tc(i) = zweigc * tc(i) + rworkc * tcthis(i)
       tw(i) = zweigw * tw(i) + rworkw * twthis(i)
       !
       tcmin(i) = tc(i) + deltat(i)
       !
       gdd0(i) = zweigc * gdd0(i) + rworkc * gdd0this(i)
       !
       gdd5(i) = zweigc * gdd5(i) + rworkc * gdd5this(i)
       !
       !
       !     Initialize this year's value of gdd0, gdd5, tc and tw to 0
       !     (climanl2 called 1st time step of the year, different from off
       !     -line IBIS)
       !
       tcthis   (i) =  100.0_r8
       twthis   (i) = - 100.0_r8
       gdd0this (i) = 0.0_r8
       gdd5this (i) = 0.0_r8

    END DO
    !
    CALL existence(TminL , &! INTENT(IN   )
         TminU , &! INTENT(IN   )
         Twarm , &! INTENT(IN   )
         GDD   , &! INTENT(IN   )
         exist , &! INTENT(OUT  )
         tcmin , &! INTENT(IN   )
         gdd5  , &! INTENT(IN   )
         gdd0  , &! INTENT(IN   )
         tw    , &! INTENT(IN   )
         npoi  , &! INTENT(IN   )
         npft    )! INTENT(IN   )
    !
    RETURN
  END SUBROUTINE climanl2
  !
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE existence(TminL    , &! INTENT(IN   )
       TminU    , &! INTENT(IN   )
       Twarm    , &! INTENT(IN   )
       GDD      , &! INTENT(IN   )
       exist    , &! INTENT(OUT  )
       tcmin    , &! INTENT(IN   )
       gdd5     , &! INTENT(IN   ) 
       gdd0     , &! INTENT(IN   ) 
       tw       , &! INTENT(IN   )
       npoi     , &! INTENT(IN   )
       npft       )! INTENT(IN   )
    ! ---------------------------------------------------------------------
    !
    ! this routine determines which plant functional types (pft's) are allowed
    ! to exist in each gridcell, based on a simple set of climatic criteria
    !
    ! the logic here is based on the biome3 model of haxeltine and prentice
    !
    ! plant functional types:
    !
    ! 1)  tropical broadleaf evergreen trees
    ! 2)  tropical broadleaf drought-deciduous trees
    ! 3)  warm-temperate broadleaf evergreen trees
    ! 4)  temperate conifer evergreen trees
    ! 5)  temperate broadleaf cold-deciduous trees
    ! 6)  boreal conifer evergreen trees
    ! 7)  boreal broadleaf cold-deciduous trees
    ! 8)  boreal conifer cold-deciduous trees
    ! 9)  evergreen shrubs
    ! 10) deciduous shrubs
    ! 11) warm (c4) grasses
    ! 12) cool (c3) grasses
    !
    !
    ! common blocks
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN   ) :: npoi            ! total number of land points
    INTEGER, INTENT(IN   ) :: npft            ! number of plant functional types

    REAL(KIND=r8), INTENT(INOUT) :: exist(npoi,npft)! probability of existence of each plant functional type in a gridcell
    REAL(KIND=r8), INTENT(IN   ) :: tcmin(npoi)     ! coldest daily temperature of current year (C)
    REAL(KIND=r8), INTENT(IN   ) :: gdd5 (npoi)     ! growing degree days > 5C
    REAL(KIND=r8), INTENT(IN   ) :: gdd0 (npoi)     ! growing degree days > 0C 
    REAL(KIND=r8), INTENT(IN   ) :: tw   (npoi)     ! warmest monthly temperature (C)
    REAL(KIND=r8), INTENT(IN   ) :: TminL(npft)     ! Absolute minimum temperature -- lower limit (upper canopy PFTs)
    REAL(KIND=r8), INTENT(IN   ) :: TminU(npft)     ! Absolute minimum temperature -- upper limit (upper canopy PFTs)
    REAL(KIND=r8), INTENT(IN   ) :: Twarm(npft)     ! Temperature of warmest month (lower canopy PFTs)
    REAL(KIND=r8), INTENT(IN   ) :: GDD  (npft)     ! minimum GDD needed (base 5 C for upper canopy PFTs, 
    ! base 0 C for lower canopy PFTs)

    !
    ! Local variables
    !
    INTEGER :: i      ! loop indice
    !
    ! ---------------------------------------------------------------------
    !
    DO i = 1, npoi
       !
       ! determine which plant types can exist in a given gridcell
       !
       exist(i,1)  = 0.0_r8
       exist(i,2)  = 0.0_r8
       exist(i,3)  = 0.0_r8
       exist(i,4)  = 0.0_r8
       exist(i,5)  = 0.0_r8
       exist(i,6)  = 0.0_r8
       exist(i,7)  = 0.0_r8
       exist(i,8)  = 0.0_r8
       exist(i,9)  = 0.0_r8
       exist(i,10) = 0.0_r8
       exist(i,11) = 0.0_r8
       exist(i,12) = 0.0_r8
       !
       ! 1) tropical broadleaf evergreen trees
       !
       !  - tcmin > 0.0
       !
       !        if (tcmin(i).gt.0.0)           exist(i,1) = 1.0
       !
       ! 2) tropical broadleaf drought-deciduous trees
       !
       !  - tcmin > 0.0
       !
       !        if (tcmin(i).gt.0.0)           exist(i,2) = 1.0
       !
       ! 3) warm-temperate broadleaf evergreen trees
       !
       !  - tcmin <   0.0 and
       !  - tcmin > -10.0
       !
       !        if ((tcmin(i).lt.0.0).and.
       !     >      (tcmin(i).gt.-10.0))       exist(i,3) = 1.0
       !
       ! 4) temperate conifer evergreen trees
       !
       !  - tcmin <    0.0 and
       !  - tcmin >  -45.0 and
       !  - gdd5  > 1200.0
       !
       !        if ((tcmin(i).lt.0.0).and.
       !     >      (tcmin(i).gt.-45.0).and.
       !     >      (gdd5(i).gt.1200.0))       exist(i,4) = 1.0
       !
       ! 5) temperate broadleaf cold-deciduous trees
       !
       !  - tcmin <    0.0 and
       !  - tcmin >  -45.0 and
       !  - gdd5  > 1200.0
       !
       !        if ((tcmin(i).lt.0.0).and.
       !     >      (tcmin(i).gt.-45.0).and.
       !     >      (gdd5(i).gt.1200.0))       exist(i,5) = 1.0
       !
       ! 6) boreal conifer evergreen trees
       !
       !  - tcmin <  -45.0 or gdd5 < 1200.0, and
       !  - tcmin >  -57.5 and
       !  - gdd5  >  350.0
       !
       !        if (((tcmin(i).lt.-45.0).or.(gdd5(i).lt.1200.0)).and.
       !     >       (tcmin(i).gt.-57.5).and.
       !     >       (gdd5(i).gt.350.0))       exist(i,6) = 1.0
       !
       ! 7) boreal broadleaf cold-deciduous trees
       !
       !  - tcmin <  -45.0 or gdd5 < 1200.0, and
       !  - tcmin >  -57.5 and
       !  - gdd5  >  350.0
       !
       !        if (((tcmin(i).lt.-45.0).or.(gdd5(i).lt.1200.0)).and.
       !     >       (tcmin(i).gt.-57.5).and.
       !     >       (gdd5(i).gt.350.0))       exist(i,7) = 1.0
       !
       ! 8) boreal conifer cold-deciduous trees
       !
       !  - tcmin <  -45.0 or gdd5 < 1200.0, and
       !  - gdd5  >  350.0
       !
       !        if (((tcmin(i).lt.-45.0).or.(gdd5(i).lt.1200.0)).and.
       !     >       (gdd5(i).gt.350.0))       exist(i,8) = 1.0
       !
       ! 9) evergreen shrubs
       !
       !  - gdd0 > 100.0
       !
       !        if (gdd0(i).gt.100.0)          exist(i,9) = 1.0
       !
       ! 10) deciduous shrubs
       !
       !  - gdd0 > 100.0
       !
       !        if (gdd0(i).gt.100.0)          exist(i,10) = 1.0
       !
       ! 11) warm (c4) grasses
       !
       !  - tw   >  22.0 and
       !  - gdd0 > 100.0
       !
       !        if ((tw(i).gt.22.0).and.
       !     >      (gdd0(i).gt.100.0))        exist(i,11) = 1.0
       !
       ! 12) cool (c3) grasses
       !
       !  - gdd0 > 100.0
       !
       !        if (gdd0(i).gt.100.0)          exist(i,12) = 1.0
       !
       !
       !*** DTP 2001/06/07: Modified version of above code reads in PFT
       !    existence criteria from external parameter file "params.veg"
       !    These are copied here for reference.... 
       !------------------------------------------------------------------
       !  TminL    TminU    Twarm    GDD    PFT
       !------------------------------------------------------------------
       !    0.0   9999.0   9999.0   9999  !   1
       !    0.0   9999.0   9999.0   9999  !   2
       !  -10.0      0.0   9999.0   9999  !   3
       !  -45.0      0.0   9999.0   1200  !   4
       !  -45.0      0.0   9999.0   1200  !   5
       !  -57.5    -45.0   9999.0    350  !   6
       !  -57.5    -45.0   9999.0    350  !   7
       ! 9999.0    -45.0   9999.0    350  !   8
       ! 9999.0   9999.0   9999.0    100  !   9
       ! 9999.0   9999.0   9999.0    100  !  10
       ! 9999.0   9999.0     22.0    100  !  11
       ! 9999.0   9999.0   9999.0    100  !  12
       !------------------------------------------------------------------

       ! 1) tropical broadleaf evergreen trees
       !
       !  - tcmin > 0.0
       !
       IF (tcmin(i).gt.TminL(1))      exist(i,1) = 1.00_r8
       !
       ! 2) tropical broadleaf drought-deciduous trees
       !
       !  - tcmin > 0.0
       !
       IF (tcmin(i).gt.TminL(2))      exist(i,2) = 1.00_r8
       !
       ! 3) warm-temperate broadleaf evergreen trees
       !
       !  - tcmin <   0.0 and
       !  - tcmin > -10.0
       !
       IF ((tcmin(i).lt.TminU(3)).and.  &
            (tcmin(i).gt.TminL(3)))    exist(i,3) = 1.00_r8
       !
       ! 4) temperate conifer evergreen trees
       !
       !  - tcmin <    0.0 and
       !  - tcmin >  -45.0 and
       !  - gdd5  > 1200.0
       !
       IF ((tcmin(i).lt.TminU(4)).and.   &
            (tcmin(i).gt.TminL(4)).and.   &
            (gdd5(i).gt.GDD(4)))       exist(i,4) = 1.00_r8
       !
       ! 5) temperate broadleaf cold-deciduous trees
       !
       !  - tcmin <    0.0 and
       !  - tcmin >  -45.0 and
       !  - gdd5  > 1200.0
       !
       IF ((tcmin(i).lt.TminU(5)).and.     &
            (tcmin(i).gt.TminL(5)).and.     & 
            (gdd5(i).gt.GDD(5)))       exist(i,5) = 1.00_r8
       !
       ! 6) boreal conifer evergreen trees
       !
       !  - tcmin <  -45.0 or gdd5 < 1200.0, and
       !  - tcmin >  -57.5 and
       !  - gdd5  >  350.0
       !
       IF (((tcmin(i).lt.TminU(6)).or.   &
            (gdd5(i).lt.GDD(4))).and.     &
            (tcmin(i).gt.TminL(6)).and.   &
            (gdd5(i).gt.GDD(6)))       exist(i,6) = 1.00_r8
       !
       ! 7) boreal broadleaf cold-deciduous trees
       !
       !  - tcmin <  -45.0 or gdd5 < 1200.0, and
       !  - tcmin >  -57.5 and
       !  - gdd5  >  350.0
       !
       IF (((tcmin(i).lt.TminU(7)).or.  &
            (gdd5(i).lt.GDD(5))).and.    & 
            (tcmin(i).gt.TminL(7)).and.  &
            (gdd5(i).gt.GDD(7)))       exist(i,7) = 1.00_r8
       !
       ! 8) boreal conifer cold-deciduous trees
       !
       !  - tcmin <  -45.0 or gdd5 < 1200.0, and
       !  - gdd5  >  350.0
       !
       IF (((tcmin(i).lt.TminU(8)).or.  &
            (gdd5(i).lt.TminL(4))).and.  &
            (gdd5(i).gt.GDD(8)))       exist(i,8) = 1.00_r8
       !
       ! 9) evergreen shrubs
       !
       !  - gdd0 > 100.0
       !
       IF (gdd0(i).gt.GDD(9))         exist(i,9) = 1.00_r8
       !
       ! 10) deciduous shrubs
       !
       !  - gdd0 > 100.0
       !
       IF (gdd0(i).gt.GDD(10))        exist(i,10) = 1.00_r8
       !
       ! 11) warm (c4) grasses
       !
       !  - tw   >  22.0 and
       !  - gdd0 > 100.0
       !
       IF ((tw(i).gt.Twarm(11)).and.  &
            (gdd0(i).gt.GDD(11)))      exist(i,11) = 1.00_r8
       !
       ! 12) cool (c3) grasses
       !
       !  - gdd0 > 100.0
       !
       IF (gdd0(i).gt.GDD(12))        exist(i,12) = 1.00_r8

    END DO
    !
    RETURN
  END SUBROUTINE existence


  !
  ! #####      #     ####    ####   ######   ####    ####   #    #  ######  #    #
  ! #    #     #    #    #  #    #  #       #    #  #    #  #    #  #       ##  ##
  ! #####      #    #    #  #       #####   #    #  #       ######  #####   # ## #
  ! #    #     #    #    #  #  ###  #       #    #  #       #    #  #       #    #
  ! #    #     #    #    #  #    #  #       #    #  #    #  #    #  #       #    #
  ! #####      #     ####    ####   ######   ####    ####   #    #  ######  #    #
  !
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE soilbgc (iyear    , &! INTENT(IN   )
       iyear0   , &! INTENT(IN   )
       imonth   , &! INTENT(IN   )
       iday     , &! INTENT(IN   )
       spin     , &! INTENT(IN   )
       spinmax  , &! INTENT(IN   )
       ayprcp   , &! INTENT(IN   )
       adfalll    , &! INTENT(IN   )
       adfallr    , &! INTENT(IN   )
       adfallw    , &! INTENT(IN   )
       falll    , &! INTENT(IN   )
       fallr    , &! INTENT(IN   )
       fallw    , &! INTENT(IN   )
       clitlm   , &! INTENT(INOUT)
       clitls   , &! INTENT(INOUT)
       clitrm   , &! INTENT(INOUT)
       clitrs   , &! INTENT(INOUT)
       clitwm   , &! INTENT(INOUT)
       clitws   , &! INTENT(INOUT)
       csoislop , &! INTENT(INOUT)
       csoislon , &! INTENT(INOUT)
       csoipas  , &! INTENT(INOUT)
       totcmic  , &! INTENT(INOUT)
       clitll   , &! INTENT(INOUT)
       clitrl   , &! INTENT(INOUT)
       clitwl   , &! INTENT(INOUT)
       decomps  , &! INTENT(IN   )
       decompl  , &! INTENT(IN   )
       tnmin    , &! INTENT(OUT  )
       totnmic  , &! INTENT(OUT  )
       totlit   , &! INTENT(OUT  )
       totalit  , &! INTENT(OUT  )
       totrlit  , &! INTENT(OUT  )
       totcsoi  , &! INTENT(OUT  )
       totfall  , &! INTENT(OUT  )
       totnlit  , &! INTENT(OUT  )
       totanlit , &! INTENT(OUT  )
       totrnlit , &! INTENT(OUT  )
       totnsoi  , &! INTENT(OUT  )
       tco2mic  , &! INTENT(OUT  )
       storedn  , &! INTENT(INOUT)
       yrleach  , &! INTENT(INOUT)
       ynleach  , &! INTENT(INOUT)
       ynleach_p ,&! INTENT(INOUT)
       tnmin_p   ,&! INTENT(OUT  )
       totnmic_p ,&! INTENT(OUT  )
       totnlit_p ,&! INTENT(OUT  )
       totanlit_p,&! INTENT(OUT  )
       totrnlit_p,&! INTENT(OUT  )
       totnsoi_p ,&! INTENT(OUT  )
       storedn_p ,&! INTENT(INOUT)
       csoi      ,&! INTENT(INOUT)  
       ta       , &! INTENT(IN   )
       vegtype0 , &! INTENT(IN   )
       hsoi     , &! INTENT(IN   )
       sand     , &! INTENT(IN   )
       clay     , &! INTENT(IN   )
       npoi     , &! INTENT(IN   )
       nsoilay  , &! INTENT(IN   )
       ndaypy     )! INTENT(IN   )
    ! --------------------------------------------------------------------------
    !
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: npoi    ! total number of land points
    INTEGER, INTENT(IN   ) :: nsoilay ! number of soil layers
    INTEGER, INTENT(IN   ) :: ndaypy  ! number of days per year
    REAL(KIND=r8), INTENT(IN   ) :: hsoi(npoi,nsoilay+1)  ! soil layer thickness (m)
    REAL(KIND=r8), INTENT(IN   ) :: sand(npoi,nsoilay)  ! percent sand of soil
    REAL(KIND=r8), INTENT(IN   ) :: clay(npoi,nsoilay)  ! percent clay of soil
    REAL(KIND=r8), INTENT(INOUT) :: storedn (npoi)    ! total storage of N in soil profile (kg_N m-2) 
    REAL(KIND=r8), INTENT(INOUT) :: yrleach (npoi)    ! annual total amount C leached from soil profile (kg_C m-2/yr)
    REAL(KIND=r8), INTENT(INOUT) :: ynleach (npoi)
    REAL(KIND=r8), INTENT(IN   ) :: adfalll(npoi)    ! day leaf litter fall                         (kg_C m-2/day)
    REAL(KIND=r8), INTENT(IN   ) :: adfallr(npoi)    ! day root litter input                        (kg_C m-2/day)
    REAL(KIND=r8), INTENT(IN   ) :: adfallw(npoi)    ! day wood litter fall                         (kg_C m-2/day)
    REAL(KIND=r8), INTENT(IN   ) :: falll   (npoi)   ! annual leaf litter fall      (kg_C m-2/year)
    REAL(KIND=r8), INTENT(IN   ) :: fallr   (npoi)   ! annual root litter input      (kg_C m-2/year)
    REAL(KIND=r8), INTENT(IN   ) :: fallw   (npoi)   ! annual wood litter fall      (kg_C m-2/year)
    REAL(KIND=r8), INTENT(INOUT) :: clitlm  (npoi)   ! carbon in leaf litter pool - metabolic       (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitls  (npoi)   ! carbon in leaf litter pool - structural      (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitrm  (npoi)   ! carbon in fine root litter pool - metabolic  (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitrs  (npoi)   ! carbon in fine root litter pool - structural (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitwm  (npoi)   ! carbon in woody litter pool - metabolic      (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitws  (npoi)   ! carbon in woody litter pool - structural     (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: csoislop(npoi)   ! carbon in soil - slow protected humus   (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: csoislon(npoi)   ! carbon in soil - slow nonprotected humus     (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: csoipas (npoi)   ! carbon in soil - passive humus    (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: totcmic (npoi)   ! total carbon residing in microbial pools     (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitll  (npoi)   ! carbon in leaf litter pool - lignin     (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitrl  (npoi)   ! carbon in fine root litter pool - lignin     (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitwl  (npoi)   ! carbon in woody litter pool - lignin     (kg_C m-2)
    REAL(KIND=r8), INTENT(IN   ) :: decomps (npoi)   ! soil organic matter decomposition factor     (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: decompl (npoi)   ! litter decomposition factor (dimensionless)
    REAL(KIND=r8), INTENT(OUT  ) :: tnmin   (npoi)   ! instantaneous nitrogen mineralization (kg_N m-2/timestep)
    REAL(KIND=r8), INTENT(OUT  ) :: totnmic (npoi)   ! total nitrogen residing in microbial pool (kg_N m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: totlit  (npoi)   ! total carbon in all litter pools (kg_C m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: totalit (npoi)   ! total standing aboveground litter (kg_C m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: totrlit (npoi)   ! total root litter carbon belowground (kg_C m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: totcsoi (npoi)   ! total carbon in all soil pools (kg_C m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: totfall (npoi)   ! total litterfall and root turnover (kg_C m-2/year)
    REAL(KIND=r8), INTENT(OUT  ) :: totnlit (npoi)   ! total nitrogen in all litter pools (kg_N m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: totanlit(npoi)   ! total standing aboveground nitrogen in litter (kg_N m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: totrnlit(npoi)   ! total root litter nitrogen belowground (kg_N m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: totnsoi (npoi)   ! total nitrogen in soil (kg_N m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: tco2mic (npoi)   ! instantaneous microbial co2 flux from soil (mol-CO2 / m-2 / second)
    REAL(KIND=r8), INTENT(IN   ) :: ayprcp  (npoi)   ! daily precitation (mm/day)

    REAL(KIND=r8), INTENT(INOUT) :: ynleach_p (npoi) ! annual total amount P leached from soil profile   (kg_P m-2/yr)
    REAL(KIND=r8), INTENT(OUT  ) :: tnmin_p   (npoi)   ! instantaneous phosphorus mineralization         (kg_P m-2/timestep)
    REAL(KIND=r8), INTENT(OUT  ) :: totnmic_p (npoi)   ! total phosphorus residing in microbial pool     (kg_P m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: totnlit_p (npoi)   ! total phosphorus in all litter pools            (kg_P m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: totanlit_p(npoi)   ! total standing aboveground phosphorus in litter (kg_P m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: totrnlit_p(npoi)   ! total root litter phosphorus belowground        (kg_P m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: totnsoi_p (npoi)   ! total phosphorus in soil                        (kg_P m-2)
    REAL(KIND=r8), INTENT(INOUT) :: storedn_p (npoi)   ! total storage of P in soil profile              (kg_P m-2) 
    REAL(KIND=r8), INTENT(INOUT) :: csoi     (npoi,nsoilay) ! global ! specific heat of soil, no pore spaces (J kg-1 deg-1)
    REAL(KIND=r8), INTENT(IN   ) :: ta       (npoi) ! INTENT(IN   )
    REAL(KIND=r8), INTENT(IN   ) :: vegtype0 (npoi) ! INTENT(IN   )
    !
    ! Arguments (input)
    !
    INTEGER, INTENT(IN   ) :: iday        ! day in month
    INTEGER, INTENT(IN   ) :: iyear     ! current year
    INTEGER, INTENT(IN   ) :: iyear0    ! initial year
    INTEGER, INTENT(IN   ) :: imonth    ! current month
    !      INTEGER, INTENT(IN   ) :: nspinsoil   ! year when soil carbon spinup stops
    INTEGER, INTENT(IN   ) :: spin        ! # of times soilbgc has been called in the current day
    INTEGER, INTENT(IN   ) :: spinmax     ! total # of times soilbgc is called per day (spinup)
    ! 
    ! local variables
    !
    INTEGER :: i            ! loop indice
    !
    REAL(KIND=r8) :: totts          ! 1/ndaypy
    REAL(KIND=r8) :: fracll         ! lignin fraction of leaves 
    REAL(KIND=r8) :: fracls         ! structural fraction of leaves 
    REAL(KIND=r8) :: fraclm         ! metabolic fraction of leaves 
    REAL(KIND=r8) :: fracrl         ! lignin fraction of roots
    REAL(KIND=r8) :: fracrs         ! structural fraction of roots 
    REAL(KIND=r8) :: fracrm         ! metabolic fraction of roots
    REAL(KIND=r8) :: fracwl         ! lignin fraction of wood
    REAL(KIND=r8) :: fracws         ! structural fraction of wood
    REAL(KIND=r8) :: fracwm         ! metabolic fraction of wood 

    REAL(KIND=r8) :: fracll_p         ! lignin fraction of leaves 
    REAL(KIND=r8) :: fracls_p         ! structural fraction of leaves 
    REAL(KIND=r8) :: fraclm_p         ! metabolic fraction of leaves 
    REAL(KIND=r8) :: fracrl_p         ! lignin fraction of roots
    REAL(KIND=r8) :: fracrs_p         ! structural fraction of roots 
    REAL(KIND=r8) :: fracrm_p         ! metabolic fraction of roots
    REAL(KIND=r8) :: fracwl_p         ! lignin fraction of wood
    REAL(KIND=r8) :: fracws_p         ! structural fraction of wood
    REAL(KIND=r8) :: fracwm_p         ! metabolic fraction of wood 

    REAL(KIND=r8) :: outclm(npoi)! c leaving leaf metabolic pool 
    REAL(KIND=r8) :: outcls(npoi)! c leaving leaf structural pool 
    REAL(KIND=r8) :: outcll(npoi)! c leaving leaf lignin pool
    REAL(KIND=r8) :: outcrm(npoi)! c leaving root metabolic pool 
    REAL(KIND=r8) :: outcrs(npoi)! c leaving root structural pool 
    REAL(KIND=r8) :: outcrl(npoi)! c leaving root lignin pool
    REAL(KIND=r8) :: outcwm(npoi)! c leaving woody metabolic carbon pool
    REAL(KIND=r8) :: outcws(npoi)! c leaving woody structural carbon pool
    REAL(KIND=r8) :: outcwl(npoi)! c leaving woody lignin carbon pool
    REAL(KIND=r8) :: outcsb(npoi)! flow of passive c to biomass
    REAL(KIND=r8) :: outcps(npoi)! flow of protected om to passive pool 
    REAL(KIND=r8) :: outcns(npoi)! flow of non-protected om to passive pool
    REAL(KIND=r8) :: outcnb(npoi)! flow of non-protected om to biomass 
    REAL(KIND=r8) :: outcpb(npoi)! flow of protected om to biomass
    REAL(KIND=r8) :: outcbp(npoi)! c leaving protected biomass pool  
    REAL(KIND=r8) :: outcbn(npoi)! c leaving non-protected biomass pool
    REAL(KIND=r8) :: totc  (npoi)! total c in soil
    !
    REAL(KIND=r8) :: dbdt    (npoi)   ! change of c in biomass pools with time 
    REAL(KIND=r8) :: dcndt   (npoi)   ! change of c in non-protected om with time
    REAL(KIND=r8) :: dcpdt   (npoi)   ! change of c in protected om with time
    REAL(KIND=r8) :: dcsdt   (npoi)   ! change of c in passive om with time
    REAL(KIND=r8) :: totmin  (npoi)   ! total nitrogen mineralization
    REAL(KIND=r8) :: totimm  (npoi)   ! total nitrogen immobilization 
    REAL(KIND=r8) :: netmin  (npoi)   ! net nitrogen mineralization

    REAL(KIND=r8) :: nbiors_p  (npoi)
    REAL(KIND=r8) :: nminrs_p  (npoi)
    REAL(KIND=r8) :: nbiols_p  (npoi)
    REAL(KIND=r8) :: nminls_p  (npoi)
    REAL(KIND=r8) :: nbiows_p  (npoi)
    REAL(KIND=r8) :: nminws_p  (npoi)
    REAL(KIND=r8) :: nbiowm_p  (npoi)
    REAL(KIND=r8) :: nminwm_p  (npoi)
    REAL(KIND=r8) :: nbiolm_p  (npoi)
    REAL(KIND=r8) :: nminlm_p  (npoi)
    REAL(KIND=r8) :: nbiorm_p  (npoi)
    REAL(KIND=r8) :: nminrm_p  (npoi)
    REAL(KIND=r8) :: nbioslon_p(npoi)
    REAL(KIND=r8) :: nminslon_p(npoi)
    REAL(KIND=r8) :: nbioslop_p(npoi)
    REAL(KIND=r8) :: nminslop_p(npoi)
    REAL(KIND=r8) :: nbiopas_p (npoi)
    REAL(KIND=r8) :: nminpas_p (npoi)
    REAL(KIND=r8) :: totimm_p  (npoi)! total phophorus immobilization 
    REAL(KIND=r8) :: totmin_p  (npoi)! total phophorus immobilization 
    REAL(KIND=r8) :: nrelps_p  (npoi)  
    REAL(KIND=r8) :: nrelns_p  (npoi) 
    REAL(KIND=r8) :: nrelbn_p  (npoi) 
    REAL(KIND=r8) :: nrelbp_p  (npoi) 
    REAL(KIND=r8) :: nrelll_p  (npoi)  
    REAL(KIND=r8) :: nrelrl_p  (npoi) 
    REAL(KIND=r8) :: nrelwl_p  (npoi) 
    REAL(KIND=r8) :: totnrel_p (npoi) 
    REAL(KIND=r8) :: netmin_p  (npoi)	! net nitrogen mineralization
    REAL(KIND=r8) :: nsoipas_p (npoi)
    REAL(KIND=r8) :: nlitlm_p  (npoi)
    REAL(KIND=r8) :: nlitls_p  (npoi)
    REAL(KIND=r8) :: nlitll_p  (npoi)
    REAL(KIND=r8) :: nlitrm_p  (npoi)
    REAL(KIND=r8) :: nlitrs_p  (npoi)
    REAL(KIND=r8) :: nlitrl_p  (npoi)
    REAL(KIND=r8) :: nlitwm_p  (npoi)
    REAL(KIND=r8) :: nlitws_p  (npoi) 
    REAL(KIND=r8) :: nlitwl_p  (npoi)



    REAL(KIND=r8) :: nbiors  (npoi)
    REAL(KIND=r8) :: nbiols  (npoi)
    REAL(KIND=r8) :: nbiows  (npoi)
    REAL(KIND=r8) :: nbiowm  (npoi)
    REAL(KIND=r8) :: nbiolm  (npoi)
    REAL(KIND=r8) :: nbiorm  (npoi)
    REAL(KIND=r8) :: nbioslon(npoi)
    REAL(KIND=r8) :: nbioslop(npoi)
    REAL(KIND=r8) :: nbiopas (npoi)
    REAL(KIND=r8) :: nminrs  (npoi)
    REAL(KIND=r8) :: nminls  (npoi)
    REAL(KIND=r8) :: nminws  (npoi)
    REAL(KIND=r8) :: nminwm  (npoi)
    REAL(KIND=r8) :: nminlm  (npoi)
    REAL(KIND=r8) :: nminrm  (npoi)
    REAL(KIND=r8) :: nminslon(npoi)
    REAL(KIND=r8) :: nminslop(npoi)
    REAL(KIND=r8) :: nminpas (npoi)
    REAL(KIND=r8) :: nrelps  (npoi)
    REAL(KIND=r8) :: nrelns  (npoi)
    REAL(KIND=r8) :: nrelbn  (npoi)
    REAL(KIND=r8) :: nrelbp  (npoi)
    REAL(KIND=r8) :: nrelll  (npoi)
    REAL(KIND=r8) :: nrelrl  (npoi)
    REAL(KIND=r8) :: nrelwl  (npoi)
    REAL(KIND=r8) :: totnrel (npoi)
    REAL(KIND=r8) :: ymintot (npoi)
    REAL(KIND=r8) :: yminmic (npoi)
    !
    ! nitrogen in litter and soil pools
    !
    REAL(KIND=r8) :: nlitlm  (npoi)
    REAL(KIND=r8) :: nlitls  (npoi)
    REAL(KIND=r8) :: nlitll  (npoi)
    REAL(KIND=r8) :: nlitrm  (npoi)
    REAL(KIND=r8) :: nlitrs  (npoi)
    REAL(KIND=r8) :: nlitrl  (npoi)
    REAL(KIND=r8) :: nlitwm  (npoi)
    REAL(KIND=r8) :: nlitws  (npoi)
    REAL(KIND=r8) :: nlitwl  (npoi)
    REAL(KIND=r8) :: nsoislop(npoi)
    REAL(KIND=r8) :: nsoipas (npoi)
    REAL(KIND=r8) :: nsoislon(npoi)
    REAL(KIND=r8) :: nsoislon_p(npoi)
    REAL(KIND=r8) :: nsoislop_p(npoi)

    !
    ! variables controlling constraints on microbial biomass 
    !
    REAL(KIND=r8) :: cmicn (npoi)
    REAL(KIND=r8) :: cmicp (npoi)
    REAL(KIND=r8) :: cmicmx(npoi)
    !
    ! variables controlling leaching, calculating co2 respiration and n deposition
    !
    REAL(KIND=r8) :: cleach   (npoi)
    REAL(KIND=r8) :: totcbegin(npoi)
    REAL(KIND=r8) :: totcend  (npoi)
    REAL(KIND=r8) :: totcin   (npoi)
    REAL(KIND=r8) :: fixsoin  (npoi)
    REAL(KIND=r8) :: deposn   (npoi)
    REAL(KIND=r8) :: depth (npoi)
    REAL(KIND=r8) :: depth2(npoi)
    REAL(KIND=r8) :: rhosoi(npoi,nsoilay)
    REAL(KIND=r8) :: zdepth
    !
    REAL(KIND=r8) :: fleach
    REAL(KIND=r8) :: h20
    !
    ! decay constants for c pools
    !
    REAL(KIND=r8) :: klm          ! leaf metabolic litter 
    REAL(KIND=r8) :: kls      ! leaf structural litter
    REAL(KIND=r8) :: kll      ! leaf lignin
    REAL(KIND=r8) :: krm      ! root metabolic litter
    REAL(KIND=r8) :: krs      ! root structural litter
    REAL(KIND=r8) :: krl      ! root lignin
    REAL(KIND=r8) :: kwm      ! woody metabolic litter
    REAL(KIND=r8) :: kws      ! woody structural litter
    REAL(KIND=r8) :: kwl      ! wood  lignin
    REAL(KIND=r8) :: kbn      ! microbial biomass --> nonprotected om 
    REAL(KIND=r8) :: kbp      ! microbial biomass --> protected om
    REAL(KIND=r8) :: knb      ! nonprotected om   --> biomass
    REAL(KIND=r8) :: kns      ! nonprotected om   --> passive c 
    REAL(KIND=r8) :: kpb      ! protected om  --> biomass
    REAL(KIND=r8) :: kps      ! protected om  --> passive c
    REAL(KIND=r8) :: ksb          ! passive c  --> biomass
    !
    ! efficiencies for microbial decomposition
    !
    REAL(KIND=r8) :: ylm       ! leaf metabolic litter decomposition 
    REAL(KIND=r8) :: yls       ! leaf structural litter decomposition
    REAL(KIND=r8) :: yll       ! leaf lignin
    REAL(KIND=r8) :: yrm       ! root metabolic litter decomposition
    REAL(KIND=r8) :: yrs       ! root structural litter decomposition
    REAL(KIND=r8) :: yrl       ! root lignin
    REAL(KIND=r8) :: ywm       ! woody metabolic litter decomposition
    REAL(KIND=r8) :: yws       ! woody structural litter decomposition
    REAL(KIND=r8) :: ywl       ! wood lignin
    REAL(KIND=r8) :: ybn       ! microbial biomass to nonprotected om
    REAL(KIND=r8) :: ybp       ! microbial biomass to protected om
    REAL(KIND=r8) :: ynb       ! nonprotected om to biomass
    REAL(KIND=r8) :: yns       ! nonprotected om to passive  c
    REAL(KIND=r8) :: ypb       ! protected om to biomass
    REAL(KIND=r8) :: yps       ! protected om to passive c
    REAL(KIND=r8) :: ysb       ! passive c to biomass
    !
    REAL(KIND=r8) :: cnr(10)   ! c:n ratios of c and litter pools
    REAL(KIND=r8) :: cnrf(10)

    REAL(KIND=r8) :: cpr(10)   ! c:p ratios of c and litter pools
    REAL(KIND=r8) :: cprf(10)

    !
    ! constants for calculating fraction of litterall in structural
    ! metabolic and lignified (resistant) fractions
    !
    REAL(KIND=r8) :: cnleaf     ! input c:n ratio of leaf litterfall 
    REAL(KIND=r8) :: cnwood     ! input c:n ratio of wood litter
    REAL(KIND=r8) :: cnroot     ! input c:n ratio of root litter turnover

    REAL(KIND=r8) :: cpleaf     ! input c:p ratio of leaf litterfall 
    REAL(KIND=r8) :: cpwood     ! input c:p ratio of wood litter
    REAL(KIND=r8) :: cproot     ! input c:p ratio of root litter turnover

    REAL(KIND=r8) :: rconst     ! value set to 1200.  from Verberne model 
    REAL(KIND=r8) :: fmax       ! maximum fraction allowed in metabolic pool
    ! 
    ! variables added to do daily time series of some values
    !
    INTEGER :: gridpt
    INTEGER :: kk,k
    !
    ! variables dealing with soil texture and algorithms
    !
    INTEGER :: msand
    INTEGER :: mclay
    INTEGER :: isoil
    !
    REAL(KIND=r8) :: fsand
    REAL(KIND=r8) :: fclay
    REAL(KIND=r8) :: cfrac
    REAL(KIND=r8) :: texfact
    REAL(KIND=r8) :: fbpom
    REAL(KIND=r8) :: fbsom
    REAL(KIND=r8) :: rdepth
    REAL(KIND=r8) :: effac
    REAL(KIND=r8) :: lig_frac
    REAL(KIND=r8) :: om_frac
    REAL(KIND=r8) :: forganic
    REAL(r8) :: om_csol      = 2.5_r8   ! heat capacity of peat soil *10^6 (J/K m3) (Farouki, 1986)
    !
    !      textcls = 1      ! sand
    !      textcls = 2      ! loamy sand
    !      textcls = 3      ! sandy loam
    !      textcls = 4      ! loam
    !      textcls = 5      ! silt loam
    !      textcls = 6      ! sandy clay loam
    !      textcls = 7      ! clay loam
    !      textcls = 8      ! silty clay loam
    !      textcls = 9      ! sandy clay
    !      textcls = 10     ! silty clay
    !      textcls = 11     ! clay
    !  Input of phosphorus in soil due the weathering of rockets
    REAL(KIND=r8) :: INput_P(1:11)=(/0.00005_r8,0.00005_r8,0.00001_r8,0.000005_r8,0.00001_r8,0.00001_r8,&
                                       0.000005_r8,0.00001_r8,0.00001_r8,0.00001_r8,0.000003_r8/) 

    gridpt = npoi   ! total number of gridpoints used
    !
    ! total timesteps (daily) used to divide litterfall into daily fractions 
    !
    totts=1.0_r8/real(ndaypy,kind=r8)
    !
    ! -------------------------------------------------------------------------------------
    ! specific maximum decay rate or growth constants; rates are per day
    ! constants are taken from Parton et al., 1987 and Verberne et al., 1990
    ! and special issue of Geoderma (comparison of 9 organic matter models) in Dec. 1997
    !
    ! leaching parameterization was changed to agree with field data,
    ! this caused a changing of the below constants.  
    !
    ! approximate factors for Verberne et al. model where efficiencies are 100%
    ! for some of the transformations: one problem was that their rate constants were
    ! based on 25C, and our modifying functions are based on 15 C...thus the rate constants
    ! are somewhat smaller compared to the Verberne et al. (1990) model parameters
    ! rates are based on a daily decomposition timestep (per day)
    ! -------------------------------------------------------------------------------------
    !
    ! leaf litter constants
    !
    klm = 0.15_r8 !dpm leaf --> microbial biomass
    kls = 0.01_r8 !spm leaf --> microbial biomass
    kll = 0.01_r8!rpm leaf --> non or protected om
    !
    ! root litter constants
    !
    krm = 0.10_r8!dpm root --> microbial biomass
    krs = 0.005_r8 !spm root --> microbial biomass
    krl = 0.005_r8!rpm root --> non or protected om 
    !
    ! woody litter constants
    !
    kwm = 0.001_r8!dpm wood --> microbial biomass
    kws = 0.001_r8 !spm wood --> microbial biomass
    kwl = 0.001_r8!rpm wood --> non or protected om 
    !
    ! biomass constants
    !
    kbn = 0.045_r8!biomass --> non protected organic matter 
    kbp = 0.005_r8!biomass --> protected organic matter
    !
    ! slow and passive c pools
    !
    knb = 0.001_r8!non protected om --> biomass
    kns = 0.000001_r8!non protected om --> stablized om
    kpb = 0.0001_r8 !protected om     --> biomass
    kps = 0.000001_r8!protected om     --> stablized om
    ksb = 8.0e-07_r8!stablized om     --> biomass
    !
    ! ---------------------------------------------------------------------
    !  yield (efficiency) with which microbes gain biomass from c source
    !  the rest is driven off as co2 respiration (microbial respiration)
    !  all of the respiration produced by microbes is assumed to leave
    !  the soil profile over the course of a year
    !  taken primarily from the models of Verberne and CENTURY
    ! ---------------------------------------------------------------------
    !
    ylm = 0.4_r8       ! metabolic material efficiencies
    yrm = 0.4_r8
    ywm = 0.4_r8
    yls = 0.3_r8       ! structural efficiencies
    yrs = 0.3_r8
    yws = 0.3_r8
    !
    yll = 1.0_r8       ! resistant fraction
    yrl = 1.0_r8 
    ywl = 1.0_r8 
    ybn = 1.0_r8       ! biomass       --> non-protected pool
    ybp = 1.0_r8       ! biomass       --> protected pool
    yps = 1.0_r8       ! protected     --> passive
    yns = 1.0_r8       ! non-protected --> passive
    !
    ysb = 0.20_r8       ! passive pool  --> biomass
    ypb = 0.20_r8       ! protected     --> biomass
    ynb = 0.25_r8       ! non-protected --> biomass
    !
    ! -------------------------------------------------------------------
    ! split of lignified litter material between protected/non-protected
    ! slow OM pools
    ! -------------------------------------------------------------------
    !
    lig_frac = 0.50_r8 
    !
    ! -------------------------------------------------------------------
    ! protected biomass as a fraction of total soil organic carbon
    ! from Verberne et al., 1990
    ! -------------------------------------------------------------------
    !
    fbsom = 0.017_r8
    !
    ! ---------------------------------------------------------------------
    ! (effac) --> efficiency of microbial biomass reincorporated
    ! into biomass pool.(from NCSOIL parameterizations; Molina et al., 1983)
    ! ---------------------------------------------------------------------
    !
    effac = 0.40_r8 
    !
    ! ---------------------------------------------------------------------
    ! define C:N ratios of substrate pools and biomass
    ! metabolic, structural, and lignin are for Leaves and roots
    ! values from Parton et al., 1987 and Whitmore and Parry, 1988
    ! index: 1 - biomass, 2 - passive pool, 3- slow protected c,
    ! 4 - slow carbon, non-protected, 5 - resistant, 6 - structural plant
    ! leaf and root litter, 7 - metabolic plant and root litter, 
    ! 8- woody biomass
    ! ---------------------------------------------------------------------
    !
    cnr(1)  = 8.0_r8       !c:n ratio of microbial biomass
    cnr(2)  = 15.0_r8      !c:n ratio of passive soil carbon
    cnr(3)  = 10.0_r8      !c:n ratio of protected slow soil carbon
    cnr(4)  = 15.0_r8      !c:n ratio of non-protected slow soil C
    cnr(5)  = 100.0_r8     !c:n ratio of resistant litter lignin
    cnr(6)  = 150.0_r8     !c:n ratio of structural plant litter
    cnr(7)  = 6.0_r8       !c:n ratio of metabolic plant litter
    cnr(8)  = 250.0_r8     !c:n Ratio of woody components

!
! ---------------------------------------------------------------------
! define C:P ratios of substrate pools and biomass
! metabolic, structural, and lignin are for Leaves and roots
! values from Parton et al., 1987 and Whitmore and Parry, 1988
! index: 1 - biomass, 2 - passive pool, 3- slow protected c,
! 4 - slow carbon, non-protected, 5 - resistant, 6 - structural plant
! leaf and root litter, 7 - metabolic plant and root litter, 
! 8- woody biomass
! ---------------------------------------------------------------------
!
       cpr(1)  = 32.0_r8       !c:p ratio of microbial biomass
       cpr(2)  = 400.0_r8      !c:p ratio of passive soil carbon
       cpr(3)  = 465.0_r8      !c:p ratio of protected slow soil carbon
       cpr(4)  = 550.0_r8      !c:p ratio of non-protected slow soil C
       cpr(5)  = 3750.0_r8     !c:p ratio of resistant litter lignin
       cpr(6)  = 3650.0_r8     !c:p ratio of structural plant litter
       cpr(7)  = 10000.0_r8    !c:p ratio of metabolic plant litter
       cpr(8)  = 7600.0_r8     !c:p Ratio of woody components
    !
    ! ---------------------------------------------------------------------
    ! calculate the fraction of wood, roots and leaves that are structural,
    ! decomposable, and resistant based on equations presented in Verberne
    ! model discussion (Geoderma, December 1997 special issue).  fmax is the
    ! maximum fraction allowed in resistant fraction, rconst is a constant
    ! defined as 1200.  The cnratio of each plant part has to be less than
    ! the value of structural defined above (i.e. 150) otherwise the equations
    ! are unstable...thus the wood litter pool value for cnr(6) is substituted
    ! with a value higher than that for cnwood (i.e. 250).  this is 
    ! insignificant for wood since 97% is structural anyways.
    !
    ! ** NOTE ******** 
    ! Would like to incorporate different C:N ratios of residue/roots for
    ! different biome types based on literature search
    ! average c:n ratio would be based on litter inputs from each pft
    ! ****************
    ! ---------------------------------------------------------------------
    !
    ! equations were changed on 1-26-99 for erratum in literature (Whitmore
    ! et al. 1997) which had an error in equations to split litterfall into
    ! the correct three fractions
    ! 
    fmax   = 0.45_r8
    rconst = 1200.0_r8

    cnleaf = 40.0_r8      ! average c:n ratio for leaf litterfall
    cnroot = 60.0_r8      ! average c:n ratio for root turnover
    cnwood = 200.0_r8     ! average c:n ratio for woody debris

    cpleaf = 408.0_r8      ! average c:p ratio for leaf litterfall
    cproot = 1170.0_r8     ! average c:p ratio for root turnover
    cpwood = 3750.0_r8     ! average c:p ratio for woody debris

    !
    ! leaf litter   [nitrogen]
    !
    fracll = fmax * (cnleaf**2)/(rconst + cnleaf**2)
    fracls = (1.0_r8/cnleaf - fracll/cnr(5) - (1.0_r8-fracll)/cnr(7))/  &
         (1.0_r8/cnr(6) - 1.0_r8/cnr(7))
    fraclm = 1.0_r8 - fracll - fracls

    !
    ! leaf litter  [phosphorus]
    !
    fracll_p = fmax * (cpleaf**2)/(rconst + cpleaf**2)                          ! lignin fraction of leaves 
    fracls_p = (1.0_r8/cpleaf - fracll_p/cpr(5) - (1.0_r8-fracll_p)/cpr(7))/  & !  structural fraction of leaves 
               (1.0_r8/cpr(6) - 1.0_r8/cpr(7))
      
    fraclm_p = 1.0_r8 - fracll_p - fracls_p                                      ! metabolic fraction of leaves 

    !
    ! root litter [nitrogen]
    !
    fracrl = fmax * (cnroot**2)/(rconst + cnroot**2)
    fracrs = (1.0_r8/cnroot - fracrl/cnr(5) - (1.0_r8-fracrl)/cnr(7))/  &
         (1.0_r8/cnr(6) - 1.0_r8/cnr(7))
    fracrm = 1.0_r8 - fracrl - fracrs

    !
    ! root litter  [phosphorus]
    !
    fracrl_p = fmax * (cproot**2)/(rconst + cproot**2)                          ! lignin fraction of roots
    fracrs_p = (1.0_r8/cproot - fracrl_p/cpr(5) - (1.0_r8-fracrl_p)/cpr(7))/  & ! structural fraction of roots 
               (1.0_r8/cpr(6) - 1.0_r8/cpr(7))

    fracrm_p = 1.0_r8 - fracrl_p - fracrs_p                                      !  metabolic fraction of roots

    !
    ! wood litter   [nitrogen]
    !
    fracwl = fmax * (cnwood**2)/(rconst + cnwood**2)
    fracws = (1.0_r8/cnwood - fracwl/cnr(5) - (1.0_r8-fracwl)/cnr(7))/  &
         (1.0_r8/cnr(8) - 1.0_r8/cnr(7))
    fracwm = 1.0_r8 - fracwl - fracws
    !
    ! wood litter [phosphorus]
    !
    fracwl_p = fmax * (cpwood**2)/(rconst + cpwood**2)                      !lignin fraction of wood
    fracws_p = (1.0_r8/cpwood - fracwl_p/cpr(5) - (1.0_r8-fracwl_p)/cpr(7))/  & ! structural fraction of wood
               (1.0_r8/cpr(8) - 1.0_r8/cpr(7)) 
    fracwm_p = 1.0_r8 - fracwl_p - fracws_p                                     ! metabolic fraction of wood 
    !
    ! ------------------------------------------------------------------------
    ! calculate the efficiency of decomposition of the material based on the
    ! C/N ratio, and the approach outlined in Modeling Plant and Soil Systems
    ! eds. Hanks and Ritchie, 1991 Article by Goodwin and Jones
    ! called the C/N ratio factor (CNRF) limit between 0.01 - 1.0
    ! this is a 3rd modifying factor to the rate of decomposition, besides
    ! the controlling factors of temperature and moisture.  Reasoning for
    ! this is to be able to account for changing C:N ratios while the models
    ! are running a simulational; although the model does not do this yet.
    ! ------------------------------------------------------------------------
    !
    ! commented out because of confusion of effect on decomposition constants
    ! just set equal to 1.0
    !
    DO i = 1,8
       !
       !        cnrf(i) = exp(-0.693 * (cnr(i) - 25.0)/25.0)
       !        if (cnrf(i) .gt. 1.0) cnrf(i) = 1.0
       !        if (cnrf(i) .le. 0.01) cnrf(i) = 0.01
       !
       cnrf(i) = 1.00_r8
       !
    END DO

    !
    ! ------------------------------------------------------------------------
    ! calculate the efficiency of decomposition of the material based on the
    ! C/P ratio, and the approach outlined in Modeling Plant and Soil Systems
    ! eds. Hanks and Ritchie, 1991 Article by Goodwin and Jones
    ! called the C/P ratio factor (CPRF) limit between 0.01 - 1.0
    ! this is a 3rd modifying factor to the rate of decomposition, besides
    ! the controlling factors of temperature and moisture.  Reasoning for
    ! this is to be able to account for changing C:P ratios while the models
    ! are running a simulational; although the model does not do this yet.
    ! ------------------------------------------------------------------------
    !
    ! commented out because of confusion of effect on decomposition constants
    ! just set equal to 1.0
    !
    DO i = 1,8
    !
    !        cnrf(i) = exp(-0.693 * (cnr(i) - 25.0)/25.0)
    !        if (cnrf(i) .gt. 1.0) cnrf(i) = 1.0
    !        if (cnrf(i) .le. 0.01) cnrf(i) = 0.01
    !
         cprf(i) = 1.00_r8   !?
    !
    END DO

    depth =0.0_r8
    depth2=0.0_r8
    DO  k = 1, nsoilay
       DO i = 1, npoi
          IF(depth(i)<=1.0_r8)THEN!m
            depth (i)=depth (i)+hsoi(i,k)
          END IF
          IF(depth2(i)<=0.30_r8)THEN!m
            depth2(i)=depth2(i)+hsoi(i,k)
          END IF
       END DO
    END DO

    !
    DO i = 1, npoi
       !
       ! ---------------------------------------------------------------------
       ! fraction of decomposing microbial biomass into protected organic
       ! matter; taken from the model of Verberne et al., 1990
       ! this is the proportion of decomposing dead microbial biomass that
       ! is transferred to a protected pool vs. a non-protected pool
       ! related to the clay content of the soil. in sandy soils, fbpom = 0.3,
       ! whereas in clay soils fbpom = 0.7.  created a linear function based
       ! on clay fraction of soil to adjust according to amount of clay in
       ! the top 1 m of profile (weighted average according to depth of each
       ! layer)
       !
       ! also take care of calculation of texfact, which is a leaching
       ! parameter based on the average sand fraction of the top 1 m of
       ! soil
       ! ---------------------------------------------------------------------
       !
       !PK         rdepth   = 1.0_r8/(hsoi(i,1) + hsoi(i,2) + hsoi(i,3) + hsoi(i,4))
       rdepth   = 1.0_r8/(depth (i))
       cfrac    = 0.0_r8
       texfact  = 0.0_r8 
       !
        zdepth=0.0_r8
        DO kk = 1, nsoilay                  ! top 1 m of soil -- 4 layers
           zdepth=zdepth+hsoi(i,kk)
           IF(zdepth<=1.0_r8)THEN
              msand    = nint(sand(i,kk)) 
              mclay    = nint(clay(i,kk)) 
              fclay    = 0.01_r8 * mclay
              fsand    = 0.01_r8 * msand 
              cfrac    = cfrac   + fclay * hsoi(i,kk)
              texfact  = texfact + fsand * hsoi(i,kk)
           END IF
        END DO
       !
       cfrac   = cfrac   * rdepth
       texfact = texfact * rdepth
       !
       ! if cfrac is greater than 0.4, set fbpom = 0.7, if cfrac is less
       ! than 0.17, set fbpom = 0.30 (sandy soil)
       !
       !        fbpom = min(max(0.3, cfrac/0.4 * 0.7),0.7)      
       fbpom = 0.50_r8
       !
       ! ------------------------------------------------------------------------
       ! total soil carbon initialized to 0 at beginning of model run
       ! used in calculation of soil co2 respiration from microbial decomposition 
       ! ------------------------------------------------------------------------
       !
       IF (iday .eq. 1 .and. imonth .eq. 1 .and. iyear .eq. iyear0) THEN
          totcbegin(i) = 0.0_r8
          storedn(i)   = 0.0_r8 
          storedn_p(i)   = 0.0_r8 
       END IF
       !
       ! ------------------------------------------------------------------------
       ! initialize yearly summation of net mineralization and co2 respiration
       ! to 0 at beginning of each year; because these quantities are usually 
       ! reported on a yearly basis, we wish to do the same in the model so we
       ! can compare easily with the data.
       ! ------------------------------------------------------------------------
       !
       IF (iday .eq. 1 .and. imonth .eq. 1) THEN
          yrleach(i) = 0.0_r8
          cleach(i)  = 0.0_r8
          ynleach(i) = 0.0_r8
          ynleach_p(i) = 0.0_r8
          ymintot(i) = 0.0_r8
          yminmic(i) = 0.0_r8
       END IF
       !       
       ! determine amount of substrate available to microbial growth
       !
       ! total timesteps (daily) used to divide litterfall into daily fractions 
       !
       !       totts=1.0_r8/float(ndaypy)
       !
       ! calculate the total amount of litterfall entering soil(C)
       !
       !               kg_C )           1               kg C
       !              -------    *   ----------    =  --------
       !               m2 * year       365 day         m2 * day
       !
       !       totcin(i) =  falll(i)*totts + fallr(i)*totts  + fallw(i)*totts
       !
       !               kg_C )            kg C
       !              -------       =  --------
       !               m2 * day         m2 * day

       totcin(i) =  adfalll(i) + adfallr(i)  + adfallw(i)

       !
       ! calculate the current total amount of carbon at each grid cell
       !             (kg_C m-2)
       totc(i) = clitlm(i) + clitls(i) + clitrm(i) + clitrs(i) +  &
            clitwm(i) + clitws(i) + csoislop(i) + csoislon(i) +  &
            csoipas(i) + totcmic(i) + clitll(i) + clitrl(i) + clitwl(i)
       !
       ! beginning amount of soil C at each timestep (used for respiration
       ! calculation)
       !
       totcbegin(i) = totc(i)
       !
       ! ------------------------------------------------------------------------
       ! split current amount of total soil microbes
       ! maximum amount of biomass is a function of the total soil C
       ! from Verberne et al., 1990
       !
       ! protected biomass as a fraction of total soil organic carbon
       ! from Verberne et al., 1990
       !
       !      fbsom = 0.017_r8
       !
       ! ------------------------------------------------------------------------
       !
       !      totcmic(i) = cmicp(i) + cmicn(i)
       cmicmx(i) = fbsom * totc(i) 
       !
       ! calculate the amount of protected and unprotected biomass
       !
       IF (totcmic(i) .ge. cmicmx(i)) THEN
          !
          cmicp(i) = cmicmx(i)
          cmicn(i) = totcmic(i) - cmicmx(i)
          !
       ELSE
          !
          cmicn(i) = 0.0_r8
          cmicp(i) = totcmic(i)
          !
       END IF
       !
       ! ---------------------------------------------------------------
       ! litter pools 
       !
       ! add in the amount of litterfall, and root turnover
       ! ---------------------------------------------------------------
       !
       msand    = nint(sand(i,1)) 
       mclay    = nint(clay(i,1)) 

       isoil =   textcls (msand,mclay)

!               kg_C )           1               kg C
!              -------    *   ----------    =  --------
!               m2 * year       365 day         m2 * day

!       clitlm(i) = clitlm(i) + (fraclm * falll(i)*totts)  + (fraclm_p * falll(i)*totts)  ! carbon in leaf litter pool - metabolic       (kg_C m-2)
!       clitls(i) = clitls(i) + (fracls * falll(i)*totts)  + (fracls_p * falll(i)*totts)  ! carbon in leaf litter pool - structural      (kg_C m-2)
!       clitll(i) = clitll(i) + (fracll * falll(i)*totts)  + (fracll_p * falll(i)*totts)  ! carbon in leaf litter pool - lignin            (kg_C m-2)

!       clitrm(i) = clitrm(i) + (fracrm * fallr(i)*totts)  + (fracrm_p * fallr(i)*totts) + (INput_P(isoil)*totts) ! carbon in fine root litter pool - metabolic  (kg_C m-2)
!       clitrs(i) = clitrs(i) + (fracrs * fallr(i)*totts)  + (fracrs_p * fallr(i)*totts)  ! carbon in fine root litter pool - structural (kg_C m-2)
!       clitrl(i) = clitrl(i) + (fracrl * fallr(i)*totts)  + (fracrl_p * fallr(i)*totts)  ! carbon in fine root litter pool - lignin     (kg_C m-2)

!       clitwm(i) = clitwm(i) + (fracwm * fallw(i)*totts)  + (fracwm_p * fallw(i)*totts)  ! carbon in woody litter pool - metabolic      (kg_C m-2)
!       clitws(i) = clitws(i) + (fracws * fallw(i)*totts)  + (fracws_p * fallw(i)*totts)  ! carbon in woody litter pool - structural     (kg_C m-2)
!       clitwl(i) = clitwl(i) + (fracwl * fallw(i)*totts)  + (fracwl_p * fallw(i)*totts)  ! carbon in woody litter pool - lignin	     (kg_C m-2)


!               kg_C )                kg C
!              -------    *      =  --------
!               m2 * year            m2 * day

       clitlm(i) = clitlm(i) + (fraclm * adfalll(i))  + (fraclm_p * adfalll(i))  ! carbon in leaf litter pool - metabolic       (kg_C m-2)
       clitls(i) = clitls(i) + (fracls * adfalll(i))  + (fracls_p * adfalll(i))  ! carbon in leaf litter pool - structural      (kg_C m-2)
       clitll(i) = clitll(i) + (fracll * adfalll(i))  + (fracll_p * adfalll(i))  ! carbon in leaf litter pool - lignin            (kg_C m-2)

       clitrm(i) = clitrm(i) + (fracrm * adfallr(i))  + (fracrm_p * adfallr(i)) + (INput_P(isoil)) ! carbon in fine root litter pool - metabolic  (kg_C m-2)
       clitrs(i) = clitrs(i) + (fracrs * adfallr(i))  + (fracrs_p * adfallr(i))  ! carbon in fine root litter pool - structural (kg_C m-2)
       clitrl(i) = clitrl(i) + (fracrl * adfallr(i))  + (fracrl_p * adfallr(i))  ! carbon in fine root litter pool - lignin     (kg_C m-2)

       clitwm(i) = clitwm(i) + (fracwm * adfallw(i))  + (fracwm_p * adfallw(i))  ! carbon in woody litter pool - metabolic      (kg_C m-2)
       clitws(i) = clitws(i) + (fracws * adfallw(i))  + (fracws_p * adfallw(i))  ! carbon in woody litter pool - structural     (kg_C m-2)
       clitwl(i) = clitwl(i) + (fracwl * adfallw(i))  + (fracwl_p * adfallw(i))  ! carbon in woody litter pool - lignin	     (kg_C m-2)

       !
       ! ---------------------------------------------------------------
       ! calculate microbial growth rates based on available C sources
       ! to microbes (substrate : litter, C in slow, passive pools)
       ! the amount of biomass added cannot be larger than the amount of
       ! available carbon from substrates and other pools at this point.
       ! ---------------------------------------------------------------
       !
       outcrs(i) = min(decomps(i) * krs * clitrs(i),clitrs(i))
       outcws(i) = min(decompl(i) * kws * clitws(i),clitws(i))
       outcls(i) = min(decompl(i) * kls * clitls(i),clitls(i))
       outclm(i) = min(decompl(i) * klm * clitlm(i),clitlm(i))
       outcrm(i) = min(decomps(i) * krm * clitrm(i),clitrm(i))
       outcwm(i) = min(decompl(i) * kwm * clitwm(i),clitwm(i))
       outcnb(i) = min(decomps(i) * knb * csoislon(i),csoislon(i))
       !
       outcpb(i) = min(decomps(i) * kpb * csoislop(i),csoislop(i))
       !
       outcsb(i) = min(decomps(i) * ksb * csoipas(i),csoipas(i))
       !
       ! ---------------------------------------------------------------
       ! calculate turnover of microbial biomass
       ! two disctinct pools: one with rapid turnover, and one with slow
       ! turnover rate
       ! ---------------------------------------------------------------
       !
       outcbp(i) = min(kbp * cmicp(i),cmicp(i))
       outcbn(i) = min(kbn * cmicn(i),cmicn(i))
       !
       ! ---------------------------------------------------------------------
       ! recycle microbes back to respective microbial pools based on effac as
       ! discussed in NCSOIL model from Molina et al., 1983
       ! ---------------------------------------------------------------------
       !
       outcbp(i) = outcbp(i) *  effac
       outcbn(i) = outcbn(i) *  effac
       !
       ! -------------------------------------------------------------------------
       ! have to adjust inputs into microbial pool for the slow
       ! and passive carbon amounts that are leaving their respective
       ! pools at an increased rate during the spinup procedure.
       ! these values should be decreased by the respective spinup factors
       ! because the microbial pools will otherwise become larger without
       ! scientific reason due to the spinup relationships used.
       ! 3 main pools: outcpb, outcnb, outcsb
       ! -------------------------------------------------------------------------
       !
       dbdt(i) =  outcrs(i) * yrs + outcws(i) * yws +&
            outcls(i) * yls + outclm(i) * ylm +&
            outcrm(i) * yrm + outcwm(i) * ywm +&
            outcnb(i) * ynb + &
            outcpb(i) * ypb +&
            outcsb(i) * ysb - outcbp(i) -&
            outcbn(i)
       ! 
       ! -------------------------------------------------------------------------
       ! change in non-protected organic matter from growth in microbial
       ! biomass, lignin input, and stablized organic matter pool
       ! the flow out of the pool from its decomposition is always less
       ! the yield--which is factored into the pool it is flowing into
       ! -------------------------------------------------------------------------
       !
       outcll(i) = min(decompl(i) * kll * clitll(i),clitll(i))
       outcrl(i) = min(decomps(i) * krl * clitrl(i),clitrl(i))
       outcwl(i) = min(decompl(i) * kwl * clitwl(i),clitwl(i))
       outcns(i) = min(decomps(i) * kns * csoislon(i),         &
            csoislon(i))
       !
       ! ------------------------------------------------------------ 
       ! the lig_frac  factor only applies to lignin content...half goes to
       ! protected slow OM, and half goes to non protected slow OM
       ! ------------------------------------------------------------
       !
       dcndt(i) =  (lig_frac * (outcll(i) * yll + outcrl(i) * yrl + &
            outcwl(i) * ywl) +                             &
            (1.0_r8 - fbpom) * (ybn * outcbn(i) +              &
            ybp * outcbp(i))) - outcnb(i) - outcns(i)

       !
       ! ------------------------------------------------------------
       ! change in protected organic matter from growth in microbial 
       ! biomass, lignin input, and stablized organic matter pool
       ! ------------------------------------------------------------
       !
       outcps(i) = min(decomps(i) * kps * csoislop(i), &
            csoislop(i))
       !
       ! ------------------------------------------------------------
       ! the lig_frac factor only applies to lignin content...half goes to
       ! protected slow OM, and half goes to non protected slow OM
       ! ------------------------------------------------------------
       !
       dcpdt(i) = (lig_frac * (outcll(i)*yll+outcrl(i)*yrl +  &
            outcwl(i) * ywl) +                          & 
            fbpom * (ybn * outcbn(i) +                  &
            ybp * outcbp(i))) - outcpb(i) - outcps(i)
       !
       ! ----------------------------------------------------------------------
       ! change in stablized organic matter (passive pool) from growth
       ! in microbial biomass, and changes in protected and unprotected
       ! SOM
       !
       ! add a loss of C due to leaching out of the profile, based
       ! on approximation of CENTURY model below 1 m in depth
       ! based on water in the profile, and texture of soil
       ! tuned to known outputs or leaching that has been measured in the field
       ! at Arlington-WI (Courtesy K. Brye, MS) and applied to the global scale
       ! on average, this calibration yields about 10-50 Kg C ha-1 yr-1 leaching
       ! depending on C in soil...will need to be tied to an amount of water
       ! flowing through the profile based upon precipitation eventually
       ! ----------------------------------------------------------------------
       !
       h20    = 0.30e-03_r8
       !
       ! h20 is a constant relating to the flow of water through the top 1 m of the
       ! profile 
       ! use texfact -- the % sand -- or texture factor effect on leaching (see Parton
       ! et al. (1991) calculated from the average sand content of top 1 m of soil
       ! in the model
       !
       fleach = h20/18.0_r8 * (0.01_r8 + 0.04_r8 * texfact)
       !
       ! --------------------------------------------------------------------
       ! change in passive organic carbon pool
       ! ---------------------------------------------------------------------
       !
       dcsdt(i) = ((yns * outcns(i)) + (yps * outcps(i))) -   &
                  outcsb(i) -  (fleach * csoipas(i))
       !
       cleach(i) = fleach * csoipas(i) + fleach * csoislop(i) +   &
                   fleach * csoislon(i)
       !
       ynleach(i) = ynleach(i) + fleach * csoipas(i)/cnr(2) +   &
                                 fleach * csoislop(i)/cnr(3) +   &
                                 fleach * csoislon(i)/cnr(4)

       ynleach_p(i) = ynleach_p(i) + fleach * csoipas (i)/cpr(2) +   & !c:p ratio of passive soil carbon
                                     fleach * csoislop(i)/cpr(3) +   & !c:p ratio of protected slow soil carbon
                                     fleach * csoislon(i)/cpr(4)       !c:p ratio of non-protected slow soil C

       !
       ! update slow pools of carbon for leaching losses
       !
       !       dcndt(i) = dcndt(i) - fleach * csoislon(i) 
       !       dcpdt(i) = dcpdt(i) - fleach * csoislop(i) 

       dcndt(i) = dcndt(i) - fleach * csoislon(i) - (  fleach * csoislon(i)/cnr(4) ) -(fleach * csoislon(i)/cpr(4))
       dcpdt(i) = dcpdt(i) - fleach * csoislop(i) - (  fleach * csoislop(i)/cnr(3) ) -(fleach * csoislop(i)/cpr(3))
       !
       IF (spin .eq. spinmax) THEN

          yrleach(i) =  cleach(i) + yrleach(i)
       !
       END IF
       !
       ! ---------------------------------------------------------------------
       ! calculate the amount of net N mineralization or immobilization
       ! ---------------------------------------------------------------------
       !
       ! uptake of n by growth of microbial biomass
       !
       ! immobilized n used for requirements of microbial growth
       ! is based on flow of carbon and the difference of C/N ratio of
       ! the microbes and their efficiency versus the C/N ratio of the
       ! material that is being decomposed 
       !
       ! ------------------------------
       ! structural root decomposition  nitrogen
       ! ------------------------------
       !
       IF (yrs/cnr(1) .gt. 1.0_r8/cnr(6)) THEN
          nbiors(i) = (1.0_r8/cnr(6) - yrs/cnr(1))   &
               * outcrs(i)
          nminrs(i) = 0.0_r8
          !
       ELSE
          nminrs(i) = (1.0_r8/cnr(6) - yrs/cnr(1))   &
               * outcrs(i)
          nbiors(i) = 0.0_r8
       END IF
       !
       ! ------------------------------
       ! structural root decomposition phosphorus
       ! ------------------------------
       !
       IF (yrs/cpr(1) .gt. 1.0_r8/cpr(6)) THEN
         nbiors_p(i) = (1.0_r8/cpr(6) - yrs/cpr(1))   &
                      * outcrs(i)
         nminrs_p(i) = 0.0_r8
       !
       ELSE
         nminrs_p(i) = (1.0_r8/cpr(6) - yrs/cpr(1))   &
                      * outcrs(i)
         nbiors_p(i) = 0.0_r8
       END IF

       !
       ! ------------------------------
       ! structural leaf decomposition  nitrogen
       ! ------------------------------
       !
       IF (yls/cnr(1) .gt. 1.0_r8/cnr(6)) THEN
          nbiols(i) = (1.0_r8/cnr(6) - yls/cnr(1))  &
               * outcls(i)
          nminls(i) = 0.0_r8
          !
       ELSE
          nminls(i) = (1.0_r8/cnr(6) - yls/cnr(1))  &
               * outcls(i)
          nbiols(i) = 0.0_r8
       END IF
       !
       ! ------------------------------
       ! structural leaf decomposition  phosphorus
       ! ------------------------------
       !
       IF (yls/cpr(1) .gt. 1.0_r8/cpr(6)) THEN
         nbiols_p(i) = (1.0_r8/cpr(6) - yls/cpr(1))  &
                      * outcls(i)
         nminls_p(i) = 0.0_r8
       !
       ELSE
         nminls_p(i) = (1.0_r8/cpr(6) - yls/cpr(1))  &
                      * outcls(i)
         nbiols_p(i) = 0.0_r8
       END IF

       !
       ! ------------------------------
       ! structural wood decomposition nitrogen
       ! ------------------------------
       !
       IF (yws/cnr(1) .gt. 1.0_r8/cnr(8)) THEN
          nbiows(i) = (1.0_r8/cnr(8) - yws/cnr(1))  &
               * outcws(i)
          nminws(i) = 0.0_r8
          !
       ELSE
          nminws(i) = (1.0_r8/cnr(8) - yws/cnr(1))  &
               * outcws(i)
          nbiows(i) = 0.0_r8
       END IF
       !
       ! ------------------------------
       ! structural wood decomposition   phosphorus
       ! ------------------------------
       ! 
       IF (yws/cpr(1) .gt. 1.0_r8/cpr(8)) THEN
         nbiows_p(i) = (1.0_r8/cpr(8) - yws/cpr(1))  &
                      * outcws(i)
         nminws_p(i) = 0.0_r8
       !
       ELSE
         nminws_p(i) = (1.0_r8/cpr(8) - yws/cpr(1))  &
                      * outcws(i)
         nbiows_p(i) = 0.0_r8
       END IF

       !
       ! ------------------------------
       ! metabolic wood decomposition    nitrogen
       ! ------------------------------
       !
       IF (ywm/cnr(1) .gt. 1.0_r8/cnr(8)) THEN
          nbiowm(i) = (1.0_r8/cnr(8) - ywm/cnr(1))  &
               * outcwm(i)
          nminwm(i) = 0.0_r8
          !
       ELSE
          nminwm(i) = (1.0_r8/cnr(8) - ywm/cnr(1))  &
               * outcwm(i)
          nbiowm(i) = 0.0_r8
       END IF
       !
       ! ------------------------------
       ! metabolic wood decomposition    phosphorus
       ! ------------------------------
       !
       IF (ywm/cpr(1) .gt. 1.0_r8/cpr(8)) THEN
         nbiowm_p(i) = (1.0_r8/cpr(8) - ywm/cpr(1))  &
                      * outcwm(i)
         nminwm_p(i) = 0.0_r8
       !
       ELSE
         nminwm_p(i) = (1.0_r8/cpr(8) - ywm/cpr(1))  &
                      * outcwm(i)
         nbiowm_p(i) = 0.0_r8
       END IF

       !
       ! ------------------------------
       ! metabolic leaf decomposition    nitrogen
       ! ------------------------------
       !
       IF (ylm/cnr(1) .gt. 1.0_r8/cnr(7)) THEN
          nbiolm(i) = (1.0_r8/cnr(7) - ylm/cnr(1))   &
               * outclm(i)
          nminlm(i) = 0.0_r8
          !
       ELSE
          nminlm(i) = (1.0_r8/cnr(7) - ylm/cnr(1))  &
               * outclm(i)
          nbiolm(i) = 0.0_r8
       END IF
       !
       ! ------------------------------
       ! metabolic leaf decomposition   phosphorus
       ! ------------------------------
       !
       IF (ylm/cpr(1) .gt. 1.0_r8/cpr(7)) THEN
         nbiolm_p(i) = (1.0_r8/cpr(7) - ylm/cpr(1))   &
                      * outclm(i)
         nminlm_p(i) = 0.0_r8
       !
       ELSE
         nminlm_p(i) = (1.0_r8/cpr(7) - ylm/cpr(1))  &
                      * outclm(i)
         nbiolm_p(i) = 0.0_r8
       END IF

       !
       ! ------------------------------
       ! metabolic root decomposition   nitrogen
       ! ------------------------------
       !
       IF (yrm/cnr(1) .gt. 1.0_r8/cnr(7)) THEN
          nbiorm(i) = (1.0_r8/cnr(7) - yrm/cnr(1))   &
               * outcrm(i)
          nminrm(i) = 0.0_r8
          !
       ELSE
          nminrm(i) = (1.0_r8/cnr(7) - yrm/cnr(1))  &
               * outcrm(i)
          nbiorm(i) = 0.0_r8
       END IF
       !
       ! ------------------------------
       ! metabolic root decomposition    phosphorus
       ! ------------------------------
       !
       IF (yrm/cpr(1) .gt. 1.0_r8/cpr(7)) THEN
         nbiorm_p(i) = (1.0_r8/cpr(7) - yrm/cpr(1))   &
                      * outcrm(i)
         nminrm_p(i) = 0.0_r8
       !
       ELSE
         nminrm_p(i) = (1.0_r8/cpr(7) - yrm/cpr(1))  &
                      * outcrm(i)
         nbiorm_p(i) = 0.0_r8
       END IF

       !
       ! ----------------------------------------------
       ! non-protected organic matter decomposition    nitrogen
       ! ----------------------------------------------
       !
       IF (ynb/cnr(1) .gt. 1.0_r8/cnr(4)) THEN
          nbioslon(i) = (1.0_r8/cnr(4) - ynb/cnr(1))  &
               * outcnb(i)
          nminslon(i) = 0.0_r8
          !
       ELSE
          nminslon(i) = (1.0_r8/cnr(4) - ynb/cnr(1))   &
               * outcnb(i)
          nbioslon(i) = 0.0_r8
       END IF


       !
       ! ----------------------------------------------
       ! non-protected organic matter decomposition     phosphorus
       ! ----------------------------------------------
       !
       IF (ynb/cpr(1) .gt. 1.0_r8/cpr(4)) THEN
         nbioslon_p(i) = (1.0_r8/cpr(4) - ynb/cpr(1))  &
                        * outcnb(i)
         nminslon_p(i) = 0.0_r8
       !
       ELSE
         nminslon_p(i) = (1.0_r8/cpr(4) - ynb/cpr(1))   &
                        * outcnb(i)
         nbioslon_p(i) = 0.0_r8
       END IF

       !
       ! ----------------------------------------------
       ! protected organic matter decomposition   nitrogen
       ! ----------------------------------------------
       !
       IF (ypb/cnr(1) .gt. 1.0_r8/cnr(3)) THEN
          nbioslop(i) = (1.0_r8/cnr(3) - ypb/cnr(1)) &
               * outcpb(i)
          nminslop(i) = 0.0_r8
          !
       ELSE
          nminslop(i) = (1.0_r8/cnr(3) - ypb/cnr(1))  &
               * outcpb(i)
          nbioslop(i) = 0.0_r8
       END IF
       !
       ! ----------------------------------------------
       ! protected organic matter decomposition     phosphorus
       ! ----------------------------------------------
       !
       IF (ypb/cpr(1) .gt. 1.0_r8/cpr(3)) THEN
         nbioslop_p(i) = (1.0_r8/cpr(3) - ypb/cpr(1)) &
                        * outcpb(i)
         nminslop_p(i) = 0.0_r8
       !
       ELSE
         nminslop_p(i) = (1.0_r8/cpr(3) - ypb/cpr(1))  &
                        * outcpb(i)
         nbioslop_p(i) = 0.0_r8
       END IF

       !
       ! ----------------------------------------------
       ! stablized organic matter decomposition   nitrogen
       ! ----------------------------------------------
       !
       IF (ysb/cnr(1) .gt. 1.0_r8/cnr(2)) THEN
          nbiopas(i) = (1.0_r8/cnr(2) - ysb/cnr(1)) &
               * outcsb(i)
          nminpas(i) = 0.0_r8
          !
       ELSE
          nminpas(i) = (1.0_r8/cnr(2) - ysb/cnr(1))      &
               * outcsb(i)
          nbiopas(i) = 0.0_r8
       END IF
       !
       ! ----------------------------------------------
       ! stablized organic matter decomposition    phosphorus
       ! ----------------------------------------------
       !
       IF (ysb/cpr(1) .gt. 1.0_r8/cpr(2)) THEN
         nbiopas_p(i) = (1.0_r8/cpr(2) - ysb/cpr(1)) &
                       * outcsb(i)
         nminpas_p(i) = 0.0_r8
       !
       ELSE
         nminpas_p(i) = (1.0_r8/cpr(2) - ysb/cpr(1))      &
                       * outcsb(i)
         nbiopas_p(i) = 0.0_r8
       END IF

       !
       ! ----------------------------------------------
       ! total immobilized N used for biomass growth
       ! ----------------------------------------------
       !
       totimm(i) = nbiors(i) + nbiols(i) + nbiows(i) + nbiowm(i)      &
                 + nbiolm(i) + nbiorm(i) + nbioslon(i) + nbioslop(i)  &
                 + nbiopas(i)
       !
       ! ----------------------------------------------
       ! total immobilized P used for biomass growth
       ! ----------------------------------------------
       !
       totimm_p(i) = nbiors_p(i) + nbiols_p(i) + nbiows_p(i) + nbiowm_p(i)      &
                   + nbiolm_p(i) + nbiorm_p(i) + nbioslon_p(i) + nbioslop_p(i)  &
                   + nbiopas_p(i)
       !
       !
       ! -----------------------------------------------------------------------------
       ! gross amount of N mineralized by decomposition of C by microbial biomass
       ! assume that N is attached to the flow of C by the C/N ratio of the substrate
       ! also assume that the amount of N attached to CO2 that is respired is also
       ! mineralized (i.e. the amount of N mineralized is related to the total outflow
       ! of carbon, and not the efficiency or yield)..see Parton et al., 1987
       ! -----------------------------------------------------------------------------
       !
       totmin(i) = nminrs(i) + nminls(i) + nminws(i) + nminwm(i)          &
            + nminlm(i) + nminrm(i) + nminslon(i) + nminslop(i)      &
            + nminpas(i)
       !
       !
       ! -----------------------------------------------------------------------------
       ! gross amount of P mineralized by decomposition of C by microbial biomass
       ! assume that P is attached to the flow of C by the C/P ratio of the substrate
       ! also assume that the amount of P attached to CO2 that is respired is also
       ! mineralized (i.e. the amount of P mineralized is related to the total outflow
       ! of carbon, and not the efficiency or yield)..see Parton et al., 1987
       ! -----------------------------------------------------------------------------
       !
       totmin_p(i) = nminrs_p(i) + nminls_p(i) + nminws_p(i) + nminwm_p(i)          &
                   + nminlm_p(i) + nminrm_p(i) + nminslon_p(i) + nminslop_p(i)      &
                   + nminpas_p(i)

       ! -----------------------------------------------------------------------------
       ! when carbon is transferred from one pool to another, each pool has a distinct
       ! C:N ratio.  In the case of pools where carbon is moving from the pool to 
       ! the microbial biomass (used for growth/assimilation), net mineralization
       ! takes place (N is released) after the requirements of building the biomass
       ! are met.  In the cases of other transformations of C, N is not conserved
       ! if it follows from one pool to another which has a different C:N ratio;
       ! either N is released or is needed to make the transformation and keep N
       ! conserved in the model. 
       !
       ! other calculations of either N release or immobilization to keep track of
       ! the budget
       !
       nrelps(i) = outcps(i) * (1.0_r8/cnr(3) - 1.0_r8/cnr(2))
       nrelns(i) = outcns(i) * (1.0_r8/cnr(4) - 1.0_r8/cnr(2))
       nrelbn(i) = (1.0_r8-fbpom) * outcbn(i) * (1.0_r8/cnr(1) - 1.0_r8/cnr(4)) +   &
            (1.0_r8-fbpom) * outcbp(i) * (1.0_r8/cnr(1) - 1.0_r8/cnr(4))
       nrelbp(i) = fbpom * outcbp(i) * (1.0_r8/cnr(1) - 1.0_r8/cnr(3)) +     &
            fbpom * outcbn(i) * (1.0_r8/cnr(1) - 1.0_r8/cnr(3))
       nrelll(i) = lig_frac * outcll(i) * (1.0_r8/cnr(5) - 1.0_r8/cnr(3)) +  &
            lig_frac * outcll(i) * (1.0_r8/cnr(5) - 1.0_r8/cnr(4))
       nrelrl(i) = lig_frac * outcrl(i) * (1.0_r8/cnr(5) - 1.0_r8/cnr(3)) +  &
            lig_frac * outcrl(i) * (1.0_r8/cnr(5) - 1.0_r8/cnr(4))
       nrelwl(i) = lig_frac * outcwl(i) * (1.0_r8/cnr(5) - 1.0_r8/cnr(3)) +  &
            lig_frac * outcwl(i) * (1.0_r8/cnr(5) - 1.0_r8/cnr(4))
       !
       totnrel(i) = nrelps(i) + nrelns(i) + nrelbn(i) +      &
            nrelbp(i) + nrelll(i) + nrelrl(i) + nrelwl(i)
       !
       ! -----------------------------------------------------------------------------
       ! when carbon is transferred from one pool to another, each pool has a distinct
       ! C:P ratio.  In the case of pools where carbon is moving from the pool to 
       ! the microbial biomass (used for growth/assimilation), net mineralization
       ! takes place (P is released) after the requirements of building the biomass
       ! are met.  In the cases of other transformations of C, P is not conserved
       ! if it follows from one pool to another which has a different C:P ratio;
       ! either P is released or is needed to make the transformation and keep N
       ! conserved in the model. 
       !
       ! other calculations of either P release or immobilization to keep track of
       ! the budget
       !
        nrelps_p(i) = outcps(i) * (1.0_r8/cpr(3) - 1.0_r8/cpr(2))
        nrelns_p(i) = outcns(i) * (1.0_r8/cpr(4) - 1.0_r8/cpr(2))

        nrelbn_p(i) = (1.0_r8-fbpom) * outcbn(i) * (1.0_r8/cpr(1) - 1.0_r8/cpr(4)) +   &
                      (1.0_r8-fbpom) * outcbp(i) * (1.0_r8/cpr(1) - 1.0_r8/cpr(4))

        nrelbp_p(i) = fbpom * outcbp(i) * (1.0_r8/cpr(1) - 1.0_r8/cpr(3)) +     &
                      fbpom * outcbn(i) * (1.0_r8/cpr(1) - 1.0_r8/cpr(3))

        nrelll_p(i) = lig_frac * outcll(i) * (1.0_r8/cpr(5) - 1.0_r8/cpr(3)) +  &
                      lig_frac * outcll(i) * (1.0_r8/cpr(5) - 1.0_r8/cpr(4))
        nrelrl_p(i) = lig_frac * outcrl(i) * (1.0_r8/cpr(5) - 1.0_r8/cpr(3)) +  &
                      lig_frac * outcrl(i) * (1.0_r8/cpr(5) - 1.0_r8/cpr(4))
        nrelwl_p(i) = lig_frac * outcwl(i) * (1.0_r8/cpr(5) - 1.0_r8/cpr(3)) +  &
                      lig_frac * outcwl(i) * (1.0_r8/cpr(5) - 1.0_r8/cpr(4))
       !
        totnrel_p(i) = nrelps_p(i) + nrelns_p(i) + nrelbn_p(i) +      &
                       nrelbp_p(i) + nrelll_p(i) + nrelrl_p(i) + nrelwl_p(i)
       !
       !
       ! -----------------------------------------------------------------------------
       ! calculate whether net mineralization or immobilization occurs
       ! on a grid cell basis -- tnmin is an instantaneous value for each time step
       ! it is passed along to stats to calculate, daily, monthly and annual totals
       ! of nitrogen mineralization
       ! this is for mineralization/immobilization that is directly related to 
       ! microbial processes (oxidation of carbon)
       !
       ! the value of totnrel(i) would need to be added to complete the budget
       ! of N in the model. Because it can add/subtract a certain amount of N
       ! from the amount of net mineralization.  However, these transformations
       ! are not directly related to microbial decomposition, so do we add them
       ! into the value or not?
       ! -----------------------------------------------------------------------------
       !
       netmin(i) = totmin(i) + totimm(i) + totnrel(i) 
       IF (netmin(i) .gt. 0.00_r8) THEN
          !
          tnmin(i) = netmin(i)
          !
       ELSE
          !  
          tnmin(i) = 0.00_r8 
          !
       END IF

       ! -----------------------------------------------------------------------------
       ! calculate whether net mineralization or immobilization occurs
       ! on a grid cell basis -- tnmin_p is an instantaneous value for each time step
       ! it is passed along to stats to calculate, daily, monthly and annual totals
       ! of nitrogen mineralization
       ! this is for mineralization/immobilization that is directly related to 
       ! microbial processes (oxidation of carbon)
       !
       ! the value of totnrel(i) would need to be added to complete the budget
       ! of N in the model. Because it can add/subtract a certain amount of N
       ! from the amount of net mineralization.  However, these transformations
       ! are not directly related to microbial decomposition, so do we add them
       ! into the value or not?
       ! -----------------------------------------------------------------------------
       !
       netmin_p(i) = totmin_p(i) + totimm_p(i) + totnrel_p(i) 
       IF (netmin_p(i) .gt. 0.00_r8) THEN
       !
           tnmin_p(i) = netmin_p(i)
       !
       ELSE
       !
           tnmin_p(i) = 0.00_r8 
       !
       END IF
       !
       !
       !
       ! convert value of tnmin of Kg-N/m2/dtime to mole-N/s
       ! based on N = .014 Kg/mole -- divide by the number of seconds in daily timestep
       !
       tnmin(i) = tnmin(i)/(86400.0_r8 * 0.0140_r8)
       !
       !
       ! convert value of tnmin of Kg-P/m2/dtime to mole-P/s
       ! based on P = .0309 Kg/mole -- divide by the number of seconds in daily timestep
       !
       tnmin_p(i) = tnmin_p(i)/(86400.0_r8 * 0.0309_r8)

       !
       ! ---------------------------------------------------
       ! update soil c pools for transformations of c and n
       ! ---------------------------------------------------
       !
       totcmic(i)  = max(totcmic(i)  + dbdt(i), 0.00_r8)
       csoislon(i) = max(csoislon(i) + dcndt(i),0.00_r8)
       csoislop(i) = max(csoislop(i) + dcpdt(i),0.00_r8)
       csoipas(i)  = max(csoipas(i)  + dcsdt(i),0.00_r8)
       clitlm(i)   = max(clitlm(i)  - outclm(i),0.00_r8)
       clitls(i)   = max(clitls(i)  - outcls(i),0.00_r8)
       clitll(i)   = max(clitll(i)  - outcll(i),0.00_r8)
       clitrm(i)   = max(clitrm(i)  - outcrm(i),0.00_r8)
       clitrs(i)   = max(clitrs(i)  - outcrs(i),0.00_r8)
       clitrl(i)   = max(clitrl(i)  - outcrl(i),0.00_r8)
       clitwm(i)   = max(clitwm(i)  - outcwm(i),0.00_r8)
       clitws(i)   = max(clitws(i)  - outcws(i),0.00_r8)
       clitwl(i)   = max(clitwl(i)  - outcwl(i),0.00_r8)
       !
       ! -----------------------------------------------------------
       ! update soil n pools based on c:n ratios of each pool
       ! this approach is assuming that the c:n ratios are remaining
       ! constant through the simulation. flow of nitrogen is attached
       ! to carbon 
       ! -----------------------------------------------------------
       !
       totnmic(i)  = totcmic(i) /cnr(1)
       nsoislon(i) = csoislon(i)/cnr(4)
       nsoislop(i) = csoislop(i)/cnr(3)
       nsoipas(i)  = csoipas(i) /cnr(2)
       nlitlm(i)   = clitlm(i)  /cnr(7)
       nlitls(i)   = clitls(i)  /cnr(6)
       nlitll(i)   = clitll(i)  /cnr(5)
       nlitrm(i)   = clitrm(i)  /cnr(7)
       nlitrs(i)   = clitrs(i)  /cnr(6)
       nlitrl(i)   = clitrl(i)  /cnr(5)
       nlitwm(i)   = clitwm(i)  /cnr(8)
       nlitws(i)   = clitws(i)  /cnr(8)
       nlitwl(i)   = clitwl(i)  /cnr(8)
       !
       ! -----------------------------------------------------------
       ! update soil P pools based on c:p ratios of each pool
       ! this approach is assuming that the c:p ratios are remaining
       ! constant through the simulation. flow of phosphorus is attached
       ! to carbon 
       ! -----------------------------------------------------------
       !
       totnmic_p(i)  = totcmic(i) /cpr(1)
       nsoislon_p(i) = csoislon(i)/cpr(4)
       nsoislop_p(i) = csoislop(i)/cpr(3)
       nsoipas_p(i)  = csoipas(i) /cpr(2)
       nlitlm_p(i)   = clitlm(i)  /cpr(7)
       nlitls_p(i)   = clitls(i)  /cpr(6)
       nlitll_p(i)   = clitll(i)  /cpr(5)
       nlitrm_p(i)   = clitrm(i)  /cpr(7)
       nlitrs_p(i)   = clitrs(i)  /cpr(6)
       nlitrl_p(i)   = clitrl(i)  /cpr(5)
       nlitwm_p(i)   = clitwm(i)  /cpr(8)
       nlitws_p(i)   = clitws(i)  /cpr(8)
       nlitwl_p(i)   = clitwl(i)  /cpr(8)

       !
       ! total above and belowground litter
       !
       totlit(i) =  clitlm(i) + clitls(i) + clitll(i) +  &
                    clitrm(i) + clitrs(i) + clitrl(i) +  &
                    clitwm(i) + clitws(i) + clitwl(i)
       !
       ! sum total aboveground litter (leaves and wood)
       !
       totalit(i) = clitlm(i) + clitls(i) + clitwm(i) +  &
                    clitll(i) + clitws(i) + clitwl(i)
       !
       ! sum total belowground litter (roots) 
       !
       totrlit(i) = clitrm(i) + clitrs(i) + clitrl(i)
       !
       ! determine total soil carbon amounts (densities are to 1 m depth; Kg/m-2)
       !
       totcsoi(i) = csoipas(i) + csoislop(i) +  &
                    totcmic(i) + csoislon(i)
       !
       ! calculate total amount of litterfall occurring (total for year)
       !
       totfall(i) = falll(i) + fallr(i) + fallw(i)
       !
       ! nitrogen 
       !
       ! total nitrogen in litter pools (above and belowground)
       !
       totnlit(i) =  nlitlm(i) + nlitls(i) + nlitrm(i) + nlitrs(i) +  &
                     nlitwm(i) + nlitws(i) + nlitll(i) + nlitrl(i) +  &
                     nlitwl(i)
       !
       ! sum total aboveground litter   (leaves and wood)
       !
       totanlit(i) = nlitlm(i) + nlitls(i) + nlitwm(i) +  &
                     nlitll(i) + nlitws(i) + nlitwl(i)
       !
       ! sum total belowground litter  (roots)
       !
       totrnlit(i) = nlitrm(i) + nlitrs(i) + nlitrl(i)
       !
       ! total soil nitrogen to 1 m depth (kg-N/m**2)
       !
       totnsoi(i) = nsoislop(i) + nsoislon(i) +  &
                    nsoipas(i)  + totnmic(i) + totnlit(i)
       !
       ! phosphorus 
       !
       ! total phosphorus in litter pools (above and belowground)
       !
       totnlit_p(i) =  nlitlm_p(i) + nlitls_p(i) + nlitrm_p(i) + nlitrs_p(i) +  &
                       nlitwm_p(i) + nlitws_p(i) + nlitll_p(i) + nlitrl_p(i) +  &
                       nlitwl_p(i)
       !
       ! sum total aboveground litter   (leaves and wood)
       !
       totanlit_p(i) = nlitlm_p(i) + nlitls_p(i) + nlitwm_p(i) +  &
                       nlitll_p(i) + nlitws_p(i) + nlitwl_p(i)
       !
       ! sum total belowground litter  (roots)
       !
       totrnlit_p(i) = nlitrm_p(i) + nlitrs_p(i) + nlitrl_p(i)
       !
       ! total soil phosphorus to 1 m depth (kg-N/m**2)
       !
       totnsoi_p(i) = nsoislop_p(i) + nsoislon_p(i) +  &
                      nsoipas_p(i)  + totnmic_p(i) + totnlit_p(i)

       !
       ! --------------------------------------------------------------------------
       ! calculate running sum of yearly net mineralization, and nitrogen in pool
       ! available to plants for uptake--during spin up period, can only count one
       ! of the cycles for each timestep--otherwise false additions will result
       ! values of yearly mineralization are in Kg/m-2
       ! --------------------------------------------------------------------------
       !
       IF (spin .eq. spinmax) THEN
          !
          storedn(i)  = storedn(i) + tnmin(i)
          storedn_p(i)  = storedn_p(i) + tnmin(i)

          !
       END IF
       !
       ! calculate total amount of carbon in soil at end of cycle
       ! this is used to help calculate the amount of carbon that is respired
       ! by decomposing microbial biomass
       !
       totcend(i) = totlit(i) + totcsoi(i)
       !
       ! --------------------------------------------------------------------------
       ! the amount of co2resp(i) is yearly value and is dependent on the amount
       ! of c input each year, the amount in each pool at beginning of the year,
       ! and the amount left in the pool at the end of the year
       ! along with the amount of root respiration contributing to the flux from
       ! calculations performed in stats.f
       ! --------------------------------------------------------------------------
       !
       IF (spin .eq. spinmax) THEN 
          !
          ! --------------------------------------------------------------------------
          ! only count the last cycle in the spin-up for co2soi
          ! when the iyear is less than the nspinsoil value...otherwise
          ! an amount of CO2 respired will be about 10 times the actual
          ! value because this routine is called articially 10 extra times
          ! each time step to spin up the soil carbon
          !
          ! add n-deposition due to rainfall once each day, and
          ! the amount of N fixed through N-fixers.  These equations
          ! are based on the annual precip input (cm) and are from
          ! the CENTURY model...Parton et al., 1987.
          ! The base equations are in units of (g) N m-2 so have to
          ! divide by 1000 to put in units of Kg.
          !
          ! the values in the equation of 0.21 and -0.18 were adjusted to reflect
          ! average daily inputs when no precipitation was falling - the original
          ! constants are for the entire year 
          ! --------------------------------------------------------------------------
          !          
          IF (iday .eq. 1 .and. imonth .eq. 1) THEN

             deposn(i)  = (0.210_r8  + 0.00280_r8 * (ayprcp(i)*0.10_r8))*1.e-3_r8
             fixsoin(i) = (-0.180_r8 + 0.140_r8  * (ayprcp(i)*0.10_r8))*1.e-3_r8
             storedn(i)   = storedn(i) + deposn(i) + fixsoin(i)
            !storedn_p(i)   = storedn_p(i) + deposn(i) + fixsoin(i)
             storedn_p(i)   = storedn_p(i) + 0.000004_r8 ! kg_P/ m2

          END IF
          !
          ! --------------------------------------------------------------------------
          ! add to the daily total of co2 flux leaving the soil from microbial
          ! respiration -- instantaneous value for each timestep
          ! since this subroutine gets called daily...instantaneous fluxes
          ! the fluxes need to be put on a per second basis, which will be dependent
          ! on the timestep.  Furthermore, because the biogeochem subroutine does
          ! not get called each timestep...an approximation for a timestep average
          ! microbial flux and nmineralization rate will be applied
          ! --------------------------------------------------------------------------
          !
          ! calculate daily co2 flux due to microbial decomposition
          !
          ! instantaneous microbial co2 flux from soil (mol-CO2 / m-2 / second)
          !
          tco2mic(i) = totcbegin(i) + totcin(i) - totcend(i) - cleach(i)
          !
          ! convert co2 flux from kg C/day  (seconds in a daily timestep) to mol-C/s
          ! based on .012 Kg C/mol
          !                      (kg_C m-2)
          tco2mic(i) = tco2mic(i)/(86400.0_r8 * 0.0120_r8) ! (kg-C/ m-2 / day) to (mol-CO2 / m-2 / second)
          !
       END IF
       !
    END DO ! fim loop npoi
    !
    ! for now, we assume that all soils have a 1% organic content -- 
    ! this is just a place holder until we couple the soil carbon
    ! dynamics to the soil physical properties
    !
    forganic = 0.010_r8
    om_frac  =forganic
    om_csol      = 2.5_r8   ! heat capacity of peat soil *10^6 (J/K m3) (Farouki, 1986)

    DO  k = 1, nsoilay
       DO i = 1, npoi
             ! density of soil material (without pores, not bulk) (kg m-3)
             ! from Campbell and Norman, 1998
             !
             rhosoi(i,k) = 2650.0_r8 * (1.0_r8 - forganic) + 1300.0_r8 * forganic 

             msand    = nint(sand(i,k)) 
             mclay    = nint(clay(i,k)) 

             !isoil =   textcls (msand,mclay)

             IF( NINT (vegtype0(i)) == 15 )THEN
                csoi(i,k) =  870.0_r8 * (1.0_r8 - forganic) + 1920.0_r8 * forganic 
!                csoi(i,k) = 152.5 + 7.122*(ta(i))!Yen (1981)
                !csoi(i,k,j) = 2.106e+3_r8  ! specific heat of ice (J deg-1 kg-1) 
             ELSE

!                csoi (i,k)   = (((1._r8-om_frac)*(2.128_r8*msand+2.385_r8*mclay) / (msand+mclay) + om_csol*om_frac)*1.e6_r8 )/rhosoi(i,k) ! J/(m3 K)  
                csoi(i,k) =  870.0_r8 * (1.0_r8 - forganic) + 1920.0_r8 * forganic
             END IF
       END DO
    END DO
  
    !
    ! return to main
    !
    RETURN
  END SUBROUTINE soilbgc
  !


 

!-------------------------------------------------------------------------
      INTEGER FUNCTION textcls (msand,mclay)
!
! adapted for ibis by cjk 01/11/01
!-------------------------------------------------------------------------
! |
! |                         T R I A N G L E
! | Main program that calls WHAT_TEXTURE, a function that classifies soil
! | in the USDA textural triangle using sand and clay %
! +-----------------------------------------------------------------------
! | Created by: aris gerakis, apr. 98 with help from brian baer
! | Modified by: aris gerakis, july 99: now all borderline cases are valid
! | Modified by: aris gerakis, 30 nov 99: moved polygon initialization to
! |              main program
! +-----------------------------------------------------------------------
! | COMMENTS
! | o Supply a data file with two columns, in free format:  1st column sand,
! |   2nd column clay %, no header.  The output is a file with the classes.
! +-----------------------------------------------------------------------
! | You may use, distribute and modify this code provided you maintain
! ! this header and give appropriate credit.
! +-----------------------------------------------------------------------
!
! code adapted for IBIS by cjk 01-11-01
!
!
      INTEGER :: msand
      INTEGER :: mclay
!
!      LOGICAL :: inpoly
!
      REAL(KIND=r8)    :: silty_loam     (1:7,1:2)
      REAL(KIND=r8)    :: sandy          (1:7,1:2)
      REAL(KIND=r8)    :: silty_clay_loam(1:7,1:2) 
      REAL(KIND=r8)    :: loam           (1:7,1:2)
      REAL(KIND=r8)    :: clay_loam      (1:7,1:2)
      REAL(KIND=r8)    :: sandy_loam     (1:7,1:2)
      REAL(KIND=r8)    :: silty_clay     (1:7,1:2)
      REAL(KIND=r8)    :: sandy_clay_loam(1:7,1:2) 
      REAL(KIND=r8)    :: loamy_sand     (1:7,1:2)
      REAL(KIND=r8)    :: clayey         (1:7,1:2)
!     REAL(KIND=r8)    :: silt           (1:7,1:2) 
      REAL(KIND=r8)    :: sandy_clay     (1:7,1:2)
!
! initalize polygon coordinates:
! each textural class reads in the sand coordinates (1,7) first, and
! then the corresponding clay coordinates (1,7)

!     data silty_loam/0, 0, 23, 50, 20, 8, 0, 12, 27, 27, 0, 0, 12, 0/
!
! because we do not have a separate silt category, have to redefine the
! polygon boundaries for the silt loam  
!
      DATA sandy           /85.0_r8, 90.0_r8, 100.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8,  &
                            10.0_r8,  0.0_r8,   0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8/
      DATA loamy_sand      /70.0_r8, 85.0_r8,  90.0_r8,85.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8,  & 
                            15.0_r8, 10.0_r8,   0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8/
      DATA sandy_loam      /50.0_r8, 43.0_r8,  52.0_r8,52.0_r8,80.0_r8,85.0_r8, 70.0_r8,  &
                             0.0_r8,  7.0_r8,   7.0_r8,20.0_r8,20.0_r8,15.0_r8,  0.0_r8/
      DATA loam            /43.0_r8, 23.0_r8,  45.0_r8,52.0_r8,52.0_r8, 0.0_r8,  0.0_r8,    &
                             7.0_r8, 27.0_r8,  27.0_r8,20.0_r8, 7.0_r8, 0.0_r8,  0.0_r8/
      DATA silty_loam      / 0.0_r8,  0.0_r8,  23.0_r8,50.0_r8, 0.0_r8, 0.0_r8,  0.0_r8, 0.0_r8,    &
                            27.0_r8, 27.0_r8,   0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8/ 
!     DATA silt            /0, 0, 8, 20, 0, 0, 0, 0, 12, 12, 0, 0, 0, 0/
      DATA sandy_clay_loam /52.0_r8, 45.0_r8, 45.0_r8, 65.0_r8, 80.0_r8, 0.0_r8, 0.0_r8,    & 
                            20.0_r8, 27.0_r8, 35.0_r8, 35.0_r8, 20.0_r8, 0.0_r8, 0.0_r8/
      DATA clay_loam       /20.0_r8, 20.0_r8, 45.0_r8, 45.0_r8, 0.0_r8, 0.0_r8, 0.0_r8,     &
                            27.0_r8, 40.0_r8, 40.0_r8, 27.0_r8, 0.0_r8, 0.0_r8, 0.0_r8/
      DATA silty_clay_loam /0.0_r8, 0.0_r8, 20.0_r8, 20.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 27.0_r8,   &
                           40.0_r8, 40.0_r8, 27.0_r8, 0.0_r8, 0.0_r8, 0.0_r8/
      DATA sandy_clay      /45.0_r8, 45.0_r8, 65.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8,      &
                            35.0_r8, 55.0_r8, 35.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8/
      DATA silty_clay      /0.0_r8, 0.0_r8, 20.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 40.0_r8,    &
                           60.0_r8, 40.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8/
      DATA clayey          /20.0_r8, 0.0_r8, 0.0_r8, 45.0_r8, 45.0_r8, 0.0_r8, 0.0_r8,      &
                            40.0_r8, 60.0_r8, 100.0_r8, 55.0_r8, 40.0_r8, 0.0_r8, 0.0_r8/
!
! polygon coordinates  
!
!     sand
!
!     >  85, 90, 100, 0, 0, 0, 0,       ! sand
!     >  70, 85, 90, 85, 0, 0, 0,       ! loamy sand
!     >  50, 43, 52, 52, 80, 85, 70,    ! sandy loam
!     >  43, 23, 45, 52, 52, 0, 0,      ! loam
!     >   0, 0, 23, 50, 0, 0, 0,        ! silt loam (combined with silt)
!     >  52, 45, 45, 65, 80, 0, 0,      ! sandy clay loam
!     >  20, 20, 45, 45, 0, 0, 0,       ! clay loam
!     >   0, 0, 20, 20, 0, 0, 0,        ! silty clay loam
!     >  45, 45, 65, 0, 0, 0, 0,        ! sandy clay
!     >   0, 0, 20, 0, 0, 0, 0,         ! silty clay 
!     >  20, 0, 0, 45, 45, 0, 0         ! clay
!
!      clay
!
!     > 0, 10, 0, 0, 0, 0, 0,           ! sand
!     > 0, 15, 10, 0, 0, 0, 0,          ! loamy sand
!     > 0, 7, 7, 20, 20, 15, 0,         ! sandy loam 
!     > 7, 27, 27, 20, 7, 0, 0,         ! loam
!     > 0, 27, 27, 0, 0, 0, 0,          ! silt loam (combined with silt)
!     > 20, 27, 35, 35, 20, 0, 0,       ! sandy clay loam
!     > 27, 40, 40, 27, 0, 0, 0,        ! clay loam
!     > 27, 40, 40, 27, 0, 0, 0,        ! silty clay loam
!     > 35, 55, 35, 0, 0, 0, 0,         ! sandy clay
!     > 40, 60, 40, 0, 0, 0, 0,         ! silty clay
!     > 40, 60, 100, 55, 40, 0, 0       ! clay
!
! +-----------------------------------------------------------------------
! | figure out what texture grid cell and layer are part of  
! | classify a soil in the triangle based on sand and clay %
! +-----------------------------------------------------------------------
! | Created by: aris gerakis, apr. 98
! | Modified by: aris gerakis, june 99.  Now check all polygons instead of
! | stopping when a right solution is found.  This to cover all borderline 
! | cases.
! +-----------------------------------------------------------------------
!
! find polygon(s) where the point is.  
!
      textcls = 0 
!
      IF (msand .gt. 0 .and. mclay .gt. 0) THEN
         IF (inpoly(sandy, 3, msand, mclay)) THEN
            textcls = 1      ! sand
         END IF
         IF (inpoly(loamy_sand, 4, msand, mclay)) THEN
            textcls = 2      ! loamy sand
         END IF
         IF (inpoly(sandy_loam, 7, msand, mclay)) THEN
            textcls = 3      ! sandy loam
         END IF
         IF (inpoly(loam, 5, msand, mclay)) THEN
            textcls = 4      ! loam
         END IF
         IF (inpoly(silty_loam, 4, msand, mclay)) THEN
            textcls = 5      ! silt loam
         END IF
         IF (inpoly(sandy_clay_loam, 5, msand, mclay)) THEN
            textcls = 6      ! sandy clay loam
         END IF
         IF (inpoly(clay_loam, 4, msand, mclay)) THEN
            textcls = 7      ! clay loam
         END IF
         IF (inpoly(silty_clay_loam, 4, msand, mclay)) THEN
            textcls = 8      ! silty clay loam
         END IF
         IF (inpoly(sandy_clay, 3, msand, mclay)) THEN
            textcls = 9      ! sandy clay
         END IF
         IF (inpoly(silty_clay, 3, msand, mclay)) THEN
            textcls = 10     ! silty clay
         END IF
         IF (inpoly(clayey, 5, msand, mclay)) THEN
            textcls = 11     ! clay
         END IF
      END IF
!
      IF (textcls .eq. 0) THEN
         textcls = 5         ! silt loam
!
!        write (*, 1000) msand, mclay
! 1000   format (/, 1x, 'Texture not found for ', f5.1, ' sand and ', f5.1, ' clay')
      END IF
!
      RETURN
      END FUNCTION textcls
!
!---------------------------------------------------------------------------
      LOGICAL FUNCTION inpoly (poly, npoints, xt, yt)
!
! adapted for ibis by cjk 01/11/01
!---------------------------------------------------------------------------
!
!                            INPOLY
!   Function to tell if a point is inside a polygon or not.
!--------------------------------------------------------------------------
!   Copyright (c) 1995-1996 Galacticomm, Inc.  Freeware source code.
!
!   Please feel free to use this source code for any purpose, commercial
!   or otherwise, as long as you don't restrict anyone else's use of
!   this source code.  Please give credit where credit is due.
!
!   Point-in-polygon algorithm, created especially for World-Wide Web
!   servers to process image maps with mouse-clickable regions.
!
!   Home for this file:  http://www.gcomm.com/develop/inpoly.c
!
!                                       6/19/95 - Bob Stein & Craig Yap
!                                       stein@gcomm.com
!                                       craig@cse.fau.edu
!--------------------------------------------------------------------------
!   Modified by:
!   Aris Gerakis, apr. 1998: 1.  translated to Fortran
!                            2.  made it work with REAL(KIND=r8) coordinates
!                            3.  now resolves the case where point falls
!                                on polygon border.
!   Aris Gerakis, nov. 1998: Fixed error caused by hardware arithmetic
!   Aris Gerakis, july 1999: Now all borderline cases are valid
!--------------------------------------------------------------------------
!   Glossary:
!   function inpoly: true=inside, false=outside (is target point inside
!                    a 2D polygon?)
!   poly(*,2):  polygon points, [0]=x, [1]=y
!   npoints: number of points in polygon
!   xt: x (horizontal) of target point
!   yt: y (vertical) of target point
!--------------------------------------------------------------------------
!
! declare arguments  
!
      INTEGER :: npoints
      INTEGER :: xt
      INTEGER :: yt 
!
      REAL(KIND=r8)    :: poly(7, 2)
!
! local variables
!
      REAL(KIND=r8)    :: xnew
      REAL(KIND=r8)    :: ynew
      REAL(KIND=r8)    :: xold
      REAL(KIND=r8)    :: yold
      REAL(KIND=r8)    :: x1
      REAL(KIND=r8)    :: y1
      REAL(KIND=r8)    :: x2
      REAL(KIND=r8)    :: y2
!
      INTEGER ::  i
!
      LOGICAL :: inside
      LOGICAL :: on_border

      inside = .false.
      on_border = .false.
!
      IF (npoints .lt. 3)  THEN
        inpoly = .false.
        RETURN
      END IF
!
      xold = poly(npoints,1)
      yold = poly(npoints,2)

      DO i = 1 , npoints
        xnew = poly(i,1)
        ynew = poly(i,2)

        IF (xnew .gt. xold)  THEN
          x1 = xold
          x2 = xnew
          y1 = yold
          y2 = ynew
        ELSE
          x1 = xnew
          x2 = xold
          y1 = ynew
          y2 = yold
        END IF

! the outer IF is the 'straddle' test and the 'vertical border' test.
! the inner IF is the 'non-vertical border' test and the 'north' test.  

! the first statement checks whether a north pointing vector crosses  
! (stradles) the straight segment.  There are two possibilities, depe-
! nding on whether xnew < xold or xnew > xold.  The '<' is because edge 
! must be "open" at left, which is necessary to keep correct count when 
! vector 'licks' a vertix of a polygon.  

        IF ((xnew .lt. xt .and. xt .le. xold)   &
           .or. (.not. xnew .lt. xt .and.       &
           .not. xt .le. xold)) THEN
!
! the test point lies on a non-vertical border:
!
          IF ((yt-y1)*(x2-x1) .eq. (y2-y1)*(xt-x1)) THEN
              
	       on_border = .true. 
!
! check if segment is north of test point.  If yes, reverse the 
! value of INSIDE.  The +0.001 was necessary to avoid errors due   
! arithmetic (e.g., when clay = 98.87 and sand = 1.13):   
!
          ELSE IF ((yt-y1)*(x2-x1) .lt. (y2-y1)*(xt-x1) + 0.001) THEN
          
	    inside = .not.inside ! cross a segment
          
	  END IF
!
! this is the rare case when test point falls on vertical border or  
! left edge of non-vertical border. The left x-coordinate must be  
! common.  The slope requirement must be met, but also point must be
! between the lower and upper y-coordinate of border segment.  There 
! are two possibilities,  depending on whether ynew < yold or ynew > 
! yold:
!
        ELSE IF ((xnew .eq. xt .or. xold .eq. xt)       &
                 .and. (yt-y1)*(x2-x1) .eq.             &
                 (y2-y1)*(xt-x1) .and. ((ynew .le. yt   &
                 .and. yt .le. yold) .or.               &
                 (.not. ynew .lt. yt .and. .not. yt .lt. yold))) THEN
       
          on_border = .true. 
       
        END IF
!
        xold = xnew
        yold = ynew
!
        END DO!  DO i = 1 , npoints  
!
! If test point is not on a border, the function result is the last state 
! of INSIDE variable.  Otherwise, INSIDE doesn't matter.  The point is
! inside the polygon if it falls on any of its borders:
!
      IF (.not. on_border) THEN
         inpoly = inside
      ELSE
         inpoly = .true.
      END IF
!
      RETURN
      END FUNCTION inpoly

 
END MODULE Sfc_Ibis_Vegetation
