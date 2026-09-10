MODULE Rad_CliRadLW
  USE Constants, ONLY :  r8

  IMPLICIT NONE
  PRIVATE

  !c-----parameters defining the size of the pre-computed tables for
  !c     transmittance using table look-up.

  !c     "nx" is the number of intervals in pressure
  !c     "no" is the number of intervals in o3 amount
  !c     "nc" is the number of intervals in co2 amount
  !c     "nh" is the number of intervals in h2o amount

  INTEGER, PARAMETER :: nx=26
  INTEGER, PARAMETER :: no=21
  INTEGER, PARAMETER :: nc=24
  INTEGER, PARAMETER :: nh=31
  INTEGER, PARAMETER :: nt=7

  REAL(KIND=r8) :: c1 (nx,nc,nt),c2 (nx,nc,nt),c3 (nx,nc,nt)
  REAL(KIND=r8) :: o1 (nx,no,nt),o2 (nx,no,nt),o3 (nx,no,nt)
  REAL(KIND=r8) :: h11(nx,nh,nt),h12(nx,nh,nt),h13(nx,nh,nt)
  REAL(KIND=r8) :: h21(nx,nh,nt),h22(nx,nh,nt),h23(nx,nh,nt)
  REAL(KIND=r8) :: h71(nx,nh,nt),h72(nx,nh,nt),h73(nx,nh,nt)


  PUBLIC ::  cliradlw,InitCliRadLW
CONTAINS
  SUBROUTINE InitCliRadLW()

    INTEGER :: ip,it,iw
    !-----copy tables to enhance the speed of co2 (band 3), o3 (band5),
    !  and h2o (bands 1, 2, and 7 only) transmission calculations
    !  using table look-up.

    LOGICAL, PARAMETER :: first=.TRUE.
    LOGICAL, PARAMETER :: high=.TRUE.

    !
    !-----------------------------------------------------------------------
    !
    !  Functions:
    !
    !-----------------------------------------------------------------------
    !
    !  REAL :: expmn
    !
    !-----------------------------------------------------------------------
    !
    !  Include files:
    !
    !-----------------------------------------------------------------------
    !
    !  include "h2o.tran3"
    !  include "co2.tran3"
    !  include "o3.tran3"

    !
    !-----------------------------------------------------------------------
    !
    !  Save variables:
    !
    !-----------------------------------------------------------------------
    !
    !  save c1,c2,c3,o1,o2,o3
    !  save h11,h12,h13,h21,h22,h23,h71,h72,h73
    !

    CALL raddata()

    IF (first) THEN

       !-----tables co2 and h2o are only used with 'high' option

       IF (high) THEN

          DO iw=1,nh
             DO ip=1,nx
                h11(ip,iw,1)=1.0_r8-h11(ip,iw,1)
                h21(ip,iw,1)=1.0_r8-h21(ip,iw,1)
                h71(ip,iw,1)=1.0_r8-h71(ip,iw,1)
             END DO
          END DO

          DO iw=1,nc
             DO ip=1,nx
                c1(ip,iw,1)=1.0_r8-c1(ip,iw,1)
             END DO
          END DO

          !-----tables are replicated to avoid memory bank conflicts

          DO it=2,nt
             DO iw=1,nc
                DO ip=1,nx
                   c1 (ip,iw,it)= c1(ip,iw,1)
                   c2 (ip,iw,it)= c2(ip,iw,1)
                   c3 (ip,iw,it)= c3(ip,iw,1)
                END DO
             END DO
             DO iw=1,nh
                DO ip=1,nx
                   h11(ip,iw,it)=h11(ip,iw,1)
                   h12(ip,iw,it)=h12(ip,iw,1)
                   h13(ip,iw,it)=h13(ip,iw,1)
                   h21(ip,iw,it)=h21(ip,iw,1)
                   h22(ip,iw,it)=h22(ip,iw,1)
                   h23(ip,iw,it)=h23(ip,iw,1)
                   h71(ip,iw,it)=h71(ip,iw,1)
                   h72(ip,iw,it)=h72(ip,iw,1)
                   h73(ip,iw,it)=h73(ip,iw,1)
                END DO
             END DO
          END DO

       END IF

       !-----always use table look-up for ozone transmittance

       DO iw=1,no
          DO ip=1,nx
             o1(ip,iw,1)=1.0_r8-o1(ip,iw,1)
          END DO
       END DO

       DO it=2,nt
          DO iw=1,no
             DO ip=1,nx
                o1 (ip,iw,it)= o1(ip,iw,1)
                o2 (ip,iw,it)= o2(ip,iw,1)
                o3 (ip,iw,it)= o3(ip,iw,1)
             END DO
          END DO
       END DO

       !       first=.FALSE.

    END IF


  END SUBROUTINE InitCliRadLW

  SUBROUTINE cliradlw  ( &
                                ! Model Info and flags
       m         , &
       np        , &
                                ! Atmospheric fields
       pl20      , &            !->   level pressure (pl)                               mb      m*(np+1)
       Pint      , &
       tl        , &            !->   layer temperature (ta)                            k       m*np
       ql        , &            !->   layer specific humidity (wa)                      g/g     m*np
       o3l       , &            !->   layer ozone mixing ratio by mass (oa)             g/g     m*np
       co2l      , &
       tgg       , &            !->   land or ocean surface temperature (tg)            k       m*ns
                                ! LW Radiation fields 
       ulwclr    , &            !->   net downward flux, clear-sky (flc)             w/m**2  m*(np+1)
       ulwtop    , &            !->   net downward flux, all-sky   (flx)             w/m**2  m*(np+1)
       atlclr    , &            !->   Heating rate in clear case (K/s) compute cooling rate profile clear-sky (flc)   
       atl       , &            !->   Heating rate (K/s) compute cooling rate profile all-sky   (flx)
       rsclr     , &            !->   Net surface flux in clear case (W/m2 )specify output fluxes for lwrad clear-sky (flc)  
       rs        , &            !->   Net surface flux (W/m2 ) specify output fluxes for lwrad all-sky   (flx)
       dlwclr    , &            !->   Downward flux at surface in clear case (W/m2 )
       dlwbot    , &
                                ! Cloud field and Microphysics
       cld       , &
       clu         ) 

    IMPLICIT NONE
    !tar
    !tarINTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
    !tar
    ! atlclr....Heating rate in clear case (K/s)
    ! atl.......Heating rate (K/s)
    ! rsclr.....Net surface flux in clear case (W/m2 )
    ! rs........Net surface flux
    ! dlwclr....Downward flux at surface in clear case (W/m2 )
    ! dlwbot....Downward flux at surface (W/m2 )

    ! Model Info and flags
    INTEGER, INTENT(in) :: m
    INTEGER, INTENT(in) :: np

    ! Atmospheric fields
    REAL(KIND=r8), INTENT(in) :: pl20(m,np)
    REAL(KIND=r8), INTENT(in) :: Pint(m,np+1)  
    REAL(KIND=r8), INTENT(in) :: tl  (m,np)
    REAL(KIND=r8), INTENT(in) :: ql  (m,np)
    REAL(KIND=r8), INTENT(in) :: o3l (m,np)
    REAL(KINd=r8), INTENT(IN) :: co2l(m,np)
    REAL(KIND=r8), INTENT(in) :: tgg (m)

    ! LW Radiation fields 
    REAL(KIND=r8), INTENT(out) :: ulwclr (m)
    REAL(KIND=r8), INTENT(out) :: ulwtop (m)
    REAL(KIND=r8), INTENT(out) :: rsclr  (m)
    REAL(KIND=r8), INTENT(out) :: rs     (m)
    REAL(KIND=r8), INTENT(out) :: dlwclr (m)
    REAL(KIND=r8), INTENT(out) :: dlwbot (m)
    REAL(KIND=r8), INTENT(out) :: atlclr (m,np)
    REAL(KIND=r8), INTENT(out) :: atl    (m,np)

    ! Cloud field and Microphysics
    REAL(KIND=r8),    INTENT(in) :: cld  (m,np)
    REAL(KIND=r8),    INTENT(in) :: clu  (m,np)

    !c-----working parameters

    INTEGER  :: np1, i,j,k 

    INTEGER, PARAMETER  :: ns=1

    REAL(KIND=r8)   :: pl(m,np+2)      !c   level pressure (pl)                               mb      m*(np+1)
    REAL(KIND=r8)   :: ta(m,np+1)      !c   layer temperature (ta)                            k       m*np
    REAL(KIND=r8)   :: wa(m,np+1)      !c   layer specific humidity (wa)                      g/g     m*np
    REAL(KIND=r8)   :: oa(m,np+1)

    REAL(KIND=r8)   :: tb   (m)        
    REAL(KIND=r8)   :: fs   (m,ns)     
    REAL(KIND=r8)   :: tg   (m,ns)     
    REAL(KIND=r8)   :: eg   (m,ns,10)  
    REAL(KIND=r8)   :: taucl(m,np+1,3) 
    REAL(KIND=r8)   :: tauc(m,np+1) 
    REAL(KIND=r8)   :: fcld (m,np+1)   
    REAL(KIND=r8)   :: flx  (m,np+2)   
    REAL(KIND=r8)   :: flc  (m,np+2) 
    REAL(KIND=r8)   :: dfdts(m,np+2) 
    REAL(KIND=r8)   :: sfcem(m)     
    REAL(KIND=r8)   :: dp   (m,np+1)   
    REAL(KIND=r8)   :: radlwin (m) 
    REAL(KIND=r8)   ::  co2(m,np+1)

    INTEGER, PARAMETER :: rlwopt=1

    !  is, k=1 is for top while k=nz is at the surface.

    !  Calculate water vapor, co2, and o3 transmittances using table
    !  look-up. rlwopt = 0, or high=.false.

    LOGICAL :: high

    IF ( rlwopt == 0 ) THEN
       high = .FALSE.
    ELSE
       high = .TRUE.
    END IF


    !c ---- initialize output fluxes
    tb   =0.0_r8
    fs   =0.0_r8
    tg   =0.0_r8
    eg   =0.0_r8
    taucl=0.0_r8
    tauc=0.0_r8
    dfdts=0.0_r8
    fcld =0.0_r8
    flx  =0.0_r8
    flc  =0.0_r8
    sfcem=0.0_r8
    dp   =0.0_r8


    ulwclr = 0.0_r8
    ulwtop = 0.0_r8
    atlclr = 0.0_r8
    atl    = 0.0_r8
    rsclr  = 0.0_r8
    rs     = 0.0_r8
    dlwclr = 0.0_r8 
    dlwbot = 0.0_r8

    !ctar     REAL(KIND=r8)      ::tck
    !ctar     REAL(KIND=r8)      ::tbk

    !c-----specify new ozone concentration

    !ctar       co2val=300.0_r8

    !co2=co2l!*1.0e-6_r8

    !c-----specify input errays in irrad

    !ctar      do k=1,np
    !ctar       do i=1,m
    !ctar         pl(i,k+1)=pl20(i,k)
    !ctar         ta(i,k)=tl(i,k)
    !ctar         wa(i,k)=ql(i,k)
    !ctar         oa(i,k)=o3l(i,k)
    !ctar       enddo
    !ctar      enddo

    DO k=1,np+1
       DO i=1,m          
          pl(i,k+1)=Pint(i,k)
       ENDDO
    ENDDO

    DO k=1,np
       DO i=1,m          
          ta(i,k+1)=tl(i,k)
          wa(i,k+1)=ql(i,k)
          oa(i,k+1)=o3l(i,k)
          co2(i,k+1) =co2l(i,k)
       ENDDO
    ENDDO
    DO i=1,m
       Co2(i,1) =co2l(i,1)
    END DO
    !c adding of 2 layers at top 0 mb --1mb ---

    np1=np+1

    !c-----high pressure level

    DO i=1,m
       !tar         pl(i,1)=1.
       pl(i,1)=0.0_r8            !mb
       pl(i,2)=pl(i,3)*0.5_r8    !mb
       ta(i,1)=tl(i,1)-20.0_r8
       !tar                           
       wa(i,1)=1.0e-12_r8        
       IF (o3l(i,1).GE.2.5e-6_r8) THEN  
          oa(i,1)=o3l(i,1)
       ELSE 
          oa(i,1)=2.5e-6_r8        
       ENDIF
       !ctar                                    
    ENDDO

    !c-----specify near surface air temperature tb as mean of lowest layer and Tg

    DO i=1,m
       tb(i)=(ta(i,np+1)+tgg(i))*0.5_r8
       !c          tb(i)=tgg(i)
    ENDDO

    !c--------specify ground properties

    DO i=1,m
       fs(i,ns)=1.0_r8
       tg(i,ns)=tgg(i)
       DO j=1,10
          eg(i,ns,j)=1.0_r8 
       ENDDO
    ENDDO

    !c-----specify cloud optical depth

    DO k=1,np1
       DO i=1,m
          dp(i,k)=pl(i,k+1)-pl(i,k)
          taucl(i,k,1)=0.0_r8
          taucl(i,k,2)=0.0_r8
          taucl(i,k,3)=0.0_r8
          fcld(i,k)=0.0_r8
       ENDDO
    ENDDO

    DO k=9,np
       DO i=1,m
          IF (cld(i,k).GT.clu(i,k).AND.cld(i,k).GT.0.1_r8) THEN 
             fcld(i,k+1)=cld(i,k)

             IF (ta(i,k+1).GT.253.0_r8) THEN
                taucl(i,k+1,2)=0.05_r8*dp(i,k+1)*1.16_r8
             ELSE
                taucl(i,k+1,1)=0.025_r8*dp(i,k+1)*1.16_r8
             ENDIF

          ENDIF
          IF (clu(i,k).GE.cld(i,k).AND.clu(i,k).GT.0.1_r8) THEN 
             fcld(i,k+1)=clu(i,k)
             taucl(i,k+1,2)=0.05_r8*dp(i,k+1)*1.16_r8
          ENDIF
       ENDDO
    ENDDO
    DO k=1,np
       DO i=1,m
          tauc(i,k)=taucl(i,k,1)+taucl(i,k,2)+taucl(i,k,3)
       END DO
    END DO
    CALL irrad( &
         m                    , &  ! INTEGER      , INTENT(IN    ) :: m
         np1                  , &  ! INTEGER      , INTENT(IN    ) :: np
         tauc    (1:m,1:np1)  , &  ! REAL(KIND=r8), INTENT(IN    ) :: taucl   (m,np)
         fcld    (1:m,1:np1)  , &  ! REAL(KIND=r8), INTENT(IN    ) :: ccld    (m,np)
         pl      (1:m,1:np1+1), &  ! REAL(KIND=r8), INTENT(IN    ) :: pl      (m,np+1)
         ta      (1:m,1:np1)  , &  ! REAL(KIND=r8), INTENT(IN    ) :: ta      (m,np)
         wa      (1:m,1:np1)  , &  ! REAL(KIND=r8), INTENT(IN    ) :: wa      (m,np)
         oa      (1:m,1:np1)  , &  ! REAL(KIND=r8), INTENT(IN    ) :: oa      (m,np)
         co2     (1:m,1:np1)  , &  ! REAL(KIND=r8), INTENT(IN    ) :: co2
         tb      (1:m)        , &  ! REAL(KIND=r8), INTENT(IN    ) :: ts      (m)
         high                 , &  ! LOGICAL      , INTENT(IN    ) :: high
         radlwin (1:m)        , &  ! REAL(KIND=r8), INTENT(OUT   ) :: radlwin (m) 
         flx     (1:m,1:np1+1), &  ! REAL(KIND=r8), INTENT(OUT   ) :: flx     (m,np+1)
         flc     (1:m,1:np1+1), &  ! REAL(KIND=r8), INTENT(OUT   ) :: flc     (m,np+1)
         dfdts   (1:m,1:np1+1), &  ! REAL(KIND=r8), INTENT(OUT   ) :: dfdts   (m,np+1)
         sfcem   (1:m)          )  ! REAL(KIND=r8), INTENT(OUT   ) :: sfcem   (m)
    !c-----compute cooling rate profile
    !ctar
    DO k=1,np
       DO i=1,m
          !         atl(i,k)=-(flx(i,k+2)-flx(i,k+1))*8.441874/  &
          atl(i,k)=-(flx(i,k+2)-flx(i,k+1))*9.7706e-5_r8/  &                  
               dp(i,k+1)
          !         atlclr(i,k)=-(flc(i,k+2)-flc(i,k+1))*8.441874/  &
          atlclr(i,k)=-(flc(i,k+2)-flc(i,k+1))*9.7706e-5_r8/  &         
               dp(i,k+1)          
       ENDDO
    ENDDO

    DO k=1,np
       DO i=1,m
          IF(ABS(atlclr(i,k)) .LT. 1.e-22_r8) atlclr(i,k)=0._r8
          IF(ABS(atl(i,k))    .LT. 1.e-22_r8) atl(i,k)   =0._r8

          !        atlclr(i,k) = atlclr(i,k) * 1.1574e-5
          !        atl(i,k) = atl(i,k) * 1.1574e-5     
       ENDDO
    ENDDO

    !c-----specify output fluxes for lwrad
    !c   emission by the surface (sfcem)                w/m**2     m
    DO i=1,m         
       ulwclr(i) = -flc(i,1)
       ulwtop(i) = -flx(i,1)
       rsclr(i)  = -flc(i,np+2)
       rs(i)     = -flx(i,np+2)
       dlwclr(i) = flc(i,np+2)-sfcem(i)
       dlwbot(i) = flx(i,np+2)-sfcem(i)
    ENDDO


  END  SUBROUTINE cliradlw


  !
  !
  !##################################################################
  !##################################################################
  !######                                                      ######
  !######                SUBROUTINE IRRAD                      ######
  !######                                                      ######
  !######                     Developed by                     ######
  !######                                                      ######
  !######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
  !######                                                      ######
  !######     Center for Analysis and Prediction of Storms     ######
  !######               University of Oklahoma                 ######
  !######                                                      ######
  !##################################################################
  !##################################################################
  !


  SUBROUTINE irrad( &
       m       , & ! INTEGER, INTENT(IN   ) :: m
       np      , & ! INTEGER, INTENT(IN   ) :: np
       taucl   , & ! REAL(KIND=r8), INTENT(IN    ) :: taucl(m,np)
       ccld    , & ! REAL(KIND=r8), INTENT(IN    ) :: ccld (m,np)
       pl      , & ! REAL(KIND=r8), INTENT(IN    ) :: pl   (m,np+1)
       ta      , & ! REAL(KIND=r8), INTENT(IN    ) :: ta   (m,np)
       wa      , & ! REAL(KIND=r8), INTENT(IN    ) :: wa   (m,np)
       oa      , & ! REAL(KIND=r8), INTENT(IN    ) :: oa   (m,np)
       co2     , & ! REAL(KIND=r8), INTENT(IN    ) :: co2
       ts      , & ! REAL(KIND=r8), INTENT(IN    ) :: ts   (m)
       high    , & ! LOGICAL      , INTENT(IN    ) :: high
       radlwin , & ! REAL(KIND=r8), INTENT(OUT   ) :: radlwin (m) 
       flx     , & ! REAL(KIND=r8), INTENT(OUT   ) :: flx     (m,np+1)
       flc     , & ! REAL(KIND=r8), INTENT(OUT   ) :: flc     (m,np+1)
       dfdts   , & ! REAL(KIND=r8), INTENT(OUT   ) :: dfdts   (m,np+1)
       st4       ) ! REAL(KIND=r8), INTENT(OUT   ) :: st4     (m)

    !
    !-----------------------------------------------------------------------
    !
    !  PURPOSE:
    !
    !  Calculate IR fluxes due to water vapor, co2, and o3. Clouds in
    !  different layers are assumed randomly overlapped.
    !
    !-----------------------------------------------------------------------
    !
    !  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
    !          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
    !              Chou (1996)
    !
    !  MODIFICATION HISTORY:
    !
    !  03/15/1996 (Yuhe Liu)
    !  Adopted the original code and formatted it in accordance with the
    !  ARPS coding standard.
    !
    !-----------------------------------------------------------------------
    !
    !  ORIGINAL COMMENTS:
    !
    !******************** CLIRAD IR1  Date: Oct. 17, 1994 ****************
    !*********************************************************************
    !
    ! This routine computes ir fluxes due to water vapor, co2, and o3.
    !   Clouds in different layers are assumed randomly overlapped.
    !
    ! This is a vectorized code.  It computes fluxes simultaneously for
    !   (m x n) soundings, which is a subset of (m x ndim) soundings.
    !   In a global climate model, m and ndim correspond to the numbers of
    !   grid boxes in the zonal and meridional directions, respectively.
    !
    ! Detailed description of the radiation routine is given in
    !   Chou and Suarez (1994).
    !
    ! There are two options for computing cooling rate profiles.
    !
    !   if high = .true., transmission functions in the co2, o3, and the
    !   three water vapor bands with strong absorption are computed using
    !   table look-up.  cooling rates are computed accurately from the
    !   surface up to 0.01 mb.
    !   if high = .false., transmission functions are computed using the
    !   k-distribution method with linear pressure scaling.  cooling rates
    !   are not calculated accurately for pressures less than 20 mb.
    !   the computation is faster with high=.false. than with high=.true.
    !
    ! The IR spectrum is divided into eight bands:
    !
    !   bnad     wavenumber (/cm)   absorber         method
    !
    !    1           0 - 340           h2o            K/T
    !    2         340 - 540           h2o            K/T
    !    3         540 - 800       h2o,cont,co2       K,S,K/T
    !    4         800 - 980       h2o,cont           K,S
    !    5         980 - 1100      h2o,cont,o3        K,S,T
    !    6        1100 - 1380      h2o,cont           K,S
    !    7        1380 - 1900          h2o            K/T
    !    8        1900 - 3000          h2o            K
    !
    ! Note : "h2o" for h2o line absorption
    !     "cont" for h2o continuum absorption
    !     "K" for k-distribution method
    !     "S" for one-parameter temperature scaling
    !     "T" for table look-up
    !
    ! The 15 micrometer region (540-800/cm) is further divided into
    !   3 sub-bands :
    !
    !   subbnad   wavenumber (/cm)
    !
    !    1          540 - 620
    !    2          620 - 720
    !    3          720 - 800
    !
    !---- Input parameters                               units    size
    !
    !   number of soundings in zonal direction (m)        n/d      1
    !   number of soundings in meridional direction (n)   n/d      1
    !   maximum number of soundings in
    !              meridional direction (ndim)            n/d      1
    !   number of atmospheric layers (np)                 n/d      1
    !   cloud optical thickness (taucl)                   n/d     m*ndim*np
    !   cloud cover (ccld)                              fraction  m*ndim*np
    !   level pressure (pl)                               mb      m*ndim*(np+1)
    !   layer temperature (ta)                            k       m*ndim*np
    !   layer specific humidity (wa)                      g/g     m*ndim*np
    !   layer ozone mixing ratio by mass (oa)             g/g     m*ndim*np
    !   surface temperature (ts)                          k       m*ndim
    !   co2 mixing ratio by volumn (co2)                  pppv     1
    !   high                                                       1
    !
    ! pre-computed tables used in table look-up for transmittance calculations:
    !
    !   c1 , c2, c3: for co2 (band 3)
    !   o1 , o2, o3: for  o3 (band 5)
    !   h11,h12,h13: for h2o (band 1)
    !   h21,h22,h23: for h2o (band 2)
    !   h71,h72,h73: for h2o (band 7)
    !
    !---- output parameters
    !
    !   net downward flux, all-sky   (flx)             w/m**2     m*ndim*(np+1)
    !   net downward flux, clear-sky (flc)             w/m**2     m*ndim*(np+1)
    !   sensitivity of net downward flux
    !    to surface temperature (dfdts)             w/m**2/k   m*ndim*(np+1)
    !   emission by the surface (st4)                  w/m**2     m*ndim
    !
    ! Notes:
    !
    !   (1)  Water vapor continuum absorption is included in 540-1380 /cm.
    !   (2)  Scattering by clouds is not included.
    !   (3)  Clouds are assumed "gray" bodies.
    !   (4)  The diffuse cloud transmission is computed to be exp(-1.66*taucl).
    !   (5)  If there are no clouds, flx=flc.
    !   (6)  plevel(1) is the pressure at the top of the model atmosphere, and
    !     plevel(np+1) is the surface pressure.
    !
    !    ARPS note: pl was replaced by pa at scalar points (layers)
    !
    !   (7)  Downward flux is positive, and upward flux is negative.
    !   (8)  dfdts is always negative because upward flux is defined as negative.
    !   (9)  For questions and coding errors, please contact with Ming-Dah Chou,
    !     Code 913, NASA/Goddard Space Flight Center, Greenbelt, MD 20771.
    !     Phone: 301-286-4012, Fax: 301-286-1759,
    !     e-mail: chou@climate.gsfc.nasa.gov
    !
    !-----parameters defining the size of the pre-computed tables for transmittance
    !  calculations using table look-up.
    !
    !  "nx" is the number of intervals in pressure
    !  "no" is the number of intervals in o3 amount
    !  "nc" is the number of intervals in co2 amount
    !  "nh" is the number of intervals in h2o amount
    !  "nt" is the number of copies to be made from the pre-computed
    !       transmittance tables to reduce "memory-bank conflict"
    !       in parallel machines and, hence, enhancing the speed of
    !       computations using table look-up.
    !       If such advantage does not exist, "nt" can be set to 1.
    !***************************************************************************
    !
    !fpp$ expand (expmn)
    !!dir$ inline always expmn
    !*$*  inline routine (expmn)
    !
    !-----------------------------------------------------------------------
    !

    !
    !-----------------------------------------------------------------------
    !
    !  Variable Declarations.
    !
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: m
    INTEGER, INTENT(IN   ) :: np



    !---- input parameters ------

    REAL(KIND=r8), INTENT(IN   ) :: taucl(m,np)
    REAL(KIND=r8), INTENT(IN   ) :: ccld (m,np)
    REAL(KIND=r8), INTENT(IN   ) :: pl   (m,np+1)
    REAL(KIND=r8), INTENT(IN   ) :: ta   (m,np)
    REAL(KIND=r8), INTENT(IN   ) :: wa   (m,np)
    REAL(KIND=r8), INTENT(IN   ) :: oa   (m,np)
    REAL(KIND=r8), INTENT(IN   ) :: co2  (m,np) 
    REAL(KIND=r8), INTENT(IN   ) :: ts   (m)

    LOGICAL      , INTENT(IN   ) :: high

    !---- output parameters ------

    REAL(KIND=r8), INTENT(OUT   ) :: radlwin (m) 
    REAL(KIND=r8), INTENT(OUT   ) :: flx     (m,np+1)
    REAL(KIND=r8), INTENT(OUT   ) :: flc     (m,np+1)
    REAL(KIND=r8), INTENT(OUT   ) :: dfdts   (m,np+1)
    REAL(KIND=r8), INTENT(OUT   ) :: st4     (m)
    !
    !-----------------------------------------------------------------------
    !
    !  Temporary arrays
    !
    !-----------------------------------------------------------------------
    !
    REAL(KIND=r8) :: fclr(m)
    REAL(KIND=r8) :: dbs(m)
    REAL(KIND=r8) :: trant(m)
    REAL(KIND=r8) :: th2o(m,6)
    REAL(KIND=r8) :: tcon(m,3)
    REAL(KIND=r8) :: tco2(m,6,2)

    REAL(KIND=r8) :: pa(m,np)
    REAL(KIND=r8) :: dt(m,np)
    REAL(KIND=r8) :: sh2o(m,np+1)
    REAL(KIND=r8) :: swpre(m,np+1)
    REAL(KIND=r8) :: swtem(m,np+1)
    REAL(KIND=r8) :: sco3(m,np+1)
    REAL(KIND=r8) :: scopre(m,np+1)
    REAL(KIND=r8) :: scotem(m,np+1)
    REAL(KIND=r8) :: dh2o(m,np)
    REAL(KIND=r8) :: dcont(m,np)
    REAL(KIND=r8) :: dco2(m,np)
    REAL(KIND=r8) :: do3(m,np)
    REAL(KIND=r8) :: flxu(m,np+1)
    REAL(KIND=r8) :: flxd(m,np+1)
    REAL(KIND=r8) :: clr(m,0:np+1)
    REAL(KIND=r8) :: blayer(m,0:np+1)

    REAL(KIND=r8) :: h2oexp(m,np,6)
    REAL(KIND=r8) :: conexp(m,np,3)

    REAL(KIND=r8) :: co2exp(m,np,6,2)
    !
    !-----------------------------------------------------------------------
    !
    !  Misc. local variables
    !
    !-----------------------------------------------------------------------
    !
    LOGICAL :: oznbnd
    LOGICAL :: co2bnd
    LOGICAL :: h2otbl
    LOGICAL :: conbnd


    REAL(KIND=r8) :: xx, w1,p1,dwe,dpe

    !---- static data -----

    !-----the following coefficients (table 2 of chou and suarez, 1995)
    !  are for computing spectrally integtrated planck fluxes of
    !  the 8 bands using eq. (22)
    !  INTEGER, PARAMETER :: mw(8) =RESHAPE(SOURCE=(/6,6,8,6,6,8,6,16/),SHAPE=(/8/))

    REAL(KIND=r8) :: cb(5,8)   =RESHAPE(SOURCE=(/&
         -2.6844E-1_r8,-8.8994E-2_r8, 1.5676E-3_r8,-2.9349E-6_r8, 2.2233E-9_r8,&
         3.7315E+1_r8,-7.4758E-1_r8, 4.6151E-3_r8,-6.3260E-6_r8, 3.5647E-9_r8,&
         3.7187E+1_r8,-3.9085E-1_r8,-6.1072E-4_r8, 1.4534E-5_r8,-1.6863E-8_r8,&
         -4.1928E+1_r8, 1.0027E+0_r8,-8.5789E-3_r8, 2.9199E-5_r8,-2.5654E-8_r8,&
         -4.9163E+1_r8, 9.8457E-1_r8,-7.0968E-3_r8, 2.0478E-5_r8,-1.5514E-8_r8,&
         -1.0345E+2_r8, 1.8636E+0_r8,-1.1753E-2_r8, 2.7864E-5_r8,-1.1998E-8_r8,&
         -6.9233E+0_r8,-1.5878E-1_r8, 3.9160E-3_r8,-2.4496E-5_r8, 4.9301E-8_r8,&
         1.1483E+2_r8,-2.2376E+0_r8, 1.6394E-2_r8,-5.3672E-5_r8, 6.6456E-8_r8/),SHAPE=(/5,8/))

    INTEGER :: i
    INTEGER :: k
    INTEGER :: ib
    INTEGER :: ik
    INTEGER :: iq
    INTEGER :: isb
    INTEGER :: k1
    INTEGER :: k2

    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !
    !  Beginning of executable code...
    !
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !
    !-----compute layer pressure (pa) and layer temperature minus 250K (dt)

    DO k=1,np
       DO i=1,m
          dt(i,k)=ta(i,k)-250.0_r8
          pa(i,k)=0.5_r8*(pl(i,k)+pl(i,k+1))
       END DO
    END DO

    !-----compute layer absorber amount

    !  dh2o : water vapor amount (g/cm**2)
    !  dcont: scaled water vapor amount for continuum absorption (g/cm**2)
    !  dco2 : co2 amount (cm-atm)stp
    !  do3  : o3 amount (cm-atm)stp
    !  the factor 1.02 is equal to 1000/980
    !  factors 789 and 476 are for unit conversion
    !  the factor 0.001618 is equal to 1.02/(.622*1013.25)
    !  the factor 6.081 is equal to 1800/296

    DO k=1,np
       DO i=1,m
          dh2o(i,k) = 1.02_r8*wa(i,k)*(pl(i,k+1)-pl(i,k))+1.e-10_r8
          dco2(i,k) = 789.0_r8*co2(i,k)*(pl(i,k+1)-pl(i,k))+1.e-10_r8
          do3 (i,k) = 476.0_r8*oa(i,k)*(pl(i,k+1)-pl(i,k))+1.e-10_r8

          !-----compute scaled water vapor amount for h2o continuum absorption
          !  following eq. (43).

          xx=pa(i,k)*0.001618_r8*wa(i,k)*wa(i,k)                       &
               *(pl(i,k+1)-pl(i,k))
          dcont(i,k) = xx*expmn(1800.0_r8/ta(i,k)-6.081_r8)+1.e-10_r8

          !-----compute effective cloud-free fraction, clr, for each layer.
          !  the cloud diffuse transmittance is approximated by using a
          !  diffusivity factor of 1.66.

          clr(i,k)=1.0_r8-(ccld(i,k)*(1.0_r8-expmn(-1.66_r8*taucl(i,k))))

       END DO
    END DO

    !-----compute column-integrated h2o amoumt, h2o-weighted pressure
    !  and temperature.  it follows eqs. (37) and (38).

    IF (high) THEN

       CALL column(m                    , &! INTEGER      , INTENT(IN   ) :: m !   number of soundings in zonal direction (m)
            np                   , &! INTEGER      , INTENT(IN   ) :: np !   number of atmospheric layers (np)
            pa   (1:m,1:np)  , &! REAL(KIND=r8), INTENT(IN   ) :: pa   (m,n,np)  ! layer pressure (pa)
            dt   (1:m,1:np)  , &! REAL(KIND=r8), INTENT(IN   ) :: dt   (m,n,np)  ! layer temperature minus 250K (dt)
            dh2o (1:m,1:np)  , &! REAL(KIND=r8), INTENT(IN   ) :: sabs0(m,n,np)  !   layer absorber amount (sabs0)
            sh2o (1:m,1:np+1), &! REAL(KIND=r8), INTENT(OUT  ) :: sabs (m,n,np+1)!   column-integrated absorber amount (sabs)
            swpre(1:m,1:np+1), &! REAL(KIND=r8), INTENT(OUT  ) :: spre (m,n,np+1)!   column absorber-weighted pressure (spre)
            swtem(1:m,1:np+1)  )! REAL(KIND=r8), INTENT(OUT  ) :: stem (m,n,np+1)!   column absorber-weighted temperature (stem)
    END IF

    !-----the surface (with an index np+1) is treated as a layer filled with
    !  black clouds.

    DO i=1,m
       clr(i,0)    = 1.0_r8
       clr(i,np+1) = 0.0_r8
       st4(i)      = 0.0_r8
    END DO

    !-----initialize fluxes

    DO k=1,np+1
       DO i=1,m
          flx(i,k)  = 0.0_r8
          flc(i,k)  = 0.0_r8
          dfdts(i,k)= 0.0_r8
          flxu(i,k) = 0.0_r8
          flxd(i,k) = 0.0_r8
       END DO
    END DO

    !-----integration over spectral bands

    DO ib=1,8

       !-----if h2otbl, compute h2o (line) transmittance using table look-up.
       !  if conbnd, compute h2o (continuum) transmittance in bands 3, 4, 5 and 6.
       !  if co2bnd, compute co2 transmittance in band 3.
       !  if oznbnd, compute  o3 transmittance in band 5.

       h2otbl=high.AND.(ib == 1 .OR. ib == 2 .OR. ib == 7)
       conbnd=ib >= 3 .AND. ib <= 6
       co2bnd=ib == 3
       oznbnd=ib == 5

       !-----blayer is the spectrally integrated planck flux of the mean layer
       !  temperature derived from eq. (22)
       !  the fitting for the planck flux is valid in the range 160-345 K.

       DO k=1,np
          DO i=1,m
             blayer(i,k)=ta(i,k)*(ta(i,k)*(ta(i,k)                 &
                  *(ta(i,k)*cb(5,ib)+cb(4,ib))+cb(3,ib))         &
                  +cb(2,ib))+cb(1,ib)
          END DO
       END DO

       !-----the earth's surface, with an index "np+1", is treated as a layer

       DO i=1,m
          blayer(i,0)   = 0.0_r8
          blayer(i,np+1)=ts(i)*(ts(i)*(ts(i)                      &
               *(ts(i)*cb(5,ib)+cb(4,ib))+cb(3,ib))          &
               +cb(2,ib))+cb(1,ib)

          !-----dbs is the derivative of the surface planck flux with respect to
          !  surface temperature (eq. 59).

          dbs(i)=ts(i)*(ts(i)*(ts(i)                              &
               *4.0_r8*cb(5,ib)+3.0_r8*cb(4,ib))+2.0_r8*cb(3,ib))+cb(2,ib)

       END DO

       !-----compute column-integrated absorber amoumt, absorber-weighted
       !  pressure and temperature for co2 (band 3) and o3 (band 5).
       !  it follows eqs. (37) and (38).

       !-----this is in the band loop to save storage

       IF( high .AND. co2bnd) THEN
          CALL column(m                     , &! INTEGER      , INTENT(IN   ) :: m !   number of soundings in zonal direction (m)
               np                    , &! INTEGER      , INTENT(IN   ) :: np !   number of atmospheric layers (np)
               pa    (1:m,1:np)  , &! REAL(KIND=r8), INTENT(IN   ) :: pa   (m,n,np)  !   layer pressure (pa)
               dt    (1:m,1:np)  , &! REAL(KIND=r8), INTENT(IN   ) :: dt   (m,n,np)  !   layer temperature minus 250K (dt)
               dco2  (1:m,1:np)  , &! REAL(KIND=r8), INTENT(IN   ) :: sabs0(m,n,np)  !   layer absorber amount (sabs0)
               sco3  (1:m,1:np+1), &! REAL(KIND=r8), INTENT(OUT  ) :: sabs (m,n,np+1)!   column-integrated absorber amount (sabs)
               scopre(1:m,1:np+1), &! REAL(KIND=r8), INTENT(OUT  ) :: spre (m,n,np+1)!   column absorber-weighted pressure (spre)
               scotem(1:m,1:np+1)  )! REAL(KIND=r8), INTENT(OUT  ) :: stem (m,n,np+1)!   column absorber-weighted temperature (stem)
       END IF

       IF(oznbnd) THEN

          CALL column(m                     , &! INTEGER      , INTENT(IN   ) :: m !   number of soundings in zonal direction (m)
               np                    , &! INTEGER      , INTENT(IN   ) :: np !   number of atmospheric layers (np)
               pa    (1:m,1:np)  , &! REAL(KIND=r8), INTENT(IN   ) :: pa   (m,n,np)  !   layer pressure (pa)
               dt    (1:m,1:np)  , &! REAL(KIND=r8), INTENT(IN   ) :: dt   (m,n,np)  !   layer temperature minus 250K (dt)
               do3   (1:m,1:np)  , &! REAL(KIND=r8), INTENT(IN   ) :: sabs0(m,n,np)  !   layer absorber amount (sabs0)
               sco3  (1:m,1:np+1), &! REAL(KIND=r8), INTENT(OUT  ) :: sabs (m,n,np+1)!   column-integrated absorber amount (sabs)
               scopre(1:m,1:np+1), &! REAL(KIND=r8), INTENT(OUT  ) :: spre (m,n,np+1)!   column absorber-weighted pressure (spre)
               scotem(1:m,1:np+1)  )! REAL(KIND=r8), INTENT(OUT  ) :: stem (m,n,np+1)!   column absorber-weighted temperature (stem)
       END IF

       !-----compute the exponential terms (eq. 32) at each layer for
       !  water vapor line absorption when k-distribution is used

       IF( .NOT. h2otbl) THEN

          CALL h2oexps(ib                       , &! INTEGER, INTENT(IN   ) :: ib!  spectral band (ib)
               m                        , &! INTEGER, INTENT(IN   ) :: m !  number of grid intervals in zonal direction (m)
               np                       , &! INTEGER, INTENT(IN   ) :: np!  number of layers (np)
               dh2o  (1:m,1:np)     , &! REAL(KIND=r8), INTENT(IN   ) :: dh2o  (m,n,np)     !  layer water vapor amount for line absorption (dh2o)
               pa    (1:m,1:np)     , &! REAL(KIND=r8), INTENT(IN   ) :: pa (m,n,np)     !  layer pressure (pa)
               dt    (1:m,1:np)     , &! REAL(KIND=r8), INTENT(IN   ) :: dt (m,n,np)     !  layer temperature minus 250K (dt)
               h2oexp(1:m,1:np,1:6)   )! REAL(KIND=r8), INTENT(INOUT) :: h2oexp(m,n,np,6)   !  6 exponentials for each layer  (h2oexp)

       END IF

       !-----compute the exponential terms (eq. 46) at each layer for
       !  water vapor continuum absorption

       IF( conbnd) THEN

          CALL conexps(ib                      ,&!INTEGER, INTENT(IN   ) :: ib!  spectral band (ib)
               m                       ,&!INTEGER, INTENT(IN   ) :: m  !  number of grid intervals in zonal direction (m)
               np                      ,&!INTEGER, INTENT(IN   ) :: np !  number of layers (np)
               dcont (1:m,1:np)    ,&!REAL(KIND=r8), INTENT(IN   ) :: dcont (m,n,np)   !  layer scaled water vapor amount for continuum absorption (dcont)
               conexp(1:m,1:np,1:3) )!REAL(KIND=r8), INTENT(INOUT) :: conexp(m,n,np,3)!  1 or 3 exponentials for each layer (conexp)

       END IF


       !-----compute the  exponential terms (eq. 32) at each layer for
       !  co2 absorption

       IF( .NOT.high .AND. co2bnd) THEN

          CALL co2exps(m                            ,& !INTEGER, INTENT(IN   ) :: m !  number of grid intervals in zonal direction (m)
               np                           ,& !INTEGER, INTENT(IN   ) :: np!  number of layers (np)
               dco2  (1:m,1:np)         ,& !REAL(KIND=r8), INTENT(IN   ) :: dco2  (m,n,np)     !  layer co2 amount (dco2)
               pa    (1:m,1:np)         ,& !REAL(KIND=r8), INTENT(IN   ) :: pa    (m,n,np)     !  layer pressure (pa)
               dt    (1:m,1:np)         ,& !REAL(KIND=r8), INTENT(IN   ) :: dt    (m,n,np)     !  layer temperature minus 250K (dt)
               co2exp(1:m,1:np,1:6,1:2)  ) !REAL(KIND=r8), INTENT(INOUT) :: co2exp(m,n,np,6,2) !  6 exponentials for each layer (co2exp)

       END IF

       !-----compute transmittances for regions between levels k1 and k2
       !  and update the fluxes at the two levels.

       DO k1=1,np

          !-----initialize fclr, th2o, tcon, and tco2

          DO i=1,m
             fclr(i)=1.0_r8
          END DO

          !-----for h2o line absorption

          IF(.NOT. h2otbl) THEN
             DO ik=1,6
                DO i=1,m
                   th2o(i,ik)=1.0_r8
                END DO
             END DO
          END IF

          !-----for h2o continuum absorption

          IF (conbnd) THEN
             DO iq=1,3
                DO i=1,m
                   tcon(i,iq)=1.0_r8
                END DO
             END DO
          END IF

          !-----for co2 absorption when using k-distribution method.
          !  band 3 is divided into 3 sub-bands, but sub-bands 3a and 3c
          !  are combined in computing the co2 transmittance.

          IF (.NOT. high .AND. co2bnd) THEN
             DO isb=1,2
                DO ik=1,6
                   DO i=1,m
                      tco2(i,ik,isb)=1.0_r8
                   END DO
                END DO
             END DO
          END IF

          !-----loop over the bottom level of the region (k2)

          DO k2=k1+1,np+1

             DO i=1,m
                trant(i)=1.0_r8
             END DO

             IF(h2otbl) THEN

                w1=-8.0_r8
                p1=-2.0_r8
                dwe=0.3_r8
                dpe=0.2_r8

                !-----compute water vapor transmittance using table look-up

                IF (ib == 1 ) THEN

                   CALL tablup(k1                   , &!INTEGER      ,INTENT(IN   ) :: k1 indices for pressure levels (k1 and k2)
                        k2                   , &!INTEGER      ,INTENT(IN   ) :: k2 indices for pressure levels (k1 and k2)
                        m                    , &!INTEGER      ,INTENT(IN   ) :: m  number of grid intervals in zonal direction (m)
                        np                   , &!INTEGER      ,INTENT(IN   ) :: np number of atmospheric layers (np)
                        nx                   , &!INTEGER      ,INTENT(IN   ) :: nx number of pressure intervals in the table (nx)
                        nh                   , &!INTEGER      ,INTENT(IN   ) :: nh number of absorber amount intervals in the table (nh)
                        nt                   , &!INTEGER      ,INTENT(IN   ) :: nt number of tables copied (nt)
                        sh2o (1:m,1:np+1)    , &!REAL(KIND=r8),INTENT(IN   ) :: sabs(m,n,np+1) column-integrated absorber amount (sabs)
                        swpre(1:m,1:np+1)    , &!REAL(KIND=r8),INTENT(IN   ) :: spre(m,n,np+1) column absorber amount-weighted pressure (spre)
                        swtem(1:m,1:np+1)    , &!REAL(KIND=r8),INTENT(IN   ) :: stem(m,n,np+1) column absorber amount-weighted temperature (stem)
                        w1                   , &!REAL(KIND=r8),INTENT(IN   ) :: w1             first value of absorber amount (log10) in the table (w1)
                        p1                   , &!REAL(KIND=r8),INTENT(IN   ) :: p1             first value of pressure (log10) in the table (p1)
                        dwe                  , &!REAL(KIND=r8),INTENT(IN   ) :: dwe             size of the interval of absorber amount (log10) in the table (dwe)
                        dpe                  , &!REAL(KIND=r8),INTENT(IN   ) :: dpe             size of the interval of pressure (log10) in the table (dpe)
                        h11  (1:nx,1:nh,1:nt), &!REAL(KIND=r8),INTENT(IN   ) :: coef1(nx,nh,nt) pre-computed coefficients (coef1, coef2, and coef3)
                        h12  (1:nx,1:nh,1:nt), &!REAL(KIND=r8),INTENT(IN   ) :: coef2(nx,nh,nt) pre-computed coefficients (coef1, coef2, and coef3)
                        h13  (1:nx,1:nh,1:nt), &!REAL(KIND=r8),INTENT(IN   ) :: coef3(nx,nh,nt) pre-computed coefficients (coef1, coef2, and coef3)
                        trant(1:m)             )!REAL(KIND=r8),INTENT(INOUT) :: tran(m,n)  updated transmittance (tran)

                END IF
                IF (ib == 2 ) THEN

                   CALL tablup(k1                    , &!INTEGER      ,INTENT(IN   ) :: k1 indices for pressure levels (k1 and k2)
                        k2                    , &!INTEGER      ,INTENT(IN   ) :: k2 indices for pressure levels (k1 and k2)
                        m                     , &!INTEGER      ,INTENT(IN   ) :: m  number of grid intervals in zonal direction (m)
                        np                    , &!INTEGER      ,INTENT(IN   ) :: np number of atmospheric layers (np)
                        nx                    , &!INTEGER      ,INTENT(IN   ) :: nx number of pressure intervals in the table (nx)
                        nh                    , &!INTEGER      ,INTENT(IN   ) :: nh number of absorber amount intervals in the table (nh)
                        nt                    , &!INTEGER      ,INTENT(IN   ) :: nt number of tables copied (nt)
                        sh2o  (1:m,1:np+1)    , &!REAL(KIND=r8),INTENT(IN   ) :: sabs(m,n,np+1) column-integrated absorber amount (sabs)
                        swpre (1:m,1:np+1)    , &!REAL(KIND=r8),INTENT(IN   ) :: spre(m,n,np+1) column absorber amount-weighted pressure (spre)
                        swtem (1:m,1:np+1)    , &!REAL(KIND=r8),INTENT(IN   ) :: stem(m,n,np+1) column absorber amount-weighted temperature (stem)
                        w1                    , &!REAL(KIND=r8),INTENT(IN   ) :: w1             first value of absorber amount (log10) in the table (w1)
                        p1                    , &!REAL(KIND=r8),INTENT(IN   ) :: p1             first value of pressure (log10) in the table (p1)
                        dwe                   , &!REAL(KIND=r8),INTENT(IN   ) :: dwe            size of the interval of absorber amount (log10) in the table (dwe)
                        dpe                   , &!REAL(KIND=r8),INTENT(IN   ) :: dpe            size of the interval of pressure (log10) in the table (dpe)
                        h21   (1:nx,1:nh,1:nt), &!REAL(KIND=r8),INTENT(IN   ) :: coef1(nx,nh,nt) pre-computed coefficients (coef1, coef2, and coef3)
                        h22   (1:nx,1:nh,1:nt), &!REAL(KIND=r8),INTENT(IN   ) :: coef2(nx,nh,nt) pre-computed coefficients (coef1, coef2, and coef3)
                        h23   (1:nx,1:nh,1:nt), &!REAL(KIND=r8),INTENT(IN   ) :: coef3(nx,nh,nt) pre-computed coefficients (coef1, coef2, and coef3)
                        trant (1:m)             )!REAL(KIND=r8),INTENT(INOUT) :: tran(m,n)  updated transmittance (tran)

                END IF
                IF (ib == 7 ) THEN

                   CALL tablup(k1                    , &!INTEGER      ,INTENT(IN   ) :: k1 indices for pressure levels (k1 and k2)
                        k2                    , &!INTEGER      ,INTENT(IN   ) :: k2 indices for pressure levels (k1 and k2)
                        m                     , &!INTEGER      ,INTENT(IN   ) :: m  number of grid intervals in zonal direction (m)
                        np                    , &!INTEGER      ,INTENT(IN   ) :: np number of atmospheric layers (np)
                        nx                    , &!INTEGER      ,INTENT(IN   ) :: nx number of pressure intervals in the table (nx)
                        nh                    , &!INTEGER      ,INTENT(IN   ) :: nh number of absorber amount intervals in the table (nh)
                        nt                    , &!INTEGER      ,INTENT(IN   ) :: nt number of tables copied (nt)
                        sh2o  (1:m,1:np+1)    , &!REAL(KIND=r8),INTENT(IN   ) :: sabs(m,n,np+1) column-integrated absorber amount (sabs)
                        swpre (1:m,1:np+1)    , &!REAL(KIND=r8),INTENT(IN   ) :: spre(m,n,np+1) column absorber amount-weighted pressure (spre)
                        swtem (1:m,1:np+1)    , &!REAL(KIND=r8),INTENT(IN   ) :: stem(m,n,np+1) column absorber amount-weighted temperature (stem)
                        w1                    , &!REAL(KIND=r8),INTENT(IN   ) :: w1                first value of absorber amount (log10) in the table (w1)
                        p1                    , &!REAL(KIND=r8),INTENT(IN   ) :: p1                first value of pressure (log10) in the table (p1)
                        dwe                   , &!REAL(KIND=r8),INTENT(IN   ) :: dwe                size of the interval of absorber amount (log10) in the table (dwe)
                        dpe                   , &!REAL(KIND=r8),INTENT(IN   ) :: dpe                size of the interval of pressure (log10) in the table (dpe)
                        h71   (1:nx,1:nh,1:nt), &!REAL(KIND=r8),INTENT(IN   ) :: coef1(nx,nh,nt) pre-computed coefficients (coef1, coef2, and coef3)
                        h72   (1:nx,1:nh,1:nt), &!REAL(KIND=r8),INTENT(IN   ) :: coef2(nx,nh,nt) pre-computed coefficients (coef1, coef2, and coef3)
                        h73   (1:nx,1:nh,1:nt), &!REAL(KIND=r8),INTENT(IN   ) :: coef3(nx,nh,nt) pre-computed coefficients (coef1, coef2, and coef3)
                        trant (1:m)             )!REAL(KIND=r8),INTENT(INOUT) :: tran(m,n)  updated transmittance (tran)


                END IF

             ELSE

                !-----compute water vapor transmittance using k-distribution.

                CALL wvkdis(ib                        ,& !INTEGER      , INTENT(IN   ) :: ib !  spectral band (ib)
                     m                         ,& !INTEGER      , INTENT(IN   ) :: m  !  number of grid intervals in zonal direction (m)
                     np                        ,& !INTEGER      , INTENT(IN   ) :: np !  number of levels (np)
                     k2-1                      ,& !INTEGER      , INTENT(IN   ) :: k  !  current level (k)
                     h2oexp(1:m,1:np,1:6)      ,& !REAL(KIND=r8), INTENT(IN   ) :: h2oexp(m,n,np,6) !  exponentials for continuum absorption (conexp)
                     conexp(1:m,1:np,1:3)      ,& !REAL(KIND=r8), INTENT(IN   ) :: conexp(m,n,np,3) !  exponentials for line absorption (h2oexp)
                     th2o  (1:m,1:6)           ,& !REAL(KIND=r8), INTENT(INOUT) :: th2o  (m,n,6)      !  transmittance between levels k1 and k2 due to water vapor line absorption (th2o)
                     tcon  (1:m,1:3)           ,& !REAL(KIND=r8), INTENT(INOUT) :: tcon  (m,n,3)      !  transmittance between levels k1 and k2 due to water vapor continuum absorption (tcon)
                     trant (1:m)                ) !REAL(KIND=r8), INTENT(INOUT) :: tran  (m,n)        !  total transmittance (tran)

             END IF

             IF(co2bnd) THEN

                IF( high ) THEN

                   !-----compute co2 transmittance using table look-up method

                   w1=-4.0_r8
                   p1=-2.0_r8
                   dwe=0.3_r8
                   dpe=0.2_r8

                   CALL tablup(k1                     , &!INTEGER      ,INTENT(IN   ) :: k1 indices for pressure levels (k1 and k2)
                        k2                     , &!INTEGER      ,INTENT(IN   ) :: k2 indices for pressure levels (k1 and k2)
                        m                      , &!INTEGER      ,INTENT(IN   ) :: m  number of grid intervals in zonal direction (m)
                        np                     , &!INTEGER      ,INTENT(IN   ) :: np number of atmospheric layers (np)
                        nx                     , &!INTEGER      ,INTENT(IN   ) :: nx number of pressure intervals in the table (nx)
                        nc                     , &!INTEGER      ,INTENT(IN   ) :: nh number of absorber amount intervals in the table (nh)
                        nt                     , &!INTEGER      ,INTENT(IN   ) :: nt number of tables copied (nt)
                        sco3   (1:m,1:np+1)    , &!REAL(KIND=r8),INTENT(IN   ) :: sabs(m,n,np+1) column-integrated absorber amount (sabs)
                        scopre (1:m,1:np+1)    , &!REAL(KIND=r8),INTENT(IN   ) :: spre(m,n,np+1) column absorber amount-weighted pressure (spre)
                        scotem (1:m,1:np+1)    , &!REAL(KIND=r8),INTENT(IN   ) :: stem(m,n,np+1) column absorber amount-weighted temperature (stem)
                        w1                     , &!REAL(KIND=r8),INTENT(IN   ) :: w1                first value of absorber amount (log10) in the table (w1)
                        p1                     , &!REAL(KIND=r8),INTENT(IN   ) :: p1                first value of pressure (log10) in the table (p1)
                        dwe                    , &!REAL(KIND=r8),INTENT(IN   ) :: dwe                 size of the interval of absorber amount (log10) in the table (dwe)
                        dpe                    , &!REAL(KIND=r8),INTENT(IN   ) :: dpe                 size of the interval of pressure (log10) in the table (dpe)
                        c1     (1:nx,1:nc,1:nt), &!REAL(KIND=r8),INTENT(IN   ) :: coef1(nx,nh,nt) pre-computed coefficients (coef1, coef2, and coef3)
                        c2     (1:nx,1:nc,1:nt), &!REAL(KIND=r8),INTENT(IN   ) :: coef2(nx,nh,nt) pre-computed coefficients (coef1, coef2, and coef3)
                        c3     (1:nx,1:nc,1:nt), &!REAL(KIND=r8),INTENT(IN   ) :: coef3(nx,nh,nt) pre-computed coefficients (coef1, coef2, and coef3)
                        trant  (1:m)       )      !REAL(KIND=r8),INTENT(INOUT) :: tran(m,n)  updated transmittance (tran)

                ELSE

                   !-----compute co2 transmittance using k-distribution method

                   CALL co2kdis(m                 , &! INTEGER          , INTENT(IN        ) :: m  !         number of grid intervals in zonal direction (m)
                        np                        , &! INTEGER          , INTENT(IN        ) :: np
                        k2-1                      , &! INTEGER          , INTENT(IN        ) :: k
                        co2exp(1:m,1:np,1:6,1:2)  , &! REAL(KIND=r8)    , INTENT(IN        ) :: co2exp(m,n,np,6,2)
                        tco2  (1:m,1:6,1:2)       , &! REAL(KIND=r8)    , INTENT(INOUT     ) :: tco2  (m,n,6,2)    !  for the various values of the absorption coefficient (tco2)
                        trant (1:m)                 )! REAL(KIND=r8)    , INTENT(INOUT     ) :: tran  (m,n)        !    total transmittance (tran)

                END IF

             END IF

             !-----compute o3 transmittance using table look-up

             IF (oznbnd) THEN

                w1=-6.0_r8
                p1=-2.0_r8
                dwe=0.3_r8
                dpe=0.2_r8

                CALL tablup(k1                      , &!INTEGER      ,INTENT(IN   ) :: k1 indices for pressure levels (k1 and k2)
                     k2                      , &!INTEGER      ,INTENT(IN   ) :: k2 indices for pressure levels (k1 and k2)
                     m                       , &!INTEGER      ,INTENT(IN   ) :: m  number of grid intervals in zonal direction (m)
                     np                      , &!INTEGER      ,INTENT(IN   ) :: np number of atmospheric layers (np)
                     nx                      , &!INTEGER      ,INTENT(IN   ) :: nx number of pressure intervals in the table (nx)
                     no                      , &!INTEGER      ,INTENT(IN   ) :: nh number of absorber amount intervals in the table (nh)
                     nt                      , &!INTEGER      ,INTENT(IN   ) :: nt number of tables copied (nt)
                     sco3    (1:m,1:np+1), &!REAL(KIND=r8),INTENT(IN   ) :: sabs(m,n,np+1) column-integrated absorber amount (sabs)
                     scopre  (1:m,1:np+1), &!REAL(KIND=r8),INTENT(IN   ) :: spre(m,n,np+1) column absorber amount-weighted pressure (spre)
                     scotem  (1:m,1:np+1), &!REAL(KIND=r8),INTENT(IN   ) :: stem(m,n,np+1) column absorber amount-weighted temperature (stem)
                     w1                      , &!REAL(KIND=r8),INTENT(IN   ) :: w1               first value of absorber amount (log10) in the table (w1)
                     p1                      , &!REAL(KIND=r8),INTENT(IN   ) :: p1               first value of pressure (log10) in the table (p1)
                     dwe                     , &!REAL(KIND=r8),INTENT(IN   ) :: dwe                 size of the interval of absorber amount (log10) in the table (dwe)
                     dpe                     , &!REAL(KIND=r8),INTENT(IN   ) :: dpe                 size of the interval of pressure (log10) in the table (dpe)
                     o1      (1:nx,1:no,1:nt), &!REAL(KIND=r8),INTENT(IN   ) :: coef1(nx,nh,nt) pre-computed coefficients (coef1, coef2, and coef3)
                     o2      (1:nx,1:no,1:nt), &!REAL(KIND=r8),INTENT(IN   ) :: coef2(nx,nh,nt) pre-computed coefficients (coef1, coef2, and coef3)
                     o3      (1:nx,1:no,1:nt), &!REAL(KIND=r8),INTENT(IN   ) :: coef3(nx,nh,nt) pre-computed coefficients (coef1, coef2, and coef3)
                     trant   (1:m)         )!REAL(KIND=r8),INTENT(INOUT) :: tran(m,n)  updated transmittance (tran)

             END IF

             !-----fclr is the clear line-of-sight between levels k1 and k2.
             !  in computing fclr, clouds are assumed randomly overlapped
             !  using eq. (10).

             DO i=1,m
                fclr(i) = fclr(i)*clr(i,k2-1)
             END DO

             !-----compute upward and downward fluxes


             !-----add "boundary" terms to the net downward flux.
             !  these are the first terms on the right-hand-side of
             !  eqs. (56a) and (56b).
             !  downward fluxes are positive.

             IF (k2 == k1+1) THEN
                DO i=1,m
                   flc(i,k1)=flc(i,k1)-blayer(i,k1)
                   flc(i,k2)=flc(i,k2)+blayer(i,k1)
                END DO
             END IF

             !-----add flux components involving the four layers above and below
             !  the levels k1 and k2.  it follows eqs. (56a) and (56b).

             DO i=1,m
                xx=trant(i)*(blayer(i,k2-1)-blayer(i,k2))
                flc(i,k1) =flc(i,k1)+xx
                xx=trant(i)*(blayer(i,k1-1)-blayer(i,k1))
                flc(i,k2) =flc(i,k2)+xx
             END DO

             !-----compute upward and downward fluxes for all-sky situation

             IF (k2 == k1+1) THEN
                DO i=1,m
                   flxu(i,k1)=flxu(i,k1)-blayer(i,k1)
                   flxd(i,k2)=flxd(i,k2)+blayer(i,k1)
                END DO
             END IF

             DO i=1,m
                xx=trant(i)*(blayer(i,k2-1)-blayer(i,k2))
                flxu(i,k1) =flxu(i,k1)+xx*fclr(i)
                xx=trant(i)*(blayer(i,k1-1)-blayer(i,k1))
                flxd(i,k2) =flxd(i,k2)+xx*fclr(i)
             END DO


          END DO

          !-----compute the partial derivative of fluxes with respect to
          !  surface temperature (eq. 59).

          DO i=1,m
             dfdts(i,k1) =dfdts(i,k1)-dbs(i)*trant(i)*fclr(i)
          END DO

       END DO

       !-----add contribution from the surface to the flux terms at the surface.

       DO i=1,m
          dfdts(i,np+1) =dfdts(i,np+1)-dbs(i)
       END DO

       DO i=1,m
          flc(i,np+1)=flc(i,np+1)-blayer(i,np+1)
          flxu(i,np+1)=flxu(i,np+1)-blayer(i,np+1)
          st4(i)=st4(i)-blayer(i,np+1)
       END DO


       !  write(7,3211) ib, flxd(1,1,52),flxu(1,1,52)
       !  write(7,3211) ib, flxd(1,1,np+1),flxu(1,1,np+1)
       !    3211 FORMAT ('ib, fluxd, fluxu=', i3,2F12.3)

    END DO

    DO k=1,np+1
       DO i=1,m
          flx(i,k)   = flxd(i,k)+flxu(i,k)
          radlwin(i) = flxd(i,np+1) 
       END DO
    END DO

    RETURN
  END SUBROUTINE irrad



  !
  !
  !##################################################################
  !##################################################################
  !######                                                      ######
  !######                SUBROUTINE COLUMN                     ######
  !######                                                      ######
  !######                     Developed by                     ######
  !######                                                      ######
  !######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
  !######                                                      ######
  !######     Center for Analysis and Prediction of Storms     ######
  !######               University of Oklahoma                 ######
  !######                                                      ######
  !##################################################################
  !##################################################################
  !


  SUBROUTINE column(m     , &! INTEGER      , INTENT(IN   ) :: m !   number of soundings in zonal direction (m)
       np    , &! INTEGER      , INTENT(IN   ) :: np !   number of atmospheric layers (np)
       pa    , &! REAL(KIND=r8), INTENT(IN   ) :: pa   (m,n,np) !   layer pressure (pa)
       dt    , &! REAL(KIND=r8), INTENT(IN   ) :: dt   (m,n,np) !   layer temperature minus 250K (dt)
       sabs0 , &! REAL(KIND=r8), INTENT(IN   ) :: sabs0(m,n,np)!   layer absorber amount (sabs0)
       sabs  , &! REAL(KIND=r8), INTENT(OUT  ) :: sabs (m,n,np+1)!   column-integrated absorber amount (sabs)
       spre  , &! REAL(KIND=r8), INTENT(OUT  ) :: spre (m,n,np+1)!   column absorber-weighted pressure (spre)
       stem    )! REAL(KIND=r8), INTENT(OUT  ) :: stem (m,n,np+1)!   column absorber-weighted temperature (stem)
    !
    !-----------------------------------------------------------------------
    !
    !  PURPOSE:
    !
    !  Calculate column-integrated (from top of the model atmosphere)
    !  absorber amount, absorber-weighted pressure and temperature.
    !
    !-----------------------------------------------------------------------
    !
    !  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
    !          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
    !              Chou (1996)
    !
    !  MODIFICATION HISTORY:
    !
    !  03/15/1996 (Yuhe Liu)
    !  Adopted the original code and formatted it in accordance with the
    !  ARPS coding standard.
    !
    !-----------------------------------------------------------------------
    !
    !  ORIGINAL COMMENTS:
    !
    !**************************************************************************
    !-----compute column-integrated (from top of the model atmosphere)
    !  absorber amount (sabs), absorber-weighted pressure (spre) and
    !  temperature (stem).
    !  computations of spre and stem follows eqs. (37) and (38).
    !
    !--- input parameters
    !   number of soundings in zonal direction (m)
    !   number of soundings in meridional direction (n)
    !   number of atmospheric layers (np)
    !   layer pressure (pa)
    !   layer temperature minus 250K (dt)
    !   layer absorber amount (sabs0)
    !
    !--- output parameters
    !   column-integrated absorber amount (sabs)
    !   column absorber-weighted pressure (spre)
    !   column absorber-weighted temperature (stem)
    !
    !--- units of pa and dt are mb and k, respectively.
    !    units of sabs are g/cm**2 for water vapor and (cm-atm)stp for co2 and o3
    !**************************************************************************
    !
    !-----------------------------------------------------------------------
    !
    !  Variable Declarations.
    !
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE

    INTEGER      , INTENT(IN   ) :: m
    INTEGER      , INTENT(IN   ) :: np

    !---- input parameters -----

    REAL(KIND=r8), INTENT(IN   ) :: pa   (m,np)
    REAL(KIND=r8), INTENT(IN   ) :: dt   (m,np)
    REAL(KIND=r8), INTENT(IN   ) :: sabs0(m,np)

    !---- output parameters -----

    REAL(KIND=r8), INTENT(OUT  ) :: sabs (m,np+1)
    REAL(KIND=r8), INTENT(OUT  ) :: spre (m,np+1)
    REAL(KIND=r8), INTENT(OUT  ) :: stem (m,np+1)


    !*********************************************************************
    INTEGER :: i
    INTEGER :: k

    !
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !
    !  Beginning of executable code...
    !
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !
    DO i=1,m
       sabs(i,1)=0.0_r8
       spre(i,1)=0.0_r8
       stem(i,1)=0.0_r8
    END DO

    DO k=1,np
       DO i=1,m
          sabs(i,k+1)=sabs(i,k)+sabs0(i,k)
          spre(i,k+1)=spre(i,k)+pa(i,k)*sabs0(i,k)
          stem(i,k+1)=stem(i,k)+dt(i,k)*sabs0(i,k)
       END DO
    END DO

    RETURN
  END SUBROUTINE column

  !
  !
  !##################################################################
  !##################################################################
  !######                                                      ######
  !######                SUBROUTINE H2OEXPS                    ######
  !######                                                      ######
  !######                     Developed by                     ######
  !######                                                      ######
  !######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
  !######                                                      ######
  !######     Center for Analysis and Prediction of Storms     ######
  !######               University of Oklahoma                 ######
  !######                                                      ######
  !##################################################################
  !##################################################################
  !


  SUBROUTINE h2oexps(ib      , &! INTEGER, INTENT(IN   ) :: ib!  spectral band (ib)
       m      , &! INTEGER, INTENT(IN   ) :: m !  number of grid intervals in zonal direction (m)
       np     , &! INTEGER, INTENT(IN   ) :: np!  number of layers (np)
       dh2o   , &! REAL(KIND=r8), INTENT(IN   ) :: dh2o(m,n,np)!  layer water vapor amount for line absorption (dh2o)
       pa     , &! REAL(KIND=r8), INTENT(IN   ) :: pa  (m,n,np)!  layer pressure (pa)
       dt     , &! REAL(KIND=r8), INTENT(IN   ) :: dt  (m,n,np)!  layer temperature minus 250K (dt)
       h2oexp   )! REAL(KIND=r8), INTENT(INOUT) :: h2oexp(m,n,np,6) !  6 exponentials for each layer  (h2oexp)
    !
    !-----------------------------------------------------------------------
    !
    !  PURPOSE:
    !
    !  Calculate exponentials for water vapor line absorption in
    !  individual layers.
    !
    !-----------------------------------------------------------------------
    !
    !
    !  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
    !          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
    !              Chou (1996)
    !
    !  MODIFICATION:
    !
    !  03/11/1996 (Yuhe Liu)
    !  Formatted code to ARPS standard format
    !
    !-----------------------------------------------------------------------
    !
    !  ORIGINAL COMMENTS:
    !
    !   compute exponentials for water vapor line absorption
    !   in individual layers.
    !
    !---- input parameters
    !  spectral band (ib)
    !  number of grid intervals in zonal direction (m)
    !  number of grid intervals in meridional direction (n)
    !  number of layers (np)
    !  layer water vapor amount for line absorption (dh2o)
    !  layer pressure (pa)
    !  layer temperature minus 250K (dt)
    !
    !---- output parameters
    !  6 exponentials for each layer  (h2oexp)
    !
    !**********************************************************************
    !
    !fpp$ expand (expmn)
    !!dir$ inline always expmn
    !*$*  inline routine (expmn)
    !
    !-----------------------------------------------------------------------
    !
    !  Variable declarations
    !
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: ib
    INTEGER, INTENT(IN   ) :: m
    INTEGER, INTENT(IN   ) :: np

    !---- input parameters ------

    REAL(KIND=r8), INTENT(IN   ) :: dh2o(m,np)
    REAL(KIND=r8), INTENT(IN   ) :: pa  (m,np)
    REAL(KIND=r8), INTENT(IN   ) :: dt  (m,np)

    !---- output parameters -----

    REAL(KIND=r8), INTENT(INOUT) :: h2oexp(m,np,6)

    !---- temporary arrays -----

    REAL(KIND=r8) :: xh

    !---- local misc. variables

    INTEGER :: i,k,ik

    !---- static data -----

    !-----xkw  are the absorption coefficients for the first
    !  k-distribution function due to water vapor line absorption
    !  (tables 4 and 7).  units are cm**2/g

    REAL(KIND=r8), PARAMETER :: xkw(8) =RESHAPE(SOURCE=(/&
         29.55_r8, 4.167E-1_r8, 1.328E-2_r8, 5.250E-4_r8, &
         5.25E-4_r8, 2.340E-3_r8, 1.320E-0_r8, 5.250E-4_r8/),SHAPE=(/8/))

    !-----mw are the ratios between neighboring absorption coefficients
    !  for water vapor line absorption (tables 4 and 7).

    INTEGER, PARAMETER :: mw(8) =RESHAPE(SOURCE=(/6,6,8,6,6,8,6,16/),SHAPE=(/8/))

    !-----aw and bw (table 3) are the coefficients for temperature scaling
    !  in eq. (25).

    REAL(KIND=r8), PARAMETER :: aw(8) =RESHAPE(SOURCE=(/&
         0.0021_r8, 0.0140_r8, 0.0167_r8, 0.0302_r8, &
         0.0307_r8, 0.0154_r8, 0.0008_r8, 0.0096_r8/),SHAPE=(/8/))

    REAL(KIND=r8), PARAMETER :: bw(8) =RESHAPE(SOURCE=(/&
         -1.01E-5_r8, 5.57E-5_r8, 8.54E-5_r8, 2.96E-4_r8,&
         2.86E-4_r8, 7.53E-5_r8,-3.52E-6_r8, 1.64E-5_r8/),SHAPE=(/8/))

    !-----expmn is an external function

    !  REAL :: expmn
    !
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !
    !  Beginning of executable code...
    !
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !

    !**********************************************************************
    !    note that the 3 sub-bands in band 3 use the same set of xkw, aw,
    !    and bw.  therefore, h2oexp for these sub-bands are identical.
    !**********************************************************************

    DO k=1,np
       DO i=1,m

          !-----xh is   the scaled water vapor amount for line absorption
          !  computed from (27).

          xh = dh2o(i,k)*(pa(i,k)*0.002_r8)                              &
               * ( 1.0_r8+(aw(ib)+bw(ib)* dt(i,k))*dt(i,k) )

          !-----h2oexp is the water vapor transmittance of the layer (k2-1)
          !  due to line absorption

          h2oexp(i,k,1) = expmn(-xh*xkw(ib))

       END DO
    END DO

    DO ik=2,6

       IF(mw(ib) == 6) THEN

          DO k=1,np
             DO i=1,m
                xh = h2oexp(i,k,ik-1)*h2oexp(i,k,ik-1)
                h2oexp(i,k,ik) = xh*xh*xh
             END DO
          END DO

       ELSE IF(mw(ib) == 8) THEN

          DO k=1,np
             DO i=1,m
                xh = h2oexp(i,k,ik-1)*h2oexp(i,k,ik-1)
                xh = xh*xh
                h2oexp(i,k,ik) = xh*xh
             END DO
          END DO

       ELSE

          DO k=1,np
             DO i=1,m
                xh = h2oexp(i,k,ik-1)*h2oexp(i,k,ik-1)
                xh = xh*xh
                xh = xh*xh
                h2oexp(i,k,ik) = xh*xh
             END DO
          END DO

       END IF
    END DO

    RETURN
  END SUBROUTINE h2oexps
  !
  !
  !##################################################################
  !##################################################################
  !######                                                      ######
  !######                SUBROUTINE CONEXPS                    ######
  !######                                                      ######
  !######                     Developed by                     ######
  !######                                                      ######
  !######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
  !######                                                      ######
  !######     Center for Analysis and Prediction of Storms     ######
  !######               University of Oklahoma                 ######
  !######                                                      ######
  !##################################################################
  !##################################################################
  !


  SUBROUTINE conexps( ib     ,&!INTEGER, INTENT(IN   ) :: ib!  spectral band (ib)
       m      ,&!INTEGER, INTENT(IN   ) :: m  !  number of grid intervals in zonal direction (m)
       np     ,&!INTEGER, INTENT(IN   ) :: np !  number of layers (np)
       dcont  ,&!REAL(KIND=r8), INTENT(IN   ) :: dcont(m,n,np) !  layer scaled water vapor amount for continuum absorption (dcont)
       conexp  )!REAL(KIND=r8), INTENT(INOUT) :: conexp(m,n,np,3)!  1 or 3 exponentials for each layer (conexp)
    !
    !-----------------------------------------------------------------------
    !
    !  PURPOSE:
    !
    !  Calculate exponentials for continuum absorption in individual
    !  layers.
    !
    !-----------------------------------------------------------------------
    !
    !  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
    !          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
    !              Chou (1996)
    !
    !  MODIFICATION HISTORY:
    !
    !  03/15/1996 (Yuhe Liu)
    !  Adopted the original code and formatted it in accordance with the
    !  ARPS coding standard.
    !
    !-----------------------------------------------------------------------
    !
    !  ORIGINAL COMMENTS:
    !
    !**********************************************************************
    !   compute exponentials for continuum absorption in individual layers.
    !
    !---- input parameters
    !  spectral band (ib)
    !  number of grid intervals in zonal direction (m)
    !  number of grid intervals in meridional direction (n)
    !  number of layers (np)
    !  layer scaled water vapor amount for continuum absorption (dcont)
    !
    !---- output parameters
    !  1 or 3 exponentials for each layer (conexp)
    !
    !**********************************************************************
    !
    !fpp$ expand (expmn)
    !!dir$ inline always expmn
    !*$*  inline routine (expmn)
    !
    !-----------------------------------------------------------------------
    !
    !  Variable Declarations.
    !
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: ib
    INTEGER, INTENT(IN   ) :: m
    INTEGER, INTENT(IN   ) :: np

    !---- input parameters ------

    REAL(KIND=r8), INTENT(IN   ) :: dcont(m,np)

    !---- updated parameters -----

    REAL(KIND=r8), INTENT(INOUT) :: conexp(m,np,3)

    !---- temporary arrays -----

    INTEGER :: i
    INTEGER :: k
    INTEGER :: iq


    !---- static data -----



    !-----xke are the absorption coefficients for the first
    !  k-distribution function due to water vapor continuum absorption
    !  (table 6).  units are cm**2/g

    REAL(KIND=r8), PARAMETER :: xke(8) =RESHAPE(SOURCE=(/&
         0.00_r8,   0.00_r8,   27.40_r8,   15.8_r8, &
         9.40_r8,   7.75_r8,     0.0_r8,    0.0_r8/),SHAPE=(/8/))


    !-----ne is the number of terms in computing water vapor
    !  continuum transmittance (Table 6).
    !  band 3 is divided into 3 sub-bands.

    INTEGER, PARAMETER :: NE(8) =RESHAPE(SOURCE=(/0,0,3,1,1,1,0,0/),SHAPE=(/8/))
    !
    !-----------------------------------------------------------------------
    !
    !  Functions:
    !
    !-----------------------------------------------------------------------
    !

    !-----expmn is an external function
    !  REAL :: expmn
    !
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !
    !  Beginning of executable code...
    !
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !
    DO k=1,np
       DO i=1,m
          conexp(i,k,1) = expmn(-dcont(i,k)*xke(ib))
       END DO
    END DO

    IF (ib == 3) THEN

       !-----the absorption coefficients for sub-bands 3b (iq=2) and 3a (iq=3)
       !  are, respectively, double and quadruple that for sub-band 3c (iq=1)
       !  (table 6).

       DO iq=2,3
          DO k=1,np
             DO i=1,m
                conexp(i,k,iq) = conexp(i,k,iq-1) *conexp(i,k,iq-1)
             END DO
          END DO
       END DO

    END IF

    RETURN
  END SUBROUTINE conexps

  !
  !
  !##################################################################
  !##################################################################
  !######                                                      ######
  !######                SUBROUTINE CO2EXPS                    ######
  !######                                                      ######
  !######                     Developed by                     ######
  !######                                                      ######
  !######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
  !######                                                      ######
  !######     Center for Analysis and Prediction of Storms     ######
  !######               University of Oklahoma                 ######
  !######                                                      ######
  !##################################################################
  !##################################################################
  !


  SUBROUTINE co2exps(m      ,& !INTEGER, INTENT(IN   ) :: m !  number of grid intervals in zonal direction (m)
       np     ,& !INTEGER, INTENT(IN   ) :: np!  number of layers (np)
       dco2   ,& !REAL(KIND=r8), INTENT(IN   ) :: dco2(m,n,np)!  layer co2 amount (dco2)
       pa     ,& !REAL(KIND=r8), INTENT(IN   ) :: pa  (m,n,np)!  layer pressure (pa)
       dt     ,& !REAL(KIND=r8), INTENT(IN   ) :: dt  (m,n,np)!  layer temperature minus 250K (dt)
       co2exp )  !REAL(KIND=r8), INTENT(INOUT) :: co2exp(m,n,np,6,2) !  6 exponentials for each layer (co2exp)
    !
    !-----------------------------------------------------------------------
    !
    !  PURPOSE:
    !
    !  Calculate co2 exponentials for individual layers.
    !
    !-----------------------------------------------------------------------
    !
    !  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
    !          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
    !              Chou (1996)
    !
    !  MODIFICATION HISTORY:
    !
    !  03/15/1996 (Yuhe Liu)
    !  Adopted the original code and formatted it in accordance with the
    !  ARPS coding standard.
    !
    !-----------------------------------------------------------------------
    !
    !  ORIGINAL COMMENTS:
    !
    !
    !**********************************************************************
    !   compute co2 exponentials for individual layers.
    !
    !---- input parameters
    !  number of grid intervals in zonal direction (m)
    !  number of grid intervals in meridional direction (n)
    !  number of layers (np)
    !  layer co2 amount (dco2)
    !  layer pressure (pa)
    !  layer temperature minus 250K (dt)
    !
    !---- output parameters
    !  6 exponentials for each layer (co2exp)
    !**********************************************************************
    !
    !fpp$ expand (expmn)
    !!dir$ inline always expmn
    !*$*  inline routine (expmn)
    !
    !-----------------------------------------------------------------------
    !
    !  Variable Declarations.
    !
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: m
    INTEGER, INTENT(IN   ) :: np

    !---- input parameters -----

    REAL(KIND=r8), INTENT(IN   ) :: dco2(m,np)
    REAL(KIND=r8), INTENT(IN   ) :: pa  (m,np)
    REAL(KIND=r8), INTENT(IN   ) :: dt  (m,np)

    !---- output parameters -----

    REAL(KIND=r8), INTENT(INOUT) :: co2exp(m,np,6,2)

    !---- temporary arrays -----

    REAL(KIND=r8) :: xc
    INTEGER :: i
    INTEGER :: k

    !---- static data -----

    !-----xkc is the absorption coefficients for the
    !  first k-distribution function due to co2 (table 7).
    !  units are 1/(cm-atm)stp.

    REAL(KIND=r8), PARAMETER :: xkc(2) =RESHAPE(SOURCE=(/2.656E-5_r8,2.656E-3_r8/),SHAPE=(/2/))

    !-----parameters (table 3) for computing the scaled co2 amount
    !  using (27).

    REAL(KIND=r8), PARAMETER :: prc(2) =RESHAPE(SOURCE=(/  300.0_r8,   30.0_r8/),SHAPE=(/2/))
    REAL(KIND=r8), PARAMETER :: pm(2)  =RESHAPE(SOURCE=(/    0.5_r8,   0.85_r8/),SHAPE=(/2/))
    REAL(KIND=r8), PARAMETER :: ac(2)  =RESHAPE(SOURCE=(/ 0.0182_r8, 0.0042_r8/),SHAPE=(/2/))
    REAL(KIND=r8), PARAMETER :: bc(2)  =RESHAPE(SOURCE=(/1.07E-4_r8,2.00E-5_r8/),SHAPE=(/2/))
    !
    !-----------------------------------------------------------------------
    !
    !  Functions:
    !
    !-----------------------------------------------------------------------
    !

    !-----function expmn

    !  REAL :: expmn
    !
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !
    !  Beginning of executable code...
    !
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !
    DO k=1,np
       DO i=1,m

          !-----compute the scaled co2 amount from eq. (27) for band-wings
          !  (sub-bands 3a and 3c).

          xc = dco2(i,k)*(pa(i,k)/prc(1))**pm(1)                      &
               *(1.0_r8+(ac(1)+bc(1)*dt(i,k))*dt(i,k))

          !-----six exponential by powers of 8 (table 7).

          co2exp(i,k,1,1)=expmn(-xc*xkc(1))

          xc=co2exp(i,k,1,1)*co2exp(i,k,1,1)
          xc=xc*xc
          co2exp(i,k,2,1)=xc*xc

          xc=co2exp(i,k,2,1)*co2exp(i,k,2,1)
          xc=xc*xc
          co2exp(i,k,3,1)=xc*xc

          xc=co2exp(i,k,3,1)*co2exp(i,k,3,1)
          xc=xc*xc
          co2exp(i,k,4,1)=xc*xc

          xc=co2exp(i,k,4,1)*co2exp(i,k,4,1)
          xc=xc*xc
          co2exp(i,k,5,1)=xc*xc

          xc=co2exp(i,k,5,1)*co2exp(i,k,5,1)
          xc=xc*xc
          co2exp(i,k,6,1)=xc*xc

          !-----compute the scaled co2 amount from eq. (27) for band-center
          !  region (sub-band 3b).

          xc = dco2(i,k)*(pa(i,k)/prc(2))**pm(2)                      &
               *(1.0_r8+(ac(2)+bc(2)*dt(i,k))*dt(i,k))

          co2exp(i,k,1,2)=expmn(-xc*xkc(2))

          xc=co2exp(i,k,1,2)*co2exp(i,k,1,2)
          xc=xc*xc
          co2exp(i,k,2,2)=xc*xc

          xc=co2exp(i,k,2,2)*co2exp(i,k,2,2)
          xc=xc*xc
          co2exp(i,k,3,2)=xc*xc

          xc=co2exp(i,k,3,2)*co2exp(i,k,3,2)
          xc=xc*xc
          co2exp(i,k,4,2)=xc*xc

          xc=co2exp(i,k,4,2)*co2exp(i,k,4,2)
          xc=xc*xc
          co2exp(i,k,5,2)=xc*xc

          xc=co2exp(i,k,5,2)*co2exp(i,k,5,2)
          xc=xc*xc
          co2exp(i,k,6,2)=xc*xc

       END DO
    END DO

    RETURN
  END SUBROUTINE co2exps


  !
  !
  !##################################################################
  !##################################################################
  !######                                                      ######
  !######                SUBROUTINE WVKDIS                     ######
  !######                                                      ######
  !######                     Developed by                     ######
  !######                                                      ######
  !######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
  !######                                                      ######
  !######     Center for Analysis and Prediction of Storms     ######
  !######               University of Oklahoma                 ######
  !######                                                      ######
  !##################################################################
  !##################################################################
  !


  SUBROUTINE wvkdis(ib     ,& !INTEGER      , INTENT(IN   ) :: ib !  spectral band (ib)
       m      ,& !INTEGER      , INTENT(IN   ) :: m  !  number of grid intervals in zonal direction (m)
       np     ,& !INTEGER      , INTENT(IN   ) :: np !  number of levels (np)
       k      ,& !INTEGER      , INTENT(IN   ) :: k  !  current level (k)
       h2oexp ,& !REAL(KIND=r8), INTENT(IN   ) :: conexp(m,n,np,3) !  exponentials for line absorption (h2oexp)
       conexp ,& !REAL(KIND=r8), INTENT(IN   ) :: h2oexp(m,n,np,6) !  exponentials for continuum absorption (conexp)
       th2o   ,& !REAL(KIND=r8), INTENT(INOUT) :: th2o(m,n,6) !  transmittance between levels k1 and k2 due to water vapor line absorption (th2o)
       tcon   ,& !REAL(KIND=r8), INTENT(INOUT) :: tcon(m,n,3) !  transmittance between levels k1 and k2 due to water vapor continuum absorption (tcon)
       tran    ) !REAL(KIND=r8), INTENT(INOUT) :: tran(m,n)   !  total transmittance (tran)
    !
    !-----------------------------------------------------------------------
    !
    !  PURPOSE:
    !
    !  Calculate water vapor transmittance using the k-distribution
    !  method.
    !

    !-----------------------------------------------------------------------
    !
    !  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
    !          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
    !              Chou (1996)
    !
    !  MODIFICATION HISTORY:
    !
    !  03/15/1996 (Yuhe Liu)
    !  Adopted the original code and formatted it in accordance with the
    !  ARPS coding standard.
    !
    !-----------------------------------------------------------------------
    !
    !  ORIGINAL COMMENTS:
    !
    !**********************************************************************
    !   compute water vapor transmittance between levels k1 and k2 for
    !   m x n soundings using the k-distribution method.
    !
    !   computations follow eqs. (34), (46), (50) and (52).
    !
    !---- input parameters
    !  spectral band (ib)
    !  number of grid intervals in zonal direction (m)
    !  number of grid intervals in meridional direction (n)
    !  number of levels (np)
    !  current level (k)
    !  exponentials for line absorption (h2oexp)
    !  exponentials for continuum absorption (conexp)
    !
    !---- updated parameters
    !  transmittance between levels k1 and k2 due to
    !    water vapor line absorption (th2o)
    !  transmittance between levels k1 and k2 due to
    !    water vapor continuum absorption (tcon)
    !  total transmittance (tran)
    !
    !**********************************************************************
    !
    !-----------------------------------------------------------------------
    !
    !  Variable declarations
    !
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: ib
    INTEGER, INTENT(IN   ) :: m
    INTEGER, INTENT(IN   ) :: np
    INTEGER, INTENT(IN   ) :: k

    !---- input parameters ------

    REAL(KIND=r8), INTENT(IN   ) :: conexp(m,np,3)
    REAL(KIND=r8), INTENT(IN   ) :: h2oexp(m,np,6)

    !---- updated parameters -----

    REAL(KIND=r8), INTENT(INOUT) :: th2o(m,6)
    REAL(KIND=r8), INTENT(INOUT) :: tcon(m,3)
    REAL(KIND=r8), INTENT(INOUT) :: tran(m)



    !---- temporary variable -----

    REAL(KIND=r8) :: trnth2o
    INTEGER :: i

    !---- static data -----

    !-----fkw is the planck-weighted k-distribution function due to h2o
    !  line absorption given in table 4 of Chou and Suarez (1995).
    !  the k-distribution function for the third band, fkw(*,3), is not used

    REAL(KIND=r8), PARAMETER :: fkw(6,8) = RESHAPE(SOURCE=(/&
         0.2747_r8,0.2717_r8,0.2752_r8,0.1177_r8,0.0352_r8,0.0255_r8, &
         0.1521_r8,0.3974_r8,0.1778_r8,0.1826_r8,0.0374_r8,0.0527_r8, &
         1.00_r8,  1.00_r8,  1.00_r8,  1.00_r8,  1.00_r8,  1.00_r8, &
         0.4654_r8,0.2991_r8,0.1343_r8,0.0646_r8,0.0226_r8,0.0140_r8, &
         0.5543_r8,0.2723_r8,0.1131_r8,0.0443_r8,0.0160_r8,0.0000_r8, &
         0.1846_r8,0.2732_r8,0.2353_r8,0.1613_r8,0.1146_r8,0.0310_r8, &
         0.0740_r8,0.1636_r8,0.4174_r8,0.1783_r8,0.1101_r8,0.0566_r8, &
         0.1437_r8,0.2197_r8,0.3185_r8,0.2351_r8,0.0647_r8,0.0183_r8/),SHAPE=(/6,8/))

    !-----gkw is the planck-weighted k-distribution function due to h2o
    !  line absorption in the 3 subbands (800-720,620-720,540-620 /cm)
    !  of band 3 given in table 7.  Note that the order of the sub-bands
    !  is reversed.

    REAL(KIND=r8), PARAMETER ::  gkw(6,3)= RESHAPE(SOURCE=(/&
         0.1782_r8,0.0593_r8,0.0215_r8,0.0068_r8,0.0022_r8,0.0000_r8, &
         0.0923_r8,0.1675_r8,0.0923_r8,0.0187_r8,0.0178_r8,0.0000_r8, &
         0.0000_r8,0.1083_r8,0.1581_r8,0.0455_r8,0.0274_r8,0.0041_r8/),SHAPE=(/6,3/))

    !-----ne is the number of terms used in each band to compute water vapor
    !  continuum transmittance (table 6).

    INTEGER, PARAMETER ::  NE(8)=RESHAPE(SOURCE=(/0,0,3,1,1,1,0,0/),SHAPE=(/8/))
    !
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !
    !  Beginning of executable code...
    !
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !

    !-----tco2 are the six exp factors between levels k1 and k2
    !  tran is the updated total transmittance between levels k1 and k2


    !-----th2o is the 6 exp factors between levels k1 and k2 due to
    !  h2o line absorption.

    !-----tcon is the 3 exp factors between levels k1 and k2 due to
    !  h2o continuum absorption.

    !-----trnth2o is the total transmittance between levels k1 and k2 due
    !  to both line and continuum absorption computed from eq. (52).

    DO i=1,m
       th2o(i,1) = th2o(i,1)*h2oexp(i,k,1)
       th2o(i,2) = th2o(i,2)*h2oexp(i,k,2)
       th2o(i,3) = th2o(i,3)*h2oexp(i,k,3)
       th2o(i,4) = th2o(i,4)*h2oexp(i,k,4)
       th2o(i,5) = th2o(i,5)*h2oexp(i,k,5)
       th2o(i,6) = th2o(i,6)*h2oexp(i,k,6)
    END DO


    IF (NE(ib) == 0) THEN


       DO i=1,m

          trnth2o      =(fkw(1,ib)*th2o(i,1)                            &
               + fkw(2,ib)*th2o(i,2)                            &
               + fkw(3,ib)*th2o(i,3)                            &
               + fkw(4,ib)*th2o(i,4)                            &
               + fkw(5,ib)*th2o(i,5)                            &
               + fkw(6,ib)*th2o(i,6))

          tran(i)=tran(i)*trnth2o

       END DO

    ELSE IF (NE(ib) == 1) THEN


       DO i=1,m

          tcon(i,1)= tcon(i,1)*conexp(i,k,1)

          trnth2o      =(fkw(1,ib)*th2o(i,1)                            &
               + fkw(2,ib)*th2o(i,2)                            &
               + fkw(3,ib)*th2o(i,3)                            &
               + fkw(4,ib)*th2o(i,4)                            &
               + fkw(5,ib)*th2o(i,5)                            &
               + fkw(6,ib)*th2o(i,6))*tcon(i,1)

          tran(i)=tran(i)*trnth2o

       END DO

    ELSE

       DO i=1,m

          tcon(i,1)= tcon(i,1)*conexp(i,k,1)
          tcon(i,2)= tcon(i,2)*conexp(i,k,2)
          tcon(i,3)= tcon(i,3)*conexp(i,k,3)

          trnth2o      = (  gkw(1,1)*th2o(i,1)                          &
               + gkw(2,1)*th2o(i,2)                          &
               + gkw(3,1)*th2o(i,3)                          &
               + gkw(4,1)*th2o(i,4)                          &
               + gkw(5,1)*th2o(i,5)                          &
               + gkw(6,1)*th2o(i,6) ) * tcon(i,1)          &
               + (  gkw(1,2)*th2o(i,1)                          &
               + gkw(2,2)*th2o(i,2)                          &
               + gkw(3,2)*th2o(i,3)                          &
               + gkw(4,2)*th2o(i,4)                          &
               + gkw(5,2)*th2o(i,5)                          &
               + gkw(6,2)*th2o(i,6) ) * tcon(i,2)          &
               + (  gkw(1,3)*th2o(i,1)                          &
               + gkw(2,3)*th2o(i,2)                          &
               + gkw(3,3)*th2o(i,3)                          &
               + gkw(4,3)*th2o(i,4)                          &
               + gkw(5,3)*th2o(i,5)                          &
               + gkw(6,3)*th2o(i,6) ) * tcon(i,3)

          tran(i)=tran(i)*trnth2o

       END DO

    END IF

    RETURN
  END SUBROUTINE wvkdis



  !
  !
  !##################################################################
  !##################################################################
  !######                                                      ######
  !######                SUBROUTINE CO2KDIS                    ######
  !######                                                      ######
  !######                     Developed by                     ######
  !######                                                      ######
  !######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
  !######                                                      ######
  !######     Center for Analysis and Prediction of Storms     ######
  !######               University of Oklahoma                 ######
  !######                                                      ######
  !##################################################################
  !##################################################################
  !


  SUBROUTINE co2kdis( &
       m       , &! INTEGER      , INTENT(IN   ) :: m  !        number of grid intervals in zonal direction (m)
       np      , &! INTEGER      , INTENT(IN   ) :: np
       k       , &! INTEGER      , INTENT(IN   ) :: k
       co2exp  , &! REAL(KIND=r8), INTENT(IN   ) :: co2exp(m,n,np,6,2)
       tco2    , &! REAL(KIND=r8), INTENT(INOUT) :: tco2(m,n,6,2)!  for the various values of the absorption coefficient (tco2)
       tran      )! REAL(KIND=r8), INTENT(INOUT) :: tran(m,n)!   total transmittance (tran)
    !
    !-----------------------------------------------------------------------
    !
    !  PURPOSE:
    !
    !  Calculate co2 transmittances using the k-distribution method
    !
    !
    !-----------------------------------------------------------------------
    !
    !  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
    !          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
    !              Chou (1996)
    !
    !  MODIFICATION HISTORY:
    !
    !  03/15/1996 (Yuhe Liu)
    !  Adopted the original code and formatted it in accordance with the
    !  ARPS coding standard.
    !
    !-----------------------------------------------------------------------
    !
    !  ORIGINAL COMMENTS:
    !
    !**********************************************************************
    !   compute co2 transmittances between levels k1 and k2 for m x n soundings
    !   using the k-distribution method with linear pressure scaling.
    !
    !   computations follow eq. (34).
    !
    !---- input parameters
    !   number of grid intervals in zonal direction (m)
    !   number of grid intervals in meridional direction (n)
    !
    !---- updated parameters
    !   transmittance between levels k1 and k2 due to co2 absorption
    !  for the various values of the absorption coefficient (tco2)
    !   total transmittance (tran)
    !
    !**********************************************************************
    !
    !-----------------------------------------------------------------------
    !
    !  Variable declarations
    !
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: m
    INTEGER, INTENT(IN   ) :: np
    INTEGER, INTENT(IN   ) :: k

    !---- input parameters -----

    REAL(KIND=r8), INTENT(IN   ) :: co2exp(m,np,6,2)

    !---- updated parameters -----

    REAL(KIND=r8), INTENT(INOUT) :: tco2(m,6,2)
    REAL(KIND=r8), INTENT(INOUT) :: tran(m)

    !---- temporary variable -----

    REAL(KIND=r8) :: xc
    INTEGER :: i

    !-----gkc is the planck-weighted co2 k-distribution function
    !  in the band-wing and band-center regions given in table 7.
    !  for computing efficiency, sub-bands 3a and 3c are combined.

    REAL(KIND=r8), PARAMETER :: gkc(6,2) = RESHAPE(SOURCE=(/&
         0.1395_r8,0.1407_r8,0.1549_r8,0.1357_r8,0.0182_r8,0.0220_r8,&
         0.0766_r8,0.1372_r8,0.1189_r8,0.0335_r8,0.0169_r8,0.0059_r8/),SHAPE=(/6,2/))


    !
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !
    !  Beginning of executable code...
    !
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !

    !-----tco2 is the 6 exp factors between levels k1 and k2.
    !  xc is the total co2 transmittance given by eq. (53).

    DO i=1,m

       !-----band-wings

       tco2(i,1,1)=tco2(i,1,1)*co2exp(i,k,1,1)
       xc=             gkc(1,1)*tco2(i,1,1)

       tco2(i,2,1)=tco2(i,2,1)*co2exp(i,k,2,1)
       xc=xc+gkc(2,1)*tco2(i,2,1)

       tco2(i,3,1)=tco2(i,3,1)*co2exp(i,k,3,1)
       xc=xc+gkc(3,1)*tco2(i,3,1)

       tco2(i,4,1)=tco2(i,4,1)*co2exp(i,k,4,1)
       xc=xc+gkc(4,1)*tco2(i,4,1)

       tco2(i,5,1)=tco2(i,5,1)*co2exp(i,k,5,1)
       xc=xc+gkc(5,1)*tco2(i,5,1)

       tco2(i,6,1)=tco2(i,6,1)*co2exp(i,k,6,1)
       xc=xc+gkc(6,1)*tco2(i,6,1)

       !-----band-center region

       tco2(i,1,2)=tco2(i,1,2)*co2exp(i,k,1,2)
       xc=xc+gkc(1,2)*tco2(i,1,2)

       tco2(i,2,2)=tco2(i,2,2)*co2exp(i,k,2,2)
       xc=xc+gkc(2,2)*tco2(i,2,2)

       tco2(i,3,2)=tco2(i,3,2)*co2exp(i,k,3,2)
       xc=xc+gkc(3,2)*tco2(i,3,2)

       tco2(i,4,2)=tco2(i,4,2)*co2exp(i,k,4,2)
       xc=xc+gkc(4,2)*tco2(i,4,2)

       tco2(i,5,2)=tco2(i,5,2)*co2exp(i,k,5,2)
       xc=xc+gkc(5,2)*tco2(i,5,2)

       tco2(i,6,2)=tco2(i,6,2)*co2exp(i,k,6,2)
       xc=xc+gkc(6,2)*tco2(i,6,2)

       tran(i)=tran(i)*xc

    END DO

    RETURN
  END SUBROUTINE co2kdis

  !
  !
  !##################################################################
  !##################################################################
  !######                                                      ######
  !######                SUBROUTINE TABLUP                     ######
  !######                                                      ######
  !######                     Developed by                     ######
  !######                                                      ######
  !######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
  !######                                                      ######
  !######     Center for Analysis and Prediction of Storms     ######
  !######               University of Oklahoma                 ######
  !######                                                      ######
  !##################################################################
  !##################################################################
  !

  SUBROUTINE tablup( &
       k1     , &!INTEGER      ,INTENT(IN   ) :: k1 indices for pressure levels (k1 and k2)
       k2     , &!INTEGER      ,INTENT(IN   ) :: k2 indices for pressure levels (k1 and k2)
       m      , &!INTEGER      ,INTENT(IN   ) :: m  number of grid intervals in zonal direction (m)
       np     , &!INTEGER      ,INTENT(IN   ) :: np number of atmospheric layers (np)
       nx     , &!INTEGER      ,INTENT(IN   ) :: nx number of pressure intervals in the table (nx)
       nh     , &!INTEGER      ,INTENT(IN   ) :: nh number of absorber amount intervals in the table (nh)
       nt     , &!INTEGER      ,INTENT(IN   ) :: nt number of tables copied (nt)
       sabs   , &!REAL(KIND=r8),INTENT(IN   ) :: sabs(m,n,np+1) column-integrated absorber amount (sabs)
       spre   , &!REAL(KIND=r8),INTENT(IN   ) :: spre(m,n,np+1) column absorber amount-weighted pressure (spre)
       stem   , &!REAL(KIND=r8),INTENT(IN   ) :: stem(m,n,np+1)  column absorber amount-weighted temperature (stem)
       w1     , &!REAL(KIND=r8),INTENT(IN   ) :: w1                first value of absorber amount (log10) in the table (w1)
       p1     , &!REAL(KIND=r8),INTENT(IN   ) :: p1                first value of pressure (log10) in the table (p1)
       dwe    , &!REAL(KIND=r8),INTENT(IN   ) :: dwe           size of the interval of absorber amount (log10) in the table (dwe)
       dpe    , &!REAL(KIND=r8),INTENT(IN   ) :: dpe           size of the interval of pressure (log10) in the table (dpe)
       coef1  , &!REAL(KIND=r8),INTENT(IN   ) :: coef1(nx,nh,nt) pre-computed coefficients (coef1, coef2, and coef3)
       coef2  , &!REAL(KIND=r8),INTENT(IN   ) :: coef2(nx,nh,nt) pre-computed coefficients (coef1, coef2, and coef3)
       coef3  , &!REAL(KIND=r8),INTENT(IN   ) :: coef3(nx,nh,nt) pre-computed coefficients (coef1, coef2, and coef3)
       tran     )!REAL(KIND=r8),INTENT(INOUT) :: tran(m,n)  updated transmittance (tran)
    !
    !-----------------------------------------------------------------------
    !
    !  PURPOSE:
    !
    !  Calculate water vapor, co2, and o3 transmittances using table
    !  look-up. rlwopt = 0, or high=.false.
    !
    !-----------------------------------------------------------------------
    !
    !  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
    !          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
    !              Chou (1996)
    !
    !  MODIFICATION HISTORY:
    !
    !  03/15/1996 (Yuhe Liu)
    !  Adopted the original code and formatted it in accordance with the
    !  ARPS coding standard.
    !
    !-----------------------------------------------------------------------
    !
    !  ORIGINAL COMMENTS:
    !
    !**********************************************************************
    !   compute water vapor, co2, and o3 transmittances between levels k1 and k2
    !   using table look-up for m x n soundings.
    !
    !   Calculations follow Eq. (40) of Chou and Suarez (1995)
    !
    !---- input ---------------------
    !  indices for pressure levels (k1 and k2)
    !  number of grid intervals in zonal direction (m)
    !  number of grid intervals in meridional direction (n)
    !  number of atmospheric layers (np)
    !  number of pressure intervals in the table (nx)
    !  number of absorber amount intervals in the table (nh)
    !  number of tables copied (nt)
    !  column-integrated absorber amount (sabs)
    !  column absorber amount-weighted pressure (spre)
    !  column absorber amount-weighted temperature (stem)
    !  first value of absorber amount (log10) in the table (w1)
    !  first value of pressure (log10) in the table (p1)
    !  size of the interval of absorber amount (log10) in the table (dwe)
    !  size of the interval of pressure (log10) in the table (dpe)
    !  pre-computed coefficients (coef1, coef2, and coef3)
    !
    !---- updated ---------------------
    !  transmittance (tran)
    !
    !  Note:
    !   (1) units of sabs are g/cm**2 for water vapor and (cm-atm)stp for co2 and o3.
    !   (2) units of spre and stem are, respectively, mb and K.
    !   (3) there are nt identical copies of the tables (coef1, coef2, and
    !    coef3).  the prupose of using the multiple copies of tables is
    !    to increase the speed in parallel (vectorized) computations.
    !    if such advantage does not exist, nt can be set to 1.
    !
    !**********************************************************************
    !
    !-----------------------------------------------------------------------
    !
    !  Variable declarations
    !
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    INTEGER     ,INTENT(IN   ) :: k1
    INTEGER     ,INTENT(IN   ) :: k2
    INTEGER     ,INTENT(IN   ) :: m
    INTEGER     ,INTENT(IN   ) :: np
    INTEGER     ,INTENT(IN   ) :: nx
    INTEGER     ,INTENT(IN   ) :: nh
    INTEGER     ,INTENT(IN   ) :: nt

    !---- input parameters -----
    REAL(KIND=r8),INTENT(IN   ) :: sabs(m,np+1)
    REAL(KIND=r8),INTENT(IN   ) :: spre(m,np+1)
    REAL(KIND=r8),INTENT(IN   ) :: stem(m,np+1)
    REAL(KIND=r8),INTENT(IN   ) :: w1
    REAL(KIND=r8),INTENT(IN   ) :: p1
    REAL(KIND=r8),INTENT(IN   ) :: dwe
    REAL(KIND=r8),INTENT(IN   ) :: dpe
    REAL(KIND=r8),INTENT(IN   ) :: coef1(nx,nh,nt)
    REAL(KIND=r8),INTENT(IN   ) :: coef2(nx,nh,nt)
    REAL(KIND=r8),INTENT(IN   ) :: coef3(nx,nh,nt)

    !---- update parameter -----

    REAL(KIND=r8), INTENT(INOUT) :: tran(m)

    !---- temporary variables -----

    REAL(KIND=r8) :: x1
    REAL(KIND=r8) :: x2
    REAL(KIND=r8) :: x3
    REAL(KIND=r8) :: we
    REAL(KIND=r8) :: pe
    REAL(KIND=r8) :: fw
    REAL(KIND=r8) :: fp
    REAL(KIND=r8) :: pa
    REAL(KIND=r8) :: pb
    REAL(KIND=r8) :: pc
    REAL(KIND=r8) :: ax
    REAL(KIND=r8) :: ba
    REAL(KIND=r8) :: bb
    REAL(KIND=r8) :: t1
    REAL(KIND=r8) :: ca
    REAL(KIND=r8) :: cb
    REAL(KIND=r8) :: t2
    INTEGER :: iw,ip,nn
    INTEGER :: i

    !
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !
    !  Beginning of executable code...
    !
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !
    DO i=1,m

       nn=MOD(i,nt)+1

       x1=sabs(i,k2)-sabs(i,k1)
       x2=(spre(i,k2)-spre(i,k1))/x1
       x3=(stem(i,k2)-stem(i,k1))/x1

       we=(LOG10(x1)-w1)/dwe
       pe=(LOG10(x2)-p1)/dpe

       we=MAX(we,w1-2.0_r8*dwe)
       pe=MAX(pe,p1)

       iw=INT(we+1.5_r8)
       ip=INT(pe+1.5_r8)

       iw=MIN(iw,nh-1)
       iw=MAX(iw, 2)

       ip=MIN(ip,nx-1)
       ip=MAX(ip, 1)

       fw=we-FLOAT(iw-1)
       fp=pe-FLOAT(ip-1)

       !-----linear interpolation in pressure

       pa = coef1(ip,iw-1,nn)*(1.0_r8-fp)+coef1(ip+1,iw-1,nn)*fp
       pb = coef1(ip,iw,  nn)*(1.0_r8-fp)+coef1(ip+1,iw,  nn)*fp
       pc = coef1(ip,iw+1,nn)*(1.0_r8-fp)+coef1(ip+1,iw+1,nn)*fp

       !-----quadratic interpolation in absorber amount for coef1

       ax = (-pa*(1.0_r8-fw)+pc*(1.0_r8+fw)) *fw*0.5_r8 + pb*(1.0_r8-fw*fw)

       !-----linear interpolation in absorber amount for coef2 and coef3

       ba = coef2(ip,iw,  nn)*(1.0_r8-fp)+coef2(ip+1,iw,  nn)*fp
       bb = coef2(ip,iw+1,nn)*(1.0_r8-fp)+coef2(ip+1,iw+1,nn)*fp
       t1 = ba*(1.0_r8-fw) + bb*fw

       ca = coef3(ip,iw,  nn)*(1.0_r8-fp)+coef3(ip+1,iw,  nn)*fp
       cb = coef3(ip,iw+1,nn)*(1.0_r8-fp)+coef3(ip+1,iw+1,nn)*fp
       t2 = ca*(1.0_r8-fw) + cb*fw

       !-----update the total transmittance between levels k1 and k2

       tran(i)= (ax + (t1+t2*x3) * x3)*tran(i)

    END DO

    RETURN
  END SUBROUTINE tablup

  !
  !##################################################################
  !##################################################################
  !######                                                      ######
  !######                FUNCTION EXPMN                        ######
  !######                                                      ######
  !######                     Developed by                     ######
  !######                                                      ######
  !######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
  !######                                                      ######
  !######     Center for Analysis and Prediction of Storms     ######
  !######               University of Oklahoma                 ######
  !######                                                      ######
  !##################################################################
  !##################################################################
  !
  REAL(KIND=r8) FUNCTION expmn(fin)
    !
    !-----------------------------------------------------------------------
    !
    !  PURPOSE:
    !
    !  Calculate exponential for arguments in the range 0> fin > -10.
    !
    !-----------------------------------------------------------------------
    !
    !  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
    !          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
    !              Chou (1996)
    !
    !  MODIFICATION HISTORY:
    !
    !  03/15/1996 (Yuhe Liu)
    !  Adopted the original code and formatted it in accordance with the
    !  ARPS coding standard.
    !
    !-----------------------------------------------------------------------
    !
    !  ORIGINAL COMMENTS:
    !
    !**************************************************************************
    ! compute exponential for arguments in the range 0> fin > -10.
    !
    !-----------------------------------------------------------------------
    !
    !  Variable declarations
    !
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE

    REAL(KIND=r8), INTENT(IN) :: fin

    REAL(KIND=r8), PARAMETER  :: one    = 1.0_r8
    REAL(KIND=r8), PARAMETER  :: expmin = -10.0_r8

    REAL(KIND=r8), PARAMETER  :: e1=1.0_r8
    REAL(KIND=r8), PARAMETER  :: e2=-2.507213E-1_r8
    REAL(KIND=r8), PARAMETER  :: e3=2.92732E-2_r8
    REAL(KIND=r8), PARAMETER  :: e4=-3.827800E-3_r8

    REAL(KIND=r8) :: tmp
    !  REAL(KIND=r8) :: expmn
    !
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !
    !  Beginning of executable code...
    !
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !
    tmp = MAX( fin, expmin )
    expmn = ((e4*tmp + e3)*tmp+e2)*tmp+e1
    expmn = expmn * expmn
    expmn = one / (expmn * expmn)

    RETURN
  END FUNCTION expmn
  !




  !
  !##################################################################
  !##################################################################
  !######                                                      ######
  !######                BLOCK DATA RADDATA                    ######
  !######                                                      ######
  !######                     Developed by                     ######
  !######                                                      ######
  !######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
  !######                                                      ######
  !######     Center for Analysis and Prediction of Storms     ######
  !######                University of Oklahoma                ######
  !######                                                      ######
  !##################################################################
  !##################################################################
  !
  SUBROUTINE  raddata()!
    !-----------------------------------------------------------------------
    !
    !  PURPOSE:
    !
    !  Initialize pre-calculated look-up tables used in radiation
    !  computation.
    !
    !-----------------------------------------------------------------------
    !
    !  AUTHOR: Yuhe Liu
    !  03/15/1996
    !
    !  Combined all files of look-up tables used in radiation computation
    !  into this file. Different tables are identified by different
    !  COMMON block names.
    !
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE

    INTEGER :: i, j,it,k
    !
    !-----------------------------------------------------------------------
    !
    !  The following DATA statements originally came from file
    !  "h2o.tran3", which define pre-computed tables used for h2o (bands
    !  1, 2, and 7 only) transmittance calculations.
    !
    !-----------------------------------------------------------------------
    !
    !  integer nx,no,nc,nh,nt
    !  parameter (nx=26,no=21,nc=24,nh=31,nt=7)
    !
    !-----------------------------------------------------------------------
    !
    ! h2o.tran3

    !  REAL :: h11(26,31,7),h12(26,31,7),h13(26,31,7)
    !  REAL :: h21(26,31,7),h22(26,31,7),h23(26,31,7)
    !  REAL :: h71(26,31,7),h72(26,31,7),h73(26,31,7)
    !  INTEGER, PARAMETER :: nx=26
    !  INTEGER, PARAMETER :: no=21
    !  INTEGER, PARAMETER :: nc=24
    !  INTEGER, PARAMETER :: nh=31
    !  INTEGER, PARAMETER :: nt=7
    !  REAL(KIND=r8) :: c1 (nx,nc,nt),c2 (nx,nc,nt),c3 (nx,nc,nt)
    !  REAL(KIND=r8) :: o1 (nx,no,nt),o2 (nx,no,nt),o3 (nx,no,nt)
    !  REAL(KIND=r8) :: h11(nx,nh,nt),h12(nx,nh,nt),h13(nx,nh,nt)
    !  REAL(KIND=r8) :: h21(nx,nh,nt),h22(nx,nh,nt),h23(nx,nh,nt)
    !  REAL(KIND=r8) :: h71(nx,nh,nt),h72(nx,nh,nt),h73(nx,nh,nt)

    REAL(KIND=r8) :: data1(nx*3*nc)
    REAL(KIND=r8) :: data2(nx*3*no)
    REAL(KIND=r8) :: data3(nx*3*nh)
    REAL(KIND=r8) :: data4(nx*3*nh)
    REAL(KIND=r8) :: data5(nx*3*nh)

    data3(1:nx*3*nh)=RESHAPE(SOURCE=(/&
                                !DATA ((h11(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 1_r8, 1)/                               &
         0.6160E-04_r8,  0.9815E-04_r8,  0.1474E-03_r8,  0.2092E-03_r8,  0.2823E-03_r8,   &
         0.3662E-03_r8,  0.4615E-03_r8,  0.5710E-03_r8,  0.6998E-03_r8,  0.8554E-03_r8,   &
         0.1049E-02_r8,  0.1295E-02_r8,  0.1612E-02_r8,  0.2021E-02_r8,  0.2548E-02_r8,   &
         0.3230E-02_r8,  0.4123E-02_r8,  0.5306E-02_r8,  0.6887E-02_r8,  0.9021E-02_r8,   &
         0.1193E-01_r8,  0.1590E-01_r8,  0.2135E-01_r8,  0.2885E-01_r8,  0.3914E-01_r8,   &
         0.5317E-01_r8,  0.7223E-01_r8,  0.9800E-01_r8,  0.1326E+00_r8,  0.1783E+00_r8,   &
         0.2373E+00_r8,&
                                ! DATA ((h12(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 1_r8, 1)_r8,&                               &
         -0.2021E-06_r8, -0.3628E-06_r8, -0.5891E-06_r8, -0.8735E-06_r8, -0.1204E-05_r8,  &
         -0.1579E-05_r8, -0.2002E-05_r8, -0.2494E-05_r8, -0.3093E-05_r8, -0.3852E-05_r8,  &
         -0.4835E-05_r8, -0.6082E-05_r8, -0.7591E-05_r8, -0.9332E-05_r8, -0.1128E-04_r8,  &
         -0.1347E-04_r8, -0.1596E-04_r8, -0.1890E-04_r8, -0.2241E-04_r8, -0.2672E-04_r8,  &
         -0.3208E-04_r8, -0.3884E-04_r8, -0.4747E-04_r8, -0.5854E-04_r8, -0.7272E-04_r8,  &
         -0.9092E-04_r8, -0.1146E-03_r8, -0.1458E-03_r8, -0.1877E-03_r8, -0.2435E-03_r8,  &
         -0.3159E-03_r8,&
                                ! DATA ((h13(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 1_r8, 1)_r8,&                               &
         0.5907E-09_r8,  0.8541E-09_r8,  0.1095E-08_r8,  0.1272E-08_r8,  0.1297E-08_r8,   &
         0.1105E-08_r8,  0.6788E-09_r8, -0.5585E-10_r8, -0.1147E-08_r8, -0.2746E-08_r8,   &
         -0.5001E-08_r8, -0.7715E-08_r8, -0.1037E-07_r8, -0.1227E-07_r8, -0.1287E-07_r8,  &
         -0.1175E-07_r8, -0.8517E-08_r8, -0.2920E-08_r8,  0.4786E-08_r8,  0.1407E-07_r8,  &
         0.2476E-07_r8,  0.3781E-07_r8,  0.5633E-07_r8,  0.8578E-07_r8,  0.1322E-06_r8,   &
         0.2013E-06_r8,  0.3006E-06_r8,  0.4409E-06_r8,  0.6343E-06_r8,  0.8896E-06_r8,   &
         0.1216E-05_r8,&
                                ! DATA ((h11(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 2_r8, 2)_r8,&                               &
         0.6166E-04_r8,  0.9828E-04_r8,  0.1477E-03_r8,  0.2097E-03_r8,  0.2833E-03_r8,   &
         0.3680E-03_r8,  0.4648E-03_r8,  0.5768E-03_r8,  0.7101E-03_r8,  0.8736E-03_r8,   &
         0.1080E-02_r8,  0.1348E-02_r8,  0.1700E-02_r8,  0.2162E-02_r8,  0.2767E-02_r8,   &
         0.3563E-02_r8,  0.4621E-02_r8,  0.6039E-02_r8,  0.7953E-02_r8,  0.1056E-01_r8,   &
         0.1412E-01_r8,  0.1901E-01_r8,  0.2574E-01_r8,  0.3498E-01_r8,  0.4763E-01_r8,   &
         0.6484E-01_r8,  0.8815E-01_r8,  0.1196E+00_r8,  0.1614E+00_r8,  0.2157E+00_r8,   &
         0.2844E+00_r8,&
                                ! DATA ((h12(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 2_r8, 2)_r8,&                               &
         -0.2017E-06_r8, -0.3620E-06_r8, -0.5878E-06_r8, -0.8713E-06_r8, -0.1201E-05_r8,  &
         -0.1572E-05_r8, -0.1991E-05_r8, -0.2476E-05_r8, -0.3063E-05_r8, -0.3808E-05_r8,  &
         -0.4776E-05_r8, -0.6011E-05_r8, -0.7516E-05_r8, -0.9272E-05_r8, -0.1127E-04_r8,  &
         -0.1355E-04_r8, -0.1620E-04_r8, -0.1936E-04_r8, -0.2321E-04_r8, -0.2797E-04_r8,  &
         -0.3399E-04_r8, -0.4171E-04_r8, -0.5172E-04_r8, -0.6471E-04_r8, -0.8150E-04_r8,  &
         -0.1034E-03_r8, -0.1321E-03_r8, -0.1705E-03_r8, -0.2217E-03_r8, -0.2889E-03_r8,  &
         -0.3726E-03_r8,&
                                ! DATA ((h13(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 2_r8, 2)_r8,&                               &
         0.5894E-09_r8,  0.8519E-09_r8,  0.1092E-08_r8,  0.1267E-08_r8,  0.1289E-08_r8,   &
         0.1093E-08_r8,  0.6601E-09_r8, -0.7831E-10_r8, -0.1167E-08_r8, -0.2732E-08_r8,   &
         -0.4864E-08_r8, -0.7334E-08_r8, -0.9581E-08_r8, -0.1097E-07_r8, -0.1094E-07_r8,  &
         -0.8999E-08_r8, -0.4669E-08_r8,  0.2391E-08_r8,  0.1215E-07_r8,  0.2424E-07_r8,  &
         0.3877E-07_r8,  0.5711E-07_r8,  0.8295E-07_r8,  0.1218E-06_r8,  0.1793E-06_r8,   &
         0.2621E-06_r8,  0.3812E-06_r8,  0.5508E-06_r8,  0.7824E-06_r8,  0.1085E-05_r8,   &
         0.1462E-05_r8,&
                                ! DATA ((h11(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 3_r8, 3)_r8,&                               &
         0.6175E-04_r8,  0.9849E-04_r8,  0.1481E-03_r8,  0.2106E-03_r8,  0.2849E-03_r8,   &
         0.3708E-03_r8,  0.4698E-03_r8,  0.5857E-03_r8,  0.7257E-03_r8,  0.9006E-03_r8,   &
         0.1126E-02_r8,  0.1425E-02_r8,  0.1823E-02_r8,  0.2353E-02_r8,  0.3059E-02_r8,   &
         0.4002E-02_r8,  0.5270E-02_r8,  0.6984E-02_r8,  0.9316E-02_r8,  0.1251E-01_r8,   &
         0.1689E-01_r8,  0.2292E-01_r8,  0.3123E-01_r8,  0.4262E-01_r8,  0.5814E-01_r8,   &
         0.7921E-01_r8,  0.1077E+00_r8,  0.1458E+00_r8,  0.1957E+00_r8,  0.2595E+00_r8,   &
         0.3380E+00_r8,&
                                ! DATA ((h12(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 3_r8, 3)_r8,&                               &
         -0.2011E-06_r8, -0.3609E-06_r8, -0.5859E-06_r8, -0.8680E-06_r8, -0.1195E-05_r8,  &
         -0.1563E-05_r8, -0.1975E-05_r8, -0.2450E-05_r8, -0.3024E-05_r8, -0.3755E-05_r8,  &
         -0.4711E-05_r8, -0.5941E-05_r8, -0.7455E-05_r8, -0.9248E-05_r8, -0.1132E-04_r8,  &
         -0.1373E-04_r8, -0.1659E-04_r8, -0.2004E-04_r8, -0.2431E-04_r8, -0.2966E-04_r8,  &
         -0.3653E-04_r8, -0.4549E-04_r8, -0.5724E-04_r8, -0.7259E-04_r8, -0.9265E-04_r8,  &
         -0.1191E-03_r8, -0.1543E-03_r8, -0.2013E-03_r8, -0.2633E-03_r8, -0.3421E-03_r8,  &
         -0.4350E-03_r8,&
                                ! DATA ((h13(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 3_r8, 3)_r8,&                               &
         0.5872E-09_r8,  0.8484E-09_r8,  0.1087E-08_r8,  0.1259E-08_r8,  0.1279E-08_r8,   &
         0.1077E-08_r8,  0.6413E-09_r8, -0.9334E-10_r8, -0.1161E-08_r8, -0.2644E-08_r8,   &
         -0.4588E-08_r8, -0.6709E-08_r8, -0.8474E-08_r8, -0.9263E-08_r8, -0.8489E-08_r8,  &
         -0.5553E-08_r8,  0.1203E-09_r8,  0.9035E-08_r8,  0.2135E-07_r8,  0.3689E-07_r8,  &
         0.5610E-07_r8,  0.8097E-07_r8,  0.1155E-06_r8,  0.1649E-06_r8,  0.2350E-06_r8,   &
         0.3353E-06_r8,  0.4806E-06_r8,  0.6858E-06_r8,  0.9617E-06_r8,  0.1315E-05_r8,   &
         0.1741E-05_r8,&
                                ! DATA ((h11(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 4_r8, 4)_r8,&                               &
         0.6189E-04_r8,  0.9882E-04_r8,  0.1488E-03_r8,  0.2119E-03_r8,  0.2873E-03_r8,   &
         0.3752E-03_r8,  0.4776E-03_r8,  0.5993E-03_r8,  0.7490E-03_r8,  0.9403E-03_r8,   &
         0.1192E-02_r8,  0.1531E-02_r8,  0.1990E-02_r8,  0.2610E-02_r8,  0.3446E-02_r8,   &
         0.4576E-02_r8,  0.6109E-02_r8,  0.8196E-02_r8,  0.1105E-01_r8,  0.1498E-01_r8,   &
         0.2039E-01_r8,  0.2785E-01_r8,  0.3809E-01_r8,  0.5209E-01_r8,  0.7112E-01_r8,   &
         0.9688E-01_r8,  0.1315E+00_r8,  0.1773E+00_r8,  0.2363E+00_r8,  0.3100E+00_r8,   &
         0.3976E+00_r8,&
                                ! DATA ((h12(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 4_r8, 4)_r8,&                               &
         -0.2001E-06_r8, -0.3592E-06_r8, -0.5829E-06_r8, -0.8631E-06_r8, -0.1187E-05_r8,  &
         -0.1549E-05_r8, -0.1953E-05_r8, -0.2415E-05_r8, -0.2975E-05_r8, -0.3694E-05_r8,  &
         -0.4645E-05_r8, -0.5882E-05_r8, -0.7425E-05_r8, -0.9279E-05_r8, -0.1147E-04_r8,  &
         -0.1406E-04_r8, -0.1717E-04_r8, -0.2100E-04_r8, -0.2580E-04_r8, -0.3191E-04_r8,  &
         -0.3989E-04_r8, -0.5042E-04_r8, -0.6432E-04_r8, -0.8261E-04_r8, -0.1068E-03_r8,  &
         -0.1389E-03_r8, -0.1820E-03_r8, -0.2391E-03_r8, -0.3127E-03_r8, -0.4021E-03_r8,  &
         -0.5002E-03_r8,&
                                ! DATA ((h13(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 4_r8, 4)_r8,&                               &
         0.5838E-09_r8,  0.8426E-09_r8,  0.1081E-08_r8,  0.1249E-08_r8,  0.1267E-08_r8,   &
         0.1062E-08_r8,  0.6313E-09_r8, -0.8241E-10_r8, -0.1094E-08_r8, -0.2436E-08_r8,   &
         -0.4100E-08_r8, -0.5786E-08_r8, -0.6992E-08_r8, -0.7083E-08_r8, -0.5405E-08_r8,  &
         -0.1259E-08_r8,  0.6099E-08_r8,  0.1732E-07_r8,  0.3276E-07_r8,  0.5256E-07_r8,  &
         0.7756E-07_r8,  0.1103E-06_r8,  0.1547E-06_r8,  0.2159E-06_r8,  0.3016E-06_r8,   &
         0.4251E-06_r8,  0.6033E-06_r8,  0.8499E-06_r8,  0.1175E-05_r8,  0.1579E-05_r8,   &
         0.2044E-05_r8,&
                                ! DATA ((h11(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 5_r8, 5)_r8,&                               &
         0.6211E-04_r8,  0.9932E-04_r8,  0.1499E-03_r8,  0.2140E-03_r8,  0.2911E-03_r8,   &
         0.3820E-03_r8,  0.4895E-03_r8,  0.6196E-03_r8,  0.7834E-03_r8,  0.9973E-03_r8,   &
         0.1285E-02_r8,  0.1678E-02_r8,  0.2216E-02_r8,  0.2951E-02_r8,  0.3953E-02_r8,   &
         0.5319E-02_r8,  0.7186E-02_r8,  0.9743E-02_r8,  0.1326E-01_r8,  0.1811E-01_r8,   &
         0.2479E-01_r8,  0.3400E-01_r8,  0.4662E-01_r8,  0.6379E-01_r8,  0.8708E-01_r8,   &
         0.1185E+00_r8,  0.1603E+00_r8,  0.2147E+00_r8,  0.2835E+00_r8,  0.3667E+00_r8,   &
         0.4620E+00_r8,&
                                ! DATA ((h12(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 5_r8, 5)_r8,&                               &
         -0.1987E-06_r8, -0.3565E-06_r8, -0.5784E-06_r8, -0.8557E-06_r8, -0.1175E-05_r8,  &
         -0.1530E-05_r8, -0.1923E-05_r8, -0.2372E-05_r8, -0.2919E-05_r8, -0.3631E-05_r8,  &
         -0.4587E-05_r8, -0.5848E-05_r8, -0.7442E-05_r8, -0.9391E-05_r8, -0.1173E-04_r8,  &
         -0.1455E-04_r8, -0.1801E-04_r8, -0.2232E-04_r8, -0.2779E-04_r8, -0.3489E-04_r8,  &
         -0.4428E-04_r8, -0.5678E-04_r8, -0.7333E-04_r8, -0.9530E-04_r8, -0.1246E-03_r8,  &
         -0.1639E-03_r8, -0.2164E-03_r8, -0.2848E-03_r8, -0.3697E-03_r8, -0.4665E-03_r8,  &
         -0.5646E-03_r8,&
                                ! DATA ((h13(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 5_r8, 5)_r8,&                               &
         0.5785E-09_r8,  0.8338E-09_r8,  0.1071E-08_r8,  0.1239E-08_r8,  0.1256E-08_r8,   &
         0.1057E-08_r8,  0.6480E-09_r8, -0.1793E-10_r8, -0.9278E-09_r8, -0.2051E-08_r8,   &
         -0.3337E-08_r8, -0.4514E-08_r8, -0.5067E-08_r8, -0.4328E-08_r8, -0.1545E-08_r8,  &
         0.4100E-08_r8,  0.1354E-07_r8,  0.2762E-07_r8,  0.4690E-07_r8,  0.7190E-07_r8,   &
         0.1040E-06_r8,  0.1459E-06_r8,  0.2014E-06_r8,  0.2764E-06_r8,  0.3824E-06_r8,   &
         0.5359E-06_r8,  0.7532E-06_r8,  0.1047E-05_r8,  0.1424E-05_r8,  0.1873E-05_r8,   &
         0.2356E-05_r8,&
                                ! DATA ((h11(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 6_r8, 6)_r8,&                               &
         0.6246E-04_r8,  0.1001E-03_r8,  0.1515E-03_r8,  0.2171E-03_r8,  0.2970E-03_r8,   &
         0.3924E-03_r8,  0.5072E-03_r8,  0.6495E-03_r8,  0.8329E-03_r8,  0.1078E-02_r8,   &
         0.1413E-02_r8,  0.1876E-02_r8,  0.2516E-02_r8,  0.3399E-02_r8,  0.4612E-02_r8,   &
         0.6276E-02_r8,  0.8562E-02_r8,  0.1171E-01_r8,  0.1605E-01_r8,  0.2205E-01_r8,   &
         0.3032E-01_r8,  0.4167E-01_r8,  0.5717E-01_r8,  0.7821E-01_r8,  0.1067E+00_r8,   &
         0.1447E+00_r8,  0.1948E+00_r8,  0.2586E+00_r8,  0.3372E+00_r8,  0.4290E+00_r8,   &
         0.5295E+00_r8,&
                                ! DATA ((h12(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 6_r8, 6)_r8,&                               &
         -0.1964E-06_r8, -0.3526E-06_r8, -0.5717E-06_r8, -0.8451E-06_r8, -0.1158E-05_r8,  &
         -0.1504E-05_r8, -0.1886E-05_r8, -0.2322E-05_r8, -0.2861E-05_r8, -0.3576E-05_r8,  &
         -0.4552E-05_r8, -0.5856E-05_r8, -0.7529E-05_r8, -0.9609E-05_r8, -0.1216E-04_r8,  &
         -0.1528E-04_r8, -0.1916E-04_r8, -0.2408E-04_r8, -0.3043E-04_r8, -0.3880E-04_r8,  &
         -0.4997E-04_r8, -0.6488E-04_r8, -0.8474E-04_r8, -0.1113E-03_r8, -0.1471E-03_r8,  &
         -0.1950E-03_r8, -0.2583E-03_r8, -0.3384E-03_r8, -0.4326E-03_r8, -0.5319E-03_r8,  &
         -0.6244E-03_r8,&
                                ! DATA ((h13(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 6_r8, 6)_r8,&                               &
         0.5713E-09_r8,  0.8263E-09_r8,  0.1060E-08_r8,  0.1226E-08_r8,  0.1252E-08_r8,   &
         0.1076E-08_r8,  0.7149E-09_r8,  0.1379E-09_r8, -0.6043E-09_r8, -0.1417E-08_r8,   &
         -0.2241E-08_r8, -0.2830E-08_r8, -0.2627E-08_r8, -0.8950E-09_r8,  0.3231E-08_r8,  &
         0.1075E-07_r8,  0.2278E-07_r8,  0.4037E-07_r8,  0.6439E-07_r8,  0.9576E-07_r8,   &
         0.1363E-06_r8,  0.1886E-06_r8,  0.2567E-06_r8,  0.3494E-06_r8,  0.4821E-06_r8,   &
         0.6719E-06_r8,  0.9343E-06_r8,  0.1280E-05_r8,  0.1705E-05_r8,  0.2184E-05_r8,   &
         0.2651E-05_r8,&
                                ! DATA ((h11(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 7_r8, 7)_r8,&                               &
         0.6299E-04_r8,  0.1013E-03_r8,  0.1541E-03_r8,  0.2220E-03_r8,  0.3058E-03_r8,   &
         0.4078E-03_r8,  0.5334E-03_r8,  0.6928E-03_r8,  0.9032E-03_r8,  0.1190E-02_r8,   &
         0.1587E-02_r8,  0.2140E-02_r8,  0.2912E-02_r8,  0.3982E-02_r8,  0.5460E-02_r8,   &
         0.7500E-02_r8,  0.1031E-01_r8,  0.1420E-01_r8,  0.1958E-01_r8,  0.2700E-01_r8,   &
         0.3721E-01_r8,  0.5118E-01_r8,  0.7019E-01_r8,  0.9593E-01_r8,  0.1305E+00_r8,   &
         0.1763E+00_r8,  0.2354E+00_r8,  0.3091E+00_r8,  0.3969E+00_r8,  0.4952E+00_r8,   &
         0.5978E+00_r8,&
                                ! DATA ((h12(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 7_r8, 7)_r8,&                               &
         -0.1932E-06_r8, -0.3467E-06_r8, -0.5623E-06_r8, -0.8306E-06_r8, -0.1136E-05_r8,  &
         -0.1472E-05_r8, -0.1842E-05_r8, -0.2269E-05_r8, -0.2807E-05_r8, -0.3539E-05_r8,  &
         -0.4553E-05_r8, -0.5925E-05_r8, -0.7710E-05_r8, -0.9968E-05_r8, -0.1278E-04_r8,  &
         -0.1629E-04_r8, -0.2073E-04_r8, -0.2644E-04_r8, -0.3392E-04_r8, -0.4390E-04_r8,  &
         -0.5727E-04_r8, -0.7516E-04_r8, -0.9916E-04_r8, -0.1315E-03_r8, -0.1752E-03_r8,  &
         -0.2333E-03_r8, -0.3082E-03_r8, -0.3988E-03_r8, -0.4982E-03_r8, -0.5947E-03_r8,  &
         -0.6764E-03_r8,&
                                ! DATA ((h13(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 7_r8, 7)_r8,&                               &
         0.5612E-09_r8,  0.8116E-09_r8,  0.1048E-08_r8,  0.1222E-08_r8,  0.1270E-08_r8,   &
         0.1141E-08_r8,  0.8732E-09_r8,  0.4336E-09_r8, -0.6548E-10_r8, -0.4774E-09_r8,   &
         -0.7556E-09_r8, -0.6577E-09_r8,  0.4377E-09_r8,  0.3359E-08_r8,  0.9159E-08_r8,  &
         0.1901E-07_r8,  0.3422E-07_r8,  0.5616E-07_r8,  0.8598E-07_r8,  0.1251E-06_r8,   &
         0.1752E-06_r8,  0.2392E-06_r8,  0.3228E-06_r8,  0.4389E-06_r8,  0.6049E-06_r8,   &
         0.8370E-06_r8,  0.1150E-05_r8,  0.1547E-05_r8,  0.2012E-05_r8,  0.2493E-05_r8,   &
         0.2913E-05_r8,&
                                ! DATA ((h11(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 8_r8, 8)_r8,&                               &
         0.6378E-04_r8,  0.1032E-03_r8,  0.1579E-03_r8,  0.2293E-03_r8,  0.3190E-03_r8,   &
         0.4305E-03_r8,  0.5712E-03_r8,  0.7541E-03_r8,  0.1001E-02_r8,  0.1342E-02_r8,   &
         0.1819E-02_r8,  0.2489E-02_r8,  0.3427E-02_r8,  0.4734E-02_r8,  0.6548E-02_r8,   &
         0.9059E-02_r8,  0.1254E-01_r8,  0.1735E-01_r8,  0.2401E-01_r8,  0.3318E-01_r8,   &
         0.4577E-01_r8,  0.6292E-01_r8,  0.8620E-01_r8,  0.1176E+00_r8,  0.1594E+00_r8,   &
         0.2139E+00_r8,  0.2827E+00_r8,  0.3660E+00_r8,  0.4614E+00_r8,  0.5634E+00_r8,   &
         0.6649E+00_r8,&
                                ! DATA ((h12(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 8_r8, 8)_r8,&                               &
         -0.1885E-06_r8, -0.3385E-06_r8, -0.5493E-06_r8, -0.8114E-06_r8, -0.1109E-05_r8,  &
         -0.1436E-05_r8, -0.1796E-05_r8, -0.2219E-05_r8, -0.2770E-05_r8, -0.3535E-05_r8,  &
         -0.4609E-05_r8, -0.6077E-05_r8, -0.8016E-05_r8, -0.1051E-04_r8, -0.1367E-04_r8,  &
         -0.1768E-04_r8, -0.2283E-04_r8, -0.2955E-04_r8, -0.3849E-04_r8, -0.5046E-04_r8,  &
         -0.6653E-04_r8, -0.8813E-04_r8, -0.1173E-03_r8, -0.1569E-03_r8, -0.2100E-03_r8,  &
         -0.2794E-03_r8, -0.3656E-03_r8, -0.4637E-03_r8, -0.5629E-03_r8, -0.6512E-03_r8,  &
         -0.7167E-03_r8,&
                                ! DATA ((h13(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 8_r8, 8)_r8,&                               &
         0.5477E-09_r8,  0.8000E-09_r8,  0.1039E-08_r8,  0.1234E-08_r8,  0.1331E-08_r8,   &
         0.1295E-08_r8,  0.1160E-08_r8,  0.9178E-09_r8,  0.7535E-09_r8,  0.8301E-09_r8,   &
         0.1184E-08_r8,  0.2082E-08_r8,  0.4253E-08_r8,  0.8646E-08_r8,  0.1650E-07_r8,   &
         0.2920E-07_r8,  0.4834E-07_r8,  0.7564E-07_r8,  0.1125E-06_r8,  0.1606E-06_r8,   &
         0.2216E-06_r8,  0.2992E-06_r8,  0.4031E-06_r8,  0.5493E-06_r8,  0.7549E-06_r8,   &
         0.1035E-05_r8,  0.1400E-05_r8,  0.1843E-05_r8,  0.2327E-05_r8,  0.2774E-05_r8,   &
         0.3143E-05_r8,&
                                ! DATA ((h11(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 9_r8, 9)_r8,&                               &
         0.6495E-04_r8,  0.1059E-03_r8,  0.1635E-03_r8,  0.2400E-03_r8,  0.3381E-03_r8,   &
         0.4631E-03_r8,  0.6246E-03_r8,  0.8394E-03_r8,  0.1134E-02_r8,  0.1546E-02_r8,   &
         0.2126E-02_r8,  0.2944E-02_r8,  0.4093E-02_r8,  0.5699E-02_r8,  0.7934E-02_r8,   &
         0.1104E-01_r8,  0.1535E-01_r8,  0.2131E-01_r8,  0.2956E-01_r8,  0.4089E-01_r8,   &
         0.5636E-01_r8,  0.7739E-01_r8,  0.1058E+00_r8,  0.1439E+00_r8,  0.1939E+00_r8,   &
         0.2578E+00_r8,  0.3364E+00_r8,  0.4283E+00_r8,  0.5290E+00_r8,  0.6314E+00_r8,   &
         0.7292E+00_r8,&
                                ! DATA ((h12(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 9_r8, 9)_r8,&                               &
         -0.1822E-06_r8, -0.3274E-06_r8, -0.5325E-06_r8, -0.7881E-06_r8, -0.1079E-05_r8,  &
         -0.1398E-05_r8, -0.1754E-05_r8, -0.2184E-05_r8, -0.2763E-05_r8, -0.3581E-05_r8,  &
         -0.4739E-05_r8, -0.6341E-05_r8, -0.8484E-05_r8, -0.1128E-04_r8, -0.1490E-04_r8,  &
         -0.1955E-04_r8, -0.2561E-04_r8, -0.3364E-04_r8, -0.4438E-04_r8, -0.5881E-04_r8,  &
         -0.7822E-04_r8, -0.1045E-03_r8, -0.1401E-03_r8, -0.1884E-03_r8, -0.2523E-03_r8,  &
         -0.3335E-03_r8, -0.4289E-03_r8, -0.5296E-03_r8, -0.6231E-03_r8, -0.6980E-03_r8,  &
         -0.7406E-03_r8,&
                                ! DATA ((h13(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 9_r8, 9)_r8,&                               &
         0.5334E-09_r8,  0.7859E-09_r8,  0.1043E-08_r8,  0.1279E-08_r8,  0.1460E-08_r8,   &
         0.1560E-08_r8,  0.1618E-08_r8,  0.1657E-08_r8,  0.1912E-08_r8,  0.2569E-08_r8,   &
         0.3654E-08_r8,  0.5509E-08_r8,  0.8964E-08_r8,  0.1518E-07_r8,  0.2560E-07_r8,   &
         0.4178E-07_r8,  0.6574E-07_r8,  0.9958E-07_r8,  0.1449E-06_r8,  0.2031E-06_r8,   &
         0.2766E-06_r8,  0.3718E-06_r8,  0.5022E-06_r8,  0.6849E-06_r8,  0.9360E-06_r8,   &
         0.1268E-05_r8,  0.1683E-05_r8,  0.2157E-05_r8,  0.2625E-05_r8,  0.3020E-05_r8,   &
         0.3364E-05_r8,&
                                ! DATA ((h11(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=10_r8,10)_r8,&                               &
         0.6661E-04_r8,  0.1098E-03_r8,  0.1716E-03_r8,  0.2554E-03_r8,  0.3653E-03_r8,   &
         0.5088E-03_r8,  0.6986E-03_r8,  0.9557E-03_r8,  0.1313E-02_r8,  0.1816E-02_r8,   &
         0.2527E-02_r8,  0.3532E-02_r8,  0.4947E-02_r8,  0.6929E-02_r8,  0.9694E-02_r8,   &
         0.1354E-01_r8,  0.1888E-01_r8,  0.2628E-01_r8,  0.3647E-01_r8,  0.5043E-01_r8,   &
         0.6941E-01_r8,  0.9514E-01_r8,  0.1297E+00_r8,  0.1755E+00_r8,  0.2347E+00_r8,   &
         0.3084E+00_r8,  0.3962E+00_r8,  0.4947E+00_r8,  0.5974E+00_r8,  0.6973E+00_r8,   &
         0.7898E+00_r8,&
                                ! DATA ((h12(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=10_r8,10)_r8,&                               &
         -0.1742E-06_r8, -0.3134E-06_r8, -0.5121E-06_r8, -0.7619E-06_r8, -0.1048E-05_r8,  &
         -0.1364E-05_r8, -0.1725E-05_r8, -0.2177E-05_r8, -0.2801E-05_r8, -0.3694E-05_r8,  &
         -0.4969E-05_r8, -0.6748E-05_r8, -0.9161E-05_r8, -0.1236E-04_r8, -0.1655E-04_r8,  &
         -0.2203E-04_r8, -0.2927E-04_r8, -0.3894E-04_r8, -0.5192E-04_r8, -0.6936E-04_r8,  &
         -0.9294E-04_r8, -0.1250E-03_r8, -0.1686E-03_r8, -0.2271E-03_r8, -0.3027E-03_r8,  &
         -0.3944E-03_r8, -0.4951E-03_r8, -0.5928E-03_r8, -0.6755E-03_r8, -0.7309E-03_r8,  &
         -0.7417E-03_r8,&
                                ! DATA ((h13(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=10_r8,10)_r8,&                               &
         0.5179E-09_r8,  0.7789E-09_r8,  0.1071E-08_r8,  0.1382E-08_r8,  0.1690E-08_r8,   &
         0.1979E-08_r8,  0.2297E-08_r8,  0.2704E-08_r8,  0.3466E-08_r8,  0.4794E-08_r8,   &
         0.6746E-08_r8,  0.9739E-08_r8,  0.1481E-07_r8,  0.2331E-07_r8,  0.3679E-07_r8,   &
         0.5726E-07_r8,  0.8716E-07_r8,  0.1289E-06_r8,  0.1837E-06_r8,  0.2534E-06_r8,   &
         0.3424E-06_r8,  0.4609E-06_r8,  0.6245E-06_r8,  0.8495E-06_r8,  0.1151E-05_r8,   &
         0.1536E-05_r8,  0.1991E-05_r8,  0.2468E-05_r8,  0.2891E-05_r8,  0.3245E-05_r8,   &
         0.3580E-05_r8,&
                                ! DATA ((h11(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=11_r8,11)_r8,&                               &
         0.6888E-04_r8,  0.1151E-03_r8,  0.1828E-03_r8,  0.2766E-03_r8,  0.4028E-03_r8,   &
         0.5713E-03_r8,  0.7987E-03_r8,  0.1111E-02_r8,  0.1548E-02_r8,  0.2167E-02_r8,   &
         0.3044E-02_r8,  0.4285E-02_r8,  0.6035E-02_r8,  0.8490E-02_r8,  0.1192E-01_r8,   &
         0.1670E-01_r8,  0.2333E-01_r8,  0.3249E-01_r8,  0.4506E-01_r8,  0.6220E-01_r8,   &
         0.8546E-01_r8,  0.1168E+00_r8,  0.1587E+00_r8,  0.2131E+00_r8,  0.2820E+00_r8,   &
         0.3653E+00_r8,  0.4609E+00_r8,  0.5630E+00_r8,  0.6645E+00_r8,  0.7599E+00_r8,   &
         0.8458E+00_r8,&
                                ! DATA ((h12(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=11_r8,11)_r8,&                               &
         -0.1647E-06_r8, -0.2974E-06_r8, -0.4900E-06_r8, -0.7358E-06_r8, -0.1022E-05_r8,  &
         -0.1344E-05_r8, -0.1721E-05_r8, -0.2212E-05_r8, -0.2901E-05_r8, -0.3896E-05_r8,  &
         -0.5327E-05_r8, -0.7342E-05_r8, -0.1011E-04_r8, -0.1382E-04_r8, -0.1875E-04_r8,  &
         -0.2530E-04_r8, -0.3403E-04_r8, -0.4573E-04_r8, -0.6145E-04_r8, -0.8264E-04_r8,  &
         -0.1114E-03_r8, -0.1507E-03_r8, -0.2039E-03_r8, -0.2737E-03_r8, -0.3607E-03_r8,  &
         -0.4599E-03_r8, -0.5604E-03_r8, -0.6497E-03_r8, -0.7161E-03_r8, -0.7443E-03_r8,  &
         -0.7133E-03_r8,&
                                ! DATA ((h13(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=11_r8,11)_r8,&                               &
         0.5073E-09_r8,  0.7906E-09_r8,  0.1134E-08_r8,  0.1560E-08_r8,  0.2046E-08_r8,   &
         0.2589E-08_r8,  0.3254E-08_r8,  0.4107E-08_r8,  0.5481E-08_r8,  0.7602E-08_r8,   &
         0.1059E-07_r8,  0.1501E-07_r8,  0.2210E-07_r8,  0.3334E-07_r8,  0.5055E-07_r8,   &
         0.7629E-07_r8,  0.1134E-06_r8,  0.1642E-06_r8,  0.2298E-06_r8,  0.3133E-06_r8,   &
         0.4225E-06_r8,  0.5709E-06_r8,  0.7739E-06_r8,  0.1047E-05_r8,  0.1401E-05_r8,   &
         0.1833E-05_r8,  0.2308E-05_r8,  0.2753E-05_r8,  0.3125E-05_r8,  0.3467E-05_r8,   &
         0.3748E-05_r8,&
                                ! DATA ((h11(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=12_r8,12)_r8,&                               &
         0.7179E-04_r8,  0.1221E-03_r8,  0.1975E-03_r8,  0.3049E-03_r8,  0.4529E-03_r8,   &
         0.6547E-03_r8,  0.9312E-03_r8,  0.1315E-02_r8,  0.1855E-02_r8,  0.2620E-02_r8,   &
         0.3705E-02_r8,  0.5243E-02_r8,  0.7414E-02_r8,  0.1046E-01_r8,  0.1472E-01_r8,   &
         0.2065E-01_r8,  0.2888E-01_r8,  0.4019E-01_r8,  0.5566E-01_r8,  0.7668E-01_r8,   &
         0.1051E+00_r8,  0.1432E+00_r8,  0.1932E+00_r8,  0.2571E+00_r8,  0.3358E+00_r8,   &
         0.4278E+00_r8,  0.5285E+00_r8,  0.6310E+00_r8,  0.7289E+00_r8,  0.8184E+00_r8,   &
         0.8954E+00_r8,&
                                ! DATA ((h12(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=12_r8,12)_r8,&                               &
         -0.1548E-06_r8, -0.2808E-06_r8, -0.4683E-06_r8, -0.7142E-06_r8, -0.1008E-05_r8,  &
         -0.1347E-05_r8, -0.1758E-05_r8, -0.2306E-05_r8, -0.3083E-05_r8, -0.4214E-05_r8,  &
         -0.5851E-05_r8, -0.8175E-05_r8, -0.1140E-04_r8, -0.1577E-04_r8, -0.2166E-04_r8,  &
         -0.2955E-04_r8, -0.4014E-04_r8, -0.5434E-04_r8, -0.7343E-04_r8, -0.9931E-04_r8,  &
         -0.1346E-03_r8, -0.1826E-03_r8, -0.2467E-03_r8, -0.3283E-03_r8, -0.4246E-03_r8,  &
         -0.5264E-03_r8, -0.6211E-03_r8, -0.6970E-03_r8, -0.7402E-03_r8, -0.7316E-03_r8,  &
         -0.6486E-03_r8,&
                                ! DATA ((h13(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=12_r8,12)_r8,&                               &
         0.5078E-09_r8,  0.8244E-09_r8,  0.1255E-08_r8,  0.1826E-08_r8,  0.2550E-08_r8,   &
         0.3438E-08_r8,  0.4532E-08_r8,  0.5949E-08_r8,  0.8041E-08_r8,  0.1110E-07_r8,   &
         0.1534E-07_r8,  0.2157E-07_r8,  0.3116E-07_r8,  0.4570E-07_r8,  0.6747E-07_r8,   &
         0.9961E-07_r8,  0.1451E-06_r8,  0.2061E-06_r8,  0.2843E-06_r8,  0.3855E-06_r8,   &
         0.5213E-06_r8,  0.7060E-06_r8,  0.9544E-06_r8,  0.1280E-05_r8,  0.1684E-05_r8,   &
         0.2148E-05_r8,  0.2609E-05_r8,  0.3002E-05_r8,  0.3349E-05_r8,  0.3670E-05_r8,   &
         0.3780E-05_r8,&
                                ! DATA ((h11(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=13_r8,13)_r8,&                               &
         0.7529E-04_r8,  0.1306E-03_r8,  0.2160E-03_r8,  0.3410E-03_r8,  0.5175E-03_r8,   &
         0.7626E-03_r8,  0.1103E-02_r8,  0.1577E-02_r8,  0.2246E-02_r8,  0.3196E-02_r8,   &
         0.4544E-02_r8,  0.6455E-02_r8,  0.9152E-02_r8,  0.1294E-01_r8,  0.1823E-01_r8,   &
         0.2560E-01_r8,  0.3577E-01_r8,  0.4971E-01_r8,  0.6870E-01_r8,  0.9443E-01_r8,   &
         0.1290E+00_r8,  0.1748E+00_r8,  0.2340E+00_r8,  0.3078E+00_r8,  0.3956E+00_r8,   &
         0.4942E+00_r8,  0.5969E+00_r8,  0.6970E+00_r8,  0.7896E+00_r8,  0.8714E+00_r8,   &
         0.9364E+00_r8,&
                                ! DATA ((h12(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=13_r8,13)_r8,&                               &
         -0.1461E-06_r8, -0.2663E-06_r8, -0.4512E-06_r8, -0.7027E-06_r8, -0.1014E-05_r8,  &
         -0.1387E-05_r8, -0.1851E-05_r8, -0.2478E-05_r8, -0.3373E-05_r8, -0.4682E-05_r8,  &
         -0.6588E-05_r8, -0.9311E-05_r8, -0.1311E-04_r8, -0.1834E-04_r8, -0.2544E-04_r8,  &
         -0.3502E-04_r8, -0.4789E-04_r8, -0.6515E-04_r8, -0.8846E-04_r8, -0.1202E-03_r8,  &
         -0.1635E-03_r8, -0.2217E-03_r8, -0.2975E-03_r8, -0.3897E-03_r8, -0.4913E-03_r8,  &
         -0.5902E-03_r8, -0.6740E-03_r8, -0.7302E-03_r8, -0.7415E-03_r8, -0.6858E-03_r8,  &
         -0.5447E-03_r8,&
                                ! DATA ((h13(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=13_r8,13)_r8,&                               &
         0.5236E-09_r8,  0.8873E-09_r8,  0.1426E-08_r8,  0.2193E-08_r8,  0.3230E-08_r8,   &
         0.4555E-08_r8,  0.6200E-08_r8,  0.8298E-08_r8,  0.1126E-07_r8,  0.1544E-07_r8,   &
         0.2130E-07_r8,  0.2978E-07_r8,  0.4239E-07_r8,  0.6096E-07_r8,  0.8829E-07_r8,   &
         0.1280E-06_r8,  0.1830E-06_r8,  0.2555E-06_r8,  0.3493E-06_r8,  0.4740E-06_r8,   &
         0.6431E-06_r8,  0.8701E-06_r8,  0.1169E-05_r8,  0.1547E-05_r8,  0.1992E-05_r8,   &
         0.2460E-05_r8,  0.2877E-05_r8,  0.3230E-05_r8,  0.3569E-05_r8,  0.3782E-05_r8,   &
         0.3591E-05_r8,&
                                ! DATA ((h11(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=14_r8,14)_r8,&                               &
         0.7911E-04_r8,  0.1403E-03_r8,  0.2378E-03_r8,  0.3847E-03_r8,  0.5973E-03_r8,   &
         0.8978E-03_r8,  0.1319E-02_r8,  0.1909E-02_r8,  0.2741E-02_r8,  0.3923E-02_r8,   &
         0.5600E-02_r8,  0.7977E-02_r8,  0.1133E-01_r8,  0.1604E-01_r8,  0.2262E-01_r8,   &
         0.3174E-01_r8,  0.4429E-01_r8,  0.6143E-01_r8,  0.8469E-01_r8,  0.1161E+00_r8,   &
         0.1579E+00_r8,  0.2124E+00_r8,  0.2813E+00_r8,  0.3647E+00_r8,  0.4603E+00_r8,   &
         0.5625E+00_r8,  0.6641E+00_r8,  0.7596E+00_r8,  0.8456E+00_r8,  0.9170E+00_r8,   &
         0.9670E+00_r8,&
                                ! DATA ((h12(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=14_r8,14)_r8,&                               &
         -0.1402E-06_r8, -0.2569E-06_r8, -0.4428E-06_r8, -0.7076E-06_r8, -0.1051E-05_r8,  &
         -0.1478E-05_r8, -0.2019E-05_r8, -0.2752E-05_r8, -0.3802E-05_r8, -0.5343E-05_r8,  &
         -0.7594E-05_r8, -0.1082E-04_r8, -0.1536E-04_r8, -0.2166E-04_r8, -0.3028E-04_r8,  &
         -0.4195E-04_r8, -0.5761E-04_r8, -0.7867E-04_r8, -0.1072E-03_r8, -0.1462E-03_r8,  &
         -0.1990E-03_r8, -0.2687E-03_r8, -0.3559E-03_r8, -0.4558E-03_r8, -0.5572E-03_r8,  &
         -0.6476E-03_r8, -0.7150E-03_r8, -0.7439E-03_r8, -0.7133E-03_r8, -0.6015E-03_r8,  &
         -0.4089E-03_r8,&
                                ! DATA ((h13(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=14_r8,14)_r8,&                               &
         0.5531E-09_r8,  0.9757E-09_r8,  0.1644E-08_r8,  0.2650E-08_r8,  0.4074E-08_r8,   &
         0.5957E-08_r8,  0.8314E-08_r8,  0.1128E-07_r8,  0.1528E-07_r8,  0.2087E-07_r8,   &
         0.2874E-07_r8,  0.4002E-07_r8,  0.5631E-07_r8,  0.7981E-07_r8,  0.1139E-06_r8,   &
         0.1621E-06_r8,  0.2275E-06_r8,  0.3136E-06_r8,  0.4280E-06_r8,  0.5829E-06_r8,   &
         0.7917E-06_r8,  0.1067E-05_r8,  0.1419E-05_r8,  0.1844E-05_r8,  0.2310E-05_r8,   &
         0.2747E-05_r8,  0.3113E-05_r8,  0.3455E-05_r8,  0.3739E-05_r8,  0.3715E-05_r8,   &
         0.3125E-05_r8,&
                                ! DATA ((h11(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=15_r8,15)_r8,&                               &
         0.8292E-04_r8,  0.1503E-03_r8,  0.2614E-03_r8,  0.4345E-03_r8,  0.6915E-03_r8,   &
         0.1061E-02_r8,  0.1584E-02_r8,  0.2320E-02_r8,  0.3359E-02_r8,  0.4831E-02_r8,   &
         0.6920E-02_r8,  0.9877E-02_r8,  0.1405E-01_r8,  0.1990E-01_r8,  0.2805E-01_r8,   &
         0.3933E-01_r8,  0.5477E-01_r8,  0.7579E-01_r8,  0.1042E+00_r8,  0.1423E+00_r8,   &
         0.1924E+00_r8,  0.2564E+00_r8,  0.3351E+00_r8,  0.4271E+00_r8,  0.5280E+00_r8,   &
         0.6306E+00_r8,  0.7286E+00_r8,  0.8182E+00_r8,  0.8952E+00_r8,  0.9530E+00_r8,   &
         0.9864E+00_r8,&
                                ! DATA ((h12(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=15_r8,15)_r8,&                               &
         -0.1378E-06_r8, -0.2542E-06_r8, -0.4461E-06_r8, -0.7333E-06_r8, -0.1125E-05_r8,  &
         -0.1630E-05_r8, -0.2281E-05_r8, -0.3159E-05_r8, -0.4410E-05_r8, -0.6246E-05_r8,  &
         -0.8933E-05_r8, -0.1280E-04_r8, -0.1826E-04_r8, -0.2589E-04_r8, -0.3639E-04_r8,  &
         -0.5059E-04_r8, -0.6970E-04_r8, -0.9552E-04_r8, -0.1307E-03_r8, -0.1784E-03_r8,  &
         -0.2422E-03_r8, -0.3237E-03_r8, -0.4203E-03_r8, -0.5227E-03_r8, -0.6184E-03_r8,  &
         -0.6953E-03_r8, -0.7395E-03_r8, -0.7315E-03_r8, -0.6487E-03_r8, -0.4799E-03_r8,  &
         -0.2625E-03_r8,&
                                ! DATA ((h13(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=15_r8,15)_r8,&                               &
         0.5891E-09_r8,  0.1074E-08_r8,  0.1885E-08_r8,  0.3167E-08_r8,  0.5051E-08_r8,   &
         0.7631E-08_r8,  0.1092E-07_r8,  0.1500E-07_r8,  0.2032E-07_r8,  0.2769E-07_r8,   &
         0.3810E-07_r8,  0.5279E-07_r8,  0.7361E-07_r8,  0.1032E-06_r8,  0.1450E-06_r8,   &
         0.2026E-06_r8,  0.2798E-06_r8,  0.3832E-06_r8,  0.5242E-06_r8,  0.7159E-06_r8,   &
         0.9706E-06_r8,  0.1299E-05_r8,  0.1701E-05_r8,  0.2159E-05_r8,  0.2612E-05_r8,   &
         0.2998E-05_r8,  0.3341E-05_r8,  0.3661E-05_r8,  0.3775E-05_r8,  0.3393E-05_r8,   &
         0.2384E-05_r8,&
                                ! DATA ((h11(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=16_r8,16)_r8,&                               &
         0.8637E-04_r8,  0.1598E-03_r8,  0.2853E-03_r8,  0.4876E-03_r8,  0.7970E-03_r8,   &
         0.1251E-02_r8,  0.1901E-02_r8,  0.2820E-02_r8,  0.4118E-02_r8,  0.5955E-02_r8,   &
         0.8557E-02_r8,  0.1224E-01_r8,  0.1742E-01_r8,  0.2467E-01_r8,  0.3476E-01_r8,   &
         0.4864E-01_r8,  0.6759E-01_r8,  0.9332E-01_r8,  0.1280E+00_r8,  0.1738E+00_r8,   &
         0.2330E+00_r8,  0.3069E+00_r8,  0.3949E+00_r8,  0.4935E+00_r8,  0.5964E+00_r8,   &
         0.6965E+00_r8,  0.7893E+00_r8,  0.8713E+00_r8,  0.9363E+00_r8,  0.9780E+00_r8,   &
         0.9961E+00_r8,&
                                ! DATA ((h12(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=16_r8,16)_r8,&                               &
         -0.1383E-06_r8, -0.2577E-06_r8, -0.4608E-06_r8, -0.7793E-06_r8, -0.1237E-05_r8,  &
         -0.1850E-05_r8, -0.2652E-05_r8, -0.3728E-05_r8, -0.5244E-05_r8, -0.7451E-05_r8,  &
         -0.1067E-04_r8, -0.1532E-04_r8, -0.2193E-04_r8, -0.3119E-04_r8, -0.4395E-04_r8,  &
         -0.6126E-04_r8, -0.8466E-04_r8, -0.1164E-03_r8, -0.1596E-03_r8, -0.2177E-03_r8,  &
         -0.2933E-03_r8, -0.3855E-03_r8, -0.4874E-03_r8, -0.5870E-03_r8, -0.6718E-03_r8,  &
         -0.7290E-03_r8, -0.7411E-03_r8, -0.6859E-03_r8, -0.5450E-03_r8, -0.3353E-03_r8,  &
         -0.1363E-03_r8,&
                                ! DATA ((h13(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=16_r8,16)_r8,&                               &
         0.6217E-09_r8,  0.1165E-08_r8,  0.2116E-08_r8,  0.3685E-08_r8,  0.6101E-08_r8,   &
         0.9523E-08_r8,  0.1400E-07_r8,  0.1959E-07_r8,  0.2668E-07_r8,  0.3629E-07_r8,   &
         0.4982E-07_r8,  0.6876E-07_r8,  0.9523E-07_r8,  0.1321E-06_r8,  0.1825E-06_r8,   &
         0.2505E-06_r8,  0.3420E-06_r8,  0.4677E-06_r8,  0.6416E-06_r8,  0.8760E-06_r8,   &
         0.1183E-05_r8,  0.1565E-05_r8,  0.2010E-05_r8,  0.2472E-05_r8,  0.2882E-05_r8,   &
         0.3229E-05_r8,  0.3564E-05_r8,  0.3777E-05_r8,  0.3589E-05_r8,  0.2786E-05_r8,   &
         0.1487E-05_r8,&
                                ! DATA ((h11(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=17_r8,17)_r8,&                               &
         0.8925E-04_r8,  0.1682E-03_r8,  0.3074E-03_r8,  0.5403E-03_r8,  0.9085E-03_r8,   &
         0.1463E-02_r8,  0.2268E-02_r8,  0.3416E-02_r8,  0.5040E-02_r8,  0.7333E-02_r8,   &
         0.1057E-01_r8,  0.1515E-01_r8,  0.2157E-01_r8,  0.3055E-01_r8,  0.4297E-01_r8,   &
         0.6001E-01_r8,  0.8323E-01_r8,  0.1146E+00_r8,  0.1565E+00_r8,  0.2111E+00_r8,   &
         0.2801E+00_r8,  0.3636E+00_r8,  0.4594E+00_r8,  0.5618E+00_r8,  0.6636E+00_r8,   &
         0.7592E+00_r8,  0.8454E+00_r8,  0.9169E+00_r8,  0.9669E+00_r8,  0.9923E+00_r8,   &
         0.9995E+00_r8,&
                                ! DATA ((h12(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=17_r8,17)_r8,&                               &
         -0.1405E-06_r8, -0.2649E-06_r8, -0.4829E-06_r8, -0.8398E-06_r8, -0.1379E-05_r8,  &
         -0.2132E-05_r8, -0.3138E-05_r8, -0.4487E-05_r8, -0.6353E-05_r8, -0.9026E-05_r8,  &
         -0.1290E-04_r8, -0.1851E-04_r8, -0.2650E-04_r8, -0.3772E-04_r8, -0.5319E-04_r8,  &
         -0.7431E-04_r8, -0.1031E-03_r8, -0.1422E-03_r8, -0.1951E-03_r8, -0.2648E-03_r8,  &
         -0.3519E-03_r8, -0.4518E-03_r8, -0.5537E-03_r8, -0.6449E-03_r8, -0.7133E-03_r8,  &
         -0.7432E-03_r8, -0.7133E-03_r8, -0.6018E-03_r8, -0.4092E-03_r8, -0.1951E-03_r8,  &
         -0.5345E-04_r8,&
                                ! DATA ((h13(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=17_r8,17)_r8,&                               &
         0.6457E-09_r8,  0.1235E-08_r8,  0.2303E-08_r8,  0.4149E-08_r8,  0.7120E-08_r8,   &
         0.1152E-07_r8,  0.1749E-07_r8,  0.2508E-07_r8,  0.3462E-07_r8,  0.4718E-07_r8,   &
         0.6452E-07_r8,  0.8874E-07_r8,  0.1222E-06_r8,  0.1675E-06_r8,  0.2276E-06_r8,   &
         0.3076E-06_r8,  0.4174E-06_r8,  0.5714E-06_r8,  0.7837E-06_r8,  0.1067E-05_r8,   &
         0.1428E-05_r8,  0.1859E-05_r8,  0.2327E-05_r8,  0.2760E-05_r8,  0.3122E-05_r8,   &
         0.3458E-05_r8,  0.3739E-05_r8,  0.3715E-05_r8,  0.3126E-05_r8,  0.1942E-05_r8,   &
         0.6977E-06_r8,&
                                ! DATA ((h11(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=18_r8,18)_r8,&                               &
         0.9151E-04_r8,  0.1750E-03_r8,  0.3265E-03_r8,  0.5891E-03_r8,  0.1020E-02_r8,   &
         0.1688E-02_r8,  0.2679E-02_r8,  0.4109E-02_r8,  0.6138E-02_r8,  0.9002E-02_r8,   &
         0.1304E-01_r8,  0.1871E-01_r8,  0.2667E-01_r8,  0.3773E-01_r8,  0.5299E-01_r8,   &
         0.7386E-01_r8,  0.1022E+00_r8,  0.1403E+00_r8,  0.1905E+00_r8,  0.2546E+00_r8,   &
         0.3335E+00_r8,  0.4258E+00_r8,  0.5269E+00_r8,  0.6297E+00_r8,  0.7280E+00_r8,   &
         0.8178E+00_r8,  0.8950E+00_r8,  0.9529E+00_r8,  0.9864E+00_r8,  0.9983E+00_r8,   &
         0.1000E+01_r8,&
                                ! DATA ((h12(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=18_r8,18)_r8,&                               &
         -0.1431E-06_r8, -0.2731E-06_r8, -0.5072E-06_r8, -0.9057E-06_r8, -0.1537E-05_r8,  &
         -0.2460E-05_r8, -0.3733E-05_r8, -0.5449E-05_r8, -0.7786E-05_r8, -0.1106E-04_r8,  &
         -0.1574E-04_r8, -0.2249E-04_r8, -0.3212E-04_r8, -0.4564E-04_r8, -0.6438E-04_r8,  &
         -0.9019E-04_r8, -0.1256E-03_r8, -0.1737E-03_r8, -0.2378E-03_r8, -0.3196E-03_r8,  &
         -0.4163E-03_r8, -0.5191E-03_r8, -0.6154E-03_r8, -0.6931E-03_r8, -0.7384E-03_r8,  &
         -0.7313E-03_r8, -0.6492E-03_r8, -0.4805E-03_r8, -0.2629E-03_r8, -0.8897E-04_r8,  &
         -0.1432E-04_r8,&
                                ! DATA ((h13(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=18_r8,18)_r8,&                               &
         0.6607E-09_r8,  0.1282E-08_r8,  0.2441E-08_r8,  0.4522E-08_r8,  0.8027E-08_r8,   &
         0.1348E-07_r8,  0.2122E-07_r8,  0.3139E-07_r8,  0.4435E-07_r8,  0.6095E-07_r8,   &
         0.8319E-07_r8,  0.1139E-06_r8,  0.1557E-06_r8,  0.2107E-06_r8,  0.2819E-06_r8,   &
         0.3773E-06_r8,  0.5107E-06_r8,  0.6982E-06_r8,  0.9542E-06_r8,  0.1290E-05_r8,   &
         0.1703E-05_r8,  0.2170E-05_r8,  0.2628E-05_r8,  0.3013E-05_r8,  0.3352E-05_r8,   &
         0.3669E-05_r8,  0.3780E-05_r8,  0.3397E-05_r8,  0.2386E-05_r8,  0.1062E-05_r8,   &
         0.2216E-06_r8,&
                                ! DATA ((h11(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=19_r8,19)_r8,&                               &
         0.9320E-04_r8,  0.1803E-03_r8,  0.3422E-03_r8,  0.6317E-03_r8,  0.1124E-02_r8,   &
         0.1915E-02_r8,  0.3121E-02_r8,  0.4890E-02_r8,  0.7421E-02_r8,  0.1100E-01_r8,   &
         0.1602E-01_r8,  0.2306E-01_r8,  0.3288E-01_r8,  0.4647E-01_r8,  0.6515E-01_r8,   &
         0.9066E-01_r8,  0.1252E+00_r8,  0.1710E+00_r8,  0.2304E+00_r8,  0.3045E+00_r8,   &
         0.3928E+00_r8,  0.4918E+00_r8,  0.5951E+00_r8,  0.6956E+00_r8,  0.7887E+00_r8,   &
         0.8709E+00_r8,  0.9361E+00_r8,  0.9780E+00_r8,  0.9961E+00_r8,  0.9999E+00_r8,   &
         0.1000E+01_r8,&
                                ! DATA ((h12(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=19_r8,19)_r8,&                               &
         -0.1454E-06_r8, -0.2805E-06_r8, -0.5296E-06_r8, -0.9685E-06_r8, -0.1695E-05_r8,  &
         -0.2812E-05_r8, -0.4412E-05_r8, -0.6606E-05_r8, -0.9573E-05_r8, -0.1363E-04_r8,  &
         -0.1932E-04_r8, -0.2743E-04_r8, -0.3897E-04_r8, -0.5520E-04_r8, -0.7787E-04_r8,  &
         -0.1094E-03_r8, -0.1529E-03_r8, -0.2117E-03_r8, -0.2880E-03_r8, -0.3809E-03_r8,  &
         -0.4834E-03_r8, -0.5836E-03_r8, -0.6692E-03_r8, -0.7275E-03_r8, -0.7408E-03_r8,  &
         -0.6865E-03_r8, -0.5459E-03_r8, -0.3360E-03_r8, -0.1365E-03_r8, -0.2935E-04_r8,  &
         -0.2173E-05_r8,&
                                ! DATA ((h13(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=19_r8,19)_r8,&                               &
         0.6693E-09_r8,  0.1312E-08_r8,  0.2538E-08_r8,  0.4802E-08_r8,  0.8778E-08_r8,   &
         0.1528E-07_r8,  0.2501E-07_r8,  0.3836E-07_r8,  0.5578E-07_r8,  0.7806E-07_r8,   &
         0.1069E-06_r8,  0.1456E-06_r8,  0.1970E-06_r8,  0.2631E-06_r8,  0.3485E-06_r8,   &
         0.4642E-06_r8,  0.6268E-06_r8,  0.8526E-06_r8,  0.1157E-05_r8,  0.1545E-05_r8,   &
         0.2002E-05_r8,  0.2478E-05_r8,  0.2897E-05_r8,  0.3245E-05_r8,  0.3578E-05_r8,   &
         0.3789E-05_r8,  0.3598E-05_r8,  0.2792E-05_r8,  0.1489E-05_r8,  0.4160E-06_r8,   &
         0.3843E-07_r8,&
                                ! DATA ((h11(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=20_r8,20)_r8,&                               &
         0.9443E-04_r8,  0.1843E-03_r8,  0.3545E-03_r8,  0.6669E-03_r8,  0.1217E-02_r8,   &
         0.2133E-02_r8,  0.3577E-02_r8,  0.5742E-02_r8,  0.8880E-02_r8,  0.1333E-01_r8,   &
         0.1958E-01_r8,  0.2830E-01_r8,  0.4039E-01_r8,  0.5705E-01_r8,  0.7988E-01_r8,   &
         0.1110E+00_r8,  0.1526E+00_r8,  0.2072E+00_r8,  0.2764E+00_r8,  0.3604E+00_r8,   &
         0.4567E+00_r8,  0.5597E+00_r8,  0.6620E+00_r8,  0.7581E+00_r8,  0.8447E+00_r8,   &
         0.9165E+00_r8,  0.9668E+00_r8,  0.9923E+00_r8,  0.9995E+00_r8,  0.1000E+01_r8,   &
         0.1000E+01_r8,&
                                ! DATA ((h12(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=20_r8,20)_r8,&                               &
         -0.1472E-06_r8, -0.2866E-06_r8, -0.5485E-06_r8, -0.1024E-05_r8, -0.1842E-05_r8,  &
         -0.3160E-05_r8, -0.5136E-05_r8, -0.7922E-05_r8, -0.1171E-04_r8, -0.1682E-04_r8,  &
         -0.2381E-04_r8, -0.3355E-04_r8, -0.4729E-04_r8, -0.6673E-04_r8, -0.9417E-04_r8,  &
         -0.1327E-03_r8, -0.1858E-03_r8, -0.2564E-03_r8, -0.3449E-03_r8, -0.4463E-03_r8,  &
         -0.5495E-03_r8, -0.6420E-03_r8, -0.7116E-03_r8, -0.7427E-03_r8, -0.7139E-03_r8,  &
         -0.6031E-03_r8, -0.4104E-03_r8, -0.1957E-03_r8, -0.5358E-04_r8, -0.6176E-05_r8,  &
         -0.1347E-06_r8,&
                                ! DATA ((h13(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=20_r8,20)_r8,&                               &
         0.6750E-09_r8,  0.1332E-08_r8,  0.2602E-08_r8,  0.5003E-08_r8,  0.9367E-08_r8,   &
         0.1684E-07_r8,  0.2863E-07_r8,  0.4566E-07_r8,  0.6865E-07_r8,  0.9861E-07_r8,   &
         0.1368E-06_r8,  0.1856E-06_r8,  0.2479E-06_r8,  0.3274E-06_r8,  0.4315E-06_r8,   &
         0.5739E-06_r8,  0.7710E-06_r8,  0.1040E-05_r8,  0.1394E-05_r8,  0.1829E-05_r8,   &
         0.2309E-05_r8,  0.2759E-05_r8,  0.3131E-05_r8,  0.3472E-05_r8,  0.3755E-05_r8,   &
         0.3730E-05_r8,  0.3138E-05_r8,  0.1948E-05_r8,  0.6994E-06_r8,  0.1022E-06_r8,   &
         0.2459E-08_r8,&
                                ! DATA ((h11(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=21_r8,21)_r8,&                               &
         0.9530E-04_r8,  0.1872E-03_r8,  0.3638E-03_r8,  0.6949E-03_r8,  0.1295E-02_r8,   &
         0.2330E-02_r8,  0.4022E-02_r8,  0.6633E-02_r8,  0.1049E-01_r8,  0.1601E-01_r8,   &
         0.2378E-01_r8,  0.3456E-01_r8,  0.4941E-01_r8,  0.6980E-01_r8,  0.9765E-01_r8,   &
         0.1353E+00_r8,  0.1852E+00_r8,  0.2494E+00_r8,  0.3287E+00_r8,  0.4216E+00_r8,   &
         0.5234E+00_r8,  0.6272E+00_r8,  0.7261E+00_r8,  0.8166E+00_r8,  0.8943E+00_r8,   &
         0.9526E+00_r8,  0.9863E+00_r8,  0.9983E+00_r8,  0.1000E+01_r8,  0.1000E+01_r8,   &
         0.1000E+01_r8,&
                                ! DATA ((h12(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=21_r8,21)_r8,&                               &
         -0.1487E-06_r8, -0.2912E-06_r8, -0.5636E-06_r8, -0.1069E-05_r8, -0.1969E-05_r8,  &
         -0.3483E-05_r8, -0.5858E-05_r8, -0.9334E-05_r8, -0.1416E-04_r8, -0.2067E-04_r8,  &
         -0.2936E-04_r8, -0.4113E-04_r8, -0.5750E-04_r8, -0.8072E-04_r8, -0.1139E-03_r8,  &
         -0.1606E-03_r8, -0.2246E-03_r8, -0.3076E-03_r8, -0.4067E-03_r8, -0.5121E-03_r8,  &
         -0.6110E-03_r8, -0.6909E-03_r8, -0.7378E-03_r8, -0.7321E-03_r8, -0.6509E-03_r8,  &
         -0.4825E-03_r8, -0.2641E-03_r8, -0.8936E-04_r8, -0.1436E-04_r8, -0.5966E-06_r8,  &
         0.0000E+00_r8,&
                                ! DATA ((h13(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=21_r8,21)_r8,&                               &
         0.6777E-09_r8,  0.1344E-08_r8,  0.2643E-08_r8,  0.5138E-08_r8,  0.9798E-08_r8,   &
         0.1809E-07_r8,  0.3185E-07_r8,  0.5285E-07_r8,  0.8249E-07_r8,  0.1222E-06_r8,   &
         0.1730E-06_r8,  0.2351E-06_r8,  0.3111E-06_r8,  0.4078E-06_r8,  0.5366E-06_r8,   &
         0.7117E-06_r8,  0.9495E-06_r8,  0.1266E-05_r8,  0.1667E-05_r8,  0.2132E-05_r8,   &
         0.2600E-05_r8,  0.3001E-05_r8,  0.3354E-05_r8,  0.3679E-05_r8,  0.3796E-05_r8,   &
         0.3414E-05_r8,  0.2399E-05_r8,  0.1067E-05_r8,  0.2222E-06_r8,  0.1075E-07_r8,   &
         0.0000E+00_r8,&
                                ! DATA ((h11(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=22_r8,22)_r8,&                               &
         0.9589E-04_r8,  0.1893E-03_r8,  0.3706E-03_r8,  0.7162E-03_r8,  0.1357E-02_r8,   &
         0.2500E-02_r8,  0.4434E-02_r8,  0.7523E-02_r8,  0.1220E-01_r8,  0.1900E-01_r8,   &
         0.2860E-01_r8,  0.4190E-01_r8,  0.6015E-01_r8,  0.8508E-01_r8,  0.1189E+00_r8,   &
         0.1643E+00_r8,  0.2233E+00_r8,  0.2976E+00_r8,  0.3865E+00_r8,  0.4865E+00_r8,   &
         0.5909E+00_r8,  0.6925E+00_r8,  0.7866E+00_r8,  0.8696E+00_r8,  0.9355E+00_r8,   &
         0.9778E+00_r8,  0.9961E+00_r8,  0.9999E+00_r8,  0.1000E+01_r8,  0.1000E+01_r8,   &
         0.1000E+01_r8,&
                                ! DATA ((h12(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=22_r8,22)_r8,&                               &
         -0.1496E-06_r8, -0.2947E-06_r8, -0.5749E-06_r8, -0.1105E-05_r8, -0.2074E-05_r8,  &
         -0.3763E-05_r8, -0.6531E-05_r8, -0.1076E-04_r8, -0.1682E-04_r8, -0.2509E-04_r8,  &
         -0.3605E-04_r8, -0.5049E-04_r8, -0.7012E-04_r8, -0.9787E-04_r8, -0.1378E-03_r8,  &
         -0.1939E-03_r8, -0.2695E-03_r8, -0.3641E-03_r8, -0.4703E-03_r8, -0.5750E-03_r8,  &
         -0.6648E-03_r8, -0.7264E-03_r8, -0.7419E-03_r8, -0.6889E-03_r8, -0.5488E-03_r8,  &
         -0.3382E-03_r8, -0.1375E-03_r8, -0.2951E-04_r8, -0.2174E-05_r8,  0.0000E+00_r8,  &
         0.0000E+00_r8,&
                                ! DATA ((h13(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=22_r8,22)_r8,&                               &
         0.6798E-09_r8,  0.1350E-08_r8,  0.2667E-08_r8,  0.5226E-08_r8,  0.1010E-07_r8,   &
         0.1903E-07_r8,  0.3455E-07_r8,  0.5951E-07_r8,  0.9658E-07_r8,  0.1479E-06_r8,   &
         0.2146E-06_r8,  0.2951E-06_r8,  0.3903E-06_r8,  0.5101E-06_r8,  0.6693E-06_r8,   &
         0.8830E-06_r8,  0.1168E-05_r8,  0.1532E-05_r8,  0.1968E-05_r8,  0.2435E-05_r8,   &
         0.2859E-05_r8,  0.3222E-05_r8,  0.3572E-05_r8,  0.3797E-05_r8,  0.3615E-05_r8,   &
         0.2811E-05_r8,  0.1500E-05_r8,  0.4185E-06_r8,  0.3850E-07_r8,  0.0000E+00_r8,   &
         0.0000E+00_r8,&
                                ! DATA ((h11(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=23_r8,23)_r8,&                               &
         0.9629E-04_r8,  0.1907E-03_r8,  0.3755E-03_r8,  0.7321E-03_r8,  0.1406E-02_r8,   &
         0.2639E-02_r8,  0.4796E-02_r8,  0.8368E-02_r8,  0.1394E-01_r8,  0.2221E-01_r8,   &
         0.3400E-01_r8,  0.5037E-01_r8,  0.7279E-01_r8,  0.1033E+00_r8,  0.1442E+00_r8,   &
         0.1983E+00_r8,  0.2673E+00_r8,  0.3516E+00_r8,  0.4489E+00_r8,  0.5533E+00_r8,   &
         0.6572E+00_r8,  0.7547E+00_r8,  0.8425E+00_r8,  0.9153E+00_r8,  0.9663E+00_r8,   &
         0.9922E+00_r8,  0.9995E+00_r8,  0.1000E+01_r8,  0.1000E+01_r8,  0.1000E+01_r8,   &
         0.1000E+01_r8,&
                                ! DATA ((h12(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=23_r8,23)_r8,&                               &
         -0.1503E-06_r8, -0.2971E-06_r8, -0.5832E-06_r8, -0.1131E-05_r8, -0.2154E-05_r8,  &
         -0.3992E-05_r8, -0.7122E-05_r8, -0.1211E-04_r8, -0.1954E-04_r8, -0.2995E-04_r8,  &
         -0.4380E-04_r8, -0.6183E-04_r8, -0.8577E-04_r8, -0.1191E-03_r8, -0.1668E-03_r8,  &
         -0.2333E-03_r8, -0.3203E-03_r8, -0.4237E-03_r8, -0.5324E-03_r8, -0.6318E-03_r8,  &
         -0.7075E-03_r8, -0.7429E-03_r8, -0.7168E-03_r8, -0.6071E-03_r8, -0.4139E-03_r8,  &
         -0.1976E-03_r8, -0.5410E-04_r8, -0.6215E-05_r8, -0.1343E-06_r8,  0.0000E+00_r8,  &
         0.0000E+00_r8,&
                                ! DATA ((h13(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=23_r8,23)_r8,&                               &
         0.6809E-09_r8,  0.1356E-08_r8,  0.2683E-08_r8,  0.5287E-08_r8,  0.1030E-07_r8,   &
         0.1971E-07_r8,  0.3665E-07_r8,  0.6528E-07_r8,  0.1100E-06_r8,  0.1744E-06_r8,   &
         0.2599E-06_r8,  0.3650E-06_r8,  0.4887E-06_r8,  0.6398E-06_r8,  0.8358E-06_r8,   &
         0.1095E-05_r8,  0.1429E-05_r8,  0.1836E-05_r8,  0.2286E-05_r8,  0.2716E-05_r8,   &
         0.3088E-05_r8,  0.3444E-05_r8,  0.3748E-05_r8,  0.3740E-05_r8,  0.3157E-05_r8,   &
         0.1966E-05_r8,  0.7064E-06_r8,  0.1030E-06_r8,  0.2456E-08_r8,  0.0000E+00_r8,   &
         0.0000E+00_r8,&
                                ! DATA ((h11(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=24_r8,24)_r8,&                               &
         0.9655E-04_r8,  0.1917E-03_r8,  0.3789E-03_r8,  0.7436E-03_r8,  0.1442E-02_r8,   &
         0.2748E-02_r8,  0.5100E-02_r8,  0.9128E-02_r8,  0.1563E-01_r8,  0.2553E-01_r8,   &
         0.3987E-01_r8,  0.5994E-01_r8,  0.8746E-01_r8,  0.1246E+00_r8,  0.1739E+00_r8,   &
         0.2376E+00_r8,  0.3170E+00_r8,  0.4107E+00_r8,  0.5141E+00_r8,  0.6198E+00_r8,   &
         0.7208E+00_r8,  0.8130E+00_r8,  0.8923E+00_r8,  0.9517E+00_r8,  0.9860E+00_r8,   &
         0.9983E+00_r8,  0.1000E+01_r8,  0.1000E+01_r8,  0.1000E+01_r8,  0.1000E+01_r8,   &
         0.1000E+01_r8,&
                                ! DATA ((h12(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=24_r8,24)_r8,&                               &
         -0.1508E-06_r8, -0.2989E-06_r8, -0.5892E-06_r8, -0.1151E-05_r8, -0.2216E-05_r8,  &
         -0.4175E-05_r8, -0.7619E-05_r8, -0.1333E-04_r8, -0.2217E-04_r8, -0.3497E-04_r8,  &
         -0.5238E-04_r8, -0.7513E-04_r8, -0.1049E-03_r8, -0.1455E-03_r8, -0.2021E-03_r8,  &
         -0.2790E-03_r8, -0.3757E-03_r8, -0.4839E-03_r8, -0.5902E-03_r8, -0.6794E-03_r8,  &
         -0.7344E-03_r8, -0.7341E-03_r8, -0.6557E-03_r8, -0.4874E-03_r8, -0.2674E-03_r8,  &
         -0.9059E-04_r8, -0.1455E-04_r8, -0.5986E-06_r8,  0.0000E+00_r8,  0.0000E+00_r8,  &
         0.0000E+00_r8,&
                                ! DATA ((h13(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=24_r8,24)_r8,&                               &
         0.6812E-09_r8,  0.1356E-08_r8,  0.2693E-08_r8,  0.5328E-08_r8,  0.1045E-07_r8,   &
         0.2021E-07_r8,  0.3826E-07_r8,  0.6994E-07_r8,  0.1218E-06_r8,  0.1997E-06_r8,   &
         0.3069E-06_r8,  0.4428E-06_r8,  0.6064E-06_r8,  0.8015E-06_r8,  0.1043E-05_r8,   &
         0.1351E-05_r8,  0.1733E-05_r8,  0.2168E-05_r8,  0.2598E-05_r8,  0.2968E-05_r8,   &
         0.3316E-05_r8,  0.3662E-05_r8,  0.3801E-05_r8,  0.3433E-05_r8,  0.2422E-05_r8,   &
         0.1081E-05_r8,  0.2256E-06_r8,  0.1082E-07_r8,  0.0000E+00_r8,  0.0000E+00_r8,   &
         0.0000E+00_r8,&
                                ! DATA ((h11(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=25_r8,25)_r8,&                               &
         0.9672E-04_r8,  0.1923E-03_r8,  0.3813E-03_r8,  0.7518E-03_r8,  0.1469E-02_r8,   &
         0.2832E-02_r8,  0.5343E-02_r8,  0.9779E-02_r8,  0.1719E-01_r8,  0.2882E-01_r8,   &
         0.4606E-01_r8,  0.7051E-01_r8,  0.1042E+00_r8,  0.1493E+00_r8,  0.2081E+00_r8,   &
         0.2824E+00_r8,  0.3720E+00_r8,  0.4736E+00_r8,  0.5803E+00_r8,  0.6844E+00_r8,   &
         0.7810E+00_r8,  0.8663E+00_r8,  0.9339E+00_r8,  0.9772E+00_r8,  0.9960E+00_r8,   &
         0.9999E+00_r8,  0.1000E+01_r8,  0.1000E+01_r8,  0.1000E+01_r8,  0.1000E+01_r8,   &
         0.1000E+01_r8,&
                                ! DATA ((h12(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=25_r8,25)_r8,&                               &
         -0.1511E-06_r8, -0.3001E-06_r8, -0.5934E-06_r8, -0.1166E-05_r8, -0.2263E-05_r8,  &
         -0.4319E-05_r8, -0.8028E-05_r8, -0.1438E-04_r8, -0.2460E-04_r8, -0.3991E-04_r8,  &
         -0.6138E-04_r8, -0.9005E-04_r8, -0.1278E-03_r8, -0.1778E-03_r8, -0.2447E-03_r8,  &
         -0.3313E-03_r8, -0.4342E-03_r8, -0.5424E-03_r8, -0.6416E-03_r8, -0.7146E-03_r8,  &
         -0.7399E-03_r8, -0.6932E-03_r8, -0.5551E-03_r8, -0.3432E-03_r8, -0.1398E-03_r8,  &
         -0.3010E-04_r8, -0.2229E-05_r8,  0.0000E+00_r8,  0.0000E+00_r8,  0.0000E+00_r8,  &
         0.0000E+00_r8,&
                                ! DATA ((h13(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=25_r8,25)_r8,&                               &
         0.6815E-09_r8,  0.1358E-08_r8,  0.2698E-08_r8,  0.5355E-08_r8,  0.1054E-07_r8,   &
         0.2056E-07_r8,  0.3942E-07_r8,  0.7349E-07_r8,  0.1315E-06_r8,  0.2226E-06_r8,   &
         0.3537E-06_r8,  0.5266E-06_r8,  0.7407E-06_r8,  0.9958E-06_r8,  0.1296E-05_r8,   &
         0.1657E-05_r8,  0.2077E-05_r8,  0.2512E-05_r8,  0.2893E-05_r8,  0.3216E-05_r8,   &
         0.3562E-05_r8,  0.3811E-05_r8,  0.3644E-05_r8,  0.2841E-05_r8,  0.1524E-05_r8,   &
         0.4276E-06_r8,  0.3960E-07_r8,  0.0000E+00_r8,  0.0000E+00_r8,  0.0000E+00_r8,   &
         0.0000E+00_r8,&
                                ! DATA ((h11(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=26_r8,26)_r8,&                               &
         0.9682E-04_r8,  0.1928E-03_r8,  0.3829E-03_r8,  0.7576E-03_r8,  0.1489E-02_r8,   &
         0.2894E-02_r8,  0.5533E-02_r8,  0.1031E-01_r8,  0.1856E-01_r8,  0.3195E-01_r8,   &
         0.5237E-01_r8,  0.8187E-01_r8,  0.1227E+00_r8,  0.1771E+00_r8,  0.2468E+00_r8,   &
         0.3321E+00_r8,  0.4312E+00_r8,  0.5384E+00_r8,  0.6455E+00_r8,  0.7463E+00_r8,   &
         0.8372E+00_r8,  0.9126E+00_r8,  0.9652E+00_r8,  0.9919E+00_r8,  0.9994E+00_r8,   &
         0.1000E+01_r8,  0.1000E+01_r8,  0.1000E+01_r8,  0.1000E+01_r8,  0.1000E+01_r8,   &
         0.1000E+01_r8,&
                                ! DATA ((h12(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=26_r8,26)_r8,&                               &
         -0.1513E-06_r8, -0.3009E-06_r8, -0.5966E-06_r8, -0.1176E-05_r8, -0.2299E-05_r8,  &
         -0.4430E-05_r8, -0.8352E-05_r8, -0.1526E-04_r8, -0.2674E-04_r8, -0.4454E-04_r8,  &
         -0.7042E-04_r8, -0.1062E-03_r8, -0.1540E-03_r8, -0.2163E-03_r8, -0.2951E-03_r8,  &
         -0.3899E-03_r8, -0.4948E-03_r8, -0.5983E-03_r8, -0.6846E-03_r8, -0.7332E-03_r8,  &
         -0.7182E-03_r8, -0.6142E-03_r8, -0.4209E-03_r8, -0.2014E-03_r8, -0.5530E-04_r8,  &
         -0.6418E-05_r8, -0.1439E-06_r8,  0.0000E+00_r8,  0.0000E+00_r8,  0.0000E+00_r8,  &
         0.0000E+00_r8,&
                                ! DATA ((h13(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=26_r8,26)_r8,&                               &
         0.6817E-09_r8,  0.1359E-08_r8,  0.2702E-08_r8,  0.5374E-08_r8,  0.1061E-07_r8,   &
         0.2079E-07_r8,  0.4022E-07_r8,  0.7610E-07_r8,  0.1392E-06_r8,  0.2428E-06_r8,   &
         0.3992E-06_r8,  0.6149E-06_r8,  0.8893E-06_r8,  0.1220E-05_r8,  0.1599E-05_r8,   &
         0.2015E-05_r8,  0.2453E-05_r8,  0.2853E-05_r8,  0.3173E-05_r8,  0.3488E-05_r8,   &
         0.3792E-05_r8,  0.3800E-05_r8,  0.3210E-05_r8,  0.2002E-05_r8,  0.7234E-06_r8,   &
         0.1068E-06_r8,  0.2646E-08_r8,  0.0000E+00_r8,  0.0000E+00_r8,  0.0000E+00_r8,   &
         0.0000E+00_r8/),SHAPE=(/nx*3*nh/))

    it=0
    DO k=1,nx
       DO i=1,3
          DO j=1,nh
             it=it+1
             IF(i==1)THEN
                !WRITE(*,'(a5,2e17.9) ' )'h11   ',data3(it),h11(k,j)
                h11(k,j,1)=data3(it)
             END IF
             IF(i==2) THEN
                !WRITE(*,'(a5,2e17.9) ' )'h12   ',data3(it),h12(k,j)
                h12(k,j,1)=data3(it)
             END IF
             IF(i==3) THEN
                !WRITE(*,'(a5,2e17.9) ' )'h13   ',data3(it),h13(k,j)
                h13(k,j,1) =data3(it)
             END IF
          END DO
       END DO
    END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    data4(1:nx*3*nh)=RESHAPE(SOURCE=(/&
                                ! DATA ((h21(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 1_r8, 1)/                               &
         0.3920E-05_r8,  0.7617E-05_r8,  0.1455E-04_r8,  0.2706E-04_r8,  0.4855E-04_r8,   &
         0.8315E-04_r8,  0.1349E-03_r8,  0.2063E-03_r8,  0.2984E-03_r8,  0.4109E-03_r8,   &
         0.5422E-03_r8,  0.6896E-03_r8,  0.8537E-03_r8,  0.1041E-02_r8,  0.1262E-02_r8,   &
         0.1534E-02_r8,  0.1870E-02_r8,  0.2286E-02_r8,  0.2803E-02_r8,  0.3444E-02_r8,   &
         0.4242E-02_r8,  0.5244E-02_r8,  0.6511E-02_r8,  0.8138E-02_r8,  0.1027E-01_r8,   &
         0.1312E-01_r8,  0.1697E-01_r8,  0.2222E-01_r8,  0.2941E-01_r8,  0.3923E-01_r8,   &
         0.5258E-01_r8, & 
                                ! DATA ((h22(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 1_r8, 1)/                               &
         -0.5622E-07_r8, -0.1071E-06_r8, -0.1983E-06_r8, -0.3533E-06_r8, -0.5991E-06_r8,  &
         -0.9592E-06_r8, -0.1444E-05_r8, -0.2049E-05_r8, -0.2764E-05_r8, -0.3577E-05_r8,  &
         -0.4469E-05_r8, -0.5467E-05_r8, -0.6654E-05_r8, -0.8137E-05_r8, -0.1002E-04_r8,  &
         -0.1237E-04_r8, -0.1528E-04_r8, -0.1884E-04_r8, -0.2310E-04_r8, -0.2809E-04_r8,  &
         -0.3396E-04_r8, -0.4098E-04_r8, -0.4960E-04_r8, -0.6058E-04_r8, -0.7506E-04_r8,  &
         -0.9451E-04_r8, -0.1207E-03_r8, -0.1558E-03_r8, -0.2026E-03_r8, -0.2648E-03_r8,  &
         -0.3468E-03_r8, & 
                                ! DATA ((h23(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 1_r8, 1)_r8, &                                &
         -0.2195E-09_r8, -0.4031E-09_r8, -0.7043E-09_r8, -0.1153E-08_r8, -0.1737E-08_r8,  &
         -0.2395E-08_r8, -0.3020E-08_r8, -0.3549E-08_r8, -0.4034E-08_r8, -0.4421E-08_r8,  &
         -0.4736E-08_r8, -0.5681E-08_r8, -0.8289E-08_r8, -0.1287E-07_r8, -0.1873E-07_r8,  &
         -0.2523E-07_r8, -0.3223E-07_r8, -0.3902E-07_r8, -0.4409E-07_r8, -0.4699E-07_r8,  &
         -0.4782E-07_r8, -0.4705E-07_r8, -0.4657E-07_r8, -0.4885E-07_r8, -0.5550E-07_r8,  &
         -0.6619E-07_r8, -0.7656E-07_r8, -0.8027E-07_r8, -0.7261E-07_r8, -0.4983E-07_r8,  &
         -0.1101E-07_r8, & 
                                ! DATA ((h21(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 2_r8, 2)_r8, &                                &
         0.3920E-05_r8,  0.7617E-05_r8,  0.1455E-04_r8,  0.2706E-04_r8,  0.4856E-04_r8,   &
         0.8318E-04_r8,  0.1349E-03_r8,  0.2065E-03_r8,  0.2986E-03_r8,  0.4114E-03_r8,   &
         0.5431E-03_r8,  0.6912E-03_r8,  0.8566E-03_r8,  0.1046E-02_r8,  0.1272E-02_r8,   &
         0.1552E-02_r8,  0.1902E-02_r8,  0.2342E-02_r8,  0.2899E-02_r8,  0.3605E-02_r8,   &
         0.4501E-02_r8,  0.5648E-02_r8,  0.7129E-02_r8,  0.9063E-02_r8,  0.1163E-01_r8,   &
         0.1509E-01_r8,  0.1980E-01_r8,  0.2626E-01_r8,  0.3510E-01_r8,  0.4717E-01_r8,   &
         0.6351E-01_r8, & 
                                ! DATA ((h22(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 2_r8, 2)_r8, &                                &
         -0.5622E-07_r8, -0.1071E-06_r8, -0.1983E-06_r8, -0.3534E-06_r8, -0.5992E-06_r8,  &
         -0.9594E-06_r8, -0.1445E-05_r8, -0.2050E-05_r8, -0.2766E-05_r8, -0.3580E-05_r8,  &
         -0.4476E-05_r8, -0.5479E-05_r8, -0.6677E-05_r8, -0.8179E-05_r8, -0.1009E-04_r8,  &
         -0.1251E-04_r8, -0.1553E-04_r8, -0.1928E-04_r8, -0.2384E-04_r8, -0.2930E-04_r8,  &
         -0.3588E-04_r8, -0.4393E-04_r8, -0.5403E-04_r8, -0.6714E-04_r8, -0.8458E-04_r8,  &
         -0.1082E-03_r8, -0.1400E-03_r8, -0.1829E-03_r8, -0.2401E-03_r8, -0.3157E-03_r8,  &
         -0.4147E-03_r8, & 
                                ! DATA ((h23(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 2_r8, 2)_r8, &                                &
         -0.2195E-09_r8, -0.4032E-09_r8, -0.7046E-09_r8, -0.1153E-08_r8, -0.1738E-08_r8,  &
         -0.2395E-08_r8, -0.3021E-08_r8, -0.3550E-08_r8, -0.4035E-08_r8, -0.4423E-08_r8,  &
         -0.4740E-08_r8, -0.5692E-08_r8, -0.8314E-08_r8, -0.1292E-07_r8, -0.1882E-07_r8,  &
         -0.2536E-07_r8, -0.3242E-07_r8, -0.3927E-07_r8, -0.4449E-07_r8, -0.4767E-07_r8,  &
         -0.4889E-07_r8, -0.4857E-07_r8, -0.4860E-07_r8, -0.5132E-07_r8, -0.5847E-07_r8,  &
         -0.6968E-07_r8, -0.8037E-07_r8, -0.8400E-07_r8, -0.7521E-07_r8, -0.4830E-07_r8,  &
         -0.7562E-09_r8, & 
                                ! DATA ((h21(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 3_r8, 3)_r8, &                                &
         0.3920E-05_r8,  0.7617E-05_r8,  0.1455E-04_r8,  0.2707E-04_r8,  0.4858E-04_r8,   &
         0.8322E-04_r8,  0.1350E-03_r8,  0.2066E-03_r8,  0.2990E-03_r8,  0.4122E-03_r8,   &
         0.5444E-03_r8,  0.6937E-03_r8,  0.8611E-03_r8,  0.1054E-02_r8,  0.1287E-02_r8,   &
         0.1579E-02_r8,  0.1949E-02_r8,  0.2424E-02_r8,  0.3037E-02_r8,  0.3829E-02_r8,   &
         0.4854E-02_r8,  0.6190E-02_r8,  0.7942E-02_r8,  0.1026E-01_r8,  0.1338E-01_r8,   &
         0.1761E-01_r8,  0.2340E-01_r8,  0.3134E-01_r8,  0.4223E-01_r8,  0.5703E-01_r8,   &
         0.7694E-01_r8, & 
                                ! DATA ((h22(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 3_r8, 3)_r8, &                                &
         -0.5622E-07_r8, -0.1071E-06_r8, -0.1983E-06_r8, -0.3535E-06_r8, -0.5994E-06_r8,  &
         -0.9599E-06_r8, -0.1446E-05_r8, -0.2052E-05_r8, -0.2769E-05_r8, -0.3586E-05_r8,  &
         -0.4487E-05_r8, -0.5499E-05_r8, -0.6712E-05_r8, -0.8244E-05_r8, -0.1021E-04_r8,  &
         -0.1272E-04_r8, -0.1591E-04_r8, -0.1992E-04_r8, -0.2489E-04_r8, -0.3097E-04_r8,  &
         -0.3845E-04_r8, -0.4782E-04_r8, -0.5982E-04_r8, -0.7558E-04_r8, -0.9674E-04_r8,  &
         -0.1254E-03_r8, -0.1644E-03_r8, -0.2167E-03_r8, -0.2863E-03_r8, -0.3777E-03_r8,  &
         -0.4959E-03_r8, & 
                                ! DATA ((h23(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 3_r8, 3)_r8, &                                &
         -0.2196E-09_r8, -0.4033E-09_r8, -0.7048E-09_r8, -0.1154E-08_r8, -0.1739E-08_r8,  &
         -0.2396E-08_r8, -0.3022E-08_r8, -0.3551E-08_r8, -0.4036E-08_r8, -0.4425E-08_r8,  &
         -0.4746E-08_r8, -0.5710E-08_r8, -0.8354E-08_r8, -0.1300E-07_r8, -0.1894E-07_r8,  &
         -0.2554E-07_r8, -0.3265E-07_r8, -0.3958E-07_r8, -0.4502E-07_r8, -0.4859E-07_r8,  &
         -0.5030E-07_r8, -0.5053E-07_r8, -0.5104E-07_r8, -0.5427E-07_r8, -0.6204E-07_r8,  &
         -0.7388E-07_r8, -0.8477E-07_r8, -0.8760E-07_r8, -0.7545E-07_r8, -0.4099E-07_r8,  &
         0.2046E-07_r8, & 
                                ! DATA ((h21(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 4_r8, 4)_r8, &                                &
         0.3919E-05_r8,  0.7618E-05_r8,  0.1455E-04_r8,  0.2708E-04_r8,  0.4860E-04_r8,   &
         0.8327E-04_r8,  0.1352E-03_r8,  0.2070E-03_r8,  0.2997E-03_r8,  0.4134E-03_r8,   &
         0.5466E-03_r8,  0.6976E-03_r8,  0.8682E-03_r8,  0.1067E-02_r8,  0.1310E-02_r8,   &
         0.1619E-02_r8,  0.2020E-02_r8,  0.2543E-02_r8,  0.3232E-02_r8,  0.4136E-02_r8,   &
         0.5328E-02_r8,  0.6906E-02_r8,  0.9004E-02_r8,  0.1182E-01_r8,  0.1562E-01_r8,   &
         0.2081E-01_r8,  0.2794E-01_r8,  0.3773E-01_r8,  0.5111E-01_r8,  0.6919E-01_r8,   &
         0.9329E-01_r8, & 
                                ! DATA ((h22(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 4_r8, 4)_r8, &                                &
         -0.5623E-07_r8, -0.1071E-06_r8, -0.1984E-06_r8, -0.3536E-06_r8, -0.5997E-06_r8,  &
         -0.9606E-06_r8, -0.1447E-05_r8, -0.2055E-05_r8, -0.2775E-05_r8, -0.3596E-05_r8,  &
         -0.4504E-05_r8, -0.5529E-05_r8, -0.6768E-05_r8, -0.8345E-05_r8, -0.1039E-04_r8,  &
         -0.1304E-04_r8, -0.1645E-04_r8, -0.2082E-04_r8, -0.2633E-04_r8, -0.3322E-04_r8,  &
         -0.4187E-04_r8, -0.5292E-04_r8, -0.6730E-04_r8, -0.8640E-04_r8, -0.1122E-03_r8,  &
         -0.1472E-03_r8, -0.1948E-03_r8, -0.2585E-03_r8, -0.3428E-03_r8, -0.4523E-03_r8,  &
         -0.5915E-03_r8, & 
                                ! DATA ((h23(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 4_r8, 4)_r8, &                                &
         -0.2196E-09_r8, -0.4034E-09_r8, -0.7050E-09_r8, -0.1154E-08_r8, -0.1740E-08_r8,  &
         -0.2398E-08_r8, -0.3024E-08_r8, -0.3552E-08_r8, -0.4037E-08_r8, -0.4428E-08_r8,  &
         -0.4756E-08_r8, -0.5741E-08_r8, -0.8418E-08_r8, -0.1310E-07_r8, -0.1910E-07_r8,  &
         -0.2575E-07_r8, -0.3293E-07_r8, -0.3998E-07_r8, -0.4572E-07_r8, -0.4980E-07_r8,  &
         -0.5211E-07_r8, -0.5287E-07_r8, -0.5390E-07_r8, -0.5782E-07_r8, -0.6650E-07_r8,  &
         -0.7892E-07_r8, -0.8940E-07_r8, -0.8980E-07_r8, -0.7119E-07_r8, -0.2452E-07_r8,  &
         0.5823E-07_r8, & 
                                ! DATA ((h21(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 5_r8, 5)_r8, &                                &
         0.3919E-05_r8,  0.7618E-05_r8,  0.1455E-04_r8,  0.2709E-04_r8,  0.4863E-04_r8,   &
         0.8337E-04_r8,  0.1354E-03_r8,  0.2075E-03_r8,  0.3007E-03_r8,  0.4152E-03_r8,   &
         0.5500E-03_r8,  0.7038E-03_r8,  0.8793E-03_r8,  0.1086E-02_r8,  0.1345E-02_r8,   &
         0.1680E-02_r8,  0.2123E-02_r8,  0.2712E-02_r8,  0.3500E-02_r8,  0.4552E-02_r8,   &
         0.5958E-02_r8,  0.7844E-02_r8,  0.1038E-01_r8,  0.1381E-01_r8,  0.1847E-01_r8,   &
         0.2487E-01_r8,  0.3366E-01_r8,  0.4572E-01_r8,  0.6209E-01_r8,  0.8406E-01_r8,   &
         0.1130E+00_r8, & 
                                ! DATA ((h22(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 5_r8, 5)_r8, &                                &
         -0.5623E-07_r8, -0.1071E-06_r8, -0.1985E-06_r8, -0.3538E-06_r8, -0.6002E-06_r8,  &
         -0.9618E-06_r8, -0.1450E-05_r8, -0.2059E-05_r8, -0.2783E-05_r8, -0.3611E-05_r8,  &
         -0.4531E-05_r8, -0.5577E-05_r8, -0.6855E-05_r8, -0.8499E-05_r8, -0.1066E-04_r8,  &
         -0.1351E-04_r8, -0.1723E-04_r8, -0.2207E-04_r8, -0.2829E-04_r8, -0.3621E-04_r8,  &
         -0.4636E-04_r8, -0.5954E-04_r8, -0.7690E-04_r8, -0.1002E-03_r8, -0.1317E-03_r8,  &
         -0.1746E-03_r8, -0.2326E-03_r8, -0.3099E-03_r8, -0.4111E-03_r8, -0.5407E-03_r8,  &
         -0.7020E-03_r8, & 
                                ! DATA ((h23(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 5_r8, 5)_r8, &                                &
         -0.2197E-09_r8, -0.4037E-09_r8, -0.7054E-09_r8, -0.1155E-08_r8, -0.1741E-08_r8,  &
         -0.2401E-08_r8, -0.3027E-08_r8, -0.3553E-08_r8, -0.4039E-08_r8, -0.4431E-08_r8,  &
         -0.4775E-08_r8, -0.5784E-08_r8, -0.8506E-08_r8, -0.1326E-07_r8, -0.1931E-07_r8,  &
         -0.2600E-07_r8, -0.3324E-07_r8, -0.4048E-07_r8, -0.4666E-07_r8, -0.5137E-07_r8,  &
         -0.5428E-07_r8, -0.5558E-07_r8, -0.5730E-07_r8, -0.6228E-07_r8, -0.7197E-07_r8,  &
         -0.8455E-07_r8, -0.9347E-07_r8, -0.8867E-07_r8, -0.5945E-07_r8,  0.5512E-08_r8,  &
         0.1209E-06_r8, & 
                                ! DATA ((h21(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 6_r8, 6)_r8, &                                &
         0.3919E-05_r8,  0.7619E-05_r8,  0.1456E-04_r8,  0.2710E-04_r8,  0.4868E-04_r8,   &
         0.8351E-04_r8,  0.1357E-03_r8,  0.2082E-03_r8,  0.3022E-03_r8,  0.4181E-03_r8,   &
         0.5554E-03_r8,  0.7134E-03_r8,  0.8963E-03_r8,  0.1117E-02_r8,  0.1397E-02_r8,   &
         0.1770E-02_r8,  0.2270E-02_r8,  0.2947E-02_r8,  0.3865E-02_r8,  0.5107E-02_r8,   &
         0.6787E-02_r8,  0.9062E-02_r8,  0.1215E-01_r8,  0.1634E-01_r8,  0.2209E-01_r8,   &
         0.2999E-01_r8,  0.4083E-01_r8,  0.5563E-01_r8,  0.7559E-01_r8,  0.1021E+00_r8,   &
         0.1364E+00_r8, & 
                                ! DATA ((h22(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 6_r8, 6)_r8, &                                &
         -0.5624E-07_r8, -0.1072E-06_r8, -0.1986E-06_r8, -0.3541E-06_r8, -0.6010E-06_r8,  &
         -0.9636E-06_r8, -0.1453E-05_r8, -0.2067E-05_r8, -0.2796E-05_r8, -0.3634E-05_r8,  &
         -0.4572E-05_r8, -0.5652E-05_r8, -0.6987E-05_r8, -0.8733E-05_r8, -0.1107E-04_r8,  &
         -0.1418E-04_r8, -0.1832E-04_r8, -0.2378E-04_r8, -0.3092E-04_r8, -0.4017E-04_r8,  &
         -0.5221E-04_r8, -0.6806E-04_r8, -0.8916E-04_r8, -0.1176E-03_r8, -0.1562E-03_r8,  &
         -0.2087E-03_r8, -0.2793E-03_r8, -0.3724E-03_r8, -0.4928E-03_r8, -0.6440E-03_r8,  &
         -0.8270E-03_r8, & 
                                ! DATA ((h23(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 6_r8, 6)_r8, &                                &
         -0.2198E-09_r8, -0.4040E-09_r8, -0.7061E-09_r8, -0.1156E-08_r8, -0.1744E-08_r8,  &
         -0.2405E-08_r8, -0.3032E-08_r8, -0.3556E-08_r8, -0.4040E-08_r8, -0.4444E-08_r8,  &
         -0.4800E-08_r8, -0.5848E-08_r8, -0.8640E-08_r8, -0.1346E-07_r8, -0.1957E-07_r8,  &
         -0.2627E-07_r8, -0.3357E-07_r8, -0.4114E-07_r8, -0.4793E-07_r8, -0.5330E-07_r8,  &
         -0.5676E-07_r8, -0.5873E-07_r8, -0.6152E-07_r8, -0.6783E-07_r8, -0.7834E-07_r8,  &
         -0.9023E-07_r8, -0.9530E-07_r8, -0.8162E-07_r8, -0.3634E-07_r8,  0.5638E-07_r8,  &
         0.2189E-06_r8, & 
                                ! DATA ((h21(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 7_r8, 7)_r8, &                                &
         0.3919E-05_r8,  0.7620E-05_r8,  0.1457E-04_r8,  0.2713E-04_r8,  0.4877E-04_r8,   &
         0.8374E-04_r8,  0.1363E-03_r8,  0.2095E-03_r8,  0.3047E-03_r8,  0.4227E-03_r8,   &
         0.5637E-03_r8,  0.7282E-03_r8,  0.9224E-03_r8,  0.1162E-02_r8,  0.1475E-02_r8,   &
         0.1898E-02_r8,  0.2475E-02_r8,  0.3267E-02_r8,  0.4354E-02_r8,  0.5839E-02_r8,   &
         0.7866E-02_r8,  0.1063E-01_r8,  0.1440E-01_r8,  0.1957E-01_r8,  0.2666E-01_r8,   &
         0.3641E-01_r8,  0.4975E-01_r8,  0.6784E-01_r8,  0.9201E-01_r8,  0.1236E+00_r8,   &
         0.1639E+00_r8, & 
                                ! DATA ((h22(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 7_r8, 7)_r8, &                                &
         -0.5626E-07_r8, -0.1072E-06_r8, -0.1987E-06_r8, -0.3545E-06_r8, -0.6022E-06_r8,  &
         -0.9665E-06_r8, -0.1460E-05_r8, -0.2078E-05_r8, -0.2817E-05_r8, -0.3671E-05_r8,  &
         -0.4637E-05_r8, -0.5767E-05_r8, -0.7188E-05_r8, -0.9080E-05_r8, -0.1165E-04_r8,  &
         -0.1513E-04_r8, -0.1981E-04_r8, -0.2609E-04_r8, -0.3441E-04_r8, -0.4534E-04_r8,  &
         -0.5978E-04_r8, -0.7897E-04_r8, -0.1047E-03_r8, -0.1396E-03_r8, -0.1870E-03_r8,  &
         -0.2510E-03_r8, -0.3363E-03_r8, -0.4475E-03_r8, -0.5888E-03_r8, -0.7621E-03_r8,  &
         -0.9647E-03_r8, & 
                                ! DATA ((h23(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 7_r8, 7)_r8, &                                &
         -0.2200E-09_r8, -0.4045E-09_r8, -0.7071E-09_r8, -0.1159E-08_r8, -0.1748E-08_r8,  &
         -0.2411E-08_r8, -0.3040E-08_r8, -0.3561E-08_r8, -0.4046E-08_r8, -0.4455E-08_r8,  &
         -0.4839E-08_r8, -0.5941E-08_r8, -0.8815E-08_r8, -0.1371E-07_r8, -0.1983E-07_r8,  &
         -0.2652E-07_r8, -0.3400E-07_r8, -0.4207E-07_r8, -0.4955E-07_r8, -0.5554E-07_r8,  &
         -0.5966E-07_r8, -0.6261E-07_r8, -0.6688E-07_r8, -0.7454E-07_r8, -0.8521E-07_r8,  &
         -0.9470E-07_r8, -0.9275E-07_r8, -0.6525E-07_r8,  0.3686E-08_r8,  0.1371E-06_r8,  &
         0.3623E-06_r8, & 
                                ! DATA ((h21(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 8_r8, 8)_r8, &                                &
         0.3919E-05_r8,  0.7622E-05_r8,  0.1458E-04_r8,  0.2717E-04_r8,  0.4890E-04_r8,   &
         0.8410E-04_r8,  0.1372E-03_r8,  0.2114E-03_r8,  0.3085E-03_r8,  0.4298E-03_r8,   &
         0.5765E-03_r8,  0.7508E-03_r8,  0.9618E-03_r8,  0.1229E-02_r8,  0.1586E-02_r8,   &
         0.2077E-02_r8,  0.2757E-02_r8,  0.3697E-02_r8,  0.5000E-02_r8,  0.6794E-02_r8,   &
         0.9260E-02_r8,  0.1264E-01_r8,  0.1728E-01_r8,  0.2365E-01_r8,  0.3242E-01_r8,   &
         0.4444E-01_r8,  0.6080E-01_r8,  0.8278E-01_r8,  0.1118E+00_r8,  0.1491E+00_r8,   &
         0.1956E+00_r8, & 
                                ! DATA ((h22(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 8_r8, 8)_r8, &                                &
         -0.5628E-07_r8, -0.1073E-06_r8, -0.1990E-06_r8, -0.3553E-06_r8, -0.6042E-06_r8,  &
         -0.9710E-06_r8, -0.1469E-05_r8, -0.2096E-05_r8, -0.2849E-05_r8, -0.3728E-05_r8,  &
         -0.4738E-05_r8, -0.5942E-05_r8, -0.7490E-05_r8, -0.9586E-05_r8, -0.1247E-04_r8,  &
         -0.1644E-04_r8, -0.2184E-04_r8, -0.2916E-04_r8, -0.3898E-04_r8, -0.5205E-04_r8,  &
         -0.6948E-04_r8, -0.9285E-04_r8, -0.1244E-03_r8, -0.1672E-03_r8, -0.2251E-03_r8,  &
         -0.3028E-03_r8, -0.4051E-03_r8, -0.5365E-03_r8, -0.6998E-03_r8, -0.8940E-03_r8,  &
         -0.1112E-02_r8, & 
                                ! DATA ((h23(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 8_r8, 8)_r8, &                                &
         -0.2204E-09_r8, -0.4052E-09_r8, -0.7088E-09_r8, -0.1162E-08_r8, -0.1755E-08_r8,  &
         -0.2422E-08_r8, -0.3053E-08_r8, -0.3572E-08_r8, -0.4052E-08_r8, -0.4474E-08_r8,  &
         -0.4898E-08_r8, -0.6082E-08_r8, -0.9046E-08_r8, -0.1400E-07_r8, -0.2009E-07_r8,  &
         -0.2683E-07_r8, -0.3463E-07_r8, -0.4334E-07_r8, -0.5153E-07_r8, -0.5811E-07_r8,  &
         -0.6305E-07_r8, -0.6749E-07_r8, -0.7346E-07_r8, -0.8208E-07_r8, -0.9173E-07_r8,  &
         -0.9603E-07_r8, -0.8264E-07_r8, -0.3505E-07_r8,  0.6878E-07_r8,  0.2586E-06_r8,  &
         0.5530E-06_r8, & 
                                ! DATA ((h21(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 9_r8, 9)_r8, &                                &
         0.3919E-05_r8,  0.7625E-05_r8,  0.1459E-04_r8,  0.2723E-04_r8,  0.4910E-04_r8,   &
         0.8465E-04_r8,  0.1385E-03_r8,  0.2143E-03_r8,  0.3143E-03_r8,  0.4407E-03_r8,   &
         0.5960E-03_r8,  0.7850E-03_r8,  0.1020E-02_r8,  0.1326E-02_r8,  0.1744E-02_r8,   &
         0.2324E-02_r8,  0.3135E-02_r8,  0.4267E-02_r8,  0.5845E-02_r8,  0.8031E-02_r8,   &
         0.1105E-01_r8,  0.1520E-01_r8,  0.2093E-01_r8,  0.2880E-01_r8,  0.3962E-01_r8,   &
         0.5440E-01_r8,  0.7435E-01_r8,  0.1009E+00_r8,  0.1353E+00_r8,  0.1787E+00_r8,   &
         0.2317E+00_r8, & 
                                ! DATA ((h22(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 9_r8, 9)_r8, &                                &
         -0.5630E-07_r8, -0.1074E-06_r8, -0.1994E-06_r8, -0.3564E-06_r8, -0.6072E-06_r8,  &
         -0.9779E-06_r8, -0.1484E-05_r8, -0.2124E-05_r8, -0.2900E-05_r8, -0.3817E-05_r8,  &
         -0.4891E-05_r8, -0.6205E-05_r8, -0.7931E-05_r8, -0.1031E-04_r8, -0.1362E-04_r8,  &
         -0.1821E-04_r8, -0.2454E-04_r8, -0.3320E-04_r8, -0.4493E-04_r8, -0.6068E-04_r8,  &
         -0.8186E-04_r8, -0.1104E-03_r8, -0.1491E-03_r8, -0.2016E-03_r8, -0.2722E-03_r8,  &
         -0.3658E-03_r8, -0.4873E-03_r8, -0.6404E-03_r8, -0.8253E-03_r8, -0.1037E-02_r8,  &
         -0.1267E-02_r8, & 
                                ! DATA ((h23(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 9_r8, 9)_r8, &                                &
         -0.2207E-09_r8, -0.4061E-09_r8, -0.7117E-09_r8, -0.1169E-08_r8, -0.1767E-08_r8,  &
         -0.2439E-08_r8, -0.3074E-08_r8, -0.3588E-08_r8, -0.4062E-08_r8, -0.4510E-08_r8,  &
         -0.4983E-08_r8, -0.6261E-08_r8, -0.9324E-08_r8, -0.1430E-07_r8, -0.2036E-07_r8,  &
         -0.2725E-07_r8, -0.3561E-07_r8, -0.4505E-07_r8, -0.5384E-07_r8, -0.6111E-07_r8,  &
         -0.6731E-07_r8, -0.7355E-07_r8, -0.8112E-07_r8, -0.8978E-07_r8, -0.9616E-07_r8,  &
         -0.9157E-07_r8, -0.6114E-07_r8,  0.1622E-07_r8,  0.1694E-06_r8,  0.4277E-06_r8,  &
         0.7751E-06_r8, & 
                                ! DATA ((h21(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=10_r8,10)_r8, &                                &
         0.3919E-05_r8,  0.7629E-05_r8,  0.1462E-04_r8,  0.2733E-04_r8,  0.4940E-04_r8,   &
         0.8549E-04_r8,  0.1406E-03_r8,  0.2188E-03_r8,  0.3232E-03_r8,  0.4572E-03_r8,   &
         0.6253E-03_r8,  0.8354E-03_r8,  0.1105E-02_r8,  0.1464E-02_r8,  0.1960E-02_r8,   &
         0.2657E-02_r8,  0.3637E-02_r8,  0.5014E-02_r8,  0.6942E-02_r8,  0.9622E-02_r8,   &
         0.1333E-01_r8,  0.1846E-01_r8,  0.2553E-01_r8,  0.3526E-01_r8,  0.4859E-01_r8,   &
         0.6667E-01_r8,  0.9084E-01_r8,  0.1225E+00_r8,  0.1629E+00_r8,  0.2127E+00_r8,   &
         0.2721E+00_r8, & 
                                ! DATA ((h22(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=10_r8,10)_r8, &                                &
         -0.5636E-07_r8, -0.1076E-06_r8, -0.2000E-06_r8, -0.3582E-06_r8, -0.6119E-06_r8,  &
         -0.9888E-06_r8, -0.1507E-05_r8, -0.2168E-05_r8, -0.2978E-05_r8, -0.3952E-05_r8,  &
         -0.5122E-05_r8, -0.6592E-05_r8, -0.8565E-05_r8, -0.1132E-04_r8, -0.1518E-04_r8,  &
         -0.2060E-04_r8, -0.2811E-04_r8, -0.3848E-04_r8, -0.5261E-04_r8, -0.7173E-04_r8,  &
         -0.9758E-04_r8, -0.1326E-03_r8, -0.1801E-03_r8, -0.2442E-03_r8, -0.3296E-03_r8,  &
         -0.4415E-03_r8, -0.5840E-03_r8, -0.7591E-03_r8, -0.9636E-03_r8, -0.1189E-02_r8,  &
         -0.1427E-02_r8, & 
                                ! DATA ((h23(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=10_r8,10)_r8, &                                &
         -0.2214E-09_r8, -0.4080E-09_r8, -0.7156E-09_r8, -0.1178E-08_r8, -0.1784E-08_r8,  &
         -0.2466E-08_r8, -0.3105E-08_r8, -0.3617E-08_r8, -0.4087E-08_r8, -0.4563E-08_r8,  &
         -0.5110E-08_r8, -0.6492E-08_r8, -0.9643E-08_r8, -0.1461E-07_r8, -0.2069E-07_r8,  &
         -0.2796E-07_r8, -0.3702E-07_r8, -0.4717E-07_r8, -0.5662E-07_r8, -0.6484E-07_r8,  &
         -0.7271E-07_r8, -0.8079E-07_r8, -0.8928E-07_r8, -0.9634E-07_r8, -0.9625E-07_r8,  &
         -0.7776E-07_r8, -0.2242E-07_r8,  0.9745E-07_r8,  0.3152E-06_r8,  0.6388E-06_r8,  &
         0.9992E-06_r8, & 
                                ! DATA ((h21(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=11_r8,11)_r8, &                                &
         0.3919E-05_r8,  0.7639E-05_r8,  0.1466E-04_r8,  0.2748E-04_r8,  0.4986E-04_r8,   &
         0.8675E-04_r8,  0.1437E-03_r8,  0.2255E-03_r8,  0.3365E-03_r8,  0.4816E-03_r8,   &
         0.6682E-03_r8,  0.9082E-03_r8,  0.1224E-02_r8,  0.1653E-02_r8,  0.2253E-02_r8,   &
         0.3099E-02_r8,  0.4296E-02_r8,  0.5984E-02_r8,  0.8352E-02_r8,  0.1165E-01_r8,   &
         0.1623E-01_r8,  0.2257E-01_r8,  0.3132E-01_r8,  0.4333E-01_r8,  0.5968E-01_r8,   &
         0.8166E-01_r8,  0.1107E+00_r8,  0.1481E+00_r8,  0.1947E+00_r8,  0.2509E+00_r8,   &
         0.3170E+00_r8, & 
                                ! DATA ((h22(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=11_r8,11)_r8, &                                &
         -0.5645E-07_r8, -0.1079E-06_r8, -0.2009E-06_r8, -0.3610E-06_r8, -0.6190E-06_r8,  &
         -0.1005E-05_r8, -0.1541E-05_r8, -0.2235E-05_r8, -0.3096E-05_r8, -0.4155E-05_r8,  &
         -0.5463E-05_r8, -0.7150E-05_r8, -0.9453E-05_r8, -0.1269E-04_r8, -0.1728E-04_r8,  &
         -0.2375E-04_r8, -0.3279E-04_r8, -0.4531E-04_r8, -0.6246E-04_r8, -0.8579E-04_r8,  &
         -0.1175E-03_r8, -0.1604E-03_r8, -0.2186E-03_r8, -0.2964E-03_r8, -0.3990E-03_r8,  &
         -0.5311E-03_r8, -0.6957E-03_r8, -0.8916E-03_r8, -0.1112E-02_r8, -0.1346E-02_r8,  &
         -0.1590E-02_r8, & 
                                ! DATA ((h23(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=11_r8,11)_r8, &                                &
         -0.2225E-09_r8, -0.4104E-09_r8, -0.7217E-09_r8, -0.1192E-08_r8, -0.1811E-08_r8,  &
         -0.2509E-08_r8, -0.3155E-08_r8, -0.3668E-08_r8, -0.4138E-08_r8, -0.4650E-08_r8,  &
         -0.5296E-08_r8, -0.6785E-08_r8, -0.9991E-08_r8, -0.1494E-07_r8, -0.2122E-07_r8,  &
         -0.2911E-07_r8, -0.3895E-07_r8, -0.4979E-07_r8, -0.6002E-07_r8, -0.6964E-07_r8,  &
         -0.7935E-07_r8, -0.8887E-07_r8, -0.9699E-07_r8, -0.9967E-07_r8, -0.8883E-07_r8,  &
         -0.4988E-07_r8,  0.4156E-07_r8,  0.2197E-06_r8,  0.5081E-06_r8,  0.8667E-06_r8,  &
         0.1212E-05_r8, & 
                                ! DATA ((h21(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=12_r8,12)_r8, &                                &
         0.3920E-05_r8,  0.7649E-05_r8,  0.1472E-04_r8,  0.2768E-04_r8,  0.5051E-04_r8,   &
         0.8856E-04_r8,  0.1481E-03_r8,  0.2352E-03_r8,  0.3557E-03_r8,  0.5169E-03_r8,   &
         0.7295E-03_r8,  0.1010E-02_r8,  0.1388E-02_r8,  0.1909E-02_r8,  0.2640E-02_r8,   &
         0.3678E-02_r8,  0.5151E-02_r8,  0.7231E-02_r8,  0.1015E-01_r8,  0.1423E-01_r8,   &
         0.1990E-01_r8,  0.2775E-01_r8,  0.3856E-01_r8,  0.5333E-01_r8,  0.7328E-01_r8,   &
         0.9980E-01_r8,  0.1343E+00_r8,  0.1778E+00_r8,  0.2308E+00_r8,  0.2936E+00_r8,   &
         0.3667E+00_r8, & 
                                ! DATA ((h22(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=12_r8,12)_r8, &                                &
         -0.5658E-07_r8, -0.1083E-06_r8, -0.2023E-06_r8, -0.3650E-06_r8, -0.6295E-06_r8,  &
         -0.1030E-05_r8, -0.1593E-05_r8, -0.2334E-05_r8, -0.3273E-05_r8, -0.4455E-05_r8,  &
         -0.5955E-05_r8, -0.7935E-05_r8, -0.1067E-04_r8, -0.1455E-04_r8, -0.2007E-04_r8,  &
         -0.2788E-04_r8, -0.3885E-04_r8, -0.5409E-04_r8, -0.7503E-04_r8, -0.1036E-03_r8,  &
         -0.1425E-03_r8, -0.1951E-03_r8, -0.2660E-03_r8, -0.3598E-03_r8, -0.4817E-03_r8,  &
         -0.6355E-03_r8, -0.8218E-03_r8, -0.1035E-02_r8, -0.1267E-02_r8, -0.1508E-02_r8,  &
         -0.1755E-02_r8, & 
                                ! DATA ((h23(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=12_r8,12)_r8, &                                &
         -0.2241E-09_r8, -0.4142E-09_r8, -0.7312E-09_r8, -0.1214E-08_r8, -0.1854E-08_r8,  &
         -0.2578E-08_r8, -0.3250E-08_r8, -0.3765E-08_r8, -0.4238E-08_r8, -0.4809E-08_r8,  &
         -0.5553E-08_r8, -0.7132E-08_r8, -0.1035E-07_r8, -0.1538E-07_r8, -0.2211E-07_r8,  &
         -0.3079E-07_r8, -0.4142E-07_r8, -0.5303E-07_r8, -0.6437E-07_r8, -0.7566E-07_r8,  &
         -0.8703E-07_r8, -0.9700E-07_r8, -0.1025E-06_r8, -0.9718E-07_r8, -0.6973E-07_r8,  &
         -0.5265E-09_r8,  0.1413E-06_r8,  0.3895E-06_r8,  0.7321E-06_r8,  0.1085E-05_r8,  &
         0.1449E-05_r8, & 
                                ! DATA ((h21(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=13_r8,13)_r8, &                                &
         0.3920E-05_r8,  0.7668E-05_r8,  0.1480E-04_r8,  0.2798E-04_r8,  0.5143E-04_r8,   &
         0.9106E-04_r8,  0.1542E-03_r8,  0.2488E-03_r8,  0.3827E-03_r8,  0.5661E-03_r8,   &
         0.8143E-03_r8,  0.1150E-02_r8,  0.1608E-02_r8,  0.2246E-02_r8,  0.3147E-02_r8,   &
         0.4428E-02_r8,  0.6248E-02_r8,  0.8823E-02_r8,  0.1244E-01_r8,  0.1750E-01_r8,   &
         0.2452E-01_r8,  0.3424E-01_r8,  0.4756E-01_r8,  0.6564E-01_r8,  0.8981E-01_r8,   &
         0.1215E+00_r8,  0.1619E+00_r8,  0.2118E+00_r8,  0.2713E+00_r8,  0.3409E+00_r8,   &
         0.4215E+00_r8, & 
                                ! DATA ((h22(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=13_r8,13)_r8, &                                &
         -0.5677E-07_r8, -0.1090E-06_r8, -0.2043E-06_r8, -0.3709E-06_r8, -0.6448E-06_r8,  &
         -0.1066E-05_r8, -0.1669E-05_r8, -0.2480E-05_r8, -0.3532E-05_r8, -0.4887E-05_r8,  &
         -0.6650E-05_r8, -0.9017E-05_r8, -0.1232E-04_r8, -0.1701E-04_r8, -0.2372E-04_r8,  &
         -0.3325E-04_r8, -0.4664E-04_r8, -0.6528E-04_r8, -0.9095E-04_r8, -0.1260E-03_r8,  &
         -0.1737E-03_r8, -0.2381E-03_r8, -0.3238E-03_r8, -0.4359E-03_r8, -0.5789E-03_r8,  &
         -0.7549E-03_r8, -0.9607E-03_r8, -0.1188E-02_r8, -0.1427E-02_r8, -0.1674E-02_r8,  &
         -0.1914E-02_r8, & 
                                ! DATA ((h23(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=13_r8,13)_r8, &                                &
         -0.2262E-09_r8, -0.4191E-09_r8, -0.7441E-09_r8, -0.1244E-08_r8, -0.1916E-08_r8,  &
         -0.2687E-08_r8, -0.3407E-08_r8, -0.3947E-08_r8, -0.4432E-08_r8, -0.5059E-08_r8,  &
         -0.5896E-08_r8, -0.7538E-08_r8, -0.1079E-07_r8, -0.1606E-07_r8, -0.2346E-07_r8,  &
         -0.3305E-07_r8, -0.4459E-07_r8, -0.5720E-07_r8, -0.6999E-07_r8, -0.8300E-07_r8,  &
         -0.9536E-07_r8, -0.1039E-06_r8, -0.1036E-06_r8, -0.8526E-07_r8, -0.3316E-07_r8,  &
         0.7909E-07_r8,  0.2865E-06_r8,  0.6013E-06_r8,  0.9580E-06_r8,  0.1303E-05_r8,   &
         0.1792E-05_r8, & 
                                ! DATA ((h21(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=14_r8,14)_r8, &                                &
         0.3923E-05_r8,  0.7691E-05_r8,  0.1491E-04_r8,  0.2836E-04_r8,  0.5261E-04_r8,   &
         0.9432E-04_r8,  0.1623E-03_r8,  0.2667E-03_r8,  0.4187E-03_r8,  0.6323E-03_r8,   &
         0.9278E-03_r8,  0.1334E-02_r8,  0.1897E-02_r8,  0.2685E-02_r8,  0.3800E-02_r8,   &
         0.5388E-02_r8,  0.7646E-02_r8,  0.1084E-01_r8,  0.1533E-01_r8,  0.2160E-01_r8,   &
         0.3031E-01_r8,  0.4231E-01_r8,  0.5866E-01_r8,  0.8064E-01_r8,  0.1097E+00_r8,   &
         0.1471E+00_r8,  0.1938E+00_r8,  0.2501E+00_r8,  0.3163E+00_r8,  0.3932E+00_r8,   &
         0.4810E+00_r8, & 
                                ! DATA ((h22(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=14_r8,14)_r8, &                                &
         -0.5703E-07_r8, -0.1098E-06_r8, -0.2071E-06_r8, -0.3788E-06_r8, -0.6657E-06_r8,  &
         -0.1116E-05_r8, -0.1776E-05_r8, -0.2687E-05_r8, -0.3898E-05_r8, -0.5493E-05_r8,  &
         -0.7607E-05_r8, -0.1048E-04_r8, -0.1450E-04_r8, -0.2024E-04_r8, -0.2845E-04_r8,  &
         -0.4014E-04_r8, -0.5658E-04_r8, -0.7947E-04_r8, -0.1110E-03_r8, -0.1541E-03_r8,  &
         -0.2125E-03_r8, -0.2907E-03_r8, -0.3936E-03_r8, -0.5261E-03_r8, -0.6912E-03_r8,  &
         -0.8880E-03_r8, -0.1109E-02_r8, -0.1346E-02_r8, -0.1591E-02_r8, -0.1837E-02_r8,  &
         -0.2054E-02_r8, & 
                                ! DATA ((h23(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=14_r8,14)_r8, &                                &
         -0.2288E-09_r8, -0.4265E-09_r8, -0.7627E-09_r8, -0.1289E-08_r8, -0.2011E-08_r8,  &
         -0.2861E-08_r8, -0.3673E-08_r8, -0.4288E-08_r8, -0.4812E-08_r8, -0.5475E-08_r8,  &
         -0.6365E-08_r8, -0.8052E-08_r8, -0.1142E-07_r8, -0.1711E-07_r8, -0.2533E-07_r8,  &
         -0.3597E-07_r8, -0.4862E-07_r8, -0.6259E-07_r8, -0.7708E-07_r8, -0.9150E-07_r8,  &
         -0.1035E-06_r8, -0.1079E-06_r8, -0.9742E-07_r8, -0.5928E-07_r8,  0.2892E-07_r8,  &
         0.1998E-06_r8,  0.4789E-06_r8,  0.8298E-06_r8,  0.1172E-05_r8,  0.1583E-05_r8,   &
         0.2329E-05_r8, & 
                                ! DATA ((h21(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=15_r8,15)_r8, &                                &
         0.3928E-05_r8,  0.7724E-05_r8,  0.1504E-04_r8,  0.2883E-04_r8,  0.5403E-04_r8,   &
         0.9826E-04_r8,  0.1722E-03_r8,  0.2891E-03_r8,  0.4644E-03_r8,  0.7172E-03_r8,   &
         0.1074E-02_r8,  0.1572E-02_r8,  0.2267E-02_r8,  0.3245E-02_r8,  0.4632E-02_r8,   &
         0.6606E-02_r8,  0.9411E-02_r8,  0.1338E-01_r8,  0.1896E-01_r8,  0.2674E-01_r8,   &
         0.3751E-01_r8,  0.5227E-01_r8,  0.7222E-01_r8,  0.9875E-01_r8,  0.1332E+00_r8,   &
         0.1768E+00_r8,  0.2299E+00_r8,  0.2927E+00_r8,  0.3660E+00_r8,  0.4504E+00_r8,   &
         0.5444E+00_r8, & 
                                ! DATA ((h22(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=15_r8,15)_r8, &                                &
         -0.5736E-07_r8, -0.1109E-06_r8, -0.2106E-06_r8, -0.3890E-06_r8, -0.6928E-06_r8,  &
         -0.1181E-05_r8, -0.1917E-05_r8, -0.2965E-05_r8, -0.4396E-05_r8, -0.6315E-05_r8,  &
         -0.8891E-05_r8, -0.1242E-04_r8, -0.1736E-04_r8, -0.2442E-04_r8, -0.3454E-04_r8,  &
         -0.4892E-04_r8, -0.6916E-04_r8, -0.9735E-04_r8, -0.1362E-03_r8, -0.1891E-03_r8,  &
         -0.2602E-03_r8, -0.3545E-03_r8, -0.4768E-03_r8, -0.6310E-03_r8, -0.8179E-03_r8,  &
         -0.1032E-02_r8, -0.1265E-02_r8, -0.1508E-02_r8, -0.1757E-02_r8, -0.1989E-02_r8,  &
         -0.2159E-02_r8, & 
                                ! DATA ((h23(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=15_r8,15)_r8, &                                &
         -0.2321E-09_r8, -0.4350E-09_r8, -0.7861E-09_r8, -0.1347E-08_r8, -0.2144E-08_r8,  &
         -0.3120E-08_r8, -0.4107E-08_r8, -0.4892E-08_r8, -0.5511E-08_r8, -0.6164E-08_r8,  &
         -0.7054E-08_r8, -0.8811E-08_r8, -0.1240E-07_r8, -0.1862E-07_r8, -0.2773E-07_r8,  &
         -0.3957E-07_r8, -0.5371E-07_r8, -0.6939E-07_r8, -0.8564E-07_r8, -0.1006E-06_r8,  &
         -0.1100E-06_r8, -0.1066E-06_r8, -0.8018E-07_r8, -0.1228E-07_r8,  0.1263E-06_r8,  &
         0.3678E-06_r8,  0.7022E-06_r8,  0.1049E-05_r8,  0.1411E-05_r8,  0.2015E-05_r8,   &
         0.3099E-05_r8, & 
                                ! DATA ((h21(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=16_r8,16)_r8, &                                &
         0.2064E-05_r8,  0.4104E-05_r8,  0.8141E-05_r8,  0.1607E-04_r8,  0.3145E-04_r8,   &
         0.6080E-04_r8,  0.1153E-03_r8,  0.2129E-03_r8,  0.3802E-03_r8,  0.6538E-03_r8,   &
         0.1078E-02_r8,  0.1699E-02_r8,  0.2576E-02_r8,  0.3805E-02_r8,  0.5540E-02_r8,   &
         0.8004E-02_r8,  0.1149E-01_r8,  0.1640E-01_r8,  0.2329E-01_r8,  0.3290E-01_r8,   &
         0.4622E-01_r8,  0.6441E-01_r8,  0.8888E-01_r8,  0.1210E+00_r8,  0.1622E+00_r8,   &
         0.2132E+00_r8,  0.2745E+00_r8,  0.3470E+00_r8,  0.4317E+00_r8,  0.5286E+00_r8,   &
         0.6335E+00_r8, & 
                                ! DATA ((h22(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=16_r8,16)_r8, &                                &
         -0.3122E-07_r8, -0.6175E-07_r8, -0.1214E-06_r8, -0.2361E-06_r8, -0.4518E-06_r8,  &
         -0.8438E-06_r8, -0.1524E-05_r8, -0.2643E-05_r8, -0.4380E-05_r8, -0.6922E-05_r8,  &
         -0.1042E-04_r8, -0.1504E-04_r8, -0.2125E-04_r8, -0.2987E-04_r8, -0.4200E-04_r8,  &
         -0.5923E-04_r8, -0.8383E-04_r8, -0.1186E-03_r8, -0.1670E-03_r8, -0.2328E-03_r8,  &
         -0.3204E-03_r8, -0.4347E-03_r8, -0.5802E-03_r8, -0.7595E-03_r8, -0.9703E-03_r8,  &
         -0.1205E-02_r8, -0.1457E-02_r8, -0.1720E-02_r8, -0.1980E-02_r8, -0.2191E-02_r8,  &
         -0.2290E-02_r8, & 
                                ! DATA ((h23(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=16_r8,16)_r8, &                                &
         -0.1376E-09_r8, -0.2699E-09_r8, -0.5220E-09_r8, -0.9897E-09_r8, -0.1819E-08_r8,  &
         -0.3186E-08_r8, -0.5224E-08_r8, -0.7896E-08_r8, -0.1090E-07_r8, -0.1349E-07_r8,  &
         -0.1443E-07_r8, -0.1374E-07_r8, -0.1386E-07_r8, -0.1673E-07_r8, -0.2237E-07_r8,  &
         -0.3248E-07_r8, -0.5050E-07_r8, -0.7743E-07_r8, -0.1097E-06_r8, -0.1369E-06_r8,  &
         -0.1463E-06_r8, -0.1268E-06_r8, -0.6424E-07_r8,  0.5941E-07_r8,  0.2742E-06_r8,  &
         0.5924E-06_r8,  0.9445E-06_r8,  0.1286E-05_r8,  0.1819E-05_r8,  0.2867E-05_r8,   &
         0.4527E-05_r8, & 
                                ! DATA ((h21(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=17_r8,17)_r8, &                                &
         0.2452E-05_r8,  0.4881E-05_r8,  0.9693E-05_r8,  0.1918E-04_r8,  0.3768E-04_r8,   &
         0.7322E-04_r8,  0.1399E-03_r8,  0.2607E-03_r8,  0.4694E-03_r8,  0.8102E-03_r8,   &
         0.1334E-02_r8,  0.2098E-02_r8,  0.3180E-02_r8,  0.4710E-02_r8,  0.6879E-02_r8,   &
         0.9952E-02_r8,  0.1429E-01_r8,  0.2038E-01_r8,  0.2894E-01_r8,  0.4085E-01_r8,   &
         0.5722E-01_r8,  0.7939E-01_r8,  0.1088E+00_r8,  0.1468E+00_r8,  0.1945E+00_r8,   &
         0.2524E+00_r8,  0.3210E+00_r8,  0.4017E+00_r8,  0.4948E+00_r8,  0.5978E+00_r8,   &
         0.7040E+00_r8, & 
                                ! DATA ((h22(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=17_r8,17)_r8, &                                &
         -0.3547E-07_r8, -0.7029E-07_r8, -0.1386E-06_r8, -0.2709E-06_r8, -0.5218E-06_r8,  &
         -0.9840E-06_r8, -0.1799E-05_r8, -0.3156E-05_r8, -0.5272E-05_r8, -0.8357E-05_r8,  &
         -0.1260E-04_r8, -0.1827E-04_r8, -0.2598E-04_r8, -0.3667E-04_r8, -0.5169E-04_r8,  &
         -0.7312E-04_r8, -0.1037E-03_r8, -0.1467E-03_r8, -0.2060E-03_r8, -0.2857E-03_r8,  &
         -0.3907E-03_r8, -0.5257E-03_r8, -0.6940E-03_r8, -0.8954E-03_r8, -0.1124E-02_r8,  &
         -0.1371E-02_r8, -0.1632E-02_r8, -0.1897E-02_r8, -0.2131E-02_r8, -0.2275E-02_r8,  &
         -0.2265E-02_r8, & 
                                ! DATA ((h23(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=17_r8,17)_r8, &                                &
         -0.1482E-09_r8, -0.2910E-09_r8, -0.5667E-09_r8, -0.1081E-08_r8, -0.2005E-08_r8,  &
         -0.3554E-08_r8, -0.5902E-08_r8, -0.8925E-08_r8, -0.1209E-07_r8, -0.1448E-07_r8,  &
         -0.1536E-07_r8, -0.1565E-07_r8, -0.1763E-07_r8, -0.2088E-07_r8, -0.2564E-07_r8,  &
         -0.3635E-07_r8, -0.5791E-07_r8, -0.8907E-07_r8, -0.1213E-06_r8, -0.1418E-06_r8,  &
         -0.1397E-06_r8, -0.1000E-06_r8, -0.4427E-08_r8,  0.1713E-06_r8,  0.4536E-06_r8,  &
         0.8086E-06_r8,  0.1153E-05_r8,  0.1588E-05_r8,  0.2437E-05_r8,  0.3905E-05_r8,   &
         0.5874E-05_r8, & 
                                ! DATA ((h21(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=18_r8,18)_r8, &                                &
         0.2878E-05_r8,  0.5733E-05_r8,  0.1140E-04_r8,  0.2258E-04_r8,  0.4446E-04_r8,   &
         0.8668E-04_r8,  0.1664E-03_r8,  0.3120E-03_r8,  0.5656E-03_r8,  0.9812E-03_r8,   &
         0.1621E-02_r8,  0.2559E-02_r8,  0.3902E-02_r8,  0.5812E-02_r8,  0.8522E-02_r8,   &
         0.1235E-01_r8,  0.1773E-01_r8,  0.2530E-01_r8,  0.3590E-01_r8,  0.5058E-01_r8,   &
         0.7059E-01_r8,  0.9736E-01_r8,  0.1323E+00_r8,  0.1767E+00_r8,  0.2312E+00_r8,   &
         0.2962E+00_r8,  0.3728E+00_r8,  0.4619E+00_r8,  0.5623E+00_r8,  0.6685E+00_r8,   &
         0.7716E+00_r8, & 
                                ! DATA ((h22(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=18_r8,18)_r8, &                                &
         -0.4064E-07_r8, -0.8066E-07_r8, -0.1593E-06_r8, -0.3124E-06_r8, -0.6049E-06_r8,  &
         -0.1148E-05_r8, -0.2118E-05_r8, -0.3751E-05_r8, -0.6314E-05_r8, -0.1006E-04_r8,  &
         -0.1526E-04_r8, -0.2232E-04_r8, -0.3196E-04_r8, -0.4525E-04_r8, -0.6394E-04_r8,  &
         -0.9058E-04_r8, -0.1284E-03_r8, -0.1812E-03_r8, -0.2533E-03_r8, -0.3493E-03_r8,  &
         -0.4740E-03_r8, -0.6315E-03_r8, -0.8228E-03_r8, -0.1044E-02_r8, -0.1286E-02_r8,  &
         -0.1544E-02_r8, -0.1810E-02_r8, -0.2061E-02_r8, -0.2243E-02_r8, -0.2291E-02_r8,  &
         -0.2152E-02_r8, & 
                                ! DATA ((h23(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=18_r8,18)_r8, &                                &
         -0.1630E-09_r8, -0.3213E-09_r8, -0.6266E-09_r8, -0.1201E-08_r8, -0.2248E-08_r8,  &
         -0.4030E-08_r8, -0.6770E-08_r8, -0.1033E-07_r8, -0.1392E-07_r8, -0.1640E-07_r8,  &
         -0.1768E-07_r8, -0.1932E-07_r8, -0.2229E-07_r8, -0.2508E-07_r8, -0.2940E-07_r8,  &
         -0.4200E-07_r8, -0.6717E-07_r8, -0.1002E-06_r8, -0.1286E-06_r8, -0.1402E-06_r8,  &
         -0.1216E-06_r8, -0.5487E-07_r8,  0.8418E-07_r8,  0.3246E-06_r8,  0.6610E-06_r8,  &
         0.1013E-05_r8,  0.1394E-05_r8,  0.2073E-05_r8,  0.3337E-05_r8,  0.5175E-05_r8,   &
         0.7255E-05_r8, & 
                                ! DATA ((h21(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=19_r8,19)_r8, &                                &
         0.3282E-05_r8,  0.6539E-05_r8,  0.1300E-04_r8,  0.2578E-04_r8,  0.5085E-04_r8,   &
         0.9936E-04_r8,  0.1914E-03_r8,  0.3608E-03_r8,  0.6586E-03_r8,  0.1153E-02_r8,   &
         0.1924E-02_r8,  0.3074E-02_r8,  0.4742E-02_r8,  0.7126E-02_r8,  0.1051E-01_r8,   &
         0.1526E-01_r8,  0.2195E-01_r8,  0.3133E-01_r8,  0.4441E-01_r8,  0.6238E-01_r8,   &
         0.8664E-01_r8,  0.1186E+00_r8,  0.1597E+00_r8,  0.2108E+00_r8,  0.2723E+00_r8,   &
         0.3450E+00_r8,  0.4300E+00_r8,  0.5272E+00_r8,  0.6324E+00_r8,  0.7378E+00_r8,   &
         0.8330E+00_r8, & 
                                ! DATA ((h22(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=19_r8,19)_r8, &                                &
         -0.4629E-07_r8, -0.9195E-07_r8, -0.1819E-06_r8, -0.3572E-06_r8, -0.6936E-06_r8,  &
         -0.1323E-05_r8, -0.2456E-05_r8, -0.4385E-05_r8, -0.7453E-05_r8, -0.1200E-04_r8,  &
         -0.1843E-04_r8, -0.2731E-04_r8, -0.3943E-04_r8, -0.5606E-04_r8, -0.7936E-04_r8,  &
         -0.1123E-03_r8, -0.1588E-03_r8, -0.2231E-03_r8, -0.3101E-03_r8, -0.4247E-03_r8,  &
         -0.5713E-03_r8, -0.7522E-03_r8, -0.9651E-03_r8, -0.1202E-02_r8, -0.1456E-02_r8,  &
         -0.1721E-02_r8, -0.1983E-02_r8, -0.2196E-02_r8, -0.2296E-02_r8, -0.2224E-02_r8,  &
         -0.1952E-02_r8, & 
                                ! DATA ((h23(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=19_r8,19)_r8, &                                &
         -0.1827E-09_r8, -0.3607E-09_r8, -0.7057E-09_r8, -0.1359E-08_r8, -0.2552E-08_r8,  &
         -0.4615E-08_r8, -0.7854E-08_r8, -0.1218E-07_r8, -0.1670E-07_r8, -0.2008E-07_r8,  &
         -0.2241E-07_r8, -0.2516E-07_r8, -0.2796E-07_r8, -0.3015E-07_r8, -0.3506E-07_r8,  &
         -0.4958E-07_r8, -0.7627E-07_r8, -0.1070E-06_r8, -0.1289E-06_r8, -0.1286E-06_r8,  &
         -0.8843E-07_r8,  0.1492E-07_r8,  0.2118E-06_r8,  0.5155E-06_r8,  0.8665E-06_r8,  &
         0.1220E-05_r8,  0.1765E-05_r8,  0.2825E-05_r8,  0.4498E-05_r8,  0.6563E-05_r8,   &
         0.8422E-05_r8, & 
                                ! DATA ((h21(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=20_r8,20)_r8, &                                &
         0.3608E-05_r8,  0.7189E-05_r8,  0.1430E-04_r8,  0.2837E-04_r8,  0.5601E-04_r8,   &
         0.1097E-03_r8,  0.2120E-03_r8,  0.4017E-03_r8,  0.7396E-03_r8,  0.1311E-02_r8,   &
         0.2226E-02_r8,  0.3624E-02_r8,  0.5685E-02_r8,  0.8652E-02_r8,  0.1286E-01_r8,   &
         0.1878E-01_r8,  0.2709E-01_r8,  0.3868E-01_r8,  0.5471E-01_r8,  0.7654E-01_r8,   &
         0.1057E+00_r8,  0.1435E+00_r8,  0.1912E+00_r8,  0.2492E+00_r8,  0.3181E+00_r8,   &
         0.3990E+00_r8,  0.4926E+00_r8,  0.5960E+00_r8,  0.7027E+00_r8,  0.8027E+00_r8,   &
         0.8852E+00_r8, & 
                                ! DATA ((h22(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=20_r8,20)_r8, &                                &
         -0.5164E-07_r8, -0.1026E-06_r8, -0.2031E-06_r8, -0.3994E-06_r8, -0.7771E-06_r8,  &
         -0.1488E-05_r8, -0.2776E-05_r8, -0.5001E-05_r8, -0.8610E-05_r8, -0.1411E-04_r8,  &
         -0.2209E-04_r8, -0.3328E-04_r8, -0.4860E-04_r8, -0.6954E-04_r8, -0.9861E-04_r8,  &
         -0.1393E-03_r8, -0.1961E-03_r8, -0.2738E-03_r8, -0.3778E-03_r8, -0.5132E-03_r8,  &
         -0.6831E-03_r8, -0.8868E-03_r8, -0.1118E-02_r8, -0.1368E-02_r8, -0.1632E-02_r8,  &
         -0.1899E-02_r8, -0.2136E-02_r8, -0.2282E-02_r8, -0.2273E-02_r8, -0.2067E-02_r8,  &
         -0.1679E-02_r8, & 
                                ! DATA ((h23(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=20_r8,20)_r8, &                                &
         -0.2058E-09_r8, -0.4066E-09_r8, -0.7967E-09_r8, -0.1539E-08_r8, -0.2904E-08_r8,  &
         -0.5293E-08_r8, -0.9116E-08_r8, -0.1447E-07_r8, -0.2058E-07_r8, -0.2608E-07_r8,  &
         -0.3053E-07_r8, -0.3418E-07_r8, -0.3619E-07_r8, -0.3766E-07_r8, -0.4313E-07_r8,  &
         -0.5817E-07_r8, -0.8299E-07_r8, -0.1072E-06_r8, -0.1195E-06_r8, -0.1031E-06_r8,  &
         -0.3275E-07_r8,  0.1215E-06_r8,  0.3835E-06_r8,  0.7220E-06_r8,  0.1062E-05_r8,  &
         0.1504E-05_r8,  0.2367E-05_r8,  0.3854E-05_r8,  0.5842E-05_r8,  0.7875E-05_r8,   &
         0.9082E-05_r8, & 
                                ! DATA ((h21(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=21_r8,21)_r8, &                                &
         0.3827E-05_r8,  0.7626E-05_r8,  0.1518E-04_r8,  0.3013E-04_r8,  0.5955E-04_r8,   &
         0.1169E-03_r8,  0.2268E-03_r8,  0.4325E-03_r8,  0.8046E-03_r8,  0.1449E-02_r8,   &
         0.2512E-02_r8,  0.4185E-02_r8,  0.6707E-02_r8,  0.1038E-01_r8,  0.1561E-01_r8,   &
         0.2298E-01_r8,  0.3328E-01_r8,  0.4754E-01_r8,  0.6708E-01_r8,  0.9340E-01_r8,   &
         0.1280E+00_r8,  0.1722E+00_r8,  0.2267E+00_r8,  0.2918E+00_r8,  0.3688E+00_r8,   &
         0.4584E+00_r8,  0.5594E+00_r8,  0.6662E+00_r8,  0.7700E+00_r8,  0.8597E+00_r8,   &
         0.9266E+00_r8, & 
                                ! DATA ((h22(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=21_r8,21)_r8, &                                &
         -0.5584E-07_r8, -0.1110E-06_r8, -0.2198E-06_r8, -0.4329E-06_r8, -0.8444E-06_r8,  &
         -0.1623E-05_r8, -0.3049E-05_r8, -0.5551E-05_r8, -0.9714E-05_r8, -0.1627E-04_r8,  &
         -0.2609E-04_r8, -0.4015E-04_r8, -0.5955E-04_r8, -0.8603E-04_r8, -0.1223E-03_r8,  &
         -0.1724E-03_r8, -0.2413E-03_r8, -0.3346E-03_r8, -0.4578E-03_r8, -0.6155E-03_r8,  &
         -0.8087E-03_r8, -0.1033E-02_r8, -0.1279E-02_r8, -0.1540E-02_r8, -0.1811E-02_r8,  &
         -0.2065E-02_r8, -0.2251E-02_r8, -0.2301E-02_r8, -0.2163E-02_r8, -0.1828E-02_r8,  &
         -0.1365E-02_r8, & 
                                ! DATA ((h23(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=21_r8,21)_r8, &                                &
         -0.2274E-09_r8, -0.4498E-09_r8, -0.8814E-09_r8, -0.1708E-08_r8, -0.3247E-08_r8,  &
         -0.5972E-08_r8, -0.1045E-07_r8, -0.1707E-07_r8, -0.2545E-07_r8, -0.3440E-07_r8,  &
         -0.4259E-07_r8, -0.4822E-07_r8, -0.5004E-07_r8, -0.5061E-07_r8, -0.5485E-07_r8,  &
         -0.6687E-07_r8, -0.8483E-07_r8, -0.9896E-07_r8, -0.9646E-07_r8, -0.5557E-07_r8,  &
         0.5765E-07_r8,  0.2752E-06_r8,  0.5870E-06_r8,  0.9188E-06_r8,  0.1291E-05_r8,   &
         0.1971E-05_r8,  0.3251E-05_r8,  0.5115E-05_r8,  0.7221E-05_r8,  0.8825E-05_r8,   &
         0.9032E-05_r8, & 
                                ! DATA ((h21(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=22_r8,22)_r8, &                                &
         0.3938E-05_r8,  0.7849E-05_r8,  0.1563E-04_r8,  0.3105E-04_r8,  0.6149E-04_r8,   &
         0.1210E-03_r8,  0.2361E-03_r8,  0.4537E-03_r8,  0.8543E-03_r8,  0.1565E-02_r8,   &
         0.2774E-02_r8,  0.4739E-02_r8,  0.7781E-02_r8,  0.1229E-01_r8,  0.1877E-01_r8,   &
         0.2791E-01_r8,  0.4064E-01_r8,  0.5812E-01_r8,  0.8180E-01_r8,  0.1132E+00_r8,   &
         0.1539E+00_r8,  0.2047E+00_r8,  0.2662E+00_r8,  0.3391E+00_r8,  0.4247E+00_r8,   &
         0.5226E+00_r8,  0.6287E+00_r8,  0.7351E+00_r8,  0.8311E+00_r8,  0.9066E+00_r8,   &
         0.9569E+00_r8, & 
                                ! DATA ((h22(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=22_r8,22)_r8, &                                &
         -0.5833E-07_r8, -0.1160E-06_r8, -0.2300E-06_r8, -0.4540E-06_r8, -0.8885E-06_r8,  &
         -0.1718E-05_r8, -0.3256E-05_r8, -0.6010E-05_r8, -0.1072E-04_r8, -0.1838E-04_r8,  &
         -0.3026E-04_r8, -0.4772E-04_r8, -0.7223E-04_r8, -0.1057E-03_r8, -0.1512E-03_r8,  &
         -0.2128E-03_r8, -0.2961E-03_r8, -0.4070E-03_r8, -0.5514E-03_r8, -0.7320E-03_r8,  &
         -0.9467E-03_r8, -0.1187E-02_r8, -0.1446E-02_r8, -0.1717E-02_r8, -0.1985E-02_r8,  &
         -0.2203E-02_r8, -0.2307E-02_r8, -0.2237E-02_r8, -0.1965E-02_r8, -0.1532E-02_r8,  &
         -0.1044E-02_r8, & 
                                ! DATA ((h23(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=22_r8,22)_r8, &                                &
         -0.2426E-09_r8, -0.4805E-09_r8, -0.9447E-09_r8, -0.1841E-08_r8, -0.3519E-08_r8,  &
         -0.6565E-08_r8, -0.1172E-07_r8, -0.1979E-07_r8, -0.3095E-07_r8, -0.4443E-07_r8,  &
         -0.5821E-07_r8, -0.6868E-07_r8, -0.7282E-07_r8, -0.7208E-07_r8, -0.7176E-07_r8,  &
         -0.7562E-07_r8, -0.8110E-07_r8, -0.7934E-07_r8, -0.5365E-07_r8,  0.2483E-07_r8,  &
         0.1959E-06_r8,  0.4731E-06_r8,  0.7954E-06_r8,  0.1123E-05_r8,  0.1652E-05_r8,   &
         0.2711E-05_r8,  0.4402E-05_r8,  0.6498E-05_r8,  0.8392E-05_r8,  0.9154E-05_r8,   &
         0.8261E-05_r8, & 
                                ! DATA ((h21(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=23_r8,23)_r8, &                                &
         0.3973E-05_r8,  0.7921E-05_r8,  0.1578E-04_r8,  0.3139E-04_r8,  0.6229E-04_r8,   &
         0.1231E-03_r8,  0.2414E-03_r8,  0.4680E-03_r8,  0.8923E-03_r8,  0.1663E-02_r8,   &
         0.3011E-02_r8,  0.5268E-02_r8,  0.8869E-02_r8,  0.1434E-01_r8,  0.2230E-01_r8,   &
         0.3359E-01_r8,  0.4924E-01_r8,  0.7057E-01_r8,  0.9908E-01_r8,  0.1363E+00_r8,   &
         0.1834E+00_r8,  0.2411E+00_r8,  0.3100E+00_r8,  0.3913E+00_r8,  0.4856E+00_r8,   &
         0.5901E+00_r8,  0.6981E+00_r8,  0.7994E+00_r8,  0.8832E+00_r8,  0.9424E+00_r8,   &
         0.9773E+00_r8, & 
                                ! DATA ((h22(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=23_r8,23)_r8, &                                &
         -0.5929E-07_r8, -0.1180E-06_r8, -0.2344E-06_r8, -0.4638E-06_r8, -0.9118E-06_r8,  &
         -0.1775E-05_r8, -0.3401E-05_r8, -0.6375E-05_r8, -0.1160E-04_r8, -0.2039E-04_r8,  &
         -0.3444E-04_r8, -0.5575E-04_r8, -0.8641E-04_r8, -0.1287E-03_r8, -0.1856E-03_r8,  &
         -0.2615E-03_r8, -0.3618E-03_r8, -0.4928E-03_r8, -0.6594E-03_r8, -0.8621E-03_r8,  &
         -0.1095E-02_r8, -0.1349E-02_r8, -0.1618E-02_r8, -0.1894E-02_r8, -0.2140E-02_r8,  &
         -0.2293E-02_r8, -0.2290E-02_r8, -0.2085E-02_r8, -0.1696E-02_r8, -0.1212E-02_r8,  &
         -0.7506E-03_r8, & 
                                ! DATA ((h23(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=23_r8,23)_r8, &                                &
         -0.2496E-09_r8, -0.4954E-09_r8, -0.9780E-09_r8, -0.1915E-08_r8, -0.3697E-08_r8,  &
         -0.6991E-08_r8, -0.1279E-07_r8, -0.2231E-07_r8, -0.3653E-07_r8, -0.5541E-07_r8,  &
         -0.7688E-07_r8, -0.9614E-07_r8, -0.1069E-06_r8, -0.1065E-06_r8, -0.9866E-07_r8,  &
         -0.8740E-07_r8, -0.7192E-07_r8, -0.4304E-07_r8,  0.1982E-07_r8,  0.1525E-06_r8,  &
         0.3873E-06_r8,  0.6947E-06_r8,  0.1000E-05_r8,  0.1409E-05_r8,  0.2253E-05_r8,   &
         0.3739E-05_r8,  0.5744E-05_r8,  0.7812E-05_r8,  0.9067E-05_r8,  0.8746E-05_r8,   &
         0.6940E-05_r8, & 
                                ! DATA ((h21(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=24_r8,24)_r8, &                                &
         0.3978E-05_r8,  0.7933E-05_r8,  0.1581E-04_r8,  0.3149E-04_r8,  0.6260E-04_r8,   &
         0.1240E-03_r8,  0.2445E-03_r8,  0.4777E-03_r8,  0.9213E-03_r8,  0.1743E-02_r8,   &
         0.3217E-02_r8,  0.5756E-02_r8,  0.9929E-02_r8,  0.1644E-01_r8,  0.2612E-01_r8,   &
         0.3995E-01_r8,  0.5910E-01_r8,  0.8497E-01_r8,  0.1190E+00_r8,  0.1626E+00_r8,   &
         0.2165E+00_r8,  0.2813E+00_r8,  0.3583E+00_r8,  0.4485E+00_r8,  0.5505E+00_r8,   &
         0.6590E+00_r8,  0.7646E+00_r8,  0.8562E+00_r8,  0.9245E+00_r8,  0.9676E+00_r8,   &
         0.9897E+00_r8, & 
                                ! DATA ((h22(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=24_r8,24)_r8, &                                &
         -0.5950E-07_r8, -0.1185E-06_r8, -0.2358E-06_r8, -0.4678E-06_r8, -0.9235E-06_r8,  &
         -0.1810E-05_r8, -0.3503E-05_r8, -0.6664E-05_r8, -0.1236E-04_r8, -0.2223E-04_r8,  &
         -0.3849E-04_r8, -0.6396E-04_r8, -0.1017E-03_r8, -0.1545E-03_r8, -0.2257E-03_r8,  &
         -0.3192E-03_r8, -0.4399E-03_r8, -0.5933E-03_r8, -0.7824E-03_r8, -0.1005E-02_r8,  &
         -0.1251E-02_r8, -0.1516E-02_r8, -0.1794E-02_r8, -0.2060E-02_r8, -0.2257E-02_r8,  &
         -0.2318E-02_r8, -0.2186E-02_r8, -0.1853E-02_r8, -0.1386E-02_r8, -0.9021E-03_r8,  &
         -0.5050E-03_r8, & 
                                ! DATA ((h23(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=24_r8,24)_r8, &                                &
         -0.2515E-09_r8, -0.5001E-09_r8, -0.9904E-09_r8, -0.1951E-08_r8, -0.3800E-08_r8,  &
         -0.7288E-08_r8, -0.1362E-07_r8, -0.2452E-07_r8, -0.4184E-07_r8, -0.6663E-07_r8,  &
         -0.9770E-07_r8, -0.1299E-06_r8, -0.1533E-06_r8, -0.1584E-06_r8, -0.1425E-06_r8,  &
         -0.1093E-06_r8, -0.5972E-07_r8,  0.1426E-07_r8,  0.1347E-06_r8,  0.3364E-06_r8,  &
         0.6209E-06_r8,  0.9169E-06_r8,  0.1243E-05_r8,  0.1883E-05_r8,  0.3145E-05_r8,   &
         0.5011E-05_r8,  0.7136E-05_r8,  0.8785E-05_r8,  0.9048E-05_r8,  0.7686E-05_r8,   &
         0.5368E-05_r8, & 
                                ! DATA ((h21(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=25_r8,25)_r8, &                                &
         0.3978E-05_r8,  0.7934E-05_r8,  0.1582E-04_r8,  0.3153E-04_r8,  0.6275E-04_r8,   &
         0.1246E-03_r8,  0.2466E-03_r8,  0.4847E-03_r8,  0.9432E-03_r8,  0.1808E-02_r8,   &
         0.3392E-02_r8,  0.6189E-02_r8,  0.1092E-01_r8,  0.1852E-01_r8,  0.3009E-01_r8,   &
         0.4683E-01_r8,  0.7008E-01_r8,  0.1012E+00_r8,  0.1416E+00_r8,  0.1920E+00_r8,   &
         0.2531E+00_r8,  0.3258E+00_r8,  0.4115E+00_r8,  0.5102E+00_r8,  0.6180E+00_r8,   &
         0.7267E+00_r8,  0.8252E+00_r8,  0.9030E+00_r8,  0.9550E+00_r8,  0.9838E+00_r8,   &
         0.9965E+00_r8, & 
                                ! DATA ((h22(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=25_r8,25)_r8, &                                &
         -0.5953E-07_r8, -0.1187E-06_r8, -0.2363E-06_r8, -0.4697E-06_r8, -0.9304E-06_r8,  &
         -0.1833E-05_r8, -0.3578E-05_r8, -0.6889E-05_r8, -0.1299E-04_r8, -0.2384E-04_r8,  &
         -0.4227E-04_r8, -0.7201E-04_r8, -0.1174E-03_r8, -0.1825E-03_r8, -0.2710E-03_r8,  &
         -0.3861E-03_r8, -0.5313E-03_r8, -0.7095E-03_r8, -0.9203E-03_r8, -0.1158E-02_r8,  &
         -0.1417E-02_r8, -0.1691E-02_r8, -0.1968E-02_r8, -0.2200E-02_r8, -0.2320E-02_r8,  &
         -0.2263E-02_r8, -0.1998E-02_r8, -0.1563E-02_r8, -0.1068E-02_r8, -0.6300E-03_r8,  &
         -0.3174E-03_r8, & 
                                ! DATA ((h23(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=25_r8,25)_r8, &                                &
         -0.2520E-09_r8, -0.5016E-09_r8, -0.9963E-09_r8, -0.1971E-08_r8, -0.3867E-08_r8,  &
         -0.7500E-08_r8, -0.1427E-07_r8, -0.2634E-07_r8, -0.4656E-07_r8, -0.7753E-07_r8,  &
         -0.1196E-06_r8, -0.1683E-06_r8, -0.2113E-06_r8, -0.2309E-06_r8, -0.2119E-06_r8,  &
         -0.1518E-06_r8, -0.5200E-07_r8,  0.9274E-07_r8,  0.2973E-06_r8,  0.5714E-06_r8,  &
         0.8687E-06_r8,  0.1152E-05_r8,  0.1621E-05_r8,  0.2634E-05_r8,  0.4310E-05_r8,   &
         0.6418E-05_r8,  0.8347E-05_r8,  0.9162E-05_r8,  0.8319E-05_r8,  0.6209E-05_r8,   &
         0.3844E-05_r8, & 
                                ! DATA ((h21(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=26_r8,26)_r8, &                                &
         0.3976E-05_r8,  0.7932E-05_r8,  0.1582E-04_r8,  0.3154E-04_r8,  0.6284E-04_r8,   &
         0.1250E-03_r8,  0.2479E-03_r8,  0.4896E-03_r8,  0.9592E-03_r8,  0.1857E-02_r8,   &
         0.3532E-02_r8,  0.6558E-02_r8,  0.1181E-01_r8,  0.2048E-01_r8,  0.3402E-01_r8,   &
         0.5399E-01_r8,  0.8188E-01_r8,  0.1190E+00_r8,  0.1664E+00_r8,  0.2243E+00_r8,   &
         0.2933E+00_r8,  0.3747E+00_r8,  0.4695E+00_r8,  0.5755E+00_r8,  0.6858E+00_r8,   &
         0.7903E+00_r8,  0.8773E+00_r8,  0.9391E+00_r8,  0.9758E+00_r8,  0.9934E+00_r8,   &
         0.9996E+00_r8, & 
                                ! DATA ((h22(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=26_r8,26)_r8, &                                &
         -0.5954E-07_r8, -0.1187E-06_r8, -0.2366E-06_r8, -0.4709E-06_r8, -0.9349E-06_r8,  &
         -0.1849E-05_r8, -0.3632E-05_r8, -0.7058E-05_r8, -0.1350E-04_r8, -0.2521E-04_r8,  &
         -0.4564E-04_r8, -0.7958E-04_r8, -0.1329E-03_r8, -0.2114E-03_r8, -0.3200E-03_r8,  &
         -0.4611E-03_r8, -0.6353E-03_r8, -0.8409E-03_r8, -0.1072E-02_r8, -0.1324E-02_r8,  &
         -0.1592E-02_r8, -0.1871E-02_r8, -0.2127E-02_r8, -0.2297E-02_r8, -0.2312E-02_r8,  &
         -0.2123E-02_r8, -0.1738E-02_r8, -0.1247E-02_r8, -0.7744E-03_r8, -0.4117E-03_r8,  &
         -0.1850E-03_r8, & 
                                ! DATA ((h23(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=26_r8,26)_r8, &                                &
         -0.2522E-09_r8, -0.5025E-09_r8, -0.9997E-09_r8, -0.1983E-08_r8, -0.3912E-08_r8,  &
         -0.7650E-08_r8, -0.1474E-07_r8, -0.2777E-07_r8, -0.5055E-07_r8, -0.8745E-07_r8,  &
         -0.1414E-06_r8, -0.2095E-06_r8, -0.2790E-06_r8, -0.3241E-06_r8, -0.3135E-06_r8,  &
         -0.2269E-06_r8, -0.5896E-07_r8,  0.1875E-06_r8,  0.4996E-06_r8,  0.8299E-06_r8,  &
         0.1115E-05_r8,  0.1467E-05_r8,  0.2236E-05_r8,  0.3672E-05_r8,  0.5668E-05_r8,   &
         0.7772E-05_r8,  0.9094E-05_r8,  0.8827E-05_r8,  0.7041E-05_r8,  0.4638E-05_r8,   &
         0.2539E-05_r8/),SHAPE=(/nx*3*nh/))

    it=0
    DO k=1,nx
       DO i=1,3
          DO j=1,nh
             it=it+1              
             IF(i==1) THEN
                !WRITE(*,'(a5,2e17.9) ' )'h21   ',data4(it),h21(k,j)
                h21(k,j,1) = data4(it)
             END IF
             IF(i==2) THEN
                !WRITE(*,'(a5,2e17.9) ' )'h22   ',data4(it),h22(k,j)
                h22(k,j,1) =data4(it)
             END IF
             IF(i==3) THEN
                !WRITE(*,'(a5,2e17.9) ' )'h23   ',data4(it),h23(k,j)
                h23(k,j,1)=data4(it)
             END IF
          END DO
       END DO
    END DO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    data5(1:nx*3*nh)=RESHAPE(SOURCE=(/&
                                ! DATA ((h71(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 1_r8, 1)/                               &
         0.1342E-04_r8,  0.2640E-04_r8,  0.5137E-04_r8,  0.9829E-04_r8,  0.1832E-03_r8,   &
         0.3284E-03_r8,  0.5585E-03_r8,  0.8907E-03_r8,  0.1327E-02_r8,  0.1856E-02_r8,   &
         0.2461E-02_r8,  0.3138E-02_r8,  0.3898E-02_r8,  0.4766E-02_r8,  0.5783E-02_r8,   &
         0.7007E-02_r8,  0.8521E-02_r8,  0.1041E-01_r8,  0.1278E-01_r8,  0.1570E-01_r8,   &
         0.1928E-01_r8,  0.2361E-01_r8,  0.2884E-01_r8,  0.3520E-01_r8,  0.4305E-01_r8,   &
         0.5286E-01_r8,  0.6531E-01_r8,  0.8127E-01_r8,  0.1019E+00_r8,  0.1287E+00_r8,   &
         0.1632E+00_r8, &
                                ! DATA ((h72(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 1_r8, 1)_r8, &                               &
         -0.5685E-08_r8, -0.1331E-07_r8, -0.3249E-07_r8, -0.8137E-07_r8, -0.2048E-06_r8,  &
         -0.4973E-06_r8, -0.1118E-05_r8, -0.2246E-05_r8, -0.3982E-05_r8, -0.6290E-05_r8,  &
         -0.9040E-05_r8, -0.1215E-04_r8, -0.1567E-04_r8, -0.1970E-04_r8, -0.2449E-04_r8,  &
         -0.3046E-04_r8, -0.3798E-04_r8, -0.4725E-04_r8, -0.5831E-04_r8, -0.7123E-04_r8,  &
         -0.8605E-04_r8, -0.1028E-03_r8, -0.1212E-03_r8, -0.1413E-03_r8, -0.1635E-03_r8,  &
         -0.1884E-03_r8, -0.2160E-03_r8, -0.2461E-03_r8, -0.2778E-03_r8, -0.3098E-03_r8,  &
         -0.3411E-03_r8, &
                                ! DATA ((h73(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 1_r8, 1)_r8, &                               &
         0.2169E-10_r8,  0.5237E-10_r8,  0.1296E-09_r8,  0.3204E-09_r8,  0.7665E-09_r8,   &
         0.1691E-08_r8,  0.3222E-08_r8,  0.5110E-08_r8,  0.6779E-08_r8,  0.7681E-08_r8,   &
         0.7378E-08_r8,  0.5836E-08_r8,  0.3191E-08_r8, -0.1491E-08_r8, -0.1022E-07_r8,   &
         -0.2359E-07_r8, -0.3957E-07_r8, -0.5553E-07_r8, -0.6927E-07_r8, -0.7849E-07_r8,  &
         -0.8139E-07_r8, -0.7853E-07_r8, -0.7368E-07_r8, -0.7220E-07_r8, -0.7780E-07_r8,  &
         -0.9091E-07_r8, -0.1038E-06_r8, -0.9929E-07_r8, -0.5422E-07_r8,  0.5379E-07_r8,  &
         0.2350E-06_r8, &
                                ! DATA ((h71(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 2_r8, 2)_r8, &                               &
         0.1342E-04_r8,  0.2640E-04_r8,  0.5137E-04_r8,  0.9829E-04_r8,  0.1832E-03_r8,   &
         0.3284E-03_r8,  0.5586E-03_r8,  0.8908E-03_r8,  0.1327E-02_r8,  0.1856E-02_r8,   &
         0.2462E-02_r8,  0.3140E-02_r8,  0.3902E-02_r8,  0.4774E-02_r8,  0.5797E-02_r8,   &
         0.7033E-02_r8,  0.8570E-02_r8,  0.1050E-01_r8,  0.1294E-01_r8,  0.1599E-01_r8,   &
         0.1979E-01_r8,  0.2448E-01_r8,  0.3024E-01_r8,  0.3738E-01_r8,  0.4633E-01_r8,   &
         0.5766E-01_r8,  0.7218E-01_r8,  0.9096E-01_r8,  0.1154E+00_r8,  0.1471E+00_r8,   &
         0.1880E+00_r8, &
                                ! DATA ((h72(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 2_r8, 2)_r8, &                               &
         -0.5684E-08_r8, -0.1331E-07_r8, -0.3248E-07_r8, -0.8133E-07_r8, -0.2047E-06_r8,  &
         -0.4971E-06_r8, -0.1117E-05_r8, -0.2245E-05_r8, -0.3981E-05_r8, -0.6287E-05_r8,  &
         -0.9035E-05_r8, -0.1215E-04_r8, -0.1565E-04_r8, -0.1967E-04_r8, -0.2444E-04_r8,  &
         -0.3036E-04_r8, -0.3780E-04_r8, -0.4694E-04_r8, -0.5779E-04_r8, -0.7042E-04_r8,  &
         -0.8491E-04_r8, -0.1013E-03_r8, -0.1196E-03_r8, -0.1399E-03_r8, -0.1625E-03_r8,  &
         -0.1879E-03_r8, -0.2163E-03_r8, -0.2474E-03_r8, -0.2803E-03_r8, -0.3140E-03_r8,  &
         -0.3478E-03_r8, &
                                ! DATA ((h73(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 2_r8, 2)_r8, &                               &
         0.2168E-10_r8,  0.5242E-10_r8,  0.1295E-09_r8,  0.3201E-09_r8,  0.7662E-09_r8,   &
         0.1690E-08_r8,  0.3220E-08_r8,  0.5106E-08_r8,  0.6776E-08_r8,  0.7673E-08_r8,   &
         0.7362E-08_r8,  0.5808E-08_r8,  0.3138E-08_r8, -0.1595E-08_r8, -0.1041E-07_r8,   &
         -0.2390E-07_r8, -0.4010E-07_r8, -0.5636E-07_r8, -0.7045E-07_r8, -0.7972E-07_r8,  &
         -0.8178E-07_r8, -0.7677E-07_r8, -0.6876E-07_r8, -0.6381E-07_r8, -0.6583E-07_r8,  &
         -0.7486E-07_r8, -0.8229E-07_r8, -0.7017E-07_r8, -0.1497E-07_r8,  0.1051E-06_r8,  &
         0.2990E-06_r8, &
                                ! DATA ((h71(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 3_r8, 3)_r8, &                               &
         0.1342E-04_r8,  0.2640E-04_r8,  0.5137E-04_r8,  0.9830E-04_r8,  0.1832E-03_r8,   &
         0.3285E-03_r8,  0.5587E-03_r8,  0.8911E-03_r8,  0.1328E-02_r8,  0.1857E-02_r8,   &
         0.2465E-02_r8,  0.3144E-02_r8,  0.3909E-02_r8,  0.4787E-02_r8,  0.5820E-02_r8,   &
         0.7075E-02_r8,  0.8646E-02_r8,  0.1064E-01_r8,  0.1319E-01_r8,  0.1643E-01_r8,   &
         0.2053E-01_r8,  0.2569E-01_r8,  0.3215E-01_r8,  0.4028E-01_r8,  0.5059E-01_r8,   &
         0.6380E-01_r8,  0.8088E-01_r8,  0.1031E+00_r8,  0.1321E+00_r8,  0.1698E+00_r8,   &
         0.2179E+00_r8, &
                                ! DATA ((h72(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 3_r8, 3)_r8, &                               &
         -0.5682E-08_r8, -0.1330E-07_r8, -0.3247E-07_r8, -0.8129E-07_r8, -0.2046E-06_r8,  &
         -0.4968E-06_r8, -0.1117E-05_r8, -0.2244E-05_r8, -0.3978E-05_r8, -0.6283E-05_r8,  &
         -0.9027E-05_r8, -0.1213E-04_r8, -0.1563E-04_r8, -0.1963E-04_r8, -0.2436E-04_r8,  &
         -0.3021E-04_r8, -0.3754E-04_r8, -0.4649E-04_r8, -0.5709E-04_r8, -0.6940E-04_r8,  &
         -0.8359E-04_r8, -0.9986E-04_r8, -0.1182E-03_r8, -0.1388E-03_r8, -0.1620E-03_r8,  &
         -0.1882E-03_r8, -0.2175E-03_r8, -0.2498E-03_r8, -0.2843E-03_r8, -0.3203E-03_r8,  &
         -0.3573E-03_r8, &
                                ! DATA ((h73(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 3_r8, 3)_r8, &                               &
         0.2167E-10_r8,  0.5238E-10_r8,  0.1294E-09_r8,  0.3198E-09_r8,  0.7656E-09_r8,   &
         0.1688E-08_r8,  0.3217E-08_r8,  0.5104E-08_r8,  0.6767E-08_r8,  0.7661E-08_r8,   &
         0.7337E-08_r8,  0.5764E-08_r8,  0.3051E-08_r8, -0.1752E-08_r8, -0.1068E-07_r8,   &
         -0.2436E-07_r8, -0.4081E-07_r8, -0.5740E-07_r8, -0.7165E-07_r8, -0.8046E-07_r8,  &
         -0.8082E-07_r8, -0.7289E-07_r8, -0.6141E-07_r8, -0.5294E-07_r8, -0.5134E-07_r8,  &
         -0.5552E-07_r8, -0.5609E-07_r8, -0.3464E-07_r8,  0.3275E-07_r8,  0.1669E-06_r8,  &
         0.3745E-06_r8, &
                                ! DATA ((h71(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 4_r8, 4)_r8, &                               &
         0.1342E-04_r8,  0.2640E-04_r8,  0.5138E-04_r8,  0.9830E-04_r8,  0.1832E-03_r8,   &
         0.3286E-03_r8,  0.5589E-03_r8,  0.8915E-03_r8,  0.1329E-02_r8,  0.1859E-02_r8,   &
         0.2468E-02_r8,  0.3150E-02_r8,  0.3920E-02_r8,  0.4806E-02_r8,  0.5855E-02_r8,   &
         0.7140E-02_r8,  0.8763E-02_r8,  0.1085E-01_r8,  0.1356E-01_r8,  0.1707E-01_r8,   &
         0.2158E-01_r8,  0.2736E-01_r8,  0.3471E-01_r8,  0.4407E-01_r8,  0.5608E-01_r8,   &
         0.7161E-01_r8,  0.9184E-01_r8,  0.1183E+00_r8,  0.1528E+00_r8,  0.1973E+00_r8,   &
         0.2537E+00_r8, &
                                ! DATA ((h72(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 4_r8, 4)_r8, &                               &
         -0.5680E-08_r8, -0.1329E-07_r8, -0.3243E-07_r8, -0.8121E-07_r8, -0.2044E-06_r8,  &
         -0.4963E-06_r8, -0.1115E-05_r8, -0.2242E-05_r8, -0.3974E-05_r8, -0.6276E-05_r8,  &
         -0.9015E-05_r8, -0.1211E-04_r8, -0.1559E-04_r8, -0.1956E-04_r8, -0.2423E-04_r8,  &
         -0.2999E-04_r8, -0.3716E-04_r8, -0.4588E-04_r8, -0.5618E-04_r8, -0.6818E-04_r8,  &
         -0.8218E-04_r8, -0.9847E-04_r8, -0.1171E-03_r8, -0.1382E-03_r8, -0.1621E-03_r8,  &
         -0.1892E-03_r8, -0.2197E-03_r8, -0.2535E-03_r8, -0.2902E-03_r8, -0.3293E-03_r8,  &
         -0.3700E-03_r8, &
                                ! DATA ((h73(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 4_r8, 4)_r8, &                               &
         0.2166E-10_r8,  0.5229E-10_r8,  0.1294E-09_r8,  0.3193E-09_r8,  0.7644E-09_r8,   &
         0.1686E-08_r8,  0.3213E-08_r8,  0.5092E-08_r8,  0.6753E-08_r8,  0.7640E-08_r8,   &
         0.7302E-08_r8,  0.5696E-08_r8,  0.2917E-08_r8, -0.1984E-08_r8, -0.1108E-07_r8,   &
         -0.2497E-07_r8, -0.4171E-07_r8, -0.5849E-07_r8, -0.7254E-07_r8, -0.8017E-07_r8,  &
         -0.7802E-07_r8, -0.6662E-07_r8, -0.5153E-07_r8, -0.3961E-07_r8, -0.3387E-07_r8,  &
         -0.3219E-07_r8, -0.2426E-07_r8,  0.8700E-08_r8,  0.9027E-07_r8,  0.2400E-06_r8,  &
         0.4623E-06_r8, &
                                ! DATA ((h71(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 5_r8, 5)_r8, &                               &
         0.1342E-04_r8,  0.2640E-04_r8,  0.5138E-04_r8,  0.9832E-04_r8,  0.1833E-03_r8,   &
         0.3287E-03_r8,  0.5591E-03_r8,  0.8922E-03_r8,  0.1330E-02_r8,  0.1862E-02_r8,   &
         0.2473E-02_r8,  0.3159E-02_r8,  0.3937E-02_r8,  0.4837E-02_r8,  0.5911E-02_r8,   &
         0.7240E-02_r8,  0.8943E-02_r8,  0.1117E-01_r8,  0.1411E-01_r8,  0.1798E-01_r8,   &
         0.2304E-01_r8,  0.2961E-01_r8,  0.3807E-01_r8,  0.4896E-01_r8,  0.6308E-01_r8,   &
         0.8148E-01_r8,  0.1056E+00_r8,  0.1371E+00_r8,  0.1780E+00_r8,  0.2303E+00_r8,   &
         0.2958E+00_r8, &
                                ! DATA ((h72(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 5_r8, 5)_r8, &                               &
         -0.5675E-08_r8, -0.1328E-07_r8, -0.3239E-07_r8, -0.8110E-07_r8, -0.2040E-06_r8,  &
         -0.4954E-06_r8, -0.1114E-05_r8, -0.2238E-05_r8, -0.3968E-05_r8, -0.6265E-05_r8,  &
         -0.8996E-05_r8, -0.1208E-04_r8, -0.1553E-04_r8, -0.1945E-04_r8, -0.2404E-04_r8,  &
         -0.2966E-04_r8, -0.3663E-04_r8, -0.4508E-04_r8, -0.5508E-04_r8, -0.6686E-04_r8,  &
         -0.8082E-04_r8, -0.9732E-04_r8, -0.1165E-03_r8, -0.1382E-03_r8, -0.1630E-03_r8,  &
         -0.1913E-03_r8, -0.2234E-03_r8, -0.2593E-03_r8, -0.2989E-03_r8, -0.3417E-03_r8,  &
         -0.3857E-03_r8, &
                                ! DATA ((h73(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 5_r8, 5)_r8, &                               &
         0.2163E-10_r8,  0.5209E-10_r8,  0.1291E-09_r8,  0.3186E-09_r8,  0.7626E-09_r8,   &
         0.1682E-08_r8,  0.3203E-08_r8,  0.5078E-08_r8,  0.6730E-08_r8,  0.7606E-08_r8,   &
         0.7246E-08_r8,  0.5592E-08_r8,  0.2735E-08_r8, -0.2325E-08_r8, -0.1162E-07_r8,   &
         -0.2576E-07_r8, -0.4268E-07_r8, -0.5938E-07_r8, -0.7262E-07_r8, -0.7827E-07_r8,  &
         -0.7297E-07_r8, -0.5786E-07_r8, -0.3930E-07_r8, -0.2373E-07_r8, -0.1295E-07_r8,  &
         -0.3728E-08_r8,  0.1465E-07_r8,  0.6114E-07_r8,  0.1590E-06_r8,  0.3257E-06_r8,  &
         0.5622E-06_r8, &
                                ! DATA ((h71(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 6_r8, 6)_r8, &                               &
         0.1342E-04_r8,  0.2640E-04_r8,  0.5139E-04_r8,  0.9834E-04_r8,  0.1833E-03_r8,   &
         0.3288E-03_r8,  0.5596E-03_r8,  0.8932E-03_r8,  0.1332E-02_r8,  0.1866E-02_r8,   &
         0.2481E-02_r8,  0.3174E-02_r8,  0.3963E-02_r8,  0.4885E-02_r8,  0.5997E-02_r8,   &
         0.7394E-02_r8,  0.9215E-02_r8,  0.1164E-01_r8,  0.1489E-01_r8,  0.1924E-01_r8,   &
         0.2501E-01_r8,  0.3258E-01_r8,  0.4243E-01_r8,  0.5523E-01_r8,  0.7195E-01_r8,   &
         0.9387E-01_r8,  0.1226E+00_r8,  0.1601E+00_r8,  0.2084E+00_r8,  0.2695E+00_r8,   &
         0.3446E+00_r8, &
                                ! DATA ((h72(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 6_r8, 6)_r8, &                               &
         -0.5671E-08_r8, -0.1326E-07_r8, -0.3234E-07_r8, -0.8091E-07_r8, -0.2035E-06_r8,  &
         -0.4941E-06_r8, -0.1111E-05_r8, -0.2232E-05_r8, -0.3958E-05_r8, -0.6247E-05_r8,  &
         -0.8966E-05_r8, -0.1202E-04_r8, -0.1544E-04_r8, -0.1929E-04_r8, -0.2377E-04_r8,  &
         -0.2921E-04_r8, -0.3593E-04_r8, -0.4409E-04_r8, -0.5385E-04_r8, -0.6555E-04_r8,  &
         -0.7965E-04_r8, -0.9656E-04_r8, -0.1163E-03_r8, -0.1390E-03_r8, -0.1649E-03_r8,  &
         -0.1947E-03_r8, -0.2288E-03_r8, -0.2675E-03_r8, -0.3109E-03_r8, -0.3575E-03_r8,  &
         -0.4039E-03_r8, &
                                ! DATA ((h73(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 6_r8, 6)_r8, &                               &
         0.2155E-10_r8,  0.5188E-10_r8,  0.1288E-09_r8,  0.3175E-09_r8,  0.7599E-09_r8,   &
         0.1675E-08_r8,  0.3190E-08_r8,  0.5059E-08_r8,  0.6699E-08_r8,  0.7551E-08_r8,   &
         0.7154E-08_r8,  0.5435E-08_r8,  0.2452E-08_r8, -0.2802E-08_r8, -0.1235E-07_r8,   &
         -0.2668E-07_r8, -0.4353E-07_r8, -0.5962E-07_r8, -0.7134E-07_r8, -0.7435E-07_r8,  &
         -0.6551E-07_r8, -0.4676E-07_r8, -0.2475E-07_r8, -0.4876E-08_r8,  0.1235E-07_r8,  &
         0.3092E-07_r8,  0.6192E-07_r8,  0.1243E-06_r8,  0.2400E-06_r8,  0.4247E-06_r8,   &
         0.6755E-06_r8, &
                                ! DATA ((h71(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 7_r8, 7)_r8, &                               &
         0.1342E-04_r8,  0.2640E-04_r8,  0.5139E-04_r8,  0.9837E-04_r8,  0.1834E-03_r8,   &
         0.3291E-03_r8,  0.5603E-03_r8,  0.8949E-03_r8,  0.1336E-02_r8,  0.1873E-02_r8,   &
         0.2494E-02_r8,  0.3198E-02_r8,  0.4005E-02_r8,  0.4960E-02_r8,  0.6130E-02_r8,   &
         0.7628E-02_r8,  0.9618E-02_r8,  0.1232E-01_r8,  0.1599E-01_r8,  0.2097E-01_r8,   &
         0.2764E-01_r8,  0.3646E-01_r8,  0.4803E-01_r8,  0.6319E-01_r8,  0.8312E-01_r8,   &
         0.1093E+00_r8,  0.1436E+00_r8,  0.1881E+00_r8,  0.2448E+00_r8,  0.3153E+00_r8,   &
         0.4001E+00_r8, &
                                ! DATA ((h72(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 7_r8, 7)_r8, &                               &
         -0.5665E-08_r8, -0.1322E-07_r8, -0.3224E-07_r8, -0.8063E-07_r8, -0.2027E-06_r8,  &
         -0.4921E-06_r8, -0.1106E-05_r8, -0.2223E-05_r8, -0.3942E-05_r8, -0.6220E-05_r8,  &
         -0.8920E-05_r8, -0.1194E-04_r8, -0.1530E-04_r8, -0.1905E-04_r8, -0.2337E-04_r8,  &
         -0.2860E-04_r8, -0.3505E-04_r8, -0.4296E-04_r8, -0.5259E-04_r8, -0.6439E-04_r8,  &
         -0.7884E-04_r8, -0.9635E-04_r8, -0.1170E-03_r8, -0.1407E-03_r8, -0.1681E-03_r8,  &
         -0.1998E-03_r8, -0.2366E-03_r8, -0.2790E-03_r8, -0.3265E-03_r8, -0.3763E-03_r8,  &
         -0.4235E-03_r8, &
                                ! DATA ((h73(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 7_r8, 7)_r8, &                               &
         0.2157E-10_r8,  0.5178E-10_r8,  0.1283E-09_r8,  0.3162E-09_r8,  0.7558E-09_r8,   &
         0.1665E-08_r8,  0.3169E-08_r8,  0.5027E-08_r8,  0.6645E-08_r8,  0.7472E-08_r8,   &
         0.7017E-08_r8,  0.5212E-08_r8,  0.2059E-08_r8, -0.3443E-08_r8, -0.1321E-07_r8,   &
         -0.2754E-07_r8, -0.4389E-07_r8, -0.5869E-07_r8, -0.6825E-07_r8, -0.6819E-07_r8,  &
         -0.5570E-07_r8, -0.3343E-07_r8, -0.7592E-08_r8,  0.1778E-07_r8,  0.4322E-07_r8,  &
         0.7330E-07_r8,  0.1190E-06_r8,  0.1993E-06_r8,  0.3348E-06_r8,  0.5376E-06_r8,   &
         0.8030E-06_r8, &
                                ! DATA ((h71(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 8_r8, 8)_r8, &                               &
         0.1342E-04_r8,  0.2641E-04_r8,  0.5141E-04_r8,  0.9842E-04_r8,  0.1836E-03_r8,   &
         0.3295E-03_r8,  0.5614E-03_r8,  0.8975E-03_r8,  0.1341E-02_r8,  0.1884E-02_r8,   &
         0.2514E-02_r8,  0.3234E-02_r8,  0.4071E-02_r8,  0.5075E-02_r8,  0.6332E-02_r8,   &
         0.7976E-02_r8,  0.1021E-01_r8,  0.1328E-01_r8,  0.1751E-01_r8,  0.2329E-01_r8,   &
         0.3109E-01_r8,  0.4146E-01_r8,  0.5517E-01_r8,  0.7325E-01_r8,  0.9710E-01_r8,   &
         0.1284E+00_r8,  0.1693E+00_r8,  0.2217E+00_r8,  0.2876E+00_r8,  0.3679E+00_r8,   &
         0.4616E+00_r8, &
                                ! DATA ((h72(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 8_r8, 8)_r8, &                               &
         -0.5652E-08_r8, -0.1318E-07_r8, -0.3210E-07_r8, -0.8018E-07_r8, -0.2014E-06_r8,  &
         -0.4888E-06_r8, -0.1099E-05_r8, -0.2210E-05_r8, -0.3918E-05_r8, -0.6179E-05_r8,  &
         -0.8849E-05_r8, -0.1182E-04_r8, -0.1509E-04_r8, -0.1871E-04_r8, -0.2284E-04_r8,  &
         -0.2782E-04_r8, -0.3403E-04_r8, -0.4177E-04_r8, -0.5145E-04_r8, -0.6354E-04_r8,  &
         -0.7853E-04_r8, -0.9681E-04_r8, -0.1185E-03_r8, -0.1437E-03_r8, -0.1729E-03_r8,  &
         -0.2072E-03_r8, -0.2475E-03_r8, -0.2942E-03_r8, -0.3457E-03_r8, -0.3973E-03_r8,  &
         -0.4434E-03_r8, &
                                ! DATA ((h73(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 8_r8, 8)_r8, &                               &
         0.2153E-10_r8,  0.5151E-10_r8,  0.1273E-09_r8,  0.3136E-09_r8,  0.7488E-09_r8,   &
         0.1649E-08_r8,  0.3142E-08_r8,  0.4980E-08_r8,  0.6559E-08_r8,  0.7346E-08_r8,   &
         0.6813E-08_r8,  0.4884E-08_r8,  0.1533E-08_r8, -0.4209E-08_r8, -0.1409E-07_r8,   &
         -0.2801E-07_r8, -0.4320E-07_r8, -0.5614E-07_r8, -0.6312E-07_r8, -0.5976E-07_r8,  &
         -0.4369E-07_r8, -0.1775E-07_r8,  0.1280E-07_r8,  0.4534E-07_r8,  0.8106E-07_r8,  &
         0.1246E-06_r8,  0.1874E-06_r8,  0.2873E-06_r8,  0.4433E-06_r8,  0.6651E-06_r8,   &
         0.9477E-06_r8, &
                                ! DATA ((h71(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 9_r8, 9)_r8, &                               &
         0.1342E-04_r8,  0.2641E-04_r8,  0.5143E-04_r8,  0.9849E-04_r8,  0.1838E-03_r8,   &
         0.3302E-03_r8,  0.5631E-03_r8,  0.9016E-03_r8,  0.1350E-02_r8,  0.1901E-02_r8,   &
         0.2546E-02_r8,  0.3292E-02_r8,  0.4171E-02_r8,  0.5250E-02_r8,  0.6634E-02_r8,   &
         0.8486E-02_r8,  0.1104E-01_r8,  0.1460E-01_r8,  0.1955E-01_r8,  0.2635E-01_r8,   &
         0.3555E-01_r8,  0.4786E-01_r8,  0.6421E-01_r8,  0.8588E-01_r8,  0.1145E+00_r8,   &
         0.1519E+00_r8,  0.2003E+00_r8,  0.2615E+00_r8,  0.3372E+00_r8,  0.4269E+00_r8,   &
         0.5281E+00_r8, &
                                ! DATA ((h72(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 9_r8, 9)_r8, &                               &
         -0.5629E-08_r8, -0.1310E-07_r8, -0.3186E-07_r8, -0.7948E-07_r8, -0.1995E-06_r8,  &
         -0.4837E-06_r8, -0.1088E-05_r8, -0.2188E-05_r8, -0.3880E-05_r8, -0.6115E-05_r8,  &
         -0.8743E-05_r8, -0.1165E-04_r8, -0.1480E-04_r8, -0.1824E-04_r8, -0.2216E-04_r8,  &
         -0.2691E-04_r8, -0.3293E-04_r8, -0.4067E-04_r8, -0.5057E-04_r8, -0.6314E-04_r8,  &
         -0.7885E-04_r8, -0.9813E-04_r8, -0.1212E-03_r8, -0.1482E-03_r8, -0.1799E-03_r8,  &
         -0.2175E-03_r8, -0.2622E-03_r8, -0.3135E-03_r8, -0.3678E-03_r8, -0.4193E-03_r8,  &
         -0.4627E-03_r8, &
                                ! DATA ((h73(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip= 9_r8, 9)_r8, &                               &
         0.2121E-10_r8,  0.5076E-10_r8,  0.1257E-09_r8,  0.3091E-09_r8,  0.7379E-09_r8,   &
         0.1623E-08_r8,  0.3097E-08_r8,  0.4904E-08_r8,  0.6453E-08_r8,  0.7168E-08_r8,   &
         0.6534E-08_r8,  0.4458E-08_r8,  0.8932E-09_r8, -0.5026E-08_r8, -0.1469E-07_r8,   &
         -0.2765E-07_r8, -0.4103E-07_r8, -0.5169E-07_r8, -0.5585E-07_r8, -0.4913E-07_r8,  &
         -0.2954E-07_r8,  0.6372E-09_r8,  0.3738E-07_r8,  0.7896E-07_r8,  0.1272E-06_r8,  &
         0.1867E-06_r8,  0.2682E-06_r8,  0.3895E-06_r8,  0.5672E-06_r8,  0.8091E-06_r8,   &
         0.1114E-05_r8, &
                                ! DATA ((h71(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=10_r8,10)_r8, &                               &
         0.1343E-04_r8,  0.2642E-04_r8,  0.5147E-04_r8,  0.9861E-04_r8,  0.1842E-03_r8,   &
         0.3312E-03_r8,  0.5658E-03_r8,  0.9080E-03_r8,  0.1364E-02_r8,  0.1928E-02_r8,   &
         0.2595E-02_r8,  0.3379E-02_r8,  0.4325E-02_r8,  0.5513E-02_r8,  0.7077E-02_r8,   &
         0.9215E-02_r8,  0.1220E-01_r8,  0.1640E-01_r8,  0.2226E-01_r8,  0.3032E-01_r8,   &
         0.4127E-01_r8,  0.5598E-01_r8,  0.7559E-01_r8,  0.1016E+00_r8,  0.1359E+00_r8,   &
         0.1804E+00_r8,  0.2372E+00_r8,  0.3081E+00_r8,  0.3935E+00_r8,  0.4916E+00_r8,   &
         0.5975E+00_r8, &
                                ! DATA ((h72(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=10_r8,10)_r8, &                               &
         -0.5597E-08_r8, -0.1300E-07_r8, -0.3148E-07_r8, -0.7838E-07_r8, -0.1964E-06_r8,  &
         -0.4759E-06_r8, -0.1071E-05_r8, -0.2155E-05_r8, -0.3822E-05_r8, -0.6019E-05_r8,  &
         -0.8586E-05_r8, -0.1139E-04_r8, -0.1439E-04_r8, -0.1764E-04_r8, -0.2134E-04_r8,  &
         -0.2591E-04_r8, -0.3188E-04_r8, -0.3978E-04_r8, -0.5011E-04_r8, -0.6334E-04_r8,  &
         -0.7998E-04_r8, -0.1006E-03_r8, -0.1253E-03_r8, -0.1547E-03_r8, -0.1895E-03_r8,  &
         -0.2315E-03_r8, -0.2811E-03_r8, -0.3363E-03_r8, -0.3917E-03_r8, -0.4413E-03_r8,  &
         -0.4809E-03_r8, &
                                ! DATA ((h73(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=10_r8,10)_r8, &                               &
         0.2109E-10_r8,  0.5017E-10_r8,  0.1235E-09_r8,  0.3021E-09_r8,  0.7217E-09_r8,   &
         0.1585E-08_r8,  0.3028E-08_r8,  0.4796E-08_r8,  0.6285E-08_r8,  0.6910E-08_r8,   &
         0.6178E-08_r8,  0.3945E-08_r8,  0.2436E-09_r8, -0.5632E-08_r8, -0.1464E-07_r8,   &
         -0.2596E-07_r8, -0.3707E-07_r8, -0.4527E-07_r8, -0.4651E-07_r8, -0.3644E-07_r8,  &
         -0.1296E-07_r8,  0.2250E-07_r8,  0.6722E-07_r8,  0.1202E-06_r8,  0.1831E-06_r8,  &
         0.2605E-06_r8,  0.3627E-06_r8,  0.5062E-06_r8,  0.7064E-06_r8,  0.9725E-06_r8,   &
         0.1304E-05_r8, &
                                ! DATA ((h71(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=11_r8,11)_r8, &                               &
         0.1343E-04_r8,  0.2644E-04_r8,  0.5152E-04_r8,  0.9878E-04_r8,  0.1847E-03_r8,   &
         0.3328E-03_r8,  0.5701E-03_r8,  0.9179E-03_r8,  0.1385E-02_r8,  0.1969E-02_r8,   &
         0.2671E-02_r8,  0.3513E-02_r8,  0.4554E-02_r8,  0.5900E-02_r8,  0.7713E-02_r8,   &
         0.1023E-01_r8,  0.1379E-01_r8,  0.1879E-01_r8,  0.2578E-01_r8,  0.3542E-01_r8,   &
         0.4854E-01_r8,  0.6622E-01_r8,  0.8983E-01_r8,  0.1211E+00_r8,  0.1619E+00_r8,   &
         0.2145E+00_r8,  0.2807E+00_r8,  0.3614E+00_r8,  0.4559E+00_r8,  0.5603E+00_r8,   &
         0.6675E+00_r8, &
                                ! DATA ((h72(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=11_r8,11)_r8, &                               &
         -0.5538E-08_r8, -0.1280E-07_r8, -0.3089E-07_r8, -0.7667E-07_r8, -0.1917E-06_r8,  &
         -0.4642E-06_r8, -0.1045E-05_r8, -0.2106E-05_r8, -0.3736E-05_r8, -0.5878E-05_r8,  &
         -0.8363E-05_r8, -0.1104E-04_r8, -0.1387E-04_r8, -0.1692E-04_r8, -0.2044E-04_r8,  &
         -0.2493E-04_r8, -0.3101E-04_r8, -0.3926E-04_r8, -0.5020E-04_r8, -0.6429E-04_r8,  &
         -0.8213E-04_r8, -0.1044E-03_r8, -0.1314E-03_r8, -0.1637E-03_r8, -0.2027E-03_r8,  &
         -0.2498E-03_r8, -0.3042E-03_r8, -0.3617E-03_r8, -0.4163E-03_r8, -0.4625E-03_r8,  &
         -0.4969E-03_r8, &
                                ! DATA ((h73(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=11_r8,11)_r8, &                               &
         0.2067E-10_r8,  0.4903E-10_r8,  0.1200E-09_r8,  0.2917E-09_r8,  0.6965E-09_r8,   &
         0.1532E-08_r8,  0.2925E-08_r8,  0.4632E-08_r8,  0.6054E-08_r8,  0.6590E-08_r8,   &
         0.5746E-08_r8,  0.3436E-08_r8, -0.2251E-09_r8, -0.5703E-08_r8, -0.1344E-07_r8,   &
         -0.2256E-07_r8, -0.3120E-07_r8, -0.3690E-07_r8, -0.3520E-07_r8, -0.2164E-07_r8,  &
         0.6510E-08_r8,  0.4895E-07_r8,  0.1037E-06_r8,  0.1702E-06_r8,  0.2502E-06_r8,   &
         0.3472E-06_r8,  0.4710E-06_r8,  0.6379E-06_r8,  0.8633E-06_r8,  0.1159E-05_r8,   &
         0.1514E-05_r8, &
                                ! DATA ((h71(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=12_r8,12)_r8, &                               &
         0.1344E-04_r8,  0.2647E-04_r8,  0.5161E-04_r8,  0.9906E-04_r8,  0.1856E-03_r8,   &
         0.3353E-03_r8,  0.5765E-03_r8,  0.9331E-03_r8,  0.1417E-02_r8,  0.2032E-02_r8,   &
         0.2785E-02_r8,  0.3712E-02_r8,  0.4892E-02_r8,  0.6456E-02_r8,  0.8604E-02_r8,   &
         0.1162E-01_r8,  0.1590E-01_r8,  0.2192E-01_r8,  0.3033E-01_r8,  0.4193E-01_r8,   &
         0.5774E-01_r8,  0.7906E-01_r8,  0.1075E+00_r8,  0.1449E+00_r8,  0.1934E+00_r8,   &
         0.2549E+00_r8,  0.3309E+00_r8,  0.4213E+00_r8,  0.5232E+00_r8,  0.6308E+00_r8,   &
         0.7348E+00_r8, &
                                ! DATA ((h72(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=12_r8,12)_r8, &                               &
         -0.5476E-08_r8, -0.1257E-07_r8, -0.3008E-07_r8, -0.7418E-07_r8, -0.1848E-06_r8,  &
         -0.4468E-06_r8, -0.1006E-05_r8, -0.2032E-05_r8, -0.3611E-05_r8, -0.5679E-05_r8,  &
         -0.8058E-05_r8, -0.1059E-04_r8, -0.1324E-04_r8, -0.1612E-04_r8, -0.1956E-04_r8,  &
         -0.2411E-04_r8, -0.3046E-04_r8, -0.3925E-04_r8, -0.5098E-04_r8, -0.6619E-04_r8,  &
         -0.8562E-04_r8, -0.1100E-03_r8, -0.1399E-03_r8, -0.1761E-03_r8, -0.2202E-03_r8,  &
         -0.2726E-03_r8, -0.3306E-03_r8, -0.3885E-03_r8, -0.4404E-03_r8, -0.4820E-03_r8,  &
         -0.5082E-03_r8, &
                                ! DATA ((h73(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=12_r8,12)_r8, &                               &
         0.2041E-10_r8,  0.4771E-10_r8,  0.1149E-09_r8,  0.2782E-09_r8,  0.6614E-09_r8,   &
         0.1451E-08_r8,  0.2778E-08_r8,  0.4401E-08_r8,  0.5736E-08_r8,  0.6189E-08_r8,   &
         0.5315E-08_r8,  0.3087E-08_r8, -0.2518E-09_r8, -0.4806E-08_r8, -0.1071E-07_r8,   &
         -0.1731E-07_r8, -0.2346E-07_r8, -0.2659E-07_r8, -0.2184E-07_r8, -0.4261E-08_r8,  &
         0.2975E-07_r8,  0.8112E-07_r8,  0.1484E-06_r8,  0.2308E-06_r8,  0.3296E-06_r8,   &
         0.4475E-06_r8,  0.5942E-06_r8,  0.7859E-06_r8,  0.1041E-05_r8,  0.1369E-05_r8,   &
         0.1726E-05_r8, &
                                ! DATA ((h71(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=13_r8,13)_r8, &                               &
         0.1345E-04_r8,  0.2650E-04_r8,  0.5173E-04_r8,  0.9947E-04_r8,  0.1868E-03_r8,   &
         0.3389E-03_r8,  0.5861E-03_r8,  0.9559E-03_r8,  0.1466E-02_r8,  0.2125E-02_r8,   &
         0.2955E-02_r8,  0.4004E-02_r8,  0.5376E-02_r8,  0.7235E-02_r8,  0.9824E-02_r8,   &
         0.1348E-01_r8,  0.1866E-01_r8,  0.2596E-01_r8,  0.3613E-01_r8,  0.5016E-01_r8,   &
         0.6929E-01_r8,  0.9504E-01_r8,  0.1292E+00_r8,  0.1738E+00_r8,  0.2309E+00_r8,   &
         0.3021E+00_r8,  0.3880E+00_r8,  0.4868E+00_r8,  0.5936E+00_r8,  0.7003E+00_r8,   &
         0.7964E+00_r8, &
                                ! DATA ((h72(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=13_r8,13)_r8, &                               &
         -0.5362E-08_r8, -0.1223E-07_r8, -0.2895E-07_r8, -0.7071E-07_r8, -0.1748E-06_r8,  &
         -0.4219E-06_r8, -0.9516E-06_r8, -0.1928E-05_r8, -0.3436E-05_r8, -0.5409E-05_r8,  &
         -0.7666E-05_r8, -0.1005E-04_r8, -0.1254E-04_r8, -0.1533E-04_r8, -0.1880E-04_r8,  &
         -0.2358E-04_r8, -0.3038E-04_r8, -0.3988E-04_r8, -0.5264E-04_r8, -0.6934E-04_r8,  &
         -0.9083E-04_r8, -0.1179E-03_r8, -0.1515E-03_r8, -0.1927E-03_r8, -0.2424E-03_r8,  &
         -0.2994E-03_r8, -0.3591E-03_r8, -0.4155E-03_r8, -0.4634E-03_r8, -0.4982E-03_r8,  &
         -0.5096E-03_r8, &
                                ! DATA ((h73(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=13_r8,13)_r8, &                               &
         0.1976E-10_r8,  0.4551E-10_r8,  0.1086E-09_r8,  0.2601E-09_r8,  0.6126E-09_r8,   &
         0.1345E-08_r8,  0.2583E-08_r8,  0.4112E-08_r8,  0.5365E-08_r8,  0.5796E-08_r8,   &
         0.5031E-08_r8,  0.3182E-08_r8,  0.5970E-09_r8, -0.2547E-08_r8, -0.6172E-08_r8,   &
         -0.1017E-07_r8, -0.1388E-07_r8, -0.1430E-07_r8, -0.6118E-08_r8,  0.1624E-07_r8,  &
         0.5791E-07_r8,  0.1205E-06_r8,  0.2025E-06_r8,  0.3032E-06_r8,  0.4225E-06_r8,   &
         0.5619E-06_r8,  0.7322E-06_r8,  0.9528E-06_r8,  0.1243E-05_r8,  0.1592E-05_r8,   &
         0.1904E-05_r8, &
                                ! DATA ((h71(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=14_r8,14)_r8, &                               &
         0.1346E-04_r8,  0.2655E-04_r8,  0.5192E-04_r8,  0.1001E-03_r8,  0.1887E-03_r8,   &
         0.3442E-03_r8,  0.6001E-03_r8,  0.9892E-03_r8,  0.1536E-02_r8,  0.2262E-02_r8,   &
         0.3200E-02_r8,  0.4420E-02_r8,  0.6053E-02_r8,  0.8300E-02_r8,  0.1146E-01_r8,   &
         0.1592E-01_r8,  0.2224E-01_r8,  0.3112E-01_r8,  0.4348E-01_r8,  0.6051E-01_r8,   &
         0.8369E-01_r8,  0.1147E+00_r8,  0.1556E+00_r8,  0.2084E+00_r8,  0.2749E+00_r8,   &
         0.3561E+00_r8,  0.4511E+00_r8,  0.5562E+00_r8,  0.6644E+00_r8,  0.7655E+00_r8,   &
         0.8499E+00_r8, &
                                ! DATA ((h72(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=14_r8,14)_r8, &                               &
         -0.5210E-08_r8, -0.1172E-07_r8, -0.2731E-07_r8, -0.6598E-07_r8, -0.1615E-06_r8,  &
         -0.3880E-06_r8, -0.8769E-06_r8, -0.1787E-05_r8, -0.3204E-05_r8, -0.5066E-05_r8,  &
         -0.7197E-05_r8, -0.9451E-05_r8, -0.1185E-04_r8, -0.1465E-04_r8, -0.1831E-04_r8,  &
         -0.2346E-04_r8, -0.3088E-04_r8, -0.4132E-04_r8, -0.5545E-04_r8, -0.7410E-04_r8,  &
         -0.9820E-04_r8, -0.1288E-03_r8, -0.1670E-03_r8, -0.2140E-03_r8, -0.2692E-03_r8,  &
         -0.3293E-03_r8, -0.3886E-03_r8, -0.4417E-03_r8, -0.4840E-03_r8, -0.5073E-03_r8,  &
         -0.4944E-03_r8, &
                                ! DATA ((h73(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=14_r8,14)_r8, &                               &
         0.1880E-10_r8,  0.4271E-10_r8,  0.9966E-10_r8,  0.2352E-09_r8,  0.5497E-09_r8,   &
         0.1205E-08_r8,  0.2334E-08_r8,  0.3765E-08_r8,  0.4993E-08_r8,  0.5532E-08_r8,   &
         0.5148E-08_r8,  0.4055E-08_r8,  0.2650E-08_r8,  0.1326E-08_r8,  0.2019E-09_r8,   &
         -0.1124E-08_r8, -0.2234E-08_r8,  0.2827E-09_r8,  0.1247E-07_r8,  0.4102E-07_r8,  &
         0.9228E-07_r8,  0.1682E-06_r8,  0.2676E-06_r8,  0.3885E-06_r8,  0.5286E-06_r8,   &
         0.6904E-06_r8,  0.8871E-06_r8,  0.1142E-05_r8,  0.1466E-05_r8,  0.1800E-05_r8,   &
         0.2004E-05_r8, &
                                ! DATA ((h71(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=15_r8,15)_r8, &                               &
         0.1348E-04_r8,  0.2663E-04_r8,  0.5216E-04_r8,  0.1009E-03_r8,  0.1912E-03_r8,   &
         0.3515E-03_r8,  0.6196E-03_r8,  0.1036E-02_r8,  0.1637E-02_r8,  0.2456E-02_r8,   &
         0.3544E-02_r8,  0.4996E-02_r8,  0.6973E-02_r8,  0.9722E-02_r8,  0.1360E-01_r8,   &
         0.1908E-01_r8,  0.2681E-01_r8,  0.3766E-01_r8,  0.5272E-01_r8,  0.7343E-01_r8,   &
         0.1015E+00_r8,  0.1388E+00_r8,  0.1874E+00_r8,  0.2492E+00_r8,  0.3257E+00_r8,   &
         0.4165E+00_r8,  0.5191E+00_r8,  0.6275E+00_r8,  0.7324E+00_r8,  0.8235E+00_r8,   &
         0.8939E+00_r8, &
                                ! DATA ((h72(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=15_r8,15)_r8, &                               &
         -0.5045E-08_r8, -0.1113E-07_r8, -0.2540E-07_r8, -0.6008E-07_r8, -0.1449E-06_r8,  &
         -0.3457E-06_r8, -0.7826E-06_r8, -0.1609E-05_r8, -0.2920E-05_r8, -0.4665E-05_r8,  &
         -0.6691E-05_r8, -0.8868E-05_r8, -0.1127E-04_r8, -0.1422E-04_r8, -0.1820E-04_r8,  &
         -0.2389E-04_r8, -0.3213E-04_r8, -0.4380E-04_r8, -0.5975E-04_r8, -0.8092E-04_r8,  &
         -0.1083E-03_r8, -0.1433E-03_r8, -0.1873E-03_r8, -0.2402E-03_r8, -0.2997E-03_r8,  &
         -0.3607E-03_r8, -0.4178E-03_r8, -0.4662E-03_r8, -0.4994E-03_r8, -0.5028E-03_r8,  &
         -0.4563E-03_r8, &
                                ! DATA ((h73(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=15_r8,15)_r8, &                               &
         0.1804E-10_r8,  0.3983E-10_r8,  0.9045E-10_r8,  0.2080E-09_r8,  0.4786E-09_r8,   &
         0.1046E-08_r8,  0.2052E-08_r8,  0.3413E-08_r8,  0.4704E-08_r8,  0.5565E-08_r8,   &
         0.5887E-08_r8,  0.5981E-08_r8,  0.6202E-08_r8,  0.6998E-08_r8,  0.8493E-08_r8,   &
         0.1002E-07_r8,  0.1184E-07_r8,  0.1780E-07_r8,  0.3483E-07_r8,  0.7122E-07_r8,   &
         0.1341E-06_r8,  0.2259E-06_r8,  0.3446E-06_r8,  0.4866E-06_r8,  0.6486E-06_r8,   &
         0.8343E-06_r8,  0.1063E-05_r8,  0.1356E-05_r8,  0.1690E-05_r8,  0.1951E-05_r8,   &
         0.2005E-05_r8, &
                                ! DATA ((h71(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=16_r8,16)_r8, &                               &
         0.1351E-04_r8,  0.2672E-04_r8,  0.5247E-04_r8,  0.1019E-03_r8,  0.1944E-03_r8,   &
         0.3609E-03_r8,  0.6451E-03_r8,  0.1098E-02_r8,  0.1772E-02_r8,  0.2718E-02_r8,   &
         0.4009E-02_r8,  0.5767E-02_r8,  0.8190E-02_r8,  0.1158E-01_r8,  0.1636E-01_r8,   &
         0.2311E-01_r8,  0.3260E-01_r8,  0.4586E-01_r8,  0.6425E-01_r8,  0.8940E-01_r8,   &
         0.1232E+00_r8,  0.1678E+00_r8,  0.2251E+00_r8,  0.2967E+00_r8,  0.3831E+00_r8,   &
         0.4825E+00_r8,  0.5901E+00_r8,  0.6976E+00_r8,  0.7946E+00_r8,  0.8725E+00_r8,   &
         0.9285E+00_r8, &
                                ! DATA ((h72(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=16_r8,16)_r8, &                               &
         -0.4850E-08_r8, -0.1045E-07_r8, -0.2334E-07_r8, -0.5367E-07_r8, -0.1265E-06_r8,  &
         -0.2980E-06_r8, -0.6750E-06_r8, -0.1406E-05_r8, -0.2601E-05_r8, -0.4239E-05_r8,  &
         -0.6201E-05_r8, -0.8389E-05_r8, -0.1091E-04_r8, -0.1413E-04_r8, -0.1859E-04_r8,  &
         -0.2500E-04_r8, -0.3432E-04_r8, -0.4761E-04_r8, -0.6595E-04_r8, -0.9030E-04_r8,  &
         -0.1219E-03_r8, -0.1624E-03_r8, -0.2126E-03_r8, -0.2708E-03_r8, -0.3327E-03_r8,  &
         -0.3926E-03_r8, -0.4458E-03_r8, -0.4871E-03_r8, -0.5045E-03_r8, -0.4777E-03_r8,  &
         -0.3954E-03_r8, &
                                ! DATA ((h73(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=16_r8,16)_r8, &                               &
         0.1717E-10_r8,  0.3723E-10_r8,  0.8093E-10_r8,  0.1817E-09_r8,  0.4100E-09_r8,   &
         0.8932E-09_r8,  0.1791E-08_r8,  0.3126E-08_r8,  0.4634E-08_r8,  0.6095E-08_r8,   &
         0.7497E-08_r8,  0.9170E-08_r8,  0.1136E-07_r8,  0.1453E-07_r8,  0.1892E-07_r8,   &
         0.2369E-07_r8,  0.2909E-07_r8,  0.3922E-07_r8,  0.6232E-07_r8,  0.1083E-06_r8,   &
         0.1847E-06_r8,  0.2943E-06_r8,  0.4336E-06_r8,  0.5970E-06_r8,  0.7815E-06_r8,   &
         0.9959E-06_r8,  0.1263E-05_r8,  0.1583E-05_r8,  0.1880E-05_r8,  0.2009E-05_r8,   &
         0.1914E-05_r8, &
                                ! DATA ((h71(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=17_r8,17)_r8, &                               &
         0.1353E-04_r8,  0.2681E-04_r8,  0.5283E-04_r8,  0.1031E-03_r8,  0.1981E-03_r8,   &
         0.3721E-03_r8,  0.6761E-03_r8,  0.1176E-02_r8,  0.1943E-02_r8,  0.3055E-02_r8,   &
         0.4612E-02_r8,  0.6766E-02_r8,  0.9760E-02_r8,  0.1396E-01_r8,  0.1986E-01_r8,   &
         0.2818E-01_r8,  0.3984E-01_r8,  0.5609E-01_r8,  0.7850E-01_r8,  0.1090E+00_r8,   &
         0.1495E+00_r8,  0.2023E+00_r8,  0.2691E+00_r8,  0.3508E+00_r8,  0.4465E+00_r8,   &
         0.5524E+00_r8,  0.6613E+00_r8,  0.7633E+00_r8,  0.8485E+00_r8,  0.9119E+00_r8,   &
         0.9543E+00_r8, &
                                ! DATA ((h72(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=17_r8,17)_r8, &                               &
         -0.4673E-08_r8, -0.9862E-08_r8, -0.2135E-07_r8, -0.4753E-07_r8, -0.1087E-06_r8,  &
         -0.2512E-06_r8, -0.5671E-06_r8, -0.1199E-05_r8, -0.2281E-05_r8, -0.3842E-05_r8,  &
         -0.5804E-05_r8, -0.8110E-05_r8, -0.1088E-04_r8, -0.1452E-04_r8, -0.1961E-04_r8,  &
         -0.2696E-04_r8, -0.3768E-04_r8, -0.5311E-04_r8, -0.7444E-04_r8, -0.1028E-03_r8,  &
         -0.1397E-03_r8, -0.1865E-03_r8, -0.2427E-03_r8, -0.3047E-03_r8, -0.3667E-03_r8,  &
         -0.4237E-03_r8, -0.4712E-03_r8, -0.5003E-03_r8, -0.4921E-03_r8, -0.4286E-03_r8,  &
         -0.3188E-03_r8, &
                                ! DATA ((h73(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=17_r8,17)_r8, &                               &
         0.1653E-10_r8,  0.3436E-10_r8,  0.7431E-10_r8,  0.1605E-09_r8,  0.3548E-09_r8,   &
         0.7723E-09_r8,  0.1595E-08_r8,  0.2966E-08_r8,  0.4849E-08_r8,  0.7169E-08_r8,   &
         0.1003E-07_r8,  0.1366E-07_r8,  0.1825E-07_r8,  0.2419E-07_r8,  0.3186E-07_r8,   &
         0.4068E-07_r8,  0.5064E-07_r8,  0.6618E-07_r8,  0.9684E-07_r8,  0.1536E-06_r8,   &
         0.2450E-06_r8,  0.3730E-06_r8,  0.5328E-06_r8,  0.7184E-06_r8,  0.9291E-06_r8,   &
         0.1180E-05_r8,  0.1484E-05_r8,  0.1798E-05_r8,  0.1992E-05_r8,  0.1968E-05_r8,   &
         0.1736E-05_r8, &
                                ! DATA ((h71(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=18_r8,18)_r8, &                               &
         0.1356E-04_r8,  0.2691E-04_r8,  0.5318E-04_r8,  0.1043E-03_r8,  0.2020E-03_r8,   &
         0.3841E-03_r8,  0.7100E-03_r8,  0.1263E-02_r8,  0.2144E-02_r8,  0.3464E-02_r8,   &
         0.5359E-02_r8,  0.8019E-02_r8,  0.1174E-01_r8,  0.1694E-01_r8,  0.2426E-01_r8,   &
         0.3451E-01_r8,  0.4883E-01_r8,  0.6870E-01_r8,  0.9593E-01_r8,  0.1326E+00_r8,   &
         0.1809E+00_r8,  0.2428E+00_r8,  0.3196E+00_r8,  0.4112E+00_r8,  0.5147E+00_r8,   &
         0.6239E+00_r8,  0.7297E+00_r8,  0.8217E+00_r8,  0.8928E+00_r8,  0.9421E+00_r8,   &
         0.9726E+00_r8, &
                                ! DATA ((h72(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=18_r8,18)_r8, &                               &
         -0.4532E-08_r8, -0.9395E-08_r8, -0.1978E-07_r8, -0.4272E-07_r8, -0.9442E-07_r8,  &
         -0.2124E-06_r8, -0.4747E-06_r8, -0.1017E-05_r8, -0.2003E-05_r8, -0.3524E-05_r8,  &
         -0.5567E-05_r8, -0.8108E-05_r8, -0.1127E-04_r8, -0.1547E-04_r8, -0.2138E-04_r8,  &
         -0.2996E-04_r8, -0.4251E-04_r8, -0.6059E-04_r8, -0.8563E-04_r8, -0.1190E-03_r8,  &
         -0.1623E-03_r8, -0.2156E-03_r8, -0.2767E-03_r8, -0.3403E-03_r8, -0.4006E-03_r8,  &
         -0.4530E-03_r8, -0.4912E-03_r8, -0.4995E-03_r8, -0.4563E-03_r8, -0.3592E-03_r8,  &
         -0.2383E-03_r8, &
                                ! DATA ((h73(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=18_r8,18)_r8, &                               &
         0.1593E-10_r8,  0.3276E-10_r8,  0.6896E-10_r8,  0.1476E-09_r8,  0.3190E-09_r8,   &
         0.6944E-09_r8,  0.1474E-08_r8,  0.2935E-08_r8,  0.5300E-08_r8,  0.8697E-08_r8,   &
         0.1336E-07_r8,  0.1946E-07_r8,  0.2707E-07_r8,  0.3637E-07_r8,  0.4800E-07_r8,   &
         0.6187E-07_r8,  0.7806E-07_r8,  0.1008E-06_r8,  0.1404E-06_r8,  0.2089E-06_r8,   &
         0.3153E-06_r8,  0.4613E-06_r8,  0.6416E-06_r8,  0.8506E-06_r8,  0.1095E-05_r8,   &
         0.1387E-05_r8,  0.1708E-05_r8,  0.1956E-05_r8,  0.2003E-05_r8,  0.1836E-05_r8,   &
         0.1483E-05_r8, &
                                ! DATA ((h71(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=19_r8,19)_r8, &                               &
         0.1358E-04_r8,  0.2700E-04_r8,  0.5350E-04_r8,  0.1054E-03_r8,  0.2057E-03_r8,   &
         0.3955E-03_r8,  0.7434E-03_r8,  0.1353E-02_r8,  0.2361E-02_r8,  0.3928E-02_r8,   &
         0.6237E-02_r8,  0.9531E-02_r8,  0.1416E-01_r8,  0.2063E-01_r8,  0.2968E-01_r8,   &
         0.4233E-01_r8,  0.5990E-01_r8,  0.8412E-01_r8,  0.1170E+00_r8,  0.1608E+00_r8,   &
         0.2177E+00_r8,  0.2894E+00_r8,  0.3764E+00_r8,  0.4768E+00_r8,  0.5855E+00_r8,   &
         0.6941E+00_r8,  0.7921E+00_r8,  0.8710E+00_r8,  0.9276E+00_r8,  0.9641E+00_r8,   &
         0.9846E+00_r8, &
                                ! DATA ((h72(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=19_r8,19)_r8, &                               &
         -0.4448E-08_r8, -0.9085E-08_r8, -0.1877E-07_r8, -0.3946E-07_r8, -0.8472E-07_r8,  &
         -0.1852E-06_r8, -0.4074E-06_r8, -0.8791E-06_r8, -0.1789E-05_r8, -0.3314E-05_r8,  &
         -0.5521E-05_r8, -0.8425E-05_r8, -0.1215E-04_r8, -0.1711E-04_r8, -0.2407E-04_r8,  &
         -0.3421E-04_r8, -0.4905E-04_r8, -0.7032E-04_r8, -0.9985E-04_r8, -0.1394E-03_r8,  &
         -0.1897E-03_r8, -0.2491E-03_r8, -0.3132E-03_r8, -0.3763E-03_r8, -0.4332E-03_r8,  &
         -0.4786E-03_r8, -0.5005E-03_r8, -0.4775E-03_r8, -0.3970E-03_r8, -0.2794E-03_r8,  &
         -0.1652E-03_r8, &
                                ! DATA ((h73(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=19_r8,19)_r8, &                               &
         0.1566E-10_r8,  0.3219E-10_r8,  0.6635E-10_r8,  0.1400E-09_r8,  0.2999E-09_r8,   &
         0.6513E-09_r8,  0.1406E-08_r8,  0.2953E-08_r8,  0.5789E-08_r8,  0.1037E-07_r8,   &
         0.1709E-07_r8,  0.2623E-07_r8,  0.3777E-07_r8,  0.5159E-07_r8,  0.6823E-07_r8,   &
         0.8864E-07_r8,  0.1134E-06_r8,  0.1461E-06_r8,  0.1960E-06_r8,  0.2761E-06_r8,   &
         0.3962E-06_r8,  0.5583E-06_r8,  0.7580E-06_r8,  0.9957E-06_r8,  0.1282E-05_r8,   &
         0.1607E-05_r8,  0.1898E-05_r8,  0.2020E-05_r8,  0.1919E-05_r8,  0.1623E-05_r8,   &
         0.1171E-05_r8, &
                                ! DATA ((h71(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=20_r8,20)_r8, &                               &
         0.1360E-04_r8,  0.2707E-04_r8,  0.5375E-04_r8,  0.1063E-03_r8,  0.2087E-03_r8,   &
         0.4052E-03_r8,  0.7731E-03_r8,  0.1437E-02_r8,  0.2575E-02_r8,  0.4414E-02_r8,   &
         0.7212E-02_r8,  0.1128E-01_r8,  0.1704E-01_r8,  0.2509E-01_r8,  0.3631E-01_r8,   &
         0.5188E-01_r8,  0.7338E-01_r8,  0.1028E+00_r8,  0.1422E+00_r8,  0.1940E+00_r8,   &
         0.2603E+00_r8,  0.3422E+00_r8,  0.4387E+00_r8,  0.5459E+00_r8,  0.6563E+00_r8,   &
         0.7597E+00_r8,  0.8461E+00_r8,  0.9105E+00_r8,  0.9536E+00_r8,  0.9791E+00_r8,   &
         0.9920E+00_r8, &
                                ! DATA ((h72(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=20_r8,20)_r8, &                               &
         -0.4403E-08_r8, -0.8896E-08_r8, -0.1818E-07_r8, -0.3751E-07_r8, -0.7880E-07_r8,  &
         -0.1683E-06_r8, -0.3640E-06_r8, -0.7852E-06_r8, -0.1640E-05_r8, -0.3191E-05_r8,  &
         -0.5634E-05_r8, -0.9046E-05_r8, -0.1355E-04_r8, -0.1953E-04_r8, -0.2786E-04_r8,  &
         -0.3995E-04_r8, -0.5752E-04_r8, -0.8256E-04_r8, -0.1174E-03_r8, -0.1638E-03_r8,  &
         -0.2211E-03_r8, -0.2854E-03_r8, -0.3507E-03_r8, -0.4116E-03_r8, -0.4633E-03_r8,  &
         -0.4966E-03_r8, -0.4921E-03_r8, -0.4309E-03_r8, -0.3215E-03_r8, -0.2016E-03_r8,  &
         -0.1061E-03_r8, &
                                ! DATA ((h73(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=20_r8,20)_r8, &                               &
         0.1551E-10_r8,  0.3147E-10_r8,  0.6419E-10_r8,  0.1356E-09_r8,  0.2860E-09_r8,   &
         0.6178E-09_r8,  0.1353E-08_r8,  0.2934E-08_r8,  0.6095E-08_r8,  0.1174E-07_r8,   &
         0.2067E-07_r8,  0.3346E-07_r8,  0.5014E-07_r8,  0.7024E-07_r8,  0.9377E-07_r8,   &
         0.1226E-06_r8,  0.1592E-06_r8,  0.2056E-06_r8,  0.2678E-06_r8,  0.3584E-06_r8,   &
         0.4892E-06_r8,  0.6651E-06_r8,  0.8859E-06_r8,  0.1160E-05_r8,  0.1488E-05_r8,   &
         0.1814E-05_r8,  0.2010E-05_r8,  0.1984E-05_r8,  0.1748E-05_r8,  0.1338E-05_r8,   &
         0.8445E-06_r8, &
                                ! DATA ((h71(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=21_r8,21)_r8, &                               &
         0.1361E-04_r8,  0.2711E-04_r8,  0.5393E-04_r8,  0.1070E-03_r8,  0.2110E-03_r8,   &
         0.4129E-03_r8,  0.7973E-03_r8,  0.1508E-02_r8,  0.2769E-02_r8,  0.4886E-02_r8,   &
         0.8229E-02_r8,  0.1322E-01_r8,  0.2037E-01_r8,  0.3038E-01_r8,  0.4427E-01_r8,   &
         0.6342E-01_r8,  0.8963E-01_r8,  0.1250E+00_r8,  0.1720E+00_r8,  0.2327E+00_r8,   &
         0.3089E+00_r8,  0.4007E+00_r8,  0.5052E+00_r8,  0.6163E+00_r8,  0.7241E+00_r8,   &
         0.8179E+00_r8,  0.8905E+00_r8,  0.9408E+00_r8,  0.9720E+00_r8,  0.9887E+00_r8,   &
         0.9961E+00_r8, &
                                ! DATA ((h72(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=21_r8,21)_r8, &                               &
         -0.4379E-08_r8, -0.8801E-08_r8, -0.1782E-07_r8, -0.3642E-07_r8, -0.7536E-07_r8,  &
         -0.1581E-06_r8, -0.3366E-06_r8, -0.7227E-06_r8, -0.1532E-05_r8, -0.3106E-05_r8,  &
         -0.5810E-05_r8, -0.9862E-05_r8, -0.1540E-04_r8, -0.2279E-04_r8, -0.3292E-04_r8,  &
         -0.4738E-04_r8, -0.6817E-04_r8, -0.9765E-04_r8, -0.1384E-03_r8, -0.1918E-03_r8,  &
         -0.2551E-03_r8, -0.3226E-03_r8, -0.3876E-03_r8, -0.4452E-03_r8, -0.4883E-03_r8,  &
         -0.5005E-03_r8, -0.4598E-03_r8, -0.3633E-03_r8, -0.2416E-03_r8, -0.1349E-03_r8,  &
         -0.6278E-04_r8, &
                                ! DATA ((h73(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=21_r8,21)_r8, &                               &
         0.1542E-10_r8,  0.3111E-10_r8,  0.6345E-10_r8,  0.1310E-09_r8,  0.2742E-09_r8,   &
         0.5902E-09_r8,  0.1289E-08_r8,  0.2826E-08_r8,  0.6103E-08_r8,  0.1250E-07_r8,   &
         0.2355E-07_r8,  0.4041E-07_r8,  0.6347E-07_r8,  0.9217E-07_r8,  0.1256E-06_r8,   &
         0.1658E-06_r8,  0.2175E-06_r8,  0.2824E-06_r8,  0.3607E-06_r8,  0.4614E-06_r8,   &
         0.6004E-06_r8,  0.7880E-06_r8,  0.1034E-05_r8,  0.1349E-05_r8,  0.1698E-05_r8,   &
         0.1965E-05_r8,  0.2021E-05_r8,  0.1857E-05_r8,  0.1500E-05_r8,  0.1015E-05_r8,   &
         0.5467E-06_r8, &
                                ! DATA ((h71(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=22_r8,22)_r8, &                               &
         0.1362E-04_r8,  0.2714E-04_r8,  0.5405E-04_r8,  0.1074E-03_r8,  0.2127E-03_r8,   &
         0.4187E-03_r8,  0.8159E-03_r8,  0.1565E-02_r8,  0.2933E-02_r8,  0.5314E-02_r8,   &
         0.9222E-02_r8,  0.1525E-01_r8,  0.2406E-01_r8,  0.3649E-01_r8,  0.5367E-01_r8,   &
         0.7717E-01_r8,  0.1090E+00_r8,  0.1514E+00_r8,  0.2067E+00_r8,  0.2770E+00_r8,   &
         0.3633E+00_r8,  0.4640E+00_r8,  0.5742E+00_r8,  0.6853E+00_r8,  0.7859E+00_r8,   &
         0.8670E+00_r8,  0.9253E+00_r8,  0.9629E+00_r8,  0.9842E+00_r8,  0.9942E+00_r8,   &
         0.9983E+00_r8, &
                                ! DATA ((h72(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=22_r8,22)_r8, &                               &
         -0.4366E-08_r8, -0.8749E-08_r8, -0.1761E-07_r8, -0.3578E-07_r8, -0.7322E-07_r8,  &
         -0.1517E-06_r8, -0.3189E-06_r8, -0.6785E-06_r8, -0.1446E-05_r8, -0.3014E-05_r8,  &
         -0.5933E-05_r8, -0.1069E-04_r8, -0.1755E-04_r8, -0.2683E-04_r8, -0.3936E-04_r8,  &
         -0.5675E-04_r8, -0.8137E-04_r8, -0.1160E-03_r8, -0.1630E-03_r8, -0.2223E-03_r8,  &
         -0.2899E-03_r8, -0.3589E-03_r8, -0.4230E-03_r8, -0.4755E-03_r8, -0.5031E-03_r8,  &
         -0.4834E-03_r8, -0.4036E-03_r8, -0.2849E-03_r8, -0.1687E-03_r8, -0.8356E-04_r8,  &
         -0.3388E-04_r8, &
                                ! DATA ((h73(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=22_r8,22)_r8, &                               &
         0.1536E-10_r8,  0.3086E-10_r8,  0.6248E-10_r8,  0.1288E-09_r8,  0.2664E-09_r8,   &
         0.5637E-09_r8,  0.1222E-08_r8,  0.2680E-08_r8,  0.5899E-08_r8,  0.1262E-07_r8,   &
         0.2527E-07_r8,  0.4621E-07_r8,  0.7678E-07_r8,  0.1165E-06_r8,  0.1640E-06_r8,   &
         0.2199E-06_r8,  0.2904E-06_r8,  0.3783E-06_r8,  0.4787E-06_r8,  0.5925E-06_r8,   &
         0.7377E-06_r8,  0.9389E-06_r8,  0.1216E-05_r8,  0.1560E-05_r8,  0.1879E-05_r8,   &
         0.2025E-05_r8,  0.1940E-05_r8,  0.1650E-05_r8,  0.1194E-05_r8,  0.6981E-06_r8,   &
         0.3103E-06_r8, &
                                ! DATA ((h71(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=23_r8,23)_r8, &                               &
         0.1362E-04_r8,  0.2716E-04_r8,  0.5413E-04_r8,  0.1077E-03_r8,  0.2139E-03_r8,   &
         0.4228E-03_r8,  0.8295E-03_r8,  0.1609E-02_r8,  0.3063E-02_r8,  0.5676E-02_r8,   &
         0.1013E-01_r8,  0.1726E-01_r8,  0.2799E-01_r8,  0.4332E-01_r8,  0.6452E-01_r8,   &
         0.9328E-01_r8,  0.1317E+00_r8,  0.1820E+00_r8,  0.2467E+00_r8,  0.3270E+00_r8,   &
         0.4228E+00_r8,  0.5308E+00_r8,  0.6434E+00_r8,  0.7499E+00_r8,  0.8395E+00_r8,   &
         0.9065E+00_r8,  0.9515E+00_r8,  0.9782E+00_r8,  0.9916E+00_r8,  0.9973E+00_r8,   &
         0.9993E+00_r8, &
                                ! DATA ((h72(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=23_r8,23)_r8, &                               &
         -0.4359E-08_r8, -0.8720E-08_r8, -0.1749E-07_r8, -0.3527E-07_r8, -0.7175E-07_r8,  &
         -0.1473E-06_r8, -0.3062E-06_r8, -0.6451E-06_r8, -0.1372E-05_r8, -0.2902E-05_r8,  &
         -0.5936E-05_r8, -0.1133E-04_r8, -0.1971E-04_r8, -0.3143E-04_r8, -0.4715E-04_r8,  &
         -0.6833E-04_r8, -0.9759E-04_r8, -0.1379E-03_r8, -0.1907E-03_r8, -0.2542E-03_r8,  &
         -0.3239E-03_r8, -0.3935E-03_r8, -0.4559E-03_r8, -0.4991E-03_r8, -0.5009E-03_r8,  &
         -0.4414E-03_r8, -0.3306E-03_r8, -0.2077E-03_r8, -0.1093E-03_r8, -0.4754E-04_r8,  &
         -0.1642E-04_r8, &
                                ! DATA ((h73(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=23_r8,23)_r8, &                               &
         0.1531E-10_r8,  0.3070E-10_r8,  0.6184E-10_r8,  0.1257E-09_r8,  0.2578E-09_r8,   &
         0.5451E-09_r8,  0.1159E-08_r8,  0.2526E-08_r8,  0.5585E-08_r8,  0.1225E-07_r8,   &
         0.2576E-07_r8,  0.5017E-07_r8,  0.8855E-07_r8,  0.1417E-06_r8,  0.2078E-06_r8,   &
         0.2858E-06_r8,  0.3802E-06_r8,  0.4946E-06_r8,  0.6226E-06_r8,  0.7572E-06_r8,   &
         0.9137E-06_r8,  0.1133E-05_r8,  0.1438E-05_r8,  0.1772E-05_r8,  0.1994E-05_r8,   &
         0.1994E-05_r8,  0.1779E-05_r8,  0.1375E-05_r8,  0.8711E-06_r8,  0.4273E-06_r8,   &
         0.1539E-06_r8, &
                                ! DATA ((h71(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=24_r8,24)_r8, &                               &
         0.1362E-04_r8,  0.2717E-04_r8,  0.5417E-04_r8,  0.1079E-03_r8,  0.2146E-03_r8,   &
         0.4256E-03_r8,  0.8393E-03_r8,  0.1641E-02_r8,  0.3163E-02_r8,  0.5968E-02_r8,   &
         0.1092E-01_r8,  0.1915E-01_r8,  0.3195E-01_r8,  0.5066E-01_r8,  0.7669E-01_r8,   &
         0.1117E+00_r8,  0.1578E+00_r8,  0.2173E+00_r8,  0.2919E+00_r8,  0.3824E+00_r8,   &
         0.4867E+00_r8,  0.5992E+00_r8,  0.7100E+00_r8,  0.8076E+00_r8,  0.8838E+00_r8,   &
         0.9370E+00_r8,  0.9702E+00_r8,  0.9880E+00_r8,  0.9959E+00_r8,  0.9988E+00_r8,   &
         0.9998E+00_r8, &
                                ! DATA ((h72(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=24_r8,24)_r8, &                               &
         -0.4354E-08_r8, -0.8703E-08_r8, -0.1742E-07_r8, -0.3499E-07_r8, -0.7074E-07_r8,  &
         -0.1441E-06_r8, -0.2971E-06_r8, -0.6195E-06_r8, -0.1309E-05_r8, -0.2780E-05_r8,  &
         -0.5823E-05_r8, -0.1165E-04_r8, -0.2152E-04_r8, -0.3616E-04_r8, -0.5604E-04_r8,  &
         -0.8230E-04_r8, -0.1173E-03_r8, -0.1635E-03_r8, -0.2211E-03_r8, -0.2868E-03_r8,  &
         -0.3567E-03_r8, -0.4260E-03_r8, -0.4844E-03_r8, -0.5097E-03_r8, -0.4750E-03_r8,  &
         -0.3779E-03_r8, -0.2522E-03_r8, -0.1409E-03_r8, -0.6540E-04_r8, -0.2449E-04_r8,  &
         -0.6948E-05_r8, &
                                ! DATA ((h73(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=24_r8,24)_r8, &                               &
         0.1529E-10_r8,  0.3060E-10_r8,  0.6142E-10_r8,  0.1241E-09_r8,  0.2535E-09_r8,   &
         0.5259E-09_r8,  0.1107E-08_r8,  0.2383E-08_r8,  0.5243E-08_r8,  0.1161E-07_r8,   &
         0.2523E-07_r8,  0.5188E-07_r8,  0.9757E-07_r8,  0.1657E-06_r8,  0.2553E-06_r8,   &
         0.3629E-06_r8,  0.4878E-06_r8,  0.6323E-06_r8,  0.7923E-06_r8,  0.9575E-06_r8,   &
         0.1139E-05_r8,  0.1381E-05_r8,  0.1687E-05_r8,  0.1952E-05_r8,  0.2029E-05_r8,   &
         0.1890E-05_r8,  0.1552E-05_r8,  0.1062E-05_r8,  0.5728E-06_r8,  0.2280E-06_r8,   &
         0.6762E-07_r8, &
                                ! DATA ((h71(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=25_r8,25)_r8, &                               &
         0.1363E-04_r8,  0.2718E-04_r8,  0.5420E-04_r8,  0.1080E-03_r8,  0.2151E-03_r8,   &
         0.4275E-03_r8,  0.8461E-03_r8,  0.1664E-02_r8,  0.3237E-02_r8,  0.6195E-02_r8,   &
         0.1156E-01_r8,  0.2080E-01_r8,  0.3573E-01_r8,  0.5818E-01_r8,  0.8981E-01_r8,   &
         0.1323E+00_r8,  0.1874E+00_r8,  0.2570E+00_r8,  0.3423E+00_r8,  0.4425E+00_r8,   &
         0.5534E+00_r8,  0.6669E+00_r8,  0.7714E+00_r8,  0.8568E+00_r8,  0.9191E+00_r8,   &
         0.9597E+00_r8,  0.9828E+00_r8,  0.9938E+00_r8,  0.9981E+00_r8,  0.9996E+00_r8,   &
         0.1000E+01_r8, &
                                ! DATA ((h72(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=25_r8,25)_r8, &                               &
         -0.4352E-08_r8, -0.8693E-08_r8, -0.1738E-07_r8, -0.3483E-07_r8, -0.7006E-07_r8,  &
         -0.1423E-06_r8, -0.2905E-06_r8, -0.6008E-06_r8, -0.1258E-05_r8, -0.2663E-05_r8,  &
         -0.5638E-05_r8, -0.1165E-04_r8, -0.2270E-04_r8, -0.4044E-04_r8, -0.6554E-04_r8,  &
         -0.9855E-04_r8, -0.1407E-03_r8, -0.1928E-03_r8, -0.2534E-03_r8, -0.3197E-03_r8,  &
         -0.3890E-03_r8, -0.4563E-03_r8, -0.5040E-03_r8, -0.4998E-03_r8, -0.4249E-03_r8,  &
         -0.3025E-03_r8, -0.1794E-03_r8, -0.8860E-04_r8, -0.3575E-04_r8, -0.1122E-04_r8,  &
         -0.2506E-05_r8, &
                                ! DATA ((h73(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=25_r8,25)_r8, &                               &
         0.1527E-10_r8,  0.3053E-10_r8,  0.6115E-10_r8,  0.1230E-09_r8,  0.2492E-09_r8,   &
         0.5149E-09_r8,  0.1068E-08_r8,  0.2268E-08_r8,  0.4932E-08_r8,  0.1089E-07_r8,   &
         0.2408E-07_r8,  0.5156E-07_r8,  0.1028E-06_r8,  0.1859E-06_r8,  0.3028E-06_r8,   &
         0.4476E-06_r8,  0.6124E-06_r8,  0.7932E-06_r8,  0.9879E-06_r8,  0.1194E-05_r8,   &
         0.1417E-05_r8,  0.1673E-05_r8,  0.1929E-05_r8,  0.2064E-05_r8,  0.1997E-05_r8,   &
         0.1725E-05_r8,  0.1267E-05_r8,  0.7464E-06_r8,  0.3312E-06_r8,  0.1066E-06_r8,   &
         0.2718E-07_r8, &
                                ! DATA ((h71(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=26_r8,26)_r8, &                               &
         0.1363E-04_r8,  0.2719E-04_r8,  0.5422E-04_r8,  0.1081E-03_r8,  0.2154E-03_r8,   &
         0.4287E-03_r8,  0.8506E-03_r8,  0.1680E-02_r8,  0.3291E-02_r8,  0.6364E-02_r8,   &
         0.1206E-01_r8,  0.2218E-01_r8,  0.3912E-01_r8,  0.6544E-01_r8,  0.1033E+00_r8,   &
         0.1544E+00_r8,  0.2198E+00_r8,  0.3008E+00_r8,  0.3970E+00_r8,  0.5060E+00_r8,   &
         0.6209E+00_r8,  0.7311E+00_r8,  0.8254E+00_r8,  0.8972E+00_r8,  0.9462E+00_r8,   &
         0.9757E+00_r8,  0.9907E+00_r8,  0.9970E+00_r8,  0.9992E+00_r8,  0.9999E+00_r8,   &
         0.1000E+01_r8, &
                                ! DATA ((h72(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=26_r8,26)_r8, &                               &
         -0.4351E-08_r8, -0.8688E-08_r8, -0.1736E-07_r8, -0.3473E-07_r8, -0.6966E-07_r8,  &
         -0.1405E-06_r8, -0.2857E-06_r8, -0.5867E-06_r8, -0.1218E-05_r8, -0.2563E-05_r8,  &
         -0.5435E-05_r8, -0.1144E-04_r8, -0.2321E-04_r8, -0.4379E-04_r8, -0.7487E-04_r8,  &
         -0.1163E-03_r8, -0.1670E-03_r8, -0.2250E-03_r8, -0.2876E-03_r8, -0.3535E-03_r8,  &
         -0.4215E-03_r8, -0.4826E-03_r8, -0.5082E-03_r8, -0.4649E-03_r8, -0.3564E-03_r8,  &
         -0.2264E-03_r8, -0.1188E-03_r8, -0.5128E-04_r8, -0.1758E-04_r8, -0.4431E-05_r8,  &
         -0.7275E-06_r8, &
                                ! DATA ((h73(ip_r8,iw_r8,1)_r8,iw=1_r8,31)_r8,ip=26_r8,26)/                               &
         0.1525E-10_r8,  0.3048E-10_r8,  0.6097E-10_r8,  0.1223E-09_r8,  0.2466E-09_r8,   &
         0.5021E-09_r8,  0.1032E-08_r8,  0.2195E-08_r8,  0.4688E-08_r8,  0.1027E-07_r8,   &
         0.2279E-07_r8,  0.4999E-07_r8,  0.1046E-06_r8,  0.2009E-06_r8,  0.3460E-06_r8,   &
         0.5335E-06_r8,  0.7478E-06_r8,  0.9767E-06_r8,  0.1216E-05_r8,  0.1469E-05_r8,   &
         0.1735E-05_r8,  0.1977E-05_r8,  0.2121E-05_r8,  0.2103E-05_r8,  0.1902E-05_r8,   &
         0.1495E-05_r8,  0.9541E-06_r8,  0.4681E-06_r8,  0.1672E-06_r8,  0.4496E-07_r8,   &
         0.9859E-08_r8/),SHAPE=(/nx*3*nh/))

    it=0
    DO k=1,nx
       DO i=1,3
          DO j=1,nh
             it=it+1              
             IF(i==1)THEN
                ! WRITE(*,'(a5,2e17.9) ' )'h81   ',data5(it),h81(k,j)
                h71(k,j,1)=data5(it)
             END IF
             IF(i==2)THEN
                ! WRITE(*,'(a5,2e17.9) ' )'h82   ',data5(it),h82(k,j)
                h72(k,j,1)=data5(it)
             END IF
             IF(i==3)THEN
                !WRITE(*,'(a5,2e17.9) ' )'h83   ',data5(it),h83(k,j)
                h73(k,j,1)=data5(it)
             END IF
          END DO
       END DO
    END DO

    !
    !-----------------------------------------------------------------------
    !
    !  The following DATA statements originally came from file
    !  "co2.tran3", which define pre-computed tables used for co2 (band
    !  3) transmittance calculations.
    !
    !-----------------------------------------------------------------------
    !
    ! co2.tran3

    !  REAL :: c1(26,24,7),c2(26,24,7),c3(26,24,7)
    !
    !  COMMON /radtab002/ c1,c2,c3
    data1(1:nx*3*nc)=RESHAPE(SOURCE=(/&
                                ! DATA ((c1(ip,iw,1),iw=1,24),ip= 1, 1)/                                &
         0.1444E-03_r8,  0.2378E-03_r8,  0.3644E-03_r8,  0.5245E-03_r8,  0.7311E-03_r8,   &
         0.1015E-02_r8,  0.1416E-02_r8,  0.1975E-02_r8,  0.2715E-02_r8,  0.3651E-02_r8,   &
         0.4813E-02_r8,  0.6264E-02_r8,  0.8084E-02_r8,  0.1037E-01_r8,  0.1327E-01_r8,   &
         0.1693E-01_r8,  0.2155E-01_r8,  0.2734E-01_r8,  0.3453E-01_r8,  0.4337E-01_r8,   &
         0.5414E-01_r8,  0.6726E-01_r8,  0.8335E-01_r8,  0.1032E+00_r8, &
                                ! DATA ((c2(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip= 1_r8, 1)_r8, &                                &
         -0.2878E-06_r8, -0.6295E-06_r8, -0.1270E-05_r8, -0.2323E-05_r8, -0.3910E-05_r8,  &
         -0.6235E-05_r8, -0.9527E-05_r8, -0.1391E-04_r8, -0.1947E-04_r8, -0.2652E-04_r8,  &
         -0.3568E-04_r8, -0.4777E-04_r8, -0.6363E-04_r8, -0.8414E-04_r8, -0.1102E-03_r8,  &
         -0.1425E-03_r8, -0.1822E-03_r8, -0.2301E-03_r8, -0.2869E-03_r8, -0.3534E-03_r8,  &
         -0.4305E-03_r8, -0.5188E-03_r8, -0.6173E-03_r8, -0.7229E-03_r8, &
                                ! DATA ((c3(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip= 1_r8, 1)_r8, &                                &
         -0.2369E-09_r8, -0.5453E-09_r8, -0.1434E-08_r8, -0.3513E-08_r8, -0.7277E-08_r8,  &
         -0.1283E-07_r8, -0.1983E-07_r8, -0.2838E-07_r8, -0.4075E-07_r8, -0.6115E-07_r8,  &
         -0.9201E-07_r8, -0.1325E-06_r8, -0.1844E-06_r8, -0.2534E-06_r8, -0.3430E-06_r8,  &
         -0.4550E-06_r8, -0.5896E-06_r8, -0.7434E-06_r8, -0.9121E-06_r8, -0.1096E-05_r8,  &
         -0.1298E-05_r8, -0.1508E-05_r8, -0.1689E-05_r8, -0.1773E-05_r8, &
                                ! DATA ((c1(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip= 2_r8, 2)_r8, &                                &
         0.1444E-03_r8,  0.2379E-03_r8,  0.3646E-03_r8,  0.5250E-03_r8,  0.7320E-03_r8,   &
         0.1016E-02_r8,  0.1419E-02_r8,  0.1981E-02_r8,  0.2726E-02_r8,  0.3669E-02_r8,   &
         0.4847E-02_r8,  0.6325E-02_r8,  0.8193E-02_r8,  0.1056E-01_r8,  0.1358E-01_r8,   &
         0.1742E-01_r8,  0.2228E-01_r8,  0.2843E-01_r8,  0.3611E-01_r8,  0.4561E-01_r8,   &
         0.5730E-01_r8,  0.7169E-01_r8,  0.8952E-01_r8,  0.1117E+00_r8, &
                                ! DATA ((c2(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip= 2_r8, 2)_r8, &                                &
         -0.2875E-06_r8, -0.6289E-06_r8, -0.1269E-05_r8, -0.2321E-05_r8, -0.3907E-05_r8,  &
         -0.6229E-05_r8, -0.9516E-05_r8, -0.1389E-04_r8, -0.1944E-04_r8, -0.2646E-04_r8,  &
         -0.3559E-04_r8, -0.4762E-04_r8, -0.6340E-04_r8, -0.8384E-04_r8, -0.1098E-03_r8,  &
         -0.1423E-03_r8, -0.1824E-03_r8, -0.2309E-03_r8, -0.2887E-03_r8, -0.3566E-03_r8,  &
         -0.4355E-03_r8, -0.5260E-03_r8, -0.6274E-03_r8, -0.7367E-03_r8, &
                                ! DATA ((c3(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip= 2_r8, 2)_r8, &                                &
         -0.2382E-09_r8, -0.5476E-09_r8, -0.1437E-08_r8, -0.3520E-08_r8, -0.7291E-08_r8,  &
         -0.1286E-07_r8, -0.1989E-07_r8, -0.2849E-07_r8, -0.4093E-07_r8, -0.6144E-07_r8,  &
         -0.9241E-07_r8, -0.1330E-06_r8, -0.1852E-06_r8, -0.2542E-06_r8, -0.3435E-06_r8,  &
         -0.4549E-06_r8, -0.5891E-06_r8, -0.7427E-06_r8, -0.9114E-06_r8, -0.1096E-05_r8,  &
         -0.1298E-05_r8, -0.1509E-05_r8, -0.1689E-05_r8, -0.1772E-05_r8, &
                                ! DATA ((c1(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip= 3_r8, 3)_r8, &                                &
         0.1445E-03_r8,  0.2381E-03_r8,  0.3650E-03_r8,  0.5258E-03_r8,  0.7335E-03_r8,   &
         0.1019E-02_r8,  0.1424E-02_r8,  0.1990E-02_r8,  0.2742E-02_r8,  0.3698E-02_r8,   &
         0.4899E-02_r8,  0.6418E-02_r8,  0.8353E-02_r8,  0.1083E-01_r8,  0.1401E-01_r8,   &
         0.1807E-01_r8,  0.2326E-01_r8,  0.2985E-01_r8,  0.3814E-01_r8,  0.4849E-01_r8,   &
         0.6134E-01_r8,  0.7733E-01_r8,  0.9732E-01_r8,  0.1223E+00_r8, &
                                ! DATA ((c2(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip= 3_r8, 3)_r8, &                                &
         -0.2871E-06_r8, -0.6279E-06_r8, -0.1267E-05_r8, -0.2318E-05_r8, -0.3901E-05_r8,  &
         -0.6219E-05_r8, -0.9500E-05_r8, -0.1386E-04_r8, -0.1939E-04_r8, -0.2639E-04_r8,  &
         -0.3547E-04_r8, -0.4743E-04_r8, -0.6313E-04_r8, -0.8354E-04_r8, -0.1096E-03_r8,  &
         -0.1425E-03_r8, -0.1831E-03_r8, -0.2325E-03_r8, -0.2916E-03_r8, -0.3613E-03_r8,  &
         -0.4425E-03_r8, -0.5361E-03_r8, -0.6412E-03_r8, -0.7550E-03_r8, &
                                ! DATA ((c3(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip= 3_r8, 3)_r8, &                                &
         -0.2404E-09_r8, -0.5516E-09_r8, -0.1443E-08_r8, -0.3531E-08_r8, -0.7310E-08_r8,  &
         -0.1290E-07_r8, -0.1997E-07_r8, -0.2865E-07_r8, -0.4121E-07_r8, -0.6185E-07_r8,  &
         -0.9295E-07_r8, -0.1337E-06_r8, -0.1860E-06_r8, -0.2549E-06_r8, -0.3438E-06_r8,  &
         -0.4547E-06_r8, -0.5888E-06_r8, -0.7424E-06_r8, -0.9114E-06_r8, -0.1097E-05_r8,  &
         -0.1300E-05_r8, -0.1511E-05_r8, -0.1691E-05_r8, -0.1771E-05_r8, &
                                ! DATA ((c1(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip= 4_r8, 4)_r8, &                                &
         0.1446E-03_r8,  0.2383E-03_r8,  0.3656E-03_r8,  0.5270E-03_r8,  0.7358E-03_r8,   &
         0.1023E-02_r8,  0.1432E-02_r8,  0.2003E-02_r8,  0.2766E-02_r8,  0.3743E-02_r8,   &
         0.4979E-02_r8,  0.6557E-02_r8,  0.8587E-02_r8,  0.1121E-01_r8,  0.1459E-01_r8,   &
         0.1895E-01_r8,  0.2455E-01_r8,  0.3170E-01_r8,  0.4076E-01_r8,  0.5217E-01_r8,   &
         0.6649E-01_r8,  0.8448E-01_r8,  0.1071E+00_r8,  0.1353E+00_r8, &
                                ! DATA ((c2(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip= 4_r8, 4)_r8, &                                &
         -0.2865E-06_r8, -0.6264E-06_r8, -0.1264E-05_r8, -0.2313E-05_r8, -0.3892E-05_r8,  &
         -0.6204E-05_r8, -0.9474E-05_r8, -0.1382E-04_r8, -0.1932E-04_r8, -0.2628E-04_r8,  &
         -0.3530E-04_r8, -0.4720E-04_r8, -0.6286E-04_r8, -0.8333E-04_r8, -0.1097E-03_r8,  &
         -0.1430E-03_r8, -0.1845E-03_r8, -0.2352E-03_r8, -0.2960E-03_r8, -0.3681E-03_r8,  &
         -0.4523E-03_r8, -0.5496E-03_r8, -0.6595E-03_r8, -0.7790E-03_r8, &
                                ! DATA ((c3(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip= 4_r8, 4)_r8, &                                &
         -0.2433E-09_r8, -0.5580E-09_r8, -0.1454E-08_r8, -0.3550E-08_r8, -0.7345E-08_r8,  &
         -0.1296E-07_r8, -0.2010E-07_r8, -0.2889E-07_r8, -0.4159E-07_r8, -0.6239E-07_r8,  &
         -0.9360E-07_r8, -0.1345E-06_r8, -0.1867E-06_r8, -0.2554E-06_r8, -0.3439E-06_r8,  &
         -0.4548E-06_r8, -0.5890E-06_r8, -0.7429E-06_r8, -0.9127E-06_r8, -0.1099E-05_r8,  &
         -0.1304E-05_r8, -0.1516E-05_r8, -0.1693E-05_r8, -0.1770E-05_r8, &
                                ! DATA ((c1(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip= 5_r8, 5)_r8, &                                &
         0.1448E-03_r8,  0.2387E-03_r8,  0.3666E-03_r8,  0.5290E-03_r8,  0.7394E-03_r8,   &
         0.1030E-02_r8,  0.1444E-02_r8,  0.2025E-02_r8,  0.2805E-02_r8,  0.3811E-02_r8,   &
         0.5099E-02_r8,  0.6761E-02_r8,  0.8919E-02_r8,  0.1173E-01_r8,  0.1538E-01_r8,   &
         0.2011E-01_r8,  0.2622E-01_r8,  0.3408E-01_r8,  0.4412E-01_r8,  0.5686E-01_r8,   &
         0.7302E-01_r8,  0.9346E-01_r8,  0.1192E+00_r8,  0.1512E+00_r8, &
                                ! DATA ((c2(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip= 5_r8, 5)_r8, &                                &
         -0.2854E-06_r8, -0.6241E-06_r8, -0.1260E-05_r8, -0.2305E-05_r8, -0.3878E-05_r8,  &
         -0.6181E-05_r8, -0.9435E-05_r8, -0.1375E-04_r8, -0.1922E-04_r8, -0.2613E-04_r8,  &
         -0.3510E-04_r8, -0.4696E-04_r8, -0.6265E-04_r8, -0.8334E-04_r8, -0.1102E-03_r8,  &
         -0.1443E-03_r8, -0.1870E-03_r8, -0.2393E-03_r8, -0.3025E-03_r8, -0.3775E-03_r8,  &
         -0.4656E-03_r8, -0.5678E-03_r8, -0.6835E-03_r8, -0.8098E-03_r8, &
                                ! DATA ((c3(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip= 5_r8, 5)_r8, &                                &
         -0.2485E-09_r8, -0.5672E-09_r8, -0.1468E-08_r8, -0.3577E-08_r8, -0.7395E-08_r8,  &
         -0.1306E-07_r8, -0.2029E-07_r8, -0.2923E-07_r8, -0.4213E-07_r8, -0.6307E-07_r8,  &
         -0.9429E-07_r8, -0.1351E-06_r8, -0.1872E-06_r8, -0.2557E-06_r8, -0.3442E-06_r8,  &
         -0.4553E-06_r8, -0.5899E-06_r8, -0.7447E-06_r8, -0.9163E-06_r8, -0.1105E-05_r8,  &
         -0.1311E-05_r8, -0.1522E-05_r8, -0.1697E-05_r8, -0.1766E-05_r8, &
                                ! DATA ((c1(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip= 6_r8, 6)_r8, &                                &
         0.1450E-03_r8,  0.2394E-03_r8,  0.3680E-03_r8,  0.5320E-03_r8,  0.7452E-03_r8,   &
         0.1040E-02_r8,  0.1462E-02_r8,  0.2058E-02_r8,  0.2863E-02_r8,  0.3914E-02_r8,   &
         0.5275E-02_r8,  0.7051E-02_r8,  0.9378E-02_r8,  0.1243E-01_r8,  0.1641E-01_r8,   &
         0.2162E-01_r8,  0.2838E-01_r8,  0.3713E-01_r8,  0.4839E-01_r8,  0.6282E-01_r8,   &
         0.8125E-01_r8,  0.1046E+00_r8,  0.1340E+00_r8,  0.1702E+00_r8, &
                                ! DATA ((c2(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip= 6_r8, 6)_r8, &                                &
         -0.2837E-06_r8, -0.6203E-06_r8, -0.1252E-05_r8, -0.2292E-05_r8, -0.3857E-05_r8,  &
         -0.6146E-05_r8, -0.9377E-05_r8, -0.1366E-04_r8, -0.1909E-04_r8, -0.2596E-04_r8,  &
         -0.3489E-04_r8, -0.4676E-04_r8, -0.6263E-04_r8, -0.8372E-04_r8, -0.1113E-03_r8,  &
         -0.1465E-03_r8, -0.1908E-03_r8, -0.2454E-03_r8, -0.3115E-03_r8, -0.3903E-03_r8,  &
         -0.4833E-03_r8, -0.5915E-03_r8, -0.7144E-03_r8, -0.8483E-03_r8, &
                                ! DATA ((c3(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip= 6_r8, 6)_r8, &                                &
         -0.2557E-09_r8, -0.5835E-09_r8, -0.1495E-08_r8, -0.3619E-08_r8, -0.7472E-08_r8,  &
         -0.1320E-07_r8, -0.2055E-07_r8, -0.2970E-07_r8, -0.4281E-07_r8, -0.6379E-07_r8,  &
         -0.9492E-07_r8, -0.1356E-06_r8, -0.1876E-06_r8, -0.2561E-06_r8, -0.3451E-06_r8,  &
         -0.4568E-06_r8, -0.5924E-06_r8, -0.7489E-06_r8, -0.9229E-06_r8, -0.1114E-05_r8,  &
         -0.1321E-05_r8, -0.1531E-05_r8, -0.1700E-05_r8, -0.1756E-05_r8, &
                                ! DATA ((c1(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip= 7_r8, 7)_r8, &                                &
         0.1454E-03_r8,  0.2404E-03_r8,  0.3704E-03_r8,  0.5368E-03_r8,  0.7542E-03_r8,   &
         0.1056E-02_r8,  0.1491E-02_r8,  0.2108E-02_r8,  0.2951E-02_r8,  0.4066E-02_r8,   &
         0.5529E-02_r8,  0.7457E-02_r8,  0.1000E-01_r8,  0.1336E-01_r8,  0.1777E-01_r8,   &
         0.2357E-01_r8,  0.3115E-01_r8,  0.4102E-01_r8,  0.5382E-01_r8,  0.7035E-01_r8,   &
         0.9154E-01_r8,  0.1184E+00_r8,  0.1519E+00_r8,  0.1926E+00_r8, &
                                ! DATA ((c2(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip= 7_r8, 7)_r8, &                                &
         -0.2811E-06_r8, -0.6147E-06_r8, -0.1241E-05_r8, -0.2272E-05_r8, -0.3825E-05_r8,  &
         -0.6093E-05_r8, -0.9294E-05_r8, -0.1354E-04_r8, -0.1893E-04_r8, -0.2577E-04_r8,  &
         -0.3472E-04_r8, -0.4672E-04_r8, -0.6292E-04_r8, -0.8464E-04_r8, -0.1132E-03_r8,  &
         -0.1500E-03_r8, -0.1964E-03_r8, -0.2539E-03_r8, -0.3237E-03_r8, -0.4074E-03_r8,  &
         -0.5065E-03_r8, -0.6222E-03_r8, -0.7532E-03_r8, -0.8955E-03_r8, &
                                ! DATA ((c3(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip= 7_r8, 7)_r8, &                                &
         -0.2673E-09_r8, -0.6022E-09_r8, -0.1531E-08_r8, -0.3682E-08_r8, -0.7588E-08_r8,  &
         -0.1342E-07_r8, -0.2092E-07_r8, -0.3030E-07_r8, -0.4357E-07_r8, -0.6448E-07_r8,  &
         -0.9533E-07_r8, -0.1359E-06_r8, -0.1879E-06_r8, -0.2570E-06_r8, -0.3470E-06_r8,  &
         -0.4598E-06_r8, -0.5971E-06_r8, -0.7563E-06_r8, -0.9334E-06_r8, -0.1127E-05_r8,  &
         -0.1336E-05_r8, -0.1543E-05_r8, -0.1701E-05_r8, -0.1738E-05_r8, &
                                ! DATA ((c1(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip= 8_r8, 8)_r8, &                                &
         0.1460E-03_r8,  0.2420E-03_r8,  0.3740E-03_r8,  0.5443E-03_r8,  0.7681E-03_r8,   &
         0.1081E-02_r8,  0.1535E-02_r8,  0.2184E-02_r8,  0.3082E-02_r8,  0.4286E-02_r8,   &
         0.5885E-02_r8,  0.8011E-02_r8,  0.1083E-01_r8,  0.1458E-01_r8,  0.1954E-01_r8,   &
         0.2608E-01_r8,  0.3468E-01_r8,  0.4595E-01_r8,  0.6068E-01_r8,  0.7979E-01_r8,   &
         0.1043E+00_r8,  0.1351E+00_r8,  0.1731E+00_r8,  0.2184E+00_r8, &
                                ! DATA ((c2(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip= 8_r8, 8)_r8, &                                &
         -0.2770E-06_r8, -0.6058E-06_r8, -0.1224E-05_r8, -0.2243E-05_r8, -0.3777E-05_r8,  &
         -0.6018E-05_r8, -0.9181E-05_r8, -0.1339E-04_r8, -0.1876E-04_r8, -0.2561E-04_r8,  &
         -0.3467E-04_r8, -0.4695E-04_r8, -0.6370E-04_r8, -0.8635E-04_r8, -0.1163E-03_r8,  &
         -0.1551E-03_r8, -0.2043E-03_r8, -0.2654E-03_r8, -0.3400E-03_r8, -0.4299E-03_r8,  &
         -0.5366E-03_r8, -0.6609E-03_r8, -0.8010E-03_r8, -0.9521E-03_r8, &
                                ! DATA ((c3(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip= 8_r8, 8)_r8, &                                &
         -0.2845E-09_r8, -0.6367E-09_r8, -0.1591E-08_r8, -0.3779E-08_r8, -0.7752E-08_r8,  &
         -0.1370E-07_r8, -0.2140E-07_r8, -0.3101E-07_r8, -0.4437E-07_r8, -0.6498E-07_r8,  &
         -0.9551E-07_r8, -0.1360E-06_r8, -0.1887E-06_r8, -0.2590E-06_r8, -0.3503E-06_r8,  &
         -0.4650E-06_r8, -0.6052E-06_r8, -0.7681E-06_r8, -0.9488E-06_r8, -0.1145E-05_r8,  &
         -0.1354E-05_r8, -0.1555E-05_r8, -0.1698E-05_r8, -0.1701E-05_r8, &
                                ! DATA ((c1(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip= 9_r8, 9)_r8, &                                &
         0.1470E-03_r8,  0.2444E-03_r8,  0.3796E-03_r8,  0.5557E-03_r8,  0.7894E-03_r8,   &
         0.1119E-02_r8,  0.1600E-02_r8,  0.2296E-02_r8,  0.3271E-02_r8,  0.4596E-02_r8,   &
         0.6373E-02_r8,  0.8752E-02_r8,  0.1193E-01_r8,  0.1617E-01_r8,  0.2181E-01_r8,   &
         0.2928E-01_r8,  0.3916E-01_r8,  0.5220E-01_r8,  0.6930E-01_r8,  0.9150E-01_r8,   &
         0.1198E+00_r8,  0.1551E+00_r8,  0.1978E+00_r8,  0.2475E+00_r8, &
                                ! DATA ((c2(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip= 9_r8, 9)_r8, &                                &
         -0.2709E-06_r8, -0.5925E-06_r8, -0.1198E-05_r8, -0.2199E-05_r8, -0.3708E-05_r8,  &
         -0.5913E-05_r8, -0.9036E-05_r8, -0.1322E-04_r8, -0.1861E-04_r8, -0.2556E-04_r8,  &
         -0.3486E-04_r8, -0.4760E-04_r8, -0.6516E-04_r8, -0.8909E-04_r8, -0.1210E-03_r8,  &
         -0.1624E-03_r8, -0.2151E-03_r8, -0.2808E-03_r8, -0.3615E-03_r8, -0.4590E-03_r8,  &
         -0.5747E-03_r8, -0.7088E-03_r8, -0.8585E-03_r8, -0.1018E-02_r8, &
                                ! DATA ((c3(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip= 9_r8, 9)_r8, &                                &
         -0.3119E-09_r8, -0.6873E-09_r8, -0.1676E-08_r8, -0.3915E-08_r8, -0.7973E-08_r8,  &
         -0.1407E-07_r8, -0.2200E-07_r8, -0.3179E-07_r8, -0.4505E-07_r8, -0.6522E-07_r8,  &
         -0.9551E-07_r8, -0.1364E-06_r8, -0.1902E-06_r8, -0.2622E-06_r8, -0.3557E-06_r8,  &
         -0.4735E-06_r8, -0.6178E-06_r8, -0.7853E-06_r8, -0.9700E-06_r8, -0.1169E-05_r8,  &
         -0.1377E-05_r8, -0.1567E-05_r8, -0.1683E-05_r8, -0.1631E-05_r8, &
                                ! DATA ((c1(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=10_r8,10)_r8, &                                &
         0.1484E-03_r8,  0.2482E-03_r8,  0.3882E-03_r8,  0.5732E-03_r8,  0.8217E-03_r8,   &
         0.1176E-02_r8,  0.1696E-02_r8,  0.2457E-02_r8,  0.3536E-02_r8,  0.5021E-02_r8,   &
         0.7028E-02_r8,  0.9728E-02_r8,  0.1336E-01_r8,  0.1821E-01_r8,  0.2470E-01_r8,   &
         0.3334E-01_r8,  0.4482E-01_r8,  0.6004E-01_r8,  0.8003E-01_r8,  0.1059E+00_r8,   &
         0.1385E+00_r8,  0.1785E+00_r8,  0.2259E+00_r8,  0.2794E+00_r8, &
                                ! DATA ((c2(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=10_r8,10)_r8, &                                &
         -0.2620E-06_r8, -0.5728E-06_r8, -0.1161E-05_r8, -0.2136E-05_r8, -0.3612E-05_r8,  &
         -0.5777E-05_r8, -0.8867E-05_r8, -0.1306E-04_r8, -0.1854E-04_r8, -0.2572E-04_r8,  &
         -0.3543E-04_r8, -0.4887E-04_r8, -0.6756E-04_r8, -0.9321E-04_r8, -0.1275E-03_r8,  &
         -0.1723E-03_r8, -0.2295E-03_r8, -0.3011E-03_r8, -0.3894E-03_r8, -0.4960E-03_r8,  &
         -0.6220E-03_r8, -0.7665E-03_r8, -0.9259E-03_r8, -0.1095E-02_r8, &
                                ! DATA ((c3(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=10_r8,10)_r8, &                                &
         -0.3497E-09_r8, -0.7575E-09_r8, -0.1791E-08_r8, -0.4091E-08_r8, -0.8254E-08_r8,  &
         -0.1452E-07_r8, -0.2265E-07_r8, -0.3255E-07_r8, -0.4557E-07_r8, -0.6530E-07_r8,  &
         -0.9557E-07_r8, -0.1374E-06_r8, -0.1930E-06_r8, -0.2675E-06_r8, -0.3644E-06_r8,  &
         -0.4866E-06_r8, -0.6362E-06_r8, -0.8090E-06_r8, -0.9981E-06_r8, -0.1199E-05_r8,  &
         -0.1402E-05_r8, -0.1574E-05_r8, -0.1643E-05_r8, -0.1506E-05_r8, &
                                ! DATA ((c1(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=11_r8,11)_r8, &                                &
         0.1506E-03_r8,  0.2538E-03_r8,  0.4010E-03_r8,  0.5992E-03_r8,  0.8693E-03_r8,   &
         0.1258E-02_r8,  0.1834E-02_r8,  0.2683E-02_r8,  0.3901E-02_r8,  0.5591E-02_r8,   &
         0.7890E-02_r8,  0.1100E-01_r8,  0.1519E-01_r8,  0.2082E-01_r8,  0.2837E-01_r8,   &
         0.3845E-01_r8,  0.5192E-01_r8,  0.6981E-01_r8,  0.9323E-01_r8,  0.1232E+00_r8,   &
         0.1606E+00_r8,  0.2055E+00_r8,  0.2571E+00_r8,  0.3133E+00_r8, &
                                ! DATA ((c2(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=11_r8,11)_r8, &                                &
         -0.2494E-06_r8, -0.5450E-06_r8, -0.1108E-05_r8, -0.2049E-05_r8, -0.3485E-05_r8,  &
         -0.5612E-05_r8, -0.8695E-05_r8, -0.1296E-04_r8, -0.1865E-04_r8, -0.2620E-04_r8,  &
         -0.3654E-04_r8, -0.5099E-04_r8, -0.7121E-04_r8, -0.9909E-04_r8, -0.1365E-03_r8,  &
         -0.1855E-03_r8, -0.2483E-03_r8, -0.3274E-03_r8, -0.4248E-03_r8, -0.5421E-03_r8,  &
         -0.6794E-03_r8, -0.8346E-03_r8, -0.1003E-02_r8, -0.1181E-02_r8, &
                                ! DATA ((c3(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=11_r8,11)_r8, &                                &
         -0.3973E-09_r8, -0.8499E-09_r8, -0.1944E-08_r8, -0.4319E-08_r8, -0.8596E-08_r8,  &
         -0.1501E-07_r8, -0.2332E-07_r8, -0.3326E-07_r8, -0.4600E-07_r8, -0.6541E-07_r8,  &
         -0.9592E-07_r8, -0.1392E-06_r8, -0.1975E-06_r8, -0.2757E-06_r8, -0.3774E-06_r8,  &
         -0.5057E-06_r8, -0.6614E-06_r8, -0.8401E-06_r8, -0.1034E-05_r8, -0.1235E-05_r8,  &
         -0.1427E-05_r8, -0.1565E-05_r8, -0.1560E-05_r8, -0.1310E-05_r8, &
                                ! DATA ((c1(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=12_r8,12)_r8, &                                &
         0.1537E-03_r8,  0.2618E-03_r8,  0.4196E-03_r8,  0.6368E-03_r8,  0.9377E-03_r8,   &
         0.1374E-02_r8,  0.2024E-02_r8,  0.2990E-02_r8,  0.4387E-02_r8,  0.6339E-02_r8,   &
         0.9008E-02_r8,  0.1263E-01_r8,  0.1752E-01_r8,  0.2410E-01_r8,  0.3296E-01_r8,   &
         0.4485E-01_r8,  0.6077E-01_r8,  0.8185E-01_r8,  0.1092E+00_r8,  0.1439E+00_r8,   &
         0.1861E+00_r8,  0.2356E+00_r8,  0.2906E+00_r8,  0.3483E+00_r8, &
                                ! DATA ((c2(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=12_r8,12)_r8, &                                &
         -0.2326E-06_r8, -0.5076E-06_r8, -0.1037E-05_r8, -0.1935E-05_r8, -0.3329E-05_r8,  &
         -0.5433E-05_r8, -0.8553E-05_r8, -0.1298E-04_r8, -0.1902E-04_r8, -0.2716E-04_r8,  &
         -0.3840E-04_r8, -0.5421E-04_r8, -0.7643E-04_r8, -0.1071E-03_r8, -0.1485E-03_r8,  &
         -0.2028E-03_r8, -0.2727E-03_r8, -0.3608E-03_r8, -0.4689E-03_r8, -0.5981E-03_r8,  &
         -0.7471E-03_r8, -0.9128E-03_r8, -0.1091E-02_r8, -0.1276E-02_r8, &
                                ! DATA ((c3(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=12_r8,12)_r8, &                                &
         -0.4576E-09_r8, -0.9621E-09_r8, -0.2122E-08_r8, -0.4569E-08_r8, -0.8941E-08_r8,  &
         -0.1550E-07_r8, -0.2398E-07_r8, -0.3404E-07_r8, -0.4663E-07_r8, -0.6577E-07_r8,  &
         -0.9677E-07_r8, -0.1421E-06_r8, -0.2042E-06_r8, -0.2878E-06_r8, -0.3962E-06_r8,  &
         -0.5318E-06_r8, -0.6946E-06_r8, -0.8798E-06_r8, -0.1078E-05_r8, -0.1275E-05_r8,  &
         -0.1445E-05_r8, -0.1523E-05_r8, -0.1413E-05_r8, -0.1030E-05_r8, &
                                ! DATA ((c1(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=13_r8,13)_r8, &                                &
         0.1579E-03_r8,  0.2729E-03_r8,  0.4454E-03_r8,  0.6893E-03_r8,  0.1033E-02_r8,   &
         0.1533E-02_r8,  0.2280E-02_r8,  0.3396E-02_r8,  0.5018E-02_r8,  0.7299E-02_r8,   &
         0.1043E-01_r8,  0.1468E-01_r8,  0.2044E-01_r8,  0.2820E-01_r8,  0.3867E-01_r8,   &
         0.5278E-01_r8,  0.7166E-01_r8,  0.9648E-01_r8,  0.1283E+00_r8,  0.1679E+00_r8,   &
         0.2150E+00_r8,  0.2685E+00_r8,  0.3256E+00_r8,  0.3836E+00_r8, &
                                ! DATA ((c2(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=13_r8,13)_r8, &                                &
         -0.2120E-06_r8, -0.4613E-06_r8, -0.9487E-06_r8, -0.1798E-05_r8, -0.3153E-05_r8,  &
         -0.5260E-05_r8, -0.8480E-05_r8, -0.1319E-04_r8, -0.1976E-04_r8, -0.2876E-04_r8,  &
         -0.4125E-04_r8, -0.5885E-04_r8, -0.8357E-04_r8, -0.1178E-03_r8, -0.1640E-03_r8,  &
         -0.2250E-03_r8, -0.3036E-03_r8, -0.4023E-03_r8, -0.5226E-03_r8, -0.6643E-03_r8,  &
         -0.8250E-03_r8, -0.1001E-02_r8, -0.1187E-02_r8, -0.1378E-02_r8, &
                                ! DATA ((c3(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=13_r8,13)_r8, &                                &
         -0.5223E-09_r8, -0.1082E-08_r8, -0.2314E-08_r8, -0.4826E-08_r8, -0.9255E-08_r8,  &
         -0.1596E-07_r8, -0.2472E-07_r8, -0.3507E-07_r8, -0.4774E-07_r8, -0.6677E-07_r8,  &
         -0.9841E-07_r8, -0.1466E-06_r8, -0.2141E-06_r8, -0.3049E-06_r8, -0.4219E-06_r8,  &
         -0.5662E-06_r8, -0.7370E-06_r8, -0.9290E-06_r8, -0.1129E-05_r8, -0.1315E-05_r8,  &
         -0.1441E-05_r8, -0.1427E-05_r8, -0.1184E-05_r8, -0.6731E-06_r8, &
                                ! DATA ((c1(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=14_r8,14)_r8, &                                &
         0.1634E-03_r8,  0.2873E-03_r8,  0.4793E-03_r8,  0.7589E-03_r8,  0.1159E-02_r8,   &
         0.1742E-02_r8,  0.2614E-02_r8,  0.3916E-02_r8,  0.5817E-02_r8,  0.8503E-02_r8,   &
         0.1220E-01_r8,  0.1723E-01_r8,  0.2403E-01_r8,  0.3324E-01_r8,  0.4572E-01_r8,   &
         0.6253E-01_r8,  0.8488E-01_r8,  0.1140E+00_r8,  0.1507E+00_r8,  0.1952E+00_r8,   &
         0.2468E+00_r8,  0.3031E+00_r8,  0.3612E+00_r8,  0.4186E+00_r8, &
                                ! DATA ((c2(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=14_r8,14)_r8, &                                &
         -0.1896E-06_r8, -0.4100E-06_r8, -0.8500E-06_r8, -0.1645E-05_r8, -0.2972E-05_r8,  &
         -0.5116E-05_r8, -0.8508E-05_r8, -0.1363E-04_r8, -0.2097E-04_r8, -0.3115E-04_r8,  &
         -0.4531E-04_r8, -0.6517E-04_r8, -0.9298E-04_r8, -0.1315E-03_r8, -0.1837E-03_r8,  &
         -0.2528E-03_r8, -0.3417E-03_r8, -0.4525E-03_r8, -0.5860E-03_r8, -0.7405E-03_r8,  &
         -0.9123E-03_r8, -0.1097E-02_r8, -0.1290E-02_r8, -0.1482E-02_r8, &
                                ! DATA ((c3(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=14_r8,14)_r8, &                                &
         -0.5738E-09_r8, -0.1185E-08_r8, -0.2465E-08_r8, -0.5015E-08_r8, -0.9502E-08_r8,  &
         -0.1640E-07_r8, -0.2567E-07_r8, -0.3675E-07_r8, -0.4990E-07_r8, -0.6893E-07_r8,  &
         -0.1014E-06_r8, -0.1533E-06_r8, -0.2278E-06_r8, -0.3278E-06_r8, -0.4553E-06_r8,  &
         -0.6098E-06_r8, -0.7899E-06_r8, -0.9876E-06_r8, -0.1184E-05_r8, -0.1342E-05_r8,  &
         -0.1393E-05_r8, -0.1254E-05_r8, -0.8703E-06_r8, -0.2593E-06_r8, &
                                ! DATA ((c1(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=15_r8,15)_r8, &                                &
         0.1697E-03_r8,  0.3043E-03_r8,  0.5205E-03_r8,  0.8458E-03_r8,  0.1319E-02_r8,   &
         0.2009E-02_r8,  0.3034E-02_r8,  0.4563E-02_r8,  0.6800E-02_r8,  0.9975E-02_r8,   &
         0.1436E-01_r8,  0.2033E-01_r8,  0.2842E-01_r8,  0.3941E-01_r8,  0.5431E-01_r8,   &
         0.7432E-01_r8,  0.1007E+00_r8,  0.1344E+00_r8,  0.1762E+00_r8,  0.2255E+00_r8,   &
         0.2806E+00_r8,  0.3385E+00_r8,  0.3964E+00_r8,  0.4532E+00_r8, &
                                ! DATA ((c2(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=15_r8,15)_r8, &                                &
         -0.1689E-06_r8, -0.3610E-06_r8, -0.7531E-06_r8, -0.1495E-05_r8, -0.2803E-05_r8,  &
         -0.5013E-05_r8, -0.8646E-05_r8, -0.1433E-04_r8, -0.2270E-04_r8, -0.3446E-04_r8,  &
         -0.5077E-04_r8, -0.7342E-04_r8, -0.1049E-03_r8, -0.1486E-03_r8, -0.2080E-03_r8,  &
         -0.2866E-03_r8, -0.3873E-03_r8, -0.5115E-03_r8, -0.6586E-03_r8, -0.8256E-03_r8,  &
         -0.1008E-02_r8, -0.1201E-02_r8, -0.1398E-02_r8, -0.1583E-02_r8, &
                                ! DATA ((c3(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=15_r8,15)_r8, &                                &
         -0.6063E-09_r8, -0.1245E-08_r8, -0.2559E-08_r8, -0.5110E-08_r8, -0.9652E-08_r8,  &
         -0.1685E-07_r8, -0.2693E-07_r8, -0.3933E-07_r8, -0.5373E-07_r8, -0.7330E-07_r8,  &
         -0.1066E-06_r8, -0.1625E-06_r8, -0.2453E-06_r8, -0.3564E-06_r8, -0.4963E-06_r8,  &
         -0.6636E-06_r8, -0.8532E-06_r8, -0.1054E-05_r8, -0.1236E-05_r8, -0.1339E-05_r8,  &
         -0.1280E-05_r8, -0.9947E-06_r8, -0.4862E-06_r8,  0.1889E-06_r8, &
                                ! DATA ((c1(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=16_r8,16)_r8, &                                &
         0.1763E-03_r8,  0.3227E-03_r8,  0.5664E-03_r8,  0.9465E-03_r8,  0.1511E-02_r8,   &
         0.2335E-02_r8,  0.3549E-02_r8,  0.5348E-02_r8,  0.7977E-02_r8,  0.1172E-01_r8,   &
         0.1692E-01_r8,  0.2402E-01_r8,  0.3369E-01_r8,  0.4685E-01_r8,  0.6465E-01_r8,   &
         0.8835E-01_r8,  0.1192E+00_r8,  0.1579E+00_r8,  0.2046E+00_r8,  0.2579E+00_r8,   &
         0.3153E+00_r8,  0.3736E+00_r8,  0.4308E+00_r8,  0.4875E+00_r8, &
                                ! DATA ((c2(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=16_r8,16)_r8, &                                &
         -0.1533E-06_r8, -0.3224E-06_r8, -0.6725E-06_r8, -0.1365E-05_r8, -0.2657E-05_r8,  &
         -0.4953E-05_r8, -0.8878E-05_r8, -0.1525E-04_r8, -0.2490E-04_r8, -0.3866E-04_r8,  &
         -0.5766E-04_r8, -0.8372E-04_r8, -0.1197E-03_r8, -0.1695E-03_r8, -0.2372E-03_r8,  &
         -0.3266E-03_r8, -0.4401E-03_r8, -0.5783E-03_r8, -0.7390E-03_r8, -0.9180E-03_r8,  &
         -0.1110E-02_r8, -0.1309E-02_r8, -0.1503E-02_r8, -0.1676E-02_r8, &
                                ! DATA ((c3(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=16_r8,16)_r8, &                                &
         -0.6109E-09_r8, -0.1253E-08_r8, -0.2555E-08_r8, -0.5100E-08_r8, -0.9703E-08_r8,  &
         -0.1729E-07_r8, -0.2846E-07_r8, -0.4291E-07_r8, -0.5988E-07_r8, -0.8135E-07_r8,  &
         -0.1158E-06_r8, -0.1750E-06_r8, -0.2657E-06_r8, -0.3890E-06_r8, -0.5438E-06_r8,  &
         -0.7262E-06_r8, -0.9256E-06_r8, -0.1122E-05_r8, -0.1269E-05_r8, -0.1284E-05_r8,  &
         -0.1086E-05_r8, -0.6553E-06_r8, -0.5588E-07_r8,  0.6660E-06_r8, &
                                ! DATA ((c1(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=17_r8,17)_r8, &                                &
         0.1824E-03_r8,  0.3403E-03_r8,  0.6126E-03_r8,  0.1054E-02_r8,  0.1727E-02_r8,   &
         0.2716E-02_r8,  0.4160E-02_r8,  0.6275E-02_r8,  0.9353E-02_r8,  0.1375E-01_r8,   &
         0.1990E-01_r8,  0.2837E-01_r8,  0.3995E-01_r8,  0.5571E-01_r8,  0.7690E-01_r8,   &
         0.1048E+00_r8,  0.1404E+00_r8,  0.1841E+00_r8,  0.2353E+00_r8,  0.2916E+00_r8,   &
         0.3500E+00_r8,  0.4078E+00_r8,  0.4646E+00_r8,  0.5215E+00_r8, &
                                ! DATA ((c2(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=17_r8,17)_r8, &                                &
         -0.1436E-06_r8, -0.2978E-06_r8, -0.6179E-06_r8, -0.1270E-05_r8, -0.2546E-05_r8,  &
         -0.4927E-05_r8, -0.9164E-05_r8, -0.1628E-04_r8, -0.2740E-04_r8, -0.4353E-04_r8,  &
         -0.6579E-04_r8, -0.9606E-04_r8, -0.1375E-03_r8, -0.1943E-03_r8, -0.2712E-03_r8,  &
         -0.3722E-03_r8, -0.4993E-03_r8, -0.6515E-03_r8, -0.8255E-03_r8, -0.1016E-02_r8,  &
         -0.1217E-02_r8, -0.1417E-02_r8, -0.1602E-02_r8, -0.1757E-02_r8, &
                                ! DATA ((c3(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=17_r8,17)_r8, &                                &
         -0.6015E-09_r8, -0.1225E-08_r8, -0.2489E-08_r8, -0.4978E-08_r8, -0.9608E-08_r8,  &
         -0.1756E-07_r8, -0.2993E-07_r8, -0.4715E-07_r8, -0.6853E-07_r8, -0.9454E-07_r8,  &
         -0.1314E-06_r8, -0.1920E-06_r8, -0.2877E-06_r8, -0.4222E-06_r8, -0.5940E-06_r8,  &
         -0.7936E-06_r8, -0.1002E-05_r8, -0.1181E-05_r8, -0.1262E-05_r8, -0.1156E-05_r8,  &
         -0.8086E-06_r8, -0.2604E-06_r8,  0.4056E-06_r8,  0.1164E-05_r8, &
                                ! DATA ((c1(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=18_r8,18)_r8, &                                &
         0.1875E-03_r8,  0.3555E-03_r8,  0.6547E-03_r8,  0.1159E-02_r8,  0.1954E-02_r8,   &
         0.3140E-02_r8,  0.4861E-02_r8,  0.7345E-02_r8,  0.1093E-01_r8,  0.1606E-01_r8,   &
         0.2330E-01_r8,  0.3339E-01_r8,  0.4729E-01_r8,  0.6613E-01_r8,  0.9118E-01_r8,   &
         0.1236E+00_r8,  0.1642E+00_r8,  0.2127E+00_r8,  0.2674E+00_r8,  0.3254E+00_r8,   &
         0.3836E+00_r8,  0.4408E+00_r8,  0.4978E+00_r8,  0.5551E+00_r8, &
                                ! DATA ((c2(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=18_r8,18)_r8, &                                &
         -0.1390E-06_r8, -0.2850E-06_r8, -0.5878E-06_r8, -0.1213E-05_r8, -0.2477E-05_r8,  &
         -0.4934E-05_r8, -0.9466E-05_r8, -0.1731E-04_r8, -0.2990E-04_r8, -0.4860E-04_r8,  &
         -0.7471E-04_r8, -0.1102E-03_r8, -0.1582E-03_r8, -0.2230E-03_r8, -0.3098E-03_r8,  &
         -0.4228E-03_r8, -0.5632E-03_r8, -0.7291E-03_r8, -0.9158E-03_r8, -0.1117E-02_r8,  &
         -0.1323E-02_r8, -0.1519E-02_r8, -0.1689E-02_r8, -0.1824E-02_r8, &
                                ! DATA ((c3(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=18_r8,18)_r8, &                                &
         -0.5875E-09_r8, -0.1187E-08_r8, -0.2392E-08_r8, -0.4777E-08_r8, -0.9321E-08_r8,  &
         -0.1745E-07_r8, -0.3091E-07_r8, -0.5121E-07_r8, -0.7879E-07_r8, -0.1130E-06_r8,  &
         -0.1555E-06_r8, -0.2161E-06_r8, -0.3110E-06_r8, -0.4527E-06_r8, -0.6400E-06_r8,  &
         -0.8571E-06_r8, -0.1069E-05_r8, -0.1209E-05_r8, -0.1193E-05_r8, -0.9472E-06_r8,  &
         -0.4712E-06_r8,  0.1597E-06_r8,  0.8892E-06_r8,  0.1671E-05_r8, &
                                ! DATA ((c1(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=19_r8,19)_r8, &                                &
         0.1914E-03_r8,  0.3677E-03_r8,  0.6900E-03_r8,  0.1252E-02_r8,  0.2172E-02_r8,   &
         0.3581E-02_r8,  0.5632E-02_r8,  0.8553E-02_r8,  0.1271E-01_r8,  0.1866E-01_r8,   &
         0.2715E-01_r8,  0.3914E-01_r8,  0.5576E-01_r8,  0.7817E-01_r8,  0.1076E+00_r8,   &
         0.1449E+00_r8,  0.1903E+00_r8,  0.2429E+00_r8,  0.2999E+00_r8,  0.3582E+00_r8,   &
         0.4157E+00_r8,  0.4727E+00_r8,  0.5303E+00_r8,  0.5878E+00_r8, &
                                ! DATA ((c2(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=19_r8,19)_r8, &                                &
         -0.1371E-06_r8, -0.2799E-06_r8, -0.5756E-06_r8, -0.1190E-05_r8, -0.2455E-05_r8,  &
         -0.4982E-05_r8, -0.9770E-05_r8, -0.1824E-04_r8, -0.3214E-04_r8, -0.5332E-04_r8,  &
         -0.8371E-04_r8, -0.1255E-03_r8, -0.1814E-03_r8, -0.2555E-03_r8, -0.3528E-03_r8,  &
         -0.4774E-03_r8, -0.6306E-03_r8, -0.8095E-03_r8, -0.1008E-02_r8, -0.1218E-02_r8,  &
         -0.1424E-02_r8, -0.1610E-02_r8, -0.1763E-02_r8, -0.1875E-02_r8, &
                                ! DATA ((c3(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=19_r8,19)_r8, &                                &
         -0.5739E-09_r8, -0.1144E-08_r8, -0.2283E-08_r8, -0.4524E-08_r8, -0.8848E-08_r8,  &
         -0.1686E-07_r8, -0.3096E-07_r8, -0.5411E-07_r8, -0.8876E-07_r8, -0.1343E-06_r8,  &
         -0.1880E-06_r8, -0.2511E-06_r8, -0.3405E-06_r8, -0.4800E-06_r8, -0.6756E-06_r8,  &
         -0.9041E-06_r8, -0.1105E-05_r8, -0.1184E-05_r8, -0.1049E-05_r8, -0.6701E-06_r8,  &
         -0.1061E-06_r8,  0.5899E-06_r8,  0.1376E-05_r8,  0.2160E-05_r8, &
                                ! DATA ((c1(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=20_r8,20)_r8, &                                &
         0.1942E-03_r8,  0.3767E-03_r8,  0.7173E-03_r8,  0.1328E-02_r8,  0.2364E-02_r8,   &
         0.4005E-02_r8,  0.6433E-02_r8,  0.9873E-02_r8,  0.1471E-01_r8,  0.2160E-01_r8,   &
         0.3150E-01_r8,  0.4565E-01_r8,  0.6537E-01_r8,  0.9183E-01_r8,  0.1260E+00_r8,   &
         0.1683E+00_r8,  0.2183E+00_r8,  0.2738E+00_r8,  0.3318E+00_r8,  0.3894E+00_r8,   &
         0.4463E+00_r8,  0.5038E+00_r8,  0.5619E+00_r8,  0.6191E+00_r8, &
                                ! DATA ((c2(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=20_r8,20)_r8, &                                &
         -0.1368E-06_r8, -0.2791E-06_r8, -0.5744E-06_r8, -0.1192E-05_r8, -0.2475E-05_r8,  &
         -0.5078E-05_r8, -0.1008E-04_r8, -0.1903E-04_r8, -0.3392E-04_r8, -0.5726E-04_r8,  &
         -0.9202E-04_r8, -0.1409E-03_r8, -0.2063E-03_r8, -0.2912E-03_r8, -0.3997E-03_r8,  &
         -0.5357E-03_r8, -0.7004E-03_r8, -0.8911E-03_r8, -0.1100E-02_r8, -0.1314E-02_r8,  &
         -0.1514E-02_r8, -0.1687E-02_r8, -0.1822E-02_r8, -0.1908E-02_r8, &
                                ! DATA ((c3(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=20_r8,20)_r8, &                                &
         -0.5617E-09_r8, -0.1109E-08_r8, -0.2184E-08_r8, -0.4264E-08_r8, -0.8251E-08_r8,  &
         -0.1582E-07_r8, -0.3000E-07_r8, -0.5536E-07_r8, -0.9656E-07_r8, -0.1549E-06_r8,  &
         -0.2254E-06_r8, -0.2999E-06_r8, -0.3855E-06_r8, -0.5114E-06_r8, -0.6982E-06_r8,  &
         -0.9171E-06_r8, -0.1082E-05_r8, -0.1078E-05_r8, -0.8280E-06_r8, -0.3561E-06_r8,  &
         0.2636E-06_r8,  0.1015E-05_r8,  0.1836E-05_r8,  0.2588E-05_r8, &
                                ! DATA ((c1(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=21_r8,21)_r8, &                                &
         0.1962E-03_r8,  0.3831E-03_r8,  0.7370E-03_r8,  0.1386E-02_r8,  0.2520E-02_r8,   &
         0.4379E-02_r8,  0.7206E-02_r8,  0.1125E-01_r8,  0.1692E-01_r8,  0.2493E-01_r8,   &
         0.3644E-01_r8,  0.5299E-01_r8,  0.7611E-01_r8,  0.1070E+00_r8,  0.1463E+00_r8,   &
         0.1936E+00_r8,  0.2474E+00_r8,  0.3046E+00_r8,  0.3621E+00_r8,  0.4188E+00_r8,   &
         0.4758E+00_r8,  0.5339E+00_r8,  0.5921E+00_r8,  0.6482E+00_r8, &
                                ! DATA ((c2(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=21_r8,21)_r8, &                                &
         -0.1371E-06_r8, -0.2803E-06_r8, -0.5786E-06_r8, -0.1206E-05_r8, -0.2522E-05_r8,  &
         -0.5211E-05_r8, -0.1040E-04_r8, -0.1969E-04_r8, -0.3522E-04_r8, -0.6021E-04_r8,  &
         -0.9890E-04_r8, -0.1552E-03_r8, -0.2312E-03_r8, -0.3285E-03_r8, -0.4496E-03_r8,  &
         -0.5973E-03_r8, -0.7728E-03_r8, -0.9733E-03_r8, -0.1188E-02_r8, -0.1399E-02_r8,  &
         -0.1590E-02_r8, -0.1748E-02_r8, -0.1862E-02_r8, -0.1922E-02_r8, &
                                ! DATA ((c3(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=21_r8,21)_r8, &                                &
         -0.5521E-09_r8, -0.1082E-08_r8, -0.2100E-08_r8, -0.4025E-08_r8, -0.7661E-08_r8,  &
         -0.1462E-07_r8, -0.2837E-07_r8, -0.5493E-07_r8, -0.1013E-06_r8, -0.1716E-06_r8,  &
         -0.2627E-06_r8, -0.3609E-06_r8, -0.4548E-06_r8, -0.5617E-06_r8, -0.7128E-06_r8,  &
         -0.8814E-06_r8, -0.9684E-06_r8, -0.8687E-06_r8, -0.5417E-06_r8, -0.3260E-07_r8,  &
         0.6298E-06_r8,  0.1417E-05_r8,  0.2233E-05_r8,  0.2884E-05_r8, &
                                ! DATA ((c1(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=22_r8,22)_r8, &                                &
         0.1974E-03_r8,  0.3872E-03_r8,  0.7505E-03_r8,  0.1427E-02_r8,  0.2638E-02_r8,   &
         0.4682E-02_r8,  0.7896E-02_r8,  0.1261E-01_r8,  0.1927E-01_r8,  0.2865E-01_r8,   &
         0.4207E-01_r8,  0.6125E-01_r8,  0.8802E-01_r8,  0.1236E+00_r8,  0.1681E+00_r8,   &
         0.2201E+00_r8,  0.2766E+00_r8,  0.3342E+00_r8,  0.3906E+00_r8,  0.4468E+00_r8,   &
         0.5043E+00_r8,  0.5628E+00_r8,  0.6204E+00_r8,  0.6743E+00_r8, &
                                ! DATA ((c2(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=22_r8,22)_r8, &                                &
         -0.1377E-06_r8, -0.2818E-06_r8, -0.5836E-06_r8, -0.1222E-05_r8, -0.2570E-05_r8,  &
         -0.5339E-05_r8, -0.1069E-04_r8, -0.2019E-04_r8, -0.3606E-04_r8, -0.6211E-04_r8,  &
         -0.1038E-03_r8, -0.1666E-03_r8, -0.2537E-03_r8, -0.3653E-03_r8, -0.5012E-03_r8,  &
         -0.6623E-03_r8, -0.8484E-03_r8, -0.1055E-02_r8, -0.1268E-02_r8, -0.1470E-02_r8,  &
         -0.1649E-02_r8, -0.1790E-02_r8, -0.1882E-02_r8, -0.1920E-02_r8, &
                                ! DATA ((c3(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=22_r8,22)_r8, &                                &
         -0.5464E-09_r8, -0.1064E-08_r8, -0.2050E-08_r8, -0.3874E-08_r8, -0.7265E-08_r8,  &
         -0.1373E-07_r8, -0.2688E-07_r8, -0.5376E-07_r8, -0.1035E-06_r8, -0.1833E-06_r8,  &
         -0.2955E-06_r8, -0.4267E-06_r8, -0.5474E-06_r8, -0.6443E-06_r8, -0.7332E-06_r8,  &
         -0.7933E-06_r8, -0.7465E-06_r8, -0.5424E-06_r8, -0.1941E-06_r8,  0.3078E-06_r8,  &
         0.9845E-06_r8,  0.1771E-05_r8,  0.2515E-05_r8,  0.2962E-05_r8, &
                                ! DATA ((c1(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=23_r8,23)_r8, &                                &
         0.1981E-03_r8,  0.3899E-03_r8,  0.7594E-03_r8,  0.1455E-02_r8,  0.2722E-02_r8,   &
         0.4913E-02_r8,  0.8469E-02_r8,  0.1387E-01_r8,  0.2168E-01_r8,  0.3275E-01_r8,   &
         0.4842E-01_r8,  0.7058E-01_r8,  0.1012E+00_r8,  0.1414E+00_r8,  0.1908E+00_r8,   &
         0.2467E+00_r8,  0.3050E+00_r8,  0.3620E+00_r8,  0.4176E+00_r8,  0.4739E+00_r8,   &
         0.5318E+00_r8,  0.5900E+00_r8,  0.6460E+00_r8,  0.6969E+00_r8, &
                                ! DATA ((c2(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=23_r8,23)_r8, &                                &
         -0.1374E-06_r8, -0.2821E-06_r8, -0.5849E-06_r8, -0.1227E-05_r8, -0.2587E-05_r8,  &
         -0.5392E-05_r8, -0.1082E-04_r8, -0.2043E-04_r8, -0.3640E-04_r8, -0.6291E-04_r8,  &
         -0.1064E-03_r8, -0.1742E-03_r8, -0.2713E-03_r8, -0.3983E-03_r8, -0.5523E-03_r8,  &
         -0.7302E-03_r8, -0.9280E-03_r8, -0.1135E-02_r8, -0.1338E-02_r8, -0.1526E-02_r8,  &
         -0.1690E-02_r8, -0.1811E-02_r8, -0.1880E-02_r8, -0.1904E-02_r8, &
                                ! DATA ((c3(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=23_r8,23)_r8, &                                &
         -0.5505E-09_r8, -0.1066E-08_r8, -0.2044E-08_r8, -0.3856E-08_r8, -0.7199E-08_r8,  &
         -0.1353E-07_r8, -0.2638E-07_r8, -0.5316E-07_r8, -0.1044E-06_r8, -0.1908E-06_r8,  &
         -0.3207E-06_r8, -0.4858E-06_r8, -0.6471E-06_r8, -0.7512E-06_r8, -0.7676E-06_r8,  &
         -0.6661E-06_r8, -0.4252E-06_r8, -0.1067E-06_r8,  0.2334E-06_r8,  0.6913E-06_r8,  &
         0.1330E-05_r8,  0.2055E-05_r8,  0.2617E-05_r8,  0.2760E-05_r8, &
                                ! DATA ((c1(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=24_r8,24)_r8, &                                &
         0.1986E-03_r8,  0.3917E-03_r8,  0.7654E-03_r8,  0.1475E-02_r8,  0.2782E-02_r8,   &
         0.5086E-02_r8,  0.8931E-02_r8,  0.1497E-01_r8,  0.2401E-01_r8,  0.3707E-01_r8,   &
         0.5549E-01_r8,  0.8107E-01_r8,  0.1157E+00_r8,  0.1604E+00_r8,  0.2138E+00_r8,   &
         0.2723E+00_r8,  0.3311E+00_r8,  0.3877E+00_r8,  0.4433E+00_r8,  0.5001E+00_r8,   &
         0.5580E+00_r8,  0.6152E+00_r8,  0.6683E+00_r8,  0.7158E+00_r8, &
                                ! DATA ((c2(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=24_r8,24)_r8, &                                &
         -0.1372E-06_r8, -0.2802E-06_r8, -0.5797E-06_r8, -0.1212E-05_r8, -0.2549E-05_r8,  &
         -0.5312E-05_r8, -0.1068E-04_r8, -0.2023E-04_r8, -0.3609E-04_r8, -0.6255E-04_r8,  &
         -0.1068E-03_r8, -0.1777E-03_r8, -0.2830E-03_r8, -0.4252E-03_r8, -0.6004E-03_r8,  &
         -0.8005E-03_r8, -0.1011E-02_r8, -0.1213E-02_r8, -0.1399E-02_r8, -0.1569E-02_r8,  &
         -0.1712E-02_r8, -0.1809E-02_r8, -0.1858E-02_r8, -0.1883E-02_r8, &
                                ! DATA ((c3(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=24_r8,24)_r8, &                                &
         -0.5589E-09_r8, -0.1071E-08_r8, -0.2083E-08_r8, -0.3965E-08_r8, -0.7465E-08_r8,  &
         -0.1408E-07_r8, -0.2720E-07_r8, -0.5414E-07_r8, -0.1061E-06_r8, -0.1958E-06_r8,  &
         -0.3365E-06_r8, -0.5263E-06_r8, -0.7236E-06_r8, -0.8420E-06_r8, -0.7920E-06_r8,  &
         -0.5166E-06_r8, -0.5606E-07_r8,  0.4032E-06_r8,  0.7697E-06_r8,  0.1160E-05_r8,  &
         0.1692E-05_r8,  0.2251E-05_r8,  0.2483E-05_r8,  0.2295E-05_r8, &
                                ! DATA ((c1(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=25_r8,25)_r8, &                                &
         0.1990E-03_r8,  0.3929E-03_r8,  0.7699E-03_r8,  0.1490E-02_r8,  0.2829E-02_r8,   &
         0.5226E-02_r8,  0.9314E-02_r8,  0.1593E-01_r8,  0.2618E-01_r8,  0.4144E-01_r8,   &
         0.6312E-01_r8,  0.9271E-01_r8,  0.1317E+00_r8,  0.1803E+00_r8,  0.2364E+00_r8,   &
         0.2956E+00_r8,  0.3541E+00_r8,  0.4111E+00_r8,  0.4680E+00_r8,  0.5255E+00_r8,   &
         0.5829E+00_r8,  0.6376E+00_r8,  0.6869E+00_r8,  0.7311E+00_r8, &
                                ! DATA ((c2(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=25_r8,25)_r8, &                                &
         -0.1366E-06_r8, -0.2772E-06_r8, -0.5687E-06_r8, -0.1179E-05_r8, -0.2458E-05_r8,  &
         -0.5093E-05_r8, -0.1024E-04_r8, -0.1950E-04_r8, -0.3503E-04_r8, -0.6107E-04_r8,  &
         -0.1053E-03_r8, -0.1782E-03_r8, -0.2902E-03_r8, -0.4470E-03_r8, -0.6455E-03_r8,  &
         -0.8715E-03_r8, -0.1095E-02_r8, -0.1290E-02_r8, -0.1458E-02_r8, -0.1604E-02_r8,  &
         -0.1716E-02_r8, -0.1783E-02_r8, -0.1821E-02_r8, -0.1862E-02_r8, &
                                ! DATA ((c3(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=25_r8,25)_r8, &                                &
         -0.5574E-09_r8, -0.1099E-08_r8, -0.2136E-08_r8, -0.4136E-08_r8, -0.7893E-08_r8,  &
         -0.1505E-07_r8, -0.2894E-07_r8, -0.5650E-07_r8, -0.1084E-06_r8, -0.1986E-06_r8,  &
         -0.3416E-06_r8, -0.5381E-06_r8, -0.7439E-06_r8, -0.8545E-06_r8, -0.7482E-06_r8,  &
         -0.3541E-06_r8,  0.2762E-06_r8,  0.9306E-06_r8,  0.1416E-05_r8,  0.1754E-05_r8,  &
         0.2105E-05_r8,  0.2343E-05_r8,  0.2113E-05_r8,  0.1689E-05_r8, &
                                ! DATA ((c1(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=26_r8,26)_r8, &                                &
         0.1993E-03_r8,  0.3940E-03_r8,  0.7739E-03_r8,  0.1503E-02_r8,  0.2871E-02_r8,   &
         0.5350E-02_r8,  0.9655E-02_r8,  0.1679E-01_r8,  0.2817E-01_r8,  0.4563E-01_r8,   &
         0.7094E-01_r8,  0.1052E+00_r8,  0.1488E+00_r8,  0.2006E+00_r8,  0.2575E+00_r8,   &
         0.3155E+00_r8,  0.3732E+00_r8,  0.4317E+00_r8,  0.4909E+00_r8,  0.5496E+00_r8,   &
         0.6058E+00_r8,  0.6569E+00_r8,  0.7019E+00_r8,  0.7436E+00_r8, &
                                ! DATA ((c2(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=26_r8,26)_r8, &                                &
         -0.1351E-06_r8, -0.2730E-06_r8, -0.5547E-06_r8, -0.1136E-05_r8, -0.2337E-05_r8,  &
         -0.4789E-05_r8, -0.9589E-05_r8, -0.1836E-04_r8, -0.3337E-04_r8, -0.5890E-04_r8,  &
         -0.1031E-03_r8, -0.1781E-03_r8, -0.2974E-03_r8, -0.4696E-03_r8, -0.6912E-03_r8,  &
         -0.9400E-03_r8, -0.1174E-02_r8, -0.1368E-02_r8, -0.1521E-02_r8, -0.1635E-02_r8,  &
         -0.1706E-02_r8, -0.1742E-02_r8, -0.1779E-02_r8, -0.1846E-02_r8, &
                                ! DATA ((c3(ip_r8,iw_r8,1)_r8,iw=1_r8,24)_r8,ip=26_r8,26)/                                &
         -0.5597E-09_r8, -0.1126E-08_r8, -0.2190E-08_r8, -0.4298E-08_r8, -0.8333E-08_r8,  &
         -0.1610E-07_r8, -0.3092E-07_r8, -0.5912E-07_r8, -0.1108E-06_r8, -0.1982E-06_r8,  &
         -0.3334E-06_r8, -0.5129E-06_r8, -0.6865E-06_r8, -0.7484E-06_r8, -0.5999E-06_r8,  &
         -0.1961E-06_r8,  0.5028E-06_r8,  0.1395E-05_r8,  0.2125E-05_r8,  0.2489E-05_r8,  &
         0.2584E-05_r8,  0.2312E-05_r8,  0.1587E-05_r8,  0.1122E-05_r8/),SHAPE=(/nx*3*nc/))

    it=0
    DO k=1,nx
       DO i=1,3
          DO j=1,nc
             it=it+1              
             IF(i==1)THEN
                ! WRITE(*,'(a5,2e16.9) ' )'c1',data1(it),c1(k,j)
                c1(k,j,1)=data1(it)
             END IF
             IF(i==2) THEN
                !WRITE(*,'(a5,2e16.9) ' )'c2',data1(it),c2(k,j)
                c2(k,j,1)=data1(it)
             END IF
             IF(i==3) THEN
                !WRITE(*,'(a5,2e16.9) ' )'c3',data1(it),c3(k,j)
                c3(k,j,1) = data1(it)
             END IF
          END DO
       END DO
    END DO

    !
    !-----------------------------------------------------------------------
    !
    !  The following DATA statements originally came from file
    !  "o3.tran3", which define pre-computed tables used for o3 (band5)
    !  transmittance calculations.
    !
    !-----------------------------------------------------------------------
    !
    ! o3.tran3

    !  REAL :: o1(26,21,7),o2(26,21,7),o3(26,21,7)

    !  COMMON /radtab003 /o1,o2,o3
    data2(1:nx*3*no)=RESHAPE(SOURCE=(/&
                                ! DATA ((o1(ip,iw,1),iw=1,21),ip= 1, 1)/                                &
         0.7117E-05_r8,  0.1419E-04_r8,  0.2828E-04_r8,  0.5629E-04_r8,  0.1117E-03_r8,   &
         0.2207E-03_r8,  0.4326E-03_r8,  0.8367E-03_r8,  0.1586E-02_r8,  0.2911E-02_r8,   &
         0.5106E-02_r8,  0.8433E-02_r8,  0.1302E-01_r8,  0.1896E-01_r8,  0.2658E-01_r8,   &
         0.3657E-01_r8,  0.4988E-01_r8,  0.6765E-01_r8,  0.9092E-01_r8,  0.1204E+00_r8,   &
         0.1559E+00_r8, &
                                ! DATA ((o2(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip= 1_r8, 1)_r8, &                                &
         0.6531E-10_r8,  0.5926E-10_r8, -0.1646E-09_r8, -0.1454E-08_r8, -0.7376E-08_r8,   &
         -0.2968E-07_r8, -0.1071E-06_r8, -0.3584E-06_r8, -0.1125E-05_r8, -0.3289E-05_r8,  &
         -0.8760E-05_r8, -0.2070E-04_r8, -0.4259E-04_r8, -0.7691E-04_r8, -0.1264E-03_r8,  &
         -0.1957E-03_r8, -0.2895E-03_r8, -0.4107E-03_r8, -0.5588E-03_r8, -0.7300E-03_r8,  &
         -0.9199E-03_r8, &
                                ! DATA ((o3(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip= 1_r8, 1)_r8, &                                &
         -0.2438E-10_r8, -0.4826E-10_r8, -0.9474E-10_r8, -0.1828E-09_r8, -0.3406E-09_r8,  &
         -0.6223E-09_r8, -0.1008E-08_r8, -0.1412E-08_r8, -0.1244E-08_r8,  0.8485E-09_r8,  &
         0.6343E-08_r8,  0.1201E-07_r8,  0.2838E-08_r8, -0.4024E-07_r8, -0.1257E-06_r8,   &
         -0.2566E-06_r8, -0.4298E-06_r8, -0.6184E-06_r8, -0.7657E-06_r8, -0.8153E-06_r8,  &
         -0.7552E-06_r8, &
                                ! DATA ((o1(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip= 2_r8, 2)_r8, &                                &
         0.7117E-05_r8,  0.1419E-04_r8,  0.2828E-04_r8,  0.5629E-04_r8,  0.1117E-03_r8,   &
         0.2207E-03_r8,  0.4326E-03_r8,  0.8367E-03_r8,  0.1586E-02_r8,  0.2912E-02_r8,   &
         0.5107E-02_r8,  0.8435E-02_r8,  0.1303E-01_r8,  0.1897E-01_r8,  0.2660E-01_r8,   &
         0.3660E-01_r8,  0.4995E-01_r8,  0.6777E-01_r8,  0.9114E-01_r8,  0.1207E+00_r8,   &
         0.1566E+00_r8, &
                                ! DATA ((o2(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip= 2_r8, 2)_r8, &                                &
         0.6193E-10_r8,  0.5262E-10_r8, -0.1774E-09_r8, -0.1478E-08_r8, -0.7416E-08_r8,   &
         -0.2985E-07_r8, -0.1071E-06_r8, -0.3584E-06_r8, -0.1124E-05_r8, -0.3287E-05_r8,  &
         -0.8753E-05_r8, -0.2069E-04_r8, -0.4256E-04_r8, -0.7686E-04_r8, -0.1264E-03_r8,  &
         -0.1956E-03_r8, -0.2893E-03_r8, -0.4103E-03_r8, -0.5580E-03_r8, -0.7285E-03_r8,  &
         -0.9171E-03_r8, &
                                ! DATA ((o3(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip= 2_r8, 2)_r8, &                                &
         -0.2436E-10_r8, -0.4822E-10_r8, -0.9466E-10_r8, -0.1827E-09_r8, -0.3404E-09_r8,  &
         -0.6220E-09_r8, -0.1008E-08_r8, -0.1414E-08_r8, -0.1247E-08_r8,  0.8360E-09_r8,  &
         0.6312E-08_r8,  0.1194E-07_r8,  0.2753E-08_r8, -0.4040E-07_r8, -0.1260E-06_r8,   &
         -0.2571E-06_r8, -0.4307E-06_r8, -0.6202E-06_r8, -0.7687E-06_r8, -0.8204E-06_r8,  &
         -0.7636E-06_r8, &
                                ! DATA ((o1(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip= 3_r8, 3)_r8, &                                &
         0.7117E-05_r8,  0.1419E-04_r8,  0.2828E-04_r8,  0.5628E-04_r8,  0.1117E-03_r8,   &
         0.2207E-03_r8,  0.4326E-03_r8,  0.8367E-03_r8,  0.1586E-02_r8,  0.2912E-02_r8,   &
         0.5109E-02_r8,  0.8439E-02_r8,  0.1303E-01_r8,  0.1899E-01_r8,  0.2664E-01_r8,   &
         0.3666E-01_r8,  0.5005E-01_r8,  0.6795E-01_r8,  0.9147E-01_r8,  0.1213E+00_r8,   &
         0.1576E+00_r8, &
                                ! DATA ((o2(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip= 3_r8, 3)_r8, &                                &
         0.5658E-10_r8,  0.4212E-10_r8, -0.1977E-09_r8, -0.1516E-08_r8, -0.7481E-08_r8,   &
         -0.2995E-07_r8, -0.1072E-06_r8, -0.3583E-06_r8, -0.1123E-05_r8, -0.3283E-05_r8,  &
         -0.8744E-05_r8, -0.2067E-04_r8, -0.4252E-04_r8, -0.7679E-04_r8, -0.1262E-03_r8,  &
         -0.1953E-03_r8, -0.2889E-03_r8, -0.4096E-03_r8, -0.5567E-03_r8, -0.7263E-03_r8,  &
         -0.9130E-03_r8, &
                                ! DATA ((o3(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip= 3_r8, 3)_r8, &                                &
         -0.2433E-10_r8, -0.4815E-10_r8, -0.9453E-10_r8, -0.1825E-09_r8, -0.3400E-09_r8,  &
         -0.6215E-09_r8, -0.1007E-08_r8, -0.1415E-08_r8, -0.1253E-08_r8,  0.8143E-09_r8,  &
         0.6269E-08_r8,  0.1186E-07_r8,  0.2604E-08_r8, -0.4067E-07_r8, -0.1264E-06_r8,   &
         -0.2579E-06_r8, -0.4321E-06_r8, -0.6229E-06_r8, -0.7732E-06_r8, -0.8277E-06_r8,  &
         -0.7752E-06_r8, &
                                ! DATA ((o1(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip= 4_r8, 4)_r8, &                                &
         0.7116E-05_r8,  0.1419E-04_r8,  0.2828E-04_r8,  0.5628E-04_r8,  0.1117E-03_r8,   &
         0.2207E-03_r8,  0.4326E-03_r8,  0.8368E-03_r8,  0.1586E-02_r8,  0.2913E-02_r8,   &
         0.5111E-02_r8,  0.8444E-02_r8,  0.1305E-01_r8,  0.1902E-01_r8,  0.2669E-01_r8,   &
         0.3676E-01_r8,  0.5022E-01_r8,  0.6825E-01_r8,  0.9199E-01_r8,  0.1222E+00_r8,   &
         0.1591E+00_r8, &
                                ! DATA ((o2(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip= 4_r8, 4)_r8, &                                &
         0.4814E-10_r8,  0.2552E-10_r8, -0.2298E-09_r8, -0.1576E-08_r8, -0.7579E-08_r8,   &
         -0.3009E-07_r8, -0.1074E-06_r8, -0.3581E-06_r8, -0.1122E-05_r8, -0.3278E-05_r8,  &
         -0.8729E-05_r8, -0.2063E-04_r8, -0.4245E-04_r8, -0.7667E-04_r8, -0.1260E-03_r8,  &
         -0.1950E-03_r8, -0.2883E-03_r8, -0.4086E-03_r8, -0.5549E-03_r8, -0.7229E-03_r8,  &
         -0.9071E-03_r8, &
                                ! DATA ((o3(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip= 4_r8, 4)_r8, &                                &
         -0.2428E-10_r8, -0.4805E-10_r8, -0.9433E-10_r8, -0.1821E-09_r8, -0.3394E-09_r8,  &
         -0.6206E-09_r8, -0.1008E-08_r8, -0.1416E-08_r8, -0.1261E-08_r8,  0.7860E-09_r8,  &
         0.6188E-08_r8,  0.1171E-07_r8,  0.2389E-08_r8, -0.4109E-07_r8, -0.1271E-06_r8,   &
         -0.2591E-06_r8, -0.4344E-06_r8, -0.6267E-06_r8, -0.7797E-06_r8, -0.8378E-06_r8,  &
         -0.7901E-06_r8, &
                                ! DATA ((o1(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip= 5_r8, 5)_r8, &                                &
         0.7116E-05_r8,  0.1419E-04_r8,  0.2827E-04_r8,  0.5627E-04_r8,  0.1117E-03_r8,   &
         0.2207E-03_r8,  0.4326E-03_r8,  0.8368E-03_r8,  0.1586E-02_r8,  0.2914E-02_r8,   &
         0.5114E-02_r8,  0.8454E-02_r8,  0.1307E-01_r8,  0.1906E-01_r8,  0.2677E-01_r8,   &
         0.3690E-01_r8,  0.5048E-01_r8,  0.6872E-01_r8,  0.9281E-01_r8,  0.1236E+00_r8,   &
         0.1615E+00_r8, &
                                ! DATA ((o2(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip= 5_r8, 5)_r8, &                                &
         0.3482E-10_r8, -0.6492E-12_r8, -0.2805E-09_r8, -0.1671E-08_r8, -0.7740E-08_r8,   &
         -0.3032E-07_r8, -0.1076E-06_r8, -0.3582E-06_r8, -0.1120E-05_r8, -0.3270E-05_r8,  &
         -0.8704E-05_r8, -0.2058E-04_r8, -0.4235E-04_r8, -0.7649E-04_r8, -0.1257E-03_r8,  &
         -0.1945E-03_r8, -0.2874E-03_r8, -0.4070E-03_r8, -0.5521E-03_r8, -0.7181E-03_r8,  &
         -0.8990E-03_r8, &
                                ! DATA ((o3(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip= 5_r8, 5)_r8, &                                &
         -0.2419E-10_r8, -0.4788E-10_r8, -0.9401E-10_r8, -0.1815E-09_r8, -0.3385E-09_r8,  &
         -0.6192E-09_r8, -0.1006E-08_r8, -0.1417E-08_r8, -0.1273E-08_r8,  0.7404E-09_r8,  &
         0.6068E-08_r8,  0.1148E-07_r8,  0.2021E-08_r8, -0.4165E-07_r8, -0.1281E-06_r8,   &
         -0.2609E-06_r8, -0.4375E-06_r8, -0.6323E-06_r8, -0.7887E-06_r8, -0.8508E-06_r8,  &
         -0.8067E-06_r8, &
                                ! DATA ((o1(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip= 6_r8, 6)_r8, &                                &
         0.7114E-05_r8,  0.1419E-04_r8,  0.2827E-04_r8,  0.5627E-04_r8,  0.1117E-03_r8,   &
         0.2207E-03_r8,  0.4325E-03_r8,  0.8369E-03_r8,  0.1587E-02_r8,  0.2916E-02_r8,   &
         0.5120E-02_r8,  0.8468E-02_r8,  0.1310E-01_r8,  0.1913E-01_r8,  0.2690E-01_r8,   &
         0.3714E-01_r8,  0.5090E-01_r8,  0.6944E-01_r8,  0.9407E-01_r8,  0.1258E+00_r8,   &
         0.1651E+00_r8, &
                                ! DATA ((o2(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip= 6_r8, 6)_r8, &                                &
         0.1388E-10_r8, -0.4180E-10_r8, -0.3601E-09_r8, -0.1820E-08_r8, -0.7993E-08_r8,   &
         -0.3068E-07_r8, -0.1081E-06_r8, -0.3580E-06_r8, -0.1117E-05_r8, -0.3257E-05_r8,  &
         -0.8667E-05_r8, -0.2049E-04_r8, -0.4218E-04_r8, -0.7620E-04_r8, -0.1253E-03_r8,  &
         -0.1937E-03_r8, -0.2860E-03_r8, -0.4047E-03_r8, -0.5481E-03_r8, -0.7115E-03_r8,  &
         -0.8885E-03_r8, &
                                ! DATA ((o3(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip= 6_r8, 6)_r8, &                                &
         -0.2406E-10_r8, -0.4762E-10_r8, -0.9351E-10_r8, -0.1806E-09_r8, -0.3370E-09_r8,  &
         -0.6170E-09_r8, -0.1004E-08_r8, -0.1417E-08_r8, -0.1297E-08_r8,  0.6738E-09_r8,  &
         0.5895E-08_r8,  0.1113E-07_r8,  0.1466E-08_r8, -0.4265E-07_r8, -0.1298E-06_r8,   &
         -0.2636E-06_r8, -0.4423E-06_r8, -0.6402E-06_r8, -0.8005E-06_r8, -0.8658E-06_r8,  &
         -0.8222E-06_r8, &
                                ! DATA ((o1(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip= 7_r8, 7)_r8, &                                &
         0.7113E-05_r8,  0.1418E-04_r8,  0.2826E-04_r8,  0.5625E-04_r8,  0.1117E-03_r8,   &
         0.2206E-03_r8,  0.4325E-03_r8,  0.8371E-03_r8,  0.1588E-02_r8,  0.2918E-02_r8,   &
         0.5128E-02_r8,  0.8491E-02_r8,  0.1315E-01_r8,  0.1923E-01_r8,  0.2710E-01_r8,   &
         0.3750E-01_r8,  0.5154E-01_r8,  0.7056E-01_r8,  0.9600E-01_r8,  0.1290E+00_r8,   &
         0.1703E+00_r8, &
                                ! DATA ((o2(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip= 7_r8, 7)_r8, &                                &
         -0.1889E-10_r8, -0.1062E-09_r8, -0.4847E-09_r8, -0.2053E-08_r8, -0.8389E-08_r8,  &
         -0.3140E-07_r8, -0.1089E-06_r8, -0.3577E-06_r8, -0.1112E-05_r8, -0.3236E-05_r8,  &
         -0.8607E-05_r8, -0.2035E-04_r8, -0.4192E-04_r8, -0.7576E-04_r8, -0.1245E-03_r8,  &
         -0.1925E-03_r8, -0.2840E-03_r8, -0.4013E-03_r8, -0.5427E-03_r8, -0.7029E-03_r8,  &
         -0.8756E-03_r8, &
                                ! DATA ((o3(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip= 7_r8, 7)_r8, &                                &
         -0.2385E-10_r8, -0.4722E-10_r8, -0.9273E-10_r8, -0.1791E-09_r8, -0.3348E-09_r8,  &
         -0.6121E-09_r8, -0.9974E-09_r8, -0.1422E-08_r8, -0.1326E-08_r8,  0.5603E-09_r8,  &
         0.5604E-08_r8,  0.1061E-07_r8,  0.6106E-09_r8, -0.4398E-07_r8, -0.1321E-06_r8,   &
         -0.2676E-06_r8, -0.4490E-06_r8, -0.6507E-06_r8, -0.8145E-06_r8, -0.8801E-06_r8,  &
         -0.8311E-06_r8, &
                                ! DATA ((o1(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip= 8_r8, 8)_r8, &                                &
         0.7110E-05_r8,  0.1418E-04_r8,  0.2825E-04_r8,  0.5623E-04_r8,  0.1116E-03_r8,   &
         0.2206E-03_r8,  0.4325E-03_r8,  0.8373E-03_r8,  0.1589E-02_r8,  0.2923E-02_r8,   &
         0.5141E-02_r8,  0.8526E-02_r8,  0.1324E-01_r8,  0.1940E-01_r8,  0.2741E-01_r8,   &
         0.3807E-01_r8,  0.5253E-01_r8,  0.7227E-01_r8,  0.9889E-01_r8,  0.1338E+00_r8,   &
         0.1777E+00_r8, &
                                ! DATA ((o2(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip= 8_r8, 8)_r8, &                                &
         -0.6983E-10_r8, -0.2063E-09_r8, -0.6785E-09_r8, -0.2416E-08_r8, -0.9000E-08_r8,  &
         -0.3243E-07_r8, -0.1100E-06_r8, -0.3574E-06_r8, -0.1104E-05_r8, -0.3205E-05_r8,  &
         -0.8516E-05_r8, -0.2014E-04_r8, -0.4151E-04_r8, -0.7508E-04_r8, -0.1234E-03_r8,  &
         -0.1907E-03_r8, -0.2811E-03_r8, -0.3966E-03_r8, -0.5355E-03_r8, -0.6924E-03_r8,  &
         -0.8613E-03_r8, &
                                ! DATA ((o3(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip= 8_r8, 8)_r8, &                                &
         -0.2353E-10_r8, -0.4659E-10_r8, -0.9153E-10_r8, -0.1769E-09_r8, -0.3313E-09_r8,  &
         -0.6054E-09_r8, -0.9899E-09_r8, -0.1430E-08_r8, -0.1375E-08_r8,  0.3874E-09_r8,  &
         0.5171E-08_r8,  0.9807E-08_r8, -0.7345E-09_r8, -0.4604E-07_r8, -0.1356E-06_r8,   &
         -0.2731E-06_r8, -0.4577E-06_r8, -0.6632E-06_r8, -0.8284E-06_r8, -0.8894E-06_r8,  &
         -0.8267E-06_r8, &
                                ! DATA ((o1(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip= 9_r8, 9)_r8, &                                &
         0.7105E-05_r8,  0.1417E-04_r8,  0.2823E-04_r8,  0.5620E-04_r8,  0.1116E-03_r8,   &
         0.2205E-03_r8,  0.4325E-03_r8,  0.8376E-03_r8,  0.1591E-02_r8,  0.2929E-02_r8,   &
         0.5162E-02_r8,  0.8581E-02_r8,  0.1336E-01_r8,  0.1966E-01_r8,  0.2790E-01_r8,   &
         0.3894E-01_r8,  0.5404E-01_r8,  0.7484E-01_r8,  0.1031E+00_r8,  0.1405E+00_r8,   &
         0.1880E+00_r8, &
                                ! DATA ((o2(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip= 9_r8, 9)_r8, &                                &
         -0.1481E-09_r8, -0.3601E-09_r8, -0.9762E-09_r8, -0.2973E-08_r8, -0.1014E-07_r8,  &
         -0.3421E-07_r8, -0.1121E-06_r8, -0.3569E-06_r8, -0.1092E-05_r8, -0.3156E-05_r8,  &
         -0.8375E-05_r8, -0.1981E-04_r8, -0.4090E-04_r8, -0.7405E-04_r8, -0.1218E-03_r8,  &
         -0.1881E-03_r8, -0.2770E-03_r8, -0.3906E-03_r8, -0.5269E-03_r8, -0.6810E-03_r8,  &
         -0.8471E-03_r8, &
                                ! DATA ((o3(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip= 9_r8, 9)_r8, &                                &
         -0.2304E-10_r8, -0.4564E-10_r8, -0.8969E-10_r8, -0.1735E-09_r8, -0.3224E-09_r8,  &
         -0.5933E-09_r8, -0.9756E-09_r8, -0.1428E-08_r8, -0.1446E-08_r8,  0.1156E-09_r8,  &
         0.4499E-08_r8,  0.8469E-08_r8, -0.2720E-08_r8, -0.4904E-07_r8, -0.1401E-06_r8,   &
         -0.2801E-06_r8, -0.4681E-06_r8, -0.6761E-06_r8, -0.8387E-06_r8, -0.8879E-06_r8,  &
         -0.8040E-06_r8, &
                                ! DATA ((o1(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=10_r8,10)_r8, &                                &
         0.7098E-05_r8,  0.1415E-04_r8,  0.2821E-04_r8,  0.5615E-04_r8,  0.1115E-03_r8,   &
         0.2204E-03_r8,  0.4325E-03_r8,  0.8382E-03_r8,  0.1593E-02_r8,  0.2940E-02_r8,   &
         0.5194E-02_r8,  0.8666E-02_r8,  0.1356E-01_r8,  0.2006E-01_r8,  0.2865E-01_r8,   &
         0.4026E-01_r8,  0.5631E-01_r8,  0.7863E-01_r8,  0.1093E+00_r8,  0.1500E+00_r8,   &
         0.2017E+00_r8, &
                                ! DATA ((o2(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=10_r8,10)_r8, &                                &
         -0.2661E-09_r8, -0.5923E-09_r8, -0.1426E-08_r8, -0.3816E-08_r8, -0.1159E-07_r8,  &
         -0.3654E-07_r8, -0.1143E-06_r8, -0.3559E-06_r8, -0.1074E-05_r8, -0.3083E-05_r8,  &
         -0.8159E-05_r8, -0.1932E-04_r8, -0.3998E-04_r8, -0.7253E-04_r8, -0.1194E-03_r8,  &
         -0.1845E-03_r8, -0.2718E-03_r8, -0.3833E-03_r8, -0.5176E-03_r8, -0.6701E-03_r8,  &
         -0.8354E-03_r8, &
                                ! DATA ((o3(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=10_r8,10)_r8, &                                &
         -0.2232E-10_r8, -0.4421E-10_r8, -0.8695E-10_r8, -0.1684E-09_r8, -0.3141E-09_r8,  &
         -0.5765E-09_r8, -0.9606E-09_r8, -0.1434E-08_r8, -0.1551E-08_r8, -0.2663E-09_r8,  &
         0.3515E-08_r8,  0.6549E-08_r8, -0.5479E-08_r8, -0.5312E-07_r8, -0.1460E-06_r8,   &
         -0.2883E-06_r8, -0.4787E-06_r8, -0.6863E-06_r8, -0.8399E-06_r8, -0.8703E-06_r8,  &
         -0.7602E-06_r8, &
                                ! DATA ((o1(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=11_r8,11)_r8, &                                &
         0.7088E-05_r8,  0.1413E-04_r8,  0.2817E-04_r8,  0.5608E-04_r8,  0.1114E-03_r8,   &
         0.2203E-03_r8,  0.4325E-03_r8,  0.8390E-03_r8,  0.1598E-02_r8,  0.2955E-02_r8,   &
         0.5242E-02_r8,  0.8796E-02_r8,  0.1386E-01_r8,  0.2067E-01_r8,  0.2978E-01_r8,   &
         0.4224E-01_r8,  0.5964E-01_r8,  0.8406E-01_r8,  0.1178E+00_r8,  0.1627E+00_r8,   &
         0.2197E+00_r8, &
                                ! DATA ((o2(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=11_r8,11)_r8, &                                &
         -0.4394E-09_r8, -0.9330E-09_r8, -0.2086E-08_r8, -0.5054E-08_r8, -0.1373E-07_r8,  &
         -0.3971E-07_r8, -0.1178E-06_r8, -0.3546E-06_r8, -0.1049E-05_r8, -0.2976E-05_r8,  &
         -0.7847E-05_r8, -0.1860E-04_r8, -0.3864E-04_r8, -0.7038E-04_r8, -0.1162E-03_r8,  &
         -0.1798E-03_r8, -0.2654E-03_r8, -0.3754E-03_r8, -0.5091E-03_r8, -0.6621E-03_r8,  &
         -0.8286E-03_r8, &
                                ! DATA ((o3(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=11_r8,11)_r8, &                                &
         -0.2127E-10_r8, -0.4216E-10_r8, -0.8300E-10_r8, -0.1611E-09_r8, -0.3019E-09_r8,  &
         -0.5597E-09_r8, -0.9431E-09_r8, -0.1450E-08_r8, -0.1694E-08_r8, -0.7913E-09_r8,  &
         0.2144E-08_r8,  0.3990E-08_r8, -0.9282E-08_r8, -0.5810E-07_r8, -0.1525E-06_r8,   &
         -0.2965E-06_r8, -0.4869E-06_r8, -0.6894E-06_r8, -0.8281E-06_r8, -0.8350E-06_r8,  &
         -0.6956E-06_r8, &
                                ! DATA ((o1(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=12_r8,12)_r8, &                                &
         0.7073E-05_r8,  0.1411E-04_r8,  0.2811E-04_r8,  0.5598E-04_r8,  0.1112E-03_r8,   &
         0.2201E-03_r8,  0.4324E-03_r8,  0.8401E-03_r8,  0.1604E-02_r8,  0.2978E-02_r8,   &
         0.5313E-02_r8,  0.8987E-02_r8,  0.1431E-01_r8,  0.2158E-01_r8,  0.3145E-01_r8,   &
         0.4512E-01_r8,  0.6440E-01_r8,  0.9161E-01_r8,  0.1293E+00_r8,  0.1793E+00_r8,   &
         0.2426E+00_r8, &
                                ! DATA ((o2(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=12_r8,12)_r8, &                                &
         -0.6829E-09_r8, -0.1412E-08_r8, -0.3014E-08_r8, -0.6799E-08_r8, -0.1675E-07_r8,  &
         -0.4450E-07_r8, -0.1235E-06_r8, -0.3538E-06_r8, -0.1014E-05_r8, -0.2827E-05_r8,  &
         -0.7407E-05_r8, -0.1759E-04_r8, -0.3676E-04_r8, -0.6744E-04_r8, -0.1120E-03_r8,  &
         -0.1742E-03_r8, -0.2585E-03_r8, -0.3683E-03_r8, -0.5034E-03_r8, -0.6594E-03_r8,  &
         -0.8290E-03_r8, &
                                ! DATA ((o3(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=12_r8,12)_r8, &                                &
         -0.1985E-10_r8, -0.3937E-10_r8, -0.7761E-10_r8, -0.1511E-09_r8, -0.2855E-09_r8,  &
         -0.5313E-09_r8, -0.9251E-09_r8, -0.1470E-08_r8, -0.1898E-08_r8, -0.1519E-08_r8,  &
         0.2914E-09_r8,  0.5675E-09_r8, -0.1405E-07_r8, -0.6359E-07_r8, -0.1584E-06_r8,   &
         -0.3020E-06_r8, -0.4893E-06_r8, -0.6821E-06_r8, -0.8021E-06_r8, -0.7834E-06_r8,  &
         -0.6105E-06_r8, &
                                ! DATA ((o1(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=13_r8,13)_r8, &                                &
         0.7053E-05_r8,  0.1407E-04_r8,  0.2804E-04_r8,  0.5584E-04_r8,  0.1110E-03_r8,   &
         0.2198E-03_r8,  0.4324E-03_r8,  0.8420E-03_r8,  0.1613E-02_r8,  0.3011E-02_r8,   &
         0.5416E-02_r8,  0.9263E-02_r8,  0.1495E-01_r8,  0.2289E-01_r8,  0.3384E-01_r8,   &
         0.4918E-01_r8,  0.7096E-01_r8,  0.1018E+00_r8,  0.1442E+00_r8,  0.2004E+00_r8,   &
         0.2708E+00_r8, &
                                ! DATA ((o2(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=13_r8,13)_r8, &                                &
         -0.1004E-08_r8, -0.2043E-08_r8, -0.4239E-08_r8, -0.9104E-08_r8, -0.2075E-07_r8,  &
         -0.5096E-07_r8, -0.1307E-06_r8, -0.3520E-06_r8, -0.9671E-06_r8, -0.2630E-05_r8,  &
         -0.6825E-05_r8, -0.1624E-04_r8, -0.3429E-04_r8, -0.6369E-04_r8, -0.1069E-03_r8,  &
         -0.1680E-03_r8, -0.2520E-03_r8, -0.3635E-03_r8, -0.5029E-03_r8, -0.6647E-03_r8,  &
         -0.8390E-03_r8, &
                                ! DATA ((o3(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=13_r8,13)_r8, &                                &
         -0.1807E-10_r8, -0.3587E-10_r8, -0.7085E-10_r8, -0.1385E-09_r8, -0.2648E-09_r8,  &
         -0.4958E-09_r8, -0.8900E-09_r8, -0.1473E-08_r8, -0.2112E-08_r8, -0.2399E-08_r8,  &
         -0.2002E-08_r8, -0.3646E-08_r8, -0.1931E-07_r8, -0.6852E-07_r8, -0.1618E-06_r8,  &
         -0.3021E-06_r8, -0.4828E-06_r8, -0.6634E-06_r8, -0.7643E-06_r8, -0.7177E-06_r8,  &
         -0.5054E-06_r8, &
                                ! DATA ((o1(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=14_r8,14)_r8, &                                &
         0.7029E-05_r8,  0.1402E-04_r8,  0.2795E-04_r8,  0.5567E-04_r8,  0.1107E-03_r8,   &
         0.2195E-03_r8,  0.4326E-03_r8,  0.8447E-03_r8,  0.1625E-02_r8,  0.3056E-02_r8,   &
         0.5554E-02_r8,  0.9638E-02_r8,  0.1584E-01_r8,  0.2470E-01_r8,  0.3713E-01_r8,   &
         0.5470E-01_r8,  0.7969E-01_r8,  0.1149E+00_r8,  0.1631E+00_r8,  0.2263E+00_r8,   &
         0.3045E+00_r8, &
                                ! DATA ((o2(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=14_r8,14)_r8, &                                &
         -0.1387E-08_r8, -0.2798E-08_r8, -0.5706E-08_r8, -0.1187E-07_r8, -0.2564E-07_r8,  &
         -0.5866E-07_r8, -0.1398E-06_r8, -0.3516E-06_r8, -0.9148E-06_r8, -0.2398E-05_r8,  &
         -0.6122E-05_r8, -0.1459E-04_r8, -0.3125E-04_r8, -0.5923E-04_r8, -0.1013E-03_r8,  &
         -0.1620E-03_r8, -0.2473E-03_r8, -0.3631E-03_r8, -0.5098E-03_r8, -0.6800E-03_r8,  &
         -0.8603E-03_r8, &
                                ! DATA ((o3(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=14_r8,14)_r8, &                                &
         -0.1610E-10_r8, -0.3200E-10_r8, -0.6337E-10_r8, -0.1245E-09_r8, -0.2408E-09_r8,  &
         -0.4533E-09_r8, -0.8405E-09_r8, -0.1464E-08_r8, -0.2337E-08_r8, -0.3341E-08_r8,  &
         -0.4467E-08_r8, -0.8154E-08_r8, -0.2436E-07_r8, -0.7128E-07_r8, -0.1604E-06_r8,  &
         -0.2945E-06_r8, -0.4666E-06_r8, -0.6357E-06_r8, -0.7187E-06_r8, -0.6419E-06_r8,  &
         -0.3795E-06_r8, &
                                ! DATA ((o1(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=15_r8,15)_r8, &                                &
         0.7003E-05_r8,  0.1397E-04_r8,  0.2785E-04_r8,  0.5550E-04_r8,  0.1104E-03_r8,   &
         0.2192E-03_r8,  0.4329E-03_r8,  0.8481E-03_r8,  0.1641E-02_r8,  0.3112E-02_r8,   &
         0.5729E-02_r8,  0.1012E-01_r8,  0.1698E-01_r8,  0.2708E-01_r8,  0.4146E-01_r8,   &
         0.6189E-01_r8,  0.9085E-01_r8,  0.1313E+00_r8,  0.1862E+00_r8,  0.2571E+00_r8,   &
         0.3433E+00_r8, &
                                ! DATA ((o2(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=15_r8,15)_r8, &                                &
         -0.1788E-08_r8, -0.3588E-08_r8, -0.7244E-08_r8, -0.1479E-07_r8, -0.3083E-07_r8,  &
         -0.6671E-07_r8, -0.1497E-06_r8, -0.3519E-06_r8, -0.8607E-06_r8, -0.2154E-05_r8,  &
         -0.5364E-05_r8, -0.1276E-04_r8, -0.2785E-04_r8, -0.5435E-04_r8, -0.9573E-04_r8,  &
         -0.1570E-03_r8, -0.2455E-03_r8, -0.3682E-03_r8, -0.5253E-03_r8, -0.7065E-03_r8,  &
         -0.8938E-03_r8, &
                                ! DATA ((o3(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=15_r8,15)_r8, &                                &
         -0.1429E-10_r8, -0.2843E-10_r8, -0.5645E-10_r8, -0.1115E-09_r8, -0.2181E-09_r8,  &
         -0.4200E-09_r8, -0.7916E-09_r8, -0.1460E-08_r8, -0.2542E-08_r8, -0.4168E-08_r8,  &
         -0.6703E-08_r8, -0.1215E-07_r8, -0.2821E-07_r8, -0.7073E-07_r8, -0.1530E-06_r8,  &
         -0.2791E-06_r8, -0.4426E-06_r8, -0.6027E-06_r8, -0.6707E-06_r8, -0.5591E-06_r8,  &
         -0.2328E-06_r8, &
                                ! DATA ((o1(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=16_r8,16)_r8, &                                &
         0.6981E-05_r8,  0.1393E-04_r8,  0.2777E-04_r8,  0.5535E-04_r8,  0.1102E-03_r8,   &
         0.2190E-03_r8,  0.4336E-03_r8,  0.8527E-03_r8,  0.1660E-02_r8,  0.3177E-02_r8,   &
         0.5930E-02_r8,  0.1068E-01_r8,  0.1836E-01_r8,  0.3000E-01_r8,  0.4686E-01_r8,   &
         0.7083E-01_r8,  0.1046E+00_r8,  0.1511E+00_r8,  0.2134E+00_r8,  0.2924E+00_r8,   &
         0.3861E+00_r8, &
                                ! DATA ((o2(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=16_r8,16)_r8, &                                &
         -0.2141E-08_r8, -0.4286E-08_r8, -0.8603E-08_r8, -0.1737E-07_r8, -0.3548E-07_r8,  &
         -0.7410E-07_r8, -0.1590E-06_r8, -0.3537E-06_r8, -0.8142E-06_r8, -0.1935E-05_r8,  &
         -0.4658E-05_r8, -0.1099E-04_r8, -0.2444E-04_r8, -0.4948E-04_r8, -0.9067E-04_r8,  &
         -0.1538E-03_r8, -0.2474E-03_r8, -0.3793E-03_r8, -0.5495E-03_r8, -0.7439E-03_r8,  &
         -0.9383E-03_r8, &
                                ! DATA ((o3(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=16_r8,16)_r8, &                                &
         -0.1295E-10_r8, -0.2581E-10_r8, -0.5136E-10_r8, -0.1019E-09_r8, -0.2011E-09_r8,  &
         -0.3916E-09_r8, -0.7585E-09_r8, -0.1439E-08_r8, -0.2648E-08_r8, -0.4747E-08_r8,  &
         -0.8301E-08_r8, -0.1499E-07_r8, -0.3024E-07_r8, -0.6702E-07_r8, -0.1399E-06_r8,  &
         -0.2564E-06_r8, -0.4117E-06_r8, -0.5669E-06_r8, -0.6239E-06_r8, -0.4748E-06_r8,  &
         -0.7013E-07_r8, &
                                ! DATA ((o1(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=17_r8,17)_r8, &                                &
         0.6967E-05_r8,  0.1390E-04_r8,  0.2772E-04_r8,  0.5527E-04_r8,  0.1101E-03_r8,   &
         0.2191E-03_r8,  0.4346E-03_r8,  0.8579E-03_r8,  0.1679E-02_r8,  0.3244E-02_r8,   &
         0.6139E-02_r8,  0.1127E-01_r8,  0.1986E-01_r8,  0.3330E-01_r8,  0.5315E-01_r8,   &
         0.8139E-01_r8,  0.1207E+00_r8,  0.1741E+00_r8,  0.2442E+00_r8,  0.3311E+00_r8,   &
         0.4312E+00_r8, &
                                ! DATA ((o2(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=17_r8,17)_r8, &                                &
         -0.2400E-08_r8, -0.4796E-08_r8, -0.9599E-08_r8, -0.1927E-07_r8, -0.3892E-07_r8,  &
         -0.7954E-07_r8, -0.1661E-06_r8, -0.3540E-06_r8, -0.7780E-06_r8, -0.1763E-05_r8,  &
         -0.4092E-05_r8, -0.9512E-05_r8, -0.2142E-04_r8, -0.4502E-04_r8, -0.8640E-04_r8,  &
         -0.1525E-03_r8, -0.2526E-03_r8, -0.3955E-03_r8, -0.5805E-03_r8, -0.7897E-03_r8,  &
         -0.9899E-03_r8, &
                                ! DATA ((o3(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=17_r8,17)_r8, &                                &
         -0.1220E-10_r8, -0.2432E-10_r8, -0.4845E-10_r8, -0.9640E-10_r8, -0.1912E-09_r8,  &
         -0.3771E-09_r8, -0.7392E-09_r8, -0.1420E-08_r8, -0.2702E-08_r8, -0.5049E-08_r8,  &
         -0.9214E-08_r8, -0.1659E-07_r8, -0.3101E-07_r8, -0.6162E-07_r8, -0.1235E-06_r8,  &
         -0.2287E-06_r8, -0.3755E-06_r8, -0.5274E-06_r8, -0.5790E-06_r8, -0.3947E-06_r8,  &
         0.1003E-06_r8, &
                                ! DATA ((o1(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=18_r8,18)_r8, &                                &
         0.6963E-05_r8,  0.1389E-04_r8,  0.2771E-04_r8,  0.5526E-04_r8,  0.1101E-03_r8,   &
         0.2193E-03_r8,  0.4359E-03_r8,  0.8629E-03_r8,  0.1698E-02_r8,  0.3305E-02_r8,   &
         0.6331E-02_r8,  0.1183E-01_r8,  0.2133E-01_r8,  0.3671E-01_r8,  0.5996E-01_r8,   &
         0.9320E-01_r8,  0.1389E+00_r8,  0.1999E+00_r8,  0.2778E+00_r8,  0.3717E+00_r8,   &
         0.4761E+00_r8, &
                                ! DATA ((o2(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=18_r8,18)_r8, &                                &
         -0.2557E-08_r8, -0.5106E-08_r8, -0.1020E-07_r8, -0.2043E-07_r8, -0.4103E-07_r8,  &
         -0.8293E-07_r8, -0.1697E-06_r8, -0.3531E-06_r8, -0.7531E-06_r8, -0.1645E-05_r8,  &
         -0.3690E-05_r8, -0.8411E-05_r8, -0.1902E-04_r8, -0.4118E-04_r8, -0.8276E-04_r8,  &
         -0.1525E-03_r8, -0.2601E-03_r8, -0.4147E-03_r8, -0.6149E-03_r8, -0.8384E-03_r8,  &
         -0.1042E-02_r8, &
                                ! DATA ((o3(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=18_r8,18)_r8, &                                &
         -0.1189E-10_r8, -0.2372E-10_r8, -0.4729E-10_r8, -0.9421E-10_r8, -0.1873E-09_r8,  &
         -0.3713E-09_r8, -0.7317E-09_r8, -0.1437E-08_r8, -0.2764E-08_r8, -0.5243E-08_r8,  &
         -0.9691E-08_r8, -0.1751E-07_r8, -0.3122E-07_r8, -0.5693E-07_r8, -0.1076E-06_r8,  &
         -0.1981E-06_r8, -0.3324E-06_r8, -0.4785E-06_r8, -0.5280E-06_r8, -0.3174E-06_r8,  &
         0.2672E-06_r8, &
                                ! DATA ((o1(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=19_r8,19)_r8, &                                &
         0.6963E-05_r8,  0.1389E-04_r8,  0.2771E-04_r8,  0.5528E-04_r8,  0.1102E-03_r8,   &
         0.2196E-03_r8,  0.4370E-03_r8,  0.8672E-03_r8,  0.1712E-02_r8,  0.3355E-02_r8,   &
         0.6488E-02_r8,  0.1230E-01_r8,  0.2262E-01_r8,  0.3989E-01_r8,  0.6677E-01_r8,   &
         0.1056E+00_r8,  0.1586E+00_r8,  0.2276E+00_r8,  0.3131E+00_r8,  0.4124E+00_r8,   &
         0.5188E+00_r8, &
                                ! DATA ((o2(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=19_r8,19)_r8, &                                &
         -0.2630E-08_r8, -0.5249E-08_r8, -0.1048E-07_r8, -0.2096E-07_r8, -0.4198E-07_r8,  &
         -0.8440E-07_r8, -0.1710E-06_r8, -0.3513E-06_r8, -0.7326E-06_r8, -0.1562E-05_r8,  &
         -0.3416E-05_r8, -0.7637E-05_r8, -0.1719E-04_r8, -0.3795E-04_r8, -0.7926E-04_r8,  &
         -0.1524E-03_r8, -0.2680E-03_r8, -0.4344E-03_r8, -0.6486E-03_r8, -0.8838E-03_r8,  &
         -0.1089E-02_r8, &
                                ! DATA ((o3(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=19_r8,19)_r8, &                                &
         -0.1188E-10_r8, -0.2369E-10_r8, -0.4725E-10_r8, -0.9417E-10_r8, -0.1875E-09_r8,  &
         -0.3725E-09_r8, -0.7365E-09_r8, -0.1445E-08_r8, -0.2814E-08_r8, -0.5384E-08_r8,  &
         -0.1008E-07_r8, -0.1816E-07_r8, -0.3179E-07_r8, -0.5453E-07_r8, -0.9500E-07_r8,  &
         -0.1679E-06_r8, -0.2819E-06_r8, -0.4109E-06_r8, -0.4555E-06_r8, -0.2283E-06_r8,  &
         0.4283E-06_r8, &
                                ! DATA ((o1(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=20_r8,20)_r8, &                                &
         0.6963E-05_r8,  0.1389E-04_r8,  0.2772E-04_r8,  0.5529E-04_r8,  0.1103E-03_r8,   &
         0.2198E-03_r8,  0.4377E-03_r8,  0.8701E-03_r8,  0.1723E-02_r8,  0.3391E-02_r8,   &
         0.6606E-02_r8,  0.1266E-01_r8,  0.2366E-01_r8,  0.4262E-01_r8,  0.7304E-01_r8,   &
         0.1179E+00_r8,  0.1789E+00_r8,  0.2563E+00_r8,  0.3487E+00_r8,  0.4516E+00_r8,   &
         0.5572E+00_r8, &
                                ! DATA ((o2(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=20_r8,20)_r8, &                                &
         -0.2651E-08_r8, -0.5291E-08_r8, -0.1056E-07_r8, -0.2110E-07_r8, -0.4221E-07_r8,  &
         -0.8462E-07_r8, -0.1705E-06_r8, -0.3466E-06_r8, -0.7155E-06_r8, -0.1501E-05_r8,  &
         -0.3223E-05_r8, -0.7079E-05_r8, -0.1581E-04_r8, -0.3517E-04_r8, -0.7553E-04_r8,  &
         -0.1510E-03_r8, -0.2746E-03_r8, -0.4528E-03_r8, -0.6789E-03_r8, -0.9214E-03_r8,  &
         -0.1124E-02_r8, &
                                ! DATA ((o3(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=20_r8,20)_r8, &                                &
         -0.1193E-10_r8, -0.2380E-10_r8, -0.4748E-10_r8, -0.9465E-10_r8, -0.1886E-09_r8,  &
         -0.3751E-09_r8, -0.7436E-09_r8, -0.1466E-08_r8, -0.2872E-08_r8, -0.5508E-08_r8,  &
         -0.1038E-07_r8, -0.1891E-07_r8, -0.3279E-07_r8, -0.5420E-07_r8, -0.8711E-07_r8,  &
         -0.1403E-06_r8, -0.2248E-06_r8, -0.3221E-06_r8, -0.3459E-06_r8, -0.1066E-06_r8,  &
         0.5938E-06_r8, &
                                ! DATA ((o1(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=21_r8,21)_r8, &                                &
         0.6964E-05_r8,  0.1389E-04_r8,  0.2772E-04_r8,  0.5530E-04_r8,  0.1103E-03_r8,   &
         0.2199E-03_r8,  0.4382E-03_r8,  0.8719E-03_r8,  0.1730E-02_r8,  0.3416E-02_r8,   &
         0.6690E-02_r8,  0.1293E-01_r8,  0.2445E-01_r8,  0.4479E-01_r8,  0.7837E-01_r8,   &
         0.1291E+00_r8,  0.1985E+00_r8,  0.2846E+00_r8,  0.3831E+00_r8,  0.4875E+00_r8,   &
         0.5902E+00_r8, &
                                ! DATA ((o2(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=21_r8,21)_r8, &                                &
         -0.2654E-08_r8, -0.5296E-08_r8, -0.1057E-07_r8, -0.2111E-07_r8, -0.4219E-07_r8,  &
         -0.8445E-07_r8, -0.1696E-06_r8, -0.3428E-06_r8, -0.7013E-06_r8, -0.1458E-05_r8,  &
         -0.3084E-05_r8, -0.6678E-05_r8, -0.1476E-04_r8, -0.3284E-04_r8, -0.7173E-04_r8,  &
         -0.1481E-03_r8, -0.2786E-03_r8, -0.4688E-03_r8, -0.7052E-03_r8, -0.9506E-03_r8,  &
         -0.1148E-02_r8, &
                                ! DATA ((o3(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=21_r8,21)_r8, &                                &
         -0.1195E-10_r8, -0.2384E-10_r8, -0.4755E-10_r8, -0.9482E-10_r8, -0.1890E-09_r8,  &
         -0.3761E-09_r8, -0.7469E-09_r8, -0.1476E-08_r8, -0.2892E-08_r8, -0.5603E-08_r8,  &
         -0.1060E-07_r8, -0.1942E-07_r8, -0.3393E-07_r8, -0.5508E-07_r8, -0.8290E-07_r8,  &
         -0.1182E-06_r8, -0.1657E-06_r8, -0.2170E-06_r8, -0.1997E-06_r8,  0.6227E-07_r8,  &
         0.7847E-06_r8, &
                                ! DATA ((o1(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=22_r8,22)_r8, &                                &
         0.6964E-05_r8,  0.1389E-04_r8,  0.2772E-04_r8,  0.5530E-04_r8,  0.1103E-03_r8,   &
         0.2200E-03_r8,  0.4385E-03_r8,  0.8731E-03_r8,  0.1735E-02_r8,  0.3433E-02_r8,   &
         0.6748E-02_r8,  0.1311E-01_r8,  0.2502E-01_r8,  0.4642E-01_r8,  0.8258E-01_r8,   &
         0.1385E+00_r8,  0.2160E+00_r8,  0.3107E+00_r8,  0.4146E+00_r8,  0.5188E+00_r8,   &
         0.6171E+00_r8, &
                                ! DATA ((o2(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=22_r8,22)_r8, &                                &
         -0.2653E-08_r8, -0.5295E-08_r8, -0.1057E-07_r8, -0.2110E-07_r8, -0.4215E-07_r8,  &
         -0.8430E-07_r8, -0.1690E-06_r8, -0.3403E-06_r8, -0.6919E-06_r8, -0.1427E-05_r8,  &
         -0.2991E-05_r8, -0.6399E-05_r8, -0.1398E-04_r8, -0.3099E-04_r8, -0.6824E-04_r8,  &
         -0.1441E-03_r8, -0.2795E-03_r8, -0.4814E-03_r8, -0.7282E-03_r8, -0.9739E-03_r8,  &
         -0.1163E-02_r8, &
                                ! DATA ((o3(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=22_r8,22)_r8, &                                &
         -0.1195E-10_r8, -0.2384E-10_r8, -0.4756E-10_r8, -0.9485E-10_r8, -0.1891E-09_r8,  &
         -0.3765E-09_r8, -0.7483E-09_r8, -0.1481E-08_r8, -0.2908E-08_r8, -0.5660E-08_r8,  &
         -0.1075E-07_r8, -0.1980E-07_r8, -0.3472E-07_r8, -0.5626E-07_r8, -0.8149E-07_r8,  &
         -0.1027E-06_r8, -0.1136E-06_r8, -0.1071E-06_r8, -0.2991E-07_r8,  0.2743E-06_r8,  &
         0.1017E-05_r8, &
                                ! DATA ((o1(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=23_r8,23)_r8, &                                &
         0.6964E-05_r8,  0.1389E-04_r8,  0.2772E-04_r8,  0.5530E-04_r8,  0.1103E-03_r8,   &
         0.2200E-03_r8,  0.4387E-03_r8,  0.8738E-03_r8,  0.1738E-02_r8,  0.3445E-02_r8,   &
         0.6786E-02_r8,  0.1324E-01_r8,  0.2541E-01_r8,  0.4757E-01_r8,  0.8567E-01_r8,   &
         0.1459E+00_r8,  0.2303E+00_r8,  0.3331E+00_r8,  0.4415E+00_r8,  0.5444E+00_r8,   &
         0.6377E+00_r8, &
                                ! DATA ((o2(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=23_r8,23)_r8, &                                &
         -0.2653E-08_r8, -0.5294E-08_r8, -0.1057E-07_r8, -0.2109E-07_r8, -0.4212E-07_r8,  &
         -0.8420E-07_r8, -0.1686E-06_r8, -0.3388E-06_r8, -0.6858E-06_r8, -0.1406E-05_r8,  &
         -0.2928E-05_r8, -0.6206E-05_r8, -0.1344E-04_r8, -0.2961E-04_r8, -0.6533E-04_r8,  &
         -0.1399E-03_r8, -0.2780E-03_r8, -0.4904E-03_r8, -0.7488E-03_r8, -0.9953E-03_r8,  &
         -0.1175E-02_r8, &
                                ! DATA ((o3(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=23_r8,23)_r8, &                                &
         -0.1195E-10_r8, -0.2384E-10_r8, -0.4756E-10_r8, -0.9485E-10_r8, -0.1891E-09_r8,  &
         -0.3767E-09_r8, -0.7492E-09_r8, -0.1485E-08_r8, -0.2924E-08_r8, -0.5671E-08_r8,  &
         -0.1084E-07_r8, -0.2009E-07_r8, -0.3549E-07_r8, -0.5773E-07_r8, -0.8208E-07_r8,  &
         -0.9394E-07_r8, -0.7270E-07_r8, -0.3947E-08_r8,  0.1456E-06_r8,  0.5083E-06_r8,  &
         0.1270E-05_r8, &
                                ! DATA ((o1(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=24_r8,24)_r8, &                                &
         0.6964E-05_r8,  0.1389E-04_r8,  0.2772E-04_r8,  0.5531E-04_r8,  0.1103E-03_r8,   &
         0.2201E-03_r8,  0.4388E-03_r8,  0.8743E-03_r8,  0.1740E-02_r8,  0.3452E-02_r8,   &
         0.6811E-02_r8,  0.1332E-01_r8,  0.2567E-01_r8,  0.4835E-01_r8,  0.8782E-01_r8,   &
         0.1511E+00_r8,  0.2412E+00_r8,  0.3507E+00_r8,  0.4627E+00_r8,  0.5638E+00_r8,   &
         0.6526E+00_r8, &
                                ! DATA ((o2(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=24_r8,24)_r8, &                                &
         -0.2653E-08_r8, -0.5294E-08_r8, -0.1056E-07_r8, -0.2109E-07_r8, -0.4210E-07_r8,  &
         -0.8413E-07_r8, -0.1684E-06_r8, -0.3379E-06_r8, -0.6820E-06_r8, -0.1393E-05_r8,  &
         -0.2889E-05_r8, -0.6080E-05_r8, -0.1307E-04_r8, -0.2861E-04_r8, -0.6310E-04_r8,  &
         -0.1363E-03_r8, -0.2758E-03_r8, -0.4969E-03_r8, -0.7681E-03_r8, -0.1017E-02_r8,  &
         -0.1186E-02_r8, &
                                ! DATA ((o3(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=24_r8,24)_r8, &                                &
         -0.1195E-10_r8, -0.2384E-10_r8, -0.4756E-10_r8, -0.9485E-10_r8, -0.1891E-09_r8,  &
         -0.3768E-09_r8, -0.7497E-09_r8, -0.1487E-08_r8, -0.2933E-08_r8, -0.5710E-08_r8,  &
         -0.1089E-07_r8, -0.2037E-07_r8, -0.3616E-07_r8, -0.5907E-07_r8, -0.8351E-07_r8,  &
         -0.8925E-07_r8, -0.4122E-07_r8,  0.8779E-07_r8,  0.3143E-06_r8,  0.7281E-06_r8,  &
         0.1500E-05_r8, &
                                ! DATA ((o1(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=25_r8,25)_r8, &                                &
         0.6964E-05_r8,  0.1389E-04_r8,  0.2772E-04_r8,  0.5531E-04_r8,  0.1103E-03_r8,   &
         0.2201E-03_r8,  0.4388E-03_r8,  0.8745E-03_r8,  0.1741E-02_r8,  0.3456E-02_r8,   &
         0.6827E-02_r8,  0.1337E-01_r8,  0.2584E-01_r8,  0.4885E-01_r8,  0.8924E-01_r8,   &
         0.1547E+00_r8,  0.2488E+00_r8,  0.3632E+00_r8,  0.4779E+00_r8,  0.5771E+00_r8,   &
         0.6627E+00_r8, &
                                ! DATA ((o2(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=25_r8,25)_r8, &                                &
         -0.2653E-08_r8, -0.5293E-08_r8, -0.1056E-07_r8, -0.2108E-07_r8, -0.4209E-07_r8,  &
         -0.8409E-07_r8, -0.1682E-06_r8, -0.3373E-06_r8, -0.6797E-06_r8, -0.1383E-05_r8,  &
         -0.2862E-05_r8, -0.5993E-05_r8, -0.1283E-04_r8, -0.2795E-04_r8, -0.6158E-04_r8,  &
         -0.1338E-03_r8, -0.2743E-03_r8, -0.5030E-03_r8, -0.7863E-03_r8, -0.1038E-02_r8,  &
         -0.1196E-02_r8, &
                                ! DATA ((o3(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=25_r8,25)_r8, &                                &
         -0.1195E-10_r8, -0.2383E-10_r8, -0.4755E-10_r8, -0.9484E-10_r8, -0.1891E-09_r8,  &
         -0.3768E-09_r8, -0.7499E-09_r8, -0.1489E-08_r8, -0.2939E-08_r8, -0.5741E-08_r8,  &
         -0.1100E-07_r8, -0.2066E-07_r8, -0.3660E-07_r8, -0.6002E-07_r8, -0.8431E-07_r8,  &
         -0.8556E-07_r8, -0.1674E-07_r8,  0.1638E-06_r8,  0.4525E-06_r8,  0.8949E-06_r8,  &
         0.1669E-05_r8, &
                                ! DATA ((o1(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=26_r8,26)_r8, &                                &
         0.6964E-05_r8,  0.1389E-04_r8,  0.2772E-04_r8,  0.5531E-04_r8,  0.1103E-03_r8,   &
         0.2201E-03_r8,  0.4389E-03_r8,  0.8747E-03_r8,  0.1741E-02_r8,  0.3458E-02_r8,   &
         0.6836E-02_r8,  0.1340E-01_r8,  0.2594E-01_r8,  0.4916E-01_r8,  0.9011E-01_r8,   &
         0.1569E+00_r8,  0.2536E+00_r8,  0.3714E+00_r8,  0.4877E+00_r8,  0.5856E+00_r8,   &
         0.6695E+00_r8, &
                                ! DATA ((o2(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=26_r8,26)_r8, &                                &
         -0.2652E-08_r8, -0.5292E-08_r8, -0.1056E-07_r8, -0.2108E-07_r8, -0.4208E-07_r8,  &
         -0.8406E-07_r8, -0.1681E-06_r8, -0.3369E-06_r8, -0.6784E-06_r8, -0.1378E-05_r8,  &
         -0.2843E-05_r8, -0.5944E-05_r8, -0.1269E-04_r8, -0.2759E-04_r8, -0.6078E-04_r8,  &
         -0.1326E-03_r8, -0.2742E-03_r8, -0.5088E-03_r8, -0.8013E-03_r8, -0.1054E-02_r8,  &
         -0.1202E-02_r8, &
                                ! DATA ((o3(ip_r8,iw_r8,1)_r8,iw=1_r8,21)_r8,ip=26_r8,26)_r8, &                                &
         -0.1194E-10_r8, -0.2383E-10_r8, -0.4754E-10_r8, -0.9482E-10_r8, -0.1891E-09_r8,  &
         -0.3768E-09_r8, -0.7499E-09_r8, -0.1489E-08_r8, -0.2941E-08_r8, -0.5752E-08_r8,  &
         -0.1104E-07_r8, -0.2069E-07_r8, -0.3661E-07_r8, -0.6012E-07_r8, -0.8399E-07_r8,  &
         -0.8183E-07_r8,  0.1930E-08_r8,  0.2167E-06_r8,  0.5434E-06_r8,  0.9990E-06_r8,  &
         0.1787E-05_r8/),SHAPE=(/nx*3*no/))

    it=0
    DO k=1,nx
       DO i=1,3
          DO j=1,no
             it=it+1
             IF(i==1) THEN
                !WRITE(*,'(a5,2e17.9) ' )'o1   ',data2(it),o1(k,j)
                o1(k,j,1)=data2(it)
             END IF
             IF(i==2)THEN
                !WRITE(*,'(a5,2e17.9) ' )'o2   ',data2(it),o2(k,j)
                o2(k,j,1)=data2(it)
             END IF
             IF(i==3) THEN
                !WRITE(*,'(a5,2e17.9) ' )'o3   ',data2(it),o3(k,j,1)
                o3(k,j,1) =data2(it)
             END IF
          END DO
       END DO
    END DO

    !
    !-----------------------------------------------------------------------
    !
    !  End of RADDATA
    !
    !-----------------------------------------------------------------------
    !

  END SUBROUTINE  raddata

END MODULE Rad_CliRadLW



!PROGRAM Main
! USE Rad_CliRadLW, ONLY:InitCliRadLW
! IMPLICIT NONE
!   REAL :: gkc1(6,2)
!  DATA gkc1/  0.1395,0.1407,0.1549,0.1357,0.0182,0.0220,                 &
!             0.0766,0.1372,0.1189,0.0335,0.0169,0.0059/
!  REAL, PARAMETER :: gkc(6,2)=RESHAPE(SOURCE=(/&
!                                           0.1395,0.1407,0.1549,0.1357,0.0182,0.0220,&
!                                           0.0766,0.1372,0.1189,0.0335,0.0169,0.0059/),SHAPE=(/6,2/))
!
! CALL InitCliRadLW(.TRUE.)
!
!-----gkc is the planck-weighted co2 k-distribution function
!  in the band-wing and band-center regions given in table 7.
!  for computing efficiency, sub-bands 3a and 3c are combined.



!END PROGRAM Main
