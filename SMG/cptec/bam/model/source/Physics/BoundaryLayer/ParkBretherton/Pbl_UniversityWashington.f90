MODULE DiffCoef
  USE PhysicalFunctions,ONLY : fpvs2es5
  IMPLICIT NONE
  SAVE
  PRIVATE

  ! Selecting Kinds
  INTEGER, PARAMETER :: r4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
  INTEGER, PARAMETER :: r16 = SELECTED_REAL_KIND(15)! Kind for 128-bits Real Numbers

  REAL (KIND=r8), PARAMETER   :: grav  =                   9.8e0_r8! gravity constant               (m/s**2)
  REAL(kind=r8),PARAMETER:: con_rd     =2.8705e+2_r8      ! gas constant air    (J/kg/K)
  REAL (KIND=r8), PARAMETER   :: con_rv =    4.6150e+2_r8 ! gas constant H2O (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_cp     =1.0046e+3_r8      ! spec heat air @p    (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_hvap   =2.5000e+6_r8      ! lat heat H2O cond   (J/kg)

  REAL (KIND=r8), PARAMETER   :: con_fvirt  =con_rv/con_rd-1.0_r8
  PUBLIC :: moninq

 CONTAINS

  SUBROUTINE moninq(ix      ,&!integer      , INTENT(IN   ) :: ix
       im      ,&!integer      , INTENT(IN   ) :: im
       km      ,&!integer      , INTENT(IN   ) :: km
       uo      ,&!real(kind=r8), INTENT(IN   ) :: uo(ix,km)
       vo      ,&!real(kind=r8), INTENT(IN   ) :: vo(ix,km)
       t1      ,&!real(kind=r8), INTENT(IN   ) :: t1(ix,km)
       q1      ,&!real(kind=r8), INTENT(IN   ) :: q1(ix,km)
       qliq    ,&!real(kind=r8), INTENT(IN   ) :: qliq(ix,km)
       bps     ,&
       swh     ,&!real(kind=r8), INTENT(IN   ) :: swh(ix,km)! swh  - real, total sky sw heating rates ( k/s )   ix,levs !
       hlw     ,&!real(kind=r8), INTENT(IN   ) :: hlw(ix,km)! hlw  - real, total sky lw heating rates ( k/s )   ix,levs !
       taux    , &!real(r8), intent(in)    :: taux(pcols)       ! Zonal wind stress at surface [ N/m2 ]
       tauy    , &!real(r8), intent(in)    :: tauy(pcols)       ! Meridional wind stress at surface [ N/m2 ]
       TSKIN   ,&!real(kind=r8), INTENT(IN   ) :: TSKIN(im)
       heat    ,&!real(kind=r8), INTENT(IN   ) :: heat(im)
       evap    ,&!real(kind=r8), INTENT(IN   ) :: evap(im)
       kpbl    ,&!integer      , INTENT(OUt  ) :: kpbl(im)
       prsi    ,&!real(kind=r8), INTENT(IN   ) :: prsi(ix,km+1)
       prsl    ,&!real(kind=r8), INTENT(IN   ) :: prsl(ix,km)
       phii    ,&!real(kind=r8), INTENT(IN   ) :: phii(ix,km+1)
       phil    ,&!real(kind=r8), INTENT(IN   ) :: phil(ix,km)
       deltim  ,&!real(kind=r8), INTENT(IN   ) :: deltim
       hpbl    ,&!real(kind=r8), INTENT(OUt  ) :: hpbl(im)
       Z0      ,&!real(kind=r8), INTENT(IN   ) :: Z0RL(im)
       mlsi    ,&!real(kind=r8), INTENT(IN   ) :: mlsi(im)!       !     slmsk      - real, sea/land/ice mask (=0/1/2)                im   !
       vfrac   ,&!real(kind=r8), INTENT(IN   ) :: vfrac(IM)!      vfrac       - real, vegetation fraction                       im   !
       dkt     ,&!real(kind=r8), INTENT(OUt  ) :: dkt (im,km+1)
       dku      )!real(kind=r8), INTENT(OUt  ) :: dku (im,km+1)
    !
    !      use machine  , only : r8
    !      use funcphys , only : fpvs
    !      use physcons, grav => con_g, rd => con_rd, cp => con_cp
    !     &,             con_hvap => con_hvap, fv => con_fvirt
    IMPLICIT NONE
    !
    !     arguments
    !
    INTEGER      , INTENT(IN   ) :: ix
    INTEGER      , INTENT(IN   ) :: im
    INTEGER      , INTENT(IN   ) :: km
    REAL(kind=r8), INTENT(IN   ) :: uo(ix,km)
    REAL(kind=r8), INTENT(IN   ) :: vo(ix,km)
    REAL(kind=r8), INTENT(IN   ) :: t1(ix,km)
    REAL(kind=r8), INTENT(IN   ) :: q1(ix,km)
    REAL(kind=r8), INTENT(IN   ) :: qliq(ix,km)
    REAL(kind=r8), INTENT(IN   ) :: bps(ix,km)
    REAL(kind=r8), INTENT(IN   ) :: swh(ix,km)!     swh      - real, total sky sw heating rates ( k/s )       ix,levs !
    REAL(kind=r8), INTENT(IN   ) :: hlw(ix,km)!     hlw      - real, total sky lw heating rates ( k/s )       ix,levs !
    REAL(kind=r8), INTENT(IN   ) :: taux  (im)!real(r8), intent(in)    :: taux(pcols)               ! Zonal wind stress at surface [ N/m2 ]
    REAL(kind=r8), INTENT(IN   ) :: tauy (im)!real(r8), intent(in)    :: tauy(pcols)               ! Meridional wind stress at surface [ N/m2 ]

    REAL(kind=r8), INTENT(IN   ) :: TSKIN(im)
    REAL(kind=r8), INTENT(IN   ) :: heat(im)
    REAL(kind=r8), INTENT(IN   ) :: evap(im)
    INTEGER      , INTENT(OUt  ) :: kpbl(im)
    REAL(kind=r8), INTENT(IN   ) :: prsi(ix,km+1)
    REAL(kind=r8), INTENT(IN   ) :: prsl(ix,km)
    REAL(kind=r8), INTENT(IN   ) :: phii(ix,km+1)
    REAL(kind=r8), INTENT(IN   ) :: phil(ix,km)
    REAL(kind=r8), INTENT(IN   ) :: deltim
    REAL(kind=r8), INTENT(OUt  ) :: hpbl(im)
    REAL(kind=r8), INTENT(IN   ) :: Z0(im)
    INTEGER(KIND=i8), INTENT(IN   ) :: mlsi(im)!     !     slmsk    - real, sea/land/ice mask (=0/1/2)                  im   !
    REAL(kind=r8), INTENT(IN   ) :: vfrac(IM)!     vfrac    - real, vegetation fraction                         im   !
    REAL(kind=r8), INTENT(OUt  ) :: dkt(im,km+1)
    REAL(kind=r8), INTENT(OUt  ) :: dku(im,km+1)


    !
    !    locals
    !
    REAL(kind=r8), PARAMETER :: xkzm_m=3.0_r8!   bkgd_vdif_m = 3.0!     xkzm_m   - real, background vertical diffusion for momentum  1    !
    REAL(kind=r8), PARAMETER :: xkzm_h=1.0_r8!   bkgd_vdif_h = 1.0!     xkzm_h   - real, background vertical diffusion for heat, q  1    !
    REAL(kind=r8), PARAMETER :: xkzm_s=0.2_r8!   bkgd_vdif_s = 0.2!     xkzm_s   - real, sigma threshold for background mom. diffusn 1    !


    REAL(kind=r8) :: SIGMAF(IM)! 


    INTEGER       :: i
 !   INTEGER       :: is
    INTEGER       :: k
    INTEGER       :: kk
    INTEGER       :: km1
    INTEGER       :: kmpbl
   ! INTEGER       :: latd
   ! INTEGER       :: lond
    INTEGER       :: lcld(im)
    INTEGER       :: icld(im)
    INTEGER       :: kcld(im)
    INTEGER       :: krad(im)
    INTEGER       :: kx1(im)
        REAL(kind=r8) :: DDVEL(IM)
    REAL(kind=r8) :: Z0RL(IM)  

    !
    REAL(kind=r8) :: hgamq(im)
    REAL(kind=r8) :: hgamt(im)
    REAL(kind=r8) :: dusfc(im)
    REAL(kind=r8) :: dvsfc(im)
    REAL(kind=r8) :: dtsfc(im)
    REAL(kind=r8) :: dqsfc(im)
    REAL(kind=r8) :: CM(IM)
    REAL(kind=r8) :: CH(IM)
    REAL(kind=r8) :: RB(IM)

    REAL(kind=r8) :: phih(im)
    REAL(kind=r8) :: phim(im)
    REAL(kind=r8) :: rbdn(im)
    REAL(kind=r8) :: rbup(im)
    REAL(kind=r8) :: beta(im)
    REAL(kind=r8) :: ustar(im)
    REAL(kind=r8) :: wscale(im)
    REAL(kind=r8) :: thermal(im)
    REAL(kind=r8) :: wstar3(im)
    !
    REAL(kind=r8) :: thvx(im,km)
    REAL(kind=r8) :: thlvx(im,km)
    REAL(kind=r8) :: qlx(im,km)
    REAL(kind=r8) :: thetae(im,km)
    REAL(kind=r8) :: qtx(im,km)
    REAL(kind=r8) :: bf(im,km-1)
    REAL(kind=r8) :: u1(im,km)
    REAL(kind=r8) :: v1(im,km)
    REAL(kind=r8) :: radx(im,km-1)
    REAL(kind=r8) :: govrth(im)
    REAL(kind=r8) :: hrad(im)
    REAL(kind=r8) :: cteit(im)
    REAL(kind=r8) :: radmin(im)
    REAL(kind=r8) :: vrad(im)
    REAL(kind=r8) :: zd(im)
    REAL(kind=r8) :: zdd(im)
    REAL(kind=r8) :: thlvx1(im)

    !
    REAL(kind=r8) :: rdzt(im,km-1)
    REAL(kind=r8) :: dktx(im,km-1)
    REAL(kind=r8) :: dkux(im,km-1)
    REAL(kind=r8) :: zi(im,km+1)
    REAL(kind=r8) :: zl(im,km)
    REAL(kind=r8) :: xkzo(im,km-1)
    REAL(kind=r8) :: xkzmo(im,km-1)
    REAL(kind=r8) :: cku(im,km-1)
    REAL(kind=r8) :: ckt(im,km-1)
   ! REAL(kind=r8) :: al(im,km-1)
   ! REAL(kind=r8) :: ad(im,km)
   ! REAL(kind=r8) :: au(im,km-1)
   ! REAL(kind=r8) :: a1(im,km)
    !REAL(kind=r8) :: a2(im,km)
    REAL(kind=r8) :: theta(im,km)
    !
    REAL(kind=r8) :: prinv(im)
    REAL(kind=r8) :: rent(im)
    INTEGER       :: kpblx(im)
    !
    REAL(kind=r8) :: hpblx(im)
    INTEGER       :: kinver(im)
 !   REAL(kind=r8) :: dv(im,km)
 !   REAL(kind=r8) :: du(im,km)
!    REAL(kind=r8) :: tau(im,km)
    REAL(kind=r8) :: xmu(im)

    REAL(kind=r8)  :: STRESS(IM)
    REAL(kind=r8)  :: FM(IM)
    REAL(kind=r8)  :: FH(IM)
    REAL(kind=r8)  :: WIND(IM)
    REAL(kind=r8) :: FM10(IM)
    REAL(kind=r8) :: FH2(IM)
    REAL(kind=r8) :: rbsoil(im)
    REAL(kind=r8) :: spd1(im)
    REAL(kind=r8) ::tsurf(im)
    LOGICAL        :: flag_iter(im)
    !
    LOGICAL       ::  pblflg(im), sfcflg(im), scuflg(im), flg(im)
    !
    REAL(kind=r8) :: bvf2
    REAL(kind=r8) :: dk
   ! REAL(kind=r8) :: dq1
   ! REAL(kind=r8) :: dsdz2
    !REAL(kind=r8) :: dsdzq
    !REAL(kind=r8) :: dsdzt
    !REAL(kind=r8) :: dsdzu
    !REAL(kind=r8) :: dsdzv
   ! REAL(kind=r8) :: dsig
    REAL(kind=r8) :: dt
  !  REAL(kind=r8) :: dthe1
  !  REAL(kind=r8) :: dtodsd
  !  REAL(kind=r8) :: dtodsu
    REAL(kind=r8) :: dw2
    REAL(kind=r8) :: hol
    REAL(kind=r8) :: hol1
    REAL(kind=r8) :: prnum
!    REAL(kind=r8) :: qtend
    REAL(kind=r8) :: rbint
    REAL(kind=r8) :: rdt
    REAL(kind=r8) :: rdz
    REAL(kind=r8) :: ri
    REAL(kind=r8) :: rl2
!    REAL(kind=r8) :: rone
!    REAL(kind=r8) :: rzero
    REAL(kind=r8) :: sflux
    REAL(kind=r8) :: shr2
    REAL(kind=r8) :: spdk2
    REAL(kind=r8) :: sri
    REAL(kind=r8) :: tem
    REAL(kind=r8) :: ti
   ! REAL(kind=r8) :: ttend
   ! REAL(kind=r8) :: tvd
   ! REAL(kind=r8) :: tvu
  !  REAL(kind=r8) :: utend
 !   REAL(kind=r8) :: vtend
    REAL(kind=r8) :: zfac
    REAL(kind=r8) :: vpert
    REAL(kind=r8) :: zk
    REAL(kind=r8) :: tem1
    REAL(kind=r8) :: tem2
   ! REAL(kind=r8) :: xkzm
   ! REAL(kind=r8) :: xkzmu
    REAL(kind=r8) :: ptem
    REAL(kind=r8) :: ptem1
    REAL(kind=r8) :: ptem2
    REAL(kind=r8) :: tx1(im)
    REAL(kind=r8) :: tx2(im)
    !
   ! REAL(kind=r8) :: u01
   ! REAL(kind=r8) :: v01
    !REAL(kind=r8) :: delu
    !REAL(kind=r8) :: delv
    REAL(kind=r8) :: cc
    REAL(kind=r8) :: ss
    REAL(kind=r8) :: ch2
    !cc
    REAL(KIND=r8), PARAMETER :: gravi=1.0_r8/grav
    REAL(KIND=r8), PARAMETER :: g=grav
   ! REAL(KIND=r8), PARAMETER :: gor=g/con_rd
   ! REAL(KIND=r8), PARAMETER :: gocp=g/con_cp
   ! REAL(KIND=r8), PARAMETER :: cont=con_cp/g
   ! REAL(KIND=r8), PARAMETER :: conq=con_hvap/g
   ! REAL(KIND=r8), PARAMETER :: conw=1.0_r8/g               ! for del in pa
!!!!!!!!     parameter(cont=1000.*con_cp/g,conq=1000.*con_hvap/g,conw=1000./g) ! for del in kpa
    REAL(KIND=r8), PARAMETER :: rlam=30.0_r8
    REAL(KIND=r8), PARAMETER :: vk=0.4_r8
   ! REAL(KIND=r8), PARAMETER :: vk2=vk*vk
    REAL(KIND=r8), PARAMETER :: prmin=0.25_r8
    REAL(KIND=r8), PARAMETER :: prmax=4.0_r8
    REAL(KIND=r8), PARAMETER :: dw2min=0.0001_r8
    REAL(KIND=r8), PARAMETER :: dkmin=0.0_r8
    REAL(KIND=r8), PARAMETER :: dkmax=1000.0_r8
    REAL(KIND=r8), PARAMETER :: rimin=-100.0_r8
    REAL(KIND=r8), PARAMETER :: rbcr=0.25_r8
    REAL(KIND=r8), PARAMETER :: wfac=7.0_r8
    REAL(KIND=r8), PARAMETER :: cfac=6.5_r8
    REAL(KIND=r8), PARAMETER :: pfac=2.0_r8
    REAL(KIND=r8), PARAMETER :: sfcfrac=0.1_r8
    !     parameter(qmin=1.e-8,xkzm=1.0,zfmin=1.e-8,aphi5=5.,aphi16=16.)
    REAL(KIND=r8), PARAMETER :: qmin=1.e-8_r8
    REAL(KIND=r8), PARAMETER :: zfmin=1.e-8_r8
    REAL(KIND=r8), PARAMETER :: aphi5=5.0_r8
    REAL(KIND=r8), PARAMETER :: aphi16=16.0_r8
    REAL(KIND=r8), PARAMETER :: tdzmin=1.e-3_r8
    REAL(KIND=r8), PARAMETER :: qlmin=1.e-12_r8
    !REAL(KIND=r8), PARAMETER :: cpert=0.25
  !  REAL(KIND=r8), PARAMETER :: sfac=5.4
    REAL(KIND=r8), PARAMETER :: h1=0.33333333_r8
 !   REAL(KIND=r8), PARAMETER :: h2=0.66666667
    REAL(KIND=r8), PARAMETER :: cldtime=500.0_r8
    REAL(KIND=r8), PARAMETER :: xkzminv=0.3_r8
    !     parameter(cldtime=500.,xkzmu=3.0,xkzminv=0.3)
    !     parameter(gamcrt=3.,gamcrq=2.e-3,rlamun=150.0)
    REAL(KIND=r8), PARAMETER :: gamcrt=3.0_r8
    REAL(KIND=r8), PARAMETER :: gamcrq=0.0_r8
    REAL(KIND=r8), PARAMETER :: rlamun=150.0_r8
    REAL(KIND=r8), PARAMETER :: rentf1=0.2_r8
    REAL(KIND=r8), PARAMETER :: rentf2=1.0_r8
    REAL(KIND=r8), PARAMETER :: radfac=0.85_r8
 !   REAL(KIND=r8), PARAMETER :: iun=84
    !
    !     parameter (zstblmax = 2500., qlcr=1.0e-5)
    !     parameter (zstblmax = 2500., qlcr=3.0e-5)
    !     parameter (zstblmax = 2500., qlcr=3.5e-5)
    !     parameter (zstblmax = 2500., qlcr=1.0e-4)
    REAL(KIND=r8), PARAMETER :: zstblmax = 2500.0_r8
    REAL(KIND=r8), PARAMETER :: qlcr=3.5e-5_r8
    !     parameter (actei = 0.23)
    REAL(KIND=r8), PARAMETER :: actei = 0.7_r8
    !
    !  ---  constant parameters:
    REAL(kind=r8),PARAMETER:: con_jcal   =4.1855E+0_r8      ! joules per calorie  ()
    REAL(kind=r8), PARAMETER :: cnwatt = -con_jcal*1.0e4_r8/60.0_r8


    !c
    !c-----------------------------------------------------------------------
    !c
!601 FORMAT(1x,' moninp lat lon step hour ',3i6,f6.1)
!602 FORMAT(1x,'    k','        z','        t','       th',         &
!         '      tvh','        q','        u','        v',              &
!         '       sp')
!603 FORMAT(1x,i5,8f9.1)
!604 FORMAT(1x,'  sfc',9x,f9.1,18x,f9.1)
!605 FORMAT(1x,'    k      zl    spd2   thekv   the1v'        &
!         ,' thermal    rbup')
!606 FORMAT(1x,i5,6f8.2)
!607 FORMAT(1x,' kpbl    hpbl      fm      fh   hgamt',         &
!         '   hgamq      ws   ustar      cd      ch')
!608 FORMAT(1x,i5,9f8.2)
!609 FORMAT(1x,' k pr dkt dku ',i5,3f8.2)
!610 FORMAT(1x,' k pr dkt dku ',i5,3f8.2,' l2 ri t2',        &
!         ' sr2  ',2f8.2,2e10.2)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !     compute preliminary variables
    !
    kpbl= 0;hpbl= 0.0_r8;dkt= 0.0_r8;dku= 0.0_r8
    lcld=0;   icld=0;    kcld=0;   krad=0;
    kx1=0;   DDVEL=0.0_r8;hgamq=0.0_r8;
    hgamt=0.0_r8;     dusfc=0.0_r8;dvsfc=0.0_r8;     dtsfc=0.0_r8;
    dqsfc=0.0_r8;     CM=0.0_r8;CH=0.0_r8;     RB=0.0_r8;

    phih=0.0_r8;phim=0.0_r8;rbdn=0.0_r8;rbup=0.0_r8;
    beta=0.0_r8;ustar=0.0_r8;wscale=0.0_r8;thermal=0.0_r8;wstar3=0.0_r8;
    
    thvx=0.0_r8;thlvx=0.0_r8;qlx=0.0_r8;thetae=0.0_r8;qtx=0.0_r8;
    bf=0.0_r8;u1=0.0_r8;v1=0.0_r8;radx=0.0_r8
     
    govrth=0.0_r8;hrad=0.0_r8;cteit=0.0_r8;radmin=0.0_r8;vrad=0.0_r8;zd=0.0_r8;
    zdd=0.0_r8;thlvx1=0.0_r8;rdzt=0.0_r8;dktx=0.0_r8;dkux=0.0_r8;zi=0.0_r8;
    zl=0.0_r8;xkzo=0.0_r8;xkzmo=0.0_r8;cku=0.0_r8;ckt=0.0_r8;theta=0.0_r8;
    prinv=0.0_r8;rent=0.0_r8;kpblx=0;hpblx=0.0_r8;kinver=0;xmu=0.0_r8;

    STRESS=0.0_r8;FM=0.0_r8;FH=0.0_r8;WIND=0.0_r8;FM10=0.0_r8;FH2=0.0_r8;
    rbsoil=0.0_r8;spd1=0.0_r8;tsurf=0.0_r8;tx1=0.0_r8;tx2=0.0_r8

    !
    !     iprt = 0
    !     if(iprt.eq.1) then
    !cc   latd = 0
    !     lond = 0
    !     else
    !cc   latd = 0
    !     lond = 0
    !     endif
    !
    dt    = 2.0_r8 * deltim
    rdt   = 1.0_r8 / dt
    km1   = km - 1
    kmpbl = km / 2
    kinver=km
    DDVEL   = 1.0_r8
    xmu    = 1.0_r8
    cc     = 0.1_r8
    ss     = 0.0_r8
    ch2     = 350.0_r8 / cnwatt
    DO i = 1, im
       IF (xmu(i) > 0.01_r8 .AND. cc > 0.01_r8) THEN
          xmu(i) = xmu(i) / cc
       ELSE
          xmu(i) = 0.0_r8
       ENDIF
    ENDDO


    !
    DO k=1,km
       DO i=1,im
          zi(i,k) = phii(i,k) !* gravi
          zl(i,k) = phil(i,k) !* gravi
          u1(i,k) = uo(i,k)
          v1(i,k) = vo(i,k)
       ENDDO
    ENDDO
    DO i=1,im
       sigmaf(i)   = MAX( vfrac(i),0.01_r8 )
       !          if (lsm == 0) sigmaf(i)   =  0.5_r8 + vfrac(i) * 0.5_r8
       Z0RL(i)  = Z0(i)  

       zi(i,km+1) = phii(i,km+1) !* gravi
       tsurf(i)=TSKIN(i)
    ENDDO
    !
    DO k = 1,km1
       DO i=1,im
          rdzt(i,k) = 1.0_r8 / (zl(i,k+1) - zl(i,k))
       ENDDO
    ENDDO
    !
    DO i=1,im
       flag_iter(i) = .TRUE.
       kx1(i) = 1
       tx1(i) = 1.0_r8 / prsi(i,1)
       tx2(i) = tx1(i)
    ENDDO


    CALL SFC_DIFF(           &
         IM                ,&!integer      , INTENT(IN   )  :: IM
         U1     (1:im,1:1) ,&!real(kind=r8), INTENT(IN   )  :: U1(IM)
         V1     (1:im,1:1) ,&!real(kind=r8), INTENT(IN   )  :: V1(IM)
         T1     (1:im,1:1) ,&!real(kind=r8), INTENT(IN   )  :: T1(IM)
         Q1     (1:im,1:1) ,&!real(kind=r8), INTENT(IN   )  :: Q1(IM)
         zi     (1:im,1:1) ,&!real(kind=r8), INTENT(IN   )  :: Z1(IM)
         taux   (1:im)      , &!real(r8), intent(in)    :: taux(pcols)       ! Zonal wind stress at surface [ N/m2 ]
         tauy   (1:im)      , &!real(r8), intent(in)    :: tauy(pcols)       ! Meridional wind stress at surface [ N/m2 ]
         TSKIN  (1:im)     ,&!real(kind=r8), INTENT(IN   )  :: TSKIN(IM)
         bps     (1:im,1:1) ,&!real(kind=r8), INTENT(IN   )  :: bps(IM)
         Z0RL   (1:im)     ,&!real(kind=r8), INTENT(INOUT)  :: Z0RL(IM)
         CM     (1:im)     ,&!real(kind=r8), INTENT(OUT  )  :: CM(IM)
         CH     (1:im)     ,&!real(kind=r8), INTENT(OUT  )  :: CH(IM)
         RB     (1:im)     ,&!real(kind=r8), INTENT(OUT  )  :: RB(IM)
         prsl   (1:im,1:1) ,&!real(kind=r8), INTENT(IN   )  :: PRSL1(IM)
         prsi   (1:im,1:1) ,&!real(kind=r8), INTENT(IN   )  :: prsi(IM)
         mlsi (1:im)     ,&!real(kind=r8), INTENT(IN   )  :: mlsi(IM)
         STRESS (1:im)     ,&!real(kind=r8), INTENT(OUT  )  :: STRESS(IM)
         FM     (1:im)     ,&!real(kind=r8), INTENT(OUT  )  :: FM(IM)
         FH     (1:im)     ,&!real(kind=r8), INTENT(OUT  )  :: FH(IM)
                                !Clu_q2m_iter [-1L/+2L]: add tsurf, flag_iter
                                !*   &                    USTAR,WIND,DDVEL,FM10,FH2)
         USTAR  (1:im)     ,&!real(kind=r8), INTENT(OUT  )  :: USTAR(IM)
         WIND   (1:im)     ,&!real(kind=r8), INTENT(OUT  )  :: WIND(IM)
         DDVEL  (1:im)     ,&!real(kind=r8), INTENT(IN   )  :: DDVEL(IM)
         FM10   (1:im)     ,&!real(kind=r8), INTENT(OUT  )  :: FM10(IM)
         FH2    (1:im)     ,&!real(kind=r8), INTENT(OUT  )  :: FH2(IM)
         SIGMAF (1:im)     ,&!real(kind=r8), INTENT(IN   )  :: SIGMAF(IM)
         tsurf  (1:im)     ,&!real(kind=r8), INTENT(IN   )  :: TSURF(IM)
         flag_iter(1:im)    )!logical      , INTENT(IN   )   :: flag_iter(im)

    rbsoil=RB
    spd1=WIND
    DO k = 1,km1
       DO i=1,im
          xkzo(i,k)  = 0.0_r8
          xkzmo(i,k) = 0.0_r8
          IF (k < kinver(i)) THEN
             !                                  vertical background diffusivity
             ptem      = prsi(i,k+1) * tx1(i)
             tem1      = 1.0_r8 - ptem
             tem1      = tem1 * tem1 * 10.0_r8
             xkzo(i,k) = xkzm_h * MIN(1.0_r8, EXP(-tem1))

             !                                  vertical background diffusivity for momentum
             IF (ptem >= xkzm_s) THEN
                xkzmo(i,k) = xkzm_m
                kx1(i)     = k + 1
             ELSE
                IF (k == kx1(i) .AND. k > 1) tx2(i) = 1.0_r8 / prsi(i,k)
                tem1 = 1.0_r8 - prsi(i,k+1) * tx2(i)
                tem1 = tem1 * tem1 * 5.0_r8
                xkzmo(i,k) = xkzm_m * MIN(1.0_r8, EXP(-tem1))
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !     if (lprnt) then
    !       print *,' xkzo=',(xkzo(ipr,k),k=1,km1)
    !       print *,' xkzmo=',(xkzmo(ipr,k),k=1,km1)
    !     endif
    !
    !  diffusivity in the inversion layer is set to be xkzminv (m^2/s)
    !
    DO k = 1,kmpbl
       DO i=1,im
          !         if(zi(i,k+1).gt.200..and.zi(i,k+1).lt.zstblmax) then
          IF(zi(i,k+1).GT.250.0_r8) THEN
             tem1 = (t1(i,k+1)-t1(i,k)) * rdzt(i,k)
             IF(tem1 .GT. 1.e-5_r8) THEN
                xkzo(i,k)  = MIN(xkzo(i,k),xkzminv)
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !
    DO i = 1,im
       dusfc(i) = 0.0_r8
       dvsfc(i) = 0.0_r8
       dtsfc(i) = 0.0_r8
       dqsfc(i) = 0.0_r8
       hgamt(i) = 0.0_r8
       hgamq(i) = 0.0_r8
       !        hgamu(i) = 0.0_r8
       !        hgamv(i) = 0.0_r8
       !        hgams(i) = 0.0_r8
       wscale(i)= 0.0_r8
       kpbl(i)  = 1
       kpblx(i) = 1
       hpbl(i)  = zi(i,1)
       hpblx(i) = zi(i,1)
       pblflg(i)= .TRUE.
       sfcflg(i)= .TRUE.
       IF(rbsoil(i).GT.0.0_r8) sfcflg(i) = .FALSE.
       scuflg(i)= .TRUE.
       IF(scuflg(i)) THEN
          radmin(i)= 0.0_r8
          cteit(i) = 0.0_r8
          rent(i)  = rentf1
          hrad(i)  = zi(i,1)
          !          hradm(i) = zi(i,1)
          krad(i)  = 1
          icld(i)  = 0
          lcld(i)  = km1
          kcld(i)  = km1
          zd(i)    = 0.0_r8
       ENDIF
    ENDDO
    !
    DO k = 1,km
       DO i = 1,im
          theta(i,k) = t1(i,k) * bps(i,k) !(i) / (i,k)
          qlx(i,k)   = MAX(qliq(i,k),qlmin)
          qtx(i,k)   = MAX(q1(i,k),qmin)+ qlx(i,k)
          ptem       = qlx(i,k)
          ptem1      = con_hvap*MAX(q1(i,k),qmin)/(con_cp*t1(i,k))
          thetae(i,k)= theta(i,k)*(1.0_r8+ptem1)
          thvx(i,k)  = theta(i,k)*(1.0_r8+con_fvirt*MAX(q1(i,k),qmin)-ptem)
          ptem2      = theta(i,k)-(con_hvap/con_cp)*ptem
          thlvx(i,k) = ptem2*(1.0_r8+con_fvirt*qtx(i,k))
       ENDDO
    ENDDO
    DO k = 1,km1
       DO i = 1,im
          dku(i,k)  = 0.0_r8
          dkt(i,k)  = 0.0_r8
          dktx(i,k) = 0.0_r8
          dkux(i,k) = 0.0_r8
          cku(i,k)  = 0.0_r8
          ckt(i,k)  = 0.0_r8
          tem       = zi(i,k+1)-zi(i,k)
          radx(i,k) = tem*(swh(i,k)*xmu(i)+hlw(i,k))
       ENDDO
    ENDDO
    !
    DO i=1,im
       flg(i)  = scuflg(i)
    ENDDO
    DO k = 1, km1
       DO i=1,im
          IF(flg(i).AND.zl(i,k).GE.zstblmax) THEN
             lcld(i)=k
             flg(i)=.FALSE.
          ENDIF
       ENDDO
    ENDDO
    !
    !  compute buoyancy flux
    !
    DO k = 1, km1
       DO i = 1, im
          bf(i,k) = (thvx(i,k+1)-thvx(i,k))*rdzt(i,k)
       ENDDO
    ENDDO
    !
    DO i = 1,im
       govrth(i) = g/theta(i,1)
    ENDDO
    !
    DO i=1,im
       beta(i)  = dt / (zi(i,2)-zi(i,1))
    ENDDO
    !
    DO i=1,im
       ustar(i) = SQRT(stress(i))
       thermal(i) = thvx(i,1)
    ENDDO
    !
    !  compute the first guess pbl height
    !
    DO i=1,im
       flg(i) = .FALSE.
       rbup(i) = rbsoil(i)
    ENDDO
    DO k = 2, kmpbl
       DO i = 1, im
          IF(.NOT.flg(i)) THEN
             rbdn(i) = rbup(i)
             spdk2   = MAX((u1(i,k)**2+v1(i,k)**2),1.0_r8)
             rbup(i) = (thvx(i,k)-thermal(i))*   &
                  (g*zl(i,k)/thvx(i,1))/spdk2
             kpbl(i) = k
             flg(i)  = rbup(i).GT.rbcr
          ENDIF
       ENDDO
    ENDDO
    DO i = 1,im
       k = kpbl(i)
       IF(rbdn(i).GE.rbcr) THEN
          rbint = 0.0_r8
       ELSEIF(rbup(i).LE.rbcr) THEN
          rbint = 1.0_r8
       ELSE
          rbint = (rbcr-rbdn(i))/(rbup(i)-rbdn(i))
       ENDIF
       hpbl(i) = zl(i,k-1) + rbint*(zl(i,k)-zl(i,k-1))
       IF(hpbl(i).LT.zi(i,kpbl(i))) kpbl(i) = kpbl(i) - 1
       hpblx(i) = hpbl(i)
       kpblx(i) = kpbl(i)
    ENDDO
    !
    !     compute similarity parameters 
    !
    DO i=1,im
       sflux = heat(i) + evap(i)*con_fvirt*theta(i,1)
       IF(sfcflg(i).AND.sflux.GT.0.0_r8) THEN
          hol = MAX(rbsoil(i)*fm(i)*fm(i)/fh(i),rimin)
          hol = MIN(hol,-zfmin)
          !
          hol1 = hol*hpbl(i)/zl(i,1)*sfcfrac
          !          phim(i) = (1.0_r8-aphi16*hol1)**(-1.0_r8/4.0_r8)
          !          phih(i) = (1.0_r8-aphi16*hol1)**(-1.0_r8/2.0_r8)
          tem     = 1.0_r8 / (1.0_r8 - aphi16*hol1)
          phih(i) = SQRT(tem)
          phim(i) = SQRT(phih(i))
          wstar3(i) = govrth(i)*sflux*hpbl(i)
          tem1      = ustar(i)**3.0_r8
          wscale(i) = (tem1+wfac*vk*wstar3(i)*sfcfrac)**h1
          !          wscale(i) = ustar(i)/phim(i)
          !          wscale(i) = min(wscale(i),ustar(i)*aphi16)
          wscale(i) = MAX(wscale(i),ustar(i)/aphi5)
       ELSE
          pblflg(i)=.FALSE.
       ENDIF
    ENDDO
    !
    ! compute counter-gradient mixing term for heat and moisture
    !
    DO i = 1,im
       IF(pblflg(i)) THEN
          hgamt(i)  = MIN(cfac*heat(i)/wscale(i),gamcrt)
          hgamq(i)  = MIN(cfac*evap(i)/wscale(i),gamcrq)
          vpert     = hgamt(i) + hgamq(i)*con_fvirt*theta(i,1)
          vpert     = MIN(vpert,gamcrt)
          thermal(i)= thermal(i)+MAX(vpert,0.0_r8)
          hgamt(i)  = MAX(hgamt(i),0.0_r8)
          hgamq(i)  = MAX(hgamq(i),0.0_r8)
       ENDIF
    ENDDO
    !
    ! compute large-scale mixing term for momentum
    !
    !     do i = 1,im
    !       flg(i) = pblflg(i)
    !       kemx(i)= 1
    !       hpbl01(i)= sfcfrac*hpbl(i)
    !     enddo
    !     do k = 1, kmpbl
    !     do i = 1, im
    !       if(flg(i).and.zl(i,k).gt.hpbl01(i)) then
    !         kemx(i) = k
    !         flg(i)  = .false.
    !       endif
    !     enddo
    !     enddo
    !     do i = 1, im
    !       if(pblflg(i)) then
    !         kk = kpbl(i)
    !         if(kemx(i).le.1) then
    !           ptem  = u1(i,1)/zl(i,1)
    !           ptem1 = v1(i,1)/zl(i,1)
    !           u01   = ptem*hpbl01(i)
    !           v01   = ptem1*hpbl01(i)
    !         else
    !           tem   = zl(i,kemx(i))-zl(i,kemx(i)-1)
    !           ptem  = (u1(i,kemx(i))-u1(i,kemx(i)-1))/tem
    !           ptem1 = (v1(i,kemx(i))-v1(i,kemx(i)-1))/tem
    !           tem1  = hpbl01(i)-zl(i,kemx(i)-1)
    !           u01   = u1(i,kemx(i)-1)+ptem*tem1
    !           v01   = v1(i,kemx(i)-1)+ptem1*tem1
    !         endif
    !         if(kk.gt.kemx(i)) then
    !           delu  = u1(i,kk)-u01
    !           delv  = v1(i,kk)-v01
    !           tem2  = sqrt(delu**2+delv**2)
    !           tem2  = max(tem2,0.1)
    !           ptem2 = -sfac*ustar(i)*ustar(i)*wstar3(i)
    !    1                /(wscale(i)**4.)
    !           hgamu(i) = ptem2*delu/tem2
    !           hgamv(i) = ptem2*delv/tem2
    !           tem  = sqrt(u1(i,kk)**2+v1(i,kk)**2)
    !           tem1 = sqrt(u01**2+v01**2)
    !           ptem = tem - tem1
    !           if(ptem.gt.0.) then
    !             hgams(i)=-sfac*vk*sfcfrac*wstar3(i)/(wscale(i)**3.)
    !           else
    !             hgams(i)=sfac*vk*sfcfrac*wstar3(i)/(wscale(i)**3.)
    !           endif
    !         else
    !           hgams(i) = 0.
    !         endif
    !       endif
    !     enddo
    !
    !  enhance the pbl height by considering the thermal excess
    !
    DO i=1,im
       flg(i)  = .TRUE.
       IF(pblflg(i)) THEN
          flg(i)  = .FALSE.
          rbup(i) = rbsoil(i)
       ENDIF
    ENDDO
    DO k = 2, kmpbl
       DO i = 1, im
          IF(.NOT.flg(i)) THEN
             rbdn(i) = rbup(i)
             spdk2   = MAX((u1(i,k)**2+v1(i,k)**2),1.0_r8)
             rbup(i) = (thvx(i,k)-thermal(i))*      &
                  (g*zl(i,k)/thvx(i,1))/spdk2
             kpbl(i) = k
             flg(i)  = rbup(i).GT.rbcr
          ENDIF
       ENDDO
    ENDDO
    DO i = 1,im
       IF(pblflg(i)) THEN
          k = kpbl(i)
          IF(rbdn(i).GE.rbcr) THEN
             rbint = 0.0_r8
          ELSEIF(rbup(i).LE.rbcr) THEN
             rbint = 1.0_r8
          ELSE
             rbint = (rbcr-rbdn(i))/(rbup(i)-rbdn(i))
          ENDIF
          hpbl(i) = zl(i,k-1) + rbint*(zl(i,k)-zl(i,k-1))
          IF(hpbl(i).LT.zi(i,kpbl(i))) kpbl(i) = kpbl(i) - 1
          IF(kpbl(i).LE.1) pblflg(i) = .FALSE.
       ENDIF
    ENDDO
    !
    !  look for stratocumulus
    !
    DO i = 1, im
       flg(i)=scuflg(i)
    ENDDO
    DO k = kmpbl,1,-1
       DO i = 1, im
          IF(flg(i).AND.k.LE.lcld(i)) THEN
             IF(qlx(i,k).GE.qlcr) THEN
                kcld(i)=k
                flg(i)=.FALSE.
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    DO i = 1, im
       IF(scuflg(i).AND.kcld(i).EQ.km1) scuflg(i)=.FALSE.
    ENDDO
    !
    DO i = 1, im
       flg(i)=scuflg(i)
    ENDDO
    DO k = kmpbl,1,-1
       DO i = 1, im
          IF(flg(i).AND.k.LE.kcld(i)) THEN
             IF(qlx(i,k).GE.qlcr) THEN
                IF(radx(i,k).LT.radmin(i)) THEN
                   radmin(i)=radx(i,k)
                   krad(i)=k
                ENDIF
             ELSE
                flg(i)=.FALSE.
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    DO i = 1, im
       IF(scuflg(i).AND.krad(i).LE.1) scuflg(i)=.FALSE.
       IF(scuflg(i).AND.radmin(i).GE.0.0_r8) scuflg(i)=.FALSE.
    ENDDO
    !
    DO i = 1, im
       flg(i)=scuflg(i)
    ENDDO
    DO k = kmpbl,2,-1
       DO i = 1, im
          IF(flg(i).AND.k.LE.krad(i)) THEN
             IF(qlx(i,k).GE.qlcr) THEN
                icld(i)=icld(i)+1
             ELSE
                flg(i)=.FALSE.
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    DO i = 1, im
       IF(scuflg(i).AND.icld(i).LT.1) scuflg(i)=.FALSE.
    ENDDO
    !
    DO i = 1, im
       IF(scuflg(i)) THEN
          hrad(i) = zi(i,krad(i)+1)
          !          hradm(i)= zl(i,krad(i))
       ENDIF
    ENDDO
    !
    DO i = 1, im
       IF(scuflg(i).AND.hrad(i).LT.zi(i,2)) scuflg(i)=.FALSE.
    ENDDO
    !
    DO i = 1, im
       IF(scuflg(i)) THEN
          k    = krad(i)
          tem  = zi(i,k+1)-zi(i,k)
          tem1 = cldtime*radmin(i)/tem
          thlvx1(i) = thlvx(i,k)+tem1
          !         if(thlvx1(i).gt.thlvx(i,k-1)) scuflg(i)=.false.
       ENDIF
    ENDDO
    ! 
    DO i = 1, im
       flg(i)=scuflg(i)
    ENDDO
    DO k = kmpbl,1,-1
       DO i = 1, im
          IF(flg(i).AND.k.LE.krad(i))THEN
             IF(thlvx1(i).LE.thlvx(i,k))THEN
                tem=zi(i,k+1)-zi(i,k)
                zd(i)=zd(i)+tem
             ELSE
                flg(i)=.FALSE.
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    DO i = 1, im
       IF(scuflg(i))THEN
          kk = MAX(1, krad(i)+1-icld(i))
          zdd(i) = hrad(i)-zi(i,kk)
       ENDIF
    ENDDO
    DO i = 1, im
       IF(scuflg(i))THEN
          zd(i) = MAX(zd(i),zdd(i))
          zd(i) = MIN(zd(i),hrad(i))
          tem   = govrth(i)*zd(i)*(-radmin(i))
          vrad(i)= tem**h1
       ENDIF
    ENDDO
    !
    !     compute inverse Prandtl number
    !
    DO i = 1, im
       IF(pblflg(i)) THEN
          tem = phih(i)/phim(i)+cfac*vk*sfcfrac
          !         prinv(i) = (1.0_r8-hgams(i))/tem
          prinv(i) = 1.0_r8 / tem
          prinv(i) = MIN(prinv(i),prmax)
          prinv(i) = MAX(prinv(i),prmin)
       ENDIF
    ENDDO
    !
    !     compute diffusion coefficients below pbl
    !
    DO k = 1, kmpbl
       DO i=1,im
          IF(pblflg(i).AND.k.LT.kpbl(i)) THEN
             !           zfac = max((1.-(zi(i,k+1)-zl(i,1))/
             !    1             (hpbl(i)-zl(i,1))), zfmin)
             zfac = MAX((1.0_r8-zi(i,k+1)/hpbl(i)), zfmin)
             tem = wscale(i)*vk*zi(i,k+1)*zfac**pfac
             !           dku(i,k) = xkzo(i,k)+wscale(i)*vk*zi(i,k+1)
             !    1                 *zfac**pfac
             dku(i,k) = xkzmo(i,k) + tem
             dkt(i,k) = xkzo(i,k)  + tem * prinv(i)
             dku(i,k) = MIN(dku(i,k),dkmax)
             !           dku(i,k) = max(dku(i,k),xkzmo(i,k))
             dkt(i,k) = MIN(dkt(i,k),dkmax)
             !           dkt(i,k) = max(dkt(i,k),xkzo(i,k))
             dktx(i,k)= dkt(i,k)
             dkux(i,k)= dku(i,k)
          ENDIF
       ENDDO
    ENDDO
    !
    ! compute diffusion coefficients based on local scheme
    !
    DO i = 1, im
       IF(.NOT.pblflg(i)) THEN
          kpbl(i) = 1
       ENDIF
    ENDDO
    DO k = 1, km1
       DO i=1,im
          IF(k.GE.kpbl(i)) THEN
             rdz  = rdzt(i,k)
             ti   = 2.0_r8/(t1(i,k)+t1(i,k+1))
             dw2  = (u1(i,k)-u1(i,k+1))**2        &
                  +(v1(i,k)-v1(i,k+1))**2
             shr2 = MAX(dw2,dw2min)*rdz*rdz
             bvf2 = g*bf(i,k)*ti
             ri   = MAX(bvf2/shr2,rimin)
             zk   = vk*zi(i,k+1)
             IF(ri.LT.0.0_r8) THEN ! unstable regime
                rl2      = zk*rlamun/(rlamun+zk)
                dk       = rl2*rl2*SQRT(shr2)
                sri      = SQRT(-ri)
                dku(i,k) = xkzmo(i,k) + dk*(1+8.0_r8*(-ri)/(1+1.746_r8*sri))
                dkt(i,k) = xkzo(i,k)  + dk*(1+8.0_r8*(-ri)/(1+1.286_r8*sri))
             ELSE             ! stable regime
                rl2      = zk*rlam/(rlam+zk)
                !!                tem      = rlam * sqrt(0.01_r8*prsi(i,k))
                !!                rl2      = zk*tem/(tem+zk)
                dk       = rl2*rl2*SQRT(shr2)
                tem1     = dk/(1+5.0_r8*ri)**2
                IF(k.GE.kpblx(i)) THEN
                   prnum = 1.0_r8 + 2.1_r8*ri
                   prnum = MIN(prnum,prmax)
                ELSE
                   prnum = 1.0_r8
                ENDIF
                dkt(i,k) = xkzo(i,k)  + tem1
                dku(i,k) = xkzmo(i,k) + tem1 * prnum
             ENDIF
             !
             dku(i,k) = MIN(dku(i,k),dkmax)
             !              dku(i,k) = max(dku(i,k),xkzmo(i,k))
             dkt(i,k) = MIN(dkt(i,k),dkmax)
             !              dkt(i,k) = max(dkt(i,k),xkzo(i,k))
             !
          ENDIF
          !
       ENDDO
    ENDDO
    !
    !  compute diffusion coefficients for cloud-top driven diffusion
    !  if the condition for cloud-top instability is met,
    !    increase entrainment flux at cloud top
    !
    DO i = 1, im
       IF(scuflg(i)) THEN
          k = krad(i)
          tem = thetae(i,k) - thetae(i,k+1)
          tem1 = qtx(i,k) - qtx(i,k+1)
          IF (tem.GT.0.0_r8.AND.tem1.GT.0.0_r8) THEN
             cteit(i)= con_cp*tem/(con_hvap*tem1)
             IF(cteit(i).GT.actei) rent(i) = rentf2
          ENDIF
       ENDIF
    ENDDO
    DO i = 1, im
       IF(scuflg(i)) THEN
          k = krad(i)
          tem1  = MAX(bf(i,k),tdzmin)
          ckt(i,k) = -rent(i)*radmin(i)/tem1
          cku(i,k) = ckt(i,k)
       ENDIF
    ENDDO
    !
    DO k = 1, kmpbl
       DO i=1,im
          IF(scuflg(i).AND.k.LT.krad(i)) THEN
             tem1=hrad(i)-zd(i)
             tem2=zi(i,k+1)-tem1
             IF(tem2.GT.0.0_r8) THEN
                ptem= tem2/zd(i)
                IF(ptem.GE.1.0_r8) ptem= 1.0_r8
                ptem= tem2*ptem*SQRT(1.0_r8-ptem)
                ckt(i,k) = radfac*vk*vrad(i)*ptem
                cku(i,k) = 0.75_r8*ckt(i,k)
                ckt(i,k) = MAX(ckt(i,k),dkmin)
                ckt(i,k) = MIN(ckt(i,k),dkmax)
                cku(i,k) = MAX(cku(i,k),dkmin)
                cku(i,k) = MIN(cku(i,k),dkmax)
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    DO k = 1, kmpbl
       DO i=1,im
          IF(scuflg(i)) THEN
             dkt(i,k) = dkt(i,k)+ckt(i,k)
             dku(i,k) = dku(i,k)+cku(i,k)
             dkt(i,k) = MIN(dkt(i,k),dkmax)
             dku(i,k) = MIN(dku(i,k),dkmax)
          ENDIF
       ENDDO
    ENDDO
    !
    !     compute tridiagonal matrix elements for heat and moisture
    !
    !PK     do i=1,im
    !PK        ad(i,1) = 1.
    !PK        a1(i,1) = t1(i,1)   + beta(i) * heat(i)
    !PK        a2(i,1) = q1(i,1,1) + beta(i) * evap(i)
    !PK     enddo

    !PK     if(ntrac.ge.2) then
    !PK       do k = 2, ntrac
    !PK         is = (k-1) * km
    !PK         do i = 1, im
    !PK           a2(i,1+is) = q1(i,1,k)
    !PK         enddo
    !PK       enddo
    !PK     endif
    !
    !PK     do k = 1,km1
    !PK       do i = 1,im
    !PK         dtodsd = dt/del(i,k)
    !PK         dtodsu = dt/del(i,k+1)
    !PK         dsig   = prsl(i,k)-prsl(i,k+1)
    !         rdz    = rdzt(i,k)*2./(t1(i,k)+t1(i,k+1))
    !PK         rdz    = rdzt(i,k)
    !PK         tem1   = dsig * dkt(i,k) * rdz

    !PK         if(pblflg(i).and.k.lt.kpbl(i)) then
    !            dsdzt = dsig*dkt(i,k)*rdz*(gocp-hgamt(i)/hpbl(i))
    !            dsdzq = dsig*dkt(i,k)*rdz*(-hgamq(i)/hpbl(i))
    !PK            ptem1 = dsig * dktx(i,k) * rdz
    !PK            tem   = 1.0 / hpbl(i)
    !PK            dsdzt = tem1 * gocp - ptem1*hgamt(i)*tem
    !PK            dsdzq = ptem1 * (-hgamq(i)*tem)
    !PK            a2(i,k)   = a2(i,k)+dtodsd*dsdzq
    !PK            a2(i,k+1) = q1(i,k+1,1)-dtodsu*dsdzq
    !PK         else
    !            dsdzt = dsig*dkt(i,k)*rdz*(gocp)
    !PK            dsdzt = tem1 * gocp
    !PK            a2(i,k+1) = q1(i,k+1,1)
    !PK         endif
    !         dsdz2 = dsig*dkt(i,k)*rdz*rdz
    !PK         dsdz2     = tem1 * rdz
    !PK         au(i,k)   = -dtodsd*dsdz2
    !PK         al(i,k)   = -dtodsu*dsdz2
    !PK         ad(i,k)   = ad(i,k)-au(i,k)
    !PK         ad(i,k+1) = 1.-al(i,k)
    !PK         a1(i,k)   = a1(i,k)+dtodsd*dsdzt
    !PK         a1(i,k+1) = t1(i,k+1)-dtodsu*dsdzt
    !PK       enddo
    !PK     enddo

    !PK     if(ntrac.ge.2) then
    !PK       do kk = 2, ntrac
    !PK         is = (kk-1) * km
    !PK         do k = 1, km1
    !PK           do i = 1, im
    !PK             a2(i,k+1+is) = q1(i,k+1,kk)
    !PK           enddo
    !PK         enddo
    !PK       enddo
    !PK     endif
    !
    !     solve tridiagonal problem for heat and moisture
    !
    !PK      call tridin(im,km,ntrac,al,ad,au,a1,a2,au,a1,a2)

    !
    !     recover tendencies of heat and moisture
    !
    !PK     do  k = 1,km
    !PK        do i = 1,im
    !PK           ttend      = (a1(i,k)-t1(i,k))*rdt
    !PK           qtend      = (a2(i,k)-q1(i,k,1))*rdt
    !PK           tau(i,k)   = tau(i,k)+ttend
    !PK           rtg(i,k,1) = rtg(i,k,1)+qtend
    !PK           dtsfc(i)   = dtsfc(i)+cont*del(i,k)*ttend
    !PK           dqsfc(i)   = dqsfc(i)+conq*del(i,k)*qtend
    !PK        enddo
    !PK     enddo
    !PK     if(ntrac.ge.2) then
    !PK       do kk = 2, ntrac
    !PK         is = (kk-1) * km
    !PK         do k = 1, km 
    !PK           do i = 1, im
    !PK             qtend = (a2(i,k+is)-q1(i,k,kk))*rdt
    !PK             rtg(i,k,kk) = rtg(i,k,kk)+qtend
    !PK           enddo
    !PK         enddo
    !PK       enddo
    !PK     endif
    !
    !     compute tridiagonal matrix elements for momentum
    !
    !PK      do i=1,im
    !PK         ad(i,1) = 1.0 + beta(i) * stress(i) / spd1(i)
    !PK         a1(i,1) = u1(i,1)
    !PK         a2(i,1) = v1(i,1)
    !PK      enddo
    !
    !PK      do k = 1,km1
    !PK        do i=1,im
    !PK          dtodsd = dt/del(i,k)
    !PK          dtodsu = dt/del(i,k+1)
    !PK          dsig   = prsl(i,k)-prsl(i,k+1)
    !PK          rdz    = rdzt(i,k)
    !PK          tem1   = dsig*dku(i,k)*rdz
    !         if(pblflg(i).and.k.lt.kpbl(i))then
    !           ptem1 = dsig*dkux(i,k)*rdz
    !           dsdzu = ptem1*(-hgamu(i)/hpbl(i))
    !           dsdzv = ptem1*(-hgamv(i)/hpbl(i))
    !           a1(i,k)   = a1(i,k)+dtodsd*dsdzu
    !           a1(i,k+1) = u1(i,k+1)-dtodsu*dsdzu
    !           a2(i,k)   = a2(i,k)+dtodsd*dsdzv
    !           a2(i,k+1) = v1(i,k+1)-dtodsu*dsdzv
    !         else
    !PK            a1(i,k+1) = u1(i,k+1)
    !PK            a2(i,k+1) = v1(i,k+1)
    !         endif
    !         dsdz2     = dsig*dku(i,k)*rdz*rdz
    !PK          dsdz2     = tem1*rdz
    !PK          au(i,k)   = -dtodsd*dsdz2
    !PK          al(i,k)   = -dtodsu*dsdz2
    !PK          ad(i,k)   = ad(i,k)-au(i,k)
    !PK          ad(i,k+1) = 1.-al(i,k)
    !PK        enddo
    !PK      enddo
    !
    !     solve tridiagonal problem for momentum
    !
    !PK      call tridi2(im,km,al,ad,au,a1,a2,au,a1,a2)
    !
    !     recover tendencies of momentum
    !
    !PK      do k = 1,km
    !PK         do i = 1,im
    !PK           utend = (a1(i,k)-u1(i,k))*rdt
    !PK            vtend = (a2(i,k)-v1(i,k))*rdt
    !PK            du(i,k)  = du(i,k)  + utend
    !PK            dv(i,k)  = dv(i,k)  + vtend
    !PK            dusfc(i) = dusfc(i) + conw*del(i,k)*utend
    !PK            dvsfc(i) = dvsfc(i) + conw*del(i,k)*vtend
    !PK         enddo
    !PK      enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  pbl height for diagnostic purpose
    !
    !      do i = 1, im
    !         hpbl(i) = hpblx(i)
    !         kpbl(i) = kpblx(i)
    !      enddo
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    RETURN
  END SUBROUTINE moninq


  SUBROUTINE SFC_DIFF( &
       IM     ,&!integer      , INTENT(IN   )  :: IM
       U1     ,&!real(kind=r8), INTENT(IN   )  :: U1(IM)
       V1     ,&!real(kind=r8), INTENT(IN   )  :: V1(IM)
       T1     ,&!real(kind=r8), INTENT(IN   )  :: T1(IM)
       Q1     ,&!real(kind=r8), INTENT(IN   )  :: Q1(IM)
       Z1     ,&!real(kind=r8), INTENT(IN   )  :: Z1(IM)
       taux   , &!real(r8), intent(in)    :: taux(pcols)       ! Zonal wind stress at surface [ N/m2 ]
       tauy   , &!real(r8), intent(in)    :: tauy(pcols)       ! Meridional wind stress at surface [ N/m2 ]
       TSKIN  ,&!real(kind=r8), INTENT(IN   )  :: TSKIN(IM)
       bps    ,&!real(kind=r8), INTENT(IN   )  :: bps(IM)
       Z0RL   ,&!real(kind=r8), INTENT(INOUT)  :: Z0RL(IM)
       CM     ,&!real(kind=r8), INTENT(OUT  )  :: CM(IM)
       CH     ,&!real(kind=r8), INTENT(OUT  )  :: CH(IM)
       RB     ,&!real(kind=r8), INTENT(OUT  )  :: RB(IM)
       PRSL1  ,&!real(kind=r8), INTENT(IN   )  :: PRSL1(IM)
       prsi   ,&!real(kind=r8), INTENT(IN   )  :: prsi(IM)
       mlsi ,&!real(kind=r8), INTENT(IN   )  :: mlsi(IM)
       STRESS ,&!real(kind=r8), INTENT(OUT  )  :: STRESS(IM)
       FM     ,&!real(kind=r8), INTENT(OUT  )  :: FM(IM)
       FH     ,&!real(kind=r8), INTENT(OUT  )  :: FH(IM)
                                !Clu_q2m_iter [-1L/+2L]: add tsurf, flag_iter
                                !*   &                    USTAR,WIND,DDVEL,FM10,FH2)
       USTAR  ,&!real(kind=r8), INTENT(OUT  )  :: USTAR(IM)
       WIND   ,&!real(kind=r8), INTENT(OUT  )  :: WIND(IM)
       DDVEL  ,&!real(kind=r8), INTENT(IN   )  :: DDVEL(IM)
       FM10   ,&!real(kind=r8), INTENT(OUT  )  :: FM10(IM)
       FH2    ,&!real(kind=r8), INTENT(OUT  )  :: FH2(IM)
       SIGMAF ,&!real(kind=r8), INTENT(IN   )  :: SIGMAF(IM)
       tsurf  ,&!real(kind=r8), INTENT(IN   )  :: TSURF(IM)
       flag_iter)!logical      , INTENT(IN   )   :: flag_iter(im)
    !
    !      USE MACHINE , ONLY : r8
    !      USE FUNCPHYS, ONLY : fpvs    
    !      USE PHYSCONS, grav => con_g, SBC => con_sbc 
    !     &,             CP => con_CP, HFUS => con_HFUS 
    !     &,             RVRDM1 => con_FVirt, RD => con_RD
    !     &,             EPS => con_eps, EPSM1 => con_epsm1

    IMPLICIT NONE
    !
    INTEGER      , INTENT(IN   )  :: IM
    REAL(kind=r8), INTENT(IN   )  :: U1(IM)
    REAL(kind=r8), INTENT(IN   )  :: V1(IM)
    REAL(kind=r8), INTENT(IN   )  :: T1(IM)
    REAL(kind=r8), INTENT(IN   )  :: Q1(IM)
    REAL(kind=r8), INTENT(IN   )  :: Z1(IM)
    REAL(kind=r8), INTENT(IN   )  :: taux (IM)!real(r8), intent(in)    :: taux(pcols)       ! Zonal wind stress at surface [ N/m2 ]
    REAL(kind=r8), INTENT(IN   )  :: tauy (IM)!real(r8), intent(in)    :: tauy(pcols)       ! Meridional wind stress at surface [ N/m2 ]
    REAL(kind=r8), INTENT(IN   )  :: TSKIN(IM)
    real(kind=r8), INTENT(IN   )  :: bps(IM)
    REAL(kind=r8), INTENT(INOUT)  :: Z0RL(IM)
    REAL(kind=r8), INTENT(OUT  )  :: CM(IM)
    REAL(kind=r8), INTENT(OUT  )  :: CH(IM)
    REAL(kind=r8), INTENT(OUT  )  :: RB(IM)
    REAL(kind=r8), INTENT(IN   )  :: PRSL1(IM)
    real(kind=r8), INTENT(IN   )  :: prsi(IM)
    INTEGER(KIND=i8)      , INTENT(IN   )  :: mlsi(IM)
    REAL(kind=r8), INTENT(OUT  )  :: STRESS(IM)
    REAL(kind=r8), INTENT(OUT  )  :: FM(IM)
    REAL(kind=r8), INTENT(OUT  )  :: FH(IM)
    REAL(kind=r8), INTENT(OUT  )  :: USTAR(IM)
    REAL(kind=r8), INTENT(OUT  )  :: WIND(IM)
    REAL(kind=r8), INTENT(IN   )  :: DDVEL(IM)
    REAL(kind=r8), INTENT(OUT  )  :: FM10(IM)
    REAL(kind=r8), INTENT(OUT  )  :: FH2(IM)
    REAL(kind=r8), INTENT(IN   )  :: SIGMAF(IM)
    REAL(kind=r8), INTENT(IN   )  :: TSURF(IM)

    !lu_q2m_iter [+1L]: add flag_iter
    LOGICAL      , INTENT(IN   )   :: flag_iter(im)

    !
    !     Locals
    !
    INTEGER              i
    !
    REAL(kind=r8) :: DTV(IM)
    REAL(kind=r8) :: HL1(IM)
    REAL(kind=r8) :: HL12(IM)
    REAL(kind=r8) :: HLINF(IM)
    REAL(kind=r8) :: PH(IM)
    REAL(kind=r8) :: PH2(IM)
    REAL(kind=r8) :: PM(IM)
    REAL(kind=r8) :: PM10(IM)
    REAL(kind=r8) :: Q0(IM)
    REAL(kind=r8) :: RAT(IM)
    REAL(kind=r8) :: THETA1(IM)
    REAL(kind=r8) :: THV1(IM)
    REAL(kind=r8) :: TV1(IM)
    REAL(kind=r8) :: TVS(IM)
    REAL(kind=r8) :: Z0(IM)
    REAL(kind=r8) :: Z0MAX(IM)
    REAL(kind=r8) :: ZTMAX(IM)
    REAL(kind=r8) :: QS1(IM)
    REAL(kind=r8) :: Z1I(IM)
    REAL(kind=r8) :: rrho(IM)

    !
    REAL(kind=r8) :: aa
    REAL(kind=r8) :: aa0
    REAL(kind=r8) :: adtv
    REAL(kind=r8) :: bb
    REAL(kind=r8) :: bb0
!    REAL(kind=r8) :: cc
    REAL(kind=r8) :: cq
    REAL(kind=r8) :: fms
    REAL(kind=r8) :: fhs
    REAL(kind=r8) :: hl0
    REAL(kind=r8) :: hl0inf
    REAL(kind=r8) :: hl110
    REAL(kind=r8) :: hlt
    REAL(kind=r8) :: hltinf
    REAL(kind=r8) :: OLINF
    REAL(kind=r8) :: restar
    !
    !c

    REAL(KIND=r8),  PARAMETER  :: CHARNOCK=.014_r8
    REAL(KIND=r8),  PARAMETER  :: CA=0.4_r8           !C CA IS THE VON KARMAN CONSTANT
    REAL(KIND=r8),  PARAMETER  :: G=grav
    REAL(KIND=r8),  PARAMETER  :: ALPHA=5.0_r8
    REAL(KIND=r8),  PARAMETER  :: A0=-3.975_r8
    REAL(KIND=r8),  PARAMETER  :: A1=12.32_r8
    REAL(KIND=r8),  PARAMETER  :: B1=-7.755_r8
    REAL(KIND=r8),  PARAMETER  :: B2=6.041_r8
    REAL(KIND=r8),  PARAMETER  :: A0P=-7.941_r8
    REAL(KIND=r8),  PARAMETER  :: A1P=24.75_r8
    REAL(KIND=r8),  PARAMETER  :: B1P=-8.705_r8
    REAL(KIND=r8),  PARAMETER  :: B2P=7.899_r8
    REAL(KIND=r8),  PARAMETER  :: VIS=1.4E-5_r8
!    REAL(KIND=r8),  PARAMETER  :: AA1=-1.076_r8
!    REAL(KIND=r8),  PARAMETER  :: BB1=0.7045_r8
!    REAL(KIND=r8),  PARAMETER  :: CC1=-0.05808_r8
!    REAL(KIND=r8),  PARAMETER  :: BB2=-0.1954_r8
!    REAL(KIND=r8),  PARAMETER  :: CC2=0.009999_r8
!    REAL(KIND=r8),  PARAMETER  :: RNU=1.51E-5_r8
!    REAL(KIND=r8),  PARAMETER  :: ARNU=0.135_r8*RNU
    REAL(R8),PARAMETER :: SHR_CONST_BOLTZ   = 1.38065e-23_R8  ! Boltzmann's constant ~ J/K/molecule
    REAL(R8),PARAMETER :: SHR_CONST_AVOGAD  = 6.02214e26_R8   ! Avogadro's number ~ molecules/kmole

    REAL(R8),PARAMETER :: SHR_CONST_RGAS    = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ       ! Universal gas constant ~ J/K/kmole
    REAL(R8),PARAMETER :: SHR_CONST_MWDAIR  = 28.966_R8       ! molecular weight dry air ~ kg/kmole

    REAL(R8),PARAMETER :: SHR_CONST_RDAIR   = SHR_CONST_RGAS/SHR_CONST_MWDAIR        ! Dry air gas constant     ~ J/K/kg
    REAL(r8),PARAMETER :: rair              = shr_const_rdair   ! Dry air gas constant     ~ J/K/kg
    REAL(r8),         PARAMETER :: ustar_min      =   0.01_r8     ! Minimum permitted value of ustar [ m/s ] 



    REAL (KIND=r8), PARAMETER   :: EPS    =    con_rd/con_rv        ! #
    REAL (KIND=r8), PARAMETER   :: EPSM1  =    con_rd/con_rv-1.0_r8 ! #
    REAL (KIND=r8), PARAMETER   :: rvrdm1=con_fvirt


    CM(IM)=0.0_r8;CH=0.0_r8;RB=0.0_r8;STRESS=0.0_r8;FM=0.0_r8;FH=0.0_r8
    USTAR=0.0_r8;WIND=0.0_r8;FM10=0.0_r8;FH2=0.0_r8;DTV=0.0_r8;HL1=0.0_r8
    HL12=0.0_r8;HLINF=0.0_r8;PH=0.0_r8;PH2=0.0_r8;PM=0.0_r8;PM10=0.0_r8
    Q0=0.0_r8;RAT=0.0_r8;THETA1=0.0_r8;THV1=0.0_r8;TV1=0.0_r8;TVS=0.0_r8
    Z0=0.0_r8;Z0MAX=0.0_r8;ZTMAX=0.0_r8;QS1=0.0_r8;Z1I=0.0_r8
    ! MBEK -- TOGA-COARE FLUX ALGORITHM
    !     PARAMETER (RNU=1.51E-5,ARNU=0.11*RNU)
    !
    !  INITIALIZE VARIABLES. ALL UNITS ARE SUPPOSEDLY M.K.S. UNLESS SPECIFIED
    !  PS IS IN PASCALS
    !  WIND IS WIND SPEED, THETA1 IS ADIABATIC SURFACE TEMP FROM LEVEL 1
    !  SURFACE ROUGHNESS LENGTH IS CONVERTED TO M FROM CM
    !
    ! Compute ustar, and kinematic surface fluxes from surface energy fluxes
    DO i = 1, IM
       rrho(i)    = rair * TSKIN(i) / prsi(i)
       ustar(i)   = MAX( SQRT( SQRT( taux(i)**2 + tauy(i)**2 ) * rrho(i) ), ustar_min )
    END DO

    DO I=1,IM
       IF(flag_iter(i)) THEN 
          !**       TSURF(I)  = TSKIN(I)                 !! <---- Clu_q2m_iter [-1L]
          WIND(I)   = SQRT(U1(I) * U1(I) + V1(I) * V1(I))         &
               + MAX(0.0_r8, MIN(DDVEL(I), 30.0_r8))
          WIND(I)   = MAX(WIND(I),1.0_r8)
          Q0(I)     = MAX(Q1(I),1.E-8_r8)
          THETA1(I) = T1(I) * bps(i) !I(I)
          TV1(I)    = T1(I) * (1.0_r8 + RVRDM1 * Q0(I))
          THV1(I)   = THETA1(I) * (1.0_r8 + RVRDM1 * Q0(I))
          !lu_q2m_iter[-1L/+2L]: TVS is computed from avg(tsurf,tskin)
          !**       TVS(I)    = TSURF(I) * (1.0_r8 + RVRDM1 * Q0(I))
          TVS(I)    = 0.5_r8 * (TSURF(I)+TSKIN(I)) * (1.0_r8 + RVRDM1 * Q0(I))
          qs1(i)    = fpvs2es5(T1(i))
          QS1(I)    = EPS * QS1(I) / (PRSL1(I) + EPSM1 * QS1(I))
          QS1(I)    = MAX(QS1(I), 1.E-8_r8)
          Q0(I)     = MIN(QS1(I),Q0(I))

          Z0(I)     = 0.01_r8 * Z0RL(i)
          !         Z1(I)     = -RD * TV1(I) * LOG(PS1(I)/PSURF(I)) / G
          Z1I(I)    = 1.0_r8 / Z1(I)
       ENDIF
    ENDDO
    !!
    !
    !  COMPUTE STABILITY DEPENDENT EXCHANGE COEFFICIENTS
    !
    !  THIS PORTION OF THE CODE IS PRESENTLY SUPPRESSED
    !
    DO I=1,IM
       IF(flag_iter(i)) THEN 
          !IF(mlsi(I).EQ.0_i8) THEN
          !   USTAR(I) = SQRT(G * Z0(I) / CHARNOCK)
          !ENDIF
          !
          !  COMPUTE STABILITY INDICES (RB AND HLINF)
          !
          Z0MAX(I) = MIN(Z0(I),1.0_r8 * Z1(I))

          !**  test xubin's new z0
          !IF (mlsi(I) .NE. 0.0_r8) THEN
          !   Z0MAX(I) = EXP( ((1.0_r8-SHDMAX(I))**2)*LOG(0.01_r8)+        &
          !        (1-((1.0_r8-SHDMAX(I))**2))*LOG(Z0MAX(I)) )
          !   IF (VEGTYPE(I) == 7) THEN
          !      Z0MAX(I) = EXP( ((1.0_r8-SHDMAX(I))**2)*LOG(0.01_r8)+     &
          !           (1-((1.0_r8-SHDMAX(I))**2))*LOG(0.07_r8) )
          !   ENDIF
          !   IF (VEGTYPE(I) == 8) THEN
          !      Z0MAX(I) = EXP( ((1.0_r8-SHDMAX(I))**2)*LOG(0.01_r8)+    &
          !           (1-((1.0_r8-SHDMAX(I))**2))*LOG(0.05_r8) )
          !   ENDIF
          !   IF (VEGTYPE(I) == 9) THEN
          !      Z0MAX(I) = EXP( ((1.0_r8-SHDMAX(I))**2)*LOG(0.01_r8)+    &
          !           (1-((1.0_r8-SHDMAX(I))**2))*LOG(0.01_r8) )
          !   ENDIF
          !   IF (VEGTYPE(I) == 11) THEN
          !      Z0MAX(I) = EXP( ((1.0_r8-SHDMAX(I))**2)*LOG(0.01_r8)+       &
          !           (1-((1.0_r8-SHDMAX(I))**2))*LOG(0.01_r8) )
          !   ENDIF
          !ENDIF

          ZTMAX(I) = Z0MAX(I)*EXP( - ((1.0_r8-SIGMAF(I))**2)       &
               *0.8_r8*CA*SQRT(USTAR(I)*0.01_r8/(1.5E-05_r8)))

          !**  test xubin's new z0

          !       ZTMAX(I) = Z0MAX(I)

          IF(mlsi(I).EQ.0_i8) THEN
             RESTAR   = USTAR(I) * Z0MAX(I) / VIS
             RESTAR   = MAX(RESTAR,0.000001_r8)
             !         RESTAR   = LOG(RESTAR)
             !         RESTAR   = MIN(RESTAR,5.)
             !         RESTAR   = MAX(RESTAR,-5.)
             !         RAT(I)   = AA1 + BB1 * RESTAR + CC1 * RESTAR ** 2
             !         RAT(I)   = RAT(I) / (1. + BB2 * RESTAR
             !    &                       + CC2 * RESTAR ** 2)
             !  Rat taken from Zeng, Zhao and Dickinson 1997
             RAT(I)   = 2.67_r8 * restar ** 0.25_r8 - 2.57_r8
             RAT(I)   = MIN(RAT(I),7.0_r8)
             ZTMAX(I) = Z0MAX(I) * EXP(-RAT(I))
          ENDIF
       ENDIF
    ENDDO
    !##DG  IF(LAT.EQ.LATD) THEN
    !##DG    PRINT *, ' z0max, ztmax, restar, RAT(I) =', 
    !##DG &   z0max, ztmax, restar, RAT(I)
    !##DG  ENDIF
    DO I = 1, IM
       IF(flag_iter(i)) THEN 
          DTV(I)   = THV1(I) - TVS(I)
          ADTV     = ABS(DTV(I))
          ADTV     = MAX(ADTV,0.001_r8)
          DTV(I)   = SIGN(1.0_r8,DTV(I)) * ADTV
          RB(I)    = G * DTV(I) * Z1(I) / (0.5_r8 * (THV1(I) + TVS(I))       &
               * WIND(I) * WIND(I))
          RB(I)    = MAX(RB(I),-5000.0_r8)
          FM(I)    = LOG((Z0MAX(I)+Z1(I)) / Z0MAX(I))
          FH(I)    = LOG((ZTMAX(I)+Z1(I)) / ZTMAX(I))
          HLINF(I) = RB(I) * FM(I) * FM(I) / FH(I)
          FM10(I)  = LOG((Z0MAX(I)+10.0_r8) / Z0MAX(I))
          FH2(I)   = LOG((ZTMAX(I)+2.0_r8) / ZTMAX(I))
       ENDIF
    ENDDO
    !##DG  IF(LAT.EQ.LATD) THEN
    !##DG    PRINT *, ' DTV, RB(I), FM(I), FH(I), HLINF =',
    !##DG &   dtv, rb, FM(I), FH(I), hlinf
    !##DG  ENDIF
    !
    !  STABLE CASE
    !
    DO I = 1, IM
       IF(flag_iter(i)) THEN 
          IF(DTV(I).GE.0.0_r8) THEN
             HL1(I) = HLINF(I)
          ENDIF
          IF(DTV(I).GE.0.0_r8.AND.HLINF(I).GT.0.25_r8) THEN
             HL0INF = Z0MAX(I) * HLINF(I) * Z1I(I)
             HLTINF = ZTMAX(I) * HLINF(I) * Z1I(I)
             AA     = SQRT(1.0_r8 + 4.0_r8 * ALPHA * HLINF(I))
             AA0    = SQRT(1.0_r8 + 4.0_r8 * ALPHA * HL0INF)
             BB     = AA
             BB0    = SQRT(1.0_r8 + 4.0_r8 * ALPHA * HLTINF)
             PM(I)  = AA0 - AA + LOG((AA + 1.0_r8) / (AA0 + 1.0_r8))
             PH(I)  = BB0 - BB + LOG((BB + 1.0_r8) / (BB0 + 1.0_r8))
             FMS    = FM(I) - PM(I)
             FHS    = FH(I) - PH(I)
             HL1(I) = FMS * FMS * RB(I) / FHS
          ENDIF
       ENDIF
    ENDDO
    !
    !  SECOND ITERATION
    !
    DO I = 1, IM
       IF(flag_iter(i)) THEN 
          IF(DTV(I).GE.0.0_r8) THEN
             HL0     = Z0MAX(I) * HL1(I) * Z1I(I)
             HLT     = ZTMAX(I) * HL1(I) * Z1I(I)
             AA      = SQRT(1.0_r8 + 4.0_r8 * ALPHA * HL1(I))
             AA0     = SQRT(1.0_r8 + 4.0_r8 * ALPHA * HL0)
             BB      = AA
             BB0     = SQRT(1.0_r8 + 4.0_r8 * ALPHA * HLT)
             PM(I)   = AA0 - AA + LOG((AA + 1.0_r8) / (AA0 + 1.0_r8))
             PH(I)   = BB0 - BB + LOG((BB + 1.0_r8) / (BB0 + 1.0_r8))
             HL110   = HL1(I) * 10.0_r8 * Z1I(I)
             AA      = SQRT(1.0_r8 + 4.0_r8 * ALPHA * HL110)
             PM10(I) = AA0 - AA + LOG((AA + 1.0_r8) / (AA0 + 1.0_r8))
             HL12(I) = HL1(I) * 2.0_r8 * Z1I(I)
             !         AA      = SQRT(1.0_r8 + 4.0_r8 * ALPHA * HL12(I))
             BB      = SQRT(1.0_r8 + 4.0_r8 * ALPHA * HL12(I))
             PH2(I)  = BB0 - BB + LOG((BB + 1.0_r8) / (BB0 + 1.0_r8))
          ENDIF
       ENDIF
    ENDDO
    !!
    !##DG  IF(LAT.EQ.LATD) THEN
    !##DG    PRINT *, ' HL1(I), PM, PH =',
    !##DG &   HL1(I),  pm, ph
    !##DG  ENDIF
    !
    !  UNSTABLE CASE
    !
    !
    !  CHECK FOR UNPHYSICAL OBUKHOV LENGTH
    !
    DO I=1,IM
       IF(flag_iter(i)) THEN 
          IF(DTV(I).LT.0.0_r8) THEN
             OLINF = Z1(I) / HLINF(I)
             IF(ABS(OLINF).LE.50.0_r8 * Z0MAX(I)) THEN
                HLINF(I) = -Z1(I) / (50.0_r8 * Z0MAX(I))
             ENDIF
          ENDIF
       ENDIF
    ENDDO
    !
    !  GET PM AND PH 
    !
    DO I = 1, IM
       IF(flag_iter(i)) THEN 
          IF(DTV(I).LT.0.0_r8.AND.HLINF(I).GE.-0.5_r8) THEN
             HL1(I)  = HLINF(I)
             PM(I)   = (A0 + A1 * HL1(I)) * HL1(I)                          &
                  / (1.0_r8 + B1 * HL1(I) + B2 * HL1(I) * HL1(I))
             PH(I)   = (A0P + A1P * HL1(I)) * HL1(I)                       &
                  / (1.0_r8 + B1P * HL1(I) + B2P * HL1(I) * HL1(I))
             HL110   = HL1(I) * 10.0_r8 * Z1I(I)
             PM10(I) = (A0 + A1 * HL110) * HL110                         &
                  / (1.0_r8 + B1 * HL110 + B2 * HL110 * HL110)
             HL12(I) = HL1(I) * 2.0_r8 * Z1I(I)
             PH2(I)  = (A0P + A1P * HL12(I)) * HL12(I)                  &
                  / (1.0_r8 + B1P * HL12(I) + B2P * HL12(I) * HL12(I))
          ENDIF
          IF(DTV(I).LT.0.AND.HLINF(I).LT.-0.5_r8) THEN
             HL1(I)  = -HLINF(I)
             PM(I)   = LOG(HL1(I)) + 2.0_r8 * HL1(I) ** (-0.25_r8) - 0.8776_r8
             PH(I)   = LOG(HL1(I)) + 0.5_r8 * HL1(I) ** (-0.5_r8) + 1.386_r8
             HL110   = HL1(I) * 10.0_r8 * Z1I(I)
             PM10(I) = LOG(HL110) + 2.0_r8 * HL110 ** (-0.25_r8) - 0.8776_r8
             HL12(I) = HL1(I) * 2.0_r8 * Z1I(I)
             PH2(I)  = LOG(HL12(I)) + 0.5_r8 * HL12(I) ** (-0.5_r8) + 1.386_r8
          ENDIF
       ENDIF
    ENDDO
    !
    !  FINISH THE EXCHANGE COEFFICIENT COMPUTATION TO PROVIDE FM AND FH
    !
    DO I = 1, IM
       IF(flag_iter(i)) THEN 
          FM(I) = FM(I) - PM(I)
          FH(I) = FH(I) - PH(I)
          FM10(I) = FM10(I) - PM10(I)
          FH2(I) = FH2(I) - PH2(I)
          CM(I) = CA * CA / (FM(I) * FM(I))
          CH(I) = CA * CA / (FM(I) * FH(I))
          CQ = CH(I)
          STRESS(I) = CM(I) * WIND(I) * WIND(I)
          !USTAR(I)  = SQRT(STRESS(I))
          !       USTAR(I) = SQRT(CM(I) * WIND(I) * WIND(I))
          STRESS(I) = (USTAR(I)**2)  

       ENDIF
    ENDDO
    !##DG  IF(LAT.EQ.LATD) THEN
    !##DG    PRINT *, ' FM, FH, CM, CH(I), USTAR =',
    !##DG &   FM, FH, CM, ch, USTAR
    !##DG  ENDIF
    !
    !  UPDATE Z0 OVER OCEAN
    !
    DO I = 1, IM
       IF(flag_iter(i)) THEN 
          IF(mlsi(I).EQ.0_i8) THEN
             !          Z0(I) = (CHARNOCK / G) * USTAR(I) ** 2
             Z0(I) = (CHARNOCK / G) * USTAR(I) * USTAR(I)
             ! MBEK -- TOGA-COARE FLUX ALGORITHM
             !         Z0(I) = (CHARNOCK / G) * USTAR(I)*USTAR(I) +  ARNU/USTAR(I)
             !  NEW IMPLEMENTATION OF Z0
             !         CC = USTAR(I) * Z0 / RNU
             !         PP = CC / (1. + CC)
             !         FF = G * ARNU / (CHARNOCK * USTAR(I) ** 3)
             !         Z0 = ARNU / (USTAR(I) * FF ** PP)
             Z0(I) = MIN(Z0(I),0.1_r8)
             Z0(I) = MAX(Z0(I),1.E-7_r8)
             Z0RL(I) = 100.0_r8 * Z0(I)
          ENDIF
       ENDIF
    ENDDO

    RETURN
  END SUBROUTINE SFC_DIFF


END MODULE DiffCoef





MODULE Pbl_UniversityWashington
  USE PhysicalFunctions,ONLY : fpvs2es5
  USE Parallelism, ONLY : myid
  USE DiffCoef,ONLY : moninq

    IMPLICIT NONE
  SAVE

  PRIVATE 
  ! vertical_diffusion_tend_
  !                         |
  !                         |____compute_tms
  !                         |
  !                         |____compute_eddy_diff_
  !                         |                      |
  !                         |                      |__trbintd_
  !                         |                      |          |
  !                         |                      |          |__sfdiag
  !                         |                      |
  !                         |                      |__austausch_atm
  !                         |                      |
  !                         |                      |__caleddy_
  !                         |                      |          |
  !                         |                      |          |__exacol
  !                         |                      |          |
  !                         |                      |          |__zisocl
  !                         |                      |
  !                         |                      |__compute_vdiff_
  !                         |                                       |
  !                         |                                       |__compute_molec_diff_
  !                         |                                       |                     |
  !                         |                                       |                     |__ubc_get_vals
  !                         |                                       |
  !                         |                                       |__vd_lu_decomp
  !                         |                                       |
  !                         |                                       |__vd_lu_solve
  !                         |                                       |
  !                         |                                       |__vd_lu_solve
  !                         |                                       |
  !                         |                                       |__vd_lu_decomp
  !                         |                                       |
  !                         |                                       |__vd_lu_solve
  !                         |                                       |
  !                         |                                       |__vd_lu_decomp
  !                         |                                       |
  !                         |                                       |__vd_lu_qdecomp
  !                         |
  !                         |____aqsat
  !                         |
  !                         |__compute_vdiff_
  !                         |                |
  !                         |                |__compute_molec_diff_
  !                         |                |                     |
  !                         |                |                     |__ubc_get_vals
  !                         |                |
  !                         |                |__vd_lu_decomp
  !                         |                |
  !                         |                |__vd_lu_solve
  !                         |                |
  !                         |                |__vd_lu_solve
  !                         |                |
  !                         |                |__vd_lu_decomp
  !                         |                |
  !                         |                |__vd_lu_solve
  !                         |                |
  !                         |                |__vd_lu_decomp
  !                         |                |
  !                         |                |__vd_lu_qdecomp
  !                         |
  !                         |__compute_vdiff_
  !                         |                |
  !                         |                |__compute_molec_diff_
  !                         |                |                     |
  !                         |                |                     |__ubc_get_vals
  !                         |                |
  !                         |                |__vd_lu_decomp
  !                         |                |
  !                         |                |__vd_lu_solve
  !                         |                |
  !                         |                |__vd_lu_solve
  !                         |                |
  !                         |                |__vd_lu_decomp
  !                         |                |
  !                         |                |__vd_lu_solve
  !                         |                |
  !                         |                |__vd_lu_decomp
  !                         |                |
  !                         |                |__vd_lu_qdecomp
  !                         |
  !                         |____positive_moisture
  !                         |
  !                         |____aqsat
  !
  !
  !
  !
  !

  !----------------------------------------------------------------------------
  ! precision/kind constants add data public
  !----------------------------------------------------------------------------
  INTEGER,PARAMETER :: R8 = SELECTED_REAL_KIND(12) ! 8 byte real
  INTEGER,PARAMETER :: R4 = SELECTED_REAL_KIND( 6) ! 4 byte real
  INTEGER,PARAMETER :: RN = KIND(1.0)              ! native real
  INTEGER,PARAMETER :: I8 = SELECTED_INT_KIND (13) ! 8 byte integer
  INTEGER,PARAMETER :: I4 = SELECTED_INT_KIND ( 6) ! 4 byte integer
  INTEGER,PARAMETER :: IN = KIND(1)                ! native integer
  INTEGER,PARAMETER :: CS = 80                     ! short char
  INTEGER,PARAMETER :: CL = 256                    ! long char
  INTEGER,PARAMETER :: CX = 512                    ! extra-long char

  CHARACTER(len=8), PARAMETER :: eddy_scheme='diag_TKE'
  INTEGER         , PARAMETER :: iulog=0
  LOGICAL                     :: do_iss=.TRUE.                            ! switch for implicit turbulent surface stress
  LOGICAL                     :: do_molec_diff = .FALSE.      ! Switch for molecular diffusion
  LOGICAL                     :: do_tms= .TRUE.                       ! Switch for turbulent mountain stress
  REAL(r8)                    :: tms_orocnst =1                 ! Converts from standard deviation to height
  LOGICAL                     :: do_pseudocon_diff = .FALSE.  ! If .true., do pseudo-conservative variables diffusion
  LOGICAL                     :: MODAL_AERO= .FALSE.     
  LOGICAL                     :: wstarent=.TRUE.                  ! .true. means use the 'wstar' entrainment closure. 
  INTEGER                     :: nturb =3                   ! Number of iteration steps for calculating eddy diffusivity [ # ]
  LOGICAL                     :: is_first_step = .TRUE.
  LOGICAL,          PARAMETER :: use_kvf        =  .FALSE.      ! .true. (.false.) : initialize kvh/kvm =  kvf ( 0. )
  LOGICAL,          PARAMETER :: use_dw_surf    =  .FALSE.       ! Used in 'zisocl'. Default is 'true'
  ! If 'true', surface interfacial energy does not contribute to the CL mean
  !            stbility functions after finishing merging.     For this case,
  !           'dl2n2_surf' is only used for a merging test based on 'l2n2'
  ! If 'false',surface interfacial enery explicitly contribute to    CL mean
  !            stability functions after finishing merging.    For this case,
  !           'dl2n2_surf' and 'dl2s2_surf' are directly used for calculating
  !            surface interfacial layer energetics

  ! --------------------------------- !
  ! PBL Parameters used in the UW PBL !
  ! --------------------------------- !

  CHARACTER,        PARAMETER :: sftype         ='z'  !'z'  !'l'           !d,l,u,z: Method for calculating saturation fraction

  CHARACTER(len=4), PARAMETER :: choice_evhc    = 'maxi'        ! 'orig',   'ramp',   'maxi'   : recommended to be used with choice_radf 
  CHARACTER(len=6), PARAMETER :: choice_radf    = 'maxi'        ! 'orig',   'ramp',   'maxi'   : recommended to be used with choice_evhc 
  CHARACTER(len=6), PARAMETER :: choice_SRCL    = 'nonamb'      ! 'origin', 'remove', 'nonamb'

  CHARACTER(len=6), PARAMETER :: choice_tunl    = 'rampsl'       !'rampcl'      ! 'origin', 'rampsl'(Sungsu), 'rampcl'(Chris)
  REAL(r8),         PARAMETER :: ctunl          =  2._r8        !  Maximum asympt leng = ctunl*tunl when choice_tunl = 'rampsl(cl)' [ no unit ]
  CHARACTER(len=6), PARAMETER :: choice_leng    = 'origin'      ! 'origin', 'takemn'
  REAL(r8),         PARAMETER :: cleng          =  3._r8        !  Order of 'leng' when choice_leng = 'origin' [ no unit ]
  CHARACTER(len=6), PARAMETER :: choice_tkes    = 'ibprod'      ! 'ibprod' (include tkes in computing bprod), 'ebprod'(exclude)

  ! Parameters for 'sedimenttaion-entrainment feedback' for liquid stratus 
  ! If .false.,  no sedimentation entrainment feedback ( i.e., use default evhc )
  LOGICAL,          PARAMETER :: ens=.TRUE.
  LOGICAL,          PARAMETER :: id_sedfact     = .TRUE.!.FALSE.
  REAL(r8),         PARAMETER :: ased           =  9._r8        !  Valid only when id_sedfact = .true.
  ! --------------------------------------------------------------------------------------------------- !
  ! Parameters governing entrainment efficiency A = a1l(i)*evhc, evhc = 1 + a2l * a3l * L * ql / jt2slv !
  ! Here, 'ql' is cloud-top LWC and 'jt2slv' is the jump in 'slv' across                                !
  ! the cloud-top entrainment zone ( across two grid layers to consider full mixture )                  !
  ! --------------------------------------------------------------------------------------------------- !

  REAL(r8),         PARAMETER :: a1l            =   0.10_r8     ! Dry entrainment efficiency for TKE closure
  ! a1l = 0.2*tunl*erat^-1.5, where erat = <e>/wstar^2 for dry CBL =  0.3.

  REAL(r8),         PARAMETER :: a1i            =   0.2_r8      ! Dry entrainment efficiency for wstar closure
  REAL(r8),         PARAMETER :: ccrit          =   0.5_r8      ! Minimum allowable sqrt(tke)/wstar. Used in solving cubic equation for 'ebrk'
  REAL(r8),         PARAMETER :: wstar3factcrit =   0.5_r8      ! 1/wstar3factcrit is the maximally allowed enhancement of 'wstar3' due to entrainment.

  REAL(r8),         PARAMETER :: a2l            =   30._r8      ! Moist entrainment enhancement param (recommended range : 10~30 )
  REAL(r8),         PARAMETER :: a3l            =   0.8_r8      ! Approximation to a complicated thermodynamic parameters

  REAL(r8),         PARAMETER :: jbumin         =   .001_r8     ! Minimum buoyancy jump at an entrainment jump, [m/s2]
  REAL(r8),         PARAMETER :: evhcmax        =   10._r8      ! Upper limit of evaporative enhancement factor

  REAL(r8),         PARAMETER :: ustar_min      =   0.01_r8     ! Minimum permitted value of ustar [ m/s ] 
  REAL(r8),         PARAMETER :: onet           =   1._r8/3._r8 ! 1/3 power in wind gradient expression [ no unit ]
  REAL(r8),         PARAMETER :: qmin(3)        =   1.0e-21_r8    ! Minimum grid-mean LWC counted as clouds [kg/kg]
  REAL(r8),         PARAMETER :: ntzero         =   1.e-12_r8   ! Not zero (small positive number used in 's2')
  REAL(r8),         PARAMETER :: b1             =   5.8_r8      ! TKE dissipation D = e^3/(b1*leng), e = b1*W.
  REAL(r8)                    :: b123                           ! b1**(2/3)
  REAL(r8),         PARAMETER :: tunl           =   0.085_r8    ! Asympt leng = tunl*(turb lay depth)
  REAL(r8),         PARAMETER :: alph1          =   0.5562_r8   ! alph1~alph5 : Galperin instability function parameters
  REAL(r8),         PARAMETER :: alph2          =  -4.3640_r8   !               These coefficients are used to calculate 
  REAL(r8),         PARAMETER :: alph3          = -34.6764_r8   !               'sh' and 'sm' from 'gh'.
  REAL(r8),         PARAMETER :: alph4          =  -6.1272_r8   !
  REAL(r8),         PARAMETER :: alph5          =   0.6986_r8   !
  REAL(r8),         PARAMETER :: ricrit         =   0.19_r8     ! Critical Richardson number for turbulence. Can be any value >= 0.19.
  REAL(r8),         PARAMETER :: ae             =   1._r8       ! TKE transport efficiency [no unit]
  REAL(r8),         PARAMETER :: rinc           =  -0.04_r8     ! Minimum W/<W> used for CL merging test 
  REAL(r8),         PARAMETER :: wpertmin       =   1.e-6_r8    ! Minimum PBL eddy vertical velocity perturbation
  REAL(r8),         PARAMETER :: wfac           =   1._r8       ! Ratio of 'wpert' to sqrt(tke) for CL.
  REAL(r8),         PARAMETER :: tfac           =   1._r8       ! Ratio of 'tpert' to (w't')/wpert for CL. Same ratio also used for q
  REAL(r8),         PARAMETER :: fak            =   8.5_r8      ! Constant in surface temperature excess for stable STL. [ no unit ]         
  REAL(r8),         PARAMETER :: rcapmin        =   0.1_r8      ! Minimum allowable e/<e> in a CL
  REAL(r8),         PARAMETER :: rcapmax        =   2.0_r8      ! Maximum allowable e/<e> in a CL
  REAL(r8),         PARAMETER :: tkemax         =  20._r8       ! TKE is capped at tkemax [m2/s2]
  REAL(r8),         PARAMETER :: tkemin         =  1e-6_r8       ! TKE is capped at tkemin [m2/s2]
  REAL(r8),         PARAMETER :: lambda         =   0.5_r8      ! Under-relaxation factor ( 0 < lambda =< 1 )
  REAL(r8),         PARAMETER :: lambdaEns         =   0.75_r8      ! Under-relaxation factor ( 0 < lambda =< 1 )


  LOGICAL,          PARAMETER :: set_qrlzero    =  .FALSE.      ! .true. ( .false.) : turning-off ( on) radiative-turbulence interaction by setting qrl = 0.

  ! ------------------------------------- !
  ! PBL Parameters not used in the UW PBL !
  ! ------------------------------------- !

  REAL(r8),         PARAMETER :: pblmaxp        =  4.e4_r8      ! PBL max depth in pressure units. 
  REAL(r8),         PARAMETER :: zkmin          =  0.01_r8      ! Minimum kneutral*f(ri). 
  REAL(r8),         PARAMETER :: betam          = 15.0_r8       ! Constant in wind gradient expression.
  REAL(r8),         PARAMETER :: betas          =  5.0_r8       ! Constant in surface layer gradient expression.
  REAL(r8),         PARAMETER :: betah          = 15.0_r8       ! Constant in temperature gradient expression.
  REAL(r8),         PARAMETER :: fakn           =  7.2_r8       ! Constant in turbulent prandtl number.
  REAL(r8),         PARAMETER :: ricr           =  0.3_r8       ! Critical richardson number.
  REAL(r8),         PARAMETER :: sffrac         =  0.1_r8       ! Surface layer fraction of boundary layer
  REAL(r8),         PARAMETER :: binm           =  betam*sffrac ! betam * sffrac
  REAL(r8),         PARAMETER :: binh           =  betah*sffrac ! betah * sffrac

  ! ------------------------------------------------------- !
  ! PBL constants set using values from other parts of code !
  ! ------------------------------------------------------- !
  REAL(R8),PARAMETER :: SHR_CONST_G       = 9.80616_R8      ! acceleration of gravity ~ m/s^2
  REAL(r8),PARAMETER :: gravit            = shr_const_g     ! gravitational acceleration (m/s**2)
  REAL(R8),PARAMETER :: SHR_CONST_KARMAN  = 0.4_R8          ! Von Karman constant
  REAL(r8),PARAMETER :: karman            = shr_const_karman     ! Von Karman constant

  REAL(R8),PARAMETER :: SHR_CONST_CPDAIR  = 1.00464e3_R8    ! specific heat of dry air   ~ J/kg/K
  REAL(r8),PARAMETER :: cpair             = shr_const_cpdair     ! specific heat of dry air (J/K/kg)
  REAL(R8),PARAMETER :: SHR_CONST_MWDAIR  = 28.966_R8       ! molecular weight dry air ~ kg/kmole
  REAL(R8),PARAMETER :: SHR_CONST_MWWV    = 18.016_R8       ! molecular weight water vapor
  REAL(R8),PARAMETER :: SHR_CONST_BOLTZ   = 1.38065e-23_R8  ! Boltzmann's constant ~ J/K/molecule
  REAL(R8),PARAMETER :: SHR_CONST_AVOGAD  = 6.02214e26_R8   ! Avogadro's number ~ molecules/kmole
  REAL(R8),PARAMETER :: SHR_CONST_RGAS    = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ       ! Universal gas constant ~ J/K/kmole
  REAL(R8),PARAMETER :: SHR_CONST_RDAIR   = SHR_CONST_RGAS/SHR_CONST_MWDAIR        ! Dry air gas constant     ~ J/K/kg
  REAL(r8),PARAMETER :: rair              = shr_const_rdair   ! Dry air gas constant     ~ J/K/kg
  REAL(r8),PARAMETER :: avogad            = shr_const_avogad     ! Avogadro's number (molecules/kmole)
  REAL(r8),PARAMETER :: boltz             = shr_const_boltz      ! Boltzman's constant (J/K/molecule)

  REAL(R8),PARAMETER :: SHR_CONST_RWV     = SHR_CONST_RGAS/SHR_CONST_MWWV          ! Water vapor gas constant ~ J/K/kg
  REAL(r8),PARAMETER :: zvir         = (shr_const_rwv/shr_const_rdair)-1.0_R8 ! (rh2o/rair) - 1

  REAL(R8),PARAMETER :: SHR_CONST_LATVAP  = 2.501e6_R8      ! latent heat of evaporation ~ J/kg
  REAL(r8),PARAMETER :: latvap      = shr_const_latvap     ! Latent heat of vaporization (J/kg)

  REAL(R8),PARAMETER :: SHR_CONST_LATICE  = 3.337e5_R8      ! latent heat of fusion      ~ J/kg
  REAL(r8),PARAMETER :: latice      = shr_const_latice     ! Latent heat of fusion (J/kg)

  REAL(R8),PARAMETER :: SHR_CONST_LATSUB  = SHR_CONST_LATICE + SHR_CONST_LATVAP     ! latent heat of sublimation ~ J/kg

  REAL(r8),PARAMETER :: latsub =SHR_CONST_LATSUB                        ! Latent heat of sublimation

  REAL(r8),PARAMETER :: g   =SHR_CONST_G                           ! Gravitational acceleration
  REAL(r8),PARAMETER :: vk  =SHR_CONST_KARMAN                           ! Von Karman's constant
  REAL(r8),PARAMETER :: ccon  =fak*sffrac*SHR_CONST_KARMAN                         ! fak * sffrac * vk
  REAL(r8),PARAMETER :: tmin  = 173.16_r8      ! min temperature (K) for table
  REAL(r8),PARAMETER :: tmax  = 375.16_r8      ! max temperature (K) for table! Maximum temperature entry in table
  REAL(r8),PARAMETER :: trice =  20.00_r8  ! Transition range from es over range to es over ice
  REAL(r8),PARAMETER :: ttrice=trice
  REAL(R8),PARAMETER :: SHR_CONST_TKFRZ   = 273.15_R8       ! freezing T of fresh water          ~ K 
  REAL(r8), PARAMETER :: tmelt       = shr_const_tkfrz      ! Freezing point of water (K)
  REAL(r8), PARAMETER :: mwdry        = shr_const_mwdair! molecular weight dry air

  INTEGER                    :: ntop_turb                      ! Top interface level to which turbulent vertical diffusion is applied ( = 1 )
  INTEGER                    :: nbot_turb                      ! Bottom interface level to which turbulent vertical diff is applied ( = pver )
  INTEGER                    :: ntop_eddy   ! Top    interface level to which eddy vertical diffusion is applied ( = 1 )
  INTEGER                    :: nbot_eddy   ! Bottom interface level to which eddy vertical diffusion is applied ( = pver )
  REAL(r8), PARAMETER   :: d0     = 1.52E20_r8         ! Diffusion factor [ m-1 s-1 ] molec sqrt(kg/kmol/K) [ unit ? ]

  REAL(r8)     , ALLOCATABLE :: ml2(:)                         ! Mixing lengths squared. Not used in the UW PBL. Used for computing free air diffusivity.
  REAL(KIND=r8), ALLOCATABLE :: hypm(:)                        ! reference pressures at midpoints
  CHARACTER*3  , ALLOCATABLE :: cnst_type(:)                   ! wet or dry mixing ratio

  ! Parameters used for Turbulent Mountain Stress

  !  real(r8), parameter :: z0fac   = 0.025_r8              ! Factor determining z_0 from orographic standard deviation
  REAL(r8), PARAMETER :: z0fac  = 0.075_r8    ! Factor determining z_0 from orographic standard deviation [ no unit ]
  !  real(r8), parameter :: z0max   = 100._r8               ! Max value of z_0 for orography
  REAL(r8), PARAMETER :: z0max  = 100._r8     ! Maximum value of z_0 for orography [ m ]

  !  real(r8), parameter :: horomin = 10._r8                ! Min value of subgrid orographic height for mountain stress
  REAL(r8), PARAMETER :: horomin= 1._r8       ! Minimum value of subgrid orographic height for mountain stress [ m ]

  REAL(r8), PARAMETER :: dv2min = 0.01_r8     ! Minimum shear squared [ m2/s2 ]
  !  real(r8), parameter :: dv2min  = 0.01_r8               ! Minimum shear squared
  REAL(r8)            :: oroconst             ! Converts from standard deviation to height [ no unit ]
  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !


  INTEGER            :: ncvmax =   -1        ! Max numbers of CLs (good to set to 'pver')
  INTEGER            :: pver   =   -1      

  INTEGER            :: nbot_molec                        ! Bottom level where molecular diffusivity is applied
  INTEGER            :: ntop_molec  ! Top    interface level to which molecular vertical diffusion is applied ( = 1 )

  INTEGER, PUBLIC    :: pcnst  = -1                      ! number of advected constituents (including water vapor)
  INTEGER            :: ntop                         ! Top interface level to which vertical diffusion is applied ( = 1 ).
  INTEGER            :: nbot                         ! Bottom interface level to which vertical diffusion is applied ( = pver ).

  ! Below stores logical array of fields to be diffused

  TYPE vdiff_selector 
     PRIVATE
     LOGICAL, POINTER, DIMENSION(:) :: fields
  END TYPE vdiff_selector

  ! Below extends .not. to operate on type vdiff_selector

  INTERFACE OPERATOR(.NOT.)
     MODULE PROCEDURE not
  END INTERFACE

  ! Below provides functionality of intrinsic any for type vdiff_selector

  INTERFACE any                           
     MODULE PROCEDURE my_any
  END INTERFACE

  TYPE(vdiff_selector) :: fieldlist_wet                ! Logical switches for moist mixing ratio diffusion
  TYPE(vdiff_selector) :: fieldlist_dry                ! Logical switches for dry mixing ratio diffusion

  INTEGER              :: ixcldice, ixcldliq           ! Constituent indices for cloud liquid and ice water
  INTEGER              :: ixnumice, ixnumliq
  CHARACTER(len=528), ALLOCATABLE :: vdiffnam(:)      ! Names of vertical diffusion tendencies
  CHARACTER(len=16)           , ALLOCATABLE :: cnst_name(:)     ! constituent names

  ! Constants for each tracer
  REAL(r8),  ALLOCATABLE   :: qmincg   (:)          ! for backward compatibility only
  REAL(r8),    ALLOCATABLE :: cnst_mw  (:)          ! molecular weight (kg/kmole)
  LOGICAL,  ALLOCATABLE    :: cnst_fixed_ubc(:) != .false.  ! upper bndy condition = fixed ?
  REAL(r8)     , ALLOCATABLE :: mw_fac(:)                      ! sqrt(1/M_q + 1/M_d) in constituent diffusivity [  unit ? ]

  !
  ! Data
  !
  INTEGER, PARAMETER:: plenest=250  ! length of saturation vapor pressure table
  !
  ! Table of saturation vapor pressure values es from tmin degrees
  ! to tmax+1 degrees k in one degree increments.  ttrice defines the
  ! transition region where es is a combination of ice & water values
  !
  REAL(r8) estbl(plenest)      ! table values of saturation vapor pressure
  LOGICAL,PARAMETER :: icephs   = .TRUE. ! false => saturation vapor press over water only
  ! Ice phase (true or false)
  REAL(r8) pcf(6)     ! polynomial coeffs -> es transition water to ice
  REAL(r8),PARAMETER :: epsilo       = shr_const_mwwv/shr_const_mwdair   ! ratio of h2o to dry air molecular weights 
  REAL(r8),PARAMETER :: epsqs=epsilo

  PUBLIC :: Init_Pbl_UniversityWashington
  PUBLIC :: Finalize_Pbl_UniversityWashington
  PUBLIC :: vertical_diffusion_tend
CONTAINS
  !
  !  Init_Pbl_UniversityWashington
  !
  SUBROUTINE Init_Pbl_UniversityWashington(pver_in,pcnst_in,ncnst,ILCON,a_hybr,b_hybr,RESTART)
    IMPLICIT NONE
    INTEGER , INTENT(IN   ) :: pver_in
    INTEGER , INTENT(IN   ) :: pcnst_in
    INTEGER , INTENT(in   )   :: ncnst          ! Number of constituents

    CHARACTER(LEN=*), INTENT(in  )  :: ILCON

    REAL(KIND=r8),    INTENT(IN) :: a_hybr(pver_in+1)
    REAL(KIND=r8),    INTENT(IN) :: b_hybr(pver_in+1)
    LOGICAL         , INTENT(IN   ) ::  RESTART

    REAL(KIND=r8)   :: ps0
    INTEGER :: k
    !    IF(TRIM(ILCON).EQ.'LSC' .OR. TRIM(ILCON).EQ.'YES' ) THEN 
    !---------------------------Local variables-----------------------------
    !
    REAL(r8) :: t          ! Temperature
    INTEGER  :: n          ! Increment counter
    INTEGER  :: lentbl     ! Calculated length of lookup table
    INTEGER  :: itype      ! Ice phase: 0 -> no ice phase
    !                        1 -> ice phase, no transition
    !                       -x -> ice phase, x degree transition
    LOGICAL  :: ip         ! Ice phase logical flag
    !
    !-----------------------------------------------------------------------
    !
    IF(RESTART) is_first_step = .FALSE.

    pcnst  = pcnst_in
    pver   = pver_in
    ncvmax =   pver ! Max numbers of CLs (good to set to 'pver')     
                    ! hypm     reference state midpoint pressures
    ALLOCATE(hypm     (pver_in));hypm=0.0_r8
    ALLOCATE(vdiffnam (pcnst)  );vdiffnam=''
    ALLOCATE(cnst_name(pcnst) );cnst_name=''
    ALLOCATE(qmincg   (pcnst) );qmincg=0.0_r8
    ALLOCATE(cnst_mw  (pcnst));cnst_mw=0.0_r8
    ALLOCATE(cnst_type(pcnst)  );cnst_type='wet'
    ALLOCATE(cnst_fixed_ubc(pcnst));cnst_fixed_ubc(1:pcnst) = .FALSE. 
    ALLOCATE(mw_fac(pcnst));mw_fac=0.0_r8

    ps0    = 1.0e5_r8            ! Base state surface pressure (pascals)
    DO k=pver_in,1,-1
!      hypm(k) =  ps0*sig(pver_in + 1 - k)
!      SB  already top to bottom 
       hypm(k) =  0.5_r8 * ( ps0 * ( b_hybr(k)+b_hybr(k+1) ) + &
                             a_hybr(k) + a_hybr(k+1) )
      !if(myid.eq.0) write(0,*) 'hypm(k) top-down',k,hypm(k)

    END DO

    qmincg=1.0e-21_r8
    cnst_mw=18.0_r8 ! Molecular weight [ kg/kmole ]
    ! ---------------------------------- !
    ! Initialize diffusion solver module !
    ! ---------------------------------- !
    CALL vertical_diffusion_init(ILCON)

    lentbl = INT(tmax-tmin+2.000001_r8)
    IF (lentbl .GT. plenest) THEN
       WRITE(0,9000) tmax, tmin, plenest
       STOP 'call endrun (GESTBL)    ! Abnormal termination'
    END IF
    !
    ! Begin building es table.
    ! Check whether ice phase requested.
    ! If so, set appropriate transition range for temperature
    !
    IF (icephs) THEN
       IF (ttrice /= 0.0_r8) THEN
          itype = -ttrice
       ELSE
          itype = 1
       END IF
    ELSE
       itype = 0
    END IF
    !
    t = tmin - 1.0_r8
    DO n=1,lentbl
       t = t + 1.0_r8
       CALL gffgch(t,estbl(n),itype)
    END DO
    !
    DO n=lentbl+1,plenest
       estbl(n) = -99999.0_r8
    END DO
    !
    ! Table complete -- Set coefficients for polynomial approximation of
    ! difference between saturation vapor press over water and saturation
    ! pressure over ice for -ttrice < t < 0 (degrees C). NOTE: polynomial
    ! is valid in the range -40 < t < 0 (degrees C).
    !
    !                  --- Degree 5 approximation ---
    !
    pcf(1) =  5.04469588506e-01_r8
    pcf(2) = -5.47288442819e+00_r8
    pcf(3) = -3.67471858735e-01_r8
    pcf(4) = -8.95963532403e-03_r8
    pcf(5) = -7.78053686625e-05_r8
    !
    !                  --- Degree 6 approximation ---
    !
    !-----pcf(1) =  7.63285250063e-02
    !-----pcf(2) = -5.86048427932e+00
    !-----pcf(3) = -4.38660831780e-01
    !-----pcf(4) = -1.37898276415e-02
    !-----pcf(5) = -2.14444472424e-04
    !-----pcf(6) = -1.36639103771e-06
    !
    RETURN
    !
9000 FORMAT('GESTBL: FATAL ERROR *********************************',/, &
            ' TMAX AND TMIN REQUIRE A LARGER DIMENSION ON THE LENGTH', &
            ' OF THE SATURATION VAPOR PRESSURE TABLE ESTBL(PLENEST)',/, &
            ' TMAX, TMIN, AND PLENEST => ', 2f7.2, i3)
    !
  END SUBROUTINE Init_Pbl_UniversityWashington



  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  SUBROUTINE vertical_diffusion_init(ILCON)
    CHARACTER(LEN=*), INTENT(IN   ) :: ILCON
    CHARACTER(128) :: errstring   ! Error status for init_vdiff
    INTEGER        :: k           ! Vertical loop index
    INTEGER        :: m
    INTEGER        :: l

    ! ----------------------------------------------------------------- !
    ! Get indices of cloud liquid and ice within the constituents array !
    ! ----------------------------------------------------------------- !
    IF(TRIM(ILCON).EQ.'LSC' .OR. TRIM(ILCON).EQ.'YES' ) THEN
        ixcldliq=1   !call cnst_get_ind( 'CLDLIQ', ixcldliq )
        ixcldice=1   !call cnst_get_ind( 'CLDICE', ixcldice )
                 !    if( microp_scheme .eq. 'MG' ) then
       ixnumliq=1   !    call cnst_get_ind( 'NUMLIQ', ixnumliq )
       ixnumice=1   !    call cnst_get_ind( 'NUMICE', ixnumi
    ELSE

       ixcldliq=3   !call cnst_get_ind( 'CLDLIQ', ixcldliq )
       ixcldice=2   !call cnst_get_ind( 'CLDICE', ixcldice )
                 !    if( microp_scheme .eq. 'MG' ) then
       ixnumliq=1   !    call cnst_get_ind( 'NUMLIQ', ixnumliq )
       ixnumice=1   !    call cnst_get_ind( 'NUMICE', ixnumice )
                 !    endif
    END IF
    !    if (masterproc) then
    !       write(iulog,*)'Initializing vertical diffusion (vertical_diffusion_init)'
    !    end if


    ! ---------------------------------------------------------------------------------------- !
    ! Initialize molecular diffusivity module                                                  !
    ! Molecular diffusion turned on above ~60 km (50 Pa) if model top is above ~90 km (.1 Pa). !
    ! Note that computing molecular diffusivities is a trivial expense, but constituent        !
    ! diffusivities depend on their molecular weights. Decomposing the diffusion matric        !
    ! for each constituent is a needless expense unless the diffusivity is significant.        !
    ! ---------------------------------------------------------------------------------------- !

    ntop_molec = 1       ! Should always be 1
    nbot_molec = 0       ! Should be set below about 70 km
    IF( hypm(1) .LT. 0.1_r8 ) THEN
       do_molec_diff = .TRUE.
       DO k = 1, pver
          IF( hypm(k) .LT. 50._r8 ) nbot_molec = k
       END DO
       CALL init_molec_diff( r8, pcnst, rair, ntop_molec, nbot_molec, mwdry, &
                              avogad, gravit, cpair, boltz )
       !call addfld( 'TTPXMLC', 'K/S', 1, 'A', 'Top interf. temp. flux: molec. viscosity', phys_decomp )
       !call add_default ( 'TTPXMLC', 1, ' ' )
       !if( masterproc ) write(iulog,fmt='(a,i3,5x,a,i3)') 'NTOP_MOLEC =', ntop_molec, 'NBOT_MOLEC =', nbot_molec
    END IF
    ! ---------------------------------- !    
    ! Initialize eddy diffusivity module !
    ! ---------------------------------- !

    ntop_eddy  = MAX(MIN(1,nbot_molec ),1)      ! No reason not to make this 1, if > 1, must be <= nbot_molec
    nbot_eddy  = pver    ! Should always be pver
    !if( masterproc ) write(iulog,fmt='(a,i3,5x,a,i3)') 'NTOP_EDDY  =', ntop_eddy, 'NBOT_EDDY  =', nbot_eddy

    SELECT CASE ( eddy_scheme )
    CASE ( 'diag_TKE' ) 
       !if( masterproc ) write(iulog,*) 'vertical_diffusion_init: eddy_diffusivity scheme: UW Moist Turbulence Scheme by Bretherton and Park'
       ! Check compatibility of eddy and shallow scheme
       !if( shallow_scheme .ne. 'UW' ) then
       !    write(iulog,*) 'ERROR: shallow convection scheme ', shallow_scheme,' is incompatible with eddy scheme ', eddy_scheme
       !    call endrun( 'convect_shallow_init: shallow_scheme and eddy_scheme are incompatible' )
       !endif
       CALL init_eddy_diff( r8, pver, gravit, cpair, rair, zvir, latvap, latice, &
            ntop_eddy, nbot_eddy, karman )
       !if( masterproc ) write(iulog,*) 'vertical_diffusion: nturb, ntop_eddy, nbot_eddy ', nturb, ntop_eddy, nbot_eddy
    CASE ( 'HB', 'HBR' )
       !if( masterproc ) write(iulog,*) 'vertical_diffusion_init: eddy_diffusivity scheme:  Holtslag and Boville'
       !call init_hb_diff( gravit, cpair, rair, zvir, ntop_eddy, nbot_eddy, karman, eddy_scheme )
    END SELECT


    ! The vertical diffusion solver must operate 
    ! over the full range of molecular and eddy diffusion

    ntop = MIN(ntop_molec,ntop_eddy)
    nbot = MAX(nbot_molec,nbot_eddy)

    ! ------------------------------------------- !
    ! Initialize turbulent mountain stress module !
    ! ------------------------------------------- !

    IF( do_tms ) THEN
       CALL init_tms( r8, tms_orocnst, karman, gravit, rair )
       !call addfld( 'TAUTMSX' ,'N/m2  ',  1,  'A',  'Zonal      turbulent mountain surface stress',  phys_decomp )
       !call addfld( 'TAUTMSY' ,'N/m2  ',  1,  'A',  'Meridional turbulent mountain surface stress',  phys_decomp )
       !call add_default( 'TAUTMSX ', 1, ' ' )
       !call add_default( 'TAUTMSY ', 1, ' ' )
       !if (masterproc) then
       !   write(iulog,*)'Using turbulent mountain stress module'
       !   write(iulog,*)'  tms_orocnst = ',tms_orocnst
       !end if
    ENDIF


    ! ---------------------------------- !
    ! Initialize diffusion solver module !
    ! ---------------------------------- !

    CALL init_vdiff( r8, pcnst, rair, gravit, fieldlist_wet, fieldlist_dry, errstring )
    IF( errstring .NE. '' ) STOP 'call endrun( errstring )'

    ! Use fieldlist_wet to select the fields which will be diffused using moist mixing ratios ( all by default )
    ! Use fieldlist_dry to select the fields which will be diffused using dry   mixing ratios.

    IF( vdiff_select( fieldlist_wet, 'u' ) .NE. '' ) STOP '!call endrun( vdiff_select( fieldlist_wet, u ) )'
    IF( vdiff_select( fieldlist_wet, 'v' ) .NE. '' ) STOP '!call endrun( vdiff_select( fieldlist_wet, v ) )'
    IF( vdiff_select( fieldlist_wet, 's' ) .NE. '' ) STOP '!call endrun( vdiff_select( fieldlist_wet, s ) )'


    DO  k = 1, pcnst
       IF (MODAL_AERO)THEN
          ! Do not diffuse droplet number - treated in dropmixnuc
          !PK if( k == ixndrop ) go to 20 
          ! Don't diffuse aerosol - treated in dropmixnuc
          !PK    do m = 1, ntot_amode
          !PK       if( k == numptr_amode(m)   ) go to 20
          !PK       do l = 1, nspec_amode(m)
          !PK         if( k == lmassptr_amode(l,m)   ) go to 20
          !PK       enddo
          !PK    enddo
       ENDIF
       IF( cnst_get_type_byind(k,pcnst) .EQ. 'wet' ) THEN
          IF( vdiff_select( fieldlist_wet, 'q', k ) .NE. '' ) STOP 'call endrun( vdiff_select( fieldlist_wet, q, k ) )'
       ELSE
          IF( vdiff_select( fieldlist_dry, 'q', k ) .NE. '' ) STOP 'call endrun( vdiff_select( fieldlist_dry, q, k ) )'
       ENDIF
20     CONTINUE
    END DO
    ! ------------------------ !
    ! Diagnostic output fields !
    ! ------------------------ !

    DO k = 1, pcnst
       vdiffnam(k) = 'VD'//cnst_name(k)
       IF( k == 1 ) vdiffnam(k) = 'VD01'    !**** compatibility with old code ****
       !call addfld( vdiffnam(k), 'kg/kg/s ', pver, 'A', 'Vertical diffusion of '//cnst_name(k), phys_decomp )
    END DO

    !  call phys_getopts( history_budget_out = history_budget )
    !  if( history_budget ) then
    !call add_default( vdiffnam(ixcldliq), 1, ' ' )
    !call add_default( vdiffnam(ixcldice), 1, ' ' )
    !  end if


  END SUBROUTINE vertical_diffusion_init
  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  !===============================================================================

  SUBROUTINE ubc_init()
    !-----------------------------------------------------------------------
    ! Initialization of time independent fields for the upper boundary condition
    ! Calls initialization routine for MSIS, TGCM and SNOE
    !-----------------------------------------------------------------------

  END SUBROUTINE ubc_init



  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !
  SUBROUTINE vertical_diffusion_tend(&
       pcols      , &!INTEGER , INTENT(IN   ) :: pcols                    ! Number of columns dimensioned
       ncol       , &!INTEGER , INTENT(IN   ) :: ncol                     !integer,  intent(in)  :: ncol! Number of columns actually used
       ncnst      , &!INTEGER , INTENT(IN   ) :: ncnst                    ! Number of constituents
       pver       , &!INTEGER , INTENT(IN   ) :: pver                     !integer,  intent(in)  :: pver ! Number of model layers
       ztodt      , &!REAL(r8), INTENT(in   ) :: ztodt                    ! 2 delta-t [ s ]
       colrad     , &! INTENT(IN   ) Cosino the colatitude [radian]
       mlsi       , &! I
       vcover     , &! I
       z0         , &
       TSK        , &! INTENT(IN   ) surface temperature
       QSFC       , &! INTENT(IN   ) surface specific temperature
       bps        , &!Exner function at layer interface   
       state_u    , &! REAL(r8), INTENT(in   ) :: state_u    (pcols,pver)  !real(r8), intent(in)  :: u(pcols,pver) ! Layer mid-point zonal wind [ m/s ]
       state_v    , &!REAL(r8), INTENT(in   ) :: state_v    (pcols,pver)  !real(r8), intent(in)  :: v(pcols,pver)  ! Layer mid-point meridional wind [ m/s ]
       state_t    , &!REAL(r8), INTENT(in   ) :: state_t    (pcols,pver)  !real(r8), intent(in)  :: t(pcols,pver)  ! Layer mid-point temperature [ K ]
       qm1        , &! initial/final constituent field
       state_qv   , &!REAL(r8), INTENT(in   ) :: state_qv   (pcols,pver) 
       state_ql   , &!REAL(r8), INTENT(in   ) :: state_ql   (pcols,pver) 
       state_qi   , &!REAL(r8), INTENT(in   ) :: state_qi   (pcols,pver) 
       state_pmid , &!REAL(r8), INTENT(in   ) :: state_pmid (pcols,pver)  !real(r8), intent(in)  :: pmid(pcols,pver)            ! Layer mid-point pressure [ Pa ]
       state_pint , &!REAL(r8), INTENT(in   ) :: state_pint (pcols,pver+1)!real(r8), intent(in)  :: pi(pcols,pver+1)      ! Interface pressure [ Pa ]
       state_exner, &!REAL(r8), INTENT(in   ) :: state_exner(pcols,pver)  !real(r8), intent(in)  :: exner(pcols,pver)            ! Layer mid-point exner function [ no unit ]
       state_zm   , &!REAL(r8), INTENT(in   ) :: state_zm   (pcols,pver)  !real(r8), intent(in)  :: zm(pcols,pver)        ! Layer mid-point height [ m ]
       state_zi   , &!REAL(r8), INTENT(in   ) :: state_zi   (pcols,pver+1)!real(r8), intent(in)  :: zi(pcols,pver+1)      ! Interface height above surface [ m ]
       state_rpdel, &!REAL(r8), INTENT(in   ) :: state_rpdel(pcols,pver)  ! 1./pdel where 'pdel' is thickness of the layer [ Pa ]
       state_pdel , &!REAL(r8), INTENT(in   ) :: state_rpdel(pcols,pver)  ! 1./pdel where 'pdel' is thickness of the layer [ Pa ]
       sgh        , &!REAL(r8), INTENT(in   ) :: sgh        (pcols)       !real(r8), intent(in)  :: sgh(pcols)                  ! Standard deviation of orography [ m ]
       landfrac   , &!REAL(r8), INTENT(in   ) :: landfrac   (pcols)       !real(r8), intent(in)  :: landfrac(pcols)       ! Land fraction [ fraction ]
       taux       , &!REAL(r8), INTENT(in   ) :: taux       (pcols)       ! x surface stress  [ N/m2 ]
       tauy       , &!REAL(r8), INTENT(in   ) :: tauy       (pcols)       ! y surface stress  [ N/m2 ]
       qrl        , &!REAL(r8), INTENT(in   ) :: qrl        (pcols,pver)  !qrl','g*W/m2',  pver,   'A',  'LW cooling rate, L', phys_decomp )
       wsedl      , &!REAL(r8), INTENT(in   ) :: wsedl      (pcols,pver)  !not used  ! Sedimentation velocity of liquid stratus cloud droplet [ m/s ]
       cldn       , &!REAL(r8), INTENT(in)    :: cldn       (pcols,pver)  !real(r8), intent(in)    :: cldn(pcols,pver)           ! Stratiform cloud fraction [ fraction ]
       shflx      , &!REAL(r8), INTENT(in)    :: shflx      (pcols)       !real(r8), intent(in)    :: shflx(pcols)           ! Sensible heat flux at surface [ unit ? ]
       cflx       , &!REAL(r8), INTENT(in)    :: cflx       (pcols,ncnst) !real(r8), intent(in)    :: qflx(pcols)           ! Water vapor flux at surface [ unit ? ]
       tauresx    , &!REAL(r8), INTENT(inout) :: tauresx    (pcols)       ! Residual stress to be added in vdiff to correct
       tauresy    , &!REAL(r8), INTENT(inout) :: tauresy    (pcols)       ! for turb stress mismatch between sfc and atm accumulated.
       kvm_in     , &!REAL(r8), INTENT(inout) :: kvm_in     (pcols,pver)  ! kvm saved from last timestep [ m2/s ]
       kvh_in     , &!REAL(r8), INTENT(inout) :: kvh_in     (pcols,pver)  ! kvh saved from last timestep [ m2/s ]

       dtv        , &!REAL(r8), INTENT(OUT  ) :: dtv   (pcols,plev)  ! temperature tendency (heating)
       dqv        , &!REAL(r8), INTENT(OUT  ) :: dqv   (pcols,plev,pcnst)  ! constituent diffusion tendency
       duv        , &!REAL(r8), INTENT(OUT  ) :: duv   (pcols,plev)   ! u-wind tendency
       dvv        , &!REAL(r8), INTENT(OUT  ) :: dvv   (pcols,plev)   ! v-wind tendency

       up1        , &!REAL(r8), INTENT(OUT  ) :: up1   (pcols,plev)! u-wind after vertical diffusion
       vp1        , &!REAL(r8), INTENT(OUT  ) :: vp1   (pcols,plev)! v-wind after vertical diffusion
       pblh       , &!REAL(r8), INTENT(OUT  ) :: pblh  (pcols)! planetary boundary layer height
       tpert      , &!REAL(r8), INTENT(OUT  ) :: tpert (pcols)! convective temperature excess
       qpert      , &!REAL(r8), INTENT(OUT  ) :: qpert (pcols,pcnst)! convective humidity and constituent excess
       tke        , &!real(r8), intent(inout) :: tke   (pcols,plev+1)
       rino       , &!REAL(r8), INTENT(INOUT) :: rino  (pcols,plev)! bulk Richardson no. from level to ref lev
       obklen     , &!REAL(r8), INTENT(OUT  ) :: obklen(pcols)! Obukhov length
       tstar      , &!REAL(r8), INTENT(OUT  ) :: tstar (pcols)
       wstar      , &!REAL(r8), INTENT(OUT  ) :: wstar (pcols)
       ustar      )
    IMPLICIT NONE
    INTEGER , INTENT(IN   ) :: pcols                    ! Number of columns dimensioned 
    INTEGER , INTENT(IN   ) :: ncol                     !integer,  intent(in)  :: ncol! Number of columns actually used
    INTEGER , INTENT(IN   ) :: pver                     !integer,  intent(in)  :: pver ! Number of model layers
    INTEGER , INTENT(IN   ) :: ncnst                    ! Number of constituents
    REAL(r8), INTENT(in   ) :: ztodt                    ! 2 delta-t [ s ]
    REAL(r8), INTENT(IN   ) :: colrad     (pcols)    
    INTEGER(KIND=i8), INTENT(IN   ) :: mlsi       (pCols)     
    REAL(r8), INTENT(IN   ) ::     vcover     (pCols)     
    REAL(r8), INTENT(IN   ) :: z0         (pCols)
    REAL(r8), INTENT(IN   ) :: TSK        (pcols)        ! INTENT(IN   ) surface temperature
    REAL(r8), INTENT(IN   ) :: QSFC       (pcols)       ! INTENT(IN   ) surface specific temperature
    REAL(r8), INTENT(in   ) :: bps        (pcols,pver)  !
    REAL(r8), INTENT(in   ) :: state_u    (pcols,pver)  !real(r8), intent(in)  :: u(pcols,pver) ! Layer mid-point zonal wind [ m/s ]
    REAL(r8), INTENT(in   ) :: state_v    (pcols,pver)  !real(r8), intent(in)  :: v(pcols,pver)  ! Layer mid-point meridional wind [ m/s ]
    REAL(r8), INTENT(in   ) :: state_t    (pcols,pver)  !real(r8), intent(in)  :: t(pcols,pver)  ! Layer mid-point temperature [ K ]
    REAL(r8), INTENT(in   ) :: state_qv   (pcols,pver) 
    REAL(r8), INTENT(in   ) :: state_ql   (pcols,pver) 
    REAL(r8), INTENT(in   ) :: state_qi   (pcols,pver) 
    REAL(r8), INTENT(in   ) :: state_pmid (pcols,pver)  !real(r8), intent(in)  :: pmid(pcols,pver)      ! Layer mid-point pressure [ Pa ]
    REAL(r8), INTENT(in   ) :: state_pint (pcols,pver+1)!real(r8), intent(in)  :: pi(pcols,pver+1)      ! Interface pressure [ Pa ]
    REAL(r8), INTENT(in   ) :: state_exner(pcols,pver)  !real(r8), intent(in)  :: exner(pcols,pver)     ! Layer mid-point exner function [ no unit ]
    REAL(r8), INTENT(in   ) :: state_zm   (pcols,pver)  !real(r8), intent(in)  :: zm(pcols,pver)        ! Layer mid-point height [ m ]
    REAL(r8), INTENT(in   ) :: state_zi   (pcols,pver+1)!real(r8), intent(in)  :: zi(pcols,pver+1)      ! Interface height above surface [ m ]
    REAL(r8), INTENT(in   ) :: state_pdel (pcols,pver)    ! layer thickness (Pa)

    REAL(r8), INTENT(in   ) :: state_rpdel(pcols,pver)  ! 1./pdel where 'pdel' is thickness of the layer [ Pa ]

    REAL(r8), INTENT(in   ) :: sgh        (pcols)       !real(r8), intent(in)  :: sgh(pcols)                  ! Standard deviation of orography [ m ]
    REAL(r8), INTENT(in   ) :: landfrac   (pcols)       !real(r8), intent(in)  :: landfrac(pcols)       ! Land fraction [ fraction ]
    REAL(r8), INTENT(in   ) :: taux       (pcols)       ! x surface stress  [ N/m2 ]
    REAL(r8), INTENT(in   ) :: tauy       (pcols)       ! y surface stress  [ N/m2 ]
    REAL(r8), INTENT(in   ) :: qrl        (pcols,pver)  !qrl','g*W/m2',  pver,   'A',  'LW cooling rate, L', phys_decomp )
    REAL(r8), INTENT(in   ) :: wsedl      (pcols,pver)  !not used  ! Sedimentation velocity of liquid stratus cloud droplet [ m/s ]
    REAL(r8), INTENT(in)    :: cldn       (pcols,pver)  !real(r8), intent(in)    :: cldn(pcols,pver)           ! Stratiform cloud fraction [ fraction ]
    REAL(r8), INTENT(in)    :: shflx      (pcols)       !real(r8), intent(in)    :: shflx(pcols)           ! Sensible heat flux at surface [ unit ? ]
    REAL(r8), INTENT(in)    :: cflx       (pcols,ncnst) !real(r8), intent(in)    :: qflx(pcols)           ! Water vapor flux at surface [ unit ? ]

    !
    ! Input/output arguments
    !
    REAL(kind=r8), INTENT(IN) :: qm1(pcols,pver,ncnst)  ! initial/final constituent field
    !
    ! Output arguments
    !


    REAL(r8), INTENT(inout) :: tauresx    (pcols)       ! Residual stress to be added in vdiff to correct
    REAL(r8), INTENT(inout) :: tauresy    (pcols)       ! for turb stress mismatch between sfc and atm accumulated.
    REAL(r8), INTENT(inout) :: kvm_in     (pcols,pver+1)  ! kvm saved from last timestep [ m2/s ]
    REAL(r8), INTENT(inout) :: kvh_in     (pcols,pver+1)  ! kvh saved from last timestep [ m2/s ]

    REAL(kind=r8), INTENT(OUT  ) :: dtv(pcols,pver)        ! temperature tendency (heating)
    REAL(kind=r8), INTENT(OUT  ) :: dqv(pcols,pver,pcnst)  ! constituent diffusion tendency
    REAL(kind=r8), INTENT(OUT  ) :: duv(pcols,pver)        ! u-wind tendency
    REAL(kind=r8), INTENT(OUT  ) :: dvv(pcols,pver)        ! v-wind tendency

    REAL(kind=r8), INTENT(OUT  ) :: up1(pcols,pver)        ! u-wind after vertical diffusion
    REAL(kind=r8), INTENT(OUT  ) :: vp1(pcols,pver)        ! v-wind after vertical diffusion
    
    REAL(kind=r8), INTENT(OUT  ) :: pblh(pcols)            ! planetary boundary layer height
    REAL(kind=r8), INTENT(OUT  ) :: tpert(pcols)           ! convective temperature excess
    REAL(kind=r8), INTENT(OUT  ) :: qpert(pcols,pcnst)     ! convective humidity and constituent excess
    REAL(kind=r8), INTENT(INOUT) :: TKE(pcols,pver+1)
    REAL(kind=r8), INTENT(INOUT) :: rino(pcols,pver)        ! bulk Richardson no. from level to ref lev
    REAL(kind=r8), INTENT(OUT  ) :: obklen(pcols)           ! Obukhov length
    REAL(kind=r8), INTENT(OUT  ) :: tstar(pcols)    
    REAL(kind=r8), INTENT(OUT  ) :: wstar(pcols)    
    REAL(r8)      , INTENT(OUT  ):: ustar(pcols)          ! Surface friction velocity [ m/s ]

    REAL(r8)   :: ptend_q(pcols,pver,ncnst) 
    REAL(r8)   :: ptend_s(pcols,pver)   
    REAL(r8)   :: ptend_u(pcols,pver)   
    REAL(r8)   :: ptend_v(pcols,pver)   

    REAL(r8)   :: state_s(pcols,pver)  !real(r8), intent(in)  :: t(pcols,pver)  ! Layer mid-point temperature [ K ]

    REAL(r8)   :: state_pmiddry(pcols,pver)  
    REAL(r8)   :: state_pintdry(pcols,pver+1)
    REAL(r8)   :: state_rpdeldry(pcols,pver)
    !REAL(r8)   :: ustar(pcols)          ! Surface friction velocity [ m/s ]
    !REAL(r8)   :: pblh(pcols)                ! PBL top height [ m ]
    !REAL(r8)   :: obklen(pcols)             ! Obukhov length [ m ]
    REAL(r8)   :: tpertPBL(pcols)
    REAL(r8)   :: qpertPBL(pcols)
    REAL(r8) :: slv_prePBL(pcols,pver)
    REAL(r8) :: slten(pcols,pver)
    REAL(r8) :: qtten(pcols,pver)
    REAL(r8) :: tem2(pcols,pver)                                    ! Saturation specific humidity and RH
    REAL(r8) :: ftem(pcols,pver)                                    ! Saturation vapor pressure before PBL
    REAL(r8) :: ftem_prePBL(pcols,pver)                             ! Saturation vapor pressure before PBL

    !real(r8)   :: kvm_out(pcols,pver+1) ! Eddy diffusivity for momentum [ m2/s ]
    !real(r8)   :: kvh_out(pcols,pver+1) ! Eddy diffusivity for heat [ m2/s ]
    !real(r8)   :: kvq(pcols,pver+1)        ! Eddy diffusivity for constituents, moisture and tracers [ m2/s ] (note not having '_out')
    REAL(r8)   :: smaw(pcols,pver+1)!real(r8), intent(out)   :: sm_aw(pcols,pver+1)   ! Normalized Galperin instability function for momentum [ no unit ]
    REAL(r8)   :: cgh(pcols,pver+1)        ! Counter-gradient term for heat [ J/kg/m ]
    REAL(r8)   :: cgs(pcols,pver+1)        ! Counter-gradient star [ cg/flux ]
    !REAL(r8)   :: tpert(pcols)          ! Convective temperature excess [ K ]
    REAL(r8)   :: qpert_loc(pcols)          ! Convective humidity excess [ kg/kg ]
    REAL(r8)   :: wpert(pcols)          ! Turbulent velocity excess [ m/s ]
    !REAL(r8)   :: tke(pcols,pver+1)        ! Turbulent kinetic energy [ m2/s2 ]
    REAL(r8)   :: bprod(pcols,pver+1)        ! Buoyancy production [ m2/s3 ] 
    REAL(r8)   :: sprod(pcols,pver+1)        ! Shear production [ m2/s3 ] 
    REAL(r8)   :: sfi(pcols,pver+1)        ! Interfacial layer saturation fraction [ fraction ]
    !integer ,external   :: fqsatd
    !integer,  external :: compute_molec_diff   ! Constituent-independent moleculuar diffusivity routine

    REAL(r8)   :: ipbl(pcols)                ! If 1, PBL is CL, while if 0, PBL is STL.
    REAL(r8)   :: kpblh(pcols)          ! Layer index containing PBL top within or at the base interface
    REAL(r8)   :: wstarPBL(pcols)        ! Convective velocity within PBL [ m/s ]
    REAL(r8)   :: turbtype(pcols,pver+1)! Turbulence type identifier at all interfaces [ no unit ]
    REAL(r8)   :: kvm (pcols,pver+1) ! Eddy diffusivity for momentum [ m2/s ]
    REAL(r8)   :: kvh (pcols,pver+1) ! Eddy diffusivity for heat [ m2/s ]
    REAL(r8)   :: kvq (pcols,pver+1)     ! Eddy diffusivity for constituents, moisture and tracers [ m2/s ] (note not having '_out')
    REAL(r8)   :: sl_prePBL(pcols,pver)
    REAL(r8)   :: qt_prePBL(pcols,pver)

    LOGICAL  :: kvinit                                              ! Tell compute_eddy_diff/ caleddy to initialize kvh, kvm (uses kvf)
    REAL(r8) :: rztodt                                              ! 1./ztodt [ 1/s ]
    REAL(r8) :: ksrftms(pcols)    !real(r8), intent(out) :: ksrf(pcols)        ! Surface drag coefficient [ kg/s/m2 ]
    REAL(r8) :: tautmsx(pcols)    !real(r8), intent(out) :: taux(pcols)        ! Surface zonal        wind stress [ N/m2 ]
    REAL(r8) :: tautmsy(pcols)    !real(r8), intent(out) :: tauy(pcols)        ! Surface meridional wind stress [ N/m2 ]
    REAL(r8) :: tautotx(pcols)    ! U component of total surface stress [ N/m2 ]
    REAL(r8) :: tautoty(pcols)    ! V component of total surface stress [ N/m2 ]
    REAL(r8) :: dtk(pcols,pver)                                     ! T tendency from KE dissipation
    REAL(r8) :: topflx(pcols)                                       ! Molecular heat flux at top interface
    REAL(r8) :: sl(pcols,pver)
    REAL(r8) :: qt(pcols,pver)
    REAL(r8) :: slv(pcols,pver)
    REAL(r8) :: slflx(pcols,pver+1)
    REAL(r8) :: qtflx(pcols,pver+1)
    REAL(r8) :: uflx(pcols,pver+1)
    REAL(r8) :: vflx(pcols,pver+1)
    REAL(r8) :: slflx_cg(pcols,pver+1)
    REAL(r8) :: qtflx_cg(pcols,pver+1)
    REAL(r8) :: uflx_cg(pcols,pver+1)
    REAL(r8) :: vflx_cg(pcols,pver+1)
    REAL(r8) :: qv_pro(pcols,pver) 
    REAL(r8) :: ql_pro(pcols,pver)
    REAL(r8) :: qi_pro(pcols,pver)
    REAL(r8) :: s_pro(pcols,pver)
    REAL(r8) :: t_pro(pcols,pver)
    REAL(r8) :: qv_aft_PBL(pcols,pver)                              ! qv after PBL diffusion
    REAL(r8) :: ql_aft_PBL(pcols,pver)                              ! ql after PBL diffusion
    REAL(r8) :: qi_aft_PBL(pcols,pver)                              ! qi after PBL diffusion
    REAL(r8) :: s_aft_PBL(pcols,pver)                               ! s after PBL diffusion
    REAL(r8) :: u_aft_PBL(pcols,pver)                               ! u after PBL diffusion
    REAL(r8) :: v_aft_PBL(pcols,pver)                               ! v after PBL diffusion
    REAL(r8) :: t_aftPBL(pcols,pver)                                ! Temperature after PBL diffusion
    REAL(r8) :: ftem_aftPBL(pcols,pver)                             ! Saturation vapor pressure after PBL
    REAL(r8) :: tten(pcols,pver)                                    ! Temperature tendency by PBL diffusion
    REAL(r8) :: rhten(pcols,pver)                                   ! RH tendency by PBL diffusion

    REAL(r8) :: rhoair


    CHARACTER(128) :: errstring   ! Error status for init_vdiff


    INTEGER  :: lchnk ,i,k  ,m                                            ! Chunk identifier
    INTEGER  :: time_index                                          ! Time level index for physics buffer access
    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !
    rhoair=0.0_r8
    DO m=1,pcnst
       DO k=1,pver
          DO i=1,ncol
             dqv(i,k,m) =0.0_r8 ! constituent diffusion tendency
             ptend_q(i,k,m)  =0.0_r8 
          END DO
       END DO
    END DO

    DO k=1,pver
       DO i=1,ncol
          dtv(i,k)=0.0_r8 ! temperature tendency (heating)
          duv(i,k)=0.0_r8 ! u-wind tendency
          dvv(i,k)=0.0_r8 ! v-wind tendency
          ptend_s(i,k)=0.0_r8 
          ptend_u(i,k)=0.0_r8 
          ptend_v(i,k)=0.0_r8 
          state_s(i,k)=0.0_r8 !real(r8), intent(in)  :: t(pcols,pver)  ! Layer mid-point temperature [ K ]
          state_pmiddry(i,k)=0.0_r8 
          up1(i,k)=0.0_r8 ! u-wind after vertical diffusion
          vp1(i,k)=0.0_r8 ! v-wind after vertical diffusion
          state_rpdeldry(i,k)=0.0_r8
          slv_prePBL(i,k) =0.0_r8
          slten(i,k) =0.0_r8
          qtten(i,k) =0.0_r8
          tem2(i,k) =0.0_r8                                    ! Saturation specific humidity and RH
          ftem(i,k) =0.0_r8                                    ! Saturation vapor pressure before PBL
          ftem_prePBL(i,k) =0.0_r8                             ! Saturation vapor pressure before PBL
          sl_prePBL(i,k) =0.0_r8
          qt_prePBL(i,k) =0.0_r8
          dtk(i,k) =0.0_r8                                     ! T tendency from KE dissipation
          sl(i,k) =0.0_r8
          qt(i,k) =0.0_r8
          slv(i,k) =0.0_r8
          qv_pro(i,k) =0.0_r8 
          ql_pro(i,k) =0.0_r8
          qi_pro(i,k) =0.0_r8
          s_pro(i,k) =0.0_r8
          t_pro(i,k) =0.0_r8
          qv_aft_PBL(i,k) =0.0_r8                              ! qv after PBL diffusion
          ql_aft_PBL(i,k) =0.0_r8                              ! ql after PBL diffusion
          qi_aft_PBL(i,k) =0.0_r8                              ! qi after PBL diffusion
          s_aft_PBL(i,k) =0.0_r8                               ! s after PBL diffusion
          u_aft_PBL(i,k) =0.0_r8                               ! u after PBL diffusion
          v_aft_PBL(i,k) =0.0_r8                               ! v after PBL diffusion
          t_aftPBL(i,k) =0.0_r8                                ! Temperature after PBL diffusion
          ftem_aftPBL(i,k) =0.0_r8                             ! Saturation vapor pressure after PBL
          tten(i,k) =0.0_r8                                    ! Temperature tendency by PBL diffusion
          rhten(i,k) =0.0_r8                                   ! RH tendency by PBL diffusion
       END DO
    END DO
    DO k=1,pver+1
       DO i=1,ncol
          state_pintdry(i,k)=0.0_r8
          smaw (i,k)=0.0_r8!real(r8), intent(out)   :: sm_aw(i,k)   ! Normalized Galperin instability function for momentum [ no unit ]
          cgh  (i,k)=0.0_r8! Counter-gradient term for heat [ J/kg/m ]
          cgs  (i,k)=0.0_r8! Counter-gradient star [ cg/flux ]
          bprod(i,k)=0.0_r8! Buoyancy production [ m2/s3 ] 
          sprod(i,k)=0.0_r8! Shear production [ m2/s3 ] 
          sfi  (i,k)=0.0_r8! Interfacial layer saturation fraction [ fraction ]
          turbtype(i,k)=0.0_r8! Turbulence type identifier at all interfaces [ no unit ]
          kvm (i,k)=0.0_r8! Eddy diffusivity for momentum [ m2/s ]
          kvh (i,k)=0.0_r8! Eddy diffusivity for heat [ m2/s ]
          kvq (i,k)=0.0_r8! Eddy diffusivity for constituents, moisture and tracers [ m2/s ] (note not having '_out')
          slflx(i,k)=0.0_r8
          qtflx(i,k)=0.0_r8
          uflx(i,k)=0.0_r8
          vflx(i,k)=0.0_r8
          slflx_cg(i,k)=0.0_r8
          qtflx_cg(i,k)=0.0_r8
          uflx_cg(i,k)=0.0_r8
          vflx_cg(i,k)=0.0_r8
       END DO
    END DO

    DO m=1,pcnst
       DO i=1,ncol
          qpert(i,m)=0.0_r8 ! convective humidity and constituent excess
       END DO
    END DO
    DO i=1,ncol
       tpertPBL(i)=0.0_r8
       qpertPBL(i)=0.0_r8
       wpert(i)=0.0_r8! Turbulent velocity excess [ m/s ]
       qpert_loc(i)=0.0_r8! Convective humidity excess [ kg/kg ]
       ipbl(i)=0.0_r8! If 1, PBL is CL, while if 0, PBL is STL.
       kpblh(i)=0.0_r8! Layer index containing PBL top within or at the base interface
       wstarPBL(i)=0.0_r8! Convective velocity within PBL [ m/s ]
       ksrftms(i)=0.0_r8!real(r8), intent(out) :: ksrf(pcols)! Surface drag coefficient [ kg/s/m2 ]
       tautmsx(i)=0.0_r8!real(r8), intent(out) :: taux(pcols)! Surface zonal wind stress [ N/m2 ]
       tautmsy(i)=0.0_r8!real(r8), intent(out) :: tauy(pcols)! Surface meridional wind stress [ N/m2 ]
       tautotx(i)=0.0_r8! U component of total surface stress [ N/m2 ]
       tautoty(i)=0.0_r8! V component of total surface stress [ N/m2 ]
       topflx(i)=0.0_r8! Molecular heat flux at top interface
       pblh(i)=0.0_r8 ! planetary boundary layer height
       tpert(i)=0.0_r8 ! convective temperature excess
       obklen(i)=0.0_r8 ! Obukhov length
       tstar(i)=0.0_r8 
       wstar(i)=0.0_r8 
    END DO

    !------
    rztodt = 1._r8 / ztodt

    
    !state_s(pcols,pver)  !real(r8), intent(in)  :: t(pcols,pver)  ! Layer mid-point temperature [ K ]
    !
    !t_aftPBL(:ncol,:pver)    = ( s_aft_PBL(:ncol,:pver) - gravit*state_zm(:ncol,:pver) ) / cpair 
    !
    !cpair * t_aftPBL(:ncol,:pver)    =  s_aft_PBL(:ncol,:pver) - gravit*state_zm(:ncol,:pver) 
    !
    !cpair * t_aftPBL(:ncol,:pver)+ gravit*state_zm(:ncol,:pver)     =  s_aft_PBL(:ncol,:pver) 
    !
    DO k=1,pver
       DO i=1,ncol
          state_s(i,k) =cpair * state_t(i,k)+ gravit*state_zm(i,k)
       END DO
    END DO

    !lchnk  = state%lchnk
    !ncol   = state%ncol

    IF( is_first_step) THEN
       ! tauresx(:ncol) = 0._r8
       ! tauresy(:ncol) = 0._r8
    ELSE
       ! tauresx(:ncol) = pbuf(tauresx_idx)%fld_ptr(1,1:ncol,1,lchnk,time_index)
       ! tauresy(:ncol) = pbuf(tauresy_idx)%fld_ptr(1,1:ncol,1,lchnk,time_index)
    ENDIF

    ! All variables are modified by vertical diffusion

    !ptend%name  = "vertical diffusion"
    !ptend%lq(:) = .TRUE.
    !ptend%ls    = .TRUE.
    !ptend%lu    = .TRUE.
    !ptend%lv    = .TRUE.

    ! ---------------------------------------- !
    ! Computation of turbulent mountain stress !
    ! ---------------------------------------- !

    ! Consistent with the computation of 'normal' drag coefficient, we are using 
    ! the raw input (u,v) to compute 'ksrftms', not the provisionally-marched 'u,v' 
    ! within the iteration loop of the PBL scheme. 

    IF( do_tms ) THEN
       CALL compute_tms( pcols      , &!integer,  intent(in)  :: pcols                 ! Number of columns dimensioned
            pver       , &!integer,  intent(in)  :: pver                  ! Number of model layers
            ncol       , &!integer,  intent(in)  :: ncol                  ! Number of columns actually used
            state_u    , &!real(r8), intent(in)  :: u(pcols,pver)         ! Layer mid-point zonal wind [ m/s ]
            state_v    , &!real(r8), intent(in)  :: v(pcols,pver)         ! Layer mid-point meridional wind [ m/s ]
            state_t    , &!real(r8), intent(in)  :: t(pcols,pver)         ! Layer mid-point temperature [ K ]
            state_pmid , &!real(r8), intent(in)  :: pmid(pcols,pver)      ! Layer mid-point pressure [ Pa ]
            state_exner, &!real(r8), intent(in)  :: exner(pcols,pver)     ! Layer mid-point exner function [ no unit ]
            state_zm   , &!real(r8), intent(in)  :: zm(pcols,pver)        ! Layer mid-point height [ m ]
            sgh        , &!real(r8), intent(in)  :: sgh(pcols)            ! Standard deviation of orography [ m ]
            ksrftms    , &!real(r8), intent(out) :: ksrf(pcols)           ! Surface drag coefficient [ kg/s/m2 ]
            tautmsx    , &!real(r8), intent(out) :: taux(pcols)           ! Surface zonal      wind stress [ N/m2 ]
            tautmsy    , &!real(r8), intent(out) :: tauy(pcols)           ! Surface meridional wind stress [ N/m2 ]
            landfrac     )!real(r8), intent(in)  :: landfrac(pcols)       ! Land fraction [ fraction ]
       ! Here, both 'taux, tautmsx' are explicit surface stresses.        
       ! Note that this 'tautotx, tautoty' are different from the total stress
       ! that has been actually added into the atmosphere. This is because both
       ! taux and tautmsx are fully implicitly treated within compute_vdiff.
       ! However, 'tautotx, tautoty' are not used in the actual numerical
       ! computation in this module.   
       DO i=1,ncol
          tautotx(i) = taux(i) + tautmsx(i)
          tautoty(i) = tauy(i) + tautmsy(i)
       END DO
    ELSE
       DO i=1,ncol
          ksrftms(i) = 0._r8
          tautotx(i) = taux(i)
          tautoty(i) = tauy(i)
       END DO
    ENDIF

    !----------------------------------------------------------------------- !
    !   Computation of eddy diffusivities - Select appropriate PBL scheme    !
    !----------------------------------------------------------------------- !

    SELECT CASE (eddy_scheme)
    CASE ( 'diag_TKE' ) 

       ! ---------------------------------------------------------------- !
       ! At first time step, have eddy_diff.F90:caleddy() use kvh=kvm=kvf !
       ! This has to be done in compute_eddy_diff after kvf is calculated !
       ! ---------------------------------------------------------------- !

       IF( is_first_step) THEN
          kvinit = .TRUE.
       ELSE
          kvinit = .FALSE.
       ENDIF

       ! ---------------------------------------------- !
       ! Get LW radiative heating out of physics buffer !
       ! ---------------------------------------------- !

       ! qrl  (pcols,pver)     !  => pbuf(pbuf_get_fld_idx('QRL'  ))%fld_ptr(1,1:pcols,1:pver,lchnk,1)
       ! wsedl(pcols,pver)     !  => pbuf(pbuf_get_fld_idx('WSEDL'))%fld_ptr(1,1:pcols,1:pver,lchnk,1)

       ! Retrieve eddy diffusivities for heat and momentum from physics buffer
       ! from last timestep ( if first timestep, has been initialized by inidat.F90 )

       !time_index      = pbuf_old_tim_idx()
       !kvm_in(:ncol,:) = pbuf(kvm_idx)%fld_ptr(1,1:ncol,1:pverp,lchnk,time_index)
       !kvh_in(:ncol,:) = pbuf(kvh_idx)%fld_ptr(1,1:ncol,1:pverp,lchnk,time_index)

       CALL compute_eddy_diff( &
            lchnk                  , &!integer,  intent(in)    :: lchnk   
            pcols                  , &!integer,  intent(in)    :: pcols                 ! Number of atmospheric columns [ # ]
            pver                   , &!integer,  intent(in)    :: pver                  ! Number of atmospheric layers  [ # ]
            ncol                   , &!integer,  intent(in)    :: ncol                  ! Number of atmospheric columns [ # ]
            mlsi                   , &
            vcover                 , &
            z0                     , & 
            TSK                    , & 
            QSFC                   , & 
            bps                    , & 
            state_t                , &!real(r8), intent(in)    :: t(pcols,pver)         ! Temperature [K]
            state_qv               , &!real(r8), intent(in)    :: qv(pcols,pver)        ! Water vapor  specific humidity [ kg/kg ]
            ztodt                  , &!real(r8), intent(in)    :: ztodt                 ! Physics integration time step 2 delta-t [ s ]
            state_ql               , &!real(r8), intent(in)    :: ql(pcols,pver)        ! Liquid water specific humidity [ kg/kg ]
            state_qi               , &!real(r8), intent(in)    :: qi(pcols,pver)        ! Ice specific humidity [ kg/kg ]
            state_rpdel            , &!real(r8), intent(in)    :: rpdel(pcols,pver)     ! 1./pdel where 'pdel' is thickness of the layer [ Pa ]
            cldn                   , &!real(r8), intent(in)    :: cldn(pcols,pver)      ! Stratiform cloud fraction [ fraction ]
            qrl                    , &!real(r8), intent(in)    :: qrl(pcols,pver)       ! LW cooling rate
            wsedl                  , &!real(r8), intent(in)    :: wsedl(pcols,pver)     ! Sedimentation velocity of liquid stratus cloud droplet [ m/s ]
            state_zm               , &!real(r8), intent(in)    :: z(pcols,pver)         ! Layer mid-point height above surface [ m ]
            state_zi               , &!real(r8), intent(in)    :: zi(pcols,pver+1)      ! Interface height above surface [ m ]
            state_pmid             , &!real(r8), intent(in)    :: pmid(pcols,pver)      ! Layer mid-point pressure [ Pa ]
            state_pint             , &!real(r8), intent(in)    :: pi(pcols,pver+1)      ! Interface pressure [ Pa ]
            state_u                , &!real(r8), intent(in)    :: u(pcols,pver)         ! Zonal velocity [ m/s ]
            state_v                , &!real(r8), intent(in)    :: v(pcols,pver)         ! Meridional velocity [ m/s ]
            taux                   , &!real(r8), intent(in)    :: taux(pcols)           ! Zonal wind stress at surface [ N/m2 ]
            tauy                   , &!real(r8), intent(in)    :: tauy(pcols)           ! Meridional wind stress at surface [ N/m2 ]
            shflx                  , &!real(r8), intent(in)    :: shflx(pcols)          ! Sensible heat flux at surface [ unit ? ]
            cflx(:,1)              , &!real(r8), intent(in)    :: qflx(pcols)           ! Water vapor flux at surface [ unit ? ]
            wstarent               , &!logical,  intent(in)    :: wstarent              ! .true. means use the 'wstar' entrainment closure. 
            nturb                  , &!integer,  intent(in)    :: nturb                 ! Number of iteration steps for calculating eddy diffusivity [ # ]
            ustar                  , &!real(r8), intent(out)   :: ustar(pcols)          ! Surface friction velocity [ m/s ]
            pblh                   , &!real(r8), intent(out)   :: pblh(pcols)           ! PBL top height [ m ]
            landfrac               , &
            kvm_in                 , &!real(r8), intent(in)    :: kvm_in(pcols,pver+1)  ! kvm saved from last timestep [ m2/s ]
            kvh_in                 , &!real(r8), intent(in)    :: kvh_in(pcols,pver+1)  ! kvh saved from last timestep [ m2/s ]
            kvm                    , &!real(r8), intent(out)   :: kvm_out(pcols,pver+1) ! Eddy diffusivity for momentum [ m2/s ]
            kvh                    , &!real(r8), intent(out)   :: kvh_out(pcols,pver+1) ! Eddy diffusivity for heat [ m2/s ]
            kvq                    , &!real(r8), intent(out)   :: kvq(pcols,pver+1)     ! Eddy diffusivity for constituents, moisture and tracers [ m2/s ] (note not having '_out')
            cgh                    , &!real(r8), intent(out)   :: cgh(pcols,pver+1)     ! Counter-gradient term for heat [ J/kg/m ]
            cgs                    , &!real(r8), intent(out)   :: cgs(pcols,pver+1)     ! Counter-gradient star [ cg/flux ]
            tpert                  , &!real(r8), intent(out)   :: tpert(pcols)          ! Convective temperature excess [ K ]
            qpert_loc              , &!real(r8), intent(out)   :: qpert(pcols)          ! Convective humidity excess [ kg/kg ]
            wpert                  , &!real(r8), intent(out)   :: wpert(pcols)          ! Turbulent velocity excess [ m/s ]
            tke                    , &!real(r8), intent(out)   :: tke(pcols,pver+1)     ! Turbulent kinetic energy [ m2/s2 ]
            bprod                  , &!real(r8), intent(out)   :: bprod(pcols,pver+1)   ! Buoyancy production [ m2/s3 ] 
            sprod                  , &!real(r8), intent(out)   :: sprod(pcols,pver+1)   ! Shear production [ m2/s3 ] 
            sfi                    , &!real(r8), intent(out)   :: sfi(pcols,pver+1)     ! Interfacial layer saturation fraction [ fraction ]
            kvinit                 , &!logical, intent(in  )   :: kvinit                ! Tell compute_eddy_diff/ caleddy to initialize kvh, kvm (uses kvf)
            tauresx                , &!real(r8), intent(inout) :: tauresx(pcols)        ! Residual stress to be added in vdiff to correct for turb
            tauresy                , &!real(r8), intent(inout) :: tauresy(pcols)        ! Stress mismatch between sfc and atm accumulated in prior timesteps
            ksrftms                , &!real(r8), intent(in)    :: ksrftms(pcols)        ! Surface drag coefficient of turbulent mountain stress [ unit ? ]
            ipbl(:)                , &!real(r8), intent(out)   :: ipbl(pcols)           ! If 1, PBL is CL, while if 0, PBL is STL.
            kpblh(:)               , &!real(r8), intent(out)   :: kpblh(pcols)          ! Layer index containing PBL top within or at the base interface
            wstarPBL(:)            , &!real(r8), intent(out)   :: wstarPBL(pcols)       ! Convective velocity within PBL [ m/s ]
            turbtype               , &!real(r8), intent(out)   :: turbtype(pcols,pver+1)! Turbulence type identifier at all interfaces [ no unit ]
            smaw                   , &!real(r8), intent(out)   :: sm_aw(pcols,pver+1)   ! Normalized Galperin instability function for momentum [ no unit ]
            rino                     )!REAL(r8), INTENT(out)   :: ri(pcols,pver)            ! Richardson number, 'n2/s2', defined at interfaces except surface [ s-2 ]

            obklen(1:ncol) = 0._r8 
            kvm_in=kvm
            kvh_in=kvh

       ! ----------------------------------------------- !       
       ! Store TKE in pbuf for use by shallow convection !
       ! ----------------------------------------------- !   
       DO i=1,ncol
          wstar(i)=wstarPBL(i)
          tpertPBL(i) = tpert(i)
          qpertPBL(i) = qpert_loc(i)
          qpert   (i,1) = qpert_loc(i)
       END DO
       !pbuf(tke_idx)%fld_ptr(1,1:ncol,1:pverp,lchnk,time_index)      = tke(:ncol,:)
       !pbuf(turbtype_idx)%fld_ptr(1,1:ncol,1:pverp,lchnk,time_index) = turbtype(:ncol,:)
       !pbuf(smaw_idx)%fld_ptr(1,1:ncol,1:pverp,lchnk,time_index)     = smaw(:ncol,:)

       ! Store updated kvh, kvm in pbuf to use here on the next timestep 

       !pbuf(kvh_idx)%fld_ptr(1,1:ncol,1:pverp,lchnk,time_index) = kvh(:ncol,:)
       !pbuf(kvm_idx)%fld_ptr(1,1:ncol,1:pverp,lchnk,time_index) = kvm(:ncol,:)
       !if( is_first_step() ) then
       !   do i = 1, pbuf_times
       !      pbuf(kvh_idx)%fld_ptr(1,1:ncol,1:pverp,lchnk,i) = kvh(:ncol,:)
       !      pbuf(kvm_idx)%fld_ptr(1,1:ncol,1:pverp,lchnk,i) = kvm(:ncol,:)
       !   enddo
       !endif
       ! Write out fields that are only used by this scheme

       !call outfld( 'BPROD   ', bprod(1,1), pcols, lchnk )
       !call outfld( 'SPROD   ', sprod(1,1), pcols, lchnk )
       !call outfld( 'SFI     ', sfi,        pcols, lchnk )



    CASE ( 'HB', 'HBR' )
       ! Modification : We may need to use 'taux' instead of 'tautotx' here, for
       !                consistency with the previous HB scheme.


    END SELECT

    !pbuf(wgustd_index)%fld_ptr(1,1:ncol,1,lchnk,1) = wpert(:ncol)
    !call outfld( 'WGUSTD' , wpert, pcols, lchnk )

    !------------------------------------ ! 
    !    Application of diffusivities     !
    !------------------------------------ !
    !qm1
    DO m=1,ncnst
       DO k=1,pver
          DO i=1,ncol
             !ptend_q(i,k,1) = state_qv(i,k)
             !ptend_q(i,k,2) = state_qi(i,k)
             !ptend_q(i,k,3) = state_ql(i,k)
             ptend_q(i,k,m) = qm1(i,k,m)
          END DO
       END DO
    END DO
 
    DO k=1,pver
       DO i=1,ncol
          ptend_s(i,k)   = state_s(i,k)
          ptend_u(i,k)   = state_u(i,k)
          ptend_v(i,k)   = state_v(i,k)
       END DO
    END DO   
    !------------------------------------------------------ !
    ! Write profile output before applying diffusion scheme !
    !------------------------------------------------------ !
    DO k=1,pver
       DO i=1,ncol
    
          sl_prePBL(i,k)  = ptend_s(i,k) -   latvap           * ptend_q(i,k,ixcldliq) &
                                         - ( latvap + latice) * ptend_q(i,k,ixcldice)
          qt_prePBL(i,k)  = ptend_q(i,k,1) + ptend_q(i,k,ixcldliq) &
                                           + ptend_q(i,k,ixcldice)
          slv_prePBL(i,k) = sl_prePBL(i,k) * ( 1._r8 + zvir*qt_prePBL(i,k) ) 

       END DO
    END DO   

    CALL aqsat( state_t, state_pmid, tem2, ftem, pcols, ncol, pver, 1, pver )
   
    ! ftem is Saturation vapor pressure before PBL
    
    DO k=1,pver
       DO i=1,pcols
          ftem_prePBL(i,k) = state_qv(i,k)/ftem(i,k)*100._r8 ! relative humidity
       END DO
    END DO

    !call outfld( 'qt_pre_PBL   ', qt_prePBL,                 pcols, lchnk )
    !call outfld( 'sl_pre_PBL   ', sl_prePBL,                 pcols, lchnk )
    !call outfld( 'slv_pre_PBL  ', slv_prePBL,                pcols, lchnk )
    !call outfld( 'u_pre_PBL    ', state%u,                   pcols, lchnk )
    !call outfld( 'v_pre_PBL    ', state%v,                   pcols, lchnk )
    !call outfld( 'qv_pre_PBL   ', state%q(:ncol,:,1),        pcols, lchnk )
    !call outfld( 'ql_pre_PBL   ', state%q(:ncol,:,ixcldliq), pcols, lchnk )
    !call outfld( 'qi_pre_PBL   ', state%q(:ncol,:,ixcldice), pcols, lchnk )
    !call outfld( 't_pre_PBL    ', state%t,                   pcols, lchnk )
    !call outfld( 'rh_pre_PBL   ', ftem_prePBL,               pcols, lchnk )

    ! --------------------------------------------------------------------------------- !
    ! Call the diffusivity solver and solve diffusion equation                          !
    ! The final two arguments are optional function references to                       !
    ! constituent-independent and constituent-dependent moleculuar diffusivity routines !
    ! --------------------------------------------------------------------------------- !

    ! Modification : We may need to output 'tautotx_im,tautoty_im' from below 'compute_vdiff' and
    !                separately print out as diagnostic output, because these are different from
    !                the explicit 'tautotx, tautoty' computed above. 
    ! Note that the output 'tauresx,tauresy' from below subroutines are fully implicit ones.

    IF( ANY(fieldlist_wet) ) THEN

       CALL compute_vdiff( &
            lchnk              , & !integer,  intent(in)    :: lchnk
            pcols              , & !integer,  intent(in)    :: pcols
            pver               , & !integer,  intent(in)    :: pver
            pcnst              , & !integer,  intent(in)    :: ncnst
            ncol               , & !integer,  intent(in)    :: ncol                      ! Number of atmospheric columns
            TSK                , & 
            QSFC               , & 
            state_pmid         , & !real(r8), intent(in)    :: pmid(pcols,pver)          ! Mid-point pressures [ Pa ]
            state_pint         , & !real(r8), intent(in)    :: pint(pcols,pver+1)        ! Interface pressures [ Pa ]
            state_rpdel        , & !real(r8), intent(in)    :: rpdel(pcols,pver)         ! 1./pdel
            state_t            , & !real(r8), intent(in)    :: t(pcols,pver)             ! Temperature [ K ]
            ztodt              , & !real(r8), intent(in)    :: ztodt                     ! 2 delta-t [ s ]
            taux               , & !real(r8), intent(in)    :: taux(pcols)               ! Surface zonal      stress. Input u-momentum per unit time per unit area into the atmosphere [ N/m2 ]
            tauy               , & !real(r8), intent(in)    :: tauy(pcols)               ! Surface meridional stress. Input v-momentum per unit time per unit area into the atmosphere [ N/m2 ]
            shflx              , & !real(r8), intent(in)    :: shflx(pcols)              ! Surface sensible heat flux [ W/m2 ]
            cflx               , & !real(r8), intent(in)    :: cflx(pcols,ncnst)         ! Surface constituent flux [ kg/m2/s ]
            ntop               , & !integer,  intent(in)    :: ntop                      ! Top    interface level to which vertical diffusion is applied ( = 1 ).
            nbot               , & !integer,  intent(in)    :: nbot                      ! Bottom interface level to which vertical diffusion is applied ( = pver ).
            kvh                , & !real(r8), intent(inout) :: kvh(pcols,pver+1)         ! Eddy diffusivity for heat [ m2/s ]
            kvm                , & !real(r8), intent(inout) :: kvm(pcols,pver+1)         ! Eddy viscosity ( Eddy diffusivity for momentum ) [ m2/s ]
            kvq                , & !real(r8), intent(inout) :: kvq(pcols,pver+1)         ! Eddy diffusivity for constituents
            cgs                , & !real(r8), intent(inout) :: cgs(pcols,pver+1)         ! Counter-gradient star [ cg/flux ]
            cgh                , & !real(r8), intent(inout) :: cgh(pcols,pver+1)         ! Counter-gradient term for heat
            state_zi           , & !real(r8), intent(in)    :: zi(pcols,pver+1)          ! Interface heights [ m ]
            ksrftms            , & !real(r8), intent(in)    :: ksrftms(pcols)            ! Surface drag coefficient for turbulent mountain stress. > 0. [ kg/s/m2 ]
            qmincg             , & !real(r8), intent(in)    :: qmincg(ncnst)             ! Minimum constituent mixing ratios from cg fluxes
            fieldlist_wet      , & !type(vdiff_selector), intent(in) :: fieldlist        ! Array of flags selecting which fields to diffuse
            ptend_u            , & !real(r8), intent(inout) :: u(pcols,pver)             ! U wind. This input is the 'raw' input wind to PBL scheme without iterative provisional update. [ m/s
            ptend_v            , & !real(r8), intent(inout) :: v(pcols,pver)             ! V wind. This input is the 'raw' input wind to PBL scheme without iterative provisional update. [ m/s
            ptend_q            , & !real(r8), intent(inout) :: q(pcols,pver,ncnst)       ! Moisture and trace constituents [ kg/kg, #/kg ? ]
            ptend_s            , & !real(r8), intent(inout) :: dse(pcols,pver)           ! Dry static energy [ J/kg ]
            tautmsx            , & !real(r8), intent(inout) :: tauresx(pcols)            ! Input  : Reserved surface stress at previous time step
            tautmsy            , & !real(r8), intent(inout) :: tauresy(pcols)            ! Output : Reserved surface stress at current  time step
            dtk                , & !real(r8), intent(out)   :: dtk(pcols,pver)           ! T tendency from KE dissipation
            topflx             , & !real(r8), intent(out)   :: topflx(pcols)             ! Molecular heat flux at the top interface
            errstring          , & !character(128), intent(out) :: errstring             ! Output status
            tauresx            , & !real(r8), intent(out)   :: tautmsx(pcols)            ! Implicit zonal      turbulent mountain surface stress [ N/m2 = kg m/s /s/m2 ]
            tauresy            , & !real(r8), intent(out)   :: tautmsy(pcols)            ! Implicit meridional turbulent mountain surface stress [ N/m2 = kg m/s /s/m2 ]
            1                    )!, & !integer,  intent(in)    :: itaures                   ! Indicator determining whether 'tauresx,tauresy' is updated (1) or non-updated (0) in this subroutine
       !do_molec_diff     , & !
       !compute_molec_diff, & !
       !compute_molec_diff ) !integer,  external, optional :: compute_molec_diff   ! Constituent-independent moleculuar diffusivity routine
       kvm_in=kvm
       kvh_in=kvh

    END IF
    IF( errstring .NE. '' ) STOP 'call endrun( errstring )'

    IF( ANY( fieldlist_dry ) ) THEN

       IF( do_molec_diff ) THEN
          errstring = "Design flaw: dry vdiff not currently supported with molecular diffusion"
          STOP 'call endrun( errstring )'
       END IF

       CALL compute_vdiff( &
            lchnk              , &!integer,  intent(in)    :: lchnk
            pcols              , &!integer,  intent(in)    :: pcols
            pver               , &!integer,  intent(in)    :: pver
            pcnst              , &!integer,  intent(in)    :: ncnst
            ncol               , &!integer,  intent(in)    :: ncol                        ! Number of atmospheric columns
            TSK                , & 
            QSFC               , & 
            state_pmid         , &!real(r8), intent(in)    :: state_pmiddry(pcols,pver)          ! Mid-point pressures [ Pa ]
            state_pint         , &!real(r8), intent(in)    :: state_pintdry(pcols,pver+1)        ! Interface pressures [ Pa ]
            state_rpdel        , &!real(r8), intent(in)    :: state_rpdeldry(pcols,pver)         ! 1./pdel
            state_t            , &!real(r8), intent(in)    :: t(pcols,pver)                ! Temperature [ K ]
            ztodt              , &!real(r8), intent(in)    :: ztodt                        ! 2 delta-t [ s ]
            taux               , &!real(r8), intent(in)    :: taux(pcols)                ! Surface zonal      stress. Input u-momentum per unit time per unit area into the atmosphere [ N/m2 ]
            tauy               , &!real(r8), intent(in)    :: tauy(pcols)                ! Surface meridional stress. Input v-momentum per unit time per unit area into the atmosphere [ N/m2 ]
            shflx              , &!real(r8), intent(in)    :: shflx(pcols)                ! Surface sensible heat flux [ W/m2 ]
            cflx               , &!real(r8), intent(in)    :: cflx(pcols,ncnst)         ! Surface constituent flux [ kg/m2/s ]
            ntop               , &!integer,  intent(in)    :: ntop                        ! Top         interface level to which vertical diffusion is applied ( = 1 ).
            nbot               , &!integer,  intent(in)    :: nbot                        ! Bottom interface level to which vertical diffusion is applied ( = pver ).
            kvh                , &!real(r8), intent(inout) :: kvh(pcols,pver+1)         ! Eddy diffusivity for heat [ m2/s ]
            kvm                , &!real(r8), intent(inout) :: kvm(pcols,pver+1)         ! Eddy viscosity ( Eddy diffusivity for momentum ) [ m2/s ]
            kvq                , &!real(r8), intent(inout) :: kvq(pcols,pver+1)         ! Eddy diffusivity for constituents
            cgs                , &!real(r8), intent(inout) :: cgs(pcols,pver+1)         ! Counter-gradient star [ cg/flux ]
            cgh                , &!real(r8), intent(inout) :: cgh(pcols,pver+1)         ! Counter-gradient term for heat
            state_zi           , &!real(r8), intent(in)    :: zi(pcols,pver+1)          ! Interface heights [ m ]
            ksrftms            , &!real(r8), intent(in)    :: ksrftms(pcols)                ! Surface drag coefficient for turbulent mountain stress. > 0. [ kg/s/m2 ]
            qmincg             , &!real(r8), intent(in)    :: qmincg(ncnst)                ! Minimum constituent mixing ratios from cg fluxes
            fieldlist_dry      , &!type(vdiff_selector)   , intent(in) :: fieldlist        ! Array of flags selecting which fields to diffuse
            ptend_u            , &!real(r8), intent(inout) :: u(pcols,pver)                ! U wind. This input is the 'raw' input wind to PBL scheme without iterative provisional update. [ m/s
            ptend_v            , &!real(r8), intent(inout) :: v(pcols,pver)                ! V wind. This input is the 'raw' input wind to PBL scheme without iterative provisional update. [ m/s
            ptend_q            , &!real(r8), intent(inout) :: q(pcols,pver,ncnst)        ! Moisture and trace constituents [ kg/kg, #/kg ? ]
            ptend_s            , &!real(r8), intent(inout) :: dse(pcols,pver)                ! Dry static energy [ J/kg ]
            tautmsx            , &!real(r8), intent(inout) :: tauresx(pcols)                ! Input  : Reserved surface stress at previous time step
            tautmsy            , &!real(r8), intent(inout) :: tauresy(pcols)                ! Output : Reserved surface stress at current  time step
            dtk                , &!real(r8), intent(out)   :: dtk(pcols,pver)                ! T tendency from KE dissipation
            topflx             , &!real(r8), intent(out)   :: topflx(pcols)                ! Molecular heat flux at the top interface
            errstring          , &!character(128), intent(out) :: errstring                ! Output status
            tauresx            , &!real(r8), intent(out)   :: tautmsx(pcols)                ! Implicit zonal      turbulent mountain surface stress [ N/m2 = kg m/s /s/m2 ]
            tauresy            , &!real(r8), intent(out)   :: tautmsy(pcols)                ! Implicit meridional turbulent mountain surface stress [ N/m2 = kg m/s /s/m2 ]
            1                  )!, &!integer,  intent(in)    :: itaures                        ! Indicator determining whether 'tauresx,tauresy' is updated (1) or non-updated (0) in this subroutine
       !do_molec_diff      , &!
       !compute_molec_diff , &!
       !compute_molec_diff )!integer,  external, optional :: compute_molec_diff        ! Constituent-independent moleculuar diffusivity routine
            kvm_in=kvm
            kvh_in=kvh

       IF( errstring .NE. '' ) STOP 'call endrun( errstring )'

    END IF

    ! Store updated tauresx, tauresy in pbuf to use here on the next timestep

    !pbuf(tauresx_idx)%fld_ptr(1,1:ncol,1,lchnk,time_index) = tauresx(:ncol)
    !pbuf(tauresy_idx)%fld_ptr(1,1:ncol,1,lchnk,time_index) = tauresy(:ncol)
    !if( is_first_step() ) then
    !    do i = 1, pbuf_times
    !       pbuf(tauresx_idx)%fld_ptr(1,1:ncol,1,lchnk,i) = tauresx(:ncol)
    !       pbuf(tauresy_idx)%fld_ptr(1,1:ncol,1,lchnk,i) = tauresy(:ncol)
    !    end do
    !end if

    IF(MODAL_AERO)THEN

       ! Add the explicit surface fluxes to the lowest layer
       ! Modification : I should check whether this explicit adding is consistent with
       !                the treatment of other tracers.

       !tmp1(:ncol) = ztodt * gravit * state_rpdel(:ncol,pver)
       !do m = 1, ntot_amode
       !   l = numptr_amode(m)
       !   ptend_q(:ncol,pver,l) = ptend_q(:ncol,pver,l) + tmp1(:ncol) * cflx(:ncol,l)
       !   do lspec = 1, nspec_amode(m)
       !      l = lmassptr_amode(lspec,m)
       !      ptend_q(:ncol,pver,l) = ptend_q(:ncol,pver,l) + tmp1(:ncol) * cflx(:ncol,l)
       !   enddo
       !enddo

    END IF
    ! -------------------------------------------------------- !
    ! Diagnostics and output writing after applying PBL scheme !
    ! -------------------------------------------------------- !
    DO k=1,pver
       DO i=1,ncol
         sl(i,k)  = ptend_s(i,k) -   latvap           * ptend_q(i,k,ixcldliq) &
                                 - ( latvap + latice) * ptend_q(i,k,ixcldice)
         qt(i,k)  = ptend_q(i,k,1) + ptend_q(i,k,ixcldliq) &
                                   + ptend_q(i,k,ixcldice)
         slv(i,k) = sl(i,k) * ( 1.0_r8 + zvir*qt(i,k) ) 
      END DO
    END DO
    DO i=1,ncol
       slflx (i,1)  = 0.0_r8
       qtflx (i,1)  = 0.0_r8
       uflx  (i,1)  = 0.0_r8
       vflx  (i,1)  = 0.0_r8

       slflx_cg(i,1) = 0.0_r8
       qtflx_cg(i,1) = 0.0_r8
       uflx_cg (i,1) = 0.0_r8
       vflx_cg (i,1) = 0.0_r8
    END DO
    DO k = 2, pver
       DO i = 1, ncol
          rhoair     = state_pint(i,k) / ( rair * ( ( 0.5*(slv(i,k)+slv(i,k-1)) - gravit*state_zi(i,k))/cpair ) )
          slflx(i,k) = kvh(i,k) * &
               ( - rhoair*(sl(i,k-1)-sl(i,k))/(state_zm(i,k-1)-state_zm(i,k)) &
               + cgh(i,k) ) 
          IF(ixcldliq >1 .OR. ixcldice>1)THEN
             qtflx(i,k) = kvh(i,k) * &
                   ( - rhoair*(qt(i,k-1)-qt(i,k))/(state_zm(i,k-1)-state_zm(i,k)) &
                   + rhoair*(cflx(i,1)+cflx(i,ixcldliq)+cflx(i,ixcldice))*cgs(i,k) )
          ELSE
             qtflx(i,k) = kvh(i,k) * &
                  ( - rhoair*(qt(i,k-1)-qt(i,k))/(state_zm(i,k-1)-state_zm(i,k)) &
                  + rhoair*(cflx(i,1))*cgs(i,k) )
 
          END IF
          uflx(i,k)  = kvm(i,k) * &
               ( - rhoair*(ptend_u(i,k-1)-ptend_u(i,k))/(state_zm(i,k-1)-state_zm(i,k)))
          vflx(i,k)  = kvm(i,k) * &
               ( - rhoair*(ptend_v(i,k-1)-ptend_v(i,k))/(state_zm(i,k-1)-state_zm(i,k)))
          slflx_cg(i,k) = kvh(i,k) * cgh(i,k)
          IF(ixcldliq >1 .OR. ixcldice>1)THEN
             qtflx_cg(i,k) = kvh(i,k) * rhoair * ( cflx(i,1) + cflx(i,ixcldliq) + cflx(i,ixcldice) ) * cgs(i,k)
          ELSE
             qtflx_cg(i,k) = kvh(i,k) * rhoair * ( cflx(i,1) ) * cgs(i,k)
          END IF
          uflx_cg(i,k)  = 0._r8
          vflx_cg(i,k)  = 0._r8
       END DO
    END DO

    ! Modification : I should check whether slflx(:ncol,pverp) is correctly computed.
    !                Note also that 'tautotx' is explicit total stress, different from
    !                the ones that have been actually added into the atmosphere.
    DO i=1,ncol
       slflx(i,pver+1) = shflx(i)
       qtflx(i,pver+1) = cflx(i,1)
       uflx (i,pver+1)  = tautotx(i)
       vflx (i,pver+1)  = tautoty(i)

       slflx_cg (i,pver+1) = 0._r8
       qtflx_cg (i,pver+1) = 0._r8
       uflx_cg  (i,pver+1)  = 0._r8
       vflx_cg  (i,pver+1)  = 0._r8
    END DO
    ! --------------------------------------------------------------- !
    ! Convert the new profiles into vertical diffusion tendencies.    !
    ! Convert KE dissipative heat change into "temperature" tendency. !
    ! --------------------------------------------------------------- !
    !qm1
    DO m=1,ncnst
       DO k=1,pver
          DO i=1,ncol
             ptend_q(i,k,m)            = ( ptend_q(i,k,m) - qm1(i,k,m)) * rztodt
          END DO
       END DO
    END DO

    DO k = 1, pver
       DO i = 1, ncol
          ptend_s(i,k)              = ( ptend_s(i,k)   -  state_s(i,k) ) * rztodt
          ptend_u(i,k)              = ( ptend_u(i,k)   -  state_u(i,k) ) * rztodt
          ptend_v(i,k)              = ( ptend_v(i,k)   -  state_v(i,k) ) * rztodt
          !ptend_q(i,k,1)            = ( ptend_q(i,k,1)        - state_qv(i,k) ) * rztodt
          !ptend_q(i,k,ixcldliq)     = ( ptend_q(i,k,ixcldliq) - state_ql(i,k) ) * rztodt
          !ptend_q(i,k,ixcldice)     = ( ptend_q(i,k,ixcldice) - state_qi(i,k) ) * rztodt
          slten(i,k)                = ( sl(i,k)       - sl_prePBL(i,k) ) * rztodt
          qtten(i,k)                = ( qt(i,k)       - qt_prePBL(i,k) ) * rztodt
       END DO
    END DO
    ! ----------------------------------------------------------- !
    ! In order to perform 'pseudo-conservative varible diffusion' !
    ! perform the following two stages:                           !
    !                                                             !
    ! I.  Re-set (1) 'qvten' by 'qtten', and 'qlten = qiten = 0'  !
    !            (2) 'sten'  by 'slten', and                      !
    !            (3) 'qlten = qiten = 0'                          !
    !                                                             !
    ! II. Apply 'positive_moisture'                               !
    !                                                             !
    ! ----------------------------------------------------------- !

    IF( eddy_scheme .EQ. 'diag_TKE' .AND. do_pseudocon_diff ) THEN
       DO k = 1, pver
          DO i = 1, ncol
             ptend_q(i,k,1) = qtten(i,k)
             ptend_s(i,k)   = slten(i,k)
             ptend_q(i,k,ixcldliq) = 0._r8
             ptend_q(i,k,ixcldice) = 0._r8
             ptend_q(i,k,ixnumliq) = 0._r8
             ptend_q(i,k,ixnumice) = 0._r8
          END DO
       END DO 
       DO i = 1, ncol
          DO k = 1, pver
             qv_pro(i,k) = state_qv(i,k)         + ptend_q(i,k,1)              * ztodt       
             ql_pro(i,k) = state_ql(i,k)         + ptend_q(i,k,ixcldliq)       * ztodt
             qi_pro(i,k) = state_qi(i,k)         + ptend_q(i,k,ixcldice)       * ztodt              
             s_pro(i,k)  = state_s(i,k)          + ptend_s(i,k)                * ztodt
             t_pro(i,k)  = state_t(i,k)          + (1.0_r8/cpair)*ptend_s(i,k) * ztodt
          END DO
       END DO
       CALL positive_moisture( cpair, latvap, latvap+latice, ncol, pver, ztodt, qmin(1), qmin(2), qmin(3),    &
            state_pdel(:ncol,pver:1:-1), qv_pro(:ncol,pver:1:-1), ql_pro(:ncol,pver:1:-1), &
            qi_pro(:ncol,pver:1:-1), t_pro(:ncol,pver:1:-1), s_pro(:ncol,pver:1:-1),       &
            ptend_q(:ncol,pver:1:-1,1), ptend_q(:ncol,pver:1:-1,ixcldliq),                 &
            ptend_q(:ncol,pver:1:-1,ixcldice), ptend_s(:ncol,pver:1:-1) )

    END IF

    ! ----------------------------------------------------------------- !
    ! Re-calculate diagnostic output variables after vertical diffusion !
    ! ----------------------------------------------------------------- !
    DO k = 1, pver
       DO i = 1, ncol
         ! qv after PBL diffusion
         qv_aft_PBL(i,k)  =   state_qv(i,k)        + ptend_q(i,k,1)        * ztodt
         ! ql after PBL diffusion
         ql_aft_PBL(i,k)  =   state_ql(i,k)        + ptend_q(i,k,ixcldliq) * ztodt
         ! qi after PBL diffusion
         qi_aft_PBL(i,k)  =   state_qi(i,k)        + ptend_q(i,k,ixcldice) * ztodt
         ! s after PBL diffusion
         s_aft_PBL (i,k)  =   state_s(i,k)         + ptend_s(i,k)          * ztodt
         ! Temperature after PBL diffusion
         !          state_s(i,k) =cpair * state_t(i,k)+ gravit*state_zm(i,k)

         t_aftPBL  (i,k)  =  (s_aft_PBL(i,k) - gravit*state_zm(i,k) ) / cpair 
         ! u after PBL diffusion
         u_aft_PBL (i,k)  =  state_u(i,k)          + ptend_u(i,k)          * ztodt
         up1(i,k)= u_aft_PBL (i,k)
         ! v after PBL diffusion
         v_aft_PBL (i,k)  =  state_v(i,k)          + ptend_v(i,k)          * ztodt
         vp1(i,k)= v_aft_PBL (i,k)

       END DO
    END DO

    CALL aqsat( t_aftPBL, state_pmid, tem2, ftem, pcols, ncol, pver, 1, pver )

    DO k = 1, pver
       DO i = 1, ncol
          ! Saturation vapor pressure after PBL
          ftem_aftPBL(i,k) = qv_aft_PBL(i,k) / ftem(i,k) * 100._r8
          ! Temperature tendency by PBL diffusion
          tten(i,k)        = ( t_aftPBL(i,k)    - state_t(i,k) )              * rztodt     
          ! RH tendency by PBL diffusion
          rhten(i,k)       = ( ftem_aftPBL(i,k) - ftem_prePBL(i,k) )          * rztodt 
       END DO
    END DO
    !----------------------------------------------------------------------
    !**********************************************************************
    !----------------------------------------------------------------------

    !
    ! Convert the diffused fields back to diffusion tendencies.
    ! Add the diffusion tendencies to the cummulative physics tendencies,
    ! except for constituents. The diffused values of the constituents
    ! replace the input values.
    !
    !rztodt = 1.0_r8/ztodt
    DO k=1,pver
       DO i=1,ncol
          duv(i,k) =  ptend_u(i,k)*SIN( colrad(i)) ! (up1(i,k)*SIN( colrad(i)) - um1(i,k)*SIN( colrad(i)))*rztodt
          dvv(i,k) =  ptend_v(i,k)*SIN( colrad(i)) ! (vp1(i,k)*SIN( colrad(i)) - vm1(i,k)*SIN( colrad(i)))*rztodt

          dtv(i,k) =  tten(i,k)*state_exner(i,k)   ! (thp(i,k) - thm(i,k))*rztodt

       END DO
       DO m=1,pcnst
          DO i=1,ncol
             dqv(i,k,m) =  ptend_q(i,k,m) !(qp1(i,k,m) - qm1(i,k,m))*rztodt
          END DO
       END DO
    END DO


    is_first_step = .FALSE.

    RETURN
  END SUBROUTINE vertical_diffusion_tend

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  SUBROUTINE positive_moisture( cp, xlv, xls, ncol, mkx, dt, qvmin, qlmin, qimin, & 
       dp, qv, ql, qi, t, s, qvten, qlten, qiten, sten )
    ! ------------------------------------------------------------------------------- !
    ! If any 'ql < qlmin, qi < qimin, qv < qvmin' are developed in any layer,         !
    ! force them to be larger than minimum value by (1) condensating water vapor      !
    ! into liquid or ice, and (2) by transporting water vapor from the very lower     !
    ! layer. '2._r8' is multiplied to the minimum values for safety.                  !
    ! Update final state variables and tendencies associated with this correction.    !
    ! If any condensation happens, update (s,t) too.                                  !
    ! Note that (qv,ql,qi,t,s) are final state variables after applying corresponding !
    ! input tendencies.                                                               !
    ! Be careful the order of k : '1': near-surface layer, 'mkx' : top layer          ! 
    ! ------------------------------------------------------------------------------- !
    IMPLICIT NONE
    INTEGER,  INTENT(in)     :: ncol, mkx
    REAL(r8), INTENT(in)     :: cp, xlv, xls
    REAL(r8), INTENT(in)     :: dt, qvmin, qlmin, qimin
    REAL(r8), INTENT(in)     :: dp(ncol,mkx)
    REAL(r8), INTENT(inout)  :: qv(ncol,mkx), ql(ncol,mkx), qi(ncol,mkx), t(ncol,mkx), s(ncol,mkx)
    REAL(r8), INTENT(inout)  :: qvten(ncol,mkx), qlten(ncol,mkx), qiten(ncol,mkx), sten(ncol,mkx)
    INTEGER   i, k
    REAL(r8)  dql, dqi, dqv, sum, aa, dum 

    ! Modification : I should check whether this is exactly same as the one used in
    !                shallow convection and cloud macrophysics.

    DO i = 1, ncol
       DO k = mkx, 1, -1    ! From the top to the 1st (lowest) layer from the surface
          dql        = MAX(0._r8,1._r8*qlmin-ql(i,k))
          dqi        = MAX(0._r8,1._r8*qimin-qi(i,k))
          qlten(i,k) = qlten(i,k) +  dql/dt
          qiten(i,k) = qiten(i,k) +  dqi/dt
          qvten(i,k) = qvten(i,k) - (dql+dqi)/dt
          sten(i,k)  = sten(i,k)  + xlv * (dql/dt) + xls * (dqi/dt)
          ql(i,k)    = ql(i,k) +  dql
          qi(i,k)    = qi(i,k) +  dqi
          qv(i,k)    = qv(i,k) -  dql - dqi
          s(i,k)     = s(i,k)  +  xlv * dql + xls * dqi
          t(i,k)     = t(i,k)  + (xlv * dql + xls * dqi)/cp
          dqv        = MAX(0._r8,1._r8*qvmin-qv(i,k))
          qvten(i,k) = qvten(i,k) + dqv/dt
          qv(i,k)    = qv(i,k)    + dqv
          IF( k .NE. 1 ) THEN 
             qv(i,k-1)    = qv(i,k-1)    - dqv*dp(i,k)/dp(i,k-1)
             qvten(i,k-1) = qvten(i,k-1) - dqv*dp(i,k)/dp(i,k-1)/dt
          ENDIF
          qv(i,k) = MAX(qv(i,k),qvmin)
          ql(i,k) = MAX(ql(i,k),qlmin)
          qi(i,k) = MAX(qi(i,k),qimin)
       END DO
       ! Extra moisture used to satisfy 'qv(i,1)=qvmin' is proportionally 
       ! extracted from all the layers that has 'qv > 2*qvmin'. This fully
       ! preserves column moisture. 
       IF( dqv .GT. 1.e-20_r8 ) THEN
          sum = 0._r8
          DO k = 1, mkx
             IF( qv(i,k) .GT. 2._r8*qvmin ) sum = sum + qv(i,k)*dp(i,k)
          ENDDO
          aa = dqv*dp(i,1)/MAX(1.e-20_r8,sum)
          IF( aa .LT. 0.5_r8 ) THEN
             DO k = 1, mkx
                IF( qv(i,k) .GT. 2._r8*qvmin ) THEN
                   dum        = aa*qv(i,k)
                   qv(i,k)    = qv(i,k) - dum
                   qvten(i,k) = qvten(i,k) - dum/dt
                ENDIF
             ENDDO
          ELSE 
             WRITE(iulog,*) 'Full positive_moisture is impossible in vertical_diffusion'
          ENDIF
       ENDIF
    END DO
    RETURN

  END SUBROUTINE positive_moisture



  SUBROUTINE init_tms( kind, oro_in, karman_in, gravit_in, rair_in )

    INTEGER,  INTENT(in) :: kind   
    REAL(r8), INTENT(in) :: oro_in, karman_in, gravit_in, rair_in

    IF( kind .NE. r8 ) THEN
       WRITE(iulog,*) 'KIND of reals passed to init_tms -- exiting.'
       STOP 'compute_tms'
    ENDIF

    oroconst = oro_in

    RETURN
  END SUBROUTINE init_tms
  !============================================================================ !
  !                                                                             !
  !============================================================================ !


  SUBROUTINE init_eddy_diff( kind, pver, gravx, cpairx, rairx, zvirx, & 
       latvapx, laticex, ntop_eddy, nbot_eddy, vkx )
    !---------------------------------------------------------------- ! 
    ! Purpose:                                                        !
    ! Initialize time independent constants/variables of PBL package. !
    !---------------------------------------------------------------- !
    !use diffusion_solver, only: init_vdiff, vdiff_select
    !use cam_history,      only: outfld, addfld, phys_decomp
    IMPLICIT NONE
    ! --------- !
    ! Arguments !
    ! --------- !
    INTEGER,  INTENT(in) :: kind       ! Kind of reals being passed in
    INTEGER,  INTENT(in) :: pver       ! Number of vertical layers
    INTEGER,  INTENT(in) :: ntop_eddy  ! Top interface level to which eddy vertical diffusivity is applied ( = 1 )
    INTEGER,  INTENT(in) :: nbot_eddy  ! Bottom interface level to which eddy vertical diffusivity is applied ( = pver )
    REAL(r8), INTENT(in) :: gravx      ! Acceleration of gravity
    REAL(r8), INTENT(in) :: cpairx     ! Specific heat of dry air
    REAL(r8), INTENT(in) :: rairx      ! Gas constant for dry air
    REAL(r8), INTENT(in) :: zvirx      ! rh2o/rair - 1
    REAL(r8), INTENT(in) :: latvapx    ! Latent heat of vaporization
    REAL(r8), INTENT(in) :: laticex    ! Latent heat of fusion
    REAL(r8), INTENT(in) :: vkx        ! Von Karman's constant

    CHARACTER(128)       :: errstring  ! Error status for init_vdiff
    INTEGER              :: k          ! Vertical loop index

    IF( kind .NE. r8 ) THEN
       WRITE(iulog,*) 'wrong KIND of reals passed to init_diffusvity -- exiting.'
       STOP 'init_eddy_diff'
    ENDIF

    ! --------------- !
    ! Basic constants !
    ! --------------- !

    ntop_turb = ntop_eddy
    nbot_turb = nbot_eddy
    b123      = b1**(2._r8/3._r8)

    ! Set the square of the mixing lengths. Only for CAM3 HB PBL scheme.
    ! Not used for UW moist PBL. Used for free air eddy diffusivity.

    ALLOCATE(ml2(pver+1));ml2=0.0_r8
    ml2(1:ntop_turb) = 0._r8
    DO k = ntop_turb + 1, nbot_turb
       ml2(k) = 30.0_r8**2
    END DO
    ml2(nbot_turb+1:pver+1) = 0._r8

    ! Initialize diffusion solver module

    CALL init_vdiff(r8, 1, rair, g, fieldlist_wet, fieldlist_dry, errstring)

    ! Select the fields which will be diffused 

    IF(vdiff_select(fieldlist_wet,'s').NE.'')   WRITE(iulog,*) 'error: ', vdiff_select(fieldlist_wet,'s')
    IF(vdiff_select(fieldlist_wet,'q',1).NE.'') WRITE(iulog,*) 'error: ', vdiff_select(fieldlist_wet,'q',1)
    IF(vdiff_select(fieldlist_wet,'u').NE.'')   WRITE(iulog,*) 'error: ', vdiff_select(fieldlist_wet,'u')
    IF(vdiff_select(fieldlist_wet,'v').NE.'')   WRITE(iulog,*) 'error: ', vdiff_select(fieldlist_wet,'v')
    
    ! ------------------------------------------------------------------- !
    ! Writing outputs for detailed analysis of UW moist turbulence scheme !
    ! ------------------------------------------------------------------- !

    !call addfld('UW_errorPBL',      'm2/s',    1,      'A',  'Error function of UW PBL',                              phys_decomp )
    !call addfld('UW_n2',            's-2',     pver,   'A',  'Buoyancy Frequency, LI',                                phys_decomp )
    !call addfld('UW_s2',            's-2',     pver,   'A',  'Shear Frequency, LI',                                   phys_decomp )
    !call addfld('UW_ri',            'no',      pver,   'A',  'Interface Richardson Number, I',                        phys_decomp )
    !call addfld('UW_sfuh',          'no',      pver,   'A',  'Upper-Half Saturation Fraction, L',                     phys_decomp )
    !call addfld('UW_sflh',          'no',      pver,   'A',  'Lower-Half Saturation Fraction, L',                     phys_decomp )
    !call addfld('UW_sfi',           'no',      pver+1, 'A',  'Interface Saturation Fraction, I',                      phys_decomp )
    !call addfld('UW_cldn',          'no',      pver,   'A',  'Cloud Fraction, L',                                     phys_decomp )
    !call addfld('UW_qrl',           'g*W/m2',  pver,   'A',  'LW cooling rate, L',                                    phys_decomp )
    !call addfld('UW_ql',            'kg/kg',   pver,   'A',  'ql(LWC), L',                                            phys_decomp )
    !call addfld('UW_chu',           'g*kg/J',  pver+1, 'A',  'Buoyancy Coefficient, chu, I',                          phys_decomp )
    !call addfld('UW_chs',           'g*kg/J',  pver+1, 'A',  'Buoyancy Coefficient, chs, I',                          phys_decomp )
    !call addfld('UW_cmu',           'g/kg/kg', pver+1, 'A',  'Buoyancy Coefficient, cmu, I',                          phys_decomp )
    !call addfld('UW_cms',           'g/kg/kg', pver+1, 'A',  'Buoyancy Coefficient, cms, I',                          phys_decomp )    
    !call addfld('UW_tke',           'm2/s2',   pver+1, 'A',  'TKE, I',                                                phys_decomp )
    !call addfld('UW_wcap',          'm2/s2',   pver+1, 'A',  'Wcap, I',                                               phys_decomp )        
    !call addfld('UW_bprod',         'm2/s3',   pver+1, 'A',  'Buoyancy production, I',                                phys_decomp )
    !call addfld('UW_sprod',         'm2/s3',   pver+1, 'A',  'Shear production, I',                                   phys_decomp )    
    !call addfld('UW_kvh',           'm2/s',    pver+1, 'A',  'Eddy diffusivity of heat, I',                           phys_decomp )
    !call addfld('UW_kvm',           'm2/s',    pver+1, 'A',  'Eddy diffusivity of uv, I',                             phys_decomp )
    !call addfld('UW_pblh',          'm',       1,      'A',  'PBLH, 1',                                               phys_decomp )
    !call addfld('UW_pblhp',         'Pa',      1,      'A',  'PBLH pressure, 1',                                      phys_decomp )
    !call addfld('UW_tpert',         'K',       1,      'A',  'Convective T excess, 1',                                phys_decomp )
    !call addfld('UW_qpert',         'kg/kg',   1,      'A',  'Convective qt excess, I',                               phys_decomp )
    !call addfld('UW_wpert',         'm/s',     1,      'A',  'Convective W excess, I',                                phys_decomp )
    !call addfld('UW_ustar',         'm/s',     1,      'A',  'Surface Frictional Velocity, 1',                        phys_decomp )
    !call addfld('UW_tkes',          'm2/s2',   1,      'A',  'Surface TKE, 1',                                        phys_decomp )
    !call addfld('UW_minpblh',       'm',       1,      'A',  'Minimum PBLH, 1',                                       phys_decomp )
    !call addfld('UW_turbtype',      'no',      pver+1, 'A',  'Interface Turbulence Type, I',                          phys_decomp )    
    !call addfld('UW_kbase_o',       'no',      ncvmax, 'A',  'Initial CL Base Exterbal Interface Index, CL',          phys_decomp )
    !call addfld('UW_ktop_o',        'no',      ncvmax, 'A',  'Initial Top Exterbal Interface Index, CL',              phys_decomp )
    !call addfld('UW_ncvfin_o',      '#',       1,      'A',  'Initial Total Number of CL regimes, CL',                phys_decomp )
    !call addfld('UW_kbase_mg',      'no',      ncvmax, 'A',  'kbase after merging, CL',                               phys_decomp )
    !call addfld('UW_ktop_mg',       'no',      ncvmax, 'A',  'ktop after merging, CL',                                phys_decomp )
    !call addfld('UW_ncvfin_mg',     '#',       1,      'A',  'ncvfin after merging, CL',                              phys_decomp )
    !call addfld('UW_kbase_f',       'no',      ncvmax, 'A',  'Final kbase with SRCL, CL',                             phys_decomp )
    !call addfld('UW_ktop_f',        'no',      ncvmax, 'A',  'Final ktop with SRCL, CL',                              phys_decomp )
    !call addfld('UW_ncvfin_f',      '#',       1,      'A',  'Final ncvfin with SRCL, CL',                            phys_decomp )
    !call addfld('UW_wet',           'm/s',     ncvmax, 'A',  'Entrainment rate at CL top, CL',                        phys_decomp )
    !call addfld('UW_web',           'm/s',     ncvmax, 'A',  'Entrainment rate at CL base, CL',                       phys_decomp )
    !call addfld('UW_jtbu',          'm/s2',    ncvmax, 'A',  'Buoyancy jump across CL top, CL',                       phys_decomp )
    !call addfld('UW_jbbu',          'm/s2',    ncvmax, 'A',  'Buoyancy jump across CL base, CL',                      phys_decomp )
    !call addfld('UW_evhc',          'no',      ncvmax, 'A',  'Evaporative enhancement factor, CL',                    phys_decomp )
    !call addfld('UW_jt2slv',        'J/kg',    ncvmax, 'A',  'slv jump for evhc, CL',                                 phys_decomp )
    !call addfld('UW_n2ht',          's-2',     ncvmax, 'A',  'n2 at just below CL top interface, CL',                 phys_decomp )
    !call addfld('UW_n2hb',          's-2',     ncvmax, 'A',  'n2 at just above CL base interface',                    phys_decomp )
    !call addfld('UW_lwp',           'kg/m2',   ncvmax, 'A',  'LWP in the CL top layer, CL',                           phys_decomp )
    !call addfld('UW_optdepth',      'no',      ncvmax, 'A',  'Optical depth of the CL top layer, CL',                 phys_decomp )
    !call addfld('UW_radfrac',       'no',      ncvmax, 'A',  'Fraction of radiative cooling confined in the CL top',  phys_decomp )
    !call addfld('UW_radf',          'm2/s3',   ncvmax, 'A',  'Buoyancy production at the CL top by radf, I',          phys_decomp )        
    !call addfld('UW_wstar',         'm/s',     ncvmax, 'A',  'Convective velocity, Wstar, CL',                        phys_decomp )
    !call addfld('UW_wstar3fact',    'no',      ncvmax, 'A',  'Enhancement of wstar3 due to entrainment, CL',          phys_decomp )
    !call addfld('UW_ebrk',          'm2/s2',   ncvmax, 'A',  'CL-averaged TKE, CL',                                   phys_decomp )
    !call addfld('UW_wbrk',          'm2/s2',   ncvmax, 'A',  'CL-averaged W, CL',                                     phys_decomp )
    !call addfld('UW_lbrk',          'm',       ncvmax, 'A',  'CL internal thickness, CL',                             phys_decomp )
    !call addfld('UW_ricl',          'no',      ncvmax, 'A',  'CL-averaged Ri, CL',                                    phys_decomp )
    !call addfld('UW_ghcl',          'no',      ncvmax, 'A',  'CL-averaged gh, CL',                                    phys_decomp )
    !call addfld('UW_shcl',          'no',      ncvmax, 'A',  'CL-averaged sh, CL',                                    phys_decomp )
    !call addfld('UW_smcl',          'no',      ncvmax, 'A',  'CL-averaged sm, CL',                                    phys_decomp )
    !call addfld('UW_gh',            'no',      pver+1, 'A',  'gh at all interfaces, I',                               phys_decomp )
    !call addfld('UW_sh',            'no',      pver+1, 'A',  'sh at all interfaces, I',                               phys_decomp )
    !call addfld('UW_sm',            'no',      pver+1, 'A',  'sm at all interfaces, I',                               phys_decomp )
    !call addfld('UW_ria',           'no',      pver+1, 'A',  'ri at all interfaces, I',                               phys_decomp )
    !call addfld('UW_leng',          'm/s',     pver+1, 'A',  'Turbulence length scale, I',                            phys_decomp )

    RETURN

  END SUBROUTINE init_eddy_diff

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  SUBROUTINE init_vdiff( kind, ncnst, rair_in, gravit_in, fieldlist_wet, fieldlist_dry, errstring )

    INTEGER,              INTENT(in)  :: kind            ! Kind used for reals
    INTEGER,              INTENT(in)  :: ncnst           ! Number of constituents
    REAL(r8),             INTENT(in)  :: rair_in         ! Input gas constant for dry air
    REAL(r8),             INTENT(in)  :: gravit_in       ! Input gravititational acceleration
    TYPE(vdiff_selector), INTENT(out) :: fieldlist_wet   ! List of fields to be diffused using moist mixing ratio
    TYPE(vdiff_selector), INTENT(out) :: fieldlist_dry   ! List of fields to be diffused using dry   mixing ratio
    CHARACTER(128),       INTENT(out) :: errstring       ! Output status

    errstring = ''
    IF( kind .NE. r8 ) THEN
       WRITE(iulog,*) 'KIND of reals passed to init_vdiff -- exiting.'
       errstring = 'init_vdiff'
       RETURN
    ENDIF

    !rair   = rair_in     
    !gravit = gravit_in 

    ALLOCATE( fieldlist_wet%fields( 3 + ncnst ) )
    fieldlist_wet%fields(:) = .FALSE.

    ALLOCATE( fieldlist_dry%fields( 3 + ncnst ) )
    fieldlist_dry%fields(:) = .FALSE.

  END SUBROUTINE init_vdiff

  ! =============================================================================== !
  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !
  SUBROUTINE compute_eddy_diff( lchnk    , &!integer,  intent(in)    :: lchnk   
       pcols    , &!integer,  intent(in)    :: pcols                     ! Number of atmospheric columns [ # ]
       pver     , &!integer,  intent(in)    :: pver                      ! Number of atmospheric layers  [ # ]
       ncol     , &!integer,  intent(in)    :: ncol                      ! Number of atmospheric columns [ # ]
       mlsi     , &
       vcover    , &
       z0       , & 
       TSK      , & 
       QSFC     , & 
       bps      , & 
       t        , &!real(r8), intent(in)    :: t(pcols,pver)             ! Temperature [K]
       qv       , &!real(r8), intent(in)    :: qv(pcols,pver)            ! Water vapor  specific humidity [ kg/kg ]
       ztodt    , &!real(r8), intent(in)    :: ztodt                     ! Physics integration time step 2 delta-t [ s ]
       ql       , &!real(r8), intent(in)    :: ql(pcols,pver)            ! Liquid water specific humidity [ kg/kg ]
       qi       , &!real(r8), intent(in)    :: qi(pcols,pver)            ! Ice specific humidity [ kg/kg ]
       rpdel    , &!real(r8), intent(in)    :: rpdel(pcols,pver)         ! 1./pdel where 'pdel' is thickness of the layer [ Pa ]
       cldn     , &!real(r8), intent(in)    :: cldn(pcols,pver)          ! Stratiform cloud fraction [ fraction ]
       qrl      , &!real(r8), intent(in)    :: qrl(pcols,pver)           ! LW cooling rate
       wsedl    , &!real(r8), intent(in)    :: wsedl(pcols,pver)         ! Sedimentation velocity of liquid stratus cloud droplet [ m/s ]
       z        , &!real(r8), intent(in)    :: z(pcols,pver)             ! Layer mid-point height above surface [ m ]
       zi       , &!real(r8), intent(in)    :: zi(pcols,pver+1)          ! Interface height above surface [ m ]
       pmid     , &!real(r8), intent(in)    :: pmid(pcols,pver)          ! Layer mid-point pressure [ Pa ]
       pi       , &!real(r8), intent(in)    :: pi(pcols,pver+1)          ! Interface pressure [ Pa ]
       u        , &!real(r8), intent(in)    :: u(pcols,pver)             ! Zonal velocity [ m/s ]
       v        , &!real(r8), intent(in)    :: v(pcols,pver)             ! Meridional velocity [ m/s ]
       taux     , &!real(r8), intent(in)    :: taux(pcols)               ! Zonal wind stress at surface [ N/m2 ]
       tauy     , &!real(r8), intent(in)    :: tauy(pcols)               ! Meridional wind stress at surface [ N/m2 ]
       shflx    , &!real(r8), intent(in)    :: shflx(pcols)              ! Sensible heat flux at surface [ unit ? ]
       qflx     , &!real(r8), intent(in)    :: qflx(pcols)               ! Water vapor flux at surface [ unit ? ]
       wstarent , &!logical,  intent(in)    :: wstarent                  ! .true. means use the 'wstar' entrainment closure. 
       nturb    , &!integer,  intent(in)    :: nturb                     ! Number of iteration steps for calculating eddy diffusivity [ # ]
       ustar    , &!real(r8), intent(out)   :: ustar(pcols)              ! Surface friction velocity [ m/s ]
       pblh     , &!real(r8), intent(out)   :: pblh(pcols)               ! PBL top height [ m ]
       landfrac   ,&
       kvm_in   , &!real(r8), intent(in)    :: kvm_in(pcols,pver+1)      ! kvm saved from last timestep [ m2/s ]
       kvh_in   , &!real(r8), intent(in)    :: kvh_in(pcols,pver+1)      ! kvh saved from last timestep [ m2/s ]
       kvm_out  , &!real(r8), intent(out)   :: kvm_out(pcols,pver+1)     ! Eddy diffusivity for momentum [ m2/s ]
       kvh_out  , &!real(r8), intent(out)   :: kvh_out(pcols,pver+1)     ! Eddy diffusivity for heat [ m2/s ]
       kvq      , &!real(r8), intent(out)   :: kvq(pcols,pver+1)         ! Eddy diffusivity for constituents, moisture and tracers [ m2/s ] (note not having '_out')
       cgh      , &!real(r8), intent(out)   :: cgh(pcols,pver+1)         ! Counter-gradient term for heat [ J/kg/m ]
       cgs      , &!real(r8), intent(out)   :: cgs(pcols,pver+1)         ! Counter-gradient star [ cg/flux ]
       tpert    , &!real(r8), intent(out)   :: tpert(pcols)              ! Convective temperature excess [ K ]
       qpert    , &!real(r8), intent(out)   :: qpert(pcols)              ! Convective humidity excess [ kg/kg ]
       wpert    , &!real(r8), intent(out)   :: wpert(pcols)              ! Turbulent velocity excess [ m/s ]
       tke2      , &!real(r8), intent(inout)   :: tke2(pcols,pver+1)         ! Turbulent kinetic energy [ m2/s2 ]
       bprod    , &!real(r8), intent(out)   :: bprod(pcols,pver+1)       ! Buoyancy production [ m2/s3 ] 
       sprod    , &!real(r8), intent(out)   :: sprod(pcols,pver+1)       ! Shear production [ m2/s3 ] 
       sfi      , &!real(r8), intent(out)   :: sfi(pcols,pver+1)         ! Interfacial layer saturation fraction [ fraction ]
       kvinit   , &!logical, intent(in  )   :: kvinit                    ! Tell compute_eddy_diff/ caleddy to initialize kvh, kvm (uses kvf)
       tauresx  , &!real(r8), intent(inout) :: tauresx(pcols)            ! Residual stress to be added in vdiff to correct for turb
       tauresy  , &!real(r8), intent(inout) :: tauresy(pcols)            ! Stress mismatch between sfc and atm accumulated in prior timesteps
       ksrftms  , &!real(r8), intent(in)    :: ksrftms(pcols)            ! Surface drag coefficient of turbulent mountain stress [ unit ? ]
       ipbl     , &!real(r8), intent(out)   :: ipbl(pcols)               ! If 1, PBL is CL, while if 0, PBL is STL.
       kpblh    , &!real(r8), intent(out)   :: kpblh(pcols)              ! Layer index containing PBL top within or at the base interface
       wstarPBL , &!real(r8), intent(out)   :: wstarPBL(pcols)           ! Convective velocity within PBL [ m/s ]
       turbtype , &!real(r8), intent(out)   :: turbtype(pcols,pver+1)    ! Turbulence type identifier at all interfaces [ no unit ]
       sm_aw    , &!real(r8), intent(out)   :: sm_aw(pcols,pver+1)       ! Normalized Galperin instability function for momentum [ no unit ]
       ri         )!REAL(r8), INTENT(out)   :: ri(pcols,pver)            ! Richardson number, 'n2/s2', defined at interfaces except surface [ s-2 ]

    !-------------------------------------------------------------------- ! 
    ! Purpose: Interface to compute eddy diffusivities.                   !
    !          Eddy diffusivities are calculated in a fully implicit way  !
    !          through iteration process.                                 !   
    ! Author:  Sungsu Park. August. 2006.                                 !
    !                       May.    2008.                                 !
    !-------------------------------------------------------------------- !

    !  use diffusion_solver, only: compute_vdiff
    !  use cam_history,      only: outfld, addfld, phys_decomp
    ! use physics_types,    only: physics_state
    !  use phys_debug_util,  only: phys_debug_col
    !  use time_manager,     only: is_first_step, get_nstep

    IMPLICIT NONE

    ! type(physics_state)     :: state                     ! Physics state variables

    ! --------------- !
    ! Input Variables !
    ! --------------- ! 

    INTEGER,  INTENT(in)    :: lchnk   
    INTEGER,  INTENT(in)    :: pcols                     ! Number of atmospheric columns [ # ]
    INTEGER,  INTENT(in)    :: pver                      ! Number of atmospheric layers  [ # ]
    INTEGER,  INTENT(in)    :: ncol                      ! Number of atmospheric columns [ # ]
    INTEGER,  INTENT(in)    :: nturb                     ! Number of iteration steps for calculating eddy diffusivity [ # ]
    LOGICAL,  INTENT(in)    :: wstarent                  ! .true. means use the 'wstar' entrainment closure. 
    LOGICAL,  INTENT(in)    :: kvinit                    ! 'true' means time step = 1 : used for initializing kvh, kvm (uses kvf or zero)
    REAL(r8), INTENT(in)    :: ztodt                     ! Physics integration time step 2 delta-t [ s ]
    REAL(r8), INTENT(in)    :: z0(pcols) 
    INTEGER(KIND=i8), INTENT(in)    :: mlsi (pcols) 
    REAL(r8), INTENT(in)    :: vcover (pcols)  
    REAL(r8), INTENT(in)    :: TSK(pcols)  
    REAL(r8), INTENT(in)    :: QSFC(pcols)  
    REAL(r8), INTENT(in)    :: bps   (pcols,pver)                   
    REAL(r8), INTENT(in)    :: t(pcols,pver)             ! Temperature [K]
    REAL(r8), INTENT(in)    :: qv(pcols,pver)            ! Water vapor  specific humidity [ kg/kg ]
    REAL(r8), INTENT(in)    :: ql(pcols,pver)            ! Liquid water specific humidity [ kg/kg ]
    REAL(r8), INTENT(in)    :: qi(pcols,pver)            ! Ice specific humidity [ kg/kg ]
    !    real(r8), intent(in)    :: s(pcols,pver)             ! Dry static energy [ J/kg ]
    REAL(r8), INTENT(in)    :: rpdel(pcols,pver)         ! 1./pdel where 'pdel' is thickness of the layer [ Pa ]
    REAL(r8), INTENT(in)    :: cldn(pcols,pver)          ! Stratiform cloud fraction [ fraction ]
    REAL(r8), INTENT(in)    :: qrl(pcols,pver)           ! LW cooling rate
    REAL(r8)    :: swh(pcols,pver)    
    REAL(r8), INTENT(in)    :: wsedl(pcols,pver)         ! Sedimentation velocity of liquid stratus cloud droplet [ m/s ]
    REAL(r8), INTENT(in)    :: z(pcols,pver)             ! Layer mid-point height above surface [ m ]
    REAL(r8), INTENT(in)    :: zi(pcols,pver+1)          ! Interface height above surface [ m ]
    REAL(r8), INTENT(in)    :: pmid(pcols,pver)          ! Layer mid-point pressure [ Pa ]
    REAL(r8), INTENT(in)    :: pi(pcols,pver+1)          ! Interface pressure [ Pa ]
    REAL(r8), INTENT(in)    :: u(pcols,pver)             ! Zonal velocity [ m/s ]
    REAL(r8), INTENT(in)    :: v(pcols,pver)             ! Meridional velocity [ m/s ]
    REAL(r8), INTENT(in)    :: taux(pcols)               ! Zonal wind stress at surface [ N/m2 ]
    REAL(r8), INTENT(in)    :: tauy(pcols)               ! Meridional wind stress at surface [ N/m2 ]
    REAL(r8), INTENT(in)    :: shflx(pcols)              ! Sensible heat flux at surface [ unit ? ]
    REAL(r8), INTENT(in)    :: qflx(pcols)               ! Water vapor flux at surface [ unit ? ]
    REAL(r8), INTENT(in)    :: kvm_in(pcols,pver+1)      ! kvm saved from last timestep [ m2/s ]
    REAL(r8), INTENT(in)    :: kvh_in(pcols,pver+1)      ! kvh saved from last timestep [ m2/s ]
    REAL(r8), INTENT(in)    :: ksrftms(pcols)            ! Surface drag coefficient of turbulent mountain stress [ unit ? ]
    REAL(r8), INTENT(in)    :: landfrac(pcols)            ! S
    ! ---------------- !
    ! Output Variables !
    ! ---------------- ! 

    REAL(r8), INTENT(out)   :: kvm_out(pcols,pver+1)     ! Eddy diffusivity for momentum [ m2/s ]
    REAL(r8), INTENT(out)   :: kvh_out(pcols,pver+1)     ! Eddy diffusivity for heat [ m2/s ]
    REAL(r8), INTENT(out)   :: kvq(pcols,pver+1)         ! Eddy diffusivity for constituents, moisture and tracers [ m2/s ] (note not having '_out')
    REAL(r8), INTENT(out)   :: ustar(pcols)              ! Surface friction velocity [ m/s ]
    REAL(r8), INTENT(out)   :: pblh(pcols)               ! PBL top height [ m ]
    REAL(r8), INTENT(out)   :: cgh(pcols,pver+1)         ! Counter-gradient term for heat [ J/kg/m ]
    REAL(r8), INTENT(out)   :: cgs(pcols,pver+1)         ! Counter-gradient star [ cg/flux ]
    REAL(r8), INTENT(out)   :: tpert(pcols)              ! Convective temperature excess [ K ]
    REAL(r8), INTENT(out)   :: qpert(pcols)              ! Convective humidity excess [ kg/kg ]
    REAL(r8), INTENT(out)   :: wpert(pcols)              ! Turbulent velocity excess [ m/s ]
    REAL(r8), INTENT(inout)   :: tke2(pcols,pver+1)         ! Turbulent kinetic energy [ m2/s2 ]
    REAL(r8), INTENT(out)   :: bprod(pcols,pver+1)       ! Buoyancy production [ m2/s3 ] 
    REAL(r8), INTENT(out)   :: sprod(pcols,pver+1)       ! Shear production [ m2/s3 ] 
    REAL(r8), INTENT(out)   :: sfi(pcols,pver+1)         ! Interfacial layer saturation fraction [ fraction ]
    REAL(r8), INTENT(out)   :: turbtype(pcols,pver+1)    ! Turbulence type identifier at all interfaces [ no unit ]
    REAL(r8), INTENT(out)   :: sm_aw(pcols,pver+1)       ! Normalized Galperin instability function for momentum [ no unit ]
    ! This is 1 when neutral condition (Ri=0), 4.964 for maximum unstable case, and 0 when Ri > Ricrit=0.19. 
    REAL(r8), INTENT(out)   :: ipbl(pcols)               ! If 1, PBL is CL, while if 0, PBL is STL.
    REAL(r8), INTENT(out)   :: kpblh(pcols)              ! Layer index containing PBL top within or at the base interface
    REAL(r8), INTENT(out)   :: wstarPBL(pcols)           ! Convective velocity within PBL [ m/s ]
    REAL(r8), INTENT(out)   :: ri(pcols,pver)            ! Richardson number, 'n2/s2', defined at interfaces except surface [ s-2 ]
    INTEGER   :: kpbl(pcols)     

    ! ---------------------- !
    ! Input-Output Variables !
    ! ---------------------- ! 

    REAL(r8), INTENT(inout) :: tauresx(pcols)            ! Residual stress to be added in vdiff to correct for turb
    REAL(r8), INTENT(inout) :: tauresy(pcols)            ! Stress mismatch between sfc and atm accumulated in prior timesteps

    ! --------------- !
    ! Local Variables !
    ! --------------- !

    INTEGER                    icol
    INTEGER                    i, k, iturb, status,m
    !integer,  external      :: qsat
    CHARACTER(128)          :: errstring                 ! Error status for compute_vdiff

    REAL(r8)                :: tautotx(pcols)            ! Total stress including tms
    REAL(r8)                :: tautoty(pcols)            ! Total stress including tms
    REAL(r8)                :: kvf(pcols,pver+1)         ! Free atmospheric eddy diffusivity [ m2/s ]
    REAL(r8)                :: kvm(pcols,pver+1)         ! Eddy diffusivity for momentum [ m2/s ]
    REAL(r8)                :: kvh(pcols,pver+1)         ! Eddy diffusivity for heat [ m2/s ]
    REAL(r8)                :: kvm_preo(pcols,pver+1)    ! Eddy diffusivity for momentum [ m2/s ]
    REAL(r8)                :: kvh_preo(pcols,pver+1)    ! Eddy diffusivity for heat [ m2/s ]
    REAL(r8)                :: kvm_pre(pcols,pver+1)     ! Eddy diffusivity for momentum [ m2/s ]
    REAL(r8)                :: kvh_pre(pcols,pver+1)     ! Eddy diffusivity for heat [ m2/s ]
    REAL(r8)                :: errorPBL(pcols)           ! Error function showing whether PBL produced convergent solution or not. [ unit ? ]
    REAL(r8)                :: s2(pcols,pver)            ! Shear squared, defined at interfaces except surface [ s-2 ]
    REAL(r8)                :: n2(pcols,pver)            ! Buoyancy frequency, defined at interfaces except surface [ s-2 ]
    REAL(r8)                :: pblhp(pcols)              ! PBL top pressure [ Pa ]
    REAL(r8)                :: minpblh(pcols)            ! Minimum PBL height based on surface stress

    REAL(r8)                :: qt(pcols,pver)            ! Total specific humidity [ kg/kg ]
    REAL(r8)                :: sfuh(pcols,pver)          ! Saturation fraction in upper half-layer [ fraction ]
    REAL(r8)                :: sflh(pcols,pver)          ! Saturation fraction in lower half-layer [ fraction ]
    REAL(r8)                :: sl(pcols,pver)            ! Liquid water static energy [ J/kg ]
    REAL(r8)                :: slv(pcols,pver)           ! Liquid water virtual static energy [ J/kg ]
    REAL(r8)                :: slslope(pcols,pver)       ! Slope of 'sl' in each layer
    REAL(r8)                :: qtslope(pcols,pver)       ! Slope of 'qt' in each layer
    REAL(r8)                :: rrho(pcols)               ! Density at the lowest layer
    REAL(r8)                :: qvfd(pcols,pver)          ! Specific humidity for diffusion [ kg/kg ]
    REAL(r8)                :: tfd(pcols,pver)           ! Temperature for diffusion [ K ]
    REAL(r8)                :: slfd(pcols,pver)          ! Liquid static energy [ J/kg ]
    REAL(r8)                :: qtfd(pcols,pver)          ! Total specific humidity [ kg/kg ] 
    REAL(r8)                :: qlfd(pcols,pver)          ! Liquid water specific humidity for diffusion [ kg/kg ]
    REAL(r8)                :: ufd(pcols,pver)           ! U-wind for diffusion [ m/s ]
    REAL(r8)                :: vfd(pcols,pver)           ! V-wind for diffusion [ m/s ]

    ! Buoyancy coefficients : w'b' = ch * w'sl' + cm * w'qt'

    REAL(r8)                :: chu(pcols,pver+1)         ! Heat buoyancy coef for dry states, defined at each interface, finally.
    REAL(r8)                :: chs(pcols,pver+1)         ! Heat buoyancy coef for sat states, defined at each interface, finally. 
    REAL(r8)                :: cmu(pcols,pver+1)         ! Moisture buoyancy coef for dry states, defined at each interface, finally.
    REAL(r8)                :: cms(pcols,pver+1)         ! Moisture buoyancy coef for sat states, defined at each interface, finally. 

    REAL(r8)                :: jnk1d(pcols)
    REAL(r8)                :: jnk2d(pcols,pver+1)  
    REAL(r8)                :: zero(pcols)
    REAL(r8)                :: zero2d(pcols,pver+1)
    REAL(r8)                :: es(1)                     ! Saturation vapor pressure
    REAL(r8)                :: qs(1)                     ! Saturation specific humidity
    REAL(r8)                :: gam(1)                    ! (L/cp)*dqs/dT
    REAL(r8)                :: ep2, templ(1), temps(1)

    ! ------------------------------- !
    ! Variables for diagnostic output !
    ! ------------------------------- !

    REAL(r8)                 :: tke(pcols,pver+1)         ! Turbulent kinetic energy [ m2/s2 ]
    REAL(r8)                :: tkes(pcols)               ! TKE at surface interface [ m2/s2 ]
    REAL(r8)                :: kbase_o(pcols,ncvmax)     ! Original external base interface index of CL from 'exacol'
    REAL(r8)                :: ktop_o(pcols,ncvmax)      ! Original external top  interface index of CL from 'exacol'
    REAL(r8)                :: ncvfin_o(pcols)           ! Original number of CLs from 'exacol'
    REAL(r8)                :: kbase_mg(pcols,ncvmax)    ! 'kbase' after extending-merging from 'zisocl'
    REAL(r8)                :: ktop_mg(pcols,ncvmax)     ! 'ktop' after extending-merging from 'zisocl'
    REAL(r8)                :: ncvfin_mg(pcols)          ! 'ncvfin' after extending-merging from 'zisocl'
    REAL(r8)                :: kbase_f(pcols,ncvmax)     ! Final 'kbase' after extending-merging & including SRCL
    REAL(r8)                :: ktop_f(pcols,ncvmax)      ! Final 'ktop' after extending-merging & including SRCL
    REAL(r8)                :: ncvfin_f(pcols)           ! Final 'ncvfin' after extending-merging & including SRCL
    REAL(r8)                :: wet(pcols,ncvmax)         ! Entrainment rate at the CL top  [ m/s ] 
    REAL(r8)                :: web(pcols,ncvmax)         ! Entrainment rate at the CL base [ m/s ]. Set to zero if CL is based at surface.
    REAL(r8)                :: jtbu(pcols,ncvmax)        ! Buoyancy jump across the CL top  [ m/s2 ]  
    REAL(r8)                :: jbbu(pcols,ncvmax)        ! Buoyancy jump across the CL base [ m/s2 ]  
    REAL(r8)                :: evhc(pcols,ncvmax)        ! Evaporative enhancement factor at the CL top
    REAL(r8)                :: jt2slv(pcols,ncvmax)      ! Jump of slv ( across two layers ) at CL top used only for evhc [ J/kg ]
    REAL(r8)                :: n2ht(pcols,ncvmax)        ! n2 defined at the CL top  interface but using sfuh(kt)   instead of sfi(kt) [ s-2 ] 
    REAL(r8)                :: n2hb(pcols,ncvmax)        ! n2 defined at the CL base interface but using sflh(kb-1) instead of sfi(kb) [ s-2 ]
    REAL(r8)                :: lwp(pcols,ncvmax)         ! LWP in the CL top layer [ kg/m2 ]
    REAL(r8)                :: opt_depth(pcols,ncvmax)   ! Optical depth of the CL top layer
    REAL(r8)                :: radinvfrac(pcols,ncvmax)  ! Fraction of radiative cooling confined in the top portion of CL top layer
    REAL(r8)                :: radf(pcols,ncvmax)        ! Buoyancy production at the CL top due to LW radiative cooling [ m2/s3 ]
    REAL(r8)                :: wstar(pcols,ncvmax)       ! Convective velocity in each CL [ m/s ]
    REAL(r8)                :: wstar3fact(pcols,ncvmax)  ! Enhancement of 'wstar3' due to entrainment (inverse) [ no unit ]
    REAL(r8)                :: ebrk(pcols,ncvmax)        ! Net mean TKE of CL including entrainment effect [ m2/s2 ]
    REAL(r8)                :: wbrk(pcols,ncvmax)        ! Net mean normalized TKE (W) of CL, 'ebrk/b1' including entrainment effect [ m2/s2 ]
    REAL(r8)                :: lbrk(pcols,ncvmax)        ! Energetic internal thickness of CL [m]
    REAL(r8)                :: ricl(pcols,ncvmax)        ! CL internal mean Richardson number
    REAL(r8)                :: ghcl(pcols,ncvmax)        ! Half of normalized buoyancy production of CL
    REAL(r8)                :: shcl(pcols,ncvmax)        ! Galperin instability function of heat-moisture of CL
    REAL(r8)                :: smcl(pcols,ncvmax)        ! Galperin instability function of mementum of CL
    REAL(r8)                :: ghi(pcols,pver+1)         ! Half of normalized buoyancy production at all interfaces
    REAL(r8)                :: shi(pcols,pver+1)         ! Galperin instability function of heat-moisture at all interfaces
    REAL(r8)                :: smi(pcols,pver+1)         ! Galperin instability function of heat-moisture at all interfaces
    REAL(r8)                :: rii(pcols,pver+1)         ! Interfacial Richardson number defined at all interfaces
    REAL(r8)                :: lengi(pcols,pver+1)       ! Turbulence length scale at all interfaces [ m ]
    REAL(r8)                :: wcap(pcols,pver+1)        ! Normalized TKE at all interfaces [ m2/s2 ]
    REAL(r8)                :: pcfm(1:pcols)
    REAL(r8)                :: pcfh(1:pcols)
    REAL(r8)                :: ycoefm  (pcols,pver+1)   !  REAL(KIND=r8), INTENT(OUT):: ycoefh (klon,klev)
    REAL(r8)                :: ycoefh  (pcols,pver+1)   !  REAL(KIND=r8), INTENT(OUT):: ycoefm (klon,klev)
    REAL(r8)                :: hpbl(1:pcols)
    !
    !  iflag_pbl doit valoir entre 6 et 9
    !      l=6, on prend  systematiquement une longueur d'equilibre
    !    iflag_pbl=6 : MY 2.0
    !    iflag_pbl=7 : MY 2.0.Fournier
    !    iflag_pbl=8 : MY 2.5
    !    iflag_pbl>=9 : MY 2.5 avec diffusion verticale
    !.......................................................................
    INTEGER, PARAMETER :: iflag_pbl=11

    ! ---------- !
    ! Initialize !
    ! ---------- !
    DO k=1,pver+1
       DO i=1,pcols
          kvm_out (i,k)= 0.0_r8! Eddy diffusivity for momentum [ m2/s ]
          kvh_out (i,k)= 0.0_r8! Eddy diffusivity for heat [ m2/s ]
          kvq     (i,k)= 0.0_r8! Eddy diffusivity for constituents, moisture and tracers [ m2/s ] (note not having '_out')
          cgh     (i,k)= 0.0_r8! Counter-gradient term for heat [ J/kg/m ]
          cgs     (i,k)= 0.0_r8! Counter-gradient star [ cg/flux ]
          tke     (i,k)= 0.0_r8! Turbulent kinetic energy [ m2/s2 ]
          bprod   (i,k)= 0.0_r8! Buoyancy production [ m2/s3 ] 
          sprod   (i,k)= 0.0_r8! Shear production [ m2/s3 ] 
          sfi     (i,k)= 0.0_r8! Interfacial layer saturation fraction [ fraction ]
          turbtype(i,k)= 0.0_r8! Turbulence type identifier at all interfaces [ no unit ]
          sm_aw   (i,k)= 0.0_r8! Normalized Galperin instability function for momentum [ no unit ]
          kvf(i,k) =0.0_r8       ! Free atmospheric eddy diffusivity [ m2/s ]
          kvm(i,k) =0.0_r8       ! Eddy diffusivity for momentum [ m2/s ]
          kvh(i,k) =0.0_r8       ! Eddy diffusivity for heat [ m2/s ]
          kvm_preo(i,k) =0.0_r8  ! Eddy diffusivity for momentum [ m2/s ]
          kvh_preo(i,k) =0.0_r8  ! Eddy diffusivity for heat [ m2/s ]
          kvm_pre(i,k) =0.0_r8   ! Eddy diffusivity for momentum [ m2/s ]
          kvh_pre(i,k) =0.0_r8   ! Eddy diffusivity for heat [ m2/s ]
          chu(i,k) =0.0_r8       ! Heat buoyancy coef for dry states, defined at each interface, finally.
          chs(i,k) =0.0_r8       ! Heat buoyancy coef for sat states, defined at each interface, finally. 
          cmu(i,k) =0.0_r8       ! Moisture buoyancy coef for dry states, defined at each interface, finally.
          cms(i,k) =0.0_r8       ! Moisture buoyancy coef for sat states, defined at each interface, finally. 
          jnk2d(i,k) =0.0_r8
          zero2d(i,k) =0.0_r8
          ghi(i,k) =0.0_r8       ! Half of normalized buoyancy production at all interfaces
          shi(i,k) =0.0_r8       ! Galperin instability function of heat-moisture at all interfaces
          smi(i,k) =0.0_r8       ! Galperin instability function of heat-moisture at all interfaces
          rii(i,k) =0.0_r8       ! Interfacial Richardson number defined at all interfaces
          lengi(i,k) =0.0_r8     ! Turbulence length scale at all interfaces [ m ]
          wcap(i,k) =0.0_r8      ! Normalized TKE at all interfaces [ m2/s2 ]

       END DO
    END DO

    DO k=1,pver
       DO i=1,pcols
          ycoefm (i,k)= 0.0_r8!
          ycoefh (i,k)= 0.0_r8!
          qt(i,k) =0.0_r8            ! Total specific humidity [ kg/kg ]
          sfuh(i,k) =0.0_r8          ! Saturation fraction in upper half-layer [ fraction ]
          sflh(i,k) =0.0_r8          ! Saturation fraction in lower half-layer [ fraction ]
          sl(i,k) =0.0_r8            ! Liquid water static energy [ J/kg ]
          slv(i,k) =0.0_r8           ! Liquid water virtual static energy [ J/kg ]
          slslope(i,k) =0.0_r8       ! Slope of 'sl' in each layer
          qtslope(i,k) =0.0_r8       ! Slope of 'qt' in each layer
          qvfd(i,k) =0.0_r8          ! Specific humidity for diffusion [ kg/kg ]
          tfd(i,k) =0.0_r8           ! Temperature for diffusion [ K ]
          slfd(i,k) =0.0_r8          ! Liquid static energy [ J/kg ]
          qtfd(i,k) =0.0_r8          ! Total specific humidity [ kg/kg ] 
          qlfd(i,k) =0.0_r8          ! Liquid water specific humidity for diffusion [ kg/kg ]
          ufd(i,k) =0.0_r8           ! U-wind for diffusion [ m/s ]
          vfd(i,k) =0.0_r8           ! V-wind for diffusion [ m/s ]
          s2(i,k) =0.0_r8            ! Shear squared, defined at interfaces except surface [ s-2 ]
          n2(i,k) =0.0_r8            ! Buoyancy frequency, defined at interfaces except surface [ s-2 ]
          ri(i,k)= 0.0_r8! Richardson number, 'n2/s2', defined at interfaces except surface [ s-2 ]
       END DO
    END DO


    DO i=1,pcols
       ustar    (i)= 0.0_r8! Surface friction velocity [ m/s ]
       pblh     (i)= 0.0_r8! PBL top height [ m ]
       tpert    (i)= 0.0_r8! Convective temperature excess [ K ]
       qpert    (i)= 0.0_r8! Convective humidity excess [ kg/kg ]
       wpert    (i)= 0.0_r8! Turbulent velocity excess [ m/s ]
       ! This is 1 when neutral condition (Ri=0), 4.964 for maximum unstable case, and 0 when Ri > Ricrit=0.19. 
       ipbl     (i)= 0.0_r8! If 1, PBL is CL, while if 0, PBL is STL.
       kpblh    (i)= 0.0_r8! Layer index containing PBL top within or at the base interface
       wstarPBL (i)= 0.0_r8! Convective velocity within PBL [ m/s ]
       rrho     (i)= 0.0_r8! Density at the lowest layer
       tautotx  (i)= 0.0_r8! Total stress including tms
       tautoty  (i)= 0.0_r8! Total stress including tms
       errorPBL (i)= 0.0_r8! Error function showing whether PBL produced convergent solution or not. [ unit ? ]
       pblhp    (i)= 0.0_r8! PBL top pressure [ Pa ]
       minpblh  (i)= 0.0_r8! Minimum PBL height based on surface stress
       jnk1d    (i)= 0.0_r8
       zero     (i)= 0.0_r8
       tkes     (i)= 0.0_r8! TKE at surface interface [ m2/s2 ]
       ncvfin_o (i)= 0.0_r8! Original number of CLs from 'exacol'
       ncvfin_mg(i)= 0.0_r8! 'ncvfin' after extending-merging from 'zisocl'
       ncvfin_f(i)= 0.0_r8! Final 'ncvfin' after extending-merging & including SRCL

    END DO

    DO m=1,ncvmax
       DO i=1,pcols
          kbase_o(i,m) = 0.0_r8    ! Original external base interface index of CL from 'exacol'
          ktop_o(i,m) = 0.0_r8     ! Original external top  interface index of CL from 'exacol'
          kbase_mg(i,m) = 0.0_r8   ! 'kbase' after extending-merging from 'zisocl'
          ktop_mg(i,m) = 0.0_r8    ! 'ktop' after extending-merging from 'zisocl'
          kbase_f(i,m) = 0.0_r8    ! Final 'kbase' after extending-merging & including SRCL
          ktop_f(i,m) = 0.0_r8     ! Final 'ktop' after extending-merging & including SRCL
          wet(i,m) = 0.0_r8        ! Entrainment rate at the CL top  [ m/s ] 
          web(i,m) = 0.0_r8        ! Entrainment rate at the CL base [ m/s ]. Set to zero if CL is based at surface.
          jtbu(i,m) = 0.0_r8       ! Buoyancy jump across the CL top  [ m/s2 ]  
          jbbu(i,m) = 0.0_r8       ! Buoyancy jump across the CL base [ m/s2 ]  
          evhc(i,m) = 0.0_r8       ! Evaporative enhancement factor at the CL top
          jt2slv(i,m) = 0.0_r8     ! Jump of slv ( across two layers ) at CL top used only for evhc [ J/kg ]
          n2ht(i,m) = 0.0_r8       ! n2 defined at the CL top  interface but using sfuh(kt)   instead of sfi(kt) [ s-2 ] 
          n2hb(i,m) = 0.0_r8       ! n2 defined at the CL base interface but using sflh(kb-1) instead of sfi(kb) [ s-2 ]
          lwp(i,m) = 0.0_r8        ! LWP in the CL top layer [ kg/m2 ]
          opt_depth(i,m) = 0.0_r8  ! Optical depth of the CL top layer
          radinvfrac(i,m) = 0.0_r8 ! Fraction of radiative cooling confined in the top portion of CL top layer
          radf(i,m) = 0.0_r8       ! Buoyancy production at the CL top due to LW radiative cooling [ m2/s3 ]
          wstar(i,m) = 0.0_r8      ! Convective velocity in each CL [ m/s ]
          wstar3fact(i,m) = 0.0_r8 ! Enhancement of 'wstar3' due to entrainment (inverse) [ no unit ]
          ebrk(i,m) = 0.0_r8       ! Net mean TKE of CL including entrainment effect [ m2/s2 ]
          wbrk(i,m) = 0.0_r8       ! Net mean normalized TKE (W) of CL, 'ebrk/b1' including entrainment effect [ m2/s2 ]
          lbrk(i,m) = 0.0_r8       ! Energetic internal thickness of CL [m]
          ricl(i,m) = 0.0_r8       ! CL internal mean Richardson number
          ghcl(i,m) = 0.0_r8       ! Half of normalized buoyancy production of CL
          shcl(i,m) = 0.0_r8       ! Galperin instability function of heat-moisture of CL
          smcl(i,m) = 0.0_r8       ! Galperin instability function of mementum of CL
       END DO
    END DO
    !##########


    ! Buoyancy coefficients : w'b' = ch * w'sl' + cm * w'qt'


    es(1)=0.0_r8! Saturation vapor pressure
    qs(1)=0.0_r8! Saturation specific humidity
    gam(1)=0.0_r8! (L/cp)*dqs/dT
    ep2=0.0_r8; templ(1)=0.0_r8; temps(1)=0.0_r8

    ! ------------------------------- !
    ! Variables for diagnostic output !
    ! ------------------------------- !


    !---------
    zero(:)     = 0._r8
    zero2d(:,:) = 0._r8

    ! ----------------------- !
    ! Main Computation Begins ! 
    ! ----------------------- !
    DO k=1,pver
       DO i=1,ncol
          ufd(i,k)  = u(i,k)
          vfd(i,k)  = v(i,k)
          tfd(i,k)  = t(i,k)
          qvfd(i,k) = qv(i,k)
          qlfd(i,k) = ql(i,k)
       END DO
    END DO
    DO iturb = 1, nturb

       ! Compute total stress by including 'tms'.
       ! Here, in computing 'tms', we can use either iteratively changed 'ufd,vfd' or the
       ! initially given 'u,v' to the PBL scheme. Note that normal stress, 'taux, tauy'
       ! are not changed by iteration. In order to treat 'tms' in a fully implicit way,
       ! I am using updated wind, here.
       DO i=1,ncol
          tautotx(i) = taux(i) - ksrftms(i) * ufd(i,pver)
          tautoty(i) = tauy(i) - ksrftms(i) * vfd(i,pver)
       END DO
       ! Calculate (qt,sl,n2,s2,ri) from a given set of (t,qv,ql,qi,u,v)

       CALL trbintd( &
            pcols    , pver    , ncol  , z       , ufd     , vfd     , tfd   , pmid    , &
            tautotx  , tautoty , ustar , rrho    , s2      , n2      , ri    , zi      , &
            pi       , cldn    , qtfd  , qvfd    , qlfd    , qi      , sfi   , sfuh    , &
            sflh     , slfd    , slv   , slslope , qtslope , chs     , chu   , cms     , &
            cmu      , minpblh ,TSK )

       ! Save initial (i.e., before iterative diffusion) profile of (qt,sl) at each iteration.         
       ! Only necessary for (qt,sl) not (u,v) because (qt,sl) are newly calculated variables. 

       IF( iturb .EQ. 1 ) THEN
          DO k=1,pver
             DO i=1,ncol
                qt(i,k) = qtfd(i,k)
                sl(i,k) = slfd(i,k)
             END DO
          END DO
       ENDIF

       ! Get free atmosphere exchange coefficients. This 'kvf' is not used in UW moist PBL scheme

       CALL austausch_atm( pcols, pver, ncol, ri, s2, kvf )

       ! Initialize kvh/kvm to send to caleddy, depending on model timestep and iteration number
       ! This is necessary for 'wstar-based' entrainment closure.

       IF( iturb .EQ. 1 ) THEN
          IF( kvinit ) THEN
             ! First iteration of first model timestep : Use free tropospheric value or zero.
             IF( use_kvf ) THEN
                DO k=1,pver+1
                   DO i=1,ncol
                      kvh(i,k) = kvf(i,k)
                      kvm(i,k) = kvf(i,k)
                   END DO
                END DO  
             ELSE
                DO k=1,pver+1
                   DO i=1,ncol
                      kvh(i,k) = 0._r8
                      kvm(i,k) = 0._r8
                   END DO
                END DO  
                !PAULO KUBOTA NEW Km Kh Kq
                DO k=1,pver+1
                   DO i=1,ncol
                      !kvm(i,k) = ycoefm  (i,k) 
                      !kvh(i,k) = ycoefh  (i,k) 
                   END DO
                END DO
             ENDIF
          ELSE
             ! First iteration on any model timestep except the first : Use value from previous timestep
             DO k=1,pver+1
                DO i=1,ncol
                   kvh(i,k) = kvh_in(i,k)
                   kvm(i,k) = kvm_in(i,k)
                END DO
             END DO  
          ENDIF
       ELSE
          ! Not the first iteration : Use from previous iteration
          DO k=1,pver+1
             DO i=1,ncol
                kvh(i,k) = kvh_out(i,k)
                kvm(i,k) = kvm_out(i,k)
             END DO
          END DO  
       ENDIF

       ! Calculate eddy diffusivity (kvh_out,kvm_out) and (tke,bprod,sprod) using
       ! a given (kvh,kvm) which are used only for initializing (bprod,sprod)  at
       ! the first part of caleddy. (bprod,sprod) are fully updated at the end of
       ! caleddy after calculating (kvh_out,kvm_out) 

       CALL caleddy( pcols     , pver      , ncol      ,                     &
            slfd      , qtfd      , qlfd      , slv      ,ufd     , &
            vfd       , pi        , z         , zi       ,          &
            qflx      , shflx     , slslope   , qtslope  ,          &
            chu       , chs       , cmu       , cms      ,sfuh    , &
            sflh      , n2        , s2        , ri       ,rrho    , &
            pblh      , ustar     ,                                 &
            kvh       , kvm       , kvh_out   , kvm_out  ,          &
            tpert     , qpert     , qrl       , kvf      , tke    , &
            wstarent  , bprod     , sprod     , minpblh  , wpert  , &
            tkes      , turbtype  , sm_aw     ,                     & 
            kbase_o   , ktop_o    , ncvfin_o  ,                     &
            kbase_mg  , ktop_mg   , ncvfin_mg ,                     &                  
            kbase_f   , ktop_f    , ncvfin_f  ,                     &                  
            wet       , web       , jtbu      , jbbu     ,          &
            evhc      , jt2slv    , n2ht      , n2hb     ,          & 
            lwp       , opt_depth , radinvfrac, radf     ,          &
            wstar     , wstar3fact,                                 &
            ebrk      , wbrk      , lbrk      , ricl     , ghcl   , & 
            shcl      , smcl      , ghi       , shi      , smi    , &
            rii       , lengi     , wcap      , pblhp    , cldn   , &
            landfrac  , ipbl      , kpblh     , wsedl    ,ycoefm  , &
            ycoefh )

  tke2=tke
  swh=0.0
!  CALL moninq(pcols                   ,&!integer      , INTENT(IN   ) :: ix
!       pcols                          ,&!integer      , INTENT(IN   ) :: im
!       pver                           ,&!integer      , INTENT(IN   ) :: km
!       u      (1:pcols,pver  :1:-1)   ,&!real(kind=r8), INTENT(IN   ) :: uo(ix,km)
 !      v      (1:pcols,pver  :1:-1)   ,&!real(kind=r8), INTENT(IN   ) :: vo(ix,km)
 !      t      (1:pcols,pver  :1:-1)   ,&!real(kind=r8), INTENT(IN   ) :: t1(ix,km)
 !      qv     (1:pcols,pver  :1:-1)   ,&!real(kind=r8), INTENT(IN   ) :: q1(ix,km)
 !      ql     (1:pcols,pver  :1:-1)   ,&!real(kind=r8), INTENT(IN   ) :: qliq(ix,km)
 !      bps    (1:pcols,pver  :1:-1)   ,&
 !      swh    (1:pcols,pver  :1:-1)   ,&!real(kind=r8), INTENT(IN   ) :: swh(ix,km)! swh  - real, total sky sw heating rates ( k/s )   ix,levs !
 !      qrl    (1:pcols,pver  :1:-1)   ,&!real(kind=r8), INTENT(IN   ) :: hlw(ix,km)! hlw  - real, total sky lw heating rates ( k/s )   ix,levs !
 !      taux   (1:pcols)               ,&!real(r8), intent(in)    :: taux(pcols)               ! Zonal wind stress at surface [ N/m2 ]
 !      tauy   (1:pcols)               ,&!real(r8), intent(in)    :: tauy(pcols)               ! Meridional wind stress at surface [ N/m2 ]
 !      TSK    (1:pcols)               ,&!real(kind=r8), INTENT(IN   ) :: TSKIN(im)
 !      shflx  (1:pcols)               ,&!real(kind=r8), INTENT(IN   ) :: heat(im)
 !      qflx   (1:pcols)               ,&!real(kind=r8), INTENT(IN   ) :: evap(im)
 !      kpbl   (1:pcols)               ,&!integer      , INTENT(OUt  ) :: kpbl(im)
 !      pi     (1:pcols,pver+1:1:-1)   ,&!real(kind=r8), INTENT(IN   ) :: prsi(ix,km+1)
 !      pmid   (1:pcols,pver  :1:-1)   ,&!real(kind=r8), INTENT(IN   ) :: prsl(ix,km)
 !      zi     (1:pcols,pver+1:1:-1)   ,&!real(kind=r8), INTENT(IN   ) :: phii(ix,km+1)
 !      z      (1:pcols,pver  :1:-1)   ,&!real(kind=r8), INTENT(IN   ) :: phil(ix,km)
 !      ztodt                          ,&!real(kind=r8), INTENT(IN   ) :: deltim
 !      hpbl   (1:pcols)               ,&!real(kind=r8), INTENT(OUt  ) :: hpbl(im)
 !      z0     (1:pcols)               ,&!real(kind=r8), INTENT(IN   ) :: Z0RL(im)
 !      mlsi   (1:pcols)               ,&!real(kind=r8), INTENT(IN   ) :: mlsi(im)!       !     slmsk      - real, sea/land/ice mask (=0/1/2)                im   !
 !      vcover (1:pcols)               ,&!real(kind=r8), INTENT(IN   ) :: vfrac(IM)!      vfrac       - real, vegetation fraction                       im   !
 !      ycoefh (1:pcols,pver+1:1:-1)   ,&!real(kind=r8), INTENT(OUt  ) :: dkt (im,km+1)
 !/      ycoefm (1:pcols,pver+1:1:-1)     )!real(kind=r8), INTENT(OUt  ) :: dku (im,km+1)
!                DO k=1,pver+1
!                   DO i=1,ncol
!                      PRINT*, k,  ycoefm (i,k),  ycoefh (i,k), hpbl   (i), kvm_out(i,k), kvh_out(i,k),pblh(i)
!                   END DO
!                END DO

       ! Calculate errorPBL to check whether PBL produced convergent solutions or not.

       IF( iturb .EQ. nturb ) THEN
          DO i = 1, ncol
             errorPBL(i) = 0._r8 
             DO k = 1, pver
                errorPBL(i) = errorPBL(i) + ( kvh(i,k) - kvh_out(i,k) )**2 
             END DO
             errorPBL(i) = SQRT(errorPBL(i)/pver)
          END DO
       END IF

       ! Eddy diffusivities which will be used for the initialization of (bprod,
       ! sprod) in 'caleddy' at the next iteration step.

       IF( iturb .GT. 1 .AND. iturb .LT. nturb ) THEN
          DO k=1,pver+1
             DO i=1,ncol
                kvm_out(i,k) = lambda * kvm_out(i,k) + ( 1._r8 - lambda ) * kvm(i,k)
                kvh_out(i,k) = lambda * kvh_out(i,k) + ( 1._r8 - lambda ) * kvh(i,k)

               ! kvm_out(i,k) = lambdaEns * kvm(i,k) + ( 2._r8*lambdaEns ) * kvm_out(i,k) +  lambdaEns *ycoefm(i,k)
               ! kvh_out(i,k) = lambdaEns * kvh(i,k) + ( 2._r8*lambdaEns ) * kvh_out(i,k) +  lambdaEns *ycoefh(i,k)

             END DO
          END DO  
       ENDIF

       ! Set nonlocal terms to zero for flux diagnostics, since not used by caleddy.
       DO k=1,pver+1
          DO i=1,ncol
             cgh(i,k) = 0._r8
             cgs(i,k) = 0._r8      
          END DO
       END DO  

       IF( iturb .LT. nturb ) THEN

          ! Each time we diffuse the original state
          DO k=1,pver
             DO i=1,ncol
                slfd(i,k)  = sl(i,k)
                qtfd(i,k)  = qt(i,k)
                ufd(i,k)   = u(i,k)
                vfd(i,k)   = v(i,k)
             END DO
          END DO  

          ! Diffuse initial profile of each time step using a given (kvh_out,kvm_out)
          ! In the below 'compute_vdiff', (slfd,qtfd,ufd,vfd) are 'inout' variables.

          CALL compute_vdiff( lchnk   ,                                   &
               pcols   , pver     , 1        , ncol         , TSK       , & 
               QSFC    , pmid     , &
               pi      , rpdel    , t        , ztodt        , taux      , &
               tauy    , shflx    , qflx     , ntop_turb    , nbot_turb , &
               kvh_out , kvm_out  , kvh_out  , cgs          , cgh       , &
               zi      , ksrftms  , zero     , fieldlist_wet,             &
               ufd     , vfd      , qtfd     , slfd         ,             &
               jnk1d   , jnk1d    , jnk2d    , jnk1d        , errstring , &
               !                               tauresx , tauresy  , 0        , .false. )
               tauresx , tauresy  , 0          )

          ! Retrieve (tfd,qvfd,qlfd) from (slfd,qtfd) in order to 
          ! use 'trbintd' at the next iteration.

          DO k = 1, pver
             DO i = 1, ncol
                ! ----------------------------------------------------- ! 
                ! Compute the condensate 'qlfd' in the updated profiles !
                ! ----------------------------------------------------- !  
                ! Option.1 : Assume grid-mean condensate is homogeneously diffused by the moist turbulence scheme.
                !            This should bs used if 'pseudodiff = .false.' in vertical_diffusion.F90.
                ! Modification : Need to be check whether below is correct in the presence of ice, qi.       
                !                I should understand why the variation of ice, qi is neglected during diffusion.
                templ(1)     = ( slfd(i,k) - g*z(i,k) ) / cpair
                status    =   fqsatd( templ(1), pmid(i,k), es(1), qs(1), gam(1), 1 )
                ep2       =  0.622_r8 
                temps(1)     =   templ(1) + ( qtfd(i,k) - qs(1) ) / ( cpair / latvap + latvap * qs(1) / ( rair * templ(1)**2 ) )
                status    =   fqsatd( temps(1), pmid(i,k), es(1), qs(1), gam(1), 1 )
                qlfd(i,k) =   MAX( qtfd(i,k) - qi(i,k) - qs(1) ,0._r8 )
                ! Option.2 : Assume condensate is not diffused by the moist turbulence scheme. 
                !            This should bs used if 'pseudodiff = .true.'  in vertical_diffusion.F90.       
                ! qlfd(i,k) = ql(i,k)
                ! ----------------------------- !
                ! Compute the other 'qvfd, tfd' ! 
                ! ----------------------------- !
                qvfd(i,k) = MAX( 0._r8, qtfd(i,k) - qi(i,k) - qlfd(i,k) )
                tfd(i,k)  = ( slfd(i,k) + latvap * qlfd(i,k) + latsub * qi(i,k) - g*z(i,k)) / cpair
             END DO
          END DO
       ENDIF

       ! Debug 
       ! icol = phys_debug_col(lchnk) 
       ! if( icol > 0 .and. get_nstep() .ge. 1 ) then
       !     write(iulog,*) ' '
       !     write(iulog,*) 'eddy_diff debug at the end of iteration' 
       !     write(iulog,*) 't,     qv,     ql,     cld,     u,     v'
       !     do k = pver-3, pver
       !        write (iulog,*) k, tfd(icol,k), qvfd(icol,k), qlfd(icol,k), cldn(icol,k), ufd(icol,k), vfd(icol,k)
       !     end do
       ! endif
       ! Debug

    END DO  ! End of 'iturb' iteration

    DO k=1,pver+1
       DO i = 1, ncol
          kvq(i,k) = kvh_out(i,k)
       END DO
    END DO
    
    ! Compute 'wstar' within the PBL for use in the future convection scheme.

    DO i = 1, ncol
       IF( ipbl(i) .EQ. 1._r8 ) THEN 
          wstarPBL(i) = MAX( 0._r8, wstar(i,1) )
       ELSE
          wstarPBL(i) = 0._r8
       ENDIF
    END DO

    ! --------------------------------------------------------------- !
    ! Writing for detailed diagnostic analysis of UW moist PBL scheme !
    ! --------------------------------------------------------------- !

    !call outfld( 'UW_errorPBL',    errorPBL,   pcols,   lchnk )

    !call outfld( 'UW_n2',          n2,         pcols,   lchnk )
    !call outfld( 'UW_s2',          s2,         pcols,   lchnk )
    !call outfld( 'UW_ri',          ri,         pcols,   lchnk )

    !call outfld( 'UW_sfuh',        sfuh,       pcols,   lchnk )
    !call outfld( 'UW_sflh',        sflh,       pcols,   lchnk )
    !call outfld( 'UW_sfi',         sfi,        pcols,   lchnk )

    !call outfld( 'UW_cldn',        cldn,       pcols,   lchnk )
    !call outfld( 'UW_qrl',         qrl,        pcols,   lchnk )
    !call outfld( 'UW_ql',          qlfd,       pcols,   lchnk )

    !call outfld( 'UW_chu',         chu,        pcols,   lchnk )
    !call outfld( 'UW_chs',         chs,        pcols,   lchnk )
    !call outfld( 'UW_cmu',         cmu,        pcols,   lchnk )
    !call outfld( 'UW_cms',         cms,        pcols,   lchnk )

    !call outfld( 'UW_tke',         tke,        pcols,   lchnk )
    !call outfld( 'UW_wcap',        wcap,       pcols,   lchnk )
    !call outfld( 'UW_bprod',       bprod,      pcols,   lchnk )
    !call outfld( 'UW_sprod',       sprod,      pcols,   lchnk )

    !call outfld( 'UW_kvh',         kvh_out,    pcols,   lchnk )
    !call outfld( 'UW_kvm',         kvm_out,    pcols,   lchnk )

    !call outfld( 'UW_pblh',        pblh,       pcols,   lchnk )
    !call outfld( 'UW_pblhp',       pblhp,      pcols,   lchnk )
    !call outfld( 'UW_tpert',       tpert,      pcols,   lchnk )
    !call outfld( 'UW_qpert',       qpert,      pcols,   lchnk )
    !call outfld( 'UW_wpert',       wpert,      pcols,   lchnk )

    !call outfld( 'UW_ustar',       ustar,      pcols,   lchnk )
    !call outfld( 'UW_tkes',        tkes,       pcols,   lchnk )
    !call outfld( 'UW_minpblh',     minpblh,    pcols,   lchnk )

    !call outfld( 'UW_turbtype',    turbtype,   pcols,   lchnk )

    !call outfld( 'UW_kbase_o',     kbase_o,    pcols,   lchnk )
    !call outfld( 'UW_ktop_o',      ktop_o,     pcols,   lchnk )
    !call outfld( 'UW_ncvfin_o',    ncvfin_o,   pcols,   lchnk )

    !call outfld( 'UW_kbase_mg',    kbase_mg,   pcols,   lchnk )
    !call outfld( 'UW_ktop_mg',     ktop_mg,    pcols,   lchnk )
    !call outfld( 'UW_ncvfin_mg',   ncvfin_mg,  pcols,   lchnk )

    !call outfld( 'UW_kbase_f',     kbase_f,    pcols,   lchnk )
    !call outfld( 'UW_ktop_f',      ktop_f,     pcols,   lchnk )
    !call outfld( 'UW_ncvfin_f',    ncvfin_f,   pcols,   lchnk ) 

    !call outfld( 'UW_wet',         wet,        pcols,   lchnk )
    !call outfld( 'UW_web',         web,        pcols,   lchnk )
    !call outfld( 'UW_jtbu',        jtbu,       pcols,   lchnk )
    !call outfld( 'UW_jbbu',        jbbu,       pcols,   lchnk )
    !call outfld( 'UW_evhc',        evhc,       pcols,   lchnk )
    !call outfld( 'UW_jt2slv',      jt2slv,     pcols,   lchnk )
    !call outfld( 'UW_n2ht',        n2ht,       pcols,   lchnk )
    !call outfld( 'UW_n2hb',        n2hb,       pcols,   lchnk )
    !call outfld( 'UW_lwp',         lwp,        pcols,   lchnk )
    !call outfld( 'UW_optdepth',    opt_depth,  pcols,   lchnk )
    !call outfld( 'UW_radfrac',     radinvfrac, pcols,   lchnk )
    !call outfld( 'UW_radf',        radf,       pcols,   lchnk )
    !call outfld( 'UW_wstar',       wstar,      pcols,   lchnk )
    !call outfld( 'UW_wstar3fact',  wstar3fact, pcols,   lchnk )
    !call outfld( 'UW_ebrk',        ebrk,       pcols,   lchnk )
    !call outfld( 'UW_wbrk',        wbrk,       pcols,   lchnk )
    !call outfld( 'UW_lbrk',        lbrk,       pcols,   lchnk )
    !call outfld( 'UW_ricl',        ricl,       pcols,   lchnk )
    !call outfld( 'UW_ghcl',        ghcl,       pcols,   lchnk )
    !call outfld( 'UW_shcl',        shcl,       pcols,   lchnk )
    !call outfld( 'UW_smcl',        smcl,       pcols,   lchnk )

    !call outfld( 'UW_gh',          ghi,        pcols,   lchnk )
    !call outfld( 'UW_sh',          shi,        pcols,   lchnk )
    !call outfld( 'UW_sm',          smi,        pcols,   lchnk )
    !call outfld( 'UW_ria',         rii,        pcols,   lchnk )
    !call outfld( 'UW_leng',        lengi,      pcols,   lchnk )

    RETURN

  END SUBROUTINE compute_eddy_diff


  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  SUBROUTINE trbintd( pcols   , pver    , ncol    ,                               &
       z       , u       , v       ,                               &
       t       , pmid    , taux    ,                               &
       tauy    , ustar   , rrho    ,                               &
       s2      , n2      , ri      ,                               &
       zi      , pi      , cld     ,                               &
       qt      , qv      , ql      , qi      , sfi     , sfuh    , &
       sflh    , sl      , slv     , slslope , qtslope ,           &
       chs     , chu     , cms     , cmu     , minpblh ,TSK )
    !----------------------------------------------------------------------- !
    ! Purpose: Calculate buoyancy coefficients at all interfaces including   !
    !          surface. Also, computes the profiles of ( sl,qt,n2,s2,ri ).   !
    !          Note that (n2,s2,ri) are defined at each interfaces except    !
    !          surface.                                                      !
    !                                                                        !
    ! Author: B. Stevens  ( Extracted from pbldiff, August, 2000 )           !
    !         Sungsu Park ( August 2006, May. 2008 )                         !
    !----------------------------------------------------------------------- !

    IMPLICIT NONE

    ! --------------- !
    ! Input arguments !
    ! --------------- !

    INTEGER,  INTENT(in)  :: pcols                            ! Number of atmospheric columns   
    INTEGER,  INTENT(in)  :: pver                             ! Number of atmospheric layers   
    INTEGER,  INTENT(in)  :: ncol                             ! Number of atmospheric columns
    REAL(r8), INTENT(in)  :: z(pcols,pver)                    ! Layer mid-point height above surface [ m ]
    REAL(r8), INTENT(in)  :: u(pcols,pver)                    ! Layer mid-point u [ m/s ]
    REAL(r8), INTENT(in)  :: v(pcols,pver)                    ! Layer mid-point v [ m/s ]
    REAL(r8), INTENT(in)  :: t(pcols,pver)                    ! Layer mid-point temperature [ K ]
    REAL(r8), INTENT(in)  :: pmid(pcols,pver)                 ! Layer mid-point pressure [ Pa ]
    REAL(r8), INTENT(in)  :: taux(pcols)                      ! Surface u stress [ N/m2 ]
    REAL(r8), INTENT(in)  :: tauy(pcols)                      ! Surface v stress [ N/m2 ]
    REAL(r8), INTENT(in)  :: zi(pcols,pver+1)                 ! Interface height [ m ]
    REAL(r8), INTENT(in)  :: pi(pcols,pver+1)                 ! Interface pressure [ Pa ]
    REAL(r8), INTENT(in)  :: cld(pcols,pver)                  ! Stratus fraction
    REAL(r8), INTENT(in)  :: qv(pcols,pver)                   ! Water vapor specific humidity [ kg/kg ]
    REAL(r8), INTENT(in)  :: ql(pcols,pver)                   ! Liquid water specific humidity [ kg/kg ]
    REAL(r8), INTENT(in)  :: qi(pcols,pver)                   ! Ice water specific humidity [ kg/kg ]
    REAL(r8), INTENT(in)  :: TSK(pcols)
    !INTEGER,  EXTERNAL    :: qsat

    ! ---------------- !
    ! Output arguments !
    ! ---------------- !

    REAL(r8), INTENT(out) :: ustar(pcols)                     ! Surface friction velocity [ m/s ]
    REAL(r8), INTENT(out) :: s2(pcols,pver)                   ! Interfacial ( except surface ) shear squared [ s-2 ]
    REAL(r8), INTENT(out) :: n2(pcols,pver)                   ! Interfacial ( except surface ) buoyancy frequency [ s-2 ]
    REAL(r8), INTENT(out) :: ri(pcols,pver)                   ! Interfacial ( except surface ) Richardson number, 'n2/s2'

    REAL(r8), INTENT(out) :: qt(pcols,pver)                   ! Total specific humidity [ kg/kg ]
    REAL(r8), INTENT(out) :: sfi(pcols,pver+1)                ! Interfacial layer saturation fraction [ fraction ]
    REAL(r8), INTENT(out) :: sfuh(pcols,pver)                 ! Saturation fraction in upper half-layer [ fraction ]
    REAL(r8), INTENT(out) :: sflh(pcols,pver)                 ! Saturation fraction in lower half-layer [ fraction ]
    REAL(r8), INTENT(out) :: sl(pcols,pver)                   ! Liquid water static energy [ J/kg ] 
    REAL(r8), INTENT(out) :: slv(pcols,pver)                  ! Liquid water virtual static energy [ J/kg ]

    REAL(r8), INTENT(out) :: chu(pcols,pver+1)                ! Heat buoyancy coef for dry states at all interfaces, finally. [ unit ? ]
    REAL(r8), INTENT(out) :: chs(pcols,pver+1)                ! heat buoyancy coef for sat states at all interfaces, finally. [ unit ? ]
    REAL(r8), INTENT(out) :: cmu(pcols,pver+1)                ! Moisture buoyancy coef for dry states at all interfaces, finally. [ unit ? ]
    REAL(r8), INTENT(out) :: cms(pcols,pver+1)                ! Moisture buoyancy coef for sat states at all interfaces, finally. [ unit ? ]
    REAL(r8), INTENT(out) :: slslope(pcols,pver)              ! Slope of 'sl' in each layer
    REAL(r8), INTENT(out) :: qtslope(pcols,pver)              ! Slope of 'qt' in each layer
    REAL(r8), INTENT(out) :: rrho(pcols)                      ! 1./bottom level density [ m3/kg ]
    REAL(r8), INTENT(out) :: minpblh(pcols)                   ! Minimum PBL height based on surface stress [ m ]

    ! --------------- !
    ! Local Variables !
    ! --------------- ! 

    INTEGER               :: i                                ! Longitude index
    INTEGER               :: k, km1                           ! Level index
    INTEGER               :: status                           ! Status returned by function calls

    REAL(r8)              :: qs(pcols,pver)                   ! Saturation specific humidity
    REAL(r8)              :: es(pcols,pver)                   ! Saturation vapor pressure
    REAL(r8)              :: gam(pcols,pver)                  ! (l/cp)*(d(qs)/dT)
    REAL(r8)              :: rdz                              ! 1 / (delta z) between midpoints
    REAL(r8)              :: dsldz                            ! 'delta sl / delta z' at interface
    REAL(r8)              :: dqtdz                            ! 'delta qt / delta z' at interface
    REAL(r8)              :: ch                               ! 'sfi' weighted ch at the interface
    REAL(r8)              :: cm                               ! 'sfi' weighted cm at the interface
    REAL(r8)              :: bfact                            ! Buoyancy factor in n2 calculations
    REAL(r8)              :: product                          ! Intermediate vars used to find slopes
    REAL(r8)              :: dsldp_a, dqtdp_a                 ! Slopes across interface above 
    REAL(r8)              :: dsldp_b(pcols), dqtdp_b(pcols)   ! Slopes across interface below
    DO k=1,pver
       DO i=1,pcols
          s2     (i,k) = 0.0_r8
          n2     (i,k) = 0.0_r8
          ri     (i,k) = 0.0_r8
          slslope(i,k) = 0.0_r8
          qtslope(i,k) = 0.0_r8
          qt     (i,k) = 0.0_r8
          sfuh   (i,k) = 0.0_r8
          sflh   (i,k) = 0.0_r8
          sl     (i,k) = 0.0_r8
          slv    (i,k) = 0.0_r8
          qs     (i,k) = 0.0_r8! Saturation specific humidity
          es     (i,k) = 0.0_r8! Saturation vapor pressure
          gam    (i,k) = 0.0_r8! (l/cp)*(d(qs)/dT)

       END DO
    END DO
    DO k=1,pver+1
       DO i=1,pcols
          sfi    (i,k) = 0.0_r8
          chu    (i,k) = 0.0_r8
          chs    (i,k) = 0.0_r8
          cmu    (i,k) = 0.0_r8
          cms    (i,k) = 0.0_r8
       END DO
    END DO
    DO i=1,pcols    
       dsldp_b(i)= 0.0_r8
       dqtdp_b(i)= 0.0_r8! Slopes across interface below

       ustar  (i)= 0.0_r8
       rrho   (i)= 0.0_r8
       minpblh(i)= 0.0_r8
    END DO
    rdz= 0.0_r8;dsldz= 0.0_r8;dqtdz= 0.0_r8;ch= 0.0_r8
    cm= 0.0_r8;bfact= 0.0_r8;product= 0.0_r8
    dsldp_a= 0.0_r8; dqtdp_a= 0.0_r8
    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !

    ! Compute ustar, and kinematic surface fluxes from surface energy fluxes

    DO i = 1, ncol
       rrho(i)    = rair * TSK(i) / pmid(i,pver)
       ustar(i)   = MAX( SQRT( SQRT( taux(i)**2 + tauy(i)**2 ) * rrho(i) ), ustar_min )
       minpblh(i) = 100.0_r8 * ustar(i)                       ! By construction, 'minpblh' is larger than 1 [m] when 'ustar_min = 0.01'. 
    END DO

    ! Calculate conservative scalars (qt,sl,slv) and buoyancy coefficients at the layer mid-points.
    ! Note that 'ntop_turb = 1', 'nbot_turb = pver'

    DO k = ntop_turb, nbot_turb
       status = fqsatd( t(1,k), pmid(1,k), es(1,k), qs(1,k), gam(1,k), ncol )
       DO i = 1, ncol
          qt(i,k)  = qv(i,k) + ql(i,k) + qi(i,k) 
          sl(i,k)  = cpair * t(i,k) + g * z(i,k) - latvap * ql(i,k) - latsub * qi(i,k)
          slv(i,k) = sl(i,k) * ( 1._r8 + zvir * qt(i,k) )
          ! Thermodynamic coefficients for buoyancy flux - in this loop these are
          ! calculated at mid-points; later,  they will be averaged to interfaces,
          ! where they will ultimately be used.  At the surface, the coefficients
          ! are taken from the lowest mid point.
          bfact    = g / ( t(i,k) * ( 1._r8 + zvir * qv(i,k) - ql(i,k) - qi(i,k) ) )
          chu(i,k) = ( 1._r8 + zvir * qt(i,k) ) * bfact / cpair
          chs(i,k) = ( ( 1._r8 + ( 1._r8 + zvir ) * gam(i,k) * cpair * t(i,k) / latvap ) / ( 1._r8 + gam(i,k) ) ) * bfact / cpair
          cmu(i,k) = zvir * bfact * t(i,k)
          cms(i,k) = latvap * chs(i,k)  -  bfact * t(i,k)
       END DO
    END DO

    DO i = 1, ncol
       chu(i,pver+1) = chu(i,pver)
       chs(i,pver+1) = chs(i,pver)
       cmu(i,pver+1) = cmu(i,pver)
       cms(i,pver+1) = cms(i,pver)
    END DO

    ! Compute slopes of conserved variables sl, qt within each layer k. 
    ! 'a' indicates the 'above' gradient from layer k-1 to layer k and 
    ! 'b' indicates the 'below' gradient from layer k   to layer k+1.
    ! We take a smaller (in absolute value)  of these gradients as the
    ! slope within layer k. If they have opposite signs,   gradient in 
    ! layer k is taken to be zero. I should re-consider whether   this
    ! profile reconstruction is the best or not.
    ! This is similar to the profile reconstruction used in the UWShCu. 

    DO i = 1, ncol
       ! Slopes at endpoints determined by extrapolation
       slslope(i,pver) = ( sl(i,pver) - sl(i,pver-1) ) / ( pmid(i,pver) - pmid(i,pver-1) )
       qtslope(i,pver) = ( qt(i,pver) - qt(i,pver-1) ) / ( pmid(i,pver) - pmid(i,pver-1) )
       slslope(i,1)    = ( sl(i,2) - sl(i,1) ) / ( pmid(i,2) - pmid(i,1) )
       qtslope(i,1)    = ( qt(i,2) - qt(i,1) ) / ( pmid(i,2) - pmid(i,1) )
       dsldp_b(i)      = slslope(i,1)
       dqtdp_b(i)      = qtslope(i,1)
    END DO

    DO k = 2, pver - 1
       DO i = 1, ncol
          dsldp_a    = dsldp_b(i)
          dqtdp_a    = dqtdp_b(i)
          dsldp_b(i) = ( sl(i,k+1) - sl(i,k) ) / ( pmid(i,k+1) - pmid(i,k) )
          dqtdp_b(i) = ( qt(i,k+1) - qt(i,k) ) / ( pmid(i,k+1) - pmid(i,k) )
          product    = dsldp_a * dsldp_b(i)
          IF( product .LE. 0._r8 ) THEN 
             slslope(i,k) = 0._r8
          ELSE IF( product .GT. 0._r8 .AND. dsldp_a .LT. 0._r8 ) THEN 
             slslope(i,k) = MAX( dsldp_a, dsldp_b(i) )
          ELSE IF( product .GT. 0._r8 .AND. dsldp_a .GT. 0._r8 ) THEN 
             slslope(i,k) = MIN( dsldp_a, dsldp_b(i) )
          END IF
          product = dqtdp_a*dqtdp_b(i)
          IF( product .LE. 0._r8 ) THEN 
             qtslope(i,k) = 0._r8
          ELSE IF( product .GT. 0._r8 .AND. dqtdp_a .LT. 0._r8 ) THEN 
             qtslope(i,k) = MAX( dqtdp_a, dqtdp_b(i) )
          ELSE IF( product .GT. 0._r8 .AND. dqtdp_a .GT. 0._r8 ) THEN 
             qtslope(i,k) = MIN( dqtdp_a, dqtdp_b(i) )
          END IF
       END DO ! i
    END DO ! k

    !  Compute saturation fraction at the interfacial layers for use in buoyancy
    !  flux computation.

    CALL sfdiag( pcols  , pver    , ncol    , qt      , ql      , sl      , & 
         pi     , pmid    , zi      , cld     , sfi     , sfuh    , &
         sflh   , slslope , qtslope  )

    ! Calculate buoyancy coefficients at all interfaces (1:pver+1) and (n2,s2,ri) 
    ! at all interfaces except surface. Note 'nbot_turb = pver', 'ntop_turb = 1'.
    ! With the previous definition of buoyancy coefficients at the surface, the 
    ! resulting buoyancy coefficients at the top and surface interfaces becomes 
    ! identical to the buoyancy coefficients at the top and bottom layers. Note 
    ! that even though the dimension of (s2,n2,ri) is 'pver',  they are defined
    ! at interfaces ( not at the layer mid-points ) except the surface. 

    DO k = nbot_turb, ntop_turb + 1, -1
       km1 = k - 1
       DO i = 1, ncol
          rdz      = 1._r8 / ( z(i,km1) - z(i,k) )
          dsldz    = ( sl(i,km1) - sl(i,k) ) * rdz
          dqtdz    = ( qt(i,km1) - qt(i,k) ) * rdz 
          chu(i,k) = ( chu(i,km1) + chu(i,k) ) * 0.5_r8
          chs(i,k) = ( chs(i,km1) + chs(i,k) ) * 0.5_r8
          cmu(i,k) = ( cmu(i,km1) + cmu(i,k) ) * 0.5_r8
          cms(i,k) = ( cms(i,km1) + cms(i,k) ) * 0.5_r8
          ch       = chu(i,k) * ( 1._r8 - sfi(i,k) ) + chs(i,k) * sfi(i,k)
          cm       = cmu(i,k) * ( 1._r8 - sfi(i,k) ) + cms(i,k) * sfi(i,k)
          n2(i,k)  = ch * dsldz +  cm * dqtdz
          s2(i,k)  = ( ( u(i,km1) - u(i,k) )**2 + ( v(i,km1) - v(i,k) )**2) * rdz**2
          s2(i,k)  = MAX( ntzero, s2(i,k) )
          ri(i,k)  = n2(i,k) / s2(i,k)
       END DO
    END DO
    DO i = 1, ncol
       n2(i,1) = n2(i,2)
       s2(i,1) = s2(i,2)
       ri(i,1) = ri(i,2)
    END DO

    RETURN

  END SUBROUTINE trbintd


  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  SUBROUTINE sfdiag( pcols   , pver    , ncol    , qt      , ql      , sl      , &
       pi      , pm      , zi      , cld     , sfi     , sfuh    , &
       sflh    , slslope , qtslope  )
    !----------------------------------------------------------------------- ! 
    !                                                                        !
    ! Purpose: Interface for calculating saturation fractions  at upper and  ! 
    !          lower-half layers, & interfaces for use by turbulence scheme  !
    !                                                                        !
    ! Method : Various but 'l' should be chosen for consistency.             !
    !                                                                        ! 
    ! Author : B. Stevens and C. Bretherton (August 2000)                    !
    !          Sungsu Park. August 2006.                                     !
    !                       May.   2008.                                     ! 
    !                                                                        !  
    ! S.Park : The computed saturation fractions are repeatedly              !
    !          used to compute buoyancy coefficients in'trbintd' & 'caleddy'.!  
    !----------------------------------------------------------------------- !

    IMPLICIT NONE       

    ! --------------- !
    ! Input arguments !
    ! --------------- !

    !INTEGER,  EXTERNAL    :: qsat

    INTEGER,  INTENT(in)  :: pcols               ! Number of atmospheric columns   
    INTEGER,  INTENT(in)  :: pver                ! Number of atmospheric layers   
    INTEGER,  INTENT(in)  :: ncol                ! Number of atmospheric columns   

    REAL(r8), INTENT(in)  :: sl(pcols,pver)      ! Liquid water static energy [ J/kg ]
    REAL(r8), INTENT(in)  :: qt(pcols,pver)      ! Total water specific humidity [ kg/kg ]
    REAL(r8), INTENT(in)  :: ql(pcols,pver)      ! Liquid water specific humidity [ kg/kg ]
    REAL(r8), INTENT(in)  :: pi(pcols,pver+1)    ! Interface pressures [ Pa ]
    REAL(r8), INTENT(in)  :: pm(pcols,pver)      ! Layer mid-point pressures [ Pa ]
    REAL(r8), INTENT(in)  :: zi(pcols,pver+1)    ! Interface heights [ m ]
    REAL(r8), INTENT(in)  :: cld(pcols,pver)     ! Stratiform cloud fraction [ fraction ]
    REAL(r8), INTENT(in)  :: slslope(pcols,pver) ! Slope of 'sl' in each layer
    REAL(r8), INTENT(in)  :: qtslope(pcols,pver) ! Slope of 'qt' in each layer

    ! ---------------- !
    ! Output arguments !
    ! ---------------- !

    REAL(r8), INTENT(out) :: sfi(pcols,pver+1)   ! Interfacial layer saturation fraction [ fraction ]
    REAL(r8), INTENT(out) :: sfuh(pcols,pver)    ! Saturation fraction in upper half-layer [ fraction ]
    REAL(r8), INTENT(out) :: sflh(pcols,pver)    ! Saturation fraction in lower half-layer [ fraction ]

    ! --------------- !
    ! Local Variables !
    ! --------------- !

    INTEGER               :: i                   ! Longitude index
    INTEGER               :: k                   ! Vertical index
    INTEGER               :: km1                 ! k-1
    INTEGER               :: status              ! Status returned by function calls
    REAL(r8)              :: sltop, slbot        ! sl at top/bot of grid layer
    REAL(r8)              :: qttop, qtbot        ! qt at top/bot of grid layer
    REAL(r8)              :: tltop(1), tlbot(1)  ! Liquid water temperature at top/bot of grid layer
    REAL(r8)              :: qxtop, qxbot        ! Sat excess at top/bot of grid layer
    REAL(r8)              :: qxm                 ! Sat excess at midpoint
    REAL(r8)              :: es(1)               ! Saturation vapor pressure
    REAL(r8)              :: qs(1)               ! Saturation spec. humidity
    REAL(r8)              :: gam(1)              ! (L/cp)*dqs/dT
    REAL(r8)              :: cldeff(pcols,pver)  ! Effective Cloud Fraction [ fraction ]

    ! ----------------------- !
    ! Main Computation Begins ! 
    ! ----------------------- !
    DO k=1,pver+1
       DO i=1,pcols
          sfi   (i,k)   = 0.0_r8
       END DO
    END DO
    DO k=1,pver
       DO i=1,pcols
          sflh  (i,k)   = 0.0_r8
          sfuh  (i,k)   = 0.0_r8
          cldeff(i,k)   = 0.0_r8
       END DO
    END DO
    
    SELECT CASE (sftype)
    CASE ('d')
       ! ----------------------------------------------------------------------- !
       ! Simply use the given stratus fraction ('horizontal' cloud partitioning) !
       ! ----------------------------------------------------------------------- !
       DO k = ntop_turb + 1, nbot_turb
          km1 = k - 1
          DO i = 1, ncol
             sfuh(i,k) = cld(i,k)
             sflh(i,k) = cld(i,k)
             sfi(i,k)  = 0.5_r8 * ( sflh(i,km1) + MIN( sflh(i,km1), sfuh(i,k) ) )
          END DO
       END DO
       DO i = 1, ncol
          sfi(i,pver+1) = sflh(i,pver) 
       END DO
    CASE ('l')
       ! ------------------------------------------ !
       ! Use modified stratus fraction partitioning !
       ! ------------------------------------------ !
       DO k = ntop_turb + 1, nbot_turb
          km1 = k - 1
          DO i = 1, ncol
             cldeff(i,k) = cld(i,k)
             sfuh(i,k)   = cld(i,k)
             sflh(i,k)   = cld(i,k)
             IF( ql(i,k) .LT. qmin(2) ) THEN
                sfuh(i,k) = 0._r8
                sflh(i,k) = 0._r8
             END IF
             ! Modification : The contribution of ice should be carefully considered.
             IF( choice_evhc .EQ. 'ramp' .OR. choice_radf .EQ. 'ramp' ) THEN 
                cldeff(i,k) = cld(i,k) * MIN( ql(i,k) / qmin(2), 1._r8 )
                sfuh(i,k)   = cldeff(i,k)
                sflh(i,k)   = cldeff(i,k)
             ELSEIF( choice_evhc .EQ. 'maxi' .OR. choice_radf .EQ. 'maxi' ) THEN 
                cldeff(i,k) = cld(i,k)
                sfuh(i,k)   = cldeff(i,k)
                sflh(i,k)   = cldeff(i,k)
             ENDIF
             ! At the stratus top, take the minimum interfacial saturation fraction
             sfi(i,k) = 0.5_r8 * ( sflh(i,km1) + MIN( sfuh(i,k), sflh(i,km1) ) )
             ! Modification : Currently sfi at the top and surface interfaces are set to be zero.
             !                Also, sfuh and sflh in the top model layer is set to be zero.
             !                However, I may need to set 
             !                         do i = 1, ncol
             !                            sfi(i,pver+1) = sflh(i,pver) 
             !                         end do
             !                for treating surface-based fog. 
             ! OK. I added below block similar to the other cases.
          END DO
       END DO
       DO i = 1, ncol
          sfi(i,pver+1) = sflh(i,pver)
       END DO
    CASE ('u')
       ! ------------------------------------------------------------------------- !
       ! Use unsaturated buoyancy - since sfi, sfuh, sflh have already been zeroed !
       ! nothing more need be done for this case.                                  !
       ! ------------------------------------------------------------------------- !
    CASE ('z')
       ! ------------------------------------------------------------------------- !
       ! Calculate saturation fraction based on whether the air just above or just !
       ! below the interface is saturated, i.e. with vertical cloud partitioning.  !
       ! The saturation fraction of the interfacial layer between mid-points k and !
       ! k+1 is computed by averaging the saturation fraction   of the half-layers !
       ! above and below the interface,  with a special provision   for cloud tops !
       ! (more cloud in the half-layer below than in the half-layer above).In each !
       ! half-layer, vertical partitioning of  cloud based on the slopes diagnosed !
       ! above is used.     Loop down through the layers, computing the saturation !
       ! fraction in each half-layer (sfuh for upper half, sflh for lower half).   !
       ! Once sfuh(i,k) is computed, use with sflh(i,k-1) to determine  saturation !
       ! fraction sfi(i,k) for interfacial layer k-0.5.                            !
       ! This is 'not' chosen for full consistent treatment of stratus fraction in !
       ! all physics schemes.                                                      !
       ! ------------------------------------------------------------------------- !
       DO k = ntop_turb + 1, nbot_turb
          km1 = k - 1
          DO i = 1, ncol
             ! Compute saturation excess at the mid-point of layer k
             sltop    = sl(i,k) + slslope(i,k) * ( pi(i,k) - pm(i,k) )      
             qttop    = qt(i,k) + qtslope(i,k) * ( pi(i,k) - pm(i,k) )
             tltop(1) = ( sltop - g * zi(i,k) ) / cpair 
             status   = fqsatd( tltop(1), pi(i,k), es(1), qs(1), gam(1), 1 )
             qxtop    = qttop - qs(1) 
             slbot    = sl(i,k) + slslope(i,k) * ( pi(i,k+1) - pm(i,k) )      
             qtbot    = qt(i,k) + qtslope(i,k) * ( pi(i,k+1) - pm(i,k) )
             tlbot(1) = ( slbot - g * zi(i,k+1) ) / cpair 
             status   = fqsatd( tlbot(1), pi(i,k+1), es(1), qs(1), gam(1), 1 )
             qxbot    = qtbot - qs(1) 
             qxm      = qxtop + ( qxbot - qxtop ) * ( pm(i,k) - pi(i,k) ) / ( pi(i,k+1) - pi(i,k) )
             ! Find the saturation fraction sfuh(i,k) of the upper half of layer k.
             IF( ( qxtop .LT. 0._r8 ) .AND. ( qxm .LT. 0._r8 ) ) THEN
                sfuh(i,k) = 0._r8 
             ELSE IF( ( qxtop .GT. 0._r8 ) .AND. ( qxm .GT. 0._r8 ) ) THEN
                sfuh(i,k) = 1._r8  
             ELSE ! Either qxm < 0 and qxtop > 0 or vice versa
                sfuh(i,k) = MAX( qxtop, qxm ) / ABS( qxtop - qxm )
             END IF
             ! Combine with sflh(i) (still for layer k-1) to get interfac layer saturation fraction
             sfi(i,k) = 0.5_r8 * ( sflh(i,k-1) + MIN( sflh(i,k-1), sfuh(i,k) ) )
             ! Update sflh to be for the lower half of layer k.             
             IF( ( qxbot .LT. 0._r8 ) .AND. ( qxm .LT. 0._r8 ) ) THEN
                sflh(i,k) = 0._r8 
             ELSE IF( ( qxbot .GT. 0._r8 ) .AND. ( qxm .GT. 0._r8 ) ) THEN
                sflh(i,k) = 1._r8 
             ELSE ! Either qxm < 0 and qxbot > 0 or vice versa
                sflh(i,k) = MAX( qxbot, qxm ) / ABS( qxbot - qxm )
             END IF
          END DO  ! i
       END DO ! k
       DO i = 1, ncol
          sfi(i,pver+1) = sflh(i,pver)  ! Saturation fraction in the lowest half-layer. 
       END DO
    END SELECT

    RETURN
  END SUBROUTINE sfdiag


  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  SUBROUTINE austausch_atm( pcols, pver, ncol, ri, s2, kvf )

    !---------------------------------------------------------------------- ! 
    !                                                                       !
    ! Purpose: Computes exchange coefficients for free turbulent flows.     !
    !          This is not used in the UW moist turbulence scheme.          !
    !                                                                       !
    ! Method:                                                               !
    !                                                                       !
    ! The free atmosphere diffusivities are based on standard mixing length !
    ! forms for the neutral diffusivity multiplied by functns of Richardson !
    ! number. K = l^2 * |dV/dz| * f(Ri). The same functions are used for    !
    ! momentum, potential temperature, and constitutents.                   !
    !                                                                       !
    ! The stable Richardson num function (Ri>0) is taken from Holtslag and  !
    ! Beljaars (1989), ECMWF proceedings. f = 1 / (1 + 10*Ri*(1 + 8*Ri))    !
    ! The unstable Richardson number function (Ri<0) is taken from  CCM1.   !
    ! f = sqrt(1 - 18*Ri)                                                   !
    !                                                                       !
    ! Author: B. Stevens (rewrite, August 2000)                             !
    !                                                                       !
    !---------------------------------------------------------------------- !
    IMPLICIT NONE

    ! --------------- ! 
    ! Input arguments !
    ! --------------- !

    INTEGER,  INTENT(in)  :: pcols                ! Number of atmospheric columns   
    INTEGER,  INTENT(in)  :: pver                 ! Number of atmospheric layers   
    INTEGER,  INTENT(in)  :: ncol                 ! Number of atmospheric columns

    REAL(r8), INTENT(in)  :: s2(pcols,pver)       ! Shear squared
    REAL(r8), INTENT(in)  :: ri(pcols,pver)       ! Richardson no

    ! ---------------- !
    ! Output arguments !
    ! ---------------- !

    REAL(r8), INTENT(out) :: kvf(pcols,pver+1)    ! Eddy diffusivity for heat and tracers

    ! --------------- !
    ! Local Variables !
    ! --------------- !

    REAL(r8)              :: fofri                ! f(ri)
    REAL(r8)              :: kvn                  ! Neutral Kv

    INTEGER               :: i                    ! Longitude index
    INTEGER               :: k                    ! Vertical index

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !
    DO k=1,pver+1
       DO i=1,ncol
          kvf(i,k)           = 0.0_r8
       END DO
    END DO
    fofri=0.0_r8;kvn=0.0_r8
    ! Compute the free atmosphere vertical diffusion coefficients: kvh = kvq = kvm. 

    DO k = ntop_turb + 1, nbot_turb
       DO i = 1, ncol
          IF( ri(i,k) < 0.0_r8 ) THEN
             fofri = SQRT( MAX( 1._r8 - 18._r8 * ri(i,k), 0._r8 ) )
          ELSE 
             fofri = 1.0_r8 / ( 1.0_r8 + 10.0_r8 * ri(i,k) * ( 1.0_r8 + 8.0_r8 * ri(i,k) ) )    
          END IF
          kvn = ml2(k) * SQRT(s2(i,k))
          kvf(i,k) = MAX( zkmin, kvn * fofri )
       END DO
    END DO

    RETURN

  END SUBROUTINE austausch_atm


  ! ---------------------------------------------------------------------------- !
  !                                                                              !
  ! The University of Washington Moist Turbulence Scheme                         !
  !                                                                              !
  ! Authors : Chris Bretherton at the University of Washington, Seattle, WA      ! 
  !           Sungsu Park at the CGD/NCAR, Boulder, CO                           !
  !                                                                              !
  ! ---------------------------------------------------------------------------- !

  SUBROUTINE caleddy( pcols        , pver         , ncol        ,                             &
       sl           , qt           , ql          , slv        , u            , &
       v            , pi           , z           , zi         ,                &
       qflx         , shflx        , slslope     , qtslope    ,                &
       chu          , chs          , cmu         , cms        , sfuh         , &
       sflh         , n2           , s2          , ri         , rrho         , &
       pblh         , ustar        ,                                           &
       kvh_in       , kvm_in       , kvh         , kvm        ,                &
       tpert        , qpert        , qrlin       , kvf        , tke          , & 
       wstarent     , bprod        , sprod       , minpblh    , wpert        , &
       tkes         , turbtype_f   , sm_aw       ,                             &
       kbase_o      , ktop_o       , ncvfin_o    ,                             & 
       kbase_mg     , ktop_mg      , ncvfin_mg   ,                             & 
       kbase_f      , ktop_f       , ncvfin_f    ,                             & 
       wet_CL       , web_CL       , jtbu_CL     , jbbu_CL    ,                &
       evhc_CL      , jt2slv_CL    , n2ht_CL     , n2hb_CL    , lwp_CL       , &
       opt_depth_CL , radinvfrac_CL, radf_CL     , wstar_CL   , wstar3fact_CL, &
       ebrk         , wbrk         , lbrk        , ricl       , ghcl         , & 
       shcl         , smcl         ,                                           &
       gh_a         , sh_a         , sm_a        , ri_a       , leng         , & 
       wcap         , pblhp        , cld         ,landfrac    , ipbl         , &
       kpblh        , wsedl        ,ycoefm  , ycoefh )

    !--------------------------------------------------------------------------------- !
    !                                                                                  !
    ! Purpose : This is a driver routine to compute eddy diffusion coefficients        !
    !           for heat (sl), momentum (u, v), moisture (qt), and other  trace        !
    !           constituents.   This scheme uses first order closure for stable        !
    !           turbulent layers (STL). For convective layers (CL), entrainment        !
    !           closure is used at the CL external interfaces, which is coupled        !
    !           to the diagnosis of a CL regime mean TKE from the instantaneous        !
    !           thermodynamic and velocity profiles.   The CLs are diagnosed by        !
    !           extending original CL layers of moist static instability   into        !
    !           adjacent weakly stably stratified interfaces,   stopping if the        !
    !           stability is too strong.   This allows a realistic depiction of        !
    !           dry convective boundary layers with a downgradient approach.           !
    !                                                                                  !   
    ! NOTE:     This routine currently assumes ntop_turb = 1, nbot_turb = pver         !
    !           ( turbulent diffusivities computed at all interior interfaces )        !
    !           and will require modification to handle a different ntop_turb.         ! 
    !                                                                                  !
    ! Authors:  Sungsu Park and Chris Bretherton. 08/2006, 05/2008.                    !
    !                                                                                  ! 
    ! For details, see                                                                 !
    !                                                                                  !
    ! 1. 'A new moist turbulence parametrization in the Community Atmosphere Model'    !
    !     by Christopher S. Bretherton & Sungsu Park. J. Climate. 22. 3422-3448. 2009. !
    !                                                                                  !
    ! 2. 'The University of Washington shallow convection and moist turbulence schemes !
    !     and their impact on climate simulations with the Community Atmosphere Model' !
    !     by Sungsu Park & Christopher S. Bretherton. J. Climate. 22. 3449-3469. 2009. !
    !                                                                                  !
    ! For questions on the scheme and code, send an email to                           !
    !     sungsup@ucar.edu or breth@washington.edu                                     !
    !                                                                                  !
    !--------------------------------------------------------------------------------- !

    ! ---------------- !
    ! Inputs variables !
    ! ---------------- !

    IMPLICIT NONE

    INTEGER,  INTENT(in) :: pcols                     ! Number of atmospheric columns   
    INTEGER,  INTENT(in) :: pver                      ! Number of atmospheric layers   
    INTEGER,  INTENT(in) :: ncol                      ! Number of atmospheric columns   
    REAL(r8), INTENT(in) :: u(pcols,pver)             ! U wind [ m/s ]
    REAL(r8), INTENT(in) :: v(pcols,pver)             ! V wind [ m/s ]
    REAL(r8), INTENT(in) :: sl(pcols,pver)            ! Liquid water static energy, cp * T + g * z - Lv * ql - Ls * qi [ J/kg ]
    REAL(r8), INTENT(in) :: slv(pcols,pver)           ! Liquid water virtual static energy, sl * ( 1 + 0.608 * qt ) [ J/kg ]
    REAL(r8), INTENT(in) :: qt(pcols,pver)            ! Total speccific humidity  qv + ql + qi [ kg/kg ] 
    REAL(r8), INTENT(in) :: ql(pcols,pver)            ! Liquid water specific humidity [ kg/kg ]
    REAL(r8), INTENT(in) :: pi(pcols,pver+1)          ! Interface pressures [ Pa ]
    REAL(r8), INTENT(in) :: z(pcols,pver)             ! Layer midpoint height above surface [ m ]
    REAL(r8), INTENT(in) :: zi(pcols,pver+1)          ! Interface height above surface, i.e., zi(pver+1) = 0 all over the globe [ m ]
    REAL(r8), INTENT(in) :: chu(pcols,pver+1)         ! Buoyancy coeffi. unsaturated sl (heat) coef. at all interfaces. [ unit ? ]
    REAL(r8), INTENT(in) :: chs(pcols,pver+1)         ! Buoyancy coeffi. saturated sl (heat) coef. at all interfaces. [ unit ? ]
    REAL(r8), INTENT(in) :: cmu(pcols,pver+1)         ! Buoyancy coeffi. unsaturated qt (moisture) coef. at all interfaces [ unit ? ]
    REAL(r8), INTENT(in) :: cms(pcols,pver+1)         ! Buoyancy coeffi. saturated qt (moisture) coef. at all interfaces [ unit ? ]
    REAL(r8), INTENT(in) :: sfuh(pcols,pver)          ! Saturation fraction in upper half-layer [ fraction ]
    REAL(r8), INTENT(in) :: sflh(pcols,pver)          ! Saturation fraction in lower half-layer [ fraction ]
    REAL(r8), INTENT(in) :: n2(pcols,pver)            ! Interfacial (except surface) moist buoyancy frequency [ s-2 ]
    REAL(r8), INTENT(in) :: s2(pcols,pver)            ! Interfacial (except surface) shear frequency [ s-2 ]
    REAL(r8), INTENT(in) :: ri(pcols,pver)            ! Interfacial (except surface) Richardson number
    REAL(r8), INTENT(in) :: qflx(pcols)               ! Kinematic surface constituent ( water vapor ) flux [ kg/m2/s ]
    REAL(r8), INTENT(in) :: shflx(pcols)              ! Kinematic surface heat flux [ unit ? ] 
    REAL(r8), INTENT(in) :: slslope(pcols,pver)       ! Slope of 'sl' in each layer [ J/kg/Pa ]
    REAL(r8), INTENT(in) :: qtslope(pcols,pver)       ! Slope of 'qt' in each layer [ kg/kg/Pa ]
    REAL(r8), INTENT(in) :: qrlin(pcols,pver)         ! Input grid-mean LW heating rate : [ K/s ] * cpair * dp = [ W/kg*Pa ]
    REAL(r8), INTENT(in) :: wsedl(pcols,pver)         ! Sedimentation velocity of liquid stratus cloud droplet [ m/s ]
    REAL(r8), INTENT(in) :: ustar(pcols)              ! Surface friction velocity [ m/s ]
    REAL(r8), INTENT(in) :: rrho(pcols)               ! 1./bottom mid-point density. Specific volume [ m3/kg ]
    REAL(r8), INTENT(in) :: kvf(pcols,pver+1)         ! Free atmosphere eddy diffusivity [ m2/s ]
    LOGICAL,  INTENT(in) :: wstarent                  ! Switch for choosing wstar3 entrainment parameterization
    REAL(r8), INTENT(in) :: minpblh(pcols)            ! Minimum PBL height based on surface stress [ m ]
    REAL(r8), INTENT(in) :: kvh_in(pcols,pver+1)      ! kvh saved from last timestep or last iterative step [ m2/s ] 
    REAL(r8), INTENT(in) :: kvm_in(pcols,pver+1)      ! kvm saved from last timestep or last iterative step [ m2/s ]
    REAL(r8), INTENT(in) :: cld(pcols,pver)           ! Stratus Cloud Fraction [ fraction ]
    REAL(r8), INTENT(in) :: landfrac(pcols)   
    REAL(r8), INTENT(in) :: ycoefm(pcols,pver+1) 
    REAL(r8), INTENT(in) :: ycoefh(pcols,pver+1) 
    ! ---------------- !
    ! Output variables !
    ! ---------------- !

    REAL(r8), INTENT(out) :: kvh(pcols,pver+1)        ! Eddy diffusivity for heat, moisture, and tracers [ m2/s ]
    REAL(r8), INTENT(out) :: kvm(pcols,pver+1)        ! Eddy diffusivity for momentum [ m2/s ]
    REAL(r8), INTENT(out) :: pblh(pcols)              ! PBL top height [ m ]
    REAL(r8), INTENT(out) :: pblhp(pcols)             ! PBL top height pressure [ Pa ]
    REAL(r8), INTENT(out) :: tpert(pcols)             ! Convective temperature excess [ K ]
    REAL(r8), INTENT(out) :: qpert(pcols)             ! Convective humidity excess [ kg/kg ]
    REAL(r8), INTENT(out) :: wpert(pcols)             ! Turbulent velocity excess [ m/s ]
    REAL(r8), INTENT(out) :: tke(pcols,pver+1)        ! Turbulent kinetic energy [ m2/s2 ], 'tkes' at surface, pver+1.
    REAL(r8), INTENT(out) :: bprod(pcols,pver+1)      ! Buoyancy production [ m2/s3 ],     'bflxs' at surface, pver+1.
    REAL(r8), INTENT(out) :: sprod(pcols,pver+1)      ! Shear production [ m2/s3 ], (ustar(i)**3)/(vk*z(i,pver))  at surface, pver+1.
    REAL(r8), INTENT(out) :: turbtype_f(pcols,pver+1) ! Turbulence type at each interface:
    ! 0. = Non turbulence interface
    ! 1. = Stable turbulence interface
    ! 2. = CL interior interface ( if bflxs > 0, surface is this )
    ! 3. = Bottom external interface of CL
    ! 4. = Top external interface of CL.
    ! 5. = Double entraining CL external interface 
    REAL(r8), INTENT(out) :: sm_aw(pcols,pver+1)      ! Galperin instability function of momentum for use in the microphysics [ no unit ]
    REAL(r8), INTENT(out) :: ipbl(pcols)              ! If 1, PBL is CL, while if 0, PBL is STL.
    REAL(r8), INTENT(out) :: kpblh(pcols)             ! Layer index containing PBL within or at the base interface

    ! --------------------------- !
    ! Diagnostic output variables !
    ! --------------------------- !

    REAL(r8) :: tkes(pcols)                           ! TKE at surface [ m2/s2 ] 
    REAL(r8) :: kbase_o(pcols,ncvmax)                 ! Original external base interface index of CL just after 'exacol'
    REAL(r8) :: ktop_o(pcols,ncvmax)                  ! Original external top  interface index of CL just after 'exacol'
    REAL(r8) :: ncvfin_o(pcols)                       ! Original number of CLs just after 'exacol'
    REAL(r8) :: kbase_mg(pcols,ncvmax)                ! kbase  just after extending-merging (after 'zisocl') but without SRCL
    REAL(r8) :: ktop_mg(pcols,ncvmax)                 ! ktop   just after extending-merging (after 'zisocl') but without SRCL
    REAL(r8) :: ncvfin_mg(pcols)                      ! ncvfin just after extending-merging (after 'zisocl') but without SRCL
    REAL(r8) :: kbase_f(pcols,ncvmax)                 ! Final kbase  after adding SRCL
    REAL(r8) :: ktop_f(pcols,ncvmax)                  ! Final ktop   after adding SRCL
    REAL(r8) :: ncvfin_f(pcols)                       ! Final ncvfin after adding SRCL
    REAL(r8) :: wet_CL(pcols,ncvmax)                  ! Entrainment rate at the CL top [ m/s ] 
    REAL(r8) :: web_CL(pcols,ncvmax)                  ! Entrainment rate at the CL base [ m/s ]
    REAL(r8) :: jtbu_CL(pcols,ncvmax)                 ! Buoyancy jump across the CL top [ m/s2 ]  
    REAL(r8) :: jbbu_CL(pcols,ncvmax)                 ! Buoyancy jump across the CL base [ m/s2 ]  
    REAL(r8) :: evhc_CL(pcols,ncvmax)                 ! Evaporative enhancement factor at the CL top
    REAL(r8) :: jt2slv_CL(pcols,ncvmax)               ! Jump of slv ( across two layers ) at CL top for use only in evhc [ J/kg ]
    REAL(r8) :: n2ht_CL(pcols,ncvmax)                 ! n2 defined at the CL top  interface but using sfuh(kt)   instead of sfi(kt) [ s-2 ]
    REAL(r8) :: n2hb_CL(pcols,ncvmax)                 ! n2 defined at the CL base interface but using sflh(kb-1) instead of sfi(kb) [ s-2 ]
    REAL(r8) :: lwp_CL(pcols,ncvmax)                  ! LWP in the CL top layer [ kg/m2 ]
    REAL(r8) :: opt_depth_CL(pcols,ncvmax)            ! Optical depth of the CL top layer
    REAL(r8) :: radinvfrac_CL(pcols,ncvmax)           ! Fraction of LW radiative cooling confined in the top portion of CL
    REAL(r8) :: radf_CL(pcols,ncvmax)                 ! Buoyancy production at the CL top due to radiative cooling [ m2/s3 ]
    REAL(r8) :: wstar_CL(pcols,ncvmax)                ! Convective velocity of CL including entrainment contribution finally [ m/s ]
    REAL(r8) :: wstar3fact_CL(pcols,ncvmax)           ! "wstar3fact" of CL. Entrainment enhancement of wstar3 (inverse)

    REAL(r8) :: gh_a(pcols,pver+1)                    ! Half of normalized buoyancy production, -l2n2/2e. [ no unit ]
    REAL(r8) :: sh_a(pcols,pver+1)                    ! Galperin instability function of heat-moisture at all interfaces [ no unit ]
    REAL(r8) :: sm_a(pcols,pver+1)                    ! Galperin instability function of momentum      at all interfaces [ no unit ]
    REAL(r8) :: ri_a(pcols,pver+1)                    ! Interfacial Richardson number                  at all interfaces [ no unit ]

    REAL(r8) :: ebrk(pcols,ncvmax)                    ! Net CL mean TKE [ m2/s2 ]
    REAL(r8) :: wbrk(pcols,ncvmax)                    ! Net CL mean normalized TKE [ m2/s2 ]
    REAL(r8) :: lbrk(pcols,ncvmax)                    ! Net energetic integral thickness of CL [ m ]
    REAL(r8) :: ricl(pcols,ncvmax)                    ! Mean Richardson number of CL ( l2n2/l2s2 )
    REAL(r8) :: ghcl(pcols,ncvmax)                    ! Half of normalized buoyancy production of CL                 
    REAL(r8) :: shcl(pcols,ncvmax)                    ! Instability function of heat and moisture of CL
    REAL(r8) :: smcl(pcols,ncvmax)                    ! Instability function of momentum of CL

    REAL(r8) :: leng(pcols,pver+1)                    ! Turbulent length scale [ m ], 0 at the surface.
    REAL(r8) :: wcap(pcols,pver+1)                    ! Normalized TKE [m2/s2], 'tkes/b1' at the surface and 'tke/b1' at
    ! the top/bottom entrainment interfaces of CL assuming no transport.
    ! ------------------------ !
    ! Local Internal Variables !
    ! ------------------------ !

    LOGICAL :: belongcv(pcols,pver+1)                 ! True for interfaces in a CL (both interior and exterior are included)
    LOGICAL :: belongst(pcols,pver+1)                 ! True for stable turbulent layer interfaces (STL)
    LOGICAL :: in_CL                                  ! True if interfaces k,k+1 both in same CL.
    LOGICAL :: extend                                 ! True when CL is extended in zisocl
    LOGICAL :: extend_up                              ! True when CL is extended upward in zisocl
    LOGICAL :: extend_dn                              ! True when CL is extended downward in zisocl

    INTEGER :: i,m                                      ! Longitude index
    INTEGER :: k                                      ! Vertical index
    INTEGER :: ks                                     ! Vertical index
    INTEGER :: ncvfin(pcols)                          ! Total number of CL in column
    INTEGER :: ncvf                                   ! Total number of CL in column prior to adding SRCL
    INTEGER :: ncv                                    ! Index of current CL
    INTEGER :: ncvnew                                 ! Index of added SRCL appended after regular CLs from 'zisocl'
    INTEGER :: ncvsurf                                ! If nonzero, CL index based on surface (usually 1, but can be > 1 when SRCL is based at sfc)
    INTEGER :: kbase(pcols,ncvmax)                    ! Vertical index of CL base interface
    INTEGER :: ktop(pcols,ncvmax)                     ! Vertical index of CL top interface
    INTEGER :: kb, kt                                 ! kbase and ktop for current CL
    INTEGER :: ktblw                                  ! ktop of the CL located at just below the current CL
    INTEGER :: turbtype(pcols,pver+1)                 ! Interface turbulence type :
    ! 0 = Non turbulence interface
    ! 1 = Stable turbulence interface
    ! 2 = CL interior interface ( if bflxs > 0, sfc is this )
    ! 3 = Bottom external interface of CL
    ! 4 = Top external interface of CL
    ! 5 = Double entraining CL external interface
    INTEGER  :: ktopbl(pcols)                         ! PBL top height or interface index 
    REAL(r8) :: bflxs(pcols)                          ! Surface buoyancy flux [ m2/s3 ]
    REAL(r8) :: rcap                                  ! 'tke/ebrk' at all interfaces of CL. Set to 1 at the CL entrainment interfaces
    REAL(r8) :: jtzm                                  ! Interface layer thickness of CL top interface [ m ]
    REAL(r8) :: jtsl                                  ! Jump of s_l across CL top interface [ J/kg ]
    REAL(r8) :: jtqt                                  ! Jump of q_t across CL top interface [ kg/kg ]
    REAL(r8) :: jtbu                                  ! Jump of buoyancy across CL top interface [ m/s2 ]
    REAL(r8) :: jtu                                   ! Jump of u across CL top interface [ m/s ]
    REAL(r8) :: jtv                                   ! Jump of v across CL top interface [ m/s ]
    REAL(r8) :: jt2slv                                ! Jump of slv ( across two layers ) at CL top for use only in evhc [ J/kg ]
    REAL(r8) :: radf                                  ! Buoyancy production at the CL top due to radiative cooling [ m2/s3 ]
    REAL(r8) :: jbzm                                  ! Interface layer thickness of CL base interface [ m ]
    REAL(r8) :: jbsl                                  ! Jump of s_l across CL base interface [ J/kg ]
    REAL(r8) :: jbqt                                  ! Jump of q_t across CL top interface [ kg/kg ]
    REAL(r8) :: jbbu                                  ! Jump of buoyancy across CL base interface [ m/s2 ]
    REAL(r8) :: jbu                                   ! Jump of u across CL base interface [ m/s ]
    REAL(r8) :: jbv                                   ! Jump of v across CL base interface [ m/s ]
    REAL(r8) :: ch                                    ! Buoyancy coefficients defined at the CL top and base interfaces using CL internal
    REAL(r8) :: cm                                    ! sfuh(kt) and sflh(kb-1) instead of sfi(kt) and sfi(kb), respectively. These are 
    ! used for entrainment calculation at CL external interfaces and SRCL identification.
    REAL(r8) :: n2ht                                  ! n2 defined at the CL top  interface but using sfuh(kt)   instead of sfi(kt) [ s-2 ]
    REAL(r8) :: n2hb                                  ! n2 defined at the CL base interface but using sflh(kb-1) instead of sfi(kb) [ s-2 ]
    REAL(r8) :: n2htSRCL                              ! n2 defined at the upper-half layer of SRCL. This is used only for identifying SRCL.
    ! n2htSRCL use SRCL internal slope sl and qt as well as sfuh(kt) instead of sfi(kt) [ s-2 ]
    REAL(r8) :: gh                                    ! Half of normalized buoyancy production ( -l2n2/2e ) [ no unit ]
    REAL(r8) :: sh                                    ! Galperin instability function for heat and moisture
    REAL(r8) :: sm                                    ! Galperin instability function for momentum
    REAL(r8) :: lbulk                                 ! Depth of turbulent layer, Master length scale (not energetic length)
    REAL(r8) :: dzht                                  ! Thickness of top    half-layer [ m ]
    REAL(r8) :: dzhb                                  ! Thickness of bottom half-layer [ m ]
    REAL(r8) :: rootp                                 ! Sqrt(net CL-mean TKE including entrainment contribution) [ m/s ]     
    REAL(r8) :: evhc                                  ! Evaporative enhancement factor: (1+E) with E = evap. cool. efficiency [ no unit ]
    REAL(r8) :: kentr                                 ! Effective entrainment diffusivity 'wet*dz', 'web*dz' [ m2/s ]
    REAL(r8) :: lwp                                   ! Liquid water path in the layer kt [ kg/m2 ]
    REAL(r8) :: opt_depth                             ! Optical depth of the layer kt [ no unit ]
    REAL(r8) :: radinvfrac                            ! Fraction of LW cooling in the layer kt concentrated at the CL top [ no unit ]
    REAL(r8) :: wet                                   ! CL top entrainment rate [ m/s ]
    REAL(r8) :: web                                   ! CL bot entrainment rate [ m/s ]. Set to zero if CL is based at surface.
    REAL(r8) :: vyt                                   ! n2ht/n2 at the CL top  interface
    REAL(r8) :: vyb                                   ! n2hb/n2 at the CL base interface
    REAL(r8) :: vut                                   ! Inverse Ri (=s2/n2) at the CL top  interface
    REAL(r8) :: vub                                   ! Inverse Ri (=s2/n2) at the CL base interface
    REAL(r8) :: fact                                  ! Factor relating TKE generation to entrainment [ no unit ]
    REAL(r8) :: trma                                  ! Intermediate variables used for solving quadratic ( for gh from ri )
    REAL(r8) :: trmb                                  ! and cubic equations ( for ebrk: the net CL mean TKE )
    REAL(r8) :: trmc                                  !
    REAL(r8) :: trmp                                  !
    REAL(r8) :: trmq                                  !
    REAL(r8) :: qq                                    ! 
    REAL(r8) :: det                                   !
    REAL(r8) :: gg                                    ! Intermediate variable used for calculating stability functions of
    ! SRCL or SBCL based at the surface with bflxs > 0.
    REAL(r8) :: dzhb5                                 ! Half thickness of the bottom-most layer of current CL regime
    REAL(r8) :: dzht5                                 ! Half thickness of the top-most layer of adjacent CL regime just below current CL
    REAL(r8) :: qrlw(pcols,pver)                      ! Local grid-mean LW heating rate : [K/s] * cpair * dp = [ W/kg*Pa ]

    REAL(r8) :: cldeff(pcols,pver)                    ! Effective stratus fraction
    REAL(r8) :: qleff                                 ! Used for computing evhc
    REAL(r8) :: tunlramp                              ! Ramping tunl
    REAL(r8) :: leng_imsi                             ! For Kv = max(Kv_STL, Kv_entrain)
    REAL(r8) :: tke_imsi                              !
    REAL(r8) :: kvh_imsi                              !
    REAL(r8) :: kvm_imsi                              !
    REAL(r8) :: alph4exs                              ! For extended stability function in the stable regime
    REAL(r8) :: ghmin                                 !   

    REAL(r8) :: sedfact                               ! For 'sedimentation-entrainment feedback' 

    ! Local variables specific for 'wstar' entrainment closure

    REAL(r8) :: cet                                   ! Proportionality coefficient between wet and wstar3
    REAL(r8) :: ceb                                   ! Proportionality coefficient between web and wstar3
    REAL(r8) :: wstar                                 ! Convective velocity for CL [ m/s ]
    REAL(r8) :: wstar3                                ! Cubed convective velocity for CL [ m3/s3 ]
    REAL(r8) :: wstar3fact                            ! 1/(relative change of wstar^3 by entrainment)
    REAL(r8) :: rmin                                  ! sqrt(p)
    REAL(r8) :: fmin                                  ! f(rmin), where f(r) = r^3 - 3*p*r - 2q
    REAL(r8) :: rcrit                                 ! ccrit*wstar
    REAL(r8) :: fcrit                                 ! f(rcrit)
    LOGICAL     noroot                                ! True if f(r) has no root r > rcrit

    !-----------------------!
    ! Start of Main Program !
    !-----------------------!
    DO i = 1, ncol
       ipbl (i)=0.0_r8! If 1, PBL is CL, while if 0, PBL is STL.
       kpblh(i)=0.0_r8! Layer index containing PBL within or at the base interface
       pblh (i)=0.0_r8! PBL top height [ m ]
       pblhp(i)=0.0_r8! PBL top height pressure [ Pa ]
       tpert(i)=0.0_r8! Convective temperature excess [ K ]
       qpert(i)=0.0_r8! Convective humidity excess [ kg/kg ]
       wpert(i)=0.0_r8! Turbulent velocity excess [ m/s ]
       ktopbl(i) =0   ! PBL top height or interface index 
       bflxs(i)=0.0_r8! Surface buoyancy flux [ m2/s3 ]
       ncvfin(i)=0 ! Total number of CL in column
       tkes(i)    =0.0_r8                    ! TKE at surface [ m2/s2 ] 
       ncvfin_o(i)=0.0_r8! Original number of CLs just after 'exacol'
       ncvfin_mg(i)=0.0_r8! ncvfin just after extending-merging (after 'zisocl') but without SRCL
       ncvfin_f(i)=0.0_r8! Final ncvfin after adding SRCL
    END DO

    DO m = 1, ncvmax
       DO i = 1, ncol   
           kbase  (i,m)   =0                 ! Vertical index of CL base interface
           ktop   (i,m)   =0                 ! Vertical index of CL top interface
           kbase_o  (i,m) =0.0_r8! Original external base interface index of CL just after 'exacol'
           ktop_o   (i,m) =0.0_r8! Original external top  interface index of CL just after 'exacol'
           kbase_mg (i,m) =0.0_r8! kbase  just after extending-merging (after 'zisocl') but without SRCL
           ktop_mg  (i,m) =0.0_r8! ktop   just after extending-merging (after 'zisocl') but without SRCL
           kbase_f  (i,m) =0.0_r8! Final kbase  after adding SRCL
           ktop_f   (i,m) =0.0_r8! Final ktop   after adding SRCL
           wet_CL   (i,m) =0.0_r8! Entrainment rate at the CL top [ m/s ] 
           web_CL   (i,m) =0.0_r8! Entrainment rate at the CL base [ m/s ]
           jtbu_CL  (i,m) =0.0_r8! Buoyancy jump across the CL top [ m/s2 ]  
           jbbu_CL  (i,m) =0.0_r8! Buoyancy jump across the CL base [ m/s2 ]  
           evhc_CL  (i,m) =0.0_r8! Evaporative enhancement factor at the CL top
           jt2slv_CL(i,m) =0.0_r8! Jump of slv ( across two layers ) at CL top for use only in evhc [ J/kg ]
           n2ht_CL  (i,m) =0.0_r8! n2 defined at the CL top  interface but using sfuh(kt)   instead of sfi(kt) [ s-2 ]
           n2hb_CL  (i,m) =0.0_r8! n2 defined at the CL base interface but using sflh(kb-1) instead of sfi(kb) [ s-2 ]
           lwp_CL   (i,m) =0.0_r8! LWP in the CL top layer [ kg/m2 ]
           opt_depth_CL(i,m) =0.0_r8! Optical depth of the CL top layer
           radinvfrac_CL(i,m) =0.0_r8! Fraction of LW radiative cooling confined in the top portion of CL
           radf_CL(i,m) =0.0_r8! Buoyancy production at the CL top due to radiative cooling [ m2/s3 ]
           wstar_CL(i,m) =0.0_r8! Convective velocity of CL including entrainment contribution finally [ m/s ]
           wstar3fact_CL(i,m) =0.0_r8! "wstar3fact" of CL. Entrainment enhancement of wstar3 (inverse)
           ebrk(i,m) =0.0_r8! Net CL mean TKE [ m2/s2 ]
           wbrk(i,m) =0.0_r8! Net CL mean normalized TKE [ m2/s2 ]
           lbrk(i,m) =0.0_r8! Net energetic integral thickness of CL [ m ]
           ricl(i,m) =0.0_r8! Mean Richardson number of CL ( l2n2/l2s2 )
           ghcl(i,m) =0.0_r8! Half of normalized buoyancy production of CL
           shcl(i,m) =0.0_r8! Instability function of heat and moisture of CL
           smcl(i,m) =0.0_r8! Instability function of momentum of CL

       END DO
    END DO

!!!!!!!!!!!!!!!!!!

    ! 0. = Non turbulence interface
    ! 1. = Stable turbulence interface
    ! 2. = CL interior interface ( if bflxs > 0, surface is this )
    ! 3. = Bottom external interface of CL
    ! 4. = Top external interface of CL.
    ! 5. = Double entraining CL external interface 
    DO k = 1, pver+1
       DO i = 1, ncol
          kvh       (i,k)=0.0_r8! Eddy diffusivity for heat, moisture, and tracers [ m2/s ]
          kvm       (i,k)=0.0_r8! Eddy diffusivity for momentum [ m2/s ]
          tke       (i,k)=0.0_r8! Turbulent kinetic energy [ m2/s2 ], 'tkes' at surface, pver+1.
          bprod     (i,k)=0.0_r8! Buoyancy production [ m2/s3 ],'bflxs' at surface, pver+1.
          sprod     (i,k)=0.0_r8! Shear production [ m2/s3 ], (ustar(i)**3)/(vk*z(i,pver))  at surface, pver+1.
          turbtype_f(i,k)=0.0_r8! Turbulence type at each interface:
          sm_aw     (i,k)=0.0_r8! Galperin instability function of momentum for use in the microphysics [ no unit ]
          turbtype  (i,k)=0     ! Interface turbulence type :
          gh_a      (i,k)=0.0_r8  ! Half of normalized buoyancy production, -l2n2/2e. [ no unit ]
          sh_a      (i,k)=0.0_r8 ! Galperin instability function of heat-moisture at all interfaces [ no unit ]
          sm_a      (i,k)=0.0_r8 ! Galperin instability function of momentum at all interfaces [ no unit ]
          ri_a      (i,k)=0.0_r8  ! Interfacial Richardson number    at all interfaces [ no unit ]
          leng      (i,k)=0.0_r8! Turbulent length scale [ m ], 0 at the surface.
          wcap      (i,k)=0.0_r8! Normalized TKE [m2/s2], 'tkes/b1' at the surface and 'tke/b1' at

       END DO
    END DO
!!!!!!!!!!!!!!!!!!
    DO k = 1, pver
       DO i = 1, ncol
          qrlw  (i,k)=0.0_r8! Local grid-mean LW heating rate : [K/s] * cpair * dp = [ W/kg*Pa ]
          cldeff(i,k)=0.0_r8! Effective stratus fraction
       END DO
    END DO


    ! the top/bottom entrainment interfaces of CL assuming no transport.
    ! ------------------------ !
    ! Local Internal Variables !
    ! ------------------------ !

   ! LOGICAL :: belongcv(pcols,pver+1)                 ! True for interfaces in a CL (both interior and exterior are included)
   ! LOGICAL :: belongst(pcols,pver+1)                 ! True for stable turbulent layer interfaces (STL)
   ! LOGICAL :: in_CL                                  ! True if interfaces k,k+1 both in same CL.
   ! LOGICAL :: extend                                 ! True when CL is extended in zisocl
   ! LOGICAL :: extend_up                              ! True when CL is extended upward in zisocl
   ! LOGICAL :: extend_dn                              ! True when CL is extended downward in zisocl
   ! INTEGER :: i                                      ! Longitude index
   ! INTEGER :: k                                      ! Vertical index
   ! INTEGER :: ks                                     ! Vertical index
   ncvf=0;ncv =0;ncvnew=0;
   ncvsurf =0;kb=0 ; kt =0;ktblw =0;
   rcap=0.0_r8;jtzm=0.0_r8;
   jtsl=0.0_r8;jtqt=0.0_r8;jtbu=0.0_r8;jtu=0.0_r8;jtv=0.0_r8;jt2slv=0.0_r8;
   radf=0.0_r8;jbzm=0.0_r8;jbsl=0.0_r8;jbqt=0.0_r8;jbbu=0.0_r8;jbu=0.0_r8;
   jbv=0.0_r8;ch=0.0_r8;cm=0.0_r8;n2ht=0.0_r8;n2hb=0.0_r8;n2htSRCL=0.0_r8;
   gh=0.0_r8;sh=0.0_r8;sm=0.0_r8;lbulk=0.0_r8;dzht=0.0_r8;dzhb=0.0_r8;rootp=0.0_r8;
   evhc=0.0_r8;kentr=0.0_r8;lwp=0.0_r8;opt_depth=0.0_r8;radinvfrac=0.0_r8;wet=0.0_r8;
   web=0.0_r8;vyt=0.0_r8;vyb=0.0_r8;vut=0.0_r8;vub=0.0_r8;fact=0.0_r8;trma=0.0_r8;
   trmb=0.0_r8;trmc=0.0_r8;trmp=0.0_r8;trmq=0.0_r8;qq=0.0_r8;det=0.0_r8;gg=0.0_r8;dzhb5=0.0_r8;
   dzht5=0.0_r8;qleff=0.0_r8;tunlramp=0.0_r8;leng_imsi=0.0_r8;tke_imsi=0.0_r8;kvh_imsi=0.0_r8;
   kvm_imsi=0.0_r8;alph4exs=0.0_r8;ghmin=0.0_r8;sedfact=0.0_r8;cet=0.0_r8;ceb=0.0_r8;wstar=0.0_r8;
   wstar3=0.0_r8;wstar3fact=0.0_r8;rmin=0.0_r8;fmin=0.0_r8;rcrit=0.0_r8;fcrit=0.0_r8;

    !----------
    !-------------------
    ! Option: Turn-off LW radiative-turbulence interaction in PBL scheme
    !         by setting qrlw = 0.  Logical parameter 'set_qrlzero'  was
    !         defined in the first part of 'eddy_diff.F90' module. 

    IF( set_qrlzero ) THEN
       qrlw(:,:) = 0._r8
    ELSE
    DO k = 1, pver
       DO i = 1, ncol
       qrlw(i,k) = qrlin(i,k)
       END DO
    END DO
    ENDIF

    ! Define effective stratus fraction using the grid-mean ql.
    ! Modification : The contribution of ice should be carefully considered.
    !                This should be done in combination with the 'qrlw' and
    !                overlapping assumption of liquid and ice stratus. 

    DO k = 1, pver
       DO i = 1, ncol
          IF( choice_evhc .EQ. 'ramp' .OR. choice_radf .EQ. 'ramp' ) THEN 
             cldeff(i,k) = cld(i,k) * MIN( ql(i,k) / qmin(2), 1._r8 )
          ELSE
             cldeff(i,k) = cld(i,k)
          ENDIF
       END DO
    END DO

    ! For an extended stability function in the stable regime, re-define
    ! alph4exe and ghmin. This is for future work.

    IF( ricrit .EQ. 0.19_r8 ) THEN
       alph4exs = alph4
       ghmin    = -3.5334_r8
    ELSEIF( ricrit .GT. 0.19_r8 ) THEN
       alph4exs = -2._r8 * b1 * alph2 / ( alph3 - 2._r8 * b1 * alph5 ) / ricrit
       ghmin    = -1.e10_r8
    ELSE
       WRITE(iulog,*) 'Error : ricrit should be larger than 0.19 in UW PBL'       
       STOP
    ENDIF

    !
    ! Initialization of Diagnostic Output
    !
    DO m = 1, ncvmax

       DO i = 1, ncol
          wet_CL(i,m)        = 0._r8
          web_CL(i,m)        = 0._r8
          jtbu_CL(i,m)       = 0._r8
          jbbu_CL(i,m)       = 0._r8
          evhc_CL(i,m)       = 0._r8
          jt2slv_CL(i,m)     = 0._r8
          n2ht_CL(i,m)       = 0._r8
          n2hb_CL(i,m)       = 0._r8                    
          lwp_CL(i,m)        = 0._r8
          opt_depth_CL(i,m)  = 0._r8
          radinvfrac_CL(i,m) = 0._r8
          radf_CL(i,m)       = 0._r8
          wstar_CL(i,m)      = 0._r8          
          wstar3fact_CL(i,m) = 0._r8
          ricl(i,m)          = 0._r8
          ghcl(i,m)          = 0._r8
          shcl(i,m)          = 0._r8
          smcl(i,m)          = 0._r8
          ebrk(i,m)          = 0._r8
          wbrk(i,m)          = 0._r8
          lbrk(i,m)          = 0._r8
       END DO
    END DO 
    DO k = 1, pver+1
       DO i = 1, ncol
          gh_a(i,k)          = 0._r8
          sh_a(i,k)          = 0._r8
          sm_a(i,k)          = 0._r8
          ri_a(i,k)          = 0._r8
          sm_aw(i,k)         = 0._r8
       END DO
    END DO 

    DO i = 1, ncol
       ipbl(i)= 0._r8
       kpblh(i)= REAL(pver,r8)
    END DO 

    ! kvh and kvm are stored over timesteps in 'vertical_diffusion.F90' and 
    ! passed in as kvh_in and kvm_in.  However,  at the first timestep they
    ! need to be computed and these are done just before calling 'caleddy'.   
    ! kvm and kvh are also stored over iterative time step in the first part
    ! of 'eddy_diff.F90'

    DO k = 1, pver + 1
       DO i = 1, ncol
          ! Initialize kvh and kvm to zero or kvf
          IF( use_kvf ) THEN
             kvh(i,k) = kvf(i,k)
             kvm(i,k) = kvf(i,k)
          ELSE
             kvh(i,k) =  0._r8!ycoefm  (i,k) 
             kvm(i,k) =  0._r8!ycoefh  (i,k) 
          END IF
          ! Zero diagnostic quantities for the new diffusion step.
          wcap(i,k) = 0._r8
          leng(i,k) = 0._r8
          tke(i,k)  = 0._r8
          turbtype(i,k) = 0
       END DO
    END DO

    ! Initialize 'bprod' [ m2/s3 ] and 'sprod' [ m2/s3 ] at all interfaces.
    ! Note this initialization is a hybrid initialization since 'n2' [s-2] and 's2' [s-2]
    ! are calculated from the given current initial profile, while 'kvh_in' [m2/s] and 
    ! 'kvm_in' [m2/s] are from the previous iteration or previous time step.
    ! This initially guessed 'bprod' and 'sprod' will be updated at the end of this 
    ! 'caleddy' subroutine for diagnostic output.
    ! This computation of 'brpod,sprod' below is necessary for wstar-based entrainment closure.

    DO k = 2, pver
       DO i = 1, ncol
          bprod(i,k) = -kvh_in(i,k) * n2(i,k)
          sprod(i,k) =  kvm_in(i,k) * s2(i,k)
       END DO
    END DO

    ! Set 'bprod' and 'sprod' at top and bottom interface.
    ! In calculating 'surface' (actually lowest half-layer) buoyancy flux,
    ! 'chu' at surface is defined to be the same as 'chu' at the mid-point
    ! of lowest model layer (pver) at the end of 'trbind'. The same is for
    ! the other buoyancy coefficients.  'sprod(i,pver+1)'  is defined in a
    ! consistent way as the definition of 'tkes' in the original code.
    ! ( Important Option ) If I want to isolate surface buoyancy flux from
    ! the other parts of CL regimes energetically even though bflxs > 0,
    ! all I should do is to re-define 'bprod(i,pver+1)=0' in the below 'do'
    ! block. Additionally for merging test of extending SBCL based on 'l2n2'
    ! in 'zisocl', I should use 'l2n2 = - wint / sh'  for similar treatment
    ! as previous code. All other parts of the code  are fully consistently
    ! treated by these change only.
    ! My future general convection scheme will use bflxs(i).

    DO i = 1, ncol
       bprod(i,1) = 0._r8 ! Top interface
       sprod(i,1) = 0._r8 ! Top interface
       ch = chu(i,pver+1) * ( 1._r8 - sflh(i,pver) ) + chs(i,pver+1) * sflh(i,pver)   
       cm = cmu(i,pver+1) * ( 1._r8 - sflh(i,pver) ) + cms(i,pver+1) * sflh(i,pver)   
       bflxs(i) = ch * shflx(i) * rrho(i) + cm * qflx(i) * rrho(i)
       IF( choice_tkes .EQ. 'ibprod' ) THEN
          bprod(i,pver+1) = bflxs(i)
       ELSE
          bprod(i,pver+1) = 0._r8
       ENDIF
       sprod(i,pver+1) = (ustar(i)**3)/(vk*z(i,pver))
    END DO

    ! Initially identify CL regimes in 'exacol'
    !    ktop  : Interface index of the CL top  external interface
    !    kbase : Interface index of the CL base external interface
    !    ncvfin: Number of total CLs
    ! Note that if surface buoyancy flux is positive ( bflxs = bprod(i,pver+1) > 0 ),
    ! surface interface is identified as an internal interface of CL. However, even
    ! though bflxs <= 0, if 'pver' interface is a CL internal interface (ri(pver)<0),
    ! surface interface is identified as an external interface of CL. If bflxs =< 0 
    ! and ri(pver) >= 0, then surface interface is identified as a stable turbulent
    ! intereface (STL) as shown at the end of 'caleddy'. Even though a 'minpblh' is
    ! passed into 'exacol', it is not used in the 'exacol'.

    CALL exacol( pcols, pver, ncol, ri, bflxs, minpblh, zi, ktop, kbase, ncvfin )

    ! Diagnostic output of CL interface indices before performing 'extending-merging'
    ! of CL regimes in 'zisocl'
    DO i = 1, ncol
       DO k = 1, ncvmax
          kbase_o(i,k) = REAL(kbase(i,k),r8)
          ktop_o(i,k)  = REAL(ktop(i,k),r8) 
          ncvfin_o(i)  = REAL(ncvfin(i),r8)
       END DO
    END DO

    ! ----------------------------------- !
    ! Perform calculation for each column !
    ! ----------------------------------- !

    DO i = 1, ncol

       ! Define Surface Interfacial Layer TKE, 'tkes'.
       ! In the current code, 'tkes' is used as representing TKE of surface interfacial
       ! layer (low half-layer of surface-based grid layer). In the code, when bflxs>0,
       ! surface interfacial layer is assumed to be energetically  coupled to the other
       ! parts of the CL regime based at the surface. In this sense, it is conceptually
       ! more reasonable to include both 'bprod' and 'sprod' in the definition of 'tkes'.
       ! Since 'tkes' cannot be negative, it is lower bounded by small positive number. 
       ! Note that inclusion of 'bprod' in the definition of 'tkes' may increase 'ebrk'
       ! and 'wstar3', and eventually, 'wet' at the CL top, especially when 'bflxs>0'.
       ! This might help to solve the problem of too shallow PBLH over the overcast Sc
       ! regime. If I want to exclude 'bprod(i,pver+1)' in calculating 'tkes' even when
       ! bflxs > 0, all I should to do is to set 'bprod(i,pver+1) = 0' in the above 
       ! initialization 'do' loop (explained above), NOT changing the formulation of
       ! tkes(i) in the below block. This is because for consistent treatment in the 
       ! other parts of the code also.

       ! tkes(i) = (b1*vk*z(i,pver)*sprod(i,pver+1))**(2._r8/3._r8)
       tkes(i) = MAX(b1*vk*z(i,pver)*(bprod(i,pver+1)+sprod(i,pver+1)), 1.e-7_r8)**(2._r8/3._r8)
       tkes(i) = MIN(tkes(i), tkemax)
       tke(i,pver+1)  = tkes(i)
       wcap(i,pver+1) = tkes(i)/b1

       ! Extend and merge the initially identified CLs, relabel the CLs, and calculate
       ! CL internal mean energetics and stability functions in 'zisocl'. 
       ! The CL nearest to the surface is CL(1) and the CL index, ncv, increases 
       ! with height. The following outputs are from 'zisocl'. Here, the dimension
       ! of below outputs are (pcols,ncvmax) (except the 'ncvfin(pcols)' and 
       ! 'belongcv(pcols,pver+1)) and 'ncv' goes from 1 to 'ncvfin'. 
       ! For 'ncv = ncvfin+1, ncvmax', below output are already initialized to be zero. 
       !      ncvfin       : Total number of CLs
       !      kbase(ncv)   : Base external interface index of CL
       !      ktop         : Top  external interface index of CL
       !      belongcv     : True if the interface (either internal or external) is CL  
       !      ricl         : Mean Richardson number of internal CL
       !      ghcl         : Normalized buoyancy production '-l2n2/2e' [no unit] of internal CL
       !      shcl         : Galperin instability function of heat-moisture of internal CL
       !      smcl         : Galperin instability function of momentum of internal CL
       !      lbrk, <l>int : Thickness of (energetically) internal CL (lint, [m])
       !      wbrk, <W>int : Mean normalized TKE of internal CL  ([m2/s2])
       !      ebrk, <e>int : Mean TKE of internal CL (b1*wbrk,[m2/s2])
       ! The ncvsurf is an identifier saying which CL regime is based at the surface.
       ! If 'ncvsurf=1', then the first CL regime is based at the surface. If surface
       ! interface is not a part of CL (neither internal nor external), 'ncvsurf = 0'.
       ! After identifying and including SRCLs into the normal CL regimes (where newly
       ! identified SRCLs are simply appended to the normal CL regimes using regime 
       ! indices of 'ncvfin+1','ncvfin+2' (as will be shown in the below SRCL part),..
       ! where 'ncvfin' is the final CL regime index produced after extending-merging 
       ! in 'zisocl' but before adding SRCLs), if any newly identified SRCL (e.g., 
       ! 'ncvfin+1') is based at surface, then 'ncvsurf = ncvfin+1'. Thus 'ncvsurf' can
       ! be 0, 1, or >1. 'ncvsurf' can be a useful diagnostic output.   

       ncvsurf = 0
       IF( ncvfin(i) .GT. 0 ) THEN 
          CALL zisocl( pcols  , pver     , i        ,           &
               z      , zi       , n2       , s2      , & 
               bprod  , sprod    , bflxs    , tkes    , landfrac,&
               ncvfin , kbase    , ktop     , belongcv, &
               ricl   , ghcl     , shcl     , smcl    , & 
               lbrk   , wbrk     , ebrk     ,           & 
               extend , extend_up, extend_dn )
          IF( kbase(i,1) .EQ. pver + 1 ) ncvsurf = 1
       ELSE
          belongcv(i,:) = .FALSE.
       ENDIF

       ! Diagnostic output after finishing extending-merging process in 'zisocl'
       ! Since we are adding SRCL additionally, we need to print out these here.

       DO k = 1, ncvmax
          kbase_mg(i,k) = REAL(kbase(i,k))
          ktop_mg(i,k)  = REAL(ktop(i,k)) 
          ncvfin_mg(i)  = REAL(ncvfin(i))
       END DO

       ! ----------------------- !
       ! Identification of SRCLs !
       ! ----------------------- !

       ! Modification : This cannot identify the 'cirrus' layer due to the condition of
       !                ql(i,k) .gt. qmin. This should be modified in future to identify
       !                a single thin cirrus layer.  
       !                Instead of ql, we may use cldn in future, including ice 
       !                contribution.

       ! ------------------------------------------------------------------------------ !
       ! Find single-layer radiatively-driven cloud-topped convective layers (SRCLs).   !
       ! SRCLs extend through a single model layer k, with entrainment at the top and   !
       ! bottom interfaces, unless bottom interface is the surface.                     !
       ! The conditions for an SRCL is identified are:                                  ! 
       !                                                                                !
       !   1. Cloud in the layer, k : ql(i,k) .gt. qmin = 1.e-5 [ kg/kg ]               !
       !   2. No cloud in the above layer (else assuming that some fraction of the LW   !
       !      flux divergence in layer k is concentrated at just below top interface    !
       !      of layer k is invalid). Then, this condition might be sensitive to the    !
       !      vertical resolution of grid.                                              !
       !   3. LW radiative cooling (SW heating is assumed uniformly distributed through !
       !      layer k, so not relevant to buoyancy production) in the layer k. However, !
       !      SW production might also contribute, which may be considered in a future. !
       !   4. Internal stratification 'n2ht' of upper-half layer should be unstable.    !
       !      The 'n2ht' is pure internal stratification of upper half layer, obtained  !
       !      using internal slopes of sl, qt in layer k (in contrast to conventional   !
       !      interfacial slope) and saturation fraction in the upper-half layer,       !
       !      sfuh(k) (in contrast to sfi(k)).                                          !
       !   5. Top and bottom interfaces not both in the same existing convective layer. !
       !      If SRCL is within the previouisly identified CL regimes, we don't define  !
       !      a new SRCL.                                                               !
       !   6. k >= ntop_turb + 1 = 2                                                    !
       !   7. Ri at the top interface > ricrit = 0.19 (otherwise turbulent mixing will  !
       !      broadly distribute the cloud top in the vertical, preventing localized    !
       !      radiative destabilization at the top interface).                          !
       !                                                                                !
       ! Note if 'k = pver', it identifies a surface-based single fog layer, possibly,  !
       ! warm advection fog. Note also the CL regime index of SRCLs itself increases    !
       ! with height similar to the regular CLs indices identified from 'zisocl'.       !
       ! ------------------------------------------------------------------------------ !

       ncv  = 1
       ncvf = ncvfin(i)

       IF( choice_SRCL .EQ. 'remove' ) GOTO 222 

       DO k = nbot_turb, ntop_turb + 1, -1 ! 'k = pver, 2, -1' is a layer index.

          IF( ql(i,k) .GT. qmin(2) .AND. ql(i,k-1) .LT. qmin(2) .AND. qrlw(i,k) .LT. 0._r8 &
               .AND. ri(i,k) .GE. ricrit ) THEN

             ! In order to avoid any confliction with the treatment of ambiguous layer,
             ! I need to impose an additional constraint that ambiguous layer cannot be
             ! SRCL. So, I added constraint that 'k+1' interface (base interface of k
             ! layer) should not be a part of previously identified CL. Since 'belongcv'
             ! is even true for external entrainment interfaces, below constraint is
             ! fully sufficient.

             IF( choice_SRCL .EQ. 'nonamb' .AND. belongcv(i,k+1) ) THEN
                go to 220 
             ENDIF

             ch = ( 1._r8 - sfuh(i,k) ) * chu(i,k) + sfuh(i,k) * chs(i,k)
             cm = ( 1._r8 - sfuh(i,k) ) * cmu(i,k) + sfuh(i,k) * cms(i,k)

             n2htSRCL = ch * slslope(i,k) + cm * qtslope(i,k)

             IF( n2htSRCL .LE. 0._r8 ) THEN

                ! Test if bottom and top interfaces are part of the pre-existing CL. 
                ! If not, find appropriate index for the new SRCL. Note that this
                ! calculation makes use of 'ncv set' obtained from 'zisocl'. The 
                ! 'in_CL' is a parameter testing whether the new SRCL is already 
                ! within the pre-existing CLs (.true.) or not (.false.). 

                in_CL = .FALSE.

                DO WHILE ( ncv .LE. ncvf )
                   IF( ktop(i,ncv) .LE. k ) THEN
                      IF( kbase(i,ncv) .GT. k ) THEN 
                         in_CL = .TRUE.
                      ENDIF
                      EXIT             ! Exit from 'do while' loop if SRCL is within the CLs.
                   ELSE
                      ncv = ncv + 1    ! Go up one CL
                   END IF
                END DO ! ncv

                IF( .NOT. in_CL ) THEN ! SRCL is not within the pre-existing CLs.

                   ! Identify a new SRCL and add it to the pre-existing CL regime group.

                   ncvfin(i)       =  ncvfin(i) + 1
                   ncvnew          =  ncvfin(i)
                   ktop(i,ncvnew)  =  k
                   kbase(i,ncvnew) =  k+1
                   belongcv(i,k)   = .TRUE.
                   belongcv(i,k+1) = .TRUE.

                   ! Calculate internal energy of SRCL. There is no internal energy if
                   ! SRCL is elevated from the surface. Also, we simply assume neutral 
                   ! stability function. Note that this assumption of neutral stability
                   ! does not influence numerical calculation- stability functions here
                   ! are just for diagnostic output. In general SRCLs other than a SRCL 
                   ! based at surface with bflxs <= 0, there is no other way but to use
                   ! neutral stability function.  However, in case of SRCL based at the
                   ! surface,  we can explicitly calculate non-zero stability functions            
                   ! in a consistent way.   Even though stability functions of SRCL are
                   ! just diagnostic outputs not influencing numerical calculations, it
                   ! would be informative to write out correct reasonable values rather
                   ! than simply assuming neutral stability. I am doing this right now.
                   ! Similar calculations were done for the SBCL and when surface inter
                   ! facial layer was merged by overlying CL in 'ziscol'.

                   IF( k .LT. pver ) THEN

                      wbrk(i,ncvnew) = 0._r8
                      ebrk(i,ncvnew) = 0._r8
                      lbrk(i,ncvnew) = 0._r8
                      ghcl(i,ncvnew) = 0._r8
                      shcl(i,ncvnew) = 0._r8
                      smcl(i,ncvnew) = 0._r8
                      ricl(i,ncvnew) = 0._r8

                   ELSE ! Surface-based fog

                      IF( bflxs(i) .GT. 0._r8 ) THEN    ! Incorporate surface TKE into CL interior energy
                         ! It is likely that this case cannot exist  since
                         ! if surface buoyancy flux is positive,  it would
                         ! have been identified as SBCL in 'zisocl' ahead. 
                         ebrk(i,ncvnew) = tkes(i)
                         lbrk(i,ncvnew) = z(i,pver)
                         wbrk(i,ncvnew) = tkes(i) / b1    

                         WRITE(iulog,*) 'Major mistake in SRCL: bflxs > 0 for surface-based SRCL'
                         WRITE(iulog,*) 'bflxs = ', bflxs(i)
                         WRITE(iulog,*) 'ncvfin_o = ', ncvfin_o(i)
                         WRITE(iulog,*) 'ncvfin_mg = ', ncvfin_mg(i)
                         DO ks = 1, ncvmax
                            WRITE(iulog,*) 'ncv =', ks, ' ', kbase_o(i,ks), ktop_o(i,ks), kbase_mg(i,ks), ktop_mg(i,ks)
                         END DO
                         STOP

                      ELSE                              ! Don't incorporate surface interfacial TKE into CL interior energy

                         ebrk(i,ncvnew) = 0._r8
                         lbrk(i,ncvnew) = 0._r8
                         wbrk(i,ncvnew) = 0._r8

                      ENDIF

                      ! Calculate stability functions (ghcl, shcl, smcl, ricl) explicitly
                      ! using an reverse procedure starting from tkes(i). Note that it is
                      ! possible to calculate stability functions even when bflxs < 0.
                      ! Previous code just assumed neutral stability functions. Note that
                      ! since alph5 = 0.7 > 0, alph3 = -35 < 0, the denominator of gh  is
                      ! always positive if bflxs > 0. However, if bflxs < 0,  denominator
                      ! can be zero. For this case, we provide a possible maximum negative
                      ! value (the most stable state) to gh. Note also tkes(i) is always a
                      ! positive value by a limiter. Also, sprod(i,pver+1) > 0 by limiter.

                      gg = 0.5_r8 * vk * z(i,pver) * bprod(i,pver+1) / ( tkes(i)**(3._r8/2._r8) )
                      IF( ABS(alph5-gg*alph3) .LE. 1.e-7_r8 ) THEN
                         ! gh = -0.28_r8
                         ! gh = -3.5334_r8
                         gh = ghmin
                      ELSE    
                         gh = gg / ( alph5 - gg * alph3 )
                      END IF
                      ! gh = min(max(gh,-0.28_r8),0.0233_r8)
                      ! gh = min(max(gh,-3.5334_r8),0.0233_r8)
                      gh = MIN(MAX(gh,ghmin),0.0233_r8)
                      ghcl(i,ncvnew) =  gh
                      shcl(i,ncvnew) =  MAX(0._r8,alph5/(1._r8+alph3*gh))
                      smcl(i,ncvnew) =  MAX(0._r8,(alph1 + alph2*gh)/(1._r8+alph3*gh)/(1._r8+alph4exs*gh))
                      ricl(i,ncvnew) = -(smcl(i,ncvnew)/shcl(i,ncvnew))*(bprod(i,pver+1)/sprod(i,pver+1))

                      ! 'ncvsurf' is CL regime index based at the surface. If there is no
                      ! such regime, then 'ncvsurf = 0'.

                      ncvsurf = ncvnew

                   END IF

                END IF

             END IF

          END IF

220       CONTINUE    

       END DO ! End of 'k' loop where 'k' is a grid layer index running from 'pver' to 2

222    CONTINUE

       ! -------------------------------------------------------------------------- !
       ! Up to this point, we identified all kinds of CL regimes :                  !
       !   1. A SBCL. By construction, 'bflxs > 0' for SBCL.                        !
       !   2. Surface-based CL with multiple layers and 'bflxs =< 0'                !
       !   3. Surface-based CL with multiple layers and 'bflxs > 0'                 !
       !   4. Regular elevated CL with two entraining interfaces                    ! 
       !   5. SRCLs. If SRCL is based at surface, it will be bflxs < 0.             !
       ! '1-4' were identified from 'zisocl' while '5' were identified separately   !
       ! after performing 'zisocl'. CL regime index of '1-4' increases with height  !
       ! ( e.g., CL = 1 is the CL regime nearest to the surface ) while CL regime   !
       ! index of SRCL is simply appended after the final index of CL regimes from  !
       ! 'zisocl'. However, CL regime indices of SRCLs itself increases with height !
       ! when there are multiple SRCLs, similar to the regular CLs from 'zisocl'.   !
       ! -------------------------------------------------------------------------- !

       ! Diagnostic output of final CL regimes indices

       DO k = 1, ncvmax
          kbase_f(i,k) = REAL(kbase(i,k))
          ktop_f(i,k)  = REAL(ktop(i,k)) 
          ncvfin_f(i)  = REAL(ncvfin(i))
       END DO

       ! ---------------------------------------- !
       ! Perform do loop for individual CL regime !
       ! ---------------------------------------- ! -------------------------------- !
       ! For individual CLs, compute                                                 !
       !   1. Entrainment rates at the CL top and (if any) base interfaces using     !
       !      appropriate entrainment closure (current code use 'wstar' closure).    !
       !   2. Net CL mean (i.e., including entrainment contribution) TKE (ebrk)      !
       !      and normalized TKE (wbrk).                                             ! 
       !   3. TKE (tke) and normalized TKE (wcap) profiles at all CL interfaces.     !
       !   4. ( kvm, kvh ) profiles at all CL interfaces.                            !
       !   5. ( bprod, sprod ) profiles at all CL interfaces.                        !
       ! Also calculate                                                              !
       !   1. PBL height as the top external interface of surface-based CL, if any.  !
       !   2. Characteristic excesses of convective 'updraft velocity (wpert)',      !
       !      'temperature (tpert)', and 'moisture (qpert)' in the surface-based CL, !
       !      if any, for use in the separate convection scheme.                     ! 
       ! If there is no surface-based CL, 'PBL height' and 'convective excesses' are !
       ! calculated later from surface-based STL (Stable Turbulent Layer) properties.!
       ! --------------------------------------------------------------------------- !

       ktblw = 0
       DO ncv = 1, ncvfin(i)

          kt = ktop(i,ncv)
          kb = kbase(i,ncv)
          ! Check whether surface interface is energetically interior or not.
          IF( kb .EQ. (pver+1) .AND. bflxs(i) .LE. 0._r8 ) THEN
             lbulk = zi(i,kt) - z(i,pver)
          ELSE
             lbulk = zi(i,kt) - zi(i,kb)
          END IF

          ! Calculate 'turbulent length scale (leng)' and 'normalized TKE (wcap)'
          ! at all CL interfaces except the surface.  Note that below 'wcap' at 
          ! external interfaces are not correct. However, it does not influence 
          ! numerical calculation and correct normalized TKE at the entraining 
          ! interfaces will be re-calculated at the end of this 'do ncv' loop. 

          DO k = MIN(kb,pver), kt, -1 
             IF( choice_tunl .EQ. 'rampcl' ) THEN
                ! In order to treat the case of 'ricl(i,ncv) >> 0' of surface-based SRCL
                ! with 'bflxs(i) < 0._r8', I changed ricl(i,ncv) -> min(0._r8,ricl(i,ncv))
                ! in the below exponential. This is necessary to prevent the model crash
                ! by too large values (e.g., 700) of ricl(i,ncv)   
                tunlramp = ctunl*tunl*(1._r8-(1._r8-1._r8/ctunl)*EXP(MIN(0._r8,ricl(i,ncv))))
                tunlramp = MIN(MAX(tunlramp,tunl),ctunl*tunl)
             ELSEIF( choice_tunl .EQ. 'rampsl' ) THEN
                tunlramp = ctunl*tunl
                ! tunlramp = 0.765_r8
             ELSE
                tunlramp = tunl
             ENDIF
             IF( choice_leng .EQ. 'origin' ) THEN
                leng(i,k) = ( (vk*zi(i,k))**(-cleng) + (tunlramp*lbulk)**(-cleng) )**(-1._r8/cleng)
                ! leng(i,k) = vk*zi(i,k) / (1._r8+vk*zi(i,k)/(tunlramp*lbulk))
             ELSE
                leng(i,k) = MIN( vk*zi(i,k), tunlramp*lbulk )              
             ENDIF
             wcap(i,k) = (leng(i,k)**2) * (-shcl(i,ncv)*n2(i,k)+smcl(i,ncv)*s2(i,k))
          END DO ! k

          ! Calculate basic cross-interface variables ( jump condition ) across the 
          ! base external interface of CL.

          IF( kb .LT. pver+1 ) THEN 

             jbzm = z(i,kb-1) - z(i,kb)                                      ! Interfacial layer thickness [m]
             jbsl = sl(i,kb-1) - sl(i,kb)                                    ! Interfacial jump of 'sl' [J/kg]
             jbqt = qt(i,kb-1) - qt(i,kb)                                    ! Interfacial jump of 'qt' [kg/kg]
             jbbu = n2(i,kb) * jbzm                                          ! Interfacial buoyancy jump [m/s2] considering saturation ( > 0 ) 
             jbbu = MAX(jbbu,jbumin)                                         ! Set minimum buoyancy jump, jbumin = 1.e-3
             jbu  = u(i,kb-1) - u(i,kb)                                      ! Interfacial jump of 'u' [m/s]
             jbv  = v(i,kb-1) - v(i,kb)                                      ! Interfacial jump of 'v' [m/s]
             ch   = (1._r8 -sflh(i,kb-1))*chu(i,kb) + sflh(i,kb-1)*chs(i,kb) ! Buoyancy coefficient just above the base interface
             cm   = (1._r8 -sflh(i,kb-1))*cmu(i,kb) + sflh(i,kb-1)*cms(i,kb) ! Buoyancy coefficient just above the base interface
             n2hb = (ch*jbsl + cm*jbqt)/jbzm                                 ! Buoyancy frequency [s-2] just above the base interface
             vyb  = n2hb*jbzm/jbbu                                           ! Ratio of 'n2hb/n2' at 'kb' interface
             vub  = MIN(1._r8,(jbu**2+jbv**2)/(jbbu*jbzm) )                  ! Ratio of 's2/n2 = 1/Ri' at 'kb' interface

          ELSE 

             ! Below setting is necessary for consistent treatment when 'kb' is at the surface.
             jbbu = 0._r8
             n2hb = 0._r8
             vyb  = 0._r8
             vub  = 0._r8
             web  = 0._r8

          END IF

          ! Calculate basic cross-interface variables ( jump condition ) across the 
          ! top external interface of CL. The meanings of variables are similar to
          ! the ones at the base interface.

          jtzm = z(i,kt-1) - z(i,kt)
          jtsl = sl(i,kt-1) - sl(i,kt)
          jtqt = qt(i,kt-1) - qt(i,kt)
          jtbu = n2(i,kt)*jtzm                                                ! Note : 'jtbu' is guaranteed positive by definition of CL top.
          jtbu = MAX(jtbu,jbumin)                                             ! But threshold it anyway to be sure.
          jtu  = u(i,kt-1) - u(i,kt)
          jtv  = v(i,kt-1) - v(i,kt)
          ch   = (1._r8 -sfuh(i,kt))*chu(i,kt) + sfuh(i,kt)*chs(i,kt) 
          cm   = (1._r8 -sfuh(i,kt))*cmu(i,kt) + sfuh(i,kt)*cms(i,kt) 
          n2ht = (ch*jtsl + cm*jtqt)/jtzm                       
          vyt  = n2ht*jtzm/jtbu                                  
          vut  = MIN(1._r8,(jtu**2+jtv**2)/(jtbu*jtzm))             

          ! Evaporative enhancement factor of entrainment rate at the CL top interface, evhc. 
          ! We take the full inversion strength to be 'jt2slv = slv(i,kt-2)-slv(i,kt)' 
          ! where 'kt-1' is in the ambiguous layer. However, for a cloud-topped CL overlain
          ! by another CL, it is possible that 'slv(i,kt-2) < slv(i,kt)'. To avoid negative
          ! or excessive evhc, we lower-bound jt2slv and upper-bound evhc.  Note 'jtslv' is
          ! used only for calculating 'evhc' : when calculating entrainment rate,   we will
          ! use normal interfacial buoyancy jump across CL top interface.

          evhc   = 1._r8
          jt2slv = 0._r8

          ! Modification : I should check whether below 'jbumin' produces reasonable limiting value.   
          !                In addition, our current formulation does not consider ice contribution. 

          IF( choice_evhc .EQ. 'orig' ) THEN

             IF( ql(i,kt) .GT. qmin(2) .AND. ql(i,kt-1) .LT. qmin(2) ) THEN 
                jt2slv = slv(i,MAX(kt-2,1)) - slv(i,kt)
                jt2slv = MAX( jt2slv, jbumin*slv(i,kt-1)/g )
                evhc   = 1._r8 + a2l * a3l * latvap * ql(i,kt) / jt2slv
                evhc   = MIN( evhc, evhcmax )
             END IF

          ELSEIF( choice_evhc .EQ. 'ramp' ) THEN

             jt2slv = slv(i,MAX(kt-2,1)) - slv(i,kt)
             jt2slv = MAX( jt2slv, jbumin*slv(i,kt-1)/g )
             evhc   = 1._r8 + MAX(cldeff(i,kt)-cldeff(i,kt-1),0._r8) * a2l * a3l * latvap * ql(i,kt) / jt2slv
             evhc   = MIN( evhc, evhcmax )

          ELSEIF( choice_evhc .EQ. 'maxi' ) THEN

             qleff  = MAX( ql(i,kt-1), ql(i,kt) ) 
             jt2slv = slv(i,MAX(kt-2,1)) - slv(i,kt)
             jt2slv = MAX( jt2slv, jbumin*slv(i,kt-1)/g )
             evhc   = 1._r8 + a2l * a3l * latvap * qleff / jt2slv
             evhc   = MIN( evhc, evhcmax )

          ENDIF

          ! Calculate cloud-top radiative cooling contribution to buoyancy production.
          ! Here,  'radf' [m2/s3] is additional buoyancy flux at the CL top interface 
          ! associated with cloud-top LW cooling being mainly concentrated near the CL
          ! top interface ( just below CL top interface ).  Contribution of SW heating
          ! within the cloud is not included in this radiative buoyancy production 
          ! since SW heating is more broadly distributed throughout the CL top layer. 

          lwp        = 0._r8
          opt_depth  = 0._r8
          radinvfrac = 0._r8 
          radf       = 0._r8

          IF( choice_radf .EQ. 'orig' ) THEN

             IF( ql(i,kt) .GT. qmin(2) .AND. ql(i,kt-1) .LT. qmin(2) ) THEN 

                lwp       = ql(i,kt) * ( pi(i,kt+1) - pi(i,kt) ) / g
                opt_depth = 156._r8 * lwp  ! Estimated LW optical depth in the CL top layer

                ! Approximate LW cooling fraction concentrated at the inversion by using
                ! polynomial approx to exact formula 1-2/opt_depth+2/(exp(opt_depth)-1))

                radinvfrac  = opt_depth * ( 4._r8 + opt_depth ) / ( 6._r8 * ( 4._r8 + opt_depth ) + opt_depth**2 )
                radf        = qrlw(i,kt) / ( pi(i,kt) - pi(i,kt+1) ) ! Cp*radiative cooling = [ W/kg ] 
                radf        = MAX( radinvfrac * radf * ( zi(i,kt) - zi(i,kt+1) ), 0._r8 ) * chs(i,kt)
                ! We can disable cloud LW cooling contribution to turbulence by uncommenting:
                ! radf = 0._r8

             END IF

          ELSEIF( choice_radf .EQ. 'ramp' ) THEN

             lwp         = ql(i,kt) * ( pi(i,kt+1) - pi(i,kt) ) / g
             opt_depth   = 156._r8 * lwp  ! Estimated LW optical depth in the CL top layer

             radinvfrac  = opt_depth * ( 4._r8 + opt_depth ) / ( 6._r8 * ( 4._r8 + opt_depth ) + opt_depth**2 )
             radinvfrac  = MAX(cldeff(i,kt)-cldeff(i,kt-1),0._r8) * radinvfrac 
             radf        = qrlw(i,kt) / ( pi(i,kt) - pi(i,kt+1) ) ! Cp*radiative cooling [W/kg] 
             radf        = MAX( radinvfrac * radf * ( zi(i,kt) - zi(i,kt+1) ), 0._r8 ) * chs(i,kt)

          ELSEIF( choice_radf .EQ. 'maxi' ) THEN

             ! Radiative flux divergence both in 'kt' and 'kt-1' layers are included 
             ! 1. From 'kt' layer
             lwp         = ql(i,kt) * ( pi(i,kt+1) - pi(i,kt) ) / g
             opt_depth   = 156._r8 * lwp  ! Estimated LW optical depth in the CL top layer

             radinvfrac  = opt_depth * ( 4._r8 + opt_depth ) / ( 6._r8 * ( 4._r8 + opt_depth ) + opt_depth**2 )
             radf        = MAX( radinvfrac * qrlw(i,kt) / ( pi(i,kt) - pi(i,kt+1) ) * ( zi(i,kt) - zi(i,kt+1) ), 0._r8 )

             ! 2. From 'kt-1' layer and add the contribution from 'kt' layer
             lwp         = ql(i,kt-1) * ( pi(i,kt) - pi(i,kt-1) ) / g
             opt_depth   = 156._r8 * lwp  ! Estimated LW optical depth in the CL top layer

             radinvfrac  = opt_depth * ( 4._r8 + opt_depth ) / ( 6._r8 * ( 4._r8 + opt_depth) + opt_depth**2 )
             radf        = radf + MAX( radinvfrac * qrlw(i,kt-1) / ( pi(i,kt-1) - pi(i,kt) ) * &
                  ( zi(i,kt-1) - zi(i,kt) ), 0.0_r8 )

             radf        = MAX( radf, 0._r8 ) * chs(i,kt) 

          ENDIF

          ! ------------------------------------------------------------------- !
          ! Calculate 'wstar3' by summing buoyancy productions within CL from   !
          !   1. Interior buoyancy production ( bprod: fcn of TKE )             !
          !   2. Cloud-top radiative cooling                                    !
          !   3. Surface buoyancy flux contribution only when bflxs > 0.        !
          !      Note that master length scale, lbulk, has already been         !
          !      corrctly defined at the first part of this 'do ncv' loop       !
          !      considering the sign of bflxs.                                 !
          ! This 'wstar3' is used for calculation of entrainment rate.          !
          ! Note that this 'wstar3' formula does not include shear production   !
          ! and the effect of drizzle, which should be included later.          !
          ! Q : Strictly speaking, in calculating interior buoyancy production, ! 
          !     the use of 'bprod' is not correct, since 'bprod' is not correct !
          !     value but initially guessed value.   More reasonably, we should ! 
          !     use '-leng(i,k)*sqrt(b1*wcap(i,k))*shcl(i,ncv)*n2(i,k)' instead !
          !     of 'bprod(i,k)', although this is still an  approximation since !
          !     tke(i,k) is not exactly 'b1*wcap(i,k)'  due to a transport term.! 
          !     However since iterative calculation will be performed after all,! 
          !     below might also be OK. But I should test this alternative.     !
          ! ------------------------------------------------------------------- !      

          dzht   = zi(i,kt)  - z(i,kt)     ! Thickness of CL top half-layer
          dzhb   = z(i,kb-1) - zi(i,kb)    ! Thickness of CL bot half-layer
          wstar3 = radf * dzht
          DO k = kt + 1, kb - 1 ! If 'kt = kb - 1', this loop will not be performed. 
             wstar3 =  wstar3 + bprod(i,k) * ( z(i,k-1) - z(i,k) )
             ! Below is an alternative which may speed up convergence.
             ! However, for interfaces merged into original CL, it can
             ! be 'wcap(i,k)<0' since 'n2(i,k)>0'.  Thus, I should use
             ! the above original one.
             ! wstar3 =  wstar3 - leng(i,k)*sqrt(b1*wcap(i,k))*shcl(i,ncv)*n2(i,k)* &
             !                    (z(i,k-1) - z(i,k))
          END DO
          IF( kb .EQ. (pver+1) .AND. bflxs(i) .GT. 0._r8 ) THEN
             wstar3 = wstar3 + bflxs(i) * dzhb
             ! wstar3 = wstar3 + bprod(i,pver+1) * dzhb
          END IF
          wstar3 = MAX( 2.5_r8 * wstar3, 0._r8 )

          ! -------------------------------------------------------------- !
          ! Below single block is for 'sedimentation-entrainment feedback' !
          ! -------------------------------------------------------------- !          

          IF( id_sedfact ) THEN
             ! wsed    = 7.8e5_r8*(ql(i,kt)/ncliq(i,kt))**(2._r8/3._r8)
             sedfact = EXP(-ased*wsedl(i,kt)/(wstar3**(1._r8/3._r8)+1.e-6))
             IF( choice_evhc .EQ. 'orig' ) THEN
                IF (ql(i,kt).GT.qmin(2) .AND. ql(i,kt-1).LT.qmin(2)) THEN
                   jt2slv = slv(i,MAX(kt-2,1)) - slv(i,kt)
                   jt2slv = MAX(jt2slv, jbumin*slv(i,kt-1)/g)
                   evhc = 1._r8+sedfact*a2l*a3l*latvap*ql(i,kt) / jt2slv
                   evhc = MIN(evhc,evhcmax)
                END IF
             ELSEIF( choice_evhc .EQ. 'ramp' ) THEN
                jt2slv = slv(i,MAX(kt-2,1)) - slv(i,kt)
                jt2slv = MAX(jt2slv, jbumin*slv(i,kt-1)/g)
                evhc = 1._r8+MAX(cldeff(i,kt)-cldeff(i,kt-1),0._r8)*sedfact*a2l*a3l*latvap*ql(i,kt) / jt2slv
                evhc = MIN(evhc,evhcmax)
             ELSEIF( choice_evhc .EQ. 'maxi' ) THEN
                qleff  = MAX(ql(i,kt-1),ql(i,kt))
                jt2slv = slv(i,MAX(kt-2,1)) - slv(i,kt)
                jt2slv = MAX(jt2slv, jbumin*slv(i,kt-1)/g)
                evhc = 1._r8+sedfact*a2l*a3l*latvap*qleff / jt2slv
                evhc = MIN(evhc,evhcmax)
             ENDIF
          ENDIF

          ! -------------------------------------------------------------------------- !
          ! Now diagnose CL top and bottom entrainment rates (and the contribution of  !
          ! top/bottom entrainments to wstar3) using entrainment closures of the form  !
          !                                                                            !        
          !                   wet = cet*wstar3, web = ceb*wstar3                       !
          !                                                                            !
          ! where cet and ceb depend on the entrainment interface jumps, ql, etc.      !
          ! No entrainment is diagnosed unless the wstar3 > 0. Note '1/wstar3fact' is  !
          ! a factor indicating the enhancement of wstar3 due to entrainment process.  !
          ! Q : Below setting of 'wstar3fact = max(..,0.5)'might prevent the possible  !
          !     case when buoyancy consumption by entrainment is  stronger than cloud  !
          !     top radiative cooling production. Is that OK ? No.  According to bulk  !
          !     modeling study, entrainment buoyancy consumption was always a certain  !
          !     fraction of other net productions, rather than a separate sum.  Thus,  !
          !     below max limit of wstar3fact is correct.   'wstar3fact = max(.,0.5)'  !
          !     prevents unreasonable enhancement of CL entrainment rate by cloud-top  !
          !     entrainment instability, CTEI.                                         !
          ! Q : Use of the same dry entrainment coefficient, 'a1i' both at the CL  top !
          !     and base interfaces may result in too small 'wstar3' and 'ebrk' below, !
          !     as was seen in my generalized bulk modeling study. This should be re-  !
          !     considered later                                                       !
          ! -------------------------------------------------------------------------- !

          IF( wstar3 .GT. 0._r8 ) THEN
             cet = a1i * evhc / ( jtbu * lbulk )
             IF( kb .EQ. pver + 1 ) THEN 
                wstar3fact = MAX( 1._r8 + 2.5_r8 * cet * n2ht * jtzm * dzht, wstar3factcrit )
             ELSE    
                ceb = a1i / ( jbbu * lbulk )
                wstar3fact = MAX( 1._r8 + 2.5_r8 * cet * n2ht * jtzm * dzht &
                     + 2.5_r8 * ceb * n2hb * jbzm * dzhb, wstar3factcrit )
             END IF
             wstar3 = wstar3 / wstar3fact       
          ELSE ! wstar3 == 0
             wstar3fact = 0._r8 ! This is just for dianostic output
             cet        = 0._r8
             ceb        = 0._r8
          END IF

          ! ---------------------------------------------------------------------------- !
          ! Calculate net CL mean TKE including entrainment contribution by solving a    !
          ! canonical cubic equation. The solution of cubic equ. is 'rootp**2 = ebrk'    !
          ! where 'ebrk' originally (before solving cubic eq.) was interior CL mean TKE, !
          ! but after solving cubic equation,  it is replaced by net CL mean TKE in the  !
          ! same variable 'ebrk'.                                                        !
          ! ---------------------------------------------------------------------------- !
          ! Solve cubic equation (canonical form for analytic solution)                  !
          !   r^3 - 3*trmp*r - 2*trmq = 0,   r = sqrt<e>                                 ! 
          ! to estimate <e> for CL, derived from layer-mean TKE balance:                 !
          !                                                                              !
          !   <e>^(3/2)/(b_1*<l>) \approx <B + S>   (*)                                  !
          !   <B+S> = (<B+S>_int * l_int + <B+S>_et * dzt + <B+S>_eb * dzb)/lbulk        !
          !   <B+S>_int = <e>^(1/2)/(b_1*<l>)*<e>_int                                    !
          !   <B+S>_et  = (-vyt+vut)*wet*jtbu + radf                                     !
          !   <B+S>_eb  = (-vyb+vub)*web*jbbu                                            !
          !                                                                              !
          ! where:                                                                       !
          !   <> denotes a vertical avg (over the whole CL unless indicated)             !
          !   l_int (called lbrk below) is aggregate thickness of interior CL layers     !
          !   dzt = zi(i,kt)-z(i,kt)   is thickness of top entrainment layer             !
          !   dzb = z(i,kb-1)-zi(i,kb) is thickness of bot entrainment layer             !
          !   <e>_int (called ebrk below) is the CL-mean TKE if only interior            !
          !                               interfaces contributed.                        !
          !   wet, web                  are top. bottom entrainment rates                !
          !                                                                              !
          ! For a single-level radiatively-driven convective layer, there are no         ! 
          ! interior interfaces so 'ebrk' = 'lbrk' = 0. If the CL goes to the            !
          ! surface, 'vyb' and 'vub' are set to zero before and 'ebrk' and 'lbrk'        !
          ! have already incorporated the surface interfacial layer contribution,        !
          ! so the same formulas still apply.                                            !
          !                                                                              !
          ! In the original formulation based on TKE,                                    !
          !    wet*jtbu = a1l*evhc*<e>^3/2/leng(i,kt)                                    ! 
          !    web*jbbu = a1l*<e>^3/2/leng(i,kt)                                         !
          !                                                                              !
          ! In the wstar formulation                                                     !
          !    wet*jtbu = a1i*evhc*wstar3/lbulk                                          !
          !    web*jbbu = a1i*wstar3/lbulk,                                              !
          ! ---------------------------------------------------------------------------- !

          fact = ( evhc * ( -vyt + vut ) * dzht + ( -vyb + vub ) * dzhb * leng(i,kb) / leng(i,kt) ) / lbulk

          IF( wstarent ) THEN

             ! (Option 1) 'wstar' entrainment formulation 
             ! Here trmq can have either sign, and will usually be nonzero even for non-
             ! cloud topped CLs.  If trmq > 0, there will be two positive roots r; we take 
             ! the larger one. Why ? If necessary, we limit entrainment and wstar to prevent
             ! a solution with r < ccrit*wstar ( Why ? ) where we take ccrit = 0.5. 

             trma = 1._r8          
             trmp = ebrk(i,ncv) * ( lbrk(i,ncv) / lbulk ) / 3._r8 + ntzero
             trmq = 0.5_r8 * b1 * ( leng(i,kt)  / lbulk ) * ( radf * dzht + a1i * fact * wstar3 )

             ! Check if there is an acceptable root with r > rcrit = ccrit*wstar. 
             ! To do this, first find local minimum fmin of the cubic f(r) at sqrt(p), 
             ! and value fcrit = f(rcrit).

             rmin  = SQRT(trmp)
             fmin  = rmin * ( rmin * rmin - 3._r8 * trmp ) - 2._r8 * trmq
             wstar = wstar3**onet
             rcrit = ccrit * wstar
             fcrit = rcrit * ( rcrit * rcrit - 3._r8 * trmp ) - 2._r8 * trmq

             ! No acceptable root exists (noroot = .true.) if either:
             !    1) rmin < rcrit (in which case cubic is monotone increasing for r > rcrit)
             !       and f(rcrit) > 0.
             ! or 2) rmin > rcrit (in which case min of f(r) in r > rcrit is at rmin)
             !       and f(rmin) > 0.  
             ! In this case, we reduce entrainment and wstar3 such that r/wstar = ccrit;
             ! this changes the coefficients of the cubic.   It might be informative to
             ! check when and how many 'noroot' cases occur,  since when 'noroot',   we
             ! will impose arbitrary limit on 'wstar3, wet, web, and ebrk' using ccrit.

             noroot  = ( ( rmin .LT. rcrit ) .AND. ( fcrit .GT. 0._r8 ) ) &
                  .OR. ( ( rmin .GE. rcrit ) .AND. ( fmin  .GT. 0._r8 ) )
             IF( noroot ) THEN ! Solve cubic for r
                trma = 1._r8 - b1 * ( leng(i,kt) / lbulk ) * a1i * fact / ccrit**3
                trma = MAX( trma, 0.5_r8 )  ! Limit entrainment enhancement of ebrk
                trmp = trmp / trma 
                trmq = 0.5_r8 * b1 * ( leng(i,kt) / lbulk ) * radf * dzht / trma
             END IF   ! noroot

             ! Solve the cubic equation

             qq = trmq**2 - trmp**3
             IF( qq .GE. 0._r8 ) THEN 
                rootp = ( trmq + SQRT(qq) )**(1._r8/3._r8) + ( MAX( trmq - SQRT(qq), 0._r8 ) )**(1._r8/3._r8)
             ELSE
                rootp = 2._r8 * SQRT(trmp) * COS( ACOS( trmq / SQRT(trmp**3) ) / 3._r8 )
             END IF

             ! Adjust 'wstar3' only if there is 'noroot'. 
             ! And calculate entrainment rates at the top and base interfaces.

             IF( noroot )  wstar3 = ( rootp / ccrit )**3     ! Adjust wstar3 
             wet = cet * wstar3                              ! Find entrainment rates
             IF( kb .LT. pver + 1 ) web = ceb * wstar3       ! When 'kb.eq.pver+1', it was set to web=0. 

          ELSE !

             ! (Option.2) wstarentr = .false. Use original entrainment formulation.
             ! trmp > 0 if there are interior interfaces in CL, trmp = 0 otherwise.
             ! trmq > 0 if there is cloudtop radiative cooling, trmq = 0 otherwise.

             trma = 1._r8 - b1 * a1l * fact
             trma = MAX( trma, 0.5_r8 )  ! Prevents runaway entrainment instability
             trmp = ebrk(i,ncv) * ( lbrk(i,ncv) / lbulk ) / ( 3._r8 * trma )
             trmq = 0.5_r8 * b1 * ( leng(i,kt)  / lbulk ) * radf * dzht / trma

             qq = trmq**2 - trmp**3
             IF( qq .GE. 0._r8 ) THEN 
                rootp = ( trmq + SQRT(qq) )**(1._r8/3._r8) + ( MAX( trmq - SQRT(qq), 0._r8 ) )**(1._r8/3._r8)
             ELSE ! Also part of case 3
                rootp = 2._r8 * SQRT(trmp) * COS( ACOS( trmq / SQRT(trmp**3) ) / 3._r8 )
             END IF   ! qq

             ! Find entrainment rates and limit them by free-entrainment values a1l*sqrt(e)

             wet = a1l * rootp * MIN( evhc * rootp**2 / ( leng(i,kt) * jtbu ), 1._r8 )   
             IF( kb .LT. pver + 1 ) web = a1l * rootp * MIN( evhc * rootp**2 / ( leng(i,kb) * jbbu ), 1._r8 )

          END IF ! wstarentr

          ! ---------------------------------------------------- !
          ! Finally, get the net CL mean TKE and normalized TKE  ! 
          ! ---------------------------------------------------- !

          ebrk(i,ncv) = MAX(MIN(rootp,tkemax),tkemin)**2
          ebrk(i,ncv) = MIN(ebrk(i,ncv),tkemax) ! Limit CL-avg TKE used for entrainment
          ebrk(i,ncv) = MAX(ebrk(i,ncv),tkemin) ! Limit CL-avg TKE used for entrainment
          wbrk(i,ncv) = ebrk(i,ncv)/b1  

          ! The only way ebrk = 0 is for SRCL which are actually radiatively cooled 
          ! at top interface. In this case, we remove 'convective' label from the 
          ! interfaces around this layer. This case should now be impossible, so 
          ! we flag it. Q: I can't understand why this case is impossible now. Maybe,
          ! due to various limiting procedures used in solving cubic equation ? 
          ! In case of SRCL, 'ebrk' should be positive due to cloud top LW radiative
          ! cooling contribution, although 'ebrk(internal)' of SRCL before including
          ! entrainment contribution (which include LW cooling contribution also) is
          ! zero. 

          IF( ebrk(i,ncv) .LE. 0._r8 ) THEN
             WRITE(iulog,*) 'CALEDDY: Warning, CL with zero TKE, i, kt, kb ', i, kt, kb
             belongcv(i,kt) = .FALSE.
             belongcv(i,kb) = .FALSE. 
          END IF

          ! ----------------------------------------------------------------------- !
          ! Calculate complete TKE profiles at all CL interfaces, capped by tkemax. !
          ! We approximate TKE = <e> at entrainment interfaces. However when CL is  !
          ! based at surface, correct 'tkes' will be inserted to tke(i,pver+1).     !
          ! Note that this approximation at CL external interfaces do not influence !
          ! numerical calculation since 'e' at external interfaces are not used  in !
          ! actual numerical calculation afterward. In addition in order to extract !
          ! correct TKE averaged over the PBL in the cumulus scheme,it is necessary !
          ! to set e = <e> at the top entrainment interface.  Since net CL mean TKE !
          ! 'ebrk' obtained by solving cubic equation already includes tkes  ( tkes !
          ! is included when bflxs > 0 but not when bflxs <= 0 into internal ebrk ),!
          ! 'tkes' should be written to tke(i,pver+1)                               !
          ! ----------------------------------------------------------------------- !

          ! 1. At internal interfaces          
          DO k = kb - 1, kt + 1, -1
             rcap = ( b1 * ae + wcap(i,k) / wbrk(i,ncv) ) / ( b1 * ae + 1._r8 )
             rcap = MIN( MAX(rcap,rcapmin), rcapmax )
             tke(i,k) = ebrk(i,ncv) * rcap
             tke(i,k) = MIN( tke(i,k), tkemax )
             kvh(i,k) = leng(i,k) * SQRT(tke(i,k)) * shcl(i,ncv)
             kvm(i,k) = leng(i,k) * SQRT(tke(i,k)) * smcl(i,ncv)
             bprod(i,k) = -kvh(i,k) * n2(i,k)
             sprod(i,k) =  kvm(i,k) * s2(i,k)
             turbtype(i,k) = 2                     ! CL interior interfaces.
             sm_aw(i,k) = smcl(i,ncv)/alph1        ! Diagnostic output for microphysics
          END DO

          ! 2. At CL top entrainment interface
          kentr = wet * jtzm
          kvh(i,kt) = kentr
          kvm(i,kt) = kentr
          bprod(i,kt) = -kentr * n2ht + radf       ! I must use 'n2ht' not 'n2'
          sprod(i,kt) =  kentr * s2(i,kt)
          turbtype(i,kt) = 4                       ! CL top entrainment interface
          trmp = -b1 * ae / ( 1._r8 + b1 * ae )
          trmq = -(bprod(i,kt)+sprod(i,kt))*b1*leng(i,kt)/(1._r8+b1*ae)/(ebrk(i,ncv)**(3._r8/2._r8))
          rcap = compute_cubic(0._r8,trmp,trmq)**2._r8
          rcap = MIN( MAX(rcap,rcapmin), rcapmax )
          tke(i,kt)  = ebrk(i,ncv) * rcap
          tke(i,kt)  = MIN( tke(i,kt), tkemax )
          sm_aw(i,kt) = smcl(i,ncv) / alph1        ! Diagnostic output for microphysics

          ! 3. At CL base entrainment interface and double entraining interfaces
          ! When current CL base is also the top interface of CL regime below,
          ! simply add the two contributions for calculating eddy diffusivity
          ! and buoyancy/shear production. Below code correctly works because
          ! we (CL regime index) always go from surface upward.

          IF( kb .LT. pver + 1 ) THEN 

             kentr = web * jbzm

             IF( kb .NE. ktblw ) THEN

                kvh(i,kb) = kentr
                kvm(i,kb) = kentr
                bprod(i,kb) = -kvh(i,kb)*n2hb     ! I must use 'n2hb' not 'n2'
                sprod(i,kb) =  kvm(i,kb)*s2(i,kb)
                turbtype(i,kb) = 3                ! CL base entrainment interface
                trmp = -b1*ae/(1._r8+b1*ae)
                trmq = -(bprod(i,kb)+sprod(i,kb))*b1*leng(i,kb)/(1._r8+b1*ae)/(ebrk(i,ncv)**(3._r8/2._r8))
                rcap = compute_cubic(0._r8,trmp,trmq)**2._r8
                rcap = MIN( MAX(rcap,rcapmin), rcapmax )
                tke(i,kb)  = ebrk(i,ncv) * rcap
                tke(i,kb)  = MIN( tke(i,kb),tkemax )

             ELSE

                kvh(i,kb) = kvh(i,kb) + kentr 
                kvm(i,kb) = kvm(i,kb) + kentr
                ! dzhb5 : Half thickness of the lowest  layer of  current CL regime
                ! dzht5 : Half thickness of the highest layer of adjacent CL regime just below current CL. 
                dzhb5 = z(i,kb-1) - zi(i,kb)
                dzht5 = zi(i,kb) - z(i,kb)
                bprod(i,kb) = ( dzht5*bprod(i,kb) - dzhb5*kentr*n2hb )     / ( dzhb5 + dzht5 )
                sprod(i,kb) = ( dzht5*sprod(i,kb) + dzhb5*kentr*s2(i,kb) ) / ( dzhb5 + dzht5 )
                trmp = -b1*ae/(1._r8+b1*ae)
                trmq = -kentr*(s2(i,kb)-n2hb)*b1*leng(i,kb)/(1._r8+b1*ae)/(ebrk(i,ncv)**(3._r8/2._r8))
                rcap = compute_cubic(0._r8,trmp,trmq)**2._r8
                rcap = MIN( MAX(rcap,rcapmin), rcapmax )
                tke_imsi = ebrk(i,ncv) * rcap
                tke_imsi = MIN( tke_imsi, tkemax )
                tke(i,kb)  = ( dzht5*tke(i,kb) + dzhb5*tke_imsi ) / ( dzhb5 + dzht5 )               
                tke(i,kb)  = MIN(tke(i,kb),tkemax)
                turbtype(i,kb) = 5                ! CL double entraining interface      

             END IF

          ELSE

             ! If CL base interface is surface, compute similarly using wcap(i,kb)=tkes/b1    
             ! Even when bflx < 0, use the same formula in order to impose consistency of
             ! tke(i,kb) at bflx = 0._r8

             rcap = (b1*ae + wcap(i,kb)/wbrk(i,ncv))/(b1*ae + 1._r8)
             rcap = MIN( MAX(rcap,rcapmin), rcapmax )
             tke(i,kb) = ebrk(i,ncv) * rcap
             tke(i,kb) = MIN( tke(i,kb),tkemax )

          END IF

          ! For double entraining interface, simply use smcl(i,ncv) of the overlying CL. 
          ! Below 'sm_aw' is a diagnostic output for use in the microphysics.
          ! When 'kb' is surface, 'sm' will be over-written later below.

          sm_aw(i,kb) = smcl(i,ncv)/alph1             

          ! Calculate wcap at all interfaces of CL. Put a  minimum threshold on TKE
          ! to prevent possible division by zero.  'wcap' at CL internal interfaces
          ! are already calculated in the first part of 'do ncv' loop correctly.
          ! When 'kb.eq.pver+1', below formula produces the identical result to the
          ! 'tkes(i)/b1' if leng(i,kb) is set to vk*z(i,pver). Note  wcap(i,pver+1)
          ! is already defined as 'tkes(i)/b1' at the first part of caleddy.

          wcap(i,kt) = (bprod(i,kt)+sprod(i,kt))*leng(i,kt)/SQRT(MAX(tke(i,kt),1.e-6_r8))
          IF( kb .LT. pver + 1 ) THEN
             wcap(i,kb) = (bprod(i,kb)+sprod(i,kb))*leng(i,kb)/SQRT(MAX(tke(i,kb),1.e-6_r8))
          END IF

          ! Save the index of upper external interface of current CL-regime in order to
          ! handle the case when this interface is also the lower external interface of 
          ! CL-regime located just above. 

          ktblw = kt 

          ! Diagnostic Output

          wet_CL(i,ncv)        = wet
          web_CL(i,ncv)        = web
          jtbu_CL(i,ncv)       = jtbu
          jbbu_CL(i,ncv)       = jbbu
          evhc_CL(i,ncv)       = evhc
          jt2slv_CL(i,ncv)     = jt2slv
          n2ht_CL(i,ncv)       = n2ht
          n2hb_CL(i,ncv)       = n2hb          
          lwp_CL(i,ncv)        = lwp
          opt_depth_CL(i,ncv)  = opt_depth
          radinvfrac_CL(i,ncv) = radinvfrac
          radf_CL(i,ncv)       = radf
          wstar_CL(i,ncv)      = wstar          
          wstar3fact_CL(i,ncv) = wstar3fact          

       END DO        ! ncv

       ! Calculate PBL height and characteristic cumulus excess for use in the
       ! cumulus convection shceme. Also define turbulence type at the surface
       ! when the lowest CL is based at the surface. These are just diagnostic
       ! outputs, not influencing numerical calculation of current PBL scheme.
       ! If the lowest CL is based at the surface, define the PBL depth as the
       ! CL top interface. The same rule is applied for all CLs including SRCL.

       IF( ncvsurf .GT. 0 ) THEN

          ktopbl(i) = ktop(i,ncvsurf)
          pblh(i)   = MAX(zi(i, ktopbl(i))-zi(i,pver+1),1.0_r8)
          pblhp(i)  = pi(i, ktopbl(i))
          wpert(i)  = MAX(wfac*SQRT(ebrk(i,ncvsurf)),wpertmin)
          tpert(i)  = MAX(ABS(shflx(i)*rrho(i)/cpair)*tfac/wpert(i),0._r8)
          qpert(i)  = MAX(ABS(qflx(i)*rrho(i))*tfac/wpert(i),0._r8)

          IF( bflxs(i) .GT. 0._r8 ) THEN
             turbtype(i,pver+1) = 2 ! CL interior interface
          ELSE
             turbtype(i,pver+1) = 3 ! CL external base interface
          ENDIF

          ipbl(i)  = 1._r8
          kpblh(i) = ktopbl(i) - 1._r8

       END IF ! End of the calculationf of te properties of surface-based CL.

       ! -------------------------------------------- !
       ! Treatment of Stable Turbulent Regime ( STL ) !
       ! -------------------------------------------- !

       ! Identify top and bottom most (internal) interfaces of STL except surface.
       ! Also, calculate 'turbulent length scale (leng)' at each STL interfaces.     

       belongst(i,1) = .FALSE.   ! k = 1 (top interface) is assumed non-turbulent
       DO k = 2, pver            ! k is an interface index
          belongst(i,k) = ( ri(i,k) .LT. ricrit ) .AND. ( .NOT. belongcv(i,k) )
          IF( belongst(i,k) .AND. ( .NOT. belongst(i,k-1) ) ) THEN
             kt = k             ! Top interface index of STL
          ELSEIF( .NOT. belongst(i,k) .AND. belongst(i,k-1) ) THEN
             kb = k - 1         ! Base interface index of STL
             lbulk = z(i,kt-1) - z(i,kb)
             DO ks = kt, kb
                IF( choice_tunl .EQ. 'rampcl' ) THEN
                   tunlramp = tunl
                ELSEIF( choice_tunl .EQ. 'rampsl' ) THEN
                   !PRINT*,'paulo kubota',ri(i,ks),tunl
                   tunlramp = MAX( 1.e-3_r8, ctunl * tunl * EXP(-LOG(ctunl)*MIN(MAX(ri(i,ks),-10.0_r8),10.0_r8)/ricrit) )
                   ! tunlramp = 0.065_r8 + 0.7_r8 * exp(-20._r8*ri(i,ks))
                ELSE
                   tunlramp = tunl
                ENDIF
                IF( choice_leng .EQ. 'origin' ) THEN
                   leng(i,ks) = ( (vk*zi(i,ks))**(-cleng) + (tunlramp*lbulk)**(-cleng) )**(-1._r8/cleng)
                   ! leng(i,ks) = vk*zi(i,ks) / (1._r8+vk*zi(i,ks)/(tunlramp*lbulk))
                ELSE
                   leng(i,ks) = MIN( vk*zi(i,ks), tunlramp*lbulk )              
                ENDIF
             END DO
          END IF
       END DO ! k

       ! Now look whether STL extends to ground.  If STL extends to surface,
       ! re-define master length scale,'lbulk' including surface interfacial
       ! layer thickness, and re-calculate turbulent length scale, 'leng' at
       ! all STL interfaces again. Note that surface interface is assumed to
       ! always be STL if it is not CL.   

       belongst(i,pver+1) = .NOT. belongcv(i,pver+1)

       IF( belongst(i,pver+1) ) THEN     ! kb = pver+1 (surface  STL)

          turbtype(i,pver+1) = 1        ! Surface is STL interface

          IF( belongst(i,pver) ) THEN   ! STL includes interior
             ! 'kt' already defined above as the top interface of STL
             lbulk = z(i,kt-1)          
          ELSE                          ! STL with no interior turbulence
             kt = pver+1
             lbulk = z(i,kt-1)
          END IF

          ! PBL height : Layer mid-point just above the highest STL interface
          ! Note in contrast to the surface based CL regime where  PBL height
          ! was defined at the top external interface, PBL height of  surface
          ! based STL is defined as the layer mid-point.

          ktopbl(i) = kt - 1
          pblh(i)   = MAX(z(i,ktopbl(i))-z(i,pver),1.0_r8)
          pblhp(i)  = 0.5_r8 * ( pi(i,ktopbl(i)) + pi(i,ktopbl(i)+1) )          

          ! Re-calculate turbulent length scale including surface interfacial
          ! layer contribution to lbulk.

          DO ks = kt, pver
             IF( choice_tunl .EQ. 'rampcl' ) THEN
                tunlramp = tunl
             ELSEIF( choice_tunl .EQ. 'rampsl' ) THEN
                tunlramp = MAX(1.e-3_r8,ctunl*tunl*EXP(-LOG(ctunl)*ri(i,ks)/ricrit))
                ! tunlramp = 0.065_r8 + 0.7_r8 * exp(-20._r8*ri(i,ks))
             ELSE
                tunlramp = tunl
             ENDIF
             IF( choice_leng .EQ. 'origin' ) THEN
                leng(i,ks) = ( (vk*zi(i,ks))**(-cleng) + (tunlramp*lbulk)**(-cleng) )**(-1._r8/cleng)
                ! leng(i,ks) = vk*zi(i,ks) / (1._r8+vk*zi(i,ks)/(tunlramp*lbulk))
             ELSE
                leng(i,ks) = MIN( vk*zi(i,ks), tunlramp*lbulk )              
             ENDIF
          END DO ! ks

          ! Characteristic cumulus excess of surface-based STL.
          ! We may be able to use ustar for wpert.

          wpert(i) = 0._r8 
          tpert(i) = MAX(shflx(i)*rrho(i)/cpair*fak/ustar(i),0._r8) ! CCM stable-layer forms
          qpert(i) = MAX(qflx(i)*rrho(i)*fak/ustar(i),0._r8)

          ipbl(i)  = 0._r8
          kpblh(i) = ktopbl(i)

       END IF

       ! Calculate stability functions and energetics at the STL interfaces
       ! except the surface. Note that tke(i,pver+1) and wcap(i,pver+1) are
       ! already calculated in the first part of 'caleddy', kvm(i,pver+1) &
       ! kvh(i,pver+1) were already initialized to be zero, bprod(i,pver+1)
       ! & sprod(i,pver+1) were direcly calculated from the bflxs and ustar.
       ! Note transport term is assumed to be negligible at STL interfaces.

       DO k = 2, pver

          IF( belongst(i,k) ) THEN

             turbtype(i,k) = 1    ! STL interfaces
             trma = alph3*alph4exs*ri(i,k) + 2._r8*b1*(alph2-alph4exs*alph5*ri(i,k))
             trmb = (alph3+alph4exs)*ri(i,k) + 2._r8*b1*(-alph5*ri(i,k)+alph1)
             trmc = ri(i,k)
             det = MAX(trmb*trmb-4._r8*trma*trmc,0._r8)
             ! Sanity Check
             IF( det .LT. 0._r8 ) THEN
                WRITE(iulog,*) 'The det < 0. for the STL in UW eddy_diff'             
                STOP
             END IF
             gh = (-trmb + SQRT(det))/(2._r8*trma)
             ! gh = min(max(gh,-0.28_r8),0.0233_r8)
             ! gh = min(max(gh,-3.5334_r8),0.0233_r8)
             gh = MIN(MAX(gh,ghmin),0.0233_r8)
             sh = MAX(0._r8,alph5/(1._r8+alph3*gh))
             sm = MAX(0._r8,(alph1 + alph2*gh)/(1._r8+alph3*gh)/(1._r8+alph4exs*gh))

             tke(i,k)   = b1*(leng(i,k)**2)*(-sh*n2(i,k)+sm*s2(i,k))
             tke(i,k)   = MIN(tke(i,k),tkemax)
             wcap(i,k)  = tke(i,k)/b1
             kvh(i,k)   = leng(i,k) * SQRT(tke(i,k)) * sh
             kvm(i,k)   = leng(i,k) * SQRT(tke(i,k)) * sm
             bprod(i,k) = -kvh(i,k) * n2(i,k)
             sprod(i,k) =  kvm(i,k) * s2(i,k)

             sm_aw(i,k) = sm/alph1     ! This is diagnostic output for use in the microphysics             

          END IF

       END DO  ! k

       ! --------------------------------------------------- !
       ! End of treatment of Stable Turbulent Regime ( STL ) !
       ! --------------------------------------------------- !

       ! --------------------------------------------------------------- !
       ! Re-computation of eddy diffusivity at the entrainment interface !
       ! assuming that it is purely STL (0<Ri<0.19). Note even Ri>0.19,  !
       ! turbulent can exist at the entrainment interface since 'Sh,Sm'  !
       ! do not necessarily go to zero even when Ri>0.19. Since Ri can   !
       ! be fairly larger than 0.19 at the entrainment interface, I      !
       ! should set minimum value of 'tke' to be 0. in order to prevent  !
       ! sqrt(tke) from being imaginary.                                 !
       ! --------------------------------------------------------------- !

       ! goto 888

       DO k = 2, pver

          IF( ( turbtype(i,k) .EQ. 3 ) .OR. ( turbtype(i,k) .EQ. 4 ) .OR. &
               ( turbtype(i,k) .EQ. 5 ) ) THEN

             trma = alph3*alph4exs*ri(i,k) + 2._r8*b1*(alph2-alph4exs*alph5*ri(i,k))
             trmb = (alph3+alph4exs)*ri(i,k) + 2._r8*b1*(-alph5*ri(i,k)+alph1)
             trmc = ri(i,k)
             det  = MAX(trmb*trmb-4._r8*trma*trmc,0._r8)
             gh   = (-trmb + SQRT(det))/(2._r8*trma)
             ! gh   = min(max(gh,-0.28_r8),0.0233_r8)
             ! gh   = min(max(gh,-3.5334_r8),0.0233_r8)
             gh   = MIN(MAX(gh,ghmin),0.0233_r8)
             sh   = MAX(0._r8,alph5/(1._r8+alph3*gh))
             sm   = MAX(0._r8,(alph1 + alph2*gh)/(1._r8+alph3*gh)/(1._r8+alph4exs*gh))

             lbulk = z(i,k-1) - z(i,k)

             IF( choice_tunl .EQ. 'rampcl' ) THEN
                tunlramp = tunl
             ELSEIF( choice_tunl .EQ. 'rampsl' ) THEN
               !tunlramp = MAX(1.e-3_r8,ctunl*tunl*EXP(-LOG(ctunl)*MIN(MAX(ri(i,k),-10.0_r8),10.0_r8)/ricrit))
                tunlramp = MAX(1.e-3_r8,ctunl*tunl*EXP(-LOG(ctunl)*MIN(MAX(ri(i,k),-10.0_r8),10.0_r8)/ricrit))

                ! tunlramp = 0.065_r8 + 0.7_r8*exp(-20._r8*ri(i,k))
             ELSE
                tunlramp = tunl
             ENDIF
             IF( choice_leng .EQ. 'origin' ) THEN
                leng_imsi = ( (vk*zi(i,k))**(-cleng) + (tunlramp*lbulk)**(-cleng) )**(-1._r8/cleng)
                ! leng_imsi = vk*zi(i,k) / (1._r8+vk*zi(i,k)/(tunlramp*lbulk))
             ELSE
                leng_imsi = MIN( vk*zi(i,k), tunlramp*lbulk )              
             ENDIF

             tke_imsi = b1*(leng_imsi**2)*(-sh*n2(i,k)+sm*s2(i,k))
             tke_imsi = MIN(MAX(tke_imsi,0._r8),tkemax)
             kvh_imsi = leng_imsi * SQRT(tke_imsi) * sh
             kvm_imsi = leng_imsi * SQRT(tke_imsi) * sm

             IF( kvh(i,k) .LT. kvh_imsi ) THEN 
                kvh(i,k)   =  kvh_imsi
                kvm(i,k)   =  kvm_imsi
                leng(i,k)  = leng_imsi
                tke(i,k)   =  tke_imsi
                wcap(i,k)  =  tke_imsi / b1
                bprod(i,k) = -kvh_imsi * n2(i,k)
                sprod(i,k) =  kvm_imsi * s2(i,k)
                sm_aw(i,k) =  sm/alph1     ! This is diagnostic output for use in the microphysics             
                turbtype(i,k) = 1          ! This was added on Dec.10.2009 for use in microphysics.
             ENDIF

          END IF

       END DO

       ! 888   continue 

       ! ------------------------------------------------------------------ !
       ! End of recomputation of eddy diffusivity at entrainment interfaces !
       ! ------------------------------------------------------------------ !

       ! As an option, we can impose a certain minimum back-ground diffusivity.

       ! do k = 1, pver+1
       !    kvh(i,k) = max(0.01_r8,kvh(i,k))
       !    kvm(i,k) = max(0.01_r8,kvm(i,k))
       ! enddo

       ! --------------------------------------------------------------------- !
       ! Diagnostic Output                                                     !
       ! Just for diagnostic purpose, calculate stability functions at  each   !
       ! interface including surface. Instead of assuming neutral stability,   !
       ! explicitly calculate stability functions using an reverse procedure   !
       ! starting from tkes(i) similar to the case of SRCL and SBCL in zisocl. !
       ! Note that it is possible to calculate stability functions even when   !
       ! bflxs < 0. Note that this inverse method allows us to define Ri even  !
       ! at the surface. Note also tkes(i) and sprod(i,pver+1) are always      !
       ! positive values by limiters (e.g., ustar_min = 0.01).                 !
       ! Dec.12.2006 : Also just for diagnostic output, re-set                 !
       ! 'bprod(i,pver+1)= bflxs(i)' here. Note that this setting does not     !
       ! influence numerical calculation at all - it is just for diagnostic    !
       ! output.                                                               !
       ! --------------------------------------------------------------------- !

       bprod(i,pver+1) = bflxs(i)

       gg = 0.5_r8*vk*z(i,pver)*bprod(i,pver+1)/(tkes(i)**(3._r8/2._r8))
       IF( ABS(alph5-gg*alph3) .LE. 1.e-7_r8 ) THEN
          ! gh = -0.28_r8
          IF( bprod(i,pver+1) .GT. 0._r8 ) THEN
             gh = -3.5334_r8
          ELSE
             gh = ghmin
          ENDIF
       ELSE    
          gh = gg/(alph5-gg*alph3)
       END IF

       ! gh = min(max(gh,-0.28_r8),0.0233_r8)
       IF( bprod(i,pver+1) .GT. 0._r8 ) THEN
          gh = MIN(MAX(gh,-3.5334_r8),0.0233_r8)
       ELSE
          gh = MIN(MAX(gh,ghmin),0.0233_r8)
       ENDIF

       gh_a(i,pver+1) = gh     
       sh_a(i,pver+1) = MAX(0._r8,alph5/(1._r8+alph3*gh))
       IF( bprod(i,pver+1) .GT. 0._r8 ) THEN       
          sm_a(i,pver+1) = MAX(0._r8,(alph1+alph2*gh)/(1._r8+alph3*gh)/(1._r8+alph4*gh))
       ELSE
          sm_a(i,pver+1) = MAX(0._r8,(alph1+alph2*gh)/(1._r8+alph3*gh)/(1._r8+alph4exs*gh))
       ENDIF
       sm_aw(i,pver+1) = sm_a(i,pver+1)/alph1
       ri_a(i,pver+1)  = -(sm_a(i,pver+1)/sh_a(i,pver+1))*(bprod(i,pver+1)/sprod(i,pver+1))

       DO k = 1, pver
          IF( ri(i,k) .LT. 0._r8 ) THEN
             trma = alph3*alph4*ri(i,k) + 2._r8*b1*(alph2-alph4*alph5*ri(i,k))
             trmb = (alph3+alph4)*ri(i,k) + 2._r8*b1*(-alph5*ri(i,k)+alph1)
             trmc = ri(i,k)
             det  = MAX(trmb*trmb-4._r8*trma*trmc,0._r8)
             gh   = (-trmb + SQRT(det))/(2._r8*trma)
             gh   = MIN(MAX(gh,-3.5334_r8),0.0233_r8)
             gh_a(i,k) = gh
             sh_a(i,k) = MAX(0._r8,alph5/(1._r8+alph3*gh))
             sm_a(i,k) = MAX(0._r8,(alph1+alph2*gh)/(1._r8+alph3*gh)/(1._r8+alph4*gh))
             ri_a(i,k) = ri(i,k)
          ELSE
             IF( ri(i,k) .GT. ricrit ) THEN
                gh_a(i,k) = ghmin
                sh_a(i,k) = 0._r8
                sm_a(i,k) = 0._r8
                ri_a(i,k) = ri(i,k)
             ELSE
                trma = alph3*alph4exs*ri(i,k) + 2._r8*b1*(alph2-alph4exs*alph5*ri(i,k))
                trmb = (alph3+alph4exs)*ri(i,k) + 2._r8*b1*(-alph5*ri(i,k)+alph1)
                trmc = ri(i,k)
                det  = MAX(trmb*trmb-4._r8*trma*trmc,0._r8)
                gh   = (-trmb + SQRT(det))/(2._r8*trma)
                gh   = MIN(MAX(gh,ghmin),0.0233_r8)
                gh_a(i,k) = gh
                sh_a(i,k) = MAX(0._r8,alph5/(1._r8+alph3*gh))
                sm_a(i,k) = MAX(0._r8,(alph1+alph2*gh)/(1._r8+alph3*gh)/(1._r8+alph4exs*gh))
                ri_a(i,k) = ri(i,k)
             ENDIF
          ENDIF

       END DO

       DO k = 1, pver + 1
          turbtype_f(i,k) = REAL(turbtype(i,k))
       END DO

    END DO   ! End of column index loop, i 

    RETURN

  END SUBROUTINE caleddy



  !
  !  exacol
  !

  !============================================================================== !
  !                                                                               !
  !============================================================================== !

  SUBROUTINE exacol( pcols, pver, ncol, ri, bflxs, minpblh, zi, ktop, kbase, ncvfin ) 

    ! ---------------------------------------------------------------------------- !
    ! Object : Find unstable CL regimes and determine the indices                  !
    !          kbase, ktop which delimit these unstable layers :                   !
    !          ri(kbase) > 0 and ri(ktop) > 0, but ri(k) < 0 for ktop < k < kbase. ! 
    ! Author : Chris  Bretherton 08/2000,                                          !
    !          Sungsu Park       08/2006, 11/2008                                  !
    !----------------------------------------------------------------------------- !

    IMPLICIT NONE

    ! --------------- !
    ! Input variables !
    ! --------------- !

    INTEGER,  INTENT(in) :: pcols                  ! Number of atmospheric columns   
    INTEGER,  INTENT(in) :: pver                   ! Number of atmospheric vertical layers   
    INTEGER,  INTENT(in) :: ncol                   ! Number of atmospheric columns   

    REAL(r8), INTENT(in) :: ri(pcols,pver)         ! Moist gradient Richardson no.
    REAL(r8), INTENT(in) :: bflxs(pcols)           ! Buoyancy flux at surface
    REAL(r8), INTENT(in) :: minpblh(pcols)         ! Minimum PBL height based on surface stress
    REAL(r8), INTENT(in) :: zi(pcols,pver+1)       ! Interface heights

    ! ---------------- !
    ! Output variables !      
    ! ---------------- !

    INTEGER, INTENT(out) :: kbase(pcols,ncvmax)    ! External interface index of CL base
    INTEGER, INTENT(out) :: ktop(pcols,ncvmax)     ! External interface index of CL top
    INTEGER, INTENT(out) :: ncvfin(pcols)          ! Total number of CLs

    ! --------------- !
    ! Local variables !
    ! --------------- !

    INTEGER              :: i
    INTEGER              :: k
    INTEGER              :: ncv
    REAL(r8)             :: rimaxentr
    REAL(r8)             :: riex(pver+1)           ! Column Ri profile extended to surface

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !
    rimaxentr=0.0_r8
    riex(1:pver+1)=0.0_r8
    DO i = 1, ncol
       ncvfin(i) = 0
    END DO

    DO ncv = 1, ncvmax
       DO i = 1, ncol
          ktop(i,ncv)  = 0
          kbase(i,ncv) = 0
       END DO
    END DO

    ! ------------------------------------------------------ !
    ! Find CL regimes starting from the surface going upward !
    ! ------------------------------------------------------ !

    rimaxentr = 0._r8   

    DO i = 1, ncol

       riex(2:pver) = ri(i,2:pver)

       ! Below allows consistent treatment of surface and other interfaces.
       ! Simply, if surface buoyancy flux is positive, Ri of surface is set to be negative.

       riex(pver+1) = rimaxentr - bflxs(i) 

       ncv = 0
       k   = pver + 1 ! Work upward from surface interface

       DO WHILE ( k .GT. ntop_turb + 2 )

          ! Below means that if 'bflxs > 0' (do not contain '=' sign), surface
          ! interface is energetically interior surface. 

          IF( riex(k) .LT. rimaxentr ) THEN 

             ! Identify a new CL

             ncv = ncv + 1

             ! First define 'kbase' as the first interface below the lower-most unstable interface
             ! Thus, Richardson number at 'kbase' is positive.

             kbase(i,ncv) = MIN(k+1,pver+1)

             ! Decrement k until top unstable level

             DO WHILE( riex(k) .LT. rimaxentr .AND. k .GT. ntop_turb + 2 )
                k = k - 1
             END DO

             ! ktop is the first interface above upper-most unstable interface
             ! Thus, Richardson number at 'ktop' is positive. 

             ktop(i,ncv) = k

          ELSE

             ! Search upward for a CL.

             k = k - 1

          END IF

       END DO ! End of CL regime finding for each atmospheric column

       ncvfin(i) = ncv    

    END DO  ! End of atmospheric column do loop

    RETURN 

  END SUBROUTINE exacol


  !
  !  zisocl
  !
  !============================================================================== !
  !                                                                               !
  !============================================================================== !

  SUBROUTINE zisocl( pcols  , pver  , long ,                                 & 
       z      , zi    , n2   ,  s2      ,                      & 
       bprod  , sprod , bflxs,  tkes    ,landfrac,                      & 
       ncvfin , kbase , ktop ,  belongcv,                      & 
       ricl   , ghcl  , shcl ,  smcl    ,                      &
       lbrk   , wbrk  , ebrk ,  extend  , extend_up, extend_dn )

    !------------------------------------------------------------------------ !
    ! Object : This 'zisocl' vertically extends original CLs identified from  !
    !          'exacol' using a merging test based on either 'wint' or 'l2n2' !
    !          and identify new CL regimes. Similar to the case of 'exacol',  !
    !          CL regime index increases with height.  After identifying new  !
    !          CL regimes ( kbase, ktop, ncvfin ),calculate CL internal mean  !
    !          energetics (lbrk : energetic thickness integral, wbrk, ebrk )  !
    !          and stability functions (ricl, ghcl, shcl, smcl) by including  !
    !          surface interfacial layer contribution when bflxs > 0.   Note  !
    !          that there are two options in the treatment of the energetics  !
    !          of surface interfacial layer (use_dw_surf= 'true' or 'false')  !
    ! Author : Sungsu Park 08/2006, 11/2008                                   !
    !------------------------------------------------------------------------ !

    IMPLICIT NONE

    ! --------------- !    
    ! Input variables !
    ! --------------- !

    INTEGER,  INTENT(in)   :: long                    ! Longitude of the column
    INTEGER,  INTENT(in)   :: pcols                   ! Number of atmospheric columns   
    INTEGER,  INTENT(in)   :: pver                    ! Number of atmospheric vertical layers   
    REAL(r8), INTENT(in)   :: z(pcols, pver)          ! Layer mid-point height [ m ]
    REAL(r8), INTENT(in)   :: zi(pcols, pver+1)       ! Interface height [ m ]
    REAL(r8), INTENT(in)   :: n2(pcols, pver)         ! Buoyancy frequency at interfaces except surface [ s-2 ]
    REAL(r8), INTENT(in)   :: s2(pcols, pver)         ! Shear frequency at interfaces except surface [ s-2 ]
    REAL(r8), INTENT(in)   :: bprod(pcols,pver+1)     ! Buoyancy production [ m2/s3 ]. bprod(i,pver+1) = bflxs 
    REAL(r8), INTENT(in)   :: sprod(pcols,pver+1)     ! Shear production [ m2/s3 ]. sprod(i,pver+1) = usta**3/(vk*z(i,pver))
    REAL(r8), INTENT(in)   :: bflxs(pcols)            ! Surface buoyancy flux [ m2/s3 ]. bprod(i,pver+1) = bflxs 
    REAL(r8), INTENT(in)   :: tkes(pcols)             ! TKE at the surface [ s2/s2 ]
    REAL(r8), INTENT(in)   :: landfrac(pcols)
    ! ---------------------- !
    ! Input/output variables !
    ! ---------------------- !

    INTEGER, INTENT(inout) :: kbase(pcols,ncvmax)     ! Base external interface index of CL
    INTEGER, INTENT(inout) :: ktop(pcols,ncvmax)      ! Top external interface index of CL
    INTEGER, INTENT(inout) :: ncvfin(pcols)           ! Total number of CLs

    ! ---------------- !
    ! Output variables !
    ! ---------------- !

    LOGICAL,  INTENT(out) :: belongcv(pcols,pver+1)   ! True if interface is in a CL ( either internal or external )
    REAL(r8), INTENT(out) :: ricl(pcols,ncvmax)       ! Mean Richardson number of internal CL
    REAL(r8), INTENT(out) :: ghcl(pcols,ncvmax)       ! Half of normalized buoyancy production of internal CL
    REAL(r8), INTENT(out) :: shcl(pcols,ncvmax)       ! Galperin instability function of heat-moisture of internal CL
    REAL(r8), INTENT(out) :: smcl(pcols,ncvmax)       ! Galperin instability function of momentum of internal CL
    REAL(r8), INTENT(out) :: lbrk(pcols,ncvmax)       ! Thickness of (energetically) internal CL ( lint, [m] )
    REAL(r8), INTENT(out) :: wbrk(pcols,ncvmax)       ! Mean normalized TKE of internal CL  [ m2/s2 ]
    REAL(r8), INTENT(out) :: ebrk(pcols,ncvmax)       ! Mean TKE of internal CL ( b1*wbrk, [m2/s2] )

    ! ------------------ !
    ! Internal variables !
    ! ------------------ !

    LOGICAL               :: extend                   ! True when CL is extended in zisocl
    LOGICAL               :: extend_up                ! True when CL is extended upward in zisocl
    LOGICAL               :: extend_dn                ! True when CL is extended downward in zisocl
    LOGICAL               :: bottom                   ! True when CL base is at surface ( kb = pver + 1 )

    INTEGER               :: i                        ! Local index for the longitude
    INTEGER               :: ncv                      ! CL Index increasing with height
    INTEGER               :: incv
    INTEGER               :: k
    INTEGER               :: kb                       ! Local index for kbase
    INTEGER               :: kt                       ! Local index for ktop
    INTEGER               :: ncvinit                  ! Value of ncv at routine entrance 
    INTEGER               :: cntu                     ! Number of merged CLs during upward   extension of individual CL
    INTEGER               :: cntd                     ! Number of merged CLs during downward extension of individual CL
    INTEGER               :: kbinc                    ! Index for incorporating underlying CL
    INTEGER               :: ktinc                    ! Index for incorporating  overlying CL

    REAL(r8)              :: wint                     ! Normalized TKE of internal CL
    REAL(r8)              :: dwinc                    ! Normalized TKE of CL external interfaces
    REAL(r8)              :: dw_surf                  ! Normalized TKE of surface interfacial layer
    REAL(r8)              :: dzinc
    REAL(r8)              :: gh
    REAL(r8)              :: sh
    REAL(r8)              :: sm
    REAL(r8)              :: gh_surf                  ! Half of normalized buoyancy production in surface interfacial layer 
    REAL(r8)              :: sh_surf                  ! Galperin instability function in surface interfacial layer  
    REAL(r8)              :: sm_surf                  ! Galperin instability function in surface interfacial layer 
    REAL(r8)              :: l2n2                     ! Vertical integral of 'l^2N^2' over CL. Include thickness product
    REAL(r8)              :: l2s2                     ! Vertical integral of 'l^2S^2' over CL. Include thickness product
    REAL(r8)              :: dl2n2                    ! Vertical integration of 'l^2*N^2' of CL external interfaces
    REAL(r8)              :: dl2s2                    ! Vertical integration of 'l^2*S^2' of CL external interfaces
    REAL(r8)              :: dl2n2_surf               ! 'dl2n2' defined in the surface interfacial layer
    REAL(r8)              :: dl2s2_surf               ! 'dl2s2' defined in the surface interfacial layer  
    REAL(r8)              :: lint                     ! Thickness of (energetically) internal CL
    REAL(r8)              :: dlint                    ! Interfacial layer thickness of CL external interfaces
    REAL(r8)              :: dlint_surf               ! Surface interfacial layer thickness 
    REAL(r8)              :: lbulk                    ! Master Length Scale : Whole CL thickness from top to base external interface
    REAL(r8)              :: lz                       ! Turbulent length scale
    REAL(r8)              :: ricll                    ! Mean Richardson number of internal CL 
    REAL(r8)              :: trma
    REAL(r8)              :: trmb
    REAL(r8)              :: trmc
    REAL(r8)              :: det
    REAL(r8)              :: zbot                     ! Height of CL base
    REAL(r8)              :: l2rat                    ! Square of ratio of actual to initial CL (not used)
    REAL(r8)              :: gg                       ! Intermediate variable used for calculating stability functions of SBCL
    REAL(r8)              :: tunlramp                 ! Ramping tunl

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- ! 

    i = long

    ! Initialize main output variables

    DO k = 1, ncvmax
       ricl(i,k) = 0._r8
       ghcl(i,k) = 0._r8
       shcl(i,k) = 0._r8
       smcl(i,k) = 0._r8
       lbrk(i,k) = 0._r8
       wbrk(i,k) = 0._r8
       ebrk(i,k) = 0._r8
    END DO
    belongcv(i,1:pver+1)=.TRUE.
    extend    = .FALSE.
    extend_up = .FALSE.
    extend_dn = .FALSE.

    ! ----------------------------------------------------------- !
    ! Loop over each CL to see if any of them need to be extended !
    ! ----------------------------------------------------------- !

    ncv = 1

    DO WHILE( ncv .LE. ncvfin(i) )

       ncvinit = ncv
       cntu    = 0
       cntd    = 0
       kb      = kbase(i,ncv) 
       kt      = ktop(i,ncv)

       ! ---------------------------------------------------------------------------- !
       ! Calculation of CL interior energetics including surface before extension     !
       ! ---------------------------------------------------------------------------- !
       ! Note that the contribution of interior interfaces (not surface) to 'wint' is !
       ! accounted by using '-sh*l2n2 + sm*l2s2' while the contribution of surface is !
       ! accounted by using 'dwsurf = tkes/b1' when bflxs > 0. This approach is fully !
       ! reasonable. Another possible alternative,  which seems to be also consistent !
       ! is to calculate 'dl2n2_surf'  and  'dl2s2_surf' of surface interfacial layer !
       ! separately, and this contribution is explicitly added by initializing 'l2n2' !
       ! 'l2s2' not by zero, but by 'dl2n2_surf' and 'ds2n2_surf' below.  At the same !
       ! time, 'dwsurf' should be excluded in 'wint' calculation below. The only diff.!
       ! between two approaches is that in case of the latter approach, contributions !
       ! of surface interfacial layer to the CL mean stability function (ri,gh,sh,sm) !
       ! are explicitly included while the first approach is not. In this sense,  the !
       ! second approach seems to be more conceptually consistent,   but currently, I !
       ! (Sungsu) will keep the first default approach. There is a switch             !
       ! 'use_dw_surf' at the first part of eddy_diff.F90 chosing one of              !
       ! these two options.                                                           !
       ! ---------------------------------------------------------------------------- !

       ! ------------------------------------------------------ !   
       ! Step 0: Calculate surface interfacial layer energetics !
       ! ------------------------------------------------------ !

       lbulk      = zi(i,kt) - zi(i,kb)
       dlint_surf = 0._r8
       dl2n2_surf = 0._r8
       dl2s2_surf = 0._r8
       dw_surf    = 0._r8
       IF( kb .EQ. pver+1 ) THEN

          IF( bflxs(i) .GT. 0._r8 ) THEN

             ! Calculate stability functions of surface interfacial layer
             ! from the given 'bprod(i,pver+1)' and 'sprod(i,pver+1)' using
             ! inverse approach. Since alph5>0 and alph3<0, denominator of
             ! gg is always positive if bprod(i,pver+1)>0.               

             gg    = 0.5_r8*vk*z(i,pver)*bprod(i,pver+1)/(tkes(i)**(3._r8/2._r8))
             gh    = gg/(alph5-gg*alph3)
             ! gh    = min(max(gh,-0.28_r8),0.0233_r8)
             gh    = MIN(MAX(gh,-3.5334_r8),0.0233_r8)
             sh    = alph5/(1._r8+alph3*gh)
             sm    = (alph1 + alph2*gh)/(1._r8+alph3*gh)/(1._r8+alph4*gh)
             ricll = MIN(-(sm/sh)*(bprod(i,pver+1)/sprod(i,pver+1)),ricrit)

             ! Calculate surface interfacial layer contribution to CL internal
             ! energetics. By construction, 'dw_surf = -dl2n2_surf + ds2n2_surf'
             ! is exactly satisfied, which corresponds to assuming turbulent
             ! length scale of surface interfacial layer = vk * z(i,pver). Note
             ! 'dl2n2_surf','dl2s2_surf','dw_surf' include thickness product.   

             dlint_surf = z(i,pver)
             dl2n2_surf = -vk*(z(i,pver)**2)*bprod(i,pver+1)/(sh*SQRT(tkes(i)))
             dl2s2_surf =  vk*(z(i,pver)**2)*sprod(i,pver+1)/(sm*SQRT(tkes(i)))
             dw_surf    = (tkes(i)/b1)*z(i,pver) 

          ELSE

             ! Note that this case can happen when surface is an external 
             ! interface of CL.
             lbulk = zi(i,kt) - z(i,pver)

          END IF

       END IF

       ! ------------------------------------------------------ !   
       ! Step 1: Include surface interfacial layer contribution !
       ! ------------------------------------------------------ !

       lint = dlint_surf
       l2n2 = dl2n2_surf
       l2s2 = dl2s2_surf          
       wint = dw_surf
       IF( use_dw_surf ) THEN
          l2n2 = 0._r8
          l2s2 = 0._r8
       ELSE
          IF(landfrac(i) > 0.5_r8 ) THEN
             l2n2 = 0._r8
             l2s2 = 0._r8
          ELSE
             l2n2 = 0._r8
             l2s2 = 0._r8
             wint = 0._r8
          END IF
       END IF

       ! --------------------------------------------------------------------------------- !
       ! Step 2. Include the contribution of 'pure internal interfaces' other than surface !
       ! --------------------------------------------------------------------------------- ! 

       IF( kt .LT. kb - 1 ) THEN ! The case of non-SBCL.

          DO k = kb - 1, kt + 1, -1       
             IF( choice_tunl .EQ. 'rampcl' ) THEN
                ! Modification : I simply used the average tunlramp between the two limits.
                tunlramp = 0.5_r8*(1._r8+ctunl)*tunl
             ELSEIF( choice_tunl .EQ. 'rampsl' ) THEN
                tunlramp = ctunl*tunl
                ! tunlramp = 0.765_r8
             ELSE
                tunlramp = tunl
             ENDIF
             IF( choice_leng .EQ. 'origin' ) THEN
                lz = ( (vk*zi(i,k))**(-cleng) + (tunlramp*lbulk)**(-cleng) )**(-1._r8/cleng)
                ! lz = vk*zi(i,k) / (1._r8+vk*zi(i,k)/(tunlramp*lbulk))
             ELSE
                lz = MIN( vk*zi(i,k), tunlramp*lbulk )              
             ENDIF
             dzinc = z(i,k-1) - z(i,k)
             l2n2  = l2n2 + lz*lz*n2(i,k)*dzinc
             l2s2  = l2s2 + lz*lz*s2(i,k)*dzinc
             lint  = lint + dzinc
          END DO

          ! Calculate initial CL stability functions (gh,sh,sm) and net
          ! internal energy of CL including surface contribution if any. 

          ! Modification : It seems that below cannot be applied when ricrit > 0.19.
          !                May need future generalization.

          ricll = MIN(l2n2/MAX(l2s2,ntzero),ricrit) ! Mean Ri of internal CL
          trma  = alph3*alph4*ricll+2._r8*b1*(alph2-alph4*alph5*ricll)
          trmb  = ricll*(alph3+alph4)+2._r8*b1*(-alph5*ricll+alph1)
          trmc  = ricll
          det   = MAX(trmb*trmb-4._r8*trma*trmc,0._r8)
          gh    = (-trmb + SQRT(det))/2._r8/trma
          ! gh    = min(max(gh,-0.28_r8),0.0233_r8)
          gh    = MIN(MAX(gh,-3.5334_r8),0.0233_r8)
          sh    = alph5/(1._r8+alph3*gh)
          sm    = (alph1 + alph2*gh)/(1._r8+alph3*gh)/(1._r8+alph4*gh)
          wint  = wint - sh*l2n2 + sm*l2s2 

       ELSE ! The case of SBCL

          ! If there is no pure internal interface, use only surface interfacial
          ! values. However, re-set surface interfacial values such  that it can
          ! be used in the merging tests (either based on 'wint' or 'l2n2')  and
          ! in such that surface interfacial energy is not double-counted.
          ! Note that regardless of the choise of 'use_dw_surf', below should be
          ! kept as it is below, for consistent merging test of extending SBCL. 

          lint = dlint_surf
          l2n2 = dl2n2_surf
          l2s2 = dl2s2_surf 
          wint = dw_surf

          ! Aug.29.2006 : Only for the purpose of merging test of extending SRCL
          ! based on 'l2n2', re-define 'l2n2' of surface interfacial layer using
          ! 'wint'. This part is designed for similar treatment of merging as in
          ! the original 'eddy_diff.F90' code,  where 'l2n2' of SBCL was defined
          ! as 'l2n2 = - wint / sh'. Note that below block is used only when (1)
          ! surface buoyancy production 'bprod(i,pver+1)' is NOT included in the
          ! calculation of surface TKE in the initialization of 'bprod(i,pver+1)'
          ! in the main subroutine ( even though bflxs > 0 ), and (2) to force 
          ! current scheme be similar to the previous scheme in the treatment of  
          ! extending-merging test of SBCL based on 'l2n2'. Otherwise below line
          ! must be commented out. Note at this stage, correct non-zero value of
          ! 'sh' has been already computed.      

          IF( choice_tkes .EQ. 'ebprod' ) THEN
             l2n2 = - wint / sh 
          ENDIF

       ENDIF

       ! Set consistent upper limits on 'l2n2' and 'l2s2'. Below limits are
       ! reasonable since l2n2 of CL interior interface is always negative.

       l2n2 = -MIN(-l2n2, tkemax*lint/(b1*sh))
       l2s2 =  MIN( l2s2, tkemax*lint/(b1*sm))

       ! Note that at this stage, ( gh, sh, sm )  are the values of surface
       ! interfacial layer if there is no pure internal interface, while if
       ! there is pure internal interface, ( gh, sh, sm ) are the values of
       ! pure CL interfaces or the values that include both the CL internal
       ! interfaces and surface interfaces, depending on the 'use_dw_surf'.       

       ! ----------------------------------------------------------------------- !
       ! Perform vertical extension-merging process                              !
       ! ----------------------------------------------------------------------- !
       ! During the merging process, we assumed ( lbulk, sh, sm ) of CL external !
       ! interfaces are the same as the ones of the original merging CL. This is !
       ! an inevitable approximation since we don't know  ( sh, sm ) of external !
       ! interfaces at this stage.     Note that current default merging test is !
       ! purely based on buoyancy production without including shear production, !
       ! since we used 'l2n2' instead of 'wint' as a merging parameter. However, !
       ! merging test based on 'wint' maybe conceptually more attractable.       !
       ! Downward CL merging process is identical to the upward merging process, !
       ! but when the base of extended CL reaches to the surface, surface inter  !
       ! facial layer contribution to the energetic of extended CL must be done  !
       ! carefully depending on the sign of surface buoyancy flux. The contribu  !
       ! tion of surface interfacial layer energetic is included to the internal !
       ! energetics of merging CL only when bflxs > 0.                           !
       ! ----------------------------------------------------------------------- !

       ! ---------------------------- !
       ! Step 1. Extend the CL upward !
       ! ---------------------------- !

       extend = .FALSE.    ! This will become .true. if CL top or base is extended

       ! Calculate contribution of potentially incorporable CL top interface

       IF( choice_tunl .EQ. 'rampcl' ) THEN
          tunlramp = 0.5_r8*(1._r8+ctunl)*tunl
       ELSEIF( choice_tunl .EQ. 'rampsl' ) THEN
          tunlramp = ctunl*tunl
          ! tunlramp = 0.765_r8
       ELSE
          tunlramp = tunl
       ENDIF
       IF( choice_leng .EQ. 'origin' ) THEN
          lz = ( (vk*zi(i,kt))**(-cleng) + (tunlramp*lbulk)**(-cleng) )**(-1._r8/cleng)
          ! lz = vk*zi(i,kt) / (1._r8+vk*zi(i,kt)/(tunlramp*lbulk))
       ELSE
          lz = MIN( vk*zi(i,kt), tunlramp*lbulk )              
       ENDIF

       dzinc = z(i,kt-1)-z(i,kt)
       dl2n2 = lz*lz*n2(i,kt)*dzinc
       dl2s2 = lz*lz*s2(i,kt)*dzinc
       dwinc = -sh*dl2n2 + sm*dl2s2

       ! ------------ !
       ! Merging Test !
       ! ------------ !

       ! do while (  dwinc .gt. ( rinc*dzinc*wint/(lint+(1._r8-rinc)*dzinc)) )  ! Merging test based on wint
       ! do while ( -dl2n2 .gt. (-rinc*dzinc*l2n2/(lint+(1._r8-rinc)*dzinc)) )  ! Merging test based on l2n2 
       DO WHILE ( -dl2n2 .GT. (-rinc*l2n2/(1._r8-rinc)) )                     ! Integral merging test

          ! Add contribution of top external interface to interior energy.
          ! Note even when we chose 'use_dw_surf='true.', the contribution
          ! of surface interfacial layer to 'l2n2' and 'l2s2' are included
          ! here. However it is not double counting of surface interfacial
          ! energy : surface interfacial layer energy is counted in 'wint'
          ! formula and 'l2n2' is just used for performing merging test in
          ! this 'do while' loop.     

          lint = lint + dzinc
          l2n2 = l2n2 + dl2n2
          l2n2 = -MIN(-l2n2, tkemax*lint/(b1*sh))
          l2s2 = l2s2 + dl2s2
          wint = wint + dwinc

          ! Extend top external interface of CL upward after merging

          kt        = kt - 1
          extend    = .TRUE.
          extend_up = .TRUE.
          IF( kt .EQ. ntop_turb ) THEN
             WRITE(iulog,*) 'zisocl: Error: Tried to extend CL to the model top'
              kt=ktop(i,1)
              EXIT
              !             STOP
          END IF

          ! If the top external interface of extending CL is the same as the 
          ! top interior interface of the overlying CL, overlying CL will be
          ! automatically merged. Then,reduce total number of CL regime by 1. 
          ! and increase 'cntu'(number of merged CLs during upward extension)
          ! by 1.

          ktinc = kbase(i,ncv+cntu+1) - 1  ! Lowest interior interface of overlying CL

          IF( kt .EQ. ktinc ) THEN

             DO k = kbase(i,ncv+cntu+1) - 1, ktop(i,ncv+cntu+1) + 1, -1

                IF( choice_tunl .EQ. 'rampcl' ) THEN
                   tunlramp = 0.5_r8*(1._r8+ctunl)*tunl
                ELSEIF( choice_tunl .EQ. 'rampsl' ) THEN
                   tunlramp = ctunl*tunl
                   ! tunlramp = 0.765_r8
                ELSE
                   tunlramp = tunl
                ENDIF
                IF( choice_leng .EQ. 'origin' ) THEN
                   lz = ( (vk*zi(i,k))**(-cleng) + (tunlramp*lbulk)**(-cleng) )**(-1._r8/cleng)
                   ! lz = vk*zi(i,k) / (1._r8+vk*zi(i,k)/(tunlramp*lbulk))
                ELSE
                   lz = MIN( vk*zi(i,k), tunlramp*lbulk )              
                ENDIF

                dzinc = z(i,k-1)-z(i,k)
                dl2n2 = lz*lz*n2(i,k)*dzinc
                dl2s2 = lz*lz*s2(i,k)*dzinc
                dwinc = -sh*dl2n2 + sm*dl2s2

                lint = lint + dzinc
                l2n2 = l2n2 + dl2n2
                l2n2 = -MIN(-l2n2, tkemax*lint/(b1*sh))
                l2s2 = l2s2 + dl2s2
                wint = wint + dwinc

             END DO

             kt        = ktop(i,ncv+cntu+1) 
             ncvfin(i) = ncvfin(i) - 1
             cntu      = cntu + 1

          END IF

          ! Again, calculate the contribution of potentially incorporatable CL
          ! top external interface of CL regime.

          IF( choice_tunl .EQ. 'rampcl' ) THEN
             tunlramp = 0.5_r8*(1._r8+ctunl)*tunl
          ELSEIF( choice_tunl .EQ. 'rampsl' ) THEN
             tunlramp = ctunl*tunl
             ! tunlramp = 0.765_r8
          ELSE
             tunlramp = tunl
          ENDIF
          IF( choice_leng .EQ. 'origin' ) THEN
             lz = ( (vk*zi(i,kt))**(-cleng) + (tunlramp*lbulk)**(-cleng) )**(-1._r8/cleng)
             ! lz = vk*zi(i,kt) / (1._r8+vk*zi(i,kt)/(tunlramp*lbulk))
          ELSE
             lz = MIN( vk*zi(i,kt), tunlramp*lbulk )              
          ENDIF

          dzinc = z(i,kt-1)-z(i,kt)
          dl2n2 = lz*lz*n2(i,kt)*dzinc
          dl2s2 = lz*lz*s2(i,kt)*dzinc
          dwinc = -sh*dl2n2 + sm*dl2s2

       END DO   ! End of upward merging test 'do while' loop

       ! Update CL interface indices appropriately if any CL was merged.
       ! Note that below only updated the interface index of merged CL,
       ! not the original merging CL.  Updates of 'kbase' and 'ktop' of 
       ! the original merging CL  will be done after finishing downward
       ! extension also later.

       IF( cntu .GT. 0 ) THEN
          DO incv = 1, ncvfin(i) - ncv
             kbase(i,ncv+incv) = kbase(i,ncv+cntu+incv)
             ktop(i,ncv+incv)  = ktop(i,ncv+cntu+incv)
          END DO
       END IF

       ! ------------------------------ !
       ! Step 2. Extend the CL downward !
       ! ------------------------------ !

       IF( kb .NE. pver + 1 ) THEN

          ! Calculate contribution of potentially incorporable CL base interface

          IF( choice_tunl .EQ. 'rampcl' ) THEN
             tunlramp = 0.5_r8*(1._r8+ctunl)*tunl
          ELSEIF( choice_tunl .EQ. 'rampsl' ) THEN
             tunlramp = ctunl*tunl
             ! tunlramp = 0.765_r8
          ELSE
             tunlramp = tunl
          ENDIF
          IF( choice_leng .EQ. 'origin' ) THEN
             lz = ( (vk*zi(i,kb))**(-cleng) + (tunlramp*lbulk)**(-cleng) )**(-1._r8/cleng)
             ! lz = vk*zi(i,kb) / (1._r8+vk*zi(i,kb)/(tunlramp*lbulk))
          ELSE
             lz = MIN( vk*zi(i,kb), tunlramp*lbulk )              
          ENDIF

          dzinc = z(i,kb-1)-z(i,kb)
          dl2n2 = lz*lz*n2(i,kb)*dzinc
          dl2s2 = lz*lz*s2(i,kb)*dzinc
          dwinc = -sh*dl2n2 + sm*dl2s2

          ! ------------ ! 
          ! Merging test !
          ! ------------ ! 

          ! In the below merging tests, I must keep '.and.(kb.ne.pver+1)',   
          ! since 'kb' is continuously updated within the 'do while' loop  
          ! whenever CL base is merged.

          ! do while( (  dwinc .gt. ( rinc*dzinc*wint/(lint+(1._r8-rinc)*dzinc)) ) &  ! Merging test based on wint
          ! do while( ( -dl2n2 .gt. (-rinc*dzinc*l2n2/(lint+(1._r8-rinc)*dzinc)) ) &  ! Merging test based on l2n2
          !             .and.(kb.ne.pver+1))
          DO WHILE( ( -dl2n2 .GT. (-rinc*l2n2/(1._r8-rinc)) ) &                     ! Integral merging test
               .AND.(kb.NE.pver+1))

             ! Add contributions from interfacial layer kb to CL interior 

             lint = lint + dzinc
             l2n2 = l2n2 + dl2n2
             l2n2 = -MIN(-l2n2, tkemax*lint/(b1*sh))
             l2s2 = l2s2 + dl2s2
             wint = wint + dwinc

             ! Extend the base external interface of CL downward after merging

             kb        =  kb + 1
             extend    = .TRUE.
             extend_dn = .TRUE.

             ! If the base external interface of extending CL is the same as the 
             ! base interior interface of the underlying CL, underlying CL  will
             ! be automatically merged. Then, reduce total number of CL by 1. 
             ! For a consistent treatment with 'upward' extension,  I should use
             ! 'kbinc = kbase(i,ncv-1) - 1' instead of 'ktop(i,ncv-1) + 1' below.
             ! However, it seems that these two methods produce the same results.
             ! Note also that in contrast to upward merging, the decrease of ncv
             ! should be performed here.
             ! Note that below formula correctly works even when upperlying CL 
             ! regime incorporates below SBCL.

             kbinc = 0
             IF( ncv .GT. 1 ) kbinc = ktop(i,ncv-1) + 1
             IF( kb .EQ. kbinc ) THEN

                DO k =  ktop(i,ncv-1) + 1, kbase(i,ncv-1) - 1

                   IF( choice_tunl .EQ. 'rampcl' ) THEN
                      tunlramp = 0.5_r8*(1._r8+ctunl)*tunl
                   ELSEIF( choice_tunl .EQ. 'rampsl' ) THEN
                      tunlramp = ctunl*tunl
                      ! tunlramp = 0.765_r8
                   ELSE
                      tunlramp = tunl
                   ENDIF
                   IF( choice_leng .EQ. 'origin' ) THEN
                      lz = ( (vk*zi(i,k))**(-cleng) + (tunlramp*lbulk)**(-cleng) )**(-1._r8/cleng)
                      ! lz = vk*zi(i,k) / (1._r8+vk*zi(i,k)/(tunlramp*lbulk))
                   ELSE
                      lz = MIN( vk*zi(i,k), tunlramp*lbulk )              
                   ENDIF

                   dzinc = z(i,k-1)-z(i,k)
                   dl2n2 = lz*lz*n2(i,k)*dzinc
                   dl2s2 = lz*lz*s2(i,k)*dzinc
                   dwinc = -sh*dl2n2 + sm*dl2s2

                   lint = lint + dzinc
                   l2n2 = l2n2 + dl2n2
                   l2n2 = -MIN(-l2n2, tkemax*lint/(b1*sh))
                   l2s2 = l2s2 + dl2s2
                   wint = wint + dwinc

                END DO

                ! We are incorporating interior of CL ncv-1, so merge
                ! this CL into the current CL.

                kb        = kbase(i,ncv-1)
                ncv       = ncv - 1
                ncvfin(i) = ncvfin(i) -1
                cntd      = cntd + 1

             END IF

             ! Calculate the contribution of potentially incorporatable CL
             ! base external interface. Calculate separately when the base
             ! of extended CL is surface and non-surface.

             IF( kb .EQ. pver + 1 ) THEN 

                IF( bflxs(i) .GT. 0._r8 ) THEN 
                   ! Calculate stability functions of surface interfacial layer
                   gg = 0.5_r8*vk*z(i,pver)*bprod(i,pver+1)/(tkes(i)**(3._r8/2._r8))
                   gh_surf = gg/(alph5-gg*alph3)
                   ! gh_surf = min(max(gh_surf,-0.28_r8),0.0233_r8)
                   gh_surf = MIN(MAX(gh_surf,-3.5334_r8),0.0233_r8)
                   sh_surf = alph5/(1._r8+alph3*gh_surf)
                   sm_surf = (alph1 + alph2*gh_surf)/(1._r8+alph3*gh_surf)/(1._r8+alph4*gh_surf)
                   ! Calculate surface interfacial layer contribution. By construction,
                   ! it exactly becomes 'dw_surf = -dl2n2_surf + ds2n2_surf'  
                   dlint_surf = z(i,pver)
                   dl2n2_surf = -vk*(z(i,pver)**2._r8)*bprod(i,pver+1)/(sh_surf*SQRT(tkes(i)))
                   dl2s2_surf =  vk*(z(i,pver)**2._r8)*sprod(i,pver+1)/(sm_surf*SQRT(tkes(i)))
                   dw_surf = (tkes(i)/b1)*z(i,pver) 
                ELSE
                   dlint_surf = 0._r8
                   dl2n2_surf = 0._r8
                   dl2s2_surf = 0._r8
                   dw_surf = 0._r8
                END IF
                ! If (kb.eq.pver+1), updating of CL internal energetics should be 
                ! performed here inside of 'do while' loop, since 'do while' loop
                ! contains the constraint of '.and.(kb.ne.pver+1)',so updating of
                ! CL internal energetics cannot be performed within this do while
                ! loop when kb.eq.pver+1. Even though I updated all 'l2n2','l2s2',
                ! 'wint' below, only the updated 'wint' is used in the following
                ! numerical calculation.                
                lint = lint + dlint_surf
                l2n2 = l2n2 + dl2n2_surf
                l2n2 = -MIN(-l2n2, tkemax*lint/(b1*sh))
                l2s2 = l2s2 + dl2s2_surf 
                wint = wint + dw_surf                

             ELSE

                IF( choice_tunl .EQ. 'rampcl' ) THEN
                   tunlramp = 0.5_r8*(1._r8+ctunl)*tunl
                ELSEIF( choice_tunl .EQ. 'rampsl' ) THEN
                   tunlramp = ctunl*tunl
                   ! tunlramp = 0.765_r8
                ELSE
                   tunlramp = tunl
                ENDIF
                IF( choice_leng .EQ. 'origin' ) THEN
                   lz = ( (vk*zi(i,kb))**(-cleng) + (tunlramp*lbulk)**(-cleng) )**(-1._r8/cleng)
                   ! lz = vk*zi(i,kb) / (1._r8+vk*zi(i,kb)/(tunlramp*lbulk))
                ELSE
                   lz = MIN( vk*zi(i,kb), tunlramp*lbulk )              
                ENDIF

                dzinc = z(i,kb-1)-z(i,kb)
                dl2n2 = lz*lz*n2(i,kb)*dzinc
                dl2s2 = lz*lz*s2(i,kb)*dzinc
                dwinc = -sh*dl2n2 + sm*dl2s2

             END IF

          END DO ! End of merging test 'do while' loop

          IF( (kb.EQ.pver+1) .AND. (ncv.NE.1) ) THEN 
             WRITE(iulog,*) 'Major mistake zisocl: the CL based at surface is not indexed 1'
             STOP
          END IF

       END IF   ! Done with bottom extension of CL 

       ! Update CL interface indices appropriately if any CL was merged.
       ! Note that below only updated the interface index of merged CL,
       ! not the original merging CL.  Updates of 'kbase' and 'ktop' of 
       ! the original merging CL  will be done later below. I should 
       ! check in detail if below index updating is correct or not.   

       IF( cntd .GT. 0 ) THEN
          DO incv = 1, ncvfin(i) - ncv
             kbase(i,ncv+incv) = kbase(i,ncvinit+incv)
             ktop(i,ncv+incv)  = ktop(i,ncvinit+incv)
          END DO
       END IF

       ! Sanity check for positive wint.

       IF( wint .LT. 0.01_r8 ) THEN
          wint = 0.01_r8
       END IF

       ! -------------------------------------------------------------------------- !
       ! Finally update CL mean internal energetics including surface contribution  !
       ! after finishing all the CL extension-merging process.  As mentioned above, !
       ! there are two possible ways in the treatment of surface interfacial layer, !
       ! either through 'dw_surf' or 'dl2n2_surf and dl2s2_surf' by setting logical !
       ! variable 'use_dw_surf' =.true. or .false.    In any cases, we should avoid !
       ! double counting of surface interfacial layer and one single consistent way !
       ! should be used throughout the program.                                     !
       ! -------------------------------------------------------------------------- !

       IF( extend ) THEN

          ktop(i,ncv)  = kt
          kbase(i,ncv) = kb

          ! ------------------------------------------------------ !   
          ! Step 1: Include surface interfacial layer contribution !
          ! ------------------------------------------------------ !        

          lbulk      = zi(i,kt) - zi(i,kb)
          dlint_surf = 0._r8
          dl2n2_surf = 0._r8
          dl2s2_surf = 0._r8
          dw_surf    = 0._r8
          IF( kb .EQ. pver + 1 ) THEN
             IF( bflxs(i) .GT. 0._r8 ) THEN
                ! Calculate stability functions of surface interfacial layer
                gg = 0.5_r8*vk*z(i,pver)*bprod(i,pver+1)/(tkes(i)**(3._r8/2._r8))
                gh = gg/(alph5-gg*alph3)
                ! gh = min(max(gh,-0.28_r8),0.0233_r8)
                gh = MIN(MAX(gh,-3.5334_r8),0.0233_r8)
                sh = alph5/(1._r8+alph3*gh)
                sm = (alph1 + alph2*gh)/(1._r8+alph3*gh)/(1._r8+alph4*gh)
                ! Calculate surface interfacial layer contribution. By construction,
                ! it exactly becomes 'dw_surf = -dl2n2_surf + ds2n2_surf'  
                dlint_surf = z(i,pver)
                dl2n2_surf = -vk*(z(i,pver)**2._r8)*bprod(i,pver+1)/(sh*SQRT(tkes(i)))
                dl2s2_surf =  vk*(z(i,pver)**2._r8)*sprod(i,pver+1)/(sm*SQRT(tkes(i)))
                dw_surf    = (tkes(i)/b1)*z(i,pver) 
             ELSE
                lbulk = zi(i,kt) - z(i,pver)
             END IF
          END IF
          lint = dlint_surf
          l2n2 = dl2n2_surf
          l2s2 = dl2s2_surf
          wint = dw_surf
          IF( use_dw_surf ) THEN
             l2n2 = 0._r8
             l2s2 = 0._r8
          ELSE
             IF(landfrac(i) > 0.5_r8 ) THEN
                l2n2 = 0._r8
                l2s2 = 0._r8
             ELSE
                l2n2 = 0._r8
                l2s2 = 0._r8
                wint = 0._r8
             END IF
          END IF

          ! -------------------------------------------------------------- !
          ! Step 2. Include the contribution of 'pure internal interfaces' !
          ! -------------------------------------------------------------- ! 

          DO k = kt + 1, kb - 1
             IF( choice_tunl .EQ. 'rampcl' ) THEN
                tunlramp = 0.5_r8*(1._r8+ctunl)*tunl
             ELSEIF( choice_tunl .EQ. 'rampsl' ) THEN
                tunlramp = ctunl*tunl
                ! tunlramp = 0.765_r8
             ELSE
                tunlramp = tunl
             ENDIF
             IF( choice_leng .EQ. 'origin' ) THEN
                lz = ( (vk*zi(i,k))**(-cleng) + (tunlramp*lbulk)**(-cleng) )**(-1._r8/cleng)
                ! lz = vk*zi(i,k) / (1._r8+vk*zi(i,k)/(tunlramp*lbulk))
             ELSE
                lz = MIN( vk*zi(i,k), tunlramp*lbulk )              
             ENDIF
             dzinc = z(i,k-1) - z(i,k)
             lint = lint + dzinc
             l2n2 = l2n2 + lz*lz*n2(i,k)*dzinc
             l2s2 = l2s2 + lz*lz*s2(i,k)*dzinc
          END DO

          ricll = MIN(l2n2/MAX(l2s2,ntzero),ricrit)
          trma = alph3*alph4*ricll+2._r8*b1*(alph2-alph4*alph5*ricll)
          trmb = ricll*(alph3+alph4)+2._r8*b1*(-alph5*ricll+alph1)
          trmc = ricll
          det = MAX(trmb*trmb-4._r8*trma*trmc,0._r8)
          gh = (-trmb + SQRT(det))/2._r8/trma
          ! gh = min(max(gh,-0.28_r8),0.0233_r8)
          gh = MIN(MAX(gh,-3.5334_r8),0.0233_r8)
          sh = alph5 / (1._r8+alph3*gh)
          sm = (alph1 + alph2*gh)/(1._r8+alph3*gh)/(1._r8+alph4*gh)
          ! Even though the 'wint' after finishing merging was positive, it is 
          ! possible that re-calculated 'wint' here is negative.  In this case,
          ! correct 'wint' to be a small positive number
          wint = MAX( wint - sh*l2n2 + sm*l2s2, 0.01_r8 )

       END IF

       ! ---------------------------------------------------------------------- !
       ! Calculate final output variables of each CL (either has merged or not) !
       ! ---------------------------------------------------------------------- !

       lbrk(i,ncv) = lint
       wbrk(i,ncv) = wint/lint
       ebrk(i,ncv) = b1*wbrk(i,ncv)
       ebrk(i,ncv) = MIN(ebrk(i,ncv),tkemax)
       ricl(i,ncv) = ricll 
       ghcl(i,ncv) = gh 
       shcl(i,ncv) = sh
       smcl(i,ncv) = sm

       ! Increment counter for next CL. I should check if the increament of 'ncv'
       ! below is reasonable or not, since whenever CL is merged during downward
       ! extension process, 'ncv' is lowered down continuously within 'do' loop.
       ! But it seems that below 'ncv = ncv + 1' is perfectly correct.

       ncv = ncv + 1

    END DO                   ! End of loop over each CL regime, ncv.

    ! ---------------------------------------------------------- !
    ! Re-initialize external interface indices which are not CLs !
    ! ---------------------------------------------------------- !

    DO ncv = ncvfin(i) + 1, ncvmax
       ktop(i,ncv)  = 0
       kbase(i,ncv) = 0
    END DO

    ! ------------------------------------------------ !
    ! Update CL interface identifiers, 'belongcv'      !
    ! CL external interfaces are also identified as CL !
    ! ------------------------------------------------ !

    DO k = 1, pver + 1
       belongcv(i,k) = .FALSE.
    END DO

    DO ncv = 1, ncvfin(i)
       DO k = ktop(i,ncv), kbase(i,ncv)
          belongcv(i,k) = .TRUE.
       END DO
    END DO

    RETURN

  END SUBROUTINE zisocl

  !
  !  compute_cubic
  !
  REAL(Kind=r8) FUNCTION compute_cubic(a,b,c)
    ! ------------------------------------------------------------------------- !
    ! Solve canonical cubic : x^3 + a*x^2 + b*x + c = 0,  x = sqrt(e)/sqrt(<e>) !
    ! Set x = max(xmin,x) at the end                                            ! 
    ! ------------------------------------------------------------------------- !
    IMPLICIT NONE
    REAL(Kind=r8), INTENT(in)     :: a
    REAL(Kind=r8), INTENT(in)     :: b
    REAL(Kind=r8), INTENT(in)     :: c
    REAL(Kind=r8)                 :: qq, rr, dd, theta, aa, bb, x1, x2, x3
    REAL(Kind=r8), PARAMETER      :: xmin = 1.e-2_r8

    qq = (a**2-3.0_r8*b)/9.0_r8 
    rr = (2.0_r8*a**3 - 9.0_r8*a*b + 27.0_r8*c)/54.0_r8

    dd = rr**2 - qq**3
    IF( dd .LE. 0._r8 ) THEN
       theta = ACOS(rr/qq**(3.0_r8/2.0_r8))
       x1 = -2.0_r8*SQRT(qq)*COS(theta/3.0_r8) - a/3.0_r8
       x2 = -2.0_r8*SQRT(qq)*COS((theta+2.0_r8*3.141592)/3.0_r8) - a/3.0_r8
       x3 = -2.0_r8*SQRT(qq)*COS((theta-2.0_r8*3.141592)/3.0_r8) - a/3.0_r8
       compute_cubic = MAX(MAX(MAX(x1,x2),x3),xmin)        
       RETURN
    ELSE
       IF( rr .GE. 0.0_r8 ) THEN
          aa = -(SQRT(rr**2-qq**3)+rr)**(1.0_r8/3.0_r8)
       ELSE
          aa =  (SQRT(rr**2-qq**3)-rr)**(1.0_r8/3.0_r8)
       ENDIF
       IF( aa .EQ. 0.0_r8 ) THEN
          bb = 0.0_r8
       ELSE
          bb = qq/aa
       ENDIF
       compute_cubic = MAX((aa+bb)-a/3.0_r8,xmin) 
       RETURN
    ENDIF

    RETURN
  END FUNCTION compute_cubic

  !===============================================================================

  SUBROUTINE ubc_get_vals (lchnk,pcols, ncol,pverp, ntop_molec, pint, zi, msis_temp, ubc_mmr)
    !-----------------------------------------------------------------------
    ! interface routine for vertical diffusion and pbl scheme
    !-----------------------------------------------------------------------

    !------------------------------Arguments--------------------------------
    INTEGER,  INTENT(in)  :: lchnk                 ! chunk identifier
    INTEGER,  INTENT(in)  :: ncol                  ! number of atmospheric columns
    INTEGER,  INTENT(in)  :: pcols
    INTEGER,  INTENT(in)  :: pverp    
    INTEGER,  INTENT(in)  :: ntop_molec            ! top of molecular diffusion region (=1)
    REAL(r8), INTENT(in)  :: pint(pcols,pverp)     ! interface pressures
    REAL(r8), INTENT(in)  :: zi(pcols,pverp)       ! interface geoptl height above sfc

    REAL(r8), INTENT(out) :: ubc_mmr(pcols,pcnst)  ! upper bndy mixing ratios (kg/kg)
    REAL(r8), INTENT(out) :: msis_temp(pcols)      ! upper bndy temperature (K)

  END SUBROUTINE ubc_get_vals

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  SUBROUTINE init_molec_diff( kind, ncnst, rair_in, ntop_molec_in, nbot_molec_in, &
                              mw_dry_in, n_avog_in, gravit_in, cpair_in, kbtz_in )
    
   ! use constituents,     only : cnst_mw
    !use upper_bc,         only : ubc_init
    
    INTEGER,  INTENT(in)  :: kind           ! Kind of reals being passed in
    INTEGER,  INTENT(in)  :: ncnst          ! Number of constituents
    INTEGER,  INTENT(in)  :: ntop_molec_in  ! Top interface level to which molecular vertical diffusion is applied ( = 1 )
    INTEGER,  INTENT(in)  :: nbot_molec_in  ! Bottom interface level to which molecular vertical diffusion is applied.
    REAL(r8), INTENT(in)  :: rair_in
    REAL(r8), INTENT(in)  :: mw_dry_in      ! Molecular weight of dry air
    REAL(r8), INTENT(in)  :: n_avog_in      ! Avogadro's number [ molec/kmol ]
    REAL(r8), INTENT(in)  :: gravit_in
    REAL(r8), INTENT(in)  :: cpair_in
    REAL(r8), INTENT(in)  :: kbtz_in        ! Boltzman constant
    INTEGER               :: m              ! Constituent index
    
    IF( kind .NE. r8 ) THEN
        WRITE(0,*) 'KIND of reals passed to init_molec_diff -- exiting.'
        STOP 'init_molec_diff'
    ENDIF
    
    !rair       = rair_in
    !mw_dry     = mw_dry_in
    !n_avog     = n_avog_in
    !gravit     = gravit_in
    !cpair      = cpair_in
    !kbtz       = kbtz_in
    !ntop_molec = ntop_molec_in
    !nbot_molec = nbot_molec_in
    
  ! Initialize upper boundary condition variables

    CALL ubc_init()

  ! Molecular weight factor in constitutent diffusivity
  ! ***** FAKE THIS FOR NOW USING MOLECULAR WEIGHT OF DRY AIR FOR ALL TRACERS ****
  ! !d0=> Diffusion factor [ m-1 s-1 ] molec sqrt(kg/kmol/K) [ unit ? ]
    ALLOCATE(mw_fac(ncnst))
    DO m = 1, ncnst
       mw_fac(m) = d0 * mw_dry_in * SQRT(1._r8/mw_dry_in + 1._r8/cnst_mw(m)) / n_avog_in
    END DO

  END SUBROUTINE init_molec_diff
  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  INTEGER FUNCTION compute_molec_diff( lchnk       ,                                                                      &
       pcols       , pver       , ncnst      , ncol     , t      , pmid   , pint        , &
       zi          , ztodt      , kvh        , kvm      , tint   , rhoi   , tmpi2       , &
       kq_scal     , ubc_t      , ubc_mmr    , dse_top  , cc_top , cd_top , cnst_mw_out , &
       cnst_fixed_ubc_out , mw_fac_out , ntop_molec_out , nbot_molec_out )

    !use upper_bc,        only : ubc_get_vals
    !use constituents,    only : cnst_mw, cnst_fixed_ubc

    ! --------------------- !
    ! Input-Output Argument !
    ! --------------------- !

    INTEGER,  INTENT(in)    :: pcols
    INTEGER,  INTENT(in)    :: pver
    INTEGER,  INTENT(in)    :: ncnst
    INTEGER,  INTENT(in)    :: ncol                      ! Number of atmospheric columns
    INTEGER,  INTENT(in)    :: lchnk                     ! Chunk identifier
    REAL(r8), INTENT(in)    :: t(pcols,pver)             ! Temperature input
    REAL(r8), INTENT(in)    :: pmid(pcols,pver)          ! Midpoint pressures
    REAL(r8), INTENT(in)    :: pint(pcols,pver+1)        ! Interface pressures
    REAL(r8), INTENT(in)    :: zi(pcols,pver+1)          ! Interface heights
    REAL(r8), INTENT(in)    :: ztodt                     ! 2 delta-t

    REAL(r8), INTENT(inout) :: kvh(pcols,pver+1)         ! Diffusivity for heat
    REAL(r8), INTENT(inout) :: kvm(pcols,pver+1)         ! Viscosity ( diffusivity for momentum )
    REAL(r8), INTENT(inout) :: tint(pcols,pver+1)        ! Interface temperature
    REAL(r8), INTENT(inout) :: rhoi(pcols,pver+1)        ! Density ( rho ) at interfaces
    REAL(r8), INTENT(inout) :: tmpi2(pcols,pver+1)       ! dt*(g*rho)**2/dp at interfaces

    REAL(r8), INTENT(out)   :: kq_scal(pcols,pver+1)     ! kq_fac*sqrt(T)*m_d/rho for molecular diffusivity
    REAL(r8), INTENT(out)   :: ubc_mmr(pcols,ncnst)      ! Upper boundary mixing ratios [ kg/kg ]
    REAL(r8), INTENT(out)   :: cnst_mw_out(ncnst)
    LOGICAL,  INTENT(out)   :: cnst_fixed_ubc_out(ncnst)
    REAL(r8), INTENT(out)   :: mw_fac_out(ncnst)
    REAL(r8), INTENT(out)   :: dse_top(pcols)            ! dse on top boundary
    REAL(r8), INTENT(out)   :: cc_top(pcols)             ! Lower diagonal at top interface
    REAL(r8), INTENT(out)   :: cd_top(pcols)             ! cc_top * dse ubc value
    INTEGER,  INTENT(out)   :: ntop_molec_out   
    INTEGER,  INTENT(out)   :: nbot_molec_out   

    ! --------------- ! 
    ! Local variables !
    ! --------------- !

    INTEGER                 :: m                          ! Constituent index
    INTEGER                 :: i                          ! Column index
    INTEGER                 :: k                          ! Level index
    REAL(r8)                :: mkvisc                     ! Molecular kinematic viscosity c*tint**(2/3)/rho
    REAL(r8)                :: ubc_t(pcols)               ! Upper boundary temperature (K)

    REAL(r8), PARAMETER   :: pwr    = 2._r8/3._r8        ! Exponentiation factor [ unit ? ]
    REAL(r8), PARAMETER   :: pr_num = 1._r8              ! Prandtl number [ no unit ]
    REAL(r8), PARAMETER   :: km_fac = 3.55E-7_r8         ! Molecular viscosity constant [ unit ? ]

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !
    kq_scal(1:pcols,1:pver+1)    =0.0_r8 ! kq_fac*sqrt(T)*m_d/rho for molecular diffusivity
    ubc_mmr(1:pcols,1:ncnst)     =0.0_r8  ! Upper boundary mixing ratios [ kg/kg ]
    cnst_mw_out(1:ncnst)       =0.0_r8 
    cnst_fixed_ubc_out(1:ncnst)=.TRUE.
    mw_fac_out(1:ncnst)=0.0_r8 
    dse_top(1:pcols)=0.0_r8 ! dse on top boundary
    cc_top(1:pcols)=0.0_r8 ! Lower diagonal at top interface
    cd_top(1:pcols)=0.0_r8 ! cc_top * dse ubc value
    ntop_molec_out=0
    nbot_molec_out=0
    mkvisc=0.0_r8
    ubc_t(1:pcols)=0.0_r8

    ! Get upper boundary values

    CALL ubc_get_vals( lchnk,pcols, ncol,pver+1, ntop_molec, pint, zi, ubc_t, ubc_mmr )

    ! Below are already computed, just need to be copied for output

    DO m=1,ncnst
       cnst_mw_out       (m)        = cnst_mw(m)
       cnst_fixed_ubc_out(m) = cnst_fixed_ubc(m)
       mw_fac_out        (m)         = mw_fac(m)
    END DO
    ntop_molec_out             = ntop_molec
    nbot_molec_out             = nbot_molec

    ! Density and related factors for moecular diffusion and ubc.
    ! Always have a fixed upper boundary T if molecular diffusion is active. Why ?
    DO i=1,ncol
       tint (i,ntop_molec) = ubc_t(i)
       rhoi (i,ntop_molec) = pint(i,ntop_molec) / ( rair * tint(i,ntop_molec) )
       tmpi2(i,ntop_molec) = ztodt * ( gravit * rhoi(i,ntop_molec))**2 &
                              / ( pmid(i,ntop_molec) - pint(i,ntop_molec) )
    END DO
    ! Compute molecular kinematic viscosity, heat diffusivity and factor for constituent diffusivity
    ! This is a key part of the code.
    DO k=1,ntop_molec-1
       DO i=1,ncol
          kq_scal(i,k) = 0._r8
       END DO
    END DO
    DO k = ntop_molec, nbot_molec
       DO i = 1, ncol
          mkvisc       = km_fac * tint(i,k)**pwr / rhoi(i,k)
          kvm(i,k)     = kvm(i,k) + mkvisc
          kvh(i,k)     = kvh(i,k) + mkvisc * pr_num
          kq_scal(i,k) = SQRT(tint(i,k)) / rhoi(i,k)
       END DO
    END DO
    DO k=nbot_molec+1,pver+1
       DO i = 1, ncol
          kq_scal(i,k) = 0._r8
       END DO
    END DO
    ! Top boundary condition for dry static energy
    DO i = 1, ncol
       dse_top(i) = cpair * tint(i,ntop_molec) + gravit * zi(i,ntop_molec)
    END DO 
    ! Top value of cc for dry static energy

    DO i = 1, ncol
       cc_top(i) = ztodt * gravit**2 * rhoi(i,ntop_molec) * km_fac * ubc_t(i)**pwr / &
            ( ( pint(i,2) - pint(i,1) ) * ( pmid(i,1) - pint(i,1) ) )
    ENDDO
    DO i = 1, ncol
       cd_top(i) = cc_top(i) * dse_top(i)
    END DO
    compute_molec_diff = 1
    RETURN
  END FUNCTION compute_molec_diff



  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  SUBROUTINE compute_vdiff( &
       lchnk               ,& !integer,  intent(in)    :: lchnk                   
       pcols               ,& !integer,  intent(in)    :: pcols
       pver                ,& !integer,  intent(in)    :: pver
       ncnst               ,& !integer,  intent(in)    :: ncnst
       ncol                ,& !integer,  intent(in)    :: ncol                      ! Number of atmospheric columns
       TSK                 , &! INTENT(IN   ) surface temperature
       QSFC                , &! INTENT(IN   ) surface specific temperature
       pmid                ,& !real(r8), intent(in)    :: pmid(pcols,pver)          ! Mid-point pressures [ Pa ]
       pint                ,& !real(r8), intent(in)    :: pint(pcols,pver+1)        ! Interface pressures [ Pa ]
       rpdel               ,& !real(r8), intent(in)    :: rpdel(pcols,pver)         ! 1./pdel
       t                   ,& !real(r8), intent(in)    :: t(pcols,pver)             ! Temperature [ K ]
       ztodt               ,& !real(r8), intent(in)    :: ztodt                     ! 2 delta-t [ s ]
       taux                ,& !real(r8), intent(in)    :: taux(pcols)               ! Surface zonal      stress. Input u-momentum per unit time per unit area into the atmosphere [ N/m2 ]
       tauy                ,& !real(r8), intent(in)    :: tauy(pcols)               ! Surface meridional stress. Input v-momentum per unit time per unit area into the atmosphere [ N/m2 ]
       shflx               ,& !real(r8), intent(in)    :: shflx(pcols)              ! Surface sensible heat flux [ W/m2 ]
       cflx                ,& !real(r8), intent(in)    :: cflx(pcols,ncnst)         ! Surface constituent flux [ kg/m2/s ]
       ntop                ,& !integer,  intent(in)    :: ntop                      ! Top    interface level to which vertical diffusion is applied ( = 1 ).
       nbot                ,& !integer,  intent(in)    :: nbot                      ! Bottom interface level to which vertical diffusion is applied ( = pver ).
       kvh                 ,& !real(r8), intent(inout) :: kvh(pcols,pver+1)         ! Eddy diffusivity for heat [ m2/s ]
       kvm                 ,& !real(r8), intent(inout) :: kvm(pcols,pver+1)         ! Eddy viscosity ( Eddy diffusivity for momentum ) [ m2/s ]
       kvq                 ,& !real(r8), intent(inout) :: kvq(pcols,pver+1)         ! Eddy diffusivity for constituents
       cgs                 ,& !real(r8), intent(inout) :: cgs(pcols,pver+1)         ! Counter-gradient star [ cg/flux ]
       cgh                 ,& !real(r8), intent(inout) :: cgh(pcols,pver+1)         ! Counter-gradient term for heat
       zi                  ,& !real(r8), intent(in)    :: zi(pcols,pver+1)          ! Interface heights [ m ]
       ksrftms             ,& !real(r8), intent(in)    :: ksrftms(pcols)            ! Surface drag coefficient for turbulent mountain stress. > 0. [ kg/s/m2 ]
       qmincg              ,& !real(r8), intent(in)    :: qmincg(ncnst)             ! Minimum constituent mixing ratios from cg fluxes
       fieldlist           ,& !type(vdiff_selector), intent(in) :: fieldlist        ! Array of flags selecting which fields to diffuse
       u                   ,& !real(r8), intent(inout) :: u(pcols,pver)             ! U wind. This input is the 'raw' input wind to PBL scheme without iterative provisional update. [ m/s ]
       v                   ,& !real(r8), intent(inout) :: v(pcols,pver)             ! V wind. This input is the 'raw' input wind to PBL scheme without iterative provisional update. [ m/s ]
       q                   ,& !real(r8), intent(inout) :: q(pcols,pver,ncnst)       ! Moisture and trace constituents [ kg/kg, #/kg ? ]
       dse                 ,& !real(r8), intent(inout) :: dse(pcols,pver)           ! Dry static energy [ J/kg ]
       tautmsx             ,& !real(r8), intent(inout) :: tauresx(pcols)            ! Input  : Reserved surface stress at previous time step
       tautmsy             ,& !real(r8), intent(inout) :: tauresy(pcols)            ! Output : Reserved surface stress at current  time step
       dtk                 ,& !real(r8), intent(out)   :: dtk(pcols,pver)           ! T tendency from KE dissipation
       topflx              ,& !real(r8), intent(out)   :: topflx(pcols)             ! Molecular heat flux at the top interface
       errstring           ,& !character(128), intent(out) :: errstring             ! Output status
       tauresx             ,& !real(r8), intent(out)   :: tautmsx(pcols)            ! Implicit zonal      turbulent mountain surface stress [ N/m2 = kg m/s /s/m2 ]
       tauresy             ,& !real(r8), intent(out)   :: tautmsy(pcols)            ! Implicit meridional turbulent mountain surface stress [ N/m2 = kg m/s /s/m2 ]
       itaures             )!,& !integer,  intent(in)    :: itaures                   ! Indicator determining whether 'tauresx,tauresy' is updated (1) or non-updated (0) in this subroutine.   
    !                           do_molec_diff       ,& !
    !                           compute_molec_diff  ,& !
    ! compute_molec_diff  ) !integer,  external, optional :: compute_molec_diff   ! Constituent-independent moleculuar diffusivity routine

    !-------------------------------------------------------------------------- !
    ! Driver routine to compute vertical diffusion of momentum, moisture, trace !
    ! constituents and dry static energy. The new temperature is computed from  !
    ! the diffused dry static energy.                                           ! 
    ! Turbulent diffusivities and boundary layer nonlocal transport terms are   !
    ! obtained from the turbulence module.                                      !
    !-------------------------------------------------------------------------- !

    !use phys_debug_util,    only : phys_debug_col
    !use time_manager,       only : is_first_step, get_nstep
    !use phys_control,       only : phys_getopts

    ! Modification : Ideally, we should diffuse 'liquid-ice static energy' (sl), not the dry static energy.
    !                Also, vertical diffusion of cloud droplet number concentration and aerosol number
    !                concentration should be done very carefully in the future version.

    ! --------------- !
    ! Input Arguments !
    ! --------------- !

    INTEGER,  INTENT(in)    :: lchnk                   
    INTEGER,  INTENT(in)    :: pcols
    INTEGER,  INTENT(in)    :: pver
    INTEGER,  INTENT(in)    :: ncnst
    INTEGER,  INTENT(in)    :: ncol                      ! Number of atmospheric columns
    INTEGER,  INTENT(in)    :: ntop                      ! Top    interface level to which vertical diffusion is applied ( = 1 ).
    INTEGER,  INTENT(in)    :: nbot                      ! Bottom interface level to which vertical diffusion is applied ( = pver ).
    INTEGER,  INTENT(in)    :: itaures                   ! Indicator determining whether 'tauresx,tauresy' is updated (1) or non-updated (0) in this subroutine.   
    REAL(r8), INTENT(in)    :: TSK  (pcols)          ! INTENT(IN   ) surface temperature
    REAL(r8), INTENT(in)    :: QSFC (pcols)          ! INTENT(IN   ) surface specific temperature
    REAL(r8), INTENT(in)    :: pmid(pcols,pver)          ! Mid-point pressures [ Pa ]
    REAL(r8), INTENT(in)    :: pint(pcols,pver+1)        ! Interface pressures [ Pa ]
    REAL(r8), INTENT(in)    :: rpdel(pcols,pver)         ! 1./pdel
    REAL(r8), INTENT(in)    :: t(pcols,pver)             ! Temperature [ K ]
    REAL(r8), INTENT(in)    :: ztodt                     ! 2 delta-t [ s ]
    REAL(r8), INTENT(in)    :: taux(pcols)               ! Surface zonal      stress. Input u-momentum per unit time per unit area into the atmosphere [ N/m2 ]
    REAL(r8), INTENT(in)    :: tauy(pcols)               ! Surface meridional stress. Input v-momentum per unit time per unit area into the atmosphere [ N/m2 ]
    REAL(r8), INTENT(in)    :: shflx(pcols)              ! Surface sensible heat flux [ W/m2 ]
    REAL(r8), INTENT(in)    :: cflx(pcols,ncnst)         ! Surface constituent flux [ kg/m2/s ]
    REAL(r8), INTENT(in)    :: zi(pcols,pver+1)          ! Interface heights [ m ]
    REAL(r8), INTENT(in)    :: ksrftms(pcols)            ! Surface drag coefficient for turbulent mountain stress. > 0. [ kg/s/m2 ]
    REAL(r8), INTENT(in)    :: qmincg(ncnst)             ! Minimum constituent mixing ratios from cg fluxes

    !logical,  intent(in)         :: do_molec_diff        ! Flag indicating multiple constituent diffusivities
    !integer,  external, optional :: compute_molec_diff   ! Constituent-independent moleculuar diffusivity routine
    TYPE(vdiff_selector), INTENT(in) :: fieldlist        ! Array of flags selecting which fields to diffuse

    ! ---------------------- !
    ! Input-Output Arguments !
    ! ---------------------- !

    REAL(r8), INTENT(inout) :: kvh(pcols,pver+1)         ! Eddy diffusivity for heat [ m2/s ]
    REAL(r8), INTENT(inout) :: kvm(pcols,pver+1)         ! Eddy viscosity ( Eddy diffusivity for momentum ) [ m2/s ]
    REAL(r8), INTENT(inout) :: kvq(pcols,pver+1)         ! Eddy diffusivity for constituents
    REAL(r8), INTENT(inout) :: cgs(pcols,pver+1)         ! Counter-gradient star [ cg/flux ]
    REAL(r8), INTENT(inout) :: cgh(pcols,pver+1)         ! Counter-gradient term for heat

    REAL(r8), INTENT(inout) :: u(pcols,pver)             ! U wind. This input is the 'raw' input wind to PBL scheme without iterative provisional update. [ m/s ]
    REAL(r8), INTENT(inout) :: v(pcols,pver)             ! V wind. This input is the 'raw' input wind to PBL scheme without iterative provisional update. [ m/s ]
    REAL(r8), INTENT(inout) :: q(pcols,pver,ncnst)       ! Moisture and trace constituents [ kg/kg, #/kg ? ]
    REAL(r8), INTENT(inout) :: dse(pcols,pver)           ! Dry static energy [ J/kg ]

    REAL(r8), INTENT(inout) :: tauresx(pcols)            ! Input  : Reserved surface stress at previous time step
    REAL(r8), INTENT(inout) :: tauresy(pcols)            ! Output : Reserved surface stress at current  time step

    ! ---------------- !
    ! Output Arguments !
    ! ---------------- !

    REAL(r8), INTENT(out)   :: dtk(pcols,pver)           ! T tendency from KE dissipation
    REAL(r8), INTENT(out)   :: tautmsx(pcols)            ! Implicit zonal      turbulent mountain surface stress [ N/m2 = kg m/s /s/m2 ]
    REAL(r8), INTENT(out)   :: tautmsy(pcols)            ! Implicit meridional turbulent mountain surface stress [ N/m2 = kg m/s /s/m2 ]
    REAL(r8), INTENT(out)   :: topflx(pcols)             ! Molecular heat flux at the top interface
    CHARACTER(128), INTENT(out) :: errstring             ! Output status

    ! --------------- !
    ! Local Variables ! 
    ! --------------- !

    INTEGER  :: i, k, m, icol                            ! Longitude, level, constituent indices
    INTEGER  :: status                                   ! Status indicator
    INTEGER  :: ntop_molec                               ! Top level where molecular diffusivity is applied
    LOGICAL  :: lqtst(pcols)                             ! Adjust vertical profiles
    LOGICAL  :: need_decomp                              ! Whether to compute a new decomposition
    LOGICAL  :: cnst_fixed_ubc(ncnst)                    ! Whether upper boundary condition is fixed
    !logical  :: do_iss                                   ! Use implicit turbulent surface stress computation

    REAL(r8) :: tmpm(pcols,pver)                         ! Potential temperature, ze term in tri-diag sol'n
    REAL(r8) :: ca(pcols,pver)                           ! - Upper diag of tri-diag matrix
    REAL(r8) :: cc(pcols,pver)                           ! - Lower diag of tri-diag matrix
    REAL(r8) :: dnom(pcols,pver)                         ! 1./(1. + ca(k) + cc(k) - cc(k)*ze(k-1))

    REAL(r8) :: tmp1(pcols)                              ! Temporary storage
    REAL(r8) :: tmpi1(pcols,pver+1)                      ! Interface KE dissipation
    REAL(r8) :: tint(pcols,pver+1)                       ! Interface temperature
    REAL(r8) :: rhoi(pcols,pver+1)                       ! rho at interfaces
    REAL(r8) :: tmpi2(pcols,pver+1)                      ! dt*(g*rho)**2/dp at interfaces
    REAL(r8) :: rrho(pcols)                              ! 1./bottom level density 

    REAL(r8) :: zero(pcols)                              ! Zero array for surface heat exchange coefficients 
    REAL(r8) :: tautotx(pcols)                           ! Total surface stress ( zonal )
    REAL(r8) :: tautoty(pcols)                           ! Total surface stress ( meridional )

    REAL(r8) :: dinp_u(pcols,pver+1)                     ! Vertical difference at interfaces, input u
    REAL(r8) :: dinp_v(pcols,pver+1)                     ! Vertical difference at interfaces, input v
    REAL(r8) :: dout_u                                   ! Vertical difference at interfaces, output u
    REAL(r8) :: dout_v                                   ! Vertical difference at interfaces, output v
    REAL(r8) :: dse_top(pcols)                           ! dse on top boundary
    REAL(r8) :: cc_top(pcols)                            ! Lower diagonal at top interface
    REAL(r8) :: cd_top(pcols)                            ! 
    REAL(r8) :: rghd(pcols,pver+1)                       ! (1/H_i - 1/H) *(g*rho)^(-1)

    REAL(r8) :: qtm(pcols,pver)                          ! Temporary copy of q
    REAL(r8) :: kq_scal(pcols,pver+1)                    ! kq_fac*sqrt(T)*m_d/rho for molecular diffusivity
    REAL(r8) :: mw_fac(ncnst)                            ! sqrt(1/M_q + 1/M_d) for this constituent
    REAL(r8) :: cnst_mw(ncnst)                           ! Molecular weight [ kg/kmole ]
    REAL(r8) :: ubc_mmr(pcols,ncnst)                     ! Upper boundary mixing ratios [ kg/kg ]
    REAL(r8) :: ubc_t(pcols)                             ! Upper boundary temperature [ K ]

    REAL(r8) :: ws(pcols)                                ! Lowest-level wind speed [ m/s ]
    REAL(r8) :: tau(pcols)                               ! Turbulent surface stress ( not including mountain stress )
    REAL(r8) :: ksrfturb(pcols)                          ! Surface drag coefficient of 'normal' stress. > 0. Virtual mass input per unit time per unit area [ kg/s/m2 ]
    REAL(r8) :: ksrf(pcols)                              ! Surface drag coefficient of 'normal' stress + Surface drag coefficient of 'tms' stress.  > 0. [ kg/s/m2 ] 
    REAL(r8) :: usum_in(pcols)                           ! Vertical integral of input u-momentum. Total zonal     momentum per unit area in column  [ sum of u*dp/g = kg m/s m-2 ]
    REAL(r8) :: vsum_in(pcols)                           ! Vertical integral of input v-momentum. Total meridional momentum per unit area in column [ sum of v*dp/g = kg m/s m-2 ]
    REAL(r8) :: usum_mid(pcols)                          ! Vertical integral of u-momentum after adding explicit residual stress
    REAL(r8) :: vsum_mid(pcols)                          ! Vertical integral of v-momentum after adding explicit residual stress
    REAL(r8) :: usum_out(pcols)                          ! Vertical integral of u-momentum after doing implicit diffusion
    REAL(r8) :: vsum_out(pcols)                          ! Vertical integral of v-momentum after doing implicit diffusion
    REAL(r8) :: tauimpx(pcols)                           ! Actual net stress added at the current step other than mountain stress
    REAL(r8) :: tauimpy(pcols)                           ! Actual net stress added at the current step other than mountain stress
    REAL(r8) :: wsmin                                    ! Minimum sfc wind speed for estimating frictional transfer velocity ksrf. [ m/s ]
    REAL(r8) :: ksrfmin                                  ! Minimum surface drag coefficient [ kg/s/m^2 ]
    REAL(r8) :: timeres                                  ! Relaxation time scale of residual stress ( >= dt ) [ s ]
    REAL(r8) :: ramda                                    ! dt/timeres [ no unit ]
    REAL(r8) :: psum
    REAL(r8) :: u_in, u_res, tauresx_in
    REAL(r8) :: v_in, v_res, tauresy_in  
    REAL(r8) :: kvh_out(pcols,pver+1) ! Eddy diffusivity for heat [ m2/s ]
    REAL(r8) :: kvm_out(pcols,pver+1) ! Eddy viscosity ( Eddy diffusivity for momentum ) [ m2/s ]
    REAL(r8) :: kvq_out(pcols,pver+1) ! Eddy diffusivity for heat [ m2/s ]
    REAL(r8) :: kvh_in(pcols,pver+1) ! Eddy diffusivity for heat [ m2/s ]
    REAL(r8) :: kvm_in(pcols,pver+1) ! Eddy viscosity ( Eddy diffusivity for momentum ) [ m2/s ]
    REAL(r8) :: kvq_in(pcols,pver+1) ! Eddy diffusivity for heat [ m2/s ]
    REAL(r8) :: ustar(pcols)   

    ! ------------------------------------------------ !
    ! Parameters for implicit surface stress treatment !
    ! ------------------------------------------------ !
    ntop_molec=0
    !wsmin    = 2.0_r8                                     ! Minimum wind speed for ksrfturb computation        [ m/s ]
    wsmin    = 1.0_r8                                     ! Minimum wind speed for ksrfturb computation        [ m/s ]
    ksrfmin  = 1.e-4_r8                                  ! Minimum surface drag coefficient                   [ kg/s/m^2 ]
    timeres  =21600._r8                                   ! Relaxation time scale of residual stress ( >= dt ) [ s ]
    !timeres  = 7200._r8                                  ! Relaxation time scale of residual stress ( >= dt ) [ s ]

    !call phys_getopts( do_iss_out = do_iss )
    dtk=0.0_r8     ! T tendency from KE dissipation
    tautmsx=0.0_r8      ! Implicit zonal      turbulent mountain surface stress [ N/m2 = kg m/s /s/m2 ]
    tautmsy=0.0_r8      ! Implicit meridional turbulent mountain surface stress [ N/m2 = kg m/s /s/m2 ]
    topflx=0.0_r8      ! Molecular heat flux at the top interface
    DO k=1,pver+1
       DO i=1,pcols
          tmpi1 (i,k)=0.0_r8! Interface KE dissipation
          tint  (i,k)=0.0_r8! Interface temperature
          rhoi  (i,k)=0.0_r8! rho at interfaces
          tmpi2 (i,k)=0.0_r8! dt*(g*rho)**2/dp at interfaces
          dinp_u(i,k)=0.0_r8! Vertical difference at interfaces, input u
          dinp_v(i,k)=0.0_r8! Vertical difference at interfaces, input v
          rghd  (i,k)=0.0_r8! (1/H_i - 1/H) *(g*rho)^(-1)
          kq_scal(i,k)=0.0_r8! kq_fac*sqrt(T)*m_d/rho for molecular diffusivity
       END DO
    END DO   

    DO k=1,pver
       DO i=1,pcols
          tmpm(i,k)=0.0_r8! Potential temperature, ze term in tri-diag sol'n
          ca  (i,k)=0.0_r8! - Upper diag of tri-diag matrix
          cc  (i,k)=0.0_r8! - Lower diag of tri-diag matrix
          dnom(i,k)=0.0_r8! 1./(1. + ca(k) + cc(k) - cc(k)*ze(k-1))
          qtm (i,k)=0.0_r8! Temporary copy of q
       END DO
    END DO
    
    DO i=1,pcols
       tmp1    (i) =0.0_r8! Temporary storage
       rrho    (i) =0.0_r8! 1./bottom level density 
       zero    (i) =0.0_r8! Zero array for surface heat exchange coefficients 
       tautotx (i) =0.0_r8! Total surface stress ( zonal )
       tautoty (i) =0.0_r8! Total surface stress ( meridional )
       dse_top (i) =0.0_r8! dse on top boundary
       cc_top  (i) =0.0_r8! Lower diagonal at top interface
       cd_top  (i) =0.0_r8! 
       ws      (i) =0.0_r8! Lowest-level wind speed [ m/s ]
       tau     (i) =0.0_r8! Turbulent surface stress ( not including mountain stress )
       ksrfturb(i) =0.0_r8! Surface drag coefficient of 'normal' stress. > 0. Virtual mass input per unit time per unit area [ kg/s/m2 ]
       ksrf    (i) =0.0_r8! Surface drag coefficient of 'normal' stress + Surface drag coefficient of 'tms' stress.  > 0. [ kg/s/m2 ] 
       usum_in (i) =0.0_r8! Vertical integral of input u-momentum. Total zonal       momentum per unit area in column  [ sum of u*dp/g = kg m/s m-2 ]
       vsum_in (i) =0.0_r8! Vertical integral of input v-momentum. Total meridional momentum per unit area in column [ sum of v*dp/g = kg m/s m-2 ]
       usum_mid(i) =0.0_r8! Vertical integral of u-momentum after adding explicit residual stress
       vsum_mid(i) =0.0_r8! Vertical integral of v-momentum after adding explicit residual stress
       usum_out(i) =0.0_r8! Vertical integral of u-momentum after doing implicit diffusion
       vsum_out(i) =0.0_r8! Vertical integral of v-momentum after doing implicit diffusion
       tauimpx (i) =0.0_r8! Actual net stress added at the current step other than mountain stress
       tauimpy (i) =0.0_r8! Actual net stress added at the current step other than mountain stress
       ubc_t   (i) =0.0_r8! Upper boundary temperature [ K ]
    END DO
    
    DO m=1,ncnst
       mw_fac  (m)  =0.0_r8                          ! sqrt(1/M_q + 1/M_d) for this constituent
       cnst_mw (m)  =0.0_r8                          ! Molecular weight [ kg/kmole ]
       DO i=1,pcols
          ubc_mmr (i,m)=0.0_r8     ! Upper boundary mixing ratios [ kg/kg ]
       END DO
    END DO
    
    ramda=0.0_r8 ! dt/timeres [ no unit ]
    psum=0.0_r8 
    u_in=0.0_r8 ; u_res=0.0_r8 ; tauresx_in=0.0_r8 
    v_in=0.0_r8 ; v_res=0.0_r8 ; tauresy_in=0.0_r8 
    dout_u=0.0_r8 ! Vertical difference at interfaces, output u
    dout_v=0.0_r8 ! Vertical difference at interfaces, output v


    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !
    errstring = ''
    IF( ( diffuse(fieldlist,'u') .OR. diffuse(fieldlist,'v') ) .AND. .NOT. diffuse(fieldlist,'s') ) THEN
       errstring = 'diffusion_solver.compute_vdiff: must diffuse s if diffusing u or v'
       RETURN
    END IF
    zero(:) = 0._r8

    ! Compute 'rho' and 'dt*(g*rho)^2/dp' at interfaces
    DO i = 1, ncol
       tint(i,1) = t(i,1)
       rhoi(i,1) = pint(i,1) / (rair*tint(i,1))
    END DO
    DO k = 2, pver
       DO i = 1, ncol
          tint(i,k)  = 0.5_r8 * ( t(i,k) + t(i,k-1) )
          rhoi(i,k)  = pint(i,k) / (rair*tint(i,k))
          tmpi2(i,k) = ztodt * ( gravit*rhoi(i,k) )**2 / ( pmid(i,k) - pmid(i,k-1) )
       END DO
    END DO
    DO i = 1, ncol
       !tint(i,pver+1) = t(i,pver)
       tint(i,pver+1) = 0.5_r8 * ( t(i,pver) + TSK(i) )
       rhoi(i,pver+1) = pint(i,pver+1) / ( rair*tint(i,pver+1) )
    END DO
    DO i = 1, ncol
       !rrho(i) = rair  * t(i,pver) / pmid(i,pver)
       rrho(i) = rair  * TSK(i) / pmid(i,pver)
       tmp1(i) = ztodt * gravit * rpdel(i,pver)
       ustar(i)   = MAX( SQRT( SQRT( taux(i)**2 + tauy(i)**2 ) * rrho(i) ), ustar_min )

    END DO
    !--------------------------------------- !
    ! Computation of Molecular Diffusivities !
    !--------------------------------------- !

    ! Modification : Why 'kvq' is not changed by molecular diffusion ? 

    IF( do_molec_diff ) THEN

       !if( (.not.present(compute_molec_diff)) ) then
       !      errstring = 'compute_vdiff: do_molec_diff true but compute_molec_diff or vd_lu_qdecomp missing'
       !      return
       !endif

       ! The next subroutine 'compute_molec_diff' :
       !     Modifies : kvh, kvm, tint, rhoi, and tmpi2
       !     Returns  : kq_scal, ubc_t, ubc_mmr, dse_top, cc_top, cd_top, cnst_mw, 
       !                cnst_fixed_ubc , mw_fac , ntop_molec 

       status = compute_molec_diff( lchnk          ,                                                                &
            pcols          , pver    , ncnst      , ncol      , t      , pmid   , pint    , &
            zi             , ztodt   , kvh        , kvm       , tint   , rhoi   , tmpi2   , &
            kq_scal        , ubc_t   , ubc_mmr    , dse_top   , cc_top , cd_top , cnst_mw , &
            cnst_fixed_ubc , mw_fac  , ntop_molec , nbot_molec )

    ELSE

       kq_scal(:,:) = 0._r8
       cd_top(:)    = 0._r8
       cc_top(:)    = 0._r8

    ENDIF
    kvh_in=kvh 
    kvm_in=kvm
    kvq_in=kvq
    IF(ens)THEN
    DO k = 2, pver  ! plev + 1
       DO i = 1, ncol
          !
          !             k=2  ****Km(k),sl*** } -----------
          !             k=3/2----si,ric,rf,km,kh,b,l -----------
          !             k=1  ****Km(k),sl*** } -----------
          !             k=1/2----si ----------------------------
          !
          !        Km(k-1) + 2*Km(k) + Km(k+1)
          ! Km = -------------------------------
          !                 4
          !
          kvm_out(i,k)=(0.25_r8*(kvm(i,k-1)+2.0_r8*kvm(i,k)+kvm(i,k+1)))
          !
          !        Kh(k-1) + 2*Kh(k) + Kh(k+1)
          ! Kh = -------------------------------
          !                 4
          !
          kvh_out(i,k)=(0.25_r8*(kvh(i,k-1)+2.0_r8*kvh(i,k)+kvh(i,k+1)))
          !
          !        Kq(k-1) + 2*Kq(k) + Kq(k+1)
          ! Kq = -------------------------------
          !                 4
          !
          kvq_out(i,k)=(0.25_r8*(kvq(i,k-1)+2.0_r8*kvq(i,k)+kvq(i,k+1)))

       END DO
    END DO
    DO i = 1, ncol
           !
           !          Km(1) + Km(2)
           ! Km = ---------------
           !               2
           !
           kvm_out(i,     1  )=     0.5_r8*  ( kvm(i,     1)+kvm(i,     2))
           !
           !          Kh(1) + Kh(2)
           ! Kh = ---------------
           !               2
           !
           kvh_out(i,     1  )=     0.5_r8*  ( kvh(i,     1)+kvm(i,     2))
           !
           !
           !          Kq(1) + Kq(2)
           ! Kq = ---------------
           !               2
           !
           kvq_out(i,     1  )=     0.5_r8*  ( kvq(i,     1)+kvq(i,     2))
           !
           !          Km(k-1) + 2*Km(k) + Km(k+1)
           ! Km = -------------------------------
           !                     4
           !
           kvm_out(i,pver+1)=(0.25_r8*(kvm(i,pver-1)+2.0_r8*kvm(i,pver+1)+kvm(i,pver  )))
           !
           !          Kh(k-1) + 2*Kh(k) + Kh(k+1)
           ! Kh = -------------------------------
           !                     4
           !
           kvh_out(i,pver+1)=(0.25_r8*(kvh(i,pver-1)+2.0_r8*kvh(i,pver+1)+kvh(i,pver  )))
           !
           !          Kq(k-1) + 2*Kq(k) + Kq(k+1)
           ! Kq = -------------------------------
           !                     4
           !
           kvq_out(i,pver+1)=(0.25_r8*(kvq(i,pver-1)+2.0_r8*kvq(i,pver+1)+kvq(i,pver  )))

    END DO
    DO k = 1, pver+1  ! plev + 1
       DO i = 1, ncol
          kvm(i,k)=kvm_out(i,k)
          kvh(i,k)=kvh_out(i,k)
          kvq(i,k)=kvq_out(i,k)
       END DO
    END DO
    END IF
    !---------------------------- !
    ! Diffuse Horizontal Momentum !
    !---------------------------- !

    IF( diffuse(fieldlist,'u') .OR. diffuse(fieldlist,'v') ) THEN

       ! Compute the vertical upward differences of the input u,v for KE dissipation
       ! at each interface.
       ! Velocity = 0 at surface, so difference at the bottom interface is -u,v(pver)
       ! These 'dinp_u, dinp_v' are computed using the non-diffused input wind.

       DO i = 1, ncol
          dinp_u(i,1) = 0._r8
          dinp_v(i,1) = 0._r8
          dinp_u(i,pver+1) = -u(i,pver)
          dinp_v(i,pver+1) = -v(i,pver)
       END DO
       DO k = 2, pver
          DO i = 1, ncol
             dinp_u(i,k) = u(i,k) - u(i,k-1)
             dinp_v(i,k) = v(i,k) - v(i,k-1)
          END DO
       END DO

       ! -------------------------------------------------------------- !
       ! Do 'Implicit Surface Stress' treatment for numerical stability !
       ! in the lowest model layer.                                     !
       ! -------------------------------------------------------------- !

       IF( do_iss ) THEN

          ! Compute surface drag coefficient for implicit diffusion 
          ! including turbulent turbulent mountain stress. 

          DO i = 1, ncol
             ws(i)       = MAX( SQRT( u(i,pver)**2._r8 + v(i,pver)**2._r8 ), wsmin )
             tau(i)      = SQRT( taux(i)**2._r8 + tauy(i)**2._r8 )
             ksrfturb(i) = MAX( tau(i) / ws(i), ksrfmin )
          END DO
          DO i = 1, ncol

             ksrf(i) = ksrfturb(i) + ksrftms(i)  ! Do all surface stress ( normal + tms ) implicitly

          END DO

          ! Vertical integration of input momentum. 
          ! This is total horizontal momentum per unit area [ kg*m/s/m2 ] in each column.
          ! Note (u,v) are the raw input to the PBL scheme, not the
          ! provisionally-marched ones within the iteration loop of the PBL scheme.  

          DO i = 1, ncol
             usum_in(i) = 0._r8
             vsum_in(i) = 0._r8
             DO k = 1, pver
                usum_in(i) = usum_in(i) + (1._r8/gravit)*u(i,k)/rpdel(i,k)
                vsum_in(i) = vsum_in(i) + (1._r8/gravit)*v(i,k)/rpdel(i,k)
             END DO
          END DO

          ! Add residual stress of previous time step explicitly into the lowest
          ! model layer with a relaxation time scale of 'timeres'.

          ramda         = ztodt / timeres
          DO i = 1, ncol
             u(i,pver) = u(i,pver) + tmp1(i)*tauresx(i)*ramda
             v(i,pver) = v(i,pver) + tmp1(i)*tauresy(i)*ramda
          END DO
          ! Vertical integration of momentum after adding explicit residual stress
          ! into the lowest model layer.

          DO i = 1, ncol
             usum_mid(i) = 0._r8
             vsum_mid(i) = 0._r8
             DO k = 1, pver
                usum_mid(i) = usum_mid(i) + (1._r8/gravit)*u(i,k)/rpdel(i,k)
                vsum_mid(i) = vsum_mid(i) + (1._r8/gravit)*v(i,k)/rpdel(i,k)
             END DO
          END DO

          ! Debug 
          ! icol = phys_debug_col(lchnk) 
          ! if ( icol > 0 .and. get_nstep() .ge. 1 ) then
          !      tauresx_in = tauresx(icol)
          !      tauresy_in = tauresy(icol)
          !      u_in  = u(icol,pver) - tmp1(icol) * tauresx(icol) * ramda
          !      v_in  = v(icol,pver) - tmp1(icol) * tauresy(icol) * ramda
          !      u_res = u(icol,pver)
          !      v_res = v(icol,pver)
          ! endif
          ! Debug

       ELSE

          ! In this case, do 'turbulent mountain stress' implicitly, 
          ! but do 'normal turbulent stress' explicitly.
          ! In this case, there is no 'redisual stress' as long as 'tms' is
          ! treated in a fully implicit wway, which is true.

          ! 1. Do 'tms' implicitly
          DO i = 1, ncol

             ksrf(i) = ksrftms(i) 
          END DO
          ! 2. Do 'normal stress' explicitly
          DO i = 1, ncol

             u(i,pver) = u(i,pver) + tmp1(i)*taux(i)
             v(i,pver) = v(i,pver) + tmp1(i)*tauy(i)
          END DO

       END IF  ! End of 'do iss' ( implicit surface stress )

       ! --------------------------------------------------------------------------------------- !
       ! Diffuse horizontal momentum implicitly using tri-diagnonal matrix.                      !
       ! The 'u,v' are input-output: the output 'u,v' are implicitly diffused winds.             !
       !    For implicit 'normal' stress : ksrf = ksrftms + ksrfturb,                            !
       !                                   u(pver) : explicitly include 'redisual normal' stress !
       !    For explicit 'normal' stress : ksrf = ksrftms                                        !
       !                                   u(pver) : explicitly include 'normal' stress          !                                              
       ! Note that in all the two cases above, 'tms' is fully implicitly treated.                !
       ! --------------------------------------------------------------------------------------- !

       CALL vd_lu_decomp( pcols , pver , ncol  ,                        &
            ksrf  , kvm  , tmpi2 , rpdel , ztodt , zero , &
            ca    , cc   , dnom  , tmpm  , ntop  , nbot )

       CALL vd_lu_solve(  pcols , pver , ncol  ,                        &
            u     , ca   , tmpm  , dnom  , ntop  , nbot , zero )

       CALL vd_lu_solve(  pcols , pver , ncol  ,                        &
            v     , ca   , tmpm  , dnom  , ntop  , nbot , zero )

       ! ---------------------------------------------------------------------- !
       ! Calculate 'total' ( tautotx ) and 'tms' ( tautmsx ) stresses that      !
       ! have been actually added into the atmosphere at the current time step. ! 
       ! Also, update residual stress, if required.                             !
       ! ---------------------------------------------------------------------- !

       DO i = 1, ncol

          ! Compute the implicit 'tms' using the updated winds.
          ! Below 'tautmsx(i),tautmsy(i)' are pure implicit mountain stresses
          ! that has been actually added into the atmosphere both for explicit
          ! and implicit approach. 

          tautmsx(i) = -ksrftms(i)*u(i,pver)
          tautmsy(i) = -ksrftms(i)*v(i,pver)

          IF( do_iss ) THEN

             ! Compute vertical integration of final horizontal momentum

             usum_out(i) = 0._r8
             vsum_out(i) = 0._r8
             DO k = 1, pver
                usum_out(i) = usum_out(i) + (1._r8/gravit)*u(i,k)/rpdel(i,k)
                vsum_out(i) = vsum_out(i) + (1._r8/gravit)*v(i,k)/rpdel(i,k)
             END DO

             ! Compute net stress added into the atmosphere at the current time step.
             ! Note that the difference between 'usum_in' and 'usum_out' are induced
             ! by 'explicit residual stress + implicit total stress' for implicit case, while
             ! by 'explicit normal   stress + implicit tms   stress' for explicit case. 
             ! Here, 'tautotx(i)' is net stress added into the air at the current time step.

             tauimpx(i) = ( usum_out(i) - usum_in(i) ) / ztodt
             tauimpy(i) = ( vsum_out(i) - vsum_in(i) ) / ztodt

             tautotx(i) = tauimpx(i) 
             tautoty(i) = tauimpy(i) 

             ! Compute redisual stress and update if required.
             ! Note that the total stress we should have added at the current step is
             ! the sum of 'taux(i) - ksrftms(i)*u(i,pver) + tauresx(i)'.

             IF( itaures .EQ. 1 ) THEN
                tauresx(i) = taux(i) + tautmsx(i) + tauresx(i) - tauimpx(i)
                tauresy(i) = tauy(i) + tautmsy(i) + tauresy(i) - tauimpy(i)
             ENDIF

          ELSE

             tautotx(i) = tautmsx(i) + taux(i)
             tautoty(i) = tautmsy(i) + tauy(i)
             tauresx(i) = 0._r8
             tauresy(i) = 0._r8

          END IF  ! End of 'do_iss' routine

       END DO ! End of 'do i = 1, ncol' routine

       ! Debug 
       ! icol = phys_debug_col(lchnk) 
       ! if ( icol > 0 .and. get_nstep() .ge. 1 ) then
       !      write(iulog,*)
       !      write(iulog,*)  'diffusion_solver debug'  
       !      write(iulog,*)
       !      write(iulog,*)  'u_in, u_res, u_out'
       !      write(iulog,*)   u_in, u_res, u(icol,pver)
       !      write(iulog,*)  'tauresx_in, tautmsx, tauimpx(actual), tauimpx(derived), tauresx_out, taux'
       !      write(iulog,*)   tauresx_in, tautmsx(icol), tauimpx(icol), -ksrf(icol)*u(icol,pver), tauresx(icol), taux(icol)
       !      write(iulog,*)
       !      write(iulog,*)  'v_in, v_res, v_out'
       !      write(iulog,*)   v_in, v_res, v(icol,pver)
       !      write(iulog,*)  'tauresy_in, tautmsy, tauimpy(actual), tauimpy(derived), tauresy_out, tauy'
       !      write(iulog,*)   tauresy_in, tautmsy(icol), tauimpy(icol), -ksrf(icol)*v(icol,pver), tauresy(icol), tauy(icol)
       !      write(iulog,*)
       !      write(iulog,*)  'itaures, ksrf, ksrfturb, ksrftms'
       !      write(iulog,*)   itaures, ksrf(icol), ksrfturb(icol), ksrftms(icol)
       !      write(iulog,*) 
       ! endif
       ! Debug

       ! ------------------------------------ !
       ! Calculate kinetic energy dissipation !
       ! ------------------------------------ !       

       ! Modification : In future, this should be set exactly same as 
       !                the ones in the convection schemes 

       ! 1. Compute dissipation term at interfaces
       !    Note that 'u,v' are already diffused wind, and 'tautotx,tautoty' are 
       !    implicit stress that has been actually added. On the other hand,
       !    'dinp_u, dinp_v' were computed using non-diffused input wind.

       ! Modification : I should check whether non-consistency between 'u' and 'dinp_u'
       !                is correctly intended approach. I think so.

       k = pver + 1
       DO i = 1, ncol
          tmpi1(i,1) = 0._r8
          tmpi1(i,k) = 0.5_r8 * ztodt * gravit * &
               ( (-u(i,k-1) + dinp_u(i,k))*tautotx(i) + (-v(i,k-1) + dinp_v(i,k))*tautoty(i) )
       END DO

       DO k = 2, pver
          DO i = 1, ncol
             dout_u = u(i,k) - u(i,k-1)
             dout_v = v(i,k) - v(i,k-1)
             tmpi1(i,k) = 0.25_r8 * tmpi2(i,k) * kvm(i,k) * &
                  ( dout_u**2 + dout_v**2 + dout_u*dinp_u(i,k) + dout_v*dinp_v(i,k) )
          END DO
       END DO

       ! 2. Compute dissipation term at midpoints, add to dry static energy

       DO k = 1, pver
          DO i = 1, ncol
             dtk(i,k) = ( tmpi1(i,k+1) + tmpi1(i,k) ) * rpdel(i,k)
             dse(i,k) = dse(i,k) + dtk(i,k)
          END DO
       END DO

    END IF ! End of diffuse horizontal momentum, diffuse(fieldlist,'u') routine

    !-------------------------- !
    ! Diffuse Dry Static Energy !
    !-------------------------- !

    ! Modification : In future, we should diffuse the fully conservative 
    !                moist static energy,not the dry static energy.

    IF( diffuse(fieldlist,'s') ) THEN

       ! Add counter-gradient to input static energy profiles

       DO k = 1, pver
          DO i=1,ncol

             dse(i,k) = dse(i,k) + ztodt * rpdel(i,k) * gravit  *                &
                  ( rhoi(i,k+1) * kvh(i,k+1) * cgh(i,k+1) &
                  - rhoi(i,k  ) * kvh(i,k  ) * cgh(i,k  ) )
           END DO
       END DO

       ! Add the explicit surface fluxes to the lowest layer
       DO i=1,ncol
          dse(i,pver) = dse(i,pver) + tmp1(i) * shflx(i)
       END DO
       ! Diffuse dry static energy

       CALL vd_lu_decomp( pcols , pver , ncol  ,                         &
            zero  , kvh  , tmpi2 , rpdel , ztodt , cc_top, &
            ca    , cc   , dnom  , tmpm  , ntop  , nbot    )

       CALL vd_lu_solve(  pcols , pver , ncol  ,                         &
            dse   , ca   , tmpm  , dnom  , ntop  , nbot  , cd_top )

       ! Calculate flux at top interface

       ! Modification : Why molecular diffusion does not work for dry static energy in all layers ?

       IF( do_molec_diff ) THEN
          DO i=1,ncol
             topflx(i) =  - kvh(i,ntop_molec) * tmpi2(i,ntop_molec) / (ztodt*gravit) * &
                  ( dse(i,ntop_molec) - dse_top(i) )
          END DO
       END IF

    ENDIF

    !---------------------------- !
    ! Diffuse Water Vapor Tracers !
    !---------------------------- !

    ! Modification : For aerosols, I need to use separate treatment 
    !                for aerosol mass and aerosol number. 

    ! Loop through constituents

    need_decomp = .TRUE.

    DO m = 1, ncnst

       IF( diffuse(fieldlist,'q',m) ) THEN

          ! Add the nonlocal transport terms to constituents in the PBL.
          ! Check for neg q's in each constituent and put the original vertical
          ! profile back if a neg value is found. A neg value implies that the
          ! quasi-equilibrium conditions assumed for the countergradient term are
          ! strongly violated.
          DO k = 1, pver
             DO i=1,ncol
                qtm(i,k) = q(i,k,m)
             END DO
          END DO     
          DO k = 1, pver
             DO i=1,ncol
             q(i,k,m) = q(i,k,m) + &
                  ztodt * rpdel(i,k) * gravit  * ( cflx(i,m) * rrho(i) ) * &
                  ( rhoi(i,k+1) * kvh(i,k+1) * cgs(i,k+1)                    &
                  - rhoi(i,k  ) * kvh(i,k  ) * cgs(i,k  ) )
             END DO
          END DO
          DO i=1,ncol
             lqtst(i) = ALL(q(i,1:pver,m) >= qmincg(m))
          END DO
          DO k = 1, pver
             DO i=1,ncol

                q(i,k,m) = MERGE( q(i,k,m), qtm(i,k), lqtst(i) )
             END DO
          END DO

          ! Add the explicit surface fluxes to the lowest layer
          ! tmp1(i) = ztodt * gravit * rpdel(i,pver)

          DO i=1,ncol
             q(i,pver,m) = q(i,pver,m) + tmp1(i) * cflx(i,m)
          END DO

          ! Diffuse constituents.

          IF( need_decomp ) THEN

             CALL vd_lu_decomp( pcols , pver , ncol  ,                         &
                  zero  , kvq  , tmpi2 , rpdel , ztodt , zero  , &
                  ca    , cc   , dnom  , tmpm  , ntop  , nbot )

             IF( do_molec_diff ) THEN

                ! Update decomposition in molecular diffusion range, include separation velocity term

                status = vd_lu_qdecomp( pcols , pver   , ncol      , cnst_fixed_ubc(m), cnst_mw(m), ubc_mmr(:,m), &
                     kvq   , kq_scal, mw_fac(m) , tmpi2            , rpdel     ,               &
                     ca    , cc     , dnom      , tmpm             , rhoi      ,               &
                     tint  , ztodt  , ntop_molec, nbot_molec       , cd_top )
             ELSE
                need_decomp =  .FALSE.
             ENDIF
          END IF

          CALL vd_lu_solve(  pcols , pver , ncol  ,                         &
               q(1:ncol,1:pver,m) , ca, tmpm  , dnom  , ntop  , nbot  , cd_top )
       END IF
    END DO
    kvh= kvh_in
    kvm= kvm_in
    kvq= kvq_in

    RETURN
  END SUBROUTINE compute_vdiff

  !                                                                             !
  !============================================================================ !

  INTEGER FUNCTION vd_lu_qdecomp( pcols , pver   , ncol       , fixed_ubc  , mw     , ubc_mmr , &
       kv    , kq_scal, mw_facm    , tmpi       , rpdel  ,           &
       ca    , cc     , dnom       , ze         , rhoi   ,           &
       tint  , ztodt  , ntop_molec , nbot_molec , cd_top )

    !------------------------------------------------------------------------------ !
    ! Add the molecular diffusivity to the turbulent diffusivity for a consitutent. !
    ! Update the superdiagonal (ca(k)), diagonal (cb(k)) and subdiagonal (cc(k))    !
    ! coefficients of the tridiagonal diffusion matrix, also ze and denominator.    !
    !------------------------------------------------------------------------------ !

    ! ---------------------- !
    ! Input-Output Arguments !
    ! ---------------------- !

    INTEGER,  INTENT(in)    :: pcols
    INTEGER,  INTENT(in)    :: pver
    INTEGER,  INTENT(in)    :: ncol                  ! Number of atmospheric columns

    INTEGER,  INTENT(in)    :: ntop_molec
    INTEGER,  INTENT(in)    :: nbot_molec

    LOGICAL,  INTENT(in)    :: fixed_ubc             ! Fixed upper boundary condition flag
    REAL(r8), INTENT(in)    :: kv(pcols,pver+1)      ! Eddy diffusivity
    REAL(r8), INTENT(in)    :: kq_scal(pcols,pver+1) ! Molecular diffusivity ( kq_fac*sqrt(T)*m_d/rho )
    REAL(r8), INTENT(in)    :: mw                    ! Molecular weight for this constituent
    REAL(r8), INTENT(in)    :: ubc_mmr(pcols)        ! Upper boundary mixing ratios [ kg/kg ]
    REAL(r8), INTENT(in)    :: mw_facm               ! sqrt(1/M_q + 1/M_d) for this constituent
    REAL(r8), INTENT(in)    :: tmpi(pcols,pver+1)    ! dt*(g/R)**2/dp*pi(k+1)/(.5*(tm(k+1)+tm(k))**2
    REAL(r8), INTENT(in)    :: rpdel(pcols,pver)     ! 1./pdel ( thickness bet interfaces )
    REAL(r8), INTENT(in)    :: rhoi(pcols,pver+1)    ! Density at interfaces [ kg/m3 ]
    REAL(r8), INTENT(in)    :: tint(pcols,pver+1)    ! Interface temperature [ K ]
    REAL(r8), INTENT(in)    :: ztodt                 ! 2 delta-t [ s ]

    REAL(r8), INTENT(inout) :: ca(pcols,pver)        ! -Upper diagonal
    REAL(r8), INTENT(inout) :: cc(pcols,pver)        ! -Lower diagonal
    REAL(r8), INTENT(inout) :: dnom(pcols,pver)      ! 1./(1. + ca(k) + cc(k) - cc(k)*ze(k-1)) , 1./(b(k) - c(k)*e(k-1))
    REAL(r8), INTENT(inout) :: ze(pcols,pver)        ! Term in tri-diag. matrix system

    REAL(r8), INTENT(out)   :: cd_top(pcols)         ! Term for updating top level with ubc

    ! --------------- !
    ! Local Variables !
    ! --------------- !

    INTEGER                 :: i                     ! Longitude index
    INTEGER                 :: k, kp1                ! Vertical indicies

    REAL(r8)                :: rghd(pcols,pver+1)    ! (1/H_i - 1/H) * (rho*g)^(-1)
    REAL(r8)                :: kmq(ncol)             ! Molecular diffusivity for constituent
    REAL(r8)                :: wrk0(ncol)            ! Work variable
    REAL(r8)                :: wrk1(ncol)            ! Work variable

    REAL(r8)                :: cb(pcols,pver)                      ! - Diagonal
    REAL(r8)                :: kvq(pcols,pver+1)                   ! Output vertical diffusion coefficient
    REAL(R8), PARAMETER     :: SHR_CONST_BOLTZ   = 1.38065e-23_R8  ! Boltzmann's constant ~ J/K/molecule
    REAL(r8), PARAMETER     :: boltz       = shr_const_boltz       ! Boltzman's constant (J/K/molecule)
    REAL(R8), PARAMETER     :: SHR_CONST_AVOGAD  = 6.02214e26_R8   ! Avogadro's number ~ molecules/kmole
    REAL(r8), PARAMETER     :: avogad      = shr_const_avogad      ! Avogadro's number (molecules/kmole)
    REAL(R8), PARAMETER     :: SHR_CONST_MWDAIR  = 28.966_R8       ! molecular weight dry air ~ kg/kmole
    REAL(r8), PARAMETER      :: mwdry        = shr_const_mwdair    ! molecular weight dry air

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !   

    ! --------------------------------------------------------------------- !
    ! Determine superdiagonal (ca(k)) and subdiagonal (cc(k)) coeffs of the !
    ! tridiagonal diffusion matrix. The diagonal elements  (cb=1+ca+cc) are !
    ! a combination of ca and cc; they are not required by the solver.      !
    !---------------------------------------------------------------------- !

    !call t_startf('vd_lu_qdecomp')

    kvq   (:,:)  = 0.0_r8
    cd_top(:)    = 0.0_r8
    rghd  (:,:)  = 0.0_r8! (1/H_i - 1/H) * (rho*g)^(-1)
    kmq   (:)    = 0.0_r8! Molecular diffusivity for constituent
    wrk0  (:)    = 0.0_r8! Work variable
    wrk1  (:)    = 0.0_r8! Work variable
    cb    (:,:)  = 0.0_r8! - Diagonal

    ! Compute difference between scale heights of constituent and dry air

    DO k = ntop_molec, nbot_molec
       DO i = 1, ncol
          rghd(i,k) = gravit / ( boltz * avogad * tint(i,k) ) * ( mw - mwdry )
          rghd(i,k) = ztodt * gravit * rhoi(i,k) * rghd(i,k) 
       ENDDO
    ENDDO

    !-------------------- !
    ! Molecular diffusion !
    !-------------------- !

    DO k = nbot_molec - 1, ntop_molec, -1
       kp1 = k + 1
       DO i = 1, ncol
          kmq(i)  = kq_scal(i,kp1) * mw_facm
          wrk0(i) = ( kv(i,kp1) + kmq(i) ) * tmpi(i,kp1)
          wrk1(i) = kmq(i) * 0.5_r8 * rghd(i,kp1)
          ! Add species separation term
          ca(i,k  )  = ( wrk0(i) - wrk1(i) ) * rpdel(i,k)
          cc(i,kp1)  = ( wrk0(i) + wrk1(i) ) * rpdel(i,kp1)
          kvq(i,kp1) = kmq(i)
       END DO
    END DO

    IF( fixed_ubc ) THEN
       DO i = 1, ncol
          cc(i,ntop_molec) = kq_scal(i,ntop_molec) * mw_facm                 &
               * ( tmpi(i,ntop_molec) + rghd(i,ntop_molec) ) &
               * rpdel(i,ntop_molec)
       END DO
    END IF

    ! Calculate diagonal elements

    DO k = nbot_molec - 1, ntop_molec + 1, -1
       kp1 = k + 1
        DO i = 1, ncol
          cb(i,k) = 1._r8 + ca(i,k) + cc(i,k)                   &
               + rpdel(i,k) * ( kvq(i,kp1) * rghd(i,kp1) &
               - kvq(i,k) * rghd(i,k) )
          kvq(i,kp1) = kv(i,kp1) + kvq(i,kp1)
       END DO
    END DO

    k   = ntop_molec
    kp1 = k + 1
    IF( fixed_ubc ) THEN
       DO i = 1, ncol
          cb(i,k) = 1.0_r8 + ca(i,k)                                 &
               + rpdel(i,k) * kvq(i,kp1) * rghd(i,kp1)   &
               + kq_scal(i,ntop_molec) * mw_facm                 &
               * ( tmpi(i,ntop_molec) - rghd(i,ntop_molec) ) &
               * rpdel(i,ntop_molec)
       END DO
    ELSE
       DO i = 1, ncol
          cb(i,k) = 1._r8 + ca(i,k) &
               + rpdel(i,k) * kvq(i,kp1) * rghd(i,kp1)
       END DO
    END IF

    k   = nbot_molec
    DO i = 1, ncol
       cb(i,k) = 1._r8 + cc(i,k) + ca(i,k) &
            - rpdel(i,k) * kvq(i,k)*rghd(i,k)
    END DO
    DO k = 1, nbot_molec + 1, -1
       DO i = 1, ncol
          cb(i,k) = 1._r8 + ca(i,k) + cc(i,k)
       END DO
    END DO

    ! Compute term for updating top level mixing ratio for ubc

    IF( fixed_ubc ) THEN
       DO i = 1, ncol
          cd_top(i) = cc(i,ntop_molec) * ubc_mmr(i)
       END DO
    END IF

    !-------------------------------------------------------- !
    ! Calculate e(k).                                         !
    ! This term is required in solution of tridiagonal matrix ! 
    ! defined by implicit diffusion equation.                 !
    !-------------------------------------------------------- !

    DO k = nbot_molec, ntop_molec + 1, -1
       DO i = 1, ncol
          dnom(i,k) = 1._r8 / ( cb(i,k) - ca(i,k) * ze(i,k+1) )
          ze(i,k)   = cc(i,k) * dnom(i,k)
       END DO
    END DO
    
    k = ntop_molec
    DO i = 1, ncol
       dnom(i,k) = 1._r8 / ( cb(i,k) - ca(i,k) * ze(i,k+1) )
    END DO
    vd_lu_qdecomp = 1
    !call t_stopf('vd_lu_qdecomp')
    RETURN

  END FUNCTION vd_lu_qdecomp

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  SUBROUTINE vd_lu_decomp( pcols, pver, ncol ,                        &
       ksrf , kv  , tmpi , rpdel, ztodt , cc_top, &
       ca   , cc  , dnom , ze   , ntop  , nbot    )
    !---------------------------------------------------------------------- !
    ! Determine superdiagonal (ca(k)) and subdiagonal (cc(k)) coeffs of the ! 
    ! tridiagonal diffusion matrix.                                         ! 
    ! The diagonal elements (1+ca(k)+cc(k)) are not required by the solver. !
    ! Also determine ze factor and denominator for ze and zf (see solver).  !
    !---------------------------------------------------------------------- !

    ! --------------------- !
    ! Input-Output Argument !
    ! --------------------- !

    INTEGER,  INTENT(in)  :: pcols                 ! Number of allocated atmospheric columns
    INTEGER,  INTENT(in)  :: pver                  ! Number of allocated atmospheric levels 
    INTEGER,  INTENT(in)  :: ncol                  ! Number of computed atmospheric columns
    INTEGER,  INTENT(in)  :: ntop                  ! Top level to operate on
    INTEGER,  INTENT(in)  :: nbot                  ! Bottom level to operate on
    REAL(r8), INTENT(in)  :: ksrf(pcols)           ! Surface "drag" coefficient [ kg/s/m2 ]
    REAL(r8), INTENT(in)  :: kv(pcols,pver+1)      ! Vertical diffusion coefficients [ m2/s ]
    REAL(r8), INTENT(in)  :: tmpi(pcols,pver+1)    ! dt*(g/R)**2/dp*pi(k+1)/(.5*(tm(k+1)+tm(k))**2
    REAL(r8), INTENT(in)  :: rpdel(pcols,pver)     ! 1./pdel  (thickness bet interfaces)
    REAL(r8), INTENT(in)  :: ztodt                 ! 2 delta-t [ s ]
    REAL(r8), INTENT(in)  :: cc_top(pcols)         ! Lower diagonal on top interface (for fixed ubc only)

    REAL(r8), INTENT(out) :: ca(pcols,pver)        ! Upper diagonal
    REAL(r8), INTENT(out) :: cc(pcols,pver)        ! Lower diagonal
    REAL(r8), INTENT(out) :: dnom(pcols,pver)      ! 1./(1. + ca(k) + cc(k) - cc(k)*ze(k-1))
    REAL(r8), INTENT(out) :: ze(pcols,pver)        ! Term in tri-diag. matrix system

    ! --------------- !
    ! Local Variables !
    ! --------------- !

    INTEGER :: i                                   ! Longitude index
    INTEGER :: k                                   ! Vertical  index

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !
    DO k = 1, pver
       DO i = 1, ncol
          ca(i,k)   =0.0_r8
          cc(i,k)   =0.0_r8
          dnom(i,k) =0.0_r8
          ze(i,k)   =0.0_r8
       END DO
    END DO

    ! Determine superdiagonal (ca(k)) and subdiagonal (cc(k)) coeffs of the 
    ! tridiagonal diffusion matrix. The diagonal elements  (cb=1+ca+cc) are
    ! a combination of ca and cc; they are not required by the solver.

    DO k = nbot - 1, ntop, -1
       DO i = 1, ncol
          ca(i,k  ) = kv(i,k+1) * tmpi(i,k+1) * rpdel(i,k  )
          cc(i,k+1) = kv(i,k+1) * tmpi(i,k+1) * rpdel(i,k+1)
       END DO
    END DO

    ! The bottom element of the upper diagonal (ca) is zero (element not used).
    ! The subdiagonal (cc) is not needed in the solver.

    DO i = 1, ncol
       ca(i,nbot) = 0._r8
    END DO

    ! Calculate e(k).  This term is 
    ! required in solution of tridiagonal matrix defined by implicit diffusion eqn.

    DO i = 1, ncol
       dnom(i,nbot) = 1._r8/(1._r8 + cc(i,nbot) + ksrf(i)*ztodt*gravit*rpdel(i,nbot))
       ze(i,nbot)   = cc(i,nbot)*dnom(i,nbot)
    END DO

    DO k = nbot - 1, ntop + 1, -1
       DO i = 1, ncol
          dnom(i,k) = 1._r8/(1._r8 + ca(i,k) + cc(i,k) - ca(i,k)*ze(i,k+1))
          ze(i,k)   = cc(i,k)*dnom(i,k)
       END DO
    END DO

    DO i = 1, ncol
       dnom(i,ntop) = 1._r8/(1._r8 + ca(i,ntop) + cc_top(i) - ca(i,ntop)*ze(i,ntop+1))
    END DO

    RETURN
  END SUBROUTINE vd_lu_decomp



  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  SUBROUTINE vd_lu_solve( pcols , pver , ncol , &
       q     , ca   , ze   , dnom , ntop , nbot , cd_top )
    !----------------------------------------------------------------------------------- !
    ! Solve the implicit vertical diffusion equation with zero flux boundary conditions. !
    ! Procedure for solution of the implicit equation follows Richtmyer and              !
    ! Morton (1967,pp 198-200).                                                          !
    !                                                                                    !
    ! The equation solved is                                                             !
    !                                                                                    !  
    !     -ca(k)*q(k+1) + cb(k)*q(k) - cc(k)*q(k-1) = d(k),                              !
    !                                                                                    !
    ! where d(k) is the input profile and q(k) is the output profile                     !
    !                                                                                    ! 
    ! The solution has the form                                                          !
    !                                                                                    !
    !     q(k) = ze(k)*q(k-1) + zf(k)                                                    !
    !                                                                                    !
    !     ze(k) = cc(k) * dnom(k)                                                        !
    !                                                                                    !  
    !     zf(k) = [d(k) + ca(k)*zf(k+1)] * dnom(k)                                       !
    !                                                                                    !
    !     dnom(k) = 1/[cb(k) - ca(k)*ze(k+1)] =  1/[1 + ca(k) + cc(k) - ca(k)*ze(k+1)]   !
    !                                                                                    !
    ! Note that the same routine is used for temperature, momentum and tracers,          !
    ! and that input variables are replaced.                                             !
    ! ---------------------------------------------------------------------------------- ! 

    ! --------------------- !
    ! Input-Output Argument !
    ! --------------------- !

    INTEGER,  INTENT(in)    :: pcols                  ! Number of allocated atmospheric columns
    INTEGER,  INTENT(in)    :: pver                   ! Number of allocated atmospheric levels 
    INTEGER,  INTENT(in)    :: ncol                   ! Number of computed atmospheric columns
    INTEGER,  INTENT(in)    :: ntop                   ! Top level to operate on
    INTEGER,  INTENT(in)    :: nbot                   ! Bottom level to operate on
    REAL(r8), INTENT(in)    :: ca(pcols,pver)         ! -Upper diag coeff.of tri-diag matrix
    REAL(r8), INTENT(in)    :: ze(pcols,pver)         ! Term in tri-diag solution
    REAL(r8), INTENT(in)    :: dnom(pcols,pver)       ! 1./(1. + ca(k) + cc(k) - ca(k)*ze(k+1))
    REAL(r8), INTENT(in)    :: cd_top(pcols)          ! cc_top * ubc value

    REAL(r8), INTENT(inout) :: q(pcols,pver)          ! Constituent field

    ! --------------- !
    ! Local Variables ! 
    ! --------------- !

    REAL(r8)                :: zf(pcols,pver)         ! Term in tri-diag solution
    INTEGER                    i, k                   ! Longitude, vertical indices

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !
    zf=0.0_r8
    ! Calculate zf(k). Terms zf(k) and ze(k) are required in solution of 
    ! tridiagonal matrix defined by implicit diffusion equation.
    ! Note that only levels ntop through nbot need be solved for.

    DO i = 1, ncol
       zf(i,nbot) = q(i,nbot)*dnom(i,nbot)
    END DO

    DO k = nbot - 1, ntop + 1, -1
       DO i = 1, ncol
          zf(i,k) = (q(i,k) + ca(i,k)*zf(i,k+1))*dnom(i,k)
       END DO
    END DO

    ! Include boundary condition on top element

    k = ntop
    DO i = 1, ncol
       zf(i,k) = (q(i,k) + cd_top(i) + ca(i,k)*zf(i,k+1))*dnom(i,k)
    END DO

    ! Perform back substitution

    DO i = 1, ncol
       q(i,ntop) = zf(i,ntop)
    END DO

    DO k = ntop + 1, nbot, +1
       DO i = 1, ncol
          q(i,k) = zf(i,k) + ze(i,k)*q(i,k-1)
       END DO
    END DO

    RETURN
  END SUBROUTINE vd_lu_solve


  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  CHARACTER(128) FUNCTION vdiff_select( fieldlist, name, qindex )
    ! --------------------------------------------------------------------- !
    ! This function sets the field with incoming name as one to be diffused !
    ! --------------------------------------------------------------------- !
    TYPE(vdiff_selector), INTENT(inout)        :: fieldlist
    CHARACTER(*),         INTENT(in)           :: name
    INTEGER,              INTENT(in), OPTIONAL :: qindex

    vdiff_select = ''
    SELECT CASE (name)
    CASE ('u','U')
       fieldlist%fields(1) = .TRUE.
    CASE ('v','V')
       fieldlist%fields(2) = .TRUE.
    CASE ('s','S')
       fieldlist%fields(3) = .TRUE.
    CASE ('q','Q')
       IF( PRESENT(qindex) ) THEN
          fieldlist%fields(3 + qindex) = .TRUE.
       ELSE
          fieldlist%fields(4) = .TRUE.
       ENDIF
    CASE default
       WRITE(vdiff_select,*) 'Bad argument to vdiff_index: ', name
    END SELECT
    RETURN

  END FUNCTION vdiff_select


  TYPE(vdiff_selector) FUNCTION NOT(a)
    ! ------------------------------------------------------------- !
    ! This function extends .not. to operate on type vdiff_selector !
    ! ------------------------------------------------------------- !    
    TYPE(vdiff_selector), INTENT(in)  :: a
    ALLOCATE(not%fields(SIZE(a%fields)))
    not%fields(:) = .NOT. a%fields(:)
  END FUNCTION not


  LOGICAL FUNCTION my_any(a)
    ! -------------------------------------------------- !
    ! This function extends the intrinsic function 'any' ! 
    ! to operate on type vdiff_selector                  ! 
    ! -------------------------------------------------- !
    TYPE(vdiff_selector), INTENT(in) :: a
    my_any = ANY(a%fields)
  END FUNCTION my_any



  LOGICAL FUNCTION diffuse(fieldlist,name,qindex)
    ! ---------------------------------------------------------------------------- !
    ! This function reports whether the field with incoming name is to be diffused !
    ! ---------------------------------------------------------------------------- !
    TYPE(vdiff_selector), INTENT(in)           :: fieldlist
    CHARACTER(*),         INTENT(in)           :: name
    INTEGER,              INTENT(in), OPTIONAL :: qindex

    SELECT CASE (name)
    CASE ('u','U')
       diffuse = fieldlist%fields(1)
    CASE ('v','V')
       diffuse = fieldlist%fields(2)
    CASE ('s','S')
       diffuse = fieldlist%fields(3)
    CASE ('q','Q')
       IF( PRESENT(qindex) ) THEN
          diffuse = fieldlist%fields(3 + qindex)
       ELSE
          diffuse = fieldlist%fields(4)
       ENDIF
    CASE default
       diffuse = .FALSE.
    END SELECT
    RETURN
  END FUNCTION diffuse


  CHARACTER*3 FUNCTION cnst_get_type_byind (ind,pcnst)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Get the type of a constituent 
    ! 
    ! Method: 
    ! <Describe the algorithm(s) used in the routine.> 
    ! <Also include any applicable external references.> 
    ! 
    ! Author:  P. J. Rasch
    ! 
    !-----------------------------Arguments---------------------------------
    !
    INTEGER, INTENT(in)   :: ind    ! global constituent index (in q array)
    INTEGER, INTENT(in)   :: pcnst
    !---------------------------Local workspace-----------------------------
    INTEGER :: m                                   ! tracer index

    !-----------------------------------------------------------------------

    IF (ind.LE.pcnst) THEN
       cnst_get_type_byind = cnst_type(ind)
    ELSE
       ! Unrecognized name
       WRITE(iulog,*) 'CNST_GET_TYPE_BYIND, ind:', ind
       STOP 'call endrun'
    ENDIF


  END FUNCTION cnst_get_type_byind

  !==============================================================================================
  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  SUBROUTINE compute_tms( &
       pcols    , &!integer,  intent(in)  :: pcols                 ! Number of columns dimensioned
       pver     , &!integer,  intent(in)  :: pver                  ! Number of model layers
       ncol     , &!integer,  intent(in)  :: ncol                  ! Number of columns actually used
       u        , &!real(r8), intent(in)  :: u(pcols,pver)         ! Layer mid-point zonal wind [ m/s ]
       v        , &!real(r8), intent(in)  :: v(pcols,pver)         ! Layer mid-point meridional wind [ m/s ]
       t        , &!real(r8), intent(in)  :: t(pcols,pver)         ! Layer mid-point temperature [ K ]
       pmid     , &!real(r8), intent(in)  :: pmid(pcols,pver)      ! Layer mid-point pressure [ Pa ]
       exner    , &!real(r8), intent(in)  :: exner(pcols,pver)     ! Layer mid-point exner function [ no unit ]
       zm       , &!real(r8), intent(in)  :: zm(pcols,pver)        ! Layer mid-point height [ m ]
       sgh      , &!real(r8), intent(in)  :: sgh(pcols)            ! Standard deviation of orography [ m ]
       ksrf     , &!real(r8), intent(out) :: ksrf(pcols)           ! Surface drag coefficient [ kg/s/m2 ]
       taux     , &!real(r8), intent(out) :: taux(pcols)           ! Surface zonal      wind stress [ N/m2 ]
       tauy     , &!real(r8), intent(out) :: tauy(pcols)           ! Surface meridional wind stress [ N/m2 ]
       landfrac   )!real(r8), intent(in)  :: landfrac(pcols)       ! Land fraction [ fraction ]

    !------------------------------------------------------------------------------ !
    ! Turbulent mountain stress parameterization                                    !  
    !                                                                               !
    ! Returns surface drag coefficient and stress associated with subgrid mountains !
    ! For points where the orographic variance is small ( including ocean ),        !
    ! the returned surface drag coefficient and stress is zero.                     !
    !                                                                               !
    ! Lastly arranged : Sungsu Park. Jan. 2010.                                     !
    !------------------------------------------------------------------------------ !

    ! ---------------------- !
    ! Input-Output Arguments ! 
    ! ---------------------- !

    INTEGER,  INTENT(in)  :: pcols                 ! Number of columns dimensioned
    INTEGER,  INTENT(in)  :: pver                  ! Number of model layers
    INTEGER,  INTENT(in)  :: ncol                  ! Number of columns actually used

    REAL(r8), INTENT(in)  :: u(pcols,pver)         ! Layer mid-point zonal wind [ m/s ]
    REAL(r8), INTENT(in)  :: v(pcols,pver)         ! Layer mid-point meridional wind [ m/s ]
    REAL(r8), INTENT(in)  :: t(pcols,pver)         ! Layer mid-point temperature [ K ]
    REAL(r8), INTENT(in)  :: pmid(pcols,pver)      ! Layer mid-point pressure [ Pa ]
    REAL(r8), INTENT(in)  :: exner(pcols,pver)     ! Layer mid-point exner function [ no unit ]
    REAL(r8), INTENT(in)  :: zm(pcols,pver)        ! Layer mid-point height [ m ]
    REAL(r8), INTENT(in)  :: sgh(pcols)            ! Standard deviation of orography [ m ]
    REAL(r8), INTENT(in)  :: landfrac(pcols)       ! Land fraction [ fraction ]

    REAL(r8), INTENT(out) :: ksrf(pcols)           ! Surface drag coefficient [ kg/s/m2 ]
    REAL(r8), INTENT(out) :: taux(pcols)           ! Surface zonal      wind stress [ N/m2 ]
    REAL(r8), INTENT(out) :: tauy(pcols)           ! Surface meridional wind stress [ N/m2 ]

    ! --------------- !
    ! Local Variables !
    ! --------------- !

    INTEGER  :: i                                  ! Loop index
    INTEGER  :: kb, kt                             ! Bottom and top of source region

    REAL(r8) :: horo                               ! Orographic height [ m ]
    REAL(r8) :: z0oro                              ! Orographic z0 for momentum [ m ]
    REAL(r8) :: dv2                                ! (delta v)**2 [ m2/s2 ]
    REAL(r8) :: ri                                 ! Richardson number [ no unit ]
    REAL(r8) :: stabfri                            ! Instability function of Richardson number [ no unit ]
    REAL(r8) :: rho                                ! Density [ kg/m3 ]
    REAL(r8) :: cd                                 ! Drag coefficient [ no unit ]
    REAL(r8) :: vmag                               ! Velocity magnitude [ m /s ]

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !
    horo    = 0._r8
    z0oro   = 0._r8
    dv2     = 0._r8
    ri      = 0._r8
    stabfri = 0._r8
    rho     = 0._r8
    cd      = 0._r8
    vmag    = 0._r8
    DO i = 1, ncol
       ksrf(i) = 0._r8
       taux(i) = 0._r8
       tauy(i) = 0._r8
    END DO
    DO i = 1, ncol

       ! determine subgrid orgraphic height ( mean to peak )

       horo = oroconst * sgh(i)

       ! No mountain stress if horo is too small

       IF( horo < horomin ) THEN

          ksrf(i) = 0._r8
          taux(i) = 0._r8
          tauy(i) = 0._r8

       ELSE

          ! Determine z0m for orography

          z0oro = MIN( z0fac * horo, z0max )

          ! Calculate neutral drag coefficient

          cd = ( karman / LOG( ( zm(i,pver) + z0oro ) / z0oro) )**2

          ! Calculate the Richardson number over the lowest 2 layers

          kt  = pver - 1
          kb  = pver
          dv2 = MAX( ( u(i,kt) - u(i,kb) )**2 + ( v(i,kt) - v(i,kb) )**2, dv2min )

          ! Modification : Below computation of Ri is wrong. Note that 'Exner' function here is
          !                inverse exner function. Here, exner function is not multiplied in
          !                the denominator. Also, we should use moist Ri not dry Ri.
          !                Also, this approach using the two lowest model layers can be potentially
          !                sensitive to the vertical resolution.  
          ! OK. I only modified the part associated with exner function.

          ri  = 2._r8 * gravit * ( t(i,kt) * exner(i,kt) - t(i,kb) * exner(i,kb) ) * ( zm(i,kt) - zm(i,kb) ) &
               / ( ( t(i,kt) * exner(i,kt) + t(i,kb) * exner(i,kb) ) * dv2 )

          ! ri  = 2._r8 * gravit * ( t(i,kt) * exner(i,kt) - t(i,kb) * exner(i,kb) ) * ( zm(i,kt) - zm(i,kb) ) &
          !                      / ( ( t(i,kt) + t(i,kb) ) * dv2 )

          ! Calculate the instability function and modify the neutral drag cofficient.
          ! We should probably follow more elegant approach like Louis et al (1982) or Bretherton and Park (2009) 
          ! but for now we use very crude approach : just 1 for ri < 0, 0 for ri > 1, and linear ramping.

          stabfri = MAX( 0._r8, MIN( 1._r8, 1._r8 - ri ) )
          cd      = cd * stabfri

          ! Compute density, velocity magnitude and stress using bottom level properties

          rho     = pmid(i,pver) / ( rair * t(i,pver) ) 
          vmag    = SQRT( u(i,pver)**2 + v(i,pver)**2 )
          ksrf(i) = rho * cd * vmag * landfrac(i)
          taux(i) = -ksrf(i) * u(i,pver)
          tauy(i) = -ksrf(i) * v(i,pver)

       END IF

    END DO

    RETURN
  END SUBROUTINE compute_tms
  !==============================================================================================

  REAL(r8) FUNCTION estblf( td )
    !
    ! Saturation vapor pressure table lookup
    !
    REAL(r8), INTENT(in) :: td         ! Temperature for saturation lookup
    !
    REAL(r8) :: ee       ! intermediate variable for es look-up
    REAL(r8) :: ai
    INTEGER  :: i
    REAL(r8), PARAMETER:: tmin  = 173.16_r8      ! min temperature (K) for table
    REAL(r8), PARAMETER:: tmax  = 375.16_r8      ! max temperature (K) for table! Maximum temperature entry in table
    REAL(r8) ttrice              ! transition range from es over H2O to es over ice
    REAL(r8), PARAMETER :: trice   =  20.00_r8       ! Transition range from es over range to es over ice
    !
    ee = MAX(MIN(td,tmax),tmin)   ! partial pressure
    i = INT(ee-tmin)+1
    ai = AINT(ee-tmin)
    estblf = (tmin+ai-ee+1.0_r8)* &
         estbl(i)-(tmin+ai-ee)* &
         estbl(i+1)
  END FUNCTION estblf
  !--xl
  !==============================================================================================

  INTEGER FUNCTION fqsatd(t    ,p    ,es    ,qs   ,gam   , len     )
    !----------------------------------------------------------------------- 
    ! Purpose: 
    ! This is merely a function interface vqsatd.
    !------------------------------Arguments--------------------------------
    ! Input arguments
    INTEGER , INTENT(in) :: len       ! vector length
    REAL(r8), INTENT(in) :: t(len)       ! temperature
    REAL(r8), INTENT(in) :: p(len)       ! pressure
    ! Output arguments
    REAL(r8), INTENT(out) :: es(len)   ! saturation vapor pressure
    REAL(r8), INTENT(out) :: qs(len)   ! saturation specific humidity
    REAL(r8), INTENT(out) :: gam(len)  ! (l/cp)*(d(qs)/dt)
    ! Call vqsatd
    es(1:len) =0.0_r8
    qs(1:len) =0.0_r8
    gam(1:len)=0.0_r8
    CALL vqsatd(t       ,p       ,es      ,qs      ,gam  , len     )
    fqsatd = 1
    RETURN
  END FUNCTION fqsatd


  SUBROUTINE vqsatd(t       ,p       ,es      ,qs      ,gam      , &
       len     )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Utility procedure to look up and return saturation vapor pressure from
    ! precomputed table, calculate and return saturation specific humidity
    ! (g/g), and calculate and return gamma (l/cp)*(d(qsat)/dT).  The same
    ! function as qsatd, but operates on vectors of temperature and pressure
    ! 
    ! Method: 
    ! 
    ! Author: J. Hack
    ! 
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    INTEGER, INTENT(in) :: len       ! vector length
    REAL(r8), INTENT(in) :: t(len)       ! temperature
    REAL(r8), INTENT(in) :: p(len)       ! pressure
    !
    ! Output arguments
    !
    REAL(r8), INTENT(out) :: es(len)   ! saturation vapor pressure
    REAL(r8), INTENT(out) :: qs(len)   ! saturation specific humidity
    REAL(r8), INTENT(out) :: gam(len)  ! (l/cp)*(d(qs)/dt)
    !
    !--------------------------Local Variables------------------------------
    !
    LOGICAL lflg   ! true if in temperature transition region
    !
    INTEGER i      ! index for vector calculations
    !
    REAL(r8) omeps     ! 1. - 0.622
    REAL(r8) trinv     ! reciprocal of ttrice (transition range)
    REAL(r8) tc        ! temperature (in degrees C)
    REAL(r8) weight    ! weight for es transition from water to ice
    REAL(r8) hltalt    ! appropriately modified hlat for T derivatives
    !
    REAL(r8) hlatsb    ! hlat weighted in transition region
    REAL(r8) hlatvp    ! hlat modified for t changes above freezing
    REAL(r8) tterm     ! account for d(es)/dT in transition region
    REAL(r8) desdt     ! d(es)/dT
    REAL(r8), PARAMETER :: trice   =  20.00_r8       ! Transition range from es over range to es over ice
    REAL(r8), PARAMETER :: ttrice=trice
    REAL(R8),PARAMETER :: SHR_CONST_TKFRZ   = 273.15_R8       ! freezing T of fresh water          ~ K 
    REAL(r8), PARAMETER :: tmelt       = shr_const_tkfrz      ! Freezing point of water (K)
    REAL(R8),PARAMETER :: SHR_CONST_LATICE  = 3.337e5_R8      ! latent heat of fusion      ~ J/kg
    REAL(R8),PARAMETER :: SHR_CONST_LATVAP  = 2.501e6_R8      ! latent heat of evaporation ~ J/kg
    REAL(r8), PARAMETER ::  hlatf     = shr_const_latice     ! Latent heat of fusion (J/kg)
    REAL(r8), PARAMETER ::  hlatv      = shr_const_latvap     ! Latent heat of vaporization (J/kg)
    REAL(R8),PARAMETER :: SHR_CONST_CPDAIR  = 1.00464e3_R8    ! specific heat of dry air   ~ J/kg/K
    REAL(r8), PARAMETER :: cp       = shr_const_cpdair     ! specific heat of dry air (J/K/kg)
    REAL(r8),PARAMETER           :: rgasv         = shr_const_rgas/shr_const_mwwv     ! Water vapor gas constant ~ J/K/kg
    REAL (KIND=r8), PARAMETER   :: rmwmd  =                0.622e0_r8! fracao molar entre a agua e o ar 

    es(1:len) =0.0_r8
    qs(1:len) =0.0_r8
    gam(1:len)=0.0_r8

    !
    !-----------------------------------------------------------------------
    !
    omeps = 1.0_r8 - epsqs
    DO i=1,len
       es(i) = fpvs2es5(t(i))  !/100.0_r8 !Pa ->mb

       !es(i) = estblf(t(i))
       !
       ! Saturation specific humidity
       !
       qs(i)=rmwmd*es(i)/MAX(((p(i)) - (1.0_r8-rmwmd)*es(i)),qmin(1))
       !qs(i) = epsqs*es(i)/(p(i) - omeps*es(i))
       !
       ! The following check is to avoid the generation of negative
       ! values that can occur in the upper stratosphere and mesosphere
       !
       qs(i) = MIN(1.0_r8,qs(i))
       !
       IF (qs(i) < 0.0_r8) THEN
          qs(i) = 1.0_r8
          es(i) = p(i)
       END IF
    END DO
    !
    ! "generalized" analytic expression for t derivative of es
    ! accurate to within 1 percent for 173.16 < t < 373.16
    !
    trinv = 0.0_r8
    IF ((.NOT. icephs) .OR. (ttrice.EQ.0.0_r8)) go to 10
    trinv = 1.0_r8/ttrice
    DO i=1,len
       !
       ! Weighting of hlat accounts for transition from water to ice
       ! polynomial expression approximates difference between es over
       ! water and es over ice from 0 to -ttrice (C) (min of ttrice is
       ! -40): required for accurate estimate of es derivative in transition
       ! range from ice to water also accounting for change of hlatv with t
       ! above freezing where const slope is given by -2369 j/(kg c) = cpv - cw
       !
       tc     = t(i) - tmelt
       lflg   = (tc >= -ttrice .AND. tc < 0.0_r8)
       weight = MIN(-tc*trinv,1.0_r8)
       hlatsb = hlatv + weight*hlatf
       hlatvp = hlatv - 2369.0_r8*tc
       IF (t(i) < tmelt) THEN
          hltalt = hlatsb
       ELSE
          hltalt = hlatvp
       END IF
       IF (lflg) THEN
          tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3) + tc*(pcf(4) + tc*pcf(5))))
       ELSE
          tterm = 0.0_r8
       END IF
       desdt  = hltalt*es(i)/(rgasv*t(i)*t(i)) + tterm*trinv
       gam(i) = hltalt*qs(i)*p(i)*desdt/(cp*es(i)*(p(i) - omeps*es(i)))
       IF (qs(i) == 1.0_r8) gam(i) = 0.0_r8
    END DO
    RETURN
    !
    ! No icephs or water to ice transition
    !
10  DO i=1,len
       !
       ! Account for change of hlatv with t above freezing where
       ! constant slope is given by -2369 j/(kg c) = cpv - cw
       !
       hlatvp = hlatv - 2369.0_r8*(t(i)-tmelt)
       IF (icephs) THEN
          hlatsb = hlatv + hlatf
       ELSE
          hlatsb = hlatv
       END IF
       IF (t(i) < tmelt) THEN
          hltalt = hlatsb
       ELSE
          hltalt = hlatvp
       END IF
       desdt  = hltalt*es(i)/(rgasv*t(i)*t(i))
       gam(i) = hltalt*qs(i)*p(i)*desdt/(cp*es(i)*(p(i) - omeps*es(i)))
       IF (qs(i) == 1.0_r8) gam(i) = 0.0_r8
    END DO
    !
    RETURN
    !
  END SUBROUTINE vqsatd


  SUBROUTINE aqsat(t       ,p       ,es      ,qs        ,ii      , &
       ILEN    ,kk      ,kstart  ,kend      )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Utility procedure to look up and return saturation vapor pressure from
    ! precomputed table, calculate and return saturation specific humidity
    ! (g/g),for input arrays of temperature and pressure (dimensioned ii,kk)
    ! This routine is useful for evaluating only a selected region in the
    ! vertical.
    ! 
    ! Method: 
    ! <Describe the algorithm(s) used in the routine.> 
    ! <Also include any applicable external references.> 
    ! 
    ! Author: J. Hack
    ! 
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    INTEGER, INTENT(in) :: ii             ! I dimension of arrays t, p, es, qs
    INTEGER, INTENT(in) :: kk             ! K dimension of arrays t, p, es, qs
    INTEGER, INTENT(in) :: ILEN           ! Length of vectors in I direction which
    INTEGER, INTENT(in) :: kstart         ! Starting location in K direction
    INTEGER, INTENT(in) :: kend           ! Ending location in K direction
    REAL(r8), INTENT(in) :: t(ii,kk)          ! Temperature
    REAL(r8), INTENT(in) :: p(ii,kk)          ! Pressure
    !
    ! Output arguments
    !
    REAL(r8), INTENT(out) :: es(ii,kk)         ! Saturation vapor pressure
    REAL(r8), INTENT(out) :: qs(ii,kk)         ! Saturation specific humidity
    !
    !---------------------------Local workspace-----------------------------
    !
    REAL(r8) omeps             ! 1 - 0.622
    INTEGER i, k           ! Indices
    REAL (KIND=r8), PARAMETER   :: rmwmd  =                0.622e0_r8! fracao molar entre a agua e o ar 
    !
    !-----------------------------------------------------------------------
    !
    omeps = 1.0_r8 - epsqs
    DO k=kstart,kend
       DO i=1,ILEN
          es(i,k) = fpvs2es5(t(i,k))  !/100.0_r8 !Pa ->mb
          !es(i,k) = estblf(t(i,k))
          !
          ! Saturation specific humidity
          !
          qs(i,k)=rmwmd*es(i,k)/MAX(((p(i,k)) - (1.0_r8-rmwmd)*es(i,k)),qmin(1))
          !qs(i,k) = epsqs*es(i,k)/(p(i,k) - omeps*es(i,k))
          !
          ! Saturation specific humidity
          !
          !
          ! The following check is to avoid the generation of negative values
          ! that can occur in the upper stratosphere and mesosphere
          !
          qs(i,k) = MIN(1.0_r8,qs(i,k))
          !
          IF (qs(i,k) < 0.0_r8) THEN
             qs(i,k) = 1.0_r8
             es(i,k) = p(i,k)
          END IF
       END DO
    END DO
    !
    RETURN
  END SUBROUTINE aqsat
  SUBROUTINE gffgch(t       ,es      ,itype   )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Computes saturation vapor pressure over water and/or over ice using
    ! Goff & Gratch (1946) relationships. 
    ! <Say what the routine does> 
    ! 
    ! Method: 
    ! T (temperature), and itype are input parameters, while es (saturation
    ! vapor pressure) is an output parameter.  The input parameter itype
    ! serves two purposes: a value of zero indicates that saturation vapor
    ! pressures over water are to be returned (regardless of temperature),
    ! while a value of one indicates that saturation vapor pressures over
    ! ice should be returned when t is less than freezing degrees.  If itype
    ! is negative, its absolute value is interpreted to define a temperature
    ! transition region below freezing in which the returned
    ! saturation vapor pressure is a weighted average of the respective ice
    ! and water value.  That is, in the temperature range 0 => -itype
    ! degrees c, the saturation vapor pressures are assumed to be a weighted
    ! average of the vapor pressure over supercooled water and ice (all
    ! water at 0 c; all ice at -itype c).  Maximum transition range => 40 c
    ! 
    ! Author: J. Hack
    ! 
    !-----------------------------------------------------------------------
    !   use shr_kind_mod, only: r8 => shr_kind_r8
    !   use physconst,    only: tmelt
    !   use abortutils,   only: endrun
    !   use cam_logfile,  only: iulog

    IMPLICIT NONE
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    REAL(r8), INTENT(in) :: t          ! Temperature
    !
    ! Output arguments
    !
    INTEGER, INTENT(inout) :: itype   ! Flag for ice phase and associated transition

    REAL(r8), INTENT(out) :: es         ! Saturation vapor pressure
    !
    !---------------------------Local variables-----------------------------
    !
    REAL(r8) e1         ! Intermediate scratch variable for es over water
    REAL(r8) e2         ! Intermediate scratch variable for es over water
    REAL(r8) eswtr      ! Saturation vapor pressure over water
    REAL(r8) f          ! Intermediate scratch variable for es over water
    REAL(r8) f1         ! Intermediate scratch variable for es over water
    REAL(r8) f2         ! Intermediate scratch variable for es over water
    REAL(r8) f3         ! Intermediate scratch variable for es over water
    REAL(r8) f4         ! Intermediate scratch variable for es over water
    REAL(r8) f5         ! Intermediate scratch variable for es over water
    REAL(r8) ps         ! Reference pressure (mb)
    REAL(r8) t0         ! Reference temperature (freezing point of water)
    REAL(r8) term1      ! Intermediate scratch variable for es over ice
    REAL(r8) term2      ! Intermediate scratch variable for es over ice
    REAL(r8) term3      ! Intermediate scratch variable for es over ice
    REAL(r8) tr         ! Transition range for es over water to es over ice
    REAL(r8) ts         ! Reference temperature (boiling point of water)
    REAL(r8) weight     ! Intermediate scratch variable for es transition
    INTEGER itypo   ! Intermediate scratch variable for holding itype

    INTEGER, PARAMETER :: iulog=0
    !
    !-----------------------------------------------------------------------
    !
    ! Check on whether there is to be a transition region for es
    !
    IF (itype < 0) THEN
       tr    = ABS(REAL(itype,r8))
       itypo = itype
       itype = 1
    ELSE
       tr    = 0.0_r8
       itypo = itype
    END IF
    IF (tr > 40.0_r8) THEN
       WRITE(iulog,900) tr
       STOP 'call endrun (GFFGCH)                ! Abnormal termination'
    END IF
    !
    IF(t < (tmelt - tr) .AND. itype == 1) go to 10
    !
    ! Water
    !
    ps = 1013.246_r8
    ts = 373.16_r8
    e1 = 11.344_r8*(1.0_r8 - t/ts)
    e2 = -3.49149_r8*(ts/t - 1.0_r8)
    f1 = -7.90298_r8*(ts/t - 1.0_r8)
    f2 = 5.02808_r8*LOG10(ts/t)
    f3 = -1.3816_r8*(10.0_r8**e1 - 1.0_r8)/10000000.0_r8
    f4 = 8.1328_r8*(10.0_r8**e2 - 1.0_r8)/1000.0_r8
    f5 = LOG10(ps)
    f  = f1 + f2 + f3 + f4 + f5
    es = (10.0_r8**f)*100.0_r8
    eswtr = es
    !
    IF(t >= tmelt .OR. itype == 0) go to 20
    !
    ! Ice
    !
10  CONTINUE
    t0    = tmelt
    term1 = 2.01889049_r8/(t0/t)
    term2 = 3.56654_r8*LOG(t0/t)
    term3 = 20.947031_r8*(t0/t)
    es    = 575.185606e10_r8*EXP(-(term1 + term2 + term3))
    !
    IF (t < (tmelt - tr)) go to 20
    !
    ! Weighted transition between water and ice
    !
    weight = MIN((tmelt - t)/tr,1.0_r8)
    es = weight*es + (1.0_r8 - weight)*eswtr
    !
20  CONTINUE
    itype = itypo
    RETURN
    !
900 FORMAT('GFFGCH: FATAL ERROR ******************************',/, &
         'TRANSITION RANGE FOR WATER TO ICE SATURATION VAPOR', &
         ' PRESSURE, TR, EXCEEDS MAXIMUM ALLOWABLE VALUE OF', &
         ' 40.0 DEGREES C',/, ' TR = ',f7.2)
    !
  END SUBROUTINE gffgch

  !
  !  Finalize_Pbl_UniversityWashington
  !
  SUBROUTINE Finalize_Pbl_UniversityWashington()
    IMPLICIT NONE

  END SUBROUTINE Finalize_Pbl_UniversityWashington

END MODULE Pbl_UniversityWashington

!PROGRAM MAIN
!  USE Pbl_UniversityWashington, ONLY :Init_Pbl_UniversityWashington, &
!                                      Finalize_Pbl_UniversityWashington
!  IMPLICIT NONE

 
! CALL Init()
! CALL Run()
! CALL Finalize()
 
 
!CONTAINS
  !
  !  Init
  !
!SUBROUTINE Init()
! IMPLICIT NONE
! INTEGER, PARAMETER :: pver =28
! INTEGER, PARAMETER :: pcnst=1
! REAL(KIND=8)       :: sig(pver)
! CALL Init_Pbl_UniversityWashington(pver,pcnst,pcnst,sig)
!END SUBROUTINE Init 
  !
  !  Run
  !
!SUBROUTINE Run()
! IMPLICIT NONE

!END SUBROUTINE Run 
  !
  !  Finalize
  !
!SUBROUTINE Finalize()
! IMPLICIT NONE
! CALL Finalize_Pbl_UniversityWashington()
!END SUBROUTINE Finalize 

!END PROGRAM MAIN
