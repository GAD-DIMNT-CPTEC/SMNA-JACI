!
!  $Author: pkubota $
!  $Date: 2006/11/13 12:56:50 $
!  $Revision: 1.4 $
!
MODULE Sfc_SeaFlux_COLA_Model

  USE Constants, ONLY :     &
       cp,            &
       hl,            &
       grav,          &
       tice,          &
       stefan,        &
       z0ice,        &
       gasr,          &
       i8  ,&      
       r8


  IMPLICIT NONE
SAVE


  PRIVATE
  PUBLIC :: InitSeaIceFlux_COLA_Model
  PUBLIC :: SeaIceFlux_COLA_Model
  INTEGER      ,    PARAMETER :: nfprt =0 
  REAL(KIND=r8),    PARAMETER :: eps   =   0.608_r8
  REAL(KIND=r8),    PARAMETER :: shrmin=   1.0e-5_r8
  REAL(KIND=r8),    PARAMETER :: facl  =   0.1_r8
  INTEGER, PARAMETER :: nitr  =   2
  REAL(KIND=r8),    PARAMETER :: gkm0  =   1.00_r8
  REAL(KIND=r8),    PARAMETER :: gkh0  =   0.10_r8
  REAL(KIND=r8),    PARAMETER :: gkm1  = 300.0_r8
  REAL(KIND=r8),    PARAMETER :: gkh1  = 300.0_r8
  REAL(KIND=r8),    PARAMETER :: vk0   =   0.4_r8
  INTEGER, PARAMETER :: kmean =   1

  REAL(KIND=r8),    PARAMETER :: epsq2=0.2_r8
  REAL(KIND=r8),    PARAMETER :: FH=1.01_r8
  REAL(KIND=r8),    PARAMETER :: a1    =   0.92_r8
  REAL(KIND=r8),    PARAMETER :: a2    =   0.74_r8
  REAL(KIND=r8),    PARAMETER :: b1    =  16.6_r8
  REAL(KIND=r8),    PARAMETER :: b2    =  10.1_r8
  REAL(KIND=r8),    PARAMETER :: c1    =   0.08_r8
  REAL(KIND=r8),    PARAMETER :: deltx =   0.0_r8


  REAL(KIND=r8), PARAMETER :: x        = 1.0_r8
  REAL(KIND=r8), PARAMETER :: xx       = x*x
  REAL(KIND=r8), PARAMETER :: g        = 1.0_r8
  REAL(KIND=r8), PARAMETER :: gravi    = 9.81_r8
  REAL(KIND=r8), PARAMETER :: agrav    = 1.0_r8/gravi
  REAL(KIND=r8), PARAMETER :: rgas     = 287.0_r8
  REAL(KIND=r8), PARAMETER :: akwnmb   = 2.5e-05_r8
  REAL(KIND=r8), PARAMETER :: lstar    = 1.0_r8/akwnmb
  REAL(KIND=r8), PARAMETER :: gocp     = gravi/1005.0_r8
  INTEGER         :: nbase
  INTEGER         :: nthin,nthinp
!  REAL(KIND=r8)   :: rbyg

  REAL(KIND=r8)               :: alfa
  REAL(KIND=r8)               :: beta
  REAL(KIND=r8)               :: gama
  REAL(KIND=r8)               :: dela
  REAL(KIND=r8)               :: r1
  REAL(KIND=r8)               :: r2
  REAL(KIND=r8)               :: r3
  REAL(KIND=r8)               :: r4
  REAL(KIND=r8)               :: s1
  REAL(KIND=r8)               :: s2
  REAL(KIND=r8)               :: rfc
  !REAL(KIND=r8), ALLOCATABLE  :: sigkiv(:)
  !REAL(KIND=r8), ALLOCATABLE  :: sigr(:)
  !REAL(KIND=r8), ALLOCATABLE  :: sigriv(:)
  !REAL(KIND=r8), ALLOCATABLE  :: a0(:)
  !REAL(KIND=r8), ALLOCATABLE  :: b0(:)
  !REAL(KIND=r8), ALLOCATABLE  :: con0(:)
  !REAL(KIND=r8), ALLOCATABLE  :: con1(:)
  !REAL(KIND=r8), ALLOCATABLE  :: con2(:)
  !REAL(KIND=r8), ALLOCATABLE  :: t0(:)
  !REAL(KIND=r8), ALLOCATABLE  :: t1(:)
  !REAL(KIND=r8)               :: c0pbl
  REAL(KIND=r8)     :: akappa
  REAL(KIND=r8)     :: gbyr

CONTAINS


  SUBROUTINE InitSeaIceFlux_COLA_Model()
    !INTEGER, INTENT(IN) :: kmax
    !REAL(KIND=r8),    INTENT(IN) :: sig(kmax)
    !REAL(KIND=r8),    INTENT(IN) :: delsig(kmax)
    !REAL(KIND=r8),    INTENT(IN) :: sigml(kmax+1)
    !INTEGER  :: k
    REAL(KIND=r8)     :: gam1
    REAL(KIND=r8)     :: gam2
    !REAL(KIND=r8)     :: sigk(kmax)
    !ALLOCATE(sigkiv(kmax))
    !ALLOCATE(sigr  (kmax)) 
    !ALLOCATE(sigriv(kmax))
    !ALLOCATE(a0    (kmax))
    !ALLOCATE(b0    (kmax))
    !ALLOCATE(con0  (kmax))
    !ALLOCATE(con1  (kmax))
    !ALLOCATE(con2  (kmax))
    !ALLOCATE(t0    (kmax))
    !ALLOCATE(t1    (kmax))
!    rbyg  =gasr/grav*delsig(1)*0.5_r8

    gam1=1.0_r8/3.0_r8-2.0_r8*a1/b1
    gam2=(b2+6.0_r8*a1)/b1
    alfa=b1*(gam1-c1)+3.0_r8*(a2+2.0_r8*a1)
    beta=b1*(gam1-c1)
    gama=a2/a1*(b1*(gam1+gam2)-3.0_r8*a1)
    dela=a2/a1* b1* gam1
    r1  =0.5_r8*gama/alfa
    r2  =    beta/gama
    r3  =2.0_r8*(2.0_r8*alfa*dela-gama*beta)/(gama*gama)
    r4  =r2*r2
    s1  =3.0_r8*a2* gam1
    s2  =3.0_r8*a2*(gam1+gam2)
    !     
    !     critical flux richardson number
    !     
    rfc =s1/s2
    akappa=gasr/cp
    !
!    DO k = 1, kmax
!       sigk  (k)=sig(k)**akappa
 !      sigkiv(k)=1.0_r8/sigk(k)
!       con0  (k)=gasr*delsig(k)/(grav*sig(k))
!    END DO
!    a0    (kmax)=0.0_r8
!    b0    (   1)=0.0_r8
!    sigr  (kmax)=0.0_r8
!    sigriv(   1)=0.0_r8
    gbyr        =(grav/gasr)**2!(m/sec**2)/(J/(Kg*K))=(m/sec**2)/((Kg*(m/sec**2)*m)/(Kg*K))
    !(m/sec**2)/((Kg*(m**2/sec**2))/(Kg*K))
    !(m/sec**2)/(m**2/sec**2*K)=K**2/m**2
!    DO k = 1, kmax-1
 !      con1  (   k)=grav*sigml(k+1)/(gasr*(sig(k)-sig(k+1)))
 !      con2  (   k)=grav*con1(k)
 !      con1  (   k)=con1(k)*con1(k)
 !      t0    (   k)=(sig(k+1)-sigml(k+1))/(sig(k+1)-sig(k))
 !      t1    (   k)=(sigml(k+1)-sig(k  ))/(sig(k+1)-sig(k))
 !      sigr  (   k)=sigk(k)*sigkiv(k+1)
 !      sigriv( k+1)=sigk(k+1)*sigkiv(k)
 !      a0    (k   )=gbyr*sigml(k+1)**2/(delsig(k  )*(sig(k)-sig(k+1)))
 !      b0    (k+1 )=gbyr*sigml(k+1)**2/(delsig(k+1)*(sig(k)-sig(k+1)))
 !   END DO
    !c0pbl=grav/(gasr*delsig(1))
  END SUBROUTINE InitSeaIceFlux_COLA_Model

  !----------------------------------------------------------------------
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !----------------------------------------------------------------------
  SUBROUTINE SeaIceFlux_COLA_Model ( &
                                ! Run Flags
       atmpbl    , &
                                ! Model information
       ncols     ,kmax           ,mskant    ,mlsi      , &
                                ! Model Geometry,&
       sinclt,&
                                ! Time info
       delt   ,&
                                ! Atmospheric fields
       gu    ,gv    ,gt    ,gq    , &
       prsi  ,prsl  ,phii  ,phil  ,&
       gps   ,qsurf ,tsurf ,tsea  , &
       speedm, &
                                ! SSIB: Total radiation absorbed at ground
       rnet , &
                                ! Turbulence fields
       Mmlen ,cu    ,ct    , &
       rmi   ,rhi   ,z0    ,zorl  ,ustar  ,tmtx  , &
       qmtx  ,umtx  ,&
                                ! Heat and Vapor Flux
       sens  ,evap  ,     &
                                ! Solving a system of linear equations by Gauss elimination 
       c3    ,b30   ,c0    , &
       b33   ,c4    ,b40   , &
       b44     )


    INTEGER      ,    INTENT(in   ) :: atmpbl
    INTEGER      ,    INTENT(in   ) :: ncols
    INTEGER      ,    INTENT(in   ) :: kmax
    INTEGER(KIND=i8), INTENT(in   ) :: mskant(ncols)
    INTEGER(KIND=i8), INTENT(inout) :: mlsi  (ncols)

    REAL(KIND=r8),    INTENT(in   ) :: sinclt(ncols)

    REAL(KIND=r8),    INTENT(in   ) :: delt

    REAL(KIND=r8),    INTENT(in   ) :: gu    (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: gv    (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: gt    (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: gq    (ncols,kmax)
    REAL(KIND=r8),    INTENT(IN   ) :: prsi  (nCols,kMax+1)  !     prsi     - real, pressure at layer interfaces             ix,levs+1  Pa
    REAL(KIND=r8),    INTENT(IN   ) :: prsl  (nCols,kMax)    !     prsl     - real, mean layer presure                       ix,levs   Pa
    REAL(KIND=r8),    INTENT(IN   ) :: phii  (nCols,kMax+1)  !===>  PHIH(K+1)  INPUT GEOPOTENTIAL @ EDGES  IN MKS units (m)
    REAL(KIND=r8),    INTENT(IN   ) :: phil  (nCols,kMax)    !===>  PHIL(K)    INPUT GEOPOTENTIAL @ LAYERS IN MKS units (m)
    REAL(KIND=r8),    INTENT(in   ) :: gps   (ncols) !Pa
    REAL(KIND=r8),    INTENT(in   ) :: qsurf (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: tsurf (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: tsea  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: speedm(ncols)

    REAL(KIND=r8),    INTENT(in   ) :: rnet  (ncols)

    REAL(KIND=r8),    INTENT(in   ) :: Mmlen (ncols)
    REAL(KIND=r8),    INTENT(inout) :: cu    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ct    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rmi   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rhi   (ncols)

    REAL(KIND=r8),    INTENT(inout) :: z0    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: zorl  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ustar (ncols)

    REAL(KIND=r8),    INTENT(inout) :: tmtx  (ncols,kmax,3)
    REAL(KIND=r8),    INTENT(inout) :: qmtx  (ncols,kmax,3)
    REAL(KIND=r8),    INTENT(inout) :: umtx  (ncols,kmax,4)
    REAL(KIND=r8),    INTENT(out  ) :: sens  (ncols)
    REAL(KIND=r8),    INTENT(out  ) :: evap  (ncols)
    REAL(KIND=r8),    INTENT(out  ) :: c3    (ncols)
    REAL(KIND=r8),    INTENT(out  ) :: b30   (ncols)
    REAL(KIND=r8),    INTENT(out  ) :: c0    (ncols)
    REAL(KIND=r8),    INTENT(out  ) :: b33   (ncols)
    REAL(KIND=r8),    INTENT(out  ) :: c4    (ncols)
    REAL(KIND=r8),    INTENT(out  ) :: b40   (ncols)
    REAL(KIND=r8),    INTENT(out  ) :: b44   (ncols)

    REAL(KIND=r8) :: gl0    (ncols)
    REAL(KIND=r8) :: rho    (ncols)
    REAL(KIND=r8) :: ah     (ncols)
    REAL(KIND=r8) :: al     (ncols)
    REAL(KIND=r8) :: dqg0   (ncols)
    REAL(KIND=r8) :: b00    (ncols)
    REAL(KIND=r8) :: b03    (ncols)
    REAL(KIND=r8) :: b04    (ncols)

    REAL(KIND=r8) :: dzm   (ncols)
    REAL(KIND=r8) :: thm   (ncols)
    REAL(KIND=r8) :: thvgm (ncols)
    REAL(KIND=r8) :: tha   (ncols)
    REAL(KIND=r8) :: cuni  (ncols)
    REAL(KIND=r8) :: cui   (ncols)
    REAL(KIND=r8) :: ctni  (ncols)
    REAL(KIND=r8) :: rib   (ncols)
    REAL(KIND=r8) :: cti   (ncols)

    REAL(KIND=r8) :: hscap
!    REAL(KIND=r8) :: gbycp
!    REAL(KIND=r8) :: gbyhl
    REAL(KIND=r8) :: dti
    REAL(KIND=r8) :: gice
!    REAL(KIND=r8) :: gr100
!    REAL(KIND=r8) :: sl1kap
    REAL(KIND=r8) :: st4

    INTEGER       :: i
    !real(r8)         :: excess(ncols)     ! Excess downward sfc latent heat flux
    !*JPB REAL(KIND=r8), PARAMETER :: dd=0.05_r8
    REAL(KIND=r8), PARAMETER :: dd=3.0_r8 ! Total depth of the ice slab (m), Using ECMWF value
    REAL(KIND=r8), PARAMETER :: tice=271.16_r8
    REAL(KIND=r8), PARAMETER :: dice=2.0_r8
    REAL(KIND=r8), PARAMETER :: hice=2.03_r8
    REAL(KIND=r8), PARAMETER :: rhoice=920.0_r8 ! Mean ice density (kg/m3)
    REAL(KIND=r8), PARAMETER :: cice=2093.0_r8  ! Heat Capacity of Ice (J/Kg)
    sens  =0.0_r8
    evap  =0.0_r8
    c3    =0.0_r8
    b30   =0.0_r8
    c0    =0.0_r8
    b33   =0.0_r8
    c4    =0.0_r8
    b40   =0.0_r8
    b44   =0.0_r8
    gl0   = Mmlen
    dzm   =0.0_r8
    thm   =0.0_r8
    thvgm =0.0_r8
    tha   =0.0_r8
    cuni  =0.0_r8
    cui   =0.0_r8
    cti   =0.0_r8

    ctni  =0.0_r8
    rib   =0.0_r8

    hscap =rhoice*cice*dd/delt
!    gr100 =gasr*0.01_r8
!    gbycp =grav/(cp*delsig(1)*100.0_r8 *sigki(1))
!    gbyhl =grav/(hl*delsig(1)*100.0_r8 )
    gice  =hice/dice ! 2.03_r8/2.0_r8
!    sl1kap=sigki(1)! sigki ! Fator de conversao de temperatura potencial sigki (k)=1.0e0_r8 / EXP(rk*LOG(sl(k)))
    st4   =stefan*4.0_r8
    dti   =1.0_r8 /delt
    CALL vntlt1 ( &
         ! Model information
         ncols , kMax  ,mskant(1:ncols), &
         mlsi(1:ncols),&
         ! Model Geometry,&
         sinclt(1:ncols),dzm(1:ncols)   , &
         ! Time info
         ! Atmospheric fields
         gu  (1:ncols,1:kMax)    ,gv (1:ncols,1:kMax)    ,gt   (1:ncols,1:kMax)    ,tsurf(1:ncols) , &
         prsi(1:ncols,1:kMax+1),prsl(1:ncols,1:kMax),phii(1:nCols,1:kMax+1),phil(1:nCols,1:kMax),&
         tsea(1:ncols)           ,thm(1:ncols)           ,thvgm(1:ncols)           ,speedm(1:ncols), &
         tha(1:ncols)            , &
         ! Turbulence fields Flux  
         cuni(1:ncols)   ,cui(1:ncols)    ,cu  (1:ncols)    ,ctni(1:ncols)   , &
         cti (1:ncols)   ,ct (1:ncols)    ,rmi (1:ncols)    ,rhi (1:ncols)   , &
         rib (1:ncols)   ,z0 (1:ncols)    ,zorl(1:ncols)    ,ustar(1:ncols)   )

    IF(atmpbl /= 1)CALL TubulenceMY2( &
         ! Model information
         ncols                 ,kmax                  , &
         ! Model Geometry,&
         sinclt(1:ncols)       , &
         ! Time info
         delt                      , &
         ! Atmospheric fields
         gu (1:ncols,1:kMax)   ,gv(1:ncols,1:kMax)    ,gt (1:ncols,1:kMax)   ,gq(1:ncols,1:kMax)        , &
         prsi(1:ncols,1:kMax+1),prsl(1:ncols,1:kMax),phii(1:nCols,1:kMax+1),phil(1:nCols,1:kMax),&
         ! Turbulence fields
         tmtx(1:ncols,1:kmax,1:3)  ,qmtx(1:ncols,1:kmax,1:3)  ,umtx(1:ncols,1:kmax,1:4) , &
         gl0(1:ncols) )

    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN
          rho  (i)=gps(i)/(gasr*gt(i,1))

          ah   (i)=(grav/(cp*((prsi(i,1)/prsi(i,1)) - (prsi(i,2)/prsi(i,1)) )* ((prsi(i,1)/(prsi(i,2)))**(gasr/cp))  ))  /  prsi(i,1)
!          ah   (i)=gbycp/gps(i)
          al   (i)=(grav/(hl*((prsi(i,1)/prsi(i,1)) - (prsi(i,2)/prsi(i,1)))  ))  /  prsi(i,1)

!          al   (i)=gbyhl/gps(i)
          dqg0 (i)=0.622_r8 *EXP(30.25353_r8 -5418.0_r8 /tsurf(i))/(tsurf(i)*tsurf(i)*(gps(i)/100.0_r8))
          IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) < tice+0.01_r8) THEN
             !
             ! Solution of sea ice
             !
             !              dTs
             ! rho*ci*Dzi* ---- = Rdown - rho*cp*ck* (Ts - Tl*Sigl) - rho*L*ck*(qs -ql) - cs*(Ts - Ti) - Sigma*Ts^4
             !              dt
             !
             !              ddTs     Rdown     rho*cp*ck* d(Ts - Tl*Sigl)   rho*L*ck*d(qs -ql)     cs*d(Ts - Ti)     Sigma*dTs^4
             ! rho*ci*Dzi* ------- = -----  - -------------------------- - --------------------- - --------------- - ----------
             !              dTs*dt   dTs                   dTs                       dTs                dTs             dTs

             !              d  dTs   Rdown                  rho*L*ck*d(qs -ql)
             ! rho*ci*Dzi* ------- = -----  - rho*cp*ck  - -------------------- -  cs - 4Sigma*dTs^3
             !              dt dTs    dTs                        dTs          
             !
             b00(i)=   hscap+cp*rho(i)*rhi(i)+hl*rho(i)*rhi(i)*dqg0(i) + gice+st4*tsurf(i)**3
             b03(i)=        -cp*rho(i)*rhi(i)*((prsi(i,1)/(prsi(i,2)))**(gasr/cp))
             b04(i)=-hl*rho(i)*rhi(i)
             ! b00 + b10 + b20 + b30 + b40    c0
             ! b01 +             b33          c3     
             ! b02                     b44    c4
             ! b03
             ! b04 
             ! Right side of eq.41 section III.A 
             ! COLA Physics Description Manual
             c0 (i)=rnet(i) -cp*rho(i)*rhi(i)*(tsurf(i)-((prsi(i,1)/(prsi(i,2)))**(gasr/cp))*gt(i,1)) &
                  -hl*rho(i)*rhi(i)*(qsurf(i)-       gq(i,1)) &
                  -gice*(tsurf(i)-tice)-stefan*tsurf(i)**4
             b30(i)=               -ah (i)*cp*rho(i)*rhi(i)
             b33(i)=tmtx(i,1,2)*dti-b30(i)*          ((prsi(i,1)/(prsi(i,2)))**(gasr/cp))
             c3 (i)=tmtx(i,1,3)    -b30(i)*(tsurf(i)-((prsi(i,1)/(prsi(i,2)))**(gasr/cp))*gt(i,1))
             b40(i)=               -al(i)*hl*rho(i)*rhi(i)* dqg0 (i)
             b44(i)=qmtx(i,1,2)*dti+al(i)*hl*rho(i)*rhi(i)
             c4 (i)=qmtx(i,1,3)    + &
                  al(i)*hl*rho(i)*rhi(i)*(qsurf(i)-gq(i,1))
             b00(i)=b00(i)-b30(i)*b03(i)/b33(i)-b40(i)*b04(i)/b44(i)
             c0 (i)=c0 (i)-c3 (i)*b03(i)/b33(i)-c4 (i)*b04(i)/b44(i)
             c0 (i)=c0 (i)/b00(i)
             IF(atmpbl == 1) THEN
                tsurf(i)=tsurf(i)+c0(i)
             END IF
             sens(i)=rho(i)*cp*  (tsurf(i)   -  gt(i,1)*(((prsi(i,1)/(prsi(i,2)))**(gasr/cp)))   )*rhi(i)
             evap(i)=rho(i)*hl*  (qsurf(i)   -  gq(i,1)            )*rhi(i)

          ELSE IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) > tice+0.01_r8) THEN
             !
             ! Solution of sea water
             !
             rho  (i) = gps(i)/((gasr)*gt(i,1))
             sens (i) = rho(i)*cp*  (tsurf(i)   -  gt(i,1)*(((prsi(i,1)/(prsi(i,2)))**(gasr/cp)))   )*rhi(i)
             evap (i) = rho(i)*hl*  (qsurf(i)   -  gq(i,1)            )*rhi(i)

          END IF
       END IF
    END DO

  END SUBROUTINE SeaIceFlux_COLA_Model



  ! vntlt1 :performs ventilation mass flux, based on deardorff, mwr, 1972?.



  SUBROUTINE vntlt1 ( &
                                ! Model information
       ncols , kMax  ,mskant,mlsi, &
                                ! Model Geometry,&
       sinclt,dzm   ,  &
                                ! Time info
                                ! Atmospheric fields
       gu    ,gv    ,gt    ,tsurf ,&
       prsi  ,prsl  ,phii  ,phil  ,&
       tsea  ,thm   ,thvgm ,speedm,&
       tha   , &
                                ! Turbulence fields Flux  
       cuni  ,cui   ,cu    ,ctni  , &
       cti   ,ct    ,rmi   ,rhi   , &
       rib   ,z0    ,zorl  ,ustar   )
    !
    !==========================================================================
    !
    ! _________________
    ! MODEL INFORMATION
    !
    ! nCols.........Number of grid points on a gaussian latitude circle
    ! kmax..........Number of sigma levels
    !
    ! __________________
    ! ATMOSPHERIC FIELDS
    !
    ! gu............(zonal      velocity)*sin(colat)
    ! gv............(meridional velocity)*sin(colat)
    ! gt............temperature
    ! tsurf.........earth's surface temperature used for radiation
    !               for the first time step when ground temperature is not yet
    !               computed (this is done by subr.tsinit ),
    ! tsea..........effective surface radiative temperature ( tgeff )
    ! speedm........speedm(i)=SQRT(gu(i)**2+gv(i)**2)*sincli, where
    !               sincli=1.0 /sinclt
    ! tha...........tha   (i)= tsurf(i)
    ! thm...........thm   (i)= gt(i)*sigki(1)
    ! dzm...........dzm   (i)=gt(i)*rbyg
    !               rbyg  =gasr/grav*delsig(1)*0.5
    ! thvgm.........thvgm (i)= tha(i)-thm(i)

    ! z0ice.........Roughness length of ice
    ! rmi...........rmi   (i)=cu(i)*ustar(i), where
    !               cu is friction  transfer coefficients
    !               ustar is surface friction velocity  (m/s)
    ! rhi...........rhi   (i)=ct(i)*ustar(i), where
    !               ct is heat transfer coefficients.
    !               ustar is surface friction velocity  (m/s)
    ! zorl..........zorl (i)= 100.0 *zgrav*speedm(i)*rhi(i)
    !               zgrav =0.032 /grav
    !
    ! ______________
    ! MODEL GEOMETRY
    !
    ! delsig........Difference between sigma levels
    ! sinclt........sinclt=SIN(colrad(latitu))
    ! sigki ........sigki (k)=1.0e0/EXP(akappa*LOG(sig(k))),  where "sig"
    !               sigma coordinate at middle of layer and akappa=gasr/cp
    ! ________________________
    ! TURBULENCE FIELDS FLUX  
    ! cuni..........neutral friction transfer  coefficients.
    ! cui...........cui   (i)=cuni(i)*EXP( aa-SQRT(aa*aa+bb*f))
    !               cui   (i)=cuni(i)*EXP(-tt+SQRT(tt*tt+ss*f))
    ! cu............Friction  transfer coefficients.
    ! ctni..........neutral heat transfer coefficients.
    ! cti...........cti   (i)=ctni(i)*EXP( qq-SQRT(qq*qq+rr*g))
    !               cti   (i)=cui (i)
    ! ct............heat transfer coefficients.
    ! rib...........bulk richardson number.
    ! z0............Roughness length
    ! ustar.........surface friction velocity  (m/s)
    !
    ! ___________________
    ! PHYSICAL CONSTANTS
    ! 
    ! gasr..........gas constant of dry air        (j/kg/k)
    ! grav..........grav   gravity constant        (m/s**2)
    !==========================================================================
    !    
    INTEGER      ,    INTENT(in   ) :: ncols
    INTEGER      ,    INTENT(in   ) :: kMax
    REAL(KIND=r8),    INTENT(in   ) :: sinclt(ncols)
    REAL(KIND=r8),    INTENT(inout) :: rmi   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rhi   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: gu    (ncols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: gv    (ncols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: gt    (ncols,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: prsi  (nCols,kMax+1)  !     prsi     - real, pressure at layer interfaces             ix,levs+1  Pa
    REAL(KIND=r8),    INTENT(IN   ) :: prsl  (nCols,kMax)    !     prsl     - real, mean layer presure                       ix,levs   Pa
    REAL(KIND=r8),    INTENT(IN   ) :: phii  (nCols,kMax+1)  !===>  PHIH(K+1)  INPUT GEOPOTENTIAL @ EDGES  IN MKS units (m)
    REAL(KIND=r8),    INTENT(IN   ) :: phil  (nCols,kMax)    !===>  PHIL(K)    INPUT GEOPOTENTIAL @ LAYERS IN MKS units (m)
    REAL(KIND=r8),    INTENT(in   ) :: tsurf (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: tsea  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: zorl  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: cuni  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: cui   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: cu    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ctni  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: cti   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ct    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: speedm(ncols)
    REAL(KIND=r8),    INTENT(inout) :: tha   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: thm   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: dzm   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: thvgm (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rib   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: z0    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ustar (ncols)
    INTEGER(KIND=i8), INTENT(in   ) :: mskant(ncols)
    INTEGER(KIND=i8), INTENT(inout) :: mlsi  (ncols)

    REAL(KIND=r8),        PARAMETER :: vkrmn=0.40_r8 
    REAL(KIND=r8),        PARAMETER :: ribc=3.05_r8
    REAL(KIND=r8),        PARAMETER :: aa=1.2270_r8
    REAL(KIND=r8),        PARAMETER :: bb=1.2642_r8
    REAL(KIND=r8),        PARAMETER :: tt=1.8900_r8
    REAL(KIND=r8),        PARAMETER :: ss=5.0519_r8
    REAL(KIND=r8),        PARAMETER :: ee=1.2743_r8
    REAL(KIND=r8),        PARAMETER :: ff=3.4805_r8
    REAL(KIND=r8),        PARAMETER :: gg=0.87581_r8
    REAL(KIND=r8),        PARAMETER :: hh=-1.5630_r8
    REAL(KIND=r8),        PARAMETER :: pp=10.815_r8
    REAL(KIND=r8),        PARAMETER :: qq=1.3462_r8
    REAL(KIND=r8),        PARAMETER :: rr=1.8380_r8
    REAL(KIND=r8),        PARAMETER :: CHARNOCK=.014_r8
    REAL(KIND=r8)                   :: sincli(ncols)
    REAL(KIND=r8)                   :: f
    REAL(KIND=r8)                   :: g
    INTEGER                :: i

    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN
          z0(i)=0.001_r8
          IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) >= 271.17_r8) THEN
             ! 
             ! Solution of sea water
             !
             z0(i)=0.01_r8*zorl(i)
          ELSE IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) < 271.17_r8) THEN
             !
             ! Solution of sea ice
             !
             z0(i)=z0ice
             mlsi(i)=2_i8
          END IF
          sincli(i)=1.0_r8 /sinclt(i)
       END IF
    END DO

    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN
          IF (tsea(i) <= 0.0_r8) THEN
             speedm(i)=SQRT(gu(i,1)**2+gv(i,1)**2)*sincli(i)
             speedm(i)=MAX(2.0_r8 ,speedm(i))

             !    rbyg  =gasr/grav*delsig(1)*0.5_r8

             !dzm   (i)=gt(i,1)*rbyg
             dzm   (i)=0.5_r8*MAX((phii(i,2) - phii(i,1)),2.0_r8)
             cuni(i)=LOG(dzm(i)/z0(i))/vkrmn*gg+hh
             ctni(i)=cuni(i)
             !
             !     stability branch based on bulk richardson number.
             !
             !thm   (i)= gt(i,1)*sigki(1)
             thm   (i)= gt(i,1)*((prsi(i,1)/(prsi(i,2)))**(gasr/cp))
             tha   (i)= tsurf(i)
             thvgm (i)= tha(i)-thm(i)
             rib   (i)=-thvgm(i)*grav*dzm(i)/ (thm(i)*speedm(i)**2)
             rib   (i)=MAX(-1.25_r8 ,rib(i))
             rib   (i)=MIN( 1.25_r8 ,rib(i))
             IF (rib(i) < 0.0_r8) THEN
                f        =LOG(1.0_r8-ee*rib(i))
                cui   (i)=cuni(i)*EXP( aa-SQRT(aa*aa+bb*f))
                g        =LOG(1.0_r8-ff*rib(i))
                cti   (i)=ctni(i)*EXP( qq-SQRT(qq*qq+rr*g))
             ELSE
                f        =LOG(1.0_r8+pp*rib(i))
                cui   (i)=cuni(i)*EXP(-tt+SQRT(tt*tt+ss*f))
                cti   (i)=cui (i)
             END IF
             cu    (i)=1.0_r8/cui(i)
             ct    (i)=1.0_r8/cti(i)
             !
             !     surface friction velocity and ventilation mass flux
             !
             ustar (i)=speedm(i)*cu(i)
             rmi   (i)=cu(i)*ustar(i)
             rhi   (i)=ct(i)*ustar(i)
          END IF
       END IF
    END DO
  END SUBROUTINE vntlt1

  !----------------------------------------------------------------------
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !----------------------------------------------------------------------
  SUBROUTINE TubulenceMY2 ( &
       ! Model information
       ncols ,kmax  , &
       ! Model Geometry,&
       sinclt ,&
       ! Time info
       delt  , &
       ! Atmospheric fields
       gu     ,gv    ,gt    ,gq    , &
       prsi,prsl,phii,phil,&
       ! Turbulence fields
       gmt   ,gmq   ,gmu   , &
       gl0)
    !==========================================================================
    !
    ! _________________
    ! MODEL INFORMATION
    !
    ! nCols.........Number of grid points on a gaussian latitude circle
    ! kmax..........Number of sigma levels
    !
    ! __________________
    ! ATMOSPHERIC FIELDS
    !
    ! gu............(zonal      velocity)*sin(colat)
    ! gv............(meridional velocity)*sin(colat)
    ! gt............temperature
    ! gq............specific humidity
    !
    ! ______________
    ! MODEL GEOMETRY
    !
    ! delsig........Difference between sigma levels
    ! sinclt........sinclt=SIN(colrad(latitu))
    !
    !_____________
    ! TIME INFO
    !...............delt   time interval
    !
    !___________________
    ! TURBULENCE FIELDS
    !
    !..gl0    maximum mixing length l0 in blackerdar's formula
    !                                  l=k0*z/(1+k0*z/l0)
    !..gmt    temperature related matrix
    !         gmt(i,k,1)*d(gt(i,k-1))/dt+gmt(i,k,2)*d(gt(i,k))/dt=gmt(i,k,3)
    !         gmt(i,1,1)=0.
    !     gmt(*,*,1)...dimensionless
    !     gmt(*,*,2)...dimensionless
    !     gmt(*,*,3)...deg/sec
    !..gmq    specific humidity related matrix
    !         gmq(i,k,1)*d(gq(i,k-1))/dt+gmq(i,k,2)*d(gq(i,k))/dt=gmq(i,k,3)
    !         gmq(i,1,1)=0.
    !     gmq(*,*,1)...dimensionless
    !     gmq(*,*,2)...dimensionless
    !     gmq(*,*,3)...kg/kg/sec
    !..gmu    wind related matrix
    !         gmu(i,k,1)*d(gu(i,k-1))/dt+gmu(i,k,2)*d(gu(i,k))/dt=gmu(i,k,3)
    !         gmu(i,k,1)*d(gv(i,k-1))/dt+gmu(i,k,2)*d(gv(i,k))/dt=gmu(i,k,4)
    !         gmu(i,1,1)=0.
    !     gmu(*,*,1)...dimensionless
    !     gmu(*,*,2)...dimensionless
    !     gmu(*,*,3)...m/sec**2
    !     gmu(*,*,4)...m/sec**2
    !
    !
    ! TubulenceMY2 :performs momentum, water vapour and sensible heat diffusion
    !         on planetary boundary layer.
    ! A. Vertical diffusion - Mellor-Yamada closure scheme
    ! The effects of mixing of heat, momentum and moisture by small scale
    ! turbulence is represented by vertical diffusion in the COLA GCM. 
    ! The mixing coefficients are calculated according to the "level 2.0" 
    ! closure scheme of Mellor and Yamada (1982). This method assumes a
    ! local balance between production and dissipation of turbulent kinetic
    ! energy. There are no explicit prognostic variables to describe 
    ! the planetary boundary layer (PBL); instead, the entire atmosphere 
    ! is represented in discrete layers which may or may not be part 
    ! of the PBL. The prognostic equations for atmospheric temperatures 
    ! and moisture are then coupled to the SSIB equations for the 
    ! ground surface and canopy, and the system of coupled equations is 
    ! solved simultaneously with vertical diffusion of heat, moisture, 
    ! and momentum as given by the Mellor and Yamada (1982) scheme.
    !
    ! In the turbulence closure method, each prognostic variable is expressed 
    ! as a sum of a large scale (resolved) part and a turbulent (sub-grid)
    ! scale part. The vertical fluxes that are expressed as quadratic terms
    ! in the turbulent quantities are assumed to be represented by vertical
    ! diffusion down the gradient of the large scale quantities, e.g.
    !-----------------------------------------------------------------------
    !
    !           input values
    !-----------------------------------------------------------------------
    !..imx.......Number of grid points on a gaussian latitude circle (ncols+ 2)
    !..ncols......Number of grid points on a gaussian latitude circle
    !..kmax......Number of grid points at vertcal  
    !..cp........Specific heat of air           (j/kg/k) 
    !..gasr......Gas constant of dry air        (j/kg/k) 
    !..grav......gravity constant               (m/s**2) 
    !..gu,gv,gt,gq is at time level t-dt
    !
    !..gu    (zonal      velocity)*sin(colat)
    !..gv    (meridional velocity)*sin(colat)
    !..gt    temperature
    !..gq    specific humidity
    !..fsen  sensible heat flux in w/m**2
    !..flat  latent   heat fulx in w/m**2
    !..fmom  momentum flux      in n/m**2
    !
    !..delsig     k=2  ****gu,gv,gt,gq,gyu,gyv,gtd,gqd,sig*** } delsig(2)
    !             k=3/2----si,ric,rf,km,kh,b,l -----------
    !             k=1  ****gu,gv,gt,gq,gyu,gyv,gtd,gqd,sig*** } delsig(1)
    !             k=1/2----si ----------------------------
    !..delt   time interval
    !..gl0    maximum mixing length l0 in blackerdar's formula
    !                                  l=k0*z/(1+k0*z/l0)
    !..csqiv  1./sin(colat)**2
    !..ncols   number of grid points on a gaussian latitude circle
    !..kpbl   number of layers pbl process is included( for u v,t )
    !..kqpbl  number of layers pbl process is included( for q     )
    !-----------------------------------------------------------------------
    !           work arrays
    !-----------------------------------------------------------------------
    !..gwrk
    !..gld
    !..gln
    !-----------------------------------------------------------------------
    !           output values
    !-----------------------------------------------------------------------
    !..gmt    temperature related matrix
    !         gmt(i,k,1)*d(gt(i,k-1))/dt+gmt(i,k,2)*d(gt(i,k))/dt=gmt(i,k,3)
    !         gmt(i,1,1)=0.
    !     gmt(*,*,1)...dimensionless
    !     gmt(*,*,2)...dimensionless
    !     gmt(*,*,3)...deg/sec
    !..gmq    specific humidity related matrix
    !         gmq(i,k,1)*d(gq(i,k-1))/dt+gmq(i,k,2)*d(gq(i,k))/dt=gmq(i,k,3)
    !         gmq(i,1,1)=0.
    !     gmq(*,*,1)...dimensionless
    !     gmq(*,*,2)...dimensionless
    !     gmq(*,*,3)...kg/kg/sec
    !..gmu    wind related matrix
    !         gmu(i,k,1)*d(gu(i,k-1))/dt+gmu(i,k,2)*d(gu(i,k))/dt=gmu(i,k,3)
    !         gmu(i,k,1)*d(gv(i,k-1))/dt+gmu(i,k,2)*d(gv(i,k))/dt=gmu(i,k,4)
    !         gmu(i,1,1)=0.
    !     gmu(*,*,1)...dimensionless
    !     gmu(*,*,2)...dimensionless
    !     gmu(*,*,3)...m/sec**2
    !     gmu(*,*,4)...m/sec**2
    !..gl0    maximum mixing length l0 in blackerdar's formula
    !         this is retained as a first guess for next time step
    !
    INTEGER, INTENT(in   ) :: ncols
    INTEGER, INTENT(in   ) :: kmax
    REAL(KIND=r8),    INTENT(in   ) :: gu    (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: gv    (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: gt    (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: gq    (ncols,kmax)
    REAL(KIND=r8),    INTENT(IN   ) :: prsi  (nCols,kMax+1)  !     prsi     - real, pressure at layer interfaces             ix,levs+1  Pa
    REAL(KIND=r8),    INTENT(IN   ) :: prsl  (nCols,kMax)    !     prsl     - real, mean layer presure                       ix,levs   Pa
    REAL(KIND=r8),    INTENT(IN   ) :: phii  (nCols,kMax+1)  !===>  PHIH(K+1)  INPUT GEOPOTENTIAL @ EDGES  IN MKS units (m)
    REAL(KIND=r8),    INTENT(IN   ) :: phil  (nCols,kMax)    !===>  PHIL(K)    INPUT GEOPOTENTIAL @ LAYERS IN MKS units (m)
    REAL(KIND=r8),    INTENT(in   ) :: delt
    REAL(KIND=r8),    INTENT(in   ) :: sinclt(ncols)
    REAL(KIND=r8),    INTENT(OUT  ) :: gmt   (ncols,kmax,3)
    REAL(KIND=r8),    INTENT(OUT  ) :: gmq   (ncols,kmax,3)
    REAL(KIND=r8),    INTENT(OUT  ) :: gmu   (ncols,kmax,4)
    REAL(KIND=r8),    INTENT(INOUT) :: gl0   (ncols)
    REAL(KIND=r8) :: PBL_CoefKm(ncols, kmax)
    REAL(KIND=r8) :: PBL_CoefKh(ncols, kmax)
    REAL(KIND=r8) :: a   (ncols,kmax)
    REAL(KIND=r8) :: b   (ncols,kmax)
    REAL(KIND=r8) :: Pbl_NRich(ncols,kmax)
    REAL(KIND=r8) :: Pbl_ATemp(ncols,kmax)
    REAL(KIND=r8) :: Pbl_ITemp(ncols,kmax)
    REAL(KIND=r8) :: Pbl_Shear(ncols,kmax)  !square of vertical wind shear
    REAL(KIND=r8) :: Pbl_Sqrtw(ncols,kmax)  !sqrt(w) or b/l
    REAL(KIND=r8) :: Pbl_MixLgh(ncols,kmax) !mixing length
    REAL(KIND=r8) :: Pbl_BRich(ncols,kmax)
    REAL(KIND=r8) :: Pbl_SmBar(ncols,kmax)
    REAL(KIND=r8) :: Pbl_ShBar(ncols,kmax)
    REAL(KIND=r8) :: Pbl_KmMixl(ncols,kmax)
    REAL(KIND=r8) :: Pbl_KhMixl(ncols,kmax)
    REAL(KIND=r8) :: Pbl_PotTep(ncols,kmax)
    REAL(KIND=r8) :: Pbl_Stabil(ncols,kmax)
    REAL(KIND=r8) :: Pbl_EddEner(ncols,kmax)
    REAL(KIND=r8) :: Pbl_HgtLyI(ncols,kmax)
    REAL(KIND=r8) :: gld (ncols)
    REAL(KIND=r8) :: gln (ncols)
    REAL(KIND=r8) :: csqiv(ncols)
    REAL(KIND=r8) :: Pbl_KM(ncols,kmax)
    REAL(KIND=r8) :: Pbl_KH(ncols,kmax)

    !     
    !     eps   ; virtual temperature correction factor
    !     shrmin; squre of minimum wind shear, this is in order to avoid
    !     large richardson number   (1.0_r8/sec**2)
    !     facl  ; appears in l0 computation
    !     nitr  ; number of iteration computing l0
    !     gkm0  ; minimum value of eddy diffusion coefficient for momentm
    !     (m/sec**2)
    !     gkh0  ; minimum value of eddy diffusion coefficient for sensible heat
    !     (m/sec**2)
    !     gkm1  ; maximum value of eddy diffusion coefficient for momentm
    !     (m/sec**2)
    !     gkh1  ; maximum value of eddy diffusion coefficient for sensible heat
    !     (m/sec**2)
    !     vk0   ; von-karman constant
    !

    INTEGER :: k 
    INTEGER :: i 
    INTEGER :: itr 
    INTEGER :: icnt(ncols)
    REAL(KIND=r8)    :: s1ms2g 
    REAL(KIND=r8)    :: x 
    REAL(KIND=r8)    :: y 
    REAL(KIND=r8)    :: fac 
    REAL(KIND=r8)    :: rfx
    !     
    !     ichk(kmax) is flag to vectorize.
    !     
    INTEGER :: ichk(ncols,kmax)
    REAL(KIND=r8)    :: aa0(ncols,kmax)
    REAL(KIND=r8)    :: bb0(ncols,kmax)
    REAL(KIND=r8)    :: tt0(ncols,kmax)
    REAL(KIND=r8)    :: tt1(ncols,kmax)
    REAL(KIND=r8)    :: sigr  (ncols,kmax)
    REAL(KIND=r8)    :: sigriv(ncols,kmax)
    REAL(KIND=r8)    :: con0  (ncols,kmax)
    REAL(KIND=r8)    :: con1  (ncols,kmax)
    REAL(KIND=r8)    :: con2  (ncols,kmax)

!    REAL(KIND=r8)    :: c(ncols)
    REAL(KIND=r8)    :: twodt
    REAL(KIND=r8)    :: twodti
    !     
    !     constants concerning mixing length scaling
    !     (mellor & yamada '82 rev.geophys.space sci. p851-874)
    !     
    IF (delt == deltx) THEN
       WRITE(UNIT=nfprt, FMT="(' ERROR - delt == 0  in TubulenceMY2' )")
       STOP "**(ymbpl0)"
    END IF
    DO i=1,ncols,1 
       ! sinclt........sinclt=SIN(colrad(latitu))
       csqiv(i)   = 1.0e0_r8/(sinclt(i))**2
    END DO

    twodt =delt
    twodti=1.0_r8/twodt
    DO k = 1, kmax
       DO i=1,nCols
          ! con0  (k)=gasr*delsig(k)/(grav*sig(k))
           con0  (i,k)=gasr*((prsi(i,k)/prsi(i,1)) - (prsi(i,k+1)/prsi(i,1)))/(grav*(prsl(i,k)/prsi(i,1)))
       END DO
    END DO

    DO i=1,nCols
       aa0    (i,kmax)=0.0_r8
       tt0   (i,kmax)=0.0_r8
       tt1   (i,kmax)=0.0_r8
       bb0    (i,   1)=0.0_r8
       sigr  (i,kmax)=0.0_r8
       sigriv(i,   1)=0.0_r8

    END DO

    DO k = 1, kmax-1
       DO i=1,nCols
       ! con0  (k)=gasr*delsig(k)/(grav*sig(k))

       ! con0  (i,k)=gasr*((prsi(i,k)/prsi(i,1)) - (prsi(i,k+1)/prsi(i,1)))/(grav*(prsl(i,k)/prsi(i,1)))

        !con1  (   k)=grav*sigml(k+1)/(gasr*(sig(k)-sig(k+1)))

        con1  (i,k) = grav*(prsi(i,k+1)/prsi(i,1))/(gasr*(((prsl(i,k)/prsi(i,1)) - (prsl(i,k+1)/prsi(i,1)))))

        !con2  (   k)=grav*con1(k)
        con2  (i,k)=grav*con1(i,k)

        !con1  (   k)=con1(k)*con1(k)
        con1  (i,k)=con1(i,k)*con1(i,k)

!       t0    (   k)=(sig(k+1)-sigml(k+1))/(sig(k+1)-sig(k))

        tt0    (i,   k)=((prsl(i,k+1)/prsi(i,1)) - (prsi  (i,k+1)/prsi  (i,1))) / &
                        ((prsl(i,k+1)/prsi(i,1)) - (prsl(i,k)/prsi(i,1)))


 !       t1    (   k)=(sigml(k+1)-sig(k  ))/(sig(k+1)-sig(k))
        tt1    (i,   k)=((prsi  (i,k+1)/prsi  (i,1)) - (prsl(i,k)/prsi(i,1)))/&
                        ((prsl  (i,k+1)/prsi  (i,1)) - (prsl(i,k)/prsi(i,1)))

       !sigr  (   k)=sigk(k)*sigkiv(k+1)
        sigr  (i,k)= ((prsl(i,k)/prsi(i,1))**akappa)*  (1.0_r8/((prsl(i,k+1)/prsi(i,1))**akappa))
       !sigriv( k+1)=sigk(k+1)*sigkiv(k)
        sigriv(i,k+1)= ((prsl(i,k+1)/prsi(i,1))**akappa )  *  (1.0_r8/((prsl(i,k)/prsi(i,1))**akappa))

       !a0    (k   )=gbyr*sigml(k+1)**2/(delsig(k  )*(sig(k)-sig(k+1)))
       !REAL(KIND=r8),    INTENT(IN   ) :: prsi  (nCols,kMax+1)  !     prsi     - real, pressure at layer interfaces             ix,levs+1  Pa
       !REAL(KIND=r8),    INTENT(IN   ) :: prsl  (nCols,kMax)    !     prsl     - real, mean layer presure                       ix,levs   Pa

       aa0    (i,k   )=gbyr*((prsi  (i,k+1)/prsi  (i,1))**2)/(((prsi(i,k)/prsi(i,1)) - (prsi(i,k+1)/prsi(i,1)))*((prsl(i,k)/prsi(i,1)) - (prsl(i,k+1)/prsi(i,1))))

       !bb0    (i,k+1 )=gbyr*sigml(k+1)**2/(delsig(k+1)*(sig(k)-sig(k+1)))

       bb0    (i,k+1 )=gbyr*((prsi  (i,k+1)/prsi  (i,1))**2)/(((prsi(i,k+1)/prsi(i,1)) - (prsi(i,k+2)/prsi(i,1)))*((prsl(i,k)/prsi(i,1)) - (prsl(i,k+1)/prsi(i,1))))

       END DO
    END DO

    DO k = 1, kmax
       DO i=1,nCols

       !                 --                  --
       !    1      1    | g       si(k+1)      |     con1(k)
       !  ------ =--- * |--- * ----------------| = ----------
       !    DZ     T    | R     sl(k) -sl(k+1) |       T
       !                 --                  --
       !                -- --
       !    DA      d  |     |
       !  ------ =---- | W'A'|
       !    DT      dZ |     |
       !                -- --
       !
       !                          ----
       !    DA      d         d  |    |
       !  ------ =---- * K * --- | A  |
       !    DT      dZ        dZ |    |
       !                          ----
       !                     
       !    DA         d      dA 
       !  ------ =K * ---- * ---- 
       !    DT         dZ     dZ 
       !                     
       !
       !
       !a0(k)  =gbyr*sigml(k+1)**2/(delsig(k  )*(sig(k)-sig(k+1)))
       !
       !                grav **2
       !              --------    * si(k+1)**2
       !                gasr**2 
       !a0(k)  = ---------------------------------
       !          ((si(k)-si(k+1))) * (sl(k)-sl(k+1)))
       !
       !         --  -- 2                             -- -- 2
       !        |   g  |       si(k+1)**2            |  1  | 
       !a0(k)  =| -----| * --------------------    = | --- |
       !        |   R  |     ((si(k)-si(k+1))  )     |  dZ | 
       !         --  --                               -- --
       !                         --   --      -- -- 2            --                   -- 2
       !                        |       |    |  1  |            |      m       kg * K   |  
       a(i,k)=twodt*aa0(i,k) !  |  2*Dt |  * | --- |    ==> s * |   ------- * --------  | 
       !                        |       |    |  dZ |            |     s**2       J      |  
       !                         --   --      -- --              --                   --
       !J = F*DX = kg m/s**2 *m = kg * m**2/s**2
       !                   --   --      -- -- 2            --                       -- 2
       !                  |       |    |  1  |            |    m         kg * K *s**2 |    K**2 * s
       !a(k)=twodt*a0(k)! |  2*Dt |  * | --- |    ==> s * | ------- * --------------- | = -----------
       !                  |       |    |  dZ |            |   s**2     kg * m**2      |    m**2  
       !                   --   --      -- --              --                       --
       !                   --   --      -- -- 2  
       !                  |       |    |  1  |         K**2 * s
       !a(k)=twodt*a0(k)! |  2*Dt |  * | --- |    ==> -----------
       !                  |       |    |  dZ |           m**2  
       !                   --   --      -- --    

       !
       !                   gbyr*sigml(k+1)**2
       !    b0(k) = -----------------------------------------------
       !               (((si(k+1)-si(k+1+1))  ) *(sig(k)-sig(k+1)))
 
       
       b(i,k)=twodt*bb0(i,k)! s * K**2/m**2 
    END DO
    END DO
    DO i=1,nCols
       !c0pbl=grav/(gasr*((prsi(k)/prsi(1)) - (prsi(k+1)/prsi(1))))

!       c(i)=twodt*grav/(gasr*((prsi(i,k)/prsi(i,1)) - (prsi(i,k+1)/prsi(i,1))))
    END DO
    !     
    !     Pbl_PotTep(3)   virtual potential temperature
    !     Pbl_ATemp (7)   temperature at the interface of two adjacent layers
    !
    !             k=2  ****gu,gv,gt,gq,gyu,gyv,gtd,gqd,sl*** } delsig(2)
    !             k=3/2----si,ric,rf,km,kh,b,l -----------
    !             k=1  ****gu,gv,gt,gq,gyu,gyv,gtd,gqd,sl*** } delsig(1)
    !             k=1/2----si ----------------------------
    !
    !                      sl(k+1)-si(k+1)
    !       t0    (   k)=-------------------
    !                      sl(k+1)-sl(k  )
    !
    !                      si(k+1)-sl(k  )
    !       t1    (   k)=-------------------
    !                      sl(k+1)-sl(k  )
    !
    !
    !                      sl(k+1)-si(k+1)              si(k+1)-sl(k  )
    !   Pbl_ATemp(i,k)  =-------------------*gt(i,k) + ------------------*gt(i,k+1)
    !                      sl(k+1)-sl(k  )              sl(k+1)-sl(k  )
    !
    !
    !
    !
    DO k = 1, kmax-1
       DO i = 1, ncols
          Pbl_ATemp(i,k)  = tt0(i,k)*gt(i,k) + tt1(i,k)*gt(i,k+1)
       END DO
    END DO
    !                      si(k+1)-sl(k  )
    !       t1    (   k)=-------------------
    !                      sl(k+1)-sl(k  )
    !
    !                          1
    !      sigkiv(k   )=--------------------
    !                      (sl(k)**akappa)
    !
    !                       --                       --  
    !                      | grav        si(k+1)       |
    !   con2  (   k)=grav* |------- * ---------------- |  
    !                      | gasr      sl(k)-sl(k+1))  |
    !                       --                       -- 
    !
    !                --                       --  2
    !               |   grav        si(k+1)     |
    !   con1  (k) = |------- * ---------------- |         
    !               |  gasr      sl(k)-sl(k+1)) |
    !                --                       -- 
    !
    ! sl.........sigma value at midpoint of                  gasr/cp
    !                                         each layer : (k=287/1005)=R/cp
    !
    !                                                                     1
    !                                             +-                   + ---
    !                                             !     k+1         k+1!  k
    !                                             !si(l)   - si(l+1)   !
    !                                     sl(l) = !--------------------!
    !                                             !(k+1) (si(l)-si(l+1)!
    !             --   --  -(R/Cp)                +-                  -+     
    !            |  P    |
    !Thetav = Tv |-------|
    !            |  P0   | 
    !             --   --
    !
    !             --   --  -(R/Cp)        --   --  -(R/Cp)
    !            |  P    |               |       |
    !sigkiv(k) = |-------|            == |sl(k)  |
    !            |  P0   |               |       | 
    !             --   --                 --   --
    !  --   --      --   -- 
    ! |  P    |    |       |
    ! |-------| =  |sl(k)  |
    ! |  P0   |    |       |
    !  --   --      --   --
    !  --          --      --           -- 
    ! |  P(k)-P(k+1) |    |               |
    ! |--------------| =  |sl(k)-sl(k+1)  |
    ! |      P0      |    |               |
    !  --          --      --           --




    DO k = 1, kmax
       DO i = 1, ncols
          Pbl_PotTep(i,k)=(1.0_r8/((prsl(i,k)/prsi(i,1))**akappa))*gt(i,k)*(1.0_r8+eps*gq(i,k))
       END DO
    END DO
    !     
    !     Pbl_Stabil(2)   stability
    !     Pbl_Shear (6)   square of vertical wind shear
    !  
    gln =1.0e-5_r8 
    DO k = 1, kmax-1
       DO i = 1, ncols
          !
          !
          !                 g         D (Theta)
          !Pbl_Stabil(2) =------- * -------------
          !                Theta      D (Z)
          !
          !
          ! P =rho*R*T and P = rho*g*Z
          !
          !                           P
          ! DP = rho*g*DZ and rho = ----
          !                          R*T
          !        P
          ! DP = ----*g*DZ
          !       R*T
          !
          !        R*T
          ! DZ = ------*DP
          !        g*P
          !
          !    1       g*P       1   
          !  ------ = ------ * ------
          !    DZ      R*T       DP  
          !
          !    T       g         P   
          !  ------ = ------ * ------
          !    DZ      R         DP  
          !
          !  --   --      --   -- 
          ! |  P    |    |       |
          ! |-------| =  |sl(k)  |
          ! |  P0   |    |       |
          !  --   --      --   --
          !
          !    T       g            si(k+1)   
          !  ------ = ------ * ----------------- =con1(k)
          !    DZ      R         sl(k) -sl(k+1)  
          !
          !                 --                  --
          !    1      1    | g       si(k+1)      |     con1(k)
          !  ------ =--- * |--- * ----------------| = ----------
          !    DZ     T    | R     sl(k) -sl(k+1) |       T
          !                 --                  --
          !                                  --                  --
          !                    g       1    | g       si(k+1)      |   
          ! Pbl_Stabil(2) = ------- * --- * |--- * ----------------|* D(Theta)
          !                  Theta     T    | R     sl(k) -sl(k+1) | 
          !                                  --                  --
          !
          !                   con2(k)*(Pbl_PotTep(i,k+1)-Pbl_PotTep(i,k))
          !Pbl_Stabil(i,k)=-----------------------------------------------------------------
          !                  ((t0(k)*Pbl_PotTep(i,k)+t1(k)*Pbl_PotTep(i,k+1))*Pbl_ATemp(i,k))
          !
          Pbl_Stabil(i,k)=con2(i,k)*(Pbl_PotTep(i,k+1)-Pbl_PotTep(i,k)) &
               /((tt0(i,k)*Pbl_PotTep(i,k)+tt1(i,k)*Pbl_PotTep(i,k+1))*Pbl_ATemp(i,k))
          Pbl_ITemp(i,k)=1.0_r8/(Pbl_ATemp(i,k)*Pbl_ATemp(i,k))
          !
          !                  --  -- 2    --  -- 2
          !                 |  dU  |    |  dV  |
          !Pbl_Shear(i,k) = | -----|  + | -----| 
          !                 |  DZ  |    |  DZ  |
          !                  --  --      --  --
          !
          !                  --  -- 2    --     -- 2
          !                 |   1  |    |         |
          !Pbl_Shear(i,k) = | -----|  * | DU + DV | 
          !                 |  DZ  |    |         |
          !                  --  --      --     --
          !
          !                         --                  -- 2    --     -- 2
          !                 1      | g       si(k+1)      |    |         |          1.0e0_r8
          !Pbl_Shear(i,k) =----- * |--- * ----------------|  * | Du + Dv |   * ---------------------
          !                 T*T    | R     sl(k) -sl(k+1) |    |         |       SIN(colrad(i))**2
          !                         --                  --      --     --
          !                  1.0e0_r8
          ! csqiv (i) = ---------------------
          !               SIN(colrad(i))**2
          !
          Pbl_Shear(i,k)=(con1(i,k)*csqiv(i))*Pbl_ITemp(i,k)*((gu(i,k)-gu(i,k+1))**2+(gv(i,k)-gv(i,k+1))**2)
          Pbl_Shear(i,k)=MAX(gln(i),Pbl_Shear(i,k))          
       END DO
    END DO
    !     
    !     Pbl_BRich(4)        richardson number
    !     Pbl_NRich(5)   flux richardson number
    !     Pbl_NRich(8)   flux richardson number
    !     
    DO k = 1,(kmax-1)
       DO i = 1, ncols
          !                g         D (Theta)
          !              ------- * -------------
          !               Theta      D (Z)
          ! Pbl_BRich = -----------------------------
          !                 --  -- 2    --  -- 2
          !                |  dU  |    |  dV  |
          !                | -----|  + | -----| 
          !                |  DZ  |    |  DZ  |
          !                 --  --      --  --
          !
          Pbl_BRich(i,k)=Pbl_Stabil(i,k)/Pbl_Shear(i,k)
          !r1 = 0.5_r8*gama/alfa
          Pbl_NRich(i,k)= r1*(Pbl_BRich(i,k)+r2 &
               -SQRT(Pbl_BRich(i,k)*(Pbl_BRich(i,k)-r3)+r4))
          Pbl_NRich(i,k)=MIN(rfc,Pbl_NRich(i,k))
          Pbl_NRich(i,k)=Pbl_NRich(i,k)
          !
          !    Pbl_SmBar and Pbl_ShBar are momentum flux and heat flux 
          !    stability parameters, respectively
          !     
          !     Pbl_ShBar(3)   shbar
          !     Pbl_SmBar(4)   smbar
          !     
          !     eliminate negative value for s1-s2*gwrk(i,1,5):
          !     gwrk(i,1,5) is s1/s2 under some circumstances
          !     which makes this expression zero.  machine roundoff
          !     can produce an unphysical negative value in this case.
          !     it is used as sqrt argument in later loop.
          !     s1  =3.0_r8*a2* gam1
          !     s2  =3.0_r8*a2*(gam1+gam2)
          !
          ! s1-s2 = 3.0*a2* gam1 - 3.0*a2*(gam1+gam2)
          ! s1-s2 = 3.0*a2* (gam1 - (gam1+gam2))
          !     
          s1ms2g=s1-s2*Pbl_NRich(i,k)
          IF (ABS(s1ms2g) < 1.0e-10_r8) s1ms2g=0.0_r8
          !
          !                     s1ms2g
          !Pbl_ShBar(i,k)=---------------------------
          !                 (1.0_r8-Pbl_NRich(i,k))
          !
          ! a1    =   0.92_r8
          ! a2    =   0.74_r8
          ! b1    =  16.6_r8
          ! b2    =  10.1_r8
          !
          !         1.0               a1
          !gam1 = -------  - 2.0  * ------
          !         3.0               b1
          !
          !       b2          a1
          !gam2=(----) + (6*------)
          !       b1          b1
          !
          !gam2=(b2  + 6.0*a1)
          !     --------------
          !          b1
          !                 --                   --
          !                | 1.0               a1  |
          !s1 = 3.0 * a2 * |------  - 2.0  * ------|
          !                | 3.0               b1  |
          !                 --                   --
          !s2  =3.0_r8*a2*(gam1+gam2)
          !
          !                   s1-s2*Pbl_NRich(i,k)
          !Pbl_ShBar(i,k)=---------------------------
          !                 (1.0_r8-Pbl_NRich(i,k))
          !
          !                        ((gam1 - (gam1+gam2)))*Pbl_NRich(i,k)
          !Pbl_ShBar(i,k)=3.0*a2* ---------------------------------------- =Sm
          !                           (1.0_r8-Pbl_NRich(i,k))

          Pbl_ShBar(i,k)=s1ms2g/(1.0_r8-Pbl_NRich(i,k))
          !     
          !     gwrk(i,1,3)=(s1-s2*gwrk(i,1,5))/(1.0_r8-gwrk(i,1,5))
          !     end of  negative sqrt argument trap
          ! c1    =   0.08_r8
          ! alfa=b1*(gam1-c1)+3.0_r8*(a2+2.0_r8*a1)
          ! beta=b1*(gam1-c1)
          !        --                              --
          !       | a2                               |
          ! gama=-|---- *(b1*(gam1+gam2)-3.0_r8*a1)  |
          !       | a1                               |
          !        --                              --
          !       --          --
          !      | a2           |
          ! dela=|---- b1* gam1 |
          !      | a1           |
          !       --          --
          !                                     ((b1*(gam1-c1))-(b1*(gam1-c1)+3.0*(a2+2.0*a1))*Pbl_NRich(i,k))
          ! Pbl_SmBar(i,k)=Pbl_ShBar(i,k) * ------------------------------------------------------------------
          !                                     --          --     --                           --
          !                                    | a2           |   | a2                            |
          !                                    |---- b1* gam1 | - |---- *(b1*(gam1+gam2)-3.0*a1)  |*Pbl_NRich(i,k))
          !                                    | a1           |   | a1                            |
          !                                     --          --     --                           --
          !
          !                                     ((b1*(gam1-c1)) - (b1*(gam1-c1) + 3.0*a2 + 6.0*a1))*Pbl_NRich(i,k))
          ! Pbl_SmBar(i,k)=Pbl_ShBar(i,k) * --------------------------------------------------------------
          !                                     --  --     --                             --
          !                                    | a2   |   |                                 |
          !                                    |----  | * |b1* gam1 -(b1*(gam1+gam2)-3.0*a1)|*Pbl_NRich(i,k))
          !                                    | a1   |   |                                 |
          !                                     --  --     --                              --

          !   
          Pbl_SmBar(i,k)=Pbl_ShBar(i,k)*(beta-alfa*Pbl_NRich(i,k))/ &
               (dela-gama*Pbl_NRich(i,k))
          !     
          !     
          !     Pbl_Sqrtw(5)   sqrt(w) or b/l
          !     Pbl_KmMixl(4)   km/l**2
          !     Pbl_KhMixl(3)   kh/l**2
          ! 
          !     The ratio of SH to SM is equal to the ratio of the turbulent
          !     flux Richardson number to the bulk (large scale) Richardson
          !     number. 
          !
          !     u^2 =  (1-2*gam1)q^2
          !     
          !     SQRT(B1*SM*(1 - Ri)*Shear)
          !             GH    SM
          !     Ri = - ---- =-----Rf
          !             GM    SH
          !             GH
          !     GM = - ----
          !             Ri   
          !                --                  --  
          !               | --  -- 2    --  -- 2 |
          !          l^2  ||  dU  |    |  dV  |  |
          !     GM =-----*|| -----|  + | -----|  |
          !          q^2  ||  DZ  |    |  DZ  |  |
          !               | --  --               |
          !                --                  --  
          ! 
          !                          --      -- 2 
          !             l^2         |  dTHETA  |  
          !     GH = - -----*beta*g*| ---------|  
          !             q^2         |  DZ      |  
          !                          --      --   

          !
          !          --                       --  1/2
          !         |      --  -- 2    --  -- 2 |
          !     q   | 1   |  dU  |    |  dV  |  |
          !    ----=|----*| -----|  + | -----|  |
          !     l   | GM  |  DZ  |    |  DZ  |  |
          !         |      --  --               |
          !          --                       --   
          !                                               
          !  1  
          ! ---- = B1*Sm*(1 - Rf)
          !  GM 
          !
          Pbl_Sqrtw (i,k)=SQRT(b1*Pbl_SmBar(i,k)*(1.0_r8-Pbl_NRich(i,k))*Pbl_Shear(i,k))
          !
          ! Km = l*q*Sm
          !  
          ! Km     q 
          !---- = ---*Sm 
          ! l^2    l
          !
          ! Kh = l*q*Sh
          !
          !  
          ! Kh     q 
          !---- = ---*Sh 
          ! l^2    l
          !
          Pbl_KmMixl(i,k)=Pbl_Sqrtw(i,k)*Pbl_SmBar(i,k)
          Pbl_KhMixl(i,k)=Pbl_Sqrtw(i,k)*Pbl_ShBar(i,k)
       END DO
    END DO
    !     
    !     Pbl_HgtLyI(1)   height at the layer interface 
    !
    !        R*T
    ! DZ = ------ * DP
    !        g*P
    !
    !        R*T      DP
    ! DZ = ------ * -----
    !         g       P
    !
    !         R      DP
    ! DZ = ------ * ----- * T
    !         g       P
    !
    !      
    ! DZ = con0(k)  * T
    !
    !                                              kg*m*m    
    !                                             --------  
    !              R      DP       (j/kg/k)        kg*K*s*s  
    ! con0(k) = ------ * ----- =  --------   =   ------------    = m/K
    !              g       P        m/s*s             m       
    !                                               -----    
    !                                                s*s          
    !            gasr        si(k) - si(k+1)
    ! con0(k)= -------- * ------------------------
    !            grav             sl(k)
    !
    ! Pbl_HgtLyI = con0(k)  * T
    !
    k=1
    DO i = 1, ncols
       Pbl_HgtLyI(i,k)=con0(i,k)*gt(i,k)
    END DO

    DO k = 2, kmax
       DO i = 1, ncols
          Pbl_HgtLyI(i,k)=Pbl_HgtLyI(i,k-1)+con0(i,k)*gt(i,k)
       END DO
    END DO
    DO itr = 1, nitr
       !     
       !     Pbl_EddEner(2)   mixing length
       !     Pbl_EddEner(2)   b     :b**2 is eddy enegy
       !     
       !...gl0   maximum mixing length l0 in blackerdar's formula
       !         this is retained as a first guess for next time step
       !                                      k0*z
       !                                l = --------
       !                                    (1 + k0*z/l0)

       gln = 0.0_r8
       gld = 0.0_r8
       Pbl_EddEner(:,1) = 0.0_r8
       DO k = 1, kmax-1
          DO i = 1, ncols
             !
             !                 vk0*gl0(i)*Z(i,k) 
             ! Pbl_EddEner = ------------------------------
             !                (gl0(i) + vk0*Z(i,k))
             !
             Pbl_EddEner(i,k+1)=vk0*gl0(i)*Pbl_HgtLyI(i,k) / (gl0(i)+vk0*Pbl_HgtLyI(i,k))
             !
             !                 vk0*gl0(i)*Z(i,k)           q          q
             ! Pbl_EddEner = ------------------------ * --------- = ---------
             !                (gl0(i) + vk0*Z(i,k))        l^2        l
             !
             Pbl_EddEner(i,k+1)=Pbl_EddEner(i,k+1)*Pbl_Sqrtw(i,k)
          END DO
       END DO
       k=1
       DO i = 1, ncols
          Pbl_EddEner(i,k)= 1.0e-3_r8
       END DO
       k=1
       DO i = 1, ncols
          x=0.5_r8*((prsi(i,k)/prsi(i,1)) - (prsi(i,k+1)/prsi(i,1)))*(Pbl_EddEner(i,k)+Pbl_EddEner(i,k+1))
          y=x*0.5_r8*Pbl_HgtLyI(i,k)
          gld(i)=gld(i)+x
          gln(i)=gln(i)+y
       END DO
       k=kmax
       DO i = 1, ncols
          x=0.5_r8*((prsi(i,k)/prsi(i,1)) - (prsi(i,k+1)/prsi(i,1)))*Pbl_EddEner(i,k)
          y=x*0.5_r8*(Pbl_HgtLyI(i,k)+Pbl_HgtLyI(i,k-1))
          gln(i)=gln(i)+y
          gld(i)=gld(i)+x
       END DO
       IF (kmax > 2) THEN
          DO k = 2, kmax-1
             DO i = 1, ncols
                x=0.5_r8*((prsi(i,k)/prsi(i,1)) - (prsi(i,k+1)/prsi(i,1)))*(Pbl_EddEner(i,k)+Pbl_EddEner(i,k+1))
                y=x*0.5_r8*(Pbl_HgtLyI(i,k)+Pbl_HgtLyI(i,k-1))
                gln(i)=gln(i)+y
                gld(i)=gld(i)+x
             END DO
          END DO
       END IF
       DO i = 1, ncols
          !                0.5 * del(k) * (q(i,k)+q(i,k+1)) * 0.5*(Z(i,k)+Z(i,k-1))
          !gl0(i)=facl * -----------------------------------------------------------
          !                        0.5 * del(k) * (q(i,k)+q(i,k+1))
          gl0(i)=facl*gln(i)/gld(i)
       END DO
       !     
       !     iteration that determines mixing length
       !     
    END DO
    !
    !     Pbl_MixLgh(5)   mixing length
    !     
    DO k = 1, kmax-1
       DO i=1,ncols
          !        vk0*gl0(i)*Z(i,k)         m^2
          ! l  = -------------------------------- = ------
          !        gl0(i) + vk0*Z(i,k)       m
          Pbl_MixLgh(i,k)=vk0*gl0(i)*Pbl_HgtLyI(i,k)/(gl0(i)+vk0*Pbl_HgtLyI(i,k))
       END DO
    END DO
    !
    !
    !     
    !     Pbl_CoefKm(1)   km = l*q*Sm
    !     Pbl_CoefKh(2)   kh = l*q*Sh
    ! 
    !     where KM and KH are the diffusion coefficients for momentum 
    !     and heat, respectively, l is the master turbulence length scale, 
    !     q2 is the turbulent kinetic energy (so q is the magnitude of 
    !     the turbulent wind velocity), and SM and SH are momentum 
    !     flux and heat flux stability parameters, respectively    
    !    
    IF (kmean == 0) THEN
       DO k = 1, kmax-1
          DO i = 1, ncols
             !
             !             Km 
             !Km = l^2 * -----
             !            l^2
             !
             Pbl_CoefKm(i,k)=MIN(gkm1,MAX(gkm0,Pbl_MixLgh(i,k)**2*Pbl_KmMixl(i,k)))
             !
             !             Kh 
             !Kh = l^2 * -----
             !            l^2
             !
             Pbl_CoefKh(i,k)=MIN(gkh1,MAX(gkh0,Pbl_MixLgh(i,k)**2*Pbl_KhMixl(i,k)))
          END DO
       END DO
    ELSE
       DO k = 1, kmax-1
          DO i = 1, ncols
             !
             !             Km 
             !Km = l^2 * -----
             !            l^2
             !
             Pbl_KM(i,k)=Pbl_MixLgh(i,k)**2*Pbl_KmMixl(i,k)
             !
             !             Kh 
             !Kh = l^2 * -----
             !            l^2
             !
             Pbl_KH(i,k)=Pbl_MixLgh(i,k)**2*Pbl_KhMixl(i,k)
          END DO
       END DO
       fac=0.25_r8
       IF (kmax >= 4) THEN
          DO k = 2, kmax-2
             DO i = 1, ncols
                !
                !             k=2  ****Km(k),sl*** } -----------
                !             k=3/2----si,ric,rf,km,kh,b,l -----------
                !             k=1  ****Km(k),sl*** } -----------
                !             k=1/2----si ----------------------------
                !
                !       Km(k-1) + 2*Km(k) + Km(k+1)
                ! Km = -------------------------------
                !                  4
                !
                Pbl_CoefKm(i,k)=fac*(Pbl_KM(i,k-1)+2.0_r8*Pbl_KM(i,k)+Pbl_KM(i,k+1))
                !
                !       Kh(k-1) + 2*Kh(k) + Kh(k+1)
                ! Kh = -------------------------------
                !                  4
                !
                Pbl_CoefKh(i,k)=fac*(Pbl_KH(i,k-1)+2.0_r8*Pbl_KH(i,k)+Pbl_KH(i,k+1))
             END DO
          END DO
       END IF
       DO i = 1, ncols
          !
          !       Km(1) + Km(2)
          ! Km = ---------------
          !            2
          !
          !
          !       Kh(1) + Kh(2)
          ! Kh = ---------------
          !            2
          !
          Pbl_CoefKm(i,     1)=0.5_r8*(Pbl_KM(i,     1)+Pbl_KM(i,     2))
          Pbl_CoefKh(i,     1)=0.5_r8*(Pbl_KH(i,     1)+Pbl_KH(i,     2))
          Pbl_CoefKm(i,kmax-1)=0.5_r8*(Pbl_KM(i,kmax-1)+Pbl_KM(i,kmax-2))
          Pbl_CoefKh(i,kmax-1)=0.5_r8*(Pbl_KH(i,kmax-1)+Pbl_KH(i,kmax-2))
       END DO
       DO k = 1, kmax-1
          DO i = 1, ncols
             Pbl_CoefKm(i,k)=MIN(gkm1,MAX(gkm0,Pbl_CoefKm(i,k))) !(m/sec**2)
             Pbl_CoefKh(i,k)=MIN(gkh1,MAX(gkh0,Pbl_CoefKh(i,k))) !(m/sec**2)
          END DO
       END DO
       rfx=rfc-0.001_r8
       IF (kmax >= 3) THEN

          icnt=0
          DO k = kmax-1, 2, -1
             DO i = 1, ncols
                IF (Pbl_NRich(i,k) > rfx .AND. Pbl_NRich(i,k-1) <= rfx) THEN
                   icnt(i)=icnt(i)+1
                   ichk(i,icnt(i))=k
                END IF
             END DO
          END DO

          DO i = 1, ncols
             IF (icnt(i) /= 0) THEN
                Pbl_CoefKm(i,ichk(i,1))=gkm0  !(m/sec**2)
                Pbl_CoefKh(i,ichk(i,1))=gkh0  !(m/sec**2)
             END IF
          END DO
       END IF
    END IF
    !                            --           --
    !      d -(U'W')       d    |        d U    | 
    !   ------------- = --------|KM *  -------- |
    !      dt              dt   |        d z    |
    !                            --           --
    !     momentum diffusion
    ! 
    CALL VDIFV(kMax,nCols,twodti,Pbl_ITemp,Pbl_CoefKm,a,b,gu,gv,gmu)
    !
    !                            --           --
    !      d -(Q'W')       d    |        d Q    | 
    !   ------------- = --------|KH *  -------- |
    !      dt              dt   |        d z    |
    !                            --           --
    !     
    !     water vapour diffusion
    !
    CALL VDIFH(kMax,nCols,twodti,Pbl_ITemp,Pbl_CoefKh,a,b,gq,gmq)       
    !     
    !                            --           --
    !      d -(T'W')       d    |        d T    | 
    !   ------------- = --------|KH *  -------- |
    !      dt              dt   |        d z    |
    !                            --           --
    !    
    !     sensible heat diffusion
    ! 
    CALL VDIFT(kMax,nCols,twodti,sigriv,sigr,Pbl_ITemp,Pbl_CoefKh,a,b,gt,gmt)

  END SUBROUTINE TubulenceMY2
  !---------------------------------------------------------------------
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !---------------------------------------------------------------------
  SUBROUTINE VDIFV(kMax,nCols,twodti,Pbl_ITemp,Pbl_CoefKm,a,b,gu,gv,gmu)
    !     ***************************************************************
    !     *                                                             *
    !     *        VERTICAL DIFFUSION OF VELOCITY COMPONENTS            *
    !     *                                                             *
    !     ***************************************************************
    !---------------------------------------------------------------------
    INTEGER      ,    INTENT(in   ) :: nCols
    INTEGER      ,    INTENT(in   ) :: kMax
    REAL(KIND=r8),    INTENT(in   ) :: twodti
    REAL(KIND=r8),    INTENT(in   ) :: a         (nCols,kmax)      !s * K**2/m**2 
    REAL(KIND=r8),    INTENT(in   ) :: b         (nCols,kmax)      !s * K**2/m**2 
    REAL(KIND=r8),    INTENT(in   ) :: Pbl_ITemp (nCols,kMax)!1/K**2
    REAL(KIND=r8),    INTENT(inout) :: Pbl_CoefKm(nCols,kMax)!(m/sec**2)
    REAL(KIND=r8),    INTENT(in   ) :: gu        (nCols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: gv        (nCols,kMax)
    REAL(KIND=r8),    INTENT(inout) :: gmu       (nCols,kMax,4)

    REAL(KIND=r8) :: Pbl_DifVzn(nCols,kMax) 
    REAL(KIND=r8) :: Pbl_DifVmd(nCols,kMax)
    REAL(KIND=r8) :: Pbl_KHbyDZ(nCols,kMax) !1/(sec * K) 
    REAL(KIND=r8) :: Pbl_KMbyDZ2(nCols,kMax) !1/(sec * K) 
    REAL(KIND=r8) :: Pbl_KMbyDZ1(nCols,kMax) !1/(sec * K) 
    REAL(KIND=r8) :: Pbl_KHbyDZ2(nCols,kMax) !1/(sec * K)
    REAL(KIND=r8) :: Pbl_KMbyDZ(nCols,kMax) !1/(sec * K)
    REAL(KIND=r8) :: Pbl_TendU(nCols,kMax) !1/(sec * K)
    REAL(KIND=r8) :: Pbl_TendV(nCols,kMax) !1/(sec * K)

    INTEGER       :: i
    INTEGER       :: k

    !     
    !     momentum diffusion
    !     
    DO k = 1, kMax-1
       DO i = 1, nCols
          Pbl_CoefKm(i,k  )=Pbl_CoefKm(i,k)*Pbl_ITemp(i,k)!(m/sec**2) * (1/K**2) = m/(sec**2 * K**2)
          !
          !                   --   --      -- -- 2   
          !                  |       |    |  1  |         K**2 * s
          !a(k)=             |  2*Dt |  * | --- |     ==>-----------
          !                  |       |    |  dZ |           m**2  
          !                   --   --      -- --     
          !                                         K**2 * s         m              1
          Pbl_KHbyDZ(i,k  )=a(i,k  )*Pbl_CoefKm(i,k)!----------- * -------------  = --------
          !                                           m**2          s**2*K**2     m * s
          !
          !                                         K**2 * s        m**2             1
          Pbl_KMbyDZ(i,k+1)=b(i,k+1)*Pbl_CoefKm(i,k)!----------- * -------------  = --------
          !                                           m**2          s**2*K**2         s
          !
          !     gwrk(1)   difference of pseudo v wind ( km is destroyed )
          !     gwrk(5)   difference of pseudo u wind ( b  is destroyed )
          !     
          Pbl_DifVzn(i,k)=gu(i,k)-gu(i,k+1)
          Pbl_DifVmd(i,k)=gv(i,k)-gv(i,k+1)
       END DO
    END DO
    DO i = 1, nCols
       Pbl_KMbyDZ2   (i,1)=0.0_r8
       Pbl_KMbyDZ1   (i,1)=1.0_r8 + Pbl_KHbyDZ(i,1)
       Pbl_KHbyDZ2   (i,1)=        -Pbl_KHbyDZ(i,1)
       !                          --             --
       ! DU       d(w'u')      d |          d U    |       m
       !------ = ------- =    ---| - Km *  ------  |  = --------
       ! Dt       dz           dz|          dz     |     s * s
       !                          --             --

       Pbl_TendU(i,1)=-twodti*Pbl_KHbyDZ(i,1)*Pbl_DifVzn(i,1)!(1/m)*(m/s)
       !
       !                     1          1           m           m
       !Pbl_TendU(i,1) = - ------ * -------- *  -------- =   --------
       !                     s          s           s         s * s
       !
       Pbl_TendV(i,1)=-twodti*Pbl_KHbyDZ(i,1)*Pbl_DifVmd(i,1)!m/s**2

       Pbl_KMbyDZ2   (i,kMax)=       - Pbl_KMbyDZ(i,kMax)
       Pbl_KMbyDZ1   (i,kMax)=1.0_r8 + Pbl_KMbyDZ(i,kMax)
       Pbl_KHbyDZ2   (i,kMax)=0.0_r8
       Pbl_TendU(i,kMax)=twodti*Pbl_KMbyDZ(i,kMax)*Pbl_DifVzn(i,kMax-1)
       Pbl_TendV(i,kMax)=twodti*Pbl_KMbyDZ(i,kMax)*Pbl_DifVmd(i,kMax-1)
    END DO
    DO k = 2, kMax-1
       DO i = 1, nCols
          Pbl_KMbyDZ2   (i,k)=      -Pbl_KMbyDZ(i,k)                   !1/(sec * K)
          Pbl_KMbyDZ1   (i,k)=1.0_r8+Pbl_KHbyDZ(i,k)+Pbl_KMbyDZ(i,k)   !1/(sec * K)
          Pbl_KHbyDZ2   (i,k)=      -Pbl_KHbyDZ(i,k)                   !1/(sec * K)
          Pbl_TendU(i,k)=(Pbl_KMbyDZ(i,k)*Pbl_DifVzn(i,k-1)  - &
               Pbl_KHbyDZ(i,k)*Pbl_DifVzn(i,k  )) * twodti
          Pbl_TendV(i,k)=(Pbl_KMbyDZ(i,k)*Pbl_DifVmd(i,k-1)  - &
               Pbl_KHbyDZ(i,k)*Pbl_DifVmd(i,k  )) * twodti
       END DO
    END DO
    DO k = kmax-1, 1, -1
       DO i = 1, ncols
          Pbl_KHbyDZ2   (i,k)=Pbl_KHbyDZ2 (i,k)/Pbl_KMbyDZ1(i,k+1)
          !
          !                                - Pbl_KMbyDZ_1(i,k) 
          !Pbl_KHbyDZ2   (i,k)=-------------------------------------------------
          !                     1.0 + Pbl_KMbyDZ_1(i,k+1) + Pbl_KMbyDZ_2(i,k+1)
          !
          Pbl_KMbyDZ1   (i,k)=Pbl_KMbyDZ1 (i,k)-Pbl_KHbyDZ2(i,k)*Pbl_KMbyDZ2(i,k+1)
          !                                                                         Pbl_KMbyDZ_1(i,k)*Pbl_KMbyDZ_2(i,k+1)   
          !Pbl_KMbyDZ1   (i,k)=1.0_r8 + Pbl_KMbyDZ_1(i,k) + Pbl_KMbyDZ_2(i,k) - -------------------------------------------------
          !                                                                      1.0 + Pbl_KMbyDZ_1(i,k+1) + Pbl_KMbyDZ_2(i,k+1)
          Pbl_TendU(i,k)=Pbl_TendU(i,k)-Pbl_KHbyDZ2(i,k)*Pbl_TendU(i,k+1)
          Pbl_TendV(i,k)=Pbl_TendV(i,k)-Pbl_KHbyDZ2(i,k)*Pbl_TendV(i,k+1)
       END DO
    END DO
    DO k = 1, kmax
       DO i = 1, ncols
          gmu(i,k,1)=Pbl_KMbyDZ2 (i,k)
          gmu(i,k,2)=Pbl_KMbyDZ1 (i,k)
          gmu(i,k,3)=Pbl_TendU   (i,k)
          gmu(i,k,4)=Pbl_TendV   (i,k)
       END DO
    END DO
  END SUBROUTINE VDIFV  
  !---------------------------------------------------------------------
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !---------------------------------------------------------------------
  SUBROUTINE VDIFH(kMax,nCols,twodti,Pbl_ITemp,Pbl_CoefKh,a,b,gq,gmq)
    !     ***************************************************************
    !     *                                                             *
    !     *         VERTICAL DIFFUSION OF MASS VARIABLES                *
    !     *                                                             *
    !     ***************************************************************
    INTEGER      ,    INTENT(in   ) :: nCols
    INTEGER      ,    INTENT(in   ) :: kMax
    REAL(KIND=r8),    INTENT(in   ) :: twodti
    REAL(KIND=r8),    INTENT(in   ) :: a         (nCols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: b         (nCols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: Pbl_ITemp (nCols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: Pbl_CoefKh(nCols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: gq        (nCols,kMax)
    REAL(KIND=r8),    INTENT(inout) :: gmq       (nCols,kMax,3)

    REAL(KIND=r8) :: Pbl_DifQms(nCols,kMax) 
    REAL(KIND=r8) :: Pbl_CoefKh2(nCols,kMax)

    REAL(KIND=r8) :: Pbl_KHbyDZ(nCols,kMax) !1/(sec * K) 
    REAL(KIND=r8) :: Pbl_KMbyDZ2(nCols,kMax) !1/(sec * K) 
    REAL(KIND=r8) :: Pbl_KMbyDZ1(nCols,kMax) !1/(sec * K) 
    REAL(KIND=r8) :: Pbl_KHbyDZ2(nCols,kMax) !1/(sec * K)
    REAL(KIND=r8) :: Pbl_KMbyDZ(nCols,kMax) !1/(sec * K)
    REAL(KIND=r8) :: Pbl_TendQ(nCols,kMax) !1/(sec * K)

    INTEGER       :: i
    INTEGER       :: k

    DO k = 1, kMax-1
       DO i = 1, nCols
          Pbl_CoefKh2(i,k  )=Pbl_CoefKh(i,k)*Pbl_ITemp(i,k)
          Pbl_KHbyDZ(i,k  )=a(i,k  )*Pbl_CoefKh2(i,k)
          Pbl_KMbyDZ(i,k+1)=b(i,k+1)*Pbl_CoefKh2(i,k)
          !     
          !     Pbl_DifQms(1)   difference of specific humidity
          !     
          Pbl_DifQms(i,k)=gq(i,k)-gq(i,k+1)  
       END DO
    END DO
    DO i = 1, ncols
       Pbl_KMbyDZ2(i,1)=0.0_r8
       Pbl_KMbyDZ1(i,1)=1.0_r8  + Pbl_KHbyDZ(i,1)
       Pbl_KHbyDZ2(i,1)=        - Pbl_KHbyDZ(i,1)
       Pbl_TendQ  (i,1)=-twodti * Pbl_KHbyDZ(i,1) * Pbl_DifQms(i,1)

       Pbl_KMbyDZ2(i,kMax)=        - Pbl_KMbyDZ(i,kMax)
       Pbl_KMbyDZ1(i,kMax)=1.0_r8  + Pbl_KMbyDZ(i,kMax)
       Pbl_KHbyDZ2(i,kMax)=0.0_r8
       Pbl_TendQ  (i,kMax)= twodti * Pbl_KMbyDZ(i,kMax) * Pbl_DifQms(i,kMax-1)
    END DO
    DO k = 2, kmax-1
       DO i = 1, ncols
          Pbl_KMbyDZ2(i,k)=       - Pbl_KMbyDZ(i,k)
          Pbl_KMbyDZ1(i,k)=1.0_r8 + Pbl_KHbyDZ(i,k)+Pbl_KMbyDZ(i,k)
          Pbl_KHbyDZ2(i,k)=       - Pbl_KHbyDZ(i,k)
          Pbl_TendQ(i,k)=(Pbl_KMbyDZ(i,k) * Pbl_DifQms(i,k-1)- &
               Pbl_KHbyDZ(i,k) * Pbl_DifQms(i,k  ))*twodti
       END DO
    END DO
    DO k = kmax-1, 1, -1
       DO i = 1, ncols
          Pbl_KHbyDZ2(i,k)=Pbl_KHbyDZ2(i,k) / Pbl_KMbyDZ1(i,k+1)
          Pbl_KMbyDZ1(i,k)=Pbl_KMbyDZ1(i,k) - Pbl_KHbyDZ2(i,k)*Pbl_KMbyDZ2(i,k+1)
          Pbl_TendQ  (i,k)=Pbl_TendQ  (i,k) - Pbl_KHbyDZ2(i,k)*Pbl_TendQ  (i,k+1)
       END DO
    END DO
    DO k = 1, kmax
       DO i = 1, ncols
          gmq(i,k,1)=Pbl_KMbyDZ2(i,k)
          gmq(i,k,2)=Pbl_KMbyDZ1(i,k)
          gmq(i,k,3)=Pbl_TendQ  (i,k)
       END DO
    END DO

  END SUBROUTINE VDIFH
  !----------------------------------------------------------------------
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !----------------------------------------------------------------------
 SUBROUTINE VDIFT(kMax,nCols,twodti,sigriv,sigr,Pbl_ITemp,Pbl_CoefKh,a,b,gt,gmt)
    !     ***************************************************************
    !     *                                                             *
    !     *         VERTICAL DIFFUSION OF MASS VARIABLES                *
    !     *                                                             *
    !     ***************************************************************
    INTEGER      ,    INTENT(in   ) :: nCols
    INTEGER      ,    INTENT(in   ) :: kMax
    REAL(KIND=r8),    INTENT(in   ) :: twodti
    REAL(KIND=r8),    INTENT(in   ) :: sigriv    (nCols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: sigr      (nCols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: a         (nCols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: b         (nCols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: Pbl_ITemp (nCols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: Pbl_CoefKh(nCols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: gt        (nCols,kMax)
    REAL(KIND=r8),    INTENT(inout) :: gmt       (nCols,kMax,3)

    REAL(KIND=r8) :: Pbl_CoefKh2(nCols,kMax)
    REAL(KIND=r8) :: Pbl_KHbyDZ(nCols,kMax) !1/(sec * K) 
    REAL(KIND=r8) :: Pbl_KMbyDZ(nCols,kMax) !1/(sec * K) 
    REAL(KIND=r8) :: Pbl_KMbyDZ2(nCols,kMax)
    REAL(KIND=r8) :: Pbl_KMbyDZ1(nCols,kMax)
    REAL(KIND=r8) :: Pbl_KHbyDZ2(nCols,kMax)
    REAL(KIND=r8) :: Pbl_TendT(nCols,kMax)
    INTEGER       :: i
    INTEGER       :: k
    !     
    !     sensible heat diffusion
    !     
    !   sigk  (   k)=       sl(k)**akappa
    !   sigkiv(   k)=1.0_r8/sl(k)**akappa
    !   sigr  (   k)=sigk(k  )*sigkiv(k+1) =(sl(k  )**akappa)*(1.0_r8/sl(k+1)**akappa)
    !   sigriv( k+1)=sigk(k+1)*sigkiv(k  ) =(sl(k+1)**akappa)*(1.0_r8/sl(k  )**akappa)
    !
    DO k = 1, kMax-1
       DO i = 1, nCols
          Pbl_CoefKh2(i,k  )=Pbl_CoefKh(i,k)*Pbl_ITemp(i,k)
          Pbl_KHbyDZ(i,k  )=a(i,k  )*Pbl_CoefKh2(i,k)
          Pbl_KMbyDZ(i,k+1)=b(i,k+1)*Pbl_CoefKh2(i,k)
       END DO
    END DO
    DO i = 1, nCols
       Pbl_KMbyDZ2(i,1)=  0.0_r8
       Pbl_KMbyDZ1(i,1)=  1.0_r8+Pbl_KHbyDZ(i,1)
       Pbl_KHbyDZ2(i,1)=-sigr(i,1)*Pbl_KHbyDZ(i,1)
       Pbl_TendT  (i,1)=-Pbl_KHbyDZ(i,1)*(gt(i,1)-sigr(i,1)*gt(i,1+1))*twodti

       Pbl_KMbyDZ2(i,kMax)=-sigriv(i,kMax)*Pbl_KMbyDZ(i,kMax)
       Pbl_KMbyDZ1(i,kMax)=    1.0_r8+Pbl_KMbyDZ(i,kMax)
       Pbl_KHbyDZ2(i,kMax)=0.0_r8
       Pbl_TendT  (i,kMax)=twodti*Pbl_KMbyDZ(i,kMax)*(sigriv(i,kMax)*gt(i,kMax-1)-gt(i,kMax))
    END DO
    DO k = 2, kMax-1
       DO i = 1, nCols
          Pbl_KMbyDZ2(i,k)=-sigriv(i,k)*Pbl_KMbyDZ(i,k)
          Pbl_KMbyDZ1(i,k)=1.0_r8+Pbl_KHbyDZ(i,k)+Pbl_KMbyDZ(i,k)
          Pbl_KHbyDZ2(i,k)=-sigr  (i,k)*Pbl_KHbyDZ(i,k)
          Pbl_TendT  (i,k)=( Pbl_KMbyDZ(i,k)*(sigriv(i,k)*gt(i,k-1) - gt(i,k))&
               -Pbl_KHbyDZ(i,k)*(gt(i,k)- sigr(i,k)*gt(i,k+1))  )*twodti
       END DO
    END DO

    DO k = kmax-1, 1, -1
       DO i = 1, ncols
          Pbl_KHbyDZ2(i,k)=Pbl_KHbyDZ2(i,k)/Pbl_KMbyDZ1(i,k+1)
          Pbl_KMbyDZ1(i,k)=Pbl_KMbyDZ1(i,k)-Pbl_KHbyDZ2(i,k  )*Pbl_KMbyDZ2(i,k+1)
          Pbl_TendT  (i,k)=Pbl_TendT  (i,k)-Pbl_KHbyDZ2(i,k  )*Pbl_TendT  (i,k+1)
       END DO
    END DO

    DO k = 1, kmax
       DO i = 1, ncols
          gmt(i,k,1)=Pbl_KMbyDZ2(i,k)
          gmt(i,k,2)=Pbl_KMbyDZ1(i,k)
          gmt(i,k,3)=Pbl_TendT  (i,k)
       END DO
    END DO
  END SUBROUTINE VDIFT

END MODULE Sfc_SeaFlux_COLA_Model
