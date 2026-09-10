!
!  $Author: pkubota $
!  $Date: 2011/04/20 14:34:46 $
!  $Revision: 1.10 $
!
MODULE FieldsDynamics

  USE Constants, Only: r8
  USE Options,   Only: &
      slagr,           &
      slhum,           &
      grid_difus,      &
      microphys,       &
      nClass,          &
      nAeros,          &
      SL_twotime_scheme

   IMPLICIT NONE
  SAVE       


  ! Spectral fields: 13 2D and 10 1D


  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qtmpp
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qtmpt
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qtmpt_si
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qrotp
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qrott
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qrott_si
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qdivp
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qdivt
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qdivt_si
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qqp
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qice
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: qvar
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qliq
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qup
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qvp
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qtphi
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qqphi
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qdiaten
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ) :: qlnpp
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ) :: qlnp_nabla4
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ) :: qlnpl
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ) :: qlnpt
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ) :: qlnpt_si
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ) :: qgzslap
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ) :: qgzs
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ) :: qplam
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ) :: qpphi
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ) :: qgzsphi
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: qozon  ! add solange 27-01-2012

  ! Gaussian fields: 28 + 2 * nscalars 3D and 12 2D

  INTEGER       :: adr_scalars
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:,:,:) :: fgpass_scalars
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgpass_fluxscalars

  ! (du/dt) tendency of zonal wind * cos(latitude)
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgyu !T  
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgyum!T-1
  ! -(dv/dt) negative of tendency of v*cos(latitude)
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgyv !T  
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgyvm!T-1
  ! specific humidity tendency (kg/kg/s)
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgqd !T  
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgqdm!T-1 
  ! specific humidity (kg/kg)
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgq  !T
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgqmm!T-1 (depois da convec)
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgqm!T-1 (antes da convec)
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgqp!T+1
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgqTot!T+1

  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgicem  
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgice  
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgicep  
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgicet  

  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:,:) :: fgvarm  
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:,:) :: fgvar  
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:,:) :: fgvarp  
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:,:) :: fgvart  

  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgliqm  
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgliq  
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgliqp  
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgliqt  
  ! virtual temperature tendency (K/s)
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgtd !T  
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgtdm!T-1
  ! virtual temperature (K)
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgtmp !T
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgtmpm!T-1
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgtmpp!T+1
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgtmpC!T+1 nilo1
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgtmpp2!T+1
  ! zonal wind (?)
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgu !T  
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgum!T-1
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fguC!T+1 nilo2

  ! meridional wind (?)
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgv !T  
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgvm!T-1
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgvC!T+1 nilo3

  ! vertical wind (?)
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgw !T  
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgwm!T-1 
  ! 
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgdiv !T  
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgdivm!T-1
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgdivC!!T+1 nilo4

  !
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgtphi !T  
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgtphim!T-1
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgtphiC!T+1 nilo5

  !
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fgpphi
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fgpphim
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fgpphiC !T+1 nilo7
  !
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgtlam
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgtlamm
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgtlamC !T+1 nilo6
  !
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fgplam
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fgplamm
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fgplamC !T+1 nilo8
  !
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgqlam
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgulam
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgvlam

  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgrot
  ! vertical omega(?) velocity (pascal/s)
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: omg!T

  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgqphi
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fguphi
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgvphi
  ! pressure at layers and interfaces (pascal) 
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgprsl
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgprsi
  ! height at layers and interfaces (meters) 
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgphil
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: fgphii

  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fgumean
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fgvmean
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fgvdlnp
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fglnps!T  
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fglnps_nabla4
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fglnpm!T-1
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fglnpsC!T+1 !nilo9
  ! surface pressure (pascal) 
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fgps !T
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fgpsp!T+1
  ! geopotential at the surface, fgzs/grav == topography
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fgzs
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fgzslam
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fgzsphi
  REAL(KIND=r8), ALLOCATABLE, TARGET, DIMENSION(:  ,:) :: fgvdlnpm








CONTAINS
  SUBROUTINE InitFields(ibMax, kMax, jbMax, mnMax, mnExtMax, &
                        kmaxloc, MNMax_si, jbMax_ext, nscalars)
    INTEGER, INTENT(IN) :: ibMax
    INTEGER, INTENT(IN) :: kMax
    INTEGER, INTENT(IN) :: kMaxloc
    INTEGER, INTENT(IN) :: jbMax
    INTEGER, INTENT(IN) :: jbMax_ext
    INTEGER, INTENT(IN) :: mnMax
    INTEGER, INTENT(IN) :: MNMax_si
    INTEGER, INTENT(IN) :: mnExtMax
    INTEGER, INTENT(IN) :: nscalars
    ALLOCATE (qtmpp(2*mnMax,kMaxloc))
    qtmpp=0.0_r8
    ALLOCATE (qtmpt(2*mnMax,kMaxloc))
    qtmpt=0.0_r8
    ALLOCATE (qtmpt_si(2*mnMax_si,kMax))
    qtmpt_si=0.0_r8
    ALLOCATE (qrotp(2*mnMax,kMaxloc))
    qrotp=0.0_r8
    ALLOCATE (qrott(2*mnMax,kMaxloc))
    qrott=0.0_r8
    ALLOCATE (qrott_si(2*mnMax_si,kMax))
    qrott_si=0.0_r8
    ALLOCATE (qdivp(2*mnMax,kMaxloc))
    qdivp=0.0_r8
    ALLOCATE (qdivt(2*mnMax,kMaxloc))
    qdivt=0.0_r8
    ALLOCATE (qdivt_si(2*mnMax_si,kMax))
    qdivt_si=0.0_r8
    ALLOCATE (qqp(2*mnMax,kMaxloc))
    qqp=0.0_r8
    IF (microphys) THEN
       IF((nClass+nAeros)>0)THEN
          ALLOCATE (qvar(2*mnMax,kMaxloc,nClass+nAeros))
          qvar=0.0_r8
       END IF
       ALLOCATE (qice(2*mnMax,kMaxloc))
       qice=0.0_r8
       ALLOCATE (qliq(2*mnMax,kMaxloc))
       qliq=0.0_r8
    ENDIF
    ALLOCATE (qdiaten(2*mnMax,kMaxloc))
    qdiaten=0.0_r8
    ALLOCATE (qup(2*mnExtMax,kMaxloc))
    qup=0.0_r8
    ALLOCATE (qvp(2*mnExtMax,kMaxloc))
    qvp=0.0_r8
    ALLOCATE (qtphi(2*mnExtMax,kMaxloc))
    qtphi=0.0_r8
    ALLOCATE (qqphi(2*mnExtMax,kMaxloc))
    qqphi=0.0_r8
    ALLOCATE (qlnpp(2*mnMax))
    qlnpp=0.0_r8
    ALLOCATE (qlnpl(2*mnMax))
    if (grid_difus) then 
       ALLOCATE (qlnp_nabla4(2*mnMax))
       qlnp_nabla4=0.0_r8
    endif
    qlnpl=0.0_r8
    ALLOCATE (qlnpt(2*mnMax))
    qlnpt=0.0_r8
    ALLOCATE (qlnpt_si(2*mnMax_si))
    qlnpt_si=0.0_r8
    ALLOCATE (qgzslap(2*mnMax))
    qgzslap=0.0_r8
    ALLOCATE (qgzs(2*mnMax))
    qgzs=0.0_r8
    ALLOCATE (qplam(2*mnMax))
    qplam=0.0_r8
    ALLOCATE (qpphi(2*mnExtMax))
    qpphi=0.0_r8
    ALLOCATE (qgzsphi(2*mnExtMax))
    qgzsphi=0.0_r8
    ALLOCATE (qozon(2*mnMax,kMaxloc))  ! solange add
    qozon=0.0_r8

    IF (nscalars .gt. 0) THEN
       ALLOCATE (fgpass_scalars(ibMax,kMax,jbmax_ext,nscalars,2))
       fgpass_scalars=0.0_r8
      ELSE
       ALLOCATE (fgpass_scalars(1,1,1,1,1))
       fgpass_scalars=0.0_r8
    ENDIF
    adr_scalars = 1
    
    ALLOCATE (fgyu(ibMax,kMax,jbMax))
    fgyu=0.0_r8
    ALLOCATE (fgyv(ibMax,kMax,jbMax))
    fgyv=0.0_r8
    ALLOCATE (fgtd(ibMax,kMax,jbMax))
    fgtd=0.0_r8
    ALLOCATE (fgqd(ibMax,kMax,jbMax_ext))
    fgqd=0.0_r8
    ALLOCATE (fgdiv(ibMax,kMax,jbMax))
    fgdiv=0.0_r8
    ALLOCATE (fgrot(ibMax,kMax,jbMax))
    fgrot=0.0_r8
    ALLOCATE (fgq(ibMax,kMax,jbMax))
    fgq=0.0_r8
    ALLOCATE (fgtmp(ibMax,kMax,jbMax))
    fgtmp=0.0_r8
    IF (SL_twotime_scheme) THEN
       ALLOCATE (fgu(ibMax,kMax,jbMax))
       fgu=0.0_r8
       ALLOCATE (fgv(ibMax,kMax,jbMax))
       fgv=0.0_r8    
       ALLOCATE (fgw(ibMax,kMax,jbMax))
       fgw=0.0_r8
       ALLOCATE (fgum(ibMax,kMax,jbMax_ext))
       fgum=0.0_r8
       ALLOCATE (fguC(ibMax,kMax,jbMax_ext))
       fguC=0.0_r8
       ALLOCATE (fgvm(ibMax,kMax,jbMax_ext))
       fgvm=0.0_r8
       ALLOCATE (fgvC(ibMax,kMax,jbMax_ext))
       fgvC=0.0_r8
       ALLOCATE (fgwm(ibMax,kMax,jbMax_ext))
       fgwm=0.0_r8
       ALLOCATE (fgtmpp(ibMax,kMax,jbMax_ext))
       fgtmpp=0.0_r8
       ALLOCATE(fgtmpC(ibMax,kMax,jbMax_ext))
       fgtmpC=0.0
       ALLOCATE(fgtmpp2(ibMax,kMax,jbMax_ext))
       fgtmpp2=0.0
       ALLOCATE (fgqp(ibMax,kMax,jbMax_ext))
       fgqp=0.0_r8
       ALLOCATE (fgqTot(ibMax,kMax,jbMax_ext))
       fgqTot=0.0_r8
       ALLOCATE (fgyum(ibMax,kMax,jbMax))
       fgyum=0.0_r8
       ALLOCATE (fgyvm(ibMax,kMax,jbMax))
       fgyvm=0.0_r8
       ALLOCATE (fgtdm(ibMax,kMax,jbMax))
       fgtdm=0.0_r8
       ALLOCATE (fgqdm(ibMax,kMax,jbMax))
       fgqdm=0.0_r8
      ELSE
       ALLOCATE (fgu(ibMax,kMax,jbMax_ext))
       fgu=0.0_r8
       ALLOCATE (fgv(ibMax,kMax,jbMax_ext))
       fgv=0.0_r8    
       ALLOCATE (fgw(ibMax,kMax,jbMax_ext))
       fgw=0.0_r8
       ALLOCATE (fgum(ibMax,kMax,jbMax))
       fgum=0.0_r8
       ALLOCATE (fguC(ibMax,kMax,jbMax))
       fguC=0.0_r8
       ALLOCATE (fgvm(ibMax,kMax,jbMax))
       fgvm=0.0_r8    
       ALLOCATE (fgvC(ibMax,kMax,jbMax))
       fgvC=0.0_r8    
       ALLOCATE (fgtmpp(ibMax,kMax,jbMax))
       fgtmpp=0.0_r8
       ALLOCATE (fgtmpC(ibMax,kMax,jbMax))
       fgtmpC=0.0_r8
       ALLOCATE (fgtmpp2(ibMax,kMax,jbMax))
       fgtmpp2=0.0_r8
       ALLOCATE (fgqp(ibMax,kMax,jbMax))
       fgqp=0.0_r8
       ALLOCATE (fgqTot(ibMax,kMax,jbMax))
       fgqTot=0.0_r8
       ALLOCATE (fgyum(ibMax,kMax,jbMax_ext))
       fgyum=0.0_r8
       ALLOCATE (fgyvm(ibMax,kMax,jbMax_ext))
       fgyvm=0.0_r8
       ALLOCATE (fgtdm(ibMax,kMax,jbMax_ext))
       fgtdm=0.0_r8
       ALLOCATE (fgqdm(ibMax,kMax,jbMax))
       fgqdm=0.0_r8
    ENDIF 
    IF((nAeros)>0)THEN
       ALLOCATE (fgpass_fluxscalars(ibMax,jbMax,nAeros))
    END IF
    IF (microphys) THEN 
       ALLOCATE (fgice(ibMax,kMax,jbMax))
       fgice=0.0_r8
       ALLOCATE (fgicep(ibMax,kMax,jbMax))
       fgicep=0.0_r8
       ALLOCATE (fgicem(ibMax,kMax,jbMax))
       fgicem=0.0_r8
       ALLOCATE (fgicet(ibMax,kMax,jbMax_ext))
       fgicet=0.0_r8
       IF((nClass+nAeros)>0)THEN
          ALLOCATE (fgvar(ibMax,kMax,jbMax,nClass+nAeros))
          fgvar=0.0_r8
          ALLOCATE (fgvarp(ibMax,kMax,jbMax,nClass+nAeros))
          fgvarp=0.0_r8
          ALLOCATE (fgvarm(ibMax,kMax,jbMax,nClass+nAeros))
          fgvarm=0.0_r8
          ALLOCATE (fgvart(ibMax,kMax,jbMax_ext,nClass+nAeros))
          fgvart=0.0_r8
       END IF 
       ALLOCATE (fgliq(ibMax,kMax,jbMax))
       fgliq=0.0_r8
       ALLOCATE (fgliqp(ibMax,kMax,jbMax))
       fgliqp=0.0_r8
       ALLOCATE (fgliqm(ibMax,kMax,jbMax))
       fgliqm=0.0_r8
       ALLOCATE (fgliqt(ibMax,kMax,jbMax_ext))
       fgliqt=0.0_r8
    ENDIF
    ALLOCATE (fgqmm(ibMax,kMax,jbMax))
    fgqmm=0.0_r8
    ALLOCATE (fgqm(ibMax,kMax,jbMax))
    fgqm=0.0_r8 
    ALLOCATE (omg(ibMax,kMax,jbMax))
    omg=0.0_r8
    ALLOCATE (fgprsl(ibMax,kMax,jbMax))
    fgprsl=0.0_r8
    ALLOCATE (fgprsi(ibMax,kMax+1,jbMax))
    fgprsi=0.0_r8
    ALLOCATE (fgphil(ibMax,kMax,jbMax))
    fgphil=0.0_r8
    ALLOCATE (fgphii(ibMax,kMax+1,jbMax))
    fgphii=0.0_r8

    ALLOCATE (fgtmpm(ibMax,kMax,jbMax))
    fgtmpm=0.0_r8
    ALLOCATE (fgdivm(ibMax,kMax,jbMax))
    fgdivm=0.0_r8
    ALLOCATE (fgdivC(ibMax,kMax,jbMax))
    fgdivC=0.0_r8
    ALLOCATE (fgtphi(ibMax,kMax,jbMax))
    fgtphi=0.0_r8
    ALLOCATE (fgqphi(ibMax,kMax,jbMax))
    fgqphi=0.0_r8
    ALLOCATE (fguphi(ibMax,kMax,jbMax))
    fguphi=0.0_r8
    ALLOCATE (fgvphi(ibMax,kMax,jbMax))
    fgvphi=0.0_r8
    ALLOCATE (fgtphim(ibMax,kMax,jbMax))
    fgtphim=0.0_r8
    ALLOCATE (fgtphiC(ibMax,kMax,jbMax))
    fgtphiC=0.0_r8
    ALLOCATE (fgtlam(ibMax,kMax,jbMax))
    fgtlam=0.0_r8
    ALLOCATE (fgqlam(ibMax,kMax,jbMax))
    fgqlam=0.0_r8
    ALLOCATE (fgulam(ibMax,kMax,jbMax))
    fgulam=0.0_r8
    ALLOCATE (fgvlam(ibMax,kMax,jbMax))
    fgvlam=0.0_r8
    ALLOCATE (fgtlamm(ibMax,kMax,jbMax))
    fgtlamm=0.0_r8
    ALLOCATE (fgtlamC(ibMax,kMax,jbMax))
    fgtlamC=0.0_r8
    IF (SL_twotime_scheme) THEN
       ALLOCATE (fgumean(ibMax,     jbMax))
       fgumean=0.0_r8
       ALLOCATE (fgvmean(ibMax,     jbMax))
       fgvmean=0.0_r8
       ALLOCATE (fgpsp(ibMax,     jbMax_ext))
       fgpsp=0.0_r8
       ALLOCATE (fgvdlnpm(ibMax,     jbMax))
       fgvdlnpm=0.0_r8
      ELSE
       ALLOCATE (fgumean(ibMax,     jbMax_ext))
       fgumean=0.0_r8
       ALLOCATE (fgvmean(ibMax,     jbMax_ext))
       fgvmean=0.0_r8
       ALLOCATE (fgpsp(ibMax,     jbMax))
       fgpsp=0.0_r8
       ALLOCATE (fgvdlnpm(ibMax,     jbMax_ext))
       fgvdlnpm=0.0_r8
    ENDIF 
    ALLOCATE (fgvdlnp(ibMax,     jbMax))
    fgvdlnp=0.0_r8
    ALLOCATE (fglnps(ibMax,     jbMax))
    fglnps=0.0_r8
    if (grid_difus) then 
       ALLOCATE (fglnps_nabla4(ibMax,     jbMax))
       fglnps_nabla4=0.0_r8
    endif
    ALLOCATE (fglnpsC(ibMax,     jbMax))
    fglnpsC=0.0_r8
    ALLOCATE (fgps(ibMax,     jbMax))
    fgps=0.0_r8
    ALLOCATE (fglnpm(ibMax,     jbMax))
    fglnpm=0.0_r8
    ALLOCATE (fgpphi(ibMax,     jbMax))
    fgpphi=0.0_r8
    ALLOCATE (fgplam(ibMax,     jbMax))
    fgplam=0.0_r8
    ALLOCATE (fgplamm(ibMax,     jbMax))
    fgplamm=0.0_r8
    ALLOCATE (fgplamC(ibMax,     jbMax))
    fgplamC=0.0_r8
    ALLOCATE (fgpphim(ibMax,     jbMax))
    fgpphim=0.0_r8
    ALLOCATE (fgpphiC(ibMax,     jbMax))
    fgpphiC=0.0_r8
    ALLOCATE (fgzs(ibMax,     jbMax))
    fgzs=0.0_r8
    ALLOCATE (fgzslam(ibMax,     jbMax))
    fgzslam=0.0_r8
    ALLOCATE (fgzsphi(ibMax,     jbMax))
    fgzsphi=0.0_r8









  END SUBROUTINE InitFields
END MODULE FieldsDynamics
