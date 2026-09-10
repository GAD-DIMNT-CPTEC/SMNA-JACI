!
!  $Author: pkubota $
!  $Date: 2009/03/03 16:36:38 $
!  $Revision: 1.16 $
!
!  Modified by Tarassova 2015
!  Climate aerosol (Kinne, 2013) coarse mode included
!  Fine aerosol mode 2000 is included
!  Modifications (3) are marked by 
!  !tar begin and  !tar end
!
MODULE NonLinearNMI
  USE Constants, ONLY: &
       gasr, grav, tref, er, pai, eriv, twomg, pscons,i8,r8,r4
  USE Utils, ONLY: &
       Rg,         &
       Tql2,       &
       Tred2,      &
       IJtoIBJB,   &
       LinearIJtoIBJB, &
       NearestIJtoIBJB, &
       SeaMaskIJtoIBJB, &
       SplineIJtoIBJB, &
       AveBoxIJtoIBJB

  USE SpecDynamics, ONLY: tm, sv, am, bm, cm, p1, p2, h1, h2, &
       Bmcm, dk, hm, snnp1_si
  USE IOLowLevel, ONLY: ReadGetNFTGZ
  USE ModTimeStep, ONLY: TimeStep, SfcGeoTrans
  USE Sizes, ONLY: kmax, nmax, mmax, imax, jmax, mymnMap, &
       ibmax ,   &
       jbmax,    &
       kmaxloc,  &
       a_hybr,   &
       b_hybr,   &
       mymnMax,  &
       mnmax,    &
       mnmax_si, &
       mymMax,   &
       lm2m,     &
       nodeHasM, &
       havesurf, &
       nsends_si,&
       nrecs_si, &
       maps_si,  &
       mapr_si,  &
       mysends_si,  &
       myrecs_si,   &
       inibr_si,    &
       inibs_si,    &
       map_four,    &
       ngroups_four,&
       nlevperg_four,&
       ThreadDecomp,&
       ThreadDecompms,&
       ibMaxPerJB

  USE FieldsDynamics, ONLY: qdivt, qdivp, &
       qrott, qrotp, &
       qtmpt, qtmpp, &
       qlnpt, qlnpp,   qlnpl, &
       qrott_si, qdivt_si, &
       qtmpt_si, qlnpt_si, &
       qqp,   qdiaten, qgzs
  USE FieldsPhysics, ONLY: &
       tg1   , &
       tg2   , &
       tg3   , &
       zorl  , &
       AlbVisDiff , &
       gtsea , &
       gco2flx,&
       gndvi , &
       geshem, &
       sheleg, &
       soilm , &
       wsib3d, &
       o3mix, &
       tracermix,&
!tar begin 
!climate aerosol parameters  coarse mode     
       aod, &
       asy, &
       ssa, &
       z_aer, &
!tar end      
!
!tar begin 
!climate aerosol parameters fine mode 2000     
       aodF, &
       asyF, &
       ssaF, &
       z_aerF
!tar end  

  USE Options, ONLY: &
       delt  ,&
       nfprt ,&
       nstep ,&
       nfin0 ,&
       nfin1 ,&
       nfnmi ,&
       ifalb ,&
       ifsst ,ifco2flx,&
       ifslm ,&
       ifsnw ,&
       ifozone, &
!tar begin       
       ifaeros, & 
!tar end      
       iftracer,&
       intsoilm,&
       sstlag,&
       intsst,&
       fint  ,&
       nftgz0,&
       nfzol ,& 
       dt    ,&
       percut,&
       jdt   ,&
       yrl   ,&
       monl  ,&
       dodyn ,&
       nfdyn ,&
       nfsst ,nfco2fx, &
       filta ,&
       filtb ,&
       ndord ,&
       istrt ,&
       ifilt ,&
       kt    ,&
       ktm   ,&
       slhum ,&
       ifndvi,&
       ifslmSib2,&
       intndvi,&
       SL_twotime_scheme,&
       reducedGrid,&
       slhum


  USE InputOutput, ONLY: &
       getsbc,&
       gread ,&
       gread4,&
       fsbc

  USE Griddynamics, ONLY:      &
       init_globconserv,       &
       init_globfluxconserv

  USE Communications, ONLY: &
       SpectoSi,            &
       SitoSpec

  USE Parallelism, ONLY: &
       MsgDump,          &
       Msgone ,          &
       DestroyParallelism, &
       FatalError,       &
       mygroup_four,     &
       maxNodes,         &
       myid_four,        &
       myid  

!  USE Dumpgraph, ONLY: &
!       dumpgra

   IMPLICIT NONE
  SAVE       

  PRIVATE
  PUBLIC :: Nlnmi
  PUBLIC :: Diaten
  PUBLIC :: Getmod

  INTEGER :: mods, niter
  INTEGER, ALLOCATABLE  :: nmodperg(:)
  INTEGER, ALLOCATABLE  :: grouphasmod(:)
  INTEGER :: modsloc
  INTEGER :: myfirstmod

  INCLUDE "mpif.h"

CONTAINS


  !
  ! diaten :computation of the diabatic terms for normal mode
  !         initialization.
  !


  SUBROUTINE Diaten(slagr,fName0,&
       ifday, tod, idate, idatec)

    LOGICAL,             INTENT(IN) :: slagr
    CHARACTER(LEN=*),  INTENT(IN) :: fName0
    INTEGER, INTENT(OUT) :: ifday
    REAL(KIND=r8)   , INTENT(OUT) :: tod
    INTEGER, INTENT(OUT), DIMENSION(4) :: idate, idatec
    !
    INTEGER :: ierr
    REAL(KIND=r8)    :: fa, fb, fb1, delta2
    REAL(KIND=r8)    :: tice=271.16e0_r8
    REAL(KIND=r8) ::   buf (iMax,jMax,4)
    REAL(KIND=r4) ::   brf (iMax,jMax)
    CHARACTER(LEN=8) :: c0
    CHARACTER(LEN=*), PARAMETER :: h="**(Diaten)**"
    INTEGER :: irec
    LOGICAL :: dohum, bckhum,dotrac

    REAL(KIND=r8)       :: qspec(2*mnmax,kmax)

    buf=0.0_r8
    fsbc=.false.
    if (myid.eq.0) write(*,*) ' in diaten '
    !
    !     start reading initial values
    !
    OPEN(UNIT=nfin0, FILE=TRIM(fName0), FORM='unformatted',ACCESS='sequential',&
         ACTION='read',STATUS='old', IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(c0,"(i8)") ierr
       CALL FatalError(h//" Open file "//TRIM(fName0)//" returned IOSTAT="//&
            TRIM(ADJUSTL(c0)))
       STOP
    END IF
    !
    CALL gread4 (nfin0, ifday, tod  , idate, idatec,qgzs  ,qlnpp , &
         qtmpp, qdivp, qrotp, qqp, a_hybr, b_hybr, dodyn, &
         nfdyn)
    REWIND nfin0
    !
    !     cold start: reset precip. to zero.
    !
    geshem = 0.0_r8
    !
    !     calculates laplacian of topography
    !
    CALL SfcGeoTrans(slagr)
    !
    !     read climatology data
    !
    buf (1:iMax,1:jMax,1:4)=0.0_r8
    irec=1
    CALL ReadGetNFTGZ(nftgz0,irec,buf(1:iMax,1:jMax,1),buf(1:iMax,1:jMax,2),buf(1:iMax,1:jMax,3))
    READ (UNIT=nfzol, REC=1) brf
    buf(1:iMax,1:jMax,4)=brf(1:iMax,1:jMax)

    IF (reducedGrid) THEN
       CALL LinearIJtoIBJB(buf(1:iMax,1:jMax,1),tg1)
        !CALL AveBoxIJtoIBJB(buf(1:iMax,1:jMax,1),tg1)
    ELSE
       CALL IJtoIBJB(buf(1:iMax,1:jMax,1) ,tg1 )
    END IF

    IF (reducedGrid) THEN
       CALL LinearIJtoIBJB(buf(1:iMax,1:jMax,2) ,tg2)
       !CALL AveBoxIJtoIBJB(buf(1:iMax,1:jMax,2),tg2)
    ELSE
       CALL IJtoIBJB(buf(1:iMax,1:jMax,2) ,tg2 )
    END IF

    IF (reducedGrid) THEN
       CALL LinearIJtoIBJB(buf(1:iMax,1:jMax,3) ,tg3)
       !CALL AveBoxIJtoIBJB(buf(1:iMax,1:jMax,3),tg3)
    ELSE
       CALL IJtoIBJB(buf(1:iMax,1:jMax,3) ,tg3 )
    END IF

    IF (reducedGrid) THEN
       CALL LinearIJtoIBJB(buf(1:iMax,1:jMax,4),zorl)
       !CALL AveBoxIJtoIBJB(buf(1:iMax,1:jMax,4),zorl)
    ELSE
       CALL IJtoIBJB(buf(1:iMax,1:jMax,4),zorl )
    END IF


    CALL getsbc (iMax ,jMax  ,kMax, AlbVisDiff,gtsea,gco2flx,gndvi,&
         soilm,sheleg,o3mix,tracermix,wsib3d,&
!tar begin
! climate aerosol parameters of coarse mode
         aod, asy, ssa, z_aer, ifaeros, & 
! climate aerosol parameters of fine mode
         aodF, asyF, ssaF, z_aerF, &
!tar end
         ifday , tod  ,idate ,idatec, &
         ifalb,ifsst,ifco2flx,ifndvi,ifslm ,ifslmSib2,&
         ifsnw,ifozone,iftracer, &
         sstlag,intsst,intndvi,intsoilm,fint ,tice  , &
         yrl  ,monl,ibMax,jbMax,ibMaxPerJB)
    !
    !     cold start (at first delt/4 ,then delt/2 )
    !
    dt= delt /4.0_r8

    ! filter arguments for first time step

    fa = 0.0_r8
    fb = 1.0_r8
    fb1 = 1.0_r8
    ifilt=0
    bckhum = .TRUE.
    dohum = .not.slhum
    dotrac = dohum
    init_globconserv = .TRUE.
    init_globfluxconserv= .TRUE.

    DO jdt=1,2
       istrt=jdt
       !
       !     calculate matrices for semi-implicit integration
       !
       delta2 = dt
       IF(slagr.and.SL_twotime_scheme) delta2 = dt/2._r8

       CALL bmcm(delta2)
       !      perform time step
       !$OMP PARALLEL
       CALL TimeStep(fb1,fa,fb,dotrac,slagr,slhum,dohum,bckhum,.FALSE., &
       .FALSE.,.FALSE.,.FALSE.,dt,jdt,ifday,tod,idatec,jdt)
       !$OMP END PARALLEL
!      call dumpgra


       ! prepare next time step, including filter arguments

       dt=dt*2.0_r8
       ktm=kt
       fb1 = 0.0_r8
       init_globconserv = .FALSE.
       init_globfluxconserv= .FALSE.
       IF (slhum) bckhum = .FALSE.
    END DO
    !
    !     smooth start
    !
    istrt=0
    !
    ! semi-implicit matrices
    !
    delta2 = dt
    IF(slagr.and.SL_twotime_scheme) delta2 = dt/2._r8
    CALL bmcm(delta2)
    !
    ! filter arguments for all remaining time steps
    !
    fa = filta
    fb = filtb
    fb1 = 0.0_r8
    IF(slagr.and.SL_twotime_scheme) THEN
       fa = 1.0_r8
       fb = 0.0_r8
    END IF
    ifilt=1
    !
    ! time step loop
    !
    DO jdt=1,nstep
       !
       !     step loop starts
       !
       ! perform time step
       !$OMP PARALLEL
       CALL TimeStep(fb1,fa,fb,dotrac,slagr,slhum,dohum,bckhum,.FALSE., &
                     .FALSE.,.TRUE.,.FALSE.,dt,jdt,ifday,tod,idatec,jdt+1)
       !$OMP END PARALLEL
!      call dumpgra
!      CALL DestroyParallelism("***MODEL EXECUTION ENDS NORMALY***")
!   stop

       ktm=kt
       !
       fb1 = fb
       IF(slagr.and.SL_twotime_scheme) fb1 = 0.0_r8
    ENDDO
    qdiaten = qdiaten / (2.0_r8*dt*nstep)

    !
  END SUBROUTINE Diaten


  ! nlnmi :  nonlinear normal mode
  !          initialization of the model,
  !          possibly diabatic. Initial
  !          values of vorticity, divergence,
  !          temperature and log of surface
  !          pressure are adjusted to supress
  !          gravity modes, according to the
  !          normal mode method of Machenhauer.



  SUBROUTINE Nlnmi(nlnminit,diabatic,slagr,fName, & 
       ifday, tod, idatec,ktm)
    LOGICAL,          INTENT(IN) :: nlnminit
    LOGICAL,          INTENT(IN) :: diabatic
    LOGICAL,          INTENT(IN) :: slagr
    CHARACTER(LEN=*), INTENT(IN) :: fName
    INTEGER, INTENT(IN) :: ifday
    REAL(KIND=r8)   , INTENT(IN) :: tod
    INTEGER, INTENT(IN), DIMENSION(4) :: idatec
    INTEGER, INTENT(INOUT) :: ktm
    !
    ! vertical mode arrays.
    ! "RealBuffer" is a scratch area that packs 
    ! all vertical mode arrays, for MPI communications and IO.
    !
    REAL(KIND=r8)                 :: eigg(kmax,kmax)
    REAL(KIND=r8)                 :: eiggt(kmax,kmax)
    REAL(KIND=r8)                 :: dotpro(kmax)
    REAL(KIND=r8)                 :: gh(kmax)
    REAL(KIND=r8)                 :: verin(kmax)

    REAL(KIND=r8),    ALLOCATABLE :: RealBuffer(:)
    !
    ! "g" contains base functions for horizontal modes,
    ! storing all simmetric and anti-simmetric base functions 
    ! for each normal mode and wave number.
    ! Since number of base functions vary with wave number,
    ! all rank 3 matrices of base functions indexed by wave number,
    ! simmetry and mode are packed into the rank 1 array "g".
    !
    ! "per" contains the period for each wave and horizontal mode.
    ! Again, all rank 3 matrices of period indexed by wave number,
    ! simmetry and mode are packed into the rank 1 array "per"
    !
    REAL(KIND=r8),    ALLOCATABLE :: g(:)
    REAL(KIND=r8),    ALLOCATABLE :: per(:)
    !
    ! "indper" points to the first element on "per" of each packed
    ! array as a function of wave number, simmetry and mode. 
    ! The array has size "jg", also indexed by wave number, simmetry and mode.
    ! Wave number and simmetry are packed in the first dimension.
    !
    INTEGER,          ALLOCATABLE :: indper(:,:)
    !
    ! "indg" points to the first element on "g" of each packed
    ! array as a function of wave number, simmetry and mode. 
    ! The array has size "jg" and "nas", also indexed by 
    ! wave number, simmetry and mode.
    ! Wave number and simmetry are packed in the first dimension.
    !
    INTEGER,          ALLOCATABLE :: indg(:,:)
    !
    ! "jg" is the first dimension of each packed array into
    ! "g" and "per".
    ! "nas" is the second dimension of each packed array into
    ! "g".
    !
    INTEGER(KIND=i8), ALLOCATABLE :: jg(:,:)
    INTEGER(KIND=i8), ALLOCATABLE :: nas(:,:)
    !
    ! "one_per" and "one_g" are scratch areas large enough to
    ! store any of the arrays packed into "per" and "g"
    !
    REAL(KIND=r8),    ALLOCATABLE :: one_per(:)
    REAL(KIND=r8),    ALLOCATABLE :: one_g(:)
    !
    ! Domain decomposition: arrays "g" and "per" are distibuted
    ! over MPI processes. Each process stores only its own wave numbers.
    ! Consequently, "jg" and "nas" store dimensions of own wave numbers.
    ! Arrays "all_jg" and "all_nas" store dimensions of all wave numbers,
    ! not only these own by each process.
    !
    INTEGER,          ALLOCATABLE :: all_jg(:,:)
    INTEGER,          ALLOCATABLE :: all_nas(:,:)
    !
    ! "mine_per_g" indexed by wave number is true iff current
    ! process stores the wave numbers (and, consequently,
    ! arrays for this wave number packed into "per" and "g")
    !
    LOGICAL                       :: mine_per_g(mMax)
    !
    ! "tag_per" and "tag_g" are MPI communication tags
    !
    INTEGER,            PARAMETER :: tag_per=43
    INTEGER,            PARAMETER :: tag_g=44
    !
    ! indices and other scratch areas
    !
    INTEGER           :: i
    INTEGER           :: ig
    INTEGER           :: ijg
    INTEGER           :: iper
    INTEGER           :: iter
    INTEGER           :: k
    INTEGER           :: j
    INTEGER           :: l
    INTEGER           :: kmod
    INTEGER           :: mode
    INTEGER           :: group
    INTEGER           :: proc
    INTEGER           :: lev
    INTEGER           :: ll
    INTEGER           :: ierr
    INTEGER           :: ip
    INTEGER           :: kp
    INTEGER           :: twice
    INTEGER           :: i0
    INTEGER           :: i1
    INTEGER           :: i2
    INTEGER           :: myms(mymmax)
    INTEGER           :: nms
    INTEGER           :: mnFirst_si
    INTEGER           :: mnLast_si
    INTEGER           :: mnRIFirst_si
    INTEGER           :: mnRILast_si
    INTEGER           :: mnRIFirst
    INTEGER           :: mnRILast
    INTEGER           :: status(MPI_STATUS_SIZE)
    INTEGER           :: ind1
    INTEGER(KIND=i8)  :: scratch1
    INTEGER(KIND=i8)  :: scratch2
    INTEGER(KIND=i8)  :: tot_g
    INTEGER(KIND=i8)  :: tot_per
    REAL(KIND=r8)     :: dt
    REAL(KIND=r8)     :: tor
    REAL(KIND=r8)     :: to
    REAL(KIND=r8)     :: fa
    REAL(KIND=r8)     :: fb
    REAL(KIND=r8)     :: fb1
    REAL(KIND=r8)     :: fcon
    REAL(KIND=r8),  ALLOCATABLE :: qgenp_si(:,:)
    REAL(KIND=r8),  ALLOCATABLE :: qgenp(:,:)
    CHARACTER(LEN=10) :: c0
    CHARACTER(LEN=10) :: c1
    CHARACTER(LEN=*), PARAMETER :: h="**(Nlnmi)**"

    REAL(KIND=r8)       :: qspec(2*mnmax,kmax)
    myms=0
    ! Read file nfnmi twice.
    ! First read to store sizes of arrays "g" and "per"
    ! into arrays "jg" and "nas".
    ! Second time to read arrays "g" and "per".

    ! allocate "RealBuffer" (buffer of vertical modes)
    ! since "RealBuffer" packs eigg(kMax,kMax), eiggt(kMax,kMax), 
    ! gh(kMax), dotpro(kMax) and tref , it has size kMax*(2*kMax+2)+1

    ALLOCATE(RealBuffer(kMax*(2*kMax+2)+1), STAT=ierr);RealBuffer=0.0_r8
    IF (ierr /= 0) THEN
       WRITE(c0,"(i10)") ierr
       CALL FatalError(h//" Allocate RealBuffer returned STAT="//&
            TRIM(ADJUSTL(c0)))
       STOP
    END IF

    ! process 0 reads vertical modes into "RealBuffer"

    IF (myId == 0) THEN
       OPEN(UNIT=nfnmi, FILE=TRIM(fName), FORM='unformatted',ACCESS='sequential',&
            ACTION='read', STATUS='old', IOSTAT=ierr) 
       IF (ierr /= 0) THEN
          WRITE(c0,"(i10)") ierr
          CALL FatalError(h//" Open file "//TRIM(fName)//" returned IOSTAT="//&
               TRIM(ADJUSTL(c0)))
          STOP
       END IF
       READ (UNIT=nfnmi, IOSTAT=ierr) RealBuffer(:)
       IF (ierr /= 0) THEN
          WRITE(c0,"(i10)") ierr
          CALL FatalError(h//" Read RealBuffer returned IOSTAT="//TRIM(ADJUSTL(c0)))
          STOP
       END IF
    END IF

    ! process 0 broadcasts "RealBuffer"

    IF (maxNodes > 1) THEN
       CALL MPI_BCAST(RealBuffer(1), SIZE(RealBuffer), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
       IF (ierr /= MPI_SUCCESS) THEN
          WRITE(c0,"(i10)") ierr
          CALL FatalError(h//" Broadcast RealBuffer returned ierr="//TRIM(ADJUSTL(c0)))
          STOP
       END IF
    END IF

    ! all processes unpack "RealBuffer"

    i0 = 0
    DO i2 = 1, kMax
       DO i1 = 1, kMax
          eigg(i1,i2) = RealBuffer(i0+i1+(i2-1)*kMax)
       END DO
    END DO
    i0 = i0 + kMax*kMax
    DO i2 = 1, kMax
       DO i1 = 1, kMax
          eiggt(i1,i2) = RealBuffer(i0+i1+(i2-1)*kMax)
       END DO
    END DO
    i0 = i0 + kMax*kMax
    DO i1 = 1, kMax
       gh(i1) = RealBuffer(i0+i1)
    END DO
    i0 = i0 + kMax
    DO i1 = 1, kMax
       dotpro(i1) = RealBuffer(i0+i1)
    END DO
    i0 = i0 + kMax
    to = RealBuffer(i0+1)

    ! how many modes

    CALL SetMods(gh)
    ALLOCATE (qgenp_si(2*mnMax_si,kmax));qgenp_si=0.0_r8
    ALLOCATE (qgenp(2*mymnMax,modsloc));qgenp=0.0_r8

    IF (mods <= 2) THEN
       niter=2
    ELSE
       niter=3
    END IF

    CALL bmcm(1.0_r8)

    ! vertical integration

    DO k=1, kmax
       verin(k)=0.0_r8
       DO j=1, kmax
          verin(k)=verin(k)-sv(j)*eigg(j,k)
       END DO
    END DO

    ! allocate buffers to continue reading

    ALLOCATE(all_jg(2*mMax,mods), STAT=ierr);all_jg=0
    IF (ierr /= 0) THEN
       WRITE(c0,"(i10)") ierr
       CALL FatalError(h//" Allocate all_jg returned STAT="//&
            TRIM(ADJUSTL(c0)))
       STOP
    END IF
    ALLOCATE(all_nas(2*mMax,mods), STAT=ierr);all_nas=0
    IF (ierr /= 0) THEN
       WRITE(c0,"(i10)") ierr
       CALL FatalError(h//" Allocate all_nas returned STAT="//&
            TRIM(ADJUSTL(c0)))
       STOP
    END IF
    ALLOCATE(one_g(3*(mMax+1)*(mMax+1)), STAT=ierr);one_g=0.0_r8
    IF (ierr /= 0) THEN
       WRITE(c0,"(i10)") ierr
       CALL FatalError(h//" Allocate one_g returned STAT="//&
            TRIM(ADJUSTL(c0)))
       STOP
    END IF
    ALLOCATE(one_per(3*(mMax+1)*(mMax+1)), STAT=ierr);one_per=0.0_r8
    IF (ierr /= 0) THEN
       WRITE(c0,"(i10)") ierr
       CALL FatalError(h//" Allocate one_per returned STAT="//&
            TRIM(ADJUSTL(c0)))
       STOP
    END IF

    ! process 0 reads array sizes and arrays (to scratch)

    IF (myId == 0) THEN

       ! capture sizes and read arrays into scratch

       DO i=1,mods 
          DO ll=1,mMax
             DO twice = 1, 2
                ind1 = 2*(ll-1)+twice
                READ(UNIT=nfnmi, IOSTAT=ierr) scratch1, scratch2
                all_jg(ind1,i) = INT(scratch1)
                all_nas(ind1,i)= INT(scratch2)
                IF (ierr /= 0) THEN
                   WRITE(c0,"(i10)") ierr
                   CALL FatalError(h//" Read all_jg, all_nas returned IOSTAT="//TRIM(ADJUSTL(c0)))
                   STOP
                ELSE IF (all_jg(ind1,i) /= 0_i8) THEN
                   CALL getperg(scratch2, one_per, one_g, scratch1)
                END IF
             END DO
          END DO
       END DO

       ! reset file position for next reading

       REWIND(UNIT=nfnmi)
    END IF

    ! broadcast all_jg, all_nas

    IF (maxNodes > 1) THEN
       CALL MPI_BCAST(all_jg(1,1), SIZE(all_jg), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
       IF (ierr /= MPI_SUCCESS) THEN
          WRITE(c0,"(i10)") ierr
          CALL FatalError(h//" Broadcast all_jg returned ierr="//TRIM(ADJUSTL(c0)))
          STOP
       END IF
       CALL MPI_BCAST(all_nas(1,1), SIZE(all_nas), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
       IF (ierr /= MPI_SUCCESS) THEN
          WRITE(c0,"(i10)") ierr
          CALL FatalError(h//" Broadcast all_nas returned ierr="//TRIM(ADJUSTL(c0)))
          STOP
       END IF
    END IF

    ! allocate indices and sizes at this process

    ALLOCATE(indper(2*mymMax,modsloc), STAT=ierr);indper=-1
    IF (ierr /= 0) THEN
       WRITE(c0,"(i10)") ierr
       CALL FatalError(h//" Allocate indper returned STAT="//&
            TRIM(ADJUSTL(c0)))
       STOP
    END IF
    ALLOCATE(indg(2*mymMax,modsloc), STAT=ierr);indg=-1
    IF (ierr /= 0) THEN
       WRITE(c0,"(i10)") ierr
       CALL FatalError(h//" Allocate indg returned STAT="//&
            TRIM(ADJUSTL(c0)))
       STOP
    END IF
    ALLOCATE(jg(2*mymMax,modsloc), STAT=ierr);jg=-1
    IF (ierr /= 0) THEN
       WRITE(c0,"(i10)") ierr
       CALL FatalError(h//" Allocate jg returned STAT="//&
            TRIM(ADJUSTL(c0)))
       STOP
    END IF
    ALLOCATE(nas(2*mymMax,modsloc), STAT=ierr);nas=-1
    IF (ierr /= 0) THEN
       WRITE(c0,"(i10)") ierr
       CALL FatalError(h//" Allocate nas returned STAT="//&
            TRIM(ADJUSTL(c0)))
       STOP
    END IF

    ! domain decomposition

    DO ll=1,mMax
       mine_per_g(ll) = nodeHasM(ll,mygroup_four) == myId_four
    END DO

    ! compute indices and sizes at this process

    ig = 1
    iper = 1
    DO i=1,modsloc
       mode = myfirstmod + i - 1
       ijg = 1
       DO ll=1,Mmax
          DO twice = 1, 2
             ind1 = 2*(ll-1)+twice
             IF (mine_per_g(ll)) THEN
                jg(ijg,i)  = all_jg(ind1,mode)
                nas(ijg,i) = all_nas(ind1,mode)
                indper(ijg,i) = iper
                iper = iper + jg(ijg,i)
                indg(ijg,i) = ig
                ig = ig + jg(ijg,i) * nas(ijg,i)
                ijg = ijg + 1
             END IF
          END DO
       END DO

       ! check correction

       IF (ijg /= 2*mymMax+1) THEN
          WRITE(c0,"(i10)") i
          CALL MsgDump(h," jg and nas not fully filled at mode "//TRIM(ADJUSTL(c0)))
          CALL FatalError(h//" jg and nas not fully filled at mode "//TRIM(ADJUSTL(c0)))
       END IF
    END DO

    ! total sizes of arrays "g" and "per"

    tot_g = ig - 1
    tot_per = iper - 1

    ! allocate arrays "g" and "per"

    ALLOCATE(per(tot_per), STAT=ierr);per=0.0_r8
    IF (ierr /= 0) THEN
       WRITE(c0,"(i10)") ierr
       WRITE(c1,"(i10)") tot_per
       CALL FatalError(h//" Allocate per of size "//&
            TRIM(ADJUSTL(c1))//" returned STAT="//&
            TRIM(ADJUSTL(c0)))
       STOP
    END IF
    ALLOCATE(g(tot_g), STAT=ierr);g=0.0_r8
    IF (ierr /= 0) THEN
       WRITE(c0,"(i10)") ierr
       WRITE(c1,"(i10)") tot_g
       CALL FatalError(h//" Allocate g of size "//&
            TRIM(ADJUSTL(c1))//" returned STAT="//&
            TRIM(ADJUSTL(c0)))
       STOP
    END IF

    !  Read file header by the second time (to scratch)

    IF (myId == 0) THEN
       READ (UNIT=nfnmi, IOSTAT=ierr) RealBuffer(:)
       IF (ierr /= 0) THEN
          WRITE(c0,"(i10)") ierr
          CALL FatalError(h//" Second read RealBuffer returned IOSTAT="//TRIM(ADJUSTL(c0)))
          STOP
       END IF
    END IF

    ! process 0 reads sizes into scratch, reads horizontal modes 
    ! into scratch and communicates horizontal modes to
    ! owner process; owner stores data into appropriate places
    !
    ! "ijg" count wave numbers and simmetry for this process

    DO i=1,mods 
       group = grouphasmod(i)
       i1 = i-myfirstmod+1
       ijg = 1
       DO ll=1,Mmax
          DO twice = 1, 2
             ind1 = 2*(ll-1)+twice

             ! process 0 reads data into one_per and one_g
             ! and sends data to owner

             IF (myId == 0) THEN
                READ(UNIT=nfnmi, IOSTAT=ierr) scratch1, scratch2
                IF (ierr /= 0) THEN
                   WRITE(c0,"(i10)") ierr
                   CALL FatalError(h//" Read scratch1, scratch2 returned IOSTAT="//TRIM(ADJUSTL(c0)))
                   STOP
                ELSE IF (scratch1 /= 0_i8) THEN
                   IF (group == mygroup_four.and.nodeHasM(ll,mygroup_four) == 0) THEN
                      CALL getperg(scratch2, per(indper(ijg,i1)), g(indg(ijg,i1)), scratch1)
                      ijg = ijg + 1
                   ELSE
                      proc = map_four(group,nodeHasM(ll,group))
                      CALL getperg(scratch2, one_per, one_g, scratch1)
                      CALL MPI_SEND(one_per(1), INT(scratch1), MPI_DOUBLE_PRECISION, &
                           proc, tag_per, MPI_COMM_WORLD, ierr)
                      IF (ierr /= MPI_SUCCESS) THEN
                         WRITE(c0,"(i10)") ierr
                         CALL FatalError(h//" Send one_per returned ierr="//TRIM(ADJUSTL(c0)))
                         STOP
                      END IF
                      CALL MPI_SEND(one_g(1), INT(scratch1*scratch2), MPI_DOUBLE_PRECISION, &
                           proc, tag_g, MPI_COMM_WORLD, ierr)
                      IF (ierr /= MPI_SUCCESS) THEN
                         WRITE(c0,"(i10)") ierr
                         CALL FatalError(h//" Send one_g returned ierr="//TRIM(ADJUSTL(c0)))
                         STOP
                      END IF
                   END IF
                END IF


             ! owner stores data at appropriate place

             ELSE IF (group == mygroup_four .AND. all_jg(ind1,i) /= 0_i8 .AND. nodeHasM(ll,mygroup_four) == myId_four) THEN
                CALL MPI_RECV(per(indper(ijg,i1)), INT(jg(ijg,i1)), MPI_DOUBLE_PRECISION, & 
                     0, tag_per, MPI_COMM_WORLD, status, ierr)
                IF (ierr /= MPI_SUCCESS) THEN
                   WRITE(c0,"(i10)") ierr
                   CALL FatalError(h//" Recv per returned ierr="//TRIM(ADJUSTL(c0)))
                   STOP
                END IF
                CALL MPI_RECV(g(indg(ijg,i1)), INT(jg(ijg,i1)*nas(ijg,i1)), MPI_DOUBLE_PRECISION, &
                     0, tag_g, MPI_COMM_WORLD, status, ierr)
                IF (ierr /= MPI_SUCCESS) THEN
                   WRITE(c0,"(i10)") ierr
                   CALL FatalError(h//" Recv one_g returned ierr="//TRIM(ADJUSTL(c0)))
                   STOP
                END IF
                ijg = ijg + 1
             END IF
          END DO
       END DO
    END DO
             

    IF (myId == 0) THEN
       CLOSE(nfnmi)
    END IF

    fb1 = 1.0_r8
    fa = 0.0_r8
    fb = 0.0_r8
    dt = 1.0_r8

    !$OMP PARALLEL PRIVATE(myms,nms,mnFirst_si,mnLast_si,mnRIFirst,mnRILast,iter,l,ip,kp,kmod,fcon,tor,lev,mnRIFirst_si,mnRILast_si)
    !
    !  Compute the iterations of the non-linear initialization
    !  -------------------------------------------------------

    CALL ThreadDecompms(mymMax, myms, nms)
    CALL ThreadDecomp(1, 2*mymnMax, mnRIFirst, mnRILast, "Nlnmi")
    CALL ThreadDecomp(1, mnMax_si, mnFirst_si, mnLast_si, "Nlnmi")
    CALL ThreadDecomp(1, 2*mnMax_si, mnRIFirst_si, mnRILast_si, "Nlnmi")
    !
    !  Compute the iterations of the non-linear initialization
    !  -------------------------------------------------------
    DO iter=1,niter
       !$OMP BARRIER
       !
       !  Compute non-linear complete tendencies
       !  --------------------------------------
       CALL TimeStep(fb1,fa,fb,.TRUE.,slagr,.FALSE.,.TRUE.,.TRUE.,.FALSE., &
                     nlnminit,.FALSE.,.FALSE.,dt,0,ifday, tod, idatec,0)
!      IF (iter.eq.1) call dumpgra
       !
       !     add diabatic heating rate to temperature tendency
       !     -------------------------------------------------
       IF (diabatic) THEN
          DO l=mnRIFirst, mnRILast
            qtmpt(l,:) = qtmpt(l,:) + qdiaten(l,:)
          END DO
       ENDIF
       !
       !    Transform fields to have all verticals in same processor
       !    --------------------------------------------------------
       !$OMP BARRIER
       !$OMP SINGLE
       CALL SpectoSi(inibs_si,inibr_si,nsends_si,nrecs_si,mysends_si,myrecs_si,&
                     maps_si,mapr_si,kmax,kmaxloc,nlevperg_four,qtmpt,qdivt,&
                     qtmpt_si,qdivt_si,qrott,qrott_si,qlnpt,qlnpt_si)
       !$OMP END SINGLE
       !
       !     create generalized pressure=(phi+r*to*q)(dot)
       !     ---------------------------------------------
       tor=tref * gasr
       DO l=1, kmax
          DO ip=mnRIFirst_si, mnRILast_si
             qgenp_si(ip,l)=tor*qlnpt_si(ip)
          END DO
          DO kp=1, kmax
             DO ip=mnRIFirst_si, mnRILast_si
                qgenp_si(ip,l)=qgenp_si(ip,l) + qtmpt_si(ip,kp)*hm(l,kp)
             END DO
          END DO
       END DO
       !$OMP BARRIER
       CALL vertic(qdivt_si,eigg,eiggt,dotpro,-1,mnFirst_si,mnLast_si)
       CALL vertic(qrott_si,eigg,eiggt,dotpro,-1,mnFirst_si,mnLast_si)
       CALL vertic(qgenp_si,eigg,eiggt,dotpro,-1,mnFirst_si,mnLast_si)
       CALL primes(qrott_si,qdivt_si,qgenp_si,gh,-1,mnFirst_si,mnLast_si)
       !$OMP BARRIER
       !$OMP SINGLE
       !
       !    Transform fields now collecting m's together, splitting modes
       !    -------------------------------------------------------------
       CALL SitoSpec(inibr_si,inibs_si,nrecs_si,nsends_si,myrecs_si,mysends_si,&
                     mapr_si,maps_si,mods,modsloc,nmodperg,.false.,qrott,qdivt,&
                     qrott_si,qdivt_si,qgenp,qgenp_si)
       !$OMP END SINGLE
       CALL horiz1(qrott,qdivt,qgenp,percut,per,g,indper,indg, &
                                           jg,nas,myms,nms)
       !
       !    Transform fields now collecting modes together, splitting m's
       !    -------------------------------------------------------------
       !$OMP BARRIER
       !$OMP SINGLE
       CALL SpectoSi(inibs_si,inibr_si,nsends_si,nrecs_si,mysends_si,myrecs_si,&
                     maps_si,mapr_si,mods,modsloc,nmodperg,qrott,qdivt,&
                     qrott_si,qdivt_si,qgenp,qgenp_si)
       !$OMP END SINGLE
       CALL primes(qrott_si,qdivt_si,qgenp_si,gh,+1,mnFirst_si,mnLast_si)
       !
       !     compute delta(q) from composite variable
       !
       DO ip=mnRIFirst_si, mnRILast_si
          qlnpt_si(ip)= 0.0_r8
       END DO
       !$OMP BARRIER
       DO kmod=1,mods
          fcon=verin(kmod)/gh(kmod)
          DO ip=mnRIFirst_si, mnRILast_si
             qlnpt_si(ip)=qlnpt_si(ip)-fcon*qgenp_si(ip,kmod)
          END DO
       END DO
       !$OMP BARRIER
       CALL vertic(qdivt_si,eigg,eiggt,dotpro,+1,mnFirst_si,mnLast_si)
       CALL vertic(qrott_si,eigg,eiggt,dotpro,+1,mnFirst_si,mnLast_si)
       CALL vertic(qgenp_si,eigg,eiggt,dotpro,+1,mnFirst_si,mnLast_si)
       !
       !     compute delta(phi) from composite variable
       !
       !$OMP BARRIER
       DO kp=1, kmax
          DO ip=mnRIFirst_si, mnRILast_si
             qgenp_si(ip,kp)=qgenp_si(ip,kp)-tor*qlnpt_si(ip)
          END DO
       END DO
       !
       !     compute delta(t) from phi
       !
       qtmpt_si(mnRIFirst_si:mnRILast_si,:) = 0.0_r8
       DO lev=1,kmax
          DO kp=1, kmax
             DO ip=mnRIFirst_si, mnRILast_si
                qtmpt_si(ip,lev)=qtmpt_si(ip,lev)+tm(lev,kp)*qgenp_si(ip,kp)
             END DO
          END DO
       END DO
       !
       !    Transform now back to normal spectral decomposition
       !    ---------------------------------------------------
       !$OMP BARRIER
       !$OMP SINGLE
       CALL SitoSpec(inibr_si,inibs_si,nrecs_si,nsends_si,myrecs_si,mysends_si,&
                     mapr_si,maps_si,kmax,kmaxloc,nlevperg_four,.false.,qrott,qdivt,&
                     qrott_si,qdivt_si,qtmpt,qtmpt_si,qlnpt,qlnpt_si)
       !$OMP END SINGLE

       IF (slagr.and.havesurf) THEN
          DO ip=mnRIFirst, mnRILast
             qlnpl(ip)=qlnpl(ip)-qlnpt(ip)
          END DO
       END IF
       DO ip=mnRIFirst, mnRILast
          qtmpp(ip,:)=qtmpp(ip,:)-qtmpt(ip,:)
          qrotp(ip,:)=qrotp(ip,:)-qrott(ip,:)
          qdivp(ip,:)=qdivp(ip,:)-qdivt(ip,:)
       END DO
       IF (havesurf) THEN
          DO ip=mnRIFirst, mnRILast
             qlnpp(ip)=qlnpp(ip)-qlnpt(ip)
          END DO
       END IF
    END DO
    !$OMP END PARALLEL
    ktm=0

    DEALLOCATE(RealBuffer)
    DEALLOCATE(qgenp)
    DEALLOCATE(qgenp_si)
    DEALLOCATE(indper)
    DEALLOCATE(indg)
    DEALLOCATE(jg)
    DEALLOCATE(nas)
    DEALLOCATE(per)
    DEALLOCATE(g)

    DEALLOCATE(one_g)
    DEALLOCATE(one_per)
    DEALLOCATE(all_jg)
    DEALLOCATE(all_nas)
  END SUBROUTINE Nlnmi




  ! vertic : performs projection of the spectral representation of all
  !          model levels of a field onto the vertical normal modes.
  !          also performs expansion of the spectral representation
  !          of a vertically projected field into a field at all model
  !          levels using the vertical normal modes.
  !
  !     input=-1 to obtain vertical mode expansion
  !     input=+1 to obtain spet. coefs. from vertical expansion



  SUBROUTINE vertic(f,eigg,eiggt,dotpro,input,mnFirst,mnLast)
    REAL(KIND=r8), INTENT(IN) :: eigg(kmax,kmax)
    REAL(KIND=r8), INTENT(IN) :: eiggt(kmax,kmax)
    REAL(KIND=r8), INTENT(IN) :: dotpro(kmax)
    REAL(KIND=r8), INTENT(INOUT) :: f(2,mnMax_si,kmax)
    INTEGER, INTENT(IN) :: input
    INTEGER, INTENT(IN) :: mnFirst
    INTEGER, INTENT(IN) :: mnLast
    REAL(KIND=r8) :: col(2,kmax)
    INTEGER :: mn, kmod, lev
    REAL(KIND=r8) :: sum1, sum2

    IF (input .lt. 0) THEN
       DO mn=mnFirst,mnLast
          DO kmod=1,mods
             sum1=0.0_r8
             sum2=0.0_r8
             DO lev=1, kmax
                sum1=sum1+eiggt(lev,kmod)* f(1,mn,lev)
                sum2=sum2+eiggt(lev,kmod)* f(2,mn,lev)
             END DO
             col(1,kmod)=dotpro(kmod)*sum1
             col(2,kmod)=dotpro(kmod)*sum2
          END DO
          DO kmod=1,mods
             f(1,mn,kmod)=col(1,kmod)
             f(2,mn,kmod)=col(2,kmod)
          END DO
       END DO
    ELSE IF (input .gt. 0) THEN
       DO mn=mnFirst,mnLast
          DO kmod=1,mods
             col(1,kmod)=f(1,mn,kmod)
             col(2,kmod)=f(2,mn,kmod)
          END DO
          DO lev=1, kmax
             f(1,mn,lev)=0.0_r8
             f(2,mn,lev)=0.0_r8
             DO kmod=1,mods
                f(1,mn,lev)=f(1,mn,lev)+eigg(lev,kmod)*col(1,kmod)
                f(2,mn,lev)=f(2,mn,lev)+eigg(lev,kmod)*col(2,kmod)
             END DO
          END DO
       END DO
    END IF
  END SUBROUTINE vertic



  ! primes : performs scaling of the vertically projected tendencies
  !          of vorticity, divergence, and composite mass variable
  !          to convert them to the correct form for the calculation
  !          of the adjustment to these fields before so doing during
  !          the current iteration of the machenauer initialization
  !          technique.
  !          also performs descaling of the vertically projected
  !          tendencies of vorticity, divergence, and composite mass
  !          variable to convert them from the correct form for the
  !          calculation of the adjustment to these fields after
  !          so doing during the current iteration of the Machenhauer
  !          initialization technique.



  SUBROUTINE primes(vord,divd,comd,gh,input,mnFirst,mnLast)
    REAL(KIND=r8), INTENT(INOUT) :: vord(2,mnMax_si,kmax)
    REAL(KIND=r8), INTENT(INOUT) :: divd(2,mnMax_si,kmax)
    REAL(KIND=r8), INTENT(INOUT) :: comd(2,mnMax_si,kmax)
    REAL(KIND=r8), INTENT(IN) :: gh(kmax)
    INTEGER, INTENT(IN) :: input
    INTEGER, INTENT(IN) :: mnFirst
    INTEGER, INTENT(IN) :: mnLast
    REAL(KIND=r8) :: w(mnMax_si), t, divdr, divdi
    INTEGER :: mn, k, m1
    IF (myid.eq.0.and.mnFirst.eq.1) THEN
       w(1) = 0.0_r8
       m1 = 2
      ELSE
       m1 = mnFirst
    ENDIF
    IF (input .lt. 0) THEN
       DO  mn=m1,mnLast
          w(mn)=1.0_r8/SQRT(REAL(snnp1_si(2*mn-1),r8))
       END DO
       DO k=1,mods
          t=1.0_r8/(er*SQRT(gh(k)))
          DO mn=mnFirst,mnLast
             vord(1,mn,k)=w(mn)*vord(1,mn,k)
             vord(2,mn,k)=w(mn)*vord(2,mn,k)
             divdr=-w(mn)*divd(2,mn,k)
             divdi= w(mn)*divd(1,mn,k)
             divd(1,mn,k)=divdr
             divd(2,mn,k)=divdi
             comd(1,mn,k)=t*comd(1,mn,k)
             comd(2,mn,k)=t*comd(2,mn,k)
          END DO
       END DO
    ELSE IF (input .gt. 0) THEN
       DO mn=m1,mnLast
          w(mn)=SQRT(REAL(snnp1_si(2*mn-1),r8))
       END DO
       DO k=1,mods
          t=er*SQRT(gh(k))
          DO mn=mnFirst,mnLast
             vord(1,mn,k)=w(mn)*vord(1,mn,k)
             vord(2,mn,k)=w(mn)*vord(2,mn,k)
             divdr= w(mn)*divd(2,mn,k)
             divdi=-w(mn)*divd(1,mn,k)
             divd(1,mn,k)=divdr
             divd(2,mn,k)=divdi
             comd(1,mn,k)=t*comd(1,mn,k)
             comd(2,mn,k)=t*comd(2,mn,k)
          END DO
       END DO
    END IF
  END SUBROUTINE primes



  ! horiz1 : arranges the vertically projected, scaled coefficients
  !          of vorticity, divergence, and composite mass variable
  !          for calculation of the adjustment to these fields by
  !          routine "horiz2". next, rearranges the calculated
  !          adjustment of vorticity, divergence, and composite
  !          variable back into the same form (and overwrites) the
  !          original input fields.



  SUBROUTINE horiz1(vord,divd,comd,v,per,g,indper,indg,jg,nas,myms,nms)
    REAL(KIND=r8), INTENT(INOUT) :: vord(2,mymnMax,modsloc)
    REAL(KIND=r8), INTENT(INOUT) :: divd(2,mymnMax,modsloc)
    REAL(KIND=r8), INTENT(INOUT) :: comd(2,mymnMax,modsloc)
    REAL(KIND=r8), INTENT(IN)    :: v
    REAL(KIND=r8), INTENT(IN)    :: per(*)
    REAL(KIND=r8), INTENT(IN)    :: g(*)
    INTEGER, INTENT(IN) :: nms
    INTEGER, INTENT(IN) :: myms(nms)
    INTEGER, INTENT(IN) :: indper(2*mymMax,modsloc)
    INTEGER(KIND=i8), INTENT(IN) :: jg(2*mymMax,modsloc)
    INTEGER, INTENT(IN) :: indg(2*mymMax,modsloc)
    INTEGER(KIND=i8), INTENT(IN) :: nas(2*mymMax,modsloc)
    !
    REAL(KIND=r8) :: sdot(2,3*(mMax+1)/2), adot(2,3*(mMax+1)/2)
    INTEGER :: nn, k, modes, l, ll, jsod, jsev, jevpod, i, mi
    INTEGER :: nends, nenda, lx, ir, nnmax, mglob
    !
    DO modes=1,modsloc
       DO mi=1,nms
          ll = myms(mi)
          i = 2*ll-1
          k=0
          mglob=lm2m(ll)
          l=mglob-1
          nnmax=mmax+1-mglob
          jsod=nnmax/2
          jsev=nnmax-jsod
          jevpod=jsev+jsod
          nends=jevpod+jsev
          nenda=jevpod+jsod
          DO nn=mglob+1,mmax,2
             lx=mymnMap(ll,nn)
             k=k+1
             DO ir=1,2
                sdot(ir,k       )=vord(ir,lx,modes)
                adot(ir,k+jsev  )=divd(ir,lx,modes)
                adot(ir,k+jevpod)=comd(ir,lx,modes)
             END DO
          END DO
          k=0
          DO nn=mglob,mmax,2
             lx=mymnMap(ll,nn)
             k=k+1
             DO ir=1,2
                adot(ir,k       )=vord(ir,lx,modes)
                sdot(ir,k+jsod  )=divd(ir,lx,modes)
                sdot(ir,k+jevpod)=comd(ir,lx,modes)
             END DO
          END DO
          IF (jg(i,modes).ne.0) &
             CALL horiz2(sdot,nas(i,modes),per(indper(i,modes)),&
                         g(indg(i,modes)),jg(i,modes),v,l)
          i = i+1
          IF (jg(i,modes).ne.0) &
             CALL horiz2(adot,nas(i,modes),per(indper(i,modes)),&
                         g(indg(i,modes)),jg(i,modes),v,l)
          k=0
          DO nn=mglob+1,mmax,2
             lx=mymnMap(ll,nn)
             k=k+1
             DO ir=1,2
                vord(ir,lx,modes)=sdot(ir,k)
                divd(ir,lx,modes)=adot(ir,k+jsev)
                comd(ir,lx,modes)=adot(ir,k+jevpod)
             END DO
          END DO
          k=0
          DO nn=mglob,mmax,2
             lx=mymnMap(ll,nn)
             k=k+1
             DO ir=1,2
                vord(ir,lx,modes)=adot(ir,k)
                divd(ir,lx,modes)=sdot(ir,k+jsod)
                comd(ir,lx,modes)=sdot(ir,k+jevpod)
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE horiz1

  SUBROUTINE getperg(nas,per,g,jg)
    INTEGER(KIND=i8), INTENT(IN) :: jg, nas
    REAL(KIND=r8), INTENT(INOUT) :: per(jg), g(jg,nas)
    !
    !     per stores periods of gravity modes,g stores eigenvectors
    !     both are read in this routine
    !

    !
    READ(UNIT=nfnmi) per,g

  END SUBROUTINE getperg

  ! horiz2 : calculates for one zonal wave number the Machenhauer
  !          adjustment to the vertically projected, properly scaled
  !          and rearranged spectral coefficients of vorticity,
  !          divergence, and composite mass variable.



  SUBROUTINE horiz2(dot,nas,per,g,jg,percut,l)
    INTEGER, INTENT(IN) :: l
    INTEGER(KIND=i8), INTENT(IN) :: jg, nas
    REAL(KIND=r8), INTENT(INOUT) :: dot(2,nas)
    REAL(KIND=r8), INTENT(IN) :: per(jg), g(jg,nas)
    REAL(KIND=r8), INTENT(IN) :: percut
    REAL(KIND=r8) :: period, dif, difcut
    REAL(KIND=r8) :: y(2,3*(mMax+1)/2), yi, yr
    INTEGER :: knit(3*(mMax+1)/2)
    INTEGER :: ndho, i, j, k, n
    !
    !     nas=vector size of sym. or  asy. tendencies stored in dot
    !     jg=nunber of gravity modes. jg=jcap for l=0 sym and asy cases
    !     for l.ne.0 jg=jcap2 for sym,jg=jcap for asy.
    !     per stores periods of gravity modes,g stores eigenvectors
    !
    dif=dk
    ndho=ndord/2
    DO j=1,jg
       knit(j)=0
       y(1,j)=0.0_r8
       y(2,j)=0.0_r8
    END DO
    n=l
    DO i=1,jg
       n=n+1
       period=ABS(per(i))
       IF (period-percut.le.0.0_r8) THEN
          !
          !     arbitray even order horizontal diffusion now used.
          !     smaller "tweek" constant used to avoid inaccuracies when
          !     reciprocal is taken.
          !     ________________________________________________________
          difcut=dif*(REAL(n*(n+1),r8)**ndho) + 1.0e-7_r8
          difcut=0.5_r8/difcut
          IF (difcut-period.ge.0.0_r8) THEN
             knit(i)=i
             DO j=1,nas
                y(1,i)=y(1,i)+g(i,j)*dot(1,j)
                y(2,i)=y(2,i)+g(i,j)*dot(2,j)
             END DO
             yr=y(1,i)
             yi=y(2,i)
             y(1,i)= per(i)*yi
             y(2,i)=-per(i)*yr
          END IF
       END IF
    END DO
    DO j=1,nas
       dot(1,j)=0.0_r8
       dot(2,j)=0.0_r8
    END DO
    DO j=1,nas
       DO k=1,jg
          IF (knit(k) .ne. 0) THEN
             dot(1,j)=dot(1,j)+g(k,j)*y(1,k)
             dot(2,j)=dot(2,j)+g(k,j)*y(2,k)
          END IF
       END DO
    END DO
  END SUBROUTINE horiz2



  !Getmod : computes vertical and horizontal
  !         modes to be used in the
  !         initialization of the model.
  !         The modes are stored in files
  !         which can be used by the model,
  !         without calling getmod again.
  !         The modes depend not only on the
  !         resolution of the model, but
  !         also on the Eulerian or
  !         Semi-Lagrangian option.



  SUBROUTINE Getmod(fName)
    CHARACTER(LEN=*), INTENT(IN) :: fName

    INTEGER, PARAMETER :: matz=1
    REAL(KIND=r8),    PARAMETER :: eps=2.0_r8**(-50)
    CHARACTER(LEN=*), PARAMETER :: h="**(GetMod)**"
    REAL(KIND=r8)    :: gh(kmax)
    INTEGER :: ierr
    
    OPEN(UNIT=nfnmi,FILE=TRIM(fName),FORM='unformatted',ACCESS='sequential',&
         ACTION='readwrite', STATUS='replace',IOSTAT=ierr)
    CALL Vermod(gh,eps,matz)

    ! how many modes

    CALL SetMods(gh)

    CALL Hormod(gh,eps)
    CLOSE(UNIT=nfnmi)
  END SUBROUTINE Getmod






  SUBROUTINE Vermod(gh,eps,matz)
    REAL(KIND=r8),    INTENT(OUT) :: gh(kmax)
    REAL(KIND=r8),    INTENT(IN ) :: eps
    INTEGER, INTENT(IN ) :: matz

    INTEGER :: j
    INTEGER :: k
    REAL(KIND=r8)    :: p
    REAL(KIND=r8)    :: siman
    REAL(KIND=r8)    :: soma
    REAL(KIND=r8)    :: er2
    REAL(KIND=r8)    :: eigg(kmax,kmax)
    REAL(KIND=r8)    :: eiggt(kmax,kmax)
    REAL(KIND=r8)    :: eigvc(kmax,kmax)
    REAL(KIND=r8)    :: dotpro(kmax)
    REAL(KIND=r8)    :: col(kmax)
    REAL(KIND=r8)    :: g(kmax,kmax)
    REAL(KIND=r8)    :: gt(kmax,kmax)

    CALL bmcm(1.0_r8)
    er2=er*er
    g = - cm * er2
    gt = TRANSPOSE (g)
    siman=-1.0_r8
    CALL Vereig(g,siman,eigvc,col,eigg,gh,eps,matz)
    siman=1.0_r8
    CALL Vereig(gt,siman,eigvc,col,eiggt,dotpro,eps,matz)

    ! dotpro=inverse dot prod. of eigenvec(g)*eigenvec(gtranspose)

    DO k=1,kmax
       soma=0.0_r8
       DO j=1,kmax
          soma=soma+eigg(j,k)*eiggt(j,k)
       END DO
       dotpro(k)=1.0_r8/soma
    END DO
    WRITE(UNIT=nfnmi)eigg,eiggt,gh,dotpro,tref
  END SUBROUTINE Vermod


  SUBROUTINE Hormod(gh,eps)
    REAL(KIND=r8),    INTENT(IN) :: eps
    REAL(KIND=r8),    INTENT(IN) :: gh(mods)

    INTEGER, PARAMETER :: ipr=-1
    INTEGER :: nxsy
    INTEGER :: nxas
    INTEGER :: k
    INTEGER :: m
    INTEGER :: nlx
    INTEGER :: nmd
    INTEGER :: lmax
    INTEGER :: klmx
    INTEGER :: mmmax
    INTEGER :: nnmax
    INTEGER(KIND=i8) :: nsy
    INTEGER(KIND=i8) :: nas
    INTEGER :: n
    INTEGER(KIND=i8) :: ncuts
    INTEGER(KIND=i8) :: ncuta
    INTEGER :: mend
    INTEGER :: modd
    INTEGER :: lend
    INTEGER :: kend
    REAL(KIND=r8)    :: alfa(mMax)
    REAL(KIND=r8)    :: beta(mMax)
    REAL(KIND=r8)    :: gama(mMax)
    REAL(KIND=r8)    :: dgl(mMax)
    REAL(KIND=r8)    :: sdg(mMax)
    REAL(KIND=r8)    :: xx(mMax,mMax)
    REAL(KIND=r8)    :: wk(3*(mMax+1)/2)
    REAL(KIND=r8)    :: ws(3*(mMax+1)/2)
    REAL(KIND=r8)    :: wa(3*(mMax+1)/2)
    REAL(KIND=r8)    :: xs(3*(mMax+1)/2*3*(mMax+1)/2)
    REAL(KIND=r8)    :: xa(3*(mMax+1)/2*3*(mMax+1)/2)
    REAL(KIND=r8)    :: es(3*(mMax+1)/2*3*(mMax+1)/2)
    REAL(KIND=r8)    :: ea(3*(mMax+1)/2*3*(mMax+1)/2)
    REAL(KIND=r8)    :: rm
    REAL(KIND=r8)    :: rn

    mend = mMax-1
    modd = MOD(mMax-1,2)
    lend = (mend+modd)/2
    kend = lend+1-modd
    nxsy = lend+2*kend
    nxas = kend+2*lend
    DO k=1,mods
       DO m=1,mMax
          nnmax=mMax-m+1
          nlx=mMax-m
          nmd=MOD(nlx,2)
          lmax=(nlx+nmd)/2
          mmmax=lmax+1-nmd
          klmx=lmax+mmmax
          nsy=lmax+2*mmmax
          nas=mmmax+2*lmax
          rm=REAL(m-1,r8)
          DO n=1,nnmax
             rn=rm+REAL(n-1,r8)
             IF (rn .EQ. 0.0_r8 ) THEN
                alfa(n)=0.0_r8
                beta(n)=0.0_r8
                gama(n)=0.0_r8
             ELSE
                alfa(n)=twomg*rm/(rn*(rn+1.0_r8))
                beta(n)=(twomg/rn)*SQRT((rn*rn-1.0_r8)*(rn*rn-rm*rm)/&
                     (4.0_r8*rn*rn-1.0_r8))
                gama(n)=eriv*SQRT(rn*(rn+1.0_r8)*gh(k))
             END IF
          END DO
          IF (m .EQ. 1) THEN

             ! symmetric case

             CALL symg0(nxsy,nsy,mMax,lmax,klmx,nmd,ipr,ncuts, &
                  eps,twomg,beta,gama, &
                  ws,sdg,dgl,xs,xx)
             CALL record(nxsy,ncuts,nsy,ws,xs,wk,es)

             ! asymmetric case

             CALL asyg0(nxas,nas,mMax,lmax,klmx,mmmax,nmd,ipr, &
                  ncuta,eps,twomg,beta,gama, &
                  wa,sdg,dgl,xa,xx)
             CALL record(nxas,ncuta,nas,wa,xa,wk,ea)
          ELSE

             ! symmetric case

             CALL symrg(nxsy,nsy,lmax,mmmax,nmd,ipr,ncuts, &
                  eps,twomg,percut,alfa,beta,gama, &
                  ws,wk,es,xs)
             CALL record(nxsy,ncuts,nsy,ws,xs,wk,es)

             ! asymmetric case

             CALL asyrg(nxas,nas,lmax,mmmax,nmd,ipr,ncuta, &
                  eps,twomg,percut,alfa,beta,gama, &
                  wa,wk,ea,xa)
             CALL record(nxas,ncuta,nas,wa,xa,wk,ea)
          END IF
       END DO
    END DO
  END SUBROUTINE Hormod






  SUBROUTINE asyg0(nxas,nas,nend1,lmax,klmx,mmax,nmd,ipr, &
       ncuta,eps,twomg,beta,gama,wa,sdg,dgl,xa,xx)
    INTEGER, INTENT(IN ) :: nxas
    INTEGER(KIND=i8), INTENT(IN ) :: nas
    INTEGER, INTENT(IN ) :: nend1
    INTEGER, INTENT(IN ) :: lmax
    INTEGER, INTENT(IN ) :: klmx
    INTEGER, INTENT(IN ) :: mmax
    INTEGER, INTENT(IN ) :: nmd
    INTEGER, INTENT(IN ) :: ipr
    INTEGER(KIND=i8), INTENT(OUT) :: ncuta
    REAL(KIND=r8),    INTENT(IN ) :: eps
    REAL(KIND=r8),    INTENT(IN ) :: twomg
    REAL(KIND=r8),    INTENT(IN ) :: beta(:)
    REAL(KIND=r8),    INTENT(IN ) :: gama(:)
    REAL(KIND=r8),    INTENT(OUT) :: wa(:)
    REAL(KIND=r8),    INTENT(OUT) :: sdg(:)
    REAL(KIND=r8),    INTENT(OUT) :: dgl(:)
    REAL(KIND=r8),    INTENT(OUT) :: xa(nxas,*)
    REAL(KIND=r8),    INTENT(OUT) :: xx(nend1,*)

    INTEGER :: n
    INTEGER :: nmx
    INTEGER :: nn
    INTEGER :: j
    INTEGER :: jj
    INTEGER :: ierr


    ! asymmetric case

    nmx=1
    n=1
    nn=2
    sdg(n)=beta(nn-1)*beta(nn)
    dgl(n)=beta(nn)*beta(nn)+beta(nn+1)*beta(nn+1)+ &
         gama(nn)*gama(nn)
    nmx=lmax-1
    DO  n=2,nmx
       nn=2*n
       sdg(n)=beta(nn-1)*beta(nn)
       dgl(n)=beta(nn)*beta(nn)+beta(nn+1)*beta(nn+1)+ &
            gama(nn)*gama(nn)
    END DO
    nmx=lmax
    n=lmax
    nn=2*n
    sdg(n)=beta(nn-1)*beta(nn)
    IF (nmd .EQ. 0) THEN
       dgl(n)=beta(nn)*beta(nn)+beta(nn+1)*beta(nn+1)+ &
            gama(nn)*gama(nn)
    ELSE
       dgl(n)=beta(nn)*beta(nn)+gama(nn)*gama(nn)
    END IF

    IF (ipr .GE. 1) THEN
       WRITE(UNIT=nfprt,FMT=*)' '
       WRITE(UNIT=nfprt,FMT=*)' sdg:'
       WRITE(UNIT=nfprt,FMT='(1p,6g12.5)')(sdg(n),n=1,nmx)
       WRITE(UNIT=nfprt,FMT=*)' dga:'
       WRITE(UNIT=nfprt,FMT='(1p,6g12.5)')(dgl(n),n=1,nmx)
    END IF

    CALL ident(nend1,nmx,xx)
    CALL tql2(nend1,nmx,dgl,sdg,xx,eps,ierr)

    DO j=1,nmx
       jj=2*j
       wa(jj-1)=-SQRT(dgl(j))
       wa(jj)=SQRT(dgl(j))
       xa(1,jj-1)=0.0_r8
       xa(1,jj)=-xa(1,jj-1)
       xa(mmax+1,jj-1)=xx(1,j)
       xa(mmax+1,jj)=xa(mmax+1,jj-1)
       xa(klmx+1,jj-1)=gama(2)*xx(1,j)/wa(jj-1)
       xa(klmx+1,jj)=-xa(klmx+1,jj-1)
       DO n=2,nmx
          nn=2*n
          xa(n,jj-1)=(beta(nn-1)*xx(n-1,j)+beta(nn)*xx(n,j))/wa(jj-1)
          xa(n,jj)=-xa(n,jj-1)
          xa(mmax+n,jj-1)=xx(n,j)
          xa(mmax+n,jj)=xa(mmax+n,jj-1)
          xa(klmx+n,jj-1)=gama(nn)*xx(n,j)/wa(jj-1)
          xa(klmx+n,jj)=-xa(klmx+n,jj-1)
       END DO
       n=nmx
       nn=2*n+1
       IF (nmd .EQ. 0) THEN
          xa(mmax,jj-1)=beta(nn)*xx(n,j)/wa(jj-1)
          xa(mmax,jj)=-xa(mmax,jj-1)
       END IF
    END DO

    ncuta=2*nmx
    CALL filter(nxas,nas,ncuta,xa,0.0_r8,eps)

    IF (ipr .GE. 1) THEN
       WRITE(UNIT=nfprt,FMT=*)' '
       WRITE(UNIT=nfprt,FMT=*)' frequency: nas=',nas,' ncuta=',ncuta
       WRITE(UNIT=nfprt,FMT='(1p,6g12.5)')(wa(n)/twomg,n=1,ncuta)
       WRITE(UNIT=nfprt,FMT=*)' period:'
       WRITE(UNIT=nfprt,FMT='(1p,6g12.5)')(1.0_r8/wa(n),n=1,ncuta)
       WRITE(UNIT=nfprt,FMT=*)' '
       WRITE(UNIT=nfprt,FMT=*)' xa: ierr=',ierr
       DO n=1,nmx
          WRITE(UNIT=nfprt,FMT='(1p,6g12.5)')(xx(n,nn),nn=1,nmx)
       END DO
       WRITE(UNIT=nfprt,FMT=*)' '
       WRITE(UNIT=nfprt,FMT=*)' xa:'
       DO n=1,nas
          WRITE(UNIT=nfprt,FMT='(1p,6g12.5)')(xa(n,nn),nn=1,MIN(6_i8,ncuta))
       END DO
    END IF
  END SUBROUTINE asyg0






  SUBROUTINE asyrg(nxas,nas,lmax,mmax,nmd,ipr,ncuta, &
       eps,twomg,percut,alfa,beta,gama,wa,wk,ea,xa)
    INTEGER, INTENT(IN) :: nxas,lmax,mmax,nmd,ipr
    INTEGER(KIND=i8), INTENT(IN) :: nas
    INTEGER(KIND=i8), INTENT(OUT) :: ncuta
    REAL(KIND=r8), INTENT(IN) :: eps,twomg,percut
    REAL(KIND=r8), INTENT(IN) ::alfa(*),beta(*),gama(*)
    REAL(KIND=r8), INTENT(OUT) :: wa(*),wk(*)
    REAL(KIND=r8), INTENT(OUT) :: ea(nxas,*),xa(nxas,*)
    !
    INTEGER :: n,nn,mm,jj,ierr
    !
    !   asymmetric case
    !
    DO nn=1,nxas
       DO mm=1,nxas
          ea(mm,nn)=0.0_r8
       END DO
       wa(nn)=0.0_r8
       wk(nn)=0.0_r8
    END DO
    !
    DO n=1,mmax
       ea(n,n)=alfa(2*n-1)
    END DO
    DO n=1,lmax
       nn=2*n
       jj=mmax+n
       ea(n,jj)=beta(nn)
       ea(jj,n)=ea(n,jj)
       IF (n.LT.lmax .OR. nmd.NE.1) THEN
          ea(n+1,jj)=beta(nn+1)
          ea(jj,n+1)=ea(n+1,jj)
       END IF
    END DO
    DO n=1,lmax
       nn=2*n
       jj=mmax+n
       ea(jj,jj)=alfa(nn)
       mm=jj+lmax
       ea(jj,mm)=gama(nn)
       ea(mm,jj)=ea(jj,mm)
    END DO
    !
    IF (ipr .GE. 3) THEN
       WRITE(UNIT=nfprt,FMT=*)' '
       WRITE(UNIT=nfprt,FMT=*)' ea:'
       DO n=1,nas
          WRITE(UNIT=nfprt,FMT='(1p,6g12.5)')(ea(n,nn),nn=1,nas)
       END DO
    END IF
    !
    CALL tred2(nxas,nas,ea,wa,wk,xa)
    CALL tql2(nxas,nas,wa,wk,xa,eps,ierr)
    !
    IF (ipr .GE. 1) THEN
       WRITE(UNIT=nfprt,FMT=*)' '
       WRITE(UNIT=nfprt,FMT=*)' wa: ierr=',ierr
       WRITE(UNIT=nfprt,FMT='(1p,6g12.5)')(wa(n),n=1,nas)
    END IF
    !
    !   reordering frequencies
    !
    CALL order(nxas,nas,wa,wk,xa,ea,percut,ncuta)
    CALL filter(nxas,nas,ncuta,xa,0.0_r8,eps)
    !
    IF (ipr .GE. 1) THEN
       WRITE(UNIT=nfprt,FMT=*)' '
       WRITE(UNIT=nfprt,FMT=*)' frequency: nas=',nas,' ncuta=',ncuta
       WRITE(UNIT=nfprt,FMT='(1p,6g12.5)')(wa(n)/twomg,n=1,ncuta)
       WRITE(UNIT=nfprt,FMT=*)' period:'
       WRITE(UNIT=nfprt,FMT='(1p,6g12.5)')(1.0_r8/wa(n),n=1,ncuta)
       WRITE(UNIT=nfprt,FMT=*)' '
       WRITE(UNIT=nfprt,FMT=*)' xa:'
       DO n=1,nas
          WRITE(UNIT=nfprt,FMT='(1p,6g12.5)')(xa(n,nn),nn=1,MIN(6_i8,ncuta))
       END DO
    END IF
  END SUBROUTINE asyrg






  SUBROUTINE filter(nm,ni,nj,a,b,c)
    !
    INTEGER, INTENT(IN) ::  nm
    INTEGER(KIND=i8), INTENT(IN) ::  ni,nj
    REAL(KIND=r8), INTENT(IN) ::  b,c
    REAL(KIND=r8), INTENT(INOUT) :: a(nm,*)
    INTEGER :: n,j
    !
    DO j=1,nj
       DO n=1,ni
          IF (ABS(a(n,j)) .LE. c) a(n,j)=b
       END DO
       !
    END DO
  END SUBROUTINE filter






  SUBROUTINE ident(nm,n,z)
    !
    !
    !   ** initialize z to identity matrix by sqrt(2)
    !
    INTEGER, INTENT(IN) ::  nm,n
    REAL(KIND=r8), INTENT(OUT) ::  z(nm,*)
    !
    INTEGER ::  i,j
    REAL(KIND=r8) :: sqrt2
    !
    sqrt2=1.0_r8/SQRT(2.0_r8)
    DO i=1,n
       DO j=1,n
          z(i,j)=0.0_r8
       END DO
       z(i,i)=sqrt2
    END DO
  END SUBROUTINE ident






  SUBROUTINE order(nm,n,fr,fw,z,zw,percut,nf)
    !
    INTEGER, INTENT(IN) :: nm
    INTEGER(KIND=i8), INTENT(IN) :: n
    REAL(KIND=r8), INTENT(IN) :: percut
    INTEGER(KIND=i8), INTENT(OUT) :: nf
    REAL(KIND=r8), INTENT(INOUT) :: fr(nm),fw(nm),z(nm,*),zw(nm,*)
    !
    INTEGER :: nm1,k,j,j1,i,jc,nc
    REAL(KIND=r8) :: chg
    !
    nm1=n-1
10  k=0
    DO j=1,nm1
       j1=j+1
       IF (ABS(fr(j)) .GT. ABS(fr(j1))) THEN
          chg=fr(j)
          DO i=1,n
             fw(i)=z(i,j)
          END DO
          fr(j)=fr(j1)
          DO i=1,n
             z(i,j)=z(i,j1)
          END DO
          fr(j1)=chg
          DO i=1,n
             z(i,j1)=fw(i)
          END DO
          k=1
       END IF
    END DO
    IF (k .NE. 0) GOTO 10
    !
    IF (percut .LE. 0.0_r8) THEN
       nf=n
       RETURN
    END IF
    !
    nc=0
    DO j=1,n
       IF (ABS(1.0_r8/fr(j)) .GT. percut) nc=j
    END DO
    nf=n-nc
    nc=nc+1
    !
    DO j=1,n
       fw(j)=fr(j)
       DO i=1,n
          zw(i,j)=z(i,j)
       END DO
    END DO
    !
    DO i=1,nm
       DO j=1,n
          z(i,j)=0.0_r8
       END DO
       fr(i)=0.0_r8
    END DO
    !
    DO jc=nc,n
       j=jc+1-nc
       fr(j)=fw(jc)
       DO i=1,n
          z(i,j)=zw(i,jc)
       END DO
    END DO
    !
    nc=nc-1
    DO jc=1,nc
       j=n+jc-nc
       fr(j)=fw(jc)
       DO i=1,n
          z(i,j)=zw(i,jc)
       END DO
    END DO
  END SUBROUTINE order






  SUBROUTINE record(nx,nc,nm,ww,xx,pp,gg)
    !
    INTEGER, INTENT(IN) :: nx
    INTEGER(KIND=i8), INTENT(IN) :: nc,nm
    REAL(KIND=r8), INTENT(IN) :: ww(nx),xx(nx,nx)
    REAL(KIND=r8), INTENT(OUT) :: pp(nc),gg(nc,nm)
    !
    INTEGER :: n,nn
    !
    WRITE(UNIT=nfnmi)nc,nm
    IF (nc .EQ. 0) RETURN
    !
    pp = 0.0_r8
    gg = 0.0_r8
    !
    DO n=1,nc
       pp(n)=1.0_r8/ww(n)
       DO nn=1,nm
          gg(n,nn)=xx(nn,n)
       END DO
    END DO
    WRITE(UNIT=nfnmi)pp,gg
  END SUBROUTINE record







  SUBROUTINE symg0(nxsy,nsy,nend1,lmax,klmx,nmd,ipr,ncuts, &
       eps,twomg,beta,gama, &
       ws,sdg,dgl,xs,xx)
    INTEGER, INTENT(IN) ::  nxsy,nend1,lmax,klmx,nmd,ipr
    INTEGER(KIND=i8), INTENT(IN) ::  nsy
    INTEGER(KIND=i8), INTENT(OUT) ::  ncuts
    REAL(KIND=r8), INTENT(IN) :: eps,twomg
    REAL(KIND=r8), INTENT(IN) :: beta(*),gama(*)
    REAL(KIND=r8), INTENT(OUT) :: ws(*),sdg(*),dgl(*)
    REAL(KIND=r8), INTENT(OUT) :: xs(nxsy,*),xx(nend1,*)
    !
    INTEGER :: n,nmx,nn,j,jj,ierr
    !
    !     symmetric case
    !
    nmx=1
    n=1
    sdg(n)=0.0_r8
    nn=2*n+1
    dgl(n)=beta(nn)*beta(nn)+beta(nn+1)*beta(nn+1)+ &
         gama(nn)*gama(nn)
    nmx=lmax-1
    DO n=2,nmx
       nn=2*n+1
       sdg(n)=beta(nn-1)*beta(nn)
       dgl(n)=beta(nn)*beta(nn)+beta(nn+1)*beta(nn+1)+ &
            gama(nn)*gama(nn)
    END DO
    IF (nmd .EQ. 0) THEN
       nmx=lmax
       n=lmax
       nn=2*n+1
       sdg(n)=beta(nn-1)*beta(nn)
       dgl(n)=beta(nn)*beta(nn)+gama(nn)*gama(nn)
    END IF
    !
    IF (ipr .GE. 1) THEN
       WRITE(UNIT=nfprt,FMT=*)' '
       WRITE(UNIT=nfprt,FMT=*)' sdg:'
       WRITE(UNIT=nfprt,FMT='(1p,6g12.5)')(sdg(n),n=1,nmx)
       WRITE(UNIT=nfprt,FMT=*)' dgs:'
       WRITE(UNIT=nfprt,FMT='(1p,6g12.5)')(dgl(n),n=1,nmx)
    END IF
    !
    CALL ident(nend1,nmx,xx)
    CALL tql2(nend1,nmx,dgl,sdg,xx,eps,ierr)
    !
    DO j=1,nmx
       jj=2*j
       ws(jj-1)=-SQRT(dgl(j))
       ws(jj)=SQRT(dgl(j))
       xs(1,jj-1)=beta(3)*xx(1,j)/ws(jj-1)
       xs(1,jj)=-xs(1,jj-1)
       xs(lmax+1,jj-1)=0.0_r8
       xs(lmax+1,jj)=xs(lmax+1,jj-1)
       xs(klmx+1,jj-1)=0.0_r8
       xs(klmx+1,jj)=-xs(klmx+1,jj-1)
       DO n=2,nmx
          nn=2*n+1
          xs(n,jj-1)=(beta(nn-1)*xx(n-1,j)+beta(nn)*xx(n,j))/ws(jj-1)
          xs(n,jj)=-xs(n,jj-1)
          xs(lmax+n,jj-1)=xx(n-1,j)
          xs(lmax+n,jj)=xs(lmax+n,jj-1)
          xs(klmx+n,jj-1)=gama(nn-2)*xx(n-1,j)/ws(jj-1)
          xs(klmx+n,jj)=-xs(klmx+n,jj-1)
       END DO
       n=lmax
       nn=2*n-1
       IF (nmd .EQ. 1) THEN
          xs(n,jj-1)=beta(nn+1)*xx(n-1,j)/ws(jj-1)
          xs(n,jj)=-xs(n,jj-1)
          xs(klmx,jj-1)=xx(n-1,j)
          xs(klmx,jj)=xs(klmx,jj-1)
          xs(nsy,jj-1)=gama(nn)*xx(n-1,j)/ws(jj-1)
          xs(nsy,jj)=-xs(nsy,jj-1)
       ELSE
          xs(klmx,jj-1)=xx(n,j)
          xs(klmx,jj)=xs(klmx,jj-1)
          xs(nsy,jj-1)=gama(nn+2)*xx(n,j)/ws(jj-1)
          xs(nsy,jj)=-xs(nsy,jj-1)
       END IF
    END DO
    !
    ncuts=2*nmx
    CALL filter(nxsy,nsy,ncuts,xs,0.0_r8,eps)
    !
    IF (ipr .GE. 1) THEN
       WRITE(UNIT=nfprt,FMT=*)' '
       WRITE(UNIT=nfprt,FMT=*)' frequency: nsy=',nsy,' ncuts=',ncuts
       WRITE(UNIT=nfprt,FMT='(1p,6g12.5)')(ws(n)/twomg,n=1,ncuts)
       WRITE(UNIT=nfprt,FMT=*)' period:'
       WRITE(UNIT=nfprt,FMT='(1p,6g12.5)')(1.0_r8/ws(n),n=1,ncuts)
       WRITE(UNIT=nfprt,FMT=*)' '
       WRITE(UNIT=nfprt,FMT=*)' xs: ierr=',ierr
       DO nn=1,nmx
          WRITE(UNIT=nfprt,FMT='(1p,6g12.5)')(xx(nn,n),n=1,nmx)
       END DO
       WRITE(UNIT=nfprt,FMT=*)' '
       WRITE(UNIT=nfprt,FMT=*)' xs:'
       DO n=1,nsy
          WRITE(UNIT=nfprt,FMT='(1p,6g12.5)')(xs(n,nn),nn=1,MIN(6_i8,ncuts))
       END DO
    END IF
  END SUBROUTINE symg0






  SUBROUTINE symrg(nxsy,nsy,lmax,mmax,nmd,ipr,ncuts, &
       eps,twomg,percut,alfa,beta,gama, &
       ws,wk,es,xs)
    !
    INTEGER, INTENT(IN) ::  nxsy,lmax,mmax,nmd,ipr
    INTEGER(KIND=i8), INTENT(IN) ::  nsy
    INTEGER(KIND=i8), INTENT(OUT) ::  ncuts
    REAL(KIND=r8), INTENT(IN) :: eps,twomg,percut
    REAL(KIND=r8), INTENT(IN) :: alfa(*),beta(*),gama(*)
    REAL(KIND=r8), INTENT(OUT) :: ws(*),wk(*)
    REAL(KIND=r8), INTENT(OUT) :: es(nxsy,*),xs(nxsy,*)
    !
    INTEGER ::  n,nn,mm,jj,ierr
    !
    !   symmetric case
    !
    DO nn=1,nxsy
       DO mm=1,nxsy
          es(mm,nn)=0.0_r8
       END DO
       ws(nn)=0.0_r8
       wk(nn)=0.0_r8
    END DO
    !
    DO n=1,lmax
       es(n,n)=alfa(2*n)
    END DO
    DO n=1,lmax
       nn=2*n
       jj=lmax+n
       es(n,jj)=beta(nn)
       es(jj,n)=es(n,jj)
       IF (n.LT.lmax .OR. nmd.NE.1) THEN
          es(n,jj+1)=beta(nn+1)
          es(jj+1,n)=es(n,jj+1)
       END IF
    END DO
    DO n=1,mmax
       nn=2*n-1
       jj=lmax+n
       es(jj,jj)=alfa(nn)
       mm=jj+mmax
       es(jj,mm)=gama(nn)
       es(mm,jj)=es(jj,mm)
    END DO
    !
    IF (ipr .GE. 3) THEN
       WRITE(UNIT=nfprt,FMT=*)' '
       WRITE(UNIT=nfprt,FMT=*)' es:'
       DO n=1,nsy
          WRITE(UNIT=nfprt,FMT='(1p,6g12.5)')(es(n,nn),nn=1,nsy)
       END DO
    END IF
    !
    CALL tred2(nxsy,nsy,es,ws,wk,xs)
    CALL tql2(nxsy,nsy,ws,wk,xs,eps,ierr)
    !
    IF (ipr .GE. 1) THEN
       WRITE(UNIT=nfprt,FMT=*)' '
       WRITE(UNIT=nfprt,FMT=*)' ws: ierr=',ierr
       WRITE(UNIT=nfprt,FMT='(1p,6g12.5)')(ws(n),n=1,nsy)
    END IF
    !
    !   reordering frequencies
    !
    CALL order(nxsy,nsy,ws,wk,xs,es,percut,ncuts)
    CALL filter(nxsy,nsy,ncuts,xs,0.0_r8,eps)
    !
    IF (ipr .GE. 1) THEN
       WRITE(UNIT=nfprt,FMT=*)' '
       WRITE(UNIT=nfprt,FMT=*)' frequency: nsy=',nsy,' ncuts=',ncuts
       WRITE(UNIT=nfprt,FMT='(1p,6g12.5)')(ws(n)/twomg,n=1,ncuts)
       WRITE(UNIT=nfprt,FMT=*)' period:'
       WRITE(UNIT=nfprt,FMT='(1p,6g12.5)')(1.0_r8/ws(n),n=1,ncuts)
       WRITE(UNIT=nfprt,FMT=*)' '
       WRITE(UNIT=nfprt,FMT=*)' xs:'
       DO n=1,nsy
          WRITE(UNIT=nfprt,FMT='(1p,6g12.5)')(xs(n,nn),nn=1,MIN(6_i8,ncuts))
       END DO
    END IF
  END SUBROUTINE symrg






  SUBROUTINE vereig(gg,siman,eigvc,col,vec,val,eps,matz)
    !
    INTEGER, INTENT(IN) :: matz
    REAL(KIND=r8), INTENT(OUT) ::  vec(kmax,kmax),val(kmax)
    REAL(KIND=r8), INTENT(OUT) ::  col(kmax)
    REAL(KIND=r8), INTENT(INOUT) ::  gg(kmax,kmax),eigvc(kmax,kmax)
    REAL(KIND=r8), INTENT(IN) ::  eps,siman
    !
    INTEGER :: kk(kmax)
    REAL(KIND=r8) :: eigvr(kmax),eigvi(kmax),wk1(kmax),wk2(kmax)
    !
    INTEGER :: k,ier,i,j,kkk
    REAL(KIND=r8) soma,rmax,e20
    !
    e20=-1.0e20_r8
    !
    CALL rg(kmax,kmax,gg,eigvr,eigvi,matz,eigvc,ier,eps,wk1,wk2)
    open(unit=98,file='eigenvalues_vert')
    write(98,*) 'vertical eigenvalues '
    do k=1,kmax
    write(98,*) eigvr(k),eigvi(k)
    enddo
    write(98,*) 'vertical eigenvectors'
    do k=1,kmax
    write(98,*) (eigvc(i,k),i=1,kmax)
    enddo
    !
    DO k=1,kmax
       kk(k)=0
       col(k)=siman*eigvr(k)
       soma=0.0_r8
       DO j=1,kmax
          soma=soma+eigvc(j,k)*eigvc(j,k)
       END DO
       !* soma=length of eigenvector k
       soma=1.0_r8/SQRT(soma)
       DO j=1,kmax
          eigvc(j,k)=soma*eigvc(j,k)
       END DO
    END DO
    !
    !   eigenvalues now have unit length.k th vector is eigvc(j,k)
    !   eigenvalues are now in col(k)
    !   next arrange in descending order
    !
    DO j=1,kmax
       rmax=e20
       kkk=0
       DO k=1,kmax
          IF(ABS(col(k)).GT.rmax) THEN
             kkk=k
             rmax=ABS(col(k))
          END IF
       END DO
       val(j)=col(kkk)
       col(kkk)=0.0_r8
       DO i=1,kmax
          vec(i,j)=eigvc(i,kkk)
       END DO
       kk(j)=kkk
    END DO
  END SUBROUTINE vereig







  SUBROUTINE w3fa03(press,temp)
    !
    REAL(KIND=r8), INTENT(IN) :: press
    REAL(KIND=r8), INTENT(OUT) ::  temp
    !
    !    this subroutine (w3fa03) computes the standard
    !    temperature (deg k) given the
    !    pressure in mb.
    !     u. s. standard atmosphere, 1962
    !    icao std atm to 20km
    !    proposed extension to 32km
    !    not valid for  pressure.lt.8.68mb
    !
    REAL(KIND=r8) grav,piso,ziso,salp,pzero,t0,alp, &
         ptrop,tstr,htrop
    REAL(KIND=r8) rovg,fkt,ar,pp0,height
    !
    DATA grav /9.80665e0_r8/, &
         piso /54.7487e0_r8/, ziso /20000.0e0_r8/, salp /-0.0010e0_r8/, &
         pzero /1013.25e0_r8/, t0 /288.15e0_r8/, alp /0.0065e0_r8/, &
         ptrop/226.321e0_r8/, tstr/216.65e0_r8/, htrop /11000.0e0_r8/
    !
    rovg=gasr/grav
    fkt=rovg*tstr
    ar=alp*rovg
    pp0=pzero**ar
    !
    IF(press.GE.piso.AND.press.LE.ptrop) THEN
       !
       !     compute isothermal cases
       !
       temp = tstr
    ELSE
       IF(press.LT.piso) THEN
          !
          !     compute lapse rate=-.0010_r8 cases
          !
          ar=salp*rovg
          pp0=piso**ar
          height=(tstr/(pp0*salp))*(pp0-(press**ar))+ziso
          temp=tstr-(height-ziso)*salp
       ELSE
          height=(t0/(pp0*alp))*(pp0-(press**ar))
          temp=t0-height*alp
       END IF
    END IF
  END SUBROUTINE w3fa03

  ! how many normal modes

  SUBROUTINE SetMods(gh)
    REAL(KIND=r8), INTENT(IN) :: gh(:)
    REAL(KIND=r8), PARAMETER :: HnCut=1000.0_r8
    CHARACTER(LEN=*), PARAMETER :: h="**(SetMods)**"
    INTEGER :: k, modperg, rest, m
    mods=0
    DO k=1,kmax
       IF (gh(k)/grav > HnCut) THEN
          mods=mods+ 1
!!$          WRITE(UNIT=nfprt,FMT='(A,I3,2(A,F8.2))') ' n = ', k, ' HnCut = ', HnCut, ' Hn = ', gh(k)/grav
       END IF
    END DO
    IF (mods == 0) THEN
       WRITE(UNIT=nfprt,FMT='(A)') ' ERROR: The Equivalent Heights of Normal Modes is Wrong'
       STOP 'SetMods  ==> (mods == 0)'
    END IF
    !
    !  split mods
    !  ----------
    IF (.NOT.ALLOCATED(nmodperg)) THEN
       ALLOCATE (nmodperg(ngroups_four))
       ALLOCATE (grouphasmod(mods))
       modperg = mods / ngroups_four
       rest = mods - modperg * ngroups_four
       m = 0
       DO k=1,ngroups_four
          IF (k.eq.mygroup_four) myfirstmod = m + 1
          IF (k.le.rest) THEN
             nmodperg(k) = modperg + 1
           ELSE
             nmodperg(k) = modperg
          END IF
          grouphasmod(m+1:m+nmodperg(k)) = k
          m = m+nmodperg(k)
       END DO
       modsloc = nmodperg(mygroup_four)
     ENDIF
  END SUBROUTINE SetMods
END MODULE NonLinearNMI
