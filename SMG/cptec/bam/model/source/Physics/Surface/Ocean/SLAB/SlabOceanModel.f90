MODULE SlabOceanModel
  USE Options, ONLY: &
       reducedGrid,nfprt
  USE Utils, ONLY: &
       IJtoIBJB, &
       LinearIJtoIBJB, &
       NearestIJtoIBJB,&
       FreqBoxIJtoIBJB

  USE FieldsPhysics, ONLY: &
   laymld,       hbath,     tdeep,sdeep
   
 IMPLICIT NONE
SAVE
  
  PRIVATE 
  ! Selecting Kinds
  INTEGER, PARAMETER :: r4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
  INTEGER, PARAMETER :: r16 = SELECTED_REAL_KIND(15)! Kind for 128-bits Real Numbers


  !  This is an example to use the ocean albedo look-up-table (ocnalbtab.bin) 
  !   to obtain the albedo at your specified spectral band, optical depth,
  !   cosine of solar zenith, wind and chlorophyll concentration.
  !  
  INTEGER, PARAMETER :: nb=24
  INTEGER, PARAMETER :: nt=16
  INTEGER, PARAMETER :: ns=15
  INTEGER, PARAMETER :: nw=7
  INTEGER, PARAMETER :: nc=5
  INTEGER, PARAMETER :: nline=8640
  REAL(KIND=r8), ALLOCATABLE  :: alb (:,:)
  REAL(KIND=r8), ALLOCATABLE  :: rflx(:,:)

  REAL(KIND=r8), PARAMETER  :: taunode(16)=(/0.0_r8, 0.05_r8, 0.1_r8, 0.16_r8, 0.24_r8,  0.35_r8, 0.5_r8, 0.7_r8, 0.99_r8, &
       1.3_r8, 1.80_r8, 2.5_r8, 5.00_r8, 9.00_r8, 15.0_r8, 25.0_r8 /)
  REAL(KIND=r8), PARAMETER  :: szanode(15)=(/0.05_r8, 0.09_r8, 0.15_r8, 0.21_r8, 0.27_r8, 0.33_r8, 0.39_r8, 0.45_r8, &
       0.52_r8, 0.60_r8, 0.68_r8, 0.76_r8, 0.84_r8, 0.92_r8, 1.0_r8 /)
  REAL(KIND=r8), PARAMETER  :: windnode(7)=(/0.0_r8, 3.0_r8, 6.0_r8, 9.0_r8, 12.0_r8, 15.0_r8, 18.0_r8 /)
  REAL(KIND=r8), PARAMETER  :: chlnode(5) =(/0.0_r8, 0.1_r8, 0.5_r8, 2.0_r8, 12.0_r8/)
  REAL(KIND=r8), PARAMETER  :: wlnode(25) =(/0.25_r8, 0.30_r8, 0.33_r8, 0.36_r8, 0.40_r8, 0.44_r8, 0.48_r8, 0.52_r8, &
       0.57_r8, 0.64_r8, 0.69_r8, 0.752_r8, 0.780_r8, 0.867_r8, 1.0_r8, 1.096_r8,&
       1.19_r8, 1.276_r8, 1.534_r8, 1.645_r8, 2.128_r8, 2.381_r8, 2.907_r8, &
       3.425_r8, 4.0_r8/)
!_________________________________________________________________________________________________
  REAL(KIND=r8)   , PARAMETER :: hslab=30.0_r8         !slab ocean depth
  INTEGER, PARAMETER :: km   =19                       !model vertical layers
  REAL(KIND=r8), PARAMETER :: capa =3950.0_r8         !heat capacity of sea water 
  REAL(KIND=r8), PARAMETER :: rhoref = 1024.438_r8    !sea water reference density, kg/m^3
  ! define model vertical structure, up to 1000 meter deep
  !layer thickness, defined within this subroutine  depth, from top down, m 
!  REAL(KIND=r8), PARAMETER :: dz(km)=(/&
!       10._r8,  5._r8,  5._r8,  5._r8,  5._r8,  5._r8,   5._r8,    5._r8,  5._r8, 10._r8,&
!       10._r8, 10._r8, 10._r8, 10._r8, 20._r8, 20._r8,  20._r8,   20._r8, 20._r8, 40._r8,&
!       40._r8, 40._r8, 60._r8, 60._r8, 60._r8,100._r8, 100._r8,  100._r8,100._r8,100._r8/)
  REAL(KIND=r8), PARAMETER :: dz(km)=(/&
        5._r8,  5._r8,  5._r8,  10._r8,  10._r8,  10._r8,   15._r8,    15._r8,   15._r8, 20._r8,&
       30._r8, 40._r8, 50._r8,  50._r8,  60._r8,  60._r8,   80._r8,    80._r8, 100._r8/)

! 5.0+  5.0+  5.0+  10.0+  10.0+  10.0+   15.0+    15.0+  15.0+ 20.0+30.0+ 40.0+ 50.0+ 50.0+ 60.0+ 60.0+  80.0+   80.0+ 1000.0
!  REAL(KIND=r8), PARAMETER :: dz(km)=(/&
!       10._r8,  10._r8,  10._r8,  20._r8,  25._r8,  25._r8,   25._r8,    25._r8,  50._r8,  50._r8,&
!       50._r8, 100._r8, 100._r8, 100._r8, 100._r8, 100._r8,  100._r8,   100._r8, 100._r8/)

  ! --- ocean depth and bottom layer of the ocean model
  INTEGER,ALLOCATABLE  :: lbottom(:,:)!   lbottom - indicates bottom of mixed layer (layer where the 
  REAL(KIND=r8) ::  zlev(km)           !sum(dz(k)), top down, m
  ! --- water optical type and coefficients for the transmission of light
!  INTEGER ,ALLOCATABLE          :: laymld(:,:)
!  REAL(KIND=r8) ,ALLOCATABLE    :: h(:,:)
!  REAL(KIND=r8) ,ALLOCATABLE    :: t(:,:,:)
!  REAL(KIND=r8) ,ALLOCATABLE    :: s(:,:,:)
  REAL(KIND=r8) ,ALLOCATABLE    :: tclim2(:,:,:,:)
  REAL(KIND=r8) ,ALLOCATABLE    :: sclim2(:,:,:,:)
  REAL(KIND=r8) ,ALLOCATABLE    ::r(:,:), d1(:,:), d2(:,:)
!_________________________________________________________________________________________________
 PUBLIC :: InitGetOceanAlb 
 PUBLIC :: GetOceanAlb
 PUBLIC :: driver
CONTAINS
  SUBROUTINE InitGetOceanAlb(imax,jMax,ibMax,jbMax,ibMaxPerJB,idatec,fNameSlabOcen,path_in,RESTART)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: iMax
    INTEGER, INTENT(IN) :: jMax
    INTEGER, INTENT(IN) :: ibMax
    INTEGER, INTENT(IN) :: jbMax
    INTEGER         , INTENT(IN   ) :: idatec(4)
    INTEGER         , INTENT(IN   ) :: ibMaxPerJB(jbMax)
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameSlabOcen
    CHARACTER(LEN=*), INTENT(IN   ) :: path_in
    LOGICAL         , INTENT(IN   ) ::  RESTART

    !    !!!Note: if the record unit in your sestem is WORDS instead of BYTE,
    !    !!!      you MUST change the "recl" (record length) in the following
    !    !!!      OPEN statement from 24*4 to 24 (recl=24).        
    INTEGER :: lrec
    INTEGER :: irec
    REAL(KIND=r4)  :: alb2 (nb)
    ALLOCATE(alb (nline,nb));alb=0.0_r8
    ALLOCATE(rflx(nline,nb));rflx=0.0_r8
    alb2=0.0_r4
    INQUIRE(IOLENGTH=lrec)alb2
    OPEN(1,file=TRIM(fNameSlabOcen),FORM='UNFORMATTED',ACCESS='DIRECT',recl=lrec)
    DO irec=1,nline
       READ(1,rec=irec)alb2
       rflx(irec,1:nb)=alb2(1:nb)
       alb (irec,1:nb)=alb2(1:nb)
    END DO
    CLOSE(1)

    CALL init_MLO  (iMax,jMAx,ibMax,jbMAx,ibMaxPerJB,idatec,path_in,RESTART)

    RETURN

  END SUBROUTINE InitGetOceanAlb

  SUBROUTINE init_MLO (iMax,jMAx,ibMax,jbMAx,ibMaxPerJB,idatec,path_in,RESTART)
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: ibMax
    INTEGER, INTENT(IN   ) :: jbMAx
    INTEGER, INTENT(IN   ) :: iMax
    INTEGER, INTENT(IN   ) :: jMAx
    INTEGER         , INTENT(IN   ) :: ibMaxPerJB(jbMax)
    CHARACTER(LEN=*), INTENT(IN   ) :: path_in
    INTEGER         , INTENT(IN   ) :: idatec(4)
    LOGICAL         , INTENT(IN   ) ::  RESTART

    REAL(KIND=r8)                   :: brf   (iMax,jMax)
    INTEGER                         :: ibrf  (iMax,jMax)
    REAL(KIND=r4)                   :: array (iMax,jMax)
    INTEGER (KIND=i8)               :: i8brf  (iMax,jMax)

    REAL(KIND=r8)                   :: bathy8(ibMax,jbMax)
    INTEGER                         :: ibr4f (ibMax,jbMax)
    INTEGER(KIND=i4)                :: br4f  (ibMax,jbMax)
    INTEGER(KIND=i8)                :: iopt  (ibMax,jbMax)
    REAL(KIND=r8)                   :: brf8  (ibMax,jbMax)
    INTEGER                         :: i,j,k,ierr,irec,im
    REAL(KIND=r8)                   :: htmp
    INTEGER :: LRecIN
    CHARACTER(LEN=5) :: c0

    ALLOCATE(lbottom(ibMax,jbMax))
    ALLOCATE( r     (ibMax,jbMax))
    ALLOCATE( d1    (ibMax,jbMax))
    ALLOCATE( d2    (ibMax,jbMax))
    ALLOCATE( tclim2(ibMax,km,12,jbMax))
    ALLOCATE( sclim2(ibMax,km,12,jbMax))
    WRITE(c0,"(i5.5)") jMax

    
    !h
    !t
    !s
    !laymld 
    !--------------------------------------------------------------------
    ! define model vertical structure, up to 1000 meter deep
    zlev(1)=dz(1)
    DO k=2,km
       zlev(k)=zlev(k-1)+dz(k)
    ENDDO
    array=0.0_r4
    ibr4f=0
    ibrf=0
    !
    !   read in ocean bathymetry (<0 over ocean; >=0 over land) 
    !   and determines the bottom of the model at each point, 
    !   lbottom indicates the first inactive layer
    INQUIRE (IOLENGTH=LRecIN) array
    OPEN(2, file=TRIM(path_in)//'/'//'ocean_depth'//'.G'//TRIM(c0), form='unformatted', &
         access='direct',recl=LRecIN, status='old', action='read', IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(path_in)//'/'//'ocean_depth'//'.G'//TRIM(c0), ierr
       STOP "**(ERROR)**"
    END IF
    READ(2,rec=1, IOSTAT=ierr)array
    IF (ierr /= 0) THEN
        WRITE(UNIT=nfprt,FMT="('**(ERROR)** Read file ',a,' returned iostat=',i4)") &
         TRIM(path_in)//'/'//'ocean_depth'//'.G'//TRIM(c0), ierr
        STOP "**(ERROR)**"
    END IF

    brf = REAL(array,KIND=r8)
    IF (reducedGrid) THEN
        CALL LinearIJtoIBJB(brf,bathy8)
    ELSE
        CALL IJtoIBJB(brf ,bathy8 )
    END IF
    CLOSE(2)

    DO j=1,jbMax 
       DO i=1, ibMaxPerJB(j)
          lbottom(i,j)=1
          DO k=1,km
             IF(-bathy8(i,j).GE.zlev(k)) lbottom(i,j)=k
          ENDDO
       ENDDO
    ENDDO


    !   global annual average optical water type from the map of 
    !   siminot and le treut (1986,jgr).
    !   water types  -   numerical value in file:
    !    land               0
    !    i                  1
    !    ii                 2
    !    iii                3
    !    ia                 4
    !    ib                 5
    INQUIRE (IOLENGTH=LRecIN) ibrf
    OPEN(2, file=TRIM(path_in)//'/'//'water_type'//'.G'//TRIM(c0), form='unformatted', &
         access='direct',recl=LRecIN, status='old', action='read', IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(path_in)//'/'//'water_type'//'.G'//TRIM(c0), ierr
       STOP "**(ERROR)**"
    END IF
    irec=0
    irec=irec+1
    READ(2,rec=1, IOSTAT=ierr)ibrf
    IF (ierr /= 0) THEN
        WRITE(UNIT=nfprt,FMT="('**(ERROR)** Read file ',a,' returned iostat=',i4)") &
         TRIM(path_in)//'/'//'water_type'//'.G'//TRIM(c0), ierr
        STOP "**(ERROR)**"
    END IF
    i8brf = INT(ibrf,KIND=i8)

    IF (reducedGrid) THEN
        CALL FreqBoxIJtoIBJB(i8brf,iopt)
    ELSE
        CALL IJtoIBJB(i8brf ,iopt )
    END IF

    CLOSE(2)
    
    !   the coefficients for the transmission of light are set by the
    !   water type (see paulson and simpson, jpo, 1977).
    DO j=1, jbMax
       DO i=1,  ibMaxPerJB(j)
          IF (iopt(i,j) .EQ. 1) THEN
             r(i,j)  = 0.58_r8
             d1(i,j) =  0.35_r8      
             d2(i,j) =  23.0_r8      
          ELSEIF (iopt(i,j) .EQ. 2) THEN
             r(i,j) = 0.77_r8
             d1(i,j) =  1.5_r8      
             d2(i,j) =  14.0_r8      
          ELSEIF (iopt(i,j) .EQ. 3) THEN
             r(i,j) = 0.78_r8
             d1(i,j) =  1.4_r8      
             d2(i,j) =  7.9_r8      
          ELSEIF (iopt(i,j) .EQ. 4) THEN
             r(i,j) = 0.62_r8
             d1(i,j) =  0.60_r8      
             d2(i,j) =  20.0_r8      
          ELSE
             !              for optical type ib => 5 then
             r(i,j) = 0.67_r8
             d1(i,j) =  1.0_r8      
             d2(i,j) =  17.0_r8      
          END IF
       END DO
    END DO

!   observed monthly ltm mean temp and salinity 
!   will be used in relaxation              
    INQUIRE (IOLENGTH=LRecIN) array
    OPEN(2, file=TRIM(path_in)//'/'//'temp_ltm_month'//'.G'//TRIM(c0), form='unformatted', &
         access='direct',recl=LRecIN, status='old', action='read', IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(path_in)//'/'//'temp_ltm_month'//'.G'//TRIM(c0), ierr
       STOP "**(ERROR)**"
    END IF
    irec=0
    DO im=1,12
       DO k=1,km
          irec=irec+1
          READ(2,rec=irec, IOSTAT=ierr)array
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Read file ',a,' returned iostat=',i4)") &
              TRIM(path_in)//'/'//'temp_ltm_month'//'.G'//TRIM(c0), ierr
             STOP "**(ERROR)**"
          END IF
          brf=REAL(array,KIND=r8)
          IF (reducedGrid) THEN
              CALL LinearIJtoIBJB(brf,brf8)
          ELSE
             CALL IJtoIBJB(brf ,brf8 )
          END IF
          DO j = 1, jbMax
             DO i = 1,  ibMaxPerJB(j) 
                   tclim2(i,k,im,j) = REAL(brf8(i,j),KIND=r8)
             END DO
          END DO
       END DO
    END DO
    CLOSE(2)

!   observed monthly ltm mean temp and salinity 
!   will be used in relaxation              
    INQUIRE (IOLENGTH=LRecIN) array
    OPEN(2, file=TRIM(path_in)//'/'//'salt_ltm_month'//'.G'//TRIM(c0), form='unformatted', &
         access='direct',recl=LRecIN, status='old', action='read', IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(path_in)//'/'//'salt_ltm_month'//'.G'//TRIM(c0), ierr
       STOP "**(ERROR)**"
    END IF
    irec=0
    DO im=1,12
       DO k=1,km
          irec=irec+1
          READ(2,rec=irec, IOSTAT=ierr)array
          IF (ierr /= 0) THEN
              WRITE(UNIT=nfprt,FMT="('**(ERROR)** Read file ',a,' returned iostat=',i4)") &
               TRIM(path_in)//'/'//'salt_ltm_month'//'.G'//TRIM(c0), ierr
              STOP "**(ERROR)**"
          END IF

          brf=REAL(array,KIND=r8)
          IF (reducedGrid) THEN
             CALL LinearIJtoIBJB(brf,brf8)
          ELSE
             CALL IJtoIBJB(brf ,brf8 )
          END IF
          DO j = 1, jbMax
             DO i = 1,  ibMaxPerJB(j)
                sclim2(i,k,im,j) = REAL(brf8(i,j),KIND=r8)*1000
             END DO
          END DO
       END DO
    END DO
    CLOSE(2)

    IF(.not.RESTART)THEN

!   determine the layer in which the mixed layer resides
!   force the mixed layer to be deeper than dz(1) and shallower
!   than zlev(lbottom)-0.5*dz(lbottom)
!   determine the layer in which the mixed layer resides
    laymld=0
    DO k=1,km
       DO j = 1, jbMax
          DO i = 1,  ibMaxPerJB(j) 
                 laymld(i,j) = 5
             !IF(zlev(k) < hslab)THEN
             !       laymld(i,j) = laymld(i,j) + 1
             !END IF
       END DO
    END DO
    END DO
    im=idatec(2)
       DO k=1,km
          DO j = 1, jbMax
             DO i = 1,  ibMaxPerJB(j)
                hbath(i,j)  =hslab
                tdeep(i,k,j)= tclim2(i,k,im,j)
                sdeep(i,k,j)= sclim2(i,k,im,j)
             END DO
          END DO
       END DO
    END IF
  END SUBROUTINE init_MLO
  ! -------------------------------------------------------------


  SUBROUTINE driver(&
                    idatec   ,&!INTEGER      , INTENT(IN   ) :: idatec(4)
                    nCols, kMax   ,&!INTEGER      , INTENT(IN   ) :: nCols
                    latco    ,&!INTEGER      , INTENT(IN   ) :: latco
                    timestep ,&!REAL(KIND=r8), INTENT(IN   ) :: timestep
                    imask    ,&!INTEGER (KIND=i8), INTENT(IN   ) :: imask (nCols)
                    xlat     ,&!REAL(KIND=r8), INTENT(IN   ) :: xlat (nCols) !  -90-90  radian
                    xevap    ,&!REAL(KIND=r8), INTENT(IN   ) :: xevap(nCols) !   xevap   - surface water flux (evaporation) (kg/m**2/s)  
                    xsw      ,&!REAL(KIND=r8), INTENT(IN   ) :: xsw  (nCols) !   xsw - net downward surface short wave radiation (w/m**2)
                    xlw      ,&!REAL(KIND=r8), INTENT(IN   ) :: xlw  (nCols) !   xlw - net upward surface long wave radiation (w/m**2)
                    xsh      ,&!REAL(KIND=r8), INTENT(IN   ) :: xsh  (nCols) !   xsh - surface sensible heat flux  (w/m**2) 
                    xlh      ,&!REAL(KIND=r8), INTENT(IN   ) :: xlh  (nCols) !   xlh - surface latent heat flux  (w/m**2)
                    fracice  ,&!REAL(KIND=r8), INTENT(IN   ) :: fracice (nCols)  !ice fraction, 0-1
                    taux     ,&!REAL(KIND=r8), INTENT(IN   ) :: taux  (nCols)  !REAL(KIND=r8),    INTENT(in) :: taux0     !surface zonal wind stress, n/m^2   
                    tauy     ,&!REAL(KIND=r8), INTENT(IN   ) :: tauy  (nCols)  !REAL(KIND=r8),    INTENT(in) :: tauy0     !surface meridional wind stress, n/m^2   
                    xprecc   ,&!REAL(KIND=r8), INTENT(IN   ) ::  xprecc(nCols)!   xprecc  - convective precipitation (mm/day)
                    xprecl   ,&!REAL(KIND=r8), INTENT(IN   ) ::  xprecl(nCols)!   xprecl  - large-scale precipitation (mm/day)
                    dump     ,&
                    sst0      )!REAL(KIND=r8), INTENT(INOUT) :: sst0(nCols)

    INTEGER      , INTENT(IN   ) :: idatec(4)
    INTEGER      , INTENT(IN   ) :: nCols,kMax
    INTEGER      , INTENT(IN   ) :: latco
    REAL(KIND=r8), INTENT(IN   ) :: timestep
    INTEGER (KIND=i8), INTENT(IN   ) :: imask (nCols)
    REAL(KIND=r8), INTENT(IN   ) :: xlat (nCols) !  -90-90  radian
    REAL(KIND=r8), INTENT(IN   ) :: xevap(nCols) !   xevap   - surface water flux (evaporation) (kg/m**2/s)  
    REAL(KIND=r8), INTENT(IN   ) :: xsw  (nCols) !   xsw     - net downward surface short wave radiation (w/m**2)
    REAL(KIND=r8), INTENT(IN   ) :: xlw  (nCols) !   xlw - net upward surface long wave radiation (w/m**2)
    REAL(KIND=r8), INTENT(IN   ) :: xsh  (nCols) !   xsh - surface sensible heat flux  (w/m**2) 
    REAL(KIND=r8), INTENT(IN   ) :: xlh  (nCols) !   xlh - surface latent heat flux  (w/m**2)
    REAL(KIND=r8), INTENT(IN   ) :: fracice (nCols)  !ice fraction, 0-1
    REAL(KIND=r8), INTENT(IN   ) :: taux  (nCols)     !REAL(KIND=r8),    INTENT(in) :: taux0     !surface zonal wind stress, n/m^2   
    REAL(KIND=r8), INTENT(IN   ) :: tauy  (nCols)     !REAL(KIND=r8),    INTENT(in) :: tauy0     !surface meridional wind stress, n/m^2   
    REAL(KIND=r8), INTENT(IN   ) ::  xprecc(nCols)!   xprecc  - convective precipitation (mm/day)
    REAL(KIND=r8), INTENT(IN   ) ::  xprecl(nCols)!   xprecl  - large-scale precipitation (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: dump (1:nCols,1:kMax ) 
    REAL(KIND=r8), INTENT(INOUT) :: sst0(nCols)


    !LOCAL 

    REAL(KIND=r8) ::   evap(nCols)!   evap   - evaporation (kg/(m**2*sec)) 
    REAL(KIND=r8) ::   precip(nCols)!   precip - precipitation (kg/(m**2*sec))
    REAL(KIND=r8) ::   qao (nCols)!   qao    = qsw - (qlw + qsh + qlh)
    REAL(KIND=r8) ::   qsw (nCols) !   qsw    - short wave radiation absorbed by the ocean (w/m**2)
    REAL(KIND=r8) ::   qcor(nCols)!   qcor   - surface heat flux correction (w/m**2) for each timestep
    REAL(KIND=r8) ::   scor(nCols)   !   scord  - daily mean surface fresh water correction (1/1000)
    REAL(KIND=r8) ::   ff(nCols)
    REAL(KIND=r8) ::  qtmp
    REAL(KIND=r8) ::  stmp,ratiot
    REAL(KIND=r8)     :: htmp(nCols)
    REAL(KIND=r8)     :: h_aux       (nCols)        !REAL(KIND=r8),    INTENT(inout) :: h     !mixed-layer depth
    REAL(KIND=r8)     :: t_aux       (nCols,1:km)   !REAL(KIND=r8),    INTENT(inout) :: t(km)   !temperature profile, K
    REAL(KIND=r8)     :: dtemp_aux       (nCols,1:km)   !REAL(KIND=r8),    INTENT(inout) :: t(km)   !temperature profile, K
    REAL(KIND=r8)     :: s_aux       (nCols,1:km)   !REAL(KIND=r8),    INTENT(inout) :: s(km)   !salinity profile, 1/1000
    INTEGER           :: laymld_aux  (nCols)        !INTEGER      ,     INTENT(inout) :: laymix    !layer in which the mixed-layer resides

    INTEGER       :: i,l,imon,k
    dtemp_aux(1:nCols,1)=0
    DO i=1,nCols
       DO k=2,km
          dtemp_aux(i,k)=(tdeep(i,k-1,latco)-tdeep(i,k,latco))/dz(k-1)
 
!           dtemp_aux(i,k)=(sst0(i,k-1,latco)-t(i,k,latco))/dz(k-1)

!           dtemp_aux(i,k)*dz(k-1)=(sst0(i,k-1,latco)-t(i,k,latco))

!           t(i,k,latco)=sst0(i,k-1,latco)-dtemp_aux(i,k)*dz(k-1)
       END DO
    ENDDO

    DO i=1,nCols
          evap(i)       =xevap(i)                                     !kg/(m**2*s) 
          qao(i)        =xsw(i)-(xlw(i)+xsh(i)+xlh(i))                !w/m**2
          qsw(i)        =xsw(i)                                       !w/m**2
          precip(i)     =(xprecc(i)+xprecl(i))                        !kg/(m**2*s)
          ff(i)         = 2.0_r8 * 7.2722e-5_r8 * SIN(xlat(i))
          qcor(i)       = 0.0
          scor(i)       = 0.0
          tdeep(i,1,latco)  = sst0(i)
    ENDDO

    DO i=1,nCols
       DO k=2,km
!           t(i,k,latco)=t(i,k-1,latco) - dtemp_aux(i,k)*dz(k-1)
       END DO
    ENDDO

    imon=idatec(2)
    !
    !
    !   for regions with depths less than hslab, treat as a hslab slab ocean
    !   for temperature; salinity is not computed and set to missing

    DO i=1,nCols
       IF(imask(i) <=0)THEN
          h_aux     (i)      =hbath  (i,latco)  
          t_aux     (i,1:km) =tdeep  (i,1:km,latco)  
          s_aux     (i,1:km) =sdeep  (i,1:km,latco)  
          laymld_aux(i)      =laymld  (i,latco)  
       IF (zlev(lbottom(i,latco)) .LE. hslab) THEN
          h_aux(i)=hslab
          t_aux(i,1) = t_aux(i,1) + timestep * (qao(i)+qcor(i))*(1-fracice(i))/(rhoref * capa * hslab)
          DO l=1, lbottom(i,latco)-1
             t_aux(i,l) = t_aux(i,1)
          END DO
       ELSE  
          CALL mlo_main(km              ,&!INTEGER      ,    INTENT(in) :: km      !=30     !model vertical layers
               timestep                 ,&!REAL(KIND=r8),    INTENT(in) :: timestep     !timestep, seconds
               dz      (1:km)           ,&!REAL(KIND=r8),    INTENT(in) :: dz0(km)     !layer thickness, m
               zlev    (1:km)           ,&!REAL(KIND=r8),    INTENT(in) :: zlev0(km)    !layer depth, m
               lbottom (i,latco)        ,&!INTEGER      ,    INTENT(in) :: lbot0     !index for the bottom of the ocean     
               r       (i,latco)        ,&!REAL(KIND=r8),    INTENT(in) :: rr0     !ocean optical properties
               d1      (i,latco)        ,&!REAL(KIND=r8),    INTENT(in) :: d10     !ocean optical properties
               d2      (i,latco)        ,&!REAL(KIND=r8),    INTENT(in) :: d20     !ocean optical properties
               ff      (i)              ,&!REAL(KIND=r8),    INTENT(in) :: coriolis0    !coriolis parameter f
               qao     (i)              ,&!REAL(KIND=r8),    INTENT(in) :: qtot0     !surface downward total heat flux, w/m^2
               qsw     (i)              ,&!REAL(KIND=r8),    INTENT(in) :: qsw0     !surface downward solar flux, w/m^2
               taux    (i)              ,&!REAL(KIND=r8),    INTENT(in) :: taux0     !surface zonal wind stress, n/m^2   
               tauy    (i)              ,&!REAL(KIND=r8),    INTENT(in) :: tauy0     !surface meridional wind stress, n/m^2   
               precip  (i)              ,&!REAL(KIND=r8),    INTENT(in) :: precip0     !precipitation rate, kg/m^2/s
               evap    (i)              ,&!REAL(KIND=r8),    INTENT(in) :: evap0     !evaporation rate, kg/m^2/s
               qcor    (i)              ,&!REAL(KIND=r8),    INTENT(in) :: qcor0     !heat flux correction term, w/m^2
               scor    (i)              ,&!REAL(KIND=r8),    INTENT(in) :: scor0     !salinity flux correction term, 1/1000
               tclim2  (i,1:km,imon,latco)   ,&!REAL(KIND=r8),    INTENT(in) :: tclim0(km)   !temperature profile for relaxation, K
               sclim2  (i,1:km,imon,latco)   ,&!REAL(KIND=r8),    INTENT(in) :: sclim0(km)   !salinity profile for relaxation, K
               fracice (i)              ,&!REAL(KIND=r8),    INTENT(in) :: fracice     !ice fraction, 0-1
               h_aux   (i)              ,&!REAL(KIND=r8),    INTENT(inout) :: h     !mixed-layer depth
               t_aux   (i,1:km)         ,&!REAL(KIND=r8),    INTENT(inout) :: t(km)     !temperature profile, K
               s_aux   (i,1:km)         ,&!REAL(KIND=r8),    INTENT(inout) :: s(km)     !salinity profile, 1/1000
               laymld_aux  (i)         )!INTEGER     ,     INTENT(inout) :: laymix    !layer in which the mixed-layer resides
       END IF

       !   calculate and accumulate heat(w/m**2) and salt (kg/m**2*sec) flux corrections
       !IF (icor .EQ. 0) THEN
       !fcord(i,j)=fcord(i,j,ijday)+(fracice(i,j)-fice(i,j))
       !
       ! -- to avoid unrealistic large heat flux in rare extreme case, 
       !    set a cap to qcord to confine sst diference within 1K in a day.
       !    note: if h=50, and (sst-t)=1.0 then qcord=2342.0 in one day
       ratiot=2342.0_r8*h_aux(i)/hslab          
       qtmp=(sst0(i)-t_aux(i,1))*(rhoref*capa*h_aux(i))/timestep
       qtmp=MAX(-ratiot,MIN(ratiot,qtmp))
       qcor(i)=qcor(i)+qtmp
       !
       IF (zlev(lbottom(i,latco)) .GT. hslab) THEN 
          stmp=((sclim2(i,1,imon,latco)-s_aux(i,1))/s_aux(i,1))*(rhoref*h_aux(i))/timestep
          scor(i) = scor(i)+stmp
       END IF
       END IF
    END DO

    !
    !
    !   for regions with depths less than hslab, treat as a hslab slab ocean
    !   for temperature; salinity is not computed and set to missing
    DO i=1,nCols
       IF(imask(i) <=0)THEN
          IF (zlev(lbottom(i,latco)) .LE. hslab) THEN
             hbath (i,latco) = hslab
             tdeep (i,1,latco)   = tdeep(i,1,latco) + timestep * (qao(i)+qcor(i))*(1-fracice(i))/(rhoref * capa * hslab)
             DO l=1, lbottom(i,latco)-1
                tdeep(i,l,latco) = tdeep(i,1,latco)
             END DO
          ELSE  
             CALL mlo_main(km              ,&!INTEGER      ,    INTENT(in) :: km      !=30     !model vertical layers
               timestep                 ,&!REAL(KIND=r8),    INTENT(in) :: timestep     !timestep, seconds
               dz      (1:km)           ,&!REAL(KIND=r8),    INTENT(in) :: dz0(km)     !layer thickness, m
               zlev    (1:km)           ,&!REAL(KIND=r8),    INTENT(in) :: zlev0(km)    !layer depth, m
               lbottom (i,latco)        ,&!INTEGER      ,    INTENT(in) :: lbot0     !index for the bottom of the ocean     
               r       (i,latco)        ,&!REAL(KIND=r8),    INTENT(in) :: rr0     !ocean optical properties
               d1      (i,latco)        ,&!REAL(KIND=r8),    INTENT(in) :: d10     !ocean optical properties
               d2      (i,latco)        ,&!REAL(KIND=r8),    INTENT(in) :: d20     !ocean optical properties
               ff      (i)              ,&!REAL(KIND=r8),    INTENT(in) :: coriolis0    !coriolis parameter f
               qao     (i)              ,&!REAL(KIND=r8),    INTENT(in) :: qtot0     !surface downward total heat flux, w/m^2
               qsw     (i)              ,&!REAL(KIND=r8),    INTENT(in) :: qsw0     !surface downward solar flux, w/m^2
               taux    (i)              ,&!REAL(KIND=r8),    INTENT(in) :: taux0     !surface zonal wind stress, n/m^2   
               tauy    (i)              ,&!REAL(KIND=r8),    INTENT(in) :: tauy0     !surface meridional wind stress, n/m^2   
               precip  (i)              ,&!REAL(KIND=r8),    INTENT(in) :: precip0     !precipitation rate, kg/m^2/s
               evap    (i)              ,&!REAL(KIND=r8),    INTENT(in) :: evap0     !evaporation rate, kg/m^2/s
               qcor    (i)              ,&!REAL(KIND=r8),    INTENT(in) :: qcor0     !heat flux correction term, w/m^2
               scor    (i)              ,&!REAL(KIND=r8),    INTENT(in) :: scor0     !salinity flux correction term, 1/1000
               tclim2  (i,1:km,imon,latco)   ,&!REAL(KIND=r8),    INTENT(in) :: tclim0(km)   !temperature profile for relaxation, K
               sclim2  (i,1:km,imon,latco)   ,&!REAL(KIND=r8),    INTENT(in) :: sclim0(km)   !salinity profile for relaxation, K
               fracice (i)              ,&!REAL(KIND=r8),    INTENT(in) :: fracice     !ice fraction, 0-1
               hbath   (i,latco)        ,&!REAL(KIND=r8),    INTENT(inout) :: h     !mixed-layer depth
               tdeep   (i,1:km,latco)   ,&!REAL(KIND=r8),    INTENT(inout) :: t(km)     !temperature profile, K
               sdeep   (i,1:km,latco)   ,&!REAL(KIND=r8),    INTENT(inout) :: s(km)     !salinity profile, 1/1000
               laymld  (i,latco)         )!INTEGER     ,     INTENT(inout) :: laymix    !layer in which the mixed-layer resides
          END IF

       END IF
       sst0(i)= MIN(MAX(tdeep(i,1,latco),sst0(i)),320.0_r8)

       !dump(i,10)= lbottom(i,latco)
       !dump(i,11)= r (i,latco)
       !dump(i,12)= d1(i,latco) 
       !dump(i,13)= d2(i,latco) 
       !dump(i,14)= hbath(i,latco)
       !dump(i,15)= qao(i) 
       !dump(i,16)= qsw(i) 
       !dump(i,17)= taux(i) 
       !dump(i,18)= tauy(i) 
       !dump(i,19)= precip(i) 
       !dump(i,20)= evap(i) 
       !dump(i,21)= qcor(i) 
       !dump(i,22)= scor(i) 
       !dump(i,23)= tclim2 (i,1,imon,latco)
       !dump(i,24)= tclim2 (i,km,imon,latco)
       !dump(i,25)= sclim2 (i,1,imon,latco)
       !dump(i,26)= sclim2 (i,km,imon,latco)
       !dump(i,27)= sst0(i)

    END DO
  END SUBROUTINE driver

  ! ----------------------------------------------------------------
  SUBROUTINE mlo_main(km          , &!INTEGER      ,    INTENT(in) :: km           !=30          !model vertical layers
                      timestep    , &!REAL(KIND=r8),    INTENT(in) :: timestep     !timestep, seconds
                      dz0         , &!REAL(KIND=r8),    INTENT(in) :: dz0(km)      !layer thickness, m
                      zlev0       , &!REAL(KIND=r8),    INTENT(in) :: zlev0(km)    !layer depth, m
                      lbot0       , &!INTEGER      ,    INTENT(in) :: lbot0        !index for the bottom of the ocean      
                      rr0         , &!REAL(KIND=r8),    INTENT(in) :: rr0          !ocean optical properties
                      d10         , &!REAL(KIND=r8),    INTENT(in) :: d10          !ocean optical properties
                      d20         , &!REAL(KIND=r8),    INTENT(in) :: d20          !ocean optical properties
                      coriolis0   , &!REAL(KIND=r8),    INTENT(in) :: coriolis0    !coriolis parameter f
                      qtot0       , &!REAL(KIND=r8),    INTENT(in) :: qtot0        !surface downward total heat flux, w/m^2
                      qsw0        , &!REAL(KIND=r8),    INTENT(in) :: qsw0         !surface downward solar flux, w/m^2
                      taux0       , &!REAL(KIND=r8),    INTENT(in) :: taux0        !surface zonal wind stress, n/m^2   
                      tauy0       , &!REAL(KIND=r8),    INTENT(in) :: tauy0        !surface meridional wind stress, n/m^2   
                      precip0     , &!REAL(KIND=r8),    INTENT(in) :: precip0      !precipitation rate, kg/m^2/s
                      evap0       , &!REAL(KIND=r8),    INTENT(in) :: evap0        !evaporation rate, kg/m^2/s
                      qcor0       , &!REAL(KIND=r8),    INTENT(in) :: qcor0        !heat flux correction term, w/m^2
                      scor0       , &!REAL(KIND=r8),    INTENT(in) :: scor0        !salinity flux correction term, 1/1000
                      tclim0      , &!REAL(KIND=r8),    INTENT(in) :: tclim0(km)   !temperature profile for relaxation, K
                      sclim0      , &!REAL(KIND=r8),    INTENT(in) :: sclim0(km)   !salinity profile for relaxation, K
                      fracice     , &!REAL(KIND=r8),    INTENT(in) :: fracice      !ice fraction, 0-1
                      h           , &!REAL(KIND=r8),    INTENT(inout) :: h         !mixed-layer depth
                      t           , &!REAL(KIND=r8),    INTENT(inout) :: t(km)     !temperature profile, K
                      s           , &!REAL(KIND=r8),    INTENT(inout) :: s(km)     !salinity profile, 1/1000
                      laymix        )!INTEGER      ,    INTENT(inout) :: laymix    !layer in which the mixed-layer resides
    ! ----------------------------------------------------------------
    ! Fanglin Yang, July 2003
    ! Purpose: 1-d mixed-layer ocean model, based on gaspar (jpo, 1988).
    ! note:
    ! 1. within the mixed layer, temperature and salinity are prognostics.
    !    mixed-layer depth is prognostics when entrainment occurs, and is
    !    diagnostics when mixed layer shoals. 
    ! 2. in layers below the mixed layer, vertical diffusion, convective 
    !    adjustment and solar radiation pentration are included. in 
    !    addition, temperature and salinity are also relaxed to observed 
    !    climatologies in a 10-year time scale.

    !
    ! vertical structure of the model
    !   --------------------------- surface
    !
    !   --------------------------- k=1, dz(1), zlev(1)
    !        
    !   --------------------------- k=2, dz(2), zlev(2)
    !              .
    !              .
    !              .
    !   --------------------------- 
    !     ^^^ h, bottom of mixed layer 
    !   --------------------------- k=laymix
    !              .
    !              .
    !              .
    !   --------------------------- k=km, dz(km), zlev(km)
    !


    IMPLICIT NONE
    !      include 'mlo.h'
    REAL(KIND=r8), PARAMETER :: rhoref = 1024.438_r8    !sea water reference density, kg/m^3
    REAL(KIND=r8), PARAMETER :: visct=1.e-5_r8          !viscocity for temperature diffusion
    REAL(KIND=r8), PARAMETER :: viscs=1.e-5_r8          !viscocity for salt diffusion
    REAL(KIND=r8), PARAMETER :: grav =9.81_r8           !gravity, kg/m/s^2
    REAL(KIND=r8), PARAMETER :: capa =3950.0_r8         !heat capacity of sea water 
    REAL(KIND=r8), PARAMETER :: smax=50.0_r8             !normal maximum salt
    REAL(KIND=r8), PARAMETER :: smin=1.0_r8             !normal minimal salt
    REAL(KIND=r8), PARAMETER :: tmin=2.68E+02_r8        !normal minimal temp
    REAL(KIND=r8), PARAMETER :: tmax=3.11E+02_r8        !normal max temp

    ! input
    INTEGER      ,    INTENT(in) :: km                !=30            !model vertical layers
    REAL(KIND=r8),    INTENT(in) :: timestep          !timestep, seconds
    REAL(KIND=r8),    INTENT(in) :: dz0(km)           !layer thickness, m
    REAL(KIND=r8),    INTENT(in) :: zlev0(km)         !layer depth, m
    INTEGER      ,    INTENT(in) :: lbot0             !index for the bottom of the ocean      
    REAL(KIND=r8),    INTENT(in) :: rr0               !ocean optical properties
    REAL(KIND=r8),    INTENT(in) :: d10               !ocean optical properties
    REAL(KIND=r8),    INTENT(in) :: d20               !ocean optical properties
    REAL(KIND=r8),    INTENT(in) :: coriolis0         !coriolis parameter f
    REAL(KIND=r8),    INTENT(in) :: qtot0             !surface downward total heat flux, w/m^2
    REAL(KIND=r8),    INTENT(in) :: qsw0              !surface downward solar flux, w/m^2
    REAL(KIND=r8),    INTENT(in) :: taux0             !surface zonal wind stress, n/m^2   
    REAL(KIND=r8),    INTENT(in) :: tauy0             !surface meridional wind stress, n/m^2   
    REAL(KIND=r8),    INTENT(in) :: precip0           !precipitation rate, kg/m^2/s
    REAL(KIND=r8),    INTENT(in) :: evap0             !evaporation rate, kg/m^2/s
    REAL(KIND=r8),    INTENT(in) :: qcor0             !heat flux correction term, w/m^2
    REAL(KIND=r8),    INTENT(in) :: scor0             !salinity flux correction term, 1/1000
    REAL(KIND=r8),    INTENT(in) :: tclim0(km)        !temperature profile for relaxation, K
    REAL(KIND=r8),    INTENT(in) :: sclim0(km)        !salinity profile for relaxation, K
    REAL(KIND=r8),    INTENT(in) :: fracice           !ice fraction, 0-1
    ! input & output
    REAL(KIND=r8),    INTENT(inout) :: h              !mixed-layer depth
    REAL(KIND=r8),    INTENT(inout) :: t(km)          !temperature profile, K
    REAL(KIND=r8),    INTENT(inout) :: s(km)          !salinity profile, 1/1000
    INTEGER      ,    INTENT(inout) :: laymix         !layer in which the mixed-layer resides
    ! local 
    REAL(KIND=r8) :: tmix, smix, taumag
    INTEGER :: k
    REAL(KIND=r8)  ::  dz(km)             !depth, from top down, m 
    REAL(KIND=r8)  :: zlev(km)                !sum(dz(k)), top down, m
    REAL(KIND=r8)  :: tclim(km)           !temperature profile for relaxation, K
    REAL(KIND=r8)  :: sclim(km)           !salinity profile for relaxation, K
    REAL(KIND=r8) :: coriolis            !coriolis parameter f  
    REAL(KIND=r8) :: d1                  !water optical coefficient
    REAL(KIND=r8) :: d2                  !water optical coefficient
    REAL(KIND=r8) :: evap                !evaporation rate, kg/m^2/s
    INTEGER ::  lbot            !level of ocean botom
    REAL(KIND=r8) :: precip              !precipitation rate, kg/m^2/s
    REAL(KIND=r8) :: qcor                !heat flux correction, w/m^2  
    REAL(KIND=r8) :: qsw                 !surface net downward solar flux, w/m^2
    REAL(KIND=r8) :: qtot                !surface total downward heat flux, w/m^2
    REAL(KIND=r8) :: ustar               !frictional velocity, m/s
    REAL(KIND=r8) :: taux                !surface wind stress, n/m^2
    REAL(KIND=r8) :: tauy                !surface wind stress, n/m^2
    REAL(KIND=r8) :: rr                  !water optical coefficient
    REAL(KIND=r8) :: scor                !water flux correction, 1/1000  

    !-----------------------------------------------------------------------------

    ! passing dummy variables that are defined in mlo.h
    DO k=1,km
       dz(k)=dz0(k)
       zlev(k)=zlev0(k)
       tclim(k)=tclim0(k)
       sclim(k)=sclim0(k)
    ENDDO
    lbot  =lbot0
    rr    =rr0
    d1    =d10
    d2    =d20
    coriolis=coriolis0
    qtot  =qtot0
    qsw   =qsw0
    precip=precip0
    evap  =evap0
    qcor  =qcor0
    scor  =scor0
    taux  =taux0
    tauy  =tauy0
    !
    ! insure a large enough u* which is obtained by assuming 
    ! the wind speed being greater than ~ 1 m/sec equivalent 
    ! to tau magnitude > 1.6e-3. ustar value under ice from 
    ! fichefet and gaspar (1988 jpo pp 184)
    taumag = MAX(SQRT(taux**2 + tauy**2),0.001625_r8)
    ustar = 0.006_r8 * fracice + SQRT(taumag/rhoref)*(1-fracice)

    tmix = t(1)
    smix = s(1)
    DO k=1, laymix-1
       t(k) = tmix
       s(k) = smix
    END DO

    ! ----------------------------------------------------------------
    !      if(iwrite.eq.1) then
    !       write(99,*) "---------------- MLO_MAIN.F -------------------"
    !       write(99,*) "year,month,day,lon,lat",iyear, imonth, iday, ilon, jlat
    !        write(99,*) "dz(k),zlev(k),tclim(k),sclim(k),t(k),s(k)"
    !       do k=1,lbot
    !        write(99,'(i4,6f10.2)') k,dz(k),zlev(k),tclim(k),sclim(k),t(k),s(k)
    !       enddo
    !       write(99,*) "laymix, h, tmix, smix, fracice"
    !       write(99,'(i4,6f10.2)')laymix,  h, tmix, smix, fracice
    !       write(99,*) "lbot, rr, d1, d2, coriolis"
    !       write(99,'(i4,3f10.2,e12.4)') lbot, rr, d1, d2, coriolis
    !       write(99,*)"qtot, qsw, taux, tauy, precip, evap"
    !       write(99,'(2f10.2, 4e12.4)')qtot, qsw, taux, tauy, precip, evap
    !       write(99,*)"ustar, qcor, scor"
    !       write(99,'(3f15.6)') ustar, qcor, scor
    !      endif
    !
    ! ----------------------------------------------------------------
    ! to represent the impact of currents and horizontal eddies
    ! in layers below the mixed layer, relaxtion is applied.
    CALL relax (km,lbot,timestep,tclim,sclim,&
         laymix,tmix, smix, t, s)

    ! convective adjustment must be called before entrainment 
    ! calculation to keep a stable stratification. this is because
    ! only when  deltab>0,  the assumption ap>0 <==> we>0 be true.      
    CALL convect(km,lbot,dz,zlev,&
         h, tmix, smix, laymix, t, s)

    ! changes in temperature, salinity and mixed-layer depth due to entrainment 
    CALL entrain(km,lbot,timestep,grav,capa,rhoref,smin  ,smax  ,tmin  ,tmax, &
         coriolis,d1,d2,evap,precip,zlev,dz,ustar,qsw,qtot,rr,&
         h, laymix, tmix, smix, t, s)
    ! changes in temperature and salinity due to
    ! surface heat and water flux, including flux corrections
    ! and penetration of solar radiation
    CALL surflux(km,lbot,capa,d1,d2,evap,precip,qcor,qsw,qtot,rhoref,&
         rr,scor,timestep,zlev,dz,&
         h, tmix, smix, laymix, t, s)
    !
    ! changes in temprature and salinity due to vertical diffusion

    CALL diffusion(h, laymix, km, lbot,   &
         visct, timestep, dz, zlev, t)

    CALL diffusion(h, laymix, km, lbot,   &
         viscs, timestep, dz, zlev, s)
    !
    tmix=t(1)
    smix=s(1)
    !      if(iwrite.eq.1) write(99,*) "---------------- MLO_MAIN.F -------------------"
    !       do k=1,laymix
    !        if(t(k).gt.tmax .or. t(k).lt.tmin) then
    !         if(iwrite.eq.1) write(99,111) iyear, imonth, iday, ilon, jlat, k, t(k)
    !         write(6,111) iyear, imonth, iday, ilon, jlat, k, t(k)
    !         endif
    !        enddo
    !        do k=1,laymix
    !         if(s(k).gt.smax .or. s(k).lt.smin) then
    !          if(iwrite.eq.1) write(99,222) iyear, imonth, iday, ilon, jlat, k, s(k)
    !          write(6,222) iyear, imonth, iday, ilon, jlat, k, s(k)
    !         endif
    !        enddo
111 FORMAT("WARN t(k)", 6i6,f10.2)
222 FORMAT("WARN s(k)", 6i6,f10.2)

  END SUBROUTINE mlo_main

  ! ------------------------------------------------------
  SUBROUTINE entrain(km,lbot,timestep,grav,capa,rhoref,smin  ,smax  ,tmin  ,tmax, &
       coriolis,d1,d2,evap,precip,zlev,dz,ustar,qsw,qtot,rr,&
       h, laymix, tmix, smix, t, s)
    ! ------------------------------------------------------
    !  purpose: update mixed-layer depth, temperature and salinity
    !           du to mixed-layer deepening/entrainment  or shoaling

    IMPLICIT NONE
    !      include 'mlo.h'

    ! input
    INTEGER, INTENT(in) :: km!, parameter :: km =30            !model vertical layers
    INTEGER , INTENT(in) :: lbot            !level of ocean botom
    REAL(KIND=r8),     INTENT(in) :: timestep!=86400.0  !timestep, seconds

    REAL(KIND=r8), INTENT(in) :: rhoref != 1024.438    !sea water reference density, kg/m^3
    REAL(KIND=r8), INTENT(in) :: smin   !=1.0             !normal minimal salt
    REAL(KIND=r8), INTENT(in) :: smax   !=50.             !normal maximum salt
    REAL(KIND=r8), INTENT(in) :: tmin   !=2.68E+02        !normal minimal temp
    REAL(KIND=r8), INTENT(in) :: tmax   !=3.11E+02        !normal max temp

    REAL(KIND=r8),    INTENT(in) :: grav !=9.81           !gravity, kg/m/s^2
    REAL(KIND=r8),    INTENT(in) :: capa !, parameter :: capa =3950.0         !heat capacity of sea water 
    REAL(KIND=r8),    INTENT(in) :: coriolis            !coriolis parameter f  
    REAL(KIND=r8),    INTENT(in) :: d1 !water optical coefficient
    REAL(KIND=r8),    INTENT(in) :: d2 !water optical coefficient
    REAL(KIND=r8),    INTENT(in) :: evap                !evaporation rate, kg/m^2/s
    REAL(KIND=r8),    INTENT(in) :: precip                !precipitation rate, kg/m^2/s
    REAL(KIND=r8),    INTENT(in) :: zlev(km)                !sum(dz(k)), top down, m
    REAL(KIND=r8),    INTENT(in)  ::  dz(km)             !depth, from top down, m 
    REAL(KIND=r8),    INTENT(in)  :: ustar               !frictional velocity, m/s
    REAL(KIND=r8),    INTENT(in)  :: qsw                 !surface net downward solar flux, w/m^2
    REAL(KIND=r8),    INTENT(in)  :: qtot                !surface total downward heat flux, w/m^2
    REAL(KIND=r8),    INTENT(in)   :: rr                  !water optical coefficient

    REAL(KIND=r8),    INTENT(inout) :: h        !mixed-layer depth, m 
    INTEGER, INTENT(inout) :: laymix   !layer in which mixed layer resides
    REAL(KIND=r8),    INTENT(inout) :: tmix     !mixed-layer temperature, K
    REAL(KIND=r8),    INTENT(inout) :: smix     !mixed-layer salinity, 1/1000
    REAL(KIND=r8),    INTENT(inout) :: t(km)    !temperature profile, K
    REAL(KIND=r8),    INTENT(inout) :: s(km)    !salinity profile, 1/1000
    ! local variables
    REAL(KIND=r8) :: tbelow   !below mixed-layer temperature, K
    REAL(KIND=r8) :: sbelow   !below mixed-layer salinity, 1/1000
    REAL(KIND=r8) :: rho(km)  !density profile, kg/m^3
    REAL(KIND=r8) :: ap, alpha, beta, we, rhomix
    REAL(KIND=r8) :: xkm1, xk, xkp1, apkm1, apk, apkp1
    REAL(KIND=r8) :: tnew(km), snew(km)
    REAL(KIND=r8) :: hnew, hmax, hmin, tmixnew, smixnew
    REAL(KIND=r8) :: allmass, allheat, allsalt
    INTEGER       :: laynew, k, count
 tbelow  =0.0_r8 !below mixed-layer temperature, K
 sbelow  =0.0_r8  !below mixed-layer salinity, 1/1000
 rho(1:km) =0.0_r8  !density profile, kg/m^3
 ap=0.0_r8 ; alpha=0.0_r8 ; beta=0.0_r8 ; we=0.0_r8 ; rhomix=0.0_r8 
 xkm1=0.0_r8 ; xk=0.0_r8 ; xkp1=0.0_r8 ; apkm1=0.0_r8 ; apk=0.0_r8 ; apkp1=0.0_r8 
 tnew(1:km)=0.0_r8 ; snew(1:km)=0.0_r8 
 hnew=0.0_r8 ; hmax=0.0_r8 ; hmin=0.0_r8 ; tmixnew=0.0_r8 ; smixnew=0.0_r8 
 allmass=0.0_r8 ; allheat=0.0_r8 ; allsalt=0.0_r8 
    ! ----------------------------------------------------------
    ! in rare case after convective adjustment mixed layer may 
    ! reside in the bottom of the ocean, then tbelow=tmix and 
    ! sbelow=smix. to prevent from deltab=0, which causes 
    ! overflowing in apwe.F, tbelow and sbelow are slightly 
    ! modified to maintain a stable stratification.    
    tbelow=MIN(t(laymix),tmix-1.e-3_r8)       !temperature immediately below mlo
    sbelow=MAX(s(laymix),smix+1.e-5_r8)       !salinity immediately below mlo

    ! initialize new sate variables
    laynew=laymix
    hnew=h
    tmixnew=tmix
    smixnew=smix
    DO k=1,lbot
       tnew(k)=t(k)
       snew(k)=s(k)
    ENDDO
    hmax=zlev(lbot)-0.5_r8*dz(lbot)         !maximum mixed-layer allowed

    !        if(iwrite.eq.1) then
    !         write(99,*)"-------------- ENTRAIN.F -------------------"
    !         write(99,*)"input: h, tmix, tbelow, smix, sbelow"
    !         write(99,'(6f12.4)') h, tmix, tbelow, smix, sbelow
    !        endif

    ! compute thermal expansion and saline contraction coefficients
    ! alpha=d(rho)/d(tmix)/rhoref, beta=d(rho)/d(smix)/rhoref
    CALL rhocoef(tmix, MAX(0.0_r8,smix), rhoref, alpha, beta)

    ! initial call to apwe, compute ap and entrainment we 
    CALL apwe(capa,grav,rhoref,coriolis,ustar,d1,d2,evap,precip,qsw,qtot,rr,&
         h, tmix, smix, tbelow, sbelow,  &
         alpha, beta, ap, 1, we)
    !      call apwe(h, tmix, smix, tbelow, sbelow, &
    !                alpha, beta, ap, 1, we)

    ! ----------------------------------------------------------
    ! ----------------------------------------------------------
    ! if ap>0, mixed-layer deepens
    IF(ap.GT.0.0_r8 .AND. h.LT.hmax ) THEN

       !      if(we.le.0.0_r8) then
       ! note: in rare cases, when initialized from observations,
       !       such as in winter over Caspin Sea, the area
       !       is covered by ice, tmix=tice(271.36), smix is smaller
       !       than 10. then alpha<0 and beta>0. moreover, both temperature 
       !       and salinity increase downward and hence rho increases downward. 
       !       convective adjustment will not adjust the profiles.  in this
       !       case, deltab<0 occurs and the condition ap>0 <==> we>0 is violated.
       !       entrainment is then skipped.
       !        if(iwrite.eq.1) then
       !        write(99,*) "inconsistence: ap>0, we<=0., skip"
       !        write(99,*) "ap, we, tmix, tbelow, smix, sbelow, alpha, beta"
       !        write(99,'(8e12.4)') ap, we, tmix, tbelow, smix, sbelow, alpha, beta
       !        endif

       IF(we>0.0_r8) THEN

          hnew=h+we*timestep
          IF(hnew.GT.hmax) THEN
             hnew=hmax
             we=(hmax-h)/timestep
          ENDIF

          ! determine new laymix after mixed layer deepens
          laynew = laymix
          DO WHILE (hnew .GT. zlev(laynew))
             laynew = laynew+1
          ENDDO

          ! determine new tbelow and sbelow using mass-weighted 
          ! mean from hnew to h
          IF(laynew.EQ.laymix) THEN         !within the same layer
             allmass= hnew-h
             allheat= t(laymix)*(hnew-h)
             allsalt= s(laymix)*(hnew-h)
          ELSEIF(laynew.EQ.laymix+1)  THEN  !cross one layer
             allmass= (zlev(laymix)-h) + (hnew-zlev(laymix))
             allheat= t(laymix)*(zlev(laymix)-h) + t(laynew)*(hnew-zlev(laymix))
             allsalt= s(laymix)*(zlev(laymix)-h) + s(laynew)*(hnew-zlev(laymix))
          ELSE                              !corss more than one layer 
             allmass= (zlev(laymix)-h)
             allheat= t(laymix)*(zlev(laymix)-h)
             allsalt= s(laymix)*(zlev(laymix)-h)
             DO k=laymix+1,laynew-1
                allmass=allmass + zlev(k)
                allheat=allheat + t(k)*zlev(k)
                allsalt=allsalt + s(k)*zlev(k)
             ENDDO
             allmass=allmass +(hnew-zlev(laynew-1))
             allheat=allheat + t(laynew)*(hnew-zlev(laynew-1))
             allsalt=allsalt + s(laynew)*(hnew-zlev(laynew-1))
          ENDIF
          IF(allmass /= 0.0_r8)THEN
             tbelow=allheat/allmass
             sbelow=allsalt/allmass
          ELSE
            tbelow=0.0_r8
            sbelow=0.0_r8
            allmass=0.0_r8
          END IF 
          ! compute new tmix and smix
          tmixnew=tmix+timestep*we*(tbelow-tmix)/hnew
          smixnew=smix+timestep*we*(sbelow-smix)/hnew

          !        if(iwrite.eq.1) then
          !         write(99,*)"-------------- ENTRAIN.F -------------------"
          !         write(99,*)"mlo deepens, we=",we 
          !         write(99,*)"laymix, laynew, h, hnew, tmix, tmixnew, smix, smixnew"
          !         write(99,'(2i4, 6f12.5)')laymix, laynew, h, hnew, tmix, tmixnew, smix, smixnew
          !         write(99,*)"tbelow, sbelow"
          !         write(99,'(6f12.5)')tbelow, sbelow
          !        endif
       ENDIF



50     CONTINUE
    ENDIF

    ! ----------------------------------------------------------
    ! ----------------------------------------------------------
    IF (ap.LE.0 .AND. h .GT. dz(1)) THEN     

       ! mixed-layer shoals. hnew satisfies ap(hnew)=0.  use 
       ! 2-point guessing method to solve ap=0.  usually the function
       ! converges after less than 5 iterations for an accuracy
       ! of abs(ap)-->1.e-18 or abs[x(k+1)-x(k)]<1.e-2 . 
       !   x(k+1)= x(k)-[x(k)-x(k-1)]*f[x(k)] / [f(x(k))-f(x(k-1))]
       ! note: 
       ! 1. newton-raphson method can be used too, but needs to 
       !    deduce d(ap)/d(h), which is nontrival either. 
       ! 2. h is set no shallower than dz(1).
       ! 3. ap has a magnitude of  about 1.e-5 to 1.0e-10
       ! 4. sometimes the change in h is unreasonably large within
       !    one timestep. therefore, hnew-h is set no larger than
       !    0.5*dz(laymix), such that h decreases in many timesteps.

       hmin=MAX(dz(1)+0.01_r8,h-0.5_r8*dz(laymix))
       count=0
       xkm1=h               !first guess
       xk  =h-0.5_r8*dz(laymix)    !second guess


       !       call apwe(xkm1, tmix, smix, tbelow, sbelow,  &
       !                 alpha, beta, apkm1, 0, we)
       CALL apwe(capa,grav,rhoref,coriolis,ustar,d1,d2,evap,precip,qsw,qtot,rr,&
            xkm1, tmix, smix, tbelow, sbelow,  &
            alpha, beta, apkm1, 0, we)

       !       call apwe(xk,   tmix, smix, tbelow, sbelow, &
       !                 alpha, beta, apk,   0, we)
       CALL apwe(capa,grav,rhoref,coriolis,ustar,d1,d2,evap,precip,qsw,qtot,rr,&
            xk,   tmix, smix, tbelow, sbelow,  &
            alpha, beta, apk,   0, we)

       DO
          xkp1=xk-(xk-xkm1)*apk/(apk-apkm1) 
          xkp1=MIN(h,MAX(dz(1), xkp1))                   !force dz(1)< hnew <h
          !         call apwe(xkp1,   tmix, smix, tbelow, sbelow, &
          !                 alpha, beta, apkp1,   0, we)
          CALL apwe(capa,grav,rhoref,coriolis,ustar,d1,d2,evap,precip,qsw,qtot,rr,&
               xkp1, tmix, smix, tbelow, sbelow,  &
               alpha, beta, apkp1, 0, we)


          count=count+1

          !          if(iwrite.eq.1) then
          !             write(99,*) "count=",count, " xkm1,xk,xkp1,apkp1",  xkm1,xk,xkp1,apkp1
          !          endif

          IF(ABS(apkp1).LE.1.0e-18_r8 .OR.xkp1.LE.hmin .OR. ABS(xkp1-xk).LT.0.1_r8)THEN
             EXIT
          ELSE 
             xkm1=xk
             xk=xkp1
             apkm1=apk
             apk=apkp1
          END IF
       END DO

       hnew=MAX(hmin,xkp1)

       ! redefine the layer (laymix) in which the mixed layer resides
       laynew = 1
       DO WHILE (hnew .GT. zlev(laynew))
          laynew = laynew+1
       ENDDO

       ! fill the retreated space with new water, assuming
       ! a linear profile. then redistribute the loss of heat 
       ! and gain of salt from the retreated space to the new 
       ! mixed layer.
       IF(laynew.EQ.laymix) THEN
          tmixnew=tmix+MAX(tmix-t(laymix),0.0_r8)*(h-hnew)/hnew
          smixnew=smix+MIN(smix-s(laymix),0.0_r8)*(h-hnew)/hnew
       ELSE
          allheat=0.0_r8
          allsalt=0.0_r8
          DO k=laynew,laymix-1
             tnew(k)=tmix+(zlev(k)-hnew)*(t(laymix)-tmix)/(zlev(laymix)-hnew)
             snew(k)=smix+(zlev(k)-hnew)*(s(laymix)-smix)/(zlev(laymix)-hnew)
             allheat=allheat+MIN(dz(k),zlev(k)-hnew)*(tmix-tnew(k))
             allsalt=allsalt+MIN(dz(k),zlev(k)-hnew)*(smix-snew(k))
          ENDDO
          tmixnew=tmix+allheat/hnew
          smixnew=smix+allsalt/hnew
       ENDIF
       tbelow=tnew(laynew)
       sbelow=snew(laynew)

       !        if(iwrite.eq.1) then
       !         write(99,*)"-------------- ENTRAIN.F -------------------"
       !         write(99,*)"mlo shoals, ap=",ap 
       !         write(99,*)"laymix, laynew, h, hnew, tmix, tmixnew, smix, smixnew"
       !         write(99,'(2i4, 6f12.5)')laymix, laynew, h, hnew, tmix, tmixnew, smix, smixnew
       !         write(99,*)"tbelow, sbelow"
       !         write(99,'(6f12.5)')tbelow, sbelow
       !        endif

    ENDIF
    ! ----------------------------------------------------------
    ! ----------------------------------------------------------

    ! update state variables 
    laymix=laynew
    h=MAX(dz(1),MIN(hmax,hnew))
    tmix=MAX(tmin,MIN(tmax,tmixnew))
    smix=MAX(smin,MIN(smax,smixnew))
    DO k=1,laymix-1
       t(k)=tmix
       s(k)=smix
    ENDDO
    DO k=laymix,lbot
       t(k)=MAX(tmin,MIN(tmax,tnew(k)))
       s(k)=MAX(smin,MIN(smax,snew(k)))
    ENDDO

  END SUBROUTINE entrain

  ! ----------------------------------------------------------------------------
  SUBROUTINE diffusion(h, laymix, km, lbot, &
       visc, dt, dz, zlev, u)
    ! ----------------------------------------------------------------------------
    !  purpose: update temperature or salinity due to vertical diffusion         

    IMPLICIT NONE

    ! input
    REAL(KIND=r8),    INTENT(in)    :: h        !mixed-layer depth
    INTEGER, INTENT(in)    :: laymix   !layer in which mlo resides
    INTEGER, INTENT(in)    :: km       !vertical layers
    INTEGER, INTENT(in)    :: lbot     !ocean bottom layer
    REAL(KIND=r8),    INTENT(in)    :: visc     !viscocity
    REAL(KIND=r8),    INTENT(in)    :: dt       !timestep in seconds 
    REAL(KIND=r8),    INTENT(in)    :: dz(km)   !layer thickness
    REAL(KIND=r8),    INTENT(in)    :: zlev(km) !layer depth
    REAL(KIND=r8),    INTENT(inout) :: u(km)    !temperature or salinity
    ! local variables
    INTEGER :: k,j,kx,kmix
    REAL(KIND=r8)    :: dx,r,aa,bb,cc,count,sum,unew(km)
    !      REAL(KIND=r8), allocatable :: utmp(:), ztmp(:)
    !      REAL(KIND=r8), allocatable :: y(:), f(:)
    !      REAL(KIND=r8), allocatable :: gama(:), alpha(:), beta(:)
    REAL(KIND=r8)    ::  utmp(0:2*km+1) 
    REAL(KIND=r8)    ::  ztmp(0:2*km+1) 
    REAL(KIND=r8)    ::     y(0:2*km) 
    REAL(KIND=r8)    ::     f(0:2*km) 
    REAL(KIND=r8)    ::  gama(0:2*km) 
    REAL(KIND=r8)    :: alpha(0:2*km) 
    REAL(KIND=r8)    ::  beta(0:2*km) 
    utmp(0:2*km+1) =0.0_r8
    ztmp(0:2*km+1) =0.0_r8
    y(0:2*km) =0.0_r8
    f(0:2*km) =0.0_r8
    gama(0:2*km) =0.0_r8
    alpha(0:2*km) =0.0_r8
    beta(0:2*km) =0.0_r8
    unew(1:km)=0.0_r8
    ! ------------------------------------------------------------
    ! ------------------------------------------------------------
    ! METHOD DESCRIPTION
    ! the parabolic first-order partial differential equation
    ! which is of the typical form d(U)/d(t)-c*[d2(U)/d(x2)]=0, 
    ! is sloved by classical implicit finite difference, 
    !  [U(m,n+1)-U(m,n)] = r*[U(m+1,n+1)-2U(m,n+1)+U(m-1,n+1)]
    ! where r=c*dt/(dx*dx), m and n represent grids in space 
    !  and time.  m=1,2,...M-1, n=1,2,....N.
    !
    ! For given 
    !   initial condition   U(m,n), m=0,1,2,....M, and
    !   boundary conditions U(0,n+1), U(M,n+1),
    ! to derive U(m,n+1), a linear tri-diagonal matix must be solved,
    ! | 1+2r -r                        |*| U(1,n+1)  | = | U(1,n)+r*U(0,n+1)  |
    ! |      -r 1+2r  -r               |*| U(2,n+1)  | = | U(2,n)             |
    ! |          -r  1+2r -r           |*| U(3,n+1)  | = | U(3,n)             |
    ! | .............................. |*| ........  | = | ......             |
    ! |                        -r 1+2r |*| U(M-1,n+1)| = | U(M-1,n)+r*U(M,n+1)|
    !
    ! Many standard routines can be used to solve this matrix. here I use the 
    ! inverse-substitution method. The quantity in question is first 
    ! projected to a uniform gird in vertical, and then restored back 
    ! to oroginal non-uniform grid after diffusion calculation.  
    ! ------------------------------------------------------------
    ! ------------------------------------------------------------
    !
    ! search the smallest depth and set as the depth for uniform grid
         dx=dz(1) 
         do k=2,lbot
          if(dz(k).lt.dx) dx=dz(k)
         enddo
    !dx=5.0_r8                             !to save cpu, skip the search
    kx=MIN(INT(zlev(lbot)/dx)-1 ,km-1)              
    IF(kx.LT.4) RETURN                  
    !
    !       allocate( utmp(0:kx+1) )
    !       allocate( ztmp(0:kx+1) )
    !       allocate(    y(0:kx) )
    !       allocate(    f(0:kx) )
    !       allocate( gama(0:kx) )
    !       allocate(alpha(0:kx) )
    !       allocate( beta(0:kx) )
    r=visc*dt/(dx*dx)
    aa=-r
    bb=1.0_r8+2.0_r8*r
    cc=-r

    ! create profile in uniform vertical grid
    DO k=0,2*km+1   
       ztmp(k)=0.5_r8*dx+k*dx           !middle depth of each layer 
    ENDDO
    kmix=MIN(INT(h/dx)-1 ,km-1  )                 
    DO k=0,kmix                    !within the miex-layer
       utmp(k)=u(1)
    ENDDO
    DO k=kmix+1,kx                 !below the miex-layer
       j=MAX(laymix-1,1)
       DO WHILE(ztmp(k).GT.zlev(j))
          j=j+1
       ENDDO
       utmp(k)=u(j)
    ENDDO

    ! set boundary conditions at k=0 and k=kx, right-hand-side matrix,
    ! determine coefficients for left-hand-side matrix for k=[1,kx-1]
    f(1)=utmp(1)+r*utmp(0)
    f(kx-1)=utmp(kx-1)+r*utmp(kx)
    DO k=2,kx-2
       f(k)=utmp(k)
    ENDDO

    alpha(1)=bb
    beta(1)=cc/alpha(1)
    DO k=2,kx-1
       gama(k)=aa
       alpha(k)=bb-gama(k)*beta(k-1)
       beta(k)=cc/alpha(k)
    ENDDO
    !
    ! forward substitute 
    y(1)=f(1)/bb
    DO k=2,kx-1
       y(k)=(f(k)-gama(k)*y(k-1))/alpha(k)
    ENDDO

    ! backward substitute
    utmp(kx-1)=y(kx-1)
    DO k=kx-2,1,-1
       utmp(k)=y(k)-beta(k)*utmp(k+1)
    ENDDO

    !-----------------------------------------------
    ! integrate and restore original non-uniform grid
    sum=0.0_r8
    count=0.0_r8
    DO k=0,kmix            !for mixed layer, derive mean u
       sum=sum+utmp(k)
       count=count+1.0_r8
    ENDDO
    DO j=1,laymix-1
       unew(j)=sum/count
    ENDDO

    k=kmix+1
    sum=0.0_r8
    count=0.0_r8
    DO j=laymix,lbot        !below mixed layer
       DO WHILE(ztmp(k).LT.zlev(j).and.k <2*km)
          sum=sum+utmp(k)
          count=count+1.0_r8
          k=k+1
       ENDDO
       IF(count.GT.0.0_r8) THEN
          unew(j)=sum/count
       ELSE
          unew(j)=u(j)          !for case [zlev(laymix)]-h<<dx, rarely happens
       ENDIF
       count=0.0_r8
       sum=0.0_r8
    ENDDO

    !      if(iwrite.eq.1) then
    !       write(99,*)"----  DIFFUSION.f: change due to diffusion  --------"
    !          write(99,*) "j, zlev(j), u(j), unew(j)-u(j)"
    !         do j=1, lbot
    !          write(99,'(i4, 5e12.4)')j, zlev(j), u(j), unew(j)-u(j)
    !         enddo
    !      endif

    DO j=1,lbot
       u(j)=unew(j)  
    ENDDO

    !       deallocate( utmp )
    !       deallocate( ztmp )
    !       deallocate(    y )
    !       deallocate(    f )
    !       deallocate( gama )
    !       deallocate(alpha )
    !       deallocate( beta )
  END SUBROUTINE diffusion
  ! ---------------------------------------------------------
  SUBROUTINE convect(km,lbot,dz,zlev,h, tmix, smix, laymix, t, s)
    ! ---------------------------------------------------------
    ! purpose:  determines if adjacent layers are convectively
    !           unstable, if they are, then mix them. mixed layer
    !           is included.

    IMPLICIT NONE
    !      include 'mlo.h'
    INTEGER, INTENT(in   ) :: km !, parameter :: km =30            !model vertical layers
    INTEGER, INTENT(in   ) ::  lbot            !level of ocean botom
    REAL(KIND=r8)   , INTENT(in   ) ::  dz(km)             !depth, from top down, m 
    REAL(KIND=r8)   , INTENT(in   )::  zlev(km)           !sum(dz(k)), top down, m

    REAL(KIND=r8)   , INTENT(inout) :: h     
    REAL(KIND=r8)   , INTENT(inout) :: tmix
    REAL(KIND=r8)   , INTENT(inout) :: smix
    INTEGER, INTENT(inout) :: laymix
    REAL(KIND=r8)   , INTENT(inout) :: t(km)
    REAL(KIND=r8)   , INTENT(inout) :: s(km)
    ! local
    !      REAL(KIND=r8), allocatable :: dztmp(:),ttmp(:),stmp(:), rhotmp(:)
    REAL(KIND=r8) :: dztmp (km+2) 
    REAL(KIND=r8) :: ttmp  (km+2) 
    REAL(KIND=r8) :: stmp  (km+2) 
    REAL(KIND=r8) :: rhotmp(km+2) 

    REAL(KIND=r8) :: sumd, sumt, sums           
    INTEGER :: k, ktmp, ltop, llow
    LOGICAL :: flag_conv
 dztmp (1:km+2) =0.0_r8
 ttmp  (1:km+2) =0.0_r8
 stmp  (1:km+2) =0.0_r8
 rhotmp(1:km+2) =0.0_r8
 sumd=0.0_r8; sumt=0.0_r8;sums=0.0_r8
    ! ----------------------------------------------
    !
    !      if(iwrite.eq.1) then
    !        write(99,*) "---------------- CONVECT.F --------------"
    !        write(99,*)" laymix, h, tmix, smix"
    !        write(99,'(i4, 3e12.4)') laymix, h, tmix, smix
    !      endif

    h=MIN(h,zlev(lbot)-0.5_r8*dz(lbot))  
    flag_conv=.FALSE.

    ! create new arrays including the mixed layer
    ktmp=lbot-laymix+2

    dztmp (1:km+2)=0.0_r8 
    ttmp  (1:km+2)=0.0_r8
    stmp  (1:km+2)=0.0_r8
    rhotmp(1:km+2)=0.0_r8

    dztmp(1)=h
    ttmp(1) =tmix
    stmp(1) =smix
    dztmp(2)=MAX(zlev(laymix)-h, 1.e-2_r8)
    ttmp(2) =t(laymix)
    stmp(2) =s(laymix)
    IF(ktmp.GE.3) THEN
       DO k=3,ktmp     
          dztmp(k)=dz(k+laymix-2)
          ttmp(k) =t(k+laymix-2)
          stmp(k) =s(k+laymix-2)
       ENDDO
    ENDIF
    DO k=1,ktmp
       CALL density(ttmp(k), MAX(stmp(k),0.0_r8), rhotmp(k))
    ENDDO

    ! convective adjustment
    DO ltop=1,ktmp-1
       DO llow=ltop+1,ktmp 
          IF (rhotmp(ltop).GT.rhotmp(llow)) THEN
             IF(ltop.EQ.1) flag_conv=.TRUE.
             sumd=0
             sumt=0.0_r8 
             sums=0.0_r8  
             DO k=ltop,llow
                sumt=sumt+ttmp(k)*dztmp(k)*rhotmp(k)
                sums=sums+stmp(k)*dztmp(k)*rhotmp(k)
                sumd=sumd+dztmp(k)*rhotmp(k)                 
             ENDDO
             DO k=ltop,llow
                ttmp(k)=sumt/sumd
                stmp(k)=sums/sumd
                CALL density(ttmp(k), MAX(stmp(k),0.0_r8), rhotmp(k))
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    !
    ! restore arrays
    tmix=ttmp(1)
    smix=stmp(1)
    DO k=1,laymix-1
       t(k)=ttmp(1)
       s(k)=stmp(1)
    ENDDO
    DO k=laymix,lbot
       t(k)=ttmp(k-laymix+2)
       s(k)=stmp(k-laymix+2)
    ENDDO
    !       
    ! if convective adjustment does occur between the mixed-layer 
    ! and layers below, find new laymix and h, which will be at 
    ! the bottom edge of a layer where temperatue and salinity jumps occur.

    IF(flag_conv) THEN
       laymix=1

       !        do k=1,lbot-1
       !         if(t(k).ne.tmix .or. s(k).ne.smix) goto 250
       !          laymix=k+1
       !          h=zlev(k)+1.0e-2_r8
       !        enddo
       ! 250  continue

       DO k=1,lbot-1
          IF(.NOT.(t(k).NE.tmix  .OR. s(k).NE.smix)) THEN
             laymix=k+1
             h=zlev(k)+1.0e-2_r8
          ENDIF
       ENDDO


    ENDIF


    !      if(iwrite.eq.1) then
    !        write(99,*)"flag_conv=", flag_conv, " laymix, h, tmix, smix"
    !        write(99,'(i4, 3f10.2)') laymix, h, tmix, smix
    !      endif

    !     deallocate( dztmp )
    !      deallocate( ttmp )
    !      deallocate( stmp )
    !      deallocate( rhotmp )

  END SUBROUTINE convect


  ! ------------------------------------------------------
  SUBROUTINE surflux(km,lbot,capa,d1,d2,evap,precip,qcor,qsw,qtot,rhoref,&
       rr,scor,timestep,zlev,dz,h, tmix, smix, laymix, t, s)
    ! ------------------------------------------------------
    !  purpose: update temperature and salinity due to 
    !    surface heat and water flux, including flux corrections

    IMPLICIT NONE
    !      include 'mlo.h'

    ! input
    INTEGER, INTENT(IN   ) :: km            ! km =30 model vertical layers
    INTEGER , INTENT(IN   ) ::  lbot            !level of ocean botom
    REAL(KIND=r8)    , INTENT(IN   ) ::capa! , parameter :: capa =3950.0         !heat capacity of sea water 
    REAL(KIND=r8)    , INTENT(IN   ) :: d1                  !water optical coefficient
    REAL(KIND=r8)    , INTENT(IN   ):: d2                  !water optical coefficient
    REAL(KIND=r8)    , INTENT(IN   ):: evap                !evaporation rate, kg/m^2/s
    REAL(KIND=r8)    , INTENT(IN   ):: precip              !precipitation rate, kg/m^2/s
    REAL(KIND=r8)    , INTENT(IN   ):: qcor                !heat flux correction, w/m^2  
    REAL(KIND=r8)    , INTENT(IN   ):: qsw                 !surface net downward solar flux, w/m^2
    REAL(KIND=r8)    , INTENT(IN   ):: qtot                !surface total downward heat flux, w/m^2
    REAL(KIND=r8)    , INTENT(IN   ):: rhoref!, parameter :: rhoref = 1024.438    !sea water reference density, kg/m^3
    REAL(KIND=r8)    , INTENT(IN   ):: rr                  !water optical coefficient
    REAL(KIND=r8)    , INTENT(IN   ):: scor                !water flux correction, 1/1000  
    REAL(KIND=r8)    , INTENT(IN   ) :: timestep!,    parameter :: timestep=86400.0  !timestep, seconds
    REAL(KIND=r8)    , INTENT(IN   ) ::  zlev(km)           !sum(dz(k)), top down, m
    REAL(KIND=r8)    , INTENT(IN   )::  dz(km)             !depth, from top down, m 

    REAL(KIND=r8),    INTENT(in)    :: h        !mixed-layer depth, m 
    INTEGER, INTENT(in)    :: laymix   !layer in which mixed layer resides
    REAL(KIND=r8),    INTENT(inout) :: tmix     !mixed-layer temperature, K
    REAL(KIND=r8),    INTENT(inout) :: smix     !mixed-layer salinity, 1/1000
    REAL(KIND=r8),    INTENT(inout) :: t(km)    !temperature profile, K
    REAL(KIND=r8),    INTENT(inout) :: s(km)    !salinity profile, 1/1000
    ! local variables
    REAL(KIND=r8) :: qswpen, qswin, qswout, exph1, exph2, dist, dtk(km)
    REAL(KIND=r8) :: tnew(km), snew(km), tmixnew, smixnew
    INTEGER :: k
    qswpen=0.0_r8; qswin=0.0_r8; qswout=0.0_r8; exph1=0.0_r8; exph2=0.0_r8; dist=0.0_r8; dtk(1:km)=0.0_r8
    tnew(1:km)=0.0_r8; snew(1:km)=0.0_r8; tmixnew=0.0_r8; smixnew=0.0_r8
    ! ----------------------------------------------------------
    ! determine solar radiation that penetrates through the mixed layer
    exph1 =EXP(MAX(-h/d1,-30.0_r8))
    exph2 =EXP(MAX(-h/d2,-30.0_r8))
    qswpen=qsw*(rr*exph1+(1.0_r8-rr)*exph2)

    ! mixed-layer temperature and salinity 
    ! changes due to surface fluxes
    tmixnew=tmix + timestep*(qtot-qswpen+qcor)/(rhoref*capa*h)
    smixnew=smix + timestep*smix*(evap-precip+scor)/(rhoref*h)

    !      if(iwrite.eq.1) then
    !       write(99,*)"---------------- SURFLUX.F ----------------"
    !       write(99,*)"Changes by Surface FLuxes: h, tmix, Dtmix, smix, Dsmix"
    !       write(99,'(6f12.5)')h, tmix, tmixnew-tmix, smix, smixnew-smix
    !      endif

    tmix=tmixnew
    smix=smixnew
    DO  k=1,laymix-1    
       t(k)=tmix            
       s(k)=smix            
    ENDDO

    ! ------------------------------------------------------------
    ! temperature changes in layers below the mixed-layer 
    ! due to the penetration of solar radiation
    DO k=1,km
       dtk(k)=0.0_r8
    ENDDO

    qswin=qswpen
    DO  k=laymix,lbot
       exph1 =EXP(MAX(-zlev(k)/d1,-30.0_r8))
       exph2 =EXP(MAX(-zlev(k)/d2,-30.0_r8))
       qswout=qsw*(rr*exph1+(1.0_r8-rr)*exph2)
       dist=dz(k)
       IF(k.EQ.laymix) dist=zlev(laymix)-h

       !       if(dist.le.1.e-2) goto 10
       !          dtk(k)=timestep*(qswin-qswout)/(rhoref*capa*dist)
       !          t(k)=t(k) + dtk(k)
       !          qswin=qswout
       ! 10   continue
       IF(dist > 1.e-2) THEN 
          dtk(k)=timestep*(qswin-qswout)/(rhoref*capa*dist)
          t(k)=t(k) + dtk(k)
          qswin=qswout
       END IF

    ENDDO

    !      if(iwrite.eq.1) then
    !       write(99,*)"---- t(k) change by solar penetration --------"
    !         write(99,*)"k, t(k), dtk(k)"
    !         do k=laymix, lbot
    !          write(99,'(i4, 2e13.5)')k, t(k)-dtk(k), dtk(k)
    !         enddo
    !      endif

  END SUBROUTINE surflux

  ! -------------------------------------------------------------
  SUBROUTINE apwe(capa,grav,rhoref,coriolis,ustar,d1,d2,evap,precip,qsw,qtot,rr,&
       h, tmix, smix, tbelow, sbelow,  &
       alpha, beta, ap, flag, we)
    ! -------------------------------------------------------------
    !  purpose: compute ap, and we if flag=1
    !  following eq.(50) of gaspar, jpo, 1988.

    IMPLICIT NONE
    !      include 'mlo.h'
    REAL(KIND=r8), INTENT(in) :: capa! =3950.0         !heat capacity of sea water 
    REAL(KIND=r8), INTENT(in) :: grav! =9.81           !gravity, kg/m/s^2
    REAL(KIND=r8), INTENT(in) :: rhoref != 1024.438    !sea water reference density, kg/m^3
    REAL(KIND=r8),    INTENT(in)  :: coriolis            !coriolis parameter f  
    REAL(KIND=r8),    INTENT(in)  :: ustar               !frictional velocity, m/s
    REAL(KIND=r8),    INTENT(in) :: d1                  !water optical coefficient
    REAL(KIND=r8),    INTENT(in) :: d2                  !water optical coefficient
    REAL(KIND=r8), INTENT(in) :: evap                !evaporation rate, kg/m^2/s
    REAL(KIND=r8), INTENT(in) :: precip              !precipitation rate, kg/m^2/s
    REAL(KIND=r8), INTENT(in)  :: qsw                 !surface net downward solar flux, w/m^2
    REAL(KIND=r8) , INTENT(in) :: qtot                !surface total downward heat flux, w/m^2
    REAL(KIND=r8), INTENT(in) :: rr                  !water optical coefficient

    REAL(KIND=r8), PARAMETER :: m1=0.45_r8, m2=2.6_r8, m3=1.9_r8, m4=2.3_r8, m5=0.6_r8
    REAL(KIND=r8), PARAMETER :: a1=0.6_r8,  a2=0.3_r8

    ! input
    REAL(KIND=r8),    INTENT(in)  :: h, tmix, smix, tbelow, sbelow
    REAL(KIND=r8),    INTENT(in)  :: alpha, beta
    INTEGER, INTENT(in) :: flag
    ! output
    REAL(KIND=r8), INTENT(out) :: ap
    REAL(KIND=r8), INTENT(out) :: we
    ! local variables
    REAL(KIND=r8) :: bh, lamada,length,hl34,hlp48,sp,cp1,cp3,c4
    REAL(KIND=r8) :: tmp1,tmp2,tmp3,tmp4, temp, deltab
    REAL(KIND=r8) :: rhomix, rhobelow
    ap=0.0_r8
    we=0.0_r8
    bh=0.0_r8; lamada=0.0_r8;length=0.0_r8;hl34=0.0_r8;hlp48=0.0_r8;sp=0.0_r8;cp1=0.0_r8;cp3=0.0_r8;c4=0.0_r8
    tmp1=0.0_r8;tmp2=0.0_r8;tmp3=0.0_r8;tmp4=0.0_r8; temp=0.0_r8; deltab=0.0_r8
    rhomix=0.0_r8; rhobelow=0.0_r8
    ! ----------------------------------------------------------
    !      if(iwrite.eq.1) then
    !       write(99,*) "-------------- APWE.F ------------------------- "
    !       write(99,*) "input: h tmix, smix, tbelow, sbelow, alpha, beta"
    !       write(99,'(5f12.6,2e12.4)') h, tmix, smix, tbelow, sbelow, alpha, beta
    !      endif

    !      call bh14(                                               bh, h, tmix, smix, alpha, beta) !bouyancy
    CALL bh14(grav,capa,rhoref,d1,d2,evap,precip,qsw,qtot,rr,bh, h, tmix, smix, alpha, beta)

    lamada=ABS(coriolis)/ustar             !inverse of ekman lengthh scale
    length=bh/(ustar*ustar*ustar)          !inverse of bulk monin-obukhov length
    temp=MAX(-30.0_r8,MIN(30.0_r8,h*length))       !provent from overflowing
    hl34=a1+a2*MAX(1.0_r8,2.5_r8*lamada*h)*EXP(temp)   !eq34
    hlp48=a1+a2*EXP(temp)                  !eq48

    cp1=( (2.0_r8-2.0_r8*m5)*(hl34/hlp48) + m4 )/6.0_r8
    cp3=( m4*(m2+m3) - (hl34/hlp48)*(m2+m3-m5*m3) )/3.0_r8
    cp3=MAX(cp3,0.1_r8)
    ap =cp3*ustar*ustar*ustar - cp1*h*bh

    ! compute we
    IF(ap.GT.0.0_r8 .AND. flag.EQ.1) THEN

       c4 =2.0_r8*m4/(m1*m1)  
       sp =(m2+m3)*ustar*ustar*ustar - 0.5_r8*h*bh
       deltab=alpha*grav*(tmix-tbelow)-beta*grav*(smix-sbelow)

       !!  add the parametrization of kim to the denominator to keep the
       !!  mixed layer from over deepening
       deltab = deltab + 9.0_r8*MAX(1.0e-4_r8,ustar**2)

       tmp1=0.5_r8*ap+cp1*sp
       tmp2=0.5_r8*ap-cp1*sp
       tmp3=-tmp1+SQRT(tmp2*tmp2 + 2.0_r8*c4*hl34*hl34*ap*sp) 
       tmp4=c4*hl34*hl34 - cp1
       we=tmp3/(tmp4*h*deltab)

    ENDIF

    !      if(iwrite.eq.1) then
    !       write(99,*) "-------------- APWE.F ------------------------- "
    !       write(99,*) "output: ap, we, tmp1, tmp2, tmp3, tmp4, deltab"
    !       write(99,'(7e12.4)') ap, we, tmp1, tmp2, tmp3, tmp4, deltab
    !      endif

  END SUBROUTINE apwe


  !------------------------------------------------------
  SUBROUTINE bh14(grav,capa,rhoref,d1,d2,evap,precip,qsw,qtot,rr,bh, h, tmix, smix, alpha, beta)
    !------------------------------------------------------
    ! purpose: coumpte term B(h), eq.(14) of gaspar, jpo, 1988.

    IMPLICIT NONE
    !      include 'mlo.h'
    ! input
    REAL(KIND=r8), INTENT(in) :: grav! =9.81           !gravity, kg/m/s^2
    REAL(KIND=r8), INTENT(in) :: capa !=3950.0         !heat capacity of sea water 
    REAL(KIND=r8), INTENT(in) :: rhoref! = 1024.438    !sea water reference density, kg/m^3

    REAL(KIND=r8), INTENT(in)  :: d1                  !water optical coefficient
    REAL(KIND=r8), INTENT(in)  :: d2                  !water optical coefficient
    REAL(KIND=r8), INTENT(in) :: evap                !evaporation rate, kg/m^2/s
    REAL(KIND=r8), INTENT(in) :: precip              !precipitation rate, kg/m^2/s
    REAL(KIND=r8), INTENT(in)  :: qsw                 !surface net downward solar flux, w/m^2
    REAL(KIND=r8) , INTENT(in) :: qtot                !surface total downward heat flux, w/m^2
    REAL(KIND=r8), INTENT(in) :: rr                  !water optical coefficient

    REAL(KIND=r8), INTENT(in) :: h      !mixed-layer depth, meter
    REAL(KIND=r8), INTENT(in) :: tmix   !mixed-layer temperature, K
    REAL(KIND=r8), INTENT(in) :: smix   !mixed-layer salinity, 1/1000
    REAL(KIND=r8), INTENT(in) :: alpha  !thermal expansion coefficient
    REAL(KIND=r8), INTENT(in) :: beta   !saline contration coefficient
    ! outpu 
    REAL(KIND=r8), INTENT(out) :: bh    !bouyancy 
    ! temporary
    REAL(KIND=r8) :: x1,x2,x3,x4,temp,sum,exph1,exph2
    bh=0.0_r8
    x1=0.0_r8
    x2=0.0_r8
    x3=0.0_r8
    x4=0.0_r8
    temp=0.0_r8
    sum=0.0_r8
    exph1=0.0_r8
    exph2=0.0_r8
    !-------------------------------------------------------
    !      if(iwrite.eq.1) then
    !       write(99,*) "------------- BH14.F ----------------"
    !       write(99,*) "input: h, tmix, smix, alpha, beta"
    !       write(99,'(3f10.2,2e12.4)')  h, tmix, smix, alpha, beta
    !      endif

    ! determine solar penetration coefficient
    exph1=EXP(MAX(-h/d1,-30.0_r8))
    exph2=EXP(MAX(-h/d2,-30.0_r8))

    temp=alpha*grav/(rhoref*capa)

    x1=beta*grav*smix*(precip-evap)/rhoref
    x2=temp*qtot
    x3=temp*qsw*(rr*exph1+(1.0_r8-rr)*exph2)         

    ! vertical integral from the surface to depth (-h) 
    sum=rr*d1*(1.0_r8-exph1)+(1.0_r8-rr)*d2*(1.0_r8-exph2)
    x4=-2.0_r8*temp*qsw*sum/h

    bh=x1+x2+x3+x4

    !      if(iwrite.eq.1) then
    !       write(99,*) "output: bh, x1, x2, x3, x4"         
    !       write(99,'(7e12.4)') bh, x1, x2, x3, x4         
    !      endif

  END SUBROUTINE bh14


  ! ----------------------------------------
  SUBROUTINE density(t, s, rho)
    ! ----------------------------------------
    IMPLICIT NONE

    ! input
    REAL(KIND=r8), INTENT(in)  :: t     !unit, K
    REAL(KIND=r8), INTENT(in)  :: s     !unit, 1/1000
    ! output
    REAL(KIND=r8), INTENT(out) :: rho   !unit, kg/m^3 
    ! local
    REAL(KIND=r8) :: tc

    ! compute density using the international equation 
    ! of state of sea water 1980, (pond and pickard, 
    ! introduction to dynamical oceanography, pp310). 
    ! compression effects are not included

    rho = 0.0_r8
    tc = t - 273.15_r8

    !  effect of temperature on density (lines 1-3)
    !  effect of temperature and salinity on density (lines 4-8)
    rho =                                                          &
         999.842594_r8                 +  6.793952e-2_r8 * tc          &
         - 9.095290e-3_r8 * tc**2        +  1.001685e-4_r8 * tc**3          &
         - 1.120083e-6_r8 * tc**4        +  6.536332e-9_r8 * tc**5          &
         + 8.24493e-1_r8 * s          -  4.0899e-3_r8 * tc * s          &
         + 7.6438e-5_r8 * tc**2 * s   -  8.2467e-7_r8 * tc**3 * s          &
         + 5.3875e-9_r8 * tc**4 * s   -  5.72466e-3_r8 * s**1.5_r8          &
         + 1.0227e-4_r8 * tc * s**1.5_r8 -  1.6546e-6_r8 * tc**2 * s**1.5_r8          &
         + 4.8314e-4_r8 * s**2

  END SUBROUTINE density
  
  ! ----------------------------------------------
  SUBROUTINE rhocoef(t, s, rhoref, alpha, beta)
    ! ----------------------------------------------

    !  compute thermal expansion coefficient (alpha) 
    !  and saline contraction coefficient (beta) using 
    !  the international equation of state of sea water 
    !  (1980). Ref: pond and pickard, introduction to 
    !  dynamical oceanography, pp310.  
    !  note: compression effects are not included

    IMPLICIT NONE
    REAL(KIND=r8), INTENT(in)  :: t, s, rhoref 
    REAL(KIND=r8), INTENT(out) :: alpha, beta  
    REAL(KIND=r8) :: tc

    tc = t - 273.15_r8

    alpha =                                                         &
         &     6.793952e-2_r8                                                   &
         &   - 2.0_r8 * 9.095290e-3_r8 * tc     +  3.0_r8 * 1.001685e-4_r8 * tc**2       &
         &   - 4.0_r8 * 1.120083e-6_r8 * tc**3  +  5.0_r8 * 6.536332e-9_r8 * tc**4       &
         &   - 4.0899e-3_r8 * s                                                 &
         &   + 2.0_r8 * 7.6438e-5_r8 * tc * s  -  3.0_r8 * 8.2467e-7_r8 * tc**2 * s      &
         &   + 4.0_r8 * 5.3875e-9_r8 * tc**3 * s                                   &
         &   + 1.0227e-4_r8 * s**1.5_r8 -  2.0_r8 * 1.6546e-6_r8 * tc * s**1.5_r8 

    !
    alpha =  -alpha/rhoref

    beta  =                                             &
         &   8.24493e-1_r8          -  4.0899e-3_r8 * tc               &
         &   + 7.6438e-5_r8 * tc**2 -  8.2467e-7_r8 * tc**3            &
         &   + 5.3875e-9_r8 * tc**4 -  1.5_r8 * 5.72466e-3_r8 * s**0.5_r8     &
         &   + 1.5_r8 * 1.0227e-4_r8 * tc * s**0.5_r8                      &
         &   -  1.5_r8 * 1.6546e-6_r8 * tc**2 * s**0.5_r8                  &
         &   + 2.0_r8 * 4.8314e-4_r8 * s

    beta = beta / rhoref

  END SUBROUTINE rhocoef


  ! ----------------------------------------------
  ! ----------------------------------------------
  SUBROUTINE relax (km         ,&
                    lbot       ,&
                    timestep   ,&
                    tclim      ,&
                    sclim      ,&
                    laymix     ,&
                    tmix       ,&
                    smix       ,&
                    t          ,&
                    s           )
    ! ----------------------------------------------

    !  relax temperature and salinity below mixed layer  to 
    !  observed climatologies with a time scale of 10 years


    IMPLICIT NONE

    !      include 'mlo.h'

    INTEGER         , INTENT(IN   ) :: km            ! km =30 model vertical layers
    INTEGER         , INTENT(IN   ) :: lbot            !level of ocean botom
    REAL(KIND=r8)   , INTENT(IN   ) :: timestep! ,    parameter :: timestep=86400.0  !timestep, seconds
    REAL(KIND=r8)   , INTENT(IN   ) :: tclim(km)!           !temperature profile for relaxation, K
    REAL(KIND=r8)   , INTENT(IN   ) :: sclim(km)!           !salinity profile for relaxation, K
    INTEGER         , INTENT(in   ) :: laymix
    REAL(KIND=r8)   , INTENT(inout) :: tmix 
    REAL(KIND=r8)   , INTENT(inout) :: smix 
    REAL(KIND=r8)   , INTENT(inout) :: t(km)
    REAL(KIND=r8)   , INTENT(inout) :: s(km)

    REAL(KIND=r8)   , PARAMETER :: rate=3.170979198376e-9_r8       !1/(10*365*86400)

    INTEGER :: k

    DO k=1, laymix-1 
       !       t(k) = t(k) + rate* (tclim(1) - t(k)) * timestep
       s(k) = s(k) + rate* (sclim(1) - s(k)) * timestep
    END DO
    !     tmix = t(1)
    smix = s(1)

    DO k=laymix, lbot 
       t(k) = t(k) + rate* (tclim(k) - t(k)) * timestep
       s(k) = s(k) + rate* (sclim(k) - s(k)) * timestep
    END DO

  END SUBROUTINE relax


  FUNCTION GetOceanAlb(ipt,tau_in,sza_in,wind_in,chl_in,wls_in,wle_in)

    IMPLICIT NONE
    INTEGER      , INTENT(IN   ) :: ipt
    REAL(KIND=r8), INTENT(IN   ) :: tau_in 
    REAL(KIND=r8), INTENT(IN   ) :: sza_in 
    REAL(KIND=r8), INTENT(IN   ) :: wind_in
    REAL(KIND=r8), INTENT(IN   ) :: chl_in 
    REAL(KIND=r8), INTENT(IN   ) :: wls_in 
    REAL(KIND=r8), INTENT(IN   ) :: wle_in 
    REAL(KIND=r8) :: GetOceanAlb
    REAL(KIND=r8) :: albb  (nb)
    REAL(KIND=r8) :: wtflx (nb)
    REAL(KIND=r8) :: rt(2)
    REAL(KIND=r8) :: rs(2)
    REAL(KIND=r8) :: rw(2)
    REAL(KIND=r8) :: rc(2)
    REAL(KIND=r8) :: tau
    REAL(KIND=r8) :: sza
    REAL(KIND=r8) :: wind
    REAL(KIND=r8) :: chl
    REAL(KIND=r8) :: wls 
    REAL(KIND=r8) :: wle 
    REAL(KIND=r8) :: twtb
    REAL(KIND=r8) :: wtb
    REAL(KIND=r8) :: wtf
    REAL(KIND=r8) :: wtt
    REAL(KIND=r8) :: albedo
    INTEGER :: ib
    INTEGER :: ib1
    INTEGER :: ib2
    INTEGER :: icc
    INTEGER :: ic
    INTEGER :: irec
    INTEGER :: irr
    INTEGER :: iss
    INTEGER :: is
    INTEGER :: itt
    INTEGER :: it
    INTEGER :: iww
    INTEGER :: iw
    !     --------------------------- INPUT ---------------------------------------
    !   
    !    specify the parameters for albedo here:

    !        tau = 0.40              !aerosol/cloud optical depth
    !        sza = 0.50              !cosine of solar zenith angle
    !        wind = 10.00            !wind speed in m/s
    !        chl = 0.10              !chlorophyll concentration in mg/m3
    !        wls = 0.69              !start wavelength (um) of your band
    !        wle = 1.19              !end wavelength (um) of your band
    !
    !   Output albedo for these input is 0.068. If not, your system may use
    !   different record unit and you need change the record length.
    !     -------------------------------------------------------------------------

    !    now find the albedo corresponding to the 4 parameters above:
    !-----------------------------------------------------------------------
    ! Computes surface albedos over ocean for Slab Ocean Model (SOM)
    !
    ! Two spectral surface albedos for direct (dir) and diffuse (dif)
    ! incident radiation are calculated. The spectral intervals are:
    !   s (shortwave)  = 0.2-0.7 micro-meters
    !   l (longwave)   = 0.7-5.0 micro-meters
    !

    !real asdir(plond)     ! Srf alb for direct  rad   0.2-0.7 micro-ms (0.20-0.69 )
    !real asdif(plond)     ! Srf alb for diffuse rad   0.2-0.7 micro-ms (1.19-2.38 )

    !real aldir(plond)     ! Srf alb for direct rad   0.7-5.0 micro-ms (0.69-1.19 )
    !real aldif(plond)     ! Srf alb for diffuse rad  0.7-5.0 micro-ms (2.38-4.0um)
    GetOceanAlb=0.8e0_r8
    albb  (1:nb)=0.0_r8
    wtflx (1:nb)=0.0_r8
    rt(1:2)=0.0_r8
    rs(1:2)=0.0_r8
    rw(1:2)=0.0_r8
    rc(1:2)=0.0_r8
    tau=0.0_r8
    sza=0.0_r8
    wind=0.0_r8
    chl=0.0_r8
    wls =0.0_r8
    wle =0.0_r8
    twtb=0.0_r8
    wtb=0.0_r8
    wtf=0.0_r8
    wtt=0.0_r8
    albedo=0.0_r8
    ib=0
    ib1=0
    ib2=0
    icc=0
    ic=0
    irec=0
    irr=0
    iss=0
    is=0
    itt=0
    it=0
    iww=0
    iw=0
    
    tau = tau_in 
    sza = sza_in 
    wind= wind_in
    chl = chl_in 
    wls = wls_in 
    wle = wle_in 
    IF(tau.LT.0.0_r8 .OR. (sza.LT.0.0_r8 .OR. sza.GT.1.0_r8) .OR. wind.LT.0.0_r8 &
         .OR. chl.LT.0.0_r8)STOP 'Err: input parameters wrong!'
    IF(tau  .GT. taunode(nt))tau=taunode(nt)
    IF(wind .GT. windnode(nw))wind=windnode(nw)
    IF(chl  .GT. 45.0_r8)chl=45.0_r8
    CALL locate(taunode,nt,tau,it)
    CALL locate(szanode,ns,sza,is)
    CALL locate(windnode,nw,wind,iw)
    CALL locate(chlnode,nc,chl,ic)
    rt(2) = (tau-taunode(it))/(taunode(it+1)-taunode(it))
    rt(1) = 1.0_r8-rt(2)
    rs(2) = (sza-szanode(is))/(szanode(is+1)-szanode(is))
    rs(1) = 1.0_r8-rs(2)
    rw(2) = (wind-windnode(iw))/(windnode(iw+1)-windnode(iw))
    rw(1) = 1.0_r8-rw(2)
    rc(2) = (chl-chlnode(ic))/(chlnode(ic+1)-chlnode(ic))
    rc(1) = 1.0_r8-rc(2)

    IF(wls .GT. wle)STOP 'Err: Start wavelength should be smaller.'
    IF(wls .LT. wlnode(1))wls=wlnode(1)
    IF(wle .GT. wlnode(25))wle=wlnode(25)
    CALL locate(wlnode,25,wls,ib1)
    CALL locate(wlnode,25,wle,ib2)

    !             ** get alb(ib1:ib2) by 4 dimensional linear interpolaton **
    DO ib=ib1,ib2
       albb(ib) = 0.0_r8
       wtflx(ib) = 0.0_r8
    ENDDO
    DO  itt=it,it+1
       DO  iss=is,is+1
          DO  iww=iw,iw+1
             DO  icc=ic,ic+1
                irec = (itt-1)*15*7*5 + (iss-1)*7*5 + (iww-1)*5 + icc
                !read(1,rec=irec) alb
                wtt = rt(itt-it+1)*rs(iss-is+1)*rw(iww-iw+1)*rc(icc-ic+1)
                DO ib=ib1,ib2
                   albb(ib) = albb(ib) + wtt*alb(irec,ib)
                ENDDO
             END DO
          END DO

          !                  *** get 24 band down flux weights ***
          irr = 8400 + (itt-1)*15 + iss
          !read(1,rec=irr) rflx

          wtf = rt(itt-it+1)*rs(iss-is+1)
          !         do ib=ib1,ib2
          DO ib=1,24
             wtflx(ib) = wtflx(ib) + wtf*rflx(irr,ib)
          ENDDO
       END DO
    END DO

    !             ** get albedo in the specified band by weighted sum **
    twtb = 0.0_r8
    albedo = 0.0_r8
    DO ib=ib1,ib2
       IF(ib .EQ. ib1)THEN
          wtb = (wlnode(ib1+1)-wls)/(wlnode(ib1+1)-wlnode(ib1))
       ELSE IF (ib .EQ. ib2)THEN
          wtb = (wle-wlnode(ib2))/(wlnode(ib2+1)-wlnode(ib2))
       ELSE
          wtb = 1.0_r8
       ENDIF
       albedo = albedo + wtb*wtflx(ib)*albb(ib)
       twtb = twtb + wtb*wtflx(ib)
    ENDDO
    IF (twtb.EQ.0.0_r8 .AND. ib1.EQ.ib2)THEN
       albedo = albb(ib1)
    ELSE
       albedo = albedo/twtb
    ENDIF
    GetOceanAlb=albedo
    !WRITE(*,'(6a6,/,6f6.2,/,a11,f6.3,/)')'Tau','COSUN','Wind', &
    !     'Chl','WL1','WL2',tau,sza,wind,chl,wls,wle,&
    !     'Albedo =',albedo

  END FUNCTION GetOceanAlb

  !=======================================================================
  SUBROUTINE locate(xx,n,x,j)
    IMPLICIT NONE
    !
    ! purpose:  given an array xx of length n, and given a value X, returns
    !           a value J such that X is between xx(j) and xx(j+1). xx must
    !           be monotonic, either increasing of decreasing. this function
    !           returns j=1 or j=n-1 if x is out of range.
    !
    ! input:
    !   xx      monitonic table
    !   n       size of xx
    !   x       single floating point value perhaps within the range of xx
    !
    INTEGER, INTENT(IN   ) :: n
    REAL(KIND=r8), INTENT(IN   ) :: x,xx(n)
    INTEGER, INTENT(OUT  ) :: j

    INTEGER :: jl,jm,ju
    j=0;jl=0;jm=0;ju=0
    
    IF(x.EQ.xx(1)) THEN
       j=1
       RETURN
    ENDIF
    IF(x.EQ.xx(n)) THEN
       j=n-1
       RETURN
    ENDIF
    jl=1
    ju=n
10  IF(ju-jl.GT.1) THEN
       jm=(ju+jl)/2
       IF((xx(n).GT.xx(1)).EQV.(x.GT.xx(jm)))THEN
          jl=jm
       ELSE
          ju=jm
       ENDIF
       GOTO 10
    ENDIF
    j=jl
    RETURN
  END SUBROUTINE locate
  !=======================================================================

END MODULE SlabOceanModel
