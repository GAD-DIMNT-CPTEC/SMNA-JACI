!  $Author: pkubota e Solange$
!  $Date: 2011/07/28 14:34:46 $
!  $Revision: 2.0 $
!
MODULE GridDump
  USE Constants, ONLY: &
       rk,  &
       i8,  &
       r8,  &
       r4,  &
       i4

  USE FieldsPhysics, ONLY:     &
      sfc

  USE Parallelism, ONLY: &
       myid,             &
       maxNodes,         &
       MsgDump,          &
       MsgOne,           &
       FatalError
       
  USE Communications, ONLY: &
       Collect_Grid_Red,    &
       Collect_Grid_Sur,    &
       Collect_Spec,        &
       Collect_Grid_Full,   &
       Collect_Grid_Sur_Print2, &   ! add solange 22-07-2011. Adpatada para escrita dos dados em formato Grade para o letkf.
       p2d

  USE Utils, ONLY:     &
       IBJBtoIJ,       &
       IJtoIBJB,       &
       SplineIBJBtoIJ, &
       LinearIBJBtoIJ, &
       SeaMaskIBJBtoIJ,&
       NearestIBJBtoIJ,&
       NearestIJtoIBJB,&
       vfirec

 IMPLICIT NONE
  SAVE       
  PRIVATE
  INTEGER,PARAMETER :: OUTSWAP  = 96    !102
  LOGICAL           :: reducedGrid
  INTEGER, ALLOCATABLE :: ibMaxPerJB(:)
  INTEGER           :: iMaxNew
  INTEGER           :: jMaxNew
  INTEGER           :: kMaxNew
  INTEGER           :: ibMax
  INTEGER           :: jbMax
  INTEGER           :: ijMaxGauQua
  INTEGER(KIND=i4)  :: version=1
  REAL(KIND=r8), ALLOCATABLE, DIMENSION(:,:) :: soil_type ! FAO/USDA soil texture

  PUBLIC :: InitGridDump
  PUBLIC :: write_GridSigma_file
!  INTERFACE WriteField2
!     MODULE PROCEDURE WriteField82Dz,WriteField81D, WriteField82D, WriteField82D_i8, WriteField81D_i8
!  END INTERFACE
 
CONTAINS
  SUBROUTINE InitGridDump(si,sl,idatec,reducedGrid_in,iMax_in,&
                          jMax_in,kMax_in,ibMax_in,jbMax_in,&
                          record_type,nfsoiltp,fNameSoilType,nfprt,ibMaxPerJB_in)


    IMPLICIT NONE
    !
    REAL(r8) ,TARGET, INTENT(INOUT) :: si(:)
    REAL(r8) ,TARGET, INTENT(INOUT) :: sl(:)    
    INTEGER, INTENT(IN   ) :: idatec(4)
    LOGICAL, INTENT(IN   ) :: reducedGrid_in
    INTEGER, INTENT(IN   ) :: iMax_in
    INTEGER, INTENT(IN   ) :: jMax_in
    INTEGER, INTENT(IN   ) :: kMax_in
    INTEGER, INTENT(IN   ) :: ibMax_in
    INTEGER, INTENT(IN   ) :: jbMax_in
    CHARACTER(LEN=*), INTENT(IN   ):: record_type
    INTEGER, INTENT(IN   ) :: nfsoiltp
    CHARACTER(LEN=*), INTENT(IN   ):: fNameSoilType
    INTEGER, INTENT(IN   ) :: nfprt
    INTEGER, INTENT(IN)    :: ibMaxPerJB_in(jbMax_in)

    REAL(KIND=r8), ALLOCATABLE :: SoilType(:,:)    ! soil texture
    INTEGER , ALLOCATABLE      :: iSoilType(:,:)   ! soil texture

    INTEGER :: ierr
    INTEGER :: IOL


    reducedGrid = reducedGrid_in
    iMaxNew=  iMax_in
    jMaxNew=  jMax_in
    kMaxNew=  kMax_in
    ibMax  =  ibMax_in
    jbMax  =  jbMax_in
    ijMaxGauQua =iMax_in*jMax_in
    ALLOCATE (ibMaxPerJB(jbMax))
    ibMaxPerJB=ibMaxPerJB_in

    ALLOCATE(SoilType(iMaxNew,jMaxNew))
    
    ALLOCATE(soil_type(ibMax,jbMax))

    !write(6,*)'In InitGridDump ',fNameSoilType,record_type(1:3),reducedGrid




    IF (record_type == 'dir' ) THEN      !sequential mode

       INQUIRE (IOLENGTH=IOL) SoilType
       OPEN(UNIT=nfsoiltp,FILE=TRIM(fNameSoilType),&
            FORM='UNFORMATTED', ACCESS='DIRECT', RECL=IOL, ACTION='READ', &
            STATUS='OLD', IOSTAT=ierr)
       IF (ierr /= 0) THEN
          WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
               TRIM(fNameSoilType), ierr
          STOP "**(ERROR)**"
       END IF
       ALLOCATE(ISoilType(iMaxNew,jMaxNew))
       READ(UNIT=nfsoiltp, REC=1) ISoilType

       SoilType=REAL(iSoilType,KIND(SoilType))
       DEALLOCATE(ISoilType)

       IF (reducedGrid) THEN
          CALL NearestIJtoIBJB(SoilType,soil_type)
       ELSE
          CALL IJtoIBJB( SoilType,soil_type)
       END IF

       CLOSE(UNIT=nfsoiltp)

    ELSE IF (record_type == 'vfm') THEN !vformat model

       OPEN(UNIT=nfsoiltp,FILE=TRIM(fNameSoilType),FORM='formatted',ACCESS='sequential',&
            ACTION='read',STATUS='old',IOSTAT=ierr)
       IF (ierr /= 0) THEN
          WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
               TRIM(fNameSoilType), ierr
          STOP "**(ERROR)**"
       END IF

       INQUIRE (IOLENGTH=IOL) SoilType
       CALL  vfirec(nfsoiltp,SoilType,iMaxNew*jMaxNew,'LIN')
       IF (reducedGrid) THEN
          CALL NearestIJtoIBJB(SoilType,soil_type)
       ELSE
          CALL IJtoIBJB( SoilType,soil_type)
       END IF

       CLOSE(UNIT=nfsoiltp)

    END IF


    
    DEALLOCATE(SoilType)
    RETURN   
  END SUBROUTINE InitGridDump
  
  
  SUBROUTINE write_GridSigma_file (&
       ifday  ,tod   ,idate  ,idatec,&
       td0    ,tg0   ,tc0    ,z0    ,&
       convc  ,convt ,convb  ,gtsea ,&
       avisb  ,avisd ,anirb  ,anird ,&
       ustar  ,sm0   ,sheleg ,lsmk  ,&
       imask  ,mlsi  ,del    ,roperm,&
       namef  ,labeli,labelf ,extw  ,&
       exdw   ,trunc ,lev    ,nfprt  )
    IMPLICIT NONE
!
!
!  Ler/Escreve do dado de superficie do NCEP na Resolucao do MCGA/CPTEC
!  Espaco Fisico: Ponto de Grade
!
!  FIXIO - routines adapted from 2001 MRF code to read and write
!  surface analysis files.
!
!  07/21/04 - First implementation by Eric
!  07/26/04 - Major rewrite by Istvan: elimination of all but two
!             subroutines and several unncesessary variables and
!             parameters
!
!  SYNOPSIS: PARA_FIXO_R Read an CPTEC/NCEP surface file (unit number
!            nread). The surface fields are first read into a
!            real*4 buffer (buffer1, buffer2, or buffer4) and then
!            copied to the the appropriate real*8 arrays used in the
!            model calculation
!            PARA_FIXIO_W Writes an CPTEC/NCEP surface file (unit number
!            nfcpt). The real*8
!            fields are first copied into a real*4 buffer
!
    INTEGER           , INTENT(IN   ) :: ifday
    REAL(KIND=r8)     , INTENT(IN   ) :: tod
    INTEGER           , INTENT(IN   ) :: idate  (:)
    INTEGER           , INTENT(IN   ) :: idatec (:)
    REAL(KIND=r8), TARGET, INTENT(IN   ) :: td0    (ibMax,jbMax)
    REAL(KIND=r8), TARGET, INTENT(IN   ) :: tg0    (ibMax,jbMax)
    REAL(KIND=r8), TARGET, INTENT(IN   ) :: tc0    (ibMax,jbMax) 
    REAL(KIND=r8), TARGET, INTENT(IN   ) :: z0     (ibMax,jbMax)  
    REAL(KIND=r8), TARGET, INTENT(IN   ) :: convc  (ibMax,jbMax)  
    REAL(KIND=r8), TARGET, INTENT(IN   ) :: convt  (ibMax,jbMax)  
    REAL(KIND=r8), TARGET, INTENT(IN   ) :: convb  (ibMax,jbMax)  
    REAL(KIND=r8), TARGET, INTENT(IN   ) :: avisb  (ibMax,jbMax)  ! Visible beam surface albedo
    REAL(KIND=r8), TARGET, INTENT(IN   ) :: avisd  (ibMax,jbMax)  ! Visible diffuse surface albedo
    REAL(KIND=r8), TARGET, INTENT(IN   ) :: anirb  (ibMax,jbMax)  ! Near-ir beam surface albedo
    REAL(KIND=r8), TARGET, INTENT(IN   ) :: anird  (ibMax,jbMax)  ! Near-ir diffuse surface albedo
    REAL(KIND=r8), TARGET, INTENT(IN   ) :: ustar  (ibMax,jbMax)
    REAL(KIND=r8), TARGET, INTENT(IN   ) :: gtsea  (ibMax,jbMax)
    REAL(KIND=r8), TARGET, INTENT(IN   ) :: sm0    (ibMax,3,jbMax)
    REAL(KIND=r8), TARGET, INTENT(IN   ) :: sheleg (ibMax,jbMax)
    REAL(KIND=r8)        , INTENT(IN   ) :: lsmk   (:)  !sem uso.
    INTEGER(KIND=i8), TARGET, INTENT(IN   ) :: imask  (ibMax,jbMax)
    INTEGER(KIND=i8), TARGET, INTENT(IN   ) :: mlsi   (ibMax,jbMax)
    REAL(KIND=r8)     , INTENT(IN   ) :: del    (:)
    CHARACTER(LEN=200), INTENT(IN   ) :: roperm
    CHARACTER(LEN=  7), INTENT(IN   ) :: namef
    CHARACTER(LEN= 10), INTENT(IN   ) :: labeli
    CHARACTER(LEN= 10), INTENT(IN   ) :: labelf
    CHARACTER(LEN=  5), INTENT(IN   ) :: extw
    CHARACTER(LEN=  5), INTENT(IN   ) :: exdw
    CHARACTER(LEN=  *), INTENT(IN   ) :: trunc
    CHARACTER(LEN=  *), INTENT(IN   ) :: lev
    INTEGER           , INTENT(IN   ) :: nfprt
    !
    !  Local variables
    !
    CHARACTER(255)          :: fn
    INTEGER                 :: ifday4
    INTEGER(KIND=i4)        :: idat4   (4)
    REAL(KIND=r8)           :: tod4
    REAL(KIND=r4)           :: xhour    
    INTEGER(KIND=i4)        :: lplsfcOut(jMaxNew/2)
    REAL(KIND=r8)           :: vcover_g(ibMax,jbMax)

    INTEGER                 :: i,j,k,ik
    CHARACTER(LEN=8)        :: labfix(4)
    INTEGER(KIND=i4)        :: iMaxNew4,jMaxNew4
    CHARACTER(LEN=*), PARAMETER :: h="**(wrprog)**"
    CHARACTER(LEN=4)        :: str1

    REAL(KIND=r8), TARGET   :: work(ibMax,jbMax,15)
    INTEGER, ALLOCATABLE    :: interp_type(:)
    TYPE(p2d), TARGET, ALLOCATABLE :: fields(:)
    INTEGER, ALLOCATABLE    :: registro_dado(:)
    INTEGER                 :: nregtot
       
    vcover_g=sfc%vcover


    IF (myid.eq.0) THEN
       CALL opnfct(ifday,tod,idate,idatec,fn,&
            roperm,namef,labeli,extw,exdw,trunc,lev)
       ifday4=ifday
       tod4=tod
       DO k=1,4
          idat4 (k) = idate (k)
          !idat4c(k) = idatec(k)
       END DO
    END IF

    IF (myid == 0) THEN
       WRITE (nfprt, FMT='(/,A)') '----------------------------------'
       WRITE (nfprt, FMT='(/,3A,I3)') ' Writing SFC  from GDAS2 NCEP File ',&
                                         TRIM(roperm)//'/'//TRIM(namef)," Unit=OUTSWAP ",OUTSWAP
    ENDIF

    iMaxNew4=iMaxNew
    jMaxNew4=jMaxNew
    lplsfcOut(1:jMaxNew/2)=iMaxNew
    xhour = REAL(tod/3600.0_r8,kind=r4)
    labfix=CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)//&
           CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)

    IF (myid.eq.0) THEN
    WRITE(nfprt, '(A)')' In PARA_FIXO_W Write an NCEP Surface File'
    WRITE(nfprt, *    )'labfix:',labfix
    WRITE(nfprt, '(A)')' xhour,idat4,iMax,Jmax,version:'
    WRITE(nfprt, *    ) xhour,idat4,iMaxNew4,jMaxNew4,version,lplsfcOut
    WRITE(nfprt, *    )
    ENDIF

    IF (myid.eq.0) THEN   !add solange 28-07-2011 
    REWIND(OUTSWAP)
    WRITE (OUTSWAP)labfix
    WRITE (OUTSWAP)xhour,idat4,iMaxNew4,jMaxNew4,version,lplsfcOut
!    DO K=3,12
!      WRITE(OUTSWAP)
!    ENDDO
    ENDIF
    ! 
    ALLOCATE(fields(26))          ! Total de campos que serao arranjados segundo definicao em "nregtot"
    ALLOCATE(interp_type(26))
    ALLOCATE(registro_dado(26+1))   ! add solange 28-07-2011 
    nregtot=4                      ! Valor maximo de registro numa matriz para a saidada num arq seuqential
    registro_dado(26+1)=0
    !
    !  Surface temperature (water and land) [K]
    !  Campo 1
    !
    DO J=1,jbMax
       DO I=1,ibMaxPerJB(j)
          IF (mlsi(I,J) == 0_i8 .or. mlsi(I,J) == 2_i8 )THEN
             work(I,J,1)=-gtsea(I,J)
          ELSE
             work(i,j,1)=REAL(gtsea(I,J),kind=r8)
          END IF
       END DO
    END DO
    fields(1)%p => work(:,:,1)
    interp_type(1) = 3
    registro_dado(1) = 1
    !
    !
    !  Wetness of surface soil moisture content - 1,2 layers
    !  sm0(:,ik,1) e sm0(:,ik,2)
    !  Campo 2
    !
    ik=1
    DO K=1, 2
       str1='SM  '
       write(str1(3:3),'(i1.1)')ik
       CALL reord (sm0(:,ik,:),     1, work(:,:,k+1), 1, imask, gtsea, str1)
       ik=ik+2
    ENDDO
    fields(2)%p => work(:,:,2)          !work_in2(:,:,1)
    interp_type(2) = 1
    fields(3)%p => work(:,:,3)           !work_in2(:,:,2)
    interp_type(3) = 1 
    registro_dado(2) = 2
    registro_dado(3) = 2
    !
    !  Snow [m]
    !  Sheleg ou snow
    !  Campo 3
    !
    fields(4)%p => sheleg
    interp_type(4) = 1 
    registro_dado(4) = 3
    !
    !  Surface soil temperature: soil temperature - 1,2 layers 
    !       tg0: ground soil temperature  (water and land) - layer 1
    !       td0: deep soil temperature  (water and land) - layer 2
    !  Campo 4
    !
    !ik=1
    !str1='TD  '
    !DO K=1, 2
    !   if ( k == 1 ) THEN
    !   CALL reord (td0,     1, work, 1, imask, gtsea, str1)
    !   ELSE
    !   CALL reord (tg0,     1, work, 1, imask, gtsea, str1)
    !   ENDIF
    !   !work_in2(:,:,K)=work(:,:)
    !   ik=ik+2
    !ENDDO

    str1='TD  '
    CALL reord (td0,     1, work(:,:,4), 1, imask, gtsea, str1)
    fields(5)%p => work(:,:,4)
    interp_type(5) = 1

    CALL reord (tg0,     1, work(:,:,5), 1, imask, gtsea, str1)
    fields(6)%p => work(:,:,5)
    interp_type(6) = 1
    registro_dado(5:6) = 4
    !
    !  Canopy temperature  (K)
    !  Campo 5
    !
    CALL reord (tc0,    1, work(:,:,6), 1, imask, gtsea, 'TD  ')
    fields(7)%p => work(:,:,6)
    interp_type(7) = 1 
    registro_dado(7) = 5
    !
    !  Surface roughness
    !  Campo 6
    !
    fields(8)%p => z0
    interp_type(8) = 1 
    registro_dado(8) = 6
    !
    !  Convective cloud cover
    !  Campo 7
    !
    fields(9)%p => convc
    interp_type(9) = 1
    registro_dado(9) = 7
    !
    !  Convective cloud base
    !  Campo 8
    !
    fields(10)%p => convb
    interp_type(10) = 1
    registro_dado(10) = 8
    !
    !  Convective cloud top
    !  Campo 9
    !
    fields(11)%p => convt
    interp_type(11) = 1
    registro_dado(11) = 9
    !
    !  Albedo surface 
    !     avisb:  albedo surface Visible Beam - Direta  ? 
    !     avisd:  albedo surface Visible Difuse - Difusa?
    !     anirb:  albedo surface  Near infrared Beam - Direta  ?
    !     anird:  albedo surface Near infrared Difuse - Difusa ? 
    !  Campo 10
    !
    fields(12)%p => avisb
    interp_type(12) = 1
    fields(13)%p => avisd
    interp_type(13) = 1
    fields(14)%p => anirb
    interp_type(14) = 1
    fields(15)%p => anird
    interp_type(15) = 1
    registro_dado(12:15) = 10
    !
    !  Mask land/sea/ice (1/0/2)   mlsi
    !  Campo 11
    !
    work(:,:,7)=REAL(mlsi,kind=r8)
    fields(16)%p => work(:,:,7)
    interp_type(16) = 1
    registro_dado(16) = 11
    !if ( myid.eq.0) then
    !   WRITE(nfprt, '(/,A)' )'Mascara do Mar no PROC 0.  Mask land/sea/ice (1/0/2)   mlsi=work'
    !   DO J=1,jbMax,9    
    !     write(nfprt, *)(work(I,J,7),i=1,ibMax,3)
    !   ENDDO
    !endif

    !
    !  VFRAC: vegetation fraction 
    !  Campo 12
    !
    CALL reord (vcover_g(:,:),    1, work(:,:,8), 1, imask, gtsea, 'W1  ')
    fields(17)%p => work(:,:,8)
    interp_type(17) = 1
    registro_dado(17) = 12
    !
    !  CANOPY: surface canopy fraction (0-1.2)
    !  Campo 13
    !
    work(:,:,9)=-99.0_r8
    fields(18)%p => work(:,:,9)
    interp_type(18) = 1
    registro_dado(18) = 13
    !
    !  F10M: surface 10 meters canopy (0-1) 
    !  0.55 to 0.95
    !  Campo 14
    !
    work(:,:,10)=-99.0_r8
    fields(19)%p => work(:,:,10)
    interp_type(19) = 1
    registro_dado(19) = 14

    !
    !  VTYPE: vegetation type (1-13) 
    !  Campo 15
    !
    work(:,:,11)=REAL(imask,kind=r8)
    fields(20)%p => work(:,:,11)
    interp_type(20) = 1
    registro_dado(20) = 15
    !
    !  STYPE: soil type (1-9)
    !
    !  No Mar de Gelo, soil_type=9.
    !  No Bioma Gelo,  soil_type=9.
    !  Campo 16
    !
    !
    work(:,:,12)=REAL(soil_type,kind=r8)
    DO j=1,jbMax
       DO i=1,ibMaxPerJB(j)
          IF (  mlsi(i,j) ==  2_i8 ) work(i,j,12)=9.0_r8
          IF ( imask(i,j) == 13_i8 .or.imask(i,j) == 15_i8 ) work(i,j,12)=9.0_r8
       END DO
    END DO
    fields(21)%p => work(:,:,12)
    interp_type(21) = 1
    registro_dado(21) = 16
    !
    !  FACSF: surface fraction of arid soil
    !  FACWF: surface fraction of vegetation
    !  Campo 17
    !
    !
    work(:,:,13)=-99.0
    fields(22)%p => work(:,:,13)
    interp_type(22) = 1
    fields(23)%p => work(:,:,13)
    interp_type(23) = 1
    registro_dado(22:23) = 17
    !
    ! surface wind stress
    !  Campo 18
    !
    fields(24)%p => ustar
    interp_type(24) = 1
    registro_dado(24) = 18
    !
    ! surface ??
    ! FFMMIn(i,j)  ???
    !  Campo 19
    !
    work(:,:,14)=-99.0_r8
    fields(25)%p => work(:,:,14)
    interp_type(25) = 1
    registro_dado(25) = 19
    !
    ! surface ??
    ! FFHHIn	   ???
    !  Campo 20
    !
    work(:,:,15)=-99.0_r8
    fields(26)%p => work(:,:,15)
    interp_type(26) = 1
    registro_dado(26) = 20
    !
    !  Collect and print fields
    !
    !CALL Collect_Grid_Sur_Print2(fields,interp_type,26,0,OUTSWAP)  !modificado por solange. Add info sobre  "registro_dado".
    CALL Collect_Grid_Sur_Print2(fields,interp_type,registro_dado,26,nregtot,0,OUTSWAP)
    

    DEALLOCATE(fields)
    DEALLOCATE(interp_type)
    DEALLOCATE(registro_dado)
    CLOSE(OUTSWAP)
    RETURN
  END SUBROUTINE write_GridSigma_file   

  SUBROUTINE opnfct( ifday, tod, idate, idatec,fn,&
       roperm,namef,labeli,extw,exdw,trunc,lev)
    INTEGER           , INTENT(IN) :: ifday
    REAL(KIND=r8)     , INTENT(IN) :: tod
    INTEGER           , INTENT(IN) :: idate(4)
    INTEGER           , INTENT(IN) :: idatec(4)
    CHARACTER(LEN=  *), INTENT(INOUT) :: fn
    CHARACTER(LEN=200), INTENT(IN) :: roperm
    CHARACTER(LEN=  7), INTENT(IN) :: namef
    CHARACTER(LEN= 10), INTENT(IN) :: labeli
    CHARACTER(LEN=  5), INTENT(IN) :: extw
    CHARACTER(LEN=  5), INTENT(IN) :: exdw
    CHARACTER(LEN=  *), INTENT(IN) :: trunc
    CHARACTER(LEN=  *), INTENT(IN) :: lev

    INTEGER :: iyi
    INTEGER :: imi
    INTEGER :: idi
    INTEGER :: ihi
    INTEGER :: iyc
    INTEGER :: imc
    INTEGER :: idc
    INTEGER :: ihc
    LOGICAL :: inic
    INTEGER :: ierr
    INTEGER,            SAVE :: icall=1
    INTEGER,            SAVE :: is
    CHARACTER(LEN= 10)       :: labelc
    CHARACTER(LEN=  3), SAVE :: ext
    CHARACTER(LEN=  6)       :: extn
    CHARACTER(LEN= 10)       :: label
    CHARACTER(LEN=8) :: c0
    CHARACTER(LEN=*), PARAMETER :: modout="/"
    CHARACTER(LEN=LEN(namef)+LEN(labeli)+LEN(labelc)+&
         LEN(extn)+LEN(trunc)+LEN(lev)) :: fNameBin
    CHARACTER(LEN=*), PARAMETER :: h="**(opnfct)**"

    inic=(ifday.EQ.0 .AND. tod.EQ.0)

    IF (icall.EQ.1 .AND. inic) THEN
       ext='icn'
    ELSEIF (icall.EQ.2 .AND. inic) THEN
       icall=3
       ext='inz'
    ELSE
       ext='fct'
    ENDIF
    IF (icall .EQ. 1) THEN
       icall=2
       is=INDEX(roperm//' ',' ')-1
    ENDIF

    iyi=idate(4)
    imi=idate(2)
    idi=idate(3)
    ihi=idate(1)
    WRITE(label,'(i4.4,3i2.2)')iyi,imi,idi,ihi
    iyc=idatec(4)
    imc=idatec(2)
    idc=idatec(3)
    ihc=idatec(1)
    WRITE(labelc,'(i4.4,3i2.2)')iyc,imc,idc,ihc
    extn(1:2)=extw(1:2)
    extn(3:5)=ext(1:3)
    extn(6:6)='.'

    CLOSE(UNIT=OUTSWAP)

    ! open binary file

    fNameBin=namef//labeli//labelc//extn//TRIM(trunc)//TRIM(lev)
    fn=roperm(1:is)//modout//TRIM(fNameBin)
    OPEN(UNIT=OUTSWAP,FILE=roperm(1:is)//modout//TRIM(fNameBin),FORM='unformatted', &
         STATUS='unknown', position='rewind', IOSTAT=ierr)

    IF (ierr /= 0) THEN
        WRITE(c0,"(i8)") ierr
        CALL FatalError(h//" Open file "//roperm(1:is)//modout//TRIM(fNameBin)//&
            " returned iostat="//TRIM(ADJUSTL(c0)))
    END IF
  END SUBROUTINE opnfct
  
  SUBROUTINE WriteField82Dz(n, field)
    INTEGER, INTENT(IN)  :: n
    REAL   (KIND=r8), INTENT(IN)  :: field(:,:,:)
    REAL   (KIND=r4) :: raux3(SIZE(field,1),SIZE(field,2), SIZE(field,3))
    
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteField82D)**"
    INTEGER :: z, k, l
    INTEGER :: d1, d2,d3
    d3=SIZE(field,3);d2=SIZE(field,2);d1=SIZE(field,1)
    DO z = 1, d3
     DO k = 1, d2
       DO l = 1, d1
          raux3(l,k,z) = REAL(field(l,k,z),kind=r4)
       END DO
    END DO
    END DO
    WRITE(UNIT=n)raux3
  END SUBROUTINE WriteField82Dz
  
  SUBROUTINE WriteField82D(n, field)
    INTEGER, INTENT(IN)  :: n
    REAL   (KIND=r8), INTENT(IN)  :: field(:,:)
    REAL   (KIND=r4) :: raux3(SIZE(field,1),SIZE(field,2))
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteField82D)**"
    INTEGER :: k, l
    INTEGER :: d1, d2
    d2=SIZE(field,2);d1=SIZE(field,1)
    DO k = 1, d2
       DO l = 1, d1
         raux3(l,k) = REAL(field(l,k),kind=r4)
       END DO
    END DO
    WRITE(UNIT=n)raux3
  END SUBROUTINE WriteField82D

  SUBROUTINE WriteField81D(n, field)
    INTEGER, INTENT(IN)  :: n
    REAL   (KIND=r8), INTENT(IN)  :: field(:)
    REAL   (KIND=r4) :: raux3(SIZE(field,1))
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteField81D)**"
    INTEGER :: l
    INTEGER :: d1
    d1=SIZE(field,1)
    DO l = 1, d1
       raux3(l) = REAL(field(l),kind=r4)
    END DO
    WRITE(UNIT=n)raux3
  END SUBROUTINE WriteField81D

  SUBROUTINE WriteField82D_i8(n, field)
    INTEGER, INTENT(IN)  :: n
    INTEGER, INTENT(IN)  :: field(:,:)
    INTEGER :: raux3(SIZE(field,1),SIZE(field,2))
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteField82D)**"
    INTEGER :: k, l
    INTEGER :: d1, d2
    d2=SIZE(field,2);d1=SIZE(field,1)
    DO k = 1, d2
       DO l = 1, d1
          raux3(l,k) = field(l,k)
       END DO
    END DO
    WRITE(UNIT=n)raux3
  END SUBROUTINE WriteField82D_i8

  SUBROUTINE WriteField81D_i8(n, field)
    INTEGER, INTENT(IN)  :: n
    INTEGER, INTENT(IN)  :: field(:)
    INTEGER :: raux3(SIZE(field,1))
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteField81D)**"
    INTEGER :: l
    INTEGER :: d1
    d1=SIZE(field,1)
    DO l = 1, d1
       raux3(l) = field(l)
    END DO
    WRITE(UNIT=n)raux3
    
  END SUBROUTINE WriteField81D_i8


  SUBROUTINE reord (datum, dim2, work, lev, imask, tsea, ittl)
    INTEGER, INTENT(IN ) :: dim2
    REAL(KIND=r8)   , INTENT(in ) :: datum(ibMax,dim2,jbMax)
    REAL(KIND=r8)   , INTENT(OUT) :: work (ibMax,jbMax)
    INTEGER         , INTENT(IN ) :: lev
    INTEGER(KIND=i8), INTENT(IN ) :: imask   (ibMax,jbMax)
    REAL(KIND=r8),    INTENT(IN ) :: tsea    (ibMax,jbMax)
    CHARACTER(LEN=4), INTENT(IN) :: ittl

    INTEGER :: j
    INTEGER :: i
    INTEGER :: ncount
    LOGICAL :: case1
    LOGICAL :: case2
    LOGICAL :: case3

    case1 = ittl == 'TD  '
    case2 = ittl == 'W1  ' .OR. ittl == 'W2  ' .OR. ittl == 'W3  '
    case3 = ittl == 'SM1 ' .OR. ittl == 'SM2 ' .OR. ittl == 'SM3 '
    DO j = 1, jbMax
       ncount=0
       DO i = 1, ibMaxPerJB(j)
          IF (imask(i,j) >= 1_i8) THEN
             ncount = ncount + 1
             work(i,j) = datum(ncount,lev,j)
          ELSE IF (case1) THEN
             work(i,j)=ABS(tsea(i,j))
          ELSE IF (case2) THEN
             work(i,j)=1.0_r8
          ELSE IF (case3) THEN
             work(i,j)=0.47_r8
          ELSE
             work(i,j)=0.0_r8
          END IF
       END DO
    END DO
  END SUBROUTINE reord

END MODULE GridDump
