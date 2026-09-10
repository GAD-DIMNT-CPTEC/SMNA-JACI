!  $Author: pkubota e Solange $
!  $Date: 2011/07/28 14:34:46 $
!  $Revision: 2.0 $
!
MODULE SpecDump
  USE Constants, ONLY: &
       rk,  &
       i8,  &
       r8,  &
       r4,  &
       pie

  USE Options, ONLY: &
      reducedGrid, nfprt

  USE Parallelism, ONLY: &
       myid,             &
       maxNodes,         &
       MsgDump,          &
       MsgOne,           &
       FatalError
  USE Utils, ONLY:     &
       colrad
       
  USE Communications, ONLY: &
       Collect_Grid_Red,    &
       Collect_Grid_Sur,    &
       Collect_Spec

  USE Sizes, ONLY:  &
       ThreadDecomp,&
       kmaxloc,     &  !add solange p nova versao do global  28-07-2011
       havesurf,    &  !add solange p nova versao do global  28-07-2011
       kMax,        &
       iMax,        &
       jMax,        &
       jbMax,       &
       mMax,        &
       mnmax,       &
       iMaxPerJ,    &
       ibMaxPerJB
       !ibPerIJ,     &
       !jbPerIJ,     &

    USE FieldsDynamics, ONLY : &
         qozon


 IMPLICIT NONE
         
  PRIVATE
  !
  !  Single-precision reals in sigma files
  !
  !flag do fortran -W1 -T101 
  integer,parameter::D=kind(1.0d0)  ! IEEE double precision or similar
  INTEGER,PARAMETER :: NTRAC=2    ! relative humidity and ozone, presumably
  !
  !  Four units are used for input and output of binary data, depending
  !  on whether byte-swapping may be required.
  !  if desired.
  !
  INTEGER,PARAMETER :: OUTSWAP  = 95  !101
  !
  !  IN and OUT correspond to the standard input and output.
  !  STDERR is the standard error.  These are intended for use with text i/o.
  !  Unit 0 is an extension of the Fortran standard, but it is commonly used
  !  for the standard error.  However, I know of one implementation where
  !  the standard error is mapped to  unit 7.  Alter as appropriate.
  !
  INTEGER,PARAMETER :: STDERR=0 
  !
  !  Logical unit to which we write a log of our actions if verbose
  !  output is requested.  In this implementation, the log file coincides
  !  with the standard output.
  !
  INTEGER,PARAMETER :: LOGFILE=0
  !
  !
  !  VERBLEV - 0:silent, 1:basics, 2:some status, 3:chatty
  !  See set_verbosity and verbosity().
  !
  INTEGER,PARAMETER::VERB_SILENT  = 0
  INTEGER,PARAMETER::VERB_TERSE   = 1
  INTEGER,PARAMETER::VERB_STATUS  = 2
  INTEGER,PARAMETER::VERB_CHATTY  = 3
  INTEGER          :: verblev = VERB_TERSE

  !Como inserir kmax aqui?   kmax, excluindo sequence em TYPE sigma_file_header
  INTEGER, PARAMETER ::  NLEVS=64
  !
  !  This header is valid for the 2001 version of the GFS code that we use.
  !  We supply default values for most components, as the only ones that
  !  are likely to change in our simulations are the forecast hour
  !  and the date.  These default values simplify the construction of the
  !  data assimilation scripts because we don't need to pass around a reference
  !  sigma file for the XI and XL coefficients.
  !
  TYPE sigma_file_header
      SEQUENCE                     ! important for working around Lahey byte-swapping deficiency
      REAL(KIND=r4) :: fhour4     = 0.0  ! forecast hour
      INTEGER       :: idate(4)   = 0    ! hh dd mm yyyy
      !
      !  XI defines values that are used to compute WSLEV and energy units - (KMAX+1)
      !
      REAL(KIND=r4) :: xi (NLEVS+1)
      !
      !  XL defines the sigma levels   (KMAX)
      !
      REAL(KIND=r4) :: xl(NLEVS)
      !
      REAL(KIND=r4) :: dummy(201-(NLEVS+1)-NLEVS)=0.0      
      REAL(KIND=r4) :: waves      =062      
      REAL(KIND=r4) :: xlayers    =NLEVS
      REAL(KIND=r4) :: trun       =1.0
      REAL(KIND=r4) :: order      =2.0
      REAL(KIND=r4) :: realform   =1.0
      REAL(KIND=r4) :: gencode    =80.0     
      REAL(KIND=r4) :: rlond      =192.0    
      REAL(KIND=r4) :: rlatd      =94.0     
      REAL(KIND=r4) :: rlonp      =192.0
      REAL(KIND=r4) :: rlatp      =94.0     
      REAL(KIND=r4) :: rlonr      =192.0
      REAL(KIND=r4) :: rlatr      =94.0     
      REAL(KIND=r4) :: tracers    =3.0
      REAL(KIND=r4) :: subcen     =0.0
      REAL(KIND=r4) :: ensemble(2)=0.0
      REAL(KIND=r4) :: ppid       =0.0
      REAL(KIND=r4) :: slid       =0.0
      REAL(KIND=r4) :: vcid       =0.0
      REAL(KIND=r4) :: vmid       =0.0
      REAL(KIND=r4) :: vtid       =21.0
      REAL(KIND=r4) :: runid4     =0.0
      REAL(KIND=r4) :: usrid4     =0.0
      !
      !  The next three fields are broken out in the 2004 GFS; the 2001 GFS
      !  has DUMMY2 as length 21.
      !
      REAL(KIND=r4)::pdryini4=0.0
      REAL(KIND=r4)::xncld=1.0
      REAL(KIND=r4)::xgf=0.0
      REAL(KIND=r4)::dummy2(18)=0.0
      
  END TYPE sigma_file_header
  

  TYPE(sigma_file_header) :: header   !acho que poderia sair daqui! solange
  !
  !
  !  The following is needed for the messy workaround to Lahey runtime 
  !  byte-swapping bugs.  SIGMA_HEADER_LENGTH is the number of default real
  !  and integer elements in the file header.
  !
  INTEGER,PARAMETER::SIGMA_HEADER_LENGTH=250
  !
  !  On Lahey, the environment variable FORT90L does not enable byte-swapping
  !  of components of derived types.  So, we create a dummy real(KIND=r4) array
  !  of length SIGMA_HEADER_LENGTH, which is the current size of the sigma file
  !  header record, read this array (which the Lahey runtime byte swaps), then 
  !  copy the bits to the dummy argument (sigh).  This code will work
  !  correctly on systems without this problem at the price of slight
  !  inefficiency, provided that REAL(KIND=r4) and INTEGER occupy the same
  !  amount of storage (which they should, assuming KIND=r4 is default real).
  !
  REAL(KIND=r4):: lahey_bug_workaround(SIGMA_HEADER_LENGTH)



  PUBLIC :: InitSpecDump
  PUBLIC :: write_sigma_file






CONTAINS
  
  

   SUBROUTINE InitSpecDump(si,sl,idatec,idate)
    IMPLICIT NONE
    !
    !  SPECTRAL_CODE_INIT - initialize the GFS spectral transforms.
    !  FILE  : if present, then we read in the header from FILE; otherwise,
    !     we use the default values supplied above.
    !  HEADER := if present, we fill in the values using the sigma file header
    !     from FILE or the default header values otherwise.
    !  COEFF := if present, the spectral coefficients from the sigma file;
    !     FILE must be specified.
    !


    INTEGER, INTENT(IN   ) :: idatec(4)
    INTEGER, INTENT(IN   ) :: idate(4)
    TYPE(sigma_file_header)::xi(kmax+1),xl(kmax), dummy(201-(kMax+1)-kMax)
    REAL(KIND=r8),    INTENT(IN) :: si(kmax+1)
    REAL(KIND=r8),    INTENT(IN) :: sl(kmax)

    header%fhour4                   = float(idatec(1)-idate(1))         ! forecast hour
    IF ( header%fhour4 < 0 ) header%fhour4=24+header%fhour4
    header%idate(1:4)              = (/idate(1),idate(2),idate(3),idate(4)/)   ! idatec=hh dd mm yyyy labelc ANL  DATE
    !
    !  XI defines values that are used to compute WSLEV and energy units
    !
    header%xi(1:(kmax+1)) = real( si(1:(kmax+1)), KIND=r4)
    !
    !  XL defines the sigma levels
    !
    header%xl(1:kmax) = real( sl(1:kmax), KIND=r4)
    !
    header%dummy(1:)                = 0.0
    
    header%waves                    = mMax-1
    header%xlayers                  = kMax
    header%trun                     = 1.0
    header%order                    = 2.0
    header%realform                 = 1.0
    header%gencode                  = 80.0
    header%rlond                    = real(iMax,KIND=r4)      !iMax
    header%rlatd                    = real(jMax,KIND=r4)      !jMax
    header%rlonp                    = real(iMax,KIND=r4)      !iMax
    header%rlatp                    = real(jMax,KIND=r4)      !jMax
    header%rlonr                    = real(iMax,KIND=r4)      !iMax
    header%rlatr                    = real(jMax,KIND=r4)      !jMax
    header%tracers                  = 3.0                     !2.0  ! relative humidity and ozone, presumably
    header%subcen                   = 0.0
    header%ensemble(1:2)            = 0.0
    header%ppid                     = 0.0
    header%slid                     = 0.0
    header%vcid                     = 0.0
    header%vmid                     = 0.0
    header%vtid                     = 21.0
    header%runid4                   = 0.0
    header%usrid4                   = 0.0
    !
    !  The next three fields are broken out in the 2004 GFS; the 2001 GFS
    !  has DUMMY2 as length 21.
    !
    header%pdryini4                 = 0.0
    header%xncld                    = 1.0
    header%xgf                      = 0.0
    header%dummy2(1:18)             = 0.0




    if (myid.eq.0) THEN

    write(nfprt,*)  "waves  mMax-1,kMax,iMax,jMax",mMax-1,kMax,iMax,jMax
    write(nfprt,*)   "SIGMA_HEADER_LENGTH: ",SIGMA_HEADER_LENGTH
    write(nfprt,*) '======================='   
    write(nfprt,*)'header%waves             ' , header%waves  
    write(nfprt,*)'header%xlayers            ' , header%xlayers 
    write(nfprt,*)'header%trun      ' , header%trun  
    write(nfprt,*)'header%order     ' , header%order
    write(nfprt,*)'header%realform            ' , header%realform
    write(nfprt,*)'header%gencode            ' , header%gencode 
    write(nfprt,*)'header%rlond             ' , header%rlond  
    write(nfprt,*)'header%rlatd             ' , header%rlatd  
    write(nfprt,*)'header%rlonp             ' , header%rlonp  
    write(nfprt,*)'header%rlatp             ' , header%rlatp  
    write(nfprt,*)'header%rlonr     ' , header%rlonr
    write(nfprt,*)'header%rlatr     ' , header%rlatr
    write(nfprt,*)'header%tracers            ' , header%tracers 
    write(nfprt,*)'header%subcen            ' , header%subcen  
    write(nfprt,*)'header%ensemble(1:2) ' , header%ensemble(1:2) 
    write(nfprt,*)'header%ppid      ' , header%ppid 
    write(nfprt,*)'header%slid      ' , header%slid 
    write(nfprt,*)'header%vcid      ' , header%vcid 
    write(nfprt,*)'header%vmid      ' , header%vmid
    write(nfprt,*)'header%vtid      ' , header%vtid  
    write(nfprt,*)'header%runid4    ' , header%runid4
    write(nfprt,*)'header%usrid4            ' , header%usrid4 
    write(nfprt,*)'header%pdryini4            ' , header%pdryini4 
    write(nfprt,*)'header%xncld     ' , header%xncld
    write(nfprt,*)'header%xgf       ' , header%xgf 
    write(nfprt,*)'header%dummy2(1:18)  ' , header%dummy2(1:18)
    write(nfprt,*)'header%xi        ' , header%xi
    write(nfprt,*)'header%xl        ' , header%xl
    write(nfprt,*)'kmax        ' , kmax, ' Size(header%xl) ',SIZE(header%xl)
    write(nfprt,*) '======================='


    endif

    
    RETURN   
  END SUBROUTINE InitSpecDump
  
  
  
  
  
  
  
  SUBROUTINE write_sigma_file (ifday ,tod   ,idate ,idatec,qrot         , &
       qdiv  ,qq    ,qlnp  ,qtmp  ,qgzs  ,o3mix,del   ,ijMaxGauQua  , &
       imax  ,jmax  ,ibMax,jbMax,kMax,roperm,namef ,labeli,labelf,extw  ,exdw  ,trunc , &
       lev, si,sl,havesurf,kmaxloc)
       !solange add 04/02/2009 lev)
       
    IMPLICIT NONE
    !  WRITE_SIGMA_FILE - write out the spectral coefficients in sigma file format.
    !  Note: the sigma file format uses single-precision values in the 2001/2004
    !  versions of the GFS; we denote the loss of precision with explicit
    !  conversions.
    !  FN  : name of the output file (any predecessor is overwritten).
    !  HEADER : sigma file header
    !  COEFF : spectral coefficients to write out.
    !
    !  Side effects:  opens and closes OUTSWAP.
    !    
    INTEGER           , INTENT(IN   ) :: ifday
    REAL(KIND=r8)     , INTENT(IN   ) :: tod
    INTEGER           , INTENT(IN   ) :: idate (:)
    INTEGER           , INTENT(IN   ) :: idatec(:)
    REAL(KIND=r8)     , INTENT(IN   ) :: qrot  (:,:)
    REAL(KIND=r8)     , INTENT(IN   ) :: qdiv  (:,:)
    REAL(KIND=r8)     , INTENT(IN   ) :: qq    (:,:)
    REAL(KIND=r8)     , INTENT(IN   ) :: qlnp  (:)
    REAL(KIND=r8)     , INTENT(IN   ) :: qtmp  (:,:)
    REAL(KIND=r8)     , INTENT(IN   ) :: qgzs  (:)
    REAL(KIND=r8)     , INTENT(IN   ) :: o3mix (:,:,:)
    REAL(KIND=r8)     , INTENT(IN   ) :: del   (:)
    INTEGER           , INTENT(IN   ) :: ijMaxGauQua
    INTEGER           , INTENT(IN   ) :: iMax
    INTEGER           , INTENT(IN   ) :: jMax
    INTEGER           , INTENT(IN   ) :: ibMax
    INTEGER           , INTENT(IN   ) :: jbMax
    INTEGER           , INTENT(IN   ) :: kMax
    INTEGER           , INTENT(IN   ) :: kMaxloc
    LOGICAL           , INTENT(IN   ) :: havesurf
    CHARACTER(LEN=200), INTENT(IN   ) :: roperm
    CHARACTER(LEN=  7), INTENT(IN   ) :: namef
    CHARACTER(LEN= 10), INTENT(IN   ) :: labeli
    CHARACTER(LEN= 10), INTENT(IN   ) :: labelf
    CHARACTER(LEN=  5), INTENT(IN   ) :: extw
    CHARACTER(LEN=  5), INTENT(IN   ) :: exdw
    CHARACTER(LEN=  *), INTENT(IN   ) :: trunc
    CHARACTER(LEN=  *), INTENT(IN   ) :: lev
    REAL(KIND=r8),    INTENT(IN) :: si(kmax+1)
    REAL(KIND=r8),    INTENT(IN) :: sl(kmax)
    !
    !  Local variables
    !
    CHARACTER(255)          :: fn
    !REAL(KIND=r8)           :: qspec1  (2*mnmax)
    REAL(KIND=r8)           :: qspec2  (2*mnmax,kMax)
    REAL(KIND=r8)           :: qspec   (2*mnmax,kMax)
    INTEGER                 :: ifday4
    INTEGER                 :: idat4   (4)
    INTEGER                 :: idat4c  (4)
    REAL(KIND=r8)           :: tod4
    !old: CHARACTER(8),PARAMETER::LAB(4)=&
    !old:     CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)
    !solange : atencao falta corrigir !!!!!!!
    CHARACTER(LEN=8),PARAMETER::LAB(4)='\0\0\0\0\0\0\0\0'
    !CHARACTER(LEN=8),PARAMETER::LAB(4)  !falta testar

    !solange REAL(KIND=r8) :: buff(2,mnmax)
    REAL(KIND=r4) :: buff(2,mnmax)
    !
    REAL(D) :: gz (mnmax,2)                     ! geopotential height
    REAL(D) :: q  (mnmax,2)                     ! surface pressure
    REAL(D) :: te (mnmax,2)                ! temperature
    REAL(D) :: di (mnmax,2)                ! divergence
    REAL(D) :: ze (mnmax,2)                ! vorticity
    REAL(D) :: rq (mnmax,2)                ! mixing ratio
    REAL(D) :: oz (mnmax,2,kMax,NTRAC)       ! ozone plus other tracers (ntrac=1->Ozone; ntrac=2->Conteudo de agua liquida)
    INTEGER :: i,j,k,ij
    !
    !  On Lahey, the environment variable FORT90L does not enable byte-swapping
    !  of components of derived types.  So, we create a dummy real(KIND=r4) array
    !  of length SIGMA_HEADER_LENGTH, which is the current size of the sigma file
    !  header record, read this array (which the Lahey runtime byte swaps), then 
    !  copy the bits to the dummy argument (sigh).  This code will work
    !  correctly on systems without this problem at the price of slight
    !  inefficiency, provided that REAL(KIND=r4) and INTEGER occupy the same
    !  amount of storage (which they should, assuming KIND=r4 is default real).
    !
    real(KIND=r4)::lahey_bug_workaround(SIGMA_HEADER_LENGTH)
    !
    INTEGER :: n

    CHARACTER(LEN=*), PARAMETER :: h="**(wrprog)**"
       

    !falta testar!
    !LAB=CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)//&
    !       CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)



    IF (myid.eq.0) THEN

       write(nfprt,*)  "mnmax,kMax,NTRAC:", mnmax,kMax,NTRAC, 'ifday:',ifday,' tod:', tod
       CALL opnfct(ifday,tod,idate,idatec,fn,&
            roperm,namef,labeli,extw,exdw,trunc,lev)

       ifday4=ifday
       tod4=tod
       DO k=1,4
          idat4 (k) = idate (k)
          idat4c(k) = idatec(k)
       END DO


!solange add 04/fev/2009 call initdump
      CALL InitSpecDump(si,sl,idatec,idate)
!fim add solange

      
       CALL dump_header(fn,LOGFILE,header,verbosity())

       !
       !  Write a NUL label at the beginning of the file.
       !
    END IF
    !
    !  Sigma file header.
    !  Lahey does not byte-swap derived types in binary i/o,
    !  which is a PITA.  So we copy the header to a sufficiently large array,
    !  which the Lahey runtime byte-swaps.  This code is portable, since
    !  the header is a sequence type, and will work correctly (with a slight
    !  loss of efficiency) on implementations that don't have Lahey's problem.
    ! 
    lahey_bug_workaround=TRANSFER(header,lahey_bug_workaround)
    !
    ! write fields at OUTSWAP
    !IF (maxnodes.gt.1) THEN
       IF (myid == 0) THEN
          WRITE(OUTSWAP)  LAB
          WRITE(OUTSWAP) lahey_bug_workaround  ! to get the byte-swapped header
       END IF
       !
       !     topography
       !
       !CALL Collect_Spec(qgzs, qspec1, 1, 0)  !old     
       IF (havesurf) CALL Collect_Spec(qgzs, qspec, 1, 1, 0)  !solange: nova versao do mcga 28-07-2001
 
       IF (myid == 0) THEN
          DO i=1,mnmax 
             gz(i,1)=real(qspec(2*i-1,1), D)
             gz(i,2)=real(  qspec(2*i,1), D)
          ENDDO
          buff=REAL(TRANSPOSE(gz(:,:)),KIND(buff))
          WRITE(OUTSWAP) buff
          !write(nfprt,*) 'Geopotential height coeff%gz:',maxval(buff),minval(buff)
       ENDIF       
       !
       !     ln surface pressure
       !
       ! CALL Collect_Spec(qlnp, qspec1, 1, 0) old
       IF (havesurf) CALL Collect_Spec(qlnp, qspec, 1, 1, 0)  !solange: nova versao do mcga 28-07-2001       
       IF (myid.eq.0) THEN          
          DO i=1,mnmax 
             q(i,1)=real(qspec(2*i-1,1), D)
             q(i,2)=real(  qspec(2*i,1), D)
          ENDDO          
          buff=REAL(TRANSPOSE(q(:,:)),KIND(buff))
          WRITE(OUTSWAP) buff
          !write(nfprt,*) 'Presao a sfc coeff%q:',maxval(buff),minval(buff),buff(1,mnmax/2),buff(2,mnmax/2)
       ENDIF
       !
       !     virtual temperature
       !
       !
       !CALL Collect_Spec(qtmp, qspec, kmax, 0)  !old
       !CALL Collect_Spec(qtmp, qspec, kmaxloc, kmaxnew, 0)  !solange: nova versao do mcga 28-07-2001
       CALL Collect_Spec(qtmp, qspec, kmaxloc, kmax, 0)  !solange: nova versao do mcga 28-07-2001
       IF (myid.eq.0) THEN
          DO k=1,kMax
             te = 0.0
             DO i=1,mnmax 
                te(i,1)=real(qspec((2*i-1),k), D)    !te(i,1)=buff(2*i-1)
                te(i,2)=real(  qspec((2*i),k), D)    !te(i,2)=buff(2*i  )
             ENDDO             
             buff=REAL(TRANSPOSE(te(:,:)),KIND(buff))
             WRITE(OUTSWAP) buff
          END DO                    
       ENDIF
       !
       !  Divergence and vorticity
       !
       !
       !     divergence
       !
       !write(nfprt,*) 'divergence'
       !CALL Collect_Spec(qdiv, qspec2, kmax, 0)  !old
       !CALL Collect_Spec(qdiv, qspec2, kmaxloc, kmaxnew, 0)  !solange: nova versao do mcga 28-07-2001
       CALL Collect_Spec(qdiv, qspec2, kmaxloc, kmax, 0)  !solange: nova versao do mcga 28-07-2001
       !
       !     vorticity
       !
       !write(nfprt,*) ' vorticity'
       !CALL Collect_Spec(qrot, qspec, kmax, 0)
       !CALL Collect_Spec(qrot, qspec, kmaxloc, kmaxnew, 0)  !solange: nova versao do mcga 28-07-2001
       CALL Collect_Spec(qrot, qspec, kmaxloc, kmax, 0)  !solange: nova versao do mcga 28-07-2001
       IF (myid.eq.0) THEN
          DO k=1,kMax
             di=0.0
             ze=0.0
              DO i=1,mnmax 
                di(i,1)=real(qspec2((2*i-1),k), D)
                di(i,2)=real(  qspec2((2*i),k), D)
                ze(i,1)=real(qspec((2*i-1),k), D)
                ze(i,2)=real(  qspec((2*i),k), D)
              ENDDO
             buff=REAL(TRANSPOSE(di(:,:)),KIND(buff))
             WRITE(OUTSWAP) buff
             buff=REAL(TRANSPOSE(ze(:,:)),KIND(buff))
             WRITE(OUTSWAP) buff
          END DO          
       END IF
       !
       !     specific humidity  -> deve ser umidade relativa, razao de mistura !!!!!
       !
       !write(nfprt,*) 'specific humidity'
       !CALL Collect_Spec(qq  , qspec, kmax, 0) !old
       !CALL Collect_Spec(qq  , qspec, kmaxloc, kmaxnew, 0)  !solange: nova versao do mcga 28-07-2001
       CALL Collect_Spec(qq  , qspec, kmaxloc, kmax, 0)  !solange: nova versao do mcga 28-07-2001
       IF (myid.eq.0) THEN
          DO k=1,kMax
             rq=0.0
             DO i=1,mnmax 
                rq(i,1)=real(qspec((2*i-1),k), D)
                rq(i,2)=real(  qspec((2*i),k), D)
             ENDDO             
             buff=REAL(TRANSPOSE(rq(:,:)),KIND(buff))
             WRITE(OUTSWAP) buff
             !write(nfprt,*) 'Relative humidity coeff%rq(:,:,k):',k,maxval(buff),minval(buff),buff(1,mnmax/2),buff(2,mnmax/2)             
          END  DO
       ENDIF
       !
       !  Ozone and other tracers
       !
       !       N=1-> Ozone, which I think is the first of two tracers in the 2004 GFS
       !
       !       N=2->Liquid water content, which is the second tracer, I think.
       !
       !write(nfprt,*) 'Ozone and other tracers'
       !CALL Collect_Spec(qozon  , qspec, kMax, 0) !old
       !CALL Collect_Spec(qozon  , qspec, kmaxloc, kmaxnew, 0)  !solange: nova versao do mcga 28-07-2001
       CALL Collect_Spec(qozon  , qspec, kmaxloc, kmax, 0)  !solange: nova versao do mcga 28-07-2001
       IF (myid.eq.0) THEN
          !solange DO n=1,NTRAC           
           DO n=1,1
             DO k=1,kMax
                oz=0.0
                DO i=1,mnmax 
                   oz(i,1,k,n)=real(qspec( (2*i-1),k), D)   
                   oz(i,2,k,n)=real(qspec(   (2*i),k), D)   
                ENDDO
                buff=REAL(TRANSPOSE(oz(:,:,k,n)),KIND(buff))
                WRITE(OUTSWAP) buff
                !write(nfprt,*)'Ozone and other tracers  coeff%oz(:,:,k,traci): traci=',n,k,maxval(buff),minval(buff),buff(1,mnmax/2),buff(2,mnmax/2)
                
             END  DO
           END DO
           oz(:,:,:,2) = 0.0
           DO n=2,NTRAC
             DO k=1,kMax
                buff=REAL(TRANSPOSE(oz(:,:,k,n)),KIND(buff))
                WRITE(OUTSWAP) buff
                !write(nfprt,*) 'Ozone and other tracers  coeff%oz(:,:,k,traci): traci=',n,k,maxval(buff),minval(buff),buff(1,mnmax/2),buff(2,mnmax/2)
             ENDDO
           END DO
           write(nfprt,*) 'Atencao !!!! oz(:,:,;,2)=0.0 FALTA INFO DE  conteudo de agua liquida'
       ENDIF
       CLOSE(OUTSWAP)
    !END IF   !descomentei 
    RETURN   
  END SUBROUTINE write_sigma_file   
  !
  !
  !
  SUBROUTINE opnfct( ifday, tod, idate, idatec,fn,&
       roperm,namef,labeli,extw,exdw,trunc,lev)

    INTEGER           , INTENT(IN) :: ifday
    REAL(KIND=r8)     , INTENT(IN) :: tod
    INTEGER           , INTENT(IN) :: idate(4)
    INTEGER           , INTENT(IN) :: idatec(4)
    CHARACTER(LEN=  *), INTENT(INOUT) :: fn
    !CHARACTER(LEN=200), INTENT(IN) :: roperm
    CHARACTER(LEN=*), INTENT(IN) :: roperm
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
!solange em 07/fev/2012
!    INTEGER,            SAVE :: icall=1
!    INTEGER,            SAVE :: is
    CHARACTER(LEN= 10)       :: labelc
!    CHARACTER(LEN=  3), SAVE :: ext
    CHARACTER(LEN=  6)       :: extn
!    CHARACTER(LEN= 23), SAVE :: modout
    CHARACTER(LEN= 10)       :: label
!    CHARACTER(LEN=*), PARAMETER :: modout="/"
!fim com solange

    INTEGER,            SAVE :: icall=1
    INTEGER,            SAVE :: is
!    CHARACTER(LEN= 10), SAVE :: labelc
    CHARACTER(LEN=  3), SAVE :: ext
!    CHARACTER(LEN=  6), SAVE :: extn
    CHARACTER(LEN=  6), SAVE :: exdn
    CHARACTER(LEN= 23), SAVE :: modout
!    CHARACTER(LEN= 10), SAVE :: label


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
    modout='/'
    
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

    fn=roperm(1:is)//TRIM(modout)//namef//labeli//labelc//extn//TRIM(trunc)//TRIM(lev)
    WRITE(UNIT=nfprt,FMT="(a,1x,i5)" ) fn,OUTSWAP


    OPEN(UNIT=OUTSWAP,FILE=fn,&
         FORM='UNFORMATTED', &
         ACCESS='sequential',ACTION='write', STATUS='replace', IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            fn, ierr
       STOP "**(ERROR)**"
    END IF


  END SUBROUTINE opnfct
  
  SUBROUTINE dump_header(fn,iunit,header,iverb)
    !  DUMP_HEADER - dump file header info for sigma file FN to the indicated unit.
    !
    CHARACTER(*),INTENT(in)            :: fn
    INTEGER     ,INTENT(in)            :: iunit
    TYPE(sigma_file_header),INTENT(in) :: header
    INTEGER     ,INTENT(in),OPTIONAL   :: iverb
    !
    INTEGER::ichat
    !
    IF(PRESENT(iverb)) THEN
       ichat=iverb
    ELSE
       ichat=verbosity()
    ENDIF

    IF(ichat.EQ.VERB_SILENT) RETURN
    !
    !  Display forcast hour and date only unless greater verbosity is requested.
    !
    IF(iverb.GT.VERB_TERSE) THEN
       WRITE(iunit,100) TRIM(fn)
       WRITE(iunit,101) 'idate: ',header%idate,'fhour4: ',header%fhour4
    ENDIF
100 FORMAT('header info for: ',a)
101 FORMAT(a,4i5,3x,a,f8.3)
    !
    IF(iverb.LE.VERB_CHATTY) RETURN
    !
    write(nfprt,*) 'iunit:',iunit,iverb,VERB_TERSE
    
    WRITE(nfprt,*) 'kMax:',kMax,' xi',header%xi(1:kMax)
    WRITE(iunit,102) 'xi',header%xi(1:kMax)
102 FORMAT(a,1x,7es11.3/(3x,7es11.3))
    WRITE(iunit,102) 'xl',header%xl
    !
    WRITE(iunit,103) 'waves: ',header%waves,'xlayers: ',header%xlayers
103 FORMAT(2(a9,1x,es11.3,1x))
    WRITE(iunit,103) 'trun: ',header%trun,'order: ',header%order
    WRITE(iunit,103) 'gencode: ',header%gencode,'tracers: ',header%tracers
    WRITE(iunit,103) 'rlond: ',header%rlond,'rlatd: ',header%rlatd
    WRITE(iunit,103) 'rlonp: ',header%rlonp,'rlatp: ',header%rlatp
    WRITE(iunit,103) 'rlonr: ',header%rlonr,'rlatr: ',header%rlatr
    WRITE(iunit,103) 'subcen: ',header%subcen,'ppid: ',header%ppid
    WRITE(iunit,103) 'slid: ',header%slid,'vcid: ',header%vcid
    WRITE(iunit,103) 'vmid: ',header%vmid,'vtid: ',header%vtid
    WRITE(iunit,103) 'runid4: ',header%runid4,'usrid4: ',header%usrid4
    WRITE(iunit,103) 'pdryini4: ',header%pdryini4,'xncld: ',header%xncld
    WRITE(iunit,103) 'xgf: ',header%xgf
    RETURN
  END SUBROUTINE dump_header
  !----------------------------------------------------------------------------
  PURE INTEGER FUNCTION verbosity()
    !
    !  VERBOSITY - return the verbosity level.
    !
    verbosity=verblev
    RETURN
  END FUNCTION verbosity
  !---------------------------------------------------------------------------
  PURE FUNCTION rmsof(f)
    REAL(KIND=r8)::rmsof
    REAL(KIND=r8),INTENT(in)::f(mnmax,2)
    !
    rmsof=SQRT(0.5*SUM(f**2))
    RETURN
  END FUNCTION rmsof

END MODULE SpecDump
