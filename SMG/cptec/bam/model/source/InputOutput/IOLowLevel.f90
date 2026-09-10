!
!  $Author: pkubota $
!  $Date: 2010/10/21 10:46:22 $
!  $Revision: 1.14 $
!
!  Modified by Tarassova 2015
!  Climate aerosol (Kinne, 2013) coarse mode included
!  Modifications (2) are marked by 
!  !tar begin and  !tar end
!
!
MODULE IOLowLevel

  USE Constants, ONLY: &
       r4,i4, r8, ndavl, ndrq, ncdg,pai

  USE Options, ONLY: &
       nfprt, nfctrl, nfsst,nfco2fx,nfslm, nferr, nfwaves, nfgauss,labelsi,labelsj, labelsi_soilm,&
       labelsj_soilm,NMSST,NMSOILM,NMCO2FX,DHFCT,CalDate,labelsi_flxco2,labelsj_flxco2,schemes,atmpbl,jull         

  USE Parallelism, ONLY: &
       myid,             &
       MsgOne,           &
       MsgOut,           &
       MsgDump,          &
       FatalError

   IMPLICIT NONE
  SAVE       


  PRIVATE

  PUBLIC :: ReadHead  
  PUBLIC :: GReadHead
  PUBLIC :: ReadField
  PUBLIC :: GReadField
  PUBLIC :: WriteHead  
  PUBLIC :: GWriteHead
  PUBLIC :: WriteField
  PUBLIC :: GWriteField
  PUBLIC :: FWriteField
  PUBLIC :: FWriteFieldb   !add solange 28-07-2011.
  PUBLIC :: WriteDir
  PUBLIC :: WriteDire
  PUBLIC :: WriteDirDhn
  PUBLIC :: ReadProgHead 
  PUBLIC :: GReadProgHead
  PUBLIC :: WriteProgHead
  PUBLIC :: ReadLandSeaMask2
  PUBLIC :: ReadVar
  PUBLIC :: ReadGetALB
  PUBLIC :: ReadGetSST
  PUBLIC :: ReadGetSST2
  PUBLIC :: ReadGetSLM
  PUBLIC :: ReadGetSLM3D
  PUBLIC :: ReadGetSNW
  PUBLIC :: ReadGetNFTGZ
  PUBLIC :: ReadMs
  PUBLIC :: WriteMs
  PUBLIC :: ReadGauss
  PUBLIC :: WriteGauss
!tar begin
!climate aerosol reading
  PUBLIC :: Read_Aeros  
!tar end
  PUBLIC :: ReadOzone !hmjb
  PUBLIC :: ReadTracer !hmjb
  PUBLIC :: InitReadWriteSpec
  PUBLIC :: WriteDiagHead
  PUBLIC :: LandSeaMask
  PUBLIC :: LandSoilmMask
  PUBLIC :: LandFlxCO2Mask
  PUBLIC :: WriteGrdHist
  PUBLIC :: WrTopoGrdHist

  INTERFACE WriteDiagHead
     MODULE PROCEDURE WriteDiagHead4, WriteDiagHead8
  END INTERFACE
  INTERFACE ReadOzone
     MODULE PROCEDURE ReadOzone8
  END INTERFACE
  INTERFACE ReadTracer
     MODULE PROCEDURE ReadTracer8,ReadTracer8s
  END INTERFACE
  INTERFACE ReadGetNFTGZ
     MODULE PROCEDURE ReadNFTGZ4, ReadNFTGZ8
  END INTERFACE
  INTERFACE ReadGetSNW 
     MODULE PROCEDURE ReadSNW4, ReadSNW8
  END INTERFACE
  INTERFACE ReadGetSLM 
     MODULE PROCEDURE ReadSLM4, ReadSLM8
  END INTERFACE
  INTERFACE ReadGetSLM3D
     MODULE PROCEDURE ReadSLM43D, ReadSLM83D
  END INTERFACE
  INTERFACE ReadGetALB
     MODULE PROCEDURE ReadAlb4, ReadAlb8
  END INTERFACE
  INTERFACE ReadGetSST 
     MODULE PROCEDURE ReadSST4, ReadSST8 
  END INTERFACE
  INTERFACE ReadGetSST2 
     MODULE PROCEDURE ReadSST4Rec, ReadSST8Rec
  END INTERFACE
  INTERFACE ReadVar 
     MODULE PROCEDURE ReadVar4, ReadVar8
  END INTERFACE
  INTERFACE GReadHead
     MODULE PROCEDURE GReadHead4, GReadHead8
  END INTERFACE
  INTERFACE ReadHead
     MODULE PROCEDURE ReadHead4, ReadHead8
  END INTERFACE
  INTERFACE WriteHead
     MODULE PROCEDURE WriteHead4, WriteHead8
  END INTERFACE
  INTERFACE GWriteHead
     MODULE PROCEDURE GWriteHead4, GWriteHead8
  END INTERFACE
  INTERFACE ReadProgHead
     MODULE PROCEDURE ReadProgHead4, ReadProgHead8
  END INTERFACE
  INTERFACE GReadProgHead
     MODULE PROCEDURE GReadProgHead4, GReadProgHead8
  END INTERFACE
  INTERFACE WriteProgHead
     MODULE PROCEDURE WriteProgHead4, WriteProgHead8
  END INTERFACE
  INTERFACE WriteDir
     MODULE PROCEDURE WriteDir8, WriteDir4
  END INTERFACE
  INTERFACE WriteDire
     MODULE PROCEDURE WriteDire4, WriteDire8 ,LogWriteDire4, LogWriteDire8 
  END INTERFACE
  INTERFACE WriteDirDhn
     MODULE PROCEDURE WriteDirDhn4, WriteDirDhn8 ,LogWriteDirDhn4, LogWriteDirDhn8 
  END INTERFACE
  INTERFACE GReadField
     MODULE PROCEDURE GReadField41D, GReadField42D, GReadField81D, GReadField82D
  END INTERFACE
  INTERFACE ReadField
     MODULE PROCEDURE ReadField41D, ReadField42D, ReadField81D, ReadField82D
  END INTERFACE
  INTERFACE WriteField
     MODULE PROCEDURE WriteField41D, WriteField42D, WriteField81D, WriteField82D
  END INTERFACE
  INTERFACE GWriteField
     MODULE PROCEDURE GWriteField41D, GWriteField42D, GWriteField81D, GWriteField82D
  END INTERFACE
  INTERFACE FWriteField
     MODULE PROCEDURE FWriteField42D, FWriteField82D
  END INTERFACE
! add solange para o letkf 28-07-2011
  INTERFACE FWriteFieldb
     MODULE PROCEDURE FWriteFieldb42D, FWriteFieldb82D, FWriteFieldb82Dz
  END INTERFACE
! fim add solange
  INTERFACE ReadLandSeaMask2
     MODULE PROCEDURE ReadLandSeaMask2_4, ReadLandSeaMask2_8
  END INTERFACE
  INTERFACE LandSeaMask
     MODULE PROCEDURE LandSeaMask4, LandSeaMask8
  END INTERFACE
  INTERFACE LandFlxCO2Mask
     MODULE PROCEDURE LandFlxCO2Mask4, LandFlxCO2Mask8
  END INTERFACE
  INTERFACE LandSoilmMask
     MODULE PROCEDURE LandSoilmMask4, LandSoilmMask8
  END INTERFACE
  INTERFACE  WriteGrdHist
     MODULE PROCEDURE WriteGrdH4, WriteGrdH8
  END INTERFACE
  INTERFACE  WrTopoGrdHist
     MODULE PROCEDURE WrTopoGrdH4,WrTopoGrdH8 
  END INTERFACE
  
  INTEGER    :: mMax 
  INTEGER    :: mnMax
  INTEGER    :: kMax 
  INTEGER    :: ijMax
  INTEGER    :: iMax 
  INTEGER    :: jMax 
  INTEGER    :: ibMax
  INTEGER    :: jbMax

CONTAINS
  SUBROUTINE InitReadWriteSpec(&
       mMax_in ,mnMax_in,kMax_in ,ijMax_in,iMax_in ,&
       jMax_in ,ibMax_in,jbMax_in)

    INTEGER          , INTENT(IN   ) :: mMax_in 
    INTEGER          , INTENT(IN   ) :: mnMax_in
    INTEGER          , INTENT(IN   ) :: kMax_in 
    INTEGER          , INTENT(IN   ) :: ijMax_in
    INTEGER          , INTENT(IN   ) :: iMax_in 
    INTEGER          , INTENT(IN   ) :: jMax_in 
    INTEGER          , INTENT(IN   ) :: ibMax_in
    INTEGER          , INTENT(IN   ) :: jbMax_in
    mMax  =    mMax_in 
    mnMax =   mnMax_in
    kMax  =    kMax_in 
    ijMax =   ijMax_in
    iMax  =    iMax_in 
    jMax  =    jMax_in 
    ibMax =   ibMax_in
    jbMax =   jbMax_in

  END SUBROUTINE InitReadWriteSpec

  SUBROUTINE ReadHead4(n, ifday, tod, idate, idatec, a, b, kMax)
    INTEGER, INTENT(IN)  :: n
    INTEGER, INTENT(IN)  :: kMax
    INTEGER(KIND=i4), INTENT(OUT) :: ifday
    REAL   (KIND=r4), INTENT(OUT) :: tod
    INTEGER(KIND=i4), INTENT(OUT) :: idate(4)
    INTEGER(KIND=i4), INTENT(OUT) :: idatec(4)
    REAL   (KIND=r4), INTENT(OUT) :: a(:)
    REAL   (KIND=r4), INTENT(OUT) :: b(:)
    REAL   (KIND=r4) :: a1(kmax+1), b1(kmax+1)
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadHead4)**"
    INTEGER :: lMax
    lMax=SIZE(a)
    IF (kMax /= lMax-1) THEN
       WRITE (UNIT=nferr, FMT='(2(A,I3))') ' kMax = ',kMax, ' is /= lMax-1', lMax-1
       STOP h
    END IF
    READ(UNIT=n)ifday, tod, idate, idatec, a1, b1
!   reverse coefficients to be from top to bottom
!
!   a(1:kmax+1) = a1(kmax+1:1:-1) / 100.
    a(1:kmax+1) = a1(kmax+1:1:-1)
    b(1:kmax+1) = b1(kmax+1:1:-1)
  END SUBROUTINE ReadHead4
  SUBROUTINE ReadHead8(n, ifday, tod, idate, idatec, a, b, kMax)
    INTEGER, INTENT(IN)  :: n
    INTEGER, INTENT(IN)  :: kMax
    INTEGER, INTENT(OUT) :: ifday
    REAL   (KIND=r8), INTENT(OUT) :: tod
    INTEGER, INTENT(OUT) :: idate(4)
    INTEGER, INTENT(OUT) :: idatec(4)
    REAL   (KIND=r8), INTENT(OUT) :: a(:)
    REAL   (KIND=r8), INTENT(OUT) :: b(:)
    INTEGER(KIND=i4) :: iaux(8),ifday4
    REAL   (KIND=r4) :: raux1(kmax+1), raux2(kmax+1),tod4
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadHead8)**"
    INTEGER :: lMax
    lMax=SIZE(a)
    IF (kMax /= lMax-1) THEN
       WRITE (UNIT=nferr, FMT='(2(A,I3))') ' kMax = ',kMax, ' is /= lMax-1', lMax-1
       STOP h
    END IF
    write(*,*) 'in ReadHead8 '
    READ(UNIT=n)ifday4,tod4,iaux, raux1, raux2
    write(*,*) 'in ReadHead8 , read ifday4', ifday4
    ifday  = ifday4
    tod    = tod4
    idate  = INT(iaux(1:4))
    idatec = INT(iaux(5:8))
    a     = REAL(raux1 , r8)
    b     = REAL(raux2 , r8)
!   reverse coefficients to be from top to bottom
!
!   a(1:kmax+1) = REAL(raux1(kmax+1:1:-1), r8) / 100.
    a(1:kmax+1) = REAL(raux1(kmax+1:1:-1), r8)
    b(1:kmax+1) = REAL(raux2(kmax+1:1:-1), r8)
  END SUBROUTINE ReadHead8
  SUBROUTINE GReadHead4(n, ifday, tod, idate, idatec, a, b, kMax)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    INTEGER(KIND=i4), INTENT(IN)  :: kMax
    INTEGER(KIND=i4), INTENT(OUT) :: ifday
    REAL   (KIND=r4), INTENT(OUT) :: tod
    INTEGER(KIND=i4), INTENT(OUT) :: idate(4)
    INTEGER(KIND=i4), INTENT(OUT) :: idatec(4)
    REAL   (KIND=r4), INTENT(OUT) :: a(:)
    REAL   (KIND=r4), INTENT(OUT) :: b(:)
    REAL   (KIND=r4) :: a1(kmax+1), b1(kmax+1)
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadHead4)**"
    INTEGER :: lMax
    lMax=SIZE(a)
    IF (kMax /= lMax-1) THEN
       WRITE (UNIT=nferr, FMT='(2(A,I3))') ' kMax = ',kMax, ' is /= lMax-1', lMax-1
       STOP h
    END IF
    READ(UNIT=n)ifday, tod, idate, idatec, a1, b1
!   reverse coefficients to be from top to bottom
!
!   a(1:kmax+1) = a1(kmax+1:1:-1) / 100.
    a(1:kmax+1) = a1(1:kmax+1:1)
    b(1:kmax+1) = b1(1:kmax+1:1)
  END SUBROUTINE GReadHead4
  SUBROUTINE GReadHead8(n, ifday, tod, idate, idatec, a, b, kMax)
    INTEGER, INTENT(IN)  :: n
    INTEGER, INTENT(IN)  :: kMax
    INTEGER, INTENT(OUT) :: ifday
    REAL   (KIND=r8), INTENT(OUT) :: tod
    INTEGER, INTENT(OUT) :: idate(4)
    INTEGER, INTENT(OUT) :: idatec(4)
    REAL   (KIND=r8), INTENT(OUT) :: a(:)
    REAL   (KIND=r8), INTENT(OUT) :: b(:)
    REAL   (KIND=r8) :: a1(kmax+1), b1(kmax+1)
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadHead8)**"
    INTEGER :: lMax
    lMax=SIZE(a)
    IF (kMax /= lMax-1) THEN
       WRITE (UNIT=nferr, FMT='(2(A,I3))') ' kMax = ',kMax, ' is /= lMax-1', lMax-1
       STOP h
    END IF
    READ(UNIT=n)ifday, tod, idate, idatec, a1, b1
!   reverse coefficients to be from top to bottom
!
!   a(1:kmax+1) = a1(kmax+1:1:-1) / 100.
    a(1:kmax+1) = a1(1:kmax+1:1)
    b(1:kmax+1) = b1(1:kmax+1:1)
  END SUBROUTINE GReadHead8

  SUBROUTINE ReadField42D(n, field)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    REAL   (KIND=r4), INTENT(OUT) :: field(:,:)
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadField42D)**"
    INTEGER :: k
    INTEGER :: d2
    d2 = SIZE(field,2)
    DO k = 1, d2
       READ(UNIT=n)field(:,k)
    END DO
  END SUBROUTINE ReadField42D
  SUBROUTINE ReadField82D(n, field)
    INTEGER, INTENT(IN)  :: n
    REAL   (KIND=r8), INTENT(OUT) :: field(:,:)
    REAL   (KIND=r4) :: raux3(SIZE(field,1))
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadField82D)**"
    INTEGER :: k
    INTEGER :: d2
    d2 = SIZE(field,2)
    DO k=1, d2
       READ(UNIT=n)raux3
       field(:,k) = REAL(raux3, r8)
    END DO
  END SUBROUTINE ReadField82D
  SUBROUTINE ReadField41D(n, field)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    REAL   (KIND=r4), INTENT(OUT) :: field(:)
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadField41D)**"
    READ(UNIT=n)field
  END SUBROUTINE ReadField41D
  SUBROUTINE ReadField81D(n, field)
    INTEGER, INTENT(IN)  :: n
    REAL   (KIND=r8), INTENT(OUT) :: field(:)
    REAL   (KIND=r4) :: raux3(SIZE(field))
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadField81D)**"
    READ(UNIT=n)raux3
    field = REAL(raux3, r8)
  END SUBROUTINE ReadField81D


  SUBROUTINE GReadField42D(n, field)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    REAL   (KIND=r4), INTENT(OUT) :: field(:,:)
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadField42D)**"
    INTEGER :: k
    INTEGER :: d2
    d2 = SIZE(field,2)
    DO k = 1, d2
       READ(UNIT=n)field(:,k)
    END DO
  END SUBROUTINE GReadField42D
  SUBROUTINE GReadField82D(n, field)
    INTEGER, INTENT(IN)  :: n
    REAL   (KIND=r8), INTENT(OUT) :: field(:,:)
    REAL   (KIND=r8) :: raux3(SIZE(field,1))
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadField82D)**"
    INTEGER :: k
    INTEGER :: d2
    d2 = SIZE(field,2)
    DO k=1, d2
       READ(UNIT=n)raux3
       field(:,k) =raux3
    END DO
  END SUBROUTINE GReadField82D
  SUBROUTINE GReadField41D(n, field)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    REAL   (KIND=r4), INTENT(OUT) :: field(:)
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadField41D)**"
    READ(UNIT=n)field
  END SUBROUTINE GReadField41D
  SUBROUTINE GReadField81D(n, field)
    INTEGER, INTENT(IN)  :: n
    REAL   (KIND=r8), INTENT(OUT) :: field(:)
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadField81D)**"
    READ(UNIT=n)field
  END SUBROUTINE GReadField81D

  SUBROUTINE WriteHead4(n, ifday, tod, idate, idatec, a, b)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    INTEGER(KIND=i4), INTENT(IN)  :: ifday
    REAL   (KIND=r4), INTENT(IN)  :: tod
    INTEGER(KIND=i4), INTENT(IN)  :: idate(4)
    INTEGER(KIND=i4), INTENT(IN)  :: idatec(4)
    REAL   (KIND=r4), INTENT(IN)  :: a(:)
    REAL   (KIND=r4), INTENT(IN)  :: b(:)
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteHead4)**"
    WRITE(UNIT=n)ifday, tod, idate, idatec, a, b
  END SUBROUTINE WriteHead4
  SUBROUTINE WriteHead8(n, ifday, tod, idate, idatec, a, b)
    INTEGER, INTENT(IN)  :: n
    INTEGER, INTENT(IN)  :: ifday
    REAL   (KIND=r8), INTENT(IN)  :: tod
    INTEGER, INTENT(IN)  :: idate(4)
    INTEGER, INTENT(IN)  :: idatec(4)
    REAL   (KIND=r8), INTENT(IN)  :: a(:)
    REAL   (KIND=r8), INTENT(IN)  :: b(:)
    INTEGER(KIND=i4) :: iaux(10)
    REAL   (KIND=r4) :: raux1(kmax+1), raux2(kmax+1)
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteHead8)**"
    iaux(  1 ) = INT (ifday , i4)
    iaux(  2 ) = INT (tod   , i4)
    iaux(3:6 ) = INT (idate , i4)
    iaux(7:10) = INT (idatec, i4)
    raux1      = REAL(a    , r4)
    raux2      = REAL(b    , r4)
    WRITE(UNIT=n)iaux, raux1, raux2
  END SUBROUTINE WriteHead8

  SUBROUTINE GWriteHead4(n, ifday, tod, idate, idatec, a, b)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    INTEGER(KIND=i4), INTENT(IN)  :: ifday
    REAL   (KIND=r4), INTENT(IN)  :: tod
    INTEGER(KIND=i4), INTENT(IN)  :: idate(4)
    INTEGER(KIND=i4), INTENT(IN)  :: idatec(4)
    REAL   (KIND=r4), INTENT(IN)  :: a(:)
    REAL   (KIND=r4), INTENT(IN)  :: b(:)
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteHead4)**"
    WRITE(UNIT=n)ifday, tod, idate, idatec, a, b
  END SUBROUTINE GWriteHead4
  SUBROUTINE GWriteHead8(n, ifday, tod, idate, idatec, a, b)
    INTEGER, INTENT(IN)  :: n
    INTEGER, INTENT(IN)  :: ifday
    REAL   (KIND=r8), INTENT(IN)  :: tod
    INTEGER, INTENT(IN)  :: idate(4)
    INTEGER, INTENT(IN)  :: idatec(4)
    REAL   (KIND=r8), INTENT(IN)  :: a(:)
    REAL   (KIND=r8), INTENT(IN)  :: b(:)
    INTEGER :: iaux(10)
    REAL   (KIND=r8) :: raux1(kmax+1), raux2(kmax+1)
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteHead8)**"
    iaux(  1 ) = INT (ifday)
    iaux(  2 ) = INT (tod)
    iaux(3:6 ) = INT (idate)
    iaux(7:10) = INT (idatec)
    raux1      = REAL(a,r8)
    raux2      = REAL(b,r8)
    WRITE(UNIT=n)iaux, raux1, raux2
  END SUBROUTINE GWriteHead8

  SUBROUTINE WriteField42D(n, field)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    REAL   (KIND=r4), INTENT(IN)  :: field(:,:)
    REAL   (KIND=r4) :: raux3(SIZE(field,1))
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteField42D)**"
    INTEGER :: k, l
    INTEGER :: d1, d2
    d2=SIZE(field,2);d1=SIZE(field,1)
    DO k = 1, d2
       DO l = 1, d1
          raux3(l) = field(l,k)
       END DO
       WRITE(UNIT=n)raux3(:)
    END DO
  END SUBROUTINE WriteField42D
  SUBROUTINE WriteField82D(n, field)
    INTEGER, INTENT(IN)  :: n
    REAL   (KIND=r8), INTENT(IN)  :: field(:,:)
    REAL   (KIND=r4) :: raux3(SIZE(field,1))
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteField82D)**"
    INTEGER :: k, l
    INTEGER :: d1, d2
    d2=SIZE(field,2);d1=SIZE(field,1)
    DO k = 1, d2
       DO l = 1, d1
          raux3(l) = REAL(field(l,k),r4)
       END DO
       WRITE(UNIT=n)raux3(:)
    END DO
  END SUBROUTINE WriteField82D


  SUBROUTINE WriteField41D(n, field)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    REAL   (KIND=r4), INTENT(IN)  :: field(:)
    REAL   (KIND=r4) :: raux3(SIZE(field,1))
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteField41D)**"
    INTEGER :: l
    INTEGER :: d1
    d1=SIZE(field,1)
    DO l = 1, d1
       raux3(l) = field(l)
    END DO
    WRITE(UNIT=n)raux3(:)
  END SUBROUTINE WriteField41D
  SUBROUTINE WriteField81D(n, field)
    INTEGER, INTENT(IN)  :: n
    REAL   (KIND=r8), INTENT(IN)  :: field(:)
    REAL   (KIND=r4) :: raux3(SIZE(field,1))
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteField81D)**"
    INTEGER :: l
    INTEGER :: d1
    d1=SIZE(field,1)
    DO l = 1, d1
       raux3(l) = REAL(field(l),r4)
    END DO
    WRITE(UNIT=n)raux3(:)
  END SUBROUTINE WriteField81D


  SUBROUTINE GWriteField42D(n, field)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    REAL   (KIND=r4), INTENT(IN)  :: field(:,:)
    REAL   (KIND=r4) :: raux3(SIZE(field,1))
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteField42D)**"
    INTEGER :: k, l
    INTEGER :: d1, d2
    d2=SIZE(field,2);d1=SIZE(field,1)
    DO k = 1, d2
       DO l = 1, d1
          raux3(l) = field(l,k)
       END DO
       WRITE(UNIT=n)raux3(:)
    END DO
  END SUBROUTINE GWriteField42D
  SUBROUTINE GWriteField82D(n, field)
    INTEGER, INTENT(IN)  :: n
    REAL   (KIND=r8), INTENT(IN)  :: field(:,:)
    REAL   (KIND=r8) :: raux3(SIZE(field,1))
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteField82D)**"
    INTEGER :: k, l
    INTEGER :: d1, d2
    d2=SIZE(field,2);d1=SIZE(field,1)
    DO k = 1, d2
       DO l = 1, d1
          raux3(l) = field(l,k)
       END DO
       WRITE(UNIT=n)raux3(:)
    END DO
  END SUBROUTINE GWriteField82D


  SUBROUTINE GWriteField41D(n, field)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    REAL   (KIND=r4), INTENT(IN)  :: field(:)
    REAL   (KIND=r4) :: raux3(SIZE(field,1))
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteField41D)**"
    INTEGER :: l
    INTEGER :: d1
    d1=SIZE(field,1)
    DO l = 1, d1
       raux3(l) = field(l)
    END DO
    WRITE(UNIT=n)raux3(:)
  END SUBROUTINE GWriteField41D
  SUBROUTINE GWriteField81D(n, field)
    INTEGER, INTENT(IN)  :: n
    REAL   (KIND=r8), INTENT(IN)  :: field(:)
    REAL   (KIND=r8) :: raux3(SIZE(field,1))
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteField81D)**"
    INTEGER :: l
    INTEGER :: d1
    d1=SIZE(field,1)
    DO l = 1, d1
       raux3(l) = field(l)
    END DO
    WRITE(UNIT=n)raux3(:)
  END SUBROUTINE GWriteField81D

  SUBROUTINE FWriteField42D(n, field)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    REAL   (KIND=r4), INTENT(IN)  :: field(:,:)
    REAL   (KIND=r4) :: raux3(SIZE(field,1)*SIZE(field,2))
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteField42D)**"
    INTEGER :: k, l
    INTEGER :: d1, d2
    d2=SIZE(field,2);d1=SIZE(field,1)
    DO k = 1, d2
       DO l = 1, d1
          raux3(l+(k-1)*d1) = field(l,k)
       END DO
    END DO
    WRITE(UNIT=n)raux3(:)
  END SUBROUTINE FWriteField42D

  SUBROUTINE FWriteField82D(n, field)
    INTEGER, INTENT(IN)  :: n
    REAL   (KIND=r8), INTENT(IN)  :: field(:,:)
    REAL   (KIND=r4) :: raux3(SIZE(field,1)*SIZE(field,2))
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteField82D)**"
    INTEGER :: k, l
    INTEGER :: d1, d2
    d2=SIZE(field,2);d1=SIZE(field,1)
    DO k = 1, d2
       DO l = 1, d1
          raux3(l+(k-1)*d1) = REAL(field(l,k),kind=r4)
       END DO
    END DO
    WRITE(UNIT=n)raux3(:)
  END SUBROUTINE FWriteField82D

! add solange 28-07-2011 para letkf
  SUBROUTINE FWriteFieldb42D(n, field)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    REAL   (KIND=r4), INTENT(IN)  :: field(:,:)
    REAL   (KIND=r4) :: raux3(SIZE(field,1)*SIZE(field,2))
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteField42D)**"
    INTEGER :: k, l
    INTEGER :: d1, d2
    d2=SIZE(field,2);d1=SIZE(field,1)
    DO k = 1, d2
       DO l = 1, d1
          raux3(l+(k-1)*d1) = field(l,k)
       END DO
    END DO
    WRITE(UNIT=n)raux3
  END SUBROUTINE FWriteFieldb42D

  SUBROUTINE FWriteFieldb82D(n, field)
    INTEGER, INTENT(IN)  :: n
    REAL   (KIND=r8), INTENT(IN)  :: field(:,:)
    REAL   (KIND=r4) :: raux3(SIZE(field,1)*SIZE(field,2))
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteField82D)**"
    INTEGER :: k, l
    INTEGER :: d1, d2
    d2=SIZE(field,2);d1=SIZE(field,1)
    DO k = 1, d2
       DO l = 1, d1
          raux3(l+(k-1)*d1) = REAL(field(l,k),KIND=r4)
       END DO
    END DO
    WRITE(UNIT=n)raux3
  END SUBROUTINE FWriteFieldb82D

  SUBROUTINE FWriteFieldb82Dz(n, field)
    INTEGER, INTENT(IN)  :: n
    REAL   (KIND=r8), INTENT(IN)  :: field(:,:,:)
    REAL   (KIND=r4) :: raux3(SIZE(field,1)*SIZE(field,2), SIZE(field,3))
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteField82D)**"
    INTEGER :: z, k, l
    INTEGER :: d1, d2,d3
    d3=SIZE(field,3);d2=SIZE(field,2);d1=SIZE(field,1)
    DO z = 1, d3
     DO k = 1, d2
       DO l = 1, d1
          raux3(l+(k-1)*d1,z) = REAL(field(l,k,z),KIND=r4)
       END DO
    END DO
    END DO
    WRITE(UNIT=n)raux3
  END SUBROUTINE FWriteFieldb82Dz
! fim add solange para letkf

 SUBROUTINE WriteDir4(n    , idate , idatec, &
                      a_hybr, b_hybr  , tod   , ifday, &
                      reqdg, combf , dodia, &
                      itcf , itavl , iavrq, &
                      nucf , lvrq  , nurq , &
                      lvcf , mxavl , icf  , &
                      opt)
   INTEGER(KIND=i4) , INTENT(IN ) :: n
   INTEGER(KIND=i4) , INTENT(IN ) :: idate(4)
   INTEGER          , INTENT(IN ) :: idatec(4)
   REAL   (KIND=r4) , INTENT(IN ) :: a_hybr(kMax+1)
   REAL   (KIND=r4) , INTENT(IN ) :: b_hybr(kMax+1)
   REAL   (KIND=r4) , INTENT(IN ) :: tod
   INTEGER(KIND=i4) , INTENT(IN ) :: ifday
   CHARACTER(len=40), INTENT(IN ) :: reqdg(1:ndrq)
   CHARACTER(len=40), INTENT(IN ) :: combf(1:ncdg)
   LOGICAL            , INTENT(IN ) :: dodia(1:ndavl)
   INTEGER(KIND=i4) , INTENT(IN ) :: itcf (1:ncdg)
   INTEGER(KIND=i4) , INTENT(IN ) :: itavl(1:ndavl)
   INTEGER(KIND=i4) , INTENT(IN ) :: iavrq(1:ndavl)
   INTEGER(KIND=i4) , INTENT(IN ) :: nucf (1:ncdg)
   INTEGER(KIND=i4) , INTENT(IN ) :: lvrq (1:ndrq )
   INTEGER(KIND=i4) , INTENT(IN ) :: nurq (1:ndrq )
   INTEGER(KIND=i4) , INTENT(IN ) :: lvcf (1:ncdg )
   INTEGER(KIND=i4) , INTENT(IN ) :: mxavl    
   INTEGER(KIND=i4) , INTENT(IN ) :: icf
   INTEGER(KIND=i4)                 :: isg(2)
   LOGICAL            , INTENT(IN ) :: opt
   CHARACTER (LEN= 5) :: imdl
   CHARACTER (LEN=150) :: jttl
   CHARACTER (LEN=21), PARAMETER :: ittl='CPTEC HIBRID VERS 2.0'
   CHARACTER (LEN= 4), PARAMETER :: nexp='0001'
   CHARACTER (LEN=40), PARAMETER :: orogra = 'TOPOGRAPHY' 
   CHARACTER (LEN=40), PARAMETER :: LSImsk = 'LAND SEA ICE MASK'
   CHARACTER (LEN=40), PARAMETER :: lnsurf = 'LN SURFACE PRESSURE'
   CHARACTER (LEN=40), PARAMETER :: divrgn = 'DIVERGENCE'
   CHARACTER (LEN=40), PARAMETER :: vortic = 'VORTICITY'
   CHARACTER (LEN=40), PARAMETER :: spechu = 'SPECIFIC HUMIDITY'
   CHARACTER (LEN=40), PARAMETER :: tempvi = 'VIRTUAL TEMPERATURE'
   CHARACTER (LEN=40), PARAMETER :: surfte = 'SURFACE TEMPERATURE'
   CHARACTER (LEN=40), PARAMETER :: srfrou = 'ROUGHNESS LENGTH'
   CHARACTER (LEN=40), PARAMETER :: deepte = 'DEEP SOIL TEMPERATURE'
   CHARACTER (LEN=40), PARAMETER :: stcnpy = 'STORAGE ON CANOPY'
   CHARACTER (LEN=40), PARAMETER :: stgrnd = 'STORAGE ON GROUND'
   CHARACTER (LEN=40), PARAMETER :: wt1soi = 'SOIL WETNESS OF SURFACE'
   CHARACTER (LEN=40), PARAMETER :: wt2soi = 'SOIL WETNESS OF ROOT ZONE'
   CHARACTER (LEN=40), PARAMETER :: wt3soi = 'SOIL WETNESS OF DRAINAGE ZONE' 
   CHARACTER (LEN=40), PARAMETER :: temp2m = 'TEMPERATURE AT 2-M FROM SURFACE' 
   CHARACTER (LEN=40), PARAMETER :: umes2m = 'SPECIFIC HUMIDITY AT 2-M FROM SURFACE' 
   CHARACTER (LEN=40), PARAMETER :: uve10m = 'ZONAL WIND AT 10-M FROM SURFACE'
   CHARACTER (LEN=40), PARAMETER :: vve10m = 'MERID WIND AT 10-M FROM SURFACE'
   CHARACTER (LEN=40), PARAMETER :: veget0 = 'MASK VEGETATION'
   CHARACTER (LEN=40), PARAMETER :: soiltx = 'MASK SOIL TEXTURE CLASSES'
   CHARACTER (LEN=40), PARAMETER :: doair0 = 'PARTIAL OXYGEN DENSITY '
   !JGERD -> fields necessary to data assimilation
   CHARACTER (LEN=40), PARAMETER :: sfcstg = 'SURFACE SOIL TEMPERATURE'
   CHARACTER (LEN=40), PARAMETER :: vegcov = 'VEGETATION COVER'
   CHARACTER (LEN=40), PARAMETER :: snodpt = 'SNOW DEPTH'
   CHARACTER (LEN=40), PARAMETER :: lqwcon = 'LIQ MIXING RATIO PROGNOSTIC'
   CHARACTER (LEN=40), PARAMETER :: icwcon = 'ICE MIXING RATIO PROGNOSTIC'
   CHARACTER (LEN=40), PARAMETER :: concld = 'CLOUD TOTAL PROGNOSTIC'
   CHARACTER (LEN=29), PARAMETER :: fmt1='(A40,2X,A4,2X,I8,3X,I4,4X,I3)'
   CHARACTER (LEN=4 ), PARAMETER :: diag='DIAG'
   INTEGER :: m
   INTEGER :: nn
   INTEGER :: ix
   INTEGER :: ihr
   INTEGER :: iday
   INTEGER :: mon
   INTEGER :: iyr
   CHARACTER(LEN=10) :: labelc
   CHARACTER(LEN=10) :: labeln
   ihr  = idatec(1)
   iday = idatec(3)
   mon  = idatec(2)
   iyr  = idatec(4)
   IF(DHFCT < 0 .and. opt )THEN
      WRITE(labelc,'(I4.4,I2.2,I2.2,I2.2)')iyr,mon,iday,ihr
      labeln=CalDate(labelc,'-')
      READ(labeln(1: 4),'(I4.4)')iyr
      READ(labeln(5: 6),'(I2.2)')mon
      READ(labeln(7: 8),'(I2.2)')iday
      READ(labeln(9:10),'(I2.2)')ihr
   END IF
   isg(1)=iMax*jMax
   isg(2)=2*mnMax

   jttl='CPTEC AGCM REVIS 1.0 2000  T    L    COLD [OPTIONS |surface=01||pbl=04|]'
   WRITE (jttl(29:32), '(i4.4)') mMax-1
   WRITE (jttl(34:36), '(i2.2)') kMax
   WRITE (jttl(61:62), '(i2.2)') schemes
   WRITE (jttl(69:70), '(i2.2)') atmpbl

   WRITE (imdl, '(A1,I4.4)') 'T', mMax-1

   WRITE (UNIT=n, FMT='(A21)')   ittl
   WRITE (UNIT=n, FMT='(A4,1X,  A4,1X, A5,1X,11I5,1X,A4)') &
                         nexp, 'SEQU', imdl, mMax, kmax+1, kmax+1, &
         ihr, iday, mon, iyr, idate, 'TAPE'
   WRITE (UNIT=n, FMT='(A)')   jttl
   WRITE (UNIT=n, FMT='(5E16.8)')   a_hybr
   WRITE (UNIT=n, FMT='(5E16.8)')   b_hybr
   WRITE (UNIT=n, FMT=fmt1) orogra, 'FIXD', 2*mnMax, 1, 10
   WRITE (UNIT=n, FMT=fmt1) LSImsk, 'FIXD', ijmax, 1, 0
   WRITE (UNIT=n, FMT=fmt1) lnsurf, 'PROG', 2*mnMax, 1, 142
   WRITE (UNIT=n, FMT=fmt1) divrgn, 'PROG', 2*mnMax, kmax, 50
   WRITE (UNIT=n, FMT=fmt1) vortic, 'PROG', 2*mnMax, kmax, 50
   WRITE (UNIT=n, FMT=fmt1) spechu, 'PROG', 2*mnMax, kmax, 0
   WRITE (UNIT=n, FMT=fmt1) tempvi, 'PROG', 2*mnMax, kmax, 40
   WRITE (UNIT=n, FMT=fmt1) srfrou, 'PROG', ijmax, 1, 10
   WRITE (UNIT=n, FMT=fmt1) surfte, 'PROG', ijmax, 1, 40
   WRITE (UNIT=n, FMT=fmt1) deepte, 'PROG', ijmax, 1, 40
   WRITE (UNIT=n, FMT=fmt1) stcnpy, 'PROG', ijmax, 1, 10
   WRITE (UNIT=n, FMT=fmt1) stgrnd, 'PROG', ijmax, 1, 10
   WRITE (UNIT=n, FMT=fmt1) wt1soi, 'PROG', ijmax, 1, 0
   WRITE (UNIT=n, FMT=fmt1) wt2soi, 'PROG', ijmax, 1, 0
   WRITE (UNIT=n, FMT=fmt1) wt3soi, 'PROG', ijmax, 1, 0
   WRITE (UNIT=n, FMT=fmt1) temp2m, 'PROG', ijmax, 1, 40
   WRITE (UNIT=n, FMT=fmt1) umes2m, 'PROG', ijmax, 1, 0
   WRITE (UNIT=n, FMT=fmt1) uve10m, 'PROG', ijmax, 1, 60
   WRITE (UNIT=n, FMT=fmt1) vve10m, 'PROG', ijmax, 1, 60
   WRITE (UNIT=n, FMT=fmt1) veget0, 'PROG', ijmax, 1, 0
   WRITE (UNIT=n, FMT=fmt1) soiltx, 'PROG', ijmax, 1, 0
   WRITE (UNIT=n, FMT=fmt1) doair0, 'PROG', ijmax, 1, 201
   !JGERD -> fields necessary to data assimilation
   WRITE (UNIT=n, FMT=fmt1) sfcstg, 'PROG', ijmax, 1, 0
   WRITE (UNIT=n, FMT=fmt1) vegcov, 'PROG', ijmax, 1, 0
   WRITE (UNIT=n, FMT=fmt1) snodpt, 'PROG', ijmax, 1, 0
   WRITE (UNIT=n, FMT=fmt1) lqwcon, 'PROG', ijmax, kmax, 0
   WRITE (UNIT=n, FMT=fmt1) icwcon, 'PROG', ijmax, kmax, 0
   WRITE (UNIT=n, FMT=fmt1) concld, 'PROG', ijmax, kmax, 0

   IF (ifday == 0 .and. tod == 0.0_r8 ) RETURN

   DO m=1,mxavl
      IF (dodia(m) .and. (iavrq(m) > 0)) THEN
         nn=iavrq(m)
         WRITE(UNIT=n,FMT=fmt1)reqdg(nn),diag,isg(itavl(m)),lvrq(nn),nurq(nn)
         IF(nfctrl(91) >= 1)WRITE(UNIT=nfprt,FMT=161) reqdg(nn),diag, &
              isg(itavl(m)),lvrq(nn),nurq(nn)
      END IF
   END DO

   IF(icf.ne.0)THEN
      DO ix=1,icf
         WRITE(UNIT=n,FMT=fmt1)combf(ix),diag,isg(itcf(ix)),lvcf(ix),nucf(ix)
         IF(nfctrl(91) >= 1)WRITE(UNIT=nfprt,FMT=161) combf(ix),diag, &
              isg(itcf(ix)),lvcf(ix),nucf(ix)
      END DO
   END IF
161 FORMAT(' ',A40,2X,A4,2X,I8,3X,I4,4X,I3)    
 END SUBROUTINE WriteDir4
  SUBROUTINE WriteDir8(n    , idate  , idatec, &
                       a_hybr, b_hybr , tod    , ifday, &
                       reqdg, combf  , dodia, &
                       itcf , itavl  , iavrq, &
                       nucf , lvrq   , nurq , &
                       lvcf , mxavl  , icf  , &
                       opt)
    INTEGER          , INTENT(IN ) :: n
    INTEGER          , INTENT(IN ) :: idate (4)
    INTEGER          , INTENT(IN ) :: idatec(4)
    REAL(KIND=r8)    , INTENT(IN ) :: a_hybr (kMax+1) 
    REAL(KIND=r8)    , INTENT(IN ) :: b_hybr (kMax+1) 
    REAL(KIND=r8)    , INTENT(IN ) :: tod
    INTEGER          , INTENT(IN ) :: ifday
    CHARACTER(len=40), INTENT(IN ) :: reqdg (1:ndrq)
    CHARACTER(len=40), INTENT(IN ) :: combf (1:ncdg)
    LOGICAL          , INTENT(IN ) :: dodia (1:ndavl)
    INTEGER          , INTENT(IN ) :: itcf  (1:ncdg)
    INTEGER          , INTENT(IN ) :: itavl(1:ndavl)
    INTEGER          , INTENT(IN ) :: iavrq(1:ndavl)
    INTEGER          , INTENT(IN ) :: nucf (1:ncdg )
    INTEGER          , INTENT(IN ) :: lvrq (1:ndrq )
    INTEGER          , INTENT(IN ) :: nurq (1:ndrq )
    INTEGER          , INTENT(IN ) :: lvcf (1:ncdg )
    INTEGER          , INTENT(IN ) :: mxavl
    INTEGER          , INTENT(IN ) :: icf
    LOGICAL          , INTENT(IN ) :: opt
    INTEGER               :: isg(2)

    CHARACTER (LEN= 5) :: imdl
    CHARACTER (LEN=150) :: jttl
    CHARACTER (LEN=21), PARAMETER :: ittl='CPTEC HIBRID VERS 2.0'
    CHARACTER (LEN= 4), PARAMETER :: nexp='0001'
    CHARACTER (LEN=40), PARAMETER :: orogra = 'TOPOGRAPHY' 
    CHARACTER (LEN=40), PARAMETER :: LSImsk = 'LAND SEA ICE MASK'
    CHARACTER (LEN=40), PARAMETER :: lnsurf = 'LN SURFACE PRESSURE'
    CHARACTER (LEN=40), PARAMETER :: divrgn = 'DIVERGENCE'
    CHARACTER (LEN=40), PARAMETER :: vortic = 'VORTICITY'
    CHARACTER (LEN=40), PARAMETER :: spechu = 'SPECIFIC HUMIDITY'
    CHARACTER (LEN=40), PARAMETER :: tempvi = 'VIRTUAL TEMPERATURE'
    CHARACTER (LEN=40), PARAMETER :: surfte = 'SURFACE TEMPERATURE'
    CHARACTER (LEN=40), PARAMETER :: srfrou = 'ROUGHNESS LENGTH'
    CHARACTER (LEN=40), PARAMETER :: deepte = 'DEEP SOIL TEMPERATURE'
    CHARACTER (LEN=40), PARAMETER :: stcnpy = 'STORAGE ON CANOPY'
    CHARACTER (LEN=40), PARAMETER :: stgrnd = 'STORAGE ON GROUND'
    CHARACTER (LEN=40), PARAMETER :: wt1soi = 'SOIL WETNESS OF SURFACE'
    CHARACTER (LEN=40), PARAMETER :: wt2soi = 'SOIL WETNESS OF ROOT ZONE'
    CHARACTER (LEN=40), PARAMETER :: wt3soi = 'SOIL WETNESS OF DRAINAGE ZONE' 
    CHARACTER (LEN=40), PARAMETER :: temp2m = 'TEMPERATURE AT 2-M FROM SURFACE' 
    CHARACTER (LEN=40), PARAMETER :: umes2m = 'SPECIFIC HUMIDITY AT 2-M FROM SURFACE' 
    CHARACTER (LEN=40), PARAMETER :: uve10m = 'ZONAL WIND AT 10-M FROM SURFACE'
    CHARACTER (LEN=40), PARAMETER :: vve10m = 'MERID WIND AT 10-M FROM SURFACE'
    CHARACTER (LEN=40), PARAMETER :: veget0 = 'MASK VEGETATION'    
    CHARACTER (LEN=40), PARAMETER :: soiltx = 'MASK SOIL TEXTURE CLASSES'
    CHARACTER (LEN=40), PARAMETER :: doair0 = 'PARTIAL OXYGEN DENSITY '
    !JGERD -> fields necessary to data assimilation
    CHARACTER (LEN=40), PARAMETER :: sfcstg = 'SURFACE SOIL TEMPERATURE'
    CHARACTER (LEN=40), PARAMETER :: vegcov = 'VEGETATION COVER'
    CHARACTER (LEN=40), PARAMETER :: snodpt = 'SNOW DEPTH'
    CHARACTER (LEN=40), PARAMETER :: lqwcon = 'LIQ MIXING RATIO PROGNOSTIC'
    CHARACTER (LEN=40), PARAMETER :: icwcon = 'ICE MIXING RATIO PROGNOSTIC'
    CHARACTER (LEN=40), PARAMETER :: concld = 'CLOUD TOTAL PROGNOSTIC'
    CHARACTER (LEN=29), PARAMETER :: fmt1='(A40,2X,A4,2X,I8,3X,I4,4X,I3)'
    CHARACTER (LEN=4 ), PARAMETER :: diag='DIAG'
    INTEGER :: m
    INTEGER :: nn
    INTEGER :: ix
    INTEGER :: ihr
    INTEGER :: iday
    INTEGER :: mon
    INTEGER :: iyr
    CHARACTER(LEN=10) :: labelc
    CHARACTER(LEN=10) :: labeln

    ihr  = idatec(1)
    iday = idatec(3)
    mon  = idatec(2)
    iyr  = idatec(4)
    IF(DHFCT < 0 .and. opt)THEN
      WRITE(labelc,'(I4.4,I2.2,I2.2,I2.2)')iyr,mon,iday,ihr
      labeln=CalDate(labelc,'-')
      READ(labeln(1: 4),'(I4.4)')iyr
      READ(labeln(5: 6),'(I2.2)')mon
      READ(labeln(7: 8),'(I2.2)')iday
      READ(labeln(9:10),'(I2.2)')ihr
    END IF
   
    isg(1)=iMax*jMax
    isg(2)=2*mnMax

    jttl='CPTEC AGCM REVIS 1.0 2000  T    L    COLD [OPTIONS |surface=01||pbl=04|]'
    WRITE (jttl(29:32), '(i4.4)') mMax-1
    WRITE (jttl(34:36), '(i2.2)') kMax
    WRITE (jttl(61:62), '(i2.2)') schemes
    WRITE (jttl(69:70), '(i2.2)') atmpbl
    WRITE (imdl, '(A1,I4.4)') 'T', mMax-1

    WRITE (UNIT=n, FMT='(A21)')   ittl
    WRITE (UNIT=n, FMT='(A4,1X,  A4,1X, A5,1X,11I5,1X,A4)') &
                         nexp, 'SEQU', imdl, mMax, kmax+1, kmax+1, &
         ihr, iday, mon, iyr, idate, 'TAPE'
    WRITE (UNIT=n, FMT='(A)')   jttl
    WRITE (UNIT=n, FMT='(5E16.8)')   a_hybr
    WRITE (UNIT=n, FMT='(5E16.8)')   b_hybr
    WRITE (UNIT=n, FMT=fmt1) orogra, 'FIXD', 2*mnMax, 1, 10
    WRITE (UNIT=n, FMT=fmt1) LSImsk, 'FIXD', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) lnsurf, 'PROG', 2*mnMax, 1, 142
    WRITE (UNIT=n, FMT=fmt1) divrgn, 'PROG', 2*mnMax, kmax, 50
    WRITE (UNIT=n, FMT=fmt1) vortic, 'PROG', 2*mnMax, kmax, 50
    WRITE (UNIT=n, FMT=fmt1) spechu, 'PROG', 2*mnMax, kmax, 0
    WRITE (UNIT=n, FMT=fmt1) tempvi, 'PROG', 2*mnMax, kmax, 40
    WRITE (UNIT=n, FMT=fmt1) srfrou, 'PROG', ijmax, 1, 10
    WRITE (UNIT=n, FMT=fmt1) surfte, 'PROG', ijmax, 1, 40
    WRITE (UNIT=n, FMT=fmt1) deepte, 'PROG', ijmax, 1, 40
    WRITE (UNIT=n, FMT=fmt1) stcnpy, 'PROG', ijmax, 1, 10
    WRITE (UNIT=n, FMT=fmt1) stgrnd, 'PROG', ijmax, 1, 10
    WRITE (UNIT=n, FMT=fmt1) wt1soi, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) wt2soi, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) wt3soi, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) temp2m, 'PROG', ijmax, 1, 40
    WRITE (UNIT=n, FMT=fmt1) umes2m, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) uve10m, 'PROG', ijmax, 1, 60
    WRITE (UNIT=n, FMT=fmt1) vve10m, 'PROG', ijmax, 1, 60
    WRITE (UNIT=n, FMT=fmt1) veget0, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) soiltx, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) doair0, 'PROG', ijmax, 1, 201
    !JGERD -> fields necessary to data assimilation
    WRITE (UNIT=n, FMT=fmt1) sfcstg, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) vegcov, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) snodpt, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) lqwcon, 'PROG', ijmax, kmax, 0
    WRITE (UNIT=n, FMT=fmt1) icwcon, 'PROG', ijmax, kmax, 0
    WRITE (UNIT=n, FMT=fmt1) concld, 'PROG', ijmax, kmax, 0

    IF (ifday == 0 .and. tod == 0.0_r8 ) RETURN

    DO m=1,mxavl
       IF (dodia(m) .and. (iavrq(m) > 0)) THEN
          nn=iavrq(m)
          WRITE(UNIT=n,FMT=fmt1)reqdg(nn),diag,isg(itavl(m)),lvrq(nn),nurq(nn)
          IF(nfctrl(91) >= 1)WRITE(UNIT=nfprt,FMT=161) reqdg(nn),diag, &
               isg(itavl(m)),lvrq(nn),nurq(nn)
       END IF
    END DO

    IF(icf.ne.0)THEN
       DO ix=1,icf
          WRITE(UNIT=n,FMT=fmt1)combf(ix),diag,isg(itcf(ix)),lvcf(ix),nucf(ix)
          IF(nfctrl(91) >= 1)WRITE(UNIT=nfprt,FMT=161) combf(ix),diag, &
               isg(itcf(ix)),lvcf(ix),nucf(ix)
       END DO
    END IF

161 FORMAT(' ',A40,2X,A4,2X,I8,3X,I4,4X,I3)
  END SUBROUTINE WriteDir8
!---------------

  SUBROUTINE WriteDirDhn4(n, idate, ihr, iday, mon, iyr,tod,a_hybr,b_hybr )
    INTEGER(KIND=i4), INTENT(IN ) :: n
    INTEGER(KIND=i4), INTENT(IN ) :: idate(4)
    INTEGER(KIND=i4), INTENT(IN ) :: ihr
    INTEGER(KIND=i4), INTENT(IN ) :: iday
    INTEGER(KIND=i4), INTENT(IN ) :: mon
    INTEGER(KIND=i4), INTENT(IN ) :: iyr
    REAL    (KIND=r4), INTENT(IN ) :: tod
    REAL(KIND=r8)     , INTENT(IN   ) :: a_hybr(:)
    REAL(KIND=r8)     , INTENT(IN   ) :: b_hybr(:)
    REAL    (KIND=r4) :: r

    CHARACTER (LEN= 5) :: imdl
    CHARACTER (LEN=150) :: jttl
    CHARACTER (LEN=21), PARAMETER :: ittl='CPTEC HIBRID VERS 2.0'
    CHARACTER (LEN= 4), PARAMETER :: nexp='0001'
    CHARACTER (LEN=40), PARAMETER :: orogra = 'TOPOGRAPHY' 
    CHARACTER (LEN=40), PARAMETER :: lseamk = 'LAND SEA MASK'
    CHARACTER (LEN=40), PARAMETER :: lnsurf = 'LN SURFACE PRESSURE'
    CHARACTER (LEN=40), PARAMETER :: divrgn = 'DIVERGENCE'
    CHARACTER (LEN=40), PARAMETER :: vortic = 'VORTICITY'
    CHARACTER (LEN=40), PARAMETER :: spechu = 'SPECIFIC HUMIDITY'
    CHARACTER (LEN=40), PARAMETER :: tempvi = 'VIRTUAL TEMPERATURE'
    CHARACTER (LEN=40), PARAMETER :: surfte = 'SURFACE TEMPERATURE'
    CHARACTER (LEN=40), PARAMETER :: srfrou = 'ROUGHNESS LENGTH'
    CHARACTER (LEN=40), PARAMETER :: deepte = 'DEEP SOIL TEMPERATURE'
    CHARACTER (LEN=40), PARAMETER :: stcnpy = 'STORAGE ON CANOPY'
    CHARACTER (LEN=40), PARAMETER :: stgrnd = 'STORAGE ON GROUND'
    CHARACTER (LEN=40), PARAMETER :: wt1soi = 'SOIL WETNESS OF SURFACE'
    CHARACTER (LEN=40), PARAMETER :: wt2soi = 'SOIL WETNESS OF ROOT ZONE'
    CHARACTER (LEN=40), PARAMETER :: wt3soi = 'SOIL WETNESS OF DRAINAGE ZONE' 
    CHARACTER (LEN=40), PARAMETER :: temp2m = 'TEMPERATURE AT 2-M FROM SURFACE' 
    CHARACTER (LEN=40), PARAMETER :: umes2m = 'SPECIFIC HUMIDITY AT 2-M FROM SURFACE' 
    CHARACTER (LEN=40), PARAMETER :: uve10m = 'ZONAL WIND AT 10-M FROM SURFACE'
    CHARACTER (LEN=40), PARAMETER :: vve10m = 'MERID WIND AT 10-M FROM SURFACE'
    CHARACTER (LEN=40), PARAMETER :: veget0 = 'MASK VEGETATION'    
    CHARACTER (LEN=40), PARAMETER :: doair0 = 'PARTIAL OXYGEN DENSITY '
    CHARACTER (LEN=40), PARAMETER :: conliq = 'LIQ MIXING RATIO PROGNOSTIC'
    CHARACTER (LEN=40), PARAMETER :: conice = 'ICE MIXING RATIO PROGNOSTIC'
    CHARACTER (LEN=40), PARAMETER :: concld = 'CLOUD TOTAL PROGNOSTIC'
    CHARACTER (LEN=29), PARAMETER :: fmt1='(A40,2X,A4,2X,I8,3X,I4,4X,I3)'

    r = tod
    jttl='CPTEC AGCM REVIS 1.0 2000  T   L    COLD'
    WRITE (jttl(29:32), '(i4.4)') mMax-1
    WRITE (jttl(33:34), '(i2.2)') kMax
    WRITE (imdl, '(A1,I4.4)') 'T', mMax-1

    WRITE (UNIT=n, FMT='(A21)')   ittl
    WRITE (UNIT=n, FMT='(A4,1X,A4,1X,A5,1X,11I5,1X,A4)') &
         nexp, 'SEQU', imdl, mMax, kmax, kmax, &
         ihr, iday, mon, iyr, idate, 'TAPE'
    WRITE (UNIT=n, FMT='(A41)')   jttl
    WRITE (UNIT=n, FMT='(5E16.8)')   a_hybr
    WRITE (UNIT=n, FMT='(5E16.8)')   b_hybr
    WRITE (UNIT=n, FMT=fmt1) orogra, 'FIXD', 2*mnMax, 1, 10
    WRITE (UNIT=n, FMT=fmt1) lseamk, 'FIXD', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) lnsurf, 'PROG', 2*mnMax, 1, 142
    WRITE (UNIT=n, FMT=fmt1) divrgn, 'PROG', 2*mnMax, kmax, 50
    WRITE (UNIT=n, FMT=fmt1) vortic, 'PROG', 2*mnMax, kmax, 50
    WRITE (UNIT=n, FMT=fmt1) spechu, 'PROG', 2*mnMax, kmax, 0
    WRITE (UNIT=n, FMT=fmt1) tempvi, 'PROG', 2*mnMax, kmax, 40
    WRITE (UNIT=n, FMT=fmt1) srfrou, 'PROG', ijmax, 1, 10
    WRITE (UNIT=n, FMT=fmt1) surfte, 'PROG', ijmax, 1, 40
    WRITE (UNIT=n, FMT=fmt1) deepte, 'PROG', ijmax, 1, 40
    WRITE (UNIT=n, FMT=fmt1) stcnpy, 'PROG', ijmax, 1, 10
    WRITE (UNIT=n, FMT=fmt1) stgrnd, 'PROG', ijmax, 1, 10
    WRITE (UNIT=n, FMT=fmt1) wt1soi, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) wt2soi, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) wt3soi, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) temp2m, 'PROG', ijmax, 1, 40
    WRITE (UNIT=n, FMT=fmt1) umes2m, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) uve10m, 'PROG', ijmax, 1, 60
    WRITE (UNIT=n, FMT=fmt1) vve10m, 'PROG', ijmax, 1, 60
    WRITE (UNIT=n, FMT=fmt1) veget0, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) doair0, 'PROG', ijmax, 1, 201
    WRITE (UNIT=n, FMT=fmt1) conliq, 'PROG', ijmax, kmax, 0
    WRITE (UNIT=n, FMT=fmt1) conice, 'PROG', ijmax, kmax, 0
    WRITE (UNIT=n, FMT=fmt1) concld, 'PROG', ijmax, kmax, 0

!160 FORMAT(A40,2X,A4,2X,I8,3X,I4,4X,I3)
!161 FORMAT(' ',A40,2X,A4,2X,I8,3X,I4,4X,I3)    
  END SUBROUTINE WriteDirDhn4

  SUBROUTINE LogWriteDirDhn4(n, idate, ihr, iday, mon, iyr,tod,a_hybr,b_hybr ,dcol,iMaxNew,jMaxNew,filename,FIELDGAUS)
    INTEGER(KIND=i4), INTENT(IN ) :: n
    INTEGER(KIND=i4), INTENT(IN ) :: iMaxNew
    INTEGER(KIND=i4), INTENT(IN ) :: jMaxNew
    INTEGER(KIND=i4), INTENT(IN ) :: idate(4)
    INTEGER(KIND=i4), INTENT(IN ) :: ihr
    INTEGER(KIND=i4), INTENT(IN ) :: iday
    INTEGER(KIND=i4), INTENT(IN ) :: mon
    INTEGER(KIND=i4), INTENT(IN ) :: iyr
    REAL    (KIND=r4), INTENT(IN ) :: tod
    REAL(KIND=r8)     , INTENT(IN   ) :: a_hybr(:)
    REAL(KIND=r8)     , INTENT(IN   ) :: b_hybr(:)
    REAL    (KIND=r4), INTENT(IN ) :: dcol(jMaxNew)
    CHARACTER(LEN=*) , INTENT(IN ) :: filename
    LOGICAL          , INTENT(IN ) :: FIELDGAUS
    REAL      (KIND=r4) :: r
    CHARACTER (LEN= 4) :: imdl
    CHARACTER (LEN=40) :: jttl
    CHARACTER (LEN=21), PARAMETER :: ittl='CPTEC HIBRID VERS 2.0'
    CHARACTER (LEN= 4), PARAMETER :: nexp='0001'
    CHARACTER (LEN=40), PARAMETER :: uve10m = 'ZONAL WIND AT 10-M FROM SURFACE         '
    CHARACTER (LEN=40), PARAMETER :: vve10m = 'MERID WIND AT 10-M FROM SURFACE         '
    CHARACTER (LEN=40), PARAMETER :: wnd100 = 'SPEED WIND AT 100-M FROM SURFACE LAYER  '
    CHARACTER (LEN=40), PARAMETER :: dir100 = 'DIR WIND AT 100-M FROM SURFACE LAYER    '
    
    CHARACTER (LEN=29), PARAMETER :: fmt1='(A40,2X,A4,2X,I8,3X,I4,4X,I3)'
    CHARACTER (LEN=3), PARAMETER :: cmth(12)=(/&
         'JAN','FEB','MAR','APR','MAY','JUN', &
         'JUL','AUG','SEP','OCT','NOV','DEC'/)
    REAL(KIND=r8)     , PARAMETER :: Undef=-2.56E+33_r8
    INTEGER :: i,j,it,nof,k



    WRITE (UNIT=n, FMT='(A)') 'DSET ^'//TRIM(filename)
    WRITE (UNIT=n, FMT='(A)') '*'
    WRITE (UNIT=n, FMT='(A)') 'OPTIONS SEQUENTIAL YREV BIG_ENDIAN'
    WRITE (UNIT=n, FMT='(A)') '*'
    WRITE (UNIT=n, FMT='(A,1PE9.2)') 'UNDEF ', Undef
    WRITE (UNIT=n, FMT='(A)') '*'
    WRITE (UNIT=n, FMT='(3A)') 'TITLE ',ittl,nexp
    WRITE (UNIT=n, FMT='(A)') '*'
    WRITE (UNIT=n, FMT='(A,I5,A,F8.3,F15.10)')'XDEF ', iMaxNew, ' LINEAR ', 0.0_r8, 360.0_r8/REAL(iMaxNew,r8)
    WRITE (UNIT=n, FMT='(A,I5,A)') 'YDEF ', jMaxNew, ' LEVELS '
    WRITE (UNIT=n, FMT='(8F25.5)') (dcol(j),j=jMaxNew,1,-1)
    WRITE (UNIT=n, FMT='(A,I5,A,10F15.10)') 'ZDEF ', kMax, ' LEVELS ',(a_hybr(k)+100000.0_r8*b_hybr(k),k=kMax,1,-1)
    it=1
    WRITE (UNIT=n, FMT='(A,I5,A,I2.2,A,I2.2,A,I4,A)') &
            'TDEF ', it, ' LINEAR ', ihr, 'Z', iday, cmth(mon), iyr, ' 6HR'
    WRITE (UNIT=n, FMT='(A)') '*'
    nof=4
    WRITE (UNIT=n, FMT='(A,I5)') 'VARS ', nof
 
    WRITE (UNIT=n, FMT='(A)') 'UVE10M    0 99 '//TRIM(uve10m)//' (M/S             )'
    WRITE (UNIT=n, FMT='(A)') 'VVE10M    0 99 '//TRIM(uve10m)//' (M/S             )'
    WRITE (UNIT=n, FMT='(A)') 'WND100    0 99 '//TRIM(wnd100)//' (M/S             )'
    WRITE (UNIT=n, FMT='(A)') 'DIR100    0 99 '//TRIM(dir100)//' (DEGREE          )'
    WRITE (UNIT=n, FMT='(A)') 'ENDVARS'

!160 FORMAT(A40,2X,A4,2X,I8,3X,I4,4X,I3)
!161 FORMAT(' ',A40,2X,A4,2X,I8,3X,I4,4X,I3)    
  END SUBROUTINE LogWriteDirDhn4


  SUBROUTINE WriteDirDhn8(n, idate, ihr, iday, mon, iyr,tod,a_hybr,b_hybr)
    INTEGER , INTENT(IN ) :: n
    INTEGER , INTENT(IN ) :: idate(4)
    INTEGER , INTENT(IN ) :: ihr
    INTEGER , INTENT(IN ) :: iday
    INTEGER , INTENT(IN ) :: mon
    INTEGER , INTENT(IN ) :: iyr
    REAL    (KIND=r8), INTENT(IN ) :: tod
    REAL(KIND=r8)     , INTENT(IN   ) :: a_hybr(:)
    REAL(KIND=r8)     , INTENT(IN   ) :: b_hybr(:)
    REAL    (KIND=r8) :: r

    CHARACTER (LEN= 5) :: imdl
    CHARACTER (LEN=41) :: jttl
    CHARACTER (LEN=21), PARAMETER :: ittl='CPTEC HIBRID VERS 2.0'
    CHARACTER (LEN= 4), PARAMETER :: nexp='0001'
    CHARACTER (LEN=40), PARAMETER :: orogra = 'TOPOGRAPHY' 
    CHARACTER (LEN=40), PARAMETER :: lseamk = 'LAND SEA MASK'
    CHARACTER (LEN=40), PARAMETER :: lnsurf = 'LN SURFACE PRESSURE'
    CHARACTER (LEN=40), PARAMETER :: divrgn = 'DIVERGENCE'
    CHARACTER (LEN=40), PARAMETER :: vortic = 'VORTICITY'
    CHARACTER (LEN=40), PARAMETER :: spechu = 'SPECIFIC HUMIDITY'
    CHARACTER (LEN=40), PARAMETER :: tempvi = 'VIRTUAL TEMPERATURE'
    CHARACTER (LEN=40), PARAMETER :: surfte = 'SURFACE TEMPERATURE'
    CHARACTER (LEN=40), PARAMETER :: srfrou = 'ROUGHNESS LENGTH'
    CHARACTER (LEN=40), PARAMETER :: deepte = 'DEEP SOIL TEMPERATURE'
    CHARACTER (LEN=40), PARAMETER :: stcnpy = 'STORAGE ON CANOPY'
    CHARACTER (LEN=40), PARAMETER :: stgrnd = 'STORAGE ON GROUND'
    CHARACTER (LEN=40), PARAMETER :: wt1soi = 'SOIL WETNESS OF SURFACE'
    CHARACTER (LEN=40), PARAMETER :: wt2soi = 'SOIL WETNESS OF ROOT ZONE'
    CHARACTER (LEN=40), PARAMETER :: wt3soi = 'SOIL WETNESS OF DRAINAGE ZONE'
    CHARACTER (LEN=40), PARAMETER :: temp2m = 'TEMPERATURE AT 2-M FROM SURFACE' 
    CHARACTER (LEN=40), PARAMETER :: umes2m = 'SPECIFIC HUMIDITY AT 2-M FROM SURFACE'     
    CHARACTER (LEN=40), PARAMETER :: uve10m = 'ZONAL WIND AT 10-M FROM SURFACE'
    CHARACTER (LEN=40), PARAMETER :: vve10m = 'MERID WIND AT 10-M FROM SURFACE'
    CHARACTER (LEN=40), PARAMETER :: veget0 = 'MASK VEGETATION'    
    CHARACTER (LEN=40), PARAMETER :: doair0 = 'PARTIAL OXYGEN DENSITY '
    CHARACTER (LEN=40), PARAMETER :: conliq = 'LIQ MIXING RATIO PROGNOSTIC'
    CHARACTER (LEN=40), PARAMETER :: conice = 'ICE MIXING RATIO PROGNOSTIC'
    CHARACTER (LEN=40), PARAMETER :: concld = 'CLOUD TOTAL PROGNOSTIC'
    CHARACTER (LEN=29), PARAMETER :: fmt1='(A40,2X,A4,2X,I8,3X,I4,4X,I3)'



    r = tod
    jttl='CPTEC AGCM REVIS 1.0 2000  T   L    COLD'
    WRITE (jttl(29:32), '(i4.4)') mMax-1
    WRITE (jttl(33:34), '(i2.2)') kMax
    WRITE (imdl, '(A1,I4.4)') 'T', mMax-1

    WRITE (UNIT=n, FMT='(A21)')   ittl
    WRITE (UNIT=n, FMT='(A4,1X,A4,1X,A5,1X,11I5,1X,A4)') &
         nexp, 'SEQU', imdl, mMax, kmax, kmax, &
         ihr, iday, mon, iyr, idate, 'TAPE'
    WRITE (UNIT=n, FMT='(A41)')   jttl
    WRITE (UNIT=n, FMT='(5E16.8)')   a_hybr
    WRITE (UNIT=n, FMT='(5E16.8)')   b_hybr
    WRITE (UNIT=n, FMT=fmt1) orogra, 'FIXD', 2*mnMax, 1, 10
    WRITE (UNIT=n, FMT=fmt1) lseamk, 'FIXD', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) lnsurf, 'PROG', 2*mnMax, 1, 142
    WRITE (UNIT=n, FMT=fmt1) divrgn, 'PROG', 2*mnMax, kmax, 50
    WRITE (UNIT=n, FMT=fmt1) vortic, 'PROG', 2*mnMax, kmax, 50
    WRITE (UNIT=n, FMT=fmt1) spechu, 'PROG', 2*mnMax, kmax, 0
    WRITE (UNIT=n, FMT=fmt1) tempvi, 'PROG', 2*mnMax, kmax, 40
    WRITE (UNIT=n, FMT=fmt1) srfrou, 'PROG', ijmax, 1, 10
    WRITE (UNIT=n, FMT=fmt1) surfte, 'PROG', ijmax, 1, 40
    WRITE (UNIT=n, FMT=fmt1) deepte, 'PROG', ijmax, 1, 40
    WRITE (UNIT=n, FMT=fmt1) stcnpy, 'PROG', ijmax, 1, 10
    WRITE (UNIT=n, FMT=fmt1) stgrnd, 'PROG', ijmax, 1, 10
    WRITE (UNIT=n, FMT=fmt1) wt1soi, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) wt2soi, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) wt3soi, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) temp2m, 'PROG', ijmax, 1, 40
    WRITE (UNIT=n, FMT=fmt1) umes2m, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) uve10m, 'PROG', ijmax, 1, 60
    WRITE (UNIT=n, FMT=fmt1) vve10m, 'PROG', ijmax, 1, 60
    WRITE (UNIT=n, FMT=fmt1) veget0, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) doair0, 'PROG', ijmax, 1, 201
    WRITE (UNIT=n, FMT=fmt1) conliq, 'PROG', ijmax, kmax, 0
    WRITE (UNIT=n, FMT=fmt1) conice, 'PROG', ijmax, kmax, 0
    WRITE (UNIT=n, FMT=fmt1) concld, 'PROG', ijmax, kmax, 0

  END SUBROUTINE WriteDirDhn8

  SUBROUTINE LogWriteDirDhn8(n, idate, ihr, iday, mon, iyr,tod,a_hybr,b_hybr,dcol,iMaxNew,jMaxNew,filename,FIELDGAUS)
    INTEGER, INTENT(IN ) :: n
    INTEGER, INTENT(IN ) :: iMaxNew
    INTEGER, INTENT(IN ) :: jMaxNew
    INTEGER, INTENT(IN ) :: idate(4)
    INTEGER, INTENT(IN ) :: ihr
    INTEGER, INTENT(IN ) :: iday
    INTEGER, INTENT(IN ) :: mon
    INTEGER, INTENT(IN ) :: iyr
    REAL   (KIND=r8), INTENT(IN ) :: tod
    REAL(KIND=r8)     , INTENT(IN   ) :: a_hybr(:)
    REAL(KIND=r8)     , INTENT(IN   ) :: b_hybr(:)
    REAL   (KIND=r8), INTENT(IN ) :: dcol(jMaxNew)
    CHARACTER(LEN=*), INTENT(IN ) :: filename
    LOGICAL         , INTENT(IN ) :: FIELDGAUS
    REAL      (KIND=r4) :: r
    CHARACTER (LEN= 4) :: imdl
    CHARACTER (LEN=40) :: jttl
    CHARACTER (LEN=21), PARAMETER :: ittl='CPTEC HIBRID VERS 2.0'
    CHARACTER (LEN= 4), PARAMETER :: nexp='0001'
    CHARACTER (LEN=40), PARAMETER :: uve10m = 'ZONAL WIND AT 10-M FROM SURFACE         '
    CHARACTER (LEN=40), PARAMETER :: vve10m = 'MERID WIND AT 10-M FROM SURFACE         '
    CHARACTER (LEN=40), PARAMETER :: wnd100 = 'SPEED WIND AT 100-M FROM SURFACE LAYER  '
    CHARACTER (LEN=40), PARAMETER :: dir100 = 'DIR WIND AT 100-M FROM SURFACE LAYER    '

    CHARACTER (LEN=29), PARAMETER :: fmt1='(A40,2X,A4,2X,I8,3X,I4,4X,I3)'
    CHARACTER (LEN=3), PARAMETER :: cmth(12)=(/&
         'JAN','FEB','MAR','APR','MAY','JUN', &
         'JUL','AUG','SEP','OCT','NOV','DEC'/)
    REAL(KIND=r8)     , PARAMETER :: Undef=-2.56E+33_r8
    INTEGER :: i,j,it,nof,k



    WRITE (UNIT=n, FMT='(A)') 'DSET ^'//TRIM(filename)
    WRITE (UNIT=n, FMT='(A)') '*'
    WRITE (UNIT=n, FMT='(A)') 'OPTIONS SEQUENTIAL YREV BIG_ENDIAN'
    WRITE (UNIT=n, FMT='(A)') '*'
    WRITE (UNIT=n, FMT='(A,1PE9.2)') 'UNDEF ', Undef
    WRITE (UNIT=n, FMT='(A)') '*'
    WRITE (UNIT=n, FMT='(3A)') 'TITLE ',ittl,nexp
    WRITE (UNIT=n, FMT='(A)') '*'
    WRITE (UNIT=n, FMT='(A,I5,A,F8.3,F15.10)')'XDEF ', iMaxNew, ' LINEAR ', 0.0_r8, 360.0_r8/REAL(iMaxNew,r8)
    WRITE (UNIT=n, FMT='(A,I5,A)') 'YDEF ', jMaxNew, ' LEVELS '
    WRITE (UNIT=n, FMT='(8F15.10)') (dcol(j),j=jMaxNew,1,-1)
    WRITE (UNIT=n, FMT='(A,I5,A)') 'ZDEF ', kMax, ' LEVELS '
    WRITE (UNIT=n, FMT='(A,I5,A,10F15.10)') 'ZDEF ', kMax, ' LEVELS ',(a_hybr(k)+100000.0_r8*b_hybr(k),k=kMax,1,-1)
    it=1
    WRITE (UNIT=n, FMT='(A,I5,A,I2.2,A,I2.2,A,I4,A)') &
            'TDEF ', it, ' LINEAR ', ihr, 'Z', iday, cmth(mon), iyr, ' 6HR'
    WRITE (UNIT=n, FMT='(A)') '*'
    nof=4
    WRITE (UNIT=n, FMT='(A,I5)') 'VARS ', nof
   
    WRITE (UNIT=n, FMT='(A)') 'UVE10M    0 99 '//TRIM(uve10m)//' (M/S             )'
    WRITE (UNIT=n, FMT='(A)') 'VVE10M    0 99 '//TRIM(uve10m)//' (M/S             )'
    WRITE (UNIT=n, FMT='(A)') 'WND100    0 99 '//TRIM(wnd100)//' (M/S             )'
    WRITE (UNIT=n, FMT='(A)') 'DIR100    0 99 '//TRIM(dir100)//' (DEGREE          )'
    WRITE (UNIT=n, FMT='(A)') 'ENDVARS'

!160 FORMAT(A40,2X,A4,2X,I8,3X,I4,4X,I3)
!161 FORMAT(' ',A40,2X,A4,2X,I8,3X,I4,4X,I3)    
  END SUBROUTINE LogWriteDirDhn8



!---------------
  SUBROUTINE WriteDire4(n, idate, ihr, iday, mon, iyr, a_hybr, b_hybr,tod)
    INTEGER(KIND=i4), INTENT(IN ) :: n
    INTEGER(KIND=i4), INTENT(IN ) :: idate(4)
    INTEGER(KIND=i4), INTENT(IN ) :: ihr
    INTEGER(KIND=i4), INTENT(IN ) :: iday
    INTEGER(KIND=i4), INTENT(IN ) :: mon
    INTEGER(KIND=i4), INTENT(IN ) :: iyr
    REAL    (KIND=r4), INTENT(IN ) :: a_hybr(kMax+1)
    REAL    (KIND=r4), INTENT(IN ) :: b_hybr(kMax+1)
    REAL    (KIND=r4), INTENT(IN ) :: tod
    REAL    (KIND=r4) :: r

    CHARACTER (LEN= 5) :: imdl
    CHARACTER (LEN=150) :: jttl
    CHARACTER (LEN=21), PARAMETER :: ittl='CPTEC HIBRID VERS 2.0'
    CHARACTER (LEN= 4), PARAMETER :: nexp='0001'
    CHARACTER (LEN=40), PARAMETER :: orogra = 'TOPOGRAPHY' 
    CHARACTER (LEN=40), PARAMETER :: LSImsk = 'LAND SEA ICE MASK'
    CHARACTER (LEN=40), PARAMETER :: lnsurf = 'LN SURFACE PRESSURE'
    CHARACTER (LEN=40), PARAMETER :: divrgn = 'DIVERGENCE'
    CHARACTER (LEN=40), PARAMETER :: vortic = 'VORTICITY'
    CHARACTER (LEN=40), PARAMETER :: spechu = 'SPECIFIC HUMIDITY'
    CHARACTER (LEN=40), PARAMETER :: tempvi = 'VIRTUAL TEMPERATURE'
    CHARACTER (LEN=40), PARAMETER :: surfte = 'SURFACE TEMPERATURE'
    CHARACTER (LEN=40), PARAMETER :: srfrou = 'ROUGHNESS LENGTH'
    CHARACTER (LEN=40), PARAMETER :: deepte = 'DEEP SOIL TEMPERATURE'
    CHARACTER (LEN=40), PARAMETER :: stcnpy = 'STORAGE ON CANOPY'
    CHARACTER (LEN=40), PARAMETER :: stgrnd = 'STORAGE ON GROUND'
    CHARACTER (LEN=40), PARAMETER :: wt1soi = 'SOIL WETNESS OF SURFACE'
    CHARACTER (LEN=40), PARAMETER :: wt2soi = 'SOIL WETNESS OF ROOT ZONE'
    CHARACTER (LEN=40), PARAMETER :: wt3soi = 'SOIL WETNESS OF DRAINAGE ZONE' 
    CHARACTER (LEN=40), PARAMETER :: temp2m = 'TEMPERATURE AT 2-M FROM SURFACE' 
    CHARACTER (LEN=40), PARAMETER :: umes2m = 'SPECIFIC HUMIDITY AT 2-M FROM SURFACE' 
    CHARACTER (LEN=40), PARAMETER :: uve10m = 'ZONAL WIND AT 10-M FROM SURFACE'
    CHARACTER (LEN=40), PARAMETER :: vve10m = 'MERID WIND AT 10-M FROM SURFACE'
    CHARACTER (LEN=40), PARAMETER :: veget0 = 'MASK VEGETATION'    
    CHARACTER (LEN=40), PARAMETER :: soiltx = 'MASK SOIL TEXTURE CLASSES'
    CHARACTER (LEN=40), PARAMETER :: doair0 = 'PARTIAL OXYGEN DENSITY '
    !JGERD -> fields necessary to data assimilation
    CHARACTER (LEN=40), PARAMETER :: sfcstg = 'SURFACE SOIL TEMPERATURE'
    CHARACTER (LEN=40), PARAMETER :: vegcov = 'VEGETATION COVER'
    CHARACTER (LEN=40), PARAMETER :: snodpt = 'SNOW DEPTH'
    CHARACTER (LEN=40), PARAMETER :: icwcon = 'ICE MIXING RATIO PROGNOSTIC'
    CHARACTER (LEN=40), PARAMETER :: lqwcon = 'LIQ MIXING RATIO PROGNOSTIC'

    CHARACTER (LEN=40), PARAMETER :: concld = 'CLOUD TOTAL PROGNOSTIC'

    CHARACTER (LEN=29), PARAMETER :: fmt1='(A40,2X,A4,2X,I8,3X,I4,4X,I3)'

    r = tod
    jttl='CPTEC AGCM REVIS 1.0 2000  T    L    COLD [OPTIONS |surface=01||pbl=04|]'
    WRITE (jttl(29:32), '(i4.4)') mMax-1
    WRITE (jttl(34:36), '(i2.2)') kMax
    WRITE (jttl(61:62), '(i2.2)') schemes
    WRITE (jttl(69:70), '(i2.2)') atmpbl
    WRITE (imdl, '(A1,I4.4)') 'T', mMax-1

    WRITE (UNIT=n, FMT='(A21)')   ittl
    WRITE (UNIT=n, FMT='(A4,1X,A4,1X,A5,1X,11I5,1X,A4)') &
         nexp, 'SEQU', imdl, mMax, kmax, kmax, &
         ihr, iday, mon, iyr, idate, 'TAPE'
    WRITE (UNIT=n, FMT='(A)')   jttl
    WRITE (UNIT=n, FMT='(5E16.8)')   a_hybr
    WRITE (UNIT=n, FMT='(5E16.8)')   b_hybr
    WRITE (UNIT=n, FMT=fmt1) orogra, 'FIXD', 2*mnMax, 1, 10
    WRITE (UNIT=n, FMT=fmt1) LSImsk, 'FIXD', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) lnsurf, 'PROG', 2*mnMax, 1, 142
    WRITE (UNIT=n, FMT=fmt1) divrgn, 'PROG', 2*mnMax, kmax, 50
    WRITE (UNIT=n, FMT=fmt1) vortic, 'PROG', 2*mnMax, kmax, 50
    WRITE (UNIT=n, FMT=fmt1) spechu, 'PROG', 2*mnMax, kmax, 0
    WRITE (UNIT=n, FMT=fmt1) tempvi, 'PROG', 2*mnMax, kmax, 40
    WRITE (UNIT=n, FMT=fmt1) srfrou, 'PROG', ijmax, 1, 10
    WRITE (UNIT=n, FMT=fmt1) surfte, 'PROG', ijmax, 1, 40
    WRITE (UNIT=n, FMT=fmt1) deepte, 'PROG', ijmax, 1, 40
    WRITE (UNIT=n, FMT=fmt1) stcnpy, 'PROG', ijmax, 1, 10
    WRITE (UNIT=n, FMT=fmt1) stgrnd, 'PROG', ijmax, 1, 10
    WRITE (UNIT=n, FMT=fmt1) wt1soi, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) wt2soi, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) wt3soi, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) temp2m, 'PROG', ijmax, 1, 40
    WRITE (UNIT=n, FMT=fmt1) umes2m, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) uve10m, 'PROG', ijmax, 1, 60
    WRITE (UNIT=n, FMT=fmt1) vve10m, 'PROG', ijmax, 1, 60
    WRITE (UNIT=n, FMT=fmt1) veget0, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) soiltx, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) doair0, 'PROG', ijmax, 1, 201
    !JGERD -> fields necessary to data assimilation
    WRITE (UNIT=n, FMT=fmt1) sfcstg, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) vegcov, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) snodpt, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) icwcon, 'PROG', ijmax, kmax, 0
    WRITE (UNIT=n, FMT=fmt1) lqwcon, 'PROG', ijmax, kmax, 0
    WRITE (UNIT=n, FMT=fmt1) concld, 'PROG', ijmax, kmax, 0

!160 FORMAT(A40,2X,A4,2X,I8,3X,I4,4X,I3)
!161 FORMAT(' ',A40,2X,A4,2X,I8,3X,I4,4X,I3)    
  END SUBROUTINE WriteDire4

  SUBROUTINE LogWriteDire4(n, idate, ihr, iday, mon, iyr, a_h, b_h,  tod,dcol,iMaxNew,jMaxNew,filename,FIELDGAUS)
    INTEGER(KIND=i4), INTENT(IN ) :: n
    INTEGER(KIND=i4), INTENT(IN ) :: iMaxNew
    INTEGER(KIND=i4), INTENT(IN ) :: jMaxNew
    INTEGER(KIND=i4), INTENT(IN ) :: idate(4)
    INTEGER(KIND=i4), INTENT(IN ) :: ihr
    INTEGER(KIND=i4), INTENT(IN ) :: iday
    INTEGER(KIND=i4), INTENT(IN ) :: mon
    INTEGER(KIND=i4), INTENT(IN ) :: iyr
    REAL    (KIND=r4), INTENT(IN ) :: a_h(kMax+1)
    REAL    (KIND=r4), INTENT(IN ) :: b_h(kMax+1)
    REAL    (KIND=r4), INTENT(IN ) :: tod

    REAL    (KIND=r4), INTENT(IN ) :: dcol(jMaxNew)
    CHARACTER(LEN=*) , INTENT(IN ) :: filename
    LOGICAL          , INTENT(IN ) :: FIELDGAUS
    REAL      (KIND=r4) :: r
    CHARACTER (LEN= 4) :: imdl
    CHARACTER (LEN=150) :: jttl
    CHARACTER (LEN=21), PARAMETER :: ittl='CPTEC HIBRID VERS 2.0'
    CHARACTER (LEN= 4), PARAMETER :: nexp='0001'
    CHARACTER (LEN=40), PARAMETER :: uve10m = 'ZONAL WIND AT 10-M FROM SURFACE'
    CHARACTER (LEN=40), PARAMETER :: vve10m = 'MERID WIND AT 10-M FROM SURFACE'
    CHARACTER (LEN=40), PARAMETER :: topo   = 'GEOPOTENTIAL                   '
    CHARACTER (LEN=40), PARAMETER :: Psnm   = 'LN SURFACE PRESSURE            '
    CHARACTER (LEN=40), PARAMETER :: uvel   = 'ZONAL WIND (U)                 '
    CHARACTER (LEN=40), PARAMETER :: vvel   = 'MERIDIONAL WIND (V)            '
    CHARACTER (LEN=40), PARAMETER :: omeg   = 'VERTICAL WIND (W)              '
    CHARACTER (LEN=40), PARAMETER :: temp   = 'ABSOLUTE TEMPERATURE           '
    CHARACTER (LEN=40), PARAMETER :: umes   = 'SPECIFIC HUMIDITY              '
    CHARACTER (LEN=40), PARAMETER :: dump   = 'DUMP FIELD                     '

    CHARACTER (LEN=29), PARAMETER :: fmt1='(A40,2X,A4,2X,I8,3X,I4,4X,I3)'
    CHARACTER (LEN=3), PARAMETER :: cmth(12)=(/&
         'JAN','FEB','MAR','APR','MAY','JUN', &
         'JUL','AUG','SEP','OCT','NOV','DEC'/)
    REAL(KIND=r8)     , PARAMETER :: Undef=-2.56E+33_r8
    INTEGER :: i,j,it,nof,k



    WRITE (UNIT=n, FMT='(A)') 'DSET ^'//TRIM(filename)
    WRITE (UNIT=n, FMT='(A)') '*'
    WRITE (UNIT=n, FMT='(A)') 'OPTIONS SEQUENTIAL YREV BIG_ENDIAN'
    WRITE (UNIT=n, FMT='(A)') '*'
    WRITE (UNIT=n, FMT='(A,1PE9.2)') 'UNDEF ', Undef
    WRITE (UNIT=n, FMT='(A)') '*'
    WRITE (UNIT=n, FMT='(3A)') 'TITLE ',ittl,nexp
    WRITE (UNIT=n, FMT='(A)') '*'
    WRITE (UNIT=n, FMT='(A,I5,A,F8.3,F15.10)')'XDEF ', iMaxNew, ' LINEAR ', 0.0_r8, 360.0_r8/REAL(iMaxNew,r8)
    WRITE (UNIT=n, FMT='(A,I5,A)') 'YDEF ', jMaxNew, ' LEVELS '
    WRITE (UNIT=n, FMT='(8F25.5)') (dcol(j),j=jMaxNew,1,-1)
    WRITE (UNIT=n, FMT='(A,I5,A,10F15.10)') 'ZDEF ', kMax, ' LEVELS ',(a_h(k)+100000.0_r8*b_h(k),k=kMax,1,-1)
    it=1
    WRITE (UNIT=n, FMT='(A,I5,A,I2.2,A,I2.2,A,I4,A)') &
            'TDEF ', it, ' LINEAR ', ihr, 'Z', iday, cmth(mon), iyr, ' 6HR'
    WRITE (UNIT=n, FMT='(A)') '*'
    nof=10
    WRITE (UNIT=n, FMT='(A,I5)') 'VARS ', nof

    WRITE (UNIT=n, FMT='(A)') 'UVE10M    0 99 '//TRIM(uve10m)//' (M/S             )'
    WRITE (UNIT=n, FMT='(A)') 'VVE10M    0 99 '//TRIM(uve10m)//' (M/S             )'
    WRITE (UNIT=n, FMT='(A)') 'geop      0 99 '//TRIM(topo)//' (M^2/s^2         )'
    WRITE (UNIT=n, FMT='(A)') 'lnsp      0 99 '//TRIM(Psnm)//' (PA              )'
    WRITE (UNIT=n, FMT='(A,I5,A)') 'uvel    ', kMax,' 99 '//TRIM(uvel)//' (M/S             )'
    WRITE (UNIT=n, FMT='(A,I5,A)') 'vvel    ', kMax,' 99 '//TRIM(vvel)//' (M/S             )'
    WRITE (UNIT=n, FMT='(A,I5,A)') 'omeg    ', kMax,' 99 '//TRIM(omeg)//' (Pa/s            )'
    WRITE (UNIT=n, FMT='(A,I5,A)') 'temp    ', kMax,' 99 '//TRIM(temp)//' (K               )'
    WRITE (UNIT=n, FMT='(A,I5,A)') 'umes    ', kMax,' 99 '//TRIM(umes)//' (KG/KG           )'
    WRITE (UNIT=n, FMT='(A,I5,A)') 'dump    ', kMax,' 99 '//TRIM(dump)//' (%               )'
    WRITE (UNIT=n, FMT='(A)') 'ENDVARS'

!160 FORMAT(A40,2X,A4,2X,I8,3X,I4,4X,I3)
!161 FORMAT(' ',A40,2X,A4,2X,I8,3X,I4,4X,I3)    
  END SUBROUTINE LogWriteDire4

  SUBROUTINE LogWriteDire8(n, idate, ihr, iday, mon, iyr,a_h,b_h,tod,dcol,iMaxNew,jMaxNew,filename,FIELDGAUS)
    INTEGER, INTENT(IN ) :: n
    INTEGER, INTENT(IN ) :: iMaxNew
    INTEGER, INTENT(IN ) :: jMaxNew
    INTEGER, INTENT(IN ) :: idate(4)
    INTEGER, INTENT(IN ) :: ihr
    INTEGER, INTENT(IN ) :: iday
    INTEGER, INTENT(IN ) :: mon
    INTEGER, INTENT(IN ) :: iyr
    REAL   (KIND=r8), INTENT(IN ) :: a_h (kMax+1)
    REAL   (KIND=r8), INTENT(IN ) :: b_h (kMax+1)
    REAL   (KIND=r8), INTENT(IN ) :: tod

    REAL   (KIND=r8), INTENT(IN ) :: dcol(jMaxNew)
    CHARACTER(LEN=*), INTENT(IN ) :: filename
    LOGICAL         , INTENT(IN ) :: FIELDGAUS
    REAL      (KIND=r4) :: r
    CHARACTER (LEN= 4) :: imdl
    CHARACTER (LEN=150) :: jttl
    CHARACTER (LEN=21), PARAMETER :: ittl='CPTEC HIBRID VERS 2.0'
    CHARACTER (LEN= 4), PARAMETER :: nexp='0001'
    CHARACTER (LEN=40), PARAMETER :: uve10m = 'ZONAL WIND AT 10-M FROM SURFACE'
    CHARACTER (LEN=40), PARAMETER :: vve10m = 'MERID WIND AT 10-M FROM SURFACE'
    CHARACTER (LEN=40), PARAMETER :: topo   = 'GEOPOTENTIAL                   '
    CHARACTER (LEN=40), PARAMETER :: Psnm   = 'LN SURFACE PRESSURE            '
    CHARACTER (LEN=40), PARAMETER :: uvel   = 'ZONAL WIND (U)                 '
    CHARACTER (LEN=40), PARAMETER :: vvel   = 'MERIDIONAL WIND (V)            '
    CHARACTER (LEN=40), PARAMETER :: omeg   = 'VERTICAL WIND (W)              '
    CHARACTER (LEN=40), PARAMETER :: temp   = 'ABSOLUTE TEMPERATURE           '
    CHARACTER (LEN=40), PARAMETER :: umes   = 'SPECIFIC HUMIDITY              '
    CHARACTER (LEN=40), PARAMETER :: dump   = 'DUMP FIELD                     '

    CHARACTER (LEN=29), PARAMETER :: fmt1='(A40,2X,A4,2X,I8,3X,I4,4X,I3)'
    CHARACTER (LEN=3), PARAMETER :: cmth(12)=(/&
         'JAN','FEB','MAR','APR','MAY','JUN', &
         'JUL','AUG','SEP','OCT','NOV','DEC'/)
    REAL(KIND=r8)     , PARAMETER :: Undef=-2.56E+33_r8
    INTEGER :: i,j,it,nof,k



    WRITE (UNIT=n, FMT='(A)') 'DSET ^'//TRIM(filename)
    WRITE (UNIT=n, FMT='(A)') '*'
    WRITE (UNIT=n, FMT='(A)') 'OPTIONS SEQUENTIAL YREV BIG_ENDIAN'
    WRITE (UNIT=n, FMT='(A)') '*'
    WRITE (UNIT=n, FMT='(A,1PE9.2)') 'UNDEF ', Undef
    WRITE (UNIT=n, FMT='(A)') '*'
    WRITE (UNIT=n, FMT='(3A)') 'TITLE ',ittl,nexp
    WRITE (UNIT=n, FMT='(A)') '*'
    WRITE (UNIT=n, FMT='(A,I5,A,F8.3,F15.10)')'XDEF ', iMaxNew, ' LINEAR ', 0.0_r8, 360.0_r8/REAL(iMaxNew,r8)
    WRITE (UNIT=n, FMT='(A,I5,A)') 'YDEF ', jMaxNew, ' LEVELS '
    WRITE (UNIT=n, FMT='(8F15.10)') (dcol(j),j=jMaxNew,1,-1)
    WRITE (UNIT=n, FMT='(A,I5,A)') 'ZDEF ', kMax, ' LEVELS '
    WRITE (UNIT=n, FMT='(8F25.5)') (a_h(k)+100000.0_r8*b_h(k),k=kMax,1,-1)
    it=1
    WRITE (UNIT=n, FMT='(A,I5,A,I2.2,A,I2.2,A,I4,A)') &
            'TDEF ', it, ' LINEAR ', ihr, 'Z', iday, cmth(mon), iyr, ' 6HR'
    WRITE (UNIT=n, FMT='(A)') '*'
    nof=10
    WRITE (UNIT=n, FMT='(A,I5)') 'VARS ', nof

    WRITE (UNIT=n, FMT='(A)') 'UVE10M    0 99 '//TRIM(uve10m)//' (M/S             )'
    WRITE (UNIT=n, FMT='(A)') 'VVE10M    0 99 '//TRIM(uve10m)//' (M/S             )'
    WRITE (UNIT=n, FMT='(A)') 'geop      0 99 '//TRIM(topo)//' (M^2/s^2         )'
    WRITE (UNIT=n, FMT='(A)') 'lnsp      0 99 '//TRIM(Psnm)//' (PA              )'
    WRITE (UNIT=n, FMT='(A,I5,A)') 'uvel    ', kMax,' 99 '//TRIM(uvel)//' (M/S             )'
    WRITE (UNIT=n, FMT='(A,I5,A)') 'vvel    ', kMax,' 99 '//TRIM(vvel)//' (M/S             )'
    WRITE (UNIT=n, FMT='(A,I5,A)') 'omeg    ', kMax,' 99 '//TRIM(omeg)//' (Pa/s            )'
    WRITE (UNIT=n, FMT='(A,I5,A)') 'temp    ', kMax,' 99 '//TRIM(temp)//' (K               )'
    WRITE (UNIT=n, FMT='(A,I5,A)') 'umes    ', kMax,' 99 '//TRIM(umes)//' (KG/KG           )'
    WRITE (UNIT=n, FMT='(A,I5,A)') 'dump    ', kMax,' 99 '//TRIM(dump)//' (%               )'

    WRITE (UNIT=n, FMT='(A)') 'ENDVARS'

!160 FORMAT(A40,2X,A4,2X,I8,3X,I4,4X,I3)
!161 FORMAT(' ',A40,2X,A4,2X,I8,3X,I4,4X,I3)    
  END SUBROUTINE LogWriteDire8


  SUBROUTINE WriteDire8(n, idate, ihr, iday, mon, iyr, a_hybr, b_hybr,tod)
    INTEGER , INTENT(IN ) :: n
    INTEGER , INTENT(IN ) :: idate(4)
    INTEGER , INTENT(IN ) :: ihr
    INTEGER , INTENT(IN ) :: iday
    INTEGER , INTENT(IN ) :: mon
    INTEGER , INTENT(IN ) :: iyr
    REAL    (KIND=r8), INTENT(IN ) :: a_hybr(kMax+1) 
    REAL    (KIND=r8), INTENT(IN ) :: b_hybr(kMax+1) 
    REAL    (KIND=r8), INTENT(IN ) :: tod
    REAL    (KIND=r8) :: r

    CHARACTER (LEN= 5) :: imdl
    CHARACTER (LEN=150) :: jttl
    CHARACTER (LEN=21), PARAMETER :: ittl='CPTEC HIBRID VERS 2.0'
    CHARACTER (LEN= 4), PARAMETER :: nexp='0001'
    CHARACTER (LEN=40), PARAMETER :: orogra = 'TOPOGRAPHY' 
    CHARACTER (LEN=40), PARAMETER :: lsimsk = 'LAND SEA ICE MASK'
    CHARACTER (LEN=40), PARAMETER :: lnsurf = 'LN SURFACE PRESSURE'
    CHARACTER (LEN=40), PARAMETER :: divrgn = 'DIVERGENCE'
    CHARACTER (LEN=40), PARAMETER :: vortic = 'VORTICITY'
    CHARACTER (LEN=40), PARAMETER :: spechu = 'SPECIFIC HUMIDITY'
    CHARACTER (LEN=40), PARAMETER :: tempvi = 'VIRTUAL TEMPERATURE'
    CHARACTER (LEN=40), PARAMETER :: surfte = 'SURFACE TEMPERATURE'
    CHARACTER (LEN=40), PARAMETER :: srfrou = 'ROUGHNESS LENGTH'
    CHARACTER (LEN=40), PARAMETER :: deepte = 'DEEP SOIL TEMPERATURE'
    CHARACTER (LEN=40), PARAMETER :: stcnpy = 'STORAGE ON CANOPY'
    CHARACTER (LEN=40), PARAMETER :: stgrnd = 'STORAGE ON GROUND'
    CHARACTER (LEN=40), PARAMETER :: wt1soi = 'SOIL WETNESS OF SURFACE'
    CHARACTER (LEN=40), PARAMETER :: wt2soi = 'SOIL WETNESS OF ROOT ZONE'
    CHARACTER (LEN=40), PARAMETER :: wt3soi = 'SOIL WETNESS OF DRAINAGE ZONE'
    CHARACTER (LEN=40), PARAMETER :: temp2m = 'TEMPERATURE AT 2-M FROM SURFACE' 
    CHARACTER (LEN=40), PARAMETER :: umes2m = 'SPECIFIC HUMIDITY AT 2-M FROM SURFACE'     
    CHARACTER (LEN=40), PARAMETER :: uve10m = 'ZONAL WIND AT 10-M FROM SURFACE'
    CHARACTER (LEN=40), PARAMETER :: vve10m = 'MERID WIND AT 10-M FROM SURFACE'
    CHARACTER (LEN=40), PARAMETER :: veget0 = 'MASK VEGETATION'    
    CHARACTER (LEN=40), PARAMETER :: soiltx = 'MASK SOIL TEXTURE CLASSES'
    CHARACTER (LEN=40), PARAMETER :: doair0 = 'PARTIAL OXYGEN DENSITY '
    !JGERD -> fields necessary to data assimilation
    CHARACTER (LEN=40), PARAMETER :: sfcstg = 'SURFACE SOIL TEMPERATURE'
    CHARACTER (LEN=40), PARAMETER :: vegcov = 'VEGETATION COVER'
    CHARACTER (LEN=40), PARAMETER :: snodpt = 'SNOW DEPTH'
    CHARACTER (LEN=40), PARAMETER :: icwcon = 'ICE MIXING RATIO PROGNOSTIC'
    CHARACTER (LEN=40), PARAMETER :: lqwcon = 'LIQ MIXING RATIO PROGNOSTIC'
    CHARACTER (LEN=40), PARAMETER :: concld = 'CLOUD TOTAL PROGNOSTIC'

    CHARACTER (LEN=29), PARAMETER :: fmt1='(A40,2X,A4,2X,I8,3X,I4,4X,I3)'



    r = tod
    jttl='CPTEC AGCM REVIS 1.0 2000  T    L    COLD [OPTIONS |surface=01||pbl=04|]'
    WRITE (jttl(29:32), '(i4.4)') mMax-1
    WRITE (jttl(34:36), '(i2.2)') kMax
    WRITE (jttl(61:62), '(i2.2)') schemes
    WRITE (jttl(69:70), '(i2.2)') atmpbl
    WRITE (imdl, '(A1,I4.4)') 'T', mMax-1

    WRITE (UNIT=n, FMT='(A21)')   ittl
    WRITE (UNIT=n, FMT='(A4,1X,A4,1X,A5,1X,11I5,1X,A4)') &
         nexp, 'SEQU', imdl, mMax, kmax, kmax, &
         ihr, iday, mon, iyr, idate, 'TAPE'
    WRITE (UNIT=n, FMT='(A41)')   jttl
    WRITE (UNIT=n, FMT='(5E16.8)')   a_hybr
    WRITE (UNIT=n, FMT='(5E16.8)')   b_hybr
    WRITE (UNIT=n, FMT=fmt1) orogra, 'FIXD', 2*mnMax, 1, 10
    WRITE (UNIT=n, FMT=fmt1) lsimsk, 'FIXD', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) lnsurf, 'PROG', 2*mnMax, 1, 142
    WRITE (UNIT=n, FMT=fmt1) divrgn, 'PROG', 2*mnMax, kmax, 50
    WRITE (UNIT=n, FMT=fmt1) vortic, 'PROG', 2*mnMax, kmax, 50
    WRITE (UNIT=n, FMT=fmt1) spechu, 'PROG', 2*mnMax, kmax, 0
    WRITE (UNIT=n, FMT=fmt1) tempvi, 'PROG', 2*mnMax, kmax, 40
    WRITE (UNIT=n, FMT=fmt1) srfrou, 'PROG', ijmax, 1, 10
    WRITE (UNIT=n, FMT=fmt1) surfte, 'PROG', ijmax, 1, 40
    WRITE (UNIT=n, FMT=fmt1) deepte, 'PROG', ijmax, 1, 40
    WRITE (UNIT=n, FMT=fmt1) stcnpy, 'PROG', ijmax, 1, 10
    WRITE (UNIT=n, FMT=fmt1) stgrnd, 'PROG', ijmax, 1, 10
    WRITE (UNIT=n, FMT=fmt1) wt1soi, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) wt2soi, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) wt3soi, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) temp2m, 'PROG', ijmax, 1, 40
    WRITE (UNIT=n, FMT=fmt1) umes2m, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) uve10m, 'PROG', ijmax, 1, 60
    WRITE (UNIT=n, FMT=fmt1) vve10m, 'PROG', ijmax, 1, 60
    WRITE (UNIT=n, FMT=fmt1) veget0, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) soiltx, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) doair0, 'PROG', ijmax, 1, 201
   !JGERD -> fields necessary to data assimilation
    WRITE (UNIT=n, FMT=fmt1) sfcstg, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) vegcov, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) snodpt, 'PROG', ijmax, 1, 0
    WRITE (UNIT=n, FMT=fmt1) icwcon, 'PROG', ijmax, kmax, 0
    WRITE (UNIT=n, FMT=fmt1) lqwcon, 'PROG', ijmax, kmax, 0
    WRITE (UNIT=n, FMT=fmt1) concld, 'PROG', ijmax, kmax, 0

  END SUBROUTINE WriteDire8
  SUBROUTINE ReadProgHead4(n, ifday, tod, idate, idatec)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    INTEGER(KIND=i4), INTENT(OUT) :: ifday
    REAL   (KIND=r4), INTENT(OUT) :: tod
    INTEGER(KIND=i4), INTENT(OUT) :: idate(4)
    INTEGER(KIND=i4), INTENT(OUT) :: idatec(4)
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadProgHead4)**"
    READ(UNIT=n)ifday, tod, idate, idatec
  END SUBROUTINE ReadProgHead4
  SUBROUTINE ReadProgHead8(n, ifday, tod, idate, idatec)
    INTEGER, INTENT(IN)  :: n
    INTEGER, INTENT(OUT) :: ifday
    REAL   (KIND=r8), INTENT(OUT) :: tod
    INTEGER, INTENT(OUT) :: idate(4)
    INTEGER, INTENT(OUT) :: idatec(4)
    INTEGER(KIND=i4) :: iaux(10)
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadProgHead8)**"
    READ(UNIT=n)iaux
    ifday  = INT(iaux(  1 ))
    tod    = INT(iaux(  2 ))
    idate  = INT(iaux(3:6 ))
    idatec = INT(iaux(7:10))
  END SUBROUTINE ReadProgHead8
  SUBROUTINE GReadProgHead4(n, ifday, tod, idate, idatec)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    INTEGER(KIND=i4), INTENT(OUT) :: ifday
    REAL   (KIND=r4), INTENT(OUT) :: tod
    INTEGER(KIND=i4), INTENT(OUT) :: idate(4)
    INTEGER(KIND=i4), INTENT(OUT) :: idatec(4)
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadProgHead4)**"
    READ(UNIT=n)ifday, tod, idate, idatec
  END SUBROUTINE GReadProgHead4
  SUBROUTINE GReadProgHead8(n, ifday, tod, idate, idatec)
    INTEGER, INTENT(IN)  :: n
    INTEGER, INTENT(OUT) :: ifday
    REAL   (KIND=r8), INTENT(OUT) :: tod
    INTEGER, INTENT(OUT) :: idate(4)
    INTEGER, INTENT(OUT) :: idatec(4)
    INTEGER :: iaux(10)
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadProgHead8)**"
    READ(UNIT=n)iaux
    ifday  = INT(iaux(  1 ))
    tod    = INT(iaux(  2 ))
    idate  = INT(iaux(3:6 ))
    idatec = INT(iaux(7:10))
  END SUBROUTINE GReadProgHead8
  SUBROUTINE WriteProgHead4(n, ifday, tod, idate, idatec)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    INTEGER(KIND=i4), INTENT(IN)  :: ifday
    REAL   (KIND=r4), INTENT(IN)  :: tod
    INTEGER(KIND=i4), INTENT(IN)  :: idate(4)
    INTEGER(KIND=i4), INTENT(IN)  :: idatec(4)
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteProgHead4)**"
    WRITE(UNIT=n)ifday, tod, idate, idatec
  END SUBROUTINE WriteProgHead4
  SUBROUTINE WriteProgHead8(n, ifday, tod, idate, idatec)
    INTEGER, INTENT(IN)  :: n
    INTEGER, INTENT(IN)  :: ifday
    REAL   (KIND=r8), INTENT(IN)  :: tod
    INTEGER, INTENT(IN)  :: idate(4)
    INTEGER, INTENT(IN)  :: idatec(4)
    INTEGER(KIND=i4) :: iaux(10)
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteProgHead8)**"
    iaux(  1 ) = INT (ifday , i4)
    iaux(  2 ) = INT (tod   , i4)
    iaux(3:6 ) = INT (idate , i4)
    iaux(7:10) = INT (idatec, i4)
    WRITE(UNIT=n)iaux
  END SUBROUTINE WriteProgHead8


  SUBROUTINE ReadMs(mPerLat,fname,jMax)
    CHARACTER(LEN=*), INTENT(IN ) :: fname
    INTEGER,          INTENT(IN ) :: jMax
    INTEGER,          INTENT(OUT) :: mPerLat(jMax)
    INTEGER(KIND=i4) :: jsize
    INTEGER(KIND=i4) :: mPerLat_Aux(jMax)
    INTEGER          ::  ierr
    CHARACTER(LEN=8) :: c0, c1
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadMs)**"

    OPEN (UNIT=nfwaves, FILE=TRIM(fname), STATUS="old", ACTION="read",&
         ACCESS="sequential", FORM="formatted", IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(c0,"(i8)") ierr
       CALL FatalError(h//" open "//TRIM(fname)//" for read returned iostat="//&
            TRIM(ADJUSTL(c0)))
       STOP
    END IF

    READ (nfwaves,*,IOSTAT=ierr) jsize
    IF (ierr /= 0) THEN
       WRITE(c0,"(i8)") ierr
       CALL FatalError(h//" read first record of "//TRIM(fname)//" returned iostat="//&
            TRIM(ADJUSTL(c0)))
       STOP
    ELSE IF (jsize /= jMax) THEN
       WRITE(c0,"(i8)") jMax
       WRITE(c1,"(i8)") jSize
       CALL FatalError(h//" latitudes ("//TRIM(ADJUSTL(c1))//") on file"//&
            TRIM(fname)//" do not match required amount ("//&
            TRIM(ADJUSTL(c0))//")")
       STOP
    END IF

    READ (nfwaves,*,IOSTAT=ierr) mPerLat_Aux(1:jMax)
    mPerLat(1:jMax)=mPerLat_Aux(1:jMax)
    IF (ierr /= 0) THEN
       WRITE(c0,"(i8)") ierr
       CALL FatalError(h//" read second record of file "//&
            TRIM(fname)//" (containing mPerLat) returned iostat="//&
            TRIM(ADJUSTL(c0)))
       STOP
    END IF

    CLOSE (UNIT=nfwaves)
  END SUBROUTINE ReadMs




  SUBROUTINE WriteMs(mPerLat,fname,jMax)
    CHARACTER(LEN=*), INTENT(IN ) :: fname
    INTEGER,          INTENT(IN ) :: jMax
    INTEGER,          INTENT(IN ) :: mPerLat(jMax)
    INTEGER(KIND=i4) :: jsize
    INTEGER(KIND=i4) :: mPerLat_Aux(jMax)
    INTEGER          :: ierr
    CHARACTER(LEN=8) :: c0
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteMs)**"

    jsize=jMax
    OPEN (UNIT=nfwaves, FILE=TRIM(fname), STATUS="replace", ACTION="write",&
         ACCESS="sequential", FORM="formatted", IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(c0,"(i8)") ierr
       CALL FatalError(h//" open "//TRIM(fname)//" for write returned iostat="//&
            TRIM(ADJUSTL(c0)))
       STOP
    ELSE
       CALL MsgOut(h," writing file "//TRIM(fname))
       CALL MsgDump(h," writing file "//TRIM(fname))
    END IF

    WRITE (nfwaves,*,IOSTAT=ierr) jsize
    IF (ierr /= 0) THEN
       WRITE(c0,"(i8)") ierr
       CALL FatalError(h//" write first record of "//TRIM(fname)//" returned iostat="//&
            TRIM(ADJUSTL(c0)))
       STOP
    END IF
    mPerLat_Aux(1:jMax)=mPerLat(1:jMax)
    WRITE (nfwaves,*,IOSTAT=ierr)mPerLat_Aux(1:jMax)
    IF (ierr /= 0) THEN
       WRITE(c0,"(i8)") ierr
       CALL FatalError(h//" write second record of file "//&
            TRIM(fname)//" (containing mPerLat) returned iostat="//&
            TRIM(ADJUSTL(c0)))
       STOP
    END IF

    CLOSE (UNIT=nfwaves)
  END SUBROUTINE WriteMs




  SUBROUTINE ReadGauss(CosGaussColat, GaussWeights, fname, MaxDegree)
    CHARACTER(LEN=*), INTENT(IN ) :: fname
    INTEGER,          INTENT(IN ) :: MaxDegree
    REAL(KIND=r8),    INTENT(OUT) :: CosGaussColat(:)
    REAL(KIND=r8),    INTENT(OUT) :: GaussWeights(:)
    INTEGER(KIND=i4) :: jsize
    INTEGER          :: MaxDegree2, ierr
    CHARACTER(LEN=8) :: c0, c1
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadGauss)**"

    MaxDegree2=MaxDegree/2
    IF (SIZE(CosGaussColat) /= MaxDegree2) THEN
       CALL FatalError(h//" size of CosGaussColat is wrong")
       STOP
    ELSE IF (SIZE(GaussWeights) /= MaxDegree) THEN
       CALL FatalError(h//" size of GaussWeights is wrong")
       STOP
    END IF
    OPEN (UNIT=nfGauss, FILE=TRIM(fname), STATUS="old", ACTION="read",&
         ACCESS="sequential", FORM="unformatted", IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(c0,"(i8)") ierr
       CALL FatalError(h//" open "//TRIM(fname)//" for read returned iostat="//&
            TRIM(ADJUSTL(c0)))
       STOP
    END IF

    READ (nfGauss,IOSTAT=ierr) jsize
    IF (ierr /= 0) THEN
       WRITE(c0,"(i8)") ierr
       CALL FatalError(h//" read first record of "//TRIM(fname)//" returned iostat="//&
            TRIM(ADJUSTL(c0)))
       STOP
    ELSE IF (jsize /= MaxDegree2) THEN
       WRITE(c0,"(i8)") MaxDegree2
       WRITE(c1,"(i8)") jSize
       CALL FatalError(h//" MaxDegree/2 ("//TRIM(ADJUSTL(c1))//") on file"//&
            TRIM(fname)//" do not match required amount ("//&
            TRIM(ADJUSTL(c0))//")")
       STOP
    END IF

    READ (nfGauss, IOSTAT=ierr) CosGaussColat
    IF (ierr /= 0) THEN
       WRITE(c0,"(i8)") ierr
       CALL FatalError(h//" read second record of "//TRIM(fname)//&
            " (containing CosGaussColat) returned iostat="//&
            TRIM(ADJUSTL(c0)))
       STOP
    END IF

    READ (nfGauss, IOSTAT=ierr) GaussWeights
    IF (ierr /= 0) THEN
       WRITE(c0,"(i8)") ierr
       CALL FatalError(h//" read third record of "//TRIM(fname)//&
            " (containing GaussWeights) returned iostat="//&
            TRIM(ADJUSTL(c0)))
       STOP
    END IF
    CLOSE (UNIT=nfGauss)
  END SUBROUTINE ReadGauss




  SUBROUTINE WriteGauss(CosGaussColat, GaussWeights, fname, MaxDegree)
    CHARACTER(LEN=*), INTENT(IN ) :: fname
    INTEGER,          INTENT(IN ) :: MaxDegree
    REAL(KIND=r8),    INTENT(IN ) :: CosGaussColat(:)
    REAL(KIND=r8),    INTENT(IN ) :: GaussWeights(:)
    INTEGER(KIND=i4) :: MaxDegree2
    INTEGER          :: ierr
    CHARACTER(LEN=8) :: c0
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteGauss)**"

    MaxDegree2=MaxDegree/2
    IF (SIZE(CosGaussColat) /= MaxDegree2) THEN
       CALL FatalError(h//" size of CosGaussColat is wrong")
       STOP
    ELSE IF (SIZE(GaussWeights) /= MaxDegree) THEN
       CALL FatalError(h//" size of GaussWeights is wrong")
       STOP
    END IF

    OPEN (UNIT=nfGauss, FILE=TRIM(fname), STATUS="replace", ACTION="write",&
         ACCESS="sequential", FORM="unformatted", IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(c0,"(i8)") ierr
       CALL FatalError(h//" open "//TRIM(fname)//" for write returned iostat="//&
            TRIM(ADJUSTL(c0)))
       STOP
    ELSE
       CALL MsgOut(h," writing file "//TRIM(fname))
       CALL MsgDump(h," writing file "//TRIM(fname))
    END IF

    WRITE (nfGauss,IOSTAT=ierr) MaxDegree2
    IF (ierr /= 0) THEN
       WRITE(c0,"(i8)") ierr
       CALL FatalError(h//" write first record of "//TRIM(fname)//" returned iostat="//&
            TRIM(ADJUSTL(c0)))
       STOP
    END IF

    WRITE (nfGauss,IOSTAT=ierr) CosGaussColat
    IF (ierr /= 0) THEN
       WRITE(c0,"(i8)") ierr
       CALL FatalError(h//" write second record of "//TRIM(fname)//&
            " (containing CosGaussColat) returned iostat="//&
            TRIM(ADJUSTL(c0)))
       STOP
    END IF

    WRITE (nfGauss,IOSTAT=ierr) GaussWeights
    IF (ierr /= 0) THEN
       WRITE(c0,"(i8)") ierr
       CALL FatalError(h//" write third record of "//TRIM(fname)//&
            " (containing GaussWeights) returned iostat="//&
            TRIM(ADJUSTL(c0)))
       STOP
    END IF

    CLOSE (UNIT=nfGauss)
  END SUBROUTINE WriteGauss


  SUBROUTINE ReadLandSeaMask2_4(fname, lsmk)
    CHARACTER(LEN=*), INTENT(IN ) :: fname
    REAL(KIND=r4),    INTENT(OUT) :: lsmk(imax*jMax)
    INTEGER(KIND=i4) :: ierr
    INTEGER          :: LRecIn
    REAL(KIND=r4) ::  bfr(imax,jMax)
    INTEGER :: i,j,ij
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadLandSeaMask2_4)**" 
    CHARACTER(LEN=256) :: line
    bfr=0.0_r4
    INQUIRE (IOLENGTH=LRecIn) bfr
    OPEN (UNIT=nfsst, FILE=TRIM(fname),FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIn, &
         ACTION='READ', STATUS='OLD', IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(line,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fname), ierr
       CALL FatalError(h//trim(line))
    END IF
    READ (UNIT=nfsst, REC=1, IOSTAT=ierr) bfr
    IF (ierr /= 0) THEN
       WRITE(line,FMT="('**(ERROR)** Read file ',a,' returned iostat=',i4)") &
            TRIM(fname), ierr
       CALL FatalError(h//trim(line))
    END IF
    ij=0
    DO j=1,jMax
       DO i=1, iMax
           ij=ij+1
          lsmk(ij)=bfr(i,j)
       END DO
    END DO
    CLOSE (UNIT=nfsst)
  END SUBROUTINE ReadLandSeaMask2_4

  SUBROUTINE ReadLandSeaMask2_8(fname, lsmk)
    CHARACTER(LEN=*  ), INTENT(IN ) :: fname
    REAL     (KIND=r8),    INTENT(OUT) :: lsmk(imax*jMax)
    INTEGER(KIND=i4) :: ierr
    INTEGER          :: i,j,ij,LRecIn
    REAL(KIND=r4) ::  bfr(imax,jMax)
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadLandSeaMask2_8)**" 
    CHARACTER(LEN=256) :: line
    bfr=0.0_r4
    INQUIRE (IOLENGTH=LRecIn) bfr
    OPEN (UNIT=nfsst, FILE=TRIM(fname),FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIn, &
         ACTION='READ', STATUS='OLD', IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(line,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fname), ierr
       CALL FatalError(h//trim(line))
    END IF
    READ (UNIT=nfsst, REC=1, IOSTAT=ierr) bfr
    IF (ierr /= 0) THEN
       WRITE(line,FMT="('**(ERROR)** Read file ',a,' returned iostat=',i4)") &
            TRIM(fname), ierr
       CALL FatalError(h//trim(line))
    END IF
    ij=0
    DO j=1,jMax
       DO i=1, iMax
          ij=ij+1
          lsmk(ij)=REAL(bfr(i,j),KIND=r8)
       END DO
    END DO
    CLOSE(UNIT=nfsst)
  END SUBROUTINE ReadLandSeaMask2_8

  SUBROUTINE  LandSeaMask4(ifsst,labeli,intsst,sstlag,fNameSSTAOI,rlsm)

    INTEGER  , INTENT(IN   )  :: ifsst
    INTEGER  , INTENT(INOUT)  :: intsst
    REAL     (KIND=r8), INTENT(OUT  )  :: sstlag
    REAL     (KIND=r4), INTENT(INOUT)  :: rlsm(:,:)
    CHARACTER(LEN=10 ), INTENT(IN   )  :: labeli
    CHARACTER(LEN=*  ), INTENT(IN   )  :: fNameSSTAOI
    INTEGER                            :: lrecl
    INTEGER(KIND=i4)                   :: nsst
    REAL     (KIND=r8)                 :: dlag
    REAL     (KIND=r4),ALLOCATABLE     :: var4(:,:)  
    INTEGER                   :: j
    INTEGER                   :: jMax
    INTEGER                   :: i  
    INTEGER                   :: imax
    INTEGER   :: ierr
    CHARACTER(LEN=*), PARAMETER :: h="**(LandSeaMask4)**" 
    CHARACTER(LEN=256) :: line

    jMax = SIZE(rlsm,2)   
    imax = SIZE(rlsm,1)
    ALLOCATE(var4(imax,jMax))

    ! Use open statement to open direct access file when ifsst .ge. 4:

    INQUIRE (IOLENGTH=lrecl) rlsm
    !lrecl=lrecl/2

    OPEN (UNIT=nfsst,file=TRIM(fNameSSTAOI),ACCESS='direct',FORM='unformatted',&
         RECL=lrecl,STATUS='old',IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(line,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fNameSSTAOI), ierr
       CALL FatalError(h//trim(line))
    END IF

    READ (UNIT=nfsst,REC=1)nsst,labelsi,labelsj
    WRITE(line,FMT='(I6,1X,a,1X,a)')nsst,labelsi,labelsj
    CALL MsgOne(h,trim(line))
    CALL daylag (labelsi,labelsj,dlag,intsst)
    intsst=NINT(dlag)
    IF (intsst > 10) intsst=-intsst
    CALL daylag (labelsi,labeli,dlag,intsst)
    sstlag=dlag
    IF ( TRIM(NMSST) == 'sstdyd' )sstlag=ABS(sstlag)
    WRITE(line,FMT='(a)')' Direct Access SST File:'
    CALL MsgOne(h,trim(line))
    WRITE(line,FMT='(1x,a4,3(1x,a2))') &
         labelsi(1:4),labelsi(5:6),labelsi(7:8),labelsi(9:10)
    CALL MsgOne(h,trim(line))
    WRITE(line,FMT='(1x,a4,3(1x,a2))') &
         labelsj(1:4),labelsj(5:6),labelsj(7:8),labelsj(9:10)
    CALL MsgOne(h,trim(line))
    WRITE(line,FMT='(i7,a)')intsst,'  Days'
    CALL MsgOne(h,trim(line))
    WRITE(line,FMT='(1x,a4,3(1x,a2))') &
         labelsi(1:4),labelsi(5:6),labelsi(7:8),labelsi(9:10)
    CALL MsgOne(h,line)
    WRITE(line,FMT='(1x,a4,3(1x,a2))') &
         labeli(1:4),labeli(5:6),labeli(7:8),labeli(9:10)
    CALL MsgOne(h,trim(line))
    IF (intsst > 0) THEN
       WRITE(line,FMT='(f7.1,a)')sstlag,'  Days'
       CALL MsgOne(h,trim(line))
    ELSE
       WRITE(line,FMT='(f7.1,a)')sstlag,'  Months'
       CALL MsgOne(h,trim(line))
    END IF
    IF (sstlag < 0) THEN
       WRITE(line, FMT=336) ifsst, sstlag
       CALL FatalError(h//trim(line))
    END IF
    READ (UNIT=nfsst,REC=2) var4
    CLOSE (UNIT=nfsst)

    DO j=1,jMax
       DO i=1,imax
          rlsm(i,j)=var4(i,j)
       END DO
    END DO
    DEALLOCATE(var4)
336 FORMAT(' FOR IFSST=',I5,' SSTLAG MUST BE SET NONNEGATIVE.  NOT ',G12.5)  
  END SUBROUTINE  LandSeaMask4
  SUBROUTINE  LandSeaMask8(ifsst,labeli,intsst,sstlag,fNameSSTAOI,rlsm)

    INTEGER  , INTENT(IN   )  :: ifsst
    INTEGER  , INTENT(INOUT)  :: intsst
    REAL     (KIND=r8), INTENT(OUT  )  :: sstlag
    REAL     (KIND=r8), INTENT(INOUT)  :: rlsm(:,:)
    CHARACTER(LEN=10 ), INTENT(IN   )  :: labeli
    CHARACTER(LEN=*  ), INTENT(IN   )  :: fNameSSTAOI
    INTEGER                            :: lrecl
    INTEGER(KIND=i4)                   :: nsst
    REAL     (KIND=r8)                 :: dlag
    REAL     (KIND=r4),ALLOCATABLE     :: var4(:,:)  
    INTEGER                   :: j
    INTEGER                   :: jMax
    INTEGER                   :: i  
    INTEGER                   :: imax
    INTEGER   :: ierr
    CHARACTER(LEN=*), PARAMETER :: h="**(LandSeaMask8)**" 
    CHARACTER(LEN=256) :: line

    jMax = SIZE(rlsm,2)   
    imax = SIZE(rlsm,1)
    ALLOCATE(var4(imax,jMax))


    ! Use open statement to open direct access file when ifsst .ge. 4:

    INQUIRE (IOLENGTH=lrecl) rlsm
    !lrecl=lrecl/2

    OPEN (UNIT=nfsst,file=TRIM(fNameSSTAOI),ACCESS='direct',FORM='unformatted',&
         RECL=lrecl,STATUS='old',IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(line,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fNameSSTAOI), ierr
       CALL FatalError(h//trim(line))
    END IF

    READ (UNIT=nfsst,REC=1)nsst,labelsi,labelsj
    WRITE(line,FMT='(I6,1X,a,1X,a)')nsst,labelsi,labelsj
    CALL MsgOne(h,trim(line))    
    CALL daylag (labelsi,labelsj,dlag,intsst)
    intsst=NINT(dlag)
    IF (intsst > 10) intsst=-intsst
    CALL daylag (labelsi,labeli,dlag,intsst)
    sstlag=dlag
    IF ( TRIM(NMSST) == 'sstdyd' )sstlag=ABS(sstlag)
    WRITE(line,FMT='(a)')' Direct Access SST File:'
    CALL MsgOne(h,trim(line))
    WRITE(line,FMT='(1x,a4,3(1x,a2))') &
         labelsi(1:4),labelsi(5:6),labelsi(7:8),labelsi(9:10)
    CALL MsgOne(h,trim(line))
    WRITE(line,FMT='(1x,a4,3(1x,a2))') &
         labelsj(1:4),labelsj(5:6),labelsj(7:8),labelsj(9:10)
    CALL MsgOne(h,trim(line))
    WRITE(line,FMT='(i7,a)')intsst,'  Days'
    CALL MsgOne(h,trim(line))
    WRITE(line,FMT='(1x,a4,3(1x,a2))') &
         labelsi(1:4),labelsi(5:6),labelsi(7:8),labelsi(9:10)
    CALL MsgOne(h,trim(line))
    WRITE(line,FMT='(1x,a4,3(1x,a2))') &
         labeli(1:4),labeli(5:6),labeli(7:8),labeli(9:10)
    CALL MsgOne(h,trim(line))
    IF (intsst > 0) THEN
       WRITE(line,FMT='(f7.1,a)')sstlag,'  Days'
       CALL MsgOne(h,trim(line))
    ELSE
       WRITE(line,FMT='(f7.1,a)')sstlag,'  Months'
       CALL MsgOne(h,trim(line))
    END IF
    IF (sstlag < 0) THEN
       WRITE(line,FMT=336) ifsst, sstlag
       CALL FatalError(h//trim(line))
    END IF
    READ (UNIT=nfsst,REC=2) var4
    CLOSE (UNIT=nfsst)

    DO j=1,jMax
       DO i=1,imax
          rlsm(i,j)=var4(i,j)
       END DO
    END DO
    DEALLOCATE(var4)
336 FORMAT(' FOR IFSST=',I5,' SSTLAG MUST BE SET NONNEGATIVE.  NOT ',G12.5)  
  END SUBROUTINE  LandSeaMask8

  SUBROUTINE  LandFlxCO2Mask4(ifco2flx,labeli,intflxco2,Flxco2lag,fNameCO2FLX,rlsm)

    INTEGER           , INTENT(IN   )  :: ifco2flx
    INTEGER           , INTENT(INOUT)  :: intflxco2
    REAL     (KIND=r8), INTENT(OUT  )  :: Flxco2lag
    REAL     (KIND=r4), INTENT(INOUT)  :: rlsm(:,:)
    CHARACTER(LEN=10 ), INTENT(IN   )  :: labeli
    CHARACTER(LEN=*  ), INTENT(IN   )  :: fNameCO2FLX
    INTEGER                            :: lrecl
    INTEGER(KIND=i4)                   :: nflxco2
    REAL     (KIND=r8)                 :: dlag
    REAL     (KIND=r4),ALLOCATABLE     :: var4(:,:)  
    INTEGER                   :: j
    INTEGER                   :: jMax
    INTEGER                   :: i  
    INTEGER                   :: imax
    INTEGER   :: ierr
    CHARACTER(LEN=*), PARAMETER :: h="**(LandFlxCO2Mask4)**" 
    CHARACTER(LEN=256) :: line

    jMax = SIZE(rlsm,2)   
    imax = SIZE(rlsm,1)
    ALLOCATE(var4(imax,jMax))

    ! Use open statement to open direct access file when ifco2flx .ge. 4:

    INQUIRE (IOLENGTH=lrecl) rlsm
    !lrecl=lrecl/2

    OPEN (UNIT=nfco2fx,file=TRIM(fNameCO2FLX),ACCESS='direct',FORM='unformatted',&
         RECL=lrecl,STATUS='old',IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(line,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fNameCO2FLX), ierr
       CALL FatalError(h//trim(line))
    END IF

    READ (UNIT=nfco2fx,REC=1)nflxco2,labelsi_flxco2,labelsj_flxco2
    WRITE(line,FMT='(I6,1X,a,1X,a)')nflxco2,labelsi_flxco2,labelsj_flxco2
    CALL MsgOne(h,trim(line))
    CALL daylag (labelsi_flxco2,labelsj_flxco2,dlag,intflxco2)
    intflxco2=NINT(dlag)
    IF (intflxco2 > 10) intflxco2=-intflxco2
    CALL daylag (labelsi_flxco2,labeli,dlag,intflxco2)
    Flxco2lag=dlag
    IF ( TRIM(NMCO2FX) == 'co2flxdyd' )Flxco2lag=ABS(Flxco2lag)
    WRITE(line,FMT='(a)')' Direct Access co2flx File:'
    CALL MsgOne(h,trim(line))
    WRITE(line,FMT='(1x,a4,3(1x,a2))') &
         labelsi_flxco2(1:4),labelsi_flxco2(5:6),labelsi_flxco2(7:8),labelsi_flxco2(9:10)
    CALL MsgOne(h,trim(line))
    WRITE(line,FMT='(1x,a4,3(1x,a2))') &
         labelsj_flxco2(1:4),labelsj_flxco2(5:6),labelsj_flxco2(7:8),labelsj_flxco2(9:10)
    CALL MsgOne(h,trim(line))
    WRITE(line,FMT='(i7,a)')intflxco2,'  Days'
    CALL MsgOne(h,trim(line))
    WRITE(line,FMT='(1x,a4,3(1x,a2))') &
         labelsi_flxco2(1:4),labelsi_flxco2(5:6),labelsi_flxco2(7:8),labelsi_flxco2(9:10)
    CALL MsgOne(h,line)
    WRITE(line,FMT='(1x,a4,3(1x,a2))') &
         labeli(1:4),labeli(5:6),labeli(7:8),labeli(9:10)
    CALL MsgOne(h,trim(line))
    IF (intflxco2 > 0) THEN
       WRITE(line,FMT='(f7.1,a)')Flxco2lag,'  Days'
       CALL MsgOne(h,trim(line))
    ELSE
       WRITE(line,FMT='(f7.1,a)')Flxco2lag,'  Months'
       CALL MsgOne(h,trim(line))
    END IF
    IF (Flxco2lag < 0) THEN
       WRITE(line, FMT=336) ifco2flx, Flxco2lag
       CALL FatalError(h//trim(line))
    END IF
    !READ (UNIT=nfco2fx,REC=2) var4
    CLOSE (UNIT=nfco2fx)

    !DO j=1,jMax
    !   DO i=1,imax
    !      rlsm(i,j)=var4(i,j)
    !   END DO
    !END DO
    DEALLOCATE(var4)
336 FORMAT(' FOR ifco2flx=',I5,' Flxco2lag MUST BE SET NONNEGATIVE.  NOT ',G12.5)  
  END SUBROUTINE  LandFlxCO2Mask4

  SUBROUTINE  LandFlxCO2Mask8(ifco2flx,labeli,intflxco2,Flxco2lag,fNameCO2FLX,rlsm)

    INTEGER           , INTENT(IN   )  :: ifco2flx
    INTEGER           , INTENT(INOUT)  :: intflxco2
    REAL     (KIND=r8), INTENT(OUT  )  :: Flxco2lag
    REAL     (KIND=r8), INTENT(INOUT)  :: rlsm(:,:)
    CHARACTER(LEN=10 ), INTENT(IN   )  :: labeli
    CHARACTER(LEN=*  ), INTENT(IN   )  :: fNameCO2FLX
    INTEGER                            :: lrecl
    INTEGER(KIND=i4)                   :: nflxco2
    REAL     (KIND=r8)                 :: dlag
    REAL     (KIND=r4),ALLOCATABLE     :: var4(:,:)  
    INTEGER                   :: j
    INTEGER                   :: jMax
    INTEGER                   :: i  
    INTEGER                   :: imax
    INTEGER   :: ierr
    CHARACTER(LEN=*), PARAMETER :: h="**(LandFlxCO2Mask4)**" 
    CHARACTER(LEN=256) :: line

    jMax = SIZE(rlsm,2)   
    imax = SIZE(rlsm,1)
    ALLOCATE(var4(imax,jMax))

    ! Use open statement to open direct access file when ifco2flx .ge. 4:

    INQUIRE (IOLENGTH=lrecl) rlsm
    !lrecl=lrecl/2

    OPEN (UNIT=nfco2fx,file=TRIM(fNameCO2FLX),ACCESS='direct',FORM='unformatted',&
         RECL=lrecl,STATUS='old',IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(line,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fNameCO2FLX), ierr
       CALL FatalError(h//trim(line))
    END IF

    READ (UNIT=nfco2fx,REC=1)nflxco2,labelsi_flxco2,labelsj_flxco2
    WRITE(line,FMT='(I6,1X,a,1X,a)')nflxco2,labelsi_flxco2,labelsj_flxco2
    CALL MsgOne(h,trim(line))
    CALL daylag (labelsi_flxco2,labelsj_flxco2,dlag,intflxco2)
    intflxco2=NINT(dlag)
    IF (intflxco2 > 10) intflxco2=-intflxco2
    CALL daylag (labelsi_flxco2,labeli,dlag,intflxco2)
    Flxco2lag=dlag
    IF ( TRIM(NMCO2FX) == 'co2flxdyd' )Flxco2lag=ABS(Flxco2lag)
    WRITE(line,FMT='(a)')' Direct Access co2flx File:'
    CALL MsgOne(h,trim(line))
    WRITE(line,FMT='(1x,a4,3(1x,a2))') &
         labelsi_flxco2(1:4),labelsi_flxco2(5:6),labelsi_flxco2(7:8),labelsi_flxco2(9:10)
    CALL MsgOne(h,trim(line))
    WRITE(line,FMT='(1x,a4,3(1x,a2))') &
         labelsj_flxco2(1:4),labelsj_flxco2(5:6),labelsj_flxco2(7:8),labelsj_flxco2(9:10)
    CALL MsgOne(h,trim(line))
    WRITE(line,FMT='(i7,a)')intflxco2,'  Days'
    CALL MsgOne(h,trim(line))
    WRITE(line,FMT='(1x,a4,3(1x,a2))') &
         labelsi_flxco2(1:4),labelsi_flxco2(5:6),labelsi_flxco2(7:8),labelsi_flxco2(9:10)
    CALL MsgOne(h,line)
    WRITE(line,FMT='(1x,a4,3(1x,a2))') &
         labeli(1:4),labeli(5:6),labeli(7:8),labeli(9:10)
    CALL MsgOne(h,trim(line))
    IF (intflxco2 > 0) THEN
       WRITE(line,FMT='(f7.1,a)')Flxco2lag,'  Days'
       CALL MsgOne(h,trim(line))
    ELSE
       WRITE(line,FMT='(f7.1,a)')Flxco2lag,'  Months'
       CALL MsgOne(h,trim(line))
    END IF
    IF (Flxco2lag < 0) THEN
       WRITE(line, FMT=336) ifco2flx, Flxco2lag
       CALL FatalError(h//trim(line))
    END IF
    !READ (UNIT=nfco2fx,REC=2) var4
    CLOSE (UNIT=nfco2fx)

    !DO j=1,jMax
    !   DO i=1,imax
    !      rlsm(i,j)=var4(i,j)
    !   END DO
    !END DO
    DEALLOCATE(var4)
336 FORMAT(' FOR ifco2flx=',I5,' Flxco2lag MUST BE SET NONNEGATIVE.  NOT ',G12.5)  
  END SUBROUTINE  LandFlxCO2Mask8


  SUBROUTINE  LandSoilmMask4(ifsoilm,labeli,intsoilm,Soilmlag,fNameSoilMoist,rlsm)

    INTEGER           , INTENT(IN   )  :: ifsoilm
    INTEGER           , INTENT(INOUT)  :: intsoilm
    REAL     (KIND=r8), INTENT(OUT  )  :: soilmlag
    REAL     (KIND=r4), INTENT(INOUT)  :: rlsm(:,:)
    CHARACTER(LEN=10 ), INTENT(IN   )  :: labeli
    CHARACTER(LEN=*  ), INTENT(IN   )  :: fNameSoilMoist
    INTEGER                            :: lrecl
    INTEGER(KIND=i4)                   :: nsoilm
    REAL     (KIND=r8)                 :: dlag
    REAL     (KIND=r4),ALLOCATABLE     :: var4(:,:)  
    INTEGER                   :: j
    INTEGER                   :: jMax
    INTEGER                   :: i  
    INTEGER                   :: imax
    INTEGER   :: ierr
    CHARACTER(LEN=*), PARAMETER :: h="**(LandSoilmMask4)**" 
    CHARACTER(LEN=256) :: line

    jMax = SIZE(rlsm,2)   
    imax = SIZE(rlsm,1)
    ALLOCATE(var4(imax,jMax))

    ! Use open statement to open direct access file when ifsoilm .ge. 4:

    INQUIRE (IOLENGTH=lrecl) rlsm
    !lrecl=lrecl/2

    OPEN (UNIT=nfslm,file=TRIM(fNameSoilMoist),ACCESS='direct',FORM='unformatted',&
         RECL=lrecl,STATUS='old',IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(line,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fNameSoilMoist), ierr
       CALL FatalError(h//trim(line))
    END IF

    READ (UNIT=nfslm,REC=1)nsoilm,labelsi_soilm,labelsj_soilm
    WRITE(line,FMT='(I6,1X,a,1X,a)')nsoilm,labelsi_soilm,labelsj_soilm
    CALL MsgOne(h,trim(line))
    CALL daylag (labelsi_soilm,labelsj_soilm,dlag,intsoilm)
    intsoilm=NINT(dlag)
    IF (intsoilm > 10) intsoilm=-intsoilm
    CALL daylag (labelsi_soilm,labeli,dlag,intsoilm)
    soilmlag=dlag
    IF ( TRIM(NMSOILM) == 'soilmdyd' )soilmlag=ABS(soilmlag)
    WRITE(line,FMT='(a)')' Direct Access SOILM File:'
    CALL MsgOne(h,trim(line))
    WRITE(line,FMT='(1x,a4,3(1x,a2))') &
         labelsi_soilm(1:4),labelsi_soilm(5:6),labelsi_soilm(7:8),labelsi_soilm(9:10)
    CALL MsgOne(h,trim(line))
    WRITE(line,FMT='(1x,a4,3(1x,a2))') &
         labelsj_soilm(1:4),labelsj_soilm(5:6),labelsj_soilm(7:8),labelsj_soilm(9:10)
    CALL MsgOne(h,trim(line))
    WRITE(line,FMT='(i7,a)')intsoilm,'  Days'
    CALL MsgOne(h,trim(line))
    WRITE(line,FMT='(1x,a4,3(1x,a2))') &
         labelsi_soilm(1:4),labelsi_soilm(5:6),labelsi_soilm(7:8),labelsi_soilm(9:10)
    CALL MsgOne(h,line)
    WRITE(line,FMT='(1x,a4,3(1x,a2))') &
         labeli(1:4),labeli(5:6),labeli(7:8),labeli(9:10)
    CALL MsgOne(h,trim(line))
    IF (intsoilm > 0) THEN
       WRITE(line,FMT='(f7.1,a)')soilmlag,'  Days'
       CALL MsgOne(h,trim(line))
    ELSE
       WRITE(line,FMT='(f7.1,a)')soilmlag,'  Months'
       CALL MsgOne(h,trim(line))
    END IF
    IF (soilmlag < 0) THEN
       WRITE(line, FMT=336) ifsoilm, soilmlag
       CALL FatalError(h//trim(line))
    END IF
    !READ (UNIT=nfslm,REC=2) var4
    CLOSE (UNIT=nfslm)

    !DO j=1,jMax
    !   DO i=1,imax
    !      rlsm(i,j)=var4(i,j)
    !   END DO
    !END DO
    DEALLOCATE(var4)
336 FORMAT(' FOR IFSOILM=',I5,' SOILMLAG MUST BE SET NONNEGATIVE.  NOT ',G12.5)  
  END SUBROUTINE  LandSoilmMask4

  SUBROUTINE  LandSoilmMask8(ifsoilm,labeli,intsoilm,soilmlag,fNameSoilMoist,rlsm)

    INTEGER           , INTENT(IN   )  :: ifsoilm
    INTEGER           , INTENT(INOUT)  :: intsoilm
    REAL     (KIND=r8), INTENT(OUT  )  :: soilmlag
    REAL     (KIND=r8), INTENT(INOUT)  :: rlsm(:,:)
    CHARACTER(LEN=10 ), INTENT(IN   )  :: labeli
    CHARACTER(LEN=*  ), INTENT(IN   )  :: fNameSoilMoist
    INTEGER                            :: lrecl
    INTEGER(KIND=i4)                   :: nsoilm
    REAL     (KIND=r8)                 :: dlag
    REAL     (KIND=r4),ALLOCATABLE     :: var4(:,:)  
    INTEGER                   :: j
    INTEGER                   :: jMax
    INTEGER                   :: i  
    INTEGER                   :: imax
    INTEGER   :: ierr
    CHARACTER(LEN=*), PARAMETER :: h="**(LandSeaMask8)**" 
    CHARACTER(LEN=256) :: line

    jMax = SIZE(rlsm,2)   
    imax = SIZE(rlsm,1)
    ALLOCATE(var4(imax,jMax))


    ! Use open statement to open direct access file when ifsst .ge. 4:

    INQUIRE (IOLENGTH=lrecl) rlsm
    !lrecl=lrecl/2

    OPEN (UNIT=nfslm,file=TRIM(fNameSoilMoist),ACCESS='direct',FORM='unformatted',&
         RECL=lrecl,STATUS='old',IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(line,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fNameSoilMoist), ierr
       CALL FatalError(h//trim(line))
    END IF

    READ (UNIT=nfslm,REC=1)nsoilm,labelsi_soilm,labelsj_soilm
    WRITE(line,FMT='(I6,1X,a,1X,a)')nsoilm,labelsi_soilm,labelsj_soilm
    CALL MsgOne(h,trim(line))
    CALL daylag (labelsi_soilm,labelsj_soilm,dlag,intsoilm)
    intsoilm=NINT(dlag)
    IF (intsoilm > 10) intsoilm=-intsoilm
    CALL daylag (labelsi_soilm,labeli,dlag,intsoilm)
    soilmlag=dlag
    IF ( TRIM(NMSOILM) == 'soilmdyd' )soilmlag=ABS(soilmlag)
    WRITE(line,FMT='(a)')' Direct Access SOILM File:'
    CALL MsgOne(h,trim(line))
    WRITE(line,FMT='(1x,a4,3(1x,a2))') &
         labelsi_soilm(1:4),labelsi_soilm(5:6),labelsi_soilm(7:8),labelsi_soilm(9:10)
    CALL MsgOne(h,trim(line))
    WRITE(line,FMT='(1x,a4,3(1x,a2))') &
         labelsj_soilm(1:4),labelsj_soilm(5:6),labelsj_soilm(7:8),labelsj_soilm(9:10)
    CALL MsgOne(h,trim(line))
    WRITE(line,FMT='(i7,a)')intsoilm,'  Days'
    CALL MsgOne(h,trim(line))
    WRITE(line,FMT='(1x,a4,3(1x,a2))') &
         labelsi_soilm(1:4),labelsi_soilm(5:6),labelsi_soilm(7:8),labelsi_soilm(9:10)
    CALL MsgOne(h,line)
    WRITE(line,FMT='(1x,a4,3(1x,a2))') &
         labeli(1:4),labeli(5:6),labeli(7:8),labeli(9:10)
    CALL MsgOne(h,trim(line))
    IF (intsoilm > 0) THEN
       WRITE(line,FMT='(f7.1,a)')soilmlag,'  Days'
       CALL MsgOne(h,trim(line))
    ELSE
       WRITE(line,FMT='(f7.1,a)')soilmlag,'  Months'
       CALL MsgOne(h,trim(line))
    END IF
    IF (soilmlag < 0) THEN
       WRITE(line, FMT=336) ifsoilm, soilmlag
       CALL FatalError(h//trim(line))
    END IF
    !READ (UNIT=nfslm,REC=2) var4
    CLOSE (UNIT=nfslm)

    !DO j=1,jMax
    !   DO i=1,imax
    !      rlsm(i,j)=var4(i,j)
    !   END DO
    !END DO
    DEALLOCATE(var4)
336 FORMAT(' FOR IFSOILM=',I5,' SOILMLAG MUST BE SET NONNEGATIVE.  NOT ',G12.5)  
  END SUBROUTINE  LandSoilmMask8

  SUBROUTINE daylag (labeli,labelf,dlag,intsst)

    CHARACTER (LEN=10 ) :: labeli
    CHARACTER (LEN=10 ) :: labelf
    REAL      (KIND=r8) :: dlag
    REAL      (KIND=r8) :: xday
    REAL      (KIND=r8) :: xday2
    REAL      (KIND=r8) :: yday
    REAL      (KIND=r8) :: tod=0.0_r8
    INTEGER , INTENT(IN   )   :: intsst
    INTEGER :: yi,mi,di,hi,yf,mf,df,hf,ndy,y,n,ndi,ndf
    INTEGER :: xyi,xmi,xdi,xhi,LenYearbyDay,nday2y
    INTEGER, DIMENSION (12) :: ndm = &
         (/31,28,31,30,31,30,31,31,30,31,30,31/)
    INTEGER, DIMENSION (12) :: ndmi = &
         (/31,28,31,30,31,30,31,31,30,31,30,31/)
    INTEGER, DIMENSION (12) :: ndmf = &
         (/31,28,31,30,31,30,31,31,30,31,30,31/)

    READ (labeli(1:4), '(i4)') yi
    READ (labeli(5:6), '(i2)') mi
    READ (labeli(7:8), '(i2)') di
    READ (labeli(9:10),'(i2)') hi
    READ (labelf(1:4), '(i4)') yf
    READ (labelf(5:6), '(i2)') mf
    READ (labelf(7:8), '(i2)') df
    READ (labelf(9:10),'(i2)') hf
    CALL jull(yi,mi,di,hi,tod,xday2,nday2y,LenYearbyDay)
    IF (LenYearbyDay == 366) ndmi(2)=29
    CALL jull(yf,mf,df,hf,tod,xday2,nday2y,LenYearbyDay)
    IF (LenYearbyDay == 366) ndmf(2)=29

    IF (intsst > 0) THEN

       ndy=0
       DO y=yi+1,yf-1
          DO n=1,12
             ndy=ndy+ndm(n)
          END DO
          xyi=y
          xmi=12
          xdi=31
          xhi=0 
          CALL jull(xyi,xmi,xdi,xhi,tod,xday2,nday2y,LenYearbyDay)
          IF (LenYearbyDay == 366) ndy=ndy+1
       END DO

       ndi=di
       DO n=1,mi-1
          ndi=ndi+ndmi(n)
       END DO
       ndf=df
       DO n=1,mf-1
          ndf=ndf+ndmf(n)
       END DO

       IF (yf .EQ. yi) THEN
          dlag=REAL(ndf-ndi,r8)+REAL(hf-hi,r8)/24.0_r8
       ELSE IF(yf .GT. yi) THEN
          ndi=365-ndi
          IF (ndmi(2) .EQ. 29) ndi=ndi+1
          dlag=REAL(ndf+ndi,r8)+REAL(hf-hi,r8)/24.0_r8+ndy
       ELSE
          dlag=-1.0_r8
       END IF

    ELSE

       IF (mf >= mi) THEN
          dlag=mf-mi+12*(yf-yi)
       ELSE
          dlag=12+mf-mi+12*(yf-yi-1)
       END IF
       xyi=yf
       xmi=12
       xdi=31
       xhi=0 
       CALL jull(xyi,xmi,xdi,xhi,tod,xday2,nday2y,LenYearbyDay)
       IF (LenYearbyDay == 366) ndmf(2)=29
       xday=REAL(df,r8)+REAL(hf,r8)/24.0_r8
       yday=1.0_r8+REAL(ndmf(mf),r8)/2.0_r8
       IF (xday <= yday) dlag=dlag-1.0_r8

    END IF

  END SUBROUTINE daylag


  SUBROUTINE ReadVar4(nfvar,irec,var,fclose)
    INTEGER(KIND=i4) , INTENT(in   ) :: nfvar
    INTEGER(KIND=i4) , INTENT(IN   ) :: irec
    REAL   (KIND=r4) , INTENT(out  ) :: var (:,:)
    INTEGER           , INTENT(IN   ) :: fclose
    REAL   (KIND=r4)  :: var4(iMax,jMax)
    INTEGER(KIND=i4)  :: i
    INTEGER(KIND=i4)  :: j
     
    READ(UNIT=nfvar,rec=irec) var4
    DO j=1,jMax
       DO i=1,iMax
          var(i,j)=var4(i,j)
       END DO
    END DO
    IF(fclose == 0)CLOSE(UNIT=nfvar,STATUS='KEEP')
  END SUBROUTINE ReadVar4
  SUBROUTINE ReadVar8(nfvar,irec,var,fclose)
    INTEGER           , INTENT(in   ) :: nfvar
    INTEGER           , INTENT(IN   ) :: irec
    REAL     (KIND=r8), INTENT(out  ) :: var (:,:)  
    INTEGER           , INTENT(IN   ) :: fclose  
    REAL     (KIND=r4) :: var8(iMax,jMax)
    INTEGER            :: i
    INTEGER            :: j 
    READ(UNIT=nfvar,rec=irec) var8
    DO j=1,jMax
       DO i=1,iMax
          var(i,j)=REAL(var8(i,j),r8)     
       END DO
    END DO
    IF(fclose == 0) CLOSE(UNIT=nfvar,STATUS='KEEP') 
  END SUBROUTINE ReadVar8

  SUBROUTINE ReadAlb4(n,irec,field)
    INTEGER(KIND=i4) , INTENT(in   ) :: n
    INTEGER(KIND=i4) , INTENT(IN   ) :: irec
    REAL   (KIND=r4) , INTENT(out  ) :: field  (:,:)
    REAL   (KIND=r8) :: raux3(SIZE(field,1),SIZE(field,2))
    INTEGER(KIND=i4) :: i,j
    INTEGER(KIND=i4) :: idim,jdim
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadAlb4)**"
    idim=SIZE(field,1) ;jdim=SIZE(field,2)
    READ(UNIT=n,rec=irec) raux3
    DO j=1,jdim
      DO i=1,idim
        field(i,j) = REAL(raux3(i,j), r4)
      END DO
    END DO  
  END SUBROUTINE ReadAlb4
  SUBROUTINE ReadAlb8(n,irec,field)
    INTEGER , INTENT(in   ) :: n
    INTEGER , INTENT(IN   ) :: irec
    REAL     (KIND=r8) , INTENT(out  ) :: field  (:,:)
    REAL     (KIND=r8) :: raux3(SIZE(field,1),SIZE(field,2))
    INTEGER(KIND=i4) :: i,j
    INTEGER(KIND=i4) :: idim,jdim
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadAlb8)**"   
    idim=SIZE(field,1) ;jdim=SIZE(field,2)
    READ(UNIT=n,rec=irec) raux3
    DO j=1,jdim
      DO i=1,idim
        field(i,j) = REAL(raux3(i,j), r8)
      END DO
    END DO  
  END SUBROUTINE ReadAlb8

  SUBROUTINE ReadSST4(n,irec,field)
    INTEGER , INTENT(in   ) :: n
    INTEGER , INTENT(IN   ) :: irec
    REAL     (KIND=r4)         , INTENT(out  ) :: field  (:,:)
    REAL     (KIND=r4)                         :: raux3(SIZE(field,1),SIZE(field,2))
    INTEGER(KIND=i4)                         :: i,j
    INTEGER(KIND=i4)                         :: idim,jdim
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadSST4)**" 
    idim=SIZE(field,1) ;jdim=SIZE(field,2)
    READ (UNIT=n, REC=irec) raux3
    DO j=1,jdim
      DO i=1,idim
        field(i,j) = REAL(raux3(i,j), r4)
      END DO
    END DO  
  END SUBROUTINE ReadSST4
  SUBROUTINE ReadSST8(n,irec,field)
    INTEGER , INTENT(IN   ) :: n
    INTEGER , INTENT(IN   ) :: irec
    REAL   (KIND=r8)          , INTENT(OUT  ) :: field(:,:)
    REAL   (KIND=r4)                          :: raux3(SIZE(field,1),SIZE(field,2))
    INTEGER(KIND=i4)                          :: i,j
    INTEGER(KIND=i4)                          :: idim,jdim
    CHARACTER(LEN=*), PARAMETER :: h="**(ReadSST8)**"   
    idim=SIZE(field,1) ;jdim=SIZE(field,2)
    READ (UNIT=n, REC=irec) raux3
    DO j=1,jdim
      DO i=1,idim
        field(i,j) = REAL(raux3(i,j), r8)
      END DO
    END DO  
  END SUBROUTINE ReadSST8

  SUBROUTINE ReadSST4Rec(n,field,irec)
    INTEGER           , INTENT(in   ) :: n
    INTEGER           , INTENT(in   ) :: irec
    REAL     (KIND=r4)         , INTENT(out  ) :: field  (:,:)
    REAL     (KIND=r4)                         :: raux3(SIZE(field,1)*SIZE(field,2))
    INTEGER(KIND=i4)                         :: i,j,ij
    INTEGER(KIND=i4)                         :: idim,jdim
    idim=SIZE(field,1) ;jdim=SIZE(field,2)
    READ(UNIT=n,rec=irec) raux3
    ij=0
    DO j=1,jdim
      DO i=1,idim
        ij=ij+1
        field(i,j) = REAL(raux3(ij), r4)
      END DO
    END DO
  END SUBROUTINE ReadSST4Rec
  SUBROUTINE ReadSST8Rec(n,field,irec)
    INTEGER          , INTENT(IN   ) :: n
    INTEGER          , INTENT(in   ) :: irec
    REAL   (KIND=r8)          , INTENT(OUT  ) :: field(:,:)
    REAL   (KIND=r4)                          :: raux3(SIZE(field,1)*SIZE(field,2))
    INTEGER(KIND=i4)                          :: i,j,ij
    INTEGER(KIND=i4)                          :: idim,jdim
    idim=SIZE(field,1) ;jdim=SIZE(field,2)
    READ(UNIT=n,rec=irec)raux3
    ij=0
    DO j=1,jdim
      DO i=1,idim
        ij=ij+1    
        field(i,j) = REAL(raux3(ij), r8)
      END DO
    END DO  
  END SUBROUTINE ReadSST8Rec

  SUBROUTINE ReadSLM4(n,irec,field)
    INTEGER(KIND=i4)  , INTENT(in   ) :: n
    INTEGER           , INTENT(in   ) :: irec
    REAL     (KIND=r4), INTENT(out  ) :: field  (:,:) 
    REAL     (KIND=r4)                :: raux3(SIZE(field,1),SIZE(field,2))
    INTEGER(KIND=i4)                  :: i,j
    INTEGER(KIND=i4)                  :: idim,jdim
    idim=SIZE(field,1) ;jdim=SIZE(field,2)
    READ(UNIT=n,rec=irec) raux3
    DO j=1,jdim
      DO i=1,idim
        field(i,j) = raux3(i,j)
      END DO
    END DO  
  END SUBROUTINE ReadSLM4
  SUBROUTINE ReadSLM8(n,irec,field)
    INTEGER           , INTENT(in   ) :: n
    INTEGER           , INTENT(in   ) :: irec
    REAL     (KIND=r8), INTENT(out  ) :: field  (:,:)
    REAL     (KIND=r4)                :: raux3(SIZE(field,1),SIZE(field,2))
    INTEGER(KIND=i4)                  :: i,j
    INTEGER(KIND=i4)                  :: idim,jdim
    idim=SIZE(field,1) ;jdim=SIZE(field,2)
    READ(UNIT=n,rec=irec) raux3
    DO j=1,jdim
      DO i=1,idim
        field(i,j) = raux3(i,j)
      END DO
    END DO
  END SUBROUTINE ReadSLM8
 SUBROUTINE ReadSLM43D(n,irec,field)
    INTEGER(KIND=i4)  , INTENT(in   ) :: n
    INTEGER           , INTENT(in   ) :: irec
    REAL     (KIND=r4), INTENT(out  ) :: field  (:,:,:) 
    REAL     (KIND=r4)                :: raux3(SIZE(field,1),SIZE(field,2),SIZE(field,3))
    INTEGER(KIND=i4)                  :: i,j,k
    INTEGER(KIND=i4)                  :: idim,jdim,kdim
    idim=SIZE(field,1) ;jdim=SIZE(field,2);kdim=SIZE(field,3)
    READ(UNIT=n,rec=irec) raux3
    DO k=1,kdim
       DO j=1,jdim
          DO i=1,idim
             field(i,j,k) = raux3(i,j,k)
          END DO
      END DO
    END DO  
  END SUBROUTINE ReadSLM43D
  SUBROUTINE ReadSLM83D(n,irec,field)
    INTEGER           , INTENT(in   ) :: n
    INTEGER           , INTENT(in   ) :: irec
    REAL     (KIND=r8), INTENT(out  ) :: field  (:,:,:)
    REAL     (KIND=r4)                :: raux3(SIZE(field,1),SIZE(field,2),SIZE(field,3))
    INTEGER(KIND=i4)                  :: i,j,k
    INTEGER(KIND=i4)                  :: idim,jdim,kdim
    idim=SIZE(field,1) ;jdim=SIZE(field,2);kdim=SIZE(field,3)
    READ(UNIT=n,rec=irec) raux3
    DO k=1,kdim
       DO j=1,jdim
          DO i=1,idim
             field(i,j,k) = raux3(i,j,k)
          END DO
      END DO
    END DO  
  END SUBROUTINE ReadSLM83D


  SUBROUTINE ReadSNW4(n,irec,field)
    INTEGER           , INTENT(IN   ) :: n
    INTEGER           , INTENT(IN   ) :: irec
    REAL     (KIND=r4)         , INTENT(out  ) :: field  (:,:)
    REAL     (KIND=r4)                         :: raux3(SIZE(field,1),SIZE(field,2))
    INTEGER(KIND=i4)                         :: i,j
    INTEGER(KIND=i4)                         :: idim,jdim
    idim=SIZE(field,1) ;jdim=SIZE(field,2)
    READ(UNIT=n,rec=irec) raux3
    DO j=1,jdim
      DO i=1,idim
        field(i,j) = REAL(raux3(i,j), r4)
      END DO
    END DO  
  END SUBROUTINE ReadSNW4
  SUBROUTINE ReadSNW8(n,irec,field)
    INTEGER         , INTENT(IN   ) :: n
    INTEGER         , INTENT(IN   ) :: irec
    REAL   (KIND=r8), INTENT(OUT  ) :: field(:,:)
    REAL   (KIND=r4)                          :: raux3(SIZE(field,1),SIZE(field,2))
    INTEGER(KIND=i4)                          :: i,j
    INTEGER(KIND=i4)                          :: idim,jdim
    idim=SIZE(field,1) ;jdim=SIZE(field,2)
    READ(UNIT=n,rec=irec)raux3
    DO j=1,jdim
      DO i=1,idim
         field(i,j) = REAL(raux3(i,j), r8)
      END DO
    END DO
  END SUBROUTINE ReadSNW8

!tar begin 
!Read climate aerosol parameters r8

  SUBROUTINE Read_Aeros(n,irec,field)
    INTEGER           , INTENT(in   ) :: n
    INTEGER           , INTENT(in   ) :: irec
    REAL     (KIND=r8), INTENT(out  ) :: field  (:,:,:)
    REAL     (KIND=r8)                :: raux3(SIZE(field,1),SIZE(field,2),SIZE(field,3))
    INTEGER(KIND=i4)                  :: i,j,k
    INTEGER(KIND=i4)                  :: idim,jdim,kdim
    idim=SIZE(field,1) ;jdim=SIZE(field,2);kdim=SIZE(field,3)
    READ(UNIT=n,rec=irec) raux3
    DO k=1,kdim
       DO j=1,jdim
          DO i=1,idim
             field(i,j,k) = raux3(i,j,k)
          END DO
      END DO
    END DO  
  END SUBROUTINE Read_Aeros

!tar end
  ! The Ozone files are written in a way that they can be read from GrADS.
  ! Therefore,  the order of the varibles inside the binary had to be x, y, z
  !   what is different from the global model x, z, y order.
  !   Moreover, the file is open in GrADS with a yrev option. So y=1 
  !   means the north pole as exepected by the global model.
  SUBROUTINE ReadOzone8(n,field,irec)
     INTEGER          , INTENT(IN   ) :: n,irec
     REAL   (KIND=r8) , INTENT(OUT  ) :: field(:,:,:)
     REAL   (KIND=r4)                 :: raux3(SIZE(field,1),SIZE(field,3),SIZE(field,2))
     CHARACTER(LEN=*), PARAMETER :: h="**(ReadOzone8)**"   
     INTEGER :: i,j,k,im,jm,km

     READ(UNIT=n,REC=irec) raux3
     ! input field is (i,k,j)
     im=size(field,1)
     jm=size(field,3)
     km=size(field,2)
     do j=1,jm
        do k=1,km
           do i=1,im
              field(i,k,j) = REAL(raux3(i,j,k), r8)
           enddo
        enddo
     enddo

  END SUBROUTINE ReadOzone8
  ! The tracer files are written in a way that they can be read from GrADS.
  ! Therefore,  the order of the varibles inside the binary had to be x, y, z
  !   what is different from the global model x, z, y order.
  !   Moreover, the file is open in GrADS with a yrev option. So y=1 
  !   means the north pole as exepected by the global model.
  SUBROUTINE ReadTracer8(n,field,irec)
     INTEGER          , INTENT(IN   ) :: n,irec
     REAL   (KIND=r8) , INTENT(OUT  ) :: field(:,:,:)
     REAL   (KIND=r4)                 :: raux3(SIZE(field,1),SIZE(field,3))
     CHARACTER(LEN=*), PARAMETER :: h="**(ReadOzone8)**"   
     INTEGER :: i,j,k,im,jm,km

     !READ(UNIT=n,REC=irec) raux3
     ! input field is (i,k,j)
     im=size(field,1)
     jm=size(field,3)
     km=size(field,2)
     do k=1,km
        READ(UNIT=n,REC=k) raux3
        PRINT*,'pkubota',k,irec
        PRINT*,raux3
        do j=1,jm
           do i=1,im
              field(i,k,j) = REAL(raux3(i,j),kind= r8)
              IF(field(i,k,j).LT.0.0e0_r8)field(i,k,j) = 1.0e-12_r8
           enddo
        enddo
     enddo

  END SUBROUTINE ReadTracer8
  SUBROUTINE ReadTracer8s(n,field)
     INTEGER          , INTENT(IN   ) :: n
     REAL   (KIND=r8) , INTENT(OUT  ) :: field(:,:,:)
     REAL   (KIND=r4)                 :: raux3(SIZE(field,1),size(field,3))
     CHARACTER(LEN=*), PARAMETER :: h="**(ReadOzone8)**"
     INTEGER :: i,j,k,im,jm,km

     !READ(UNIT=n,REC=irec) raux3
     ! input field is (i,k,j)
     im=size(field,1)
     jm=size(field,3)
     km=size(field,2)
     do k=1,km
        READ(UNIT=n) raux3
        do j=1,jm
           do i=1,im
              field(i,k,j) = REAL(raux3(i,j),kind= r8)
              IF(field(i,k,j).LT.0.0e0_r8)field(i,k,j) = 1.0e-12_r8
           enddo
        enddo
     enddo

  END SUBROUTINE ReadTracer8s
  SUBROUTINE ReadNFTGZ4(n,irec,field1,field2,field3)
    INTEGER, INTENT(in   ) :: n 
    INTEGER, INTENT(in   ) :: irec
    REAL   (KIND=r8), INTENT(out  ) :: field1 (:,:)
    REAL   (KIND=r8), INTENT(out  ) :: field2 (:,:)
    REAL   (KIND=r8), INTENT(out  ) :: field3 (:,:)
    REAL   (KIND=r4)                :: raux1(SIZE(field1,1),SIZE(field1,2))
    REAL   (KIND=r4)                :: raux2(SIZE(field2,1),SIZE(field2,2))
    REAL   (KIND=r4)                :: raux3(SIZE(field3,1),SIZE(field3,2))

    READ (UNIT=n, Rec=irec) raux1
    READ (UNIT=n, Rec=irec+1) raux2
    READ (UNIT=n, Rec=irec+2) raux3

    field1 = REAL(raux1, r8)
    field2 = REAL(raux2, r8)
    field3 = REAL(raux3, r8)

    !REWIND(UNIT=n)
  END SUBROUTINE ReadNFTGZ4
  SUBROUTINE ReadNFTGZ8(n,irec,field1,field2,field3)
    INTEGER(KIND=i4), INTENT(in   ) :: n 
    INTEGER         , INTENT(in   ) :: irec
    REAL   (KIND=r4), INTENT(out  ) :: field1 (:,:)
    REAL   (KIND=r4), INTENT(out  ) :: field2 (:,:)
    REAL   (KIND=r4), INTENT(out  ) :: field3 (:,:)
    REAL   (KIND=r4)                :: raux1(SIZE(field1,1),SIZE(field1,2))
    REAL   (KIND=r4)                :: raux2(SIZE(field2,1),SIZE(field2,2))
    REAL   (KIND=r4)                :: raux3(SIZE(field3,1),SIZE(field3,2))

    READ (UNIT=n, Rec=irec ) raux1
    READ (UNIT=n, Rec=irec+1) raux2
    READ (UNIT=n, Rec=irec+2) raux3

    field1 = REAL(raux1, r4)
    field2 = REAL(raux2, r4)
    field3 = REAL(raux3, r4)

    REWIND(UNIT=n)
  END SUBROUTINE ReadNFTGZ8

  SUBROUTINE WriteDiagHead4(n, ifday, tod)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    INTEGER(KIND=i4), INTENT(IN)  :: ifday
    REAL   (KIND=r4), INTENT(IN)  :: tod
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteProgHead4)**"
    WRITE(UNIT=n)ifday, tod
  END SUBROUTINE WriteDiagHead4
  SUBROUTINE WriteDiagHead8(n, ifday, tod)
    INTEGER, INTENT(IN)  :: n
    INTEGER, INTENT(IN)  :: ifday
    REAL   (KIND=r8), INTENT(IN)  :: tod
    INTEGER(KIND=i4) :: iaux(2)
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteProgHead8)**"
    iaux(  1 ) = INT (ifday , i4)
    iaux(  2 ) = INT (tod   , i4)
    WRITE(UNIT=n)iaux
  END SUBROUTINE WriteDiagHead8

  SUBROUTINE WriteGrdH4(n,field1,field2)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    REAL   (KIND=r4), INTENT(IN)  :: field1(:)
    REAL   (KIND=r4), INTENT(IN)  :: field2(:,:)
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteGrdH4)**"
    WRITE(UNIT=n)field1
    WRITE(UNIT=n)field2
  END SUBROUTINE WriteGrdH4
  SUBROUTINE WriteGrdH8(n,field1,field2)
    INTEGER, INTENT(IN)  :: n
    REAL   (KIND=r8), INTENT(IN)  :: field1(:)
    REAL   (KIND=r8), INTENT(IN)  :: field2(:,:)
    REAL   (KIND=r4) :: raux1(SIZE(field1,1))
    REAL   (KIND=r4) :: raux2(SIZE(field2,1),SIZE(field2,2))   
    CHARACTER(LEN=*), PARAMETER :: h="**(WriteGrdH8)**"
    raux1 = REAL(field1, r4)
    raux2 = REAL(field2, r4)
    WRITE(UNIT=n)raux1
    WRITE(UNIT=n)raux2
  END SUBROUTINE WriteGrdH8  
  SUBROUTINE WrTopoGrdH4(n,field2)
    INTEGER(KIND=i4), INTENT(IN)  :: n
    REAL   (KIND=r4), INTENT(IN)  :: field2(:)
    CHARACTER(LEN=*), PARAMETER :: h="**(WrTopoGrdH4)**"
    WRITE(UNIT=n)field2
  END SUBROUTINE WrTopoGrdH4
  SUBROUTINE WrTopoGrdH8(n,field2)
    INTEGER, INTENT(IN)  :: n
    REAL   (KIND=r8), INTENT(IN)  :: field2(:)
    REAL   (KIND=r4) :: raux2(SIZE(field2,1))   
    CHARACTER(LEN=*), PARAMETER :: h="**(WrTopoGrdH8)**"
    raux2 = REAL(field2, r4)
    WRITE(UNIT=n)raux2
  END SUBROUTINE WrTopoGrdH8
END MODULE IOLowLevel
