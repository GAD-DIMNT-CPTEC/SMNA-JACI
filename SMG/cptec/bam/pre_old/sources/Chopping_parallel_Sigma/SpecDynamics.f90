!
!  $Author: pkubota $
!  $Date: 2011/04/07 16:00:31 $
!  $Revision: 1.15 $
!
MODULE SpecDynamics

  USE Utils, ONLY : &
       Epslon                ! intent(in)

  USE Parallelism, ONLY : &
       mygroup_four, & ! intent(in)
       maxnodes,     & ! intent(in)
       myid

  USE InputParameters,    ONLY : &
       EMRad1,       & ! intent(in)
       EMRad,        & ! intent(in)
       r8              ! intent(in)

  USE Sizes,  ONLY : &
       mMax,         & ! intent(in)
       nMax,         & ! intent(in)
       nExtMax,      & ! intent(in)
       mnMax,        & ! intent(in)
       mnExtMax,     & ! intent(in)
       mnExtMap,     & ! intent(in)
       mnMap,        & ! intent(in)
       nMap,         & ! intent(in)
       mMap,         & ! intent(in)
       lm2m,         & ! intent(in)
       kMax,         & ! intent(in)
       kMaxloc,      & ! intent(in)
       myfirstlev,   & ! intent(in)
       mylastlev,    & ! intent(in)
       mymMax,       & ! intent(in)
       mymnMax,      & ! intent(in)
       mymnExtMax,   & ! intent(in)
       mymnExtMap,   & ! intent(in)
       mymnMap,      & ! intent(in)
       mynMap,       & ! intent(in)
       mymextMap,    & ! intent(in)
       mynextMap,    & ! intent(in)
       mymMap,       & ! intent(in)
       haveM1,       & ! intent(in)
       havesurf,     & ! intent(in)
       nodehasm,     & ! intent(in)
       ngroups_four, & ! intent(in)
       nlevperg_four,& ! intent(in)
       rpi,          & ! intent(in)
       del,          & ! intent(in)
       ci              ! intent(in) 

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: InitDztouv
  PUBLIC :: dztouv
  PUBLIC :: InitUvtodz
  PUBLIC :: uvtodz

  REAL(KIND=r8) :: er, eriv
  !  Module exports two routines:
  !     InitDZtoUV:        Should be invoked once, before any other module 
  !                        routine; sets up local constants and mappings;
  !     DZtoUV:            Velocity fields from Divergence and Vorticity; 
  !                        use values computed by InitDZtoUV

  !  Module require values from modules Sizes, AssocLegFunc and Constants

  REAL(KIND=r8), ALLOCATABLE :: alfa_dz(:)        ! er*Epslon(m,n)/n  for m<n<=nExtMax; 0 otherwise
  REAL(KIND=r8), ALLOCATABLE :: alfa_dzNp1(:)     ! alfa_dz(m,n+1) 
  REAL(KIND=r8), ALLOCATABLE :: beta_dz(:)        ! m*er/(n*(n+1)) for m/=0 and m<=n<=nMax;
  !                                          er/(n+1)     for m=0;



  ! Observe, in the relation to be computed, that u and v are defined for 
  ! 1<=mn<=mnExtMax, while Div and Vor for 1<=mn<=mnMax. Consequently, a 
  ! mapping from 1:mnExtMax to 1:mnMax
  ! has to be computed.

  ! This mapping will have faults, since there is no Div or Vor at (m,nExtMax).

  ! Furthermore, the relation requires mappings from (m,nExt) to (m,n), 
  ! (m,n+1) and (m,n-1)

  ! Mapping function mnp1_dz(1:2*mnExtMax) gives index of (m,nExt) on (m,n+1). 
  ! It is faulty on (*,nMax:nExtMax). Since it is only used by the last term,
  ! faulty values have to be multipied by 0 (on alfa_dzNp1)

  ! Mapping function mnm1_dz(1:2*mnExtMax) gives index of (m,nExt) on (m,n-1). 
  ! It is faulty on (m,m) for all m. Since it is only used by the second term,
  ! faulty values have to be multipied by 0 (on alfa_dz)

  ! Mapping function mnir_dz(1:2*mnExtMax) gives index of (m,nExt) on (m,n-1) and
  ! multiplies by i (trading imaginary by real and correcting sign). It is
  ! faulty on (m,nExtMax). To correct the fault, beta_dz(m,nExtMax) is set to 0.

  INTEGER, ALLOCATABLE :: mnir_dz(:)
  INTEGER, ALLOCATABLE :: mnm1_dz(:)
  INTEGER, ALLOCATABLE :: mnp1_dz(:)


  ! DIVERGENCE AND VORTICITY FROM VELOCITY FIELDS
  !
  ! Implements the following relations:
  !
  !    m              m        m       m        m         m   m       m   m
  ! Div   =CMPLX(-Alfa * Imag(U ), Alfa * Real(U )) + Beta * V  - Gama * V
  !    n              n        n       n        n        n+1 n+1      n  n-1
  !
  !    m              m        m       m        m         m   m       m   m
  ! Vor   =CMPLX(-Alfa * Imag(V ), Alfa * Real(V )) + Beta * U  - Gama * U
  !    n              n        n       n        n        n+1 n+1      n  n-1
  !
  ! for 0<=m<=mMax, m<=n<=nMax, where
  !
  !  m   m
  ! U = V = 0 for n < m
  !  n   n

  !  Module exports two routines:
  !     InitUvtodz:        Should be invoked once, before any other module 
  !                        routine; sets up local constants and mappings;
  !     Uvtodz:            Divergence and Vorticity from Velocity fields; 
  !                        use values computed by InitUvtodz

  !  Module require values from modules Sizes, AssocLegFunc and Constants

  REAL(KIND=r8),    ALLOCATABLE :: alfa_uv(:)        ! eriv*m 
  !                                           for m<=n<=nMax
  REAL(KIND=r8),    ALLOCATABLE :: beta_uv(:)        ! eriv*n*Epslon(m,n+1) 
  !                                           for m<=n<=nMax;
  REAL(KIND=r8),    ALLOCATABLE :: gama_uv(:)        ! eriv*(n+1)*Epslon(m,n) 
  !                                           for m<=n<=nMax;

  ! Observe, in the relation to be computed, that u and v are defined for 
  ! 1<=mn<=mnExtMax, while Div and Vor for 1<=mn<=mnMax. Consequently, a 
  ! mapping from 1:mnMax to 1:mnExtMax has to be computed.

  ! In fact, the relation requires 3 mappings:
  !    1) (m,n) to (m,nExt)    implemented by index array mnir(mn);
  !    2) (m,n) to (m,nExt+1)  implemented by index array mnp1(mn);
  !    3) (m,n) to (m,nExt-1)  implemented by index array mnm1(mn);
  ! for m<=n<=nMax

  ! Mapping functions (1) and (2) are easily computed; mapping function (3)
  ! will be faulty for n=m; since it is used only at the last term,
  ! faulty values have to be multipied by 0 (on gama_uv(m,m))

  ! Mapping function mnir(1:2*mnMax) gives index of (m,n) on (m,nExt) and
  ! multiplies by i (trading immaginary by real). Other mappings keep in
  ! place the real and immaginary components.

  INTEGER, ALLOCATABLE :: mnir_uv(:)
  INTEGER, ALLOCATABLE :: mnm1_uv(:)
  INTEGER, ALLOCATABLE :: mnp1_uv(:)

  INTEGER :: kGlob
  REAL(KIND=r8)    :: alphaGlob
  REAL(KIND=r8)    :: betaGlob
  INTEGER, ALLOCATABLE :: ncrit(:)

CONTAINS


  ! InitDZtoUV: Mapping Functions and Local Constants


  SUBROUTINE InitDZtoUV()
    INTEGER :: m, mglob, n, mn, mn2, indir, indnp1, indnm1
    REAL(KIND=r8) :: aux

    ! mapping mnir_dz
    er = EMRad
    eriv = EMRad1
    ALLOCATE(mnir_dz(2*mymnExtMax))
    !$OMP PARALLEL DO PRIVATE(mglob,n,mn,indir)
    DO m = 1, mymMax
       mglob = lm2m(m)
       DO n = mglob, nMax
          mn = mymnExtMap(m,n)
          indir = mymnMap(m,n)
          mnir_dz(2*mn-1) = 2*indir
          mnir_dz(2*mn  ) = 2*indir-1
       END DO
       mn = mymnExtMap(m,nExtMax)
       mnir_dz(2*mn-1) = 1     ! faulty mapping # 1
       mnir_dz(2*mn  ) = 1     ! faulty mapping # 1
    END DO
    !$OMP END PARALLEL DO

    ! mapping mnm1_dz

    ALLOCATE(mnm1_dz  (2*mymnExtMax))
    !$OMP PARALLEL DO PRIVATE(mglob,n,mn,indnm1)
    DO m = 1, mymMax
       mglob = lm2m(m)
       mn = mymnExtMap(m,mglob)
       indnm1 = mymnMap(m,mglob)       
       mnm1_dz(2*mn-1) = 1    ! faulty mapping # 2
       mnm1_dz(2*mn  ) = 1    ! faulty mapping # 2
       DO n = mglob+1, nExtMax
          mn = mymnExtMap(m,n)
          indnm1 = mymnMap(m,n-1)
          mnm1_dz(2*mn-1) = 2*indnm1-1
          mnm1_dz(2*mn  ) = 2*indnm1  
       END DO
    END DO
    !$OMP END PARALLEL DO

    ! mapping mnp1_dz

    ALLOCATE(mnp1_dz  (2*mymnExtMax))
    mnp1_dz = 0
    !$OMP PARALLEL DO PRIVATE(mglob,n,mn,indnp1)
    DO m = 1, mymMax
       mglob = lm2m(m)
       DO n = mglob, nMax-1
          mn = mymnExtMap(m,n)
          indnp1 = mymnMap(m,n+1)
          mnp1_dz(2*mn-1) = 2*indnp1-1
          mnp1_dz(2*mn  ) = 2*indnp1  
       END DO
       DO n = nMax, nExtMax
          mn = mymnExtMap(m,n)
          mnp1_dz(2*mn-1) = 1 ! faulty mapping # 3
          mnp1_dz(2*mn  ) = 1 ! faulty mapping # 3
       END DO
    END DO
    !$OMP END PARALLEL DO

    ! constant beta_dz

    ALLOCATE(beta_dz(2*mymnExtMax))
    !$OMP PARALLEL DO PRIVATE(mglob,n,mn,aux)
    DO m = 1, mymMax
       mglob = lm2m(m)
       aux = er/REAL(mglob,r8)
       mn = mymnExtMap(m,mglob)
       beta_dz(2*mn-1) = aux
       beta_dz(2*mn  ) = -aux
       DO n = mglob+1, nMax
          aux = REAL(mglob-1,r8)*er/REAL((n-1)*n,r8)
          mn = mymnExtMap(m,n)
          beta_dz(2*mn-1) = aux
          beta_dz(2*mn  ) = -aux
       END DO
       mn = mymnExtMap(m,nExtMax)
       beta_dz(2*mn-1) = 0.0_r8           ! corrects faulty mapping # 1
       beta_dz(2*mn  ) = 0.0_r8           ! corrects faulty mapping # 1
    END DO
    !$OMP END PARALLEL DO

    ! constant alfa_dz

    ALLOCATE(alfa_dz    (2*mymnExtMax))
    !$OMP PARALLEL DO PRIVATE(mglob,n,mn,aux)
    DO m = 1, mymMax
       mglob = lm2m(m)
       mn = mymnExtMap(m,mglob)
       alfa_dz(2*mn-1) = 0.0_r8           ! corrects faulty mapping # 2
       alfa_dz(2*mn  ) = 0.0_r8           ! corrects faulty mapping # 2
       DO n = mglob+1, nExtMax
          mn = mymnExtMap(m,n)
          aux = er * Epslon(mn) / REAL(n-1,r8)
          alfa_dz(2*mn-1) = aux
          alfa_dz(2*mn  ) = aux
       END DO
    END DO
    !$OMP END PARALLEL DO

    ! constant alfa_dz mapped to n-1

    ALLOCATE(alfa_dzNp1 (2*mymnExtMax))
    !$OMP PARALLEL DO PRIVATE(mglob,n,mn,mn2,aux)
    DO m = 1, mymMax
       mglob = lm2m(m)
       DO n = mglob, nMax - 1
          mn  = mymnExtMap(m,n)
          mn2 = mymnExtMap(m,n+1)
          aux  = er * Epslon(mn2) / REAL(n,r8)
          alfa_dzNp1(2*mn-1) = aux
          alfa_dzNp1(2*mn  ) = aux
       END DO
       DO n = nMax, nExtMax
          mn = mymnExtMap(m,n)
          alfa_dzNp1(2*mn-1) = 0.0_r8     ! corrects faulty mapping # 3
          alfa_dzNp1(2*mn  ) = 0.0_r8     ! corrects faulty mapping # 3
       END DO
    END DO
    !$OMP END PARALLEL DO

  END SUBROUTINE InitDZtoUV


  ! DZtoUV: Velocity fields from divergence and vorticity
  !
  ! Implements the following relations:
  !
  !  m            m          m       m         m         m     m       m     m
  ! U  =CMPLX(Beta * Imag(Div), -Beta *Real(Div))  - alfa * Vor  + alfa * Vor
  !  n            n          n       n         n         n    n-1     n+1   n+1
  !
  !  m            m          m       m         m         m     m       m     m
  ! V  =CMPLX(Beta * Imag(Vor), -Beta *Real(Vor))  + alfa * Div  - alfa * Div
  !  n            n          n       n         n         n    n-1     n+1   n+1
  !
  ! for 0<=m<=mMax, m<=n<=nExtMax, where
  !
  !    m     m
  ! Div = Vor = 0 for n > nMax or n < m
  !    n     n



  SUBROUTINE DZtoUV(qdivp, qrotp, qup, qvp, mnRIExtFirst, mnRIExtLast)
    REAL(KIND=r8),    INTENT(IN ) :: qdivp(2*mymnMax,kMaxloc)
    REAL(KIND=r8),    INTENT(IN ) :: qrotp(2*mymnMax,kMaxloc)
    REAL(KIND=r8),    INTENT(OUT) :: qup(2*mymnExtMax,kMaxloc)
    REAL(KIND=r8),    INTENT(OUT) :: qvp(2*mymnExtMax,kMaxloc)
    INTEGER, INTENT(IN ) :: mnRIExtFirst
    INTEGER, INTENT(IN ) :: mnRIExtLast
    INTEGER :: mn, k   

    DO k = 1, kMaxloc
       DO mn = mnRIExtFirst, mnRIExtLast
          qup(mn,k) =                          - &
               alfa_dz   (mn) * qrotp(mnm1_dz(mn),k) + &
               alfa_dzNp1(mn) * qrotp(mnp1_dz(mn),k) + &
               beta_dz   (mn) * qdivp(mnir_dz(mn),k) 
          qvp(mn,k) =                          + &
               alfa_dz   (mn) * qdivp(mnm1_dz(mn),k) - &
               alfa_dzNp1(mn) * qdivp(mnp1_dz(mn),k) + &
               beta_dz   (mn) * qrotp(mnir_dz(mn),k) 
       END DO
    END DO
  END SUBROUTINE DZtoUV

  ! InitUvtodz: Mapping Functions and Local Constants


  SUBROUTINE InitUvtodz()
    INTEGER :: m, mglob, n, mn, mnExt

    ! mapping mnir_uv

    er = EMRad
    eriv = EMRad1
    ALLOCATE(mnir_uv(2*mymnMax))
    !$OMP PARALLEL DO PRIVATE(mglob,n,mn,mnExt)
    DO m = 1, mymMax
       mglob = lm2m(m)
       DO n = mglob, nMax
          mn    = mymnMap(m,n)
          mnExt = mymnExtMap(m,n)
          mnir_uv(2*mn-1) = 2*mnExt
          mnir_uv(2*mn  ) = 2*mnExt-1
       END DO
    END DO
    !$OMP END PARALLEL DO

    ! mapping mnm1_uv

    ALLOCATE(mnm1_uv(2*mymnMax))
    !$OMP PARALLEL DO PRIVATE(mglob,n,mn,mnExt)
    DO m = 1, mymMax
       mglob = lm2m(m)
       mn    = mymnMap(m,mglob)
       mnExt = mymnExtMap(m,mglob)
       mnm1_uv(2*mn-1) = 2*mnExt-1   ! faulty mapping
       mnm1_uv(2*mn  ) = 2*mnExt     ! faulty mapping
       DO n = mglob+1, nMax
          mn    = mymnMap(m,n)
          mnExt = mymnExtMap(m,n-1)
          mnm1_uv(2*mn-1) = 2*mnExt-1
          mnm1_uv(2*mn  ) = 2*mnExt
       END DO
    END DO
    !$OMP END PARALLEL DO

    ! mapping mnp1_uv

    ALLOCATE(mnp1_uv(2*mymnMax))
    !$OMP PARALLEL DO PRIVATE(mglob,n,mn,mnExt)
    DO m = 1, mymMax
       mglob = lm2m(m)
       DO n = mglob, nMax
          mn    = mymnMap(m,n)
          mnExt = mymnExtMap(m,n+1)
          mnp1_uv(2*mn-1) = 2*mnExt-1
          mnp1_uv(2*mn  ) = 2*mnExt
       END DO
    END DO
    !$OMP END PARALLEL DO

    ! constant alfa_uv

    ALLOCATE(alfa_uv(2*mymnMax))
    !$OMP PARALLEL DO PRIVATE(mglob,n,mn)
    DO m = 1, mymMax
       mglob = lm2m(m)
       DO n = mglob, nMax
          mn = mymnMap(m,n)
          alfa_uv(2*mn-1) = -REAL(mglob-1,r8)
          alfa_uv(2*mn  ) = REAL(mglob-1,r8)
       END DO
    END DO
    !$OMP END PARALLEL DO

    ! constant beta_uv

    ALLOCATE(beta_uv(2*mymnMax))
    !$OMP PARALLEL DO PRIVATE(mglob,n,mn,mnExt)
    DO m = 1, mymMax
       mglob = lm2m(m)
       DO n = mglob, nMax
          mn    = mymnMap(m,n)
          mnExt = mymnExtMap(m,n+1)
          beta_uv(2*mn-1) = REAL(n-1,r8)*Epslon(mnExt)
          beta_uv(2*mn  ) = REAL(n-1,r8)*Epslon(mnExt)
       END DO
    END DO
    !$OMP END PARALLEL DO

    ! constant gama_uv

    ALLOCATE(gama_uv(2*mymnExtMax))
    !$OMP PARALLEL DO PRIVATE(mglob,n,mn,mnExt)
    DO m = 1, mymMax
       mglob = lm2m(m)
       mn = mymnMap(m,mglob)
       gama_uv(2*mn-1) = 0.0_r8     ! corrects faulty mapping
       gama_uv(2*mn  ) = 0.0_r8     ! corrects faulty mapping
       DO n = mglob+1, nMax
          mn    = mymnMap(m,n)
          mnExt = mymnExtMap(m,n)
          gama_uv(2*mn-1) = REAL(n,r8) * Epslon(mnExt)
          gama_uv(2*mn  ) = REAL(n,r8) * Epslon(mnExt)
       END DO
    END DO
    !$OMP END PARALLEL DO

  END SUBROUTINE InitUvtodz


  ! Uvtodz: Divergence and Vorticity from Velocity fields 


  SUBROUTINE Uvtodz(qup, qvp, qdivt, qrott, mnRIFirst, mnRILast)
    REAL(KIND=r8),    INTENT(IN ) :: qup(2*mymnExtMax,kMaxloc)
    REAL(KIND=r8),    INTENT(IN ) :: qvp(2*mymnExtMax,kMaxloc)
    REAL(KIND=r8),    INTENT(OUT) :: qdivt(2*mymnMax,kMaxloc)
    REAL(KIND=r8),    INTENT(OUT) :: qrott(2*mymnMax,kMaxloc)
    INTEGER, INTENT(IN ) :: mnRIFirst
    INTEGER, INTENT(IN ) :: mnRILast
    INTEGER :: mn, k

    DO k = 1, kMaxloc
       DO mn = mnRIFirst, mnRILast
          qdivt(mn,k) = &
               alfa_uv(mn) * qup(mnir_uv(mn),k) + &
               beta_uv(mn) * qvp(mnp1_uv(mn),k) - &
               gama_uv(mn) * qvp(mnm1_uv(mn),k)
          qdivt(mn,k) = eriv * qdivt(mn,k)
          qrott(mn,k) = &
               alfa_uv(mn) * qvp(mnir_uv(mn),k) - &
               beta_uv(mn) * qup(mnp1_uv(mn),k) + &
               gama_uv(mn) * qup(mnm1_uv(mn),k)
          qrott(mn,k) = eriv * qrott(mn,k)
       END DO
    END DO
  END SUBROUTINE Uvtodz
END MODULE SpecDynamics
