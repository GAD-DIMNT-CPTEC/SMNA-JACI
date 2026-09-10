!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
MODULE FastFourierTransform

  USE InputParameters, ONLY : r8, nferr, Imax

  IMPLICIT NONE

  ! DESCRIPTION
  !
  ! Computes direct and inverse FFT transforms of sequences of 
  !          n real numbers for restricted values of n. Input
  !          data for the direct transform and output data
  !          from the inverse transform is supposed to represent
  !          a periodic function with period n+1.
  !
  ! n should be in the form: n=2*m
  ! where m =2**i * 3**j * 5**k, with i>=1 and j,k>=0
  !
  ! The direct FFT of a sequence of n real numbers fIn(1:n) is
  ! a sequence of n+1 real numbers fOut(1:n+1) such that
  !   fOut(2*k+1) = 1/m * SUM(fIn(j+1)*COS(PI*j*k/m)), k=0,...,m
  !   fOut(2*k+2) = 1/m * SUM(fIn(j+1)*SIN(PI*j*k/m)), k=0,...,m-1
  ! where both summations are taken for j=0,...,n-1
  ! 
  ! The inverse FFT of a sequence of n+1 real numbers fOut(1:n+1) is
  ! a sequence of n real numbers fIn(1:n) such that
  !   fIn(j+1) = 0.5_r8*(fOut(1)*COS(0.0r) + fOut(n+1)*COS(PI*j)) +
  !                 SUM(fOut(2*k+1)*COS(PI*j*k/m)) +
  !                 SUM(fOut(2*k+2)*SIN(PI*j*k/m))
  ! where both summations are taken for k=1,...,m-1
  !
  ! USAGE
  !
  ! Procedure CreateFFT(n) creates hidden data structures
  ! to compute sets of FFTs of length n. Should be invoked
  ! prior to the first FFT of that size. 
  !
  ! Procedures DirFFT and InvFFT compute multiple direct and inverse
  ! FFTs of size n. Successive invocations of these procedures
  ! use the same hidden data structures created by a single call
  ! to CreateFFT.
  !
  ! DestroyFFT destoys internal data structures. Usefull if
  ! one has to change the size of FFTs. In that case, 
  ! CreateFFT should be invoked with the new size, after the
  ! FFT of the previous size was destroied.
  !
  ! This configuration don't allow simultaneous computing of
  ! FFT of multiple sizes. A data type is required for that
  ! case.
  ! 
  !
  ! ALGORITHM
  !
  ! The fundamental FFT algorithm uses 2, 3, 4 and 5 as bases and
  ! operates over m complex numbers, expressed as pairs of real numbers.
  ! The algorithm requires m = 2**i * 3**j * 5**k, 
  ! with i >= 1 and at j, k >= 0
  !
  ! The direct and the inverse procedures use the fundamental FFT 
  ! algorithm over complexes builded from even and odd functions 
  ! of the input data set. Properties of even and odd functions are
  ! used to unscramble the desired results (in the direct case) 
  ! or to scramble input (in the inverse case) from/to the fundamental
  ! algorithm.
  !
  ! Periodicity is required to compute even/odd functions.
  !
  ! REFERENCES
  !
  ! Rader, C. M., "Discrete Fourier Transforms When the Number of
  ! Data Samples Is Prime", Proceedings of the IEEE, June 1968, 
  ! pp 1107-1008.
  ! 
  ! Singleton, R. C., "Algol Procedures for the Fast Fourier Transform",
  ! CACM, November 1969, pp 773-775.
  !
  ! Brigham, E. O., "The Fast Fourier Transform", Prentice-Hall
  !
  ! CODE HISTORY
  !
  ! Converted from ancient Fortran 77 1-D code by CPTEC in 2001.
  ! The conversion procedure transposes input data sets to achieve
  ! higher efficiency on memory accesses. Performance increases
  ! with the number of transforms computed simultaneously.

  PUBLIC :: CreateFFT, DestroyFFT, DirFFT, InvFFT

  INTERFACE InvFFT
     MODULE PROCEDURE InvD1, InvD2, InvD3
  END INTERFACE

  INTERFACE DirFFT
     MODULE PROCEDURE DirD1, DirD2, DirD3
  END INTERFACE

  PRIVATE

  LOGICAL :: created=.FALSE.

  INTEGER :: nGiven

  INTEGER, PARAMETER :: nBase=4
  INTEGER, PARAMETER :: Base(nBase) = (/ 4, 2, 3, 5 /)
  INTEGER, PARAMETER :: Permutation(nBase) = (/ 2, 3, 1, 4 /)

  REAL (KIND=r8) :: radi
  REAL (KIND=r8) :: sin60
  REAL (KIND=r8) :: sin36
  REAL (KIND=r8) :: sin72
  REAL (KIND=r8) :: cos36
  REAL (KIND=r8) :: cos72

  INTEGER, ALLOCATABLE :: Factors(:)

  REAL (KIND=r8), ALLOCATABLE :: Trigs(:)


CONTAINS


  ! CreateFFT: Initializes internal values for FFTs of size nIn


  SUBROUTINE CreateFFT()

    IMPLICIT NONE

    CHARACTER (LEN=15), PARAMETER :: h="**(CreateFFT)**" ! header
    IF (created) THEN
       WRITE (UNIT=nferr, FMT='(2A)') h, &
             ' invoked without destroying previous FFT'
       STOP
    END IF
    created=.TRUE.
    nGiven=Imax
    CALL Factorize (nGiven)
    CALL TrigFactors (nGiven)

  END SUBROUTINE CreateFFT


  SUBROUTINE Factorize (nIn)

    IMPLICIT NONE

    ! Factorize: Factorizes nIn/2 in powers of 4, 3, 2, 5, if possible.
    !            Otherwise, stops with error message

    INTEGER, INTENT(IN) :: nIn

    CHARACTER(LEN=15), PARAMETER :: h="**(Factorize)**" ! header
    CHARACTER(LEN=15) :: charInt ! Character representation of integer

    INTEGER :: Powers(nBase)
    INTEGER :: nOut
    INTEGER :: sumPowers
    INTEGER :: ifac
    INTEGER :: i
    INTEGER :: j

    CALL NextPossibleSize (nIn, nOut, Powers)

    IF (nIn /= nOut) THEN
       WRITE (charInt, FMT='(I15)') nIn
       WRITE (UNIT=nferr, FMT='(4A)') h, &
             ' FFT size = ', TRIM(ADJUSTL(charInt)), ' not factorizable '
       WRITE (charInt, FMT='(I15)') nOut
       WRITE (UNIT=nferr, FMT='(3A)') h, &
             ' Next factorizable FFT size is ', TRIM(ADJUSTL(charInt))
       STOP
    END IF

    sumPowers=SUM(Powers)
    ALLOCATE (Factors(sumPowers+1))
    Factors(1)=sumPowers
    ifac=1
    DO i=1,nBase
       j=Permutation(i)
       Factors(ifac+1:ifac+Powers(j))=Base(j)
       ifac=ifac+Powers(j)
    END DO

  END SUBROUTINE Factorize


  SUBROUTINE NextPossibleSize(nIn, nOut, Powers)

    IMPLICIT NONE

    ! NextPossibleSize: Next possible fft size and its factors

    INTEGER, INTENT(IN ) :: nIn
    INTEGER, INTENT(OUT) :: nOut
    INTEGER, INTENT(OUT) :: Powers(nBase)

    CHARACTER (LEN=22), PARAMETER :: h="**(NextPossibleSize)**" ! header
    CHARACTER (LEN=15) :: charNIn ! Character representation of nIn

    INTEGER :: i 
    INTEGER :: left ! portion of nOut/2 yet to be factorized

    INTEGER, PARAMETER :: limit = HUGE(nIn)-1   ! Maximum representable integer

    ! nOut = positive even 

    IF (nIn <= 0) THEN
       WRITE (charNIn, FMT='(I15)') nIn
       WRITE (UNIT=nferr, FMT='(3A)') h, &
             ' Meaningless FFT size = ', TRIM(ADJUSTL(charNIn))
       STOP
    ELSE IF (MOD(nIn,2) == 0) THEN
       nOut=nIn
    ELSE
       nOut=nIn+1
    END IF

    ! Loop over evens, starting from nOut, looking for
    ! next factorizable even/2

    DO
       left=nOut/2
       Powers=0

       ! factorize nOut/2

       DO i=1,nBase 
          DO
             IF (MOD(left, Base(i)) == 0) THEN
                Powers(i)=Powers(i)+1
                left=left/Base(i)
             ELSE
                EXIT
             END IF
          END DO
       END DO

       IF (left == 1) THEN
          EXIT
       ELSE IF (nOut < limit) THEN
          nOut=nOut+2
       ELSE
          WRITE (charNIn, FMT='(i15)') nIn
          WRITE (UNIT=nferr, FMT='(4A)') h, ' Next factorizable FFT size > ', &
                 TRIM(ADJUSTL(charNIn)), ' is not representable in this machine'
          STOP
       END IF
    END DO

  END SUBROUTINE NextPossibleSize


  SUBROUTINE TrigFactors (nIn)

    IMPLICIT NONE

    ! TrigFactors: Sin and Cos required to compute FFT of size nIn

    INTEGER, INTENT(IN) :: nIn

    INTEGER :: nn, nh, i
    REAL (KIND=r8) :: radi, pi, del, angle

    radi=ATAN(1.0_r8)/45.0_r8
    sin60=SIN(60.0_r8*radi)
    sin36=SIN(36.0_r8*radi)
    sin72=SIN(72.0_r8*radi)
    cos36=COS(36.0_r8*radi)
    cos72=COS(72.0_r8*radi)

    nn=nIn/2
    nh=(nn+1)/2
    ALLOCATE (Trigs(2*(nn+nh)))

    pi=2.0_r8*ASIN(1.0_r8)
    del=(2.0_r8*pi)/REAL(nn,r8)

    DO i=1,2*nn,2
       angle=0.5_r8*REAL(i-1,r8)*del
       Trigs(i)=COS(angle)
       Trigs(i+1)=SIN(angle)
    END DO

    del=0.5_r8*del
    DO i=1,2*nh,2
       angle=0.5_r8*REAL(i-1,r8)*del
       Trigs(2*nn+i)=COS(angle)
       Trigs(2*nn+i+1)=SIN(angle)
    END DO

  END SUBROUTINE TrigFactors


  SUBROUTINE DestroyFFT()

    IMPLICIT NONE

    ! DestroyFFT: Destrois internal data structure

    CHARACTER (LEN=16), PARAMETER :: h="**(DestroyFFT)**" ! header

    IF (.NOT. created) THEN
       WRITE (UNIT=nferr, FMT='(2A)') h, ' there is no FFT to destroy'
       STOP
    END IF

    created=.FALSE.
    DEALLOCATE (Factors)
    DEALLOCATE (Trigs)

  END SUBROUTINE DestroyFFT


  SUBROUTINE SplitFour (fin, a, b, ldin, n, nh, lot)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ldin
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: nh
    INTEGER, INTENT(IN) :: lot

    REAL (KIND=r8), INTENT(IN) :: fin(ldin,lot)
    REAL (KIND=r8), INTENT(OUT) :: a(lot,nh)
    REAL (KIND=r8), INTENT(OUT) :: b(lot,nh)

    INTEGER :: i, j

    REAL (KIND=r8) :: c, s

!cdir nodep
    DO i=1,lot
       a(i,1)=fin(1,i)+fin(n+1,i)
       b(i,1)=fin(1,i)-fin(n+1,i)
    END DO

    DO j=2,(nh+1)/2
       c=Trigs(n+2*j-1)
       s=Trigs(n+2*j)
!cdir nodep
       DO i=1,lot
          a(i,j)=(fin(2*j-1,i)+fin(n+3-2*j,i)) &
                 - (s*(fin(2*j-1,i)-fin(n+3-2*j,i)) &
                 + c*(fin(2*j,i)+fin(n+4-2*j,i)))
          a(i,nh+2-j)=(fin(2*j-1,i)+fin(n+3-2*j,i)) &
                         + (s*(fin(2*j-1,i)-fin(n+3-2*j,i)) &
                         + c*(fin(2*j,i)+fin(n+4-2*j,i)))
          b(i,j)=(c*(fin(2*j-1,i)-fin(n+3-2*j,i)) &
                 - s*(fin(2*j,i)+fin(n+4-2*j,i)))&
                 + (fin(2*j,i)-fin(n+4-2*j,i))
          b(i,nh+2-j)=(c*(fin(2*j-1,i)-fin(n+3-2*j,i)) &
                         - s*(fin(2*j,i)+fin(n+4-2*j,i)))&
                         - (fin(2*j,i)-fin(n+4-2*j,i))
       END DO
    END DO
    IF ( (nh>=2) .AND. (MOD(nh,2)==0) ) THEN
!cdir nodep
       DO i=1,lot
          a(i,nh/2+1)=2.0_r8*fin(nh+1,i)
          b(i,nh/2+1)=-2.0_r8*fin(nh+2,i)
       END DO
    END IF

  END SUBROUTINE SplitFour


  SUBROUTINE JoinFour (a, b, fout, ldout, n, nh, lot)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ldout
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: nh
    INTEGER, INTENT(IN) :: lot

    REAL (KIND=r8), INTENT(OUT) :: fout(ldout,lot)
    REAL (KIND=r8), INTENT(IN) :: a(lot,nh)
    REAL (KIND=r8), INTENT(IN) :: b(lot,nh)

    INTEGER :: i, j

    REAL (KIND=r8) :: scale, scalh
    REAL (KIND=r8) :: c, s

    scale=1.0_r8/REAL(n,r8)
    scalh=0.5_r8*scale

!cdir nodep
    DO i=1,lot
       fout(1,i)=scale*(a(i,1)+b(i,1))
       fout(n+1,i)=scale*(a(i,1)-b(i,1))
       fout(2,i)=0.0_r8
    END DO

    DO j=2,(nh+1)/2
       c=Trigs(n+2*j-1)
       s=Trigs(n+2*j)
!cdir nodep
       DO i=1,lot
          fout(2*j-1,i)=scalh*((a(i,j)+a(i,nh+2-j)) &
                              + (c*(b(i,j)+b(i,nh+2-j))&
                              + s*(a(i,j)-a(i,nh+2-j))))
          fout(n+3-2*j,i)=scalh*((a(i,j)+a(i,nh+2-j)) &
                                - (c*(b(i,j)+b(i,nh+2-j))&
                                + s*(a(i,j)-a(i,nh+2-j))))
          fout(2*j,i)=scalh*((c*(a(i,j)-a(i,nh+2-j))&
                         - s*(b(i,j)+b(i,nh+2-j)))&
                         + (b(i,nh+2-j)-b(i,j)))
          fout(n+4-2*j,i)=scalh*((c*(a(i,j)-a(i,nh+2-j))&
                                - s*(b(i,j)+b(i,nh+2-j)))&
                                - (b(i,nh+2-j)-b(i,j)))
       END DO
    END DO

    IF ((nh>=2) .AND. (MOD(nh,2)==0)) THEN
!cdir nodep
       DO i=1,lot
          fout(nh+1,i)= scale*a(i,nh/2+1)
          fout(nh+2,i)=-scale*b(i,nh/2+1)
       END DO
    END IF

    fout(n+2:ldout,:)=0.0_r8

  END SUBROUTINE JoinFour


  SUBROUTINE SplitGaus (fin, a, b, ldin, nh, lot)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ldin
    INTEGER, INTENT(IN) :: nh
    INTEGER, INTENT(IN) :: lot

    REAL (KIND=r8), INTENT(IN) :: fin(ldin,lot)
    REAL (KIND=r8), INTENT(OUT) :: a(lot,nh)
    REAL (KIND=r8), INTENT(OUT) :: b(lot,nh)

    INTEGER :: i, j

    DO j=1,nh
       DO i=1,lot
          a(i,j)=fin(2*j-1,i)
          b(i,j)=fin(2*j,i)
       END DO
    END DO

  END SUBROUTINE SplitGaus


  SUBROUTINE JoinGaus (a, b, fout, ldout, nh, lot)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ldout
    INTEGER, INTENT(IN) :: nh
    INTEGER, INTENT(IN) :: lot

    REAL (KIND=r8), INTENT(OUT) :: fout(ldout,lot)
    REAL (KIND=r8), INTENT(IN) :: a(lot,nh)
    REAL (KIND=r8), INTENT(IN) :: b(lot,nh)

    INTEGER :: i, j

    DO j=1,nh
       DO i=1,lot
          fout(2*j-1,i)=a(i,j)
          fout(2*j,i)=b(i,j)
       END DO
    END DO

    fout(2*nh+1:ldout,:)=0.0_r8

  END SUBROUTINE JoinGaus


  SUBROUTINE OnePass (a, b, c, d, lot, nh, ifac, la)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: lot
    INTEGER, INTENT(IN) :: nh         ! = PROD(factor(1:K))
    INTEGER, INTENT(IN) :: ifac       ! = factor(k)
    INTEGER, INTENT(IN) :: la         ! = PROD(factor(1:k-1))

    REAL (KIND=r8), INTENT(IN) :: a(lot,nh)
    REAL (KIND=r8), INTENT(IN) :: b(lot,nh)
    REAL (KIND=r8), INTENT(OUT) :: c(lot,nh)
    REAL (KIND=r8), INTENT(OUT) :: d(lot,nh)

    INTEGER :: m
    INTEGER :: jump
    INTEGER :: i, j, k
    INTEGER :: ia, ja
    INTEGER :: ib, jb, kb
    INTEGER :: ic, jc, kc
    INTEGER :: id, jd, kd
    INTEGER :: ie, je, ke

    REAL (KIND=r8) :: c1, s1
    REAL (KIND=r8) :: c2, s2
    REAL (KIND=r8) :: c3, s3
    REAL (KIND=r8) :: c4, s4
    REAL (KIND=r8) :: wka, wkb
    REAL (KIND=r8) :: wksina, wksinb
    REAL (KIND=r8) :: wkaacp, wkbacp
    REAL (KIND=r8) :: wkaacm, wkbacm

    m=nh/ifac
    jump=(ifac-1)*la

    ia=0
    ib=m
    ic=2*m
    id=3*m
    ie=4*m

    ja=0
    jb=la
    jc=2*la
    jd=3*la
    je=4*la

    IF (ifac == 2) THEN
       DO j=1,la
!cdir nodep
          DO i=1,lot
             c(i,j+ja)=a(i,j+ia)+a(i,j+ib)
             c(i,j+jb)=a(i,j+ia)-a(i,j+ib)
             d(i,j+ja)=b(i,j+ia)+b(i,j+ib)
             d(i,j+jb)=b(i,j+ia)-b(i,j+ib)
          END DO
       END DO
       DO k=la,m-1,la
          kb=k+k
          c1=Trigs(kb+1)
          s1=Trigs(kb+2)
          ja=ja+jump
          jb=jb+jump
          DO j=k+1,k+la
!cdir nodep
             DO i=1,lot
                wka=a(i,j+ia)-a(i,j+ib)
                c(i,j+ja)=a(i,j+ia)+a(i,j+ib)
                wkb=b(i,j+ia)-b(i,j+ib)
                d(i,j+ja)=b(i,j+ia)+b(i,j+ib)
                c(i,j+jb)=c1*wka-s1*wkb
                d(i,j+jb)=s1*wka+c1*wkb
             END DO
          END DO
       END DO
    ELSEIF (ifac == 3) THEN
       DO j=1,la
!cdir  nodep
          DO i=1,lot
             wka=a(i,j+ib)+a(i,j+ic)
             wksina=sin60*(a(i,j+ib)-a(i,j+ic))
             wkb=b(i,j+ib)+b(i,j+ic)
             wksinb=sin60*(b(i,j+ib)-b(i,j+ic))
             c(i,j+ja)=a(i,j+ia)+wka
             c(i,j+jb)=(a(i,j+ia)-0.5_r8*wka)-wksinb
             c(i,j+jc)=(a(i,j+ia)-0.5_r8*wka)+wksinb
             d(i,j+ja)=b(i,j+ia)+wkb
             d(i,j+jb)=(b(i,j+ia)-0.5_r8*wkb)+wksina
             d(i,j+jc)=(b(i,j+ia)-0.5_r8*wkb)-wksina
          END DO
       END DO
       DO k=la,m-1,la
          kb=k+k
          kc=kb+kb
          c1=Trigs(kb+1)
          s1=Trigs(kb+2)
          c2=Trigs(kc+1)
          s2=Trigs(kc+2)
          ja=ja+jump
          jb=jb+jump
          jc=jc+jump
          DO j=k+1,k+la
!cdir nodep
             DO i=1,lot
                wka=a(i,j+ib)+a(i,j+ic)
                wksina=sin60*(a(i,j+ib)-a(i,j+ic))
                wkb=b(i,j+ib)+b(i,j+ic)
                wksinb=sin60*(b(i,j+ib)-b(i,j+ic))
                c(i,j+ja)=a(i,j+ia)+wka
                d(i,j+ja)=b(i,j+ia)+wkb
                c(i,j+jb)=c1*((a(i,j+ia)-0.5_r8*wka)-wksinb) &
                          - s1*((b(i,j+ia)-0.5_r8*wkb)+wksina)
                d(i,j+jb)=s1*((a(i,j+ia)-0.5_r8*wka)-wksinb) &
                          + c1*((b(i,j+ia)-0.5_r8*wkb)+wksina)
                c(i,j+jc)=c2*((a(i,j+ia)-0.5_r8*wka)+wksinb) &
                          - s2*((b(i,j+ia)-0.5_r8*wkb)-wksina)
                d(i,j+jc)=s2*((a(i,j+ia)-0.5_r8*wka)+wksinb) &
                          + c2*((b(i,j+ia)-0.5_r8*wkb)-wksina)
             END DO
          END DO
       END DO
    ELSEIF (ifac == 4) THEN
       DO j=1,la
!cdir nodep
          DO i=1,lot
             wkaacp=a(i,j+ia)+a(i,j+ic)
             wkaacm=a(i,j+ia)-a(i,j+ic)
             wkbacp=b(i,j+ia)+b(i,j+ic)
             wkbacm=b(i,j+ia)-b(i,j+ic)
             c(i,j+ja)=wkaacp+(a(i,j+ib)+a(i,j+id))
             c(i,j+jc)=wkaacp-(a(i,j+ib)+a(i,j+id))
             d(i,j+jb)=wkbacm+(a(i,j+ib)-a(i,j+id))
             d(i,j+jd)=wkbacm-(a(i,j+ib)-a(i,j+id))
             d(i,j+ja)=wkbacp+(b(i,j+ib)+b(i,j+id))
             d(i,j+jc)=wkbacp-(b(i,j+ib)+b(i,j+id))
             c(i,j+jb)=wkaacm-(b(i,j+ib)-b(i,j+id))
             c(i,j+jd)=wkaacm+(b(i,j+ib)-b(i,j+id))
          END DO
       END DO
       DO k=la,m-1,la
          kb=k+k
          kc=kb+kb
          kd=kc+kb
          c1=Trigs(kb+1)
          s1=Trigs(kb+2)
          c2=Trigs(kc+1)
          s2=Trigs(kc+2)
          c3=Trigs(kd+1)
          s3=Trigs(kd+2)
          ja=ja+jump
          jb=jb+jump
          jc=jc+jump
          jd=jd+jump
          DO j=k+1,k+la
!cdir nodep
             DO i=1,lot
                wkaacp=a(i,j+ia)+a(i,j+ic)
                wkbacp=b(i,j+ia)+b(i,j+ic)
                wkaacm=a(i,j+ia)-a(i,j+ic)
                wkbacm=b(i,j+ia)-b(i,j+ic)
                c(i,j+ja)=wkaacp+(a(i,j+ib)+a(i,j+id))
                d(i,j+ja)=wkbacp+(b(i,j+ib)+b(i,j+id))
                c(i,j+jc)=c2*(wkaacp-(a(i,j+ib)+a(i,j+id))) &
                          - s2*(wkbacp-(b(i,j+ib)+b(i,j+id))) 
                d(i,j+jc)=s2*(wkaacp-(a(i,j+ib)+a(i,j+id))) &
                          + c2*(wkbacp-(b(i,j+ib)+b(i,j+id)))
                c(i,j+jb)=c1*(wkaacm-(b(i,j+ib)-b(i,j+id))) &
                          - s1*(wkbacm+(a(i,j+ib)-a(i,j+id)))
                d(i,j+jb)=s1*(wkaacm-(b(i,j+ib)-b(i,j+id))) &
                          + c1*(wkbacm+(a(i,j+ib)-a(i,j+id)))
                c(i,j+jd)=c3*(wkaacm+(b(i,j+ib)-b(i,j+id))) &
                          - s3*(wkbacm-(a(i,j+ib)-a(i,j+id)))
                d(i,j+jd)=s3*(wkaacm+(b(i,j+ib)-b(i,j+id))) &
                          + c3*(wkbacm-(a(i,j+ib)-a(i,j+id)))
             END DO
          END DO
       END DO
    ELSEIF (ifac == 5) THEN
       DO j=1,la
!cdir nodep
        DO i=1,lot
             c(i,j+ja)=a(i,j+ia)+(a(i,j+ib)+a(i,j+ie))+(a(i,j+ic)+a(i,j+id))
             d(i,j+ja)=b(i,j+ia)&
                       + (b(i,j+ib)+b(i,j+ie))&
                       + (b(i,j+ic)+b(i,j+id))
             c(i,j+jb)=(a(i,j+ia)&
                       + cos72*(a(i,j+ib)+a(i,j+ie))&
                       - cos36*(a(i,j+ic)+a(i,j+id)))&
                       -(sin72*(b(i,j+ib)-b(i,j+ie))&
                       + sin36*(b(i,j+ic)-b(i,j+id)))
             c(i,j+je)=(a(i,j+ia)&
                       + cos72*(a(i,j+ib)+a(i,j+ie))&
                       - cos36*(a(i,j+ic)+a(i,j+id)))&
                       +(sin72*(b(i,j+ib)-b(i,j+ie))&
                       + sin36*(b(i,j+ic)-b(i,j+id)))
             d(i,j+jb)=(b(i,j+ia)&
                       + cos72*(b(i,j+ib)+b(i,j+ie))&
                       - cos36*(b(i,j+ic)+b(i,j+id)))&
                       +(sin72*(a(i,j+ib)-a(i,j+ie))&
                       + sin36*(a(i,j+ic)-a(i,j+id)))
             d(i,j+je)=(b(i,j+ia)&
                       + cos72*(b(i,j+ib)+b(i,j+ie))&
                       - cos36*(b(i,j+ic)+b(i,j+id)))&
                       -(sin72*(a(i,j+ib)-a(i,j+ie))&
                       + sin36*(a(i,j+ic)-a(i,j+id)))
             c(i,j+jc)=(a(i,j+ia)&
                       - cos36*(a(i,j+ib)+a(i,j+ie))&
                       + cos72*(a(i,j+ic)+a(i,j+id)))&
                       -(sin36*(b(i,j+ib)-b(i,j+ie))&
                       - sin72*(b(i,j+ic)-b(i,j+id)))
             c(i,j+jd)=(a(i,j+ia)&
                       - cos36*(a(i,j+ib)+a(i,j+ie))&
                       + cos72*(a(i,j+ic)+a(i,j+id)))&
                       +(sin36*(b(i,j+ib)-b(i,j+ie))&
                       - sin72*(b(i,j+ic)-b(i,j+id)))
             d(i,j+jc)=(b(i,j+ia)&
                       - cos36*(b(i,j+ib)+b(i,j+ie))&
                       + cos72*(b(i,j+ic)+b(i,j+id)))&
                       +(sin36*(a(i,j+ib)-a(i,j+ie))&
                       - sin72*(a(i,j+ic)-a(i,j+id)))
             d(i,j+jd)=(b(i,j+ia)&
                       - cos36*(b(i,j+ib)+b(i,j+ie))&
                       + cos72*(b(i,j+ic)+b(i,j+id)))&
                       -(sin36*(a(i,j+ib)-a(i,j+ie))&
                       - sin72*(a(i,j+ic)-a(i,j+id)))
          END DO
       END DO
       DO k=la,m-1,la
          kb=k+k
          kc=kb+kb
          kd=kc+kb
          ke=kd+kb
          c1=Trigs(kb+1)
          s1=Trigs(kb+2)
          c2=Trigs(kc+1)
          s2=Trigs(kc+2)
          c3=Trigs(kd+1)
          s3=Trigs(kd+2)
          c4=Trigs(ke+1)
          s4=Trigs(ke+2)
          ja=ja+jump
          jb=jb+jump
          jc=jc+jump
          jd=jd+jump
          je=je+jump
          DO j=k+1,k+la
!cdir nodep
             DO i=1,lot
                c(i,j+ja)=a(i,j+ia)&
                          + (a(i,j+ib)+a(i,j+ie))&
                          + (a(i,j+ic)+a(i,j+id))
                d(i,j+ja)=b(i,j+ia)&
                          + (b(i,j+ib)+b(i,j+ie))&
                          + (b(i,j+ic)+b(i,j+id))
                c(i,j+jb)=c1*((a(i,j+ia)&
                          + cos72*(a(i,j+ib)+a(i,j+ie))&
                          - cos36*(a(i,j+ic)+a(i,j+id)))&
                          - (sin72*(b(i,j+ib)-b(i,j+ie))&
                          + sin36*(b(i,j+ic)-b(i,j+id))))&
                          - s1*((b(i,j+ia)&
                          + cos72*(b(i,j+ib)+b(i,j+ie))&
                          - cos36*(b(i,j+ic)+b(i,j+id)))&
                          + (sin72*(a(i,j+ib)-a(i,j+ie))&
                          + sin36*(a(i,j+ic)-a(i,j+id))))
                d(i,j+jb)=s1*((a(i,j+ia)&
                          + cos72*(a(i,j+ib)+a(i,j+ie))&
                          - cos36*(a(i,j+ic)+a(i,j+id)))&
                          - (sin72*(b(i,j+ib)-b(i,j+ie))&
                          + sin36*(b(i,j+ic)-b(i,j+id))))&
                          + c1*((b(i,j+ia)&
                          + cos72*(b(i,j+ib)+b(i,j+ie))&
                          - cos36*(b(i,j+ic)+b(i,j+id)))&
                          + (sin72*(a(i,j+ib)-a(i,j+ie)) &
                          + sin36*(a(i,j+ic)-a(i,j+id))))
                c(i,j+je)=c4*((a(i,j+ia)&
                          + cos72*(a(i,j+ib)+a(i,j+ie))&
                          - cos36*(a(i,j+ic)+a(i,j+id)))&
                          + (sin72*(b(i,j+ib)-b(i,j+ie)) &
                          + sin36*(b(i,j+ic)-b(i,j+id)))) &
                          - s4*((b(i,j+ia)&
                          + cos72*(b(i,j+ib)+b(i,j+ie))&
                          - cos36*(b(i,j+ic)+b(i,j+id)))&
                          - (sin72*(a(i,j+ib)-a(i,j+ie))&
                          +  sin36*(a(i,j+ic)-a(i,j+id))))
                d(i,j+je)=s4*((a(i,j+ia)&
                          + cos72*(a(i,j+ib)+a(i,j+ie))&
                          - cos36*(a(i,j+ic)+a(i,j+id)))&
                          + (sin72*(b(i,j+ib)-b(i,j+ie))&
                          + sin36*(b(i,j+ic)-b(i,j+id))))&
                          + c4*((b(i,j+ia)&
                          + cos72*(b(i,j+ib)+b(i,j+ie))&
                          - cos36*(b(i,j+ic)+b(i,j+id))) &
                          - (sin72*(a(i,j+ib)-a(i,j+ie))&
                          + sin36*(a(i,j+ic)-a(i,j+id))))
                c(i,j+jc)=c2*((a(i,j+ia)&
                          - cos36*(a(i,j+ib)+a(i,j+ie))&
                          + cos72*(a(i,j+ic)+a(i,j+id)))&
                          - (sin36*(b(i,j+ib)-b(i,j+ie))&
                          - sin72*(b(i,j+ic)-b(i,j+id))))&
                          - s2*((b(i,j+ia)&
                          - cos36*(b(i,j+ib)+b(i,j+ie))&
                          + cos72*(b(i,j+ic)+b(i,j+id)))&
                          + (sin36*(a(i,j+ib)-a(i,j+ie))&
                          - sin72*(a(i,j+ic)-a(i,j+id)))) 
                d(i,j+jc)=s2*((a(i,j+ia)&
                          - cos36*(a(i,j+ib)+a(i,j+ie))&
                          + cos72*(a(i,j+ic)+a(i,j+id)))&
                          - (sin36*(b(i,j+ib)-b(i,j+ie))&
                          - sin72*(b(i,j+ic)-b(i,j+id))))&
                          + c2*((b(i,j+ia)&
                          - cos36*(b(i,j+ib)+b(i,j+ie))&
                          + cos72*(b(i,j+ic)+b(i,j+id)))&
                          + (sin36*(a(i,j+ib)-a(i,j+ie))&
                          - sin72*(a(i,j+ic)-a(i,j+id))))
                c(i,j+jd)=c3*((a(i,j+ia)&
                          - cos36*(a(i,j+ib)+a(i,j+ie))&
                          + cos72*(a(i,j+ic)+a(i,j+id)))&
                          + (sin36*(b(i,j+ib)-b(i,j+ie))&
                          - sin72*(b(i,j+ic)-b(i,j+id))))&
                          - s3*((b(i,j+ia)&
                          - cos36*(b(i,j+ib)+b(i,j+ie))&
                          + cos72*(b(i,j+ic)+b(i,j+id)))&
                          - (sin36*(a(i,j+ib)-a(i,j+ie))&
                          - sin72*(a(i,j+ic)-a(i,j+id))))
                d(i,j+jd)=s3*((a(i,j+ia)&
                          - cos36*(a(i,j+ib)+a(i,j+ie))&
                          + cos72*(a(i,j+ic)+a(i,j+id)))&
                          + (sin36*(b(i,j+ib)-b(i,j+ie))&
                          - sin72*(b(i,j+ic)-b(i,j+id))))&
                          + c3*((b(i,j+ia)&
                          - cos36*(b(i,j+ib)+b(i,j+ie))&
                          + cos72*(b(i,j+ic)+b(i,j+id)))&
                          - (sin36*(a(i,j+ib)-a(i,j+ie))&
                          - sin72*(a(i,j+ic)-a(i,j+id))))
             END DO
          END DO
       END DO
    ENDIF

  END SUBROUTINE OnePass


  SUBROUTINE Dir (fin, fout, ldin, ldout, n, lot)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ldin
    INTEGER, INTENT(IN) :: ldout
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: lot

    REAL (KIND=r8), INTENT(IN) :: fin (ldin,lot)
    REAL (KIND=r8), INTENT(OUT) :: fout(ldout,lot)

    INTEGER :: nh
    INTEGER :: nfax
    INTEGER :: la
    INTEGER :: k

    LOGICAL :: ab2cd

    CHARACTER (LEN=12), PARAMETER :: h="**(Dir)**"
    CHARACTER (LEN=15) :: charInt1, charInt2

    REAL (KIND=r8), DIMENSION (lot,n/2) :: a, b, c, d

    IF (.NOT. created) THEN
       WRITE (UNIT=nferr, FMT='(2A)') h, ' FFT was not created'
       STOP
    ELSE IF (n /= nGiven) THEN
       WRITE (charInt1, FMT='(I15)') n
       WRITE (charInt2, FMT='(I15)') nGiven
       WRITE (UNIT=nferr, FMT='(4A)') h, &
             ' FFT invoked with size ', TRIM(ADJUSTL(charInt1)), &
             ' but created with size ', TRIM(ADJUSTL(charInt2))
       STOP
    ELSE IF (ldout < n+1) THEN
       WRITE (charInt1, FMT='(I15)') ldout
       WRITE (charInt2, FMT='(I15)') n+1
       WRITE (UNIT=nferr, FMT='(4A)') h, &
             ' Output field has first dimension ', TRIM(ADJUSTL(charInt1)), &
             '; should be at least ', TRIM(ADJUSTL(charInt2))
       STOP
    ELSE IF (ldin < n) THEN
       WRITE (charInt1, FMT='(I15)') ldin
       WRITE (charInt2, FMT='(I15)') n
       WRITE (UNIT=nferr, FMT='(4A)') h, &
             ' Input field has first dimension ', TRIM(ADJUSTL(charInt1)), &
             '; should be at least ', TRIM(ADJUSTL(charInt2))
       STOP
    END IF

    nfax=Factors(1)
    nh=n/2

    CALL SplitGaus (fin, a, b, ldin, nh, lot)

    la=1
    ab2cd=.TRUE.
    DO k=1,nfax
       IF (ab2cd) THEN
          CALL OnePass (a, b, c, d, lot, nh, Factors(k+1), la)
          ab2cd=.FALSE.
       ELSE
          CALL OnePass (c, d, a, b, lot, nh, Factors(k+1), la)
          ab2cd=.TRUE.
       END IF
       la=la*Factors(k+1)
    END DO

    IF (ab2cd) THEN
       CALL JoinFour (a, b, fout, ldout, n, nh, lot)
    ELSE
       CALL JoinFour (c, d, fout, ldout, n, nh, lot)
    END IF

  END SUBROUTINE Dir


  SUBROUTINE DirD1 (fin, fout)

    IMPLICIT NONE

    REAL (KIND=r8), INTENT(IN) :: fin(:)
    REAL (KIND=r8), INTENT(OUT) :: fout(:)

    INTEGER :: din, dout

    din=SIZE(fin)
    dout=SIZE(fout)

    CALL Dir (fin, fout, din, dout, nGiven, 1)

  END SUBROUTINE DirD1


  SUBROUTINE DirD2 (fin, fout)

    REAL (KIND=r8), INTENT(IN) :: fin(:,:)
    REAL (KIND=r8), INTENT(OUT) :: fout(:,:)

    INTEGER :: d2in, d2out, din, dout

    CHARACTER (LEN=11), PARAMETER :: h="**(DirD2)**"

    d2in=SIZE(fin,2)
    d2out=SIZE(fout,2)

    IF (d2in /= d2out) THEN
       WRITE (UNIT=nferr, FMT='(2A,2I10)') h, &
             ' Error: dim 2 of fin, fout, differ: ', d2in, d2out
       STOP
    END IF

    din=SIZE(fin,1)
    dout=SIZE(fout,1)

    CALL Dir (fin, fout, din, dout, nGiven, d2in)

  END SUBROUTINE DirD2


  SUBROUTINE DirD3 (fin, fout)

    IMPLICIT NONE

    REAL (KIND=r8), INTENT(IN) :: fin(:,:,:)
    REAL (KIND=r8), INTENT(OUT) :: fout(:,:,:)

    INTEGER :: d2in, d2out, d3in, d3out, din, dout, dio

    CHARACTER (LEN=11), PARAMETER :: h="**(DirD3)**"

    d2in=SIZE(fin,2)
    d2out=SIZE(fout,2)
    IF (d2in /= d2out) THEN
       WRITE (UNIT=nferr, FMT='(2A,2I10)') h, &
             ' Error: dim 2 of fin, fout, differ: ', d2in, d2out
       STOP
    END IF

    d3in=SIZE(fin,3)
    d3out=SIZE(fout,3)
    IF (d3in /= d3out) THEN
       WRITE (UNIT=nferr, FMT='(2A,2I10)') h, &
             ' Error: dim 3 of fin, fout, differ: ', d3in, d3out
       STOP
    END IF

    din=SIZE(fin,1)
    dout=SIZE(fout,1)
    dio=d2in*d3in

    CALL Dir (fin, fout, din, dout, nGiven, dio)

  END SUBROUTINE DirD3


  SUBROUTINE Inv (fin, fout, ldin, ldout, n, lot)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ldin
    INTEGER, INTENT(IN) :: ldout
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: lot

    REAL (KIND=r8), INTENT(IN) :: fin (ldin,lot)
    REAL (KIND=r8), INTENT(OUT) :: fout(ldout,lot)


    INTEGER :: nh
    INTEGER :: nfax
    INTEGER :: la
    INTEGER :: k

    LOGICAL :: ab2cd

    CHARACTER (LEN=12), PARAMETER :: h="**(Inv)**"
    CHARACTER (LEN=15) :: charInt1, charInt2

    REAL (KIND=r8), DIMENSION(lot,n/2) :: a, b, c, d

    IF (.NOT. created) THEN
       WRITE (UNIT=nferr, FMT='(4A)') h, ' FFT was not created'
       STOP
    ELSE IF (n /= nGiven) THEN
       WRITE (charInt1, FMT='(I15)') n
       WRITE (charInt2, FMT='(I15)') nGiven
       WRITE (UNIT=nferr, FMT='(4A)') h, &
             ' FFT invoked with size ', TRIM(ADJUSTL(charInt1)), &
             ' but created with size ', TRIM(ADJUSTL(charInt2))
       STOP
    ELSE IF (ldin < n+1) THEN
       WRITE (charInt1, FMT='(I15)') ldin
       WRITE (charInt2, FMT='(I15)') n+1
       WRITE (UNIT=nferr, FMT='(4A)') h, &
             ' Input field has first dimension ', TRIM(ADJUSTL(charInt1)), &
             '; should be at least ', TRIM(ADJUSTL(charInt2))
       STOP
    ELSE IF (ldout < n) THEN
       WRITE (charInt1, FMT='(I15)') ldout
       WRITE (charInt2, FMT='(I15)') n
       WRITE (UNIT=nferr, FMT='(4A)') h, &
             ' Output field has first dimension ', TRIM(ADJUSTL(charInt1)), &
             '; should be at least ', TRIM(ADJUSTL(charInt2))
       STOP
    END IF

    nfax=Factors(1)
    nh=n/2

    CALL SplitFour (fin, a, b, ldin, n, nh, lot)

    la=1
    ab2cd=.TRUE.
    DO k=1,nfax
       IF (ab2cd) THEN
          CALL OnePass (a, b, c, d, lot, nh, Factors(k+1), la)
          ab2cd=.FALSE.
       ELSE
          CALL OnePass (c, d, a, b, lot, nh, Factors(k+1), la)
          ab2cd=.TRUE.
       END IF
       la=la*Factors(k+1)
    END DO

    IF (ab2cd) THEN
       CALL JoinGaus (a, b, fout, ldout, nh, lot)
    ELSE
       CALL JoinGaus (c, d, fout, ldout, nh, lot)
    END IF

  END SUBROUTINE Inv


  SUBROUTINE InvD1 (fin, fout)

    IMPLICIT NONE

    REAL (KIND=r8), INTENT(IN) :: fin(:)
    REAL (KIND=r8), INTENT(OUT) :: fout(:)

    INTEGER :: din, dout

    din=SIZE(fin)
    dout=SIZE(fout)

    CALL Inv (fin, fout, din, dout, nGiven, 1)

  END SUBROUTINE InvD1


  SUBROUTINE InvD2 (fin, fout)

    IMPLICIT NONE

    REAL (KIND=r8), INTENT(IN ) :: fin(:,:)
    REAL (KIND=r8), INTENT(OUT) :: fout(:,:)

    INTEGER :: d2in, d2out, din, dout

    CHARACTER (LEN=11), PARAMETER :: h="**(InvD2)**"

    d2in=SIZE(fin,2)
    d2out=SIZE(fout,2)

    IF (d2in /= d2out) THEN
       WRITE (UNIT=nferr, FMT='(2A,2I10)') h, &
             ' Error: dim 2 of fin, fout, differ: ', d2in, d2out
       STOP
    END IF

    din=SIZE(fin,1)
    dout=SIZE(fout,1)

    CALL Inv (fin, fout, din, dout, nGiven, d2in)

  END SUBROUTINE InvD2


  SUBROUTINE InvD3 (fin, fout)

    IMPLICIT NONE

    REAL (KIND=r8), INTENT(IN ) :: fin(:,:,:)
    REAL (KIND=r8), INTENT(OUT) :: fout(:,:,:)

    INTEGER :: d2in, d2out, d3in, d3out, din, dout, dio

    CHARACTER (LEN=11), PARAMETER :: h="**(InvD3)**"

    d2in=SIZE(fin,2)
    d2out=SIZE(fout,2)
    IF (d2in /= d2out) THEN
       WRITE (UNIT=nferr, FMT='(2A,2I10)') h, &
             ' Error: dim 2 of fin, fout, differ: ', d2in, d2out
       STOP
    END IF

    d3in=SIZE(fin,3)
    d3out=SIZE(fout,3)
    IF (d3in /= d3out) THEN
       WRITE (UNIT=nferr, FMT='(2A,2I10)') h, &
             ' Error: dim 3 of fin, fout, differ: ', d3in, d3out
       STOP
    END IF

    din=SIZE(fin,1)
    dout=SIZE(fout,1)
    dio=d2in*d3in

    CALL Inv (fin, fout, din, dout, nGiven, dio)

  END SUBROUTINE InvD3


END MODULE FastFourierTransform
