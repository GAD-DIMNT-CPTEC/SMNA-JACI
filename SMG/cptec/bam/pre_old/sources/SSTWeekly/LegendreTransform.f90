!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
MODULE LegendreTransform

  USE InputParameters, ONLY : r8, EMRad1, EMRad12, nferr, &
                              Jmax, Jmaxhf, Mend1, Mend2, &
                              Mnwv0, Mnwv1, Mnwv2, Mnwv3

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: CreateLegTrans, CreateSpectralRep, CreateGaussRep, &
            transs, Spec2Four, Four2Spec, SplitTrans

  INTEGER, ALLOCATABLE, DIMENSION (:,:) :: la0, la1

  REAL (KIND=r8), ALLOCATABLE, DIMENSION (:) :: eps

  REAL (KIND=r8), ALLOCATABLE, DIMENSION (:) :: ColRad, rCs2, Wgt

  REAL (KIND=r8), ALLOCATABLE, DIMENSION (:,:) :: LegS2F, LegExtS2F, LegDerS2F

  REAL (KIND=r8), ALLOCATABLE, DIMENSION (:,:) :: LegF2S, LegDerNS, LegDerEW

  INTEGER, ALLOCATABLE, DIMENSION (:) :: LenDiag, LenDiagExt, &
                                         LastPrevDiag, LastPrevDiagExt

  LOGICAL :: created=.FALSE.

  INTERFACE SplitTrans
     MODULE PROCEDURE SplitTrans2D, SplitTrans3D
  END INTERFACE

  INTERFACE Spec2Four
     MODULE PROCEDURE Spec2Four1D, Spec2Four2D
  END INTERFACE

  INTERFACE Four2Spec
     MODULE PROCEDURE Four2Spec1D, Four2Spec2D
  END INTERFACE


CONTAINS


  SUBROUTINE CreateSpectralRep ()

    IMPLICIT NONE

    INTEGER :: l, mm, nn

    REAL (KIND=r8) :: am, an

    ALLOCATE (la0(Mend1,Mend1))
    ALLOCATE (la1(Mend1,Mend2))
    ALLOCATE (eps(Mnwv1))

    l=0
    DO nn=1,Mend1
       DO mm=1,Mend2-nn
          l=l+1
          la0(mm,nn)=l
       END DO
    END DO
    l=0
    DO mm=1,Mend1
       l=l+1
       la1(mm,1)=l
    END DO
    DO nn=2,Mend2
       DO mm=1,Mend1+2-nn
          l=l+1
          la1(mm,nn)=l
       END DO
    END DO

    DO l=1,Mend1
       eps(l)=0.0_r8
    END DO
    l=Mend1
    DO nn=2,Mend2
       DO mm=1,Mend1+2-nn
          l=l+1
          am=mm-1
          an=mm+nn-2
          eps(l)=SQRT((an*an-am*am)/(4.0_r8*an*an-1.0_r8))
       END DO
    END DO

  END SUBROUTINE CreateSpectralRep


  SUBROUTINE CreateGaussRep ()

    IMPLICIT NONE

    INTEGER :: j
    REAL (KIND=r8) :: rd

    ALLOCATE (ColRad(Jmaxhf))
    ALLOCATE (rCs2(Jmaxhf))
    ALLOCATE (Wgt(Jmaxhf))

    CALL gLats ()

  END SUBROUTINE CreateGaussRep


  SUBROUTINE gLats ()

    IMPLICIT NONE

    ! Calculates Gaussian Latitudes and Gaussian Weights 
    ! for Use in Grid-Spectral and Spectral-Grid Transforms

    INTEGER :: j
    REAL (KIND=r8) :: Epsil
    REAL (KIND=r8) :: scale, dgColIn, gCol, dgCol, p2, p1

    Epsil=EPSILON(1.0_r8)*100.0_r8
    scale=2.0_r8/(REAL(Jmax,r8)*REAL(Jmax,r8))
    dgColIn=ATAN(1.0_r8)/REAL(Jmax,r8)
    gCol=0.0_r8
    DO j=1,Jmaxhf
       dgCol=dgColIn
       DO
          CALL LegendrePolynomial (Jmax, gCol, p2)
          DO
             p1=p2
             gCol=gCol+dgCol
             CALL LegendrePolynomial (Jmax, gCol, p2)
             IF (SIGN(1.0_r8,p1) /= SIGN(1.0_r8,p2)) EXIT
          END DO
          IF (dgCol <= Epsil) EXIT
          gCol=gCol-dgCol
          dgCol=dgCol*0.25_r8
       END DO
       ColRad(j)=gCol
       CALL LegendrePolynomial (Jmax-1, gCol, p1)
       Wgt(j)=scale*(1.0_r8-COS(gCol)*COS(gCol))/(p1*p1)
       rCs2(j)=1.0_r8/(SIN(gCol)*SIN(gCol))
    END DO

  END SUBROUTINE gLats


  SUBROUTINE LegendrePolynomial (N, Colatitude, Pln)

    IMPLICIT NONE

    ! Calculates the Value of the Ordinary Legendre 
    ! Function of Given Order at a Specified Colatitude.  
    ! Used to Determine Gaussian Latitudes.

    INTEGER, INTENT(IN ) :: N
    REAL (KIND=r8), INTENT(IN ) :: Colatitude
    REAL (KIND=r8), INTENT(OUT) :: Pln

    INTEGER :: i
    REAL (KIND=r8) :: x, y1, y2, y3, g

    x=COS(Colatitude)
    y1=1.0_r8
    y2=x
    DO i=2,N
       g=x*y2
       y3=g-y1+g-(g-y1)/REAL(i,r8)
       y1=y2
       y2=y3
    END DO
    Pln=y3

  END SUBROUTINE LegendrePolynomial


  SUBROUTINE CreateLegTrans ()

    IMPLICIT NONE

    CHARACTER (LEN=20) :: h="**(CreateLegTrans)**"

    INTEGER :: diag

    IF (created) THEN
       WRITE (UNIT=nferr, FMT='(2A)') h, ' already created'
       STOP
    ELSE
       created=.TRUE.
    END IF

    ! Associated Legendre Functions

    ALLOCATE (LegS2F   (Mnwv2, jMaxHf))
    ALLOCATE (LegDerS2F(Mnwv2, jMaxHf))
    ALLOCATE (LegExtS2F(Mnwv3, jMaxHf))
    ALLOCATE (LegF2S   (Mnwv2, jMaxHf))
    ALLOCATE (LegDerNS (Mnwv2, jMaxHf))
    ALLOCATE (LegDerEW (Mnwv2, jMaxHf))

    CALL LegPols()

    ! diagonal length

    ALLOCATE (LenDiag(Mend1))
    ALLOCATE (LenDiagExt(Mend2))
    DO diag=1,Mend1
       LenDiag(diag)=2*(Mend1+1-diag)
    END DO
    LenDiagExt(1) = 2*Mend1
    DO diag=2,Mend2
       LenDiagExt(diag)=2*(Mend1+2-diag)
    END DO

    ! last element previous diagonal

    ALLOCATE (LastPrevDiag(Mend1))
    ALLOCATE (LastPrevDiagExt(Mend2))
    DO diag=1,Mend1
       LastPrevDiag(diag)=(diag-1)*(2*Mend1+2-diag)
    END DO
    LastPrevDiagExt(1)=0
    DO diag=2,Mend2
       LastPrevDiagExt(diag)=(diag-1)*(2*Mend1+4-diag)-2
    END DO

  END SUBROUTINE CreateLegTrans


  SUBROUTINE LegPols()

   IMPLICIT NONE

    INTEGER :: j, l, nn, mm, mn, lx

    REAL (KIND=r8) :: rd
    REAL (KIND=r8) :: pln(Mnwv1)
    REAL (KIND=r8) :: dpln(Mnwv0)
    REAL (KIND=r8) :: der(Mnwv0)
    REAL (KIND=r8) :: plnwcs(Mnwv0)

    rd=45.0_r8/ATAN(1.0_r8)
    DO j=1,jMaxHf
       CALL pln2 (pln, colrad, j, eps, la1)
       l=0
       DO nn=1,Mend1
          DO mm=1,Mend2-nn
             l=l+1
             lx=la1(mm,nn)
             LegS2F(2*l-1,j)=pln(lx)
             LegS2F(2*l,j)=pln(lx)
          END DO
       END DO
       DO mn=1,Mnwv2
          LegF2S(mn,j)=LegS2F(mn,j)*wgt(j)
       END DO
       DO mn=1,Mnwv1
          LegExtS2F(2*mn-1,j)=pln(mn)
          LegExtS2F(2*mn,j)=pln(mn)
       END DO
       CALL plnder (pln, dpln, der, plnwcs, rcs2(j), wgt(j), eps, la1)
       DO mn=1,Mnwv0
          LegDerS2F(2*mn-1,j)=dpln(mn)
          LegDerS2F(2*mn,j)=dpln(mn)
          LegDerNS(2*mn-1,j)=der(mn)
          LegDerNS(2*mn,j)=der(mn)
          LegDerEW(2*mn-1,j)=plnwcs(mn)
          LegDerEW(2*mn,j)=plnwcs(mn)
       END DO
    END DO

  END SUBROUTINE LegPols


  SUBROUTINE pln2 (sln, colrad, lat, eps, la1)

    IMPLICIT NONE

    !      pln2: calculates the associated legendre functions
    !            at one specified latitude.

    REAL (KIND=r8), INTENT(OUT) :: sln(Mnwv1)
    REAL (KIND=r8), INTENT(IN) :: colrad(jMaxHf)
    REAL (KIND=r8), INTENT(IN) :: eps(Mnwv1)

    INTEGER, INTENT(IN) :: lat
    INTEGER, INTENT(IN) :: la1(Mend1,Mend2)

    INTEGER :: mm, nn, lx, ly, lz

    REAL (KIND=r8) :: colr, sinlat, coslat, prod

    LOGICAL, SAVE :: first = .TRUE.

    REAL (KIND=r8), ALLOCATABLE, DIMENSION (:), SAVE :: x, y
    REAL (KIND=r8), SAVE :: rthf

    IF (first) THEN
       ALLOCATE (x(Mend1))
       ALLOCATE (y(Mend1))
       first=.FALSE.
       DO mm=1,Mend1
          x(mm)=SQRT(2.0_r8*mm+1.0_r8)
          y(mm)=SQRT(1.0_r8+0.5_r8/REAL(mm,r8))
       END DO
       rthf=SQRT(0.5_r8)
    ENDIF
    colr=colrad(lat)
    sinlat=COS(colr)
    coslat=SIN(colr)
    prod=1.0_r8
    DO mm=1,Mend1
       sln(mm)=rthf*prod
       !     line below should only be used where exponent range is limted
       !     if(prod < flim) prod=0.0_r8
       prod=prod*coslat*y(mm)
    END DO

    DO mm=1,Mend1
       sln(mm+Mend1)=x(mm)*sinlat*sln(mm)
    END DO
    DO nn=3,Mend2
       DO mm=1,Mend1+2-nn
          lx=la1(mm,nn)
          ly=la1(mm,nn-1)
          lz=la1(mm,nn-2)
          sln(lx)=(sinlat*sln(ly)-eps(ly)*sln(lz))/eps(lx)
       END DO
    END DO

  END SUBROUTINE pln2


  SUBROUTINE plnder (pln, dpln, der, plnwcs, rcs2l, wgtl, eps, la1)

    IMPLICIT NONE

    !     calculates zonal and meridional pseudo-derivatives as
    !     well as laplacians of the associated legendre functions.

    REAL (KIND=r8), INTENT(INOUT) :: pln(Mnwv1)
    REAL (KIND=r8), INTENT(OUT) :: dpln(Mnwv0)
    REAL (KIND=r8), INTENT(OUT) :: der(Mnwv0)
    REAL (KIND=r8), INTENT(OUT) :: plnwcs(Mnwv0)
    REAL (KIND=r8), INTENT(IN) :: rcs2l
    REAL (KIND=r8), INTENT(IN) :: wgtl
    REAL (KIND=r8), INTENT(IN) :: eps(Mnwv1)

    INTEGER, INTENT(IN) :: la1(Mend1,Mend2)

    INTEGER :: n, l, nn, mm, mn, lm, l0, lp

    REAL (KIND=r8) :: raa, wcsa
    REAL (KIND=r8) :: x(Mnwv1)

    LOGICAL, SAVE :: first=.TRUE.

    REAL (KIND=r8), ALLOCATABLE, SAVE :: an(:)
    ! 
    !     compute pln derivatives
    ! 
    IF (first) THEN
       ALLOCATE (an(Mend2))
       DO n=1,Mend2
          an(n)=REAL(n-1,r8)
       END DO
       first=.FALSE.
    END IF
    raa=wgtl*EMRad12
    wcsa=rcs2l*wgtl*EMRad1
    l=0
    DO mm=1,Mend1
       l=l+1
       x(l)=an(mm)
    END DO
    DO nn=2,Mend2
       DO mm=1,Mend1+2-nn
          l=l+1
          x(l)=an(mm+nn-1)
       END DO
    END DO
    l=Mend1
    DO nn=2,Mend1
       DO mm=1,Mend2-nn
          l=l+1
          lm=la1(mm,nn-1)
          l0=la1(mm,nn)
          lp=la1(mm,nn+1)
          der(l)=x(lp)*eps(l0)*pln(lm)-x(l0)*eps(lp)*pln(lp)
       END DO
    END DO
    DO mm=1,Mend1
       der(mm)=-x(mm)*eps(mm+Mend1)*pln(mm+Mend1)
    END DO
    DO mn=1,Mnwv0
       dpln(mn)=der(mn)
       der(mn)=wcsa*der(mn)
    END DO
    l=0
    DO nn=1,Mend1
       DO mm=1,Mend2-nn
          l=l+1
          l0=la1(mm,nn)
          plnwcs(l)=an(mm)*pln(l0)
       END DO
    END DO
    DO mn=1,Mnwv0
       plnwcs(mn)=wcsa*plnwcs(mn)
    END DO
    DO nn=1,Mend1
       DO mm=1,Mend2-nn
          l0=la1(mm,nn)
          lp=la1(mm,nn+1)
          pln(l0)=x(l0)*x(lp)*raa*pln(l0)
       END DO
    END DO

  END SUBROUTINE plnder


  SUBROUTINE transs (Ldim, isign, a)

    IMPLICIT NONE

    !     transp: after input, transposes scalar arrays
    !             of spectral coefficients by swapping
    !             the order of the subscripts
    !             representing the degree and order
    !             of the associated legendre functions.
    ! 
    !     argument(dimensions)        description
    ! 
    !     Ldim                 input: number of layers.
    !     a(Mnwv2,Ldim)        input: spectral representation of a
    !                                 global field at "n" levels.
    !                                 isign=+1 diagonalwise storage
    !                                 isign=-1 coluMnwise   storage
    !                         output: spectral representation of a
    !                                 global field at "n" levels.
    !                                 isign=+1 coluMnwise   storage
    !                                 isign=-1 diagonalwise storage

    INTEGER, INTENT(IN) :: Ldim
    INTEGER, INTENT(IN) :: isign

    REAL (KIND=r8), INTENT(INOUT) :: a(Mnwv2,Ldim)

    REAL (KIND=r8) :: w(Mnwv2)

    INTEGER :: k
    INTEGER :: l
    INTEGER :: lx
    INTEGER :: mn
    INTEGER :: mm
    INTEGER :: nlast
    INTEGER :: nn

    IF (isign == 1) THEN
       DO k=1,Ldim
          l=0
          DO mm=1,Mend1
             nlast=Mend2-mm
             DO nn=1,nlast
                l=l+1
                lx=la0(mm,nn)
                w(2*l-1)=a(2*lx-1,k)
                w(2*l)=a(2*lx,k)
             END DO
          END DO
          DO mn=1,Mnwv2
             a(mn,k)=w(mn)
          END DO
       END DO
    ELSE
       DO k=1,Ldim
          l=0
          DO mm=1,Mend1
             nlast=Mend2-mm
             DO nn=1,nlast
                l=l+1
                lx=la0(mm,nn)
                w(2*lx-1)=a(2*l-1,k)
                w(2*lx)=a(2*l,k)
             END DO
          END DO
          DO mn=1,Mnwv2
             a(mn,k)=w(mn)
          END DO
       END DO
    END IF

  END SUBROUTINE transs


  SUBROUTINE SumSpec (Nmax, Mmax, Mnwv, Imax, Jmax, Jmaxhf, Kmax, &
                      Spec, Leg, Four, Len, LastPrev)

    IMPLICIT NONE

    ! fourier representation from spectral representation

    INTEGER, INTENT(IN) :: Nmax           ! # assoc leg func degrees (=trunc+1 or +2)
    INTEGER, INTENT(IN) :: Mmax           ! # assoc leg func waves (=trunc+1)
    INTEGER, INTENT(IN) :: Mnwv           ! # spectral coef, real+imag 
    INTEGER, INTENT(IN) :: Imax           ! # longitudes * 2 (real+imag four coef)
    INTEGER, INTENT(IN) :: Jmax           ! # latitudes (full sphere)
    INTEGER, INTENT(IN) :: Jmaxhf         ! # latitudes (hemisphere)
    INTEGER, INTENT(IN) :: Kmax           ! # verticals
    INTEGER, INTENT(IN) :: Len(Nmax)      ! diagonal length (real+imag)
    INTEGER, INTENT(IN) :: LastPrev(Nmax) ! last element previous diagonal (real+imag)

    REAL (KIND=r8), INTENT(IN) :: Spec(Mnwv,Kmax)       ! spectral field
    REAL (KIND=r8), INTENT(IN) :: Leg (Mnwv,Jmaxhf)     ! associated legendre function
    REAL (KIND=r8), INTENT(OUT) :: Four(Imax,Jmax,Kmax) ! full fourier field

    INTEGER :: j, k, ele, diag

    REAL (KIND=r8) :: OddDiag(2*Mmax,Jmaxhf,Kmax)
    REAL (KIND=r8) :: EvenDiag(2*Mmax,Jmaxhf,Kmax)

    ! initialize diagonals

    OddDiag=0.0_r8
    EvenDiag=0.0_r8

    ! sum odd diagonals (n+m even)

    DO k=1,Kmax
       DO j=1,Jmaxhf
          DO diag=1,Nmax,2
             DO ele=1,Len(diag)
                OddDiag(ele,j,k)=OddDiag(ele,j,k)+ &
                     Leg(ele+LastPrev(diag),j)*Spec(ele+LastPrev(diag),k)
             END DO
          END DO
       END DO
    END DO

    ! sum even diagonals (n+m odd)

    DO k=1,Kmax
       DO j=1,Jmaxhf
          DO diag=2,Nmax,2
             DO ele=1,Len(diag)
                EvenDiag(ele,j,k)=EvenDiag(ele,j,k)+ &
                     Leg(ele+LastPrev(diag),j)*Spec(ele+LastPrev(diag),k)
             END DO
          END DO
       END DO
    END DO

    ! use even-odd properties

    !$cdir nodep
    DO k=1,Kmax
       DO j=1,Jmaxhf
          DO ele=1,2*Mmax
             Four(ele,j,k)=OddDiag(ele,j,k)+EvenDiag(ele,j,k)
             Four(ele,Jmax+1-j,k)=OddDiag(ele,j,k)-EvenDiag(ele,j,k)
          END DO
       END DO
    END DO
    Four(2*Mmax+1:Imax,:,:) = 0.0_r8

  END SUBROUTINE SumSpec


  SUBROUTINE SumFour (Nmax, Mmax, Mnwv, Imax, Jmax, Jmaxhf, Kmax, &
                      Spec, Leg, Four, Len, LastPrev)

    IMPLICIT NONE

    ! spectral representation from fourier representation

    INTEGER, INTENT(IN) :: Nmax           ! # assoc leg func degrees (=trunc+1 or +2)
    INTEGER, INTENT(IN) :: Mmax           ! # assoc leg func waves (=trunc+1)
    INTEGER, INTENT(IN) :: Mnwv           ! # spectral coef, real+imag 
    INTEGER, INTENT(IN) :: Imax           ! # longitudes * 2 (real+imag four coef)
    INTEGER, INTENT(IN) :: Jmax           ! # latitudes (full sphere)
    INTEGER, INTENT(IN) :: Jmaxhf         ! # latitudes (hemisphere)
    INTEGER, INTENT(IN) :: Kmax           ! # verticals
    INTEGER, INTENT(IN) :: Len(Nmax)      ! diagonal length (real+imag)
    INTEGER, INTENT(IN) :: LastPrev(Nmax) ! last element previous diagonal (real+imag)

    REAL (KIND=r8), INTENT(OUT) :: Spec(Mnwv,Kmax)     ! spectral field
    REAL (KIND=r8), INTENT(IN) :: Leg (Mnwv,Jmaxhf)    ! associated legendre function
    REAL (KIND=r8), INTENT(IN) :: Four(Imax,Jmax,Kmax) ! full fourier field

    INTEGER :: j, jj, k, ele, diag

    REAL (KIND=r8), DIMENSION(2*Mmax,Jmaxhf,Kmax) :: FourEven, FourOdd

    ! initialize result

    Spec=0.0_r8

    ! use even-odd properties

    DO k=1,Kmax
       DO j=1,Jmaxhf
          jj=Jmax-j+1
          DO ele=1,2*Mmax
             FourEven(ele,j,k)=Four(ele,j,k)+Four(ele,jj,k)
             FourOdd (ele,j,k)=Four(ele,j,k)-Four(ele,jj,k)
          END DO
       END DO
    END DO

    ! sum odd diagonals (n+m even)

    DO k=1,Kmax
       DO j=1,Jmaxhf
          DO diag=1,Nmax,2
             DO ele=1,Len(diag)
                Spec(ele+LastPrev(diag),k)=Spec(ele+LastPrev(diag),k)+ &
                     FourEven(ele,j,k)*Leg(ele+LastPrev(diag),j)
             END DO
          END DO
       END DO
    END DO

    ! sum even diagonals (n+m odd)

    DO k=1,Kmax
       DO j=1,Jmaxhf
          DO diag=2,Nmax,2
             DO ele=1,Len(diag)
                Spec(ele+LastPrev(diag),k)=Spec(ele+LastPrev(diag),k)+ &
                     FourOdd(ele,j,k)*Leg(ele+LastPrev(diag),j)
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE SumFour


  SUBROUTINE Spec2Four2D (Spec, Four, der)

    IMPLICIT NONE

    REAL (KIND=r8), INTENT(IN) :: Spec(:,:)
    REAL (KIND=r8), INTENT(OUT) :: Four(:,:,:)

    LOGICAL, INTENT(IN), OPTIONAL :: der

    INTEGER :: s1, s2, f1, f2, f3

    LOGICAL :: extended, derivate

    CHARACTER (LEN=17) :: h="**(Spec2Four2D)**"

    IF (.NOT. created) THEN
       WRITE (UNIT=nferr, FMT='(2A)') h, &
             ' Module not created; invoke InitLegTrans prior to this call'
       STOP
    END IF

    s1=SIZE(Spec,1); s2=SIZE(Spec,2)
    f1=SIZE(Four,1); f2=SIZE(Four,2); f3=SIZE(Four,3)

    IF (s1 == Mnwv2) THEN
       extended=.FALSE.
    ELSE IF (s1 == Mnwv3) THEN
       extended=.TRUE.
    ELSE
       WRITE (UNIT=nferr, FMT='(2A,I10)') h, &
             ' wrong first dim of spec: ', s1
    END IF

    IF (s2 /= f3) THEN
       WRITE (UNIT=nferr, FMT= '(2A,2I10)') h, &
             ' vertical layers of spec and four dissagre :', s2, f3
       STOP
    END IF

    IF (f1 < 2*Mend1) THEN
       WRITE (UNIT=nferr, FMT= '(2A,2I10)') h, &
             ' first dimension of four too small: ', f1, 2*Mend1
       STOP
    END IF

    IF (f2 /= 2*jMaxHf) THEN
       WRITE (UNIT=nferr, FMT= '(2A,I10,A,I10)') h, &
             ' second dimension of four is ', f2, '; should be ', 2*jMaxHf
       STOP
    END IF

    IF (PRESENT(der)) THEN
       derivate=der
    ELSE
       derivate=.FALSE.
    END IF

    IF (derivate .AND. extended) THEN
       WRITE (UNIT=nferr, FMT= '(2A)') h, &
             ' derivative cannot be applied to extended gaussian field'
       STOP
    END IF

    IF (extended) THEN
       CALL SumSpec (Mend2, Mend1, Mnwv3, f1, f2, jMaxHf, f3, &
                     Spec, LegExtS2F, Four, LenDiagExt, LastPrevDiagExt)
    ELSE IF (derivate) THEN
       CALL SumSpec (Mend1, Mend1, Mnwv2, f1, f2, jMaxHf, f3, &
                     Spec, LegDerS2F, Four, LenDiag, LastPrevDiag)
    ELSE
       CALL SumSpec (Mend1, Mend1, Mnwv2, f1, f2, jMaxHf, f3, &
                     Spec, LegS2F, Four, LenDiag, LastPrevDiag)
    END IF

  END SUBROUTINE Spec2Four2D


  SUBROUTINE Spec2Four1D (Spec, Four, der)

    IMPLICIT NONE

    REAL (KIND=r8), INTENT(IN) :: Spec(:)
    REAL (KIND=r8), INTENT(OUT) :: Four(:,:)

    LOGICAL, INTENT(IN), OPTIONAL :: der

    INTEGER :: s1, f1, f2, f3

    LOGICAL :: extended, derivate

    CHARACTER (LEN=15) :: h="**(Spec2Four1D)**"

    IF (.NOT. created) THEN
       WRITE (UNIT=nferr, FMT='(2A)') h, &
             ' Module not created; invoke InitLegTrans prior to this call'
       STOP
    END IF

    s1=SIZE(Spec,1)
    f1=SIZE(Four,1); f2=SIZE(Four,2)

    IF (s1 == Mnwv2) THEN
       extended=.FALSE.
    ELSE IF (s1 == Mnwv3) THEN
       extended=.TRUE.
    ELSE
       WRITE (UNIT=nferr, FMT='(2A,I10)') h, &
             ' wrong first dim of spec : ', s1
    END IF

    IF (f1 < 2*Mend1) THEN
       WRITE (UNIT=nferr, FMT='(2A,2I10)') h, &
             ' first dimension of four too small: ', f1, 2*Mend1
       STOP
    END IF

    IF (f2 /= 2*jMaxHf) THEN
       WRITE (UNIT=nferr, FMT='(2A,I10,A,I10)') h, &
             ' second dimension of four is ', f2, '; should be ', 2*jMaxHf
       STOP
    END IF

    IF (PRESENT(der)) THEN
       derivate = der
    ELSE
       derivate = .FALSE.
    END IF

    IF (derivate .AND. extended) THEN
       WRITE (UNIT=nferr, FMT='(2A)') h, &
             ' derivative cannot be applied to extended gaussian field'
       STOP
    END IF

    f3=1
    IF (extended) THEN
       CALL SumSpec (Mend2, Mend1, Mnwv3, f1, f2, jMaxHf, f3, &
                     Spec, LegExtS2F, Four, LenDiagExt, LastPrevDiagExt)
    ELSE IF (derivate) THEN
       CALL SumSpec (Mend1, Mend1, Mnwv2, f1, f2, jMaxHf, f3, &
                     Spec, LegDerS2F, Four, LenDiag, LastPrevDiag)
    ELSE
       CALL SumSpec (Mend1, Mend1, Mnwv2, f1, f2, jMaxHf, f3, &
                     Spec, LegS2F, Four, LenDiag, LastPrevDiag)
    END IF

  END SUBROUTINE Spec2Four1D


  SUBROUTINE Four2Spec2D (Four, Spec)

    IMPLICIT NONE

    REAL (KIND=r8), INTENT(OUT) :: Spec(:,:)
    REAL (KIND=r8), INTENT(IN ) :: Four(:,:,:)

    INTEGER :: s1, s2, f1, f2, f3

    CHARACTER (LEN=17) :: h="**(Four2Spec2D)**"

    IF (.NOT. created) THEN
       WRITE (UNIT=nferr, FMT= '(2A)') h, &
             ' Module not created; invoke InitLegTrans prior to this call'
       STOP
    END IF

    s1=SIZE(Spec,1); s2=SIZE(Spec,2)
    f1=SIZE(Four,1); f2=SIZE(Four,2); f3=SIZE(Four,3)

    IF (s1 /= Mnwv2) THEN
       WRITE (UNIT=nferr, FMT='(2A,I10)') h, &
             ' wrong first dim of spec: ', s1
    END IF

    IF (s2 /= f3) THEN
       WRITE (UNIT=nferr, FMT='(2A,2I10)') h, &
             ' vertical layers of spec and four dissagre: ', s2, f3
       STOP
    END IF

    IF (f1 < 2*Mend1) THEN
       WRITE (UNIT=nferr, FMT='(2A,2I10)') h, &
             ' first dimension of four too small: ', f1, 2*Mend1
       STOP
    END IF

    IF (f2 /= 2*jMaxHf) THEN
       WRITE (UNIT=nferr, FMT='(2A,I10,A,I10)') h, &
             ' second dimension of four is ', f2, '; should be ', 2*jMaxHf
       STOP
    END IF

    CALL SumFour (Mend1, Mend1, Mnwv2, f1, f2, jMaxHf, f3, &
                  Spec, LegF2S, Four, LenDiag, LastPrevDiag)

  END SUBROUTINE Four2Spec2D


  SUBROUTINE Four2Spec1D (Four, Spec)

    IMPLICIT NONE

    REAL (KIND=r8), INTENT(OUT) :: Spec(:)
    REAL (KIND=r8), INTENT(IN) :: Four(:,:)

    INTEGER :: s1, f1, f2, f3

    CHARACTER (LEN=17) :: h="**(Four2Spec1D)**"

    IF (.NOT. created) THEN
       WRITE (UNIT=nferr, FMT='(2A)') h, &
             ' Module not created; invoke InitLegTrans prior to this call'
       STOP
    END IF

    s1=SIZE(Spec,1)
    f1=SIZE(Four,1); f2=SIZE(Four,2)

    IF (s1 /= Mnwv2) THEN
       WRITE (UNIT=nferr, FMT='(2A,I10)') h, &
             ' wrong first dim of spec: ', s1
    END IF

    IF (f1 < 2*Mend1) THEN
       WRITE (UNIT=nferr, FMT='(2A,2I10)') h, &
             ' first dimension of four too small: ', f1, 2*Mend1
       STOP
    END IF

    IF (f2 /= 2*jMaxHf) THEN
       WRITE (UNIT=nferr, FMT='(2A,I10,A,I10)') h, &
             ' second dimension of four is ', f2, '; should be ', 2*jMaxHf
       STOP
    END IF

    f3=1
    CALL SumFour (Mend1, Mend1, Mnwv2, f1, f2, jMaxHf, f3, &
                  Spec, LegF2S, Four, LenDiag, LastPrevDiag)

  END SUBROUTINE Four2Spec1D


  SUBROUTINE SplitTrans3D (full, north, south)

    IMPLICIT NONE

    REAL (KIND=r8), INTENT(IN) :: full (:,:,:)
    REAL (KIND=r8), INTENT(OUT) :: north(:,:,:)
    REAL (KIND=r8), INTENT(OUT) :: south(:,:,:)

    INTEGER :: if1, in1, is1
    INTEGER :: if2, in2, is2
    INTEGER :: if3, in3, is3
    INTEGER :: i, j, k

    CHARACTER (LEN=18), PARAMETER :: h="**(SplitTrans3D)**"

    if1=SIZE(full,1); in1=SIZE(north,1); is1=SIZE(south,1)
    if2=SIZE(full,2); in2=SIZE(north,2); is2=SIZE(south,2)
    if3=SIZE(full,3); in3=SIZE(north,3); is3=SIZE(south,3)

    IF (in1 /= is1) THEN
       WRITE (UNIT=nferr, FMT='(2A,2I10)') h, &
             ' dim 1 of north and south dissagree: ', in1, is1
       STOP
    END IF
    IF (in2 /= is2) THEN
       WRITE (UNIT=nferr, FMT='(2A,2I10)') h, &
             ' dim 2 of north and south dissagree: ', in2, is2
       STOP
    END IF
    IF (in3 /= is3) THEN
       WRITE (UNIT=nferr, FMT='(2A,2I10)') h, &
             ' dim 3 of north and south dissagree: ', in3, is3
       STOP
    END IF

    IF (in1 < if1) THEN
       WRITE (UNIT=nferr, FMT='(2A,2I10)') h, &
             ' first dimension of north too small: ', in1, if1
       STOP
    END IF
    IF (if2 /= 2*in3) THEN
       WRITE (UNIT=nferr, FMT='(2A,2I10)') h, &
             ' second dimension of full /= 2*third dimension of north: ', if2, 2*in3
       STOP
    END IF
    IF (if3 /= in2) THEN
       WRITE (UNIT=nferr, FMT='(2A,2I10)') h, &
             ' second dimension of north and third dimension of full dissagree: ', in2, if3
       STOP
    END IF

    north=0.0_r8
    south=0.0_r8
    DO k=1,in2
       DO j=1,in3
          DO i=1,if1
             north(i,k,j)=full(i,j,k)
             south(i,k,j)=full(i,if2-j+1,k)
          END DO
       END DO
    END DO

  END SUBROUTINE SplitTrans3D


  SUBROUTINE SplitTrans2D (full, north, south)

    IMPLICIT NONE

    REAL (KIND=r8), INTENT(IN) :: full (:,:)
    REAL (KIND=r8), INTENT(OUT) :: north(:,:)
    REAL (KIND=r8), INTENT(OUT) :: south(:,:)

    CHARACTER (LEN=18), PARAMETER :: h="**(SplitTrans2D)**"

    INTEGER :: if1, in1, is1
    INTEGER :: if2, in2, is2
    INTEGER :: i, j

    if1=SIZE(full,1); in1=SIZE(north,1); is1=SIZE(south,1)
    if2=SIZE(full,2); in2=SIZE(north,2); is2=SIZE(south,2)

    IF (in1 /= is1) THEN
       WRITE (UNIT=nferr, FMT='(2A,2I10)') h, &
             ' dim 1 of north and south dissagree: ', in1, is1
       STOP
    END IF
    IF (in2 /= is2) THEN
       WRITE (UNIT=nferr, FMT='(2A,2I10)') h, &
             ' dim 2 of north and south dissagree: ', in2, is2
       STOP
    END IF

    IF (in1 < if1) THEN
       WRITE (UNIT=nferr, FMT='(2A,2I10)') h, &
             ' first dimension of north too small: ', in1, if1
       STOP
    END IF
    IF (if2 /= 2*in2) THEN
       WRITE (UNIT=nferr, FMT='(2A,2I10)') h, &
             ' second dimension of full /= 2*second dimension of north: ', if2, 2*in2
       STOP
    END IF

    north=0.0_r8
    south=0.0_r8
    DO j=1,in2
       DO i=1,if1
          north(i,j)=full(i,j)
          south(i,j)=full(i,if2-j+1)
       END DO
    END DO

  END SUBROUTINE SplitTrans2D


END MODULE LegendreTransform
