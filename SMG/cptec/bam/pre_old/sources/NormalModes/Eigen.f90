!
!  $Author: bonatti $
!  $Date: 2009/03/20 18:00:00 $
!  $Revision: 1.1.1.1 $
!
MODULE Eigen

   USE InputParameters, ONLY: r8, nferr

   IMPLICIT NONE

   PRIVATE

   !REAL (KIND=r8) :: TolH, TolW, TolZ
   REAL (KIND=r8) :: TolW, TolZ

   PUBLIC :: GetEigenRGM, tql2, tred2


CONTAINS


SUBROUTINE GetEigenRGM (nm, n, a, wr, wi, MakeVec, z, iErr, Eps)
 
   ! Calculates the Eigenvalues and/or Eigenvectors of a
   ! Real General Matrix

   ! Arguments:

   ! nm - Row dimension of matrix A at the calling routine: Input
   !      Integer variable

   ! n - Current dimension of matrix A: Input
   !     Integer variable; must be <= nm

   ! a - Real matrix (destroyed, modified at Output): Input
   !     Real array with dimensions a(nm,n)

   ! wr - Real part of the eigenvalues: Output
   !      Real vector with dimensions wr(n)

   ! wi - Imaginary part of the eigenvalues: Output
   !      Real vector with dimensions wi(n)

   ! OBS: There is no ordenation for the eigenvalues, except
   !      that for the conjugate complex pairs are put together
   !      and the pair with real positive imaginary part comes
   !      in first place.

   ! MakeVec - Integer variable: Input
   !         = 0 - only eigenvalues non-filtered
   !         = 1 - eigenvalues and eigenvectors normalized and filtered
   !         = 2 - eigenvalues and eigenvectors non-norm. and non-filt
   !         = 3 - eigenvalues and eigenvectors normalized, non-filt
   !               and without zeroes for <= TolX (see bellow and zNorma)

   ! z - Eigenvectors: real and imaginary parts, so that:
   !     a) for a real j-th eigenvalue wr(j) /= 0 .AND. wi(j) == 0
   !        the j-th eigenvector is (z(i,j),i=1,n)
   !     b) for a imaginary j-th eigenvalue with wi(j) /= 0,
   !        the (j+1)-th eigenvalue is its conjugate complex,
   !        the j-th eigenvector has real part (z(i,j),i=1,n) and
   !        imaginary part (z(i,j+1),i=1,n), and the (j+1)-th
   !        eigenvector has real part (z(i,j),i=1,n) and
   !        imaginary part (-z(i,j+1),i=1,n).
   !        Real array with dimensions z(nm,n): Output

   ! iErr - is a integer variable: Output, indicating:
   !        - if n > nm, then the routine stop calculations
   !          and returns with iErr=10*n
   !        - if 50 iteractions is exceeded for the j-th eigenvalue
   !          computation, then the routine stop calculations
   !          and returns with iErr=j and the j+1, j+2, ..., n
   !          eigenvalues are computed, but none eigenvector is computed
   !        - for a normal termination iErr is set 0.

   ! Eps - is a machine dependent parameter specifying the
   !       relative precision of the floating point arithmetic.
   !       It must be recomputed for the specific machine in use.
   !       Eps=EPISILON(Eps)
   !       Real variable: Input

   ! TolH - Tolerance value to filter the Hessemberg matrix
   !        Real variable: Local

   ! TolW - Tolerance value to filter the Eigenvalues
   !        Real variable: Local

   ! TolZ - Tolerance value to filter the Eigenvectors
   !        Real variable: Local

   IMPLICIT NONE

   INTEGER, INTENT (IN) :: nm
   INTEGER, INTENT (IN) :: n
   INTEGER, INTENT (IN) :: MakeVec

   REAL (KIND=r8), INTENT (IN) :: Eps

   REAL (KIND=r8), INTENT (IN OUT) :: a(nm,n)

   REAL (KIND=r8), INTENT (OUT) :: wr(n)
   REAL (KIND=r8), INTENT (OUT) :: wi(n)
   REAL (KIND=r8), INTENT (OUT) :: z(nm,n)

   INTEGER, INTENT (OUT) :: iErr

   INTEGER :: Low, High

   REAL (KIND=r8) :: Scale(n), Ort(n)

   !TolH=Eps
   TolW=1.0E-12_r8
   TolZ=1.0E-12_r8
   iErr=0
   IF (n > nm) THEN
      iErr=10*n
      RETURN
   END IF

   ! Performing the Balance of the input real general matrix in place.
   CALL Balanc (nm, n, a, Low, High, Scale)

   ! Performing the redution of the Balanced matrix 
   ! in place to the Hessemberg superior form.
   ! It is used similarity Orthogonal transformations.
   CALL Orthes (nm, n, Low, High, a, Ort)

   IF (MakeVec /= 0) THEN
      ! Saving the transformations above for eigenvector computations.
      CALL Ortran (nm, n, Low, High, a, Ort, z)
   END IF

   ! Computing the eigenvalues/eigenvectors of the 
   ! Hessemberg matrix using the QR method.
   CALL HQRValue (nm, n, Low, High, a, wr, wi, z, iErr, MakeVec, Eps)

   IF (MakeVec /= 0 .AND. iErr == 0) THEN

     ! Back-transforming the eigenvectors of the Hessembeg matrix 
     ! to the eigenvectors of the original input matrix.
      CALL BalBak (nm, n, Low, High, Scale, n, z)

      ! Normalizing and filtering the eigenvectors.
      ! See MakeVec above and comments inside zNorma routine.
      CALL zNorma (nm, n, wr, wi, z, MakeVec)

   END IF

   IF (iErr /= 0) WRITE (UNIT=nferr, FMT='(/,A,I5,A,/)') &
      ' ***** The ', iErr, '-th Eigenvalue Did Not Converge *****'

END SUBROUTINE GetEigenRGM


SUBROUTINE Balanc (nm, n, a, Low, High, Scale)
 
   ! Balances a real general matrix, and
   ! isolates eigenvalues whenever possible.

   IMPLICIT NONE

   INTEGER, INTENT (IN) :: nm
   INTEGER, INTENT (IN) :: n

   REAL (KIND=r8), INTENT (IN OUT) :: a(nm,*)

   INTEGER, INTENT (OUT) :: Low
   INTEGER, INTENT (OUT) :: High

   REAL (KIND=r8), INTENT (OUT) :: Scale(*)

   INTEGER :: i, j, k, l, m, jj, iexc

   LOGICAL :: NoConv

   REAL (KIND=r8) :: c, f, g, r, s, b2, Radix

   ! Radix is a machine dependent parameter specifying
   !       the base of the machine floating pont representation.

   Radix=2.0_r8

   b2=Radix*Radix
   k=1
   l=n
   GO TO 100

   ! In-line procedure for row and column exchange.

   20 Scale(m)=REAL(j,r8)
   IF (j == m) GO TO 50

   DO i=1,l
      f=a(i,j)
      a(i,j)=a(i,m)
      a(i,m)=f
   END DO

   DO i=k,n
      f=a(j,i)
      a(j,i)=a(m,i)
      a(m,i)=f
   END DO

   50 SELECT CASE (iexc)
      CASE (1)
         GO TO 80
      CASE (2)
         GO TO 130
   END SELECT

   ! Search for rows isolating an eigenvalue and push them down.

   80 IF (l == 1) GO TO 280
   l=l-1

   ! For j = l  step -1 until 1 DO -- .

   100 CONTINUE
   loop120: DO jj=1,l
      j=l+1-jj

      DO i=1,l
         IF (i == j) CYCLE
         IF (a(j,i) /= 0.0_r8) CYCLE loop120
      END DO

      m=l
      iexc=1
      GO TO 20
   END DO loop120

   GO TO 140

   ! Search for columns isolating an eigenvalue and push them left.

   130 k=k+1

   140 CONTINUE
   loop170: DO j=k,l

      DO i=k,l
         IF (i == j) CYCLE
         IF (a(i,j) /= 0.0_r8) CYCLE loop170
      END DO

      m=k
      iexc=2
      GO TO 20
   END DO loop170

   ! Now Balance the submatrix in rows k to l.

   DO i=k,l
      Scale(i)=1.0_r8
   END DO

   ! Interative loop for norm reduction.

   190 NoConv=.FALSE.

   DO i=k,l
      c=0.0_r8
      r=0.0_r8

      DO j=k,l
         IF (j == i) CYCLE
         c=c+ABS(a(j,i))
         r=r+ABS(a(i,j))
      END DO

      ! Guard against zero c or r due to underflow.

      IF (c == 0.0_r8 .OR. r == 0.0_r8) CYCLE
      g=r/Radix
      f=1.0_r8
      s=c+r
      DO WHILE (c < g)
         f=f*Radix
         c=c*b2
      END DO
      g=r*Radix
      DO WHILE (c >= g)
         f=f/Radix
         c=c/b2
      END DO

      ! Now Balance

      IF ((c+r)/f >= 0.95_r8*s) CYCLE
      g=1.0_r8/f
      Scale(i)=Scale(i)*f
      NoConv=.TRUE.

      DO j=k,n
         a(i,j)=a(i,j)*g
      END DO

      DO j=1,l
         a(j,i)=a(j,i)*f
      END DO

   END DO

   IF (NoConv) GO TO 190

   280 Low=k
   High=l

END SUBROUTINE Balanc


SUBROUTINE OrtHes (nm, n, Low, High, a, Ort)
 
! Reduces a Real General Matrix to Upper
! Hessemberg Form using Orthogonal Similarity

   IMPLICIT NONE

   INTEGER, INTENT (IN) :: nm
   INTEGER, INTENT (IN) :: n
   INTEGER, INTENT (IN) :: Low
   INTEGER, INTENT (IN) :: High

   REAL (KIND=r8), INTENT (IN OUT) :: a(nm,n)

   REAL (KIND=r8), INTENT (OUT) :: Ort(High)

   INTEGER :: i, j, m, ii, jj, la, mp, kp1

   REAL (KIND=r8) :: f, g, h, Scal

   la=High-1
   kp1=Low+1
   IF (la < kp1) RETURN

   DO m=kp1,la
      h=0.0_r8
      Ort(m)=0.0_r8
      Scal=0.0_r8

      ! Scale column

      DO i=m,High
         Scal=Scal+ABS(a(i,m-1))
      END DO

      IF (Scal == 0.0_r8) CYCLE
      mp=m+High

      ! For i = High step -1 until m DO --

      DO ii=m,High
         i=mp-ii
         Ort(i)=a(i,m-1)/Scal
         h=h+Ort(i)*Ort(i)
      END DO

      g=-SIGN(SQRT(h),Ort(m))
      h=h-Ort(m)*g
      Ort(m)=Ort(m)-g

      ! Form (I-(u*ut)/h) * a

      DO j=m,n
         f=0.0_r8

         ! For i = High step -1 until m DO --

         DO ii=m,High
            i=mp-ii
            f=f+Ort(i)*a(i,j)
         END DO

         f=f/h

         DO i=m,High
            a(i,j)=a(i,j)-f*Ort(i)
            !IF (ABS(a(i,j)) < TolH) a(i,j)=TolH
         END DO

      END DO

      ! Form (I-(u*ut)/h) * a * (I-(u*ut)/h)

      DO i=1,High
         f=0.0_r8

         ! For j = High step -1 until m DO --

         DO jj=m,High
            j=mp-jj
            f=f+Ort(j)*a(i,j)
         END DO

         f=f/h

         DO j=m,High
            a(i,j)=a(i,j)-f*Ort(j)
            !IF (ABS(a(i,j)) < TolH) a(i,j)=TolH
         END DO

      END DO

      Ort(m)=Scal*Ort(m)
      a(m,m-1)=Scal*g
      !IF (ABS(a(m,m-1)) < TolH) a(m,m-1)=TolH
   END DO

END SUBROUTINE OrtHes


SUBROUTINE Ortran (nm, n, Low, High, a, Ort, z)
 
   ! Accumulates the Orthogonal Similarity Tranformations
   ! used in the Reduction of a Real General Matrix to
   ! Upper Hessemberg Form by OrtHes.

   IMPLICIT NONE

   INTEGER, INTENT (IN) :: nm
   INTEGER, INTENT (IN) :: n
   INTEGER, INTENT (IN) :: Low
   INTEGER, INTENT (IN) :: High

   REAL (KIND=r8), INTENT (IN) :: a(nm,n)

   REAL (KIND=r8), INTENT (OUT) :: Ort(n)

   REAL (KIND=r8), INTENT (OUT) :: z(nm,n)

   INTEGER :: i, j, kl, mm, mp, mp1

   REAL (KIND=r8) :: g

   ! Initialize z to identity matrix

   DO i=1,n
      DO j=1,n
         z(i,j)=0.0_r8
      END DO
      z(i,i)=1.0_r8
   END DO

   kl=High-Low-1
   IF (kl < 1) RETURN

   ! For mp = High-1 step -1 until Low+1 DO --

   DO mm=1,kl
      mp=High-mm
      IF (a(mp,mp-1) == 0.0_r8) CYCLE
      mp1=mp+1

      DO i=mp1,High
         Ort(i)=a(i,mp-1)
      END DO

      DO j=mp,High
         g=0.0_r8

         DO i=mp,High
            g=g+Ort(i)*z(i,j)
         END DO

         ! Divisor below is negative of h formed in OrtHes
         ! double division avoids possible underflow.

         g=(g/Ort(mp))/a(mp,mp-1)

         DO i=mp,High
            z(i,j)=z(i,j)+g*Ort(i)
            !IF (ABS(z(i,j)) < TolH) z(i,j)=TolH
         END DO

      END DO

   END DO

END SUBROUTINE Ortran


SUBROUTINE HQRValue (nm, n, Low, High, h, wr, wi, z, &
                     iErr, MakeVec, Eps)
 
   ! Computes the Eigenvalues and/or Eigenvectors of a
   ! Real Upper Hessemberg Matrix using the QR Method.

   ! Eps is a machine dependent parameter specifying the
   ! relative precision of the floating point arithmetic.
   ! It must be recomputed and replaced for the specific
   ! machine in use.

   IMPLICIT NONE

   INTEGER, INTENT (IN) :: nm
   INTEGER, INTENT (IN) :: n
   INTEGER, INTENT (IN) :: Low
   INTEGER, INTENT (IN) :: High
   INTEGER, INTENT (IN) :: MakeVec

   REAL (KIND=r8), INTENT (IN) :: Eps

   REAL (KIND=r8), INTENT (IN OUT) :: h(nm,n)
   REAL (KIND=r8), INTENT (IN OUT) :: z(nm,n)

   INTEGER, INTENT (OUT) :: iErr

   REAL (KIND=r8), INTENT (OUT) :: wr(n)
   REAL (KIND=r8), INTENT (OUT) :: wi(n)

   INTEGER :: i, j, k, l, m, en, ll, mm, na, its, mp2, enm2

   REAL (KIND=r8) :: p, q, r, s, t, x, w, y, zz, Norm

   LOGICAL :: notlas

   iErr=0
   Norm=0.0_r8
   k=1

   ! Store Roots Isolated by Balanc and Compute Matrix Norm

   DO i=1,n

      DO j=k,n
         !IF (ABS(h(i,j)) < TolH) h(i,j)=TolH
         Norm=Norm+ABS(h(i,j))
      END DO

      k=i
      IF (i >= Low .AND. i <= High) CYCLE
      wr(i)=h(i,i)
      wi(i)=0.0_r8
   END DO

   en=High
   t=0.0_r8

   ! Search for Next Eigenvalues.

   60 IF (en < Low) GO TO 340
   its=0
   na=en-1
   enm2=na-1

   ! Look for Single Small Sub-Diagonal Element
   ! For l=en Step -1 Until Low DO --

   70 DO ll=Low,en
      l=en+Low-ll
      IF (l == Low) EXIT
      s=ABS(h(l-1,l-1))+ABS(h(l,l))
      IF (s == 0.0_r8) s=Norm
      IF (ABS(h(l,l-1)) <= Eps*s) EXIT
   END DO

   ! Form Shift

   100 x=h(en,en)
   IF (l == en) GO TO 270
   y=h(na,na)
   w=h(en,na)*h(na,en)
   IF (l == na) GO TO 280

   IF (its == 30) THEN
      ! Set Error - No Convergence to an Eigenvalue after 30 Iterations
      iErr=en
      RETURN
   END IF
   IF (its /= 10 .AND. its /= 20) GO TO 130

   ! Form Exceptional Shift

   t=t+x

   DO i=Low,en
      h(i,i)=h(i,i)-x
   END DO

   s=ABS(h(en,na))+ABS(h(na,enm2))
   x=0.75_r8*s
   y=x
   w=-0.4375_r8*s*s
   130 its=its+1

   ! Look for Two Consecutive Small Sub-diagonal Elements
   ! For m=en-2 step -1 until l DO --

   DO mm=l,enm2
      m=enm2+l-mm
      zz=h(m,m)
      r=x-zz
      s=y-zz
      p=(r*s-w)/h(m+1,m)+h(m,m+1)
      q=h(m+1,m+1)-zz-r-s
      r=h(m+2,m+1)
      s=ABS(p)+ABS(q)+ABS(r)
      p=p/s
      q=q/s
      r=r/s
      IF (m == l) EXIT
      IF (ABS(h(m,m-1))*(ABS(q)+ABS(r)) <= Eps*ABS(p)*  &
         (ABS(h(m-1,m-1))+ABS(zz)+ABS(h(m+1,m+1)))) EXIT
   END DO

   150 mp2=m+2

      DO i=mp2,en
      h(i,i-2)=0.0_r8
      IF (i == mp2) CYCLE
      h(i,i-3)=0.0_r8
   END DO

   ! Double QR Step Involving Rows l to End and Columns m to en

   DO k=m,na
      notlas=k /= na
      IF (k /= m) THEN
         p=h(k,k-1)
         q=h(k+1,k-1)
         r=0.0_r8
         IF(notlas) r=h(k+2,k-1)
         x=ABS(p)+ABS(q)+ABS(r)
         IF (x == 0.0_r8) CYCLE
         p=p/x
         q=q/x
         r=r/x
      END IF
      s=SIGN(SQRT(p*p+q*q+r*r),p)
      IF (k /= m) THEN
         h(k,k-1)=-s*x
      ELSE
         IF (l /= m) h(k,k-1)=-h(k,k-1)
      END IF
      p=p+s
      x=p/s
      y=q/s
      zz=r/s
      q=q/p
      r=r/p

      ! Row Modification

      DO j=k,n
         p=h(k,j)+q*h(k+1,j)
         IF (notlas) THEN
            p=p+r*h(k+2,j)
            h(k+2,j)=h(k+2,j)-p*zz
         END IF
         h(k+1,j)=h(k+1,j)-p*y
         h(k,j)=h(k,j)-p*x
      END DO

      j=MIN(en,k+3)

      ! Column Modification

      DO i=1,j
         p=x*h(i,k)+y*h(i,k+1)
         IF (notlas) THEN
            p=p+zz*h(i,k+2)
            h(i,k+2)=h(i,k+2)-p*r
         END IF
         h(i,k+1)=h(i,k+1)-p*q
         h(i,k)=h(i,k)-p
      END DO

      IF (MakeVec == 0) CYCLE

      ! Accumulate Transformations

      DO i=Low,High
         p=x*z(i,k)+y*z(i,k+1)
         IF (notlas) THEN
            p=p+zz*z(i,k+2)
            z(i,k+2)=z(i,k+2)-p*r
         END IF
         !IF (ABS(p) < TolH) p=TolH
         z(i,k+1)=z(i,k+1)-p*q
         z(i,k)=z(i,k)-p
      END DO

   END DO

   GO TO 70

   ! One Root Found

   270 h(en,en)=x+t
   wr(en)=h(en,en)
   wi(en)=0.0_r8
   en=na
   GO TO 60

   ! Two Roots Found

   280 p=(y-x)*0.5_r8
   q=p*p+w
   zz=SQRT(ABS(q))
   h(en,en)=x+t
   x=h(en,en)
   h(na,na)=y+t
   IF (q < 0.0_r8) GO TO 320

   ! Real Pair

   zz=p+SIGN(zz,p)
   wr(na)=x+zz
   wr(en)=wr(na)
   IF (zz /= 0.0_r8) wr(en)=x-w/zz
   wi(na)=0.0_r8
   wi(en)=0.0_r8

   IF (MakeVec == 0) GO TO 330

   x=h(en,na)
   s=ABS(x)+ABS(zz)
   p=x/s
   q=zz/s
   r=SQRT(p*p+q*q)
   p=p/r
   q=q/r

   ! Row Modification

   DO j=na,n
      zz=h(na,j)
      h(na,j)=q*zz+p*h(en,j)
      h(en,j)=q*h(en,j)-p*zz
   END DO

   ! Column Modification

   DO i=1,en
      zz=h(i,na)
      h(i,na)=q*zz+p*h(i,en)
      h(i,en)=q*h(i,en)-p*zz
   END DO

   ! Accumulate Transformations

   DO i=Low,High
      zz=z(i,na)
      z(i,na)=q*zz+p*z(i,en)
      z(i,en)=q*z(i,en)-p*zz
   END DO

   GO TO 330

   ! Complex Pair

   320 wr(na)=x+p
   wr(en)=x+p
   wi(na)=zz
   wi(en)=-zz

   330 en=enm2

   GO TO 60

   340 IF (MakeVec == 0) RETURN

   ! All Roots Found

   ! Backsubstitute to Find Vectors of Upper Triangular Form

   IF (Norm /= 0.0_r8) &
      CALL HQRVector (nm, n, Low, High, h, wr, wi, z, Eps, Norm)

END SUBROUTINE HQRValue


SUBROUTINE HQRVector (nm, n, Low, High, h, wr, wi, z, Eps, Norm)
 
   ! Backsubstitutes to Find Vectors of Upper Triangular Form

   IMPLICIT NONE

   INTEGER, INTENT (IN) :: nm
   INTEGER, INTENT (IN) :: n
   INTEGER, INTENT (IN) :: Low
   INTEGER, INTENT (IN) :: High

   REAL (KIND=r8), INTENT (IN) :: Eps
   REAL (KIND=r8), INTENT (IN) :: Norm
   REAL (KIND=r8), INTENT (IN) :: wr(n)
   REAL (KIND=r8), INTENT (IN) :: wi(n)

   REAL (KIND=r8), INTENT (IN OUT) :: h(nm,n)

   REAL (KIND=r8), INTENT (OUT) :: z(nm,n)

   INTEGER :: i, j, k, m, en, ii, jj, na, nn, enm2

   REAL (KIND=r8) :: p, q, r, s, t, x, w, y, &
                     ar, ai, br, bi, ra, sa, zz

   ! For en=n Step -1 Until 1 DO --

   DO nn=1,n

      en=n+1-nn
      p=wr(en)
      q=wi(en)
      na=en-1

      IF (q == 0.0_r8) THEN

         ! Real Vector

         m=en
         h(en,en)=1.0_r8
         IF (na == 0) CYCLE

         ! For i=en-1 Step -1 Until 1 DO -- .

         DO ii=1,na
            i=en-ii
            w=h(i,i)-p
            r=h(i,en)

            IF (m <= na) THEN
               DO j=m,na
                  r=r+h(i,j)*h(j,en)
               END DO
            END IF

            IF (wi(i) < 0.0_r8) THEN
               zz=w
               s=r
               CYCLE
            END IF
            m=i
            IF (wi(i) == 0.0_r8) THEN
               t=w
               IF (w == 0.0_r8) t=Eps*Norm
               h(i,en)=-r/t
               CYCLE
            END IF

            ! Solve Real Equations

            x=h(i,i+1)
            y=h(i+1,i)
            q=(wr(i)-p)*(wr(i)-p)+wi(i)*wi(i)
            t=(x*s-zz*r)/q
            h(i,en)=t
            IF (ABS(x) > ABS(zz)) THEN
               h(i+1,en)=(-r-w*t)/x
               CYCLE
            END IF
            h(i+1,en)=(-s-y*t)/zz

         END DO

         ! End Real Vector

         CYCLE

      ELSE IF (q < 0.0_r8) THEN

         ! Complex Vector

         m=na

         ! Last Vector Component Chosen Imaginary
         ! so that Eigenvector Matrix is Triangular

         IF (ABS(h(en,na)) > ABS(h(na,en))) THEN
            h(na,na)=q/h(en,na)
            h(na,en)=-(h(en,en)-p)/h(en,na)
         ELSE
            h(na,na)=CmplxDiv(0.0_r8,-h(na,en),h(na,na)-p,q)
            h(na,en)=CmplxDiv(-h(na,en),0.0_r8,h(na,na)-p,q)
         END IF
         h(en,na)=0.0_r8
         h(en,en)=1.0_r8
         enm2=na-1
         IF (enm2 == 0) CYCLE

         ! For i=en-2 step -1 Until 1 DO --

         DO ii=1,enm2
            i=na-ii
            w=h(i,i)-p
            ra=0.0_r8
            sa=h(i,en)

            DO j=m,na
               ra=ra+h(i,j)*h(j,na)
               sa=sa+h(i,j)*h(j,en)
            END DO

            IF (wi(i) < 0.0_r8) THEN
               zz=w
               r=ra
               s=sa
               CYCLE
            END IF
            m=i
            IF (wi(i) == 0.0_r8) THEN
               h(i,na)=CmplxDiv(-ra,-sa,w,q)
               h(i,en)=CmplxDiv(-sa,ra,w,q)
            CYCLE
            END IF

            ! Solve Complex Equations

            x=h(i,i+1)
            y=h(i+1,i)
            ar=x*r-zz*ra+q*sa
            ai=x*s-zz*sa-q*ra
            br=(wr(i)-p)*(wr(i)-p)+wi(i)*wi(i)-q*q
            bi=(wr(i)-p)*2.0_r8*q
            IF (br == 0.0_r8 .AND. bi == 0.0_r8) br=Eps*Norm*  &
               (ABS(w)+ABS(q)+ABS(x)+ABS(y)+ABS(zz))
            h(i,na)=CmplxDiv(ar,ai,br,bi)
            h(i,en)=CmplxDiv(ai,-ar,br,bi)
            IF (ABS(x) > (ABS(zz)+ABS(q))) THEN
               h(i+1,na)=(-ra-w*h(i,na)+q*h(i,en))/x
               h(i+1,en)=(-sa-w*h(i,en)-q*h(i,na))/x
               CYCLE
            END IF
            h(i+1,na)=CmplxDiv(-r-y*h(i,na),-s-y*h(i,en),zz,q)
            h(i+1,en)=CmplxDiv(-s-y*h(i,en),r+y*h(i,na),zz,q)
         END DO

         ! End Complex Vector

      ELSE

         CYCLE

      END IF

   END DO

   ! End Back Substitution

   ! Vectors of Isolated Roots

   DO i=1,n
      IF (i >= Low .AND. i <= High) CYCLE
      DO j=i,n
        z(i,j)=h(i,j)
      END DO
   END DO

   ! Multiply by Transformations Matrix to give
   ! Vectors of Original Full Matrix
   ! For j=n Step -1 Until Low DO --

   DO jj=Low,n
      j=n+Low-jj
      m=MIN(j,High)
      DO i=Low,High
         zz=0.0_r8
         DO k=Low,m
            zz=zz+z(i,k)*h(k,j)
         END DO
         z(i,j)=zz
      END DO
   END DO

END SUBROUTINE HQRVector


REAL (KIND=r8) FUNCTION CmplxDiv (z1, z2, z3, z4)

   IMPLICIT NONE

   REAL (KIND=r8), INTENT (IN) :: z1, z2, z3, z4
 
   CmplxDiv=(z1*z3+z2*z4)/(z3*z3+z4*z4)

END FUNCTION CmplxDiv


SUBROUTINE BalBak (nm, n, Low, High, Scale, m, z)
 
   ! Forms the Eigenvectors of a Real General Matrix from
   ! the Eigenvectors of that Matrix Transformed by Balanc.

   IMPLICIT NONE

   INTEGER, INTENT (IN) :: nm
   INTEGER, INTENT (IN) :: n
   INTEGER, INTENT (IN) :: Low
   INTEGER, INTENT (IN) :: High
   INTEGER, INTENT (IN) :: m

   REAL (KIND=r8), INTENT (IN) :: Scale(n)

   REAL (KIND=r8), INTENT (IN OUT) :: z(nm,n)

   INTEGER :: i ,j ,k ,ii

   REAL (KIND=r8) :: s

   IF (m == 0) RETURN

   IF (High /= Low) THEN

      DO i=Low,High
         s=Scale(i)

         ! Left hand eigenvectors are back transformed if the
         ! foregoing statment is replaced by s=1.0/Scale(i)

         DO j=1,m
            z(i,j)=z(i,j)*s
         END DO

      END DO

   END IF

   ! For I=Low-1  step -1 until 1,
   !       High+1 step  1 until n DO --

   DO ii=1,n
      i=ii
      IF (i >= Low .AND. i <= High) CYCLE
      IF (i < Low) i=Low-ii
      k=INT(Scale(i))
      IF (k == i) CYCLE

      DO j=1,m
         s=z(i,j)
         z(i,j)=z(k,j)
         z(k,j)=s
      END DO

   END DO

END SUBROUTINE BalBak


SUBROUTINE zNorma (nm, n, wr, wi, z, MakeVec)
 
   ! Normalizes and filters the eigenvectors and 
   ! filters the eigenvalues.

   ! It sets zz=a+b*i, corresponding to the maximum
   ! absolute value of the eigenvector, to:

   ! a) 1         - if b == 0
   ! b) 1+(b/a)*i - if  ABS(a) >= ABS(b)
   ! c) 1+(a/b)*i - if  ABS(a) < ABS(b)

   IMPLICIT NONE

   INTEGER, INTENT (IN) :: nm
   INTEGER, INTENT (IN) :: n
   INTEGER, INTENT (IN) :: MakeVec

   REAL (KIND=r8), INTENT (IN OUT) :: wr(n)
   REAL (KIND=r8), INTENT (IN OUT) :: wi(n)
   REAL (KIND=r8), INTENT (IN OUT) :: z(nm,n)

   INTEGER :: i, j, ic, j1

   REAL (KIND=r8) :: zz, div

   REAL (KIND=r8) :: h(nm,n)

   IF (MakeVec == 2) RETURN

   DO j=1,n
      DO i=1,n
         h(i,j)=z(i,j)
      END DO
   END DO

   DO j=1,n

      IF (wi(j) == 0.0_r8) THEN

         zz=0.0_r8
         DO i=1,n
            zz=MAX(zz,ABS(h(i,j)))
            IF (zz == ABS(h(i,j))) ic=i
         END DO

         DO i=1,n
            z(i,j)=h(i,j)/h(ic,j)
         END DO

      ELSE IF (wi(j) > 0.0_r8) THEN

         zz=0.0_r8
         j1=j+1
         DO i=1,n
            div=h(i,j)*h(i,j)+h(i,j1)*h(i,j1)
            zz=MAX(zz,div)
            IF (zz == div) ic=i
         END DO
         IF (ABS(h(ic,j)) < ABS(h(ic,j1))) THEN
            div=1.0_r8/h(ic,j1)
         ELSE
            div=1.0_r8/h(ic,j)
         END IF
         IF (div /= 0.0_r8) THEN
            DO i=1,n
               z(i,j)=h(i,j)*div
               z(i,j1)=h(i,j1)*div
            END DO
         END IF

      END IF

   END DO

   IF (MakeVec == 3) RETURN

   div=0.0_r8
   DO j=1,n
      zz=SQRT(wr(j)*wr(j)+wi(j)*wi(j))
      div=MAX(div,zz)
   END DO
   IF (div <= 0.0_r8) div=1.0_r8

   DO j=1,n
      IF (ABS(wr(j)/div) < TolW) wr(j)=0.0_r8
      IF (ABS(wi(j)/div) < TolW) wi(j)=0.0_r8
      DO i=1,n
         IF (ABS(z(i,j)) < TolZ) z(i,j)=0.0_r8
      END DO
   END DO

END SUBROUTINE zNorma


SUBROUTINE tql2 (nm, n, d, e, z, Eps, iErr)
 
   ! Computes the Eigenvalues and Eigenvectors of a Real
   ! Symmetric Tridiagonal Matrix Using the QL Method.

   ! Eps is a Machine Dependent Parameter Specifying the
   ! Relative Precision of the Floating Point Arithmetic.

   IMPLICIT NONE

   INTEGER, INTENT (IN) :: nm
   INTEGER, INTENT (IN) :: n

   REAL (KIND=r8), INTENT (IN) :: Eps

   REAL (KIND=r8), INTENT (IN OUT) :: e(n)
   REAL (KIND=r8), INTENT (IN OUT) :: d(n)
   REAL (KIND=r8), INTENT (IN OUT) :: z(nm,n)

   INTEGER, INTENT (OUT) :: iErr

   INTEGER :: i, j, k, l, m, ii, l1, mml

   REAL (KIND=r8) :: b, c, f, g, h, p, r, s

   INTEGER :: MaxIter=50

   iErr=0

   DO i=2,n
      e(i-1)=e(i)
   END DO

   f=0.0_r8
   b=0.0_r8
   e(n)=0.0_r8

   DO l=1,n
      j=0
      h=Eps*(ABS(d(l))+ABS(e(l)))
      IF (b < h) b=h

      ! Look for Small Sub-diagonal Element
      DO m=l,n
         IF (ABS(e(m)) <= b) EXIT
         ! e(n) is Always 0, so there is No Exit
         !      Through the Bottom of the Loop
      END DO

      120 IF (m == l) GO TO 220
      130 IF (j == MaxIter) THEN

         ! No Convergence to an Eigenvalue after MaxIter Iterations
         iErr=l
         WRITE (UNIT=nferr, FMT='(/,A,I5,A,/)') &
               ' *** The ', l, '-th Eigenvalue Did Not Converge ***'
         RETURN
      END IF

      j=j+1

      ! Form Shift

      l1=l+1
      g=d(l)
      p=(d(l1)-g)/(2.0_r8*e(l))
      r=SQRT(p*p+1.0_r8)
      d(l)=e(l)/(p+SIGN(r,p))
      h=g-d(l)

      DO i=l1,n
         d(i)=d(i)-h
      END DO

      f=f+h

      ! QL Transformation

      p=d(m)
      c=1.0_r8
      s=0.0_r8
      mml=m-l

      ! For i=m-1 Step -1 Until l DO --

      DO ii=1,mml
         i=m-ii
         g=c*e(i)
         h=c*p

         IF (ABS(p) < ABS(e(i))) THEN
            c=p/e(i)
            r=SQRT(c*c+1.0_r8)
            e(i+1)=s*e(i)*r
            s=1.0_r8/r
            c=c*s
         ELSE
            c=e(i)/p
            r=SQRT(c*c+1.0_r8)
            e(i+1)=s*p*r
            s=c/r
            c=1.0_r8/r
         END IF

         p=c*d(i)-s*g
         d(i+1)=h+s*(c*g+s*d(i))

         ! Form Vector

         DO k=1,n
            h=z(k,i+1)
            z(k,i+1)=s*z(k,i)+c*h
            z(k,i)=c*z(k,i)-s*h
         END DO

      END DO

      e(l)=s*p
      d(l)=c*p
      IF (ABS(e(l)) > b) GO TO 130
      220 d(l)=d(l)+f
   END DO

   ! Order Eigenvalues and Eigenvectors

   DO ii=2,n
      i=ii-1
      k=i
      p=d(i)

      DO j=ii,n
         IF (d(j) < p) THEN
            k=j
            p=d(j)
         END IF
      END DO

      IF (k /= i) THEN
         d(k)=d(i)
         d(i)=p

         DO j=1,n
            p=z(j,i)
            z(j,i)=z(j,k)
            z(j,k)=p
         END DO
      END IF

   END DO

END SUBROUTINE tql2


SUBROUTINE tred2 (nm, n, a, d, e, z)
 
   IMPLICIT NONE

   INTEGER, INTENT (IN) :: nm
   INTEGER, INTENT (IN) :: n

   REAL (KIND=r8), INTENT (IN) :: a(nm,n)

   REAL (KIND=r8), INTENT (OUT) :: d(n)
   REAL (KIND=r8), INTENT (OUT) :: e(n)
   REAL (KIND=r8), INTENT (OUT) :: z(nm,n)

   INTEGER :: i, j, k, l,  ii,  jp1

   REAL (KIND=r8) :: f, g, h, hh, Scal

   DO i=1,n
      DO j=1,i
         z(i,j)=a(i,j)
      END DO
   END DO

   IF (n == 1) GO TO 320

   ! For i=n Step -1 Until 2 DO

   DO ii=2,n
      i=n+2-ii
      l=i-1
      h=0.0_r8
      Scal=0.0_r8
      IF (l < 2) GO TO 130

      ! Scale Row
      DO k=1,l
         Scal=Scal+ABS(z(i,k))
      END DO

      IF (Scal /= 0.0_r8) GO TO 140

      130 e(i)=z(i,l)
      GO TO 290

      140 DO k=1,l
         z(i,k)=z(i,k)/Scal
         h=h+z(i,k)*z(i,k)
      END DO

      f=z(i,l)
      g=-SIGN(SQRT(h),f)
      e(i)=Scal*g
      h=h-f*g
      z(i,l)=f-g
      f=0.0_r8

      DO j=1,l
         z(j,i)=z(i,j)/h
         g=0.0_r8

         ! Form Element of a*u
         DO k=1,j
            g=g+z(j,k)*z(i,k)
         END DO

         jp1=j+1
         IF (l < jp1) GO TO 220

         DO k=jp1,l
            g=g+z(k,j)*z(i,k)
         END DO

         ! Form Element of p
         220 e(j)=g/h
         f=f+e(j)*z(i,j)

      END DO

      hh=f/(h+h)
      ! Form Reduced a
      DO j=1,l
         f=z(i,j)
         g=e(j)-hh*f
         e(j)=g

         DO k=1,j
            z(j,k)=z(j,k)-f*e(k)-g*z(i,k)
         END DO
      END DO

      290 d(i)=h
   END DO

   320 d(1)=0.0_r8
   e(1)=0.0_r8

   ! Accumulation of Transformation Matrices

   DO i=1,n
      l=i-1
      IF (d(i) == 0.0_r8) GO TO 380

      DO j=1,l
         g=0.0_r8

         DO k=1,l
            g=g+z(i,k)*z(k,j)
         END DO

         DO k=1,l
            z(k,j)=z(k,j)-g*z(k,i)
         END DO
      END DO

      380 d(i)=z(i,i)
      z(i,i)=1.0_r8
      IF (l < 1) CYCLE

      DO j=1,l
         z(i,j)=0.0_r8
         z(j,i)=0.0_r8
      END DO

   END DO

END SUBROUTINE tred2


END MODULE Eigen
