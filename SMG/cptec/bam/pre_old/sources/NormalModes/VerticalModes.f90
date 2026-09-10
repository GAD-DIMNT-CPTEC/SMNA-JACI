!
!  $Author: bonatti $
!  $Date: 2009/03/20 18:00:00 $
!  $Revision: 1.1.1.1 $
!
MODULE VerticalModes

USE InputParameters, ONLY: r8, Kmax, KmaxP, KmaxM, &
                           go, Rd, Cp, RE, Ps, SqRt2, Eps, &
                           nfprt, nfmod, MakeVec, PrtOut

USE Eigen , ONLY: GetEigenRGM

IMPLICIT NONE

PRIVATE

PUBLIC :: GetVerticalModes


CONTAINS


SUBROUTINE GetVerticalModes (gh, del)

   IMPLICIT NONE
 
   REAL (KIND=r8), INTENT (IN) :: del(Kmax)

   REAL (KIND=r8), INTENT (IN OUT) :: gh(Kmax)

   INTEGER :: k, i, j

   REAL (KIND=r8) :: p, siman, sum, dt, RE2

   REAL (KIND=r8) :: ci(KmaxP), si(KmaxP), sl(Kmax), cl(Kmax), &
                     rpi(KmaxM), sv(Kmax), p1(Kmax), p2(Kmax), &
                     h1(Kmax), h2(Kmax), dotpro(Kmax), tov(Kmax)

   REAL (KIND=r8) :: am(Kmax,Kmax), hm(Kmax,Kmax), tm(Kmax,Kmax), &
                     bm(Kmax,Kmax), cm(Kmax,Kmax), g(Kmax,Kmax), &
                     gt(Kmax,Kmax), eigg(Kmax,Kmax), eiggt(Kmax,Kmax)

   CALL SetSigma (ci, si, del, sl, cl, rpi)

   DO k=1,Kmax
      p=sl(k)*Ps
      CALL StandardAtm (p, tov(k))
   END DO

   CALL amhmtm (del, rpi, sv, p1, p2, am, hm, tm)

   dt=1.0_r8
   CALL bmcm (tov, p1, p2, h1, h2, del, ci, bm, cm, dt, sv, am)

   ! cm=g if a=1 and dt=1

   RE2=RE*RE
   DO i=1,Kmax
      DO j=1,Kmax
         cm(i,j)=Rd*tov(i)*sv(j)
         DO k=1,Kmax
            cm(i,j)=cm(i,j)+RE2*am(i,k)*bm(k,j)
         END DO
      END DO
   END DO

   IF (PrtOut) THEN
      WRITE (UNIT=nfprt, FMT='(/,A,/)') ' Matrix cm:'
      DO i=1,Kmax
         WRITE (UNIT=nfprt, FMT='(A,I6)') ' i = ', i
         WRITE (UNIT=nfprt, FMT='(1P6G12.5)') cm(i,:)
      END DO
   END IF

   DO i=1,Kmax
      DO j=1,Kmax
         g(i,j)=cm(i,j)
         gt(i,j)=cm(j,i)
      END DO
   END DO

   siman=-1.0_r8
   CALL VerticalEigen (g, siman, am, cl, eigg, gh)

   siman=1.0_r8
   CALL VerticalEigen (gt, siman, am, cl, eiggt, dotpro)

   ! dotpro=inverse dot prod. of eigenvec(g)*eigenvec(gtranspose)

   DO k=1,Kmax
      sum=0.0_r8
      DO j=1,Kmax
         sum=sum+eigg(j,k)*eiggt(j,k)
      END DO
      dotpro(k)=1.0_r8/sum
   END DO

   DO k=1,Kmax
      WRITE (UNIT=nfprt, FMT='(/,A,1P,2G12.5)') &
            ' g Eigenvalue and DotProd = ', gh(k), dotpro(k)
      WRITE (UNIT=nfprt, FMT='(/,A)') ' g Eigenvectors :'
      WRITE (UNIT=nfprt, FMT='(1X,1P,6G12.5)') (eigg(j,k),j=1,Kmax)
      WRITE (UNIT=nfprt, FMT='(/,A)') ' gTranpose Eigenvectors :'
      WRITE (UNIT=nfprt, FMT='(1X,1P,6G12.5)') (eiggt(j,k),j=1,Kmax)
   END DO

   WRITE (UNIT=nfmod) REAL(eigg,KIND=r8), REAL(eiggt,KIND=r8), &
                      REAL(gh,KIND=r8), REAL(dotpro,KIND=r8), REAL(tov,KIND=r8)

END SUBROUTINE GetVerticalModes


SUBROUTINE SetSigma (ci, si, del, sl, cl, rpi)
 
   ! Calculates Quantities Related to the 
   ! Discretization of the Sigma Coordinate

   ! ci(KmaxP)  Output : Sigma value at each level
   ! si(KmaxP)  Output : si(l)=1.0-ci(l)
   ! del(Kmax)  Output : Sigma spacing for each layer
   ! sl(Kmax)   Output : Sigma value at midpoint of
   !                     each layer (K=Rd/Cp):

   !                                                     1
   !                             +-                   + ---
   !                             !     K+1         K+1!  K
   !                             !si(l)   - si(l+1)   !
   !                     sl(l) = !--------------------!
   !                             !(K+1) (si(l)-si(l+1)!
   !                             +-                  -+

   ! cl(Kmax)   Output : cl(l)=1.0-sl(l)
   ! rpi(KmaxM) Output : Ratios of "pi" at adjacent layers:

   !                              +-     -+ K
   !                              !sl(l+1)!
   !                     rpi(l) = !-------!
   !                              ! sl(l) !
   !                              +-     -+

   IMPLICIT NONE

   REAL (KIND=r8), INTENT (IN) :: del(Kmax)

   REAL (KIND=r8), INTENT (OUT) :: ci(KmaxP)
   REAL (KIND=r8), INTENT (OUT) :: si(KmaxP)
   REAL (KIND=r8), INTENT (OUT) :: sl(Kmax)
   REAL (KIND=r8), INTENT (OUT) :: cl(Kmax)
   REAL (KIND=r8), INTENT (OUT) :: rpi(KmaxM)

   INTEGER :: k

   REAL (KIND=r8) :: Kp, Kp1, dif, siKp, siKp1

   WRITE (UNIT=nfprt, FMT='(/,3(A,I4),/)') &
         ' Begin SetSigma  Kmax = ', Kmax, &
         '  KmaxP = ', KmaxP, '  KmaxM = ', KmaxM

!cdir novector
   ci(1)=0.0_r8
   DO k=1,Kmax
      ci(k+1)=ci(k)+del(k)
   END DO
   ci(KmaxP)=1.0_r8

   DO k=1,KmaxP
      si(k)=1.0_r8-ci(k)
   END DO

   Kp=Rd/Cp
   Kp1=Kp+1.0_r8
   DO k=1,Kmax
      ! dif=si(k)**Kp1-si(k+1)**Kp1
      siKp=EXP(Kp1*LOG(si(k)))
      IF (k <= Kmax-1) THEN
         siKp1=EXP(Kp1*LOG(si(k+1)))
      ELSE
         siKp1=0.0_r8
      END IF
      dif=(siKp-siKp1)/(Kp1*(si(k)-si(k+1)))
      ! sl(k)=dif**(1.0_r8/Kp)
      sl(k)=EXP(LOG(dif)/Kp)
      cl(k)=1.0_r8-sl(k)
   END DO

   ! Compute pi Ratios for Temperature Matrix
   DO k=1,KmaxM
      ! rpi(k)=(sl(k+1)/sl(k))**Kp
      rpi(k)=EXP(Kp*LOG(sl(k+1)/sl(k)))
   END DO

   DO k=1,KmaxP
     WRITE (UNIT=nfprt, FMT='(A,I4,2(A,F10.5))') &
           ' Interface = ', k, '  ci = ', ci(k), '  si = ', si(k)
   END DO

   WRITE (UNIT=nfprt, FMT='(A)') ' '
   DO k=1,Kmax
     WRITE (UNIT=nfprt, FMT='(A,I4,3(A,F10.5))') &
           '     Layer = ', k, '  cl = ', cl(k), &
           '  sl = ', sl(k), '  del = ', del(k)
   END DO

   WRITE (UNIT=nfprt, FMT='(/,A)') ' rpi: '
   WRITE (UNIT=nfprt, FMT='(6(F10.5,2X))') rpi
   WRITE (UNIT=nfprt, FMT='(/,A,/)') ' End SetSigma'

END SUBROUTINE SetSigma


SUBROUTINE StandardAtm (Press, Temp)

   ! Computes the standard Height (meters), Temperature (K) and
   ! potential Temperature (deg k) given the Pressure in hPa
   ! U. S. Standard Atmosphere, 1962
   ! ICAO SStandard Atmosphere up to 20KM
   ! Proposed Extension up to 32km
   ! Not Valid for Height > 32km or Pressure < 8.68 hPa
 
   IMPLICIT NONE

   REAL (KIND=r8), INTENT (IN) :: Press

   REAL (KIND=r8), INTENT (OUT) :: Temp

   REAL (KIND=r8) :: Height, Theta, RdbyCp, Rdbygo, fkt, ar, pp0

   REAL (KIND=r8) :: Piso=54.7487_r8
   REAL (KIND=r8) :: Ziso=20000.0_r8
   REAL (KIND=r8) :: Salp=-0.0010_r8
   REAL (KIND=r8) :: Pzero=1013.25_r8
   REAL (KIND=r8) :: T0=288.15_r8
   REAL (KIND=r8) :: Alp=0.0065_r8
   REAL (KIND=r8) :: Ptrop=226.321_r8
   REAL (KIND=r8) :: Tstr=216.65_r8
   REAL (KIND=r8) :: Htrop=11000.0_r8

   RdbyCp=Rd/Cp
   Rdbygo=Rd/go

   IF (Press < Piso) THEN
      ! Compute Lapse Rate=-0.0010 Cases
      ar=Salp*Rdbygo
      pp0=Piso**ar
      Height=(Tstr/(pp0*Salp))*(pp0-(Press**ar))+Ziso
      Temp=Tstr-(Height-Ziso)*Salp

   ELSE IF (Press > Ptrop) THEN
      ! Compute Lapse Rate=0.0065 Cases
      ar=Alp*Rdbygo
      pp0=Pzero**ar
      Height=(T0/(pp0*Alp))*(pp0-(Press**ar))
      Temp=T0-Height*Alp

   ELSE
      ! Compute Isothermal Cases
      fkt=Rdbygo*Tstr
      Height=Htrop+fkt*LOG(Ptrop/Press)
      Temp=Tstr

   END IF

   Theta=Temp*((Ps/Press)**RdbyCp)

END SUBROUTINE StandardAtm


SUBROUTINE amhmtm (del, rpi, sv, p1, p2, am, hm, tm)
 
   ! Calculates the matrix associated with the 
   ! Hydrostatic Equation according to the Arakawa 
   ! vertical finite differencing scheme.
   ! Also calculates the "pi" ratios.

   ! del(Kmax)     Input : Sigma spacing for each layer
   ! rpi(KmaxM)    Input : Ratios of "pi" at adjacent layers

   ! sv(Kmax)      Output : sv(l)=-del(l)
   ! p1(Kmax)      Output : p1(l)=1/rpi(l); p1(Kmax) = 0
   ! p2(Kmax)      Output : p2(l+1)=rpi(l); p2(1) = 0
   ! am(Kmax,Kmax) Output : "hm" matrix divided by square
   !                        of earth radius for laplacian
   ! hm(Kmax,Kmax) Output : Matrix relating geopotential
   !                        to temperature
   ! tm(Kmax,Kmax) Output : Inverse of "hm" matrix

   IMPLICIT NONE

   REAL (KIND=r8), INTENT (IN) :: del(Kmax)
   REAL (KIND=r8), INTENT (IN) :: rpi(KmaxM)

   REAL (KIND=r8), INTENT (OUT) :: sv(Kmax)
   REAL (KIND=r8), INTENT (OUT) :: p1(Kmax)
   REAL (KIND=r8), INTENT (OUT) :: p2(Kmax)
   REAL (KIND=r8), INTENT (OUT) :: am(Kmax,Kmax)
   REAL (KIND=r8), INTENT (OUT) :: hm(Kmax,Kmax)
   REAL (KIND=r8), INTENT (OUT) :: tm(Kmax,Kmax)

   INTEGER :: k, i, j

   REAL (KIND=r8) :: RESqInv

   RESqInv=1.0_r8/(RE*RE)

   DO k=1,Kmax
    DO i=1,kMax
      hm(k,i)=0.0_r8
      tm(k,i)=0.0_r8
    end do
   END DO

!cdir novector
   DO k=1,KmaxM
      hm(k,k)=1.0_r8
      tm(k,k)=0.5_r8*Cp*(rpi(k)-1.0_r8)
   END DO

   DO k=1,KmaxM
      hm(k,k+1)=-1.0_r8
      tm(k,k+1)=0.5_r8*Cp*(1.0_r8-1.0_r8/rpi(k))
   END DO
   DO k=1,Kmax
      hm(Kmax,k)=del(k)
      tm(Kmax,k)=Rd*del(k)
   END DO

   CALL InvertMatrix (hm, Kmax)

   DO i=1,Kmax
      DO j=1,Kmax
         am(i,j)=0.0_r8
         DO k=1,Kmax
            am(i,j)=am(i,j)+hm(i,k)*tm(k,j)
         END DO
      END DO
   END DO

   ! Divide by a**2 for Laplacian and 
   ! store am in tm and divide am

!cdir vector
   DO k=1,Kmax
     DO i=1,kMax
      tm(k,i)=am(k,i)
      hm(k,i)=am(k,i)
      am(k,i)=am(k,i)*RESqInv
     enddo
   END DO

   CALL InvertMatrix (tm, Kmax)

!cdir novector
   DO k=1,Kmax
      sv(k)=-del(k)
   END DO
   DO k=1,KmaxM
      p1(k)=1.0_r8/rpi(k)
      p2(k+1)=rpi(k)
   END DO
   p1(Kmax)=0.0_r8
   p2(1)=0.0_r8

END SUBROUTINE amhmtm


SUBROUTINE bmcm (tov, p1, p2, h1, h2, del, ci, &
                 bm, cm, dt, sv, am)
 
   ! Calculates arrays used in the semi-implicit integration
   ! scheme due to terms arising from the thermodynamic and
   ! divergence equations.

   ! dt            Input : time step
   ! tov(kmax)     Input : Temperature of Standard Atmosphere (K)
   ! p1(kmax)      Input : Reciprocals of "pi" ratios
   ! p2(kmax)      Input : "pi" ratios
   ! del(kmax)     Input : Sigma spacing for each layer
   ! ci(kmaxp)     Input : Sigma value for each level
   ! sv(kmax)      Input : Negative of "del"
   ! am(kmax,kmax) Input : Matrix relating geopotential

   ! h1(kmax)      Output : h1(l)=p1(l)*tov(l+1)-tov(l); h1(kmax)=0
   ! h2(kmax)      Output : h2(l)=tov(l)-p2(l)*tov(l-1); h2(1)=0
   ! bm(kmax,kmax) Output : Matrix relating temperature
   !                        tendency to the divergence
   ! cm(kmax,kmax) Output : Matrix resulting from product of
   !                        "am" and "bm" and further addition
   !                        of "tov" terms to temperature

   IMPLICIT NONE

   REAL (KIND=r8), INTENT (IN) :: dt

   REAL (KIND=r8), INTENT (IN) :: tov(kmax)
   REAL (KIND=r8), INTENT (IN) :: p1(kmax )
   REAL (KIND=r8), INTENT (IN) :: p2(kmax)
   REAL (KIND=r8), INTENT (IN) :: del(kmax)
   REAL (KIND=r8), INTENT (IN) :: ci(kmaxp)
   REAL (KIND=r8), INTENT (IN) :: sv(kmax)
   REAL (KIND=r8), INTENT (IN) :: am(kmax,kmax)

   REAL (KIND=r8), INTENT (OUT) :: h1(kmax)
   REAL (KIND=r8), INTENT (OUT) :: h2(kmax)
   REAL (KIND=r8), INTENT (OUT) :: bm(kmax,kmax)
   REAL (KIND=r8), INTENT (OUT) :: cm(kmax,kmax)

   INTEGER :: k, i, j

   REAL (KIND=r8) :: Kp, RdByaa

   REAL (KIND=r8) :: x1(kmax), x2(kmax)

   Kp=Rd/Cp

!cdir novector
   DO k=1,kmax-1
      h1(k)=p1(k)*tov(k+1)-tov(k)
   END DO
   h1(kmax)=0.0_r8
   h2(1)=0.0_r8
   DO k=2,kmax
      h2(k)=tov(k)-p2(k)*tov(k-1)
   END DO

   DO k=1,kmax
      x1(k)=Kp*tov(k)+0.5_r8*(ci(k+1)*h1(k)+ci(k)*h2(k))/del(k)
      x2(k)=0.5_r8*(h1(k)+h2(k))/del(k)
   END DO

   DO j=1,kmax
      DO k=1,kmax
         bm(k,j)=-x1(k)*del(j)
      END DO
   END DO
   DO k=1,kmax
      DO j=1,k
         bm(k,j)=bm(k,j)+x2(k)*del(j)
      END DO
   END DO
   DO k=1, kmax
      bm(k,k)=bm(k,k)-0.5_r8*h2(k)
   END DO

   RdByaa=Rd/(RE*RE)
   DO i=1,kmax
      DO j=1,kmax
         cm(i,j)=0.0_r8
         DO k=1,kmax
            cm(i,j)=cm(i,j)+am(i,k)*bm(k,j)
         END DO
         cm(i,j)=(cm(i,j)+RdByaa*tov(i)*sv(j))*dt*dt
      END DO
   END DO

END SUBROUTINE bmcm


SUBROUTINE VerticalEigen (gg, siman, eigvc, col, vec, val)

   IMPLICIT NONE
 
   REAL (KIND=r8), INTENT (IN) :: siman

   REAL (KIND=r8), INTENT (IN OUT) :: gg(kmax,kmax)

   REAL (KIND=r8), INTENT (OUT) :: eigvc(kmax,kmax)
   REAL (KIND=r8), INTENT (OUT) :: col(kmax)
   REAL (KIND=r8), INTENT (OUT) :: vec(kmax,kmax)
   REAL (KIND=r8), INTENT (OUT) :: val(kmax)

   INTEGER :: iErr, i, j, k, kkk

   REAL (KIND=r8) :: sum, rmax

   INTEGER :: kk(kmax)

   REAL (KIND=r8) :: eigvr(kmax), eigvi(kmax)

   CALL GetEigenRGM (kmax, kmax, gg, eigvr, eigvi, MakeVec, eigvc, iErr, Eps)

   WRITE (UNIT=nfprt, FMT='(/,A,I4,A)') &
         ' iErr = ', iErr, ' Eigenvalues and Vectors Follow'
   DO k=1,kmax
      WRITE (UNIT=nfprt, FMT='(/,1P,2G12.5)') eigvr(k), eigvi(k)
      WRITE (UNIT=nfprt, FMT='(1X,1P,6G12.5)') (eigvc(i,k),i=1,kmax)
   END DO
   WRITE (UNIT=nfprt, FMT='(/,A,/)') ' End Eigenvalues and vectors'

   DO k=1,kmax
      kk(k)=0
      col(k)=siman*eigvr(k)
      sum=0.0_r8
      DO j=1,kmax
         sum=sum+eigvc(j,k)*eigvc(j,k)
      END DO
      ! sum=length of eigenvector k
      sum=1.0_r8/SQRT(sum)
      DO j=1,kmax
         eigvc(j,k)=sum*eigvc(j,k)
      END DO
   END DO

   ! Eigenvalues now have unit length
   ! k-th vector is eigvc(j,k)
   ! Eigenvalues are now in col(k)
   ! Next arrange in descending order

   WRITE (UNIT=nfprt, FMT='(1X,1P,6G12.5)') (col(k),k=1,kmax)
   DO j=1,kmax
      IF (col(j) == 0.0_r8) WRITE (UNIT=nfprt, FMT='(/,A,/)') &
                                  ' Zero Eigenvalue'
   END DO

   WRITE (UNIT=nfprt, FMT='(A)') ' '
   WRITE (UNIT=nfprt, FMT='(1X,20I4)') (kk(j),j=1,kmax)
   DO j=1,kmax
      rmax=-1.0E20_r8
      kkk=0
      DO k=1,kmax
         IF (ABS(col(k)) > rmax) kkk=k
         IF (ABS(col(k)) > rmax) rmax=ABS(col(k))
      END DO
      val(j)=col(kkk)
      col(kkk)=0.0_r8
      DO i=1,kmax
         vec(i,j)=eigvc(i,kkk)
      END DO
      kk(j)=kkk
   END DO
   WRITE (UNIT=nfprt, FMT='(A)') ' '
   WRITE (UNIT=nfprt, FMT='(1X,20I4)') (kk(j),j=1,kmax)
   WRITE (UNIT=nfprt, FMT='(1X,1P,6G12.5)') (val(k),k=1,kmax)

END SUBROUTINE VerticalEigen


SUBROUTINE InvertMatrix (a, n)
 
   ! Inverts a Matrix in Place using Gauss-Jordan Technique.
   ! The Determinant is also calculated.
   ! Input Matrix is destroyed.

   IMPLICIT NONE

   INTEGER, INTENT (IN) :: n ! Order of Matrix a

   REAL (KIND=r8), INTENT (IN OUT) :: a(n*n) ! Input: Matrix to be Inverted
                                             ! Output: Inverse of Input Matrix

   REAL (KIND=r8) :: d ! Determinant of Input Matrix
   REAL (KIND=r8) :: biga, hold

   INTEGER :: nk, k, kk, j, iz, i, ij, ki, ji, jp, jk, ik, kj, jq, jr

   INTEGER :: l(n), m(n)

   ! Search for Largest Element

   d=1.0_r8
   nk=-n

!cdir novector
   DO k=1,n

      nk=nk+n
      l(k)=k
      m(k)=k
      kk=nk+k
      biga=a(kk)

      DO j=k,n
         iz=n*(j-1)
         DO i=k,n
            ij=iz+i
            IF (ABS(biga)- ABS(a(ij)) < 0.0_r8) THEN
               biga=a(ij)
               l(k)=i
               m(k)=j
            END IF
         END DO
      END DO
 
      ! Interchange Rows
  
      j=l(k)
      IF (j > k) THEN
         ki=k-n
!cdir nodep
         DO i=1,n
            ki=ki+n
            hold=-a(ki)
            ji=ki-k+j
            a(ki)=a(ji)
            a(ji)=hold
         END DO
      END IF
  
      ! Interchange Columns
  
      i=m(k)
      IF (i > k) THEN
         jp=n*(i-1)
!cdir nodep
         DO j=1,n
            jk=nk+j
            ji=jp+j
            hold=-a(jk)
            a(jk)=a(ji)
            a(ji) =hold
         END DO
      END IF
  
      ! Divide Column by Minus pivot 
      ! Value of pivot Element is Contained in biga)
  
      IF (biga == 0.0_r8) THEN
         d=0.0_r8
         RETURN
      END IF

      DO i=1,n
         IF (i /= k) THEN
            ik=nk+i
            a(ik)=a(ik)/(-biga)
         END IF
      END DO
  
      ! Reduce Matrix
  
      DO i=1,n
         ik=nk+i
         ij=i-n
!cdir nodep
         DO j=1,n
            ij=ij+n
            IF (i == k) CYCLE
            IF (j == k) CYCLE
            kj=ij-i+k
            a(ij)=a(ik)*a(kj)+a(ij)
         END DO
      END DO
  
      ! Divide Row by pivot
  
      kj=k-n
      DO j=1,n
         kj=kj+n
         IF (j == k) CYCLE
         a(kj)=a(kj)/biga
      END DO
  
      ! Product of pivots
  
      d=d*biga
  
      ! Replace pivot by Reciprocal
  
     a(kk)=1.0_r8/biga

   END DO

   ! Final Row and Column Interchange

   k=n
   DO
      k=k-1
      IF (k <= 0) RETURN

      i=l(k)
      IF (i > k) THEN
         jq=n*(k-1)
         jr=n*(i-1)
!cdir nodep
         DO j=1,n
            jk=jq+j
            hold=a(jk)
            ji=jr+j
            a(jk)=-a(ji)
            a(ji)=hold
         END DO
      END IF

      j=m(k)
      IF (j > k) THEN
         ki=k-n
!cdir nodep
         DO i=1,n
            ki=ki+n
            hold=a(ki)
            ji=ki-k+j
            a(ki)=-a(ji)
            a(ji)=hold
         END DO
      END IF
   END DO

END SUBROUTINE InvertMatrix


END MODULE VerticalModes
