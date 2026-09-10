!
!  $Author: bonatti $
!  $Date: 2009/03/20 18:00:00 $
!  $Revision: 1.1.1.1 $
!
MODULE HorizontalModes

USE InputParameters, ONLY: r8, i8, Mend, Mend1, NxAsy, NxSym, &
                           nfprt, nfmod, PrtOut, &
                           go, RE, twoOmega, PerCut, SqRt2, Eps

USE Eigen, ONLY: tql2, tred2

IMPLICIT NONE

REAL (KIND=r8) :: Tolx

PRIVATE

PUBLIC :: GetHorizontalModes


CONTAINS


SUBROUTINE GetHorizontalModes (Mods, gh)

   IMPLICIT NONE
 
   ! Triangular Truncation Horizontal Modes

   INTEGER, INTENT (IN) :: Mods

   REAL (KIND=r8), INTENT (IN) :: gh(Mods)

   INTEGER :: k, m, nlx, nmd, Nmax, Lmax, Mmax, KLmx, &
              nSym, nAsy, n, nCutSym, nCutAsy

   REAL (KIND=r8) :: REi, rm, rn

   REAL (KIND=r8) :: Alfa(Mend1), Beta(Mend1), Gamma(Mend1), &
                     wSym(NxSym), wAsy(NxAsy)
   REAL (KIND=r8) :: xSym(NxSym,NxSym), xAsy(NxAsy,NxAsy)

   !Tolx=1.0E-09_r8
   Tolx=1.0E-12_r8
   REi=1.0_r8/RE

   DO k=1,Mods

      WRITE (UNIT=nfprt, FMT='(/,2(A6,I5))') &
            ' Mend = ', Mend, ' Mode = ', k
      WRITE (UNIT=nfprt, FMT='(/,1P,3(A7,G12.5))') &
            ' 2*Omega = ', twoOmega, '  Hn = ', gh(k)/go, ' PerCut = ', PerCut

      DO m=1,Mend1

         Nmax=Mend1-m+1
         nlx=Mend1-m
         nmd=MOD(nlx,2)
         Lmax=(nlx+nmd)/2
         Mmax=Lmax+1-nmd
         KLmx=Lmax+Mmax
         nSym=Lmax+2*Mmax
         nAsy=Mmax+2*Lmax
         WRITE (UNIT=nfprt, FMT='(/,4(A,I5))') &
               ' m = ', m-1,' Nmax = ', Nmax-1, &
               ' Lmax = ', Lmax, ' Mmax = ', Mmax

         rm=REAL(m-1,r8)
         DO n=1,Nmax
            rn=rm+REAL(n-1,r8)
            IF (rn == 0.0_r8) THEN
               Alfa(n)=0.0_r8
               Beta(n)=0.0_r8
               Gamma(n)=0.0_r8
            ELSE
               Alfa(n)=twoOmega*rm/(rn*(rn+1.0_r8))
               Beta(n)=(twoOmega/rn)*SQRT((rn*rn-1.0_r8)* &
                       (rn*rn-rm*rm)/(4.0_r8*rn*rn-1.0_r8))
               Gamma(n)=REi*SQRT(rn*(rn+1.0_r8)*gh(k))
            END IF
         END DO

         IF (PrtOut) THEN
            WRITE (UNIT=nfprt, FMT='(/,A)') ' Alfa :'
            WRITE (UNIT=nfprt, FMT='(1P,6G12.5)') (Alfa(n),n=1,Nmax)
            WRITE (UNIT=nfprt, FMT='(A)')' Beta :'
            WRITE (UNIT=nfprt, FMT='(1P,6G12.5)') (Beta(n),n=1,Nmax)
            WRITE (UNIT=nfprt, FMT='(A)') ' Gamma :'
            WRITE (UNIT=nfprt, FMT='(1P,6G12.5)') (Gamma(n),n=1,Nmax)
         END IF

         IF (m == 1) THEN

            ! Symmetric Case: s = 0
            CALL SymSeq0 (nSym, Lmax, KLmx, nmd, nCutSym,  &
                          Beta, Gamma, wSym, xSym)
            CALL Record (NxSym, nCutSym, nSym, wSym, xSym)

            ! Asymmetric Case: s = 0
            CALL AsySeq0 (nAsy, Lmax, KLmx, Mmax, nmd, nCutAsy, &
                          Beta, Gamma, wAsy, xAsy)
            CALL Record (NxAsy, nCutAsy, nAsy, wAsy, xAsy)

         ELSE

            ! Symmetric Case: s /= 0
            CALL SymSne0 (nSym, Lmax, Mmax, nmd, nCutSym, &
                          Alfa, Beta, Gamma, wSym, xSym)
            CALL Record (NxSym, nCutSym, nSym, wSym, xSym)

            ! Asymmetric Case: s /= 0
            CALL AsySne0 (nAsy, Lmax, Mmax, nmd, nCutAsy, &
                          Alfa, Beta, Gamma, wAsy, xAsy)
            CALL Record (NxAsy, nCutAsy, nAsy, wAsy, xAsy)

         END IF

      END DO
   END DO

END SUBROUTINE GetHorizontalModes


SUBROUTINE SymSeq0 (nSym, Lmax, KLmx, nmd, nCutSym, &
                    Beta, Gamma, wSym, xSym)

   IMPLICIT NONE

   INTEGER, INTENT (IN) :: nSym
   INTEGER, INTENT (IN) :: Lmax
   INTEGER, INTENT (IN) :: KLmx
   INTEGER, INTENT (IN) :: nmd

   REAL (KIND=r8), INTENT (IN) :: Beta(Mend1)
   REAL (KIND=r8), INTENT (IN) :: Gamma(Mend1)

   INTEGER, INTENT (OUT) :: nCutSym

   REAL (KIND=r8), INTENT (OUT) :: wSym(NxSym)
   REAL (KIND=r8), INTENT (OUT) :: xSym(NxSym,NxSym)

   INTEGER :: n, nn, j, jj, iErr, Lmax1

   REAL (KIND=r8) :: SubDiagSym(NxSym), DiagSym(NxSym)
   REAL (KIND=r8) :: ySym(Mend1,Lmax)

   ! Symmetric Case: s = 0

   Lmax1=1
   n=1
   SubDiagSym(n)=0.0_r8
   nn=2*n+1
   DiagSym(n)=Beta(nn)*Beta(nn)+Beta(nn+1)*Beta(nn+1)+ Gamma(nn)*Gamma(nn)
   Lmax1=Lmax-1
   DO n=2,Lmax1
      nn=2*n+1
      SubDiagSym(n)=Beta(nn-1)*Beta(nn)
      DiagSym(n)=Beta(nn)*Beta(nn)+Beta(nn+1)*Beta(nn+1)+ Gamma(nn)*Gamma(nn)
   END DO
   IF (nmd == 0) THEN
      Lmax1=Lmax
      n=Lmax
      nn=2*n+1
      SubDiagSym(n)=Beta(nn-1)*Beta(nn)
      DiagSym(n)=Beta(nn)*Beta(nn)+Gamma(nn)*Gamma(nn)
   END IF

   IF (PrtOut) THEN
      WRITE (UNIT=nfprt, FMT='(/,A)') ' SubDiagSym :'
      WRITE (UNIT=nfprt, FMT='(1P,6G12.5)') (SubDiagSym(n),n=1,Lmax1)
      WRITE (UNIT=nfprt, FMT='(A)') ' dga :'
      WRITE (UNIT=nfprt, FMT='(1P,6G12.5)') (DiagSym(n),n=1,Lmax1)
   END IF

   ySym=0.0_r8
   DO j=1,Lmax1
      ySym(j,j)=SqRt2
   END DO
   CALL tql2 (Mend1, Lmax1, DiagSym, SubDiagSym, ySym, Eps, iErr)

   DO j=1,Lmax1
      jj=2*j
      wSym(jj-1)=-SQRT(DiagSym(j))
      wSym(jj)=SQRT(DiagSym(j))
      xSym(1,jj-1)=Beta(3)*ySym(1,j)/wSym(jj-1)
      xSym(1,jj)=-xSym(1,jj-1)
      xSym(Lmax+1,jj-1)=0.0_r8
      xSym(Lmax+1,jj)=xSym(Lmax+1,jj-1)
      xSym(KLmx+1,jj-1)=0.0_r8
      xSym(KLmx+1,jj)=-xSym(KLmx+1,jj-1)
      DO n=2,Lmax1
         nn=2*n+1
         xSym(n,jj-1)=(Beta(nn-1)*ySym(n-1,j)+Beta(nn)*ySym(n,j))/wSym(jj-1)
         xSym(n,jj)=-xSym(n,jj-1)
         xSym(Lmax+n,jj-1)=ySym(n-1,j)
         xSym(Lmax+n,jj)=xSym(Lmax+n,jj-1)
         xSym(KLmx+n,jj-1)=Gamma(nn-2)*ySym(n-1,j)/wSym(jj-1)
         xSym(KLmx+n,jj)=-xSym(KLmx+n,jj-1)
      END DO
      n=Lmax
      nn=2*n-1
      IF (nmd == 1) THEN
         xSym(n,jj-1)=Beta(nn+1)*ySym(n-1,j)/wSym(jj-1)
         xSym(n,jj)=-xSym(n,jj-1)
         xSym(KLmx,jj-1)=ySym(n-1,j)
         xSym(KLmx,jj)=xSym(KLmx,jj-1)
         xSym(nSym,jj-1)=Gamma(nn)*ySym(n-1,j)/wSym(jj-1)
         xSym(nSym,jj)=-xSym(nSym,jj-1)
      ELSE
         xSym(KLmx,jj-1)=ySym(n,j)
         xSym(KLmx,jj)=xSym(KLmx,jj-1)
         xSym(nSym,jj-1)=Gamma(nn+2)*ySym(n,j)/wSym(jj-1)
         xSym(nSym,jj)=-xSym(nSym,jj-1)
      END IF
   END DO

   nCutSym=2*Lmax1
   WHERE (ABS(xSym) <= Tolx) xSym=0.0_r8

   WRITE (UNIT=nfprt, FMT='(/,2(A,I4))') &
         ' Frequency: nSym = ', nSym, ' nCutSym = ', nCutSym
   WRITE (UNIT=nfprt, FMT='(1P,6G12.5)') (wSym(n)/twoOmega,n=1,nCutSym)
   WRITE (UNIT=nfprt, FMT='(A)')' Period :'
   WRITE (UNIT=nfprt, FMT='(1P,6G12.5)') (1.0_r8/wSym(n),n=1,nCutSym)

   IF (PrtOut) THEN
      WRITE (UNIT=nfprt, FMT='(/,A,I4)') ' xSym : iErr = ', iErr
      DO n=1,Lmax1
         WRITE (UNIT=nfprt, FMT='(1P,6G12.5)') (ySym(n,nn),nn=1,Lmax1)
      END DO
      WRITE (UNIT=nfprt, FMT='(/,A)') ' xSym :'
      DO n=1,nSym
         WRITE (UNIT=nfprt, FMT='(1P,6G12.5)') (xSym(n,nn),nn=1,MIN(6,nCutSym))
      END DO
   END IF

END SUBROUTINE SymSeq0


SUBROUTINE AsySeq0 (nAsy, Lmax, KLmx, Mmax, nmd, nCutAsy, &
                    Beta, Gamma, wAsy, xAsy)
 
   IMPLICIT NONE

   INTEGER,  INTENT (IN) :: nAsy
   INTEGER,  INTENT (IN) :: Lmax
   INTEGER,  INTENT (IN) :: KLmx
   INTEGER,  INTENT (IN) :: Mmax
   INTEGER,  INTENT (IN) :: nmd

   REAL (KIND=r8),  INTENT (IN) :: Beta(Mend1)
   REAL (KIND=r8),  INTENT (IN) :: Gamma(Mend1)

   INTEGER,  INTENT (OUT) :: nCutAsy

   REAL (KIND=r8),  INTENT (OUT) :: wAsy(NxAsy)
   REAL (KIND=r8),  INTENT (OUT) :: xAsy(NxAsy,NxAsy)

   INTEGER :: n, nn, j, jj, iErr, Lmax1

   REAL (KIND=r8) :: SubDiagAsy(NxAsy), DiagAsy(NxAsy)
   REAL (KIND=r8) :: yAsy(Mend1,Lmax)

   ! Asymmetric Case: s = 0

   Lmax1=1
   n=1
   SubDiagAsy(n)=0.0_r8
   nn=2*n
   SubDiagAsy(n)=Beta(nn-1)*Beta(nn)
   DiagAsy(n)=Beta(nn)*Beta(nn)+Beta(nn+1)*Beta(nn+1)+ Gamma(nn)*Gamma(nn)
   Lmax1=Lmax-1
   DO n=2,Lmax1
      nn=2*n
      SubDiagAsy(n)=Beta(nn-1)*Beta(nn)
      DiagAsy(n)=Beta(nn)*Beta(nn)+Beta(nn+1)*Beta(nn+1)+ Gamma(nn)*Gamma(nn)
   END DO
   Lmax1=Lmax
   n=Lmax
   nn=2*n
   SubDiagAsy(n)=Beta(nn-1)*Beta(nn)
   IF (nmd == 0) THEN
      DiagAsy(n)=Beta(nn)*Beta(nn)+Beta(nn+1)*Beta(nn+1)+ Gamma(nn)*Gamma(nn)
   ELSE
      DiagAsy(n)=Beta(nn)*Beta(nn)+Gamma(nn)*Gamma(nn)
   END IF

   IF (PrtOut) THEN
      WRITE (UNIT=nfprt, FMT='(/,A)') ' SubDiagAsy :'
      WRITE (UNIT=nfprt, FMT='(1P,6G12.5)') (SubDiagAsy(n),n=1,Lmax1)
      WRITE (UNIT=nfprt, FMT='(A)') ' dga :'
      WRITE (UNIT=nfprt, FMT='(1P,6G12.5)') (DiagAsy(n),n=1,Lmax1)
   END IF

   yAsy=0.0_r8
   DO j=1,Lmax1
      yAsy(j,j)=SqRt2
   END DO
   CALL tql2 (Mend1, Lmax1, DiagAsy, SubDiagAsy, yAsy, Eps, iErr)

   DO j=1,Lmax1
      jj=2*j
      wAsy(jj-1)=-SQRT(DiagAsy(j))
      wAsy(jj)=SQRT(DiagAsy(j))
      xAsy(1,jj-1)=0.0_r8
      xAsy(1,jj)=-xAsy(1,jj-1)
      xAsy(Mmax+1,jj-1)=yAsy(1,j)
      xAsy(Mmax+1,jj)=xAsy(Mmax+1,jj-1)
      xAsy(KLmx+1,jj-1)=Gamma(2)*yAsy(1,j)/wAsy(jj-1)
      xAsy(KLmx+1,jj)=-xAsy(KLmx+1,jj-1)
      DO n=2,Lmax1
         nn=2*n
         xAsy(n,jj-1)=(Beta(nn-1)*yAsy(n-1,j)+Beta(nn)*yAsy(n,j))/wAsy(jj-1)
         xAsy(n,jj)=-xAsy(n,jj-1)
         xAsy(Mmax+n,jj-1)=yAsy(n,j)
         xAsy(Mmax+n,jj)=xAsy(Mmax+n,jj-1)
         xAsy(KLmx+n,jj-1)=Gamma(nn)*yAsy(n,j)/wAsy(jj-1)
         xAsy(KLmx+n,jj)=-xAsy(KLmx+n,jj-1)
      END DO
      n=Lmax1
      nn=2*n+1
      IF (nmd == 0) THEN
         xAsy(Mmax,jj-1)=Beta(nn)*yAsy(n,j)/wAsy(jj-1)
         xAsy(Mmax,jj)=-xAsy(Mmax,jj-1)
      END IF
   END DO

   nCutAsy=2*Lmax1
   WHERE (ABS(xAsy) <= Tolx) xAsy=0.0_r8

   WRITE (UNIT=nfprt, FMT='(/,2(A,I4))') &
         ' Frequency: nAsy = ', nAsy, ' nCutAsy = ', nCutAsy
   WRITE (UNIT=nfprt, FMT='(1P,6G12.5)') (wAsy(n)/twoOmega,n=1,nCutAsy)
   WRITE (UNIT=nfprt, FMT='(A)')' Period :'
   WRITE (UNIT=nfprt, FMT='(1P,6G12.5)') (1.0_r8/wAsy(n),n=1,nCutAsy)

   IF (PrtOut) THEN
      WRITE (UNIT=nfprt, FMT='(/,A,I4)') ' xAsy : iErr = ', iErr
      DO n=1,Lmax1
         WRITE (UNIT=nfprt, FMT='(1P,6G12.5)') (yAsy(n,nn),nn=1,Lmax1)
      END DO
      WRITE (UNIT=nfprt, FMT='(/,A)') ' xAsy :'
      DO n=1,nAsy
         WRITE (UNIT=nfprt, FMT='(1P,6G12.5)') (xAsy(n,nn),nn=1,MIN(6,nCutAsy))
      END DO
   END IF

END SUBROUTINE AsySeq0


SUBROUTINE SymSne0 (nSym, Lmax, Mmax, nmd, nCutSym,  &
                    Alfa, Beta, Gamma, wSym, xSym)
 
   IMPLICIT NONE

   INTEGER, INTENT (IN) :: nSym
   INTEGER, INTENT (IN) :: Lmax
   INTEGER, INTENT (IN) :: Mmax
   INTEGER, INTENT (IN) :: nmd

   REAL (KIND=r8), INTENT (IN) :: Alfa(Mend1)
   REAL (KIND=r8), INTENT (IN) :: Beta(Mend1)
   REAL (KIND=r8), INTENT (IN) :: Gamma(Mend1)

   INTEGER, INTENT (OUT) :: nCutSym

   REAL (KIND=r8), INTENT (OUT) :: wSym(NxSym)
   REAL (KIND=r8), INTENT (OUT) :: xSym(NxSym,NxSym)

   INTEGER :: n, nn, mm, jj, iErr

   REAL (KIND=r8) :: wkSym(NxSym)
   REAL (KIND=r8) :: xkSym(NxSym,NxSym)

   ! Symmetric Case: s /= 0

   wSym=0.0_r8
   wkSym=0.0_r8
   xkSym=0.0_r8

   DO n=1,Lmax
      xkSym(n,n)=Alfa(2*n)
   END DO
   DO n=1,Lmax
      nn=2*n
      jj=Lmax+n
      xkSym(n,jj)=Beta(nn)
      xkSym(jj,n)=xkSym(n,jj)
      IF (n == Lmax .AND. nmd == 1) EXIT
      xkSym(n,jj+1)=Beta(nn+1)
      xkSym(jj+1,n)=xkSym(n,jj+1)
   END DO
   DO n=1,Mmax
      nn=2*n-1
      jj=Lmax+n
      xkSym(jj,jj)=Alfa(nn)
      mm=jj+Mmax
      xkSym(jj,mm)=Gamma(nn)
      xkSym(mm,jj)=xkSym(jj,mm)
   END DO

   IF (PrtOut) THEN
      WRITE (UNIT=nfprt, FMT='(/,A)') ' xkSym :'
      DO n=1,nSym
         WRITE (UNIT=nfprt, FMT='(1P,6G12.5)') (xkSym(n,nn),nn=1,nSym)
      END DO
   END IF

   CALL tred2 (NxSym, nSym, xkSym, wSym, wkSym, xSym)
   CALL tql2 (NxSym, nSym, wSym, wkSym, xSym, Eps, iErr)

   IF (PrtOut) THEN
      WRITE (UNIT=nfprt, FMT='(/,A,I4)') ' wSym : iErr = ', iErr
      WRITE (UNIT=nfprt, FMT='(1P,6G12.5)') (wSym(n),n=1,nSym)
   END IF

   ! Reordering Frequencies

   CALL Order (NxSym, nSym, wSym, xSym, PerCut, nCutSym)
   WHERE (ABS(xSym) <= Tolx) xSym=0.0_r8

   WRITE (UNIT=nfprt, FMT='(/,2(A,I4))') &
         ' Frequency: nSym = ', nSym, ' nCutSym = ', nCutSym
   WRITE (UNIT=nfprt, FMT='(1P,6G12.5)') (wSym(n)/twoOmega,n=1,nCutSym)
   WRITE (UNIT=nfprt, FMT='(A)')' Period :'
   WRITE (UNIT=nfprt, FMT='(1P,6G12.5)') (1.0_r8/wSym(n),n=1,nCutSym)

   IF (PrtOUt) THEN
      WRITE (UNIT=nfprt, FMT='(/,A)') ' xSym :'
      DO n=1,nSym
         WRITE (UNIT=nfprt, FMT='(1P,6G12.5)') (xSym(n,nn),nn=1,MIN(6,nCutSym))
      END DO
   END IF

END SUBROUTINE SymSne0


SUBROUTINE AsySne0 (nAsy, Lmax, Mmax, nmd, nCutAsy,  &
                    Alfa, Beta, Gamma, wAsy, xAsy)
 
   IMPLICIT NONE

   INTEGER, INTENT (IN) :: nAsy
   INTEGER, INTENT (IN) :: Lmax
   INTEGER, INTENT (IN) :: Mmax
   INTEGER, INTENT (IN) :: nmd

   REAL (KIND=r8), INTENT (IN) :: Alfa(Mend1)
   REAL (KIND=r8), INTENT (IN) :: Beta(Mend1)
   REAL (KIND=r8), INTENT (IN) :: Gamma(Mend1)

   INTEGER, INTENT (OUT) :: nCutAsy

   REAL (KIND=r8), INTENT (OUT) :: wAsy(NxAsy)
   REAL (KIND=r8), INTENT (OUT) :: xAsy(NxAsy,NxAsy)

   INTEGER :: n, nn, mm, jj, iErr

   REAL (KIND=r8) :: wkAsy(NxAsy)
   REAL (KIND=r8) :: xkAsy(NxAsy,NxAsy)

   ! Asymmetric Case: s /= 0

   wAsy=0.0_r8
   wkAsy=0.0_r8
   xkAsy=0.0_r8

   DO n=1,Mmax
      xkAsy(n,n)=Alfa(2*n-1)
   END DO
   DO n=1,Lmax
      nn=2*n
      jj=Mmax+n
      xkAsy(n,jj)=Beta(nn)
      xkAsy(jj,n)=xkAsy(n,jj)
      IF (n == Lmax .AND. nmd == 1) EXIT
      xkAsy(n+1,jj)=Beta(nn+1)
      xkAsy(jj,n+1)=xkAsy(n+1,jj)
   END DO
   DO n=1,Lmax
      nn=2*n
      jj=Mmax+n
      xkAsy(jj,jj)=Alfa(nn)
      mm=jj+Lmax
      xkAsy(jj,mm)=Gamma(nn)
      xkAsy(mm,jj)=xkAsy(jj,mm)
   END DO

   IF (PrtOut) THEN
      WRITE (UNIT=nfprt, FMT='(/,A)') ' xkAsy :'
      DO n=1,nAsy
         WRITE (UNIT=nfprt, FMT='(1P,6G12.5)') (xkAsy(n,nn),nn=1,nAsy)
      END DO
   END IF

   CALL tred2 (NxAsy, nAsy, xkAsy, wAsy, wkAsy, xAsy)
   CALL tql2 (NxAsy, nAsy, wAsy, wkAsy, xAsy, Eps, iErr)

   IF (PrtOut) THEN
      WRITE (UNIT=nfprt, FMT='(/,A,I4)') ' wAsy : iErr = ', iErr
      WRITE (UNIT=nfprt, FMT='(1P,6G12.5)') (wAsy(n),n=1,nAsy)
   END IF

   ! Reordering Frequencies

   CALL Order (NxAsy, nAsy, wAsy, xAsy, PerCut, nCutAsy)
   WHERE (ABS(xAsy) <= Tolx) xAsy=0.0_r8

   WRITE (UNIT=nfprt, FMT='(/,2(A,I4))') &
         ' Frequency: nAsy = ', nAsy, ' nCutAsy = ', nCutAsy
   WRITE (UNIT=nfprt, FMT='(1P,6G12.5)') (wAsy(n)/twoOmega,n=1,nCutAsy)
   WRITE (UNIT=nfprt, FMT='(A)')' Period :'
   WRITE (UNIT=nfprt, FMT='(1P,6G12.5)') (1.0_r8/wAsy(n),n=1,nCutAsy)

   IF (PrtOUt) THEN
      WRITE (UNIT=nfprt, FMT='(/,A)') ' xAsy :'
      DO n=1,nAsy
         WRITE (UNIT=nfprt, FMT='(1P,6G12.5)') (xAsy(n,nn),nn=1,MIN(6,nCutAsy))
      END DO
   END IF

END SUBROUTINE AsySne0


SUBROUTINE Record (NxSymAsy, nCutSymAsy, nSymAsy, wSymAsy, xSymAsy)
 
   IMPLICIT NONE

   INTEGER, INTENT (IN) :: NxSymAsy
   INTEGER, INTENT (IN) :: nCutSymAsy
   INTEGER, INTENT (IN) :: nSymAsy

   REAL (KIND=r8), INTENT (IN) :: wSymAsy(NxSymAsy)
   REAL (KIND=r8), INTENT (IN) :: xSymAsy(NxSymAsy,NxSymAsy)

   INTEGER :: n, nn

   REAL (KIND=r8) :: Period(nCutSymAsy)
   REAL (KIND=r8) :: Vector(nCutSymAsy,nSymAsy)

   WRITE (UNIT=nFPRT, FMT='(A,3I4)') &
         ' From Record: NxSymAsy, nCutSymAsy, nSymAsy : ', &
           NxSymAsy, nCutSymAsy, nSymAsy

   WRITE (UNIT=nfmod) INT(nCutSymAsy,KIND=i8), INT(nSymAsy,KIND=i8)

   IF (nCutSymAsy == 0) RETURN

   Period=0.0_r8
   Vector=0.0_r8

   DO n=1,nCutSymAsy
      Period(n)=1.0_r8/wSymAsy(n)
      DO nn=1,nSymAsy
         Vector(n,nn)=xSymAsy(nn,n)
      END DO
   END DO

   WRITE (UNIT=nfmod) REAL(Period,KIND=r8), REAL(Vector,KIND=r8)

END SUBROUTINE Record


SUBROUTINE Order (NxSymAsy, nSymAsy, wSymAsy, xSymAsy, PerCut, nCutSymAsy)
 
   IMPLICIT NONE

   INTEGER, INTENT (IN) :: NxSymAsy
   INTEGER, INTENT (IN) :: nSymAsy

   REAL (KIND=r8), INTENT (IN) :: PerCut

   REAL (KIND=r8), INTENT (IN OUT) :: wSymAsy(NxSymAsy)
   REAL (KIND=r8), INTENT (IN OUT) :: xSymAsy(NxSymAsy,nSymAsy)

   INTEGER, INTENT (OUT) :: nCutSymAsy

   INTEGER :: k, j, j1, i, jCut, nCut, NxSymAsy1

   REAL (KIND=r8) :: wSymAsyChang

   REAL (KIND=r8) :: xwSymAsyChang(NxSymAsy)

   REAL (KIND=r8) :: xSymAsyChang(NxSymAsy,nSymAsy)

   NxSymAsy1=nSymAsy-1
   10 k=0
   DO j=1,NxSymAsy1
      j1=j+1
      IF (ABS(wSymAsy(j)) > ABS(wSymAsy(j1))) THEN
         wSymAsyChang=wSymAsy(j)
         DO i=1,nSymAsy
            xwSymAsyChang(i)=xSymAsy(i,j)
         END DO
         wSymAsy(j)=wSymAsy(j1)
         DO i=1,nSymAsy
            xSymAsy(i,j)=xSymAsy(i,j1)
         END DO
         wSymAsy(j1)=wSymAsyChang
         DO i=1,nSymAsy
            xSymAsy(i,j1)=xwSymAsyChang(i)
         END DO
         k=1
      END IF
   END DO
   IF (k /= 0) GO TO 10

   IF (PerCut <= 0.0_r8) THEN
      nCutSymAsy=nSymAsy
      RETURN
   END IF

   nCut=0
   DO j=1,nSymAsy
      IF (ABS(1.0_r8/wSymAsy(j)) > PerCut) nCut=j
   END DO
   nCutSymAsy=nSymAsy-nCut
   nCut=nCut+1

   DO j=1,nSymAsy
      xwSymAsyChang(j)=wSymAsy(j)
      DO i=1,nSymAsy
         xSymAsyChang(i,j)=xSymAsy(i,j)
      END DO
   END DO

   wSymAsy=0.0_r8
   xSymAsy=0.0_r8

   DO jCut=nCut,nSymAsy
      j=jCut+1-nCut
      wSymAsy(j)=xwSymAsyChang(jCut)
      DO i=1,nSymAsy
         xSymAsy(i,j)=xSymAsyChang(i,jCut)
      END DO
   END DO

   nCut=nCut-1
   DO jCut=1,nCut
      j=nSymAsy+jCut-nCut
      wSymAsy(j)=xwSymAsyChang(jCut)
      DO i=1,nSymAsy
         xSymAsy(i,j)=xSymAsyChang(i,jCut)
      END DO
   END DO

END SUBROUTINE Order


END MODULE HorizontalModes
