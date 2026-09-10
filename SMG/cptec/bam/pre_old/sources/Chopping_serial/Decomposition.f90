!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
MODULE Decomposition

  USE InputParameters, ONLY : r8, &
                              Mend1 => Mend1Out, Mend2 => Mend2Out, Mnwv2 => Mnwv2Out, &
                              Imax => ImaxOut, Jmax => JmaxOut, Kmax => KmaxOut, &
                              Imx => ImxOut, Jmaxhf => JmaxhfOut, ICaseRec

  USE Fourier, ONLY : fft991

  USE Legendre, ONLY : SymAsy, FourierToSpecCoef, &
                       GetConv, GetVort, la0, colrad


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: DectoSpherHarm, UVtoDivgVort, TransSpherHarm
 

CONTAINS


SUBROUTINE DectoSpherHarm (f, qf, Ldim)

  IMPLICIT NONE

  INTEGER, INTENT (IN) :: Ldim

  REAL (KIND=r8), DIMENSION (Imax,Jmax,Ldim), INTENT (IN)  :: f
  REAL (KIND=r8), DIMENSION (Mnwv2,Ldim),     INTENT (OUT) :: qf

  INTEGER :: i, j, k, inc

  REAL (KIND=r8), DIMENSION (Imx,Ldim) :: fn, fs, gn, gs

  qf=0.0_r8
  inc=1

  DO j=1,Jmaxhf
     DO k=1,Ldim
        DO i=1,Imax
           fn(i,k)=f(i,j,k)
           fs(i,k)=f(i,Jmax+1-j,k)
        ENDDO
        DO i=Imax+1,Imx
           fn(i,k)=0.0_r8
           fs(i,k)=0.0_r8
        ENDDO
     ENDDO

     CALL fft991 (fn(1:,1), gn(1:,1), inc, Imx, Imax, Ldim, ICaseRec)
     CALL fft991 (fs(1:,1), gs(1:,1), inc, Imx, Imax, Ldim, ICaseRec)
     CALL SymAsy (fn, fs, Ldim)
     CALL FourierToSpecCoef (fn, fs, qf, Ldim, j)

  ENDDO

END SUBROUTINE DectoSpherHarm


SUBROUTINE UVtoDivgVort (u, v, qdiv, qrot)

  IMPLICIT NONE

  REAL (KIND=r8), DIMENSION (Imax,Jmax,Kmax), INTENT (IN)  :: u, v
  REAL (KIND=r8), DIMENSION (Mnwv2,Kmax)    , INTENT (OUT) :: qrot, qdiv

  INTEGER :: i, j, k, inc

  REAL (KIND=r8), DIMENSION (Jmaxhf) :: coslat

  REAL (KIND=r8), DIMENSION (Imx,Kmax) :: un, us, vn, vs, gw

  qrot=0.0_r8
  qdiv=0.0_r8
  inc=1

  coslat=SIN(colrad)
  DO j=1,Jmaxhf
     DO k=1,Kmax
        DO i=1,Imax
           un(i,k)=u(i,j,k)*coslat(j)
           us(i,k)=u(i,Jmax+1-j,k)*coslat(j)
           vn(i,k)=v(i,j,k)*coslat(j)
           vs(i,k)=v(i,Jmax+1-j,k)*coslat(j)
        ENDDO
        DO i=Imax+1,Imx
           un(i,k)=0.0_r8
           us(i,k)=0.0_r8
           vn(i,k)=0.0_r8
           vs(i,k)=0.0_r8
        ENDDO
     ENDDO
     CALL fft991 (un(1:,1), gw(1:,1), inc, Imx, Imax, Kmax, ICaseRec)
     CALL fft991 (us(1:,1), gw(1:,1), inc, Imx, Imax, Kmax, ICaseRec)
     CALL fft991 (vn(1:,1), gw(1:,1), inc, Imx, Imax, Kmax, ICaseRec)
     CALL fft991 (vs(1:,1), gw(1:,1), inc, Imx, Imax, Kmax, ICaseRec)
     CALL SymAsy (un, us, Kmax)
     CALL SymAsy (vn, vs, Kmax)
     CALL GetConv  (us, un, vs, vn, qdiv, j)
     CALL GetVort  (us, un, vs, vn, qrot, j)
  ENDDO
 
! Get_Conv computes convergence, 
!          signal of qdiv must be changed to get divergence

  qdiv=-qdiv

END SUBROUTINE UVtoDivgVort


SUBROUTINE TransSpherHarm (Ldim, qCoef, ICase)

  IMPLICIT NONE

  INTEGER, INTENT (IN) :: Ldim, ICase

  REAL (KIND=r8), DIMENSION (Mnwv2,Ldim), INTENT (IN OUT) :: qCoef

  INTEGER :: k, l, lx, mn, mm, Nmax, nn

  REAL (KIND=r8), DIMENSION (Mnwv2) :: qWork

  qWork=0.0_r8
  IF (ICase == 1) THEN
     DO k=1,Ldim
        l=0
        DO mm=1,Mend1
           Nmax=Mend2-mm
           DO nn=1,Nmax
              l=l+1
              lx=la0(mm,nn)
              qWork(2*l-1)=qCoef(2*lx-1,k)
              qWork(2*l)=qCoef(2*lx,k)
           ENDDO
        ENDDO
        DO mn=1,Mnwv2
           qCoef(mn,k)=qWork(mn)
        ENDDO
     ENDDO
  ELSE
     DO k=1,Ldim
        l=0
        DO mm=1,Mend1
           Nmax=Mend2-mm
           DO nn=1,Nmax
              l=l+1
              lx=la0(mm,nn)
              qWork(2*lx-1)=qCoef(2*l-1,k)
              qWork(2*lx)=qCoef(2*l,k)
           ENDDO
        ENDDO
        DO mn=1,Mnwv2
           qCoef(mn,k)=qWork(mn)
        ENDDO
     ENDDO
  ENDIF

END SUBROUTINE TransSpherHarm


END MODULE Decomposition
