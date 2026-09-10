!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
MODULE Fourier

  USE InputParameters, ONLY : r8, Imax => ImaxOut, &
                              Ifmx => MFactorFourier,  &
                              Imxt => MTrigsFourier

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: InitFFT, ClsMemFFT, fft991

  INTEGER, DIMENSION (:), ALLOCATABLE :: ifax

  REAL (KIND=r8), DIMENSION (:), ALLOCATABLE :: trigs


CONTAINS


SUBROUTINE InitFFT

  IMPLICIT NONE

  INTEGER, SAVE :: ifp=1

  IF (ifp == 1) THEN
     CALL SetMemFFT
     ifp=0
  END IF

  CALL fax    (Imax, 2)
  CALL fftrig (Imax, 2)

END SUBROUTINE InitFFT


SUBROUTINE SetMemFFT

  IMPLICIT NONE

  ALLOCATE (ifax (Ifmx))
  ALLOCATE (trigs(Imxt))

END SUBROUTINE SetMemFFT


SUBROUTINE ClsMemFFT

  IMPLICIT NONE

  DEALLOCATE (ifax)
  DEALLOCATE (trigs)

END SUBROUTINE ClsMemFFT


SUBROUTINE fax (n, mode)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n
  INTEGER, INTENT(IN) :: mode

  INTEGER :: i, l, k, ii, nn, inc, item, nfax, istop

  nn=n
  IF (ABS(mode) == 1) GO TO 10
  IF (ABS(mode) == 8) GO TO 10
  nn=n/2
  IF ((nn+nn) == n) GO TO 10
  ifax(1)=-99
  RETURN

10 CONTINUE
  k=1
20 CONTINUE
  IF (MOD(nn,4) /= 0) GO TO 30
  k=k+1
  ifax(k)=4
  nn=nn/4
  IF (nn == 1) GO TO 80
  GO TO 20
30 CONTINUE
  IF (MOD(nn,2) /= 0) GO TO 40
  k=k+1
  ifax(k)=2
  nn=nn/2
  IF (nn == 1) GO TO 80
40 CONTINUE
  IF (MOD(nn,3) /= 0) GO TO 50
  k=k+1
  ifax(k)=3
  nn=nn/3
  IF (nn == 1) GO TO 80
  GO TO 40
50 CONTINUE
  l=5
  inc=2
60 CONTINUE
  IF (MOD(nn,l) /= 0) GO TO 70
  k=k+1
  ifax(k)=l
  nn=nn/l
  IF (nn == 1) GO TO 80
  GO TO 60
70 CONTINUE
  l=l+inc
  inc=6-inc
  GO TO 60
80 CONTINUE
  ifax(1)=k-1
  nfax=ifax(1)
  IF (nfax /= 1) THEN
     DO ii=2,nfax
        istop=nfax+2-ii
        DO i=2,istop
           IF (ifax(i+1) >= ifax(i)) CYCLE
           item=ifax(i)
           ifax(i)=ifax(i+1)
           ifax(i+1)=item
        END DO
     END DO
  END IF

END SUBROUTINE fax


SUBROUTINE fftrig (n, mode)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n
  INTEGER, INTENT(IN) :: mode

  INTEGER :: i, l, la, nh, nn, imode

  REAL (KIND=r8) :: pi, del, angle

  pi=2.0_r8*ASIN(1.0_r8)
  imode=ABS(mode)
  nn=n
  IF (imode > 1 .AND. imode < 6) nn=n/2
  del=(pi+pi)/REAL(nn,r8)
  l=nn+nn
  DO i=1,l,2
     angle=0.5_r8*REAL(i-1,r8)*del
     trigs(i)=COS(angle)
     trigs(i+1)=SIN(angle)
  END DO
  IF (imode == 1) RETURN
  IF (imode == 8) RETURN
  del=0.5_r8*del
  nh=(nn+1)/2
  l=nh+nh
  la=nn+nn
  DO i=1,l,2
     angle=0.5_r8*REAL(i-1,r8)*del
     trigs(la+i)=COS(angle)
     trigs(la+i+1)=SIN(angle)
  END DO
  IF (imode <= 3) RETURN
  del=0.5_r8*del
  la=la+nn
  IF (mode /= 5) THEN
     DO i=2,nn
        angle=REAL(i-1,r8)*del
        trigs(la+i)=2.0_r8*SIN(angle)
     END DO
  ELSE
     del=0.5_r8*del
     DO i=2,n
        angle=REAL(i-1,r8)*del
        trigs(la+i)=SIN(angle)
     END DO
  END IF

END SUBROUTINE fftrig


SUBROUTINE fft991 (a, work, inc, jump, n, lot, isign)

  ! fft991 :transforms grid point representations of fields into their
  !         fourier representations at all model levels and vice-versa.
  !
  !        i=sqrt(-1)
  !        0<x<2*pai
  !        max=n/2
  !     f(x)=f(0)+f(  1)*exp(-i*  1*x)+f(-  1)*exp(+i*  1*x)
  !              +f(  2)*exp(-i*  2*x)+f(-  2)*exp(+i*  2*x)
  !              +........................
  !              +f(  m)*exp(-i*  m*x)+f(-  m)*exp(+i*  m*x)
  !              +........................
  !              +f(max)*exp(-i*max*x)+f(-max)*exp(+i*max*x)
  !-----------------------------------------------------------------------
  !..  a    (jump,lot)   isign=+1
  !                 input array  for m=0,1,2,....,n/2
  !                      a(2*m+1)=real(f(m))
  !                      a(2*m+2)=imag(f(m))
  !                      a(  1)=cos(m= 0)   ,a(2)=0.
  !                      a(  3)=cos(m= 1)   ,a(4)=sin(m=1)
  !                      a(  5)=cos(m= 2)   ,a(6)=sin(m=2)
  !                      .................................................
  !                      a(n-1)=cos(m=n/2-1),a(n)=sin(m=n/2-1)
  !                      a(n+1)=cos(m=n/2)
  !                 output array
  !                      grid point values
  !                      a(n+1)=0.,......a(jump)=0.
  !..  a    (jump,lot)  isign=-1
  !                 input array
  !                      grid point values
  !                      a(n+1)=0.,......a(jump)=0.
  !                 output array  for m=0,1,2,....,n/2
  !                      a(2*m+1)=real(f(m))
  !                      a(2*m+2)=imag(f(m))
  !                      a(  1)=cos(m= 0)   ,a(2)=0.
  !                      a(  3)=cos(m= 1)   ,a(4)=sin(m= 1)
  !                      a(  5)=cos(m= 2)   ,a(6)=sin(m= 2)
  !                      .................................................
  !                      a(n-1)=cos(m=n/2-1),a(n)=sin(m=n/2-1)
  !                      a(n+1)=cos(m=n/2)
  !..  work (jump,lot)     work  array
  !..  trigs(.gt.1.5*Imax) sin,cos coefficient
  !..  ifax (*)            factors of Imax
  !..  inc                 increment of data spacing
  !..  jump                jump.ge.n+1
  !..  n                   number of grids to be trasformed
  !..  lot                 number of lines transformed simultaneously
  !..  isign=+1            wave to grid transform
  !          -1            grid to wave transform
  !-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: inc, jump, n, lot, isign

  REAL (KIND=r8), DIMENSION(:), INTENT(INOUT) :: a
  REAL (KIND=r8), DIMENSION(:), INTENT(INOUT) :: work

  INTEGER :: i, j, k, l, m, ia, ib, la, nx, nh
  INTEGER :: ink, igo, nfax, ibase, jbase

  nfax=ifax(1)
  nx=n+1
  nh=n/2
  ink=inc+inc
  IF (isign == 1) THEN
     CALL fft99a (a, work, inc, jump, n, lot)
     igo=60
  ELSE
     igo=50
     IF (MOD(nfax,2) /= 1)THEN
        ibase=1
        jbase=1
        DO l=1,lot
           i=ibase
           j=jbase
!cdir nodep
           DO m=1,n
              work(j)=a(i)
              i=i+inc
              j=j+1
           END DO
           ibase=ibase+jump
           jbase=jbase+nx
        END DO
        igo=60
     END IF
  END IF
  ia=1
  la=1
  DO k=1,nfax
     IF (igo == 60) THEN
        CALL vpassm (work(1:), work(2:), a(ia:), a(ia+inc:), &
                     2, ink, nx, jump, lot, nh, ifax(k+1), la)
        igo=50
     ELSE
        CALL vpassm (a(ia:), a(ia+inc:), work(1:), work(2:), &
                     ink, 2, jump, nx, lot, nh, ifax(k+1), la)
        igo=60
     END IF
     la=la*ifax(k+1)
  END DO
  IF (isign == -1) THEN
     CALL fft99b (work, a, inc, jump, n, lot)
  ELSE
     IF (MOD(nfax,2) /= 1)THEN
        ibase=1
        jbase=1
        DO l=1,lot
           i=ibase
           j=jbase
!cdir nodep
           DO m=1,n
              a(j)=work(i)
              i=i+1
              j=j+inc
           END DO
           ibase=ibase+nx
           jbase=jbase+jump
        END DO
     END IF
     ib=n*inc+1
!cdir nodep
     DO l=1,lot
        a(ib)=0.0_r8
        a(ib+inc)=0.0_r8
        ib=ib+jump
     END DO
  END IF

END SUBROUTINE fft991


SUBROUTINE fft99a (a, work, inc, jump, n, lot)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: inc, jump, n, lot

  REAL (KIND=r8), DIMENSION(:), INTENT(IN)  :: a
  REAL (KIND=r8), DIMENSION(:), INTENT(OUT) :: work

  INTEGER :: l, k, ia, ib, ja, jb, nh, nx, ink
  INTEGER :: iabase, ibbase, jabase, jbbase

  REAL (KIND=r8) :: c, s

  nh=n/2
  nx=n+1
  ink=inc+inc
  ia=1
  ib=n*inc+1
  ja=1
  jb=2
!cdir nodep
  DO l=1,lot
     work(ja)=a(ia)+a(ib)
     work(jb)=a(ia)-a(ib)
     ia=ia+jump
     ib=ib+jump
     ja=ja+nx
     jb=jb+nx
  END DO
  iabase=2*inc+1
  ibbase=(n-2)*inc+1
  jabase=3
  jbbase=n-1
  DO k=3,nh,2
     ia=iabase
     ib=ibbase
     ja=jabase
     jb=jbbase
     c=trigs(n+k)
     s=trigs(n+k+1)
!cdir nodep
     DO l=1,lot
        work(ja)=(a(ia)+a(ib))- &
                 (s*(a(ia)-a(ib))+c*(a(ia+inc)+a(ib+inc)))
        work(jb)=(a(ia)+a(ib))+ &
                 (s*(a(ia)-a(ib))+c*(a(ia+inc)+a(ib+inc)))
        work(ja+1)=(c*(a(ia)-a(ib))-s*(a(ia+inc)+a(ib+inc)))+ &
                      (a(ia+inc)-a(ib+inc))
        work(jb+1)=(c*(a(ia)-a(ib))-s*(a(ia+inc)+a(ib+inc)))- &
                      (a(ia+inc)-a(ib+inc))
        ia=ia+jump
        ib=ib+jump
        ja=ja+nx
        jb=jb+nx
     END DO
     iabase=iabase+ink
     ibbase=ibbase-ink
     jabase=jabase+2
     jbbase=jbbase-2
  END DO
  IF (iabase == ibbase) THEN
     ia=iabase
     ja=jabase
!cdir nodep
     DO l=1,lot
        work(ja)=2.0_r8*a(ia)
        work(ja+1)=-2.0_r8*a(ia+inc)
        ia=ia+jump
        ja=ja+nx
     END DO
  END IF

END SUBROUTINE fft99a


SUBROUTINE fft99b (work, a, inc, jump, n, lot)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: inc, jump, n, lot

  REAL (KIND=r8), DIMENSION(:), INTENT(IN)  :: work
  REAL (KIND=r8), DIMENSION(:), INTENT(OUT) :: a

  INTEGER :: l, k, ia, ib, ja, jb, nh, nx, ink
  INTEGER :: iabase, ibbase, jabase, jbbase

  REAL (KIND=r8) :: c, s, scale

  nh=n/2
  nx=n+1
  ink=inc+inc
  scale=1.0_r8/REAL(n,r8)
  ia=1
  ib=2
  ja=1
  jb=n*inc+1
!cdir nodep
  DO l=1,lot
     a(ja)=scale*(work(ia)+work(ib))
     a(jb)=scale*(work(ia)-work(ib))
     a(ja+inc)=0.0_r8
     a(jb+inc)=0.0_r8
     ia=ia+nx
     ib=ib+nx
     ja=ja+jump
     jb=jb+jump
  END DO
  scale=0.5_r8*scale
  iabase=3
  ibbase=n-1
  jabase=2*inc+1
  jbbase=(n-2)*inc+1
  DO k=3,nh,2
     ia=iabase
     ib=ibbase
     ja=jabase
     jb=jbbase
     c=trigs(n+k)
     s=trigs(n+k+1)
!cdir nodep
     DO l=1,lot
        a(ja)=scale*((work(ia)+work(ib)) &
             +(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
        a(jb)=scale*((work(ia)+work(ib)) &
             -(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
        a(ja+inc)=scale*((c*(work(ia)-work(ib))-s*(work(ia+1)+ &
             work(ib+1)))+(work(ib+1)-work(ia+1)))
        a(jb+inc)=scale*((c*(work(ia)-work(ib))-s*(work(ia+1)+ &
             work(ib+1)))-(work(ib+1)-work(ia+1)))
        ia=ia+nx
        ib=ib+nx
        ja=ja+jump
        jb=jb+jump
     END DO
     iabase=iabase+2
     ibbase=ibbase-2
     jabase=jabase+ink
     jbbase=jbbase-ink
  END DO
  IF (iabase == ibbase)THEN
     ia=iabase
     ja=jabase
     scale=2.0_r8*scale
!cdir nodep
     DO l=1,lot
        a(ja)=scale*work(ia)
        a(ja+inc)=-scale*work(ia+1)
        ia=ia+nx
        ja=ja+jump
     END DO
  END IF

END SUBROUTINE fft99b


SUBROUTINE vpassm (a, b, c, d, inc1, inc2, inc3, inc4, lot, n, ifac, la)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: inc1, inc2, inc3, inc4, lot, n, ifac, la

  REAL (KIND=r8), DIMENSION (:), INTENT(IN)  :: a, b
  REAL (KIND=r8), DIMENSION (:), INTENT(OUT) :: c, d

  INTEGER :: i, j, k, l, m
  INTEGER :: ia, ja, ib, jb, kb, ic, jc, kc, id, jd, kd, ie, je, ke
  INTEGER :: igo, ijk, la1
  INTEGER :: iink, jink, jump, ibase, jbase

  REAL (KIND=r8) :: c1, s1,c2, s2, c3, s3, c4, s4
  REAL (KIND=r8) :: wka, wkb, wksina, wksinb, wkaacp, wkbacp, wkaacm, wkbacm

  INTEGER, SAVE :: ifr=0

  REAL (KIND=r8), SAVE :: radi,sin60,sin36,sin72,cos36,cos72

  IF (ifr == 0) THEN
     radi=ATAN(1.0_r8)/45.0_r8
     sin60=SIN(60.0_r8*radi)
     sin36=SIN(36.0_r8*radi)
     sin72=SIN(72.0_r8*radi)
     cos36=COS(36.0_r8*radi)
     cos72=COS(72.0_r8*radi)
     ifr=1
  END IF
  m=n/ifac
  iink=m*inc1
  jink=la*inc2
  jump=(ifac-1)*jink
  ibase=0
  jbase=0
  igo=ifac-1
  IF (igo == 1) THEN
     ia=1
     ja=1
     ib=ia+iink
     jb=ja+jink
     DO l=1,la
        i=ibase
        j=jbase
!cdir nodep
        DO ijk=1,lot
           c(ja+j)=a(ia+i)+a(ib+i)
           d(ja+j)=b(ia+i)+b(ib+i)
           c(jb+j)=a(ia+i)-a(ib+i)
           d(jb+j)=b(ia+i)-b(ib+i)
           i=i+inc3
           j=j+inc4
        END DO
        ibase=ibase+inc1
        jbase=jbase+inc2
     END DO
     IF (la == m) RETURN
     la1=la+1
     jbase=jbase+jump
     DO k=la1,m,la
        kb=k+k-2
        c1=trigs(kb+1)
        s1=trigs(kb+2)
        DO l=1,la
           i=ibase
           j=jbase
!cdir nodep
           DO ijk=1,lot
              wka=a(ia+i)-a(ib+i)
              wkb=b(ia+i)-b(ib+i)
              c(ja+j)=a(ia+i)+a(ib+i)
              d(ja+j)=b(ia+i)+b(ib+i)
              c(jb+j)=c1*wka-s1*wkb
              d(jb+j)=s1*wka+c1*wkb
              i=i+inc3
              j=j+inc4
           END DO
           ibase=ibase+inc1
           jbase=jbase+inc2
        END DO
        jbase=jbase+jump
     END DO
  ELSEIF (igo == 2) THEN
     ia=1
     ja=1
     ib=ia+iink
     jb=ja+jink
     ic=ib+iink
     jc=jb+jink
     DO l=1,la
        i=ibase
        j=jbase
!cdir  nodep
        DO ijk=1,lot
           wka=a(ib+i)+a(ic+i)
           wkb=b(ib+i)+b(ic+i)
           wksina=sin60*(a(ib+i)-a(ic+i))
           wksinb=sin60*(b(ib+i)-b(ic+i))
           c(ja+j)=a(ia+i)+wka
           d(ja+j)=b(ia+i)+wkb
           c(jb+j)=(a(ia+i)-0.5_r8*wka)-wksinb
           c(jc+j)=(a(ia+i)-0.5_r8*wka)+wksinb
           d(jb+j)=(b(ia+i)-0.5_r8*wkb)+wksina
           d(jc+j)=(b(ia+i)-0.5_r8*wkb)-wksina
           i=i+inc3
           j=j+inc4
        END DO
        ibase=ibase+inc1
        jbase=jbase+inc2
     END DO
     IF (la == m) RETURN
     la1=la+1
     jbase=jbase+jump
     DO k=la1,m,la
        kb=k+k-2
        kc=kb+kb
        c1=trigs(kb+1)
        s1=trigs(kb+2)
        c2=trigs(kc+1)
        s2=trigs(kc+2)
        DO l=1,la
           i=ibase
           j=jbase
!cdir nodep
           DO ijk=1,lot
              wka=a(ib+i)+a(ic+i)
              wkb=b(ib+i)+b(ic+i)
              wksina=sin60*(a(ib+i)-a(ic+i))
              wksinb=sin60*(b(ib+i)-b(ic+i))
              c(ja+j)=a(ia+i)+wka
              d(ja+j)=b(ia+i)+wkb
              c(jb+j)=c1*((a(ia+i)-0.5_r8*wka)-wksinb) &
                     -s1*((b(ia+i)-0.5_r8*wkb)+wksina)
              d(jb+j)=s1*((a(ia+i)-0.5_r8*wka)-wksinb) &
                     +c1*((b(ia+i)-0.5_r8*wkb)+wksina)
              c(jc+j)=c2*((a(ia+i)-0.5_r8*wka)+wksinb) &
                     -s2*((b(ia+i)-0.5_r8*wkb)-wksina)
              d(jc+j)=s2*((a(ia+i)-0.5_r8*wka)+wksinb) &
                     +c2*((b(ia+i)-0.5_r8*wkb)-wksina)
              i=i+inc3
              j=j+inc4
           END DO
           ibase=ibase+inc1
           jbase=jbase+inc2
        END DO
        jbase=jbase+jump
     END DO
  ELSEIF (igo == 3) THEN
     ia=1
     ja=1
     ib=ia+iink
     jb=ja+jink
     ic=ib+iink
     jc=jb+jink
     id=ic+iink
     jd=jc+jink
     DO l=1,la
        i=ibase
        j=jbase
!cdir nodep
        DO ijk=1,lot
           wkaacp=a(ia+i)+a(ic+i)
           wkbacp=b(ia+i)+b(ic+i)
           wkaacm=a(ia+i)-a(ic+i)
           wkbacm=b(ia+i)-b(ic+i)
           c(ja+j)=wkaacp+(a(ib+i)+a(id+i))
           c(jc+j)=wkaacp-(a(ib+i)+a(id+i))
           d(ja+j)=wkbacp+(b(ib+i)+b(id+i))
           d(jc+j)=wkbacp-(b(ib+i)+b(id+i))
           c(jb+j)=wkaacm-(b(ib+i)-b(id+i))
           c(jd+j)=wkaacm+(b(ib+i)-b(id+i))
           d(jb+j)=wkbacm+(a(ib+i)-a(id+i))
           d(jd+j)=wkbacm-(a(ib+i)-a(id+i))
           i=i+inc3
           j=j+inc4
        END DO
        ibase=ibase+inc1
        jbase=jbase+inc2
     END DO
     IF (la == m) RETURN
     la1=la+1
     jbase=jbase+jump
     DO k=la1,m,la
        kb=k+k-2
        kc=kb+kb
        kd=kc+kb
        c1=trigs(kb+1)
        s1=trigs(kb+2)
        c2=trigs(kc+1)
        s2=trigs(kc+2)
        c3=trigs(kd+1)
        s3=trigs(kd+2)
        DO l=1,la
           i=ibase
           j=jbase
!cdir nodep
           DO ijk=1,lot
              wkaacp=a(ia+i)+a(ic+i)
              wkbacp=b(ia+i)+b(ic+i)
              wkaacm=a(ia+i)-a(ic+i)
              wkbacm=b(ia+i)-b(ic+i)
              c(ja+j)=wkaacp+(a(ib+i)+a(id+i))
              d(ja+j)=wkbacp+(b(ib+i)+b(id+i))
              c(jc+j)= c2*(wkaacp-(a(ib+i)+a(id+i))) &
                      -s2*(wkbacp-(b(ib+i)+b(id+i))) 
              d(jc+j)= s2*(wkaacp-(a(ib+i)+a(id+i))) &
                      +c2*(wkbacp-(b(ib+i)+b(id+i)))
              c(jb+j)= c1*(wkaacm-(b(ib+i)-b(id+i))) &
                      -s1*(wkbacm+(a(ib+i)-a(id+i)))
              d(jb+j)= s1*(wkaacm-(b(ib+i)-b(id+i))) &
                      +c1*(wkbacm+(a(ib+i)-a(id+i)))
              c(jd+j)= c3*(wkaacm+(b(ib+i)-b(id+i))) &
                      -s3*(wkbacm-(a(ib+i)-a(id+i)))
              d(jd+j)= s3*(wkaacm+(b(ib+i)-b(id+i))) &
                      +c3*(wkbacm-(a(ib+i)-a(id+i)))
              i=i+inc3
              j=j+inc4
           END DO
           ibase=ibase+inc1
           jbase=jbase+inc2
        END DO
        jbase=jbase+jump
     END DO
  ELSEIF (igo == 4) THEN
     ia=1
     ja=1
     ib=ia+iink
     jb=ja+jink
     ic=ib+iink
     jc=jb+jink
     id=ic+iink
     jd=jc+jink
     ie=id+iink
     je=jd+jink
     DO l=1,la
        i=ibase
        j=jbase
!cdir nodep
        DO ijk=1,lot
           c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
           d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
           c(jb+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+ &
                a(id+i)))-(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
           c(je+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+ &
                a(id+i)))+(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
           d(jb+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+ &
                b(id+i)))+(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
           d(je+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+ &
                b(id+i)))-(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
           c(jc+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+ &
                a(id+i)))-(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
           c(jd+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+ &
                a(id+i)))+(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
           d(jc+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+ &
                b(id+i)))+(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
           d(jd+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+ &
                b(id+i)))-(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
           i=i+inc3
           j=j+inc4
        END DO
        ibase=ibase+inc1
        jbase=jbase+inc2
     END DO
     IF (la == m) RETURN
     la1=la+1
     jbase=jbase+jump
     DO k=la1,m,la
        kb=k+k-2
        kc=kb+kb
        kd=kc+kb
        ke=kd+kb
        c1=trigs(kb+1)
        s1=trigs(kb+2)
        c2=trigs(kc+1)
        s2=trigs(kc+2)
        c3=trigs(kd+1)
        s3=trigs(kd+2)
        c4=trigs(ke+1)
        s4=trigs(ke+2)
        DO l=1,la
           i=ibase
           j=jbase
!cdir nodep
           DO ijk=1,lot
              c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
              d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
              c(jb+j)=c1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36* &
                   (a(ic+i)+a(id+i)))-(sin72*(b(ib+i)-b(ie+i)) &
                   +sin36*(b(ic+i)-b(id+i)))) &
                   -s1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36* &
                   (b(ic+i)+b(id+i)))+(sin72*(a(ib+i)-a(ie+i)) &
                   +sin36*(a(ic+i)-a(id+i))))
              d(jb+j)=s1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36* &
                   (a(ic+i)+a(id+i)))-(sin72*(b(ib+i)-b(ie+i)) &
                   +sin36*(b(ic+i)-b(id+i)))) &
                   +c1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36* &
                   (b(ic+i)+b(id+i)))+(sin72*(a(ib+i)-a(ie+i)) &
                   +sin36*(a(ic+i)-a(id+i))))
              c(je+j)=c4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36* &
                   (a(ic+i)+a(id+i)))+(sin72*(b(ib+i)-b(ie+i)) &
                   +sin36*(b(ic+i)-b(id+i)))) &
                   -s4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36* &
                   (b(ic+i)+b(id+i)))-(sin72*(a(ib+i)-a(ie+i))+ &
                   sin36*(a(ic+i)-a(id+i))))
              d(je+j)=s4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36* &
                   (a(ic+i)+a(id+i)))+(sin72*(b(ib+i)-b(ie+i))+ &
                   sin36*(b(ic+i)-b(id+i))))+c4*((b(ia+i)+ &
                   cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
                   -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
              c(jc+j)=c2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72* &
                   (a(ic+i)+a(id+i)))-(sin36*(b(ib+i)-b(ie+i))- &
                   sin72*(b(ic+i)-b(id+i))))-s2*((b(ia+i)-cos36* &
                   (b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))+(sin36 &
                   *(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))) 
              d(jc+j)=s2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72* &
                   (a(ic+i)+a(id+i)))-(sin36*(b(ib+i)-b(ie+i))- &
                   sin72*(b(ic+i)-b(id+i))))+c2*((b(ia+i)-cos36* &
                   (b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))+(sin36 &
                   *(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
              c(jd+j)=c3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72* &
                   (a(ic+i)+a(id+i)))+(sin36*(b(ib+i)-b(ie+i))- &
                   sin72*(b(ic+i)-b(id+i))))-s3*((b(ia+i)-cos36* &
                   (b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))-(sin36 &
                   *(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
              d(jd+j)=s3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72* &
                   (a(ic+i)+a(id+i)))+(sin36*(b(ib+i)-b(ie+i))-sin72 &
                   *(b(ic+i)-b(id+i))))+c3*((b(ia+i)-cos36*(b(ib+i)+ &
                   b(ie+i))+cos72*(b(ic+i)+b(id+i)))-(sin36* &
                   (a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
              i=i+inc3
              j=j+inc4
           END DO
           ibase=ibase+inc1
           jbase=jbase+inc2
        END DO
        jbase=jbase+jump
     END DO
  ENDIF

END SUBROUTINE vpassm


END MODULE Fourier
