!
!  $Author: pkubota $
!  $Date: 2006/10/30 18:41:06 $
!  $Revision: 1.3 $
!
MODULE SigmaToPressure

  USE Constants, ONLY : r8, nferr,ExtrapoAdiabatica,nfprt
  USE Sizes, ONLY : Ibmax, Jbmax, Kmax, Lmax, ibmaxperjb
  USE GaussSigma, ONLY : a_hybr, b_hybr
  USE Parallelism, ONLY : myid

  IMPLICIT NONE

  !     subprogram  documentation  block
  ! 
  !     subprogram: sigtop         sigma to pressure calculation
  !     author: j. sela            org: w/nmc42    date: 18 nov 83
  ! 
  !     abstract: interpolates values on the gaussian grid from the
  !     sigma coordinate system to the mandatory pressure levels.
  !     assumes that relative humidity, temperature, horizontal
  !     wind components and vertical velocity vary in the vertical
  !     with the log of pressure. obtains height at pressure levels
  !     by integrating the hydrostatic equation.
  ! 
  !     usage:  call sigtop (ts, psmb, top, tvp, zp, si, alnpmd)
  ! 
  !     input variables:
  ! 
  !     names       meaning/content/purpose/units/type        interface
  !     -----       ----------------------------------        ---------
  ! 
  !     ts          real array of absolute tvirt (k)          arg list
  !                 in sigma layers.
  !     psmb        real array of surface pressure            arg list
  !                 (milibars)
  !     top         real array of terrain height (m)          arg list
  !     si          real array of dimensionless values        arg list
  !                 used to define p at interface of
  !                 model layers.
  !     alnpmd      real array of dimensionless values        arg list
  !                 of log(pmand).
  ! 
  !     output variables:
  ! 
  !     names       meaning/content/purpose/units/type        interface
  !     -----       ----------------------------------        ---------
  ! 
  !     tvp         real array of virtual temperature on      arg list
  !                 mandatory pressure levels. (k)
  ! 
  !     zp          real array of geopotential height on      arg list
  !                 mandatory pressure levels. (m)
  ! 
  !     attributes: language: fortran

  PRIVATE

  PUBLIC :: sig2po, sig2pz, sigtop, gavint


CONTAINS


  SUBROUTINE sig2po (psmb, alnpmd, bfi1, bfo1, bfi2, bfo2, bfi3, bfo3)

    USE Constants, ONLY : Po, Pt, RdByCp, CpByRd

    IMPLICIT NONE

    REAL (KIND=r8), INTENT(IN ) :: psmb   (Ibmax,Jbmax)
    REAL (KIND=r8), INTENT(IN ) :: alnpmd (Lmax)
    REAL (KIND=r8), INTENT(IN ) :: bfi1   (Ibmax,Kmax,Jbmax)
    REAL (KIND=r8), INTENT(OUT) :: bfo1   (Ibmax,Lmax,Jbmax)

    REAL (KIND=r8), INTENT(IN ), OPTIONAL :: bfi2(Ibmax,Kmax,Jbmax)
    REAL (KIND=r8), INTENT(OUT), OPTIONAL :: bfo2(Ibmax,Lmax,Jbmax)
    REAL (KIND=r8), INTENT(IN ), OPTIONAL :: bfi3(Ibmax,Kmax,Jbmax)
    REAL (KIND=r8), INTENT(OUT), OPTIONAL :: bfo3(Ibmax,Lmax,Jbmax)

    REAL (KIND=r8) :: wtl

    INTEGER :: nb, i, j, k, kp
    INTEGER :: ks(Ibmax,Lmax)

    REAL (KIND=r8) :: ddp   (Ibmax,Lmax)
    REAL (KIND=r8) :: pim   (Ibmax,Kmax)
    REAL (KIND=r8) :: alnpm (Ibmax,Kmax)
    REAL (KIND=r8) :: alnp  (Ibmax,Kmax+1)
    REAL (KIND=r8) :: pi    (Ibmax,Kmax+1)
    REAL (KIND=r8) :: p     (Ibmax,Kmax+1)
    REAL (KIND=r8) :: rdelp (Ibmax,Kmax-1)
    REAL (KIND=r8) :: b1    (Ibmax,Kmax-1)
    REAL (KIND=r8) :: b2    (Ibmax,Kmax-1)
    REAL (KIND=r8) :: b3    (Ibmax,Kmax-1)

    ! initialization

    nb=1
    IF (PRESENT(bfi2)) THEN
       IF (PRESENT(bfo2)) THEN
          nb=2
       ELSE
          WRITE (UNIT=nferr, FMT='(A)') ' sig2po: bfo2 required '
          STOP 9001
       END IF
       IF (PRESENT(bfi3)) THEN
          IF (PRESENT(bfo3)) THEN
             nb=3
          ELSE
             WRITE (UNIT=nferr, FMT='(A)') ' sig2po: bfo3 required '
             STOP 9001
          END IF
       END IF
    END IF

    wtl=LOG(Po)

    ! compute pressure on each sigma surface 
    ! set the highest to Pt (See Constants)

    DO j=1,Jbmax
       DO k=1,Kmax
          DO i=1,Ibmaxperjb(j)
                 ! transform from Pa to mbar
             p(i,k)= 1.e-2_r8 * a_hybr(k) + b_hybr(k)*psmb(i,j)
          END DO
       END DO
       p(:,Kmax+1)=Pt

    ! log of pressure and exner function on each sigma surface.

       DO k=1,Kmax+1
          DO i=1,Ibmaxperjb(j)
             alnp(i,k)=LOG(MAX(p(i,k),0.00000001_r8))
             pi(i,k)=EXP(RdByCp*(alnp(i,k)-wtl))
          END DO
       END DO
!       IF(myid ==0)THEN
!          WRITE (UNIT=nfprt, FMT=*) '  exp(alnpmd(k)),  alnpmd(k)'
!          DO k=1,Lmax
!             WRITE (UNIT=nfprt, FMT=*) k,  exp(alnpmd(k)),  alnpmd(k)
!          END DO
!
!          WRITE (UNIT=nfprt, FMT=*)   'exp(alnp(k)),  alnp(k)'
!          DO k=1,kmax
!             WRITE (UNIT=nfprt, FMT=*) k,  exp(alnp(1,k)),  alnp(1,k),psmb(1,j)
!          END DO
!        END IF
    ! mean value of the exner function,
    ! log of pim and 
    ! log of the mean value of pressure 
    ! in each layer

       DO k=1,Kmax
          DO i=1,Ibmaxperjb(j)
             pim(i,k)=(pi(i,k+1)-pi(i,k))/ &
                  (RdByCp*(alnp(i,k+1)-alnp(i,k)))
             alnpm(i,k)=wtl+CpByRd*LOG(pim(i,k))
          END DO
       END DO

    ! interpolated t absolute. ts replaced by tm

       IF (nb == 3) THEN
          DO k=1,Kmax-1
             DO i=1,Ibmaxperjb(j)
                rdelp(i,k)=1.0_r8/(alnpm(i,k+1)-alnpm(i,k))
                b1(i,k)=(bfi1(i,k+1,j)-bfi1(i,k,j))*rdelp(i,k)
                b2(i,k)=(bfi2(i,k+1,j)-bfi2(i,k,j))*rdelp(i,k)
                b3(i,k)=(bfi3(i,k+1,j)-bfi3(i,k,j))*rdelp(i,k)
             END DO
          END DO
       ELSE IF (nb == 2) THEN
          DO k=1,Kmax-1
             DO i=1,Ibmaxperjb(j)
                rdelp(i,k)=1.0_r8/(alnpm(i,k+1)-alnpm(i,k))
                b1(i,k)=(bfi1(i,k+1,j)-bfi1(i,k,j))*rdelp(i,k)
                b2(i,k)=(bfi2(i,k+1,j)-bfi2(i,k,j))*rdelp(i,k)
             END DO
          END DO
       ELSE 
          DO k=1,Kmax-1
             DO i=1,Ibmaxperjb(j)
                rdelp(i,k)=1.0_r8/(alnpm(i,k+1)-alnpm(i,k))
                b1(i,k)=(bfi1(i,k+1,j)-bfi1(i,k,j))*rdelp(i,k)
             END DO
          END DO
       END IF

    ! heights to pressure surfaces  -  hydrostatic interpolation
    ! winds and temps by linear interpolation with ln(p)
    ! search for middle of sigma layer above kp

       ks(:,:)=Kmax
       DO kp=1,Lmax
          DO k=1,Kmax-1
             DO i=1,Ibmaxperjb(j)
                IF ((ks(i,kp) == Kmax) .AND. (alnpmd(kp) > alnpm(i,k))) THEN
                   ks(i,kp)=k
                END IF
             END DO
          END DO
       END DO

    ! find values and slopes for upward or downward extrapolation
    ! as well as for interpolation away from tropopause

       IF (nb == 3) THEN
          DO kp=1,Lmax
             DO i=1,Ibmaxperjb(j)
                IF (ks(i,kp) /= 1) THEN
                   ddp(i,kp)=alnpmd(kp)-alnpm(i,ks(i,kp))
                   bfo1(i,kp,j)=bfi1(i,ks(i,kp),j)+b1(i,ks(i,kp)-1)*ddp(i,kp)
                   bfo2(i,kp,j)=bfi2(i,ks(i,kp),j)+b2(i,ks(i,kp)-1)*ddp(i,kp)
                   bfo3(i,kp,j)=bfi3(i,ks(i,kp),j)+b3(i,ks(i,kp)-1)*ddp(i,kp)
                ELSE
                   bfo1(i,kp,j)=bfi1(i,1,j)
                   bfo2(i,kp,j)=bfi2(i,1,j)
                   bfo3(i,kp,j)=bfi3(i,1,j)
                END IF
             END DO
          END DO
       ELSE IF (nb == 2) THEN
          DO kp=1,Lmax
             DO i=1,Ibmaxperjb(j)
                IF (ks(i,kp) /= 1) THEN
                   ddp(i,kp)=alnpmd(kp)-alnpm(i,ks(i,kp))
                   bfo1(i,kp,j)=bfi1(i,ks(i,kp),j)+b1(i,ks(i,kp)-1)*ddp(i,kp)
                   bfo2(i,kp,j)=bfi2(i,ks(i,kp),j)+b2(i,ks(i,kp)-1)*ddp(i,kp)
                ELSE
                   bfo1(i,kp,j)=bfi1(i,1,j)
                   bfo2(i,kp,j)=bfi2(i,1,j)
                END IF
             END DO
          END DO
       ELSE
          DO kp=1,Lmax
             DO i=1,Ibmaxperjb(j)
                IF (ks(i,kp) /= 1) THEN
                   ddp(i,kp)=alnpmd(kp)-alnpm(i,ks(i,kp))
                   bfo1(i,kp,j)=bfi1(i,ks(i,kp),j)+b1(i,ks(i,kp)-1)*ddp(i,kp)
                ELSE
                   bfo1(i,kp,j)=bfi1(i,1,j)
                END IF
             END DO
          END DO
       END IF
    END DO

  END SUBROUTINE sig2po


  SUBROUTINE sig2pz (ts, psmb, top, tvp, zp, alnpmd)

    USE Constants, ONLY : Po, Pt, RdByCp, RdByGrav, CpByRd, &
                          Tref, Zref, TVVTa, TVVTb

    IMPLICIT NONE

    REAL (KIND=r8), INTENT(IN) :: ts(Ibmax,Kmax,Jbmax)
    REAL (KIND=r8), INTENT(IN) :: psmb(Ibmax,Jbmax)
    REAL (KIND=r8), INTENT(IN) :: top(Ibmax,Jbmax)
    REAL (KIND=r8), INTENT(OUT) :: tvp(Ibmax,Lmax,Jbmax)
    REAL (KIND=r8), INTENT(OUT) :: zp(Ibmax,Lmax,Jbmax)
    REAL (KIND=r8), INTENT(IN) :: alnpmd(Lmax)

    INTEGER :: i, j, k, kp
    INTEGER :: ks(Ibmax,Lmax)

    REAL (KIND=r8) :: wtl, ddp, hmh, part

    REAL (KIND=r8), DIMENSION(Ibmax) :: gamma, tmsl, tstr, zlay
    REAL (KIND=r8), DIMENSION(Ibmax,Kmax) :: pim, alnpm, dh, bz, zt
    REAL (KIND=r8), DIMENSION(Ibmax,Kmax+1) :: alnp, pi, p, z, tlev, balnp
    REAL (KIND=r8), DIMENSION(Ibmax,Kmax-1) :: rdelp, btv

    ! constants

    wtl=LOG(Po)

    ! compute pressure on each sigma surface and
    ! set the top of the model at Pt (See Constants)

    DO j=1,Jbmax
       DO k=1,Kmax
          DO i=1,Ibmaxperjb(j)
             ! transform from Pa to mbar
             p(i,k) = 1.e-2_r8 * a_hybr(k) + b_hybr(k) * psmb(i,j)
          END DO
       END DO
       p(:,Kmax+1)=Pt

    ! compute the log of pressure and exner function on each sigma surface

       DO k=1,Kmax+1
          DO i=1,Ibmaxperjb(j)
             alnp(i,k)=LOG(MAX(p(i,k),0.00000001_r8))
             pi(i,k)=EXP(RdByCp*(alnp(i,k)-wtl))
          END DO
       END DO

    ! compute the mean value of the exner function and
    ! the log of the mean value of pressure in each layer.

       DO k=1,Kmax
          DO i=1,Ibmaxperjb(j)
             pim(i,k)=(pi(i,k+1)-pi(i,k))/ &
                  (RdByCp*(alnp(i,k+1)-alnp(i,k)))
             alnpm(i,k)=wtl+CpByRd*LOG(pim(i,k))
          END DO
       END DO

    ! compute the height of the sigma surfaces.

       DO i=1,Ibmaxperjb(j)
          z(i,1)=top(i,j)
       END DO
       DO k=1,Kmax
          DO i=1,Ibmaxperjb(j)
             z(i,k+1)=z(i,k)+RdByGrav*ts(i,k,j)* &
                                    (alnp(i,k)-alnp(i,k+1))
          END DO
       END DO

    ! compute the temperature on sigma surfaces. model the temp
    ! profile as linear in log of pressure. extrapolate to determine
    ! the value of temp at the highest and lowest sigma surface.

       DO k=2,Kmax
          DO i=1,Ibmaxperjb(j)
             tlev(i,k)=((alnp(i,k)-alnpm(i,k))*ts(i,k-1,j)+ &
                          (alnpm(i,k-1)-alnp(i,k))*ts(i,k,j))/ &
                          (alnpm(i,k-1)-alnpm(i,k))
          END DO
       END DO
       DO i=1,Ibmaxperjb(j)
          tlev(i,1)=((alnp(i,1)-alnpm(i,2))*ts(i,1,j)+ &
                          (alnpm(i,1)-alnp(i,1))*ts(i,2,j))/ &
                          (alnpm(i,1)-alnpm(i,2))
          tlev(i,Kmax+1)= &
               ((alnp(i,Kmax+1)-alnpm(i,Kmax))*ts(i,Kmax-1,j)+ &
               (alnpm(i,Kmax-1)-alnp(i,Kmax+1))*ts(i,Kmax,j))/ &
               (alnpm(i,Kmax-1)-alnpm(i,Kmax))
       END DO
       DO i=1,Ibmaxperjb(j)
          zlay(i)=0.5_r8*(z(i,1)+z(i,2))
          tstr(i)=ts(i,1,j)+TVVTa*(zlay(i)-z(i,1))
          tmsl(i)=ts(i,1,j)+TVVTa*zlay(i)
          gamma(i)=0.0_r8
          IF (tmsl(i) > Tref) THEN
             tmsl(i)=Tref
             IF (tstr(i) > Tref) tmsl(i)=Tref- &
                  TVVTb*(tstr(i)-Tref)*(tstr(i)-Tref)
          END IF
          IF (z(i,1) > Zref) gamma(i)=(tstr(i)-tmsl(i))/z(i,1)
       END DO
       DO k=1,Kmax
          DO i=1,Ibmaxperjb(j)
             dh(i,k)=alnp(i,k+1)-alnp(i,k)
             bz(i,k)=RdByGrav*(tlev(i,k+1)-tlev(i,k))/dh(i,k)
             zt(i,k)=0.5_r8*(z(i,k+1)+z(i,k))+ &
                  0.125_r8*bz(i,k)*dh(i,k)*dh(i,k)
             balnp(i,k)=0.5_r8*(alnp(i,k+1)+alnp(i,k))
          END DO
       END DO
       DO k=1,Kmax-1
          DO i=1,Ibmaxperjb(j)
             rdelp(i,k)=1.0_r8/(alnpm(i,k+1)-alnpm(i,k))
             btv(i,k)=(ts(i,k+1,j)-ts(i,k,j))*rdelp(i,k)
          END DO
       END DO

    ! heights to pressure surfaces  -  hydrostatic interpolation
    ! winds and temps by linear interpolation with ln(p).
    ! search for middle of sigma layer above kp
    ! if fall thru must extrapolate up

       ks(:,:) = Kmax
       DO kp=1,Lmax
          DO k=1,Kmax-1
             DO i=1,Ibmaxperjb(j)
                IF ((ks(i,kp) == Kmax) .AND. (alnpmd(kp) > alnpm(i,k))) THEN
                   ks(i,kp) = k
                END IF
             END DO
          END DO
       END DO

    ! find values and slopes for upward or downward extrapolation
    ! as well as for interpolation away from tropopause

       DO kp=1,Lmax
          DO i=1,Ibmaxperjb(j)
             IF (ks(i,kp) /= 1) THEN
                ddp=alnpmd(kp)-alnpm(i,ks(i,kp))
                tvp(i,kp,j)=ts(i,ks(i,kp),j)+btv(i,ks(i,kp)-1)*ddp
             ELSE
                tvp(i,kp,j)=ts(i,1,j)+btv(i,1)*(alnpmd(kp)-alnpm(i,1))
             END IF
          END DO
       END DO

    ! start with another search thru the sigma levels
    ! to find the level (ks)  above the desires mandatory level (kp)

       ks(:,:) = Kmax+1
       DO kp=1,Lmax
          DO k=1,Kmax
             DO i=1,Ibmaxperjb(j)
                IF ((ks(i,kp) == Kmax+1) .AND. (alnpmd(kp) > alnp(i,k))) THEN
                   ks(i,kp) = k
                END IF
             END DO
          END DO
       END DO
       DO kp=1,Lmax
          DO i=1,Ibmaxperjb(j)
             IF (ks(i,kp) == 1) THEN
                ! zp(i,kp) is below ground.
                part=RdByGrav*(alnp(i,1)-alnpmd(kp))
                zp(i,kp,j)=z(i,1)+tstr(i)*part/(1.0_r8-0.5_r8*part*gamma(i))
             ELSE
                ! zp references a pressure level in the free air
                hmh=alnpmd(kp)-balnp(i,ks(i,kp)-1)
                zp(i,kp,j)=zt(i,ks(i,kp)-1)- &
                     (RdByGrav*ts(i,ks(i,kp)-1,j)+ &
                     0.5_r8*bz(i,ks(i,kp)-1)*hmh)*hmh
             END IF
          END DO
       END DO
    END DO

  END SUBROUTINE sig2pz


  SUBROUTINE sigtop (tm, gts, gsh, gss, psmb, tg, rg,rq, pmand, alnpmd)

    USE Constants, ONLY : Po, Pt, CTv, a, b, To, Eo, eps, eps1, &
                          RdByCp, RdByGrav, CpByRd, PRHcut, RHmin, RHmax, SHmin

    IMPLICIT NONE

    REAL (KIND=r8), INTENT(OUT) :: tm(Ibmax,Kmax,Jbmax)
    REAL (KIND=r8), INTENT(IN) :: gts(Ibmax,Kmax,Jbmax)
    REAL (KIND=r8), INTENT(IN) :: gsh(Ibmax,Kmax,Jbmax)
    REAL (KIND=r8), INTENT(OUT) :: gss(Ibmax,Kmax,Jbmax)
    REAL (KIND=r8), INTENT(IN) :: psmb(Ibmax,Jbmax)
    REAL (KIND=r8), INTENT(OUT) :: tg(Ibmax,Lmax,Jbmax)
    REAL (KIND=r8), INTENT(OUT) :: rg(Ibmax,Lmax,Jbmax)
    REAL (KIND=r8), INTENT(OUT) :: rq(Ibmax,Lmax,Jbmax)
    REAL (KIND=r8), INTENT(IN) :: pmand (Lmax)
    REAL (KIND=r8), INTENT(IN) :: alnpmd(Lmax)

    INTEGER :: i, j, k
    INTEGER :: kp
    INTEGER :: ks(Ibmax,Lmax)

    REAL (KIND=r8) :: ddp(Ibmax,Lmax)
    REAL (KIND=r8) :: pim(Ibmax,Kmax)
    REAL (KIND=r8) :: alnpm(Ibmax,Kmax)
    REAL (KIND=r8) :: alnp(Ibmax,Kmax+1)
    REAL (KIND=r8) :: pi(Ibmax,Kmax+1)
    REAL (KIND=r8) :: p(Ibmax,Kmax+1)
    REAL (KIND=r8) :: rdelp(Ibmax,Kmax-1)
    REAL (KIND=r8) :: brh(Ibmax,Kmax-1)
    REAL (KIND=r8) :: bqq(Ibmax,Kmax-1)
    REAL (KIND=r8) :: bt(Ibmax,Kmax-1)

    REAL (KIND=r8) :: es
    REAL (KIND=r8) :: ee
    REAL (KIND=r8) :: wtl

    wtl=LOG(Po)

    ! compute absolute temperature t=tv/(1.0_r8+0.61_r8*q)

    DO j=1,Jbmax
       DO k=1,Kmax
          DO i=1,Ibmaxperjb(j)
             tm(i,k,j)=gts(i,k,j)/(1.0_r8+CTv*gsh(i,k,j))
          END DO
       END DO

    ! compute pressure on each sigma surface 

       DO k=1,Kmax+1
          DO i=1,Ibmaxperjb(j)
            ! transform from Pa to mbar
             p(i,k) = 1.e-2_r8 * a_hybr(k) + b_hybr(k) * psmb(i,j)
          END DO
       END DO

    ! vapour pressure form specific humidity

       DO k=1,Kmax
          DO i=1,Ibmaxperjb(j)
             gss(i,k,j)=MAX(gsh(i,k,j),SHmin)
             es=Eo*EXP(a*(tm(i,k,j)-To)/(tm(i,k,j)-b))
             ee=0.5_r8*(p(i,k)+p(i,k+1))*gss(i,k,j)/(eps+eps1*gss(i,k,j))
             gss(i,k,j)=ee/es
             gss(i,k,j)=MIN(gss(i,k,j),1.0_r8)
          END DO
       END DO

    ! set the top of the model at Pt (See Constants).

       p(:,Kmax+1)=Pt

    ! log of pressure and exner function on each sigma surface

       DO k=1,Kmax+1
          DO i=1,Ibmaxperjb(j)
             alnp(i,k)=LOG(MAX(p(i,k),0.00000001_r8))
             pi(i,k)=EXP(RdByCp*(alnp(i,k)-wtl))
          END DO
       END DO

!       IF(myid ==0)THEN
!          WRITE (UNIT=nfprt, FMT=*) '  exp(alnpmd(k)),  alnpmd(k)'
!          DO k=1,Lmax
!             WRITE (UNIT=nfprt, FMT=*) k,  exp(alnpmd(k)),  alnpmd(k)
!          END DO
!
!          WRITE (UNIT=nfprt, FMT=*)   'exp(alnp(k)),  alnp(k)'
!          DO k=1,kmax
!             WRITE (UNIT=nfprt, FMT=*) k,  exp(alnp(1,k)),  alnp(1,k),psmb(1,j)
!          END DO
!        END IF

    ! mean value of the exner function, log of pim and
    ! the log of the mean value of pressure in each layer

       DO k=1,Kmax
          DO i=1,Ibmaxperjb(j)
             pim(i,k)=(pi(i,k+1)-pi(i,k))/ &
                  (RdByCp*(alnp(i,k+1)-alnp(i,k)))
             alnpm(i,k)=wtl+CpByRd*LOG(pim(i,k))
          END DO
       END DO

    ! return interpolated t absolute. ts replaced by tm

       DO k=1,Kmax-1
          DO i=1,Ibmaxperjb(j)
             rdelp(i,k)=1.0_r8/(alnpm(i,k+1)-alnpm(i,k))
             bt(i,k)=(tm(i,k+1,j)-tm(i,k,j))*rdelp(i,k)
          END DO
       END DO
       DO k=1,Kmax-1
          DO i=1,Ibmaxperjb(j)
             brh(i,k)=(gss(i,k+1,j)-gss(i,k,j))*rdelp(i,k)
             bqq(i,k)=(gsh(i,k+1,j)-gsh(i,k,j))*rdelp(i,k)
          END DO
       END DO

    ! heights to pressure surfaces  -  hydrostatic interpolation
    ! winds and temps by linear interpolation with ln(p)
    ! search for middle of sigma layer above kp

       ks(:,:)=Kmax
       DO kp=1,Lmax
          DO k=1,Kmax-1
             DO i=1,Ibmaxperjb(j)
                IF ((ks(i,kp) == Kmax) .AND. (alnpmd(kp) > alnpm(i,k))) THEN
                   ks(i,kp)=k
                END IF
             END DO
          END DO
       END DO

    ! find values and slopes for upward or downward extrapolation
    ! as well as for interpolation away from tropopause
    ! ts to tm for absolute t

       DO kp=1,Lmax
          DO i=1,Ibmaxperjb(j)
             IF (ks(i,kp) /= 1) THEN
                IF (pmand(kp) < PRHcut) THEN
                   rg(i,kp,j)=0.0_r8
                   rq(i,kp,j)=1.0e-12_r8
                ELSE
                   IF (ks(i,kp) <= Kmax) THEN
                      rg(i,kp,j)=gss(i,ks(i,kp),j)+brh(i,ks(i,kp)-1)*(alnpmd(kp)-alnpm(i,ks(i,kp)))
                      rq(i,kp,j)=gsh(i,ks(i,kp),j)+bqq(i,ks(i,kp)-1)*(alnpmd(kp)-alnpm(i,ks(i,kp)))
                   ELSE
                      rg(i,kp,j)=gss(i,Kmax,j)
                      rq(i,kp,j)=gsh(i,Kmax,j)
                   END IF
                END IF
                ddp(i,kp)=alnpmd(kp)-alnpm(i,ks(i,kp))
                tg(i,kp,j)=tm(i,ks(i,kp),j)+bt(i,ks(i,kp)-1)*ddp(i,kp)
             ELSE
                IF(.NOT.ExtrapoAdiabatica)THEN
                   rg(i,kp,j)=gss(i,1,j)
                   rq(i,kp,j)=gsh(i,1,j)
                   tg(i,kp,j)=tm(i,1,j)
                ELSE
                   rg(i,kp,j)=gss(i,1,j)
                   rq(i,kp,j)=gsh(i,1,j)*EXP(0.5_r8*RdByCp*ABS(alnpm(i,ks(i,kp))-alnpmd(kp)))
                   tg(i,kp,j)=tm (i,1,j)*EXP(0.5_r8*RdByCp*ABS(alnpm(i,ks(i,kp))-alnpmd(kp)))  
                END IF
             END IF
             rg(i,kp,j)=MAX(rg(i,kp,j),RHmin)
             rg(i,kp,j)=MIN(rg(i,kp,j),RHmax)
          END DO
       END DO
    END DO

  END SUBROUTINE sigtop


  SUBROUTINE gavint (nlevs, nlevr, gausin, gauout, psmb, pmand)

    IMPLICIT NONE

    ! vertical interpolation of gaussian grid fields 

    INTEGER, INTENT(IN) :: nlevs
    INTEGER, INTENT(IN) :: nlevr

    REAL (KIND=r8), INTENT(INOUT) :: gausin(Ibmax,Kmax,Jbmax)
    REAL (KIND=r8), INTENT(OUT) :: gauout(Ibmax,nlevr,Jbmax)
    REAL (KIND=r8), INTENT(IN) :: psmb(Ibmax,Jbmax)
    REAL (KIND=r8), INTENT(IN) :: pmand(Lmax)

    INTEGER :: k
    INTEGER :: j
    INTEGER :: i
    INTEGER :: lat
    INTEGER :: int

    REAL (KIND=r8) :: deltap
    REAL (KIND=r8) :: df
    REAL (KIND=r8) :: dp
    REAL (KIND=r8) :: work(Ibmax,Lmax)
    REAL (KIND=r8) :: pin(Ibmax,Kmax)

    LOGICAL :: Above(Ibmax,Lmax)

    gauout=0.0_r8
    IF (nlevs > 1 .AND. nlevs < Kmax) THEN
       DO j=1,Jbmax
          DO k=nlevs+1,Kmax
             DO i=1,Ibmaxperjb(j)
                gausin(i,k,j)=0.0_r8
             END DO
          END DO
       END DO
    END IF

    DO lat=1,Jbmax

       ! single level inputs

       IF (nlevs ==  1) THEN
          DO i=1,Ibmaxperjb(lat)
             work(i,1)=gausin(i,1,lat)
          END DO

          ! vertical interpolation

       ELSE
          DO k=1,Kmax
             DO i=1,Ibmaxperjb(lat)
                ! transform from Pa to mbar
                pin(i,k)=( 1.e-2_r8 * (0.5_r8*(a_hybr(k)+a_hybr(k+1))) + &
                                (0.5_r8*(b_hybr(k)+b_hybr(k+1)))*psmb(i,lat))
             END DO
          END DO

          Above = .FALSE.
          DO k=1,Lmax
             DO i=1,Ibmaxperjb(lat)
                IF (pmand(k) > pin(i,1)) THEN
                   work(i,k)=gausin(i,1,lat)
                   Above(i,k)=.TRUE.
                END IF
             END DO
          END DO

          DO k=1,Lmax
             DO i=1,Ibmaxperjb(lat)
                IF (.NOT. Above(i,k) .AND. pmand(k) <= pin(i,Kmax)) THEN
                   work(i,k)=gausin(i,Kmax,lat)
                   Above(i,k)=.TRUE.
                END IF
             END DO
          END DO

          DO k=1,Lmax
             DO int=1,Kmax-1
                DO i=1,Ibmaxperjb(lat)
                   IF (.NOT. Above(i,k) .AND. &
                        pin(i,int)  >= pmand(k) .AND. &
                        pin(i,int+1) < pmand(k)) THEN

                      ! interpolation linear in p

                      deltap=pmand(k)-pin(i,int)
                      df=gausin(i,int+1,lat)-gausin(i,int,lat)
                      dp=pin(i,int+1)-pin(i,int)
                      work(i,k)=gausin(i,int,lat)+(df/dp)*deltap
                   END IF
                END DO
             END DO
          END DO
       END IF

       DO k=1,nlevr
          DO i=1,Ibmaxperjb(lat)
             gauout(i,k,lat)=work(i,k)
          END DO
       END DO

    END DO

  END SUBROUTINE gavint


END MODULE SigmaToPressure
