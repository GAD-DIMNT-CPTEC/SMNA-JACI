!
!  $Author: pkubota $
!  $Date: 2006/10/30 18:37:46 $
!  $Revision: 1.2 $
!
MODULE GaussPressure

  USE Constants, ONLY : r8, nferr,ExtrapoAdiabatica

  USE Sizes, ONLY : Ibmax, Jbmax, Lmax, ibmaxperjb

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: getsh, getslp, getth, lowtmp


  LOGICAL :: first_getslp=.TRUE.
  INTEGER :: l5=0
  INTEGER :: l10=0

  LOGICAL :: first_lowtmp=.TRUE.
  INTEGER :: mlvl, nlvlm, nlvlp, nlvlpx, nlvlpp
  REAL (KIND=r8), ALLOCATABLE :: at (:,:)
  REAL (KIND=r8), ALLOCATABLE :: ata(:,:)
  REAL (KIND=r8), ALLOCATABLE :: dlg(:)
  REAL (KIND=r8), ALLOCATABLE :: ac (:)
  REAL (KIND=r8), ALLOCATABLE :: bc (:)

CONTAINS


  SUBROUTINE getsh (sg, rg,rq, pmand, tg, psmb)

    USE Constants, ONLY : a, b, To, Eo, Eps, Eps1, CTv, SHmin

    IMPLICIT NONE

    REAL (KIND=r8), INTENT(INOUT) :: sg(Ibmax,Lmax,Jbmax) ! tv on input, q on output
    REAL (KIND=r8), INTENT(INOUT) :: rg(Ibmax,Lmax,Jbmax) ! r.h. on input
    REAL (KIND=r8), INTENT(IN   ) :: rq(Ibmax,Lmax,Jbmax)
    REAL (KIND=r8), INTENT(IN   ) :: pmand(Lmax)
    REAL (KIND=r8), INTENT(INOUT) :: tg(Ibmax,Lmax,Jbmax)
    REAL (KIND=r8), INTENT(IN   ) :: psmb(Ibmax,Jbmax)

    REAL (KIND=r8), PARAMETER :: Epsl=1.0E-7_r8

    INTEGER :: k, j, i, it
    REAL (KIND=r8) :: es, ee, q1, ue, pq, ff, df, des
    REAL (KIND=r8) :: Q00(Ibmax,Lmax)
    REAL (KIND=r8) :: T00(Ibmax,Lmax)
    LOGICAL :: Above(Ibmax,Lmax)
    LOGICAL :: NotDone(Ibmax,Lmax)

    ! Above ground or not

    DO j=1,Jbmax
       DO k=1,Lmax
          DO i=1,Ibmaxperjb(j)
             Above(i,k) = pmand(k) <= psmb(i,j)
          END DO
       END DO
       NotDone = .NOT. Above

    ! Above ground case

       DO k=1,Lmax
          DO i=1,Ibmaxperjb(j)
             IF(Above(i,k))THEN
                q00(i,k) = rq(i,k,j)
                rg(i,k,j)=MIN(MAX(rg(i,k,j),0.0_r8),1.0_r8)
                es=Eo*EXP(a*(tg(i,k,j)-To)/(tg(i,k,j)-b))
                ee=rg(i,k,j)*es
                sg(i,k,j)=Eps*ee/(pmand(k)-Eps1*ee)
             ELSE
                q00(i,k) = 0.1_r8
             END IF
          END DO
       END DO

    ! below ground case
    ! sg contains tvg from postg4

      ! q0=0.1_r8
       DO it=1,200
          DO k=1,Lmax
             DO i=1,Ibmaxperjb(j)
                IF(NotDone(i,k) )THEN
                   ee=Eo*EXP(a*&
                        (sg(i,k,j)-To*(1.0_r8+CTv*q00(i,k)))/&
                        (sg(i,k,j)- b*(1.0_r8+Ctv*q00(i,k))))
                   ue=rg(i,k,j)*ee
                   pq=(pmand(k)-Eps1*ue)
                   ff=q00(i,k)-Eps*ue/pq
                   des=(-a*Ctv*To)/&
                        (sg(i,k,j)-b*(1.0_r8+CTv*q00(i,k)))+b*CTv/(&
                        (sg(i,k,j)-b*(1.0_r8+CTv*q00(i,k)))*     &
                        (sg(i,k,j)-b*(1.0_r8+CTv*q00(i,k))))
                   df=1.0_r8-Eps*pmand(k)*rg(i,k,j)*des/(pq*pq)
                   q1=ff/df
                   q00(i,k)=q00(i,k)-q1
                   NotDone(i,k) = ABS(q1) > Epsl
                END IF
             END DO
          END DO
       END DO

       DO k=1,Lmax
          DO i=1,Ibmaxperjb(j)
             IF (NotDone(i,k)) THEN
                WRITE (UNIT=nferr, FMT=*)' q0=',q00(i,k),' q1=',q1,' i=',i,' j=',j,' k=',k
                WRITE (UNIT=nferr, FMT=*)' pm=',pmand(k),' ts=',sg(i,k,j),' rs=',rg(i,k,j)
                WRITE (UNIT=nferr, FMT=*)' ***   Disaster in getsh   ***'
                STOP ' ***   Disaster in getsh   ***'
             END IF
          END DO
       END DO

    ! set specific humdity and temperature to values 
    ! consistent with relative humidity and virtual temperature

       DO k=1,Lmax
          DO i=1,Ibmaxperjb(j)
             IF(.NOT. Above(i,k))THEN
                q00(i,k)=MAX(rq(i,k,j),q00(i,k))
                t00(i,k)=sg(i,k,j)/(1.0_r8+CTv*q00(i,k))
                tg(i,k,j)=Max(t00(i,k),tg(i,k,j))
                sg(i,k,j)=q00(i,k)
                es=Eo*EXP(a*(tg(i,k,j)-To)/(tg(i,k,j)-b))
                ee=pmand(k)*q00(i,k)/(Eps+Eps1*q00(i,k))
                rg(i,k,j)=MIN(MAX(ee/es,0.0_r8),1.0_r8)
             END IF
             sg(i,k,j)=MAX(sg(i,k,j),SHmin)
          END DO
       END DO
    END DO

  END SUBROUTINE getsh


  SUBROUTINE getslp (slpg, zg, pmand)

    USE Constants, ONLY : Po, P5

    IMPLICIT NONE

    ! calculates sea level pressure from
    ! empirical nmc formula based upon 1000 mb height
    ! and 1000 mb - 500 mb thickness.

    REAL (KIND=r8), INTENT(OUT) :: slpg(Ibmax,Jbmax)
    REAL (KIND=r8), INTENT(IN ) :: zg(Ibmax,Lmax,Jbmax)
    REAL (KIND=r8), INTENT(IN ) :: pmand(Lmax)

    REAL (KIND=r8), PARAMETER :: dp=0.01_r8
    REAL (KIND=r8), PARAMETER :: ct=1.5422885_r8
!    REAL (KIND=r8), PARAMETER :: ct=2.5422885_r8
!    REAL (KIND=r8), PARAMETER :: ct=3.5422885_r8
!    REAL (KIND=r8), PARAMETER :: ct=4.5422885_r8

!    LOGICAL, SAVE :: first_getslp=.TRUE.
!    INTEGER, SAVE :: l5=0
!    INTEGER, SAVE :: l10=0

    INTEGER :: l, j, i

    IF (first_getslp) THEN

       DO l=1,Lmax
          IF (ABS(pmand(l)-P5) < dp) THEN
             l5=l
             EXIT
          END IF
       END DO

       IF (l5 == 0) THEN
          WRITE (UNIT=nferr, FMT='(A)') &
               ' Unable to Find 500 Mb Level in Mandatory Pressure Array'
          STOP 20
       END IF

       DO l=1,Lmax
          IF (ABS(pmand(l)-Po) < dp) THEN
             l10=l
             EXIT
          END IF
       END DO

       IF (l10 == 0) THEN
          WRITE (UNIT=nferr, FMT='(A)') &
               ' Unable to Find 1000 Mb Level in Mandatory Pressure Array'
          STOP 40
       END IF

       first_getslp=.FALSE.
    END IF

    DO j=1,Jbmax
       DO i=1,Ibmaxperjb(j)
          slpg(i,j)=Po*(EXP(zg(i,l10,j)/(ct*(zg(i,l5,j)-zg(i,l10,j))))-1.0_r8)
       END DO
    END DO

  END SUBROUTINE getslp


  SUBROUTINE getth (thg, pmand, tg)

    USE Constants, ONLY : Po, RdByCp

    IMPLICIT NONE

    REAL (KIND=r8), INTENT(OUT) :: thg(Ibmax,Lmax,Jbmax)
    REAL (KIND=r8), INTENT(IN ) :: tg(Ibmax,Lmax,Jbmax)
    REAL (KIND=r8), INTENT(IN ) :: pmand(Lmax)

    INTEGER :: k, j, i
    REAL (KIND=r8) :: wtl, exner

    wtl=LOG(Po)
    DO k=1,Lmax
       exner=RdByCp*(wtl-LOG(pmand(k)))
       exner=EXP(exner)
       DO j=1,Jbmax
          DO i=1,Ibmaxperjb(j)
             thg(i,k,j)=tg(i,k,j)*exner
          END DO
       END DO
    END DO

  END SUBROUTINE getth


  SUBROUTINE lowtmp (zsv, tsv, pmand)

    USE Constants, ONLY : GravByRd

    IMPLICIT NONE

    !     author : maj t flattery                 date: long ago
    ! 
    !     abstract:
    !     given as input a set of heights at the mandatory pressures,
    !     constructs a set of temperatures at the mandatory levels.
    !     the method involves the construction of thickness
    !     temperatures for the layers bounded by the mandatory
    !     pressure levels.  using these thickness temperatures
    !     a set of equations is developed to serve as constraints
    !     on the level temperatures.  the solution of the set of
    !     equations is obtained by matrix methods which simulate
    !     a least square approximation to the level temperatures
    !     that minimize the error with which the constraining
    !     equations are satisfied.
    ! 
    !     input:
    !     zsv(Lmax)   : heights at mandatory pressures levels (meters)
    !     pmand(Lmax) : mandatory pressures levels (mb)
    ! 
    !     output:
    !     tsv(Lmax) : calculated temperatures at mandatory levels
    ! 
    !     subroutines called: iminv

    REAL (KIND=r8), INTENT(IN ) :: zsv(Ibmax,Lmax,Jbmax)
    REAL (KIND=r8), INTENT(OUT) :: tsv(Ibmax,Lmax,Jbmax)
    REAL (KIND=r8), INTENT(IN ) :: pmand(Lmax)

    INTEGER :: i, j, k, l, l1
    INTEGER :: ir(Lmax)
    INTEGER :: jc(Lmax)

    REAL (KIND=r8) :: det
    REAL (KIND=r8) :: tsvxix(Ibmax,Lmax)
    REAL (KIND=r8) :: aunk(Ibmax,Lmax)
    REAL (KIND=r8) :: rpl(Lmax)
    REAL (KIND=r8) :: unk(Ibmax,2*Lmax)
    REAL (KIND=r8) :: tbar(Ibmax,Lmax)

!    LOGICAL, SAVE :: first_lowtmp=.TRUE.
!    INTEGER, SAVE :: mlvl, nlvlm, nlvlp, nlvlpx, nlvlpp

!    REAL (KIND=r8), ALLOCATABLE, SAVE :: at (:,:)
!    REAL (KIND=r8), ALLOCATABLE, SAVE :: ata(:,:)
!    REAL (KIND=r8), ALLOCATABLE, SAVE :: dlg(:)
!    REAL (KIND=r8), ALLOCATABLE, SAVE :: ac (:)
!    REAL (KIND=r8), ALLOCATABLE, SAVE :: bc (:)

    IF (Lmax < 3) THEN
      WRITE (UNIT=nferr, FMT='(/,A,/,A,/)') &
           ' Lmax = nmand must be at least 3 to use lowtmp ', &
           ' See NameList for further details'
      STOP ' Lmax = nmand < 3 in lowtmp : see NameList '
    END IF

    IF (first_lowtmp) THEN

       ALLOCATE (at(Lmax,2*Lmax))
       ALLOCATE (ata(Lmax,Lmax))
       ALLOCATE (dlg(Lmax))
       ALLOCATE (ac(Lmax))
       ALLOCATE (bc(Lmax))
       mlvl=Lmax
       rpl=LOG(pmand)
       DO i=2,mlvl
          dlg(i)=rpl(i)-rpl(i-1)
          dlg(i)=1.0_r8/dlg(i)
       END DO
       nlvlm=mlvl-1
       ac(1)=0.0_r8
       bc(1)=ac(1)
       ac(mlvl)=ac(1)
       bc(mlvl)=ac(1)
       DO k=2,nlvlm
          ac(k)=(rpl(k+1)-rpl(k))/(rpl(k+1)-rpl(k-1))
          bc(k)=(rpl(k)-rpl(k-1))/(rpl(k+1)-rpl(k-1))
       END DO
       nlvlp=nlvlm+mlvl-1
       nlvlpp=nlvlp+1
       at=0.0_r8
       DO i=1,nlvlm
          at(i,i)=0.5_r8
          at(i+1,i)=at(i,i)
       END DO
       DO i=1,mlvl
          at(i,i+nlvlm)=1.0_r8
       END DO
       DO i=1,mlvl
          DO j=1,mlvl
             ata(i,j)=0.0_r8
             DO k=1,nlvlpp
                ata(i,j)=ata(i,j)+at(i,k)*at(j,k)
             END DO
          END DO
       END DO
       CALL iminv (ata, Lmax, det, ir, jc)
       nlvlpx=nlvlp-1
       first_lowtmp=.FALSE.
    END IF

    ! deltz in m.

    DO j=1,Jbmax
       DO k=2,mlvl
          DO i=1,Ibmaxperjb(j)
             tbar(i,k)=(zsv(i,k,j)-zsv(i,k-1,j))*dlg(k)
          END DO
       END DO
       DO l=1,nlvlm
          DO i=1,Ibmaxperjb(j)
             unk(i,l)=tbar(i,l+1)
          END DO
       END DO
       DO l=mlvl,nlvlpx
          DO i=1,Ibmaxperjb(j)
             unk(i,l+1)=ac(l-mlvl+2)*tbar(i,l-mlvl+2)+ &
                             bc(l-mlvl+2)*tbar(i,l-mlvl+3)
          END DO
       END DO
       DO i=1,Ibmaxperjb(j)
          unk(i,mlvl)=-unk(i,mlvl+1)+2.0_r8*tbar(i,2)
          unk(i,nlvlpp)=-unk(i,nlvlp)+2.0_r8*tbar(i,mlvl)
       END DO
       aunk=0.0_r8
       DO l=1,mlvl
          DO l1=1,nlvlpp
             DO i=1,Ibmaxperjb(j)
                aunk(i,l)=aunk(i,l)+at(l,l1)*unk(i,l1)
             END DO
          END DO
       END DO
       tsvxix=0.0_r8
       DO l=1,mlvl
          DO l1=1,mlvl
             DO i=1,Ibmaxperjb(j)
                tsvxix(i,l)=tsvxix(i,l)+ata(l,l1)*aunk(i,l1)
             END DO
          END DO
       END DO
       tsv(:,:,j)=tsvxix*GravByRd
    END DO

    ! tsv in deg celsius.


  END SUBROUTINE lowtmp


  SUBROUTINE iminv (a, n, d, l, m)

    IMPLICIT NONE

    !        purpose
    !           invert a matrix
    ! 
    !        description of parameters
    !           a - input matrix, destroyed in computation and replaced by
    !               resultant inverse.
    !           n - order of matrix a
    !           d - resultant determinant
    !           l - work vector of length n
    !           m - work vector of length n
    ! 
    !        remarks
    !           matrix a must be a general matrix
    ! 
    !        subroutines and function subprograms required
    !           none
    ! 
    !        method
    !           the standard gauss-jordan method is used. the determinant
    !           is also calculated. a determinant of zero indicates that
    !           the matrix is singular.
    ! 
    !     additional remarks:
    ! 
    !        if a double precision version of this routine is desired, the
    !        c in column 1 should be removed from the double precision
    !        statement which follows.
    ! 
    !     double precision a, d, biga, hold
    ! 
    !        the c must also be removed from double precision statements
    !        appearing in other routines used in conjunction with this
    !        routine.
    ! 
    !        the double precision version of this subroutine must also
    !        contain double precision fortran functions.  abs in statement
    !        10 must be changed to dabs.

    REAL (KIND=r8), INTENT(INOUT) :: a(*)
    INTEGER, INTENT(IN) :: n
    REAL (KIND=r8), INTENT(OUT) :: d
    INTEGER, INTENT(OUT) :: l(*)
    INTEGER, INTENT(OUT) :: m(*)

    INTEGER :: nk, k, kk, j, iz, i, ij, ki
    INTEGER :: ji, jp, jk, ik, kj, jq, jr
    REAL (KIND=r8) :: biga, hold

    ! search for largest element

    d=1.0_r8
    nk=-n
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
             IF(ABS(biga) < ABS(a(ij))) THEN
                biga=a(ij)
                l(k)=i
                m(k)=j
             END IF
          END DO
       END DO

       ! interchange rows

       j=l(k)
       IF (j > k) THEN
          ki=k-n
          DO i=1,n
             ki=ki+n
             hold=-a(ki)
             ji=ki-k+j
             a(ki)=a(ji)
             a(ji) =hold
          END DO
       END IF

       ! interchange columns

       i=m(k)
       IF (i > k) THEN
          jp=n*(i-1)
          DO j=1,n
             jk=nk+j
             ji=jp+j
             hold=-a(jk)
             a(jk)=a(ji)
             a(ji) =hold
          END DO
       END IF

       ! divide column by minus pivot (pivot in biga)

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

       ! reduce matrix

       DO i=1,n
          ik=nk+i
          ij=i-n
          DO j=1,n
             ij=ij+n
             IF( (i /= k) .AND. (j /= k) ) THEN
                kj=ij-i+k
                a(ij)=a(ik)*a(kj)+a(ij)
             END IF
          END DO
       END DO

       ! divide row by pivot

       kj=k-n
       DO j=1,n
          kj=kj+n
          IF (j /= k) THEN
             a(kj)=a(kj)/biga
          END IF
       END DO

       ! product of pivots

       d=d*biga

       ! replace pivot by reciprocal

       a(kk)=1.0_r8/biga
    END DO

    ! final row and column interchange

    DO k=n-1,1,-1
       i=l(k)
       IF (i > k) THEN
          jq=n*(k-1)
          jr=n*(i-1)
          DO j=1,n
             jk=jq+j
             hold=a(jk)
             ji=jr+j
             a(jk)=-a(ji)
             a(ji) =hold
          END DO
       END IF
       j=m(k)
       IF (j > k) THEN
          ki=k-n
          DO i=1,n
             ki=ki+n
             hold=a(ki)
             ji=ki-k+j
             a(ki)=-a(ji)
             a(ji) =hold
          END DO
       END IF
    END DO

  END SUBROUTINE iminv


END MODULE GaussPressure
