!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
MODULE Legendre

  USE InputParameters, ONLY: r8, nfprt, &
                             Mend  => MendOut,  Mend1 => Mend1Out, &
                             Mend2 => Mend2Out, Mend3 => Mend3Out, &
                             Mnwv0 => Mnwv0Out, Mnwv1 => Mnwv1Out, &
                             Mnwv2 => Mnwv2Out, Mnwv3 => Mnwv3Out, &
                             Imax => ImaxOut, Jmax => JmaxOut, Kmax => KmaxOut, &
                             Imx  => ImxOut,  Jmaxhf => JmaxhfOut, &
                             EMRad, EMRad1, EMRad12, EMRad2

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: InitLegendre, ClsMemLegendre, &
            FourierToSpecCoef, SymAsy, GetConv, GetVort, &
            SpecCoefToFourierScal, SpecCoefToFourierVect

  INTEGER, DIMENSION (:,:), ALLOCATABLE, PUBLIC :: la0, la1

  REAL (KIND=r8), DIMENSION (:), ALLOCATABLE, PUBLIC :: glat, colrad, rcs2, snnp1, eps

  REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE, PUBLIC :: qln, qdln, qlnw, qder, qlnwcs, qlnv

  REAL (KIND=r8), DIMENSION (:), ALLOCATABLE :: wgt, pln, dpln, der, plnwcs

CONTAINS


SUBROUTINE InitLegendre

  IMPLICIT NONE

  INTEGER :: j, l, nn, Mmax, mm, lx, mn

  REAL (KIND=r8) :: rd, sn

  INTEGER, SAVE :: ifp=1

  IF (ifp == 1) THEN
     CALL SetMemLegendre
     ifp=0
  END IF

  l=0
  DO nn=1,Mend1
     Mmax=Mend2-nn
     DO mm=1,Mmax
        l=l+1
        la0(mm,nn)=l
        sn=REAL(mm+nn-2,r8)
        IF (sn /= 0.0_r8) THEN
           sn=-EMRad2/(sn*(sn+1.0_r8))
           snnp1(2*l-1)=sn
           snnp1(2*l)=sn
        ELSE
           snnp1(2*l-1)=0.0_r8
           snnp1(2*l)=0.0_r8
        END IF
     END DO
  END DO
  l=0
  DO mm=1,Mend1
     l=l+1
     la1(mm,1)=l
  END DO
  DO nn=2,Mend2
     Mmax=Mend3-nn
     DO mm=1,Mmax
        l=l+1
        la1(mm,nn)=l
     END DO
  END DO
  CALL epslon
  CALL glats
  rd=45.0_r8/ATAN(1.0_r8)
  DO j=1,Jmaxhf
     glat(j)=90.0_r8-colrad(j)*rd
     glat(Jmax-j+1)=-glat(j)
     CALL pln2 (j)
     l=0
     DO nn=1,Mend1
        Mmax=Mend2-nn
        DO mm=1,Mmax
           l=l+1
           lx=la1(mm,nn)
           qln(2*l-1,j)=pln(lx)
           qln(2*l,j)=pln(lx)
        END DO
     END DO
     DO mn=1,Mnwv2
        qlnw(mn,j)=qln(mn,j)*wgt(j)
     END DO
     DO mn=1,Mnwv1
        qlnv(2*mn-1,j)=pln(mn)
        qlnv(2*mn,j)=pln(mn)
     END DO
     CALL plnder (rcs2(j), wgt(j))
     DO mn=1,Mnwv0
        qdln(2*mn-1,j)=dpln(mn)
        qdln(2*mn,j)=dpln(mn)
        qder(2*mn-1,j)=der(mn)
        qder(2*mn,j)=der(mn)
        qlnwcs(2*mn-1,j)=plnwcs(mn)
        qlnwcs(2*mn,j)=plnwcs(mn)
     END DO
  END DO
  WRITE (UNIT=nfprt, FMT='(/,A)') ' Gaussian Latitudes:'
  WRITE (UNIT=nfprt, FMT='(6F10.3)') glat(1:Jmaxhf)

END SUBROUTINE InitLegendre


SUBROUTINE SetMemLegendre

  IMPLICIT NONE

   ALLOCATE (la0(Mend1,Mend1))
   ALLOCATE (la1(Mend1,Mend2))

   ALLOCATE (glat(Jmax))
   ALLOCATE (colrad(Jmaxhf))
   ALLOCATE (rcs2(Jmaxhf))
   ALLOCATE (snnp1(Mnwv2))
   ALLOCATE (eps(Mnwv1))

   ALLOCATE (qln(Mnwv2,Jmaxhf))
   ALLOCATE (qdln(Mnwv2,Jmaxhf))
   ALLOCATE (qlnw(Mnwv2,Jmaxhf))
   ALLOCATE (qder(Mnwv2,Jmaxhf))
   ALLOCATE (qlnwcs(Mnwv2,Jmaxhf))
   ALLOCATE (qlnv(Mnwv3,Jmaxhf))

   ALLOCATE (wgt(Jmaxhf))
   ALLOCATE (pln(Mnwv1))
   ALLOCATE (dpln(Mnwv0))
   ALLOCATE (der(Mnwv0))
   ALLOCATE (plnwcs(Mnwv0))

END SUBROUTINE SetMemLegendre


SUBROUTINE ClsMemLegendre

  IMPLICIT NONE

   DEALLOCATE (la0)
   DEALLOCATE (la1)

   DEALLOCATE (glat)
   DEALLOCATE (colrad)
   DEALLOCATE (rcs2)
   DEALLOCATE (snnp1)
   DEALLOCATE (eps)

   DEALLOCATE (qln)
   DEALLOCATE (qdln)
   DEALLOCATE (qlnw)
   DEALLOCATE (qder)
   DEALLOCATE (qlnwcs)
   DEALLOCATE (qlnv)

   DEALLOCATE (wgt)
   DEALLOCATE (pln)
   DEALLOCATE (dpln)
   DEALLOCATE (der)
   DEALLOCATE (plnwcs)

END SUBROUTINE ClsMemLegendre


SUBROUTINE epslon
 
  !     epslon : calculates an array which is used in a 
  !              recursion relation to calculate the 
  !              associated legendre functions.
  ! 
  !     argument(dimensions)          description
  ! 
  !     eps(Mnwv1)           output : array used in various
  !                                   calculations involving
  !                                   spherical harmonics.
  !                                   values at the end of each
  !                                   column are stored from Mnwv0+1
  !      Mend                 input : truncation of the triangular
  !                                   spectral domain.
  !      Mnwv1                input : number of wave coefficients

  IMPLICIT NONE
 
  INTEGER :: l, nn, Mmax, mm

  REAL (KIND=r8) :: am, an

  DO l=1,Mend1
     eps(l)=0.0_r8
  END DO
  l=Mend1
  DO nn=2,Mend2
     Mmax=Mend3-nn
     DO mm=1,Mmax
        l=l+1
        am=REAL(mm-1,r8)
        an=REAL(mm+nn-2,r8)
        eps(l)=SQRT((an*an-am*am)/(4.0_r8*an*an-1.0_r8))
     END DO
  END DO

END SUBROUTINE epslon


SUBROUTINE glats
 
  !     glats: calculates gaussian latitudes and 
  !            gaussian weights for use in grid-spectral 
  !            and spectral-grid transforms.
  ! 
  !      glats calls the subroutine poly
  ! 
  !     argument(dimensions)        description
  ! 
  !     colrad(Jmaxhf)      output: co-latitudes for gaussian
  !                                 latitudes in one hemisphere.
  !     wgt(Jmaxhf)         output: gaussian weights.
  !     rcs2(Jmaxhf)        output: 1.0/cos(gaussian latitude)**2
  !     Jmaxhf               input: number of gaussian latitudes
  !                                 in one hemisphere.
 
  IMPLICIT NONE

  INTEGER :: j

  REAL (KIND=r8) :: prec, scale, dradz, rad, drad, p2, p1

  prec=100.0_r8*EPSILON(prec)
  scale=2.0_r8/(REAL(Jmax,r8)*REAL(Jmax,r8))
  dradz=ATAN(1.0_r8)/REAL(Jmax,r8)
  rad=0.0_r8
  DO j=1,Jmaxhf
     drad=dradz
     DO
        CALL poly (Jmax, rad, p2)
        DO
           p1=p2
           rad=rad+drad
           CALL poly (Jmax, rad, p2)
           IF (SIGN(1.0_r8,p1) /= SIGN(1.0_r8,p2)) EXIT
        END DO
        IF (drad <= prec) EXIT
        rad=rad-drad
        drad=drad*0.25_r8
     END DO
     colrad(j)=rad
     CALL poly (Jmax-1, rad, p1)
     wgt(j)=scale*(1.0_r8-COS(rad)*COS(rad))/(p1*p1)
     rcs2(j)=1.0_r8/(SIN(rad)*SIN(rad))
  END DO

END SUBROUTINE glats


SUBROUTINE poly (n, rad, p)

  ! poly : calculates the value of the ordinary legendre function
  !        of given order at a specified latitude.  used to
  !        determine gaussian latitudes.
  !
  !***********************************************************************
  !
  ! poly is called by the subroutine glats.
  !
  ! poly calls no subroutines.
  !
  !***********************************************************************
  !
  ! argument(dimensions)                       description
  !
  !             n                   input : order of the ordinary legendre
  !                                         function whose value is to be
  !                                         calculated. set in routine
  !                                         "glats".
  !            rad                  input : colatitude (in radians) at
  !                                         which the value of the ordinar
  !                                         legendre function is to be
  !                                         calculated. set in routine
  !                                         "glats".
  !             p                  output : value of the ordinary legendre
  !                                         function of order "n" at
  !                                         colatitude "rad".
  !
  !***********************************************************************
  !

  IMPLICIT NONE

  INTEGER, INTENT (IN ) :: n
  REAL    (KIND=r8), INTENT (IN ) :: rad
  REAL    (KIND=r8), INTENT (OUT) :: p

  INTEGER :: i

  REAL (KIND=r8) :: x, y1, y2, y3, g

  x=COS(rad)
  y1=1.0_r8
  y2=x
  DO i=2,n
     g=x*y2
     y3=g-y1+g-(g-y1)/REAL(i,r8)
     y1=y2
     y2=y3
  END DO
  p=y3

END SUBROUTINE poly


SUBROUTINE pln2 (lat)
 
  !      pln2: calculates the associated legendre functions
  !            at one specified latitude.
  ! 
  !     argument(dimensions)        description
  ! 
  !     pln (Mnwv1)         output: values of associated legendre
  !                                 functions at one gaussian
  !                                 latitude specified by the
  !                                 argument "lat".
  ! 
  !      colrad(Jmaxhf)      input: colatitudes of gaussian grid
  !                                 (in radians). calculated
  !                                 in routine "glats".
  !      lat                 input: gaussian latitude index. set
  !                                 by calling routine.
  !      eps                 input: factor that appears in recusio
  !                                 formula of a.l.f.
  !      Mend                input: triangular truncation wave number
  !      Mnwv1               input: number of elements
  !      la1(Mend1,Mend1+1)  input: numbering array of pln1
  ! 
  !      x(Mend1)            local: work area
  !      y(Mend1)            local: work area
 
  IMPLICIT NONE

  INTEGER, INTENT (IN ) :: lat

  INTEGER :: mm, nn, Mmax, lx, ly, lz

  REAL (KIND=r8) :: colr, sinlat, coslat, prod

  INTEGER, SAVE :: ifp=1

  REAL (KIND=r8), SAVE :: rthf

  REAL (KIND=r8), DIMENSION (:), ALLOCATABLE, SAVE :: x, y

  IF (ifp == 1) THEN
     ALLOCATE (x(Mend1))
     ALLOCATE (y(Mend1))
     ifp=0
     DO mm=1,Mend1
        x(mm)=SQRT(2.0_r8*mm+1.0_r8)
        y(mm)=SQRT(1.0_r8+0.5_r8/REAL(mm,r8))
     END DO
     rthf=SQRT(0.5_r8)
  END IF
  colr=colrad(lat)
  sinlat=COS(colr)
  coslat=SIN(colr)
  prod=1.0_r8
  DO mm=1,Mend1
     pln(mm)=rthf*prod
     !     line below should only be used where exponent range is limted
     !     IF (prod < flim) prod=0.0_r8
     prod=prod*coslat*y(mm)
  END DO

  DO mm=1,Mend1
     pln(mm+Mend1)=x(mm)*sinlat*pln(mm)
  END DO
  DO nn=3,Mend2
     Mmax=Mend3-nn
     DO mm=1,Mmax
        lx=la1(mm,nn)
        ly=la1(mm,nn-1)
        lz=la1(mm,nn-2)
        pln(lx)=(sinlat*pln(ly)-eps(ly)*pln(lz))/eps(lx)
     END DO
  END DO

END SUBROUTINE pln2


SUBROUTINE plnder (rcs2l, wgtl)
 
  !     calculates zonal and meridional pseudo-derivatives as
  !     well as laplacians of the associated legendre functions.
  ! 
  !     argument(dimensions)             description
  ! 
  !     pln   (Mnwv1)            input : associated legendre function
  !                                      values at gaussian latitude
  !                             output : pln(l,n)=
  !                                      pln(l,n)*(l+n-2)*(l+n-1)/a**2.
  !                                      ("a" denotes earth's radius)
  !                                      used in calculating the
  !                                      laplacian of global fields
  !                                      in spectral form.
  !     pdln  (Mnwv0)           output : derivatives of
  !                                      associated legendre functions
  !                                      at gaussian latitude
  !     der   (Mnwv0)           output : cosine-weighted derivatives of
  !                                      associated legendre functions
  !                                      at gaussian latitude
  !     plnwcs(Mnwv0)           output : plnwcs(l,n)=
  !                                      pln(l,n)*(l-1)/sin(latitude)**2.
  !     rcs2l                    input : 1.0/sin(latitude)**2 at
  !                                      gaussian latitude
  !     wgtl                     input : gaussian weight, at gaussian
  !                                      latitude
  !     eps   (Mnwv1)            input : array of constants used to
  !                                      calculate "der" from "pln".
  !                                      computed in routine "epslon".
  !     la1(Mend1,Mend2)         input : numbering array for pln
 
  IMPLICIT NONE

  REAL (KIND=r8), INTENT (IN) :: rcs2l, wgtl

  INTEGER :: n, l, nn, Mmax, mm, mn, lm, l0, lp

  REAL (KIND=r8) :: raa, wcsa

  REAL (KIND=r8), DIMENSION (Mnwv1) :: x

  INTEGER, SAVE :: ifp=1

  REAL (KIND=r8), DIMENSION (:), ALLOCATABLE, SAVE :: an

! Compute Pln Derivatives
 
  IF (ifp == 1) THEN
     ALLOCATE (an(Mend2))
     DO n=1,Mend2
        an(n)=REAL(n-1,r8)
     END DO
     ifp=0
  END IF

  raa=wgtl*EMRad12
  wcsa=rcs2l*wgtl*EMRad1
  l=0
  DO mm=1,Mend1
     l=l+1
     x(l)=an(mm)
  END DO
  DO nn=2,Mend2
     Mmax=Mend3-nn
     DO mm=1,Mmax
        l=l+1
        x(l)=an(mm+nn-1)
     END DO
  END DO
  l=Mend1
  DO nn=2,Mend1
     Mmax=Mend2-nn
     DO mm=1,Mmax
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
     Mmax=Mend2-nn
     DO mm=1,Mmax
        l=l+1
        l0=la1(mm,nn)
        plnwcs(l)=an(mm)*pln(l0)
     END DO
  END DO
  DO mn=1,Mnwv0
     plnwcs(mn)=wcsa*plnwcs(mn)
  END DO
  DO nn=1,Mend1
     Mmax=Mend2-nn
     DO mm=1,Mmax
        l0=la1(mm,nn)
        lp=la1(mm,nn+1)
        pln(l0)=x(l0)*x(lp)*raa*pln(l0)
     END DO
  END DO

END SUBROUTINE plnder


SUBROUTINE FourierToSpecCoef (fp, fm, fln, Ldim, lat)
 
  !     calculates spectral representations of
  !     global fields from fourier representations
  !     of symmetric and anti-symmetric portions
  !     of global fields.
  ! 
  !     argument(dimensions)        description
  ! 
  !     fp(Imax+2,Ldim)      input: fourier representation of
  !                                 symmetric portion of a
  !                                 global field at one gaussian
  !                                 latitude.
  !     fm(Imax+2,Ldim)      input: fourier representation of
  !                                 anti-symmetric portion of a
  !                                 global field at one gaussian
  !                                 latitude.
  !     fln(Mnwv2,Ldim)      input: spectral representation of
  !                                 (the laplacian of) a global
  !                                 field. includes contributions
  !                                 from gaussian latitudes up to
  !                                 but not including current
  !                                 iteration of gaussian loop
  !                                 in calling routine.
  !                         output: spectral representation of
  !                                 (the laplacian of) a global
  !                                 field. includes contributions
  !                                 from gaussian latitudes up to
  !                                 and including current
  !                                 iteration of gaussian loop
  !                                 in calling routine.
  !     Ldim                 input: number of vertical layers.
  !     lat                  input: current index of gaussian
  !                                 loop in calling routine.

  IMPLICIT NONE

  INTEGER, INTENT (IN) :: Ldim, lat

  REAL (KIND=r8), DIMENSION (Imx,Ldim),   INTENT (IN ) :: fp, fm
  REAL (KIND=r8), DIMENSION (Mnwv2,Ldim), INTENT (INOUT) :: fln

  INTEGER :: k, l, nn, Mmax, mm, mn

  REAL (KIND=r8), DIMENSION (Mnwv2) :: s

  DO k=1,Ldim
     l=0
     DO nn=1,Mend1
        Mmax=2*(Mend2-nn)
        IF (MOD(nn-1,2) == 0) THEN
           DO mm=1,Mmax
              l=l+1
              s(l)=fp(mm,k)
           END DO
        ELSE
           DO mm=1,Mmax
              l=l+1
              s(l)=fm(mm,k)
           END DO
        END IF
     END DO
     DO mn=1,Mnwv2
        fln(mn,k)=fln(mn,k)+s(mn)*qlnw(mn,lat)
     END DO
  END DO

END SUBROUTINE FourierToSpecCoef


SUBROUTINE SymAsy (a, b, Ldim)

  !     converts the fourier representations of a field at two
  !     parallels at the same latitude in the northern and
  !     southern hemispheres into the fourier representations
  !     of the symmetric and anti-symmetric portions of that
  !     field at the same distance from the equator as the
  !     input latitude circles.
  ! 
  !     argument(dimensions)        description
  ! 
  !     a(Imx,Ldim)          input: fourier representation of one
  !                                 latitude circle of a field
  !                                 from the northern hemisphere
  !                                 at "n" levels in the vertical.
  !                         output: fourier representation of the
  !                                 symmetric portion of a field
  !                                 at the same latitude as the
  !                                 input, at "n" levels in 
  !                                 the vertical.
  !     b(Imx,Ldim)          input: fourier representation of one
  !                                 latitude circle of a field
  !                                 from the southern hemisphere
  !                                 at "n" levels in the vertical.
  !                         output: fourier representation of the
  !                                 anti-symmetric portion of a
  !                                 field at the same latitude as
  !                                 the input, at "n" levels in
  !                                 the vertical.
  !     t(Imx,Ldim)                 temporary storage
  ! 
  !     Ldim                 input: number of layers.
 
  IMPLICIT NONE

  INTEGER, INTENT (IN) :: Ldim

  REAL (KIND=r8), DIMENSION (Imx,Ldim), INTENT (IN OUT) :: a, b

  INTEGER :: i, k

  REAL (KIND=r8), DIMENSION (Imx,Ldim) :: t

  DO k=1,Ldim
     DO i=1,Imax
        t(i,k)=a(i,k)
        a(i,k)=t(i,k)+b(i,k)
        b(i,k)=t(i,k)-b(i,k)
     END DO
  END DO

END SUBROUTINE SymAsy


SUBROUTINE GetConv (am, ap, bm, bp, fln, lat)

  !     calculates the spectral representations of the horizontal
  !     convergence of pseudo-vector fields from the fourier
  !     representations of the symmetric and anti-symmetric
  !     portions of the two individual fields.
  ! 
  !     argument(dimensions)         description
  ! 
  !     am(Imax+2,Kmax)       input: fourier representation of
  !                                  anti-symmetric portion of
  !                                  zonal pseudo-wind field at
  !                                  one gaussian latitude.
  !     ap(Imax+2,Kmax)       input: fourier representation of
  !                                  symmetric portion of zonal
  !                                  pseudo-wind field at one
  !                                  gaussian latitude.
  !     bm(Imax+2,Kmax)       input: fourier representation of
  !                                  anti-symmetric portion of
  !                                  meridional pseudo-wind field
  !                                  at one gaussian latitude.
  !     bp(Imax+2,Kmax)       input: fourier representation of
  !                                  symmetric portion of
  !                                  meridional pseudo-wind field
  !                                  at one gaussian latitude.
  !     fln(Mnwv2,Kmax)       input: spectral representation of
  !                                  the divergence of the global
  !                                  wind field. includes
  !                                  contributions from gaussian
  !                                  latitudes up to but not
  !                                  including current iteration
  !                                  of gaussian loop in calling
  !                                  routine.
  !                          output: spectral representation of
  !                                  the divergence of the global
  !                                  wind field. includes
  !                                  contributions from gaussian
  !                                  latitudes up to and
  !                                  including current iteration
  !                                  of gaussian loop in calling
  !                                  routine.
  !     lat                   input: current index of gaussian
  !                                  loop in calling routine.
 
  IMPLICIT NONE

  INTEGER, INTENT (IN) :: lat

  REAL (KIND=r8), DIMENSION (Imx,Kmax),   INTENT (IN)  :: am, ap, bm, bp
  REAL (KIND=r8), DIMENSION (Mnwv2,Kmax), INTENT (INOUT) :: fln

  INTEGER :: k, l, nn, Mmax, mm, mn

  REAL (KIND=r8), DIMENSION (Mnwv2) :: s 

  DO k=1,Kmax
     l=0
     DO nn=1,Mend1
        Mmax=2*(Mend2-nn)
        IF (MOD(nn-1,2) == 0) THEN
           DO mm=1,Mmax
              l=l+1
              s(l)=bm(mm,k)
           END DO
        ELSE
           DO mm=1,Mmax
              l=l+1
              s(l)=bp(mm,k)
           END DO
        END IF
     END DO
     DO mn=1,Mnwv2
        fln(mn,k)=fln(mn,k)+s(mn)*qder(mn,lat)
     END DO
     l=0
     DO nn=1,Mend1
        Mmax=2*(Mend2-nn)
        IF (MOD(nn-1,2) == 0) THEN
           DO mm=1,Mmax,2
              l=l+1
              s(2*l-1)= ap(mm+1,k)
              s(2*l)=-ap(mm,k)
           END DO
        ELSE
           DO mm=1,Mmax,2
              l=l+1
              s(2*l-1)= am(mm+1,k)
              s(2*l)=-am(mm,k)
           END DO
        END IF
     END DO
     DO mn=1,Mnwv2
        fln(mn,k)=fln(mn,k)+s(mn)*qlnwcs(mn,lat)
     END DO
  END DO

END SUBROUTINE GetConv


SUBROUTINE GetVort (am, ap, bm, bp, fln, lat)

  !     calculates the spectral representations of the horizontal
  !     vorticity   of pseudo-vector fields from the fourier
  !     representations of the symmetric and anti-symmetric
  !     portions of the two individual fields.
  ! 
  !     argument(dimensions)         description
  ! 
  !     am(Imax+2,Kmax)       input: fourier representation of
  !                                  anti-symmetric portion of
  !                                  zonal pseudo-wind field at
  !                                  one gaussian latitude.
  !     ap(Imax+2,Kmax)       input: fourier representation of
  !                                  symmetric portion of zonal
  !                                  pseudo-wind field at one
  !                                  gaussian latitude.
  !     bm(Imax+2,Kmax)       input: fourier representation of
  !                                  anti-symmetric portion of
  !                                  meridional pseudo-wind field
  !                                  at one gaussian latitude.
  !     bp(Imax+2,Kmax)       input: fourier representation of
  !                                  symmetric portion of
  !                                  meridional pseudo-wind field
  !                                  at one gaussian latitude.
  !     fln(Mnwv2,Kmax)       input: spectral representation of
  !                                  the vorticity  of the global
  !                                  wind field. includes
  !                                  contributions from gaussian
  !                                  latitudes up to but not
  !                                  including current iteration
  !                                  of gaussian loop in calling
  !                                  routine.
  !                          output: spectral representation of
  !                                  the vorticity  of the global
  !                                  wind field. includes
  !                                  contributions from gaussian
  !                                  latitudes up to and
  !                                  including current iteration
  !                                  of gaussian loop in calling
  !                                  routine.
  !     lat                   input: current index of gaussian
  !                                  loop in calling routine.
 
  IMPLICIT NONE

  INTEGER, INTENT (IN) :: lat

  REAL (KIND=r8), DIMENSION (Imx,Kmax),   INTENT (IN)  :: am, ap, bm, bp
  REAL (KIND=r8), DIMENSION (Mnwv2,Kmax), INTENT (INOUT) :: fln

  INTEGER :: k, l, nn, Mmax, mm, mn

  REAL (KIND=r8), DIMENSION (Mnwv2) :: s 
 
  DO k=1,Kmax
     l=0
     DO nn=1,Mend1
        Mmax=2*(Mend2-nn)
        IF (MOD(nn-1,2) == 0) THEN
           DO mm=1,Mmax
              l=l+1
              s(l)=am(mm,k)
           END DO
        ELSE
           DO mm=1,Mmax
              l=l+1
              s(l)=ap(mm,k)
           END DO
        END IF
     END DO
     DO mn=1,Mnwv2
        fln(mn,k)=fln(mn,k)+s(mn)*qder(mn,lat)
     END DO
     l=0
     DO nn=1,Mend1
        Mmax=2*(Mend2-nn)
        IF (MOD(nn-1,2) == 0) THEN
           DO mm=1,Mmax,2
              l=l+1
              s(2*l-1)=-bp(mm+1,k)
              s(2*l)=+bp(mm,k)
           END DO
        ELSE
           DO mm=1,Mmax,2
              l=l+1
              s(2*l-1)=-bm(mm+1,k)
              s(2*l)=+bm(mm,k)
           END DO
        END IF
     END DO
     DO mn=1,Mnwv2
        fln(mn,k)=fln(mn,k)+s(mn)*qlnwcs(mn,lat)
     END DO
  END DO

END SUBROUTINE GetVort


SUBROUTINE SpecCoefToFourierScal (fln, ap, am, Ldim, lat)

  !     calculates the fourier representation of a field at a
  !     pair of latitudes symmetrically located about the
  !     equator. the calculation is made using the spectral
  !     representation of the field and the values of the
  !     associated legendre functions at that latitude.
  !     it is designed to triangular truncation only and
  !     for scalar fields.
  !
  !     argument(dimensions)            description
  !
  !     fln(Mnwv2,Ldim)        input: spectral representation of a
  !                                   global field.
  !     ap(Imax+2,Ldim)       output: fourier representation of
  !                                   a global field at the
  !                                   latitude in the northern
  !                                   hemisphere at which the
  !                                   associated legendre functions
  !                                   have been defined.
  !     am(Imax+2,Ldim)       output: fourier representation of
  !                                   a global field at the
  !                                   latitude in the southern
  !                                   hemisphere at which the
  !                                   associated legendre functions
  !                                   have been defined.
  !     Ldim                   input: number of vertical levels.
  !     lat                    input: current index of gaussian
  !                                   loop in calling routine.

  IMPLICIT NONE

  INTEGER, INTENT (IN) :: Ldim, lat

  REAL (KIND=r8), DIMENSION (Mnwv2,Ldim), INTENT (IN)  :: fln
  REAL (KIND=r8), DIMENSION (Imx,Ldim),   INTENT (OUT) :: ap, am

  INTEGER :: l, k, mm, mn, nn, Mmax, Mstr, Mmax1, Mend1d

  REAL (KIND=r8), DIMENSION (Mnwv2) :: s

  DO k=1,Ldim
     Mend1d=2*Mend1
     Mmax1=2*Mend
     l=Mend1d+Mmax1
     DO mn=1,Mnwv2
        s(mn)=qln(mn,lat)*fln(mn,k)
     END DO
     DO nn=3,Mend1
        Mmax=2*(Mend2-nn)
        IF (MOD(nn-1,2) == 0) THEN
           Mstr=0
        ELSE
           Mstr=Mend1d
        END IF
!cdir nodep
        DO mm=1,Mmax
           l=l+1
           s(mm+Mstr)=s(mm+Mstr)+s(l)
        END DO
     END DO
     DO mm=1,Mend1d
        ap(mm,k)=s(mm)
        am(mm,k)=s(mm)
     END DO
     DO mm=1,Mmax1
        ap(mm,k)=ap(mm,k)+s(mm+Mend1d)
        am(mm,k)=am(mm,k)-s(mm+Mend1d)
     END DO
  END DO

END SUBROUTINE SpecCoefToFourierScal


SUBROUTINE SpecCoefToFourierVect (fln, ap, am, Ldim, lat)

  !     calculates the fourier representation of a field at a
  !     pair of latitudes symmetrically located about the
  !     equator. the calculation is made using the spectral
  !     representation of the field and the values of the
  !     associated legendre functions at that latitude.
  !     it is designed to triangular truncation only and
  !     for wind fields.
  !
  !     argument(dimensions)            description
  !
  !     fln(Mnwv3,Ldim)        input: spectral representation of a
  !                                   global field.
  !     ap(Imax+2,Ldim)       output: fourier representation of
  !                                   a global field at the
  !                                   latitude in the northern
  !                                   hemisphere at which the
  !                                   associated legendre functions
  !                                   have been defined.
  !     am(Imax+2,Ldim)       output: fourier representation of
  !                                   a global field at the
  !                                   latitude in the southern
  !                                   hemisphere at which the
  !                                   associated legendre functions
  !                                   have been defined.
  !     Ldim                   input: number of vertical levels.
  !     lat                    input: current index of gaussian
  !                                   loop in calling routine.

  IMPLICIT NONE

  INTEGER, INTENT (IN) :: Ldim, lat

  REAL (KIND=r8), DIMENSION (Mnwv3,Ldim), INTENT (IN)  :: fln
  REAL (KIND=r8), DIMENSION (Imx,Ldim),   INTENT (OUT) :: ap, am

  INTEGER :: l, k, mm, mn, nn, Mmax, Mstr, Mmax1, Mend1d

  REAL (KIND=r8), DIMENSION (Mnwv3) :: s

  DO k=1,Ldim
     Mend1d=2*Mend1
     Mmax1=2*Mend1
     l=Mend1d+Mmax1
     DO mn=1,Mnwv3
        s(mn)=qlnv(mn,lat)*fln(mn,k)
     END DO
     DO nn=3,Mend2
        Mmax=2*(Mend3-nn)
        IF (MOD(nn-1,2) == 0) THEN
           Mstr=0
        ELSE
           Mstr=Mend1d
        END IF
!cdir nodep
        DO mm=1,Mmax
           l=l+1
           s(mm+Mstr)=s(mm+Mstr)+s(l)
        END DO
     END DO
     DO mm=1,Mend1d
        ap(mm,k)=s(mm)
        am(mm,k)=s(mm)
     END DO
     DO mm=1,Mmax1
        ap(mm,k)=ap(mm,k)+s(mm+Mend1d)
        am(mm,k)=am(mm,k)-s(mm+Mend1d)
     END DO
  END DO

END SUBROUTINE SpecCoefToFourierVect


END MODULE Legendre
