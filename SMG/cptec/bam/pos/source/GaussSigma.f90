!
!  $Author: pkubota $
!  $Date: 2006/10/30 18:38:08 $
!  $Revision: 1.2 $
!
MODULE GaussSigma

  USE Constants, ONLY : r8, nfprt
  USE Parallelism, ONLY : myid
  USE Sizes, ONLY : Kmax, ibmax, jbmax, ibmaxperjb

  IMPLICIT NONE

  PRIVATE

  REAL (KIND=r8), ALLOCATABLE, DIMENSION(:), PUBLIC :: a_hybr, b_hybr, c_hybr
  REAL (KIND=r8), ALLOCATABLE, DIMENSION(:), PUBLIC :: a_hybr_cb, c_hybr_cb

  PUBLIC :: CreateHybrCoor, Omegas, pWater


CONTAINS


  SUBROUTINE CreateHybrCoor()

    IMPLICIT NONE

    ALLOCATE (a_hybr(Kmax+1), b_hybr(kMax+1), c_hybr(kMax))
    ALLOCATE (a_hybr_cb(Kmax+1), c_hybr_cb(kMax))

  END SUBROUTINE CreateHybrCoor


  SUBROUTINE Omegas (dphi, dlam, ug, vg, dg, rcl, omg, psmb)

    IMPLICIT NONE

    REAL (KIND=r8), INTENT(IN   ) :: dphi(Ibmax,Jbmax)
    REAL (KIND=r8), INTENT(IN   ) :: dlam(Ibmax,Jbmax)
    REAL (KIND=r8), INTENT(IN   ) :: ug  (Ibmax,Kmax,Jbmax)
    REAL (KIND=r8), INTENT(IN   ) :: vg  (Ibmax,Kmax,Jbmax)
    REAL (KIND=r8), INTENT(IN   ) :: dg  (Ibmax,Kmax,Jbmax)
    REAL (KIND=r8), INTENT(IN   ) :: rcl (Ibmax,Jbmax)
    REAL (KIND=r8), INTENT(OUT  ) :: omg (Ibmax,Kmax,Jbmax)
    REAL (KIND=r8), INTENT(IN   ) :: psmb(Ibmax,Jbmax)

    INTEGER :: i,j,k

    REAL (KIND=r8) :: phalf (Ibmax,Kmax+1,Jbmax)
    REAL (KIND=r8) :: delp  (Ibmax,Kmax,Jbmax)
    REAL (KIND=r8) :: delb  (Kmax)
    REAL (KIND=r8) :: rpi   (Ibmax,Kmax,Jbmax)
    REAL (KIND=r8) :: alpha (Ibmax,Kmax,Jbmax)
    REAL (KIND=r8) :: adveps(Ibmax,Kmax,Jbmax)
    REAL (KIND=r8) :: psint (Ibmax,Kmax,Jbmax)
    REAL (KIND=r8) :: divint(Ibmax,Kmax,Jbmax)
    REAL (KIND=r8) :: ps    (Ibmax,Jbmax)
    REAL (KIND=r8) :: c_hybr_l (Kmax)
    REAL (KIND=r8) :: c_hybr_Pa_l(kMax)
    ! Compute Omega in Layers ( in Pa / s )

    DO k=1,Kmax
        !delb(k)=b_hybr(k+1)-b_hybr(k)
       delb     (k) = b_hybr(k) - b_hybr(k+1)
       !c_hybr(k)   =  a_hybr   (k+1) * b_hybr(k) - a_hybr   (k) * b_hybr(k+1)

       c_hybr_cb(k) =  a_hybr_cb(k) * b_hybr(k+1) - a_hybr_cb(k+1) * b_hybr(k) 
       c_hybr   (k) = 10.0_r8 * c_hybr_cb(k)

       c_hybr_Pa_l(k) =  a_hybr(k) * b_hybr(k+1) - a_hybr(k+1) * b_hybr(k) 
       c_hybr_l   (k) =  c_hybr_Pa_l(k)
    ENDDO

    DO j=1,Jbmax
       DO i=1,Ibmaxperjb(j)
          ps(i,j)= 100.0_r8*psmb(i,j) 
         !phalf(i,kmax+1) = ps(i)
          phalf(i,1,j) = ps(i,j)
       ENDDO
       DO k=1,kMax   
          DO i=1,Ibmaxperjb(j)
             phalf (i,k+1,j) =  a_hybr(k+1) + b_hybr(k+1) * ps(i,j)
            ! delp(i,k)    = phalf(i,k+1)-phalf(i,k)
             delp  (i,k,j)   = phalf(i,k,j)-phalf(i,k+1,j)
            !rcl(i,j)  AuxGaussColat = 1.0_r16/(SinGaussColat*SinGaussColat)

            ! adveps(i,k) = rcl(i)*(u(i,k) * plam(i) + v(i,k) * pphi(i) )
             adveps(i,k,j)   = rcl(i,j)*(ug(i,k,j)*dlam(i,j)+vg(i,k,j)*dphi(i,j))
          ENDDO
       ENDDO
       !alpha(i,1) = log(2.0_r8)
       k=kMax
       DO i=1,Ibmaxperjb(j)
          alpha(i,k,j) = log(2.0_r8)
          !psint(i,k) = delb(k) * adveps(i,k)
          psint(i,k,j) = delb(k) * adveps(i,k,j)
         ! divint(i,k) = delp(i,k) * div(i,k)
          divint(i,k,j) = delp(i,k,j) * dg(i,k,j)
       ENDDO

       DO k=kMax-1,1,-1
          DO i=1,Ibmaxperjb(j)
             !psint(i,k) = psint(i,k-1) + delb(k) * adveps(i,k)
             psint(i,k,j) = psint(i,k+1,j) + delb(k) * adveps(i,k,j)
             !divint(i,k) = divint(i,k-1) + delp(i,k) * div(i,k)
             divint(i,k,j) = divint(i,k+1,j) + delp(i,k,j) * dg(i,k,j)
             !rpi(i,k) = log (phalf(i,k+1)/phalf(i,k))
             rpi(i,k,j) = log (phalf(i,k,j)/phalf(i,k+1,j))
             !alpha(i,k) = 1.0_r8 - phalf(i,k) * rpi(i,k) / delp(i,k) 
             alpha(i,k,j) = 1.0_r8 - phalf(i,k,j) * rpi(i,k,j) / delp(i,k,j)
             omg(i,k,j) = rpi(i,k,j)*(psint(i,k+1,j)*ps(i,j) + divint(i,k+1,j)) + &
                          alpha(i,k,j) * (delp(i,k,j)*dg(i,k,j)+&
                          adveps(i,k,j)*ps(i,j)*delb(k) ) - ps(i,j) *(delb(k)+ &
                          c_hybr_l(k) * rpi(i,k,j) / delp(i,k,j) ) * adveps(i,k,j)
          ENDDO
       ENDDO
       k = kMax
       DO i=1,Ibmaxperjb(j)
          omg(i,k,j) = alpha(i,k,j) * (delp(i,k,j)*dg(i,k,j) + &
                    adveps(i,k,j)*ps(i,j)*delb(k)) - ps(i,j)*delb(k)*adveps(i,k,j)
       ENDDO
       DO k=1,kMax   
          DO i=1,Ibmaxperjb(j)
             omg(i,k,j) =-0.01*omg(i,k,j)*0.5_r8*(phalf(i,k,j)+phalf(i,k+1,j))/delp(i,k,j)
          ENDDO
       ENDDO
    ENDDO

   
  END SUBROUTINE Omegas


  SUBROUTINE pWater (jjsh, Pw, Psmb)

    USE Constants, ONLY : CvMbPa, Grav

    IMPLICIT NONE

    ! Multiply Matrix jjsh by Vector Del and Scales Results by Psmb

    REAL (KIND=r8), INTENT(IN)  :: jjsh(Ibmax,Kmax,Jbmax)
    REAL (KIND=r8), INTENT(OUT) :: Pw(Ibmax,Jbmax)
    REAL (KIND=r8), INTENT(IN)  :: Psmb(Ibmax,Jbmax)

    INTEGER :: i, k, j, ibLim

    REAL (KIND=r8) :: Fac
    REAL (KIND=r8) :: phalf(Ibmax,Kmax+1)
    REAL (KIND=r8) :: delp(Ibmax,Kmax)

    Fac=CvMbPa/Grav
    Pw=0.0_r8
    DO j=1,Jbmax
       ibLim = Ibmaxperjb(j)
       DO i=1,ibLim
          phalf(i,1) = Psmb(i,j)
       ENDDO
       DO k=1,kMax   
          DO i=1,ibLim
             phalf(i,k+1) = 10._r8 * a_hybr_cb(k+1) + b_hybr(k+1) * Psmb(i,j)
             delp(i,k) = phalf(i,k)-phalf(i,k+1)
             Pw(i,j)=Pw(i,j)+jjsh(i,k,j)*delp(i,k)
          END DO
       END DO
       DO i=1,ibLim
          Pw(i,j)=Pw(i,j)*Fac
       END DO
    END DO

  END SUBROUTINE pWater

END MODULE GaussSigma
