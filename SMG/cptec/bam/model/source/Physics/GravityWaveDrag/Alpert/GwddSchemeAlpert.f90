!
!  $Author: pkubota $
!  $Date: 2008/04/09 12:40:01 $
!  $Revision: 1.1 $
!
MODULE GwddSchemeAlpert

  USE Constants, ONLY :     &
       cp,            &
       grav,          &
       gasr,          &
       r8

  USE Options, ONLY :  &
       nferr, nfprt


    IMPLICIT NONE
  SAVE


  PRIVATE
  PUBLIC :: InitGwddSchAlpert
  PUBLIC :: GwddSchAlpert
  REAL(KIND=r8), PARAMETER :: gravi    = 9.81_r8
  REAL(KIND=r8), PARAMETER :: rgas     = 287.0_r8 ! gas constant of dry air        (j/kg/k)
  REAL(KIND=r8), PARAMETER :: agrav    = 1.0_r8/gravi
  REAL(KIND=r8), PARAMETER :: akwnmb   = 2.5e-05_r8      ! (j/kg) 
  REAL(KIND=r8), PARAMETER :: lstar    = 1.0_r8/akwnmb  ! (kg/J) 
  REAL(KIND=r8), PARAMETER :: g        = 1.0_r8
  REAL(KIND=r8), PARAMETER :: gocp     = gravi/1005.0_r8!g/cp

CONTAINS

  SUBROUTINE InitGwddSchAlpert()
    IMPLICIT NONE

  END SUBROUTINE InitGwddSchAlpert
  !----------------------------------------------------------------------
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !----------------------------------------------------------------------
  SUBROUTINE GwddSchAlpert( prsi ,prsl  ,phii ,phil    ,&
       u, v, t, chug, chvg, xdrag, ydrag,xtens,&
       ytens,var, varcut, ncols, kmax)
    ! The parameterization of gravity wave drag (GWD) is based on
    ! simplified theoretical concepts and observational evidence.
    ! The parameterization consists of determines the drag due
    ! to gravity waves at the surface and its vertical variation.
    ! We have adopted the scheme of Alpert et al. (1988) as
    ! implemented by Kirtman et al. (1992). The parameterization
    ! of the drag at the surface is based on the formulation
    ! described by Pierrehumbert (1987) and Helfand et al. (1987).

    !
    !
    ! GwddSchAlpert   :change on gwdd by cptec on 29 july 1994 to improve
    !         vectorization performance on gwdd: dr. j.p. bonatti
    !==========================================================================
    !  ncols......Number of grid points on a gaussian latitude circle
    !  jmax......Number of gaussian latitudes
    !  kmax......Number of sigma levels
    !  latco.....latitude
    !  xdrag.....gravity wave drag surface zonal stress
    !  ydrag.....gravity wave drag surface meridional stress
    !  t.........Temperature
    !  u.........(zonal      velocity)*sin(colat)
    !  v.........(meridional velocity)*sin(colat)
    !  chug......gravity wave drag zonal momentum change
    !  chvg......gravity wave drag meridional momentum change
    !  psfc......surface pressure
    !  si........si(l)=1.0-ci(l).
    !  ci........sigma value at each level.
    !  sl........sigma value at midpoint of
    !                                         each layer : (k=287/1005)
    !
    !                                                                     1
    !                                             +-                   + ---
    !                                             !     k+1         k+1!  k
    !                                             !si(l)   - si(l+1)   !
    !                                     sl(l) = !--------------------!
    !                                             !(k+1) (si(l)-si(l+1)!
    !                                             +-                  -+
    !  del.......sigma spacing for each layer computed in routine "setsig".
    !  varcut....cut off height variance in m**2 for gravity wave drag
    !  var.......Surface height variance
    !  nfprt.....standard print out unit
    !            0 no print, 1 less detail, 2 more detail, 3 most detail
    !==========================================================================

    INTEGER, INTENT(in   ) :: ncols
    INTEGER, INTENT(in   ) :: kmax

    REAL(KIND=r8),    INTENT(inout  ) :: xdrag(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: ydrag(ncols)

    REAL(KIND=r8),    INTENT(in   ) :: prsi   (ncols,kMax+1)  !     prsi     - real, pressure at layer interfaces [Pa]
    REAL(KIND=r8),    INTENT(in   ) :: prsl   (ncols,kMax)    !     prsl     - real, mean layer presure [Pa]
    REAL(KIND=r8),    INTENT(in   ) :: phii   (nCols,kMax+1) !===>  PHIH(K+1)  INPUT GEOPOTENTIAL @ EDGES  IN MKS units (m)
    REAL(KIND=r8),    INTENT(in   ) :: phil   (nCols,kMax)   !===>  PHIL(K) INPUT GEOPOTENTIAL @ LAYERS IN MKS units (m)
    REAL(KIND=r8),    INTENT(in   ) :: t    (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: u    (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: v    (ncols,kmax)

    REAL(KIND=r8),    INTENT(inout) :: chug (ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: chvg (ncols,kmax)

    REAL(KIND=r8),    INTENT(in   ) :: varcut
    REAL(KIND=r8),    INTENT(inout) :: var   (ncols)
    REAL(KIND=r8),    INTENT(out  ) :: xtens (ncols,kmax+1)
    REAL(KIND=r8),    INTENT(out  ) :: ytens (ncols,kmax+1)

    REAL(KIND=r8)    :: psfc (ncols)    !mb
    REAL(KIND=r8)    :: ro    (ncols,kmax)
    REAL(KIND=r8)    :: pp    (ncols,kmax)
    REAL(KIND=r8)    :: tensio(ncols,kmax+1)
    REAL(KIND=r8)    :: dz    (ncols,kmax)
    REAL(KIND=r8)    :: ppp   (ncols,kmax+1)
    REAL(KIND=r8)    :: dragsf(ncols)
    REAL(KIND=r8)    :: nbar  (ncols)
    REAL(KIND=r8)    :: bv    (ncols,kmax)
    REAL(KIND=r8)    :: robar (ncols)
    REAL(KIND=r8)    :: ubar  (ncols)
    REAL(KIND=r8)    :: vbar  (ncols)
    INTEGER :: icrilv(ncols)
    REAL(KIND=r8)    :: speeds(ncols)
    REAL(KIND=r8)    :: ang   (ncols)
    REAL(KIND=r8)    :: coef  (ncols)
    REAL(KIND=r8)    :: DeltaP    (ncols,kmax)



    INTEGER :: k,kk
    INTEGER :: i
    REAL(KIND=r8)    :: vai1
    REAL(KIND=r8)    :: cte
    REAL(KIND=r8)    :: Fr
    REAL(KIND=r8)    :: gstar
    REAL(KIND=r8)    :: roave
    REAL(KIND=r8)    :: velco
    REAL(KIND=r8)    :: fro2
    REAL(KIND=r8)    :: delve2
    REAL(KIND=r8)    :: richsn
    REAL(KIND=r8)    :: crifro
    REAL(KIND=r8)    :: vsqua
    REAL(KIND=r8)    :: dsigma
    REAL(KIND=r8)    :: rbyg
    INTEGER, PARAMETER  :: nsurf=1
    INTEGER          :: nthin(ncols)
    INTEGER          :: nbase(ncols)
    nbase=0
    nthin=0
    xtens=0.0_r8
    ytens=0.0_r8
    DO k=1,kmax
       DO i = 1,ncols
         DeltaP(i,k) = ((prsi(i,k)) - (prsi(i,k+1)))/prsi(i,1)
       END DO
    END DO
    !
    !     surface pressure and constrain the variance
    !
    !PK rbyg=gasr/grav*DeltaP(i,1)*0.5e0_r8
    DO i=1,ncols
       icrilv(i)=0
       !PK psfc(i)=EXP(psfc(i))*10.0_r8 !mb
       psfc(i)=prsi   (i,nsurf)/100.0_r8 
       var(i) =MIN(varcut,var(i))
       var(i) =MAX(0.0_r8,var(i))
       coef(i)=gravi/psfc(i)*0.01_r8
       !
       !     end of input and elementary computations
       !
       !     surface and base layer stress
       !
       !     base layer stress is defined in terms of a vertical ave.
       !
       robar(i)=0.0_r8
       ubar(i)=0.0_r8
       vbar(i)=0.0_r8
       nbar(i)=0.0_r8
       !
       !     compute pressure at every level
       !
       !PK ppp(i,kmax+1)=si(kmax+1)*psfc(i)!mb
       ppp(i,kmax+1)=prsi   (i,kmax+1)/100.0_r8 !mb
       !PK ppp(i,1)=si(1)*psfc(i)
       ppp(i,1)=prsi   (i,nsurf)/100.0_r8 !mb

       !PK pp(i,1)=sl(1)*psfc(i)           !mb
       pp(i,1)=prsl   (i,nsurf)/100.0_r8 !mb

       ro(i,1)=pp(i,1)/(rgas*t(i,1))
       rbyg=gasr/grav*DeltaP(i,1)*0.5e0_r8
       dz(i,1)=MAX((rbyg * t(i,1)),0.5_r8)      
    END DO
    !
    !     compute pressure at every layer, density
    !     at every layer
    !
    DO k=2,kmax
       DO i=1,ncols
          ppp(i,k)=prsi   (i,k)/100.0_r8!si(k)*psfc(i)
          pp (i,k)=prsl   (i,k)/100.0_r8!sl(k)*psfc(i)
          !
          !         100*N        kg*K
          ! ro = ----------- * ---------
          !           m*m        J * K
          !
          !         100*N        kg*K
          ! ro = ----------- * ---------
          !           m*m        N*m* K
          ! 
          !         100          kg
          ! ro = ----------- * ---------
          !           m*m        m

          !         100 kg  
          ! ro = -----------
          !           m*m*  

          !
          !
          ro(i,k)=pp(i,k)/(rgas*t(i,k))! (mb/( (J/(kg*K))*K
          !
          !     compute dz at every level from 2 to kmax
          !
          !       1      1      p
          !dz =-------*----- *---- 
          !       g      ro 
          !
          !       1      Rgas *K     p
          !dz =-------*---------- *---- 
          !       g      P 

          !      Rgas *K  
          !dz =--------- 
          !       g       


          !      J        S*S*K  
          !dz =--------*--------- 
          !      kg*K      m       
          !
          !      N*m       S*S  
          !dz =--------*--------- 
          !      kg       m       

          !      kg*m*m       S*S  
          !dz =--------*--------- 
          !      kg*S*S       m       

          !
          !dz = m
          !              

          dz(i,k)=MAX(0.5_r8*agrav*(1.0_r8/ro(i,k-1) + 1.0_r8/ro(i,k))*(pp(i,k-1)-pp(i,k)),0.01_r8)
          !dz(i,k)=0.5_r8*gasr*(t(i,k-1)+t(i,k))* LOG(pp(i,k-1)/pp(i,k))/grav
       END DO
    END DO

    !
    !     nthin is the number of low layers strapped together
    !
    DO i=1,ncols
       DO kk=1,(kmax+1)
          IF((prsi(i,kk)/prsi(i,nsurf)).LT.0.6667_r8) EXIT
       END DO
       nbase(i)=kk
       IF(nbase(i).LE.1)THEN
          WRITE(UNIT=nfprt,FMT=9976)nbase(i),prsi(i,kk)/prsi(i,nsurf)
          WRITE(UNIT=nferr,FMT=9976)nbase(i),prsi(i,kk)/prsi(i,nsurf)
          STOP 9976
       END IF
       IF(0.6667_r8-prsi(i,nbase(i))/prsi(i,nsurf).GT.prsi(i,nbase(i)-1)/prsi(i,nsurf)-0.6667_r8)nbase(i)=nbase(i)-1
       IF(nbase(i).LE.1)THEN
          WRITE(UNIT=nfprt,FMT=9976)nbase(i),prsi(i,kk)/prsi(i,nsurf)
          WRITE(UNIT=nferr,FMT=9976)nbase(i),prsi(i,kk)/prsi(i,nsurf)
          STOP 9976
       END IF
    END DO

    !
    !     vaisala frequency
    ! (j/kg/k)
    !       g 
    !gocp=-------
    !       cp    

    !       m       kg*K
    !gocp=-------*-------
    !      s*s      J
    !
    !       m       kg*K
    !gocp=-------*-------
    !      s*s      N*m

    !       m       kg*K*s*s
    !gocp=-------*-------
    !      s*s      kg*m

    !       K
    !gocp=-------
    !       m

    !
    !DO k=2,nbase
    DO k=1,kMax
       DO i=1,ncols
          IF(k>=2 .and. k <=nbase(i))THEN
             !
             !           K
             !vai1 = --------
             !           m
             vai1=MAX(0.0_r8,(t(i,k)-t(i,k-1))/dz(i,k)+gocp)

             !        K     m        1 
             !bv = --------*-------*----
             !        m     s*s      K

             !       1 
             !bv = ------
             !      s*s

             bv(i,k)=SQRT(vai1*2.0_r8*gravi/(t(i,k)+t(i,k-1)))
          END IF
       END DO
    END DO
    !
    !     mass weighted veritcal average of density, velocity
    !
    !DO k=1,nbase-1
    DO k=1,kMax
       DO i=1,ncols
          IF(k>=1 .and. k <=nbase(i)-1)THEN
             !
             !            kg*K
             ! robar = ----------
             !           J * K
             !
             !           kg*K
             ! robar = ----------
             !           N*m* K
             ! 
             !             kg
             ! robar =  ----------
             !            N*m

             !            kg  
             ! robar = -----------
             !           kg*m*m  
             robar(i)=robar(i) + ro(i,k)*(ppp(i,k)-ppp(i,k+1))
             ubar(i)=ubar(i)+u(i,k)*(ppp(i,k)-ppp(i,k+1))
             vbar(i)=vbar(i)+v(i,k)*(ppp(i,k)-ppp(i,k+1))
          END IF
       END DO
    END DO
    !
    !     vertical mass weighted average of the brunt-vaisiala freq.
    !
    !DO k=2,nbase
    DO k=1,kMax
       DO i=1,ncols
          IF(k >= 2 .and. k <= nbase(i))THEN
             !          1 
             !bv = ------
             !         s*s

             !          1     1 m*m
             !nbar = ------*--------
             !         s*s   100*N
          
             !          1     1 m*m *s*s
             !nbar = ------*--------------
             !         s*s   100*(kg*m*m)

             !          1
             !nbar =---------
             !        100*kg

             nbar(i)=nbar(i)+bv(i,k)*(pp(i,k-1)-pp(i,k))
          END IF
       END DO
    END DO
    DO i=1,ncols
       cte=1.0_r8/max(ppp(i,1)-ppp(i,nbase(i)),0.000001_r8)
       !            kg           m*m
       ! robar = -----------*-------
       !           kg*m*m         N

       !             1
       ! robar = --------
       !              N

       robar(i)=robar(i)*cte*100.0_r8
       ubar(i)=ubar(i)*cte
       vbar(i)=vbar(i)*cte
       !          1 
       !nbar = ------
       !         s*s

       nbar(i)=nbar(i)/max(pp(i,1)-pp(i,nbase(i)),0.000001_r8)
       !
       !     end vertical average
       !     definition of surface wind vector
       !
       !speeds(i)=SQRT(ubar(i)*ubar(i)+vbar(i)*vbar(i))
       speeds(i)=Max(1.0e-12_r8,SQRT(ubar(i)*ubar(i)+vbar(i)*vbar(i)))
       ang(i)=ATAN2(vbar(i),ubar(i))
       !
       !     stress at the surface level lev=1
       !
       !      Fr = N*h/U
       !     "Fr" "Fr" is the Froude number
       !     " var" "h"  is the amplitude of the orographic perturbation
       !     "nbar" "N"  is the Brunt-Vaisala frequency
       !     "U "   "speeds" is the surface wind speed
       !
       !     tgw = <rho*U*w> It is the momentum flux due to the gravity
       !                     waves averaged over a grid box is written as
       !
       !     Pierrehumbert (1987) argued that convective instability and
       !     wave breaking at the Earths surface are so prevalent that
       !     nonlinear effects become very common and cannot be neglected.
       !     Based on scaling arguments and the results from numerical
       !     experiments, Pierrehumbert (1987) obtained a formula for the
       !     surface wind stress due to gravity waves. He concluded that

       !              (rho*U**3)      (Fr**2)
       !     |tgw| = ------------*----------------
       !                (N*l)       (1 + Fr**2)
       !
       !     "tgw"  "tensio"    is in the direction of the surface wind.
       !


       !
       !     stress at the surface level lev=1
       !
       IF(speeds(i) == 0.0_r8 .OR. nbar(i) == 0.0_r8)THEN
          tensio(i,1)=0.0_r8
       ELSE
          Fr=nbar(i)*SQRT(var(i))/speeds(i)
          !
          !     use non linear weighting function
          !
          gstar=g*Fr*Fr/(Fr*Fr+1.0_r8)!m/s**2

          !          1 
          !nbar = ------
          !         s*s
          !             1
          ! robar = --------
          !             N

          !            1    m*m*m    s*s      kg
          ! tensio = *---* -------*------* -------
          !            N    s*s*s     1      N*m

          !            m*m     kg*s*s
          ! tensio = -------- -------
          !            N*s     kg*m 

          !            m        s
          ! tensio = -------- -------
          !            N      

          !                         m*m*m    
          ! tensio = --------*---* -------*- -------
          !                   kg      s    


          tensio(i,1)=gstar*(robar(i)*speeds(i)*speeds(i)*speeds(i))/ &
               (nbar(i)*lstar)
       ENDIF
       xtens(i,1)=COS(ang(i))*tensio(i,1)
       ytens(i,1)=SIN(ang(i))*tensio(i,1)
       !
       !     save surface values
       !
       dragsf(i)=tensio(i,1)
       xdrag(i)=xtens(i,1)
       ydrag(i)=ytens(i,1)
    END DO
    IF(ANY(speeds == 0.0_r8))WRITE(UNIT=nfprt,FMT=222)
    !
    !     nthin is the number of low layers strapped together
    !
    DO i=1,ncols
       DO kk=1,(kmax+1)
          IF( prsi(i,kk)/prsi(i,nsurf).GT.0.025_r8)EXIT
       END DO
       nthin(i)=kk
    END DO


!    IF(nthin.GT.1)THEN
       DO k=1,(kmax+1)
          DO i=1,ncols
             IF(nthin(i)>1 .and. k >=1 .and. k <= nthin(i))THEN
                tensio(i,k)=tensio(i,1)
                xtens(i,k)=xtens(i,1)
                ytens(i,k)=ytens(i,1)
             END IF
          END DO
       END DO
 !   ENDIF
    !
    !     scalar product of lower wind vector and surface wind
    !
    DO k=1,(kmax+1)
!    DO k=nthin+1,nbase
       DO i=1,ncols
          IF(k >=nthin(i)+1 .and. k <= nbase(i))THEN
             !
             !     *100 to convert to newton/m2
             !     velocity component paralell to surface velocity
             !
             velco=0.5_r8*((u(i,k-1)+u(i,k))*ubar(i)+ &
                  (v(i,k-1)+v(i,k))*vbar(i))/speeds(i)
             !
             !     tau doesn't change in the base layer because of a
             !     critical level i.e. velco < 0.0_r8
             !
             IF(velco.LE.0.0_r8)THEN
                tensio(i,k)=tensio(i,k-1)
                !
                !     froude number squared
                !
             ELSE
                roave=50.0_r8*(ro(i,k-1)+ro(i,k))
                fro2=bv(i,k)/(akwnmb*roave*velco*velco*velco)* &
                     tensio(i,k-1)
                !
                !     denominator of richardson number
                !
                delve2=(u(i,k)-u(i,k-1))*(u(i,k)-u(i,k-1))+ &
                      (v(i,k)-v(i,k-1))*(v(i,k)-v(i,k-1))
                !
                !     richardson number
                !
                IF(delve2.NE.0.0_r8)THEN
                   richsn=dz(i,k)*bv(i,k)*dz(i,k)*bv(i,k)/delve2
                ELSE
                   richsn=99999.0_r8
                ENDIF
                !
                !     tau in the base layer does not change because of the
                !     richardson criterion
                !
                IF(richsn.LE.0.25_r8)THEN
                   tensio(i,k)=tensio(i,k-1)
                   !
                   !     tau in the base layer does change if the local froude
                   !     excedes the critical froude number... the so called
                   !     froude number reduction.
                   !
                ELSE
                   crifro=(1.0_r8-0.25_r8/richsn)**2
                   IF(k.EQ.2)crifro= MIN (0.7_r8,crifro)
                   IF(fro2.GT.crifro)THEN
                      tensio(i,k)=crifro/fro2*tensio(i,k-1)
                   ELSE
                      tensio(i,k)=tensio(i,k-1)
                   ENDIF
                ENDIF
             ENDIF
             xtens(i,k)=tensio(i,k)*COS(ang(i))
             ytens(i,k)=tensio(i,k)*SIN(ang(i))
          END IF
       END DO
    END DO
    !
    !     stress from base level to top level
    !
    !DO k=nbase+1,(kmax+1)
    DO k=1,kMax
       DO i=1,ncols
          IF(k >=nbase(i)+1 .and. k <= (kmax+1))THEN
             !
             !     the stress is always initialized to zero
             !
             tensio(i,k)=0.0_r8
             IF(icrilv(i).NE.1.AND.k.NE.(kmax+1))THEN
                !
                !     vaisala frequency
                !
                vai1=(t(i,k)-t(i,k-1))/dz(i,k)+gocp
                vsqua=vai1*2.0_r8*gravi/(t(i,k)+t(i,k-1))
                !
                !     velocity component paralell to surface velocity
                !     scalar product of upper and surface wind vector
                !
                velco=0.5_r8*((u(i,k-1)+u(i,k))*ubar(i)+ &
                     (v(i,k-1)+v(i,k))*vbar(i))/speeds(i)
                !
                !     froude number squared
                !     fro2=vaisd/(akwnmb*roave*velco*velco*velco)*tensio(i,k-1)
                !     denominator of richardson number
                !
                delve2=(u(i,k)-u(i,k-1))*(u(i,k)-u(i,k-1))+ &
                       (v(i,k)-v(i,k-1))*(v(i,k)-v(i,k-1))
                !
                !     richardson number
                !
                IF(delve2.NE.0.0_r8)THEN
                   richsn=dz(i,k)*dz(i,k)*vsqua/delve2
                ELSE
                   richsn=99999.0_r8
                ENDIF
               IF(vai1.GE.0.0_r8.AND.velco.GE.0.0_r8.AND.richsn.GT.0.25_r8) THEN
                   !
                   !     *100 to convert to newton/m2
                   !
                   roave=50.0_r8*(ro(i,k-1)+ro(i,k))
                   !
                   !     vaisala frequency
                   !
                   fro2=SQRT(vsqua)/(akwnmb*roave*velco*velco*velco)* &
                         tensio(i,k-1)
                   !
                   !     critical froude number
                   !
                   crifro=(1.0_r8-0.25_r8/richsn)**2
                   !
                   !     end critical froude number
                   !
                   IF(fro2.GE.crifro)THEN
                      tensio(i,k)=crifro/fro2*tensio(i,k-1)
                   ELSE
                      tensio(i,k)=tensio(i,k-1)
                   ENDIF
                   xtens(i,k)=tensio(i,k)*COS(ang(i))
                   ytens(i,k)=tensio(i,k)*SIN(ang(i))
                ELSE
                   icrilv(i)   = 1
                   tensio(i,k) = 0.0_r8
                   xtens(i,k)  = 0.0_r8
                  ytens(i,k)  = 0.0_r8
               ENDIF
             ELSE
                tensio(i,k) = 0.0_r8
                xtens(i,k)  = 0.0_r8
                ytens(i,k)  = 0.0_r8
             ENDIF
          END IF
       END DO
    END DO
    !
    !     end stress
    !
    !     momentum change for free atmosphere
    !
    !DO k=nthin+1,kmax
    DO k=1,kmax
       DO i=1,ncols
          IF(k >=nthin(i)+1 .and. k <= (kmax))THEN
             chug(i,k)=-coef(i)/DeltaP(i,k)*(xtens(i,k+1)-xtens(i,k))
             chvg(i,k)=-coef(i)/DeltaP(i,k)*(ytens(i,k+1)-ytens(i,k))
          END IF
       END DO
    END DO
    !
    !     momentum change near the surface
    !     if lowest layer is very thin, it is strapped to next layer
    !
    !dsigma=si(nthin+1)-si(1)
    DO i=1,ncols
       dsigma=(prsi(i,nthin(i)+1)   - prsi(i,1))/prsi(i,nsurf)
       chug(i,1)=coef(i)/dsigma*(xtens(i,nthin(i)+1)-xtens(i,1))
       chvg(i,1)=coef(i)/dsigma*(ytens(i,nthin(i)+1)-ytens(i,1))
    END DO
    !IF(nthin.GT.1)THEN
       !DO k=2,nthin
       DO k=1,kmax
          DO i=1,ncols
             IF(nthin(i)> 1 .and. k >=2 .and. k <= nthin(i))THEN
                chug(i,k)=chug(i,1)
                chvg(i,k)=chvg(i,1)
             END IF 
          END DO
       END DO
    !ENDIF
222 FORMAT('SPEEDS EQ ZERO IN GWD')
9976 FORMAT(' MODEL TOO COARSE FOR MULTILAYER BASE LAYER.'/, &
         ' CHANGE THE MODEL VERTICAL STRUCTURE OR RUN WITH IGWD SET TO NO.'/, &
         ' NBASE=',I5,' SI='/(' ',5G16.8/))

  END SUBROUTINE GwddSchAlpert
END MODULE GwddSchemeAlpert
