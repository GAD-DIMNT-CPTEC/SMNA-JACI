!
!  $Author: pkubota $
!  $Date: 2009/07/15 18:51:26 $
!  $Revision: 1.7 $
!
MODULE PostLoop

  USE Constants, ONLY : r8, nferr, nfprt, nfpos, nfctl, nfdir, nfrfd, &
                        ndv, mdf, kdv, nvv, iclcd, nudv, nureq, chrdv, chreq, &
                        SHmin, Undef, Binary, EMRad1, res,RunRecort
  USE Sizes,    ONLY : Imax, Jmax, Kmax, Lmax, Mnmax, mymnmax, ibmax,haveM1, &
                       jbmax, kmaxloc, lmaxloc, ReshapeVerticalGroups, &
                       iPerIJB, jPerIJB, &
                       ibmaxperjb, havesurf, mymnextmax, mymmax, lm2m, mymnmap
  USE Parallelism, ONLY : myid
  USE PrblSize, ONLY : ngaus, pmand, alnpmd
  USE FileAccess, ONLY : ReadHeader, WriteField, opnpos, CloseFiles, skipf, &
                         ReadFieldG, ReadFieldSp, ReadField
  USE FileGrib, ONLY : GDSPDSSETION, WriteGrbField
  USE Communications, ONLY : Collect_Grid_Full
  USE Conversion, ONLY : GiveUnit, cnvout 
  USE Utils, ONLY : scase, lati, cosiv, rcl, coslat
  USE InputArrays, ONLY : qlnpp, qlnppmb,qlnppPA,qplam, qpphi, qpphiPa,qgzs, wd, pscb, lnpscb,lnpsPa, &
                          qup, qvp, qrotp, qdivp, qtmpp, qqp, &
                          qug, qvg, qrotg, qdivg, wk, wl, wlphi, &
                          psmb, top, lsmk, fgplam, fgplamPa,fgpphi, fgpphiPa,fgomega, fgdivq, &
                          og, ga_l, fgu, fgv, fgdiv, fgq, fgtmp, qa, ga,&
                          z0rl,tsfc,dest,stcp,stso,wssf,wrzo,wdzo,t2mt,q2mt,td2m
  USE SigmaToPressure, ONLY : sigtop, gavint, sig2pz, sig2po
  USE GaussPressure, ONLY : getth, getsh, getslp, lowtmp
  USE GaussSigma, ONLY : pwater, omegas, a_hybr, b_hybr, c_hybr, a_hybr_cb
  USE RegInterp, ONLY : Idim=>IdimOut, Jdim=>JdimOut, mgaus, gLats, &
                        RegInt, DoAreaInterpolation,DoAreaGausInterpolation
  USE tables, ONLY: table1,table2,table3,size_tb
  USE Transform, ONLY: Rectrg, DepositGridToSpec, DoGridToSpec, &
                       DestroyGridToSpec, CreateGridToSpec, &
                       CreateSpecToGrid, DepositSpecToGrid, &
                       DepositSpecToGridAndDelLamGrid, DoSpecToGrid, &
                       DestroySpecToGrid, DepositSpecToDelLamGrid
  USE SpecDynamics, ONLY : dztouv, gozrim, Uvtodz, snnp1
  USE infgdspds, only : hh_fct,dd_fct,mm_fct,yy_fct,hh_anl, dd_anl, &
                        mm_anl, yy_anl
  USE PhysicalFunctions, ONLY: CalcDewPoint

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: InitPostLoop
  PUBLIC :: postgl
  INTEGER, PARAMETER :: nopf=36   !PK

  INTEGER, PARAMETER :: ndi=250
  INTEGER, PARAMETER :: ndp=ndi+ndv!207

  INTEGER :: nharm(ndi) ! either ngaus or 2*mnmax
  INTEGER :: nlevs(ndi) ! number of levels for this field
  INTEGER :: ifldcd(ndi)! ????
  INTEGER :: mkdir(ndi)
  INTEGER :: nuc(ndp)
  INTEGER :: mko(ndp)
  INTEGER :: ife(ndp)
  INTEGER :: nvo(ndp)
  INTEGER :: nuco(ndp)
  INTEGER :: nfe(ndp)
  INTEGER :: jfe(ndp)
  INTEGER :: mkdv(ndv)  !pk
  INTEGER :: lif(ndv)   !pk

  CHARACTER (LEN=4)  :: alias(ndp)
  CHARACTER (LEN=4)  :: aliop(ndp)
  CHARACTER (LEN=4)  :: prodia(ndi) ! either PROG or DIAG
  CHARACTER (LEN=40)  :: chrdsc(ndi) ! input file field name
  CHARACTER (LEN=40)  :: chrdo(ndp)
  CHARACTER (LEN=40)  :: chrop(ndp)

  CHARACTER (LEN=4) :: rcode='N000' !irrelevant
  CHARACTER (LEN=40) :: title='NEW POSTP RHOMB/TRIANG OUTPUT TEST      '

  LOGICAL  :: first
  LOGICAL  :: second2
  LOGICAL  :: dopf(nopf,2)  !PK

  CHARACTER (LEN=40), DIMENSION (5,2) :: inf
  CHARACTER (LEN=40), DIMENSION (nopf,2) :: opf  !PK

  CHARACTER (LEN=40) :: t2ps='TIME MEAN SURFACE PRESSURE              '

  INCLUDE 'mpif.h'

CONTAINS
  SUBROUTINE InitPostLoop()
    IMPLICIT NONE
    first=.TRUE.
    second2=.FALSE.
  END SUBROUTINE InitPostLoop


  SUBROUTINE postgl (nFFrs, nFBeg, nFEnd, nFile)

    IMPLICIT NONE

    ! reads spectral forecast coefficients of topography,
    ! log of sfc pressure, temperature, divergence, vorticity
    ! and humidity in sigma layers. converts these values
    ! to selected mandatory pressure levels.

    INTEGER, INTENT(IN) :: nFFrs
    INTEGER, INTENT(IN) :: nFBeg
    INTEGER, INTENT(IN) :: nFEnd
    INTEGER, INTENT(IN) :: nFile

    INTEGER :: indate(4)
    INTEGER :: idate(4)
    INTEGER :: idatec(4)
    INTEGER :: nfld
    INTEGER :: nflp
    INTEGER :: nof
    INTEGER :: kp
    INTEGER :: ihr
    INTEGER :: iday
    INTEGER :: month
    INTEGER :: iyr
    INTEGER :: kfld
    INTEGER :: kflo
    INTEGER :: ks
    INTEGER :: ifday
    INTEGER :: i,ii
    INTEGER :: j
    INTEGER :: kpds(200),kgds(200)
    INTEGER :: ierr
    INTEGER :: kfl_save

    REAL (KIND=r8) :: fhour
    REAL (KIND=r8) :: hr
    REAL (KIND=r8) :: tod

    REAL (KIND=r8) :: qg(Ibmax,Jbmax) !intent(out) rwrite
    REAL (KIND=r8) :: topreg(Idim,Jdim), lsmkreg(Idim,Jdim)

    LOGICAL :: mean
    LOGICAL :: newday

    CHARACTER (LEN=4) :: dtin !irrelevant
    CHARACTER (LEN=10) :: labelp
    CHARACTER (LEN=40) :: specal !irrelevant
    CHARACTER (LEN=256) :: fname
    CHARACTER (LEN=4), PARAMETER :: diag='DIAG'
!    CHARACTER (LEN=4), SAVE :: rcode='N000' !irrelevant
    CHARACTER (LEN=20), PARAMETER :: type='PRESSURE HISTORY    '
!    CHARACTER (LEN=40), SAVE :: title='NEW POSTP RHOMB/TRIANG OUTPUT TEST      '
    CHARACTER (LEN=3), PARAMETER :: cmth(12)=(/&
         'JAN','FEB','MAR','APR','MAY','JUN', &
         'JUL','AUG','SEP','OCT','NOV','DEC'/)

    ! read input sigma file

    IF (nFBeg > nFile) GOTO 4500
    CALL opnpos (labelp)
    CALL recon (nfld, nflp, nof, indate, title, &
         specal, rcode, dtin, nFile)

    IF (nfld == 0) GOTO 4500

    CALL ReadHeader (ifday, tod, idate, idatec)
    fhour=REAL( ifday*24,r8) + tod/3600.0_r8
    ihr=idatec(1)
    month=idatec(2)
    iday=idatec(3)
    iyr=idatec(4)
    hr=REAL(ihr,r8)+MOD(tod,3600.0_r8)/3600.0_r8
    IF (Myid.ne.0) GOTO 3000
    WRITE (UNIT=*, FMT='(A,I6,A,F8.1)') ' ForecastDay =', ifday, '  TimeOfDay =', tod
    WRITE (UNIT=*, FMT='(A,3I3,I5)')    ' InitialDate =', idate
    WRITE (UNIT=*, FMT='(A,3I3,I5)')    ' CurrentDate =', idatec

    if (binary)then
       CALL GeraBinCtl(12,fname,title,labelp,nof,ndp,nvo,aliop,chrop,nuco)
    else
       CALL GeraGribCtl(12,fname,title,labelp,nof,ndp,nvo,aliop)
    endif

    WRITE (UNIT=nfprt, FMT='(2(A,I5))') ' nFFrs = ', nFFrs,' nFEnd =', nFEnd
    WRITE (UNIT=nfprt, FMT='(2(A,I5))') ' nFile = ', nFile,' nFBeg =', nFBeg


    WRITE (UNIT=*, FMT='(/,A,F6.2,2X,8I5,/)') &
         '  Hour = ', hr, indate, ihr, month, iday, iyr

3000 CONTINUE
    !
    !   for  Infgdspds module (common for GDSPDSSETION)
    !
    hh_fct=ihr
    dd_fct=iday
    mm_fct=month
    yy_fct=iyr
    hh_anl=idate(1)
    dd_anl=idate(3)
    mm_anl=idate(2)
    yy_anl=idate(4)

    CALL GDSPDSSETION (kgds, kpds )

    IF (myid.eq.0) THEN
       print*,'Dimensions: nx = ',Idim,' ny = ',Jdim
       print*,'LATITUDES = ',LATI(1),LATI(Jdim)
    ENDIF

    ! read spectral coefficients of orography
    ! and land sea mask

    CALL ReshapeVerticalGroups(kmax, kmaxloc)
    CALL ReadFieldSp (2*Mnmax, 1, qgzs)
    CALL rectrg (qgzs, top, 1, 1)
    CALL ReadFieldG (ngaus,1,lsmk)
    CALL Collect_Grid_Full(top, ga(:,:,1), 1, 0)
    CALL Collect_Grid_Full(lsmk, ga(:,:,2), 1, 0)
    DO j=1,SIZE(ga,2)
       DO i=1,SIZE(ga,1)
          IF(ga(i,j,1) < 0.0_r8)THEN
             ga(i,j,1)=0.0_r8
          ENDIF
       END DO
    END DO
    IF(myid.ne.0) GOTO 3100

    IF (RegInt) THEN
       IF (Binary) THEN
          !
          ! write binary format
          !
          CALL DoAreaInterpolation (ga(:,:,1), topreg)
          CALL WriteField (mgaus, topreg)
          CALL DoAreaInterpolation (ga(:,:,2), lsmkreg)
          CALL WriteField (mgaus, lsmkreg)
          !
       ELSE
          !
          ! write grib format
          !
          CALL DoAreaInterpolation (ga(:,:,1), topreg)
          CALL WriteGrbField ('TOPO',mgaus, kgds, kpds, topreg,1)
          CALL DoAreaInterpolation (ga(:,:,2), lsmkreg)
          CALL WriteGrbField ('LSMK',mgaus, kgds, kpds, lsmkreg,1)
          !       
       END IF
    ELSE
       IF (Binary) THEN
          !
          ! write binary format
          !
          CALL DoAreaGausInterpolation(ga(:,:,1), topreg)
          CALL WriteField (mgaus, topreg)
          CALL DoAreaGausInterpolation(ga(:,:,2), lsmkreg)
          CALL WriteField (mgaus, lsmkreg)
          !
       ELSE
          !
          ! write grib format
          !
          CALL DoAreaGausInterpolation(ga(:,:,1), topreg)
          CALL WriteGrbField ('TOPO',mgaus, kgds, kpds, topreg,1)
          CALL DoAreaGausInterpolation(ga(:,:,2), lsmkreg)
          CALL WriteGrbField ('LSMK',mgaus, kgds, kpds, lsmkreg,1)
          !
       ENDIF
    END IF

    ! compute sig to p on ggrid
    mean=.FALSE.
    if(myid.eq.0) WRITE (UNIT=nfprt, FMT='(A,I5,2I3,F7.2)') ' Processing: ', iyr, month, iday, hr
    if(myid.eq.0) WRITE (UNIT=nfprt, FMT='(A,L6)') ' Starting postgg: Mean = ', mean

3100 CONTINUE
    kfld=1
    kflo=1
    mean=.FALSE.
    newday=.TRUE.
    CALL postgg (psmb, pmand, alnpmd, mean, kfld, kflo, newday, top,nfld)
    IF(myid.eq.0) WRITE (UNIT=nfprt, FMT='(/,A,/)') ' postgg Finished'
    ! skip PROG fields
    DO kp=1,nfld
       IF (prodia(kp) == diag) EXIT
    ENDDO
    kp=kp-1
    ks=kfld
    ! loop over remaining PROG fields
    DO kfld=ks,kp
       CALL ReshapeVerticalGroups(kmax, kmaxloc)
       IF(myid.eq.0) WRITE (UNIT=nfprt, FMT='(A, I3)') ' Processing Input Field Number: ', kfld
       IF (mkdir(kfld)<1 .OR. mkdir(kfld)==100) THEN
          IF(myid.eq.0) WRITE (UNIT=nfprt, FMT='(A,I5)') ' Prog skipf(nlevs(kfld)), nLevs = ',nlevs(kfld)
          CALL skipf (nlevs(kfld))
       ELSEIF (nharm(kfld) == ngaus) THEN
          kfl_save = kflo
          IF (nlevs(kfld).gt.1) THEN
             CALL ReadFieldG(ngaus, nlevs(kfld), ga_l)
             IF (MOD(mkdir(kfld),2) == 1) THEN
                CALL gavint(nlevs(kfld),nvo(kflo),ga_l,og,psmb,pmand)
             ELSEIF (mkdir(kfld) > 1) THEN
                CALL scase (kflo,nfe,iclcd,ga_l)
                CALL gavint (nlevs(kfld),nvo(kflo),ga_l,og,psmb,pmand)
             ENDIF
             IF (myid.eq.0)PRINT*,'pkubota4',kfld,mkdir(kfld),chrop(kflo),nlevs(kflo)
             CALL rwrite (nvo(kflo),kflo,og,ga,.true.)
             kflo=kflo+1
           ELSE
             IF (myid.eq.0) THEN
                CALL ReadField (ngaus, nlevs(kfld), ga)
                IF (MOD(mkdir(kfld),2) == 1) THEN
                      PRINT*,'pkubota1',kfld,mkdir(kfld),chrop(kflo)

                   CALL rwrite (nvo(kflo),kflo,og,ga,.false.)
                   kflo=kflo+1
                ELSEIF (mkdir(kfld) > 1) THEN
                   CALL scase (kflo,nfe,iclcd,ga)
                      PRINT*,'pkubota2',kfld,mkdir(kfld),chrop(kflo)
                   CALL rwrite (nvo(kflo),kflo,og,ga,.false.)
                   kflo=kflo+1
                ENDIF
              ELSE
                CALL skipf (nlevs(kfld))
                kflo=kflo+1
             ENDIF
          ENDIF
             IF(myid.eq.0.and.kfl_save.eq.kflo) WRITE (UNIT=nfprt, FMT='(A,I5)') ' Field Skipped, nLevs = ', nlevs(kfld)
       ELSEIF (nharm(kfld) == 2*mnmax) THEN
          CALL ReadFieldSp (2*mnmax, nlevs(kfld), wk)
          IF(myid.eq.0) WRITE (UNIT=nfprt, FMT='(A,L6)') ' Mean   = ', mean
          CALL stog (psmb, alnpmd, nlevs(kfld), mean, kflo, nof, kfld, wk, og, qg)
          IF (MOD(mkdir(kfld),2) == 1) THEN
             !og USED
             CALL rwrite (nvo(kflo), kflo, og, ga, .true.)
             kflo=kflo+1
          ENDIF
          IF (mkdir(kfld) > 1) THEN
             !qg USED
             CALL rwrite (nvo(kflo), kflo, qg, qa, .true.) 
             kflo=kflo+1
          ENDIF
       ENDIF
    END DO
    kfld=kp+1
    IF (kfld > nfld) THEN
       IF(myid.eq.0) WRITE (UNIT=nfprt, FMT='(/,A,I5,/)') ' Progs Finished, Diags Skipped for nFile = ', nFile
       GOTO 4500
    ENDIF
    IF (nFile <= 0) THEN
       IF(myid.eq.0) WRITE (UNIT=nfprt, FMT='(/,A,I5,/)') ' Progs Finished, Diags Skipped for nFile = ', nFile
       GOTO 4500
    ENDIF
    IF(myid.eq.0) WRITE (UNIT=nfprt, FMT='(/,A,/)') ' Progs Finished, Diags Starting ...'

    !     process diagnostic fields

    mean=.TRUE.
    IF(myid.eq.0) WRITE (UNIT=nfprt, FMT='(A,L6)') ' Starting postgg: Mean = ', mean
    CALL postgg2 (psmb, pmand, alnpmd, mean, kfld, kflo, newday, top)
    IF(myid.eq.0) WRITE (UNIT=nfprt, FMT='(/,A,/)') ' postgg Finished'
    ks=kfld
    !loop over remaining DIAG fields
    DO kfld=ks,nfld
       CALL ReshapeVerticalGroups(kmax, kmaxloc)
       IF(myid.eq.0) WRITE (UNIT=nfprt, FMT='(A, I3)') ' Processing Input Field Number: ', kfld
       IF (mkdir(kfld)<1 .OR. mkdir(kfld)==100) THEN
          IF(myid.eq.0) WRITE (UNIT=nfprt, FMT='(A,I5)') ' Prog skipf(nlevs(kfld)), nLevs = ',nlevs(kfld)
          CALL skipf(nlevs(kfld))
       ELSEIF (nharm(kfld) == ngaus) THEN
          kfl_save = kflo
 
          IF(chrop(kflo) == '2 METRE DEWPOINT TEMPERATURE            ') THEN
             !IF (myid.eq.0)PRINT*,'pkubota3',kfld,mkdir(kfld),chrop(kflo),nlevs(kflo)
             CALL rwrite (nvo(kflo), kflo, td2m, qa, .true.)
             kflo=kflo+1
          ENDIF
             IF (nlevs(kfld).gt.1) THEN
                CALL ReadFieldG(ngaus, nlevs(kfld), ga_l)

                IF (MOD(mkdir(kfld),2) == 1) THEN
                   CALL gavint (nlevs(kfld), nvo(kflo), ga_l, og, psmb, pmand)
                ENDIF

                IF (mkdir(kfld) > 1) THEN
                   CALL scase (kflo, nfe, iclcd, ga_l)
                   CALL gavint (nlevs(kfld), nvo(kflo), ga_l, og, psmb, pmand)
                ENDIF
               !IF (myid.eq.0)PRINT*,'pkubota4',kfld,mkdir(kfld),chrop(kflo),nlevs(kflo)

                CALL rwrite (nvo(kflo),kflo,og,ga,.true.)
                kflo=kflo+1
             ELSE
                IF (myid.eq.0) THEN
                   CALL ReadField (ngaus, nlevs(kfld), ga)
                   IF (MOD(mkdir(kfld),2) == 1) THEN
                      !PRINT*,'pkubota1',kfld,mkdir(kfld),chrop(kflo)

                      CALL rwrite (nvo(kflo),kflo,og,ga,.false.)
                      kflo=kflo+1
                   ENDIF

                   IF (mkdir(kfld) > 1) THEN
                      !PRINT*,'pkubota2',kfld,mkdir(kfld),chrop(kflo)
                       CALL scase (kflo,nfe,iclcd,ga)
                       CALL rwrite (nvo(kflo),kflo,og,ga,.false.)
                       kflo=kflo+1
                   ENDIF
                ELSE
                   CALL skipf (nlevs(kfld))
                   kflo=kflo+1
                ENDIF
             ENDIF

       ELSEIF (nharm(kfld) == 2*mnmax) THEN
          CALL ReadFieldSp (2*mnmax, nlevs(kfld), wk)
          IF(myid.eq.0) WRITE (UNIT=nfprt, FMT='(A,L6)') ' Mean   = ',mean
          CALL stog (psmb, alnpmd, nlevs(kfld), mean, kflo, nof, kfld, wk, og, qg)
          IF (MOD(mkdir(kfld),2) == 1) THEN
             CALL rwrite (nvo(kflo), kflo, og, ga, .true.)
             kflo=kflo+1
          ENDIF
          IF (mkdir(kfld) > 1) THEN
             CALL rwrite (nvo(kflo), kflo, qg, qa, .true.)
             kflo=kflo+1
          ENDIF
       ENDIF
       CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
    END DO
    IF(myid.eq.0) WRITE (UNIT=nfprt, FMT='(/,A,/)') ' Diagnostics Completed ...'

    ! end of diagnostics

4500 CONTINUE

   ! PK IF (myid.eq.0) THEN
       CALL CloseFiles()

       IF(.NOT. Binary)THEN
          CALL BACLOSE(51,ierr)
          CALL W3TAGE('GCMPOST ')
       ENDIF
    !PK ENDIF

  END SUBROUTINE postgl


  SUBROUTINE postgg (psmb, pmand, alnpmd, mean, kfld, kflo, newday, top,nfld)

    IMPLICIT NONE

    INTEGER, INTENT(INOUT) :: kfld
    INTEGER, INTENT(INOUT) :: kflo
    INTEGER, INTENT(INOUT) :: nfld

    REAL (KIND=r8), INTENT(IN) :: pmand(Lmax)
    REAL (KIND=r8), INTENT(IN) :: alnpmd(Lmax)
    REAL (KIND=r8), INTENT(IN) :: top(Ibmax,Jbmax)
    REAL (KIND=r8), INTENT(OUT) :: psmb(Ibmax,Jbmax)

    LOGICAL, INTENT(IN) :: mean
    LOGICAL, INTENT(INOUT) :: newday

    INTEGER, PARAMETER :: nopf=36  !PK
    INTEGER :: id
    INTEGER :: i
    INTEGER :: kt
    INTEGER :: mm
    INTEGER :: ii
    INTEGER :: j
    INTEGER :: ib,jb,kp
    INTEGER :: m, mg, mglob, n, mn, mf
    INTEGER :: l
    INTEGER :: ll
    INTEGER :: nn
    INTEGER :: nmax,nCols

    REAL (KIND=r8) :: ug(Ibmax,Lmax,Jbmax)
    REAL (KIND=r8) :: vg(Ibmax,Lmax,Jbmax)
    REAL (KIND=r8) :: rq(Ibmax,Lmax,Jbmax)
    REAL (KIND=r8) :: og(Ibmax,Lmax,Jbmax)

    REAL (KIND=r8) :: root2  
    CHARACTER (LEN=4), PARAMETER :: prog='DIAG'

    LOGICAL :: timemean
    !    LOGICAL, SAVE :: first=.TRUE.
    !    LOGICAL, SAVE :: second=.FALSE.
    !    LOGICAL, SAVE :: dopf(nopf,2)

    !    CHARACTER (LEN=40), DIMENSION (5,2), SAVE :: inf
    !    CHARACTER (LEN=40), DIMENSION (nopf,2), SAVE :: opf

    !    CHARACTER (LEN=40), SAVE :: t2ps='TIME MEAN SURFACE PRESSURE              '

    inf(:,1) = (/ &
         'LN SURFACE PRESSURE                     ', &
         'DIVERGENCE                              ', &
         'VORTICITY                               ', &
         'SPECIFIC HUMIDITY                       ', &
         'VIRTUAL TEMPERATURE                     ' /)
    inf(:,2) = (/ &
         'TIME MEAN LN SURFACE PRESSURE           ', &
         'TIME MEAN DIVERGENCE                    ', &
         'TIME MEAN VORTICITY                     ', &
         'TIME MEAN SPECIFIC HUMIDITY             ', &
         'TIME MEAN VIRTUAL TEMPERATURE           ' /)

    opf(:,1) = (/ &
         'SURFACE PRESSURE                        ', &! 1
         'MASK                                    ', &! 2
         'SURFACE ZONAL WIND (U)                  ', &! 3
         'ZONAL WIND (U)                          ', &! 4
         'SURFACE MERIDIONAL WIND (V)             ', &! 5
         'MERIDIONAL WIND (V)                     ', &! 6
         'OMEGA                                   ', &! 7
         'DIVERGENCE                              ', &! 8
         'VORTICITY                               ', &! 9
         'STREAM FUNCTION                         ', &!10
         'ZONAL WIND PSI                          ', &!11
         'MERIDIONAL WIND PSI                     ', &!12
         'VELOCITY POTENTIAL                      ', &!13
         'ZONAL WIND CHI                          ', &!14
         'MERIDIONAL WIND CHI                     ', &!15
         'VIRTUAL TEMPERATURE                     ', &!16
         'GEOPOTENTIAL HEIGHT                     ', &!17
         'SEA LEVEL PRESSURE                      ', &!18
         'SURFACE ABSOLUTE TEMPERATURE            ', &!19
         'ABSOLUTE TEMPERATURE                    ', &!20
         'SURFACE RELATIVE HUMIDITY               ', &!21
         'RELATIVE HUMIDITY                       ', &!22
         'SPECIFIC HUMIDITY                       ', &!23
         'INST. PRECIP. WATER                     ', &!24
         'POTENTIAL TEMPERATURE                   ', &!25
         'ROUGHNESS LENGTH                        ', &!26 !PK---
         'SURFACE TEMPERATURE                     ', &!27 
         'DEEP SOIL TEMPERATURE                   ', &!28 
         'STORAGE ON CANOPY                       ', &!29 
         'STORAGE ON GROUND                       ', &!30 
         'SOIL WETNESS OF SURFACE                 ', &!31 
         'SOIL WETNESS OF ROOT ZONE               ', &!32 
         'SOIL WETNESS OF DRAINAGE ZONE           ', &!33 
         'TEMPERATURE AT 2-M FROM SURFACE         ', &!34 
         '2 METRE DEWPOINT TEMPERATURE            ', &!35
         'SPECIFIC HUMIDITY AT 2-M FROM SURFACE   '/) !35 !PK---

    opf(:,2) = (/ &
         'TIME MEAN SURFACE PRESSURE              ', &! 1
         'TIME MEAN MASK                          ', &! 2
         'TIME MEAN SURFACE ZONAL WIND (U)        ', &! 3
         'TIME MEAN ZONAL WIND (U)                ', &! 4
         'TIME MEAN SURFACE MERIDIONAL WIND (V)   ', &! 5
         'TIME MEAN MERIDIONAL WIND (V)           ', &! 6
         'TIME MEAN DERIVED OMEGA                 ', &! 7
         'TIME MEAN DIVERGENCE                    ', &! 8
         'TIME MEAN VORTICITY                     ', &! 9
         'TIME MEAN STREAM FUNCTION               ', &!10
         'TIME MEAN ZONAL WIND PSI                ', &!11
         'TIME MEAN MERIDIONAL WIND PSI           ', &!12
         'TIME MEAN VELOCITY POTENTIAL            ', &!13
         'TIME MEAN ZONAL WIND CHI                ', &!14
         'TIME MEAN MERIDIONAL WIND CHI           ', &!15
         'TIME MEAN VIRTUAL TEMPERATURE           ', &!16
         'TIME MEAN GEOPOTENTIAL HEIGHT           ', &!17
         'TIME MEAN SEA LEVEL PRESSURE            ', &!18
         'TIME MEAN SURFACE ABSOLUTE TEMPERATURE  ', &!19
         'TIME MEAN ABSOLUTE TEMPERATURE          ', &!20
         'TIME MEAN SURFACE RELATIVE HUMIDITY     ', &!21
         'TIME MEAN RELATIVE HUMIDITY             ', &!22
         'TIME MEAN SPECIFIC HUMIDITY             ', &!23
         'TIME MEAN PRECIP. WATER                 ', &!24
         'TIME MEAN POTENTIAL TEMPERATURE         ', &!25
         'TIME MEAN_ROUGHNESS LENGTH              ', &!26 !PK---
         'TIME MEAN_SURFACE TEMPERATURE           ', &!27
         'TIME MEAN_DEEP SOIL TEMPERATURE         ', &!28
         'TIME MEAN_STORAGE ON CANOPY             ', &!29
         'TIME MEAN_STORAGE ON GROUND             ', &!30
         'TIME MEAN_SOIL WETNESS OF SURFACE       ', &!31
         'TIME MEAN_SOIL WETNESS OF ROOT ZONE     ', &!32
         'TIME MEAN_SOIL WETNESS OF DRAINAGE ZONE ', &!33
         'TIME MEAN_TEMPERATURE AT 2M FROM SURFACE', &!34
         'TIME MEAN 2 METRE DEWPOINT TEMPERATURE  ', &!35 !PK
         'TIME MEAN_SPEC-HUMID AT 2-M FROM SURFACE'/) !35 !PK---

    timemean=.TRUE.
    IF (first) THEN

       IF (mean) THEN
          IF (myid.eq.0) THEN
             WRITE (UNIT=nferr, FMT='(A3L4)') &
                ' First, Second, Mean Inconsistent:',&
                  first, second2, mean
          END IF
          STOP 8100
       END IF
       id=1

       ! input field has to be complete and well ordered,

       kt=kfld
       DO mm=1,5
          IF (inf(mm,1) == chrdsc(kt)) THEN
             kt=kt+1
          ELSE
             IF (myid.eq.0) THEN
                WRITE (UNIT=nferr, FMT='(A,L2,2(A,I4),A40,/,40X,A40)') &
                   ' At Mean = ', mean, ' mm = ', mm, ' kt = ', kt, &
                   ' inf = ', inf(mm,id), ' chrdsc = ', chrdsc(kt)
             END IF
             STOP 9100
          END IF
       END DO

       ! output field has to be well ordered.
       ! ignore undesired output fields.

       kt=kflo
       DO ii=1,nopf
            IF (myid.eq.0) WRITE (UNIT=nfprt, FMT='(" opf(ii,1) = ",A)') opf(ii,1)
            IF (myid.eq.0) WRITE (UNIT=nfprt, FMT='(" chrop(kt) = ",A)') chrop(kt)

          IF (opf(ii,1) == chrop(kt)) THEN
!            IF (myid.eq.0) WRITE (UNIT=nfprt, FMT='(" opf(ii,1) = ",A)') opf(ii,1)
!            IF (myid.eq.0) WRITE (UNIT=nfprt, FMT='(" chrop(kt) = ",A)') chrop(kt)
             dopf(ii,1)=.TRUE.
             kt=kt+1
          ELSE
             dopf(ii,1)=.FALSE.
          END IF
       END DO
       first=.FALSE.
       second2=.TRUE.
       IF (myid.eq.0) &
       WRITE (UNIT=nfprt, FMT='(" dopf(,1) = ",50L2)') (dopf(i,1),i=1,nopf)

    ELSE IF (second2) THEN
       IF (.NOT. newday) THEN
          IF (.NOT.mean) THEN
             IF (myid.eq.0) &
             WRITE (UNIT=nferr, FMT='(A3L4)') &
                   ' First, Second, Mean Inconsistent:',&
                     first, second2, mean
             STOP 8100
          END IF

          id=2
          inf(1,2)=t2ps

          ! check that the time mean input field is complete and well ordered

          kt=kfld
          DO mm=1,5
             IF (inf(mm,2) == chrdsc(kt)) THEN
                kt=kt+1
             ELSE
                IF (myid.eq.0) &
                WRITE (UNIT=nfprt, FMT='(A)') ' Time Mean is not Available, Ignore it ...'
                timemean=.FALSE.
             END IF
          END DO

          ! check that the time mean output field is well ordered.
          ! ignore undesired output fields.

          kt=kflo
          DO ii=1,nopf
             IF (opf(ii,2) == chrop(kt)) THEN
                dopf(ii,2)=.TRUE.
                kt=kt+1
             ELSE
                dopf(ii,2)=.FALSE.
             END IF
          END DO
          IF (myid.eq.0) &
          WRITE (UNIT=nfprt, FMT='(" dopf(,2) = ",50L2)') (dopf(i,2),i=1,nopf)

       END IF
    END IF
    root2  = SQRT(2.0_r8)
    IF (timemean) THEN

       ! surface pressure

       IF (.NOT.mean) THEN

          ! case full field:
          ! read spectral ln surface pressure;
          ! legandre transform to grid;
          ! get surface pressure and convert from centibar to milibar

          IF (myid.eq.0) &
          WRITE (nfprt,*) ' lendo qlnpp spectral '
          id=1
          CALL ReadFieldSp (nharm(kfld), nlevs(kfld), qlnppPA)
          CALL Rectrg(qlnppPA, lnpsPa, 1, 1)
          qlnppmb=qlnppPA
          qlnpp=qlnppPA
          !      transform surface pressure from Pascal to cbar
          qlnpp(1) = qlnpp(1) - log(1000._r8) * root2
          !      transform surface pressure from Pascal to mbar
          qlnppmb(1) = qlnppmb(1) - log(100._r8) * root2
       ELSE

          ! case time mean:
          ! read grid surface pressure;
          ! convert from centibar to milibar and get log
          ! legendre transform from spectral to grid

          id=2
          pscb = 0.
          CALL ReadFieldG(nharm(kfld), nlevs(kfld), pscb)
          DO j=1,Jbmax
             DO i=1,Ibmaxperjb(j)
                psmb(i,j)=pscb(i,j)*10.0_r8
                lnpscb(i,j)=LOG(pscb(i,j))
             END DO
          END DO
 !PK         CALL Rectrg(qlnpp, lnpscb, 1, 1)
          CALL Rectrg(qlnppPA, lnpsPa, 1, 1)
          !      transform surface pressure from Pascal to cbar
          DO j=1,Jbmax
             DO i=1,Ibmaxperjb(j)
                 lnpscb(i,j) = EXP(lnpsPa(i,j))/1000.0_r8
                 lnpscb(i,j) = log(lnpscb(i,j))
             END DO
          END DO

       END IF

       ! divergence
       kfld=kfld+1
       CALL ReadFieldSp (nharm(kfld), nlevs(kfld), qdivp)
       ! vorticity
       kfld=kfld+1
       CALL ReadFieldSp (nharm(kfld), nlevs(kfld), qrotp)
       ! humidity 
       kfld=kfld+1
       CALL ReadFieldSp (nharm(kfld), nlevs(kfld), qqp)
       ! temperature
       kfld=kfld+1
       CALL ReadFieldSp (nharm(kfld), nlevs(kfld), qtmpp)
       kfld=kfld+1
       !PK ----begin
       ! ROUGHNESS LENGTH 
       CALL ReadFieldG(nharm(kfld), nlevs(kfld), z0rl) !(ga_l(ibmax,jbmax))
       kfld=kfld+1
       ! SURFACE TEMPERATURE
       CALL ReadFieldG(nharm(kfld), nlevs(kfld), tsfc) !(ga_l(ibmax,jbmax))
       tsfc=ABS(tsfc)
       kfld=kfld+1
       ! DEEP SOIL TEMPERATURE
       CALL ReadFieldG(nharm(kfld), nlevs(kfld), dest) !(ga_l(ibmax,jbmax))
       kfld=kfld+1
       ! STORAGE ON CANOPY
       CALL ReadFieldG(nharm(kfld), nlevs(kfld), stcp) !(ga_l(ibmax,jbmax))
       kfld=kfld+1
       ! STORAGE ON GROUND 
       CALL ReadFieldG(nharm(kfld), nlevs(kfld), stso) !(ga_l(ibmax,jbmax))
       kfld=kfld+1
       ! SOIL WETNESS OF SURFACE 
       CALL ReadFieldG(nharm(kfld), nlevs(kfld), wssf) !(ga_l(ibmax,jbmax))
       kfld=kfld+1
       ! SOIL WETNESS OF ROOT ZONE 
       CALL ReadFieldG(nharm(kfld), nlevs(kfld), wrzo) !(ga_l(ibmax,jbmax))
       kfld=kfld+1
       ! SOIL WETNESS OF DRAINAGE ZONE
       CALL ReadFieldG(nharm(kfld), nlevs(kfld), wdzo) !(ga_l(ibmax,jbmax))
       kfld=kfld+1
       ! TEMPERATURE AT 2-M FROM SURFACE 
       CALL ReadFieldG(nharm(kfld), nlevs(kfld), t2mt) !(ga_l(ibmax,jbmax))
       kfld=kfld+1
       ! SPECIFIC HUMIDITY AT 2-M FROM SURFACE
       CALL ReadFieldG(nharm(kfld), nlevs(kfld), q2mt) !(ga_l(ibmax,jbmax))
       kfld=kfld+1
       !PK----end

!
       !  Transform all fields to grid space (kmax levels)
!
       IF (havesurf) CALL gozrim(qlnpp, qpphi, 1, 2*mymnExtMax)
       IF (havesurf) CALL gozrim(qlnppPa, qpphiPa, 1, 2*mymnExtMax)

       CALL dztouv(qdivp, qrotp, qup, qvp, 1, 2*mymnExtMax)
!
       CALL CreateSpecToGrid(5, 4, 5, 6)
       CALL DepositSpecToGridAndDelLamGrid(qlnpp, lnpscb, fgplam)
       CALL DepositSpecToGridAndDelLamGrid(qlnppPa, lnpsPa, fgplamPa)
!
       CALL DepositSpecToGrid(qpphi, fgpphi)
       CALL DepositSpecToGrid(qpphiPa, fgpphiPa)
!
       CALL DepositSpecToGrid(qup,fgu)
       CALL DepositSpecToGrid(qvp,fgv)
       CALL DepositSpecToGrid(qdivp, fgdiv)
       CALL DepositSpecToGrid(qtmpp, fgtmp)
       CALL DepositSpecToGrid(qqp, fgq)
!
       CALL DoSpecToGrid()
       CALL DestroySpecToGrid()
!
       IF (.NOT.mean) THEN
          DO j=1,Jbmax
             DO i=1,Ibmaxperjb(j)
!PK                pscb(i,j)=EXP(lnpscb(i,j))/100.0_r8
                pscb(i,j)=EXP(lnpsPa(i,j))/1000.0_r8
                psmb(i,j)=pscb(i,j)*10.0_r8
             END DO
          END DO
       ENDIF
       IF (dopf(1,id) .OR. dopf(2,id)) THEN
       IF(myid.eq.0) WRITE (UNIT=nfprt, FMT='(A)') 'postgg: ps and masK'

          ! generates mask array containing
          ! ones above surface and -ones below surface.

          DO j=1,Jbmax
             DO l=1,Lmax
                DO i=1,Ibmaxperjb(j)
                   IF (psmb(i,j) < pmand(l)) THEN
                      og(i,l,j)=-pmand(l)/psmb(i,j)
                   ELSE
                      og(i,l,j)=pmand(l)/psmb(i,j)
                   END IF
                END DO
             END DO
          END DO

          ! surface pressure

          IF (dopf(1,id)) THEN
             CALL rwrite (nvo(kflo), kflo, pscb, qa, .true.)
             kflo=kflo+1
          END IF

          ! mask for vertical extrapolation inside terrain

          IF (dopf(2,id)) THEN
             CALL rwrite (nvo(kflo), kflo, og, ga, .true.)
             kflo=kflo+1
          END IF

       END IF

       IF (dopf(3,id) .OR. dopf(4,id) .OR. dopf(5,id) .OR. &
            dopf(6,id) .OR. dopf(7,id) .OR. dopf(8,id) .OR. &
            dopf(10,id) .OR. dopf(11,id) .OR. dopf(12,id) .OR. &
            dopf(13,id) .OR. dopf(14,id) .OR. dopf(15,id)) THEN
          IF (myid.eq.0) WRITE (UNIT=nfprt, FMT='(A)') ' Winds: us, up, vs, vp and omega'
          CALL Winds (psmb, alnpmd, kflo, dopf, id, nvo, ndp, &
                      ug, vg, og)
       END IF

       ! this is a subset of the above condition with the exception of
       ! dopf(9,id)
       IF (dopf(8,id) .OR. dopf(9,id) .OR. &
            dopf(10,id) .OR. dopf(11,id) .OR. dopf(12,id) .OR. &
            dopf(13,id) .OR. dopf(14,id) .OR. dopf(15,id)) THEN
          IF (myid.eq.0) WRITE (UNIT=nfprt, FMT='(A)') ' uvtodz: div and vort'
          ! Recompute Divergence and Vorticy using U an V at p-levels.
          CALL ReshapeVerticalGroups(lmax, kmaxloc)
          CALL CreateGridToSpec(2, 0)
          CALL DepositGridToSpec(qug, ug)
          CALL DepositGridToSpec(qvg, vg)
          CALL DoGridToSpec()
          CALL DestroyGridToSpec()
          CALL Uvtodz(qug, qvg, qdivg, qrotg, 1, 2*mymnmax)
          CALL CreateSpecToGrid(2, 0, 2, 0)
          CALL DepositSpecToGrid(qdivg, ug)
          CALL DepositSpecToGrid(qrotg, vg)
          CALL DoSpecToGrid()
          CALL DestroySpecToGrid()
       END IF

       ! horizontal divergence at p levels

       IF (dopf(8,id)) THEN
          CALL rwrite (nvo(kflo), kflo, ug, ga, .true.)
          kflo=kflo+1
       END IF

       ! vertical componente of vorticity at p levels

       IF (dopf(9,id)) THEN
          CALL rwrite (nvo(kflo), kflo, vg, ga, .true.)
          kflo=kflo+1
       END IF

       ! stream function at p levels

       IF (dopf(10,id) .OR. dopf(11,id) .OR. dopf(12,id)) THEN
          DO l=1,Lmaxloc
             mf = 1
             IF (haveM1) THEN 
                wk(1,l)=0.0_r8
                wk(2,l)=0.0_r8
                mf = 3
             ENDIF
             DO m=mf,2*mymnmax
                wk(m,l)=qrotg(m,l)*snnp1(m)
             END DO
          END DO
          IF (dopf(10,id)) THEN
             if(myid.eq.0) WRITE (UNIT=nfprt, FMT='(A)') ' postgg: psi'
             CALL rectrg (wk, ug, Lmaxloc, Lmax)
             CALL rwrite (nvo(kflo), kflo, ug, ga, .true.)
             kflo=kflo+1
          END IF
       END IF

       ! zonal wind psi at p levels

       IF (dopf(11,id)) THEN
          if(myid.eq.0) WRITE (UNIT=nfprt, FMT='(A)') ' postgg: upsi'
          CALL CreateSpecToGrid(2, 0, 2, 0)
          CALL DepositSpecToDelLamGrid(wk,vg)
          CALL gozrim(wk, wlphi, 1, 2*mymnextmax)
          CALL DepositSpecToGrid(wlphi, ug)
          CALL DoSpecToGrid()
          CALL DestroySpecToGrid()
          DO j=1,Jbmax
             DO l=1,Lmax
                DO i=1,Ibmaxperjb(j)
                   ug(i,l,j)=-ug(i,l,j)*cosiv(i,j)
                END DO
             END DO
          END DO
          CALL rwrite (nvo(kflo), kflo, ug, ga, .true.)
          kflo=kflo+1
       END IF

       ! meridional wind psi at p levels

       IF (dopf(12,id)) THEN
          if(myid.eq.0) WRITE (UNIT=nfprt, FMT='(A)') ' postgg: vpsi'
          DO j=1,Jbmax
             DO l=1,Lmax
                DO i=1,Ibmaxperjb(j)
                   vg(i,l,j)=vg(i,l,j)*cosiv(i,j)
                END DO
             END DO
          END DO
          CALL rwrite (nvo(kflo), kflo, vg, ga, .true.)
          kflo=kflo+1
       END IF

       ! velocity potential at p levels

       IF (dopf(13,id) .OR. dopf(14,id) .OR. dopf(15,id)) THEN
          DO l=1,Lmaxloc
             mf = 1
             IF (haveM1) THEN 
                wk(1,l)=0.0_r8
                wk(2,l)=0.0_r8
                mf = 3
             ENDIF
             DO m=mf,2*mymnmax
                wk(m,l)=qdivg(m,l)*snnp1(m)
             END DO
          END DO
          IF (dopf(13,id)) THEN
             if(myid.eq.0) WRITE (UNIT=nfprt, FMT='(A)') ' postgg: chi'
             CALL rectrg (wk, vg, Lmaxloc, Lmax)
             CALL rwrite (nvo(kflo), kflo, vg, ga, .true.)
             kflo=kflo+1
          END IF
       END IF

       ! zonal wind chi at p levels

       IF (dopf(14,id)) THEN
          if(myid.eq.0) WRITE (UNIT=nfprt, FMT='(A)') ' postgg: uchi'
          CALL CreateSpecToGrid(2, 0, 2, 0)
          CALL DepositSpecToDelLamGrid(wk,ug)
          CALL gozrim(wk, wlphi, 1, 2*mymnextmax)
          CALL DepositSpecToGrid(wlphi, vg)
          CALL DoSpecToGrid()
          CALL DestroySpecToGrid()
          DO j=1,Jbmax
             DO l=1,Lmax
                DO i=1,Ibmaxperjb(j)
                   ug(i,l,j)=ug(i,l,j)*cosiv(i,j)
                END DO
             END DO
          END DO
          CALL rwrite (nvo(kflo), kflo, ug, ga, .true.)
          kflo=kflo+1
       END IF

       ! meridional wind chi at p levels

       IF (dopf(15,id)) THEN
          if(myid.eq.0) WRITE (UNIT=nfprt, FMT='(A)') ' postgg: vchi'
          DO j=1,Jbmax
             DO l=1,Lmax
                DO i=1,Ibmaxperjb(j)
                   vg(i,l,j)=vg(i,l,j)*cosiv(i,j)
                END DO
             END DO
          END DO
          CALL rwrite (nvo(kflo), kflo, vg, ga, .true.)
          kflo=kflo+1
       END IF


       ! it is necessary to call Heights in the case of
       ! absolute temperature (dopf(20,id)),
       ! relative humidity (dopf(22,id)),
       ! specific humidity (dopf(23,id)) and
       ! to save the virtual temperature in  og
       ! that will be used in getsh called at Humidity

       IF (dopf(16,id) .OR. dopf(17,id) .OR. dopf(18,id) .OR. &
           dopf(20,id) .OR. dopf(22,id) .OR. dopf(23,id)) THEN
          if(myid.eq.0) WRITE (UNIT=nfprt, FMT='(A)') ' Heights: Tv, zp and slp'
          CALL Heights (psmb, pmand, alnpmd, kflo, dopf, id, &
               nvo, ndp, top, ug, vg, og)
       END IF

       ! the above includes 16, 17, 18
       ! the one below includes 19, 21, 24
       ! both include 20, 22, 23
       IF (dopf(19,id) .OR. dopf(20,id) .OR. dopf(21,id) .OR. &
           dopf(22,id) .OR. dopf(23,id) .OR. dopf(24,id)) THEN
          if(myid.eq.0) WRITE (UNIT=nfprt, FMT='(A)') ' Humidity: Ts, Tp, rhs, rhp, sh, pw and theta'
          CALL Humidity (psmb, pmand, alnpmd, kflo, dopf, id, &
               nvo, ndp, ug, vg,rq, og)
       END IF

       DO j=1,Jbmax
          nCols=Ibmaxperjb(j)
          td2m(1:nCols,j)=CalcDewPoint(t2mt(1:nCols,j),q2mt(1:nCols,j),psmb(1:nCols,j),nCols)
       END DO



       ! ROUGHNESS LENGTH 
         IF (dopf(26,id)) THEN
            CALL scase (kflo, nfe, iclcd, z0rl)
            CALL rwrite (nvo(kflo), kflo, z0rl, qa, .true.)
            kflo=kflo+1
         END IF
       ! SURFACE TEMPERATURE
         IF (dopf(27,id)) THEN
            CALL scase (kflo, nfe, iclcd, tsfc)
            CALL rwrite (nvo(kflo), kflo, tsfc, qa, .true.)
            kflo=kflo+1
         END IF
       ! DEEP SOIL TEMPERATURE
         IF (dopf(28,id)) THEN
            CALL scase (kflo, nfe, iclcd, dest)
            CALL rwrite (nvo(kflo), kflo, dest, qa, .true.)
            kflo=kflo+1
         END IF
       ! STORAGE ON CANOPY
         IF (dopf(29,id)) THEN
            CALL scase (kflo, nfe, iclcd, stcp)
            CALL rwrite (nvo(kflo), kflo, stcp, qa, .true.)
            kflo=kflo+1
         END IF
       ! STORAGE ON GROUND 
         IF (dopf(30,id)) THEN
            CALL scase (kflo, nfe, iclcd, stso)
            CALL rwrite (nvo(kflo), kflo, stso, qa, .true.)
            kflo=kflo+1
         END IF
       ! SOIL WETNESS OF SURFACE 
         IF (dopf(31,id)) THEN
            CALL scase (kflo, nfe, iclcd, wssf)
            CALL rwrite (nvo(kflo), kflo, wssf, qa, .true.)
            kflo=kflo+1
         END IF
       ! SOIL WETNESS OF ROOT ZONE 
         IF (dopf(32,id)) THEN
            CALL scase (kflo, nfe, iclcd, wrzo)
            CALL rwrite (nvo(kflo), kflo, wrzo, qa, .true.)
            kflo=kflo+1
         END IF
       ! SOIL WETNESS OF DRAINAGE ZONE
         IF (dopf(33,id)) THEN
            CALL scase (kflo, nfe, iclcd, wdzo)
            CALL rwrite (nvo(kflo), kflo, wdzo, qa, .true.)
            kflo=kflo+1
         END IF 
       ! TEMPERATURE AT 2-M FROM SURFACE 
         IF (dopf(34,id)) THEN
            CALL scase (kflo, nfe, iclcd, t2mt)
            CALL rwrite (nvo(kflo), kflo, t2mt, qa, .true.)
            kflo=kflo+1
         END IF
         kp=0
         DO i=1,nfld
            IF (prodia(i) == prog)kp=kp+1
         ENDDO
         IF(kp == 0)THEN
            ! 2 METRE DEWPOINT TEMPERATURE 
            IF (dopf(35,id)) THEN
               CALL rwrite (nvo(kflo), kflo, td2m, qa, .true.)
               kflo=kflo+1
            END IF
         END IF
       ! SPECIFIC HUMIDITY AT 2-M FROM SURFACE
         IF (dopf(36,id)) THEN
            CALL rwrite (nvo(kflo), kflo, q2mt, qa, .true.)
            kflo=kflo+1
         END IF

       newday=.FALSE.

    END IF

  END SUBROUTINE postgg


  SUBROUTINE postgg2 (psmb, pmand, alnpmd, mean, kfld, kflo, newday, top)

    IMPLICIT NONE

    INTEGER, INTENT(INOUT) :: kfld
    INTEGER, INTENT(INOUT) :: kflo

    REAL (KIND=r8), INTENT(IN) :: pmand(Lmax)
    REAL (KIND=r8), INTENT(IN) :: alnpmd(Lmax)
    REAL (KIND=r8), INTENT(IN) :: top(Ibmax,Jbmax)
    REAL (KIND=r8), INTENT(OUT) :: psmb(Ibmax,Jbmax)

    LOGICAL, INTENT(IN) :: mean
    LOGICAL, INTENT(INOUT) :: newday

    INTEGER, PARAMETER :: nopf=36  !PK
    INTEGER :: id
    INTEGER :: i
    INTEGER :: kt
    INTEGER :: mm
    INTEGER :: ii
    INTEGER :: j
    INTEGER :: ib,jb
    INTEGER :: m, mg, mglob, n, mn, mf
    INTEGER :: l
    INTEGER :: ll
    INTEGER :: nn
    INTEGER :: nmax,nCols

    REAL (KIND=r8) :: ug(Ibmax,Lmax,Jbmax)
    REAL (KIND=r8) :: vg(Ibmax,Lmax,Jbmax)
    REAL (KIND=r8) :: rq(Ibmax,Lmax,Jbmax)
    REAL (KIND=r8) :: og(Ibmax,Lmax,Jbmax)
    REAL (KIND=r8) :: root2  

    LOGICAL :: timemean
    !    LOGICAL, SAVE :: first=.TRUE.
    !    LOGICAL, SAVE :: second=.FALSE.
    !    LOGICAL, SAVE :: dopf(nopf,2)

    !    CHARACTER (LEN=40), DIMENSION (5,2), SAVE :: inf
    !    CHARACTER (LEN=40), DIMENSION (nopf,2), SAVE :: opf

    !    CHARACTER (LEN=40), SAVE :: t2ps='TIME MEAN SURFACE PRESSURE              '

    inf(:,1) = (/ &
         'LN SURFACE PRESSURE                     ', &
         'DIVERGENCE                              ', &
         'VORTICITY                               ', &
         'SPECIFIC HUMIDITY                       ', &
         'VIRTUAL TEMPERATURE                     ' /)
    inf(:,2) = (/ &
         'TIME MEAN LN SURFACE PRESSURE           ', &
         'TIME MEAN DIVERGENCE                    ', &
         'TIME MEAN VORTICITY                     ', &
         'TIME MEAN SPECIFIC HUMIDITY             ', &
         'TIME MEAN VIRTUAL TEMPERATURE           ' /)

    opf(:,1) = (/ &
         'SURFACE PRESSURE                        ', &! 1
         'MASK                                    ', &! 2
         'SURFACE ZONAL WIND (U)                  ', &! 3
         'ZONAL WIND (U)                          ', &! 4
         'SURFACE MERIDIONAL WIND (V)             ', &! 5
         'MERIDIONAL WIND (V)                     ', &! 6
         'OMEGA                                   ', &! 7
         'DIVERGENCE                              ', &! 8
         'VORTICITY                               ', &! 9
         'STREAM FUNCTION                         ', &!10
         'ZONAL WIND PSI                          ', &!11
         'MERIDIONAL WIND PSI                     ', &!12
         'VELOCITY POTENTIAL                      ', &!13
         'ZONAL WIND CHI                          ', &!14
         'MERIDIONAL WIND CHI                     ', &!15
         'VIRTUAL TEMPERATURE                     ', &!16
         'GEOPOTENTIAL HEIGHT                     ', &!17
         'SEA LEVEL PRESSURE                      ', &!18
         'SURFACE ABSOLUTE TEMPERATURE            ', &!19
         'ABSOLUTE TEMPERATURE                    ', &!20
         'SURFACE RELATIVE HUMIDITY               ', &!21
         'RELATIVE HUMIDITY                       ', &!22
         'SPECIFIC HUMIDITY                       ', &!23
         'INST. PRECIP. WATER                     ', &!24
         'POTENTIAL TEMPERATURE                   ', &!25
         'ROUGHNESS LENGTH                        ', &!26 !PK---
         'SURFACE TEMPERATURE                     ', &!27 
         'DEEP SOIL TEMPERATURE                   ', &!28 
         'STORAGE ON CANOPY                       ', &!29 
         'STORAGE ON GROUND                       ', &!30 
         'SOIL WETNESS OF SURFACE                 ', &!31 
         'SOIL WETNESS OF ROOT ZONE               ', &!32 
         'SOIL WETNESS OF DRAINAGE ZONE           ', &!33 
         'TEMPERATURE AT 2-M FROM SURFACE         ', &!34 
         'SPECIFIC HUMIDITY AT 2-M FROM SURFACE   ', &!36 !PK---
         '2 METRE DEWPOINT TEMPERATURE            '/) !36 !PK---

    opf(:,2) = (/ &
         'TIME MEAN SURFACE PRESSURE              ', &! 1
         'TIME MEAN MASK                          ', &! 2
         'TIME MEAN SURFACE ZONAL WIND (U)        ', &! 3
         'TIME MEAN ZONAL WIND (U)                ', &! 4
         'TIME MEAN SURFACE MERIDIONAL WIND (V)   ', &! 5
         'TIME MEAN MERIDIONAL WIND (V)           ', &! 6
         'TIME MEAN DERIVED OMEGA                 ', &! 7
         'TIME MEAN DIVERGENCE                    ', &! 8
         'TIME MEAN VORTICITY                     ', &! 9
         'TIME MEAN STREAM FUNCTION               ', &!10
         'TIME MEAN ZONAL WIND PSI                ', &!11
         'TIME MEAN MERIDIONAL WIND PSI           ', &!12
         'TIME MEAN VELOCITY POTENTIAL            ', &!13
         'TIME MEAN ZONAL WIND CHI                ', &!14
         'TIME MEAN MERIDIONAL WIND CHI           ', &!15
         'TIME MEAN VIRTUAL TEMPERATURE           ', &!16
         'TIME MEAN GEOPOTENTIAL HEIGHT           ', &!17
         'TIME MEAN SEA LEVEL PRESSURE            ', &!18
         'TIME MEAN SURFACE ABSOLUTE TEMPERATURE  ', &!19
         'TIME MEAN ABSOLUTE TEMPERATURE          ', &!20
         'TIME MEAN SURFACE RELATIVE HUMIDITY     ', &!21
         'TIME MEAN RELATIVE HUMIDITY             ', &!22
         'TIME MEAN SPECIFIC HUMIDITY             ', &!23
         'TIME MEAN PRECIP. WATER                 ', &!24
         'TIME MEAN POTENTIAL TEMPERATURE         ', &!25
         'TIME MEAN_ROUGHNESS LENGTH              ', &!26 !PK---
         'TIME MEAN_SURFACE TEMPERATURE           ', &!27
         'TIME MEAN_DEEP SOIL TEMPERATURE         ', &!28
         'TIME MEAN_STORAGE ON CANOPY             ', &!29
         'TIME MEAN_STORAGE ON GROUND             ', &!30
         'TIME MEAN_SOIL WETNESS OF SURFACE       ', &!31
         'TIME MEAN_SOIL WETNESS OF ROOT ZONE     ', &!32
         'TIME MEAN_SOIL WETNESS OF DRAINAGE ZONE ', &!33
         'TIME MEAN_TEMPERATURE AT 2M FROM SURFACE', &!34
         'TIME MEAN_SPEC-HUMID AT 2-M FROM SURFACE', &!36 !PK
         'TIME MEAN 2 METRE DEWPOINT TEMPERATURE  '/) !36 !PK---

    timemean=.TRUE.
    IF (first) THEN

       IF (mean) THEN
          IF (myid.eq.0) THEN
             WRITE (UNIT=nferr, FMT='(A3L4)') &
                ' First, Second, Mean Inconsistent:',&
                  first, second2, mean
          END IF
          STOP 8100
       END IF
       id=1

       ! input field has to be complete and well ordered,

       kt=kfld
       DO mm=1,5
          IF (inf(mm,1) == chrdsc(kt)) THEN
             kt=kt+1
          ELSE
             IF (myid.eq.0) THEN
                WRITE (UNIT=nferr, FMT='(A,L2,2(A,I4),A40,/,40X,A40)') &
                   ' At Mean = ', mean, ' mm = ', mm, ' kt = ', kt, &
                   ' inf = ', inf(mm,id), ' chrdsc = ', chrdsc(kt)
             END IF
             STOP 9100
          END IF
       END DO

       ! output field has to be well ordered.
       ! ignore undesired output fields.

       kt=kflo
       DO ii=1,nopf
            IF (myid.eq.0) WRITE (UNIT=nfprt, FMT='(" opf(ii,1) = ",A)') opf(ii,1)
            IF (myid.eq.0) WRITE (UNIT=nfprt, FMT='(" chrop(kt) = ",A)') chrop(kt)

          IF (opf(ii,1) == chrop(kt)) THEN
!            IF (myid.eq.0) WRITE (UNIT=nfprt, FMT='(" opf(ii,1) = ",A)') opf(ii,1)
!            IF (myid.eq.0) WRITE (UNIT=nfprt, FMT='(" chrop(kt) = ",A)') chrop(kt)
             dopf(ii,1)=.TRUE.
             kt=kt+1
          ELSE
             dopf(ii,1)=.FALSE.
          END IF
       END DO
       first=.FALSE.
       second2=.TRUE.
       IF (myid.eq.0) &
       WRITE (UNIT=nfprt, FMT='(" dopf(,1) = ",50L2)') (dopf(i,1),i=1,nopf)

    ELSE IF (second2) THEN
       IF (.NOT. newday) THEN
          IF (.NOT.mean) THEN
             IF (myid.eq.0) &
             WRITE (UNIT=nferr, FMT='(A3L4)') &
                   ' First, Second, Mean Inconsistent:',&
                     first, second2, mean
             STOP 8100
          END IF

          id=2
          inf(1,2)=t2ps

          ! check that the time mean input field is complete and well ordered

          kt=kfld
          DO mm=1,5
             IF (inf(mm,2) == chrdsc(kt)) THEN
                kt=kt+1
             ELSE
                IF (myid.eq.0) &
                WRITE (UNIT=nfprt, FMT='(A)') ' Time Mean is not Available, Ignore it ...'
                timemean=.FALSE.
             END IF
          END DO

          ! check that the time mean output field is well ordered.
          ! ignore undesired output fields.

          kt=kflo
          DO ii=1,nopf
             IF (opf(ii,2) == chrop(kt)) THEN
                dopf(ii,2)=.TRUE.
                kt=kt+1
             ELSE
                dopf(ii,2)=.FALSE.
             END IF
          END DO
          IF (myid.eq.0) &
          WRITE (UNIT=nfprt, FMT='(" dopf(,2) = ",50L2)') (dopf(i,2),i=1,nopf)

       END IF
    END IF
    root2  = SQRT(2.0_r8)
    IF (timemean) THEN

       ! surface pressure

       IF (.NOT.mean) THEN

          ! case full field:
          ! read spectral ln surface pressure;
          ! legandre transform to grid;
          ! get surface pressure and convert from centibar to milibar

          IF (myid.eq.0) &
          WRITE (nfprt,*) ' lendo qlnpp spectral '
          id=1
          CALL ReadFieldSp (nharm(kfld), nlevs(kfld), qlnppPA)
          CALL Rectrg(qlnppPA, lnpsPa, 1, 1)
          qlnppmb=qlnppPA
          qlnpp=qlnppPA
          !      transform surface pressure from Pascal to cbar
          qlnpp(1) = qlnpp(1) - log(1000._r8) * root2
          !      transform surface pressure from Pascal to mbar
          qlnppmb(1) = qlnppmb(1) - log(100._r8) * root2
       ELSE

          ! case time mean:
          ! read grid surface pressure;
          ! convert from centibar to milibar and get log
          ! legendre transform from spectral to grid

          id=2
          pscb = 0.
          CALL ReadFieldG(nharm(kfld), nlevs(kfld), pscb)
          DO j=1,Jbmax
             DO i=1,Ibmaxperjb(j)
                psmb(i,j)=pscb(i,j)*10.0_r8
                lnpscb(i,j)=LOG(pscb(i,j))
             END DO
          END DO
 !PK         CALL Rectrg(qlnpp, lnpscb, 1, 1)
          CALL Rectrg(qlnppPA, lnpsPa, 1, 1)
          !      transform surface pressure from Pascal to cbar
          DO j=1,Jbmax
             DO i=1,Ibmaxperjb(j)
                 lnpscb(i,j) = EXP(lnpsPa(i,j))/1000.0_r8
                 lnpscb(i,j) = log(lnpscb(i,j))
             END DO
          END DO

       END IF

       ! divergence
       kfld=kfld+1
       CALL ReadFieldSp (nharm(kfld), nlevs(kfld), qdivp)
       ! vorticity
       kfld=kfld+1
       CALL ReadFieldSp (nharm(kfld), nlevs(kfld), qrotp)
       ! humidity 
       kfld=kfld+1
       CALL ReadFieldSp (nharm(kfld), nlevs(kfld), qqp)
       ! temperature
       kfld=kfld+1
       CALL ReadFieldSp (nharm(kfld), nlevs(kfld), qtmpp)
       kfld=kfld+1
!
       !  Transform all fields to grid space (kmax levels)
!
       IF (havesurf) CALL gozrim(qlnpp, qpphi, 1, 2*mymnExtMax)
       IF (havesurf) CALL gozrim(qlnppPa, qpphiPa, 1, 2*mymnExtMax)

       CALL dztouv(qdivp, qrotp, qup, qvp, 1, 2*mymnExtMax)
!
       CALL CreateSpecToGrid(5, 4, 5, 6)
       CALL DepositSpecToGridAndDelLamGrid(qlnpp, lnpscb, fgplam)
       CALL DepositSpecToGridAndDelLamGrid(qlnppPa, lnpsPa, fgplamPa)
!
       CALL DepositSpecToGrid(qpphi, fgpphi)
       CALL DepositSpecToGrid(qpphiPa, fgpphiPa)
!
       CALL DepositSpecToGrid(qup,fgu)
       CALL DepositSpecToGrid(qvp,fgv)
       CALL DepositSpecToGrid(qdivp, fgdiv)
       CALL DepositSpecToGrid(qtmpp, fgtmp)
       CALL DepositSpecToGrid(qqp, fgq)
!
       CALL DoSpecToGrid()
       CALL DestroySpecToGrid()
!
       IF (.NOT.mean) THEN
          DO j=1,Jbmax
             DO i=1,Ibmaxperjb(j)
!PK                pscb(i,j)=EXP(lnpscb(i,j))/100.0_r8
                pscb(i,j)=EXP(lnpsPa(i,j))/1000.0_r8
                psmb(i,j)=pscb(i,j)*10.0_r8
             END DO
          END DO
       ENDIF
       IF (dopf(1,id) .OR. dopf(2,id)) THEN
       IF(myid.eq.0) WRITE (UNIT=nfprt, FMT='(A)') 'postgg: ps and masK'

          ! generates mask array containing
          ! ones above surface and -ones below surface.

          DO j=1,Jbmax
             DO l=1,Lmax
                DO i=1,Ibmaxperjb(j)
                   IF (psmb(i,j) < pmand(l)) THEN
                      og(i,l,j)=-pmand(l)/psmb(i,j)
                   ELSE
                      og(i,l,j)=pmand(l)/psmb(i,j)
                   END IF
                END DO
             END DO
          END DO

          ! surface pressure

          IF (dopf(1,id)) THEN
             CALL rwrite (nvo(kflo), kflo, pscb, qa, .true.)
             kflo=kflo+1
          END IF

          ! mask for vertical extrapolation inside terrain

          IF (dopf(2,id)) THEN
             CALL rwrite (nvo(kflo), kflo, og, ga, .true.)
             kflo=kflo+1
          END IF

       END IF

       IF (dopf(3,id) .OR. dopf(4,id) .OR. dopf(5,id) .OR. &
            dopf(6,id) .OR. dopf(7,id) .OR. dopf(8,id) .OR. &
            dopf(10,id) .OR. dopf(11,id) .OR. dopf(12,id) .OR. &
            dopf(13,id) .OR. dopf(14,id) .OR. dopf(15,id)) THEN
          IF (myid.eq.0) WRITE (UNIT=nfprt, FMT='(A)') ' Winds: us, up, vs, vp and omega'
          CALL Winds (psmb, alnpmd, kflo, dopf, id, nvo, ndp, &
                      ug, vg, og)
       END IF

       ! this is a subset of the above condition with the exception of
       ! dopf(9,id)
       IF (dopf(8,id) .OR. dopf(9,id) .OR. &
            dopf(10,id) .OR. dopf(11,id) .OR. dopf(12,id) .OR. &
            dopf(13,id) .OR. dopf(14,id) .OR. dopf(15,id)) THEN
          IF (myid.eq.0) WRITE (UNIT=nfprt, FMT='(A)') ' uvtodz: div and vort'
          ! Recompute Divergence and Vorticy using U an V at p-levels.
          CALL ReshapeVerticalGroups(lmax, kmaxloc)
          CALL CreateGridToSpec(2, 0)
          CALL DepositGridToSpec(qug, ug)
          CALL DepositGridToSpec(qvg, vg)
          CALL DoGridToSpec()
          CALL DestroyGridToSpec()
          CALL Uvtodz(qug, qvg, qdivg, qrotg, 1, 2*mymnmax)
          CALL CreateSpecToGrid(2, 0, 2, 0)
          CALL DepositSpecToGrid(qdivg, ug)
          CALL DepositSpecToGrid(qrotg, vg)
          CALL DoSpecToGrid()
          CALL DestroySpecToGrid()
       END IF

       ! horizontal divergence at p levels

       IF (dopf(8,id)) THEN
          CALL rwrite (nvo(kflo), kflo, ug, ga, .true.)
          kflo=kflo+1
       END IF

       ! vertical componente of vorticity at p levels

       IF (dopf(9,id)) THEN
          CALL rwrite (nvo(kflo), kflo, vg, ga, .true.)
          kflo=kflo+1
       END IF

       ! stream function at p levels

       IF (dopf(10,id) .OR. dopf(11,id) .OR. dopf(12,id)) THEN
          DO l=1,Lmaxloc
             mf = 1
             IF (haveM1) THEN 
                wk(1,l)=0.0_r8
                wk(2,l)=0.0_r8
                mf = 3
             ENDIF
             DO m=mf,2*mymnmax
                wk(m,l)=qrotg(m,l)*snnp1(m)
             END DO
          END DO
          IF (dopf(10,id)) THEN
             if(myid.eq.0) WRITE (UNIT=nfprt, FMT='(A)') ' postgg: psi'
             CALL rectrg (wk, ug, Lmaxloc, Lmax)
             CALL rwrite (nvo(kflo), kflo, ug, ga, .true.)
             kflo=kflo+1
          END IF
       END IF

       ! zonal wind psi at p levels

       IF (dopf(11,id)) THEN
          if(myid.eq.0) WRITE (UNIT=nfprt, FMT='(A)') ' postgg: upsi'
          CALL CreateSpecToGrid(2, 0, 2, 0)
          CALL DepositSpecToDelLamGrid(wk,vg)
          CALL gozrim(wk, wlphi, 1, 2*mymnextmax)
          CALL DepositSpecToGrid(wlphi, ug)
          CALL DoSpecToGrid()
          CALL DestroySpecToGrid()
          DO j=1,Jbmax
             DO l=1,Lmax
                DO i=1,Ibmaxperjb(j)
                   ug(i,l,j)=-ug(i,l,j)*cosiv(i,j)
                END DO
             END DO
          END DO
          CALL rwrite (nvo(kflo), kflo, ug, ga, .true.)
          kflo=kflo+1
       END IF

       ! meridional wind psi at p levels

       IF (dopf(12,id)) THEN
          if(myid.eq.0) WRITE (UNIT=nfprt, FMT='(A)') ' postgg: vpsi'
          DO j=1,Jbmax
             DO l=1,Lmax
                DO i=1,Ibmaxperjb(j)
                   vg(i,l,j)=vg(i,l,j)*cosiv(i,j)
                END DO
             END DO
          END DO
          CALL rwrite (nvo(kflo), kflo, vg, ga, .true.)
          kflo=kflo+1
       END IF

       ! velocity potential at p levels

       IF (dopf(13,id) .OR. dopf(14,id) .OR. dopf(15,id)) THEN
          DO l=1,Lmaxloc
             mf = 1
             IF (haveM1) THEN 
                wk(1,l)=0.0_r8
                wk(2,l)=0.0_r8
                mf = 3
             ENDIF
             DO m=mf,2*mymnmax
                wk(m,l)=qdivg(m,l)*snnp1(m)
             END DO
          END DO
          IF (dopf(13,id)) THEN
             if(myid.eq.0) WRITE (UNIT=nfprt, FMT='(A)') ' postgg: chi'
             CALL rectrg (wk, vg, Lmaxloc, Lmax)
             CALL rwrite (nvo(kflo), kflo, vg, ga, .true.)
             kflo=kflo+1
          END IF
       END IF

       ! zonal wind chi at p levels

       IF (dopf(14,id)) THEN
          if(myid.eq.0) WRITE (UNIT=nfprt, FMT='(A)') ' postgg: uchi'
          CALL CreateSpecToGrid(2, 0, 2, 0)
          CALL DepositSpecToDelLamGrid(wk,ug)
          CALL gozrim(wk, wlphi, 1, 2*mymnextmax)
          CALL DepositSpecToGrid(wlphi, vg)
          CALL DoSpecToGrid()
          CALL DestroySpecToGrid()
          DO j=1,Jbmax
             DO l=1,Lmax
                DO i=1,Ibmaxperjb(j)
                   ug(i,l,j)=ug(i,l,j)*cosiv(i,j)
                END DO
             END DO
          END DO
          CALL rwrite (nvo(kflo), kflo, ug, ga, .true.)
          kflo=kflo+1
       END IF

       ! meridional wind chi at p levels

       IF (dopf(15,id)) THEN
          if(myid.eq.0) WRITE (UNIT=nfprt, FMT='(A)') ' postgg: vchi'
          DO j=1,Jbmax
             DO l=1,Lmax
                DO i=1,Ibmaxperjb(j)
                   vg(i,l,j)=vg(i,l,j)*cosiv(i,j)
                END DO
             END DO
          END DO
          CALL rwrite (nvo(kflo), kflo, vg, ga, .true.)
          kflo=kflo+1
       END IF


       ! it is necessary to call Heights in the case of
       ! absolute temperature (dopf(20,id)),
       ! relative humidity (dopf(22,id)),
       ! specific humidity (dopf(23,id)) and
       ! to save the virtual temperature in  og
       ! that will be used in getsh called at Humidity

       IF (dopf(16,id) .OR. dopf(17,id) .OR. dopf(18,id) .OR. &
           dopf(20,id) .OR. dopf(22,id) .OR. dopf(23,id)) THEN
          if(myid.eq.0) WRITE (UNIT=nfprt, FMT='(A)') ' Heights: Tv, zp and slp'
          CALL Heights (psmb, pmand, alnpmd, kflo, dopf, id, &
               nvo, ndp, top, ug, vg, og)
       END IF

       ! the above includes 16, 17, 18
       ! the one below includes 19, 21, 24
       ! both include 20, 22, 23
       IF (dopf(19,id) .OR. dopf(20,id) .OR. dopf(21,id) .OR. &
           dopf(22,id) .OR. dopf(23,id) .OR. dopf(24,id)) THEN
          if(myid.eq.0) WRITE (UNIT=nfprt, FMT='(A)') ' Humidity: Ts, Tp, rhs, rhp, sh, pw and theta'
          CALL Humidity (psmb, pmand, alnpmd, kflo, dopf, id, &
               nvo, ndp, ug, vg,rq, og)
       END IF

       newday=.FALSE.

    END IF

  END SUBROUTINE postgg2
  
  
  SUBROUTINE Winds (psmb, alnpmd, kflo, dopf, id, nvo, ndp, ug, vg, og)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ndp
    INTEGER, INTENT(IN) :: id
    INTEGER, INTENT(INOUT) :: kflo
    INTEGER, INTENT(IN) :: nvo(ndp)

    REAL (KIND=r8), INTENT(IN) :: psmb(Ibmax,Jbmax)
    REAL (KIND=r8), INTENT(IN) :: alnpmd(Lmax)
    REAL (KIND=r8), INTENT(OUT) :: ug(Ibmax,Lmax,Jbmax)
    REAL (KIND=r8), INTENT(OUT) :: vg(Ibmax,Lmax,Jbmax)
    REAL (KIND=r8), INTENT(OUT) :: og(Ibmax,Lmax,Jbmax)

    REAL (KIND=r8) :: us(Ibmax,Jbmax)
    REAL (KIND=r8) :: vs(Ibmax,Jbmax) 

    LOGICAL, INTENT(IN) :: dopf(:,:)

    INTEGER :: j
    INTEGER :: l
    INTEGER :: i

    !cdir nodep
    DO j=1,jbmax
       DO i=1, ibmaxperjb(j)
          us(i,j)=fgu(i,1,j)/coslat(i,j)
          vs(i,j)=fgv(i,1,j)/coslat(i,j)
       END DO
    END DO

    ! omegas calculates vertical motions.

    CALL omegas (fgpphiPa, fgplamPa, fgu, fgv, fgdiv, rcl, fgomega, psmb)



    ! sig2po performs vertical interpolation

    !ug,vg and og SET
    CALL sig2po (psmb, alnpmd, fgu, ug, fgv, vg, fgomega, og)

    ! convert pseudo-wind to wind field

    DO j=1,Jbmax
       DO l=1,Lmax
          DO i=1,Ibmaxperjb(j)
             ug(i,l,j)=ug(i,l,j)/coslat(i,j)
             vg(i,l,j)=vg(i,l,j)/coslat(i,j)
          END DO
       END DO
    END DO

    ! zonal wind at first sigma layer

    IF (dopf(3,id)) THEN
       CALL rwrite (nvo(kflo), kflo, us, ga, .true.)
       kflo=kflo+1
    END IF

    ! zonal wind at p levels

    IF (dopf(4,id)) THEN
       CALL rwrite (nvo(kflo), kflo, ug, ga, .true.)
       kflo=kflo+1
    END IF

    ! meridional wind at first sigma layer

    IF (dopf(5,id)) THEN
       CALL rwrite (nvo(kflo), kflo, vs, ga, .true.)
       kflo=kflo+1
    END IF

    ! meridional wind at p levels

    IF (dopf(6,id)) THEN
       CALL rwrite (nvo(kflo), kflo, vg, ga, .true.)
       kflo=kflo+1
    END IF

    ! vertical p-velocity at p levels

    IF (dopf(7,id)) THEN
       CALL rwrite (nvo(kflo), kflo, og, ga, .true.)
       kflo=kflo+1
    END IF

    ! prepare ug and vg for computing divergence and vorticity

    DO j=1,Jbmax
       DO l=1,Lmax
          DO i=1,Ibmaxperjb(j)
             ug(i,l,j)=ug(i,l,j)/coslat(i,j)
             vg(i,l,j)=vg(i,l,j)/coslat(i,j)
          END DO
       END DO
    END DO

  END SUBROUTINE Winds


  SUBROUTINE Heights (psmb, pmand, alnpmd, kflo, dopf, id, &
       nvo, ndp, top, zp, tv, tvsav)


    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ndp
    INTEGER, INTENT(IN) :: id
    INTEGER, INTENT(INOUT) :: kflo
    INTEGER, INTENT(IN) :: nvo(ndp)

    REAL (KIND=r8), INTENT(IN) :: psmb(Ibmax,Jbmax)
    REAL (KIND=r8), INTENT(IN) :: pmand(Lmax)
    REAL (KIND=r8), INTENT(IN) :: alnpmd(Lmax)
    REAL (KIND=r8), INTENT(IN) :: top(Ibmax,Jbmax)
    REAL (KIND=r8), INTENT(OUT) :: zp(Ibmax,Lmax,Jbmax)
    REAL (KIND=r8), INTENT(OUT) :: tv(Ibmax,Lmax,Jbmax)
    REAL (KIND=r8), INTENT(OUT) :: tvsav(Ibmax,Lmax,Jbmax)

    LOGICAL, INTENT(IN) :: dopf(:,:)

    INTEGER :: i, j, l

    REAL (KIND=r8) :: tc(Ibmax,Lmax,Jbmax)
    REAL (KIND=r8) :: ts(Ibmax,Lmax,Jbmax)
    REAL (KIND=r8) :: slpg(Ibmax,Jbmax)

    zp=0.0_r8
    tv=0.0_r8
    slpg=0.0_r8

    ! legendre transform from spectral to fourier
    ! fourier transform from fourier to grid


    CALL sig2pz (fgtmp, psmb, top, tv, zp, alnpmd)

    ! calculation for below ground virtual temp. based on height
    ! and save virtual temperature at p levels for next routine

    CALL lowtmp (zp, tc, pmand)
    DO j=1,Jbmax
       DO l=1,Lmax
          DO i=1,Ibmaxperjb(j)
             IF (pmand(l) > psmb(i,j)) tv(i,l,j)=tc(i,l,j)
             tvsav(i,l,j)=tv(i,l,j)
          ENDDO
       ENDDO
    ENDDO

    ! virtual temperature at p levels

    IF (dopf(16,id)) THEN
       CALL rwrite (nvo(kflo), kflo, tv, ga, .true.)
       kflo=kflo+1
    ENDIF

    ! geopotential height at p levels

    IF (dopf(17,id)) THEN
       CALL rwrite (nvo(kflo), kflo, zp, ga, .true.)
       kflo=kflo+1
    ENDIF

    ! reduced sea level pressure

    IF (dopf(18,id)) THEN
       CALL getslp (slpg, zp, pmand)
       CALL rwrite (nvo(kflo), kflo, slpg, ga, .true.)
       kflo=kflo+1
    ENDIF

  END SUBROUTINE Heights


  SUBROUTINE Humidity (psmb, pmand, alnpmd, kflo, dopf, id, &
       nvo, ndp, ta, rh,rq, tv)


    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ndp
    INTEGER, INTENT(IN) :: id
    INTEGER, INTENT(INOUT) :: kflo
    INTEGER, INTENT(IN) :: nvo(ndp)

    REAL (KIND=r8), INTENT(IN) :: psmb(Ibmax,Jbmax)
    REAL (KIND=r8), INTENT(IN) :: pmand(Lmax)
    REAL (KIND=r8), INTENT(IN) :: alnpmd(Lmax)
    REAL (KIND=r8), INTENT(OUT) :: ta(Ibmax,Lmax,Jbmax)
    REAL (KIND=r8), INTENT(OUT) :: rh(Ibmax,Lmax,Jbmax)
    REAL (KIND=r8), INTENT(OUT) :: rq(Ibmax,Lmax,Jbmax)
    REAL (KIND=r8), INTENT(INOUT) :: tv(Ibmax,Lmax,Jbmax)

    LOGICAL, INTENT(IN) :: dopf(:,:)

    INTEGER :: j, i

    REAL (KIND=r8) :: gta(Ibmax,Kmax,Jbmax)
    REAL (KIND=r8) :: gss(Ibmax,Kmax,Jbmax)
    REAL (KIND=r8) :: ts(Ibmax,Jbmax)
    REAL (KIND=r8) :: rs(Ibmax,Jbmax)
    REAL (KIND=r8) :: pw(Ibmax,Jbmax)

    ta=0.0_r8
    rh=0.0_r8
    ts=0.0_r8
    rs=0.0_r8
    pw=0.0_r8

    DO j=1,jbMax
       DO i = 1,ibMaxperjb(j)
          fgq(i,:,j)=MAX(fgq(i,:,j),SHmin)
       ENDDO
    ENDDO

    ! compute precipitable water from specific humidity

    CALL pwater (fgq, pw, psmb)

    !     vertical interpolation.

    CALL sigtop (gta, fgtmp,fgq, gss, psmb, ta, rh,rq, pmand, alnpmd)

    DO j=1,Jbmax
       DO i=1,Ibmaxperjb(j)
          ts(i,j)=gta(i,1,j)
          rs(i,j)=gss(i,1,j)
       END DO
    END DO

    ! getsh calculates specific humidity from relative humidity
    ! and virtual temperature saved in tv (tvsav) at subroutine Heights.

    IF (dopf(20,id) .OR. dopf(22,id) .OR. dopf(23,id)) THEN
       CALL getsh (tv, rh,rq, pmand, ta, psmb)
    END IF

    ! absolute temperature of the first sigma layer

    IF (dopf(19,id)) THEN
       CALL rwrite (nvo(kflo), kflo, ts, ga, .true.)
       kflo=kflo+1
    END IF

    ! absolute temperature at p levels

    IF (dopf(20,id)) THEN
       CALL rwrite (nvo(kflo), kflo, ta, ga, .true.)
       kflo=kflo+1
    END IF

    ! relative humidity at first sigma layer

    IF (dopf(21,id)) THEN
       CALL rwrite (nvo(kflo), kflo, rs, ga, .true.)
       kflo=kflo+1
    END IF

    ! relative humidity at p levels

    IF (dopf(22,id)) THEN
       CALL rwrite (nvo(kflo), kflo, rh, ga, .true.)
       kflo=kflo+1
    END IF

    ! specific humidity at p levels

    IF (dopf(23,id)) THEN
       CALL rwrite (nvo(kflo), kflo, tv, ga, .true.)
       kflo=kflo+1
    END IF

    ! precipitable water integrated over sigma layers

    IF (dopf(24,id)) THEN
       CALL rwrite (nvo(kflo), kflo, pw, ga, .true.)
       kflo=kflo+1
    END IF

    ! potential temperature at p levels

    IF (dopf(25,id)) THEN
       CALL getth (tv, pmand, ta)
       CALL rwrite (nvo(kflo), kflo, tv, ga, .true.)
       kflo=kflo+1
    END IF

  END SUBROUTINE Humidity


  SUBROUTINE stog (psmb, alnpmd, Ldim, mean, kflo, nof, kfld, di, og, qg)

    IMPLICIT NONE

    ! general purpose single field spectral sigma to regular pressure
    ! conversion subroutine.  input array assumed to be in 1st slot

    INTEGER, INTENT(IN) :: Ldim
    INTEGER, INTENT(IN) :: kflo
    INTEGER, INTENT(IN) :: nof
    INTEGER, INTENT(IN) :: kfld

    REAL (KIND=r8), INTENT(IN) :: psmb(Ibmax,Jbmax)
    REAL (KIND=r8), INTENT(IN) :: alnpmd(Lmax)
    REAL (KIND=r8), INTENT(IN) :: di(2*mymnmax,Kmax)
    REAL (KIND=r8), INTENT(OUT) :: og(Ibmax,Lmax,Jbmax)
    REAL (KIND=r8), INTENT(OUT) :: qg(Ibmax,Jbmax)

    LOGICAL, INTENT(IN) :: mean

    INTEGER :: lout, j, l, i, ijKmax
    INTEGER :: jjj, mmm, nnn, mm1, nn1, jj1, nnx, k, ier

    IF (Ldim == 1) THEN
       lout=1
    ELSE
       lout=Lmax
    END IF

    qg=0.0_r8
    og(:,1:lout,:)=0.0_r8

    jjj=jfe(kflo)
    mmm=mko(jjj)
    nnn=ife(jfe(kflo))
    mm1=0
    nn1=1
    IF (mmm == 1 .AND. kflo < nof) THEN
       jj1=jfe(kflo+1)
       mm1=mko(jj1)
       nn1=ife(jj1)
       IF (kfld /= lif(nn1)) mm1=0
    END IF

    IF (lout.eq.1) THEN
       CALL CreateSpecToGrid(0, 1, 0, 1)
       CALL DepositSpecToGrid(di(:,1), fgpphi)
     ELSE
       CALL CreateSpecToGrid(1, 0, 1, 0)
       CALL DepositSpecToGrid(di, fgdiv)
    ENDIF
    CALL DoSpecToGrid()
    CALL DestroySpecToGrid()

    IF (Ldim==Kmax .AND. mean) THEN

       IF (  (mmm==2 .AND. iclcd(nnn)==3)  .OR.&
            (mm1==2 .AND. iclcd(nn1)==3)) THEN
          IF (mmm==2) nnx=nnn
          IF (mm1==2) nnx=nn1
          IF (ifldcd(kfld) /= nureq(nnx,2)) THEN
             ier=0
             ijKmax = ibmax*jbmax*kmax
             CALL cnvout (ijKmax, ifldcd(kfld), nureq(nnx,2), fgdiv, fgq, ier)
             CALL pwater (fgq, fgplam, psmb)
             IF (ier > 0) WRITE (UNIT=nfprt, FMT='(6(A,I4))') &
                 ' Conversion Error at kflo = ', kflo, ' kfld = ', kfld, &
                 ' Error = ', ier, ' nfe = ', nnx, &
                 ' ifldcd = ', ifldcd(kfld), ' nureq = ', nureq(nnx,2)
          ELSE
             CALL pwater (fgdiv, fgplam, psmb)
          END IF
       END IF

       CALL sig2po (psmb, alnpmd, fgdiv, fgdivq)
    ELSE
       DO j=1,Jbmax
          DO i=1,Ibmaxperjb(j)
             fgdivq(i,1,j)=fgpphi(i,j)
          END DO
       END DO
    END IF


    IF ((mmm==2 .AND. iclcd(nnn)==4) .OR. &
         (mm1==2 .AND. iclcd(nn1)==4)) THEN
       DO j=1,Jbmax
          DO i=1,Ibmaxperjb(j)
             fgplam(i,j)=EXP(fgdivq(i,1,j))
          END DO
       END DO
    END IF

    IF (mmm /= 2) THEN
       DO j=1,Jbmax
          DO k=1,Lmax
             DO i=1,Ibmaxperjb(j)
                og(i,k,j)=fgdivq(i,k,j)
             END DO
          END DO
       END DO
    END IF

    IF ((mmm==2 .AND. iclcd(nnn)>=3) .OR. &
         (mm1==2 .AND. iclcd(nn1)>=3)) THEN
       DO j=1,Jbmax
          DO i=1,Ibmax
             qg(i,j)=fgplam(i,j)
          END DO
       END DO
    END IF
    IF (ALL(qg == 0.0_r8) .AND. ALL(og(:,1:lout,:) == 0.0_r8)) &
       WRITE (UNIT=nferr, FMT=*) " stog: og and qg NOT SET"

  END SUBROUTINE stog


  SUBROUTINE recon (nfld, nflp, nof, indate, title, &
       specal, rcode, dtin, nFile)

    IMPLICIT NONE

    ! It is used to reconcile the differences between the
    ! input file directory and the requested field directory and
    ! produce a final output directory.  Input datasets are:
    ! 1.  input file directory on dataset unit nfdir
    !     order of fields in directory must correspond to fields
    !     in the input file
    ! 2.  requested field directory on dataset unit nfrfd
    !     order of fields is independent
    !
    ! Derived field table is at Constants Module
    ! order of fields to be derived must correspond to the
    ! order they are processed in the code.  order of fields
    ! required for deriving a field is independent


    INTEGER, INTENT(IN) :: nFile
    INTEGER, INTENT(OUT) :: nfld
    INTEGER, INTENT(OUT) :: nflp
    INTEGER, INTENT(OUT) :: nof
    INTEGER, INTENT(OUT) :: indate(4)

    CHARACTER (LEN=4), INTENT(OUT) :: rcode !irrelevant
    CHARACTER (LEN=4), INTENT(OUT) :: dtin !irrelevant
    CHARACTER (LEN=40), INTENT(OUT) :: title
    CHARACTER (LEN=40), INTENT(OUT) :: specal

    INTEGER :: nwn !irrelevant
    INTEGER :: Mend1 ! irrelevant
    INTEGER :: levin ! # of vertical layers
    INTEGER :: levqin !irrelevant
    INTEGER :: in
    INTEGER :: kk,k
    INTEGER :: mm
    INTEGER :: ll
    INTEGER :: nn
    INTEGER :: mxdir
    INTEGER :: ii
    INTEGER :: jj
    INTEGER :: ifex
    INTEGER :: idate(4)

    CHARACTER (LEN=5) :: trunc !irrelevant
    CHARACTER (LEN=4) :: sdain !irrelevant
    CHARACTER (LEN=20) :: type
    CHARACTER (LEN=20) :: type1='COLA SIGMA HISTORY 4'
    CHARACTER (LEN=20) :: type2='CPTEC SIGMA VERS 2.0'

    ! read in input file directory

    READ (UNIT=nfdir, FMT='(A20)', END=4000) TYPE
    IF(myid.eq.0) WRITE (UNIT=nfprt, FMT='(1X,A)') TYPE
    IF(myid.eq.0) WRITE (UNIT=nfprt, FMT='(1X,A)') type1
    IF(myid.eq.0) WRITE (UNIT=nfprt, FMT='(1X,A)') type2

    !only uses levin and indate
!    WRITE (UNIT=n   , FMT='(A4,1X,A4,1X,A5,1X,11I5,1X,A4   )'                  ) nexp , 'SEQU', imdl , mMax , kmax , kmax  ,ihr,iday,mon,iyr, idate, 'TAPE'
    READ (UNIT=nfdir, FMT='(A4,1X,A4,1X,A5,1X,11I5,1X,A4   )', END=4000        ) rcode, sdain , trunc, Mend1, levin, levqin, idate          , indate, dtin
    IF(myid.eq.0) WRITE (UNIT=nfprt, FMT='(1X,A4,1X,A4,1X,A5,1X,11I5,1X,A4   )') rcode, sdain , trunc, Mend1, levin, levqin, idate          , indate, dtin
    READ (UNIT=nfdir, FMT='(2A40)', END=4000) title, specal
    IF(myid.eq.0) WRITE (UNIT=nfprt, FMT='(1X,2A40)') title, specal
    READ (UNIT=nfdir, FMT='(5E16.8)', END=4000) a_hybr_cb
    DO  k=1,levin
        a_hybr(k)=a_hybr_cb(levin-k+1)  !  Pa
    END DO
    READ (UNIT=nfdir, FMT='(5E16.8)', END=4000) a_hybr_cb
    DO  k=1,levin
        b_hybr(k)=a_hybr_cb(levin-k+1)
    END DO
    a_hybr_cb=0.0_r8
    IF(myid.eq.0) WRITE (UNIT=nfprt,FMT=*) ' writing a_hybr , levin = ',levin
    IF(myid.eq.0) WRITE (UNIT=nfprt, FMT='(1X,5E16.8)') (a_hybr(in),in=1,levin)
    IF(myid.eq.0) WRITE (UNIT=nfprt,FMT=*) ' writing b_hybr , levin = ',levin
    IF(myid.eq.0) WRITE (UNIT=nfprt, FMT='(1X,5E16.8)') (b_hybr(in),in=1,levin)
    !a_hybr_cb(k) = a_hybr(k) / 100._r4  ! transform from Pa to mbar
    a_hybr_cb = a_hybr / 1000._r8        ! transform from Pa to cbar
    kk=0
    nflp=0
    ! read in nfld field descriptions, skipping those where
    ! prodia == FIXD
    DO nfld=1,ndi
       DO
          READ (UNIT=nfdir, FMT='(A40,2X,A4,2X,I8,3X,I4,4X,I3)', END=30) &
                chrdsc(nfld), prodia(nfld), nharm(nfld), nlevs(nfld), ifldcd(nfld)
          IF (prodia(nfld) /= 'FIXD') EXIT
       END DO
    END DO
    IF(myid.eq.0) WRITE (UNIT=nfprt, FMT='(3(A,I5),A,I3)') ' Directory Table Overflow.  nfld = ', &
           nfld, ' nflp = ', nflp, ' ndv = ',ndv, ' kk = ', kk
    STOP 3100
30  CONTINUE
    nfld=nfld-1
    mkdir(1:nfld)=0

    ! read in requested field directory

    DO nflp=1,ndp
       READ (UNIT=nfrfd, FMT='(A40,I5,1X,A4)', END=40) &
                   chrdo(nflp), nuc(nflp), alias(nflp)
    END DO
    IF(myid.eq.0) WRITE (UNIT=nfprt, FMT='(3(A,I5),A,I3)') ' Directory Table Overflow.  nfld = ', &
           nfld, ' nflp = ', nflp, ' ndv = ',ndv, ' kk = ', kk
    STOP 3100
40  CONTINUE
    nflp=nflp-1
    mko(1:nflp)=0
    ife(1:nflp)=0
    mkdv(1:ndv)=0
    lif(1:ndv)=0

    ! determine direct availability of requested fields

    DO mm=1,nfld
       ! see if field is to be defered
       DO ll=1,ndv
          IF (chrdsc(mm) == chrdv(ll)) THEN
             ! mark input field as to be defered (derived from itself?)
             ! mark derived field as pending derivation from input field
             mkdir(mm)=100
             mkdv(ll)=2
             EXIT
          END IF
       END DO
       IF (ll <= ndv) CYCLE 

       DO nn=1,nflp
          IF (chrdsc(mm) == chrdo(nn)) THEN
             ! mark input field as needed directly
             ! mark requested field as being available directly from inputfield
             ! save input field number
             mkdir(mm)=1
             mko(nn)=1
             ife(nn)=mm
             EXIT
          END IF
       END DO
    END DO

    ! examine derived field table
    DO nn=1,nflp
       ! ignore directly available fields
       IF (mko(nn) == 1) CYCLE
       DO ll=1,ndv
          ! ignore undesired derivable fields
          IF (chrdv(ll) /= chrdo(nn)) CYCLE
          kk=kdv(ll)
          DO jj=1,kk
             DO mm=1,nfld
                ! verify availability of required fields for derivation
                ! matches required? if not, skip this input
                IF (chrdsc(mm) /= chreq(ll,jj)) CYCLE
                ! matches target? if so, dont save it's number
                IF (chrdsc(mm) /= chrdv(ll)) EXIT
                ! save input file field number
                IF (mkdv(ll) == 2) ife(nn)=mm
                EXIT
             END DO
             ! all inputs were inspected and no match occurred
             ! some required field is missing from the input
             IF (mm > nfld) EXIT
          END DO
          ! all inputs were inspected and no match occurred
          ! some required field is missing from the input
          IF (jj <= kk) EXIT
          ! ife(nn) contains last input required for this derivation
          IF (mkdv(ll) == 2) THEN
             ! mark requested field as defered if listed in input file directory
             ! save derivation rule number on "high" part of ife
             mko(nn)=3
             ife(nn)=ife(nn)+1000*ll
          ELSE 
             ! mark requested field as available by derivation
             ! save derived field number corresponding to requested field
             mko(nn)=2
             ife(nn)=ll
          ENDIF
          ! mark derived field as computable from input fields
          mkdv(ll)=1
          EXIT
       END DO
    END DO

    ! review input file directory
    DO ll=1,ndv
       ! ignore underivable or undesired fields in derived field table
       IF (mkdv(ll) /= 1) CYCLE
       kk=kdv(ll)
       DO  mm=1,nfld
          ! this is the same test as the first
          ! successfully derivable fields have mkdir == 200
          IF (chrdv(ll) == chrdsc(mm)) mkdir(mm)=mkdir(mm)+100
          DO jj=1,kk
             ! match current input field with a required field
             IF (chreq(ll,jj) /= chrdsc(mm)) CYCLE
             ! mark input field as needed for computing derived fields
             ! if it is not the derived field itself
             IF (chrdv(ll) /= chrdsc(mm)) mkdir(mm)=mkdir(mm)+2
             mxdir=mm
             EXIT
          END DO
       END DO
       ! save last required field input file sequence number
       lif(ll)=mxdir
    END DO

    ! create output file directory
    nof=0
    ll=1
    DO mm=1,nfld
       IF (nFile <= 0) THEN
          IF (prodia(mm) == 'DIAG') mkdir(mm)=-1
       ENDIF
       IF (mkdir(mm) < 1) CYCLE
       ! if field is directly available
       IF (MOD(mkdir(mm),2) == 1) THEN
          nof=nof+1
          ! add directly requested field to output directory
          chrop(nof)=chrdsc(mm)
          nvo(nof)=Lmax
          IF (nlevs(mm) == 1) nvo(nof)=1
       ENDIF
       ! any derivation rules left?
       IF (ll > ndv) CYCLE
       DO
          IF (mkdv(ll) == 1) THEN
             ! mkdir >= 1 means is required for some computation
             IF (mkdir(mm)>1 .AND. lif(ll)<=mm) THEN
                ! add requested derived field to output directory after last
                ! required input field has been processed
                nof=nof+1
                chrop(nof)=chrdv(ll)
                IF (nvv(ll) == 1) nvo(nof)=1
                IF (nvv(ll) == 2) nvo(nof)=Lmax
             ELSE
                EXIT
             ENDIF
          END IF
          ll=ll+1
          IF (ll > ndv) EXIT
       END DO
    END DO
    DO ii=1,nof
       nfe(ii)=0
       jfe(ii)=0
       DO nn=1,nflp
          IF (chrop(ii) == chrdo(nn)) THEN
             ifex=MOD(ife(nn),1000)
             IF (nuc(nn) == -1) THEN
                IF (mko(nn) == 1) nuco(ii)=ifldcd(ifex)
                IF (mko(nn) == 2) nuco(ii)=nudv(ifex)
                IF (mko(nn) == 3) nuco(ii)=ifldcd(ifex)
                ! if (mko(nn) == 3) nuco(ii)=nudv(ifex)
             ELSE
                nuco(ii)=nuc(nn)
             END IF
             IF ( mko(nn) == 3) THEN
                nfe(ii)=Min(ife(nn)/1000,ndv )
             ELSE
                nfe(ii)=Min(ifex,ndv )
             END IF
             jfe(ii)=nn
             ife(nn)=ifex
             EXIT
          ENDIF
       END DO
    END DO
    ! print input file directory anotated by field requirements
    IF(myid.eq.0) THEN
     WRITE (UNIT=nfprt, FMT='(/,A,/,A,3I3,I5,A,3I3,I5,/,&
                         &A,T22,A,T51,A,T60,A,T67,A,T73,A,T78,A)') &
          ' I n p u t   F i l e   D i r e c t o r y', &
          '   idate = ', idate, ' indate = ', indate, &
          ' number', 'chrdsc', 'prodia', 'len', 'lev', 'code', 'mkdir'
    WRITE (UNIT=nfprt, FMT='(1X,I3,T9,A40,3X,A4,2X,I8,I4,2X,I5,3X,I3)') &
          (mm, chrdsc(mm), prodia(mm), nharm(mm), nlevs(mm), &
           ifldcd(mm), mkdir(mm), mm=1,nfld)
    ! print requested file directory anotated by field availability
    WRITE (UNIT=nfprt, FMT='(/,A,/,A,T24,A,T56,A,T64,A,T71,A)') &
          ' R e q u e s t e d   F i e l d   D i r e c t o r y', &
          ' number', 'chrdo', 'nuc', 'mko', 'ife'
    WRITE (UNIT=nfprt, FMT='(1X,T3,I3,T11,A40,3X,I5,4X,I3,4X,I5,1X,A4)') &
          (nn, chrdo(nn), nuc(nn), mko(nn), ife(nn), alias(nn),nn=1,nflp)
    ! print derived field table anotated by field availability
    WRITE (UNIT=nfprt, FMT='(/,A,/)') ' D e r i v e d   F i e l d   T a b l e'
    DO ll=1,ndv
       WRITE (UNIT=nfprt, FMT='(A,I3,A,A40,/,6(A,I3))') &
             ' number', ll, ' chrdv:  ', chrdv(ll), ' kdv', kdv(ll), &
             ' nvv', nvv(ll), ' mkdv', mkdv(ll), ' lif', lif(ll), &
             ' iclcd', iclcd(ll), ' nudv', nudv(ll)
       WRITE (UNIT=nfprt, FMT='(1X,T18,A40,2X,I3)') &
             (chreq(ll,kk), nureq(ll,kk), kk=1,kdv(ll))
    END DO
    ENDIF
    DO ii=1,nof
       aliop(ii)=alias(jfe(ii))
    ENDDO
    ! print output file directory
    IF(myid.eq.0) THEN
    WRITE (UNIT=nfprt, FMT='(/,A,/,A,T22,A,T54,A,T60,A,T67,A,T73,A)') &
          ' O u t p u t   F i l e   D i r e c t o r y', &
          ' number', 'chrop', 'nvo', 'nuco', 'nfe', 'jfe'
    WRITE (UNIT=nfprt, FMT='(1X,I3,5X,A40,3X,I3,2X,I5,4X,I3,3X,I3,1X,A4)') &
          (ii, chrop(ii), nvo(ii), nuco(ii), nfe(ii), &
           jfe(ii), aliop(ii), ii=1,nof)
    ENDIF

    REWIND (UNIT=nfrfd)

    RETURN
    4000 WRITE (UNIT=nferr, FMT='(A)') &
               ' Unexpected End of File in Input File Directory.'
    STOP 4100

  END SUBROUTINE recon


  SUBROUTINE rwrite (lev, kflo, field, bfr, collect)

    IMPLICIT NONE

    !     write gaussian grid to output file

    INTEGER, INTENT(IN) :: lev
    INTEGER, INTENT(IN) :: kflo
    LOGICAL, INTENT(IN) :: collect

    REAL (KIND=r8), INTENT(IN) :: field(Ibmax,lev,Jbmax)
    REAL (KIND=r8), INTENT(INOUT) :: bfr(Imax,Jmax,lev)

    INTEGER :: l,ier,j,i

    LOGICAL :: t1, t2, t3, t4

    REAL (KIND=r8) :: bbfr(Imax,Jmax), breg(Idim,Jdim)

    INTEGER :: kpds(200),kgds(200)
    
    IF (collect) CALL Collect_Grid_Full(field, bfr, lev, 0)

    IF (myid.ne.0) return

    WRITE (UNIT=nfprt, FMT='(2(A,I5))') ' Lev = ', lev, ' kflo = ', kflo

    t1 = mko(jfe(kflo)) == 2
    t2 = mko(jfe(kflo)) == 0
    t3 = nuco(kflo)     == nudv(nfe(kflo)) 
    t4 = nuco(kflo)     == ifldcd(ife(jfe(kflo))) 

    !    
    CALL GDSPDSSETION (kgds, kpds )

    IF (.NOT. (t1 .OR. t2 .OR. t4)) THEN
       DO l=1,lev
          ier=0
          CALL cnvout (ngaus, ifldcd(ife(jfe(kflo))), nuco(kflo), &
               bfr(1,1,l), bbfr, ier)
          IF (ier > 0) WRITE (UNIT=nfprt, FMT='(6(A,I4))') &
             ' Conversion Error at kflo = ', kflo, ' Error = ', ier, &
             ' jfe = ', jfe(kflo), ' ife = ', ife(jfe(kflo)), &
             ' ifldcd = ', ifldcd(ife(jfe(kflo))), ' nuco = ',nuco(kflo)
          IF (RegInt) THEN
             IF (Binary) THEN
                !
                ! binary format
                !
                CALL DoAreaInterpolation (bbfr, breg)
                CALL WriteField (mgaus, breg)
             ELSE
                !
                ! grib format
                !
                CALL DoAreaInterpolation (bbfr, breg)
                CALL WriteGrbField (aliop(kflo),mgaus,kgds,kpds,breg,l)
             END IF
          ELSE
             IF (Binary) THEN
                !
                ! binary format
                !
                CALL DoAreaGausInterpolation(bbfr, breg)
                CALL WriteField (mgaus, breg)
             ELSE             
                !
                ! grib format
                !
                CALL DoAreaGausInterpolation(bbfr, breg)
                CALL WriteGrbField (aliop(kflo),mgaus,kgds,kpds,breg,l)
             ENDIF
          END IF
       END DO
    ELSE IF (t1 .AND. (.NOT. t3)) THEN
       DO l=1,lev
          ier=0
          CALL cnvout (ngaus, nudv(nfe(kflo)), nuco(kflo), &
               bfr(1,1,l), bbfr, ier)
          IF (ier > 0) WRITE (UNIT=nfprt, FMT='(5(A,I4))') &
             ' Conversion Error at kflo = ', kflo, ' Error = ', ier, &
             ' nfe = ', nfe(kflo), ' nudv = ', nudv(nfe(kflo)), ' nuco = ',nuco(kflo)
          IF (RegInt) THEN
            IF (Binary) THEN
                !
                ! binary format
                !
                CALL DoAreaInterpolation (bbfr, breg)
                CALL WriteField (mgaus, breg)
            ELSE
                !
                ! grib format
                !
                CALL DoAreaInterpolation (bbfr, breg)
                CALL WriteGrbField (aliop(kflo),mgaus,kgds,kpds, breg,l)
            END IF 
          ELSE
             IF (Binary) THEN
                !
                ! binary format
                !
                CALL DoAreaGausInterpolation(bbfr, breg)
                CALL WriteField (mgaus, breg)
                !
             ELSE
                !
                ! grib format
                !
                CALL DoAreaGausInterpolation(bbfr, breg)
                CALL WriteGrbField (aliop(kflo),mgaus,kgds,kpds, breg,l)
                !
             ENDIF
          END IF
       END DO
    ELSE IF ((t1 .AND. t3) .OR. t2 .OR. &
         ((.NOT. t1) .AND. (.NOT. t2) .AND. t4)) THEN
       DO l=1,lev
          DO j=1,Jmax
             DO i=1,Imax
                bbfr(i,j)=bfr(i,j,l)
             END DO
          END DO
          IF (RegInt) THEN
             IF (Binary) THEN
                !
                ! binary format
                !
                CALL DoAreaInterpolation (bbfr, breg)
                CALL WriteField (mgaus, breg)
             ELSE
                !
                ! grib format
                !
                CALL DoAreaInterpolation (bbfr, breg)
                CALL WriteGrbField (aliop(kflo),mgaus,kgds,kpds,breg,l)
             END IF
          ELSE
             IF (Binary) THEN
                !
                ! binary format
                !
                CALL DoAreaGausInterpolation(bbfr, breg)
                CALL WriteField (mgaus, breg)
                !
             ELSE
                !
                ! grib format
                !
                CALL DoAreaGausInterpolation(bbfr, breg)
                CALL WriteGrbField (aliop(kflo),mgaus,kgds,kpds,breg,l)
                !
             ENDIF
          END IF
       END DO
    END IF

  END SUBROUTINE rwrite
  
  SUBROUTINE GeraGribCtl(unt,fname,title,labelp,nof,ndp,nvo,aliop)
    USE Constants, ONLY : res,RunRecort,RecLat,RecLon,newlat0,newlat1,newlon0,newlon1
    USE RegInterp, ONLY : gLats,glond, glatd

   IMPLICIT NONE
   INTEGER    ,INTENT(IN   ) :: unt
   CHARACTER(LEN=256),INTENT(INOUT) :: fname
   CHARACTER(LEN=40) ,INTENT(IN   ) :: title
   CHARACTER(LEN=10) ,INTENT(IN   ) :: labelp
   INTEGER  ,INTENT(IN   ) :: nof
   INTEGER  ,INTENT(IN   ) :: ndp
   INTEGER  ,INTENT(IN   ) :: nvo  (ndp)
   CHARACTER(LEN=4),INTENT(IN   ) :: aliop(ndp)

   !
   INTEGER :: ifna
   INTEGER :: ifnb
   INTEGER :: it,iy, im, id, ih,ii,j,k,inv,itypelev,iiplev,i
   CHARACTER (LEN=20), PARAMETER :: type='PRESSURE HISTORY    '
   CHARACTER (LEN=3), PARAMETER :: cmth(12)=(/&
         'JAN','FEB','MAR','APR','MAY','JUN', &
         'JUL','AUG','SEP','OCT','NOV','DEC'/)

   !
   !  write output directory
   !
   INQUIRE (UNIT=unt, NAME=fname)
   ifnb=INDEX(fname//' ',' ')-5
   ifna=ifnb+1
   DO
      ifna=ifna-1
      IF (fname(ifna:ifna) == '/') EXIT
   END DO
   ifna=ifna+1
   PRINT*,' OUTPUT FILE: '//fname,ifna,ifna
   WRITE (UNIT=*, FMT='(A)') ' OUTPUT FILE: '//fname(ifna:ifnb)//'.grb'
   WRITE (UNIT=unt, FMT='(A)') 'dset ^'//fname(ifna:ifnb)//'.grb'
   WRITE (UNIT=unt, FMT='(A)') '*'
   WRITE (UNIT=unt, FMT='(A)') 'index ^'//fname(ifna:ifnb)//'.idx'    
   WRITE (UNIT=unt, FMT='(A)') '*'
   IF (RegInt) THEN
      WRITE (UNIT=unt, FMT='(A)') 'undef -2.56E+33'
   ELSE
      WRITE (UNIT=unt, FMT='(A)') 'undef 9.999E+20'
   END IF
   WRITE (UNIT=unt, FMT='(A)') '*'
   WRITE (UNIT=unt, FMT='(3A)') 'title ',type,title
   WRITE (UNIT=unt, FMT='(A)') '*'
   WRITE (UNIT=unt, FMT='(A,I6)') 'dtype grib',table1(size_tb(1))%id 
   WRITE (UNIT=unt, FMT='(A)') '*'
   WRITE (UNIT=unt, FMT='(A)') 'options yrev'
   WRITE (UNIT=unt, FMT='(A)') '*'
   IF(RunRecort) THEN
       IF (RegInt) THEN
          WRITE (UNIT=unt, FMT='(A,I5,A,F8.3,F15.10)') &
               'xdef ', newlon1-newlon0+1, ' linear ', glond(newlon0), 360.0_r8/REAL(Idim,r8)
          WRITE (UNIT=unt, FMT='(A,I5,A,F8.3,F15.10)') &
               'ydef ', newlat1-newlat0+1, ' linear ', glatd(newlat1), 180.0_r8/REAL(Jdim,r8)
       ELSE
          WRITE (UNIT=unt, FMT='(A,I5,A,F8.3,F15.10)') &
               'xdef ', newlon1-newlon0+1, ' linear ', glond(newlon0), 360.0_r8/REAL(Idim,r8)
          WRITE (UNIT=unt, FMT='(A,I5,A)') 'ydef ', newlat1-newlat0+1, ' levels '
          IF(res<=0)THEN
             WRITE (UNIT=12, FMT='(8F10.5)') (glatd(j),j=newlat1,newlat0,-1)
          ELSE
             WRITE (UNIT=12, FMT='(8F10.5)') (glatd(j),j=newlat1,newlat0,-1)
          END IF
       END IF
    ELSE
      IF (RegInt) THEN
         WRITE (UNIT=unt, FMT='(A,I5,A,F8.3,F15.10)') &
               'xdef ', Idim, ' linear ', 0.0_r8, 360.0_r8/REAL(Idim,r8)
         WRITE (UNIT=unt, FMT='(A,I5,A,F8.3,F15.10)') &
              'ydef ', Jdim, ' linear ',  -90.0_r8, 180.0_r8/REAL(Jdim-1,r8)
      ELSE
         WRITE (UNIT=unt, FMT='(A,I5,A,F8.3,F15.10)') &
               'xdef ', Idim, ' linear ', 0.0_r8, 360.0_r8/REAL(Idim,r8)
         WRITE (UNIT=unt, FMT='(A,I5,A)') 'ydef ', Jdim, ' levels '
         IF(res<=0)THEN
           WRITE (UNIT=12, FMT='(8F10.5)') (lati(j),j=jMax,1,-1)
         ELSE
           WRITE (UNIT=12, FMT='(8F10.5)') (gLats(j),j=Jdim,1,-1)
         END IF
      END IF
    END IF

   it=1
   READ (labelp, FMT='(I4,3I2)') iy, im, id, ih
   WRITE (UNIT=unt, FMT='(A,I5,A,I2.2,A,I2.2,A,I4,A)') &
        'tdef ', it, ' linear ', ih, 'Z', id, cmth(im), iy, ' 6hr'
   WRITE (UNIT=unt, FMT='(A)') '*'
   IF (Lmax <= 10) THEN
      WRITE (UNIT=unt, FMT='(A,I5,A,10I5)') 'zdef ', Lmax,'  levels ', &
           (NINT(pmand(k)),k=1,Lmax)
   ELSE
      WRITE (UNIT=unt, FMT='(A,I5,A,10I5)') 'zdef ', Lmax, ' levels ', &
           (NINT(pmand(k)),k=1,10)
      WRITE (UNIT=unt, FMT='((16X,10I5))') (NINT(pmand(k)),k=11,Lmax)
   END IF
   WRITE (UNIT=unt, FMT='(A,I5)') 'vars ', nof+2
   WRITE (UNIT=unt, FMT='(A)') 'topo  0 132,1,0 '// &
        '** surface TOPOGRAPHY [m]'
   WRITE (UNIT=unt, FMT='(A)') 'lsmk  0  81,1,0 '// &
        '** surface LAND SEA MASK [0,1]'
   DO ii=1,nof
       inv=nvo(ii)
       IF (inv == 1) inv=0
       DO i=1,size_tb(1)
           IF(aliop(ii) == table1(i)%name) THEN
                itypelev=100
                iiplev  = 0
               DO k=1,size_tb(2)
                  IF(TRIM(table1(i)%level) == TRIM(table2(k)%level_type))THEN 
                     itypelev=table2(k)%default 
                     iiplev  =table2(k)%p2
                  END IF 
               END DO   
              WRITE (UNIT=unt, FMT='(A,I5,I5,A1,I5,A1,I5,1X,A)')table1(i)%name, inv,table1(i)%id,',',&
                     itypelev,',',iiplev,' ** '//table1(i)%level//table1(i)%title//'('//table1(i)%unit//')'
           END IF
       END DO
    ENDDO
    WRITE (UNIT=unt, FMT='(A)') 'endvars'  
  END SUBROUTINE GeraGribCtl
  
  
  SUBROUTINE GeraBinCtl(unt,fname,title,labelp,nof,ndp,nvo,aliop,chrop,nuco)
    USE Constants, ONLY : Undef,res,RunRecort,RecLat,RecLon,newlat0,newlat1,newlon0,newlon1
    USE RegInterp, ONLY : gLats,glond, glatd

   IMPLICIT NONE
   INTEGER    ,INTENT(IN   ) :: unt
   CHARACTER(LEN=256),INTENT(INOUT) :: fname
   CHARACTER(LEN=40) ,INTENT(IN   ) :: title
   CHARACTER(LEN=10) ,INTENT(IN   ) :: labelp
   INTEGER  ,INTENT(IN   ) :: nof
   INTEGER  ,INTENT(IN   ) :: ndp
   INTEGER  ,INTENT(IN   ) :: nvo  (ndp)
   CHARACTER(LEN=4),INTENT(IN   ) :: aliop(ndp)
   CHARACTER (LEN=40),INTENT(IN   ):: chrop(ndp)
   INTEGER,INTENT(IN   ) :: nuco(ndp)

   !
   INTEGER :: ifna
   INTEGER :: ifnb
   INTEGER :: it,iy, im, id, ih,ii,j,k,inv
   CHARACTER (LEN=20), PARAMETER :: type='PRESSURE HISTORY    '
   CHARACTER (LEN=3), PARAMETER :: cmth(12)=(/&
         'JAN','FEB','MAR','APR','MAY','JUN', &
         'JUL','AUG','SEP','OCT','NOV','DEC'/)

   !
   !  write output directory
   !
       INQUIRE (UNIT=11, NAME=fname)
       ifnb=INDEX(fname//' ',' ')-1
       ifna=ifnb+1
       DO
          ifna=ifna-1
          IF (fname(ifna:ifna) == '/') EXIT
       END DO
       ifna=ifna+1
       WRITE (UNIT=*, FMT='(/,A,/)') ' OUTPUT FILE: '//fname(ifna:ifnb)
       WRITE (UNIT=12, FMT='(A)') 'DSET ^'//fname(ifna:ifnb)
       WRITE (UNIT=12, FMT='(A)') '*'
       WRITE (UNIT=12, FMT='(A)') 'OPTIONS SEQUENTIAL YREV BIG_ENDIAN'
       WRITE (UNIT=12, FMT='(A)') '*'
       WRITE (UNIT=12, FMT='(A,1PE9.2)') 'UNDEF ', Undef
       WRITE (UNIT=12, FMT='(A)') '*'
       WRITE (UNIT=12, FMT='(3A)') 'TITLE ',type,title
       WRITE (UNIT=12, FMT='(A)') '*'

       IF(RunRecort) THEN
          IF (RegInt) THEN
             WRITE (UNIT=unt, FMT='(A,I5,A,F8.3,F15.10)') &
                  'xdef ', newlon1-newlon0+1, ' linear ', glond(newlon0), 360.0_r8/REAL(Idim,r8)
             WRITE (UNIT=unt, FMT='(A,I5,A,F8.3,F15.10)') &
                  'ydef ', newlat1-newlat0+1, ' linear ', glatd(newlat1), 180.0_r8/REAL(Jdim,r8)
          ELSE
             WRITE (UNIT=unt, FMT='(A,I5,A,F8.3,F15.10)') &
                  'xdef ', newlon1-newlon0+1, ' linear ', glond(newlon0), 360.0_r8/REAL(Idim,r8)
             WRITE (UNIT=unt, FMT='(A,I5,A)') 'ydef ', newlat1-newlat0+1, ' levels '
             IF(res<=0)THEN
                WRITE (UNIT=12, FMT='(8F10.5)') (glatd(j),j=newlat1,newlat0,-1)
             ELSE
                WRITE (UNIT=12, FMT='(8F10.5)') (glatd(j),j=newlat1,newlat0,-1)
             END IF
          END IF

       ELSE

          IF (RegInt) THEN
             WRITE (UNIT=12, FMT='(A,I5,A,F8.3,F15.10)') &
               'XDEF ', Idim, ' LINEAR ', 0.0_r8, 360.0_r8/REAL(Idim,r8)
             WRITE (UNIT=12, FMT='(A,I5,A,F8.3,F15.10)') &
                  'YDEF ', Jdim, ' LINEAR ', -90.0_r8, 180.0_r8/REAL(Jdim-1,r8)
          ELSE
             WRITE (UNIT=12, FMT='(A,I5,A,F8.3,F15.10)') &
               'XDEF ', Idim, ' LINEAR ', 0.0_r8, 360.0_r8/REAL(Idim,r8)
             WRITE (UNIT=12, FMT='(A,I5,A)') 'YDEF ', Jdim, ' LEVELS '
             IF(res<=0)THEN
               WRITE (UNIT=12, FMT='(8F10.5)') (lati(j),j=jMax,1,-1)
             ELSE
               WRITE (UNIT=12, FMT='(8F10.5)') (gLats(j),j=Jdim,1,-1)
             END IF
          END IF
       END IF

       IF (Lmax <= 10) THEN
          WRITE (UNIT=12, FMT='(A,I5,A,10I5)') 'ZDEF ', Lmax,'  LEVELS ', &
               (NINT(pmand(k)),k=1,Lmax)
       ELSE
          WRITE (UNIT=12, FMT='(A,I5,A,10I5)') 'ZDEF ', Lmax, ' LEVELS ', &
               (NINT(pmand(k)),k=1,10)
          WRITE (UNIT=12, FMT='((16X,10I5))') (NINT(pmand(k)),k=11,Lmax)
       END IF
       it=1
       READ (labelp, FMT='(I4,3I2)') iy, im, id, ih
       WRITE (UNIT=12, FMT='(A,I5,A,I2.2,A,I2.2,A,I4,A)') &
            'TDEF ', it, ' LINEAR ', ih, 'Z', id, cmth(im), iy, ' 6HR'
       WRITE (UNIT=12, FMT='(A)') '*'
       WRITE (UNIT=12, FMT='(A,I5)') 'VARS ', nof+2
       WRITE (UNIT=12, FMT='(A)') 'TOPO    0 99 '// &
            'TOPOGRAPHY                              (M               )'
       WRITE (UNIT=12, FMT='(A)') 'LSMK    0 99 '// &
            'LAND SEA MASK                           (NO DIM          )'
       DO ii=1,nof
          inv=nvo(ii)
          IF (inv == 1) inv=0
          WRITE (UNIT=12, FMT='(A,I5,I3,1X,A)') aliop(ii), inv, 99, &
               chrop(ii)//'('//GiveUnit(nuco(ii))//')'
       ENDDO
       WRITE (UNIT=12, FMT='(A)') 'ENDVARS'
  END SUBROUTINE GeraBinCtl
END MODULE PostLoop
