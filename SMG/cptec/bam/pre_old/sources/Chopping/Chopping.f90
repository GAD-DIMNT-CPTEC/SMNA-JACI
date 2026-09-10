!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
PROGRAM Chopping

  USE InputParameters, ONLY: InitParameters, r8, nferr, nfinp, nfprt, nftop, &
                             nfsig, nfnmc, nfcpt, nfozw, nftrw, nficr, nfozr, &
                             nftrr, nficw, nfozg, nftrg, nfgrd, nfctl, &
                             MendInp, KmaxInp, MendOut, KmaxOut, MendMin, &
                             MendCut, Iter, SmthPerCut, GetOzone, GetTracers, &
                             GrADS, GrADSOnly, GDASOnly, SmoothTopo, &
                             ImaxOut, JmaxOut, TruncInp, &
                             Mnwv2Inp, Mnwv2Out, Mnwv3Inp, Mnwv3Out, &
                             NTracers, Kdim, ICaseDec, &
                             ForecastDay, TimeOfDay, cTv, MonChar, &
                             DataCPT, DataInp, DataOut, DataTop, DataSig, DataSigInp, &
                             GDASInp, OzonInp, TracInp, OzonOut, TracOut, &
                             DirMain, DirInp, DirOut, DirTop, &
                             DirSig, DirGrd, DGDInp

  USE InputArrays, ONLY: GetArrays, ClsArrays, GetSpHuTracers, ClsSpHuTracers, &
                         DateInitial, DateCurrent, DelSInp, SigIInp, SigLInp, SigIOut, &
                         SigLOut, qWorkInp, qWorkOut, DelSigmaInp, SigInterInp, &
                         SigLayerInp, DelSigmaOut, SigInterOut, SigLayerOut, &
                         qTopoInp, qLnPsInp, qTopoOut, qLnPsOut, qTopoRec, &
                         qDivgInp, qVortInp, qTvirInp, qDivgOut, qVortOut, &
                         qTvirOut, qUvelInp, qVvelInp, qSpHuInp, qSpHuOut, &
                         gWorkOut, gTopoInp, gTopoOut, gTopoDel, gLnPsInp, &
                         gPsfcInp, gLnPsOut, gPsfcOut, gUvelInp, gVvelInp, &
                         gTvirInp,  gDivgInp, gVortInp, gPresInp, gUvelOut, &
                         gVvelOut, gTvirOut, gPresOut, gSpHuInp, gSpHuOut

  USE Fourier, ONLY: InitFFT, ClsMemFFT

  USE Legendre, ONLY: InitLegendre, ClsMemLegendre

  USE VerticalInterpolation, ONLY: VertSigmaInter

  USE Recomposition, ONLY: InitRecomposition, ClsMemRecomposition, &
                           RecompositionScalar, RecompositionVector, &
                           DivgVortToUV, glat, coslat

  USE Utils, ONLY: NewSigma, SigmaInp, NewPs

  USE Decomposition, ONLY: DectoSpherHarm, UVtoDivgVort, TransSpherHarm

  IMPLICIT NONE

  INTEGER :: ios, nRec, IOL

  INTEGER :: i, j, k, nt

  LOGICAL :: GetNewTop=.FALSE., GetNewSig=.FALSE., ExistGDAS=.FALSE., &
             ExistGANLCPT=.FALSE., ExistGANLSMT=.FALSE., ExistGANL=.FALSE., &
             VerticalInterp=.TRUE.

  CHARACTER (LEN=12) :: Tdef='  z         '

  REAL (KIND=r8), DIMENSION (:,:,:), ALLOCATABLE :: gOzonOut

  REAL (KIND=r8), DIMENSION (:,:,:,:), ALLOCATABLE :: gTracOut

  CALL InitParameters

  CALL GetArrays

  INQUIRE (FILE=TRIM(DGDInp)//TRIM(GDASInp), EXIST=ExistGDAS)
  INQUIRE (FILE=TRIM(DirInp)//TRIM(DataCPT), EXIST=ExistGANLCPT)
  INQUIRE (FILE=TRIM(DirInp)//TRIM(DataInp), EXIST=ExistGANLSMT)

  WRITE (UNIT=nfprt, FMT='(A,L6)') TRIM(DGDInp)//TRIM(GDASInp), ExistGDAS
  WRITE (UNIT=nfprt, FMT='(A,L6)') TRIM(DirInp)//TRIM(DataCPT), ExistGANLCPT
  WRITE (UNIT=nfprt, FMT='(A,L6)') TRIM(DirInp)//TRIM(DataInp), ExistGANLSMT

  IF (.NOT.ExistGDAS .AND. .NOT.ExistGANLCPT .AND. .NOT.ExistGANLSMT) THEN
     WRITE (UNIT=nfprt, FMT='(/,3(A,/))') &
           ' The NCEP Input File Does Not Exist and', &
           ' The CPTEC Input File Does Not Exist Also and', &
           ' The CPTEC Topo-Smoothed Input File Does Not Exist Also'
     STOP 'No NCEP or GANL File'
  END IF

  ExistGANL = ExistGANLCPT .OR. ExistGANLSMT
  IF (ExistGANLSMT) THEN
     IF (SmoothTopo) SmoothTopo=.FALSE.
  ELSE
     IF (ExistGANLCPT) THEN
        DataInp=DataCPT
        IF (.NOT.SmoothTopo .AND. MendOut > MendMin) SmoothTopo=.TRUE.
     END IF
  END IF

  IF (ExistGDAS .AND. .NOT.ExistGANL) THEN

    CALL GDAStoGANL
    CALL GetSpHuTracers

  ELSE

    WRITE (UNIT=nfprt, FMT='(/,A)') ' GANL File Already Exists'

    INQUIRE (FILE=TRIM(DirInp)//TRIM(OzonInp), EXIST=GetOzone)
    INQUIRE (FILE=TRIM(DirInp)//TRIM(TracInp), EXIST=GetTracers)
    IF (GetOzone) THEN
      NTracers=NTracers+1
      IF (GetTracers) THEN
        WRITE (UNIT=nfprt, FMT='(/,A)') &
              ' Considering Just One Other Tracer Than Ozone'
        NTracers=NTracers+1
      END IF
    END IF
    CALL GetSpHuTracers

    IF (GetOzone) THEN
      WRITE (UNIT=nfprt, FMT='(/,A)') ' Get Ozone from '//TRIM(OzonInp)
      OPEN (UNIT=nfozr, FILE=TRIM(DirInp)//TRIM(OzonInp), FORM='UNFORMATTED', &
            ACCESS='SEQUENTIAL', ACTION='READ', STATUS='OLD', IOSTAT=ios)
      IF (ios /= 0) THEN
        WRITE (UNIT=nferr, FMT='(3A,I4)') ' ** (Error) ** Open file ', &
                                          TRIM(TRIM(DirInp)//TRIM(OzonInp)), &
                                          ' returned IOStat = ', ios
        STOP ' ** (Error) **'
      END IF
      nt=2
      DO k=1,KmaxInp
        READ (UNIT=nfozr) qWorkInp
        qSpHuInp(:,k,nt)=qWorkInp
      END DO
      CLOSE(UNIT=nfozr)

      IF (GetTracers) THEN
        WRITE (UNIT=nfprt, FMT='(/,A)') ' Get Tracers from '//TRIM(TracInp)
        OPEN (UNIT=nftrr, FILE=TRIM(DirInp)//TRIM(TracInp), FORM='UNFORMATTED', &
              ACCESS='SEQUENTIAL', ACTION='READ', STATUS='OLD', IOSTAT=ios)
        IF (ios /= 0) THEN
          WRITE (UNIT=nferr, FMT='(3A,I4)') ' ** (Error) ** Open file ', &
                                            TRIM(TRIM(DirInp)//TRIM(TracInp)), &
                                            ' returned IOStat = ', ios
          STOP ' ** (Error) **'
        END IF
        nt=3
        DO k=1,KmaxInp
          READ (UNIT=nftrr) qWorkInp
          qSpHuInp(:,k,nt)=qWorkInp
        END DO
        CLOSE(UNIT=nftrr)
      ELSE
        WRITE (UNIT=nfprt, FMT='(/,A)') ' Other Tracers File Does Not Exist'
      END IF

    ELSE

      WRITE (UNIT=nfprt, FMT='(/,A)') &
            ' Ozone file Does Not Exist. Ignore Other Tracers'
      GetTracers=.FALSE.
      NTracers=1

    END IF

  END IF
  WRITE (UNIT=nfprt, FMT='(/,A,I5)') ' NTracers = ', NTracers

  CALL ICRead

  qTopoOut=0.0_r8
  gTopoDel=0.0_r8
  INQUIRE (FILE=TRIM(DirTop)//TRIM(DataTop), EXIST=GetNewTop)
  IF (GetNewTop) THEN
    WRITE (UNIT=nfprt, FMT='(/,A)') ' Getting New Topography'
    OPEN (UNIT=nftop, FILE=TRIM(DirTop)//TRIM(DataTop), FORM='UNFORMATTED', &
          ACCESS='SEQUENTIAL', ACTION='READ', STATUS='OLD', IOSTAT=ios)
    IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') ' ** (Error) ** Open file ', &
                                        TRIM(TRIM(DirTop)//TRIM(DataTop)), &
                                        ' returned IOStat = ', ios
      STOP ' ** (Error) **'
    END IF
    READ  (UNIT=nftop) qWorkOut
    CLOSE (UNIT=nftop)
    qTopoOut=qWorkOut
    WRITE (UNIT=nfprt, FMT='(/,A)') ' TopoOut for New Topography:'
    WRITE (UNIT=nfprt, FMT='(1P3G12.5)') qTopoOut(1), MINVAL(qTopoOut(2:)), &
                                         MAXVAL(qTopoOut(2:))
  ELSE
    IF (SmoothTopo) THEN
      WRITE (UNIT=nfprt, FMT='(/,A)') ' Chopping Old Topography for Smoothing'
      CALL Chop (qTopoInp(:), qTopoOut(:), MendInp, MendOut, Kdim)
      WRITE (UNIT=nfprt, FMT='(/,A)') ' TopoOut after Chopping for Smoothing:'
      WRITE (UNIT=nfprt, FMT='(1P3G12.5)') qTopoOut(1), MINVAL(qTopoOut(2:)), &
                                           MAXVAL(qTopoOut(2:))
    END IF
  END IF
  IF (SmoothTopo) THEN
    WRITE (UNIT=nfprt, FMT='(/,A)') ' Smoothing Topography'
    CALL SmoothCoef (qTopoOut, MendOut, Kdim, MendCut)
    WRITE (UNIT=nfprt, FMT='(/,A)') ' TopoOut after Smoothing:'
    WRITE (UNIT=nfprt, FMT='(1P3G12.5)') qTopoOut(1), MINVAL(qTopoOut(2:)), &
                                         MAXVAL(qTopoOut(2:))
  END IF

  INQUIRE (FILE=TRIM(DirSig)//TRIM(DataSig), EXIST=GetNewSig)
  IF (GetNewSig) THEN
    WRITE (UNIT=nfprt, FMT='(/,A)') ' Getting New Delta Sigma'
    OPEN  (UNIT=nfsig, FILE=TRIM(DirSig)//TRIM(DataSig), FORM='FORMATTED', &
          ACCESS='SEQUENTIAL', ACTION='READ', STATUS='OLD', IOSTAT=ios)
    IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') ' ** (Error) ** Open file ', &
                                        TRIM(TRIM(DirSIg)//TRIM(DataSig)), &
                                        ' returned IOStat = ', ios
      STOP ' ** (Error) **'
    END IF
    READ  (UNIT=nfsig, FMT='(5F9.6)') DelSigmaOut
    CLOSE (UNIT=nfsig)
    IF (KmaxOut == KmaxInp .AND. MAXVAL(ABS(DelSigmaOut-DelSigmaInp)) < 1.0E-04_r8) THEN
      WRITE (UNIT=nfprt, FMT='(/,A)') ' KmaxOut = KmaxInp And DelSima Is Quite The Same'
      WRITE (UNIT=nfprt, FMT='(A,1PG12.5,/)') ' MAXVAL(ABS(DelSigmaOut-DelSigmaInp)) = ', &
                                                MAXVAL(ABS(DelSigmaOut-DelSigmaInp))
      SigInterOut=SigInterInp
      SigLayerOut=SigLayerInp
      DelSigmaOut=DelSigmaInp
    ELSE
      CALL NewSigma
    END IF
  ELSE
  IF (KmaxOut /= KmaxInp) THEN
    WRITE (UNIT=nfprt, FMT='(A,/)')   ' Error In Getting New Sigma: KmaxInp /= KmaxOut'
    WRITE (UNIT=nfprt, FMT='(2(A,I5))') ' KmaxInp = ', KmaxInp, ' KmaxOut = ', KmaxOut
    STOP
  END IF
  SigInterOut=SigInterInp
  SigLayerOut=SigLayerInp
  DelSigmaOut=DelSigmaInp
  END IF

  IF (KmaxOut == KmaxInp .AND. &
     (.NOT.GetNewTop .AND. .NOT.SmoothTopo)) VerticalInterp = .FALSE.

  IF (VerticalInterp .OR. GrADS) THEN

    CALL ICRecomposition

    IF (GrADS) THEN
      WRITE (Tdef(1: 2),'(I2.2)') DateCurrent(1)
      WRITE (Tdef(4: 5),'(I2.2)') DateCurrent(3)
      WRITE (Tdef(6: 8),'(A3)')   MonChar(DateCurrent(2))
      WRITE (Tdef(9:12),'(I4.4)') DateCurrent(4)
      CALL GetGrADSInp
      IF (GrADSOnly) STOP ' GrADS Only'
    END IF

    CALL ClsMemRecomposition

  END IF

  IF (VerticalInterp) THEN

  WRITE (UNIT=nfprt, FMT='(/,A)') 'Doing Vertical Interpolation'

  DO k=1,KmaxInp
    DO j=1,JmaxOut
      DO I=1,ImaxOut
        gPresInp(i,j,k)=gPsfcInp(i,j)*SigLayerInp(k)
      END DO
    END DO
  END DO

  CALL NewPs

  DO k=1,KmaxOut
    DO j=1,JmaxOut
      DO I=1,ImaxOut
        gPresOut(i,j,k)=gPsfcOut(i,j)*SigLayerOut(k)
      END DO
    END DO
  END DO
  gTvirInp=gTvirInp/(1.0_r8+cTv*gSpHuInp(:,:,:,1))
  write (*,*) ' pInp: ', minval(gPresInp(:,:,:)), maxval(gPresInp(:,:,:)) 
  write (*,*) ' pOut: ', minval(gPresOut(:,:,:)), maxval(gPresOut(:,:,:)) 
  CALL VertSigmaInter (ImaxOut, JmaxOut, &
                       KmaxInp, KmaxOut, NTracers, &
                       gPresInp(:,:,:), gUvelInp(:,:,:), gVvelInp(:,:,:), &
                       gTvirInp(:,:,:), gSpHuInp(:,:,:,:), &
                       gPresOut(:,:,:), gUvelOut(:,:,:), gVvelOut(:,:,:), &
                       gTvirOut(:,:,:), gSpHuOut(:,:,:,:))

  gTvirInp=gTvirInp*(1.0_r8+cTv*gSpHuInp(:,:,:,1))
  gTvirOut=gTvirOut*(1.0_r8+cTv*gSpHuOut(:,:,:,1))

  CALL ICDecomposition

  ELSE

  WRITE (UNIT=nfprt, FMT='(A,/)') 'Doing Chopping'

  IF (.NOT.GetNewTop .AND. .NOT.SmoothTopo) &
  CALL Chop (qTopoInp(:), qTopoOut(:), MendInp, MendOut, Kdim)
  CALL Chop (qLnPsInp(:), qLnPsOut(:), MendInp, MendOut, Kdim)
  CALL Chop (qTvirInp, qTvirOut, MendInp, MendOut, KmaxOut)
  CALL Chop (qVortInp, qVortOut, MendInp, MendOut, KmaxOut)
  CALL Chop (qDivgInp, qDivgOut, MendInp, MendOut, KmaxOut)
  DO nt=1,NTracers
    CALL  Chop (qSpHuInp(:,:,nt), qSpHuOut(:,:,nt), MendInp, MendOut, KmaxOut)
  END DO

  END IF

  CALL ICWrite

  CALL ClsArrays

  CALL ClsSpHuTracers

PRINT *, "*** Chopping ENDS NORMALLY ***"

CONTAINS


SUBROUTINE GDAStoGANL

  IMPLICIT NONE

  CHARACTER (LEN=1), DIMENSION (32) :: Descriptor

  WRITE (UNIT=nfprt, FMT='(/,A)') ' Getting GANL from GDAS NCEP File'

  OPEN (UNIT=nfnmc, FILE=TRIM(DGDInp)//TRIM(GDASInp), FORM='UNFORMATTED', &
        ACCESS='SEQUENTIAL', ACTION='READ', STATUS='OLD', IOSTAT=ios)
  IF (ios /= 0) THEN
    WRITE (UNIT=nferr, FMT='(3A,I4)') ' ** (Error) ** Open file ', &
                                      TRIM(TRIM(DGDInp)//TRIM(GDASInp)), &
                                      ' returned IOStat = ', ios
    STOP ' ** (Error) **'
  END IF

  OPEN (UNIT=nfcpt, FILE=TRIM(DirInp)//TRIM(DataInp), FORM='UNFORMATTED', &
        ACCESS='SEQUENTIAL', ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
  IF (ios /= 0) THEN
    WRITE (UNIT=nferr, FMT='(3A,I4)') ' ** (Error) ** Open file ', &
                                      TRIM(TRIM(DirInp)//TRIM(DataInp)), &
                                      ' returned IOStat = ', ios
    STOP ' ** (Error) **'
  END IF

  ! Descriptor Label (See NMC Office Note 85) (Not Used at CPTEC)
  READ (UNIT=nfnmc) Descriptor

  ! Forecast Time     - TimeOfDay
  ! Initial Hour      - DateInitial(1)
  ! Initial Month     - DateInitial(2)
  ! Initial Day       - DateInitial(3)
  ! Initial Year      - DateInitial(4)
  ! Sigma Interfaces  - SigInterInp
  ! Sigma Layers      - SigLayerInp
  READ (UNIT=nfnmc) TimeOfDay, DateInitial, SigIInp, SigLInp

  IF (ANY(SigIInp < 0.0_r8 .OR. SigIInp > 1.0_r8)) THEN
    WRITE (UNIT=nferr, FMT='(/,A)') ' SigI and SIgLi will be recalculated based on DelSInp'
    INQUIRE (FILE=TRIM(DirSig)//TRIM(DataSigInp), EXIST=GetNewSig)
    IF (GetNewSig) THEN
      WRITE (UNIT=nfprt, FMT='(/,A)') ' Getting New Delta Sigma'
      OPEN  (UNIT=nfsig, FILE=TRIM(DirSig)//TRIM(DataSigInp), FORM='FORMATTED', &
            ACCESS='SEQUENTIAL', ACTION='READ', STATUS='OLD', IOSTAT=ios)
      IF (ios /= 0) THEN
        WRITE (UNIT=nferr, FMT='(3A,I4)') ' ** (Error) ** Open file ', &
                                          TRIM(TRIM(DirSIg)//TRIM(DataSigInp)), &
                                          ' returned IOStat = ', ios
        STOP ' ** (Error) **'
      END IF
      READ  (UNIT=nfsig, FMT='(5F9.6)') DelSInp
      CLOSE (UNIT=nfsig)
      CALL SigmaInp
    ELSE
      WRITE (UNIT=nferr, FMT='(A)') ' There is no file : '//TRIM(DirSig)//TRIM(DataSigInp)
    END IF
  END IF

  DateCurrent=DateInitial
  ! Forecast Day      - ForecastDay
  ! Time of Day       - TimeOfDay
  ! Initial Date      - DateInitial
  ! Current Date      - DateCurrent
  WRITE (UNIT=nfcpt) ForecastDay, TimeOfDay, DateInitial, &
                     DateCurrent, SigIInp, SigLInp

  ! Spectral Coefficients of Orography (m)
  READ  (UNIT=nfnmc) qWorkInp
  WRITE (UNIT=nfcpt) qWorkInp

  ! Spectral coefficients of ln(Ps) (ln(hPa)/10)
  READ  (UNIT=nfnmc) qWorkInp
  WRITE (UNIT=nfcpt) qWorkInp

  ! Spectral Coefficients of Virtual Temp (K)
  DO k=1,KmaxInp
    READ  (UNIT=nfnmc) qWorkInp
    WRITE (UNIT=nfcpt) qWorkInp
  END DO

  ! Spectral Coefficients of Divergence and Vorticity (1/seg)
  DO k=1,KmaxInp
    READ  (UNIT=nfnmc) qWorkInp
    WRITE (UNIT=nfcpt) qWorkInp
    READ  (UNIT=nfnmc) qWorkInp
    WRITE (UNIT=nfcpt) qWorkInp
  END DO

  ! Spectral Coefficients of Specific Humidity (g/g)
  DO k=1,KmaxInp
    READ  (UNIT=nfnmc) qWorkInp
    WRITE (UNIT=nfcpt) qWorkInp
  END DO

  CLOSE(UNIT=nfcpt)

  IF (GetOzone) THEN

    ! Spectral Coefficients of Ozone (?)
    OPEN (UNIT=nfozw, FILE=TRIM(DirInp)//TRIM(OzonInp), FORM='UNFORMATTED', &
          ACCESS='SEQUENTIAL', ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
    IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') ' ** (Error) ** Open file ', &
                                        TRIM(TRIM(DirInp)//TRIM(OzonInp)), &
                                        ' returned IOStat = ', ios
      STOP ' ** (Error) **'
    END IF

    DO k=1,KmaxInp
      READ  (UNIT=nfnmc) qWorkInp
      WRITE (UNIT=nfozw) qWorkInp
    END DO
    CLOSE(UNIT=nfozw)
    NTracers=NTracers+1

    IF (GetTracers) THEN
      ! Spectral Coefficients of Tracers (?)
      OPEN (UNIT=nftrw, FILE=TRIM(DirInp)//TRIM(TracInp), FORM='UNFORMATTED', &
            ACCESS='SEQUENTIAL', ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
        WRITE (UNIT=nferr, FMT='(3A,I4)') ' ** (Error) ** Open file ', &
                                          TRIM(TRIM(DirInp)//TRIM(TracInp)), &
                                          ' returned IOStat = ', ios
        STOP ' ** (Error) **'
      END IF

      ios=0
      Tracer: DO
        DO k=1,KmaxInp
          READ (UNIT=nfnmc, IOSTAT=ios) qWorkInp
          IF (ios /= 0) THEN
            IF (ios == -1) THEN
              WRITE (UNIT=nfprt, FMT='(/,A,I5,A)') ' End of File Found - NTracers = ', &
                                                     NTracers, '  in:'
            ELSE
              WRITE (UNIT=nfprt, FMT='(/,A,I5,A)') ' Reading Error - ios = ', ios, '  in:'
            END IF
            WRITE (UNIT=nfprt, FMT='(1X,A,/)') TRIM(DGDInp)//TRIM(GDASInp)
            EXIT Tracer
          END IF
          WRITE (UNIT=nftrw) qWorkInp
        END DO
        NTracers=NTracers+1
      END DO Tracer
      CLOSE(UNIT=nftrw)
    END IF

  END IF

  CLOSE(UNIT=nfnmc)

  IF (GDASOnly) STOP ' GDASOnly = .TRUE. '

END SUBROUTINE GDAStoGANL


SUBROUTINE ICRead

  IMPLICIT NONE

  OPEN (UNIT=nficr, FILE=TRIM(DirInp)//TRIM(DataInp), FORM='UNFORMATTED', &
        ACCESS='SEQUENTIAL', ACTION='READ', STATUS='OLD', IOSTAT=ios)
  IF (ios /= 0) THEN
    WRITE (UNIT=nferr, FMT='(3A,I4)') ' ** (Error) ** Open file ', &
                                      TRIM(TRIM(DirInp)//TRIM(DataInp)), &
                                      ' returned IOStat = ', ios
    STOP ' ** (Error) **'
  END IF

  READ (UNIT=nficr) ForecastDay, TimeOfDay, DateInitial, &
                    DateCurrent, SigIInp, SigLInp

  SigInterInp=SigIInp
  SigLayerInp=SigLInp
  DO k=1,KmaxInp
     DelSigmaInp(k)=SigInterInp(k)-SigInterInp(k+1)
  END DO

  WRITE (UNIT=nfprt, FMT='(/,A,I5,A,F15.4)') ' ForecastDay = ', ForecastDay, &
                                             ' TimeOfDay = ', TimeOfDay
  WRITE (UNIT=nfprt, FMT='(/,A,4I5)') ' DateInitial = ', DateInitial
  WRITE (UNIT=nfprt, FMT='(/,A,4I5)') ' DateCurrent = ', DateCurrent
  WRITE (UNIT=nfprt, FMT='(/,A)')  ' DelSigmaInp:'
  WRITE (UNIT=nfprt, FMT='(7F10.6)') DelSigmaInp
  WRITE (UNIT=nfprt, FMT='(/,A)')  ' SigInterInp:'
  WRITE (UNIT=nfprt, FMT='(7F10.6)') SigInterInp
  WRITE (UNIT=nfprt, FMT='(/,A)')  ' SigLayerInp:'
  WRITE (UNIT=nfprt, FMT='(7F10.6)') SigLayerInp

  READ (UNIT=nficr) qWorkInp
  qTopoInp(:)=qWorkInp
  WRITE (UNIT=nfprt, FMT='(/,A)') ' TopoInp:'
  WRITE (UNIT=nfprt, FMT='(I5,1P3G12.5)') 0, qTopoInp(1), &
                                          MINVAL(qTopoInp(2:)), &
                                          MAXVAL(qTopoInp(2:))

  READ (UNIT=nficr) qWorkInp
  qLnPsInp(:)=qWorkInp
  WRITE (UNIT=nfprt, FMT='(/,A)') ' LnPsInp:'
  WRITE (UNIT=nfprt, FMT='(I5,1P3G12.5)') 0, qLnPsInp(1), &
                                         MINVAL(qLnPsInp(2:)), &
                                         MAXVAL(qLnPsInp(2:))

  WRITE (UNIT=nfprt, FMT='(/,A)') ' TvirInp:'
  DO k=1,KmaxInp
    READ (UNIT=nficr) qWorkInp
    qTvirInp(:,k)=qWorkInp
    WRITE (UNIT=nfprt, FMT='(I5,1P3G12.5)') k, qTvirInp(1,k), &
                                            MINVAL(qTvirInp(2:,k)), &
                                            MAXVAL(qTvirInp(2:,k))
  END DO

  WRITE (UNIT=nfprt, FMT='(/,A)') ' DivgInp - VortInp:'
  DO k=1,KmaxInp
    READ (UNIT=nficr) qWorkInp
    qDivgInp(:,k)=qWorkInp
    WRITE (UNIT=nfprt, FMT='(I5,1P3G12.5)') k, qDivgInp(1,k), &
                                            MINVAL(qDivgInp(2:,k)), &
                                            MAXVAL(qDivgInp(2:,k))
    READ (UNIT=nficr) qWorkInp
    qVortInp(:,k)=qWorkInp
    WRITE (UNIT=nfprt, FMT='(I5,1P3G12.5)') k, qVortInp(1,k), &
                                            MINVAL(qVortInp(2:,k)), &
                                            MAXVAL(qVortInp(2:,k))
  END DO

  WRITE (UNIT=nfprt, FMT='(/,A)') ' SpHuInp:'
  DO k=1,KmaxInp
    READ (UNIT=nficr) qWorkInp
    qSpHuInp(:,k,1)=qWorkInp
    WRITE (UNIT=nfprt, FMT='(I5,1P3G12.5)') k, qSpHuInp(1,k,1), &
                                            MINVAL(qSpHuInp(2:,k,1)), &
                                            MAXVAL(qSpHuInp(2:,k,1))
  END DO

  CLOSE (UNIT=nficr)

  IF (GetOzone) THEN

    OPEN (UNIT=nfozr, FILE=TRIM(DirInp)//TRIM(OzonInp), FORM='UNFORMATTED', &
          ACCESS='SEQUENTIAL', ACTION='READ', STATUS='OLD', IOSTAT=ios)
    IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') ' ** (Error) ** Open file ', &
                                        TRIM(TRIM(DirInp)//TRIM(OzonInp)), &
                                        ' returned IOStat = ', ios
      STOP ' ** (Error) **'
    END IF

    nt=2
    DO k=1,KmaxInp
      READ (UNIT=nfozr) qWorkInp
      qSpHuInp(:,k,nt)=qWorkInp
    END DO
    CLOSE(UNIT=nfozr)

    IF (GetTracers) THEN
      OPEN (UNIT=nftrr, FILE=TRIM(DirInp)//TRIM(TracInp), FORM='UNFORMATTED', &
            ACCESS='SEQUENTIAL', ACTION='READ', STATUS='OLD', IOSTAT=ios)
      IF (ios /= 0) THEN
        WRITE (UNIT=nferr, FMT='(3A,I4)') ' ** (Error) ** Open file ', &
                                          TRIM(TRIM(DirInp)//TRIM(TracInp)), &
                                          ' returned IOStat = ', ios
        STOP ' ** (Error) **'
      END IF
      Do nt=3,NTracers
        DO k=1,KmaxInp
          READ (UNIT=nftrr) qWorkInp
          qSpHuInp(:,k,nt)=qWorkInp
        END DO
      END DO
      CLOSE(UNIT=nftrr)
    END IF

  END IF

END SUBROUTINE ICRead


SUBROUTINE ICWrite

  IMPLICIT NONE

  OPEN (UNIT=nficw, FILE=TRIM(DirOut)//TRIM(DataOut), FORM='UNFORMATTED', &
        ACCESS='SEQUENTIAL', ACTION='WRITE', STATUS='UNKNOWN', IOSTAT=ios)
  IF (ios /= 0) THEN
    WRITE (UNIT=nferr, FMT='(3A,I4)') ' ** (Error) ** Open file ', &
                                      TRIM(TRIM(DirOut)//TRIM(DataOut)), &
                                      ' returned IOStat = ', ios
    STOP ' ** (Error) **'
  END IF

  SigIOut=SigInterOut
  SigLOut=SigLayerOut

  WRITE (UNIT=nficw) ForecastDay, TimeOfDay, DateInitial, &
                     DateCurrent, SigIOut, SigLOut

  WRITE (UNIT=nfprt, FMT='(/,A,I5,A,F15.4)') ' ForecastDay = ', ForecastDay, &
                                             ' TimeOfDay = ', TimeOfDay
  WRITE (UNIT=nfprt, FMT='(/,A,4I5)') ' DateInitial = ', DateInitial
  WRITE (UNIT=nfprt, FMT='(/,A,4I5)') ' DateCurrent = ', DateCurrent
  WRITE (UNIT=nfprt, FMT='(/,A)')  ' DelSigmaOut:'
  WRITE (UNIT=nfprt, FMT='(7F10.6)') DelSigmaOut
  WRITE (UNIT=nfprt, FMT='(/,A)')  ' SigInterOut:'
  WRITE (UNIT=nfprt, FMT='(7F10.6)') SigInterOut
  WRITE (UNIT=nfprt, FMT='(/,A)')  ' SigLayerOut:'
  WRITE (UNIT=nfprt, FMT='(7F10.6)') SigLayerOut

  qWorkOut=qTopoOut(:)
  WRITE (UNIT=nficw) qWorkOut
  WRITE (UNIT=nfprt, FMT='(/,A)') ' TopoOut:'
  WRITE (UNIT=nfprt, FMT='(I5, 1P3G12.5)') 0, qTopoOut(1), &
                                           MINVAL(qTopoOut(2:)), &
                                           MAXVAL(qTopoOut(2:))

  qWorkOut=qLnPsOut(:)
  WRITE (UNIT=nficw) qWorkOut
  WRITE (UNIT=nfprt, FMT='(/,A)') ' LnPsOut:'
  WRITE (UNIT=nfprt, FMT='(I5,1P3G12.5)') 0, qLnPsOut(1), &
                                          MINVAL(qLnPsOut(2:)), &
                                          MAXVAL(qLnPsOut(2:))

  WRITE (UNIT=nfprt, FMT='(/,A)') ' TvirOut:'
  DO k=1,KmaxOut
    qWorkOut=qTvirOut(:,k)
    WRITE (UNIT=nficw) qWorkOut
    WRITE (UNIT=nfprt, FMT='(I5,1P3G12.5)') k, qTvirOut(1,k), &
                                            MINVAL(qTvirOut(2:,k)), &
                                            MAXVAL(qTvirOut(2:,k))
  END DO

  WRITE (UNIT=nfprt, FMT='(/,A)') ' DivgOut - VortOut:'
  DO k=1,KmaxOut
    qWorkOut=qDivgOut(:,k)
    WRITE (UNIT=nficw) qWorkOut
    WRITE (UNIT=nfprt, FMT='(I5,1P3G12.5)') k, qDivgOut(1,k), &
                                            MINVAL(qDivgOut(2:,k)), &
                                            MAXVAL(qDivgOut(2:,k))
    qWorkOut=qVortOut(:,k)
    WRITE (UNIT=nficw) qWorkOut
    WRITE (UNIT=nfprt, FMT='(I5,1P3G12.5)') k, qVortOut(1,k), &
                                            MINVAL(qVortOut(2:,k)), &
                                            MAXVAL(qVortOut(2:,k))
  END DO

  WRITE (UNIT=nfprt, FMT='(/,A)') ' SpHuOut:'
  DO k=1,KmaxOut
    qWorkOut=qSpHuOut(:,k,1)
    WRITE (UNIT=nficw) qWorkOut
    WRITE (UNIT=nfprt, FMT='(I5,1P3G12.5)') k, qSpHuOut(1,k,1), &
                                            MINVAL(qSpHuOut(2:,k,1)), &
                                            MAXVAL(qSpHuOut(2:,k,1))
  END DO
  WRITE (UNIT=nfprt, FMT='(A)') ' '

  CLOSE(UNIT=nficw)

  IF (GetOzone) THEN

    CALL InitRecomposition (MendOut, ImaxOut, JmaxOut, KmaxOut)

    ALLOCATE (gOzonOut (ImaxOut,JmaxOut,KmaxOut))
    nt=2
    CALL RecompositionScalar (KmaxOut, qSpHuOut(:,:,nt), gOzonOut(:,:,:))

    INQUIRE (IOLENGTH=IOL) gWorkOut
    OPEN (UNIT=nfozg, FILE=TRIM(DirOut)//TRIM(OzonOut), FORM='UNFORMATTED', &
          ACCESS='DIRECT', RECL=IOL, ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
    IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') ' ** (Error) ** Open file ', &
                                        TRIM(TRIM(DirOut)//TRIM(OzonOut)), &
                                        ' returned IOStat = ', ios
      STOP ' ** (Error) **'
    END IF
    DO k=1,KmaxOut
      gWorkOut=gOzonOut(:,:,k)
      WRITE (UNIT=nfozg, REC=k) gWorkOut
    END DO
    CLOSE(UNIT=nfozg)
    DEALLOCATE (gOzonOut)

    IF (GetTracers) THEN

      ALLOCATE (gTracOut (ImaxOut,JmaxOut,KmaxOut,NTracers-2))
      DO nt=1,Ntracers-2
        CALL RecompositionScalar (KmaxOut, qSpHuOut(:,:,nt+2), gTracOut(:,:,:,nt))
      END DO

      OPEN (UNIT=nftrg, FILE=TRIM(DirOut)//TRIM(TracOut), FORM='UNFORMATTED', &
          ACCESS='SEQUENTIAL', ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
    IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') ' ** (Error) ** Open file ', &
                                        TRIM(TRIM(DirOut)//TRIM(TracOut)), &
                                        ' returned IOStat = ', ios
      STOP ' ** (Error) **'
    END IF
      DO nt=1,Ntracers-2
        DO k=1,KmaxOut
          gWorkOut=gTracOut(:,:,k,nt)
          WRITE (UNIT=nftrg) gWorkOut
        END DO
      END DO
      CLOSE(UNIT=nftrg)
      DEALLOCATE (gTracOut)
    END IF

    CALL ClsMemRecomposition

  END IF

END SUBROUTINE ICWrite


SUBROUTINE ICRecomposition

  IMPLICIT NONE

  CALL InitRecomposition (MendInp, ImaxOut, JmaxOut, KmaxInp)

  CALL RecompositionScalar (Kdim, qTopoInp(:), gTopoInp(:,:))
  IF (GetNewTop .OR. SmoothTopo) THEN
    WRITE (UNIT=nfprt, FMT='(/,A)') ' TopoOut before Chopping:'
    WRITE (UNIT=nfprt, FMT='(1P3G12.5)') qTopoOut(1), MINVAL(qTopoOut(2:)), &
                                         MAXVAL(qTopoOut(2:))
    CALL Chop (qTopoOut(:), qTopoRec(:), MendOut, MendInp, Kdim)
    WRITE (UNIT=nfprt, FMT='(/,A)') ' TopoRec after Chopping:'
    WRITE (UNIT=nfprt, FMT='(1P3G12.5)') qTopoRec(1), MINVAL(qTopoRec(2:)), &
                                         MAXVAL(qTopoRec(2:))
    CALL RecompositionScalar (Kdim, qTopoRec(:), gTopoOut(:,:))
    gTopoDel=gTopoOut-gTopoInp
    !gWorkOut=gTopoInp(:,:)
    !WRITE (UNIT=77) gWorkOut
    !gWorkOut=gTopoOut(:,:)
    !WRITE (UNIT=77) gWorkOut
    !gWorkOut=gTopoDel(:,:)
    !WRITE (UNIT=77) gWorkOut
  ELSE
    gTopoOut=gTopoInp
  END IF

  CALL RecompositionScalar (Kdim, qLnPsInp(:), gLnPsInp(:,:))
  gPsfcInp=10.0_r8*EXP(gLnPsInp)

  CALL RecompositionScalar (KmaxInp, qTvirInp, gTvirInp)

  CALL RecompositionScalar (KmaxInp, qDivgInp, gDivgInp)

  CALL RecompositionScalar (KmaxInp, qVortInp, gVortInp)

  DO nt=1,Ntracers
    CALL RecompositionScalar (KmaxInp, qSpHuInp(:,:,nt), gSpHuInp(:,:,:,nt))
  END DO

  CALL DivgVortToUV (qDivgInp, qVortInp, qUvelInp, qVvelInp)
  CALL RecompositionVector (KmaxInp, qUvelInp, gUvelInp)
  DO k=1,KmaxInp
    DO j=1,JmaxOut
      DO i=1,ImaxOut
        gUvelInp(i,j,k)=gUvelInp(i,j,k)/coslat(j)
      END DO
    END DO
  END DO
  CALL RecompositionVector (KmaxInp, qVvelInp, gVvelInp)
  DO k=1,KmaxInp
    DO j=1,JmaxOut
      DO i=1,ImaxOut
        gVvelInp(i,j,k)=gVvelInp(i,j,k)/coslat(j)
      END DO
    END DO
  END DO

END SUBROUTINE ICRecomposition


SUBROUTINE ICDecomposition

  IMPLICIT NONE

  CALL InitFFT
  CALL InitLegendre

  WRITE (UNIT=nfprt, FMT='(/,A)') ' TopoOut before Dec:'
  WRITE (UNIT=nfprt, FMT='(1P3G12.5)') qTopoOut(1), MINVAL(qTopoOut(2:)), &
                                       MAXVAL(qTopoOut(2:))
  CALL DectoSpherHarm (gTopoOut, qTopoOut, Kdim)
  WRITE (UNIT=nfprt, FMT='(/,A)') ' TopoOut after Dec:'
  WRITE (UNIT=nfprt, FMT='(1P3G12.5)') qTopoOut(1), MINVAL(qTopoOut(2:)), &
                                       MAXVAL(qTopoOut(2:))
  CALL TransSpherHarm (Kdim,    qTopoOut, ICaseDec)
  WRITE (UNIT=nfprt, FMT='(/,A)') ' TopoOut after Trans:'
  WRITE (UNIT=nfprt, FMT='(1P3G12.5)') qTopoOut(1), MINVAL(qTopoOut(2:)), &
                                       MAXVAL(qTopoOut(2:))

  CALL DectoSpherHarm (gLnPsOut, qLnPsOut, Kdim)
  CALL TransSpherHarm (Kdim,    qLnPsOut, ICaseDec)

  CALL DectoSpherHarm (gTvirOut, qTvirOut, KmaxOut)
  CALL TransSpherHarm (KmaxOut, qTvirOut, ICaseDec)

  DO nt=1,NTracers
    CALL DectoSpherHarm (gSpHuOut(:,:,:,nt), qSpHuOut(:,:,nt), KmaxOut)
    CALL TransSpherHarm (KmaxOut, qSpHuOut(:,:,nt), ICaseDec)
  END DO

  CALL UVtoDivgVort (gUvelOut, gVvelOut, qDivgOut, qVortOut)
  CALL TransSpherHarm (KmaxOut, qDivgOut, ICaseDec)
  CALL TransSpherHarm (KmaxOut, qVortOut, ICaseDec)

  CALL ClsMemFFT
  CALL ClsMemLegendre

END SUBROUTINE ICDecomposition


SUBROUTINE Chop (qInp, qOut, MendI, MendO, Kmax)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: MendI, MendO, Kmax

  REAL (KIND=r8), DIMENSION ((MendI+1)*(MendI+2), Kmax), INTENT(IN) :: qInp

  REAL (KIND=r8), DIMENSION ((MendO+1)*(MendO+2), Kmax), INTENT(Out) :: qOut

  INTEGER :: n, ma, mb, mai, mbi, maf, mbf, maif

  IF (Kmax /= 1 .AND. KmaxInp /= KmaxOut) THEN
    WRITE (UNIT=nfprt, FMT='(A,/)')   ' Error In SUBROUTINE Chop: KmaxInp /= KmaxOut'
    WRITE (UNIT=nfprt, FMT='(2(A,I5))') ' KmaxInp = ', KmaxInp, ' KmaxOut = ', KmaxOut
    STOP
  END IF

  qOut=0.0_r8

  IF (MendO == MendI) THEN

    WRITE (UNIT=nfprt, FMT='(/,A,1PG12.5)') &
          ' SUBROUTINE Chop: MendO == MendI, qInp(1,1) = ', qInp(1,1)
    qOut=qInp

  ELSE
 
    IF (MendO < MendI) THEN
      WRITE (UNIT=nfprt, FMT='(/,A,1PG12.5)') &
            ' SUBROUTINE Chop: MendO < MendI, qInp(1,1) = ', qInp(1,1)
      DO k=1,Kmax
        mbf=0
        maif=0
        DO n=1,MendO+1
          mb=2*(MendO+2-n)
          ma=2*(MendI+2-n)
          mbi=mbf+1
          mbf=mbi+mb-1
          mai=maif+1
          maif=mai+ma-1
          maf=maif-(ma-mb)
          qOut(mbi:mbf,k)=qInp(mai:maf,k)
        END DO
      END DO
    ELSE
      WRITE (UNIT=nfprt, FMT='(/,A,1PG12.5)') &
            ' SUBROUTINE Chop: MendO > MendI, qInp(1,1) = ', qInp(1,1)
      DO k=1,Kmax
        mbf=0
        maif=0
        DO n=1,MendI+1
          ma=2*(MendO+2-n)
          mb=2*(MendI+2-n)
          mai=maif+1
          maif=mai+ma-1
          maf=maif-(ma-mb)
          mbi=mbf+1
          mbf=mbi+mb-1
          qOut(mai:maf,k)=qInp(mbi:mbf,k)
        END DO
      END DO
    END IF

  END IF

END SUBROUTINE Chop


SUBROUTINE SmoothCoef (qOut, MendO, Kmax, MendC)

! Smoothes Spherical Harmonics Coefficients
!          Using Hoskin's(?) Filter

  IMPLICIT NONE

  INTEGER, INTENT (IN) :: MendO, Kmax, MendC

  REAL (KIND=r8), DIMENSION ((MendO+1)*(MendO+2),Kmax), INTENT (IN OUT) :: qOut

  INTEGER :: k, mn, nn, mm, Mmax, mns

  REAL (KIND=r8) :: rm, rn, cx, red, rmn, ck

  REAL (KIND=r8), DIMENSION ((MendO+1)*(MendO+2)) :: Weight

  INTEGER :: icall=0
  SAVE :: icall

  WRITE (UNIT=nfprt, FMT='(/,A,I5)') ' Doing Smoothing : ', icall

  IF (icall == 0) THEN
    rm=REAL(2*MendC,r8)*REAL(2*MendC+1,r8)
    rn=REAL(MendC-1,r8)*REAL(MendC,r8)
    red=(SmthPerCut)**(-(rm*rm)/(rn*rn)/Iter)
    WRITE (UNIT=nfprt, FMT='(/,A,F10.4)') ' From SmoothCoef: Red = ', red
    cx=-LOG(red)/(rm*rm)
    mn=0
    DO nn=1,MendO+1
      Mmax=MendO-nn+2
      DO mm=1,Mmax
        mn=mn+1
        mns=mm+nn-1
        rmn=REAL(mns-1,r8)*REAL(mns,r8)
        ck=(EXP(cx*rmn*rmn))**Iter
        Weight(2*mn-1)=ck
        Weight(2*mn)=ck
      END DO
    END DO
    icall=icall+1
  END IF

  DO k=1,Kmax
    DO mn=1,(MendO+1)*(MendO+2)
      qOut(mn,k)=qOut(mn,k)*Weight(mn)
    END DO
  END DO

END SUBROUTINE SmoothCoef


SUBROUTINE GetGrADSInp

  IMPLICIT NONE

  WRITE (UNIT=nfprt, FMT='(/,A,L6,/)') ' GrADS     = ', GrADS
  WRITE (UNIT=nfprt, FMT='(/,A,L6,/)') ' GrADSOnly = ', GrADSOnly

  INQUIRE (IOLENGTH=IOL) gWorkOut
  OPEN (UNIT=nfgrd, FILE=TRIM(DirGrd)//TRIM(DataOut)//'.GrADS', FORM='UNFORMATTED', &
        ACCESS='DIRECT', RECL=IOL, ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
  IF (ios /= 0) THEN
    WRITE (UNIT=nferr, FMT='(3A,I4)') ' ** (Error) ** Open file ', &
                                      TRIM(TRIM(DirGrd)//TRIM(DataOut))//'.GrADS', &
                                      ' returned IOStat = ', ios
    STOP ' ** (Error) **'
  END IF
  gWorkOut=gTopoOut(:,:)
  nRec=1
  WRITE (UNIT=nfgrd, REC=nRec) gWorkOut
  gWorkOut=gPsfcInp(:,:)
  nRec=2
  WRITE (UNIT=nfgrd, REC=nRec) gWorkOut
  DO k=1,KmaxInp
    gWorkOut=gTvirInp(:,:,k)
    nRec=2+INT(k)
    WRITE (UNIT=nfgrd, REC=nRec) gWorkOut
  END DO
  DO k=1,KmaxInp
    gWorkOut=gDivgInp(:,:,k)
    nRec=2+INT(k+KmaxInp)
    WRITE (UNIT=nfgrd, REC=nRec) gWorkOut
  END DO
  DO k=1,KmaxInp
    gWorkOut=gVortInp(:,:,k)
    nRec=2+INT(k+2*KmaxInp)
    WRITE (UNIT=nfgrd, REC=nRec) gWorkOut
  END DO
  DO k=1,KmaxInp
    gWorkOut=gSpHuInp(:,:,k,1)
    nRec=2+INT(k+3*KmaxInp)
    WRITE (UNIT=nfgrd, REC=nRec) gWorkOut
  END DO
  DO k=1,KmaxInp
    gWorkOut=gUvelInp(:,:,k)
    nRec=2+INT(k+4*KmaxInp)
    WRITE (UNIT=nfgrd, REC=nRec) gWorkOut
  END DO
  DO k=1,KmaxInp
    gWorkOut=gVvelInp(:,:,k)
    nRec=2+INT(k+5*KmaxInp)
    WRITE (UNIT=nfgrd, REC=nRec) gWorkOut
  END DO
  IF (GetOzone) THEN
    nt=2
    DO k=1,KmaxInp
      gWorkOut=gSpHuInp(:,:,k,nt)
      nRec=2+INT(k+6*KmaxInp)
      WRITE (UNIT=nfgrd, REC=nRec) gWorkOut
    END DO
  END IF
  IF (GetTracers) THEN
    DO nt=3,NTracers
      DO k=1,KmaxInp
        gWorkOut=gSpHuInp(:,:,k,nt)
        nRec=2+INT(k+(4+nt)*KmaxInp)
        WRITE (UNIT=nfgrd, REC=nRec) gWorkOut
      END DO
    END DO
  END IF
  CLOSE (UNIT=nfgrd)

  OPEN (UNIT=nfctl, FILE=TRIM(DirGrd)//TRIM(DataOut)//'.GrADS.ctl', FORM='FORMATTED', &
        ACCESS='SEQUENTIAL', ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
  IF (ios /= 0) THEN
    WRITE (UNIT=nferr, FMT='(3A,I4)') ' ** (Error) ** Open file ', &
                                      TRIM(TRIM(DirGrd)//TRIM(DataOut))//'.GrADS.ctl', &
                                      ' returned IOStat = ', ios
    STOP ' ** (Error) **'
  END IF
  WRITE (UNIT=nfctl, FMT='(A)') 'dset ^'//TRIM(DataOut)//'.GrADS'
  WRITE (UNIT=nfctl, FMT='(A)') 'options yrev big_endian'
  WRITE (UNIT=nfctl, FMT='(A)') 'undef 1e20'
  WRITE (UNIT=nfctl, FMT='(A,I5,A,2F12.6)') 'xdef ',ImaxOut,' linear ', &
                                            0.0_r8, 360.0_r8/REAL(ImaxOut,r8)
  WRITE (UNIT=nfctl, FMT='(A,I5,A)') 'ydef ',JmaxOut,' levels '
  WRITE (UNIT=nfctl, FMT='(6F10.3)') glat(JmaxOut:1:-1)
  WRITE (UNIT=nfctl, FMT='(A,I5,A)') 'zdef ',KmaxInp,' levels '
  WRITE (UNIT=nfctl, FMT='(6F10.3)') 1000.0_r8*SigLayerInp
  WRITE (UNIT=nfctl, FMT='(3A)') 'tdef 1 linear ', Tdef,' 1dy'
  IF (NTracers == 1) THEN
    WRITE (UNIT=nfctl, FMT='(A)') 'vars 8'
  ELSE
    WRITE (UNIT=nfctl, FMT='(A,I5)') 'vars ',7+NTracers
  END IF
  WRITE (UNIT=nfctl, FMT='(A)') 'topo   0 99 Topography        '//TruncInp//' (m)'
  WRITE (UNIT=nfctl, FMT='(A)') 'pslc   0 99 Surface Pressure  '//TruncInp//' (hPa)'
  WRITE (UNIT=nfctl, FMT='(A,I3,A)') 'tvir ',KmaxInp,' 99 Virt Temperature  '// &
                                     TruncInp//' (K)'
  WRITE (UNIT=nfctl, FMT='(A,I3,A)') 'divg ',KmaxInp,' 99 Divergence        '// &
                                     TruncInp//' (1/s)'
  WRITE (UNIT=nfctl, FMT='(A,I3,A)') 'vort ',KmaxInp,' 99 Vorticity         '// &
                                     TruncInp//' (1/s)'
  WRITE (UNIT=nfctl, FMT='(A,I3,A)') 'umes ',KmaxInp,' 99 Specific Humidity '// &
                                     TruncInp//' (kg/kg)'
  WRITE (UNIT=nfctl, FMT='(A,I3,A)') 'uvel ',KmaxInp,' 99 Zonal Wind        '// &
                                     TruncInp//' (m/s)'
  WRITE (UNIT=nfctl, FMT='(A,I3,A)') 'vvel ',KmaxInp,' 99 Meridional Wind   '// &
                                     TruncInp//' (m/s)'
  IF (GetOzone) THEN
    WRITE (UNIT=nfctl, FMT='(A,I3,A)') 'ozon ',KmaxInp,' 99 Ozone             '// &
                                       TruncInp//' (?)'
  END IF
  IF (GetTracers) THEN
    DO nt=3,NTracers
      WRITE (UNIT=nfctl, FMT='(A,I1,A,I3,A)') 'trc',nt-2,' ',KmaxInp, &
                             ' 99 Tracer            '//TruncInp//' (?)'
    END DO
  END IF
  WRITE (UNIT=nfctl, FMT='(A)') 'endvars'
  CLOSE (UNIT=nfctl)

END SUBROUTINE GetGrADSInp


END PROGRAM Chopping
