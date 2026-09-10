!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
PROGRAM PorceClayMaskIBIS

   USE InputParameters, ONLY: r4, r8, Undef, UndefG, &
                              nferr, nfprt, nflsm, nfcvm, &
                              nfvgm, nflsv, nfout, nfctl, &
                              Imax, Jmax, Idim, Jdim,Layer, &
                              nLats, mskfmt, NameLSM, VarName, &
                              NameLSMSSiB, VarNameVeg, &
                              DirPreOut, DirModelIn, DirMain, &
                              GrADS, &
                              InitInputParameters

   USE AreaIntegerInterp, ONLY: gLats, InitAreaIntegerInterp, &
                                        DoAreaIntegerInterp

   IMPLICIT NONE
   ! Horizontal Areal Interpolator
   ! Interpolate Regular To Gaussian
   ! Regular Input Data is Assumed to be Oriented with
   ! the North Pole and Greenwich as the First Point

   ! Set Undefined Value for Input Data at Locations which
   ! are not to be Included in Interpolation

   INTEGER :: LRecIn, LRecOut, LRecGad, ios, &
              i, j, k,im, ip, jm, jp, kk, &
              k1, k2, k3, k4, k5, k6, k7, k8

   INTEGER, DIMENSION (:,:,:), ALLOCATABLE :: PorceClayMaskIBISIn
   INTEGER, DIMENSION (:,:,:), ALLOCATABLE :: PorceClayMaskIBISOut, VegMaskSave
   INTEGER, DIMENSION (:,:)  , ALLOCATABLE :: LandSeaMask, LSMaskSave
   INTEGER, DIMENSION (:,:,:), ALLOCATABLE :: PorceClayMaskIBISInput
   INTEGER, DIMENSION (:,:,:), ALLOCATABLE :: PorceClayMaskIBISOutput
   LOGICAL, PARAMETER :: dumpLocal=.FALSE.

   REAL (KIND=r4), DIMENSION (:,:,:), ALLOCATABLE :: Gad

   CALL InitInputParameters ()

   CALL InitAreaIntegerInterp ()

   ALLOCATE (PorceClayMaskIBISInput(Idim,Jdim,Layer))
   ALLOCATE (PorceClayMaskIBISOutput(Imax,Jmax,Layer))
   ALLOCATE (PorceClayMaskIBISIn(Idim,Jdim,Layer))
   ALLOCATE (PorceClayMaskIBISOut(Imax,Jmax,Layer), VegMaskSave(Imax,Jmax,Layer))
   ALLOCATE (LandSeaMask(Imax,Jmax), LSMaskSave(Imax,Jmax))
   ALLOCATE (Gad(Imax,Jmax,Layer))

   ! Land Sea Mask : Input
   OPEN (UNIT=nflsm, FILE=TRIM(DirMain)//DirPreOut//NameLSM//nLats, &
         FORM='FORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='READ', STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirPreOut//NameLSM//nLats, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   READ  (UNIT=nflsm, FMT=mskfmt) LandSeaMask
   LSMaskSave=LandSeaMask
   CLOSE (UNIT=nflsm)

   ! Read In Input PorceClayMaskIBIS
   INQUIRE (IOLENGTH=LRecIn) PorceClayMaskIBISIn
   OPEN (UNIT=nfcvm, FILE=TRIM(DirMain)//DirPreOut//TRIM(VarName)//'.dat', &
         FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIn, ACTION='READ', &
         STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirPreOut//TRIM(VarName)//'.dat', &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   READ  (UNIT=nfcvm, REC=1) PorceClayMaskIBISIn
   CLOSE (UNIT=nfcvm)

   ! Interpolate Input Regular PorceClayMaskIBIS To Gaussian Grid Output
   PorceClayMaskIBISInput=INT(PorceClayMaskIBISIn)
   DO k=1,Layer
      CALL DoAreaIntegerInterp (PorceClayMaskIBISInput(:,:,k), PorceClayMaskIBISOutput(:,:,k))
   END DO
   DO k=1,Layer
      DO j=1,jMax
         DO i=1,iMax
            IF(LSMaskSave(i,j) /= 0 .AND.PorceClayMaskIBISOutput(i,j,k) ==0.0 )PorceClayMaskIBISOutput(i,j,k)=60
         END DO
      END DO
   END DO

   VegMaskSave=INT(PorceClayMaskIBISOutput)
   
   ! Fix Problems on Vegetation Mask and Land Sea Mask
   CALL VegetationMaskCheck ()

   PorceClayMaskIBISOut=INT(PorceClayMaskIBISOutput)
 
   ! Write Out Adjusted Interpolated PorceClayMaskIBIS
   INQUIRE (IOLENGTH=LRecOut) PorceClayMaskIBISOut
   OPEN (UNIT=nfvgm, FILE=TRIM(DirMain)//DirModelIn//TRIM(VarNameVeg)//nLats, &
         FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecOut, ACTION='WRITE', &
         STATUS='REPLACE', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirModelIn//TRIM(VarNameVeg)//nLats, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   WRITE (UNIT=nfvgm, REC=1) PorceClayMaskIBISOut
   CLOSE (UNIT=nfvgm)

   ! Land Sea Mask : Output
   OPEN (UNIT=nflsv, FILE=TRIM(DirMain)//DirPreOut//NameLSMSSiB//nLats, &
         FORM='FORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirPreOut//NameLSMSSiB//nLats, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   WRITE (UNIT=nflsv, FMT=mskfmt) LandSeaMask
   CLOSE (UNIT=nflsv)

   IF (GrADS) THEN
      INQUIRE (IOLENGTH=LRecGad) Gad(:,:,1)
      OPEN (UNIT=nfout, FILE=TRIM(DirMain)//DirPreOut//TRIM(VarNameVeg)//nLats, &
            FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecGad, ACTION='WRITE', &
            STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE (UNIT=nferr, FMT='(3A,I4)') &
               ' ** (Error) ** Open file ', &
                 TRIM(DirMain)//DirPreOut//TRIM(VarNameVeg)//nLats, &
               ' returned IOStat = ', ios
         STOP  ' ** (Error) **'
      END IF
      Gad(:,:,1)=REAL(LSMaskSave,r4)
      WRITE (UNIT=nfout, REC=1) Gad(:,:,1)
      Gad(:,:,1)=REAL(LandSeaMask,r4)
      WRITE (UNIT=nfout, REC=2) Gad(:,:,1)
      DO k=1,Layer
         Gad(:,:,k)=REAL(VegMaskSave(:,:,k),r4)
         WRITE (UNIT=nfout, REC=2+k) Gad(:,:,k)
      END DO 
      DO k=1,Layer        
         Gad(:,:,k)=REAL(PorceClayMaskIBISOut(:,:,k),r4)
         WRITE (UNIT=nfout, REC=2+Layer+k) Gad(:,:,k)
      END DO 
      CLOSE (UNIT=nfout)
      ! Write GrADS Control File
      OPEN (UNIT=nfctl, FILE=TRIM(DirMain)//DirPreOut//TRIM(VarNameVeg)//nLats//'.ctl', &
            FORM='FORMATTED', ACCESS='SEQUENTIAL', &
            ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE (UNIT=nferr, FMT='(3A,I4)') &
               ' ** (Error) ** Open file ', &
                 TRIM(DirMain)//DirPreOut//TRIM(VarNameVeg)//nLats//'.ctl', &
               ' returned IOStat = ', ios
         STOP  ' ** (Error) **'
      END IF
      WRITE (UNIT=nfctl, FMT='(A)') 'DSET '// &
             TRIM(DirMain)//DirPreOut//TRIM(VarNameVeg)//nLats
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'OPTIONS YREV BIG_ENDIAN'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A,1PG12.5)') 'UNDEF ', UndefG
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'TITLE Vegetation Mask on a Gaussian Grid'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A,I5,A,F8.3,F15.10)') &
                          'XDEF ',Imax,' LINEAR ',0.0_r8,360.0_r8/REAL(Imax,r8)
      WRITE (UNIT=nfctl, FMT='(A,I5,A)') 'YDEF ',Jmax,' LEVELS '
      WRITE (UNIT=nfctl, FMT='(8F10.5)') gLats(Jmax:1:-1)
      WRITE (UNIT=nfctl, FMT='(A)') 'ZDEF  6 LEVELS 0 1 2 3 4 5'
      WRITE (UNIT=nfctl, FMT='(A)') 'TDEF  1 LINEAR JAN2005 1MO'
      WRITE (UNIT=nfctl, FMT='(A)') 'VARS  4'
      WRITE (UNIT=nfctl, FMT='(A)') 'LSMO  0 99 Land Sea Mask Before Fix  [0-Sea 1-Land]'
      WRITE (UNIT=nfctl, FMT='(A)') 'LSMK  0 99 Land Sea Mask For Model   [0-Sea 1-Land]'
      WRITE (UNIT=nfctl, FMT='(A)') 'VGMO  6 99 Vegetation Mask Before Fix [0 to 13]'
      WRITE (UNIT=nfctl, FMT='(A)') 'VEGM  6 99 Vegetation Mask For Model  [0 to 13]'
      WRITE (UNIT=nfctl, FMT='(A)') 'ENDVARS'
      CLOSE (UNIT=nfctl)
   END IF

PRINT *, "*** PorceClayMaskIBIS ENDS NORMALLY ***"

CONTAINS


SUBROUTINE VegetationMaskCheck ()

   IMPLICIT NONE
   DO k=1,Layer
      im=0
      ip=0
      jm=0
      jp=0
      k1=0
      k2=0
      k3=0
      k4=0
      k5=0
      k6=0
      k7=0
      k8=0
      kk=0
   10 CONTINUE
      DO j=1,Jmax
         DO i=1,Imax
            IF (LandSeaMask(i,j) /= 0 .AND. PorceClayMaskIBISOutput(i,j,k) == Undef) THEN
               IF (im == 0) THEN
                  im=Imax
               ELSE
                  im=i-1
               END IF
               IF (ip > Imax) THEN
                  ip=1
               ELSE
                  ip=i+1
               END IF
               IF (jm == 0) THEN
                  jm=1
               ELSE
               jm=j-1
               END IF
               IF (jp > Jmax) THEN
                  jp=Jmax
               ELSE
                  jp=j+1
               END IF
               im=max(min(Imax,im),1)
               jm=max(min(Jmax,jm),1)
               IF (LandSeaMask(im,jm) /= 0) THEN
                  k1=PorceClayMaskIBISOutput(im,jm,k)
               ELSE
                  k1=0
               END IF
               IF (LandSeaMask(i,jm) /= 0) THEN
                  k2=PorceClayMaskIBISOutput(i,jm,k)
               ELSE
                  k2=0
               END IF
               ip=max(min(Imax,ip),1)
               IF (LandSeaMask(ip,jm) /= 0) THEN
                  k3=PorceClayMaskIBISOutput(ip,jm,k)
               ELSE
                  k3=0
               END IF
               IF (LandSeaMask(im,j) /= 0) THEN
                  k4=PorceClayMaskIBISOutput(im,j,k)
               ELSE
                  k4=0
               END IF
               IF (LandSeaMask(ip,j) /= 0) THEN
                  k5=PorceClayMaskIBISOutput(ip,j,k)
               ELSE
                  k5=0
               END IF
               jp=max(min(Jmax,jp),1)
               IF (LandSeaMask(im,jp) /= 0) THEN
                  k6=PorceClayMaskIBISOutput(im,jp,k)
               ELSE
                  k6=0
               END IF
               IF (LandSeaMask(i,jp) /= 0) THEN
                  k7=PorceClayMaskIBISOutput(i,jp,k)
               ELSE
                  k7=0
               END IF
               IF (LandSeaMask(ip,jp) /= 0) THEN
                  k8=PorceClayMaskIBISOutput(ip,jp,k)
               ELSE
                  k8=0
               END IF

               IF (k1+k2+k3+k4+k5+k6+k7+k8 == 0) THEN
                  CALL VegetationMaskFix ()
                  GO TO 10
               END IF

               IF (k1 /= 0) THEN
                  kk=k1
               ELSE
                  kk=-1
               END IF
               IF (k2 /= 0 .AND. kk == -1) THEN
                  kk=k2
               ELSE IF (k2 /= 0 .AND. kk /= k2) THEN
                  CALL VegetationMaskFix ()
                  GO TO 10
               END IF

               IF (k3 /= 0 .AND. kk == -1) THEN
                  kk=k3
               ELSE IF (k3 /= 0 .AND. kk /= k3) THEN
                  CALL VegetationMaskFix ()
                 GO TO 10
               END IF

               IF (k4 /= 0 .AND. kk == -1) THEN
                  kk=k4
               ELSE IF (k4 /= 0 .AND. kk /= k4) THEN
                  CALL VegetationMaskFix ()
                  GO TO 10
               END IF

               IF (k5 /= 0 .AND. kk == -1) THEN
                  kk=k5
               ELSE IF (k5 /= 0 .AND. kk /= k5) THEN
                  CALL VegetationMaskFix ()
                  GO TO 10
               END IF

               IF (k6 /= 0 .AND. kk == -1) THEN
                  kk=k6
               ELSE IF (k6 /= 0 .AND. kk /= k6) THEN
                  CALL VegetationMaskFix ()
                  GO TO 10
               END IF

               IF (k7 /= 0 .AND. kk == -1) THEN
                  kk=k7
               ELSE IF (k7 /= 0 .AND. kk /= k7) THEN
                  CALL VegetationMaskFix ()
                  GO TO 10
               END IF

               IF (k8 /= 0 .AND. kk == -1) THEN
                  kk=k8
               ELSE IF (k8 /= 0 .AND. kk /= k8) THEN
                  CALL VegetationMaskFix ()
                  GO TO 10
               END IF

               PorceClayMaskIBISOutput(i,j,k)=kk
               IF(dumpLocal)WRITE (UNIT=*, FMT='(A,2I4,A,I2)') &
                     ' Undefined Location i, j = ', i, j, &
                     ' Filled with Nearby Value = ', kk
             END IF

             IF (LandSeaMask(i,j) == 0) PorceClayMaskIBISOutput(i,j,k)=Undef

          END DO
       END DO

    END DO
END SUBROUTINE VegetationMaskCheck


SUBROUTINE VegetationMaskFix ()

   IMPLICIT NONE
 
   WRITE (UNIT=*, FMT='(A,3I8,/,A,/,(11X,3I8))') &
         ' At Land Point Value Undefined at i, j = ', i, j,k, &
         ' With no Unambigous Nearby Values.  Local Area:', &
           k1, k2, k3, k4, PorceClayMaskIBISOutput(i,j,k), k5, k6, k7, k8
   WRITE (UNIT=*, FMT='(9(A31,I8))') &
           'k1                            =',k1, &
           'k2                            =',k2, &
           'k3                            =',k3, &
           'k4                            =',k4, &
           'PorceClayMaskIBISOutput(i,j,k)=',PorceClayMaskIBISOutput(i,j,k),  &
           'k5                            =',k5, &
           'k6                            =',k6, &
           'k7                            =',k7, &
           'k8                            =',k8


   IF (LandSeaMask(i,j) == 1) THEN
     LandSeaMask(i,j)=0
     IF(dumpLocal)WRITE (UNIT=*, FMT='(A)') ' Land Sea Mask Changed from 1 to 0'
   ELSE
     LandSeaMask(i,j)=1
     IF(dumpLocal)WRITE (UNIT=*, FMT='(A)') ' Land Sea Mask Changed from 0 to 1'
   END IF

END SUBROUTINE VegetationMaskFix


END PROGRAM PorceClayMaskIBIS
