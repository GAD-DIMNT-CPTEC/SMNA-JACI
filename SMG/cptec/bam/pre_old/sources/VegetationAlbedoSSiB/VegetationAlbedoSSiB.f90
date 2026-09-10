!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
PROGRAM VegetationAlbedoSSiB

   IMPLICIT NONE

   INTEGER, PARAMETER :: &
            r4 = SELECTED_REAL_KIND(6) ! Kind for 32-bits Real Numbers

   INTEGER :: ios

   CHARACTER (LEN=528) :: DirMain

   CHARACTER (LEN=13) :: DirModelIn='model/datain/'

   CHARACTER (LEN=12) :: DirBCs='pre/databcs/'

   CHARACTER (LEN=11) :: AlbedoBCs='sibalb.form'

   CHARACTER (LEN=10) :: AlbedoOut='AlbedoSSiB'

   CHARACTER (LEN=11) :: VegetationBCs='sibveg.form'

   CHARACTER (LEN=14) :: VegetationOut='VegetationSSiB'

   CHARACTER (LEN=24) :: NameNML='VegetationAlbedoSSiB.nml'

   INTEGER :: nferr=0    ! Standard Error Print Out
   INTEGER :: nfinp=5    ! Standard Read In
   INTEGER :: nfprt=6    ! Standard Print Out
   INTEGER :: nfvgi=10   ! To Read Formatted Veg-Morpho-Physio-SSiB Data
   INTEGER :: nfvgo=20   ! To Write Unformatted Veg-Morpho-Physio-SSiB Data
   INTEGER :: nfali=30   ! To Read Formatted Albedo SSiB Data
   INTEGER :: nfalo=40   ! To Write Unformatted Albedo SSiB Data

   NAMELIST /InputDim/ DirMain

   DirMain='./ '

   OPEN (UNIT=nfinp, FILE='./'//NameNML, &
         FORM='FORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='READ', STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              './'//NameNML, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   READ  (UNIT=nfinp, NML=InputDim)
   CLOSE (UNIT=nfinp)

   WRITE (UNIT=nfprt, FMT='(/,A)')  ' &InputDim'
   WRITE (UNIT=nfprt, FMT='(A)')    '  DirMain = '//TRIM(DirMain)
   WRITE (UNIT=nfprt, FMT='(A,/)')  ' /'

   CALL VegMorphoPhysioSSiBData ()

   CALL AlbedoSSiBData ()

PRINT *, "*** VegetationAlbedoSSiB ENDS NORMALLY ***"

CONTAINS


SUBROUTINE VegMorphoPhysioSSiBData ()
 
   ! Vegetation Morphological and Physiological SSiB Data

   IMPLICIT NONE

   INTEGER, PARAMETER :: ityp=13, imon=12, icg=2, &
                         iwv=3, ild=2, idp=3, ibd=2

   INTEGER :: i, ntyp, iv, im, k, mmm

   ! Vegetation and Soil Parameters

   REAL (KIND=r4) tran(ityp,icg,iwv,ild), &
                  ref(ityp,icg,iwv,ild), &
                  rstpar(ityp,icg,iwv), &
                  soref(ityp,iwv), &
                  chil(ityp,icg), &
                  topt(ityp,icg), &
                  tll(ityp,icg), &
                  tu(ityp,icg), &
                  defac(ityp,icg), &
                  ph1(ityp,icg), &
                  ph2(ityp,icg), &
                  rootd(ityp,icg), &
                  rlmax(ityp,icg), &
                  rootca(ityp,icg), &
                  rplant(ityp,icg), &
                  rdres(ityp,icg), &
                  bee(ityp), &
                  phsat(ityp), &
                  satco(ityp), &
                  poros(ityp), &
                  slope(ityp), &
                  zdepth(ityp,idp), &
                  green(ityp,imon,icg), &
                  xcover(ityp,imon,icg), &
                  zlt(ityp,imon,icg), &
                  x0x(ityp,imon),&
                  xd(ityp,imon), &
                  z2(ityp,imon), &
                  z1(ityp,imon), &
                  xdc(ityp,imon), &
                  xbc(ityp,imon), &
                  rootl(ityp,imon,icg)

   WRITE(UNIT=nfprt, FMT='(/,A,/)') ' From VegMorphoPhysioSSiBData:'

   OPEN (UNIT=nfvgi, FILE=TRIM(DirMain)//DirBCs//VegetationBCs, &
         FORM='FORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='READ', STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirBCs//VegetationBCs, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   READ (UNIT=nfvgi, FMT='(A)')
   DO
      READ  (UNIT=nfvgi, FMT='(I3)') ntyp
      WRITE (UNIT=nfprt, FMT='(A,I3,A)') ' Type = ', ntyp, ' Read In'
      READ  (UNIT=nfvgi, FMT=*) ((rstpar(ntyp,iv,im),im=1,3),iv=1,2), &
                                (chil(ntyp,iv),iv=1,2), &
                                (topt(ntyp,iv),iv=1,2), &
                                (tll(ntyp,iv),iv=1,2), &
                                (tu(ntyp,iv),iv=1,2), &
                                (defac(ntyp,iv),iv=1,2), &
                                (ph1(ntyp,iv),iv=1,2), &
                                (ph2(ntyp,iv),iv=1,2), &
                                (rootd (ntyp,iv),iv=1,2), &
                                bee(ntyp),phsat(ntyp), &
                                satco(ntyp), poros(ntyp), &
                                (zdepth(ntyp,k),k=1,3)
      DO i=1,12
         READ (UNIT=nfvgi, FMT=*) mmm, (zlt(ntyp,i,iv),iv=1,2), &
                                  (green (ntyp,i,iv),iv=1,2), &
                                  z2(ntyp,i), z1(ntyp,i), &
                                  (xcover(ntyp,i,iv),iv=1,2), &
                                  x0x(ntyp,i), xd(ntyp,i), &
                                  xbc(ntyp,i),xdc(ntyp,i)
      END DO
      IF(ntyp == ityp) EXIT
   END DO
   CLOSE (UNIT=nfvgi)

   OPEN (UNIT=nfvgo, FILE=TRIM(DirMain)//DirModelIn//VegetationOut, &
         FORM='UNFORMATTED', ACCESS='SEQUENTIAL', ACTION='WRITE', &
         STATUS='REPLACE', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirModelIn//VegetationOut, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   WRITE (UNIT=nfvgo) rstpar, chil, topt, tll, tu, defac, ph1, ph2, &
                      rootd, bee, phsat, satco, poros, zdepth
   WRITE (UNIT=nfvgo) green, xcover, zlt, x0x, xd, z2, z1, xdc, xbc
   CLOSE (UNIT=nfvgo)

   DO ntyp=1,ityp
      WRITE (UNIT=nfprt, FMT='(I3)') ntyp
      WRITE (UNIT=nfprt, FMT='(6F12.4,/,2F8.2,/,6F8.1,/,2F10.4,4F8.2,/,&
                              &2F6.1,/,2F10.4,/,E12.5,F10.4,/,3E12.5)') &
                         ((rstpar(ntyp,iv,im),im=1,3),iv=1,2), &
                         (chil(ntyp,iv),iv=1,2), &
                         (topt(ntyp,iv),iv=1,2), &
                         (tll(ntyp,iv),iv=1,2), &
                         (tu(ntyp,iv),iv=1,2), &
                         (defac(ntyp,iv),iv=1,2), &
                         (ph1(ntyp,iv),iv=1,2), &
                         (ph2(ntyp,iv),iv=1,2), &
                         (rootd(ntyp,iv),iv=1,2), &
                         bee(ntyp), phsat(ntyp), &
                         satco(ntyp), poros(ntyp), &
                         (zdepth(ntyp,k),k=1,3)
      DO i=1,12
         WRITE (UNIT=nfprt, FMT='(I3,/,2F10.5,2F10.5,2F8.4,/,&
                                 &2F10.5,2F10.5,/,2F10.2)') &
                            i, (zlt(ntyp,i,iv),iv=1,2), &
                            (green (ntyp,i,iv),iv=1,2), &
                            z2(ntyp,i), z1(ntyp,i), &
                            (xcover(ntyp,i,iv),iv=1,2), &
                            x0x(ntyp,i), xd(ntyp,i), &
                            xbc(ntyp,i), xdc(ntyp,i)
      END DO
   END DO

   WRITE(UNIT=nfprt, FMT='(/)')

END SUBROUTINE VegMorphoPhysioSSiBData


SUBROUTINE AlbedoSSiBData ()

   IMPLICIT NONE

   INTEGER, PARAMETER :: ityp=13, imon=12, icg=2, njj=6, &
                         ild=2, ibd=2, nj=9, nk=3

   ! Albedo Parameters

   REAL (KIND=r4) :: cedfu(ityp,imon,nj), &
                     cedir(ityp,imon,nj,3), &
                     cedfu1(2,ityp,imon,njj,3), &
                     cedir1(2,ityp,imon,njj,nk,3), &
                     cedfu2(2,ityp,imon,njj,3), &
                     cedir2(2,ityp,imon,njj,nk,3), &
                     cledfu(ityp,imon,nj), &
                     cledir(ityp,imon,nj,3), &
                     cether(ityp,imon,2), &
                     xmiu(imon,nk), &
                     xmiw(imon,nk)

   INTEGER :: itp, l, j, i, k, ml, mon

   WRITE(UNIT=nfprt, FMT='(/,A,/)') ' From AlbedoSSiBData:'

   OPEN (UNIT=nfali, FILE=TRIM(DirMain)//DirBCs//AlbedoBCs, &
         FORM='FORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='READ', STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirBCs//AlbedoBCs, &
            ' returned IOStat = ', ios 
      STOP  ' ** (Error) **'
   END IF
   DO itp=1,12
      DO l=1,2
         READ (UNIT=nfali, FMT='(6E12.5)') (cether(itp,mon,l),mon=1,12)
      END DO
      DO j=1,nj
         READ (UNIT=nfali, FMT='(6E12.5)') (cledfu(itp,mon,j),mon=1,12)
         DO l=1,3
            READ (UNIT=nfali, FMT='(6E12.5)') (cledir(itp,mon,j,l),mon=1,12)
         END DO
      END DO
      DO i=1,nk
         READ (UNIT=nfali, FMT='(6E12.5)') (xmiu(mon,i),mon=1,12)
      END DO
      DO j=1,nj
         READ (UNIT=nfali, FMT='(6E12.5)') (cedfu(itp,mon,j),mon=1,12)
         DO l=1,3
            READ (UNIT=nfali, FMT='(6E12.5)') (cedir(itp,mon,j,l),mon=1,12)
         END DO
      END DO
      DO ml=1,2
         DO j=1,njj
            DO l=1,3
               READ (UNIT=nfali, FMT='(6E12.5)') (cedfu1(ml,itp,mon,j,l),mon=1,12)
               READ (UNIT=nfali, FMT='(6E12.5)') (cedfu2(ml,itp,mon,j,l),mon=1,12)
            END DO
         END DO
         DO j=1,njj
            DO k=1,nk
               DO l=1,3
                  READ (UNIT=nfali, FMT='(6E12.5)') (cedir1(ml,itp,mon,j,k,l),mon=1,12)
                  READ (UNIT=nfali, FMT='(6E12.5)') (cedir2(ml,itp,mon,j,k,l),mon=1,12)
               END DO
            END DO
         END DO
      END DO
   END DO
   DO l=1,2
      READ (UNIT=nfali, FMT='(7E12.5)') (cether(13,mon,l),mon=1,7)
   END DO
   DO j=1,nj
      READ (UNIT=nfali, FMT='(7E12.5)') (cledfu(13,mon,j),mon=1,7)
      DO l=1,3
         READ (UNIT=nfali, FMT='(7E12.5)') (cledir(13,mon,j,l),mon=1,7)
      END DO
   END DO
   DO i=1,nk
      READ (UNIT=nfali, FMT='(7E12.5)') (xmiw(mon,i),mon=1,7)
   END DO
   DO i=1,nk
      WRITE(UNIT=nfprt, FMT='(7E12.5)') (xmiw(mon,i),mon=1,7)
   END DO
   DO j=1,nj
      READ (UNIT=nfali, FMT='(7E12.5)') (cedfu(13,mon,j),mon=1,7)
      DO l=1,3
         READ (UNIT=nfali, FMT='(7E12.5)') (cedir(13,mon,j,l),mon=1,7)
      END DO
   END DO
   DO ml=1,2
      DO j=1,njj
         DO l=1,3
            READ (UNIT=nfali, FMT='(7E12.5)') (cedfu1(ml,13,mon,j,l),mon=1,7)
            READ (UNIT=nfali, FMT='(7E12.5)') (cedfu2(ml,13,mon,j,l),mon=1,7)
         END DO
      END DO
      DO j=1,njj
         DO k=1,nk
            DO l=1,3
               READ (UNIT=nfali, FMT='(7E12.5)') (cedir1(ml,13,mon,j,k,l),mon=1,7)
               READ (UNIT=nfali, FMT='(7E12.5)') (cedir2(ml,13,mon,j,k,l),mon=1,7)
            END DO
         END DO
      END DO
   END DO
   CLOSE (UNIT=nfali)

   OPEN (UNIT=nfalo, FILE=TRIM(DirMain)//DirModelIn//AlbedoOut, &
         FORM='UNFORMATTED', ACCESS='SEQUENTIAL', ACTION='WRITE', &
         STATUS='REPLACE', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirModelIn//AlbedoOut, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   WRITE (UNIT=nfalo) cedfu, cedir, cedfu1, cedir1, cedfu2, cedir2, &
                      cledfu, cledir, xmiu, cether, xmiw
   CLOSE (UNIT=nfalo)

   WRITE(UNIT=nfprt, FMT='(/)')

END SUBROUTINE AlbedoSSiBData


END PROGRAM VegetationAlbedoSSiB
