!
!  $Author: pkubota $
!  $Date: 2007/10/18 18:39:59 $
!  $Revision: 1.7 $
!
MODULE FileAccess

  USE Constants, ONLY : r4, i4, r8, nferr, nfprt, nffct, nfpos, &
                        nfctl, nfmdf, nfppf, nfdir, nfrfd,rfd,&
                        RunRecort,RecLat,RecLon,newlat0,newlat1,newlon0,newlon1
  USE Utils, ONLY : getpoint
  USE Sizes, ONLY : imax, jmax, ibmax, jbmax, ibMaxPerJB, iPerIJB, jPerIJB, &
                    mnMax, mymnMax, mymmax, msinproc, mnmap, mymnmap, mmax, &
                    myfirstlev, mylastlev

  USE RegInterp, ONLY : gLats,glond, glatd

  USE Parallelism, ONLY : myid, myid_four

  IMPLICIT NONE

  PRIVATE

  CHARACTER (LEN=128) :: dirctl=' ' ! dir. archived outputs 
                                    ! (used only on a write statement)
  CHARACTER (LEN=7  ) :: namep =' ' ! post processed file prefix
  CHARACTER (LEN=5  ) :: extp  =' ' ! post processed file extension
  CHARACTER (LEN=5  ) :: exdp  =' ' ! post processed directives file ext.
  CHARACTER (LEN=23 ) :: modout=' ' ! model output subdirectory
  CHARACTER (LEN=33 ) :: moddir=' ' ! preffix name of model output arquive with file names
  CHARACTER (LEN=400) :: modfls=' ' ! total name of model output arquive with file names

  INTEGER :: id = 1 ! LEN_TRIM(dirctl) (used only on a write statement)
  INTEGER :: im = 1 ! LEN_TRIM(modfls)


  PUBLIC :: InitFiles, opnpos, CloseFiles, ReadHeader, ReadField, &
            ReadFieldSp, ReadFieldG, skipf, WriteField

CONTAINS

SUBROUTINE InitFiles ()

  USE Constants, ONLY : nFBeg, nFFrs, nFEnd, trunc, lev, labeli, labelf, &
                        prefx,prefy, req, datain,datalib,dataout

  IMPLICIT NONE

  INTEGER :: ios
  INTEGER :: nr
  INTEGER :: n

  CHARACTER (LEN=1)  :: skip  
  CHARACTER (LEN=33) :: posctl
!  CHARACTER (LEN=12) :: posinp

  LOGICAL :: lex

  IF(RunRecort)THEN
     IF ( req == "e" ) THEN
       dirctl='/forecast/ETA'
       namep='GPOSETA'
     ELSE
       dirctl='/forecast/'//TRIM(prefx)
       namep='G'//TRIM(prefy)//TRIM(prefx)
     ENDIF
  ELSE
     IF ( req == "e" ) THEN
        dirctl='/forecast/ETA'
        namep='GPOSETA'
     ELSE
        dirctl='/forecast/'//TRIM(prefx)
        namep='GPOS'//TRIM(prefx)
     END IF
  END IF
  IF (req == 'p') THEN
    rfd='rfd.pnt'
    extp='P.unf'
    exdp='P.ctl'
  ELSE IF (req == 's') THEN
    rfd='rfd.sfc'
    extp='S.unf'
    exdp='S.ctl'
  ELSE IF (req == 'c') THEN
    rfd='rfd.clm'
    extp='C.unf'
    exdp='C.ctl'
  ELSE IF (req == 'e') THEN
    rfd='rfd.eta'
    extp='E.unf'
    exdp='E.ctl'
  ELSE IF (req == 'g') THEN
       rfd='rfd.ens'
       extp='P.unf'
       exdp='P.ctl'
  ELSE
    rfd='rfd    '
    extp='D.unf'
    exdp='D.ctl'
  END IF

  id=MAX(1_i4, LEN_TRIM(dirctl))

  moddir='GFCT'//prefx//labeli//labelf//'F.dir.'
  modfls=TRIM(datain)//'/'//moddir//TRIM(trunc)//TRIM(lev)//'.files'
  im=MAX(1_i4, LEN_TRIM(modfls))

  INQUIRE(file=modfls(1:im),EXIST=lex)
  IF (.NOT.lex) THEN
     WRITE(*,'(A)') 'ERROR:: InitFiles(): Listing of forecast files does not exist!'
     WRITE(*,'(A)') 'FILE='//modfls(1:im)
     STOP
  ENDIF
  OPEN (UNIT=nfmdf, FILE=TRIM(modfls(1:im)), FORM='FORMATTED', ACCESS='SEQUENTIAL', &
        ACTION='READ', STATUS='OLD', IOSTAT=ios)

  IF (ios /= 0) THEN
     WRITE (UNIT=nferr, FMT='(3A,I4)') ' ** (Error) ** Open file ', &
                                       TRIM(modfls), &
                                       ' returned IOStat = ', ios 
     STOP ' ** (Error) **'
  END IF


  ! Unit nfmdf will be permanently connected to a file that contains
  ! a list of pairs of forecast output files, sorted
  ! in ascending forecasted time. Some of those
  ! will be input to the post-processing.

  ! Skip pairs of files that won't be used,
  ! corresponding to undesired forecasted times.

   nr=2*(nFBeg-nFFrs)
   IF (nr > 0) THEN
      DO n=1,nr
         READ (UNIT=nfmdf, FMT='(A)') skip
      END DO
   END IF

   IF (nFBeg <= 0 .AND. nFEnd <= 0) THEN
      posctl=namep//labeli//labelf//exdp(1:2)//'icn.'
   ELSE
      posctl=namep//labeli//labelf//exdp(1:2)//'fct.'
   END IF

  ! Unit nfppf will be permanently connected to an output file
  ! that will contain the names of post-processing output
  ! files in the archiving machine.

   OPEN (UNIT=nfppf, FILE=TRIM(dataout)//'/'//TRIM(posctl)//TRIM(trunc)//TRIM(lev)//'.lst', &
         FORM='FORMATTED', ACCESS='SEQUENTIAL', ACTION='WRITE', &
         STATUS='REPLACE', IOSTAT=ios)

   IF (ios /= 0) THEN
      IF(myid.eq.0) WRITE (UNIT=nferr, FMT='(3A,I4)') ' ** (Error) ** Open file ', &
                                        TRIM(dataout)//'/'//TRIM(posctl)//TRIM(trunc)//TRIM(lev)//'.lst', &
                                        ' returned IOStat = ', ios
      STOP ' ** (Error) **'
   END IF

   ! Unit nfrfd -> Requested Post-processed Fields

   OPEN (UNIT=nfrfd, FILE=TRIM(datalib)//'/'//rfd, FORM='FORMATTED', &
         ACCESS='SEQUENTIAL', ACTION='READ', STATUS='OLD', IOSTAT=ios)

   IF (ios /= 0) THEN
      IF(myid.eq.0) WRITE (UNIT=nferr, FMT='(3A,I4)') ' ** (Error) ** Open file ', &
                                        TRIM(datalib)//'/'//rfd, &
                                        ' returned IOStat = ', ios
      STOP ' ** (Error) **'
   END IF


END SUBROUTINE InitFiles


SUBROUTINE opnpos (labelp)

  USE Constants, ONLY : trunc, lev, labeli, Binary,datain,datalib,dataout

  IMPLICIT NONE

  ! Opnpos: Compose filenames and connect files to units:
  ! Two files are connected once and used throughout execution:
  !     Unit nfmdf, to read forecast output filenames;
  !     Unit nfppf, to write archieved post-processing output filenames.
  ! Four file connections will change at each invocation, according
  ! to the forecasted output time beeing post-processed:
  !     Unit nfdir, to read forecast output file structure;
  !     Unit nffct, to read forecast output data;
  !     Unit nfpos, to write post-processed data;
  !     Unit nfctl, to write post-processed geometry;

  CHARACTER (LEN=10), INTENT(OUT) :: labelp

  INTEGER :: ios
  INTEGER :: in
  CHARACTER (LEN=256) :: namfct1
  CHARACTER (LEN=256) :: namfct2
  INTEGER (kind=i4) :: ierr

  READ (UNIT=nfmdf, FMT='(A)', IOSTAT=ierr) namfct1
  IF(ierr /= 0) THEN
      IF(myid.eq.0) WRITE (nferr,*) TRIM(modfls)
      IF(myid.eq.0) WRITE (nferr,*)'<- THE FILE LIST WAS FINALIZED ->'
      STOP
  END IF

  ! Unit nfdir will be connected to a file that contains 
  ! the structure of the forecast output file for the current time
  ! Opened for reading.


  OPEN (UNIT=nfdir, FILE=TRIM(namfct1), FORM='FORMATTED', &
        ACCESS='SEQUENTIAL', ACTION='READ', STATUS='OLD', IOSTAT=ios)

  IF (ios /= 0) THEN
     IF(myid.eq.0) WRITE (UNIT=nferr, FMT='(3A,I4)') ' ** (Error) ** Open file ', &
                                       TRIM(namfct1), &
                                       ' returned IOStat = ', ios
      STOP ' ** (Error) **'
  END IF



  READ (UNIT=nfmdf, FMT='(A)') namfct2

  ! Unit nffct will be connected to the forecast output file for 
  ! the current time
  ! Opened for reading.


  OPEN (UNIT=nffct, FILE=TRIM(namfct2), FORM='UNFORMATTED', &
        ACCESS='SEQUENTIAL', ACTION='READ', STATUS='OLD', IOSTAT=ios)

  IF (ios /= 0) THEN
     WRITE (UNIT=nferr, FMT='(3A,I4)') ' ** (Error) ** Open file ', &
                                       TRIM(namfct2), &
                                       ' returned IOStat = ', ios
      STOP ' ** (Error) **'
  END IF

  ! Compose parts of output filename

  in=MAX(1,LEN_TRIM(namfct1))
  IF (namfct1(in-12:in-11) /= 'ir') exdp(4:5)=namfct1(in-12:in-11)
  extp(3:5)=namfct2(in-13:in-11)
  labelp=namfct2(in-25:in-16)
  
  IF (Binary) then
  
  IF(myid.eq.0) WRITE (UNIT=nfprt, FMT='(A)') TRIM(dataout)//'/'//namep// &
                                labeli//labelp//extp//'.'//TRIM(trunc)//TRIM(lev)

  ! Unit nfpos will be opened for writing orography and
  ! land sea mask for the current time.

  OPEN (UNIT=nfpos, FILE=TRIM(dataout)//'/'//namep// &
                      labeli//labelp//extp//'.'//TRIM(trunc)//TRIM(lev), &
        FORM='UNFORMATTED', ACCESS='SEQUENTIAL', ACTION='WRITE', &
        STATUS='REPLACE', IOSTAT=ios)

  IF (ios /= 0) THEN
     WRITE (UNIT=nferr, FMT='(3A,I4)') ' ** (Error) ** Open file ', &
                                       TRIM(dataout)//'/'//namep// &
                                       labeli//labelp//extp//'.'//TRIM(trunc)//TRIM(lev), &
                                       ' returned IOStat = ', ios
      STOP ' ** (Error) **'
  END IF

  IF(myid.eq.0) WRITE (UNIT=nfprt, FMT='(A)') TRIM(dataout)//'/'//namep//labeli// &
                                labelp//extp//'.'//TRIM(trunc)//TRIM(lev)//'.ctl'

  ! Unit nfctl will be opened for writing control
  ! values for the current time.

  OPEN (UNIT=nfctl, FILE=TRIM(dataout)//'/'//namep//labeli// &
                 labelp//extp//'.'//TRIM(trunc)//TRIM(lev)//'.ctl', &
        FORM='FORMATTED', ACCESS='SEQUENTIAL', ACTION='WRITE', &
        STATUS='REPLACE', IOSTAT=ios)
  IF (ios /= 0) THEN
     WRITE (UNIT=nferr, FMT='(3A,I4)') ' ** (Error) ** Open file ', &
                                       TRIM(dataout)//'/'//namep//labeli// &
                                       labelp//extp//'.'//TRIM(trunc)//TRIM(lev)//'.ctl', &
                                       ' returned IOStat = ', ios
      STOP ' ** (Error) **'
  END IF

  ! Dump names of post-processing output
  ! files in the archiving machine on unit nfppf

  IF(myid.eq.0) WRITE (UNIT=nfppf, FMT='(A)') TRIM(dirctl)//'/'//namep//labeli// &
                                labelp//extp//'.'//TRIM(trunc)//TRIM(lev)//'.ctl'
  IF(myid.eq.0) WRITE (UNIT=nfppf, FMT='(A)') TRIM(dirctl)//'/'//namep//labeli// &
                                labelp//extp//'.'//TRIM(trunc)//TRIM(lev)

  ELSE
  ! Unit 12 will be opened for writing control
  ! values for the current time.

  IF(myid.eq.0) WRITE (UNIT=nferr, FMT='(A)') TRIM(dataout)//'/'//namep//labeli// &
                            labelp//extp//'.'//TRIM(trunc)//TRIM(lev)//'.ctl'
  OPEN (UNIT=nfctl, FILE=TRIM(dataout)//'/'//namep//labeli// &
        labelp//extp//'.'//TRIM(trunc)//TRIM(lev)//'.ctl', STATUS='REPLACE')

  ! Dump names of post-processing output
  ! files in the archiving machine on unit 14.

  IF(myid.eq.0) WRITE (UNIT=nfppf, FMT='(A)') TRIM(dirctl)//'/'//namep//labeli// &
                                labelp//extp//'.'//TRIM(trunc)//TRIM(lev)//'.ctl'
  IF(myid.eq.0) WRITE (UNIT=nfppf, FMT='(A)') TRIM(dirctl)//'/'//namep//labeli// &
                                labelp//extp//'.'//TRIM(trunc)//TRIM(lev)//'.grb'

!
! open Grib file
!
  IF(myid.eq.0) WRITE (UNIT=nferr, FMT='(A)') TRIM(dataout)//'/'//namep// &
                            labeli//labelp//extp//'.'//TRIM(trunc)//TRIM(lev)//'.grb'

  CALL BAOPEN(51,TRIM(dataout)//'/'//namep//labeli// &
                 labelp//extp//'.'//TRIM(trunc)//TRIM(lev)//'.grb',ierr)

  if(ierr.ne.0)then
    stop ' Error in BAOPEN'
  endif

  END IF






END SUBROUTINE opnpos

SUBROUTINE ReadField (mdim, ldim, bufa)

  IMPLICIT NONE

  INTEGER, INTENT(IN ) :: mdim
  INTEGER, INTENT(IN ) :: ldim
  REAL (KIND=r8),    INTENT(OUT) :: bufa(mdim,ldim)

  INTEGER :: l

  REAL (KIND=r4) :: bufb(mdim)

  DO l=1,ldim
     READ (UNIT=nffct) bufb
     bufa(:,l)=REAL(bufb,r8)
  END DO

END SUBROUTINE ReadField


SUBROUTINE ReadFieldG (mdim, ldim, bufa)

  IMPLICIT NONE

  INTEGER, INTENT(IN ) :: mdim
  INTEGER, INTENT(IN ) :: ldim
  REAL (KIND=r8),    INTENT(OUT) :: bufa(ibmax,ldim,jbmax)

  INTEGER :: i,j,l,ib,jb

  REAL (KIND=r4) :: bufb(mdim)

  IF (mdim.ne.imax*jmax) THEN
     IF (myid.eq.0) THEN
        write(*,*) ' Inconsistence in reading files in ReadFieldG '
        write(*,*) ' mdim ',mdim, 'imax*jmax ',imax*jmax
     ENDIF
     STOP
  ENDIF

  DO l=1,ldim
     READ (UNIT=nffct) bufb(:)
     DO jb = 1, jbMax
        DO ib = 1,ibMaxPerJB(jb)
           i = iPerIJB(ib,jb)
           j = jPerIJB(ib,jb)
           bufa(ib,l,jb)=REAL(bufb(i+(j-1)*imax),r8)
        END DO
     END DO
  END DO


END SUBROUTINE ReadFieldG

SUBROUTINE ReadFieldSp (mdim,klev,bufa)


  IMPLICIT NONE

  INTEGER,           INTENT(IN)  :: klev
  INTEGER,           INTENT(IN)  :: mdim
  REAL (KIND=r8),    INTENT(OUT) :: bufa(2*mymnMax,klev)

  INTEGER :: mm, nn, i1, i2, m, k

  REAL (KIND=r4) :: bufb(2*mnMax,klev)

  IF (mdim.ne.2*mnmax) then
     IF (myid.eq.0) THEN
        WRITE (nferr,*) ' trying to read spectral field of',&
                   ' wrong dimension ',mdim, ' mnmax = ',mnmax
     END IF
     STOP
  ENDIF

  DO k=1,klev
     READ (UNIT=nffct) bufb(:,k)
  END DO
  IF (klev.eq.1) THEN
       DO mm=1,mymmax
          m = msinproc(mm,myid_four)
          i1 = 2*mnmap(m,m)-1
          i2 = 2*mymnmap(mm,m)-1
          DO nn=0,2*(mmax-m)+1
             bufa(i2+nn,1) = REAL(bufb(i1+nn,1),r8)
          ENDDO
       ENDDO
     ELSE
       DO k=myfirstlev,mylastlev
          DO mm=1,mymmax
             m = msinproc(mm,myid_four)
             i1 = 2*mnmap(m,m)-1
             i2 = 2*mymnmap(mm,m)-1
             DO nn=0,2*(mmax-m)+1
                bufa(i2+nn,k+1-myfirstlev) = REAL(bufb(i1+nn,k),r8)
             ENDDO
          ENDDO
       ENDDO
   ENDIF

END SUBROUTINE ReadFieldSp

SUBROUTINE ReadHeader (ifday, tod, idate, idatec)

  IMPLICIT NONE

  INTEGER, INTENT (OUT) :: ifday
  REAL (KIND=r8),INTENT (OUT) :: tod
  INTEGER, INTENT (OUT) :: idate(4)
  INTEGER, INTENT (OUT) :: idatec(4)

  INTEGER (KIND=i4) :: ifday4
  REAL    (KIND=r4) :: tod4
  INTEGER (KIND=i4) :: idate4(4)
  INTEGER (KIND=i4) :: idatec4(4)

  READ (UNIT=nffct) ifday4, tod4, idate4, idatec4

  ifday=INT(ifday4)
  tod=REAL(tod4,r8)
  idate=INT(idate4)
  idatec=INT(idatec4)

END SUBROUTINE ReadHeader

SUBROUTINE WriteField (ndim, bfr)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ndim
  REAL (KIND=r8), INTENT(IN) :: bfr(:,:)
  INTEGER        :: ijk,j,i
  REAL (KIND=r8) :: bfrg(SIZE(bfr,1),SIZE(bfr,2))
  REAL (KIND=r4) :: bfr4(ndim)

    bfr4=0.0_r4
    bfrg=0.0_r8
    IF(RunRecort)THEN
       IF(RecLon(1) < 0.0_r8 .or. RecLon(2) < 0.0_r8)THEN
          bfrg=CSHIFT (bfr, SHIFT=SIZE(bfr,1)/2,DIM = 1) 
       ELSE
          bfrg=bfr
       END IF
       ijk=0
       DO j=newlat0,newlat1
          DO i=newlon0,newlon1
             ijk=ijk+1
             bfr4(ijk)=REAL(bfrg(i,j),r4)
          END DO
       END DO
    ELSE
       ijk=0
       DO j=1,SIZE(bfr,2)
          DO i=1,SIZE(bfr,1)
             ijk=ijk+1
             bfr4(ijk)=REAL(bfr(i,j),r4)
          END DO
       END DO
    END IF

    WRITE (UNIT=nfpos) bfr4(1:ijk)

END SUBROUTINE WriteField



SUBROUTINE skipf (lev)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: lev

  INTEGER :: l

  DO l=1,lev
     READ (UNIT=nffct)
  END DO

END SUBROUTINE skipf


SUBROUTINE CloseFiles ()

   IMPLICIT NONE

   CLOSE (UNIT=nffct)
   CLOSE (UNIT=nfpos)
   CLOSE (UNIT=nfctl)
   CLOSE (UNIT=nfdir)
   CLOSE (UNIT=51)
END SUBROUTINE CloseFiles


END MODULE FileAccess
