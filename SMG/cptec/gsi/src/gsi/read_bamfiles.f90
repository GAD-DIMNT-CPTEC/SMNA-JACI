subroutine read_bamfiles(mype)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    read_files       get info about atm & sfc guess files
!   prgmmr: derber           org: np23                date: 2002-11-14
!
! abstract:  This routine determines how many global atmospheric and
!            surface guess files are present.  The valid time for each
!            guess file is determine.  The time are then sorted in
!            ascending order.  This information is broadcast to all
!            mpi tasks.
!
! program history log:
!   2002-11-14  derber
!   2004-06-16  treadon - update documentation
!   2004-08-02  treadon - add only to module use, add intent in/out
!   2004-12-02  treadon - replace mpe_ibcast (IBM extension) with
!                         standard mpi_bcast
!   2005-01-27  treadon - make use of sfcio module
!   2005-02-18  todling - no need to read entire sfc file; only head needed
!   2005-03-30  treadon - clean up formatting of write statements
!   2006-01-09  treadon - use sigio to read gfs spectral coefficient file header
!   2007-03-01  tremolet - measure time from beginning of assimilation window
!   2007-04-17  todling  - getting nhr_assimilation from gsi_4dvar
!   2008-05-27  safford - rm unused vars
!   2009-01-07  todling - considerable revamp (no pre-assigned dims)
!   2009-08-26  li      - add variables to handle nst guess files
!   2010-04-20  jing    - set hrdifsig_all and hrdifsfc_all for non-ESMF cases.
!   2010-12-06  Huang   - make use of nemsio_module to check whether atm and sfc files
!                         are in NEMSIO format and get header informaion 'lpl'
!   2010-12-16  li      - (1) set nfldnst, ntguesnst, ifilenst, hrdifnst, hrdifnst_all for nst files. 
!                         (2) add zero initialization of nfldsfc
!   2011-04-04  Huang   - change looping over 0,99 to find existed sigf and sfcf
!                         twice. Use fcst_hr_sig and fcst_hr_sfc to store info of
!                         files found in first loop of 0,99. Use nfldsig and nfldsfc
!                         to access needed sigf and sfcf w/ fcst_hr_sig and *_sfc.
!
!   input argument list:
!     mype     - mpi task id
!
!   output argument list:
!
!   comments:
!     The difference of time Info between operational GFS IO (gfshead%, sfc_head%),
!      analysis time (iadate), and NEMSIO (idate=)
!
!       gfshead & sfc_head            NEMSIO Header           Analysis time (obsmod)
!       ===================   ============================  ==========================
!         %idate(1)  Hour     idate(1)  Year                iadate(1)  Year
!         %idate(2)  Month    idate(2)  Month               iadate(2)  Month
!         %idate(3)  Day      idate(3)  Day                 iadate(3)  Day
!         %idate(4)  Year     idate(4)  Hour                iadate(4)  Hour
!                             idate(5)  Minute              iadate(5)  Minute
!                             idate(6)  Scaled seconds
!                             idate(7)  Seconds multiplier
!
!     The difference of header forecasting hour Info bewteen operational GFS IO
!      (gfshead%, sfc_head%) and NEMSIO
!
!           gfshead & sfc_head                NEMSIO Header
!       ==========================     ============================
!       %fhour  FCST Hour (r_kind)     nfhour     FCST Hour (i_kind)
!                                      nfminute   FCST Mins (i_kind)
!                                      nfsecondn  FCST Secs (i_kind) numerator
!                                      nfsecondd  FCST Secs (i_kind) denominator
!
!       %fhour = float(nfhour) + float(nfminute)/r60 + float(nfsecondn)/float(nfsecondd)/r3600

! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_kind,r_single,i_kind,i_llong
  use mpimod, only: mpi_rtype,mpi_comm_world,ierror,npe,mpi_itype

  use guess_grids, only: nfldsig,nfldsfc,nfldnst,ntguessig,ntguessfc,ntguesnst,&
                         ifilesig,ifilesfc,ifilenst,hrdifsig,hrdifsfc,hrdifnst,&
                         create_gesfinfo

  use guess_grids, only: hrdifsig_all,hrdifsfc_all,hrdifnst_all

  use gsi_4dvar, only: l4dvar, iwinbgn, winlen, nhr_assimilation
  
  use gridmod, only: nlat_sfc,nlon_sfc,lpl_gfs,dx_gfs
  use constants, only: zero,r60inv,r60,r3600,i_missing
  use obsmod, only: iadate
  use gsi_nstcouplermod, only: nst_gsi

  use read_obsmod, only: gsi_inquire

  !
  ! BAM model sigio
  !
  use sigioBAMMod, only : BAMFile
!  use sigio_BAMMod, only : BAM_Open, BAM_Close, BAMFile, &
!                           BAM_GetTimeInfo, BAM_GetDims, &
!                           BAM_GetOneDim,                &
!                           BAM_GetAvailUnit
  implicit none

! Declare passed variables
  integer(i_kind),intent(in   ) :: mype

! Declare local parameters
  integer(i_kind),parameter:: lunsfc=11
  integer(i_kind),parameter:: lunatm=12
  integer(i_kind),parameter:: lunnst=13
  integer(i_kind),parameter:: num_lpl=2000
  integer(i_kind),parameter:: max_file = 100
  real(r_kind),parameter:: r0_001=0.001_r_kind

! Declare local variables
  logical(4) fexist
  character(10) filename
  character(32) fileDir, fileFct
  integer(i_kind) i,j,iwan,npem1,iret
  integer(i_kind) nhr_half
  integer(i_kind) iamana(3)
  integer(i_kind) nminanl,nmings,nming2,ndiff
  integer(i_kind),dimension(4):: idateg
  integer(i_kind),dimension(2):: i_ges
  integer(i_kind),allocatable,dimension(:):: nst_ges
  integer(i_kind),dimension(5):: idate5
  integer(i_kind),dimension(num_lpl):: lpl_dum
  integer(i_kind),dimension(7):: idate
  integer(i_kind) :: nfhour, nfminute, nfsecondn, nfsecondd
  integer(i_kind),dimension(:),allocatable:: irec
  integer(i_llong) :: lenbytes
  real(r_single) hourg4
  real(r_kind) hourg,t4dv
  real(r_kind),allocatable,dimension(:,:):: time_atm
  real(r_kind),allocatable,dimension(:,:):: time_sfc
  real(r_kind),allocatable,dimension(:,:):: time_nst

!  type(sfcio_head):: sfc_head
!  type(sigio_head):: sigatm_head
!  type(nstio_head):: nst_head
!  type(nemsio_gfile) :: gfile_atm,gfile_sfc,gfile_nst

  type(BAMFile) :: BAM
  integer :: ierr
  integer :: lonb, latb

  integer(i_kind),dimension(5):: itime
  integer(i_kind),dimension(5):: ftime
  real(r_kind)    :: ijday, fjday
  integer(i_kind), external :: iw3jdn


  character(len=64), parameter :: myname_= 'read_bamfiles( )'
!-----------------------------------------------------------------------------
! Initialize variables
  nhr_half=nhr_assimilation/2
  if(nhr_half*2 < nhr_assimilation) nhr_half=nhr_half+1
  npem1=npe-1

  fexist=.true.
  nfldsig=0
  nfldsfc=0
  nfldnst=0
  iamana=0

! Check for non-zero length BAM files on single task
! ATM and SFC fields are at same file 
!
  if(mype==npem1)then

     allocate( irec(max_file) )
     irec=i_missing
     do i=0,max_file-1
        write(filename,'(''BAM.fct.'',i2.2)')i
        call gsi_inquire(lenbytes,fexist,filename,mype)

        ! counting and save time existing files
        if(fexist .and. lenbytes>0) then
           ! counting
           nfldsig = nfldsig+1
           nfldsfc = nfldsfc+1
           ! save time
           irec(nfldsig) = i
        end if
     enddo

     if(nfldsig==0) then
        write(6,*)trim(myname_),' ***ERROR*** NO valid BAM fields; aborting ', trim(filename)
        call stop2(169)
     end if

     
     !----------------------------------------------------!
     allocate(time_atm(nfldsig,2),time_sfc(nfldsfc,2))
     !
     ! Esta duas linhas foram adicionadas para corrigir
     ! um possivel bug caso o horario dos arquivos que 
     ! foram encontrados nao estajam na mesma janela de 
     ! de tempo do procedimento de AD.
     ! Se os arquivos nao estiverem na mesma janela de
     ! tempo estas duas matrizes ficam com valor zero
     ! e durante a leitura dos arquivos de first gues
     ! o gsi tentara ler os arquivos com a extencao '*.00'
     ! o que nao esta correto.
     !
     time_atm = -1
     time_sfc = -1
     !
     !----------------------------------------------------!

     write(6,'(A5,1x,I3.1,1x,A10)')'Found',nfldsig,'BAM files:'
     do i=1,nfldsig
        write(*,'(A4,1x,I3.1,1x,A8,I2.2)')'File',i,'BAM.dir.',irec(i)
     enddo

! Let a single task query the guess files.

!    Convert analysis time to minutes relative to fixed date
     call w3fs21(iadate,nminanl)
     write(6,*)trim(myname_),' :  analysis date,minutes ',iadate,nminanl

!    Check for consistency of times from atmospheric guess files.
     iwan=0
     do i=1,nfldsig

        write(fileFct,'("BAM.fct.",I2.2)') irec(i)
        write(fileDir,'("BAM.dir.",I2.2)') irec(i)

        call BAM%Open(trim(fileDir), istat=ierr)
        if (ierr.ne.0)then
           write(6,'(3(1x,A))')trim(myname_),' Problem to open/read BAM file: ',trim(fileDir)
           call stop2(99)
        endif

        !
        ! Get Time information
        !
        itime(:) = 0
        itime(1) = BAM%GetTimeInfo('iyr') ! Year of Century
        itime(2) = BAM%GetTimeInfo('imo') ! Month of year
        itime(3) = BAM%GetTimeInfo('idy') ! Day of month
        itime(4) = BAM%GetTimeInfo('ihr') ! Hour of day
        ijday = real(iw3jdn(itime(1),itime(2),itime(3)),r_kind)+&
                real(itime(4),r_kind)/24.0

        ftime(:) = 0
        ftime(1) = BAM%GetTimeInfo('fyr') ! Year of Century
        ftime(2) = BAM%GetTimeInfo('fmo') ! Month of year
        ftime(3) = BAM%GetTimeInfo('fdy') ! Day of month
        ftime(4) = BAM%GetTimeInfo('fhr') ! Hour of day
        fjday = real(iw3jdn(ftime(1),ftime(2),ftime(3)),r_kind)+&
                real(ftime(4),r_kind)/24.0

        hourg = abs(fjday-ijday)*24.0

        !
        ! Grid information
        !

        lonb = BAM%GetOneDim('IMax')
        latb = BAM%GetOneDim('JMax')

        i_ges(1) = lonb
        i_ges(2) = latb + 2

        !
        ! lons per lat
        !
        
        lpl_dum = 0
        lpl_dum(1:latb/2) = lonb

        call BAM%Close(  )

        !
        ! Convert time to minutes relative to fixed date
        !        
        call w3fs21(itime,nmings)
        ! nmings -> Integer number od minutes since 1 jan 1978
        ! itime(1) Year
        ! itime(2) Month
        ! itime(3) Day
        ! itime(4) Hour
        ! itime(5) Minute

        nming2 = nmings+60*hourg
        write(6,*)trim(myname_),' :  bam guess file, nming2 ',hourg,ftime,nming2
        t4dv = real((nming2-iwinbgn),r_kind)*r60inv
        if (l4dvar) then
           if (t4dv<zero .OR. t4dv>winlen) cycle
        else
           ndiff=nming2-nminanl
           if(abs(ndiff) > 60*nhr_half ) cycle
        endif

        iwan = iwan+1
        if(nminanl==nming2) then
           iamana(1) = iwan ! atm
           iamana(2) = iwan ! sfc
        endif

        time_atm(iwan,1) = t4dv
        time_atm(iwan,2) = irec(i)+r0_001

        time_sfc(iwan,1) = t4dv
        time_sfc(iwan,2) = irec(i)+r0_001
     end do

     !
     ! Caso as datas dos arquivos do BAM nao estejam na mesma janela
     ! de tempo do GSI o procedimento de AD sera abortado
     !
     if (sum(time_atm(2,:)).lt. 0.0 .or. &
         sum(time_sfc(2,:)).lt.0 )then
         write(6,*)trim(myname_),' ***ERROR*** No valid BAM file time; aborting'
         write(6,*)trim(myname_),' ***ERROR*** Verify BAM file times !'
         write(6,*)trim(myname_),' ***ERROR*** GSI time window !', winlen
        call stop2(169)
     endif

     deallocate( irec )
  end if

! Broadcast guess file information to all tasks
  call mpi_bcast(nfldsig,1,mpi_itype,npem1,mpi_comm_world,ierror)
  call mpi_bcast(nfldsfc,1,mpi_itype,npem1,mpi_comm_world,ierror)

!  if (nst_gsi > 0) call mpi_bcast(nfldnst,1,mpi_itype,npem1,mpi_comm_world,ierror)
  if(.not.allocated(time_atm)) allocate(time_atm(nfldsig,2))
  if(.not.allocated(time_sfc)) allocate(time_sfc(nfldsfc,2))
  if(.not.allocated(time_nst)) allocate(time_nst(nfldnst,2))

  call mpi_bcast(time_atm,2*nfldsig,mpi_rtype,npem1,mpi_comm_world,ierror)
  call mpi_bcast(time_sfc,2*nfldsfc,mpi_rtype,npem1,mpi_comm_world,ierror)

  if (nst_gsi > 0) call mpi_bcast(time_nst,2*nfldnst,mpi_rtype,npem1,mpi_comm_world,ierror)

  call mpi_bcast(iamana,3,mpi_rtype,npem1,mpi_comm_world,ierror)
  call mpi_bcast(i_ges,2,mpi_itype,npem1,mpi_comm_world,ierror)

  nlon_sfc=i_ges(1)
  nlat_sfc=i_ges(2)

  call mpi_bcast(lpl_dum,num_lpl,mpi_itype,npem1,mpi_comm_world,ierror)

  allocate(lpl_gfs(nlat_sfc/2))
  allocate(dx_gfs(nlat_sfc/2))

  lpl_gfs(1)=1  ! singularity at pole
  dx_gfs(1) = 360._r_kind / lpl_gfs(1)
  do j=2,nlat_sfc/2
     lpl_gfs(j)=lpl_dum(j-1)
     dx_gfs(j) = 360._r_kind / lpl_gfs(j)
  enddo

! Allocate space for guess information files
  call create_gesfinfo

! Load time information for atm guess field sinfo into output arrays
  ntguessig = iamana(1)
  do i=1,nfldsig
     hrdifsig(i) = time_atm(i,1)
     ifilesig(i) = nint(time_atm(i,2))
     hrdifsig_all(i) = hrdifsig(i)
  end do
  if(mype == 0) write(6,*)trim(myname_),' :  bam fcst files used in analysis  :  ',&
       (ifilesig(i),i=1,nfldsig),(hrdifsig(i),i=1,nfldsig),ntguessig

! Load time information for surface guess field info into output arrays
  ntguessfc = iamana(2)
  do i=1,nfldsfc
     hrdifsfc(i) = time_sfc(i,1)
     ifilesfc(i) = nint(time_sfc(i,2))
     hrdifsfc_all(i) = hrdifsfc(i)
  end do

  if(mype == 0) write(6,*)trim(myname_),' :  sfc fcst files used in analysis:  ',&
       (ifilesfc(i),i=1,nfldsfc),(hrdifsfc(i),i=1,nfldsfc),ntguessfc
  
  deallocate(time_atm,time_sfc)
  
!! Load time information for nst guess field info into output arrays
!  ntguesnst = iamana(3)
!  if ( nst_gsi > 0 ) then
!    do i=1,nfldnst
!       hrdifnst(i) = time_nst(i,1)
!       ifilenst(i) = nint(time_nst(i,2))
!       hrdifnst_all(i) = hrdifnst(i)
!    end do
!    if(mype == 0) write(6,*)trim(myname_),' :  nst fcst files used in analysis:  ',&
!         (ifilenst(i),i=1,nfldnst),(hrdifnst(i),i=1,nfldnst),ntguesnst
!    deallocate(time_nst)
!  endif

! End of routine
  return
end subroutine read_bamfiles

