!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: gesinfo_cptec
!
! !DESCRIPTION:
!   This subroutine reads metadata from CPTEC BAM guess files and sets
!   vertical coordinate information, guess time, and analysis time structure
!   used by the GSI system.
!
! !INTERFACE:
subroutine gesinfo_cptec(nhr_assimilation, hourg, idate4)

  use kinds,        only: i_kind, r_kind
  use gridmod,      only: nsig, jcap_b, ak5, bk5, ck5, tref5, idvc5, idsl5, &
                          idlsm, nlat, nlon
  use mpimod,       only: npe,mype
  use obsmod,       only: iadatemn
  use constants,    only: zero, h300
  use sigioBAMMod,  only: BAMFile
  use m_chars,      only: uppercase
  implicit none
!
! !INPUT/OUTPUT PARAMETERS:
  integer(i_kind), intent(in)    :: nhr_assimilation
  real(r_kind),    intent(out)   :: hourg
  integer(i_kind), intent(out)   :: idate4(5)
!
! !REVISION HISTORY:
!   09 Jul 2025 - J. Gerd - Initial modularization from GSI `gesinfo`
!
! !EXAMPLE:
!   call gesinfo_cptec(nhr_assimilation, hourg, idate4)
!
!EOP
!-------------------------------------------------------------------------------
  
  type(BAMFile) :: BAM
  integer       :: IMax, JMax, KMax, Mend
  integer       :: nvcoord, ierr, k
  integer       :: mype_out 
  character(len=32) :: fileDir, fileFct
  real(r_kind)    :: ijday, fjday
  integer(i_kind), dimension(4) :: itime, ftime

  integer(i_kind), external :: iw3jdn
  external stop2

  ! File naming
  write(fileFct,'("BAM.fct.",I2.2)') nhr_assimilation
  write(fileDir,'("BAM.dir.",I2.2)') nhr_assimilation

  inquire(file=fileDir, exist=ierr)
  if (.not. ierr) then
     write(6,*)'GESINFO_CPTEC:  BAM guess file not available: ', trim(fileDir)
     call stop2(99)
  end if

  ! Open BAM file
  call BAM%Open(trim(fileDir), istat=ierr)
  if (ierr /= 0) then
     write(6,*)'GESINFO_CPTEC: Failed to open BAM file: ', trim(fileDir), ' Status: ', ierr
     call stop2(99)
  end if

  call BAM%GetDims(IMax, JMax, KMax, Mend)

  if (KMax /= nsig .or. Mend /= jcap_b) then
     write(6,'(A,4I6)') 'GESINFO_CPTEC: BAM dims (IMax,JMax,KMax,Mend):', IMax, JMax, KMax, Mend
     write(6,'(A,4I6)') 'GESINFO_CPTEC: Namelist (nlon,nlat,nsig,jcap_b):', nlon, nlat, nsig, jcap_b
     call stop2(99)
  end if

  ! Get time info
  itime = (/ BAM%GetTimeInfo('iyr'), BAM%GetTimeInfo('imo'), &
             BAM%GetTimeInfo('idy'), BAM%GetTimeInfo('ihr') /)
  ftime = (/ BAM%GetTimeInfo('fyr'), BAM%GetTimeInfo('fmo'), &
             BAM%GetTimeInfo('fdy'), BAM%GetTimeInfo('fhr') /)

  ijday = real(iw3jdn(itime(1),itime(2),itime(3)), r_kind) + real(itime(4), r_kind)/24.0
  fjday = real(iw3jdn(ftime(1),ftime(2),ftime(3)), r_kind) + real(ftime(4), r_kind)/24.0
  hourg = abs(fjday - ijday) * 24.0

  idate4 = (/ itime(4), itime(2), itime(3), itime(1), 0 /) ! hour, month, day, year, minute
  iadatemn = (/ itime(1), itime(2), itime(3), itime(4), 0 /)

  ! Surface model info
  idlsm = BAM%getPhysics('lsm')

  !Load vertical coordinate structure

  ! vertical coordinate id (idvc5)
  !      1 for sigma, 2 for ec-hybrid, 3 for ncep hybrid
  !
  ! number of vcoord profiles
  !   for idvc=1, nvcoord=1: sigma interface
  !   for idvc=2, nvcoord=2: hybrid interface a and b
  !
  ! type of sigma structure (idsl) for sigma case
  !   1 for phillips or 2 for mean
  !
  
  if (BAM%isHybrid)then
  ! BAM hybrid vertical coordenate
     idvc5   = 2
     nvcoord = 2
  else
  ! BAM sigma-p vertical coordenate
     idvc5   = 1
     idsl5   = 1
     nvcoord = 1
  endif

  ! Load surface pressure and thermodynamic variable ids!
  !   Initializing coefficients for generalized
  !   vertical coordinate
  !

  select case (nvcoord)
    case (1) ! Sigma
       call BAM%GetVerticalCoord('si', bk5(1:nsig+1))

    case (2) !Hybrid -> Sigma-P
       call BAM%GetVerticalCoord('ak', ak5(1:nsig+1))
       call BAM%GetVerticalCoord('bk', bk5(1:nsig+1))
       ak5(1:nsig+1) = ak5(1:nsig+1) * 0.001_r_kind

    case default
       write(6,*)'GESINFO_CPTEC: Invalid nvcoord value:', nvcoord
       call stop2(85)
  end select
  
  !
  !    Load reference temperature array (used by general coordinate)        
  !

  tref5(1:nsig) = h300

  !    Echo select header information to stdout
  mype_out = npe / 2
  if(mype==mype_out) then
     write(6,120) Mend,KMax,JMax,IMax,&
                  idlsm, nvcoord
120  format('GESINFO:  jcap_b=',i4,', levs=',i3,', latb=',i5,&
          ', lonb=',i5,', idlsm=', i3,&
          ', nvcoord=',i3)

     do k=1,nsig
        write(6,130) k,ak5(k),bk5(k),ck5(k),tref5(k)
     end do
     k=nsig+1
     write(6,130) k,ak5(k),bk5(k),ck5(k)
130  format(3x,'k,ak,bk,ck,tref=',i3,1x,4(g19.12,1x))
  endif

  !
  ! Close BAM header file
  !
  call BAM%Close( )

  if (mype==mype_out) &
       write(6,*)'GESINFO:  READ CPTEC DIR FILE OK'

end subroutine gesinfo_cptec

