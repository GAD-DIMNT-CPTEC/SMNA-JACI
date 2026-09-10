module m_stdio
  use EndianUtility
  use TypeKinds
!
! A F90 module defines std. I/O parameters
!
! Description
!   Define system dependent I/O parameters.
!
! History
!   * 23 Mar 2011 - J.G. de Mattos - Initial Code
!
  implicit none
  private

  public   :: stdin    ! a unit linked to UNIX stdin
  public   :: stdout   ! a unit linked to UNIX stdout
  public   :: stderr   ! a unit linked to UNIX stderr
  public   :: openIEEE
  public   :: readIEEE
  public   :: getAvailUnit

  ! Defines standar i/o units.

  integer, parameter :: stdin  = 5
  integer, parameter :: stdout = 6

  ! Generic setting for UNIX other than HP-UX

  integer, parameter :: stderr = 0

  contains

  function openIEEE(unit, file, status, action)result(iret)
     integer(Long),    intent(in) :: unit
     character(len=*), intent(in) :: file
     character(len=*), intent(in) :: status
     character(len=*), intent(in) :: action
     integer(Long)                :: iret

    open( unit   = unit,          &
          file   = trim(file),    &
          status = trim(status),  &
          form   = 'unformatted', &
          action = trim(action),  &
          access = 'stream',      &
          iostat = iret           &
     )

  end function

!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: SRead1D_ - read and return a 1D sigle precision BAM field.
!
! 
! !DESCRIPTION: Esta rotina lê e retorna um campo modelo BAM em um vetor 1D em
!               precisão simples
!
! !INTERFACE:
!   
   subroutine readIEEE(funit, fld, byteSwap, idx, istat)
      implicit none
!
! !INPUT PARAMETERS:
! 
      integer(Long),            intent(in   ) :: funit
      integer(LLong), optional, intent(in   ) :: idx
      logical,                  intent(in   ) :: byteSwap ! .true. machine is big Endian
!
! !OUTPUT PARAMETERS:
! 
      real(Single),              intent(  out) :: fld(:)
      integer(Long),   optional, intent(  out) :: istat
!
! !REVISION HISTORY: 
!
!  11 Oct 2016 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
      character(len=*), parameter          :: myname_ = ':: readIEEE( ... )'

      integer :: iret
      character(len=64) :: msg

      if(present(istat)) istat = 0

      if(present(idx))then
         read(unit = funit, POS=idx, iostat = iret, iomsg= msg) fld
      else
         read(unit = funit, iostat = iret) fld
      endif
      if (iret.ne.0)then
         write(stdout,*)trim(myname_),': error to read field, ',iret
         write(stdout,*)trim(myname_),':',trim(msg)
         if(present(istat)) istat = iret
      endif

      ! if machine is little endian do byteSwap
      if (.not.byteSwap) fld = swap_endian(fld)

   end subroutine
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !FUNCTION: GetAvailUnit
!
! !DESCRIPTON: function to return next available logical unit
!
!
!             
!                 
! !INTERFACE:
!
   function getAvailUnit( exclude ) result( lu )

      implicit none
!
! !INPUT PARAMETERS:
!
      ! Skip this logical unit      
      integer, optional :: exclude(:)
!
! !RETURN VALUE:
!
      integer :: lu
!
! !REVISION HISTORY: 
!  11 Oct 2016 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
      character(len=*), parameter :: myname_ = ':: getAvailUnit( ... )'

      integer,          parameter :: MaxLogicalUnit = 254
      integer :: iunit
      integer :: ios
      logical :: isopen
      integer :: i

      

      ! start open loop for lun search
      find_unit:do iunit = 10, MaxLogicalUnit

         if (present(exclude))then
            do i=1,size(exclude)
               if (iunit.eq.exclude(i)) cycle find_unit
            enddo
         endif

         inquire (Unit = iunit, opened = isopen, iostat = ios )
         
         if (.not.isopen.and.ios.eq.0)then
            lu = iunit
            return
         endif

         if(iunit .eq. MaxLogicalUnit)then
             write(stderr,'(4A)')'Error at ',trim(myname_),':','Units from 10 to 254 are already in use!'
            stop
         endif

      end do find_unit

   end function GetAvailUnit


end module m_stdio
