!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !MODULE: cptecbam_io -- Implements a CPTEC BAM interface to Read/Write
!                         atmospheric and surface at GSI.
!
! !DESCRIPTON:  
!              
!             
!                 
!\\
!\\
! !INTERFACE:
!
module cptecbam_io

! GSI kinds
   use kinds, only: r_kind

   implicit none
   private
!
! !PUBLIC MEMBER FUNCTIONS:
!

   public :: read_bam    ! read bam atmopheric guess fields
   public :: read_bamsfc ! read bam surface guess fields
   public :: write_bam   ! write bam atmospheric analysis fields
!   public :: write_bamsfc! write bam surface analysis fields

!
! !PRIVATE MEMBER VARIABLES:
!
   !
   ! Transfers Variables
   !

   real(r_kind), allocatable :: g_z  ( :,: ) ! Orography (m)
   real(r_kind), allocatable :: g_ps ( :,: ) ! Surface Pressure (kPa)
   real(r_kind), allocatable :: g_tv (:,:,:) ! Virtural Temperature (K)
   real(r_kind), allocatable :: g_vor(:,:,:) ! Vorticity (m s^-1)
   real(r_kind), allocatable :: g_div(:,:,:) ! Divergence (m s^-1)
   real(r_kind), allocatable :: g_u  (:,:,:) ! Zonal Wind (m/s)
   real(r_kind), allocatable :: g_v  (:,:,:) ! Meridional Wind (m/s)
   real(r_kind), allocatable :: g_q  (:,:,:) ! Specific Humidy (kg/kg)
   real(r_kind), allocatable :: g_ql (:,:,:) ! Liq mixing ratio prognostic (kg/kg)
   real(r_kind), allocatable :: g_qi (:,:,:) ! Ice mixing ratio prognostic (kg/kg)
   real(r_kind), allocatable :: g_cw (:,:,:) ! Total Cloud Water Content (kg/kg)
   real(r_kind), allocatable :: g_oz (:,:,:) ! Ozone

!
! write out units
!

   integer, parameter :: stdinp = 5 ! standard input
   integer, parameter :: stdout = 6 ! standard output
   integer, parameter :: stderr = 0 ! standard error
!
! Undef value
!

   integer, parameter :: iUdef = -9999
   real,    parameter :: Udef  = -9.99E9


!
! !REVISION HISTORY:
!
!   28 Apr 2016 - J. G. de Mattos -  Initial code.
!   07 Mar 2021 - J. G. de Mattos -  Change Name of microphycis Variables 
!
!EOP
!-----------------------------------------------------------------------------!

   character(len=64),parameter :: myname='cptecbam_io'

   contains

!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !ROUTINE: GenReadBAM
!              
! !DESCRIPTON: rotine to get atmospheric guess fields from BAM and return to
!              GSI
!             
!                 
!\\
!\\
! !INTERFACE:
!

   subroutine read_bam( mype )
!
! !USES:
!
      ! GSI kinds
      use kinds, only: r_kind

      !GSI bundle
      use gsi_metguess_mod, only: gsi_metguess_bundle
      use gsi_bundlemod, only: gsi_bundlegetpointer, gsi_bundleprint

      use guess_grids, only: nfldsig 

      implicit none
!
! !INPUT PARAMETERS:
!
      integer, intent(in   ) :: mype ! mpi task id

! !REVISION HISTORY:
!
!   28 Apr 2016 - J. G. de Mattos -  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      character(len=64),parameter :: myname_=trim(myname)//' :: read_bam( )'

      !
      ! GSI bungle guess fields
      !

      real(r_kind), pointer, dimension(:,:  ) :: ges_z_it    => NULL()
      real(r_kind), pointer, dimension(:,:  ) :: ges_ps_it   => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: ges_u_it    => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: ges_v_it    => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: ges_vor_it  => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: ges_div_it  => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: ges_tv_it   => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: ges_q_it    => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: ges_oz_it   => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: ges_cw_it   => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: ges_ql_it   => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: ges_qi_it   => NULL()
  
      !
      ! Auxiliary vars
      !

      integer :: istatus
      integer :: it, i

#ifdef DEBUG
    WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif

      !
      ! Loop over all time gues fields
      !

      do it = 1, nfldsig

         !
         ! Allocating transfer fields
         !

         call AllocateAtmFields( )

 
         !
         ! Read BAM fields
         !

         call GenReadBAM ( mype, it )
         !
         ! Put fields at GSI Bundle
         !

         ! Orography

         call gsi_bundlegetpointer (gsi_metguess_bundle(it), 'z', ges_z_it, istatus) 
         if(istatus.eq.0) ges_z_it = g_z
         

         ! Surface Pressure

         call gsi_bundlegetpointer (gsi_metguess_bundle(it), 'ps', ges_ps_it, istatus) 
         if(istatus.eq.0) ges_ps_it = g_ps 

         ! Zonal wind component

         call gsi_bundlegetpointer (gsi_metguess_bundle(it), 'u', ges_u_it, istatus) 
         if(istatus.eq.0) ges_u_it = g_u

         ! Meridional wind component

         call gsi_bundlegetpointer (gsi_metguess_bundle(it), 'v', ges_v_it, istatus) 
         if(istatus.eq.0) ges_v_it = g_v
      
         ! Vorticity

         call gsi_bundlegetpointer (gsi_metguess_bundle(it), 'vor', ges_vor_it, istatus) 
         if(istatus.eq.0) ges_vor_it = g_vor

         ! Divergency

         call gsi_bundlegetpointer (gsi_metguess_bundle(it), 'div', ges_div_it, istatus) 
         if(istatus.eq.0) ges_div_it = g_div

         ! Virtual Temperature

         call gsi_bundlegetpointer (gsi_metguess_bundle(it), 'tv', ges_tv_it, istatus) 
         if(istatus.eq.0) ges_tv_it = g_tv

         ! Specific Humidity

         call gsi_bundlegetpointer (gsi_metguess_bundle(it), 'q', ges_q_it, istatus) 
         if(istatus.eq.0) ges_q_it = g_q

         ! Ozone

!         call gsi_bundlegetpointer (gsi_metguess_bundle(it), 'oz', ges_oz_it, istatus) 
!         if(istatus.eq.0) ges_oz_it = g_oz

         ! Total Cloud Water Content

         call gsi_bundlegetpointer (gsi_metguess_bundle(it), 'cw', ges_cw_it, istatus) 
         if(istatus.eq.0) ges_cw_it = g_ql + g_qi

         
         ! Liq mixing ratio prognostic

         call gsi_bundlegetpointer (gsi_metguess_bundle(it), 'ql',ges_ql_it, istatus) 
         if(istatus.eq.0) ges_ql_it = g_ql

         ! Ice mixing ratio prognostic

         call gsi_bundlegetpointer (gsi_metguess_bundle(it), 'qi', ges_qi_it, istatus)           
         if(istatus.eq.0) ges_qi_it = g_qi

         !
         ! Clear ges transfer fields
         !

         call DeallocateAtmFields ( )

     enddo

  end subroutine
!EOC
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !ROUTINE: GenReadBAM
!
! !DESCRIPTON: rotine to read atmospheric guess fields
!
!             
!                 
!\\
!\\
! !INTERFACE:
!

   subroutine GenReadBAM( mype, it )
!
! !USES:
!
      ! Model read 
      use sigioBAMMod, only : BAMFile
      use MiscMod, only: r4, r8

      ! GSI grid informations
      use gridmod, only : nlon, nlat, lat2, lon2, jcap, jcap_b, nsig
      use gridmod, only : grd=>grd_a

      ! GSI file ges information
      use guess_grids, only: ifilesig 

      !GSI kinds
      use kinds, only : r_kind

      !MPI
      use mpimod, only : NPe, mpi_comm_world, mpi_itype


      implicit none
!
! !INPUT PARAMETERS:
!
      integer, intent(in   ) :: mype ! mpi task id
      integer, intent(in   ) :: it   ! gues field number
!
! !PARAMTERS:
!

      integer, parameter                  :: NVars = 8
      character(len=40), dimension(NVars) :: VName = [                             &
                                                     'TOPOGRAPHY                 ',& !1
                                                     'LN SURFACE PRESSURE        ',& !2
                                                     'VIRTUAL TEMPERATURE        ',& !3
                                                     'DIVERGENCE                 ',& !4
                                                     'VORTICITY                  ',& !5
                                                     'SPECIFIC HUMIDITY          ',& !6
                                                     'LIQ MIXING RATIO PROGNOSTIC',& !7
                                                     'ICE MIXING RATIO PROGNOSTIC' & !8
                                                     ] 
  

! !REVISION HISTORY:
!
!   28 Apr 2016 - J. G. de Mattos -  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC      
      character(len=64), parameter :: myname_=trim(myname)//' :: GenReadBAM( )'

      type(BAMFile) :: BAM  ! BAM files data type

      real(r8),     allocatable :: grid_in(:,:)
      real(r8),     allocatable :: grid_b(:,:)
      real(r8),     allocatable :: grid(:,:)
      real(r_kind) :: work(grd%itotsub)

      !
      !  variables in spectral space
      !
  
      real(r8), dimension(:),   allocatable :: divq
      real(r8), dimension(:),   allocatable :: vorq
      real(r8), dimension(:),   allocatable :: uveq
      real(r8), dimension(:),   allocatable :: vveq
  
      !
      !  variables in physical space
      !
  
      real(r8), dimension(:,:), allocatable :: gu, gu_
      real(r8), dimension(:,:), allocatable :: gv, gv_
      real(r8) :: guW(grd%itotsub)
      real(r8) :: gvW(grd%itotsub)


      !
      ! Auxiliary Variables
      !

      real(r4), dimension( : ), allocatable :: clat

      integer :: istat

      integer, allocatable :: WVar(:)
      integer, allocatable :: WLev(:)

      integer :: Nflds
      integer :: imax, imax_b
      integer :: jmax, jmax_b
      integer :: kmax
      integer :: Mend
      integer :: MnWv2
      integer :: MnWv3
      integer, allocatable :: nlevs(:)

      integer :: i, j, k, ij
      integer :: ii,jj,kk
      integer :: icount
      integer :: ivar
      integer :: iret
      integer :: ilev
      integer :: iPe
      integer :: LastPeUsed

      integer :: WrkPe
      character(len=80) :: linha
      character(len=80) :: fileFct
      character(len=80) :: fileDir

      
      !
      ! Define files names
      ! fct => sequential binary ieee file (spectral and grid points field)
      ! dir => ascii file with field information
      !

      !
      ! Files to read and get first-guess
      ! !! hardwire files! not a best way !!
      !

#ifdef DEBUG
      WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif

      write(fileFct,'("BAM.fct.",I2.2)') ifilesig(it)
      write(fileDir,'("BAM.dir.",I2.2)') ifilesig(it)

      !
      ! All tasks open and read header
      !
      Nflds      = NVars * grd%nsig
      LastPeUsed = Nflds - 1
      WrkPe      = 0

      if (MyPe .eq. 0) WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)

      if (MyPe <= LastPeUsed) then

         call BAM%Open(trim(fileDir), trim(fileFct), istat=istat)

         if (istat.ne.0)then

            write(stdout,'(2(A,1x),I4,1x,A,1x,I4)')trim(myname_),'Problem to open/read BAM files at ',MyPe, 'rank! istat,', istat

            return
         endif

      endif

      ! Get some info to broadcast for all Pe's

      allocate(NLevs(NVars))

      if (MyPe .eq. WrkPe)then

           do ivar = 1, Nvars
              NLevs(ivar) = BAM%GetNlevels(trim(VName(ivar)),iret)
              if(iret .ne. 0)then
                 write(6,'(2A,1x,I4)')trim(myname_),':: *** Error ***: BAM file should be openned at', MyPe
                 stop
              endif         
           enddo

          call BAM%GetDims(imax_b, jmax_b, kmax, Mend)

          !
          ! consistency check
          !

          MnWv2 = (Mend+1)*(Mend+2)
          MnWv3 = (MnWv2+2)*(Mend+1)

          ! from namelist and gridmod
          imax  = grd%nlon
          jmax  = grd%nlat-2

          if ( Mend .ne.  jcap_b      .or. &
               KMax .ne.  grd%nsig   .and. &
               MyPe .eq. 0                 &
             )then

             write (6,'(3(1x,A))') trim(myname_),':: *** ERROR *** reading ', trim(BAM%fBin)
             write (6,'(2(1x,A,1x,I4))')'<TRC>',Mend,'<->.',jcap_b
             write (6,'(2(1x,A,1x,I4))')'<KMax>',KMax,'<->.',grd%nsig

             return
          endif

      endif
      ! Broadcast some info to all tasks
      call mpi_bcast(imax_b,    1,mpi_itype,WrkPe,mpi_comm_world,iret)
      call mpi_bcast(jmax_b,    1,mpi_itype,WrkPe,mpi_comm_world,iret)
      call mpi_bcast(  imax,    1,mpi_itype,WrkPe,mpi_comm_world,iret)
      call mpi_bcast(  jmax,    1,mpi_itype,WrkPe,mpi_comm_world,iret)
      call mpi_bcast(  kmax,    1,mpi_itype,WrkPe,mpi_comm_world,iret)
      call mpi_bcast(  Mend,    1,mpi_itype,WrkPe,mpi_comm_world,iret)
      call mpi_bcast( MnWv2,    1,mpi_itype,WrkPe,mpi_comm_world,iret)
      call mpi_bcast( MnWv3,    1,mpi_itype,WrkPe,mpi_comm_world,iret)
      call mpi_bcast( NLevs,NVars,mpi_itype,WrkPe,mpi_comm_world,iret)

      !
      ! Process guess fields according to type of input file.   BAM_SIGIO files
      ! are spectral coefficient files and need to be transformed to the grid.
      ! Once on the grid, fields need to be scattered from the full domain to 
      ! sub-domains.
      !

      allocate(wvar(NPe))
      allocate(wlev(NPe))

      iCount = 1

      do ivar = 1, NVars
         do ilev = 1, NLevs(ivar)
         
           iPe = iCount - 1

            !----------------------------------------
            ! This is used to scatter between all Pe's
            ! Como cada Pe ira pegar um campo
            ! em um nivel, estas duas variaveis dirao
            ! qual variavel é lida em cada Pe
            !
            wvar(iCount) = ivar ! What Variable
            wlev(iCount) = ilev ! What Level
            !----------------------------------------

            if ( MyPe .eq. iPe )then

               if(jcap.ne.jcap_b)then

                  allocate(grid_b(imax_b,jmax_b))

                  call BAM%GetField(trim(VName(ivar)), ilev, grid_b, istat=istat)

                  allocate(grid(imax,jmax))

                  call lterp(grid_b, jcap_b, jcap, grid)
                  
                  deallocate(grid_b)

               else

                  allocate(grid(imax,jmax))

                  call BAM%GetField(trim(VName(ivar)), ilev, grid, istat=istat)
               
               endif

               !
               ! some necessary adjusts to Specific Humidity field
               !

               if (ivar .eq. 6 )then
                  do j=1,jmax
                     do i=1,imax
                        grid(i,j)=MAX(1.0e-12_r8,grid(i,j))
                     enddo
                  enddo
               endif
               

               !
               ! reorganize grid to be used by GSI
               !
               !  - adds a southern and northern latitude row to the input grid,
               !    this is necessary by GSI to make some interpolations
               !
               !  - reorder the output array so that it is a one-dimensional 
               !    array read in an order consistently with that assumed for total
               !    domain gsi grids.
               !

               call GenFill_ns(grid,work)

               deallocate(grid)

            endif

            if ( iCount .eq. NPe .or. ( ivar .eq. NVars .and. ilev .eq. NLevs(ivar)) )then

               !
               ! Transfer contents of 2-d array global to 3-d subdomain array
               ! Fields are scattered from the full domain to sub-domains.
               !
               call GenReload( work,  & ! Input Field read by Pe
                               WVar,  & !       Position in Var list
                               WLev,  & !       Level
                               icount,& !       MaxCount until here
                               g_z,   & ! OutPut Topography
                               g_ps,  & !        Surface Pressure
                               g_tv,  & !        Virtural Temperature
                               g_vor, & !        Vorticity
                               g_div, & !        Divergency
                               g_q,   & !        Specific Humidity
                               g_ql,  & !        Liq mixing ratio prognostic
                               g_qi   & !        Ice mixing ratio prognostic
                               )

               !
               ! after transfer and scatter reset all counters
               !
               
               WVar   = 0
               WLev   = 0
               icount = 0

            endif

            icount = icount + 1
         enddo

      enddo

      deallocate(wvar)
      deallocate(wlev)


      ! Convert Surface Pressure from ln(pslc) to pslc.
      ! NOTA: GSI use pslc in millibar but this conversion
      !       is made by internals subroutines. So here we 
      !       need only convert from ln. At this point
      !       pslc need be in centibar.

     
      if (bam%isHybrid)then
!$omp parallel do  schedule(dynamic,1) private(i,j)      
         do j = 1, lon2
            do i = 1, lat2
               g_ps(i,j) =  exp (g_ps(i,j))*1e-3_r4
            enddo
         enddo
!$omp end parallel do

      else
!$omp parallel do  schedule(dynamic,1) private(i,j)      
         do j = 1, lon2
            do i = 1, lat2
               g_ps(i,j) =  exp (g_ps(i,j))
            enddo
         enddo
!$omp end parallel do
      endif
      !
      ! Compute u and v from Vorticity and Divergence
      ! One level per Pe!
      ! convert from spectral to grid and scatter from the
      ! full domain to sub-domains.
      !
      guW = udef
      gvW = udef
      
      allocate(wlev(nPe))
      wlev   = 0
      icount = 1

      do ilev = 1, kmax

         iPe = icount - 1

         wlev(icount) = ilev

         if (MyPe .eq. iPe)then

            allocate(gu(imax_b,jmax_b))
            allocate(gv(imax_b,jmax_b))

            call BAM%getUV(ilev, gu, gv)

            allocate(clat(JMax_b))

            call BAM%GetWCoord('clat', clat, iret)

            do j=1,jmax_b
               gu(1:imax_b,j) = gu(1:imax_b,j)/clat(j)
               gv(1:imax_b,j) = gv(1:imax_b,j)/clat(j)
            enddo

            deallocate(clat)

            !
            ! Interpola se o background nao for igual a analise
            !
            if(jcap .ne. jcap_b)then

               allocate(gu_(imax,jmax))
               call lterp(gu, jcap_b, jcap, gu_)
               deallocate(gu)

               allocate(gv_(imax,jmax))
               call lterp(gv, jcap_b, jcap, gv_)
               deallocate(gv)

               call GenFillUV_ns(gu_, gv_, guW, gvW)
               deallocate (gu_)
               deallocate (gv_)

            else
               call GenFillUV_ns(gu, gv, guW, gvW)
               deallocate (gu)
               deallocate (gv)
            endif

         endif

         if (iCount .eq. NPe .or. iLev .eq. kmax)then
            ! Transfer contents of 2-d array global to 3-d subdomain array
            ! Fields are scattered from the full domain to sub-domains.
            !
            
 
            call GenReload_uv ( guW,   & ! Input U component
                                gvW,   & !       V component
                                wlev,  & !       Level for each Pe
                                icount,& !       MaxCount until this step
                                g_u,   & ! OutPut U component
                                g_v    & !        V component
                              )

            iCount = 0
            wlev   = 0

         endif
         icount = icount + 1
      enddo
      deallocate(wlev)
      !
      ! Close BAM Files
      !
      if (MyPe <= LastPeUsed) call BAM%Close(iret)

      return
   end subroutine
!EOC
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !ROUTINE: GenReload
!
! !DESCRIPTON: rotine to allocate ges transfer fields
!
!             
!                 
!\\
!\\
! !INTERFACE:
!

   subroutine AllocateAtmFields( )

      use gridmod, only: lat2, lon2, nsig
      implicit none

! !REVISION HISTORY:
!
!   28 Apr 2016 - J. G. de Mattos -  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC  
      character(len=64), parameter :: myname_=trim(myname)//' :: AllocateAtmFields( )'


#ifdef DEBUG
    WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif

      allocate( g_z  (lat2, lon2) )
      allocate( g_ps (lat2, lon2) )
      allocate( g_tv (lat2, lon2, nsig) )
      allocate( g_vor(lat2, lon2, nsig) )
      allocate( g_div(lat2, lon2, nsig) )
      allocate( g_u  (lat2, lon2, nsig) )
      allocate( g_v  (lat2, lon2, nsig) )
      allocate( g_q  (lat2, lon2, nsig) )
!      allocate( g_oz (lat2, lon2, nsig) )
      allocate( g_cw (lat2, lon2, nsig) )
      allocate( g_ql (lat2, lon2, nsig) )
      allocate( g_qi (lat2, lon2, nsig) )

   end subroutine

!EOC
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !ROUTINE: GenReload
!
! !DESCRIPTON: rotine to Deallocate ges transfer fields
!
!             
!                 
!\\
!\\
! !INTERFACE:
!

   subroutine DeallocateAtmFields( )

      implicit none

! !REVISION HISTORY:
!
!   28 Apr 2016 - J. G. de Mattos -  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC      
      character(len=64), parameter :: myname_=trim(myname)//' :: DeAllocateAtmFields( )'


#ifdef DEBUG
    WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif

      deallocate( g_z  )
      deallocate( g_ps )
      deallocate( g_tv )
      deallocate( g_vor)
      deallocate( g_div)
      deallocate( g_u  )
      deallocate( g_v  )
      deallocate( g_q  )
!      deallocate( g_oz )
      deallocate( g_cw )
      deallocate( g_ql )
      deallocate( g_qi )

   end subroutine

!EOC
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !ROUTINE: GenReload
!
! !DESCRIPTON: rotine to transfer contents of 2-d array global to 3-d subdomain array
!              and scatter from the full domain to sub-domains.
!             
!                 
!\\
!\\
! !INTERFACE:
!

subroutine GenReload( work, wvar, wlev, icount, &
                      g_z,g_ps,g_tv,g_vor,g_div,g_q,g_ql,g_qi)

! !USES:

  use kinds, only: r_kind,i_kind
  use mpimod, only: npe,mpi_comm_world,ierror,mpi_rtype
  use gridmod, only: grd => grd_a
  implicit none

! !INPUT PARAMETERS:

  real(r_kind),dimension(grd%itotsub)   ,intent(in   ) :: work
  integer(i_kind),dimension(npe)        ,intent(in   ) :: wvar
  integer(i_kind),dimension(npe)        ,intent(in   ) :: wlev
  integer(i_kind)                       ,intent(in   ) :: icount

! !OUTPUT PARAMETERS:

  real(r_kind),dimension(:,:),  intent(  out) :: g_z
  real(r_kind),dimension(:,:),  intent(  out) :: g_ps
  real(r_kind),dimension(:,:,:),intent(  out) :: g_tv
  real(r_kind),dimension(:,:,:),intent(  out) :: g_vor
  real(r_kind),dimension(:,:,:),intent(  out) :: g_div
  real(r_kind),dimension(:,:,:),intent(  out) :: g_q
  real(r_kind),dimension(:,:,:),intent(  out) :: g_ql
  real(r_kind),dimension(:,:,:),intent(  out) :: g_qi


! !DESCRIPTION: Transfer contents of 2-d array global to 3-d subdomain array
!
! !REVISION HISTORY:
!   2004-05-14  treadon
!   2004-07-15  todling, protex-compliant prologue
!   2014-12-03  derber     - introduce vdflag and optimize routines
!   2016-05-31  de Mattos  - adapt to BAM model
!
! !REMARKS:
!
!   language: f90
!   machine:  ibm rs/6000 sp; sgi origin 2000; compaq/hp
!
! !AUTHOR:
!   treadon          org: np23                date: 2004-05-14
!
!EOP
!-------------------------------------------------------------------------
!BOC
      character(len=64), parameter :: myname_=trim(myname)//' :: GenReload( )'

  integer(i_kind) :: i,j,k,ij
  integer(i_kind) :: Pe
  integer(i_kind) :: var
  integer(i_kind) :: lev
  real(r_kind),dimension(grd%lat2*grd%lon2,npe):: sub


#ifdef DEBUG
    WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif

  call mpi_alltoallv(work,          & 
                     grd%ijn_s,     &
                     grd%displs_s,  &
                     mpi_rtype,     &
                     sub,           &
                     grd%irc_s,     &
                     grd%ird_s,     &
                     mpi_rtype,     &
                     mpi_comm_world,&
                     ierror)

!$omp parallel do  schedule(dynamic,1) private(Pe,i,j,ij,var,lev)
  do Pe = 1, icount

     var = wvar(Pe)
     lev = wlev(Pe)

     select case ( var )

        case (1) ! Topography

           ij=0
           do j=1,grd%lon2
              do i=1,grd%lat2
                 ij=ij+1
                 g_z(i,j)=sub(ij,Pe)
              end do
           end do

        case (2) ! Surface Pressure

           ij=0
           do j=1,grd%lon2
              do i=1,grd%lat2
                 ij=ij+1
                 g_ps(i,j)=sub(ij,Pe)
              end do
           end do

        case (3) ! Virtual Temperature

           ij=0
           do j=1,grd%lon2
              do i=1,grd%lat2
                 ij=ij+1
                 g_tv(i,j,lev)=sub(ij,Pe)
              end do
           end do

        case (4) ! Divergence

          ij=0
          do j=1,grd%lon2
             do i=1,grd%lat2
                ij=ij+1
                g_div(i,j,lev)=sub(ij,Pe)
             end do
          end do

        case (5) ! Vorticity

          ij=0
          do j=1,grd%lon2
             do i=1,grd%lat2
                ij=ij+1
                g_vor(i,j,lev)=sub(ij,Pe)
             end do
          end do

        case (6) ! Specific Humidity

           ij=0
           do j=1,grd%lon2
              do i=1,grd%lat2
                 ij=ij+1
                 g_q(i,j,lev)=sub(ij,Pe)
              end do
           end do


        case (7) ! Liq mixing ratio prognostic

           ij=0
           do j=1,grd%lon2
              do i=1,grd%lat2
                 ij=ij+1
                 g_ql(i,j,lev)=sub(ij,Pe)
              end do
           end do

        case (8) ! Ice mixing ratio prognostic
        
           ij=0
           do j=1,grd%lon2
              do i=1,grd%lat2
                 ij=ij+1
                 g_qi(i,j,lev)=sub(ij,Pe)
              end do
           end do

     end select
  enddo


  return
end subroutine GenReload

subroutine GenReload_uv( UGrdIn, VGrdIn, wlev, icount, UGrdOut, VGrdOut)

! !USES:

  use kinds, only: r_kind,i_kind
  use mpimod, only: npe,mpi_comm_world,ierror,mpi_rtype
  use gridmod, only: grd => grd_a
  implicit none

! !INPUT PARAMETERS:

  real(r_kind),    intent(in) :: UGrdIn(grd%itotsub)
  real(r_kind),    intent(in) :: VGrdIn(grd%itotsub)
  integer(i_kind), intent(in) :: wlev(npe)
  integer(i_kind), intent(in) :: icount

! !OUTPUT PARAMETERS:

  real(r_kind),intent(out) :: UGrdOut(:,:,:)
  real(r_kind),intent(out) :: VGrdOut(:,:,:)


! !DESCRIPTION: Transfer contents of 2-d array global to 3-d subdomain array
!
! !REVISION HISTORY:
!   2004-05-14  treadon
!   2004-07-15  todling, protex-compliant prologue
!   2014-12-03  derber     - introduce vdflag and optimize routines
!   2016-05-31  de Mattos  - adapt to BAM model
!
! !REMARKS:
!
!   language: f90
!   machine:  ibm rs/6000 sp; sgi origin 2000; compaq/hp
!
! !AUTHOR:
!   treadon          org: np23                date: 2004-05-14
!
!EOP
!-------------------------------------------------------------------------
!BOC
      character(len=64), parameter :: myname_=trim(myname)//' :: GenReload_UV( )'

  integer(i_kind) :: i,j,k,ij
  integer(i_kind) :: Pe
  integer(i_kind) :: var
  integer(i_kind) :: klev
  real(r_kind),allocatable :: subU(:,:)
  real(r_kind),allocatable :: subV(:,:)


#ifdef DEBUG
    WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif


  allocate (subU(grd%lat2*grd%lon2,npe))
  allocate (subV(grd%lat2*grd%lon2,npe))

  call mpi_alltoallv(UGrdIn,        & 
                     grd%ijn_s,     &
                     grd%displs_s,  &
                     mpi_rtype,     &
                     subU,          &
                     grd%irc_s,     &
                     grd%ird_s,     &
                     mpi_rtype,     &
                     mpi_comm_world,&
                     ierror)

  call mpi_alltoallv(VGrdIn,        & 
                     grd%ijn_s,     &
                     grd%displs_s,  &
                     mpi_rtype,     &
                     subV,          &
                     grd%irc_s,     &
                     grd%ird_s,     &
                     mpi_rtype,     &
                     mpi_comm_world,&
                     ierror)


!$omp parallel do  schedule(dynamic,1) private(k,i,j,ij,klev)

     ! each Pe has only one field. so kmax .eq. used Pe to read
     do pe = 1, icount

        k  = wlev(pe)

        ij = 0
        do j=1,grd%lon2
           do i=1,grd%lat2
              ij = ij+1
              UGrdOut(i,j,k) = subU(ij,pe)
           end do
        end do
        
        ij = 0
        do j=1,grd%lon2
           do i=1,grd%lat2
              ij = ij+1
              VGrdOut(i,j,k) = subV(ij,pe)
           end do
        end do

     enddo

!$omp end parallel do



  deallocate (subU)
  deallocate (subV)

  return
end subroutine GenReload_uv


 subroutine GenFill_ns(GridIn,GridOut)

! !USES:

   use MiscMod, only: r8
   use kinds, only: r_kind,i_kind
   use constants, only: zero,one
   use general_sub2grid_mod, only: sub2grid_info
   use gridmod, only: grd => grd_a

   implicit none

! !INPUT PARAMETERS:

   real(r8),            intent(in   ) :: GridIn (:,:) ! input grid <imax,jmax>(nlon,nlat-2)
   real(r_kind),        intent(  out) :: GridOut( : )! output grid <itotsub>

! !DESCRIPTION: This routine adds a southern and northern latitude
!               row to the input grid.  The southern row contains
!               the longitudinal mean of the adjacent latitude row.
!               The northern row contains the longitudinal mean of
!               the adjacent northern row.
!
!               The added rows correpsond to the south and north poles.
!
!               In addition to adding latitude rows corresponding to the
!               south and north poles, the routine reorder the output
!               array so that it is a one-dimensional array read in
!               an order consisten with that assumed for total domain
!               gsi grids.
!
!               The assumed order for the input grid is longitude as
!               the first dimension with array index increasing from
!               east to west.  The second dimension is latitude with
!               the index increasing from north to south.  This ordering
!               differs from that used in the GSI.
!
!               The GSI ordering is latitude first with the index
!               increasing from south to north.  The second dimension is
!               longitude with the index increasing from east to west.
!
!               Thus, the code below also rearranges the indexing and
!               order of the dimensions to make the output grid
!               consistent with that which is expected in the rest of
!               gsi.
!
!
! !REVISION HISTORY:
!   2004-08-27  treadon
!
! !REMARKS:
!   language: f90
!   machine:  ibm rs/6000
!
! !AUTHOR:
!   treadon          org: np23                date: 2004-08-27
!
!EOP
!-------------------------------------------------------------------------
!BOC
      character(len=64), parameter :: myname_=trim(myname)//' :: GenFill_ns( )'


!  Declare local variables
   integer(i_kind) i,j,k,nlatm2
   real(r_kind) rnlon,sumn,sums

#ifdef DEBUG
    WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif

!  Compute mean along southern and northern latitudes

   sumn   = zero
   sums   = zero
   nlatm2 = grd%nlat-2

   do i = 1, grd%nlon

      sumn = sumn + real(GridIn(i,1),r_kind)
      sums = sums + real(GridIn(i,nlatm2),r_kind)

   end do

   rnlon = one/float(grd%nlon)
   sumn  = sumn*rnlon
   sums  = sums*rnlon

!  Transfer local work array to output grid
   do k = 1, grd%itotsub

      j = grd%nlat-grd%ltosi_s(k)

      if(j .eq. grd%nlat-1) then

         GridOut(k) = sums

      else if(j .eq. 0) then

         GridOut(k) = sumn

      else

        i          = grd%ltosj_s(k)
        GridOut(k) = real(GridIn(i,j),r_kind)

      end if

   end do

   return
 end subroutine GenFill_ns

 subroutine GenFillUV_ns(UGrdIn, VGrdIn, UGrdOut, VGrdOut)

! !USES:

   use sigioBAMMod, only: BAMFile
   use MiscMod, only: r4, r8, rd
   use MiscMod, only: GetLongitudes, GetGaussianLatitudes

   use kinds, only: r_kind, i_kind
   use constants, only: zero
   use gridmod, only : grd=>grd_a

   implicit none

! !INPUT PARAMETERS:

   real(r8),            intent(in   ) :: UGrdIn(:,:)
   real(r8),            intent(in   ) :: VGrdIn(:,:)
   real(r_kind),        intent(  out) :: UGrdOut(:)
   real(r_kind),        intent(  out) :: VGrdOut(:)

! !DESCRIPTION: This routine adds a southern and northern latitude
!               row to the input grid.  The southern row contains
!               the longitudinal mean of the adjacent latitude row.
!               The northern row contains the longitudinal mean of
!               the adjacent northern row.
!
!               The added rows correpsond to the south and north poles.
!
!               In addition to adding latitude rows corresponding to the
!               south and north poles, the routine reorder the output
!               array so that it is a one-dimensional array read in
!               an order consisten with that assumed for total domain
!               gsi grids.
!
!               The assumed order for the input grid is longitude as
!               the first dimension with array index increasing from
!               east to west.  The second dimension is latitude with
!               the index increasing from north to south.  This ordering
!               differs from that used in the GSI.
!
!               The GSI ordering is latitude first with the index
!               increasing from south to north.  The second dimension is
!               longitude with the index increasing from east to west.
!
!               Thus, the code below also rearranges the indexing and
!               order of the dimensions to make the output grid
!               consistent with that which is expected in the rest of
!               gsi.
!
!
! !REVISION HISTORY:
!   2004-08-27  treadon
!
! !REMARKS:
!   language: f90
!   machine:  ibm rs/6000
!
! !AUTHOR:
!   treadon          org: np23                date: 2004-08-27
!
!EOP
!-------------------------------------------------------------------------
!  Declare local variables
   integer(i_kind) i,j,k,nlatm2
   real(r_kind) polnu, polnv, polsu, polsv
   real(r8), pointer :: rlon(:)
   real(r4), pointer :: slon(:)
   real(r4), pointer :: clon(:)

   character(len=64), parameter :: myname_=trim(myname)//' :: GenFillUV_ns( )'


#ifdef DEBUG
    WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif

   !
   ! Get BAM grid Informations
   !

   allocate(rlon(grd%nlon))
   allocate(slon(grd%nlon))
   allocate(clon(grd%nlon))

   call GetLongitudes(grd%nlon, 0.0_r8, rlon)

   do i=1,grd%nlon
      slon(i) = sin(rlon(i)/rd)
      clon(i) = cos(rlon(i)/rd)
   enddo

!  Compute mean along southern and northern latitudes
   polnu=zero
   polnv=zero
   polsu=zero
   polsv=zero
   nlatm2=grd%nlat-2
   do i=1,grd%nlon
      polnu=polnu+UGrdIn(i,1     )*clon(i)-VGrdIn(i,1     )*slon(i)
      polnv=polnv+UGrdIn(i,1     )*slon(i)+VGrdIn(i,1     )*clon(i)
      polsu=polsu+UGrdIn(i,nlatm2)*clon(i)+VGrdIn(i,nlatm2)*slon(i)
      polsv=polsv+UGrdIn(i,nlatm2)*slon(i)-VGrdIn(i,nlatm2)*clon(i)
   end do
   polnu=polnu/float(grd%nlon)
   polnv=polnv/float(grd%nlon)
   polsu=polsu/float(grd%nlon)
   polsv=polsv/float(grd%nlon)

!  Transfer local work array to output grid
   do k=1,grd%itotsub
      j=grd%nlat-grd%ltosi_s(k)
      i=grd%ltosj_s(k)
      if(j == grd%nlat-1)then
        UGrdOut(k) = polsu*clon(i)+polsv*slon(i)
        VGrdOut(k) = polsu*slon(i)-polsv*clon(i)
      else if(j == 0) then
        UGrdOut(k) = polnu*clon(i)+polnv*slon(i)
        VGrdOut(k) = -polnu*slon(i)+polnv*clon(i)
      else
        UGrdOut(k)=UGrdIn(i,j)
        VGrdOut(k)=VGrdIn(i,j)
      end if
   end do

   deallocate(rlon, slon, clon)

   return
 end subroutine GenFillUV_ns
!EOC
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !FUNCTION: GenReadBAMSFC
!
! !DESCRIPTON: subroutine to get BAM sfc fields and tranfer to GSI
!
!             
!                 
!\\
!\\
! !INTERFACE:
!
   subroutine Read_BAMSFC(   iope,  mype,  &
                             fact10,       &! 10-meter wind factor
                             td0,          &! surface temperature [k]
                             sheleg,       &! snow depth [m]
                             vtype,        &! vegetation type  [1-13]
                             vcover,       &! vegetation cover [-]
                             stype,        &! soil type [1-9]
                             tg0,          &! surface soil temperature [k]
                             w0,           &! surface soil moisture [fraction]
                             lsimsk,       &! land sea mask
                             Z0,           &! surface roughness length [cm]
                             topo,         &! Orography [ m ]                   
                             use_sfc_any   &
                             )  
!
!  !USES:
!
      ! GSI guess number actual count of sfc in-cache time slots 
      use guess_grids, only: nfldsfc

      ! GSI grid - SFC grid sizes
      use gridmod, only: nlat_sfc, nlon_sfc

      ! GSI kinds
      use kinds, only: r_kind,i_kind

      ! MPI GSI 
      use mpimod, only: mpi_itype,mpi_rtype,mpi_comm_world

      implicit none

!
! !INPUT PARAMETERS:
!
    integer(i_kind), intent(in   ) :: iope        ! mpi task handling i/o
    integer(i_kind), intent(in   ) :: mype        ! mpi task id
    logical,         intent(in   ) :: use_sfc_any !

!
! !OUTPUT PARAMETERS:
!

    integer(i_kind), dimension(:,:)  , intent(inout) :: lsimsk ! land sea mask 
    real(r_kind)   , dimension(:,:)  , intent(inout) :: vtype  ! vegetation type  [1-13]
    real(r_kind)   , dimension(:,:)  , intent(inout) :: stype  ! soil type [1-9]
    real(r_kind)   , dimension(:,:)  , intent(inout) :: topo   ! Orography [ m ]
    real(r_kind)   , dimension(:,:,:), intent(inout) :: fact10 ! 10-meter wind factor
    real(r_kind)   , dimension(:,:,:), intent(inout) :: td0    ! surface temperature [k]
    real(r_kind)   , dimension(:,:,:), intent(inout) :: sheleg ! snow depth [m]
    real(r_kind)   , dimension(:,:,:), intent(inout) :: vcover ! vegetation cover [-]
    real(r_kind)   , dimension(:,:,:), intent(inout) :: tg0    ! surface soil temperature [k]
    real(r_kind)   , dimension(:,:,:), intent(inout) :: w0     ! surface soil moisture [fraction]
    real(r_kind)   , dimension(:,:,:), intent(inout) :: Z0     ! surface roughness length [cm]   


! !REVISION HISTORY:
!
!   06 Jul 2016 - J. G. de Mattos -  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC

      character(len=100), parameter :: myname_=trim(myname)//':: Read_BAMSFC( )'

      !
      ! Auxiliary vars
      !

      integer :: istat
      integer :: it
      integer :: npts
      integer :: nptsall


#ifdef DEBUG
      WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif

      npts    = nlat_sfc * nlon_sfc
      nptsall = npts * nfldsfc

      !
      ! Loop over all time sfc gues fields
      !
      do it = 1, nfldsfc

         !
         ! Read BAM fields
         !

         call GenReadBAMSFC  (   iope,  mype,     &
                                 lsimsk (:,:),    &
                                 vtype  (:,:),    &
                                 stype  (:,:),    &
                                 topo   (:,:),    &
                                 fact10 (:,:,it), &
                                 td0    (:,:,it), &
                                 sheleg (:,:,it), &
                                 vcover (:,:,it), &
                                 tg0    (:,:,it), &
                                 w0     (:,:,it), &
                                 Z0     (:,:,it), &
                                 use_sfc_any,     &
                                 it               &
                             )
      enddo

   end subroutine

!EOC
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !FUNCTION: GenReadBAMSFC
!
! !DESCRIPTON: subroutine to read sfc files from BAM fct files
!
!             
!                 
!\\
!\\
! !INTERFACE:
!
   subroutine GenReadBAMSFC(   iope,  mype,               &
                               lsimsk, vtype, stype, topo,& 
                               fact10,td0, sheleg,vcover, &
                               tg0, w0, Z0,               &
                               use_sfc_any, it )
       ! Model read 
       use sigioBAMMod, only : BAMFile
       use MiscMod, only: r4, r8

       ! GSI kinds
       use kinds, only: r_kind,i_kind

       ! GSI guess files names
       use guess_grids, only: ifilesfc

      ! GSI grid informations
      use gridmod, only : grd=>grd_a, nlat_sfc, nlon_sfc

      ! MPI GSI 
      use mpimod, only: mpi_itype,mpi_rtype,mpi_comm_world

       implicit none

!
! !INPUT PARAMETERS:
!
    integer(i_kind), intent(in   ) :: iope        ! mpi task handling i/o
    integer(i_kind), intent(in   ) :: mype        ! mpi task id
    logical,         intent(in   ) :: use_sfc_any !
    integer(i_kind), intent(in   ) :: it
!
! !OUTPUT PARAMETERS:
!

    integer(i_kind), dimension(:,:), intent(inout) :: lsimsk ! land sea mask 
    real(r_kind)   , dimension(:,:), intent(inout) :: vtype  ! vegetation type  [1-13]
    real(r_kind)   , dimension(:,:), intent(inout) :: stype  ! soil type [1-9]
    real(r_kind)   , dimension(:,:), intent(inout) :: topo   ! Orography [ m ]
    real(r_kind)   , dimension(:,:), intent(inout) :: fact10 ! 10-meter wind factor
    real(r_kind)   , dimension(:,:), intent(inout) :: td0    ! surface temperature [k]
    real(r_kind)   , dimension(:,:), intent(inout) :: sheleg ! snow depth [m]
    real(r_kind)   , dimension(:,:), intent(inout) :: vcover ! vegetation cover [-]
    real(r_kind)   , dimension(:,:), intent(inout) :: tg0    ! surface soil temperature [k]
    real(r_kind)   , dimension(:,:), intent(inout) :: w0     ! surface soil moisture [fraction]
    real(r_kind)   , dimension(:,:), intent(inout) :: Z0     ! surface roughness length [cm]

!
! !PARAMTERS:
!

       integer, parameter                  :: NVars = 11
       character(len=40), dimension(NVars) :: VName = [                             &
                                                       'TOPOGRAPHY                ',& ! 01 
                                                       'LAND SEA ICE MASK         ',& ! 02
                                                       'ROUGHNESS LENGTH          ',& ! 03
                                                       'SURFACE TEMPERATURE       ',& ! 04
                                                       'SNOW DEPTH                ',& ! 05
                                                       'SOIL WETNESS OF SURFACE   ',& ! 06
                                                       'SURFACE SOIL TEMPERATURE  ',& ! 07
                                                       'VEGETATION COVER          ',& ! 08
                                                       'MASK VEGETATION           ',& ! 09
                                                       'MASK SOIL TEXTURE CLASSES ',& ! 10
                                                       '10-meter WIND FACTOR      ' & ! 11
                                                      ]

! !REVISION HISTORY:
!
!   06 Jul 2016 - J. G. de Mattos -  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      character(len=100), parameter :: myname_=trim(myname)//' :: GenReadSFCBAM( )'

      type(BAMFile) :: BAM  ! BAM files data type

      real(r8), allocatable :: OutGrid(:,:)
      real(r8), allocatable :: TmpGrid(:,:)

      integer :: imax
      integer :: jmax
      integer :: kmax
      integer :: Mend
      integer :: MnWv2
      integer :: MnWv3

      integer :: ivar
      integer :: istat
      integer :: iret
      integer :: j

      integer :: npts

      !
      !  variables in spectral space
      !
  
      real(r8), dimension(:),   allocatable :: divq
      real(r8), dimension(:),   allocatable :: vorq
      real(r8), dimension(:),   allocatable :: uveq
      real(r8), dimension(:),   allocatable :: vveq
  
      !
      !  variables in physical space
      !
  
      real(r8), dimension(:,:), allocatable :: gu
      real(r8), dimension(:,:), allocatable :: gv
      real(r8), dimension(:,:), allocatable :: u10m
      real(r8), dimension(:,:), allocatable :: v10m

      !
      ! Auxiliary Variables
      !

      real(r4), dimension( : ), allocatable :: clat
      real(r8), dimension(:,:), allocatable :: speed
      real(r8), dimension(:,:), allocatable :: speed10m

      character(len=80) :: fileFct
      character(len=80) :: fileDir

#ifdef DEBUG
      WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif

      !
      ! Define files names
      ! fct => sequential binary ieee file (spectral and grid points field)
      ! dir => ascii file with field information
      !

      !
      ! Files to read and get first-guess
      ! !! hardwire files names! not a best way !!
      !

      write(fileFct,'("BAM.fct.",I2.2)') ifilesfc(it)
      write(FileDir,'("BAM.dir.",I2.2)') ifilesfc(it)

      !
      ! Open BAM File
      !
      call BAM%Open(trim(fileDir), trim(fileFct), istat=istat)

      if (istat.ne.0)then

         write(6,*)trim(myname_),'Problem to open/read BAM files at ',MyPe, 'rank!'

         return
      endif

      !
      ! consistency check
      !
      
      call BAM%GetDims(imax, jmax, kmax, Mend)
      npts = nlat_sfc*nlon_sfc

      if(MyPe.eq.0)then
         write (6,'(3(1x,A))') trim(myname_),':: *** BAM_SFC *** reading', trim(BAM%fBin)
         write (6,'(2(1x,A,1x,I4))')'<IMax>',IMax,'---',grd%nlon
         write (6,'(2(1x,A,1x,I4))')'<JMax>',JMax,'---',grd%nlat-2
         write (6,'(2(1x,A,1x,I4))')'<KMax>',KMax,'---.',grd%nsig
      endif

      do ivar = 1, NVars
         if (myPe .ne. (ivar-1) )cycle
         if (                    &
             ( ivar .eq.  6 .or. & ! SOIL WETNESS OF SURFACE
               ivar .eq.  7 .or. & ! SURFACE SOIL TEMPERATURE
               ivar .eq.  8 .or. & ! VEGETATION COVER
               ivar .eq.  9 .or. & ! MASK VEGETATION
               ivar .eq. 10      & ! MASK SOIL TEXTURE CLASSES
             ) .and. .not. use_sfc_any &
            ) cycle

         allocate ( OutGrid( IMax, JMax ) )

         if ( ivar .eq. 11 )then
            call BAM%GetDims(imax, jmax, kmax, Mend)

            allocate(gu(imax,jmax))
            allocate(gv(imax,jmax))
            allocate(u10m(imax,jmax))
            allocate(v10m(imax,jmax))

            iret = 0
            call BAM%GetField('ZONAL WIND AT 10-M FROM SURFACE', 1, u10m, istat=istat)
            iret = iret + istat
            call BAM%GetField('MERID WIND AT 10-M FROM SURFACE', 1, v10m, istat=istat)
            iret = iret + istat

            call BAM%getUV(1,gu,gv)

            if ( iret .eq. 0 )then

               allocate(clat(JMax))

               call BAM%GetWCoord('clat', clat, iret)

               do j=1,jmax
                  gu(1:imax,j) = gu(1:imax,j)/clat(j)
                  gv(1:imax,j) = gv(1:imax,j)/clat(j)
               enddo

               deallocate(clat)

               allocate (speed(imax,jmax))

               speed = sqrt(gu*gu + gv*gv)

               deallocate (gu)
               deallocate (gv)

               allocate(speed10m(imax,jmax))

               speed10m = sqrt(u10m*u10m + v10m*v10m)

               OutGrid = speed10m/speed

               deallocate(speed)
               deallocate(speed10m)

            else
               deallocate(gu)
               deallocate(gv)
               deallocate(u10m)
               deallocate(v10m) 

            endif


         else
            call BAM%GetField(trim(VName(ivar)), 1, OutGrid, istat=istat)
         endif

         select case ( ivar )
            case ( 1) ! Topography [ m ]
               if (istat .eq. 0 )then
                  call SFCTrans(OutGrid,topo)
                else
                  topo = Udef
               endif


            case ( 2) ! Land Sea Ice Mask [L=1,S=0,I=2]

               if (istat .eq. 0 )then
                  allocate(TmpGrid(JMax+2,IMax))
                  call SFCTrans(OutGrid,TmpGrid)

                  lsimsk = nint(abs(TmpGrid))

                  deallocate(TmpGrid)
               else
                  lsimsk = iUdef
               endif


            case ( 3) ! Roughness Length [ cm ]
               
               if (istat .eq. 0 )then
                  !convert from m to cm
                  OutGrid = OutGrid * 100.0_r8

                  call SFCTrans(OutGrid,Z0)
               else
                  Z0 = Udef
               endif


            case ( 4) ! Surface Temperature [ k ]

               if (istat .eq. 0 )then
                  OutGrid = abs(OutGrid)
                  call SFCTrans(OutGrid,td0)
               else
                  td0 = Udef
               endif

            case ( 5) ! Snow Depth [ m ]

               if (istat .eq. 0 )then
                  !convert from [mm] to [m]
                  OutGrid = OutGrid * 0.001_r8

                  call SFCTrans(OutGrid,sheleg)
               else
                  sheleg = Udef
               endif

            case ( 6) ! Soil Wetness of Surface [ fraction ]

               if (istat .eq. 0 )then
                  call SFCTrans(OutGrid,w0)
               else
                  w0 = Udef
               endif

            case ( 7) ! Surface Soil Temperature [ k ]

               if (istat .eq. 0 )then
                  call SFCTrans(OutGrid,tg0)
               else
                  tg0 = Udef
               endif


            case ( 8) ! Vegetation Cover [ - ]

               if (istat .eq. 0 )then
                  call SFCTrans(OutGrid,vcover)
               else
                  vcover = Udef
               endif

!                                          SSIB    :    IBIS
            case ( 9) ! Vegetation Type [ 1 - 13 ] : [ 0 - 15 ]

               if (istat .eq. 0 )then
                  call SFCTrans(OutGrid,vtype)
               else
                  vtype = Udef
               endif


            case (10) ! Soil Type

               if (istat .eq. 0 )then
                  call SFCTrans(OutGrid,stype)
               else
                  stype = Udef
               endif


            case (11) ! 10-meter Wind Factor [ - ]

               if (istat .eq. 0 )then
                  call SFCTrans(OutGrid,fact10)
               else
                  fact10 = 0.5!Udef
               endif


         end select

         deallocate(OutGrid)

      enddo

     call mpi_bcast(  topo, npts, mpi_rtype,  0,  mpi_comm_world, istat)
     call mpi_bcast(lsimsk, npts, mpi_itype,  1,  mpi_comm_world, istat)
     call mpi_bcast(    Z0, npts, mpi_rtype,  2,  mpi_comm_world, istat)
     call mpi_bcast(   td0, npts, mpi_rtype,  3,  mpi_comm_world, istat)
     call mpi_bcast(sheleg, npts, mpi_rtype,  4,  mpi_comm_world, istat)
     if (use_sfc_any)then
        call mpi_bcast(    w0, npts, mpi_rtype,  5,  mpi_comm_world, istat)
        call mpi_bcast(   tg0, npts, mpi_rtype,  6,  mpi_comm_world, istat)
        call mpi_bcast(vcover, npts, mpi_rtype,  7,  mpi_comm_world, istat)
        call mpi_bcast (vtype, npts, mpi_rtype,  8,  mpi_comm_world, istat)
        call mpi_bcast( stype, npts, mpi_rtype,  9,  mpi_comm_world, istat)
        call mpi_bcast(fact10, npts, mpi_rtype, 10,  mpi_comm_world, istat)
     endif


      !
      ! Close BAM Files
      !

      call BAM%Close( )

      return
   end subroutine

!EOC
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !FUNCTION: LUAvail
!
! !DESCRIPTON: function to return next available logical unit
!
!
!             
!                 
!\\
!\\
! !INTERFACE:
!
   subroutine SFCTrans( GridIn, GridOut )
! !USES:
      use MiscMod, only: r8
      use kinds, only: r_kind, i_kind
      use constants, only: zero,one

      implicit none
! !IMPUT PARAMETERS:

      real(r8),            intent(in   ) :: GridIn (:,:) ! input grid  ->  IMax,JMax
      real(r_kind),        intent(inout) :: GridOut(:,:) ! output grid ->  JMax,IMax

! !DESCRIPTION: This routine adds a southern and northern latitude
!               row to the input grid.  The southern row contains
!               the longitudinal mean of the adjacent latitude row.
!               The northern row contains the longitudinal mean of
!               the adjacent northern row.
!
!               The added rows correpsond to the south and north poles.
!
!               In addition to adding latitude rows corresponding to the
!               south and north poles, the routine reorder the output
!               array so that it is a one-dimensional array read in
!               an order consisten with that assumed for total domain
!               gsi grids.
!
!               The assumed order for the input grid is longitude as
!               the first dimension with array index increasing from
!               east to west.  The second dimension is latitude with
!               the index increasing from north to south.  This ordering
!               differs from that used in the GSI.
!
!               The GSI ordering is latitude first with the index
!               increasing from south to north.  The second dimension is
!               longitude with the index increasing from east to west.
!
!               Thus, the code below also rearranges the indexing and
!               order of the dimensions to make the output grid
!               consistent with that which is expected in the rest of
!               gsi.
!
!
! !REVISION HISTORY:
!   06 Jul 2016 - J. G. de Mattos -  Initial code.
!
!
!EOP
!-------------------------------------------------------------------------
!BOC
      integer         :: IMax, JMax
      integer(i_kind) :: i, j, k, nlatm2
      real(r_kind)    :: rnlon, sumn, sums

      character(len=64), parameter :: myname_=trim(myname)//' :: SFCTrans( )'

#ifdef DEBUG
    WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif
!  Compute mean along southern and northern latitudes


      IMax = size(GridIn,1) ! lon
      JMax = size(GridIn,2) ! lat

      sumn   = zero
      sums   = zero

      do i = 1, IMax

         sumn = sumn + real( GridIn ( i,    1 ), r_kind)
         sums = sums + real( GridIn ( i, JMax ), r_kind)

      end do

      rnlon = one/float(IMax)
      sumn  = sumn * rnlon
      sums  = sums * rnlon

!  Transfer local work array to output grid

      do j = 1, IMax

         GridOut(1,j) = sums

         do i = 2, JMax+1
            GridOut(i,j) = GridIn( j, JMax+2-i )
         end do

         GridOut(JMax+2, j) = sumn

      end do

   return


   end subroutine
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !ROUTINE: Write_BAM
!              
! !DESCRIPTON: rotine to write atmospheric analysis fields from GSI to BAM 
!
!             
!                 
!\\
!\\
! !INTERFACE:
!

   subroutine write_bam(increment, mype, mype_atm, mype_sfc )
!
! !USES:
!
      ! Model read 

      ! GSI kinds
      use kinds, only: r_kind, i_kind

      !GSI bundle
      use gsi_metguess_mod, only: gsi_metguess_bundle
      use gsi_bundlemod, only: gsi_bundlegetpointer
      
      !GSI Guess informations
      use guess_grids, only: nfldsig 
      use guess_grids, only: ntguessig,ntguessfc,ifilesig,nfldsig

      use gsi_4dvar, only: lwrite4danl

      implicit none
!
! !INPUT PARAMETERS:
!
      integer(i_kind), intent(in   ) :: increment
      integer(i_kind), intent(in   ) :: mype
      integer(i_kind), intent(in   ) :: mype_atm
      integer(i_kind), intent(in   ) :: mype_sfc

! !REVISION HISTORY:
!
!   06 Jul 2016 - J. G. de Mattos -  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      character(len=64), parameter :: myname_=trim(myname)//' :: write_bam( )'


      !
      ! GSI bungle guess fields
      !

      real(r_kind), pointer, dimension(:,:  ) :: anl_z_it    => NULL()
      real(r_kind), pointer, dimension(:,:  ) :: anl_ps_it   => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: anl_u_it    => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: anl_v_it    => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: anl_vor_it  => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: anl_div_it  => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: anl_tv_it   => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: anl_q_it    => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: anl_oz_it   => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: anl_cw_it   => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: anl_ql_it   => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: anl_qi_it   => NULL()
  
      !
      ! Auxiliary vars
      !

      integer :: istatus
      integer :: ntlevs
      integer :: it
      integer :: itout
      integer :: i

      character(len=80) :: fileAnl


#ifdef DEBUG
      WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif

      !
      ! Write atmospheric analysis file
      ! Loop over all time analysis fields
      !

      if (.not.lwrite4danl) then
         ntlevs = 1
      else
         ntlevs = nfldsig
      end if
      


      do it = 1, ntlevs

         !
         ! define anl file name
         !

         if (increment>0) then
            fileAnl = 'BAM.inc'
            itout   = increment
            if(mype.eq.0) write(6,*) 'WRITE_BAM: writing time slot ', itout
         else if (.not.lwrite4danl) then
            fileAnl = 'BAM.anl'
            itout    = ntguessig
            if(mype.eq.0) write(6,*) 'WRITE_BAM: writing single analysis state for F ', itout
         else
            write(fileAnl,'("BAM.anl.",I2.2)') ifilesig(it)
            itout = it
            if(mype.eq.0) write(6,*) 'WRITE_BAM: writing full analysis state for F ', itout
         endif


         !
         ! Allocating anl transfer fields
         !

         call AllocateAtmFields( )

         !
         ! Get fields from GSI Bundle
         !

         ! Orography
         call gsi_bundlegetpointer (gsi_metguess_bundle(itout), 'z', anl_z_it, istatus)
         if(istatus.eq.0) g_z = anl_z_it

         ! Surface Pressure
         call gsi_bundlegetpointer (gsi_metguess_bundle(itout), 'ps', anl_ps_it, istatus)
         if(istatus.eq.0) g_ps = anl_ps_it

         ! Zonal wind component
         call gsi_bundlegetpointer (gsi_metguess_bundle(itout), 'u', anl_u_it, istatus)
         if(istatus.eq.0) g_u = anl_u_it

         ! Meridional wind component
         call gsi_bundlegetpointer (gsi_metguess_bundle(itout), 'v', anl_v_it, istatus)
         if(istatus.eq.0) g_v = anl_v_it
      
         ! Vorticity
         call gsi_bundlegetpointer (gsi_metguess_bundle(itout), 'vor', anl_vor_it, istatus)
         if(istatus.eq.0) g_vor = anl_vor_it

         ! Divergency
         call gsi_bundlegetpointer (gsi_metguess_bundle(itout), 'div', anl_div_it, istatus)
         if(istatus.eq.0) g_div = anl_div_it

         ! Virtual Temperature
         call gsi_bundlegetpointer (gsi_metguess_bundle(itout), 'tv', anl_tv_it, istatus)
         if(istatus.eq.0) g_tv = anl_tv_it

         ! Specific Humidity
         call gsi_bundlegetpointer (gsi_metguess_bundle(itout), 'q', anl_q_it, istatus)
         if(istatus.eq.0) g_q = anl_q_it

         ! Ozone
         call gsi_bundlegetpointer (gsi_metguess_bundle(itout), 'oz', anl_oz_it, istatus)
         if(istatus.eq.0) g_oz = anl_oz_it

         ! Total Cloud Water Content
         call gsi_bundlegetpointer (gsi_metguess_bundle(itout), 'cw', anl_cw_it, istatus)
         if(istatus.eq.0) g_cw = anl_cw_it

         ! Liq mixing ratio prognostic

         call gsi_bundlegetpointer (gsi_metguess_bundle(itout), 'ql',anl_ql_it, istatus)
         if(istatus.eq.0) g_ql = anl_ql_it

         ! Ice mixing ratio prognostic

         call gsi_bundlegetpointer (gsi_metguess_bundle(itout), 'qi', anl_qi_it, istatus)
         if(istatus.eq.0) g_qi = anl_qi_it
         
         !
         ! Write BAM fields
         !

         call GenWriteBAM (fileAnl, mype, mype_atm, itout )


         !
         ! Clear anl transfer fields
         !

         call DeallocateAtmFields ( )

     enddo

  end subroutine
!EOC
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !ROUTINE: GenWriteBAM
!
! !DESCRIPTON: rotine to read atmospheric guess fields
!
!             
!                 
!\\
!\\
! !INTERFACE:
!

   subroutine GenWriteBAM( fileAnl, MyPe, MyPe_Out, it )
!
! !USES:
!
      ! Model read 
      use sigioBAMMod, only : BAMFile
      use sigioBAMMod, only : BAM_WriteField, BAM_SendField
      use MiscMod, only: i4, i8, r4, r8
      
      ! GSI grid informations
      use gridmod, only : nlon, nlat, lat2, lon2, nsig, jcap, jcap_b
      use gridmod, only : grd=>grd_a

      ! GSI file guess information
      use guess_grids, only: ntguessig, ifilesig

      ! GSI analysis time
      use obsmod, only: iadate

      ! GSI 4dvar case
      use gsi_4dvar, only: ibdate, nhr_obsbin, lwrite4danl

      !GSI kinds
      use kinds, only : r_kind, i_kind

      !GSI constants
      use constants, only:zero 

      !MPI
      use mpimod, only : NPe, mpi_comm_world, mpi_itype


      implicit none
!
! !INPUT PARAMETERS:
!
      character(len=*), intent(in   ) :: fileAnl  ! name of anlFile
      integer,          intent(in   ) :: MyPe     ! current mpi task id
      integer,          intent(in   ) :: MyPe_Out ! mpi task id to write
      integer,          intent(in   ) :: it

!
! !PARAMTERS:
!

      integer, parameter                  :: NVars = 6
      character(len=40), dimension(NVars) :: VName = [                            &
                                                     'TOPOGRAPHY                ',& !1
                                                     'LN SURFACE PRESSURE       ',& !2
                                                     'VIRTUAL TEMPERATURE       ',& !3
                                                     'DIVERGENCE                ',& !4
                                                     'VORTICITY                 ',& !5
                                                     'SPECIFIC HUMIDITY         ' & !6
!                                                     'CLOUD LIQUID WATER CONTENT',& !7
!                                                     'CLOUD ICE WATER CONTENT   ' & !8
                                                     ] 

! !REVISION HISTORY:
!
!   28 Apr 2016 - J. G. de Mattos -  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC      
      character(len=100), parameter :: myname_=':: GenWriteBAM()'

      type(BAMFile) :: bam, bamOut  ! BAM files data type

      !
      ! Variables to write GFCT spectral file
      !

      integer (i4)                          :: ifday, rc, nymd, nhms
      real    (r4)                          :: tod
      integer (r4), dimension(4)            :: idate, idatec

      real(r8),     dimension(:),   allocatable :: workq ! Spectral Space
      real(r_kind), dimension(:),   allocatable :: work1 ! Physical Space
      real(r8),     dimension(:,:), allocatable :: work2 ! Physical Space
      real(r8),     dimension(:,:), allocatable :: work3 !
      real(r8),     dimension(:,:), allocatable :: work4 !
      !
      !  variables in physical space
      !
  

      !
      ! Auxiliary Variables
      !

      integer(i_kind),dimension(5) :: mydate
      integer(i_kind),dimension(8) :: ida,jda
      real(r_kind),   dimension(5) :: fha
      real(r_kind) :: psaux

      integer :: istat

      integer :: WVar(NPe)
      integer :: WLev(NPe)


      integer :: Nflds
      integer :: imax, imax_b
      integer :: jmax, jmax_b
      integer :: kmax
      integer :: Mend
      integer :: mnwv2
      integer, allocatable :: nlevs(:)

      integer :: i, j, k, ij
      integer :: icount
      integer :: ivar
      integer :: iret
      integer :: ilev
      integer :: iPe
      integer :: WrkPe

      integer :: iret_write

      character(len=80) :: fileFct
      character(len=80) :: fileDir
      character(len=80) :: fAnlDir

#ifdef DEBUG
      WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif

      ! configure some necessary info

      Nflds = NVars * grd%nsig
      WrkPe = NPe-1

      !
      ! Set fct and dir guess file name
      !


      write(fAnlDir,'("BAM.dir.",I2.2)') ifilesig(ntguessig-1)

      !
      ! all tasks open dir files to get model info
      ! but only the output task will open BAM anl file

      if (MyPe .eq. MyPe_Out)then
         ! Only one task should open output file for write
         !
         !  * fileAnl was get from input argument
         !
         call bamOut%Open(trim(fAnlDir), trim(fileAnl), mode='w', ftype='anl', istat=istat)
         if (istat.ne.0)then
            write(6,*)trim(myname_),'Problem to open BAM anl file to write at ',MyPe, 'rank!', istat
            stop
         endif

      else

         call bamOut%Open(trim(fAnlDir), istat=istat)
         if (istat.ne.0)then
            write(6,*)trim(myname_),'Problem to open BAM dir file to get model info at ',MyPe, 'rank!', istat
            stop
         endif

      endif

      ! get vars info at all tasks
      allocate(NLevs(NVars))
      do ivar = 1, Nvars
         NLevs(ivar) = bamOut%GetNlevels(trim(VName(ivar)),iret)
         if(iret .ne. 0)then
            write(stdout,'(2A,1x,I4)')trim(myname_),':: *** Error ***: BAM file should be openned at', MyPe
            stop
         endif         
      enddo

      call bamOut%GetDims(imax_b, jmax_b, kmax, Mend)
      MnWv2 = (Mend+1) * (Mend+2)
      imax  = grd%nlon
      jmax  = grd%nlat-2

      !
      ! only output task will open BAM anl file
      !

      if(MyPe .eq. MyPe_Out)then

         write (stdout,'(3(1x,A))') trim(myname_),':: *** BAM *** writing', trim(bamOut%fBin)
         write (stdout,'(2(1x,A,1x,I4))')'<JCAP>',jcap_b,'<->',jcap
         write (stdout,'(2(1x,A,1x,I4))')'<IMax>',IMax_b,'<->',imax
         write (stdout,'(2(1x,A,1x,I4))')'<JMax>',JMax_b,'<->',jmax
         write (stdout,'(2(1x,A,1x,I4))')'<KMax>',KMax,'<->',grd%nsig

         ! Load date
         if (.not.lwrite4danl) then
           mydate = iadate
         else
         !  increment mydate ...
            mydate = ibdate
            fha(:) = zero ; ida=0; jda=0
            fha(2) = real(nhr_obsbin*(it-1))  ! relative time interval in hours
            ida(1) = mydate(1) ! year
            ida(2) = mydate(2) ! month
            ida(3) = mydate(3) ! day
            ida(4) = 0         ! time zone
            ida(5) = mydate(4) ! hour
   
      ! Move date-time forward by nhr_assimilation hours
            call w3movdat(fha,ida,jda)
            mydate(1) = jda(1)
            mydate(2) = jda(2)
            mydate(3) = jda(3)
            mydate(4) = jda(5)
         end if
   
         call bamOut%WriteAnlHeader(mydate, iret)

      endif
      

      !
      ! Some Ajustments
      !
!$omp parallel do  schedule(dynamic,1) private(i,j)      
      ! Surface Pressure is ln(pslc) in kiloPascal [kPa] or centibar [cb]
      do j=1, lon2
         do i=1,lat2
            g_ps(i,j) = log(g_ps(i,j))
         enddo
      enddo
!$omp end parallel do

      !
      !------------------------------------------!
      !
      wvar   = 0
      ilev   = 0
      iCount = 1

      do ivar = 1, NVars
         do ilev = 1, NLevs(ivar)

            iPe = iCount - 1

            !----------------------------------------
            ! This is used to scatter between all Pe's
            ! Como cada Pe ira pegar um campo
            ! em um nivel, estas duas variaveis dirao
            ! qual variavel é lida em cada Pe
            !
            wvar(iCount) = ivar ! What Variable
            wlev(iCount) = ilev ! What Level
            !----------------------------------------

            if ( icount .eq. NPe .or. ( ivar .eq. NVars .and. ilev .eq. NLevs(NVars)) )then

               !
               ! Transfer contents of 3-d subdomain array to 2-d array global
               !

               allocate(work1(IMax*(JMax+2)))

               call GenGatherBAM( g_z,   & ! Input  Topography
                                  g_ps,  & !        Surface Pressure
                                  g_tv,  & !        Virtural Temperature
                                  g_vor, & !        Vorticity
                                  g_div, & !        Divergency
                                  g_q,   & !        Specific Humidity
                                  g_ql,  & !        Liq mixing ratio prognostic
                                  g_qi,  & !        Ice mixing ratio prognostic
                                  WVar,  & !        Position in Var list
                                  WLev,  & !        Level
                                  icount,& !        MaxCount until here
                                  MyPe,  &
                                  VName, &
                                  work1  & ! OutPut Field read by Pe
                                )

               if (MyPe .lt. iCount)then

                  allocate(work2(IMax,JMax))

                  call load_grid(work1, work2)

                  deallocate(work1)

                  !------------------------------------------------------------------------!
                  !
                  ! Convert from grid to spec
                  !


                  if(jcap.ne.jcap_b)then

                     !
                     ! if jcap ne. jcap_b we need get original orography from guess
                     ! to avoid some short waves (need investigate better this)
                     !

                     if( wvar(MyPe+1) .eq. 1 )then ! orography

                        write(fileFct,'("BAM.fct.",I2.2)') ifilesig(ntguessig)
                        write(fileDir,'("BAM.dir.",I2.2)') ifilesig(ntguessig)

                        call bam%Open(trim(fileDir), trim(fileFct), istat=istat)
                        if (istat.ne.0)then
                           write(6,*)trim(myname_),'Problem to open BAM fct file to get orography at ',MyPe, 'rank!'
                           stop
                        endif

                        allocate(workq(mnwv2))

                        call bam%GetField(trim(VName(1)), 1, workq, istat=istat)

                        call bam%close( )

                     else

                        allocate(work3(IMax_b,JMax_b))

                        call lterp(work2, jcap, jcap_b, work3)

                        deallocate(work2)

                        allocate(workq(mnwv2))

                        call bamOut%Grid2Spec(work3, workq)

                        deallocate(work3)

                     endif

                  else

                     allocate(workq(mnwv2))

                     call bamOut%Grid2Spec(work2, workq)

                     deallocate(work2)
 
                  endif

                  if (MyPe .eq. MyPe_Out)then

                     call BAM_WriteField(bamOut%uBin, workq, MyPe_Out, icount)

                  else   ! send to MyPe_out/

                     call BAM_SendField(MyPe, MyPe_Out, workq)

                  endif

                  if(allocated(work2)) deallocate(work2)
                  if(allocated(workq)) deallocate(workq)

               endif
               !
               ! reset all counters
               !

               if(allocated(work1))deallocate(work1)
               
               WVar   = 0
               WLev   = 0
               icount = 0

            endif

            !
            ! Transform from Grid to Spec
            !

            icount = icount + 1
         enddo

      enddo

      !
      ! Close BAM Files
      !

      call bamOut%Close( )

      return
   end subroutine
!EOC
!-------------------------------------------------------------------------
!    NOAA/NCEP, National Centers for Environmental Prediction GSI        !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  load_grid --- strip off south/north latitude rows
!
! !INTERFACE:
!
 subroutine load_grid(grid_in,grid_out)

! !USES:
   use gridmod, only: grd => grd_a
   use kinds, only : r_kind, i_kind
   use MiscMod, only: r8
   implicit none

! !INPUT PARAMETERS:

   real(r_kind),dimension(max(grd%iglobal,grd%itotsub)), intent(in   ) :: grid_in  ! input grid
   real(r8),    dimension(grd%nlon,grd%nlat-2),          intent(  out) :: grid_out ! output grid

! !DESCRIPTION: This routine prepares grids for use in splib
!               grid to spectral tranforms.  This preparation
!               entails to two steps
!                  1) reorder indexing of the latitude direction.
!                     The GSI ordering is south to north.  The 
!                     ordering assumed in splib routines is north
!                     to south.
!                  2) The global GSI adds two latitude rows, one
!                     for each pole.  These latitude rows are not
!                     needed in the grid to spectral transforms of
!                     splib.  The code below strips off these
!                     "pole rows"
!
! !REVISION HISTORY:
!   2004-08-27  treadon
!   2013-10-25  todling - move from gridmod to this module
!
! !REMARKS:
!   language: f90
!   machine:  ibm rs/6000
!
! !AUTHOR:
!   treadon          org: np23                date: 2004-08-27
!
!EOP
!-------------------------------------------------------------------------
   integer(i_kind) i,j,k,nlatm1,jj
   real(r_kind),dimension(grd%nlon,grd%nlat):: grid

!  Transfer input grid from 1d to 2d local array.  As loading
!  local array, reverse direction of latitude index.  Coming
!  into the routine the order is south --> north.  On exit
!  the order is north --> south
   do k=1,grd%iglobal
      i=grd%nlat-grd%ltosi(k)+1
      j=grd%ltosj(k)
      grid(j,i)=grid_in(k)
   end do
   
!  Transfer contents of local array to output array.
   nlatm1=grd%nlat-1
   do j=2,nlatm1
      jj=j-1
      do i=1,grd%nlon
         grid_out(i,jj)=real(grid(i,j),r8)
      end do
   end do
   
   return
 end subroutine load_grid

!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !ROUTINE: GenGatherBAM
!              
! !DESCRIPTION: Transfer contents of 3d subdomains to 2d work arrays over pes
!
!             
!                 
!\\
!\\
! !INTERFACE:
!

   subroutine GenGatherBAM( g_z, g_ps, g_tv, g_vor, g_div, g_q, g_ql, g_qi,&
                            wvar, wlev, icount, MyPe, VNAme, work  )
!
! !USES:
!

      use kinds, only: r_kind,i_kind
      use mpimod, only: npe,mpi_comm_world,ierror,mpi_rtype
      use gridmod, only: strip
      use gridmod, only: grd => grd_a

      implicit none
!
! !INPUT PARAMETERS:
!
  real(r_kind),dimension(:,:),    intent(in   ) :: g_z
  real(r_kind),dimension(:,:),    intent(in   ) :: g_ps
  real(r_kind),dimension(:,:,:),  intent(in   ) :: g_tv
  real(r_kind),dimension(:,:,:),  intent(in   ) :: g_vor
  real(r_kind),dimension(:,:,:),  intent(in   ) :: g_div
  real(r_kind),dimension(:,:,:),  intent(in   ) :: g_q
  real(r_kind),dimension(:,:,:),  intent(in   ) :: g_ql
  real(r_kind),dimension(:,:,:),  intent(in   ) :: g_qi
  integer(i_kind),dimension(npe), intent(in   ) :: wvar
  integer(i_kind),dimension(npe), intent(in   ) :: wlev
  integer(i_kind),                intent(in   ) :: icount
  integer(i_kind),                intent(in   ) :: MyPe
  character(len=*),dimension(:), intent(in) :: VName


!
! !OUTPUT PARAMETERS:
!

  real(r_kind), dimension(:),     intent(inout) :: work

!
! !REVISION HISTORY:
!   2013-06-19  treadon
!   2013-10-24  todling   - update interface to strip
!   2016-08-04  de Mattos - Adpat to BAM model
!
!EOP
!-------------------------------------------------------------------------
!BOC

  real(r_kind),dimension(grd%lat1*grd%lon1,npe) :: sub

  integer :: ivar
  integer :: ilev
  integer :: Pe
  integer :: var



!$omp parallel do  schedule(dynamic,1) private(Pe,ivar,ilev)

      do Pe = 1, icount

         ivar = wvar(Pe)
         ilev = wlev(Pe)

         select case ( ivar )

            case (1) ! Topography

               call strip ( g_z(:,:), sub(:,Pe) )

            case (2) ! ln Surface Pressure

               call strip ( g_ps(:,:), sub(:,Pe) )
           
            case (3) ! Virtual Temperature

               call strip ( g_tv(:,:,ilev), sub(:,Pe) )

            case (4) ! Divergence

               call strip ( g_div(:,:,ilev), sub(:,Pe) )

            case (5) ! Vorticity

               call strip ( g_vor(:,:,ilev), sub(:,Pe) )

            case (6) ! Specific Humidity

               call strip ( g_q(:,:,ilev), sub(:,Pe) )

            case (7) ! Liq mixing ratio prognostic

               call strip ( g_ql(:,:,ilev), sub(:,Pe) )

            case (8) ! Ice mixing ratio prognostic
        
               call strip ( g_qi(:,:,ilev), sub(:,Pe) )

         end select
      enddo

      call mpi_alltoallv(sub,           &
                         grd%isc_g,     &
                         grd%isd_g,     &
                         mpi_rtype,     &
                         work,          &
                         grd%ijn,       &
                         grd%displs_g,  &
                         mpi_rtype,     &
                         mpi_comm_world,&
                         ierror         &
                         )
   end subroutine
!EOC

subroutine lterp(iField, iMend, oMend, oField)
   use MiscMod, only: GetImaxJmax, GetLongitudes, GetGaussianLatitudes, r8
   use coord_compute,   only: compute_grid_coord!, r8

   real(r8), intent(in   ) :: iField(:,:)
   integer,  intent(in   ) :: iMend
   integer,  intent(in   ) :: oMend
   real(r8), intent(inout) :: oField(:,:)

   ! Local variables

   real    :: y1, y2, y3, y4
   real    :: t, u, d
   integer :: iIMax, iJMax
   integer :: oIMax, oJMax
   integer :: i, i1, ii
   integer :: j, j1, jj
   integer :: iret
   
   real,              dimension(200)     :: gDesci             ! Grid description parameters of input field
   real,              dimension(200)     :: gDesco             ! Grid description parameters of output field
   real(r8), allocatable, dimension(:)   :: olat               ! latitudes in degrees of output field
   real(r8), allocatable, dimension(:)   :: olon               ! longitudes in degrees of output field
   real(r8), allocatable, dimension(:)   :: ilat               ! latitudes in degrees of input field
   real(r8), allocatable, dimension(:)   :: ilon               ! longitudes in degrees of input field

   real,     allocatable, dimension(:)   :: rlat               ! latitudes in degrees of input field
   real,     allocatable, dimension(:)   :: rlon               ! longitudes in degrees of input field

   real,     allocatable, dimension(:)   :: xpts
   real,     allocatable, dimension(:)   :: ypts
   real(r8),     allocatable, dimension(:,:) :: ya
   real(r8),     allocatable, dimension(:)   :: x1a
   real(r8),     allocatable, dimension(:)   :: x2a

!   real, parameter    :: udef  = 9.999E20

   !
   ! Get input info
   !
   call GetImaxJmax(iMend, iIMax, iJMax)
   if(size(iField,1).ne.iImax.or.size(iField,2).ne.iJMax)then
      write(stdout,*)'iField error:'
      write(stdout,*)'isize in  :',size(iField,1),size(iField,2)
      write(stdout,*)'isize jcap:',iIMax,iJMax
      stop
   endif

   allocate(ilon(iIMax))
   call GetLongitudes (iIMax, 0.0_r8, ilon)
   allocate(ilat(iJMax))
   call GetGaussianLatitudes(iJMax,ilat)
   gDesci    = 0
   gDesci(1) = 4
   gDesci(2) = iIMax
   gDesci(3) = iJMax
   gDesci(4) = maxval(ilat)
   gDesci(5) = minval(ilon)
   gDesci(6) = 128
   gDesci(7) = minval(ilat)
   gDesci(8) = maxval(ilon)
   gDesci(9) = abs(ilon(2)-ilon(1))
   gDesci(10)= iJMax/2
   gDesci(11)= 64
   gDesci(20)= 0

   DeAllocate(ilat)
   DeAllocate(ilon)

   !
   ! Get output info
   !
   call GetImaxJmax(oMend, oIMax, oJMax)
   if(size(oField,1).ne.oImax.or.size(oField,2).ne.oJMax)then
      write(stdout,*)'oField error:'
      write(stdout,*)'osize in  :',size(oField,1),size(oField,2)
      write(stdout,*)'osize jcap:',oIMax,oJMax
      stop
   endif

   allocate(olon(oIMax))
   call GetLongitudes (oIMax, 0.0_r8, olon)
   allocate(olat(oJMax))
   call GetGaussianLatitudes(oJMax,olat)

   Allocate(xpts(oIMax))
   Allocate(ypts(oJMax))

   call compute_grid_coord(gDesci, real(olon,4), real(olat,4), udef, xpts, ypts, iret)

!   write(stdout,*)'xpts:',minval(xpts),maxval(xpts)
   ypts = ypts + 1 ! <- to use with a ghost zone

   !
   ! Create a ghost zone at north and south poles
   !

   Allocate(ya(iIMax,iJMax+2))
   d               = 1/float(iIMax)
   ya(:,        1) = sum(iField(:,1))*d
   ya(:,2:iJMax+1) = iField
   ya(:,  iJMax+2) = sum(iField(:,iJMax))*d

   !
   ! Create input Grid 
   !

   Allocate(x1a(iIMax+1))
   x1a(1:iIMax+1) = (/(i,i=1,iIMax+1)/) !  <- Global field last point + 1 = 1st point
   Allocate(x2a(iJMax+2))
   x2a(1:iJMax+2) =(/(i,i=1,iJMax+2)/)

   do jj=1,oJMax
      do ii=1,oIMax

         if(xpts(ii).eq.udef .or. ypts(jj).eq. udef)then
            write(stdout,'(A1,x,2I6)')'wrong x,y points at:',ii,jj
            oField(ii,jj) = udef
            cycle
         endif

         i  = int(xpts(ii))
         i1 = mod(i,iIMax) + 1 ! <- Global field last point + 1 = 1st point
         j  = int(ypts(jj))
         j1 = j + 1


         y1 = ya( i, j)
         y2 = ya(i1, j)
         y3 = ya(i1,j1)
         y4 = ya( i,j1)

         t  = (real(xpts(ii),r8) - x1a(i))/(x1a(i+1)-x1a(i))
         u  = (real(ypts(jj),r8) - x2a(j))/(x2a(j+1)-x2a(j))

         oField(ii,jj) = (1-t)*(1-u)*y1 + t*(1-u)*y2 + t*u*y3 + (1-t)*u*y4

      enddo
   enddo

   DeAllocate(xpts)
   DeAllocate(ypts)
   DeAllocate(ya)
   DeAllocate(x1a)
   DeAllocate(x2a)

end subroutine

!!-----------------------------------------------------------------------------!
!!             Modeling and Development Division - DMD/CPTEC/INPE              !
!!-----------------------------------------------------------------------------!
!!
!!BOP
!!
!! !IROUTINE: WriteField_MPI - write fields of BAM files from a MPI Pe.
!!
!! 
!! !DESCRIPTION: Esta rotina escreve um campo do modelo BAM 
!!               
!!
!! !INTERFACE:
!!   
!  subroutine WriteField(OutUnit, field, OutPe, iCount, istat)
!
!     implicit none
!!
!! !INPUT PARAMETERS:
!! 
!     ! Output logical unit
!     integer(i4),           intent(in   ) :: OutUnit
!
!     ! Field to be writed
!     real(r8),              intent(in   ) :: field(:)
!
!     ! How many Pe's are working
!     integer(i4),           intent(in   ) :: iCount
!
!     ! What Pe will write
!     integer(i4),           intent(in   ) :: OutPe
!!
!! !OUTPUT PARAMETERS:
!! 
!
!     integer(i4), optional, intent(  out) :: istat
!
!!
!! !REVISION HISTORY: 
!!
!!  11 Oct 2016 - J. G. de Mattos - Initial Version
!!
!!
!!EOP
!!-----------------------------------------------------------------------------!
!!BOC
!!
!     character(len=100), parameter :: myname_=':: WriteFields_MPI_( ... )'
!
!     real(r8), allocatable :: buff(:)
!     integer :: sizebuff
!     integer :: Pe
!     integer :: iret
!   
!     if(present(istat)) istat = 0
!
!     sizebuff = size(field)
!
!     allocate(buff(sizebuff))
!      
!     do Pe = 0, iCount-1
!
!        if(Pe .eq. OutPe)then
!        
!           call WriteField_Serial_(OutUnit, field, iret)
!
!           if(iret .ne. 0)then
!              write(stdout,*)trim(myname_),': error to write field, ',iret
!              if(present(istat)) istat = iret
!              return
!           endif
!
!        else
!
!           call mpi_recv(buff, sizebuff, MPI_DOUBLE, Pe, MPITag, MPI_COMM_WORLD, status, iret)
!
!           call WriteField_Serial_(OutUnit, buff, iret)
!           if(iret .ne. 0)then
!              if(present(istat)) istat = iret
!              return
!           endif
!        endif
!
!     end do
!      
!     deallocate(buff)
!
!  end subroutine
!!
!!EOC
!!
!!-----------------------------------------------------------------------------!
!!             Modeling and Development Division - DMD/CPTEC/INPE              !
!!-----------------------------------------------------------------------------!
!!
!!BOP
!!
!! !IROUTINE: SendField -
!!
!! 
!! !DESCRIPTION:
!!
!!
!! !INTERFACE:
!!   
!  subroutine SendField(MyPe, toPe, Field, istat)
!
!     implicit none
!!
!! !INPUT PARAMETERS:
!!
!     ! Source Pe
!     integer(i4),           intent(in   ) :: MyPe
!
!     ! Target Pe
!     integer(i4),           intent(in   ) :: toPe
!
!     ! Field to be send
!     real(r8),              intent(in   ) :: field(:)
!!
!! !OUTPUT PARAMETERS:
!! 
!     integer(i4), optional, intent(  out) :: istat
!!
!! !REVISION HISTORY: 
!!
!!  11 Oct 2016 - J. G. de Mattos - Initial Version
!!
!!
!!EOP
!!-----------------------------------------------------------------------------!
!!BOC
!!
!     character(len=100), parameter :: myname_=':: BAM_SendField_( ... )'
!
!     integer :: iret
!     integer :: sizefield
!
!     if(present(istat)) istat = 0
!
!     sizefield = size(field)
!      
!     call mpi_send(field, sizefield, MPI_DOUBLE, toPe, MPITag, MPI_COMM_WORLD, iret)
!
!     if(iret .ne. 0)then
!        write(*,'(2A,I5,1x,A,1x,I5)')trim(myname_),': ERROR to send field from',MyPe,'to',toPe
!        if(present(istat)) istat = iret
!        return
!     endif
!
!  end subroutine
!!
!!EOC
!!-----------------------------------------------------------------------------!

end module
