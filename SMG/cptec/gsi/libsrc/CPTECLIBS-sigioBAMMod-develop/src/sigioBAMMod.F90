!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOI
!
! !TITLE: Input/Output BAM interface Documentation \\ Version 1.0.0
!
! !AUTHORS: João Gerd Zell de Mattos
!
! !AFFILIATION: Modeling and Development Division, CPTEC/INPE
!
! !DATE: October 11, 2016
!
! !INTRODUCTION:
!      Input/Output BAM interface (sigio\_BAMMod) is a Fortran 90 collection of 
!      routines/functions for accessing the forecasting data files of the
!      Brazilian Global Atmospheric Model (BAM) in spectral format.
!
! \subsection{BAM data files}
!
! Os arquivos de previsão do modelo BAM são constituídos por dois arquivos em
! formatos distintos: 
!   \begin{enumerate}
!      \item Um arquivo no formato ASCII, denominado arquivo {\bf dir} contém um 
!            cabeçalho descrevendo algumas informações da simulação, tais
!            como: (a) data inicial e final da simulação, (b) número de níveis
!            verticais e (c) delta Z entre os níveis verticais. Este arquivo também
!            possui uma tabela contendo as variáveis disponíveis no arquivo de
!            previsão, identificando o tipo de cada variável
!            (Diagnostica/Prognóstica), a quantidade de níveis verticais e se é
!            uma variável em ponto de grade ou em formato spectral;
!      \item Um arquivo no formato `IEEE', denominado arquivo {\bf fct} que contém um 
!            cabeçalho com a data da simulacão, seguido pelos campos diagnosticos e 
!            prognosticos do modelo BAM. Estes campos estão na mesma sequência em que 
!            aparecem na tabela de variáveis do aquivo {\bf dir} \\
!            {\bf NOTA:} um arquivo no formato `IEEE' é refenciado como um
!            arquivo fortran não formatado (``unformatted'') com acesso sequencial
!            e com o comprimento dos registros variável. No sistema Linux/Unix comum, é 
!            apenas um arquivo com registros utilizando palavras de 4 bytes e um 
!            ``cabecalho'' de 4 bytes indicando o tamanho do registro em bytes.
!   \end{enumerate}
!   
! \newpage
! \subsection{Principais Rotinas/Funções}
!
! \begin{verbatim}
!  ------------------------+-----------------------------------------
!     Routine/Function     |             Description
!  ------------------------+-----------------------------------------
!                          |
!   Open                   | Open a BAM file
!   Close                  | Close a BAM file
!   GetField               | Return a BAM field from a file
!   GetUV                  ! Return fields of zonal and meridional wind 
!   GetOneDim              | Return just one of 4 BAM field dimensions
!   GetDims                | Return BAM field dimensions
!   GetOFS                 !
!   GetNlevels             | Return the number of leves of a field
!   GetVarNames            ! Return a table with all variables names (one by row)
!   GetVerticalCoord       | Return information about vertical coordinate
!   GetWCoord              | Return BAM world coordinates
!   GetTimeInfo            | Return time info
!   GetPhysics             | Return Physics Scheme
!   isHybrid               | Return .true. or .false. (.false. model use Sigma)
!   WriteAnlHeader         | Write the Anl file header
!   WriteField             | Write a BAM field (mpi or serial)
!                          |
!  ------------------------+-----------------------------------------
! \end{verbatim}
!
! \subsection{Exemplo de uso}
!
! O primeiro passo é carregar este módulo no programa fortran e definir uma estrutura 
! de dados contendo as informações do BAM.
! 
! \begin{enumerate}
!
!    \item Defina no início do programa fortran o uso do módulo sigio\_BAMMod:
!
!    \begin{verbatim}
!       use sigioBAMMod, only: BAMFile
!    \end{verbatim}
!
!    \item Defina uma variável que conterá a estrutura de dados:
!
!    \begin{verbatim}
!       type(BAMFile) :: bam
!    \end{verbatim}
!
!    \item Defina algumas variáveis para auxiliar na leitura das informacões
!
!    \begin{verbatim}
!        character(len=256) :: FileFct
!        character(len=256) :: FileDir
!        integer            :: iMax
!        integer            :: jMax
!        integer            :: kMax
!        integer            :: Mend
!        integer            :: MnvWv2
!        integer            :: iret
!        real, allocatable  :: grid(:)
!    \end{verbatim}
!
!    \item Defina qual o arquivo deverá ser lido (fct/anl,dir). Os arquivos do bam são
!    compostos por um arquivo de cabecalho (dir) e um arquivo binario (fct/anl):
! 
!    \begin{verbatim}
!      FileDir = 'GFCTNMC20131231002013123106F.dir.TQ0062L028'
!      FileFct = 'GFCTNMC20131231002013123106F.fct.TQ0062L028'
!    \end{verbatim}
! 
!    \item Utilize a rotina específica para a abertura do arquivo:
! 
!    \begin{verbatim}
!  
!       call bam%Open(FileFct, FileDir, iret)
!   
!    \end{verbatim}
!
!    \item Obtenha as informacões da grade:
!
!    \begin{verbatim}
!
!    iMax=bam%GetOneDim('imax')
!    jMax=bam%GetOneDim('jmax')
!    kMax=bam%GetOneDim('kmax')
!    Mend=bam%GetOneDim('mend')
!
!    MnWv2 = (Mend+1)*(Mend+2)
!
!    \end{verbatim}
! 
! 
!    \item Faça a leitura do campo disponível no BAM, defina antes o nome da 
!          variável, o nível vertical e aloque um vetor do tamanho necessário para
!          retornar o campo solicitado:
! 
!    \begin{verbatim}
!    
!       Allocate(grid(iMax*jMax))
!       
!       VName = 'VIRTUAL TEMPERATUE'
!       ilev  = 1
!       
!       call bam%GetField(trim(VName(ivar)), ilev, grid, iret)
!   
!    \end{verbatim}
! 
!    \item Depois que forem lidos todos os campos necessários do modelo BAM feche o arquivo:
!    
!    \begin{verbatim}
!        call bam%close(iret)
!    \end{verbatim}
! 
! \end{enumerate}
! 
!Veja os prólogos na próxima seção para obter detalhes adicionais.
!
!EOI
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!
Module sigioBAMMod

#ifdef useMPI
   use MPI
#endif
   use ModConstants
   use TransformTools, only : specTrans
   implicit none
   Private

#ifdef useMPI
   ! MPI Parameters
   integer :: MPITag = 100
   integer :: status(MPI_Status_size)
#endif

   

   ! BAM header structure
   type, extends(specTrans):: BAM_head
!      private
      character (len = 21) :: ittl
      character (len = 40) :: specal
      character (len =150) :: jttl
      character (len =  4) :: sdain
      character (len =  4) :: dtin
      character (len =  5) :: trunc
      integer (I4)         :: nexp
      integer (I4)         :: mMax
      integer (I4)         :: iosize
      integer (I4)         :: kmax
      integer (I4)         :: lsm
      integer (I4)         :: pbl
      integer (I4)         :: ihr, idy, imo, iyr
      integer (I4)         :: fhr, fdy, fmo, fyr

      !  - ID of Vertical Coordinate
      !    = 1 Sigma (0 for old files)
      !    = 2 ec-Hybrid
      !    = 3 NCEP General Hybrid
      real (I4) :: idvc

      ! vertical coordinate info
      !  - ID of Sigma Structure (IDSL)
      !    = 1 Phillips
      !    = 2 Mean
      real (I4) :: idsl

      ! sigma levels
      real (r8), dimension (:), pointer :: delSig  => null()
      real (r8), dimension (:), pointer :: SigInt  => null()
      real (r8), dimension (:), pointer :: SigMid  => null()

      ! Hybrid-Sigma levels
      real (r8), dimension (:), pointer :: ak  => null()
      real (r8), dimension (:), pointer :: bk  => null()
      real (r8), dimension (:), pointer :: ck  => null()


      ! lat/lon iformations
      real(r4), dimension(:), pointer :: rlon => null()! real longitude
      real(r4), dimension(:), pointer :: clon => null()! rlon cosine
      real(r4), dimension(:), pointer :: slon => null()! rlon sine
      real(r4), dimension(:), pointer :: rlat => null()! real latitude
      real(r4), dimension(:), pointer :: clat => null()! rlat cosine
      real(r4), dimension(:), pointer :: slat => null()! rlat sine

   end type

   ! BAM fields structure
   type, public:: BAM_fld
      character (len = 40)   :: Name  ! Long field Name
      character (len =  4)   :: ProDia! Prognostig / Diagnostic 
      integer (I4)           :: nharm ! variable array size
      integer (I8)           :: pmark ! Initial Position Marker (bytes)
      integer (I4)           :: nlevs ! variable # of levs
      integer (I4)           :: fldunit
      real (r4), allocatable :: grid(:,:,:)
      type(BAM_fld), pointer :: next => null()
   end type

   ! BAM data type structure
   type, extends(BAM_head), public :: BAMFile
!      private
      character(len=strlen)  :: fBin  ! Anl or Fct file name
      character(len=strlen)  :: fHead ! dir, dic or din file name (header informations)
      character(len=strlen)  :: ftype ! type of binary file (anl or fct)
      character(len=strlen)  :: oMode ! Open mode (w: write, r: read)
      logical, public        :: isHybrid ! .true. or .false. (if .false. vertical coodinate is sigma)
      integer (I4), pointer  :: uBin => null()  ! Anl of Fct file opened unit
      integer (I4), pointer  :: uHead => null() ! dir, dic or din file opened unit
!      type(BAM_head)         :: head  ! header informations (see BAM_Head type)
      type(BAM_fld), pointer :: root => null()
      type(BAM_fld), pointer :: flds => null()
      integer (I4)           :: fcount
      contains
         procedure, public :: open => bamOpen_
         procedure, public :: close => bamClose_
         procedure, public :: getOneDim => bamGetOneDim_
         procedure, public :: getTimeInfo => bamGetTimeInfo_
         procedure, public :: getDims => bamGetDims_
         procedure, public :: getOFS => bamGetOFS_
         procedure, public :: getNLevels => bamGetNLevels_
         procedure, public :: getWCoord => bamGetWCoord_
         procedure, public :: getVerticalCoord => bamGetVerticalCoord_
         procedure, public :: getPhysics => bamGetPhysics_
         procedure, public :: getSigValues => bamGetSigValues_
         procedure, public :: writeAnlHeader => bamWriteAnlHeader_
         procedure, public :: getVarNames => bamGetVarNames_

         procedure         :: SGetFLD_1d, SGetFLD_2d,&
                              DGetFLD_1d, DGetFLD_2d
         generic,   public :: getField => SGetFLD_1d, SGetFLD_2d,&
                                          DGetFLD_1d, DGetFLD_2d

         procedure         :: getUV1D_, getUV2D_
         generic,   public :: getUV => getUV1D_, getUV2D_
   end type
!
! !PUBLIC MEMBER FUNCTIONS:
!

   public :: BAM_WriteField     ! Write a BAM field (mpi or serial)
#ifdef useMPI
   public :: BAM_SendField      ! if mpi code is used to send a field to write
#endif
   public :: BAM_GetAvailUnit


   interface BAM_GetHeadInfo
      module procedure ReadHead_
   end interface

   interface ReadField
      module procedure SRead1D_, SRead2D_, SRead3D_, &
                       DRead1D_, DRead2D_, DRead3D_, &
                       SkipRead_
   end interface

#ifdef useMPI

   interface BAM_SendField
      module procedure   BAM_SendField_, BAM_SendFieldr4_
   end interface

   interface BAM_WriteField
      module procedure   WriteField_MPI_,   WriteField_Serial_, &
                       WriteField_MPIr4_, WriteField_Serialr4_
   end interface

#else

   interface BAM_WriteField
      module procedure   WriteField_Serial_, WriteField_Serialr4_
   end interface

#endif

   interface perr
      module procedure perr1_, perr2_
   end interface

   interface BAM_GetAvailUnit
      module procedure GetAvailUnit
   end interface
   
!
! !REVISION HISTORY:
!
! 	11 Oct 2016 - J. G. Z. de Mattos - Initial Version
!
!-----------------------------------------------------------------------------!


   contains
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: BAM_Open - routine used to open BAM type files for read and
!                       write purposes.
!
! 
! !DESCRIPTION: Esta rotina é uma interface para abrir os arquivos
!               do modelo BAM. Há dois modos para a abertura dos aquivos do BAM, 
!               um somente leitura e outro para escrita. No modo somente
!               leitura, caso não seja definido explicitamente qual o arquivo
!               que deve ser lido, abre-se por padrão os arquivos necessários
!               para a leitura do arquivo de previsão, ou seja, abre o arquivo
!               DIR e o arquivo FCT. O arquivo DIR serve como um arquivo
!               descritor que indica qual a posicão de cada variável dentro do
!               arquivo FCT. Já no modo de escrita, caso não seja especificado,
!               abre-se somente o arquivo de condicão inicial do BAM. No modo de
!               escrita, caso o arquivo já exista ele é sobrescrito pelo novo
!               arquivo.
!
!
! !INTERFACE:
!
   subroutine bamOpen_(self, hFile, bFile, mode, ftype, istat)

      implicit none
!
! !INPUT PARAMETERS:
!
      ! BAM file structure
      class(BAMFile),             intent(inout) :: self

      ! BAM binary file (Analisys or Forecast)
      character(len=*),           intent(in   ) :: hFile

      ! BAM header file
      character(len=*), optional, intent(in   ) :: bFile

      ! Open mode: R - read-only; W - Write
      character(len=*), optional, intent(in   ) :: mode

      ! File type: dir, fct, anl
      character(len=*), optional, intent(in   ) :: ftype
      ! Can be:
      !        1. fct: forecast BAM file
      !        2. anl: initial condition BAM file
      !
      ! At Read-only mode default ftype is to open dir and fct together
      ! At Write mode default ftype is to open anl file only
      !
!
! !OUTPUT PARAMETERS:
!
      integer,          optional, intent(  out) :: istat

!
! !SEE ALSO:
!
!   Open_( ) interface to open each type of BAM file
!
! !REVISION HISTORY: 
!  11 Oct 2016 - J. G. de Mattos - Initial Version
!  29 sep 2020 - J. G. de Mattos - modify to class
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
      character(len=100), parameter :: myname_='bamOpen_( ... )'

      integer          :: iret
      logical          :: ok

      if(present(istat))istat = 0

      !----------------------------------------------------------------!
      ! Some Sanity Check
      if(present(ftype))then
         self%ftype = lcase(ftype)
      else
         if (present(bFile))then
            self%ftype = 'fct'
         else
            self%ftype = 'dir'
         endif
      endif

      if(present(mode))then
         self%oMode = trim(ucase(mode))
      else
         self%oMode = 'R'
      endif

      if (trim(self%ftype) .ne. 'anl' .and. &
          trim(self%ftype) .ne. 'fct' .and. &
          trim(self%ftype) .ne. 'dir')then
         write(stdout,*)trim(myname_),': ftype should be `anl`, `fct` .or. `dir`'
         write(stdout,*) 'ftype was set to : ', trim(self%ftype)
         stop
      endif

      if((trim(self%ftype) .eq. 'anl' .or. trim(self%ftype) .eq. 'fct').and.&
          .not.present(bFile)) then

            write(stdout,*)trim(myname_),': Missing binary file name!'
            stop
            
      endif
      !----------------------------------------------------------------!

      self%fHead = trim(hFile)

      call openHeader_(self, iret)
      if (iret .ne. 0)then
         call perr(trim(myname_),'problem to OPEN header file '//trim(hFile), iret)
         if (present(istat))then
            istat = iret
            return
         endif
         stop
      endif

      call ReadHead_(self,iret)
      if (iret .ne. 0)then
         call perr(trim(myname_),'problem to READ header file '//trim(hFile), iret)
         if (present(istat))then
            istat = iret
            return
         endif
         stop
      endif

      ! if ftype is equal a dir, dont open binary file
      if (trim(self%ftype) .eq. 'dir') return

      self%fBin  = trim(bFile)
      call openBin_(self, iret)
      if (iret .ne. 0)then
         call perr(trim(myname_),'problem to OPEN binary file '//trim(bFile), iret)
         if (present(istat))then
            istat = iret
            return
         endif
         stop
      endif


   end subroutine
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: OpenBin_ - routine that call low level routines to open FCT type 
!                    of BAM files
!
! 
! !DESCRIPTION: Esta rotina é uma interface para abrir os arquivos binários do
!               do modelo BAM. Há dois modos para a abertura dos aquivos: 
!               (1) somente leitura e (2) escrita.  É importante salientar que no 
!               modo de escrita, caso o arquivo já exista ele é sobrescrito pelo 
!               novo arquivo.
!
!
! !INTERFACE:
!

   subroutine OpenBin_( BFile, istat )
      implicit none
!
! !INPUT PARAMETERS:
!
      ! BAM file structure
      type(BAMFile),     intent(inout)  :: BFile
      
!
! !OUTPUT PARAMETERS:
!
      integer, optional, intent(  out)  :: istat
      !
      ! -80 : File not found
      !

! !REVISION HISTORY: 
!  11 Oct 2016 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
      character(len=*), parameter :: myname_ = 'OpenBin_( ... )'

      logical          :: existe
      integer          :: iret

#ifdef DEBUG
      WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif

      if(present(istat)) istat = 0

      inquire(file=trim(BFile%fBin),exist=existe)

      !Get next available logical unit
      allocate(BFile%uBin)
      BFile%uBin = BAM_GetAvailUnit(  )

      select case (trim(BFile%oMode))

         case ('R')

            if(.not.existe)then
               iret = -80
               if(present(istat))istat=iret
               nullify(BFile%uBin)
               return
            endif

            open( unit   = BFile%uBin,         &
                  file   = trim(BFile%fBin),   &
                  status = 'old',              &
                  form   = 'unformatted',      &
                  action = 'read',             &
                  access = 'stream',           &
                  iostat = iret                &
                )

         case ('W')

            if(existe)then
               write(6,'(A)')'Overwriting existing file!'
            endif

            open( unit   = BFile%uBin,         &
                  file   = trim(BFile%fBin),   &
                  status = 'unknown',          &
                  form   = 'unformatted',      &
                  action = 'write',            &
                  iostat = iret                &
                )

      end select

      if (iret.ne.0)then
!         write(6,'(A,1x,A,1x,A,1x,I3)')'BAM error =>',trim(myname_),&
!                                       'Error to open model file: '//trim(BFile%ffct),iret
         if(present(istat))istat=iret
         return
      endif

   endsubroutine
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: OpenHeader_ - routine that call low level routines to open Header 
!                    of BAM files (dir, dic or din)
!
! 
! !DESCRIPTION: Esta rotina é uma interface para abrir os arquivos
!               descritores (DIR) do modelo BAM. Há dois modos para a abertura 
!               dos aquivos: (1) somente leitura e (2) escrita.  É importante 
!               salientar que no modo de escrita, caso o arquivo já exista ele
!               é sobrescrito pelo novo arquivo.
!
!
! !INTERFACE:
!

   subroutine OpenHeader_( BFile, istat )
      implicit none
!
! !INPUT PARAMETERS:
!
      ! BAM file structure
      type(BAMFile),     intent(inout)  :: BFile
      
! !OUTPUT PARAMETERS:
!
      integer, optional, intent(  out)  :: istat
      !
      ! -80 : File not found
      !

! !REVISION HISTORY: 
!  11 Oct 2016 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
      character(len=*), parameter :: myname_ = 'OpenDir_( ... )'

      logical          :: existe
      integer          :: iret

#ifdef DEBUG
      WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif

      if(present(istat)) istat = 0

      inquire(file=trim(BFile%fHead),exist=existe)

      !Get next available logical unit
      allocate(BFile%uHead)
      BFile%uHead = BAM_GetAvailUnit(  )

!      select case (trim(md))
!
!         case ('R')

            if(.not.existe)then
               iret = -80
               if(present(istat))istat=iret
               nullify(BFile%uHead)
               return
            endif

            open(unit   = BFile%uHead,       &
                 file   = trim(BFile%fHead), &
                 status = 'old',             &
                 action = 'read',            &
                 form   = 'formatted',       &
                 iostat = iret               &
                )

!         case ('W')
!
!            if(existe)then
!               write(stdout,'(A26)')'Overwriting existing file!'
!            endif
!
!            open(unit   = BFile%uHead,       &
!                 file   = trim(BFile%fHead), &
!                 status = 'new',             &
!                 action = 'write',           &
!                 form   = 'formatted',       &
!                 iostat = iret               &
!                )
!
!      end select

      if (iret.ne.0)then
!         write(6,'(A,1x,A,1x,A,1x,A,1x,I3)')'BAM error =>',trim(myname_),&
!                                       'Error to open model file: ',trim(BFile%fdir),iret
         if(present(istat))istat=iret
         return
      endif

   endsubroutine
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: ReadHead_ - routine that read informations from dir file and
!                        initialize all BAM structure used to read and write BAM
!                        files.
!
! 
! !DESCRIPTION: Esta rotina lê as informacões contidas no arquivo descritor do
!               modelo BAM (arquivo DIR) e então inicializa a estrutura de dados
!               (BAMFile) usada neste módulo para ler e escrever os arquivos do
!               BAM.
!
!
! !INTERFACE:
!

   subroutine ReadHead_(BFile, istat )
   
      USE MiscMod, ONLY: GetLongitudes, GetGaussianLatitudes, GetImaxJmax
      implicit none
!
! !INPUT PARAMETERS:
!
      ! BAM file structure
      type(BAMFile), intent(inout)   :: BFile

!
! !OUTPUT PARAMETERS:
!
      integer, optional, intent(out) :: istat
!
! !SEE ALSO:
!
!     setsig( ) - calcula niveis sigma
!     InitParameters( ) - Inicializa algumas informacoes do modelo BAM
!
! !REVISION HISTORY: 
!  11 Oct 2016 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
      character(len=*), parameter :: myname_ = 'ReadHead_( ... )'

      integer, parameter :: iosize = 4 ! single precision file

      integer :: k
      integer :: iret
      integer :: icol, isum
      integer :: pmark
      integer :: lastSize
      integer :: Mend
      real(kind=r8), allocatable :: glats(:)
      real(kind=r8), allocatable :: glons(:)

#ifdef DEBUG
      WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif

      rewind(BFile%uHead)

      read (BFile%uHead,'(A21)') BFile%ittl
      
      BFile%isHybrid = check('hibrid',BFile%ittl)

      read (BFile%uHead,'(A4,1X,A4,1X,A5,1X,11I5,1X,A4)') &
                        BFile%nexp, BFile%sdain, BFile%trunc, &
                        BFile%mMax, BFile%kmax, BFile%kmax, &
                        BFile%fhr, BFile%fdy, BFile%fmo, BFile%fyr, &
                        BFile%ihr, BFile%imo, BFile%idy, BFile%iyr, &
                        BFile%dtin


      if (BFile%isHybrid)then
         read (BFile%uHead,'(2A150)') BFile%jttl
         read (BFile%jttl(61:62),*)BFile%LSM
         read (BFile%jttl(69:70),*)BFile%PBL

         allocate(BFile%ak(BFile%kmax))
         allocate(BFile%bk(BFile%kmax))

         ! reverse coefficients to be from top to bottom
         ! This reverse order is just at model environment
         ! GANL and GFCT files aren't in reverse order.
         ! But at output dir file these coeffients are in 
         ! reverse order to. So we need, by now, reverse to be
         ! in conform with output fields.
         ! 
         icol = 5
         isum = 0
         do while(isum .lt. BFile%kmax)

            icol = min(icol, (BFile%kmax-isum))
            !read (BFile%uHead,*)(BFile%ak(k+isum),k=1,icol)
            read (BFile%uHead,*)(BFile%ak(BFile%kmax-(k+isum)+1),k=1,icol)
            isum = isum + icol

         enddo

         icol = 5
         isum = 0
         do while(isum .lt. BFile%kmax)

            icol = min(icol, (BFile%kmax-isum))
            !read (BFile%uHead,*)(BFile%bk(k+isum),k=1,icol)
            read (BFile%uHead,*)(BFile%bk(BFile%kmax-(k+isum)+1),k=1,icol)
            isum = isum + icol

         enddo
         
         ! This is only necessary for BAM Hybrid model 
         ! Parameter kmax from DIR file have the total numbers of
         ! the hybrid levels so we need subtract one to get the 
         ! BAM vertical levels

         BFile%kmax = BFile%kmax - 1

         !
         ! vertical coordinate info
         !

         BFile%idvc = 2
         BFile%idsl = 0

      else

         read (BFile%uHead,'(2A41)') BFile%jttl, BFile%specal
         !
         ! Na versão sigma do modelo esses parametros não estão
         ! no arquivo DIR. Além disso o GSI so foi testado com
         ! o modelo SSiB.
         !
         BFile%LSM = 1
         BFile%PBL = -1

         allocate(BFile%delSig(BFile%kmax),stat=iret)
!--------------------------------------------------------------------------------------!
!
!      read (BFile%udir,'(5E16.8)') (BFile%head%delSig(k),k=1,BFile%head%kmax)

         icol = 5
         isum = 0
         do while(isum .lt. BFile%kmax)
            icol = min(icol, (BFile%kmax-isum))
            READ (BFile%uHead,*) (BFile%delSig(k+isum),k=1,icol)
            isum = isum + icol
         enddo
!--------------------------------------------------------------------------------------!

         allocate(BFile%SigMid(BFile%kmax),stat=iret)
         allocate(BFile%SigInt(BFile%kmax+1),stat=iret)
   
         call setsig(&
                     BFile%SigInt, & ! sigma value at each interface.
                     BFile%SigMid, & ! sigma value at midpoint of each layer
                     BFile%DelSig  & ! sigma spacing for each layer.
                    )
         !
         ! vertical coordinate info
         !

         BFile%idvc = 1
         BFile%idsl = 1
         
      endif

      allocate(BFile%root,stat=iret)
      BFile%flds => BFile%root
      BFile%fcount = 0

      !
      ! pmark and LastSize are used to obtain start point to read a forecast
      ! BAM file. We assume that BAM file has a header of size iosize + 10*iosize + iosize bytes 
      ! (at a single precision iosize = 4). We start read before header!
      !

      BFile%iosize = iosize

      pmark     = 0
      LastSize  = iosize + 10*iosize + iosize
      
      if (trim(BFile%ftype) .eq. 'anl')then
         if (BFile%isHybrid)then
            LastSize = LastSize +  BFile%iosize * ( 2 * (BFile%kmax+1) )
         else
            LastSize = LastSize + BFile%iosize * ( 2* BFile%kmax + 1 )
         endif
            
      endif
      
      do
         read (BFile%uHead,'(A40,2X,A4,2X,I8,3X,I4,4X,I3)',&
               END=35,ERR=34,iostat=iret)           &
                                         BFile%flds%Name,   &
                                         BFile%flds%ProDia, &
                                         BFile%flds%nharm,  &
                                         BFile%flds%nlevs,  &
                                         BFile%flds%fldunit
         !
         ! Obtain Position marker to read directly a BAM forescast file.
         ! We assume that the integer value occupies a 32-bit (4-byte)
         ! for a single precision file and 64-bit (8-byte) for a double
         ! precion
         !
         
         pmark            = pmark + LastSize
         BFile%flds%pmark = pmark + (iosize + 1)
         
         LastSize         = (BFile%flds%nharm*BFile%iosize+BFile%iosize*2) * BFile%flds%nlevs
         
         ! allocate next field information
         allocate(BFile%flds%next,stat=iret)
         BFile%flds => BFile%flds%next

         BFile%fcount = BFile%fcount + 1
      enddo

34    write(*,'(A)') '!!! wrong header file '//trim(BFile%fHead),iret
      if(present(istat))then
         istat = iret
         return
      endif

35    continue

      !
      ! Init Grid and Spectral parameters
      !

      Mend = BFile%mMax-1
      call BFile%initTransform(Mend)

      !
      ! Lat/Lon informations
      !
      allocate(glats(BFile%yMax))

      call GetGaussianLatitudes(BFile%yMax,glats)

      allocate(BFile%rlat(BFile%yMax))

      BFile%rlat ( BFile%yMax:1:-1 ) = gLats ( 1:BFile%yMax )

      deallocate(glats)

      allocate(glons(BFile%xMax))
      call GetLongitudes(BFile%xMax,0.0_r8, glons)
      
      allocate(BFile%rlon(BFile%xMax))

      BFile%rlon ( 1:BFile%xMax ) = glons( 1:BFile%xMax )
      
      deallocate(glons)

      return 

40    write (*, '(A)') &
              ' Unexpected End of File in Input File Directory.'
      stop 4100

   end subroutine
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: BAM_Close_ - routine that close all BAM files opened and reset all
!                         informations in BAMFile structure.
!
! 
! !DESCRIPTION: Esta rotina fecha todos os arquivos do BAM (ANL, FCT, DIR) que
!               estejam abertos e também reinicia todas as informacões contidas
!               na estrutura de dados BAMFile.
!
!
! !INTERFACE:
!

   subroutine bamClose_(BFile, istat)

      implicit none
!
! !INPUT PARAMETERS:
!      
      ! BAM file structure
      class(BAMFile),     intent(inout) :: BFile
!
! !OUTPUT PARAMETERS:
!
      integer, optional, intent(  out) :: istat
! !SEE ALSO:
!
!   ResetHead_( ) reset all information from BAMFile structure
!
!
! !REVISION HISTORY: 
!  11 Oct 2016 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
      character(len=100), parameter :: myname_=':: BAM_Close( ... )'

      type(BAM_fld), pointer :: tmp => null()
      integer       :: iret
      logical       :: isopen

      if(present(istat)) istat = 0

      !
      ! close ANL file if is opened
      !

      if (associated(BFile%uBin))then
         inquire(BFile%uBin, opened=isopen)
      else
         isopen = .false.
      endif

      if(isopen) then 
         close(BFile%uBin, STATUS='KEEP', iostat=iret)
         if(iret.ne.0)then
            write(*,'(3(1x,A))')trim(myname_),': ERROR to close file',trim(BFile%fBin)
            if(present(istat)) istat = iret
         endif
         nullify(BFile%uBin)
         BFile%fBin = ''
      endif

      !
      ! close DIR file if is opened
      !
      if (associated(BFile%uHead))then
         inquire(BFile%uHead, opened=isopen)
      else
         isopen = .false.
      endif
      if(isopen) then 
         close(BFile%uHead, STATUS='KEEP', iostat=iret)
         if(iret.ne.0)then
            write(*,'(3(1x,A))')trim(myname_),': ERROR to close file',trim(BFile%fHead)
            if(present(istat)) istat = iret
         endif

         nullify(BFile%uHead)
         BFile%fHead = ''

         !
         ! reset all header data
         !

         call ResetHead_(BFile)

         !
         ! reset transfom data
         !

         call BFile%destroyTransform( )

         !
         ! reset vartable
         !

         BFile%fcount = -1

         tmp => BFile%root%next
         do
            deallocate(BFile%root)
            if(.not.associated(tmp)) exit
            BFile%root => tmp
            tmp => tmp%next         
         enddo

      endif

      return

   end subroutine
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: ResetHead_ - routine that reset all informations in BAMFile 
!                         structure.
!
! 
! !DESCRIPTION: Esta rotina reinicia todas as informacões contidas
!               na estrutura de dados BAMFile.
!
!
! !INTERFACE:
!

   subroutine ResetHead_(BFile)
      implicit none
      
!
! !INPUT PARAMETERS:
!      
      ! BAM file structure
      type(BAMFile), intent(inout) :: BFile

! !REVISION HISTORY: 
!  11 Oct 2016 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
      character(len=100), parameter :: myname_=':: ResetHead_()'
      integer :: iret

      BFile%ittl  = ''
      BFile%jttl  = ''
      BFile%sdain = ''
      BFile%dtin  = ''
      BFile%trunc = ''
      BFile%nexp  = -1
      BFile%mMax  = -1
      BFile%kmax  = -1
      BFile%ihr   = -1
      BFile%idy   = -1
      BFile%imo   = -1
      BFile%iyr   = -1
      BFile%fhr   = -1
      BFile%fdy   = -1
      BFile%fmo   = -1
      BFile%fyr   = -1

      if(associated(BFile%delsig))&
                deallocate(BFile%delsig)

      if(associated(BFile%SigInt))&
                deallocate(BFile%SigInt)

      if(associated(BFile%SigMid))&
                deallocate(BFile%SigMid)

      if(associated(BFile%rlon))&
                deallocate(BFile%rlon)

      if(associated(BFile%clon))&
                deallocate(BFile%clon)

      if(associated(BFile%slon))&
                deallocate(BFile%slon)

      if(associated(BFile%rlat))&
                deallocate(BFile%rlat)

      if(associated(BFile%clat))&
                deallocate(BFile%clat)

      if(associated(BFile%slat))&
                deallocate(BFile%slat)

      return

   end subroutine
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: BAM_GetOneDim - return information about one dimension of BAM file.
!
!
! 
! !DESCRIPTION: Rotina que retorma informacão sobre um das dimensões do arquivo
!               do modelo BAM. Podem ser retornadas informacões sobre:
!               
!               1. IMax : Número de pontos na direcão I (longitude)
!               2. JMax : Número de pontos na direcão J (latitude)
!               3. KMax : Número de pontos na direcão K (níveis verticais)
!               4. Mend : comprimento de onda (spectral)
!
! !INTERFACE:
!

   function bamGetOneDim_(BFile, DName) result(dim)
      
      implicit none
!
! !INPUT PARAMETERS:
!      
      ! BAM file structure
      class(BAMFile),     intent(in   ) :: BFile

      !BAM dimension name
      character(len=*),  intent(in   ) :: DName
      !
      ! Can be:
      !   1. imax - points in i direction (longitude)
      !   2. jmax - points in j direction (latitude)
      !   3. kmax - points in k direction (levels)
      !   4. mend - wave number (spectral)
      !
!
! !OUTPUT PARAMETERS:
!
      integer                          :: dim
!
! !REVISION HISTORY: 
!  11 Oct 2016 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!

      character(len=100), parameter :: myname_=':: BAM_GetOneDim( ... )'

      character(len=5) :: DimName

      DimName = Trim(Adjustl(lcase(DName)))

      select case (DimName)
         case ('imax','xmax')
            dim = BFile%xMax
         case ('jmax','ymax')
            dim = BFile%yMax
         case ('kmax')
            dim = BFile%KMax
         case ('mend')
            dim = BFile%mMax - 1
         case ('mnwv2')
            dim = BFile%MnWv2
         case ('mnwv3')
            dim = BFile%MnWv3
         case default
            write(*,'(4A)')trim(myname_),': wrong dimension <',trim(DimName),'>'
            write(*,'( A)')'try: imax, jmax, kmax, mend'
            dim = -1
      end select

   end function
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: BAM_GetTimeInfo - routine to return some information about the time
!                              of the BAM file 
!
! 
! !DESCRIPTION: Esta rotina retorna informacões sobre a data dos arquivos do
!               modelo BAM. Podem ser retornadas as seguintes informacões:
!               
!               1. ihr: Hora da condicão inicial utilizada para a simulacão
!               2. iyr: Ano da condicão inicial utilizada para a simulacão
!               3. idy: dia da condicão inicial utilizada para a simulacão
!               4. imo: mês da condicão inicial utilizada para a simulacão
!               5. fhr: hora da previsão da simulacão
!               6. fyr: ano da previsão da simulacão
!               7. fdy: dia da previsão da simulacão
!               8. fmo: mês da previsão da simulacão
!
! !INTERFACE:
!

   function bamGetTimeInfo_(BFile, DName) result(dt)
      
      implicit none
!
! !INPUT PARAMETERS:
!      
      ! BAM file structure
      class(BAMFile),     intent(in   ) :: BFile

      ! BAM time request
      character(len=*),  intent(in   ) :: DName
      !
      ! Can be:
      !   1. ihr: request hour of initial condition
      !   2. iyr: request year of initial condition
      !   3. idy: request day of initial condition
      !   4. imo: request month of initial condition
      !   5. fhr: request hour of forecast
      !   6. fyr: request year of forecast
      !   7. fdy: request day of forecast
      !   8. fmo: request month of forecast
!
! !OUTPUT PARAMETERS:
!
      ! time of simulation
      integer                          :: dt
!
! !REVISION HISTORY: 
!  11 Oct 2016 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
      character(len=100), parameter :: myname_=':: BAM_GetTimeInfo( ... )'

      character(len=4) :: DimName

      DimName = Trim(Adjustl(lcase(DName)))

      select case (DimName)
         case ('iyr')
            dt = BFile%iyr
         case ('imo')
            dt = BFile%imo
         case ('idy')
            dt = BFile%idy
         case ('ihr')
            dt = BFile%ihr
         case ('fyr')
            dt = BFile%fyr
         case ('fmo')
            dt = BFile%fmo
         case ('fdy')
            dt = BFile%fdy
         case ('fhr')
            dt = BFile%fhr

         case default
            write(*,'(4A)')trim(myname_),': wrong dimension <',trim(DimName),'>'
            write(*,'( A)')'try:'
            write(*,'( A)')' * ihr, idy, imo, iyr: Initial Time'
            write(*,'( A)')' * fhr, fdy, fmo, fyr: Forecast Time'
            write(*,'( A)')''
            dt = -1
      end select

   end function
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: BAM_GetDims - routine to return informations about all BAM file
!                          dimensions. 
!
! 
! !DESCRIPTION: Esta rotina retorna informacões sobre todas as dimensões do
!               modelo BAM. São retornadas as seguintes informacões:
!
!               1. IMax : Número de pontos na direcão I (longitude)
!               2. JMax : Número de pontos na direcão J (latitude)
!               3. KMax : Número de pontos na direcão K (níveis verticais)
!               4. Mend : comprimento de onda (spectral)

!
! !INTERFACE:
!

   subroutine bamGetDims_(BFile, Imax, Jmax, Kmax, Mend, istat)

      implicit none
!
! !INPUT PARAMETERS:
!      
      ! BAM file structure
      class(BAMFile),     intent(in   ) :: BFile

      !   1. imax - points in i direction (longitude)
      integer,           intent(  out) :: IMax

      !   2. jmax - points in j direction (latitude)
      integer,           intent(  out) :: JMax

      !   3. kmax - points in k direction (levels)
      integer,           intent(  out) :: KMax

      !   4. mend - wave number (spectral)
      integer,           intent(  out) :: Mend
!
! !OUTPUT PARAMETERS:
!
      integer, optional, intent(  out) :: istat
!
! !REVISION HISTORY: 
!  11 Oct 2016 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!

      character(len=100), parameter :: myname_=':: BAM_GetDims( ... )'
     
      if(present(istat)) istat = 0

      IMax  = BFile%xMax
      JMax  = BFile%yMax
      KMax  = BFile%KMax
      Mend  = BFile%mMax - 1

   end subroutine
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IFUNCTION: BAM_GetNLevels - return how many vertical levels has one
!                             give variable
!
! 
! !DESCRIPTION: Esta funcão retorna o número de níveis verticais que uma dada
!               variável possui.
!
!
!
! !INTERFACE:
!

   function bamGetNlevels_(BFile,VName,istat) result(nlevs)

      implicit none

!
! !INPUT PARAMETERS:
!      
      ! BAM file structure
      class(BAMFile),     intent(in   ) :: BFile

      ! Name of a BAM variable
      character(len=*),  intent(in   ) :: VName
      
!
! !OUTPUT PARAMETERS:
!
      integer, optional, intent(  out) :: istat
!
! !RETURN VALUE:
!
      integer                          :: NLevs
!
! !REVISION HISTORY: 
!  11 Oct 2016 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
      character(len=100), parameter :: myname_=':: BAM_GetNLevels( ... )'
     

      type(BAM_fld), pointer :: fld
      integer                :: f
      character(len=strlen)  :: InqName
      character(len=strlen)  :: FLDName

      if(present(istat)) istat = 0


      InqName = trim(adjustl(lcase(VName)))

      NLevs = -1

      fld => BFile%root
      do f=1,BFile%fcount

         FLDName = trim(adjustl(lcase(fld%name)))

         if(trim(InqName) .eq. trim(FLDName))then
            NLevs = fld%nlevs
            if(present(istat)) istat = 0
            return
         endif

         fld => fld%next
      enddo

      if(NLevs .eq. -1 )then
         write(*,'(4A)')trim(myname_),': variable not found ! <',trim(InqName),trim(FLDName)
         if(present(istat)) istat = -1
      endif

      return

   end function
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: GetOFS_ - return Original field Size.
!
!
! 
! !DESCRIPTION: Rotina que retorma o tamanho original do campo do BAM.
!               Se o campo estiver no espaço espectral retornará 
!               (Mend+1)*(Mend+2), caso esteja no espaço físico, retornará
!               IMax*JMax.
!               
!
! !INTERFACE:
!

   function bamGetOFS_(BFile,VName,istat) result(OFSize)

      implicit none

!
! !INPUT PARAMETERS:
!      
      ! BAM file structure
      class(BAMFile),     intent(in   ) :: BFile

      ! Name of a BAM variable
      character(len=*),  intent(in   ) :: VName
      
!
! !OUTPUT PARAMETERS:
!
      integer, optional, intent(  out) :: istat
!
! !RETURN VALUE:
!
      integer                          :: OFSize
!
! !REVISION HISTORY: 
!  15 Jun 2018 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
      character(len=100), parameter :: myname_=':: GetOFS_( ... )'
     

      type(BAM_fld), pointer :: fld
      integer                :: f
      character(len=strlen)  :: InqName
      character(len=strlen)  :: FLDName

      if(present(istat)) istat = 0


      InqName = trim(adjustl(lcase(VName)))

      OFSize = -1

      fld => BFile%root
      do f=1,BFile%fcount

         FLDName = trim(adjustl(lcase(fld%name)))

         if(trim(InqName) .eq. trim(FLDName))then
            OFSize = fld%nharm
            if(present(istat)) istat = 0
            return
         endif

         fld => fld%next
      enddo

      if(OFSize .eq. -1 )then
         write(*,'(4A)')trim(myname_),': variable not found ! <',trim(InqName),trim(FLDName)
         if(present(istat)) istat = -1
      endif

      return

   end function
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: bamGetPhysics - return physical scheme used by BAM model
!
! 
! !DESCRIPTION: Esta rotina retorna um valor inteiro referente à opção física
!               utilizada pelo modelo BAM. As seguintes informacões são disponíveis:
!
!              1. lsm : esquema de superfície utillizado (1-SSiB, 2-SiB2, 3-IBIS)
!              2. pbl : esquema de camada linite utilizao ()
!
! !INTERFACE:
!

   function bamGetPhysics_(BFile, what)result(scheme)
!
! !INPUT PARAMETERS:
!      
      ! BAM file structure
      class(BAMFile)   :: BFile

      ! Name of inquired physic type scheme
      character(len=*) :: what
      ! 
      ! Can be:     
      !       1. lsm : esquema de superfície utillizado
      !       2. pbl : esquema de camada linite utilizao

!
! !OUTPUT PARAMETERS:
!
      integer(i4) :: scheme
!
! !REVISION HISTORY: 
!  25 Mar 2021 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
      character(len=100) :: inquired

      character(len=100), parameter :: myname_=':: bamGetPhysics_( ... )'

      inquired = lcase(what)

      select case(trim(inquired))
         case('pbl')
            scheme = BFile%pbl
         case('lsm')
            scheme = BFile%lsm
         case('idvc')
            scheme = BFile%idvc
         case('idsl')
            scheme = BFile%idsl
         case default
            write(*,'(4A)')trim(myname_),': physics type not yeat implemented ! <',trim(Inquired),'>'
            scheme = -1
      end select
   end function
!
!EOC
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: BAM_GetWCoord - routine to return informations world coordinates
!                            of BAM model file
!
! 
! !DESCRIPTION: Esta rotina retorna um vetor contendo o valor de uma determinada
!               informacão geográfica do modelo BAM. Podem ser retornadas as 
!               seguintes informacões:
!
!              1. rlon : longitude real
!              2. slon : seno da longitude real
!              3. clon : cosseno da longitude real
!              4. rlat : latitude real
!              5. slat : seno da latitude real
!              6. clat : cosseno da latitude real
!
! !INTERFACE:
!

   subroutine  bamGetWCoord_(BFile, what, coord, istat)
      implicit none
!
! !INPUT PARAMETERS:
!      
      ! BAM file structure
      class(BAMFile),     intent(in   ) :: BFile

      ! Name of inquired Word Coordinate
      character(len=*),  intent(in   ) :: what
      ! 
      ! Can be:
      !        1. rlon : real longitude
      !        2. slon : sine of real longitude
      !        3. clon : cosine of real longitude
      !        4. rlat : real latitude
      !        5. slat : sine of real latitude
      !        6. clat : cosine of real latitude
      !
!
! !OUTPUT PARAMETERS:
!
      real(r4),          intent(inout) :: coord(:)
      integer, optional, intent(  out) :: istat
!
! !REVISION HISTORY: 
!  11 Oct 2016 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!

      character(len=100), parameter :: myname_=':: BAM_GetWCoord( ... )'

      integer :: lenIn
      integer :: lenOu

      select case (trim(what))
         case('rlon') ! real longitude

            !check consistensy
            lenIn = size(coord)
            lenOu = BFile%xMax
            if (lenIn.ne.lenOu) then
               write (stdout,'(2(1x,A))')trim(myname_),':ERROR, wrong size of coordinate array:'
               write (stdout,'(I4,1x,A,1x,I4)')lenIn,' .ne. ',lenOu
               return
            endif

            coord = BFile%rlon

         case('slon') ! sine of longitude

            !check consistensy
            lenIn = size(coord)
            lenOu = BFile%xMax
            if (lenIn.ne.lenOu) then
               write (stdout,'(2(1x,A))')trim(myname_),':ERROR, wrong size of coordinate array:'
               write (stdout,'(I4,1x,A,1x,I4)')lenIn,' .ne. ',lenOu
               return
            endif

            coord = sin(BFile%rlon/rd)

         case('clon') !cosine of longitude

            !check consistensy
            lenIn = size(coord)
            lenOu = BFile%xMax
            if (lenIn.ne.lenOu) then
               write (stdout,'(2(1x,A))')trim(myname_),':ERROR, wrong size of coordinate array:'
               write (stdout,'(I4,1x,A,1x,I4)')lenIn,' .ne. ',lenOu
               return
            endif

            coord = cos(BFile%rlon/rd)

         case('rlat') ! real latitude

            !check consistensy
            lenIn = size(coord)
            lenOu = BFile%yMax
            if (lenIn.ne.lenOu) then
               write (stdout,'(2(1x,A))')trim(myname_),':ERROR, wrong size of coordinate array:'
               write (stdout,'(I4,1x,A,1x,I4)')lenIn,' .ne. ',lenOu
               return
            endif

            coord = BFile%rlat

         case('slat') ! sine of latitude

            !check consistensy
            lenIn = size(coord)
            lenOu = BFile%yMax
            if (lenIn.ne.lenOu) then
               write (stdout,'(2(1x,A))')trim(myname_),':ERROR, wrong size of coordinate array:'
               write (stdout,'(I4,1x,A,1x,I4)')lenIn,' .ne. ',lenOu
               return
            endif

            coord = sin(BFile%rlat/rd)

         case('clat') ! cosine of latitude

            !check consistensy
            lenIn = size(coord)
            lenOu = BFile%yMax
            if (lenIn.ne.lenOu) then
               write (stdout,'(2(1x,A))')trim(myname_),':ERROR, wrong size of coordinate array:'
               write (stdout,'(I4,1x,A,1x,I4)')lenIn,' .ne. ',lenOu
               return
            endif

            coord = cos(BFile%rlat/rd)

         case default
            write(stdout,'(2(1x,A))')trim(myname_),': ERROR'
            write(stdout,'(2(1x,A))')'Wrong coordinate inquired :',trim(what)
            write(stdout,'(2(1x,A))')'validy coordinates:'
            write(stdout,'(2(1x,A))')' '
            write(stdout,'(2(1x,A))')'rlon : real longitude'
            write(stdout,'(2(1x,A))')'clon : cosine of real longitude'
            write(stdout,'(2(1x,A))')'slon : sine of real longitude'
            write(stdout,'(2(1x,A))')' '
            write(stdout,'(2(1x,A))')'rlat : real latitude'
            write(stdout,'(2(1x,A))')'clat : cosine of real latitude'
            write(stdout,'(2(1x,A))')'slat : sine of real latitude'

            if(present(istat))istat = -1
            return 

      end select

   endsubroutine
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: BAM_GetVerticalCoord - return informations vertical coordinates
!                            of BAM model.
!
! 
! !DESCRIPTION: Esta rotina retorna um vetor contendo informacões sobre a
!               coordenada vertical do modelo BAM.
!               Podem ser retornadas as seguintes informacões:
!              1. Ak :
!              2. Bk :
!              3. Ck :
!              4. SI :
!              5. SL :
!              6. DL :
!
!
! !INTERFACE:
!

   subroutine bamGetVerticalCoord_(BFile, what, vcoord, istat)
      implicit none
!
! !INPUT PARAMETERS:
!      
      ! BAM file structure
      class(BAMFile),     intent(in   ) :: BFile

      ! Name of inquired Vertical Coordinate parameter
      character(len=*),  intent(in   ) :: what
!
! !OUTPUT PARAMETERS:
!   
      real(r8),          intent(inout) :: vcoord(:)
      integer, optional, intent(  out) :: istat

! !REMARKS:
!
!  Pressure is defined as:
!       \begin{equation}
!          p_{(i,j,k)} = A_{k}P_{0}+B_{k}P_{s}(i,j)
!       \end{equation}
!
!  where $p$ is the pressure at a given level and latitude, longitude grid point. 
!  The coefficients $A$, $B$ and $P_{0}$ are constants. $P_{s}$ is the model's current 
!  surface pressure. $P_{0}$ is set in the model code. The input model initial 
!  conditions dataset sets $A$ and $B$ through the variables hyam, hyai, hybm, and hybi. 
!  The subscript "i" refers to interface levels, and "m" refers to the mid-point levels. 
!  "hyam" then refers to Hybrid level "A" coefficient on the interfaces. 
!
!  More details on the theoretical nature of the vertical coordinate system can 
!  be found in Collins et al. [2004]. 
!
!  Collins, W. D., P. J. Rasch, and Others, Description of the NCAR Community 
!  Atmosphere Model (CAM 3.0), Technical Report NCAR/TN-464+STR, National Center 
!  for Atmospheric Research, Boulder, Colorado, 210 pp., 2004. 
!
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

      character(len=100), parameter :: myname_=':: BAM_GetWCoord( ... )'

      character(len=2) :: vcn
   
      if(present(istat)) istat = 0

      vcn = ucase(trim(what))

      !Sanity Check

      select case (trim(vcn) )
         case ('AK', 'BK', 'CK')
            if ( .not. BFile%isHybrid)then
               write(stdout,'(2(1x,A))')trim(myname_),': ak, bk or ck are just for Hybrid-sigma vertical coordinate'
               if(present(istat))istat = -2
               vcoord = undef
               return
            endif

         case ('SI', 'SL', 'DL')
            if (BFile%isHybrid)then
               write(stdout,'(2(1x,A))')trim(myname_),': si, sl or dl are just for Sigma vertical coordinate'
               if(present(istat))istat = -2
               vcoord = undef
               return
            endif
         case default
            write (stdout,'(3(1x,A))')trim(myname_),': unknown requested parameter,:',trim(vcn)
            if(present(istat))istat = -1
      end select

      select case (trim(vcn))
         case ( 'AK' )
            vcoord = BFile%ak
         case ( 'BK' )
            vcoord = BFile%bk
         case ( 'CK' )
            vcoord = BFile%ck
         case ( 'SI' )
            vcoord = BFile%SigInt
         case ( 'SL' )
            vcoord = BFile%SigMid
         case ( 'DL' )
            vcoord = BFile%DelSig
         case default
            write (stdout,'(3(1x,A))')trim(myname_),': unknown requested parameter,:',trim(vcn)
            if(present(istat))istat = -1
      end select

      return
   end subroutine
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: BAM_GetSigValues - return vertical sigma parameters
!
! 
! !DESCRIPTION: Esta rotina retorna os valores dos seguintes parametros da
!               coordenada vertical sigma:
!
!              1. DelSig : diferenca entre os niveis sigma
!              2. SigInt : interface entre os niveis sigma
!              3. SigMid : ponto médio entre os níveis sigma
!
! !INTERFACE:
!

   subroutine bamGetSigValues_(BFile, DelSig, SigInt, SigMid, istat)
      implicit none
!
! !INPUT PARAMETERS:
!      
      ! BAM file structure
      class(BAMFile),  intent(in   ) :: BFile
!
! !OUTPUT PARAMETERS:
!
      real(r8),          intent(inout) :: DelSig(:)
      real(r8),          intent(inout) :: SigInt(:)
      real(r8),          intent(inout) :: SigMid(:)
      integer, optional, intent(  out) :: istat
      
! !REVISION HISTORY: 
!
!  11 Oct 2016 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
      character(len=100), parameter :: myname_=':: BAM_GetSigmaValues( ... )'

      if(present(istat)) istat = 0
      
      DelSig = BFile%DelSig
      SigInt = BFile%SigInt
      SigMid = BFile%SigMid


   end subroutine
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IFUNCTION: bamGetVarNames_ - return a table with all variable names
!
! 
! !DESCRIPTION: Esta funcão retorna uma tabela contendo todas as variáveis
!               disponíveis no arquivo do modelo bam
!
!
! !INTERFACE:
!
subroutine  bamGetVarNames_(self,table)

!
! !INPUT PARAMETERS:
!
   class(BAMFile), intent(in) :: self

!
! ! OUTPUT PARAMETERS:
!
   character(len=*), allocatable :: table(:)
!
! !REVISION HISTORY: 
!  28 Apr 2021 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
   character(len=100), parameter :: myname_=':: bamGetVarNames_( ... )'

   type(BAM_fld), pointer :: fields => null()
   integer       :: i

   allocate(table(self%fcount))
   fields => self%root
   do i=1,self%fcount
      table(i) = trim(fields%name)
      fields => fields%next
   enddo
   return

end subroutine
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IFUNCTION: FindFLD_ - find position of a forecasted field in the BAM file
!
! 
! !DESCRIPTION: Esta funcão retorna a posicão de um determinado campo previsto
!               pelo modelo BAM indicado pelo arquivo DIR e assim poder acessa-lo
!               no arquivo FCT.
!
!
! !INTERFACE:
!

   function FindFLD_(BFile, wfld, wlev, iosize, fosize) result(idx)
      implicit none
!
! !INPUT PARAMETERS:
!     

      type(BAMFile),     intent(in   ) :: BFile ! file information
      character (len=*), intent(in   ) :: wfld  ! what field name ?
      integer,           intent(in   ) :: wlev  ! what field level ? 
      integer,           intent(in   ) :: iosize! size of record markers
!
! !RETURN VALUE:
!
      integer, intent(out) :: fosize
      integer :: idx ! Field position
!
! !REVISION HISTORY: 
!  11 Oct 2016 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
      character(len=100), parameter :: myname_=':: FindFLD_( ... )'
      
      character(len=strlen)  :: InqName
      character(len=strlen)  :: FldName

      type(BAM_fld), pointer :: fld
      integer :: ilev
      integer :: ifld

!      real, allocatable :: tst(:)
!      integer           :: isize, osize, i, tmpsize

      !
      ! Get root file information
      !

      fld => BFile%root

      !
      ! BAM fct files contain spectral an grid point data.
      ! This data are stored in a sequential unformatted 
      ! format, so we need to access (through) the whole file!
      !

      InqName = Trim(adjustl(lcase(wfld)))

!      idx  = 1
      ! Here we need skip header from BAM fct file
      ! to get initial position of first field!
      ! Head have a 10 positions!
      !
      ! NOTA: Fortran typically puts a size header at the 
      !       start and end of the data, in this case it 
      !       is using 4 byte record markers. So we need
      !       add 8 bytes to idx
      !       Add iosize is the size of record markers

      idx = 10*BFile%iosize + BFile%iosize*2
      if (trim(BFile%ftype) .eq. 'anl')then
         idx = idx + BFile%iosize * (2*BFile%kmax + 1)
      endif
      !
      ! ----------------------------------------------
      ilev = 0
          
      do ifld=1,BFile%fcount

         FldName = Trim(Adjustl(lcase(fld%name)))

         if ( trim(InqName) .eq. trim(FldName) )then

            if(wlev .gt. fld%nlevs)then
               write(stdout,*)trim(myname_),': error, wlev > nlevs :', wlev, fld%nlevs
               return
            endif


             idx = idx + (wlev-1)*(fld%nharm*BFile%iosize + BFile%iosize*2)
             fosize = fld%nharm

            exit

         else

            idx = idx + (fld%nharm*BFile%iosize+BFile%iosize*2)*fld%nlevs

         endif

         fld => fld%next
      enddo

      idx = idx + (BFile%iosize + 1)

   end function
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: SGetFLD_1D - return a single precision 1D array with a requested
!                         forecast BAM field.
!
! 
! !DESCRIPTION: Esta rotina retorna um vetor 1D em precisão simples contendo o
!               campo previsto pelo modelo BAM que foi requisitado.
!
! !INTERFACE:
!

   subroutine SGetFLD_1D(self, wfld, wlev, grd, istat)

      implicit none
!
! !INPUT PARAMETERS:
! 
      class(BAMFile),    intent(in   ) :: self ! file information
      character (len=*), intent(in   ) :: wfld  ! what field name ?
      integer,           intent(in   ) :: wlev  ! what field level ? 
!
! !OUTPUT PARAMETERS:
! 
      real(r4),          intent(inout) :: grd(:)
      integer, optional, intent(  out) :: istat
!
! !SEE ALSO:
!
!    DGetFLD_1D( ) - return a double precision 1D BAM Field
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
      character(len=*), parameter      :: myname_ = ':: SGetFLD_1D( ... )'

      real(r8) :: TmpGrd(size(grd,1))
      integer  :: iret

      call self%DGetFLD_1D(wfld, wlev, TmpGrd, iret)
      if(present(istat))istat = iret

      grd = real(TmpGrd,r4)

      return

   end subroutine
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: SGetFLD_2D - return a single precision 2D array with a requested
!                         forecast BAM field.
!
! 
! !DESCRIPTION: Esta rotina retorna um vetor 2D em precisão simples contendo o
!               campo previsto pelo modelo BAM que foi requisitado.
!
! !INTERFACE:
!

   subroutine SGetFLD_2D(self, wfld, wlev, grd, istat)

      implicit none
!
! !INPUT PARAMETERS:
! 
      class(BAMFile),    intent(in   ) :: self  ! file information
      character (len=*), intent(in   ) :: wfld  ! what field name ?
      integer,           intent(in   ) :: wlev  ! what field level ? 
!
! !OUTPUT PARAMETERS:
! 
      real(r4),          intent(inout) :: grd(:,:)
      integer, optional, intent(  out) :: istat
!
! !SEE ALSO:
!
!    DGetFLD_2D( ) - return a double precision 2D BAM Field
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
      character(len=*), parameter      :: myname_ = ':: SGetFLD_2D( ... )'

      
      real(r8) :: TmpGrd(size(grd,1),size(grd,2))
      integer  :: iret

      call self%DGetFLD_2D(wfld, wlev, TmpGrd, iret)
      if(present(istat))istat = iret

      grd = real(TmpGrd,r4)

      return

   end subroutine
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: DGetFLD_2D - return a double precision 2D array with a requested
!                         forecast BAM field.
!
! 
! !DESCRIPTION: Esta rotina retorna um vetor 2D em precisão dupla contendo o
!               campo previsto pelo modelo BAM que foi requisitado.
!
! !INTERFACE:
!

   subroutine DGetFLD_2D(self, wfld, wlev, grd, istat)

      implicit none
!
! !INPUT PARAMETERS:
! 
      class(BAMFile),    intent(in   ) :: self  ! file information
      character (len=*), intent(in   ) :: wfld  ! what field name ?
      integer,           intent(in   ) :: wlev  ! what field level ? 
!
! !OUTPUT PARAMETERS:
! 
      real(r8),          intent(inout) :: grd(:,:)
      integer, optional, intent(  out) :: istat
!
! !SEE ALSO:
!
!    DGetFLD_2D( ) - return a double precision 1D BAM Field
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
      character(len=*), parameter      :: myname_ = ':: DGetFLD_2D( ... )'

      real(r8), allocatable :: fld (:)
      integer :: i, j, k
      integer :: IMax, JMax
      integer :: npts
      integer :: iret

      npts = size(grd)

      allocate(fld(npts))

      call self%DGetFLD_1D(wfld, wlev, fld, iret)
      if(present(istat)) istat = iret

      IMax = size(grd,1)
      JMax = size(grd,2)

      k = 1
      do j = 1,JMax
         do i = 1, IMax
            grd(i,j) = fld(k)
            k = k + 1
         enddo
      enddo
      
   end subroutine
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: SGetFLD_1D - return a double precision 1D array with a requested
!                         forecast BAM field.
!
! 
! !DESCRIPTION: Esta rotina retorna um vetor 1D em precisão dupla contendo o
!               campo previsto pelo modelo BAM que foi requisitado.
!
! !INTERFACE:
!     
   subroutine DGetFLD_1D(self, wfld, wlev, grd, istat)

      implicit none
!
! !INPUT PARAMETERS:
! 
      class(BAMFile),    intent(in   ) :: self  ! file information
      character (len=*), intent(in   ) :: wfld  ! what field name ?
      integer,           intent(in   ) :: wlev  ! what field level ? 
!
! !OUTPUT PARAMETERS:
! 
      real(r8),          intent(inout) :: grd(:)
      integer, optional, intent(  out) :: istat
!
! !SEE ALSO:
!
!    BAM_GetDims( ) - return BAM field dimensions
!    FindFLD_( ) - return BAM field position
!    ReadField( ) - read a BAM field in double or sigle precision
!    InitRecomposition( ) - initialize spectral2grid routines
!    RecompositionScalar( ) - convert spectral2grid fields
!    ClsMemRecomposition( ) - finalize spectral2grid routines
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
      character(len=*), parameter      :: myname_ = ':: DGetFLD_1D( ... )'

      integer :: ifld
      integer :: iret

      type(BAM_fld), pointer :: fld

      real (R8), dimension (:),   allocatable :: bufx
      real (R4), dimension (:),   allocatable :: bufs
      real (R8), dimension (:,:), allocatable :: buft
      real (R8), dimension (:,:), allocatable :: bufw

      integer :: iMax
      integer :: jMax
      integer :: Mend
      integer :: Mnwv2

      integer :: i, j, k
      integer(I8) :: idx
      integer :: npts

      character(len=strlen)  :: InqName
      character(len=strlen)  :: FldName

      character(len=80) :: linha
      character(len=1024) :: acc, fmt

      if(present(istat)) istat = 0

      if (.not. associated(self%uBin))then
         write(stdout,*)''
         write(stdout,*)'!!! WARNING !!!'
         write(stdout,*)'Binary file not opened!'
         write(stdout,*)'You can not get a field!'
         write(stdout,*)''
         write(stdout,*)'Review the options of the model OPEN statment!'
         write(stdout,*)''
         write(stdout,*)'Binary File: <empty>'
         write(stdout,*)'Header File:',trim(self%fHead)
         write(stdout,*)'Open mode  :',trim(self%oMode)
         write(stdout,*)''


         if(present(istat)) istat = -1
         return
      endif

      !
      ! get BAM file dimensions
      !

      Mend  = self%GetOneDim('Mend')
      iMax  = self%GetOneDim('iMax')
      jMax  = self%GetOneDim('jMax')
      MnWv2 = self%GetOneDim('MnWv2')

      !
      ! back to begin of file
      !

      rewind(self%uBin) 


      !
      ! Find BAM Field and get bit position and real size
      !
      
      idx = -1
      InqName = Trim(adjustl(lcase(wfld)))

      fld => self%root
      do ifld = 1, self%Fcount

         FldName = Trim(Adjustl(lcase(fld%name)))

         if ( trim(InqName) .eq. trim(FldName) )then

            if(wlev .gt. fld%nlevs)then
               write(stdout,*)trim(myname_),': error, wlev > nlevs :', wlev, fld%nlevs
               return
            endif

            idx = fld%pmark + ( (wlev-1)*(fld%nharm*R4 + R4*2))
            
            exit

         endif

         fld => fld%next
      enddo

      if (idx .lt. 0)then
         call perr(trim(myname_),'Field Not Found ...['//trim(self%fBin)//' ]')
         if(present(istat))istat = idx
         return
      endif

      !
      ! Get Field
      !
      allocate(bufs(fld%nharm))
      call ReadField (self%uBin, bufs, idx, iret)

      if (iret.ne.0)then
         if(present(istat)) istat = iret
         write(stdout,*)trim(myname_),': error to read file ::','['//trim(self%fBin)//' ], ',iret
      endif

      npts = size(grd)

      if ( ( npts .eq. (IMax*JMax) ) .and. (fld%nharm .eq. mnwv2 )) then ! rever este if

         !
         ! field is spectral but need return grid point
         !

         allocate (bufx (mnwv2))

         bufx(1:mnwv2) = real(bufs(1:mnwv2),R8)

         deallocate(bufs)

         allocate (bufw (imax,jmax))
         call self%Spec2Grid(bufx,bufw)
         
         k = 1
         do j = 1, JMax
            do i = 1, IMax
               grd(k) = bufw(i,j)
               k = k + 1
            enddo
         enddo

         deallocate(bufw)

      else

         grd = real(bufs,r8)

      endif

      return

   end subroutine
!
!EOC
!
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
   subroutine SRead1D_(funit, fld, idx, istat)
      implicit none
!
! !INPUT PARAMETERS:
! 
      integer(i4),           intent(in   ) :: funit
      integer(i8), optional, intent(in   ) :: idx
!
! !OUTPUT PARAMETERS:
! 
      real(r4),              intent(  out) :: fld(:)
      integer(i4), optional, intent(  out) :: istat
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
      character(len=*), parameter          :: myname_ = ':: SRead1D_( ... )'

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

   end subroutine
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: SRead2D_ - read and return a 2D sigle precision BAM field.
!
! 
! !DESCRIPTION: Esta rotina lê e retorna um campo modelo BAM em um matriz 2D em
!               precisão simples
!
! !INTERFACE:
!   

   subroutine SRead2D_(funit, fld, idx, istat)
      implicit none
!
! !INPUT PARAMETERS:
! 
      integer(i4),           intent(in   ) :: funit
      integer(i8), optional, intent(in   ) :: idx
!
! !OUTPUT PARAMETERS:
! 
      real(r4),              intent(  out) :: fld(:,:)
      integer(i4), optional, intent(  out) :: istat
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
      character(len=*), parameter      :: myname_ = ':: SRead2D_( ... )'

      integer :: iret

      if(present(istat))istat=0

      if(present(idx))then
         read(unit = funit, POS=idx, iostat = iret) fld
      else
         read(unit = funit, iostat = iret) fld
      endif

      if (iret.ne.0)then
         write(stdout,*)trim(myname_),': error to read field, ',iret
         if(present(istat)) istat = iret
      endif
      
   end subroutine
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: SRead3D_ - read and return a 3D sigle precision BAM field.
!
! 
! !DESCRIPTION: Esta rotina lê e retorna um campo modelo BAM em um vetor 3D em
!               precisão simples
!
! !INTERFACE:
!   
   subroutine SRead3D_(funit, fld, idx, istat)
      implicit none
!
! !INPUT PARAMETERS:
! 
      integer(i4),           intent(in   ) :: funit
      integer(i8), optional, intent(in   ) :: idx
!
! !OUTPUT PARAMETERS:
! 
      real(r4),              intent(  out) :: fld(:,:,:)
      integer(i4), optional, intent(  out) :: istat
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
      character(len=*), parameter      :: myname_ = ':: SRead3D_( ... )'

      integer :: iret

      if(present(istat))istat=0

      if(present(idx))then
         read(unit = funit, POS=idx, iostat = iret) fld
      else
         read(unit = funit, iostat = iret) fld
      endif
      if (iret.ne.0)then
         write(stdout,*)trim(myname_),': error to read file, ',iret
         if(present(istat)) istat = iret
      endif

   end subroutine
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: DRead1D_ - read and return a 1D double precision BAM field.
!
! 
! !DESCRIPTION: Esta rotina lê e retorna um campo modelo BAM em um vetor 1D em
!               precisão dupla.
!
! !INTERFACE:
!   
   subroutine DRead1D_(funit, fld, idx, istat)
      implicit none
!
! !INPUT PARAMETERS:
! 
      integer(i4),           intent(in   ) :: funit
      integer(i8), optional, intent(in   ) :: idx
!
! !OUTPUT PARAMETERS:
! 
      real(r8),              intent(  out) :: fld(:)
      integer(i4), optional, intent(  out) :: istat
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
      character(len=*), parameter          :: myname_ = ':: SRead1D_( ... )'

      integer :: iret

      if(present(istat)) istat = 0

      if(present(idx))then
         read(unit = funit, POS=idx, iostat = iret) fld
      else
         read(unit = funit, iostat = iret) fld
      endif

      if (iret.ne.0)then
         write(stdout,*)trim(myname_),': error to read file, ',iret
         if(present(istat)) istat = iret
      endif

   end subroutine
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: DRead2D_ - read and return a 2D double precision BAM field.
!
! 
! !DESCRIPTION: Esta rotina lê e retorna um campo modelo BAM em um matriz 2D em
!               precisão dupla.
!
! !INTERFACE:
!   

   subroutine DRead2D_(funit, fld, idx, istat)
      implicit none
!
! !INPUT PARAMETERS:
! 
      integer(i4),           intent(in   ) :: funit
      integer(i8), optional, intent(in   ) :: idx

!
! !OUTPUT PARAMETERS:
! 
      real(r8),              intent(  out) :: fld(:,:)
      integer(i4), optional, intent(  out) :: istat
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
      character(len=*), parameter      :: myname_ = ':: DRead2D_( ... )'

      integer :: iret

      if(present(istat))istat=0

      if(present(idx))then
         read(unit = funit, POS=idx, iostat = iret) fld
      else
         read(unit = funit, iostat = iret) fld
      endif

      if (iret.ne.0)then
         write(stdout,*)trim(myname_),': error to read file, ',iret
         if(present(istat)) istat = iret
      endif
      
   end subroutine
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: DRead3D_ - read and return a 3D double precision BAM field.
!
! 
! !DESCRIPTION: Esta rotina lê e retorna um campo modelo BAM em um matriz 3D em
!               precisão dupla.
!
! !INTERFACE:
!   
   subroutine DRead3D_(funit, fld, idx, istat)
      implicit none
!
! !INPUT PARAMETERS:
! 
      integer(i4),           intent(in   ) :: funit
      integer(i8), optional, intent(in   ) :: idx
!
! !OUTPUT PARAMETERS:
! 
      real(r8),              intent(  out) :: fld(:,:,:)
      integer(i4), optional, intent(  out) :: istat
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
      character(len=*), parameter      :: myname_ = ':: DRead3D_( ... )'

      integer :: iret

      if(present(istat))istat=0

      if(present(idx))then
         read(unit = funit, POS=idx, iostat = iret) fld
      else
         read(unit = funit, iostat = iret) fld
      endif

      if (iret.ne.0)then
         write(stdout,*)trim(myname_),': error to read file, ',iret
         if(present(istat)) istat = iret
      endif

   end subroutine
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: SkipRead_ - skip a sequential read in a BAM field.
!
! 
! !DESCRIPTION: Esta rotina pula a leitura sequencial de um campo modelo BAM.
!
!
! !INTERFACE:
!   

   subroutine SkipRead_(funit, istat)
      implicit none
!
! !INPUT PARAMETERS:
! 
      integer(i4),           intent(in   ) :: funit
!
! !OUTPUT PARAMETERS:
! 
      integer(i4), optional, intent(  out) :: istat
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
      character(len=*), parameter      :: myname_ = 'SkipRead_( ... )'

      integer :: iret

      if(present(istat))istat=0

      read(unit = funit, iostat = iret)
      if (iret.ne.0)then
         write(stdout,*)trim(myname_),': error to skip read file, ',iret
         if(present(istat)) istat = iret
      endif

   end subroutine
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: BAM_WriteAnlHeader - write a Header information in a Initial
!                                 condition BAM file
!
!
! !DESCRIPTION: Esta rotina escreve o cabecalho do arquivo de condicão inicial
!               do modelo BAM.
!
! !INTERFACE:

  subroutine bamWriteAnlHeader_(self, mydate, istat)

     implicit none
!
! !INPUT PARAMETERS:
!
     class(BAMFile),    intent(in   ) :: self      ! BAM files data type
     integer,           intent(in   ) :: mydate(:) ! Anl date
     
!
! !OUTPUT PARAMETERS:
!
     integer, optional, intent(  out) :: istat ! code error status

! !REVISION HISTORY:
! 	03 Aug 2016 - J. G. de Mattos - Initial code
!
!EOP
!-----------------------------------------------------------------------
!BOC
     character(len=100), parameter :: myname_=':: BAM_WriteAnlHeader( ... )'

     integer :: iret

     integer (i4) :: idate (4)
     integer (i4) :: idatec(4)
     integer (i4) :: ifday
     real (r4)    :: tod

     integer (i4) :: KMax
     integer (i4) :: KMaxP

     integer(i4) :: idsl
     integer(i4) :: idvc 


     real (r8), allocatable :: del(:)
     real (r8), allocatable :: ak(:)
     real (r8), allocatable :: bk(:)

!     integer :: k
!     real :: sumdel

    ! Get forecast day, time of day, dates

    ifday = 0
    tod   = 0.0

    ! -------------------
    ! idate(1) - hour
    ! idate(2) - month
    ! idate(3) - day
    ! idate(4) - year

    idate(4) = mydate(1)
    idate(2) = mydate(2)
    idate(3) = mydate(3)
    idate(1) = mydate(4)
    idatec   = idate

    !
    ! Get sigma levels
    !
    KMax  = self%GetOneDim('KMax')
    KMaxP = KMax + 1

    if (self%isHybrid)then
       allocate( ak (KMaxP) )
       allocate( bk (KMaxP) )
       call self%getVerticalCoord('ak',ak)
       call self%getVerticalCoord('bk',bk)
    else
!       allocate( del( KMax) )
       allocate( ak (KMaxP) )
       allocate( bk ( KMax) )
!       call self%GetSigValues (del, ak, bk, iret)
       call self%getVerticalCoord('si',ak)
       call self%getVerticalCoord('sl',bk)
    endif

    idvc = self%getPhysics('idvc')
    idsl = self%getPhysics('idsl')


!    do k=1,kmaxp
!       write (6,'(a,i2,1(a,f10.6))') &
!            ' level = ', k, ' si  = ', real(si(k),r4)
!    enddo
!    write (6,'(a)')' '
!    sumdel = 0.0
!    do k=1,kmax
!       sumdel = sumdel + del(k)
!       write (6,'(a,i2,2(a,f10.6))') &
!            ' layer = ', k, ' sl = ', real(sl(k),r4), ' del = ', del(k)
!    enddo
!    write (6,'(/,a,i3,a,f12.8,/)') ' kmax = ', kmax, ' sum del = ', sumdel

!    sl    = sl(kmax:1:-1)
!    si    = si(kmaxp:1:-1)
!    si(1) = 0.000001_r8 ! to avoid floating point exceptions

    
    Write (self%uBin, iostat = iret) ifday, tod, idate, idatec, real(ak,r4), real(bk,r4), int(idvc, i8), int(idsl,i8)
    if(present(istat)) istat = iret

!    deallocate(del)
    deallocate(ak)
    deallocate(bk)

  end subroutine
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: WriteField_MPI - write fields of BAM files from a MPI Pe.
!
! 
! !DESCRIPTION: Esta rotina escreve um campo do modelo BAM 
!               
!
! !INTERFACE:
!   
#ifdef useMPI
  subroutine WriteField_MPI_(OutUnit, field, OutPe, iCount, istat)

     implicit none
!
! !INPUT PARAMETERS:
! 
     ! Output logical unit
     integer(i4),           intent(in   ) :: OutUnit

     ! Field to be writed
     real(r8),              intent(in   ) :: field(:)

     ! How many Pe's are working
     integer(i4),           intent(in   ) :: iCount

     ! What Pe will write
     integer(i4),           intent(in   ) :: OutPe
!
! !OUTPUT PARAMETERS:
! 

     integer(i4), optional, intent(  out) :: istat

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
     character(len=100), parameter :: myname_=':: WriteFields_MPI_( ... )'

     real(r8), allocatable :: buff(:)
     integer :: sizebuff
     integer :: Pe
     integer :: iret
   
     if(present(istat)) istat = 0

     sizebuff = size(field)

     allocate(buff(sizebuff))
      
     do Pe = 0, iCount-1

        if(Pe .eq. OutPe)then
        
           call WriteField_Serial_(OutUnit, field, iret)

           if(iret .ne. 0)then
              write(stdout,*)trim(myname_),': error to write field, ',iret
              if(present(istat)) istat = iret
              return
           endif

        else

           call mpi_recv(buff, sizebuff, MPI_DOUBLE, Pe, MPITag, MPI_COMM_WORLD, status, iret)

           call WriteField_Serial_(OutUnit, buff, iret)
           if(iret .ne. 0)then
              if(present(istat)) istat = iret
              return
           endif
        endif

     end do
      
     deallocate(buff)

  end subroutine
#endif
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: WriteField_Serial_ - 
!
! 
! !DESCRIPTION: 
!              
!
! !INTERFACE:
!   

  subroutine WriteField_Serial_(OutUnit, Field, istat)

     implicit none
!
! !INPUT PARAMETERS:
!
     ! Output logical unit
     integer(i4),           intent(in   ) :: OutUnit

     ! Field to be writed
     real(r8),              intent(in   ) :: Field(:)
!
! !OUTPUT PARAMETERS:
! 

     integer(i4), optional, intent(  out) :: istat

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
     character(len=100), parameter :: myname_=':: WriteFields_Serial_( ... )'

     real(r4), allocatable :: tmp(:)
     integer :: iret
     integer :: siz

     siz=size(Field)
     allocate(tmp(siz))
     tmp = real(Field,r8)

     if(present(istat)) istat = 0

!     write( OutUnit, iostat = iret ) Field
     write( OutUnit, iostat = iret ) tmp

     if(iret .ne. 0)then
        write(*,'(1A)')trim(myname_),': ERROR to write BAM file'
        if(present(istat)) istat = iret
        return
     endif

     return

  end subroutine
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: BAM_SendField_ -
!
! 
! !DESCRIPTION:
!
!
! !INTERFACE:
!   
#ifdef useMPI
  subroutine BAM_SendField_(MyPe, toPe, Field, istat)

     implicit none
!
! !INPUT PARAMETERS:
!
     ! Source Pe
     integer(i4),           intent(in   ) :: MyPe

     ! Target Pe
     integer(i4),           intent(in   ) :: toPe

     ! Field to be send
     real(r8),              intent(in   ) :: field(:)
!
! !OUTPUT PARAMETERS:
! 
     integer(i4), optional, intent(  out) :: istat
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
     character(len=100), parameter :: myname_=':: BAM_SendField_( ... )'

     integer :: iret
     integer :: sizefield

     if(present(istat)) istat = 0

     sizefield = size(field)
      
     call mpi_send(field, sizefield, MPI_DOUBLE, toPe, MPITag, MPI_COMM_WORLD, iret)

     if(iret .ne. 0)then
        write(*,'(2A,I5,1x,A,1x,I5)')trim(myname_),': ERROR to send field from',MyPe,'to',toPe
        if(present(istat)) istat = iret
        return
     endif

  end subroutine
#endif
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: WriteField_MPIr4 - write fields of BAM files from a MPI Pe.
!
! 
! !DESCRIPTION: Esta rotina escreve um campo do modelo BAM 
!               
!
! !INTERFACE:
!   
#ifdef useMPI
  subroutine WriteField_MPIr4_(OutUnit, field, OutPe, iCount, istat)

     implicit none
!
! !INPUT PARAMETERS:
! 
     ! Output logical unit
     integer(i4),           intent(in   ) :: OutUnit

     ! Field to be writed
     real(r4),              intent(in   ) :: field(:)

     ! How many Pe's are working
     integer(i4),           intent(in   ) :: iCount

     ! What Pe will write
     integer(i4),           intent(in   ) :: OutPe
!
! !OUTPUT PARAMETERS:
! 

     integer(i4), optional, intent(  out) :: istat

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
     character(len=100), parameter :: myname_=':: WriteFields_MPIr4_( ... )'

     real(r4), allocatable :: buff(:)
     integer :: sizebuff
     integer :: Pe
     integer :: iret
   
     if(present(istat)) istat = 0

     sizebuff = size(field)

     allocate(buff(sizebuff))
      
     do Pe = 0, iCount-1

        if(Pe .eq. OutPe)then
        
           call WriteField_Serialr4_(OutUnit, field, iret)

           if(iret .ne. 0)then
              write(stdout,*)trim(myname_),': error to write field, ',iret
              if(present(istat)) istat = iret
              return
           endif

        else

           call mpi_recv(buff, sizebuff, MPI_FLOAT, Pe, MPITag, MPI_COMM_WORLD, status, iret)

           call WriteField_Serialr4_(OutUnit, buff, iret)
           if(iret .ne. 0)then
              if(present(istat)) istat = iret
              return
           endif
        endif

     end do
      
     deallocate(buff)

  end subroutine
#endif
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: WriteField_Serialr4 - 
!
! 
! !DESCRIPTION: 
!              
!
! !INTERFACE:
!   

  subroutine WriteField_Serialr4_(OutUnit, Field, istat)

     implicit none
!
! !INPUT PARAMETERS:
!
     ! Output logical unit
     integer(i4),           intent(in   ) :: OutUnit

     ! Field to be writed
     real(r4),              intent(in   ) :: Field(:)
!
! !OUTPUT PARAMETERS:
! 

     integer(i4), optional, intent(  out) :: istat

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
     character(len=100), parameter :: myname_=':: WriteFields_Serialr4_( ... )'

     real(r4), allocatable :: tmp(:)
     integer :: iret
     integer :: siz

     siz=size(Field)
     allocate(tmp(siz))

     if(present(istat)) istat = 0

     write( OutUnit, iostat = iret ) Field
!     write( OutUnit, iostat = iret ) tmp

     if(iret .ne. 0)then
        write(*,'(1A)')trim(myname_),': ERROR to write BAM file'
        if(present(istat)) istat = iret
        return
     endif

     return

  end subroutine
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: BAM_SendField_ -
!
! 
! !DESCRIPTION:
!
!
! !INTERFACE:
!   
#ifdef useMPI
  subroutine BAM_SendFieldr4_(MyPe, toPe, Field, istat)

     implicit none
!
! !INPUT PARAMETERS:
!
     ! Source Pe
     integer(i4),           intent(in   ) :: MyPe

     ! Target Pe
     integer(i4),           intent(in   ) :: toPe

     ! Field to be send
     real(r4),              intent(in   ) :: field(:)
!
! !OUTPUT PARAMETERS:
! 
     integer(i4), optional, intent(  out) :: istat
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
     character(len=100), parameter :: myname_=':: BAM_SendField_( ... )'

     integer :: iret
     integer :: sizefield

     if(present(istat)) istat = 0

     sizefield = size(field)
      
     call mpi_send(field, sizefield, MPI_FLOAT, toPe, MPITag, MPI_COMM_WORLD, iret)

     if(iret .ne. 0)then
        write(*,'(2A,I5,1x,A,1x,I5)')trim(myname_),': ERROR to send field from',MyPe,'to',toPe
        if(present(istat)) istat = iret
        return
     endif

  end subroutine
#endif

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: upper_case - convert lowercase letters to uppercase.
!
! !DESCRIPTION:
!
! !INTERFACE:

  function ucase(str) result(ustr)
    implicit none
  character(len=*), intent(in) :: str
  character(len=len(str))      :: ustr

! !REVISION HISTORY:
! 	13Aug96 - J. Guo	- (to do)
!EOP
!-----------------------------------------------------------------------
    integer i
    integer,parameter :: il2u=ichar('A')-ichar('a')

    ustr=str
    do i=1,len_trim(str)
      if(str(i:i).ge.'a'.and.str(i:i).le.'z')&
         ustr(i:i)=char(ichar(str(i:i))+il2u)
    end do
  end function ucase

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: lower_case - convert uppercase letters to lowercase.
!
! !DESCRIPTION:
!
! !INTERFACE:

  function lcase(str) result(lstr)
    implicit none
    character(len=*), intent(in) :: str
    character(len=len(str))      :: lstr

! !REVISION HISTORY:
! 	13Aug96 - J. Guo	- (to do)
!EOP
!-----------------------------------------------------------------------
    integer i
    integer,parameter :: iu2l=ichar('a')-ichar('A')

    lstr=str
    do i=1,len_trim(str)
      if(str(i:i).ge.'A'.and.str(i:i).le.'Z')&
         lstr(i:i)=char(ichar(str(i:i))+iu2l)
    end do
  end function lcase
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: check - check if the substring is contained in the string.
!
! !DESCRIPTION:
!
! !INTERFACE:

  function check(substr, str) result(found)
     character(len=*) :: substr
     character(len=*) :: str
     logical          :: found
! 	25Mar21 - J. Gerd	- (to do)
!EOP
!-----------------------------------------------------------------------

     integer :: i, istart, iend
     integer :: strLen
     integer :: subStrLen

     character(len=1024) :: string, subString 

     strLen    = len_trim(str)
     subStrLen = len_trim(substr)

     string    = lcase(str)
     subString = lcase(subStr)

     found = .false.
     do i = 1, strLen-subStrLen
        istart = i
        iend   = i+subStrLen-1
        if (trim(substring) == string(istart:iend)) found = .true.
     enddo

     
  end function

  subroutine perr1_(where, message)
     implicit none
     character(len=*), intent(in) :: where
     character(len=*), intent(in) :: message

     write(stderr,'(4A)')'Error at ',trim(where),' : ', trim(message)

  end subroutine

  subroutine perr2_(where, message, cod)
     implicit none
     character(len=*), intent(in) :: where
     character(len=*), intent(in) :: message
     integer,          intent(in) :: cod

     write(stderr,'(4A,1x,I4)')'Error at ',trim(where),' : ', trim(message),cod

  end subroutine
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
   function GetAvailUnit( exclude ) result( lu )

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
      character(len=*), parameter :: myname_ = ':: GetAvailUnit( ... )'

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
            call perr(trim(myname_),'Units from 10 to 254 are already in use!')
            stop
         endif

      end do find_unit

   end function GetAvailUnit


  subroutine setsig (si, sl, del)

    !     calculates the 
    !     sigma value at each interface and the 
    !     sigma value at midpoint of each layer given the
    !     sigma spacing for each layer.
    !
    !     argument(dimensions)          description
    !
    !     del(kmax)            input  : sigma spacing for each layer.
    !     si(kmaxp)            output : sigma value at each interface.
    !     sl(kmax)             output : sigma value at midpoint of
    !                                   each layer : (k=287/1005)
    !
    !                                                              1
    !                                   +-                      + ---
    !                                   |      k+1          k+1 |  k
    !                                   | si(l)    - si(l+1)    |
    !                           sl(l) = | --------------------- |
    !                                   | (k+1) (si(l)-si(l+1)) |
    !                                   +-                     -+
    !
    !     ci(kmaxp)            local  : ci(l)=1.0-si(l).

    implicit none

    real (kind=r8), intent (in   ) :: del(:)
    real (kind=r8), intent (inout) :: sl(:)
    real (kind=r8), intent (inout) :: si(:)

    integer (kind=i4) :: kmax
    integer (kind=i4) :: kmaxp
    integer (kind=i4) :: kmaxm, k
    real    (kind=r8) :: sumdel, rk, rk1, sirk, sirk1, dif

    real (kind=r8), allocatable  :: ci(:)

    kmax  = size(del)
    kmaxp = kmax + 1

    allocate(ci(kmaxp))

    ci(1)  = 0.0
    sumdel = 0.0

    do k = 1, kmax
       sumdel  = sumdel+del(k)
       ci(k+1) = ci(k)+del(k)
    enddo

    ci(kmaxp) = 1.0
    rk        = 287.05/1005.0
    rk1       = rk+1.0

    do k=1,kmaxp
       si( k ) = 1.0 - ci( k )
    enddo

    kmaxm = kmax-1
    do k=1,kmax
       !
       !     dif=si(k)**rk1-si(k+1)**rk1
       !
       sirk=exp(rk1*log(si(k)))
       if (k .le. kmaxm) then
          sirk1=exp(rk1*log(si(k+1)))
       else
          sirk1=0.0
       endif
       dif=sirk-sirk1
       dif=dif/(rk1*(si(k)-si(k+1)))
       !
       !     sl(k)=dif**(one/rk)
       !
       sl(k)=exp(log(dif)/rk)
    enddo


!    write (6,'(/,a,/)') ' from setsig: '
!    do k=1,kmaxp
!       write (6,'(a,i2,2(a,f10.6))') &
!            ' level = ', k, ' ci = ', ci(k), ' si  = ', si(k)
!    enddo
!    write (6,'(a)')' '
!    do k=1,kmax
!       write (6,'(a,i2,2(a,f10.6))') &
!            ' layer = ', k, ' sl = ', sl(k), ' del = ', del(k)
!    enddo
!    write (6,'(/,a,i3,a,f12.8,/)') ' kmax = ', kmax, ' sum del = ', sumdel

    deallocate(ci)

  end subroutine setsig
  !EOC
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: getUV - routine to get U and V components of wind.
!
! 
! !DESCRIPTION: Esta rotina é utilizada para obter os campos de u e v do modelo 
!               BAM.
!               Nesta rotina são lidas as variáveis diverfencia e vorticidade,
!               e, logo após, estes campos são convetidos nas componentes zonal
!               e meridional do vento.
!
! !INTERFACE:

  subroutine getUV2D_(self, wlev, uvel, vvel, istat)
      implicit none
!
! !INPUT PARAMETERS:
! 
      class(BAMFile),    intent(in   ) :: self  ! file information
      integer,           intent(in   ) :: wlev  ! what field level ? 

!
! !OUTPUT PARAMETERS:
! 
      real(r8),          intent(  out) :: uvel(:,:)
      real(r8),          intent(  out) :: vvel(:,:)
      integer, optional, intent(  out) :: istat
!
! !SEE ALSO:
!
!
!
!
! !REVISION HISTORY: 
!
!  20 Aug 2019 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
     real(r8), allocatable :: uwind(:)
     real(r8), allocatable :: vwind(:)
     integer :: ipt, jpt, iret

     ipt = size(uvel,1)
     jpt = size(uvel,2)

     allocate(uwind(ipt*jpt))
     allocate(vwind(ipt*jpt))

     call self%getUV1D_(wlev, uwind, vwind, iret)

     uvel = reshape(uwind,[ipt,jpt])
     vvel = reshape(vwind,[ipt,jpt])

     if(present(istat)) istat = iret

  end subroutine
!EOC
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: getUV - routine to get U and V components of wind.
!
! 
! !DESCRIPTION: Esta rotina é utilizada para obter os campos de u e v do modelo 
!               BAM.
!               Nesta rotina são lidas as variáveis diverfencia e vorticidade,
!               e, logo após, estes campos são convetidos nas componentes zonal
!               e meridional do vento.
!
! !INTERFACE:

  subroutine getUV1D_(self, wlev, uvel, vvel, istat)
      implicit none
!
! !INPUT PARAMETERS:
! 
      class(BAMFile),    intent(in   ) :: self ! file information
      integer,           intent(in   ) :: wlev ! what field level ? 

!
! !OUTPUT PARAMETERS:
! 
      real(r8),          intent(  out) :: uvel(:)
      real(r8),          intent(  out) :: vvel(:)
      integer, optional, intent(  out) :: istat
!
! !SEE ALSO:
!
!
!
!
! !REVISION HISTORY: 
!
!  14 Aug 2019 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
!     type(legendre) :: leg
     integer :: iMax
     integer :: jMax
     integer :: MnWv2
     integer :: MnWv3
     integer :: iret

     real(r8), dimension(:),   allocatable :: divq
     real(r8), dimension(:),   allocatable :: vorq
     real(r8), dimension(:),   allocatable :: uveq
     real(r8), dimension(:),   allocatable :: vveq
     real(r8), dimension(:,:), allocatable :: buft
     real(r8), dimension(:,:), allocatable :: bufw

     !
     ! get BAM file dimensions
     !

     iMax  = self%GetOneDim('iMax')
     jMax  = self%GetOneDim('jMax')
     MnWv2 = self%GetOneDim('MnWv2')
     MnWv3 = self%GetOneDim('MnWv3')

     !
     ! Get Divergence and Vorticity
     !

     allocate(divq(MnWv2))
     allocate(vorq(MnWv2))
     call self%getField('DIVERGENCE', wlev, divq, iret)
     if(present(istat)) istat = iret
     call self%getField('VORTICITY' , wlev, vorq, iret)
     if(present(istat)) istat = istat + iret


     !
     ! Transform to uvel and vvel
     !

     allocate(uveq(MnWv3))
     allocate(vveq(MnWv3))

     call self%DivgVortToUV(divq, vorq, uveq, vveq)
     deallocate(divq)
     deallocate(vorq)


     if(size(uvel) .eq. (iMax*jMax))then

!        call createFFT (iMax)
        
        !
        ! UVEL
        !

        allocate (bufw (imax,jmax))
        call self%Spec2Grid(uveq,bufw)
        deallocate(uveq)

        uvel = reshape(bufw,[size(uvel)])
        deallocate(bufw)

        !
        ! VVEL
        !

        allocate (bufw (imax,jmax))
        call self%Spec2Grid(vveq,bufw)
        deallocate(vveq)

        vvel = reshape(bufw,[size(vvel)])
        deallocate(bufw)

     else

        uvel = uveq
        vvel = vveq

        deallocate(uveq)
        deallocate(vveq)

     endif

  end subroutine
!
!EOC
!
!-----------------------------------------------------------------------------!


end module sigiobammod
