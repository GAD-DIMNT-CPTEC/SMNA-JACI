!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOI
! !TITLE: 
!      spc2grd Conversor de arquivo fct spectral para ponto de grade
!
! !AUTHORS: Luiz F. Sapucci based on fct2anl_trans.f90 from João Gerd Zell de Mattos
!
! !AFFILIATION: Modeling and Development Division, CPTEC/INPE
!
! !DESCRIPTION:
!      Este programa faz a conversao de um arquivo de previsao do modelo global
!      BAM (arquivo fct) no espaco espectral para ponto de grade onde o GSI trabalha a assimilacao.
!      Uma potencial aplicacao e na avaliacao dos campos do modelos usados como background 
!      no processo de assimilacao. 
!      Sao necessarios dois arquivos:
!         1 - arquivo scestral (previsao ou analise)
!         2 - arquivo dir      (previsao ou analise)
!
! ! INOUT FILE:
!  input argument list:
!     spc2grd.config                     - argument list with number and list of variable required   
!     BAM.fct.hh or BAM.anl = [nameFile] - spectral format file of forecast or analysis from BAM model (hh=06 hour for background)
!     BAM.dir.hh or BAM.dir.anl.gsi      - directives of the forecast file from BAM model (hh=06 hour for background)
!     BAM.bin.ctl.template               - template ctl file for binary file created
!
!   output
!     [nameFile].bin      - binary file of the BAM analyis or forecast of hh with field required  
!     [nameFile].bin.ctl  - ctl file for read the binry file created.
!
! !REMARKS:
!     (1) Resolucoes aptas nessa rotina TQ0062L28 TQ0299L64 incluir template para outras resolucoes
!
! !HISTORY REVISION:
!     03 Abr 2017 - Luiz Sapucci    - Initial Code only forecast file
!     07 Abr 2017 - Luiz Sapucci    - Version with analysis and forecast file TQ0062L28
!     18 Jan 2019 - Luiz Sapucci    - Adaptacoes no config para rodar a versão TQ0299L64 funcionando no Eval
!     12 Jul 2021 - J.G.Z de Mattos - Add new version of sigio_BAMMod 
!
! !USES:
!     * sigio_BAMMod - Modulo contendo as rotinas necessarias para a leitura e 
!                      escrita dos arquivos do BAM
!EOI
!-----------------------------------------------------------------------------!
!


program spc2grd
   use sigio_BAMMod, only : BAMFile
   use ModConstants, only : r8,      & ! Kind for 64-bits Real Numbers
                            r4         ! Kind for 64-bits Real Numbers


   implicit none

   integer :: istat

   integer :: IMax, i, j, N
   integer :: JMax
   integer :: KMax
   integer :: Mend
   integer :: Nvar

   real    :: mwave2
   integer :: mnwv2, mnwv3
   
   integer :: mydate(5) 

   integer :: icount
   integer :: ivar
   integer :: ilev
   integer :: Nlevs
   integer :: NVars

   type(BAMFile) :: BAM

   real(r8), allocatable :: field(:),umimax(:),umimin(:)
   real(r4), allocatable :: grid(:,:)

   integer, parameter :: stderr = 0
   integer, parameter :: stdinp = 5
   integer, parameter :: stdout = 6

   character(len=20) :: myname_='spc2grd'
   character(len=20) :: spcFILEname,dirFILEname
   character(len=31) :: ctlTemplate
   character(len=20) :: binName,ctlName
   character(len=82) :: mensa = '#'
   character(len=40), dimension(30) :: VName,Vlable

   write (*,*)
   write (*,*) "-------------------------------------"
   write (*,*) " GSI util rotine: spc2grd:  "
   write (*,*) "-------------------------------------"
   write (*,*)

   ! Open files config and template
      
   OPEN (UNIT=9,FILE='spc2grd.config',STATUS='old')

   !
   ! lendo o arquivo de configuracao e carregando os parametros
   !

   do while (mensa(1:1) == '#' .OR. mensa(1:1) == '$' .OR. mensa(1:1) == ' ')
     read(9,'(A40)')mensa
!     write(*,*)'dentro loop',mensa
   enddo 
 
   backspace(9)
   
   read(9,'(A13,A15)')mensa,spcFILEname ! Nome do arquivo de entrada espectral previsao ou analise
   read(9,'(A13,A15)')mensa,dirFILEname ! Nome do arquivo de directivas do arquivo espectral de previsao ou analise
   read(9,'(A13,A31)')mensa,ctlTemplate ! Nome do arquivo ctl template para a resolucao em uso do tipo Q0299L064
   
   read(9,'(A6,I4)')mensa,NVars

   do Nvar = 1, NVars
     read(9,'(A4,2x,A26)')Vlable(Nvar),VName(Nvar)
   enddo

   ! Make output file name and open them  

   N = 1
   do while (spcFILEname(N:N) /= ' ')
      N=N+1
   enddo 

   binName=spcFILEname(1:N-1)//'.bin'
   ctlName=spcFILEname(1:N-1)//'.bin.ctl'

   OPEN (UNIT=13,FILE=binName,status='unknown',form='unformatted')
   OPEN (UNIT=14,FILE=ctlTemplate,status='old')
   OPEN (UNIT=15,FILE=ctlName,status='unknown')

! Generating the ctl file for this bin file

   write(15,'(A6,A20)')'DSET ^',binName
      
   do while (mensa(1:4) /= 'TDEF')
     read(14,'(A82)')mensa
     write(15,'(A82)')mensa
   ENDDO
 
   write(15,'(A4,1X,I2)')'VARS',NVars

   !
   ! Open and read header
   !

   call BAM%Open(dirFILEname, spcFILEname, istat=istat)
   if (istat.ne.0)then
      write(*,*)'Problem to read BAM files!'
   endif

   ! Get BAM dimensions
   call BAM%getDims(imax, jmax, kmax, Mend, istat)

   mnwv2 = (Mend+1)*(Mend+2)
   mnwv3 = (mnWV2+2)*(Mend+1)
   write(stdout,*)'IMAx',IMax
   write(stdout,*)'JMax',JMax
   write(stdout,*)'KMax',KMax
   write(stdout,*)'Mend',Mend
   write(stdout,*)'mnWV2',mnwv2
   write(stdout,*)

   ! Get atmospheric fields necessary to IC file
   ! - Nota: Os campos são lidos e escritos no formato
   !         espectral. Em caso de duvida veja a rotina 
   !         BAM_GetField no modulo sigio_BAMMod.F90

   icount = 1
   ivar   = 1
   do while (ivar .le. nvars)
      NLevs = BAM%getNlevels(trim(VName(ivar)),istat)
   
      allocate(umimax(NLevs))
      allocate(umimin(NLevs))
      
      do ilev = 1, NLevs

  
            allocate(grid(imax,jmax))
                     
!DEBUG      write(stdout,'(A10,A20,1x,A9,1x,I4)')'Reading: ',trim(VName(ivar)),'at ilevel', ilev
            call BAM%getField(trim(VName(ivar)), ilev, grid, istat)

            umimax(ilev)=MAXVAL(grid)
            umimin(ilev)=MINVAL(grid)      
            write(13) grid

            deallocate(grid)

         icount = icount + 1
   
      enddo

      write(*,*)
      write(*,*)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
      write(stdout,'(A20,A20,A25)')'%%%%%%%%%%%%%%     ',trim(VName(ivar)),'          %%%%%%%%%%%%%%%'
      write(stdout,*)"         level       Maximum values             Minimum  values "

      do ilev = 1, NLevs 
        write(stdout,*) ilev,umimax(ilev),umimin(ilev)      
      enddo

      write(stdout,*)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
      write(stdout,*)

      write(15,'(A4,I3,I3,2X,A26)')Vlable(ivar),Nlevs,99,VName(ivar) ! Variable information in the ctl file 

      deallocate(umimax)
      deallocate(umimin)


      ivar = ivar + 1

   enddo
 
   write(15,'(A7)')'ENDVARS'  ! Finalizing the ctl file 
   call BAM%close(istat)

   close(13)
   close(14)   
   close(15)
   
   write(stdout,'(A)')'Program ends normaly'
   
end program
