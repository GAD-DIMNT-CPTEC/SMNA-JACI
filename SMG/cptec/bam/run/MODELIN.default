!namelist
!############################### Change Log ##################################
! 1.0.0.0
!
!  $Author: pkubota $
!  $Date: 2009/03/03 16:36:38 $
!  $Revision: 1.13 $
!
!
!#############################################################################
!
! trunc    =021,vert     =09,dt       =1800.0,
! trunc    =062,vert     =28,dt       =1200.0,
! trunc    =126,vert     =28,dt       =600.0
! trunc    =170,vert     =42,dt       =450.0
! trunc    =213,vert     =42,dt       =360.0
! trunc    =254,vert     =64,dt       =300.0
! trunc    =341,vert     =64,dt       =200.0  !!!teste=225.0 
!############################################################################
&MODEL_RES
 trunc    = #TRUNC#,		         ! TRC   : three-digit triangular truncation
 vert     = #NLEV#,		         ! LV    : two-digit number of vertical sigma-layers
 dt       = #DELT#,	            ! dt    : delta t
 IDATEI   = #LABELI#,            ! LABELI: initial forecasting label
 IDATEW   = #LABELW#,            ! LABELC: final forecasting label for cold
 IDATEF   = #LABELF#,            ! LABELF: final forecasting label for warm
                                 !         if cold, then LABELF=LABELC
                                 !
 NMSST    ='#NMSST#',            ! NMSST : sst file name
                                 !    sstaoi : for optimum interpolated SST
                                 !             climatology (original 1x1 degree)
                                 !    sstanp : for optimum interpolated SST
                                 !             climatology plus persisted SST anomaly
                                 !    sstwkl : for weekly run mean of SST for
                                 !             the week finished at initial day
                                 !             minus one (original 1x1 degree)
                                 !    sstwkd : for weekly run mean of SST
                                 !             direct access file (original 1x1 degree)
                                 !    sstmtd : for monthly run mean of SST
                                 !             direct access file (original 1x1 degree)
                                 !    sstdyd : for daily run mean of SST
                                 !             direct access file (original 1x1 degree)

 NMSNOW   ='snowaoi',            ! NMSNOW : snow file name
 NMSOILM  ='#NMSOILM#',          ! NMSOILM='soilmaoi' NMSOILM='soilmwkl'
                                 !          snowwkl :-- > daily (NCEP file)
                                 !                        for weekly run mean of snow for
                                 !                        the week finished at initial day
                                 !                        minus one (original 1x1 degree)
                                 !          snowaoi :-- > for optimum interpolated SST
                                 !                       climatology 
 DHFCT     = #DHFCT#,            ! DHFCT : > 0 interval in hours to output diagnostics
                                 !            = 0 to use default list of interval in hours 
                                 !            < 0 interval in months to output diagnostics 
 DHRES     = #DHRES#,            ! DHRES : interval in hours to output restart,         
                                 !        -365 (anual), -24 (mensal) ,6 (6 em 6) ,24 (1dia)
                                 !        equal zero to use default list        
 DHDHN     = 0,                  ! DHDHN : interval in hours to output DHN diagnostic
 NHDHN     = 0,                  ! NHDHN : time in hours to stop DHN diagnostics,         
                                 !            equal zero to not execute DHN diagnostics         
 DHEXT     = 0,                  ! DHEXT : interval in hours to output extra diagnostics
 NHEXT     = 00000,              ! NHEXT : time in hours to stop extra diagnostics,         
                                 !            equal zero to not execute extra diagnostics  
 DOGRH     = .FALSE.,            ! DOGRH : logical (T or F) to do grid history for         
                                 !            selected points
 DOPRC     = .FALSE.,            ! DOPRC : logical (T or F) to do time step output
                                 !            of global precipitation          
 StrFormat = 'new',              ! Get NCEP Analysis and dumping to 'old' format or 'new' format
 FIELDGAUS = .FALSE.,
 FILELIST  = .FALSE.,
 PREFX     = '#PREFX#',          ! PREFX : preffix for name of output files
 PREFY     = '#PREFY#', 	      ! PREFY : preffix for name of input files
 TABLE     ='p',                 ! TABLE : indicator of the desire table:          
                                 !                     n => for desirtable (default)
                                 !                     p => for desirtable.pnt
                                 !                     c => for desirtable.clm
                                 !                      LABELx: hhddmmyyyy
                                 !                              yyyy = four digit year
                                 !                                mm = two digit month
                                 !                                dd = two digit day
                                 !                                hh = two digit hour
 path_in        = '#PATHIN#', 
 dirfNameOutput = '#PATHOU#', 
 dirRstInput    = '#RSTIN#', 
 dirRstOutput   = '#RSTOU#', 
/
&MODEL_IN
 slagr              =#SLAGR#,       ! Semi-Lagrangian option (.FALSE. for Eulerian Model) (TRUE for TQ666)
 slhum              =.TRUE.,       ! Humidit should be lagrangian  ? (Eulerian Model)
grid_difus          =.TRUE.,       !
 microphys          =.TRUE.,       ! Microphysics requires slhum to be true
 nlnminit           =.TRUE.,       ! do normal mode non linear initialization
 diabatic           =.TRUE.,       ! diabatic or not initialization (.FALSE. if nlnminit=.FALSE.)
 eigeninit          =#EIGENINIT#,  ! eigenInit  --> .FALSE.
 mgiven             =#MGIVEN#,     ! mgiven  --> .FALSE.
 gaussgiven         =#GAUSSGIVEN#, ! gaussgiven  --> .FALSE.
 reducedGrid        =.TRUE.,       ! reduced    --> .TRUE. 
 linearGrid         =.FALSE.,      ! lineargrid --> .FALSE.
 GenAssFiles        =.FALSE.,      ! .TRUE. flag to generate  assimilation files
 GenRestFiles       =#GENRST#,     ! .TRUE. to generate restart files.
 rmRestFiles        =#RMRST#,      ! .TRUE. to remove restart files after read them.
 MasCon             =.TRUE.,     ! .TRUE. to do Mass Conservation: ln(ps)
 MasCon_ps          =.FALSE.,    ! .TRUE. to do Mass Conservation: ps
 FluxCon_co2        =.FALSE.
 nscalars           = 0,         ! Number of passive variables
 tamBlock           = 64,        ! quantidade de fft enviadas por bloco   
 ibdim_size         = 6,         ! tamanho de bloco basico (ibmax)     
 nproc_vert         = 28
/
&PHYSPROC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ixxxx=yes  the physical process included
!     ixxxx=no   the physical process excluded
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 UNIFIED = #UNIFIED#,  ! Enable Unified Physics
 !rootmode='MILENA'    ! rootmode='JACKSON' ! Jackson et al., 1996: or rootmode='MILENA' Milena et al., 2020:

 ISWRAD  = 'CRD'       ! iswrad = NON: don't do sw radiation, LCH: lacis & hansen
                       !          CRD: clirad                UKM: ukmet office
                       !          RRTMG: cptec               CRDTF: clirad tarasova 
 ILWRAD  = 'CRD'       ! ilwrad = NON: don't do lwradiation, HRS: Harshvardhan
                       !          CRD: clirad                UKM: ukmet office
                       !          RRTMG: cptec               CRDTF: clirad tarasova
                       !
 ICCON   = 'ARA',      ! iccon=KUO:cumulus convection(kuo)
                       ! iccon=ARA:cumulus convection(ARAKAWA)
                       ! iccon=RAS:cumulus convection(RELAXED ARAKAWA
                       ! iccon=GRE:cumulus convection(grell-nilo)
                       ! iccon=GEC:cumulus convection(grell-cptec)
                       ! iccon=ZMC:cumulus convection(Zhang-McFarlane)
                       ! iccon=NON:without cumulus convection
                       !
 ISCON   = 'TIED',     ! iscon=TIED:shallow convection this process follows cumulus convection (tiedke)
                       ! iscon=MFLX:shallow convection this process follows cumulus convection (tiedke)
                       !
 ILCON   = 'HUMO',     ! ilcon=yes:large scale condensation
                       !  'NON' 
                       !  'YES' 
                       !  'LSC' 
                       !  'MIC'
                       !  'HWRF'
                       !  'HGFS'
                       !  'UKMO'
                       !  'MORR'
                       !  'HUMO'
                       !
 IGWD    = 'YES',      ! igwd =NON:gravity wave 
                       ! igwd =YES:gravity wave 
                       ! igwd =CAM:gravity wave 
                       ! igwd =USS:gravity wave 
                       ! igwd =GMB:gravity wave
                       ! isimp ="YES "         ! isimp=yes:simplified physics version. 
 CRDCLD = 6,           ! cloud scheme =1 (old) =4 (ccm3) =5 (cam5)  =6 (gfs)
 ISIMP   = 'NO ',      ! isimp=yes:simplified physics version.
 thermcell=1,          !
 atmpbl = 4,           !
                       ! 
 schemes = 3 ,         ! 1 ssib; 2 sib2; 3 ibis
 isimco2 = 0,          !
 gbaco2  = 0.0e-10,    !-4.1109e-09  !-3.01522e-09, !      
 OCFLUX  = 'UKME',     ! OCFLUX  = 'COLA',
                       ! OCFLUX  = 'UKME',
                       ! OCFLUX  = 'WGFS',  
 SLABOCEAN  = 'SLAB', ! SLABOCEAN = 'COLA'
                       ! SLABOCEAN = 'SLAB'
 ICEMODEL   = 'SSIB',  ! ICEMODEL = 'COLA'
                       ! ICEMODEL = 'SSIB'
 nAeros     = 0
 /
&PHYSCS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   the following are vertical resolution dependent cloud parameters
!   used in cldgen.  correct settings for these parameters need to be
!   determined experimentally.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 swint   = #SWINT#, ! short wave radiation call interval in hours swint=1.0e0
 trint   = #TRINT#, ! ir subr. call interval in hours  trint=3.0e0
                    ! long  wave radiation call interval in hours
                    ! physical constants for simple physics options
 co2val  =370.0,    ! co2val is wgne standard value in ppm 
                    ! co2val=345.0 old value 
                    ! co2val=370.0 new value
 /

&COMCON
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     files
!     ifxxx=0    xxx is not processed
!     ifxxx=1    xxx is set to month=idatec(2) in the first call,
!                but not processed from the subsequent calls.
!                ifxxx is set to zero after interpolation
!     ifxxx=2    xxx is interpolated to current day and time every fint
!                hours synchronized to 00z regardless of initial time.
!                interpolation is continuous (every time step) if fint<0.
!     ifxxx=3    xxx is interpolated to current day and time when ifday=0
!                and tod=0.0 but not processed otherwise
!                ( appropriate only when xxx is predicted )
!
!                the following are for sst only (fint applies as in
!                ifxxx=2):
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 initlz  =#INITLZ# , ! initlz =2 diabatic initialization and normal mode initialization
                   !         =1 diabatic with no normal mode initialization
                   !         =0 adiabatic with no normal mode initialization
                   !         <0 same as >0 with sib variables read in instead of 
                   !            initialized
 fint    =6,       ! surface boundary calling interval in hours
 ifsst   =-1,      ! ifsst=4  sst is linearly interpolated from continuous 
                   !              direct access data set to current day and 
                   !              time.data set is assumed to be spaced every 
                   !              intsst days or every calendar month is 
                   !              intsst < 0.
                   ! ifsst=5  sst is expanded from piecewise cubic 
                   !              coefficients in        direct access data set to 
                   !              current day and time. data set
                   !                is assumed to be spaced every intsst days.
 filta   =0.92,
 alfa    =1.0  !0.533,
 ifsnw   =3,
 ifalb   =0,
 ifslm   =3,  
 ifco2   =0,       ! 0  co2val Using a single global value for co2 concentration
                   !-1  "with variable co2val following 2nd order fit to Mauna Loa data (1958-04)")
 ifozone =1,       !ifozone =0 "OZONE field is constant";ifozone =1 analisys advec
 iftracer=0,       !"TRACER field"
                   !tar begin
                   !  ifaeros=0  ....original aerosol: fixed continental and oceanic aerosol types
                   !  ifaeros=2 .... Climate aerosol of Kinne, 2013
                   !
 ifaeros=0         !
                   !tar end
 dk=8.0E+15,
 tk=6.0E+15,
Gamma_Factor=2.15e-3 !GammaFactor =1.15e-3  --> TQ0062  non-dimensional empirical factors
                     !GammaFactor =3.25e-2  --> TQ0126  non-dimensional empirical factors
                     !GammaFactor =3.25e-1  --> TQ0299  non-dimensional empirical factors
                     !GammaFactor =6.25e-0  --> TQ0299  non-dimensional empirical factors
 sflyr=3,            !

/
&NumberOutPutForecast
  maxtfm2= 17
  cthrfx( 1: 4) =   6.0, 12.0,  18.0,  24.0,
  cthrfx( 5: 8) =  30.0, 36.0,  42.0,  48.0,
  cthrfx( 9:12) =  54.0, 60.0,  66.0,  72.0,
  cthrfx(13:16) =  84.0, 96.0, 120.0, 144.0,
  cthrfx(17:17) = 168.0,
/
