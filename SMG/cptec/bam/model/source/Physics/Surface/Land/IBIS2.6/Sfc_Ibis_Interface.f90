MODULE Sfc_Ibis_Interface
  !
  ! InitFieldsIbis__vegin
  !              |                
  !              |__RD_PARAM
  !              |
  !              |__readit__ibismap__cellbox
  !              |                |
  !              |                |__cellbox 
  !              |
  !              |__climanl__existence
  !              |
  !              |__initial__coldstart
  !              |        |
  !              |        |__restart__existence
  !              |        |
  !              |        |__inisurf
  !              |        |
  !              |        |__inisnow
  !              |        |
  !              |        |__inisoil
  !              |        |    
  !              |        |__iniveg
  !              |        | 
  !              |        |__inisum
  !
  ! IbisDrv__Ibis__pheno
  !       |     |
  !       |     |__soilbgc
  !       |     |
  !       |     |__lsxmainn_______setsoi
  !       |     |              |
  !       |     |              |__fwetcal
  !       |     |              |
  !       |     |              |__solset
  !       |     |              |
  !       |     |              |__solsur
  !       |     |              |
  !       |     |              |__solalb___twostr__twoset
  !       |     |              |        |
  !       |     |              |        |__twostr__twoset
  !       |     |              |
  !       |     |              |__solarf
  !       |     |              |
  !       |     |              |__irrad
  !       |     |              |
  !       |     |              |__cascade___mix
  !       |     |              |         | 
  !       |     |              |         |__mix
  !       |     |              |         |
  !       |     |              |         |__steph2o__mix
  !       |     |              |         |         |
  !       |     |              |         |         |__mix
  !       |     |              |         |
  !       |     |              |         |__steph2o__mix
  !       |     |              |         |         |
  !       |     |              |         |         |__mix
  !       |     |              |         |
  !       |     |              |         |__mix
  !       |     |              |         |
  !       |     |              |         |__mix
  !       |     |              |         |
  !       |     |              |         |__steph2o__mix
  !       |     |              |         |         |
  !       |     |              |         |         |__mix  
  !       |     |              |         |
  !       |     |              |         |
  !       |     |              |         |__mix
  !       |     |              |         |
  !       |     |              |         |__mix
  !       |     |              | 
  !       |     |              |__fwetcal
  !       |     |              |
  !       |     |              |__canopy__canini
  !       |     |              |       |
  !       |     |              |       |__drystress
  !       |     |              |       |
  !       |     |              |       |__turcof__fstrat
  !       |     |              |       |       |
  !       |     |              |       |       |__fstrat
  !       |     |              |       |
  !       |     |              |       |__stomata
  !       |     |              |       |
  !       |     |              |       |__turvap___impexp
  !       |     |              |       |       |
  !       |     |              |       |       |__impexp
  !       |     |              |       |       |
  !       |     |              |       |       |__impexp
  !       |     |              |       |       |
  !       |     |              |       |       |__impexp
  !       |     |              |       |       |
  !       |     |              |       |       |__impexp2
  !       |     |              |       |       |
  !       |     |              |       |       |__linsolve
  !       |     |              |       |
  !       |     |              |       |__tscreen
  !       |     |              |       |
  !       |     |              |       |__tscreen
  !       |     |              |
  !       |     |              |__cascad2__steph2o2
  !       |     |              |        |
  !       |     |              |        |__steph2o2
  !       |     |              |        |
  !       |     |              |        |__steph2o2
  !       |     |              |
  !       |     |              |__noveg
  !       |     |              |
  !       |     |              |__snow__snowheat___tridia
  !       |     |              |     |
  !       |     |              |     |__vadapt
  !       |     |              |     |
  !       |     |              |     |__MixSnow
  !       |     |              |
  !       |     |              |__soilctl__soilh2o__tridia
  !       |     |                       |
  !       |     |                       |__soilheat__tridia
  !       |     |                       |
  !       |     |                       |__wadjust
  !       |     |                       |
  !       |     |                       |__wadjust
  !       |     |
  !       |     |
  !       |     |__sumnow
  !       |     |
  !       |     |__sumday
  !       |     |
  !       |     |__summonth
  !       |     |
  !       |     |__sumyear
  !       |     |
  !       |     |__solset
  !       |     |
  !       |     |__solsur
  !       |     |
  !       |     |__solalb___twostr__twoset
  !       |     |       |
  !       |     |       |___twostr__twoset
  !       |     |
  !       |     |__solarf
  !       |
  !       |__co2
  !       |
  !       |__dynaveg2__fire
  !       |         |
  !       |         |__vegmap
  !       |
  !       |__climanl2__existence
  !
  !
  ! Albedo_IBIS__fwetcal
  !           |
  !           |__solset
  !           |
  !           |__solsur
  !           |
  !           |__solalb___twostr__twoset
  !                   |
  !                   |___twostr__twoset

  USE Constants, ONLY :     &
       oceald   ,     &
       gasr     ,     &
       icealn   ,     &
       icealv   ,     &
       r8,i8,r4,i4,MinFlux,MaxFlux

  USE Sfc_Ibis_Fiels , ONLY : spinmax ,isimco2 ,doalb,  &
       isimveg ,isimco2,isimfire,nband, nsoilay, nsnolay,nlpoints,&
       pi,stef,vonk,grav,tmelt,hfus,hvap,hsub,ch2o,npft, &
       cice,cair,cvap,rair,rvap,cappa,rhow,epsilon,&
       ginvap  ,gsuvap  ,gtrans,gtransu,gtransl,grunof,& 
       gdrain  ,gadjust ,a10td,a10ancub,a10ancuc,a10ancls, & 
       a10ancl3,a10ancl4,a10scalparamu,a10daylightu,a10scalparaml,&
       a10daylightl,vmax_pft,tau15,kc15,ko15,cimax,gammaub,alpha3,&
       theta3,beta3,coefmub,coefbub,gsubmin,gammauc,coefmuc,coefbuc ,&
       gsucmin,gammals,coefmls,coefbls,gslsmin,gammal3,coefml3, &
       coefbl3,gsl3min,gammal4,alpha4,theta4,beta4,coefml4,coefbl4, &
       gsl4min, wliqum,wliqu,wliqu0,wliqumax,wsnoum,wsnou,wsnou0,wsnoumax,tum,tu,tu0,wliqsm,wliqs,wliqs0,wliqsmax,& 
       wsnosm,wsnos,wsnos0,wsnosmax,tsm,ts,ts0,wliqlm,wliql,wliql0,wliqlmax,wsnolm,wsnol,wsnol0,wsnolmax,tlm,tl,tl0,topparu,&
       topparl,fl,fu,lai,sai,rhoveg,tauveg,orieh,oriev,wliqmin, &
       wsnomin,t12,tdripu,tblowu,tdrips,tblows,t34,tdripl,tblowl,&
       ztop,alaiml,zbot,alaimu,froot,q34,q12,su,cleaf,dleaf,ss        ,& 
       cstem,dstem,sl,cgrass,ciub,ciuc,exist,csub,gsub,csuc,gsuc,&
       agcub,agcuc,ancub,ancuc,totcondub,totconduc,cils,cil3,&         
       cil4,csls,gsls,csl3,gsl3,csl4,gsl4,agcls,agcl4,agcl3,ancls,&         
       ancl4,ancl3,totcondls,totcondl3,totcondl4,chu,chs,chl,frac,& 
       tlsub,z0sno,rhos,consno,hsnotop,hsnomin,fimin,fimax,fi,&     
       tsnom,tsno,tsno0,hsno,sand,clay,poros,wsoim,wsoi,wsoi0,wisoim,wisoi,wisoi0,consoi,zwpmax, wpud,&        
       wipud,wpudmax, qglif ,tsoim ,tsoi,tsoi0,hvasug,hvasui,albsav,albsan,&
       tg,tim,ti,ti0,z0soi,swilt,sfield,stressl,stressu,stresstl,stresstu,&
       csoi,rhosoi,hsoi,suction,bex,upsoiu,upsoil,heatg,heati,&
       hydraul,porosflo,ibex,bperm,hflo,o2conc,co2conc,cbiow,&
       sapfrac,cbior ,adrain,adsnow,adaet,adtrunoff,& 
       adsrunoff,addrainage,adrh,adsnod,adsnof,adwsoi,adtsoi, &
       adwisoi,adtlaysoi,adwlaysoi,adwsoic,adtsoic,adco2mic,adco2root, & 
       adco2soi,adco2ratio,adnmintot,decompl,decomps,tnmin,amrain,&     
       amsnow,amaet,amtrunoff,amsrunoff,amdrainage,amtemp,&    
       amqa,amsolar,amirup,amirdown,amsens,amlatent,amlaiu,& 
       amlail,amtsoi,amwsoi,amwisoi,amvwc,amawc,amsnod,amsnof,amnpp,adnpp,&  
       amnpptot,amco2mic,amco2root,amco2soi,amco2ratio,amneetot, &
       amnmintot,amts2,amtransu,amtransl,amsuvap,aminvap,amalbedo,&
       amtsoil,amwsoil,amwisoil,tnpptot,&
       aysolar,ayirup,ayirdown,aysens,aylatent,ayprcp,&    
       ayaet,aytrans,aytrunoff,aysrunoff,aydrainage,aydwtot,aywsoi, &
       aywisoi,aytsoi,ayvwc,ayawc,aystresstu,aystresstl,aygpp,aygpptot,&  
       aynpp,aynpptot,ayco2mic,ayco2root,ayco2soi,ayneetot,ayrootbio, & 
       aynmintot,ayalit,ayblit,aycsoi,aycmic,ayanlit,aybnlit,aynsoi,ayalbedo, &    
       totalit,totrlit,totcsoi,totcmic,totanlit,totrnlit,totnsoi,&   
       totnmic,totlit,totfall,totnlit,firefac,wtot,storedn,yrleach,& 
       ynleach,adfalll,adfallr,adfallw,falll,fallr,fallw,clitlm,clitls,clitrm,clitrs,clitwm,&  
       clitws,csoislop,csoislon,csoipas,clitll,clitrl,clitwl,tw,&
       tc,agddu,tempu,agddl,templ,dropu,dropls,dropl4,dropl3,plai,adplai,&
       deltat,gdd0,gdd0this,tcthis,twthis,tcmin,gdd5,gdd5this,TminL,&        
       TminU,Twarm,GDD,aleaf,awood,cbiol,aroot,specla,td ,vzero,&  
       biomass,totlaiu,totlail,totbiou,totbiol,woodnorm,vegtype0,&
       tauwood0,tauwood,tauleaf,tauroot,xminlai,cdisturb,ayanpp,&
       ayanpptot,asurd,asuri,ndaypm,idateprev,iyear0,&
       ndtimes,nmtimes,nytimes,nppdummy,tco2root,&
       tneetot,tco2mic,iMaskIBIS,nVegClass,&
       ynleach_p ,tnmin_p   ,totnmic_p ,totnlit_p ,totanlit_p, &
       totrnlit_p,totnsoi_p ,storedn_p  ,adcbiol,adcbior,adcbiow,beta1,beta2,stressfac,avmuir_factor

  USE Sfc_Ibis_LsxMain  , ONLY : co2,lsxmain, fwetcal,solset,solsur ,solalb,solarf

  USE Sfc_Ibis_Vegetation, ONLY : dynaveg1,dynaveg2,climanl2,&
       pheno ,sumnow,sumday,summonth,sumyear,soilbgc,DailyDynaVeg

  USE FieldsPhysics, ONLY: Tsfc0,Qsfc0,Tsfcm,Qsfcm,w0,wm,capac0,capacm,td0,tdm,tcm,tc0,tgm,tg0,qm0,tm0,sheleg,gco2flx,cflxm

  USE Sfc_SeaFlux_Interface   , Only :  seasfc

  USE Utils  , ONLY:       &
       totflux

  USE Parallelism, ONLY: &
       MsgOne,           &
       myId,             &
       maxNodes
  USE Diagnostics, ONLY: &
       updia, &
       dodia, &
       nDiag_biomau, &
       nDiag_biomal, &
       nDiag_tlaiup, & 
       nDiag_tlailw, & 
       nDiag_tstnsp, &
       nDiag_wsttot, &
       nDiag_lidecf, &
       nDiag_somdfa, &
       nDiag_facuca, &
       nDiag_fsfclc, &
       nDiag_frsnow, &
       nDiag_insnpp, & ! instantaneous npp (mol-CO2 / m-2 / second)
       nDiag_insnee, & ! instantaneous net ecosystem exchange of co2 per timestep (kg_C m-2/timestep)
       nDiag_grbdy0, & ! annual total growing degree days for current year > 0C
       nDiag_grbdy5, & ! annual total growing degree days for current year > 5C  
       nDiag_avet2m, & ! monthly average 2-m surface-air temperature 
       nDiag_monnpp, & ! monthly total npp for ecosystem (kg-C/m**2/month)
       nDiag_monnee, & ! monthly total net ecosystem exchange of CO2 (kg-C/m**2/month)
       nDiag_yeanpp, & ! annual total npp for ecosystem (kg-c/m**2/yr)
       nDiag_yeanee, & ! annual total NEE for ecosystem (kg-C/m**2/yr) 
       nDiag_upclai, & ! upper canopy single-sided leaf area index (area leaf/area veg)
       nDiag_lwclai, & ! lower canopy single-sided leaf area index (area leaf/area veg)
       nDiag_pfts01, & ! pft tropical broadleaf evergreen trees
       nDiag_pfts02, & ! pft tropical broadleaf drought-deciduous trees
       nDiag_pfts03, & ! pft warm-temperate broadleaf evergreen trees
       nDiag_pfts04, & ! pft temperate conifer evergreen trees
       nDiag_pfts05, & ! pft temperate broadleaf cold-deciduous trees
       nDiag_pfts06, & ! pft boreal conifer evergreen trees
       nDiag_pfts07, & ! pft boreal broadleaf cold-deciduous trees
       nDiag_pfts08, & ! pft boreal conifer cold-deciduous trees
       nDiag_pfts09, & ! pft evergreen shrubs
       nDiag_pfts10, & ! pft cold-deciduous shrubs
       nDiag_pfts11, & ! pft warm (c4) grasses
       nDiag_pfts12, & ! pft cool (c3) grasses
       nDiag_biol01, & ! cbiol tropical broadleaf evergreen trees
       nDiag_biol02, & ! cbiol tropical broadleaf drought-deciduous trees
       nDiag_biol03, & ! cbiol warm-temperate broadleaf evergreen trees
       nDiag_biol04, & ! cbiol temperate conifer evergreen trees
       nDiag_biol05, & ! cbiol temperate broadleaf cold-deciduous trees
       nDiag_biol06, & ! cbiol boreal conifer evergreen trees
       nDiag_biol07, & ! cbiol boreal broadleaf cold-deciduous trees
       nDiag_biol08, & ! cbiol boreal conifer cold-deciduous trees
       nDiag_biol09, & ! cbiol evergreen shrubs
       nDiag_biol10, & ! cbiol cold-deciduous shrubs
       nDiag_biol11, & ! cbiol warm (c4) grasses
       nDiag_biol12, & ! cbiol cool (c3) grasses
       nDiag_ynpp01, & ! ynpp tropical broadleaf evergreen trees
       nDiag_ynpp02, & ! ynpp tropical broadleaf drought-deciduous trees
       nDiag_ynpp03, & ! ynpp warm-temperate broadleaf evergreen trees
       nDiag_ynpp04, & ! ynpp temperate conifer evergreen trees
       nDiag_ynpp05, & ! ynpp temperate broadleaf cold-deciduous trees
       nDiag_ynpp06, & ! ynpp boreal conifer evergreen trees
       nDiag_ynpp07, & ! ynpp boreal broadleaf cold-deciduous trees
       nDiag_ynpp08, & ! ynpp boreal conifer cold-deciduous trees
       nDiag_ynpp09, & ! ynpp evergreen shrubs
       nDiag_ynpp10, & ! ynpp cold-deciduous shrubs
       nDiag_ynpp11, & ! ynpp warm (c4) grasses
       nDiag_ynpp12, & ! ynpp cool (c3) grasses
       nDiag_cmontp, & !coldest monthly temperature                         (C)
       nDiag_wmontp, & !warmest monthly temperature                         (C)
       nDiag_atogpp, & !annual total gpp for ecosystem                             (kg-c/m**2/yr)
       nDiag_toigpp, & !instantaneous gpp                           (mol-CO2 / m-2 / second)
       nDiag_fxcsol, & !instantaneous fine co2 flux from soil                           (mol-CO2 / m-2 / second)
       nDiag_mcsoil, & !instantaneous microbial co2 flux from soil         (mol-CO2 / m-2 / second)
       nDiag_cagcub, & !canopy average gross photosynthesis rate - broadleaf  (mol_co2 m-2 s-1)
       nDiag_cagcuc, & !canopy average gross photosynthesis rate - conifer    (mol_co2 m-2 s-1)
       nDiag_cagcls, & !canopy average gross photosynthesis rate - shrubs     (mol_co2 m-2 s-1)
       nDiag_cagcl4, & !canopy average gross photosynthesis rate - c4 grasses (mol_co2 m-2 s-1)
       nDiag_cagcl3, & !canopy average gross photosynthesis rate - c3 grasses (mol_co2 m-2 s-1)
       nDiag_cancub, & !canopy average net photosynthesis rate - broadleaf    (mol_co2 m-2 s-1)
       nDiag_cancuc, & !canopy average net photosynthesis rate - conifer      (mol_co2 m-2 s-1)
       nDiag_cancls, & !canopy average net photosynthesis rate - shrubs         (mol_co2 m-2 s-1)
       nDiag_cancl4, & !canopy average net photosynthesis rate - c4 grasses   (mol_co2 m-2 s-1)
       nDiag_cancl3, & !canopy average net photosynthesis rate - c3 grasses   (mol_co2 m-2 s-1)
       nDiag_cicoub, & !intercellular co2 concentration - broadleaf        (mol_co2/mol_air)
       nDiag_cicouc, & !intercellular co2 concentration - conifer        (mol_co2/mol_air)
       nDiag_cscoub, & !leaf boundary layer co2 concentration - broadleaf     (mol_co2/mol_air)
       nDiag_gscoub, & !upper canopy stomatal conductance - broadleaf                (mol_co2 m-2 s-1)
       nDiag_cscouc, & !leaf boundary layer co2 concentration - conifer         (mol_co2/mol_air)
       nDiag_gscouc, & !upper canopy stomatal conductance - conifer        (mol_co2 m-2 s-1)
       nDiag_cicols, & !intercellular co2 concentration - shrubs        (mol_co2/mol_air)
       nDiag_cicol3, & !intercellular co2 concentration - c3 plants        (mol_co2/mol_air)
       nDiag_cicol4, & !intercellular co2 concentration - c4 plants        (mol_co2/mol_air)
       nDiag_cscols, & !leaf boundary layer co2 concentration - shrubs          (mol_co2/mol_air)
       nDiag_gscols, & !lower canopy stomatal conductance - shrubs        (mol_co2 m-2 s-1)
       nDiag_cscol3, & !leaf boundary layer co2 concentration - c3 plants     (mol_co2/mol_air)
       nDiag_gscol3, & !lower canopy stomatal conductance - c3 grasses          (mol_co2 m-2 s-1)
       nDiag_cscol4, & !leaf boundary layer co2 concentration - c4 plants     (mol_co2/mol_air)
       nDiag_gscol4, & !lower canopy stomatal conductance - c4 grasses          (mol_co2 m-2 s-1)
       nDiag_tcthis, & !coldest monthly temperature of current year          (C)
       nDiag_twthis    !warmest monthly temperature of current year          (C)

  USE Options, ONLY: OCFLUX,rootmode,omlmodel,oml_hml0,SLABOCEAN,co2val,atmpbl,ICEMODEL,nClass,nAeros,indexchem,typechem,gbaco2,jull
  USE SlabOceanModel  , Only : GetOceanAlb
  USE Sfc_SeaIceFlux_WRF_Model  , Only : GetIceOceanAlb

  IMPLICIT NONE
SAVE

  PRIVATE

  REAL (KIND=r8), PARAMETER   :: z0ice  =       0.001e0_r8! 
  REAL(KIND=r8)   , PARAMETER   :: rgas   = 287.0_r8!  dry air gas constant (J deg^-1 kg^-1)
  REAL(KIND=r8)   , PARAMETER   :: kapa   = 0.2861328125_r8!  rgas/cp (unitless)
  REAL(KIND=r8)   , PARAMETER   :: cp     = rgas/kapa!  specific of heat of dry air at constant pressure (J deg^-1 kg^-1)
  REAL(KIND=r8)   , PARAMETER   :: hltm   = 2.52e6_r8            !  latent heat of vaporization (J kg^-1)
  REAL(KIND=r8)   , PARAMETER   :: stefan = 5.67e-8_r8           !  Stefan-Boltzmann constant
  REAL(KIND=r8), PARAMETER, PUBLIC :: MAPL_AIRMW  = 28.97_r8                  ! kg/Kmole
  REAL(KIND=r8), PARAMETER, PUBLIC :: MAPL_H2OMW  = 18.01_r8                  ! kg/Kmole

  REAL(KIND=r8), PARAMETER, PUBLIC :: MAPL_VIREPS = MAPL_AIRMW/MAPL_H2OMW-1.0_r8   ! --


  PUBLIC :: Ibis_Interface,Albedo_IBIS 

CONTAINS

  SUBROUTINE Ibis_Interface(intg                ,istrt             , &
       jdt                ,jb                   ,dtc3x             ,nCols             , &
       ktm                ,initlz               , &
       kt                 ,iswrad               ,ilwrad            ,kMax              , &
       tod                ,idatec               ,filta             ,epsflt            , &
       gt                 ,gq                   ,gu                ,gv                , &
       prsi               ,prsl                 ,phii              ,phil              , &
       gps                ,tmtx                 ,qmtx              ,umtx              , &
       zenith             ,colrad               ,fira2             ,xvisb             , &
       xvisd              ,xnirb                ,xnird             ,ppli              , &
       ppci               ,snow                 ,SwSfcUp, &
       tseam              ,tsea                 ,mskant            ,speedm            , &
       slrad              ,tsurf                ,qsurf             ,zorl              , &
       taux               ,tauy                 ,sens              ,evap              , &
       umom               ,vmom                 ,rmi               ,rhi               , &
       z0                 ,ustar                ,hc                ,hg                , &
       ec                 ,eg                   ,ts2               ,qs2               , &
       qsfc               ,tsfc                 ,z0sea             ,d                 , &
       cu                 ,imask                ,Ustarm            ,tgrd              , &
       roff               ,ect                  ,eci               ,egt               , &
       egi                ,egs                  ,rho               ,bstar             , &
       HML                ,HUML                 ,HVML              ,TSK               , &
       cldtot             ,ySwSfcNet            ,LwSfcNet          ,pblh              , &
       QCF                ,QCL                  ,sm0               ,mlsi              , &
       LwSfcDown          ,month                ,Mmlen             ,co2m              , &
       cflx               ,topog                ,dump             )

    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: jdt
    REAL(KIND=r8)   , INTENT(IN   ) :: tod
    INTEGER, INTENT(IN   ) :: idatec (4)  
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: jb
    INTEGER, INTENT(IN   ) :: ktm
    INTEGER, INTENT(IN   ) :: initlz
    INTEGER, INTENT(IN   ) :: kt
    INTEGER, INTENT(IN   ) :: kMax
    REAL(KIND=r8)   , INTENT(IN   ) :: dtc3x         ! model timestep (seconds)
    REAL(KIND=r8)   , INTENT(INOUT) :: gt     (nCols,kMax)          ! air temperature (K)
    REAL(KIND=r8)   , INTENT(INOUT) :: gq     (nCols,kMax)          ! specific humidity (kg_h2o/kg_air)
    REAL(KIND=r8)   , INTENT(IN   ) :: gu     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: gv     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: prsi   (nCols,kMax+1)  !     prsi     - real, pressure at layer interfaces             ix,levs+1  Pa
    REAL(KIND=r8)   , INTENT(IN   ) :: prsl   (nCols,kMax)    !     prsl     - real, mean layer presure                       ix,levs   Pa
    REAL(KIND=r8)   , INTENT(IN   ) :: phii   (nCols,kMax+1)  !===>  PHIH(K+1)  INPUT GEOPOTENTIAL @ EDGES  IN MKS units (m)
    REAL(KIND=r8)   , INTENT(IN   ) :: phil   (nCols,kMax)    !===>  PHIL(K)    INPUT GEOPOTENTIAL @ LAYERS IN MKS units (m)
    REAL(KIND=r8)   , INTENT(IN   ) :: gps    (nCols)          ! surface pressure (nPa)
    REAL(KIND=r8)   , INTENT(IN   ) :: xvisb  (nCols)!solad(i,1).Downward Surface shortwave fluxe visible beam    (cloudy)
    REAL(KIND=r8)   , INTENT(IN   ) :: xvisd  (nCols)!solai(i,1) !.Downward Surface shortwave fluxe visible diffuse (cloudy)
    REAL(KIND=r8)   , INTENT(IN   ) :: xnirb  (nCols)!solad(i,2).Downward Surface shortwave fluxe Near-IR beam    (cloudy)
    REAL(KIND=r8)   , INTENT(IN   ) :: xnird  (nCols)!solai(1,2) !.Downward Surface shortwave fluxe Near-IR diffuse (cloudy)
    REAL(KIND=r8)   , INTENT(IN   ) :: fira2  (nCols)          ! incoming ir flux (W m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: zenith (nCols)          ! cosine of solar zenith angle
    REAL(KIND=r8)   , INTENT(IN   ) :: colrad (nCols)           
    REAL(KIND=r8)   , INTENT(IN   ) :: ppli   (nCols)! Precipitation rate ( large scale )       (mm/s)
    REAL(KIND=r8)   , INTENT(IN   ) :: ppci   (nCols)! Precipitation rate ( cumulus )           (mm/s)
    REAL(KIND=r8)   , INTENT(IN   ) :: snow   (nCols)! snowfall rate (mm/s or kg m-2 s-1 of water)
    REAL(KIND=r8)   , INTENT(IN   ) :: SwSfcUp(nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: tseam  (nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: tsea   (nCols)

    INTEGER(KIND=i8), INTENT(IN   ) :: mskant(ncols)

    REAL(KIND=r8)   , INTENT(INOUT) :: speedm (nCols) ! wind speed (m s-1)
    REAL(KIND=r8)   , INTENT(INOUT) :: tmtx   (nCols,kMax,3)
    REAL(KIND=r8)   , INTENT(INOUT) :: qmtx   (nCols,kMax,3)
    REAL(KIND=r8)   , INTENT(INOUT) :: umtx   (nCols,kMax,4)
    REAL(KIND=r8)   , INTENT(IN   ) :: slrad  (nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: tsurf  (nCols) 
    REAL(KIND=r8)   , INTENT(IN   ) :: qsurf  (nCols) 
    REAL(KIND=r8)   , INTENT(INOUT) :: zorl   (nCols) 
    REAL(KIND=r8)   , INTENT(OUT  ) :: taux   (nCols)
    REAL(KIND=r8)   , INTENT(OUT  ) :: tauy   (nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: ts2    (nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: qs2    (nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: filta
    REAL(KIND=r8)   , INTENT(IN   ) :: epsflt
    INTEGER         , INTENT(IN   ) :: intg
    INTEGER         , INTENT(IN   ) :: istrt
    CHARACTER(len=*), INTENT(IN   ) :: iswrad
    CHARACTER(len=*), INTENT(IN   ) :: ilwrad
    REAL(KIND=r8)   , INTENT(INOUT) :: sens(1:nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: evap(1:nCols) 
    REAL(KIND=r8)   , INTENT(INOUT) :: umom  (1:nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: vmom  (1:nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: rmi   (1:nCols) 
    REAL(KIND=r8)   , INTENT(INOUT) :: rhi   (1:nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: tsfc  (1:nCols) 
    REAL(KIND=r8)   , INTENT(INOUT) :: qsfc  (1:nCols) 
    REAL(KIND=r8)   , INTENT(INOUT) :: z0    (1:nCols) 
    REAL(KIND=r8)   , INTENT(INOUT) :: ustar (1:nCols )
    REAL(KIND=r8)   , INTENT(OUT  ) :: hc (1:nCols)
    REAL(KIND=r8)   , INTENT(OUT  ) :: hg (1:nCols)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ec (1:nCols)
    REAL(KIND=r8)   , INTENT(OUT  ) :: eg (1:nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: z0sea (1:nCols)
    REAL(KIND=r8)   , INTENT(OUT  ) :: d(1:nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: cu(1:nCols)
    INTEGER(KIND=i8), INTENT(IN   ) :: imask    (nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: Ustarm(nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: tgrd(nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: roff(nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: ect(nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: eci(nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: egt(nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: egi(nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: egs(nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: rho(nCols)
    REAL(KIND=r8)   , INTENT(OUT  ) :: bstar(nCols)
    REAL(KIND=r8)    ,INTENT(IN OUT) ::dump(1:nCols,1:kMax )

    REAL(KIND=r8),    INTENT(INOUT) :: HML  (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: HUML (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: HVML (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: TSK  (ncols)
    REAL(KIND=r8)   , INTENT(IN   ) :: cldtot (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: ySwSfcNet (ncols)
    REAL(KIND=r8)   , INTENT(IN   ) :: LwSfcNet (ncols)
    REAL(KIND=r8)   , INTENT(IN   ) :: pblh (ncols)
    REAL(KIND=r8)   , INTENT(IN   ) :: QCF(nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: QCL(nCols,kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: sm0   (ncols,3)
    INTEGER(KIND=i8), INTENT(INOUT) :: mlsi  (ncols)
    REAL   (KIND=r8), INTENT(IN   ) :: LwSfcDown(1:nCols)
    INTEGER         , INTENT(IN   ) :: month(1:nCols)
    REAL   (KIND=r8), INTENT(IN   ) :: Mmlen (1:nCols)
    REAL   (KIND=r8), INTENT(IN   ) :: co2m  (1:nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: cflx (1:nCols,nClass+nAeros)
    REAL   (KIND=r8), INTENT(IN   ) :: topog (1:nCols)

    !  idatec(1)....hour(00/12)
    !  idatec(2)....month
    !  idatec(3)....day of month
    !  idatec(4)....year
    !   INCLUDE 'comatm.h'
    REAL(KIND=r8)   :: cond  (1:nCols)
    REAL(KIND=r8)   :: stor  (1:nCols)
    REAL(KIND=r8)         :: ta     (nCols)        ! air temperature (K)
    REAL(KIND=r8)         :: qa     (nCols)        ! specific humidity (kg_h2o/kg_air)
    REAL(KIND=r8)         :: ux     (nCols)
    REAL(KIND=r8)         :: uy     (nCols)
    REAL(KIND=r8)         :: psurf  (nCols)        ! surface pressure (Pa)
    REAL(KIND=r8)         :: solad  (nCols,nband) ! direct downward solar flux (W m-2)
    REAL(KIND=r8)         :: solai  (nCols,nband) ! diffuse downward solar flux (W m-2)
    REAL(KIND=r8)         :: fira   (nCols)        ! incoming ir flux (W m-2)
    REAL(KIND=r8)         :: coszen (nCols)        ! cosine of solar zenith angle
    REAL(KIND=r8)         :: raina  (nCols)        ! rainfall rate (mm/s or kg m-2 s-1)
    REAL(KIND=r8)         :: snowa  (nCols)        ! snowfall rate (mm/s or kg m-2 s-1 of water)
    REAL(KIND=r8)         :: ua     (nCols)        ! wind speed (m s-1)
    REAL(KIND=r8)         :: tmin   (nCols)        ! minimun soil temperature (K)
    REAL(KIND=r8)         :: tmax   (nCols)        ! maximun soil temperature (K)
    INTEGER, PARAMETER :: ndaypy=365  ! number of days per year(365or366)
    REAL(KIND=r8)         :: tgpptot(nCols) 
    REAL(KIND=r8)         :: disturbf(nCols)    ! local annual fire disturbance regime (m2/m2/yr)
    REAL(KIND=r8)         :: disturbo(nCols)    ! local fraction of biomass pool lost every year to disturbances other than fire
    REAL(KIND=r8)   :: fvapa    (nCols)         ! local ! downward h2o vapor flux between za & z12 at za (kg m-2 s-1)
    REAL(KIND=r8)   :: fsena    (nCols)         ! local ! downward sensible heat flux between za & z12 at za (W m-2)
    REAL(KIND=r8)   :: xsea       (1:nCols) 
    REAL(KIND=r8)   :: bstar1       (1:nCols)
    REAL(KIND=r8)   :: diag       (1:nCols)
    REAL(KIND=r8)   :: firb     (nCols)          ! local ! net upward ir radiation at reference
    REAL(KIND=r8)   :: GSW (nCols)
    REAL(KIND=r8)   :: GLW (nCols)
    REAL(kind=r8)   :: dlwflx (nCols)!    - real, total sky sfc downward lw flux ( w/m**2 )   im   !
    REAL(kind=r8)   :: sfcnsw (nCols)!     sfcnsw   - real, total sky sfc netsw flx into ground(w/m**2) im   !
    REAL(kind=r8)   :: sfcnlw (nCols)!     sfcnlw   - real, total sky sfc netlw flx into ground(w/m**2) im   !

    REAL(kind=r8)   :: sfcdsw (nCols)!     sfcdsw   - real, total sky sfc downward sw flux ( w/m**2 )   im   !
    REAL(kind=r8)   :: zorl_landice(nCols)
    REAL(KIND=r8)   :: tprcp_landice(nCols)
    REAL(kind=r8)   :: tsoi_landice   (1:nCols,1:nsoilay)! ,    INTENT(in	)  Var_StepP    , &
    REAL(kind=r8)   :: sheleg_landice(nCols)
    REAL(kind=r8)   :: tsea_landice(nCols)
    REAL(kind=r8)   :: qsurf_landice(nCols)
    REAL(kind=r8)   :: evap_landice(nCols)
    REAL(kind=r8)   :: hflx_landice(nCols)
    REAL(kind=r8)   :: bstar_landice(nCols)
    REAL(kind=r8)   ::taux_landice (nCols)
    REAL(kind=r8)   ::tauy_landice (nCols)
    REAL(kind=r8)   :: ustar_landice(nCols)
    REAL(kind=r8)   :: co2initm(nCols) !mol/mol
    REAL(kind=r8)   :: co2diag(nCols) !mol/mol
    REAL(kind=r8)   :: za(nCols)
    REAL(kind=r8)   :: bps(nCols)
    ! atmospheric level za (W m-2)

    !
    ! Arguments (input)
    !
    LOGICAL :: InitMod
    !REAL(KIND=r8)    :: calday     ! current julian day (1-365.99)
    !  INTEGER :: iday     ! current day in month (1-31 or 1-30, or 1-28)
    INTEGER :: iday         ! day number  (passed in)
    !INTEGER :: idayprev   ! day in month of previous timestep

    !   integer:: imonth         ! current month (1 - 12)
    INTEGER :: imonth         ! month number (passed in)   
    INTEGER :: imonthprev ! month of previous timestep 

    INTEGER :: iyear
    !INTEGER :: iyearprev  ! year of previous timestep 

    !INTEGER :: idayout         ! write out daily output
    !INTEGER :: imonthout! write out daily output
    !INTEGER :: iyearout 
    INTEGER :: nLndPts
    INTEGER :: i,k,ind,nint,IntSib,itr,itrac,kk
    REAL (KIND=r8), PARAMETER   :: tice   =      271.16e0_r8! constant tice
    REAL(r8), PARAMETER :: amdc  = 0.658114_r8! Molecular weight of dry air / carbon dioxide
    hc    =0.0_r8
    hg    =0.0_r8
    ec    =0.0_r8
    eg    =0.0_r8
    d     =0.0_r8
    bstar =0.0_r8
    co2diag=0.0_r8
    !  idatec(1)....hour(00/12)
    !  idatec(2)....month
    !  idatec(3)....day of month
    !  idatec(4)....year
    iday      = idatec   (3) ! current day in month (1-31 or 1-30, or 1-28)
    imonth    = idatec   (2) ! current month (1 - 12)  
    iyear     = idatec   (4) ! current year 
    !co2init = co2val*0.000001_r8

    imonthprev= idateprev(2) ! month of previous timestep 
    !idayprev  = idateprev(3) ! day in month of previous timestep
    !iyearprev = idateprev(4) ! year of previous timestep 

    !calday = julday (imonth,iday,iyear,tod)
    !CALL MsgOne(h," IBIS: IBIS Tendencies")
    !IF(myId ==0)WRITE(*,*)imonthprev
    nLndPts=0
    DO i=1,nCols
       GSW(I) = xvisb  (i)+xvisd  (i)+xnirb  (i)+xnird  (i)
       GLW(I) = fira2  (i)
       sfcnsw (i  )  = (xvisb  (i)+xvisd  (i)+xnirb  (i)+xnird  (i) ) - SwSfcUp(i) 
       sfcnlw (i  )  = LwSfcDown(i)                                   - (stefan*tsurf(i)**4)
       IF (iMask(i) >= 1_i8) THEN
          nLndPts=nLndPts+1 
          mlsi(i) = 1_i8   !add solange 13-11-2012
          !  ---  parameter constants for gas volume mixing ratioes 
          ! Gases other than N2, O2, Ar, and H2O are present in the atmosphere at extremely low 
          ! concentrations and are called trace gases. Despite their low concentrations, these 
          ! trace gases can be of critical importance for the greenhouse effect, the ozone layer, 
          ! smog, and other environmental issues. Mixing ratios of trace gases are commonly given
          ! in units of parts per million volume ( ppmv or simply ppm), parts per billion volume ( ppbv or ppb), 
          ! or parts per trillion volume ( pptv or ppt); 1 ppmv = 1x10-6 mol/mol, 1 ppbv = 1x10-9 mol/mol, 
          ! and 1 pptv = 1x10-12 mol/mol. For example, the present-day CO2 concentration is 365 ppmv (365x10-6 mol/mol).
          !
          !
          !            1 gas volume X
          !  1 ppm = --------------------= 1e-6 mol/mol
          !             1e6 air volumes
          !  
          !  CH4
          !   
          !  1 mole C= 12.01 g
          !  4 moles H= 4.04 g
          !  Therefore, 1 mole of CH4 = 12.01g + 4.04g = 16.05g/mol
          !
          !  mwdry       = 28.966_R8       ! molecular weight dry air ~ kg/kmole ~ g/mol


          ! mol/mol               kg/kg    
          IF(isimco2 == 0)THEN
             co2initm  (nLndPts  )= co2val*0.000001_r8  !mol/mol
          ELSE
             co2initm  (nLndPts  )= co2m(i)*amdc    !convert kg/kg  to mol/mol
          END IF
          !bps   (ncount)=sigki(1)
          bps   (nLndPts)    = (prsi(i,1)/(prsi(i,2)))**(gasr/cp)

          ta     (nLndPts  ) = gt(i,1)
          qa     (nLndPts  ) = MAX(gq(i,1),0.00000012_r8) 
          za     (nLndPts  ) = MAX((phii(i,2) - phii(i,1))*0.5_r8,0.5_r8)  
          !IF(topog(i) > )
          ux     (nLndPts  ) = gu(i,1) /SIN(colrad(i))
          IF( ux(nLndPts  ) <= 0_r8 .and. ux(nLndPts) > -0.1_r8)ux(nLndPts  )=-0.1_r8
          IF( ux(nLndPts  ) >= 0_r8 .and. ux(nLndPts) <  0.1_r8)ux(nLndPts  )= 0.1_r8
          uy     (nLndPts  ) = gv(i,1) /SIN(colrad(i))
          IF( uy(nLndPts  ) <= 0_r8 .and. uy(nLndPts) >-0.1_r8)uy(nLndPts  )=-0.1_r8
          IF( uy(nLndPts  ) >= 0_r8 .and. uy(nLndPts) < 0.1_r8)uy(nLndPts  )= 0.1_r8
          speedm (nLndPts  ) = SQRT(( ux(nLndPts))**2  + ( uy(nLndPts))**2)
          ua     (nLndPts  ) = speedm (nLndPts) 

          psurf  (nLndPts  ) = gps    (i)*100.00_r8          ! surface pressure hPa -> Pa
          solad  (nLndPts,1) = xvisb  (i)  !! number of solar radiation wavebands : vis, nir
          solai  (nLndPts,1) = xvisd  (i)
          solad  (nLndPts,2) = xnirb  (i)
          solai  (nLndPts,2) = xnird  (i)
          fira   (nLndPts  ) = fira2  (i)      
          coszen (nLndPts  ) = zenith (i)
          raina  (nLndPts  ) = (ppli(i) + ppci(i)-snow   (i) )*(1.0_r8/dtc3x)  !convert mm/s to m/s
          tprcp_landice(nLndPts) =(ppli(i) + ppci(i) -snow(i)) !mm
          snowa  (nLndPts  ) = snow   (i) *(1.0_r8/dtc3x)  !convert mm/s to m/s
          q34    (nLndPts,jb) =  MAX(q34(nLndPts,jb),0.00000012_r8) 
          q12    (nLndPts,jb) =  MAX(q34(nLndPts,jb),0.00000012_r8) 
          dlwflx (nLndPts  )  = LwSfcDown(i)
          sfcnsw (nLndPts  )  = ySwSfcNet(i)
          sfcdsw (nLndPts  )  =xvisb  (i)+xvisd  (i)+xnirb  (i)+xnird  (i)
          zorl_landice (nLndPts  )  =zorl (i  )
          sheleg_landice (nLndPts  )  = sheleg(i,jb)
       END IF
    END DO
    InitMod = (initlz >= 0 .AND. ktm == -1 .AND. kt == 0 )
    IF(InitMod)THEN
       nLndPts=0
       DO i=1,nCols
          !HML  (i) = oml_hml0 + 2.0*SQRT((gu(i,1) /SIN(colrad(i)))**2  + (gv(i,1) /SIN(colrad(i)))**2) 
          IF( iMask(i) >= 1_i8)THEN 
             nLndPts=nLndPts+1
             tu0(nLndPts,jb) = ta     (nLndPts  )
             tu (nLndPts,jb) = ta     (nLndPts  )
             tum(nLndPts,jb) = ta     (nLndPts  )
             ts0 (nLndPts,jb) = ta     (nLndPts  )
             ts (nLndPts,jb) = ta     (nLndPts  )
             tsm (nLndPts,jb) = ta     (nLndPts  )
             tl0 (nLndPts,jb) = ta     (nLndPts  )
             tl (nLndPts,jb) = ta     (nLndPts  )
             tlm (nLndPts,jb) = ta     (nLndPts  )

             t12(nLndPts,jb) = ta     (nLndPts  )

             t34(nLndPts,jb) = ta     (nLndPts  )       

             tg (nLndPts,jb) = tgm     (nLndPts,jb) 
             wliqu(nLndPts,jb) = wliqum  (nLndPts,jb) 
             wliqs(nLndPts,jb) = wliqsm   (nLndPts,jb) 
             wliql(nLndPts,jb) = wliqlm   (nLndPts,jb) 
             wsnou(nLndPts,jb) = wsnoum   (nLndPts,jb)          
             wsnos(nLndPts,jb) = wsnosm   (nLndPts,jb) 
             wsnol(nLndPts,jb) = wsnolm   (nLndPts,jb) 
             q12(nLndPts,jb) = qa     (nLndPts  )       
             q34(nLndPts,jb) = qa     (nLndPts  )       
             DO k=1,nsoilay 
                tsoi0(nLndPts,k,jb) =tsoim (nLndPts,k,jb)
                tsoim(nLndPts,k,jb) =tsoim (nLndPts,k,jb)
                tsoi (nLndPts,k,jb) =tsoim (nLndPts,k,jb)
             END DO
             DO k=1,nsoilay
                wsoi0(nLndPts,k,jb)  = wsoim (nLndPts,k,jb)
                wsoim(nLndPts,k,jb)  = wsoim (nLndPts,k,jb)
                wsoi (nLndPts,k,jb)  = wsoim (nLndPts,k,jb)
             END DO
             DO k=1,nsoilay
                wisoi0(nLndPts,k,jb)  = wisoi0 (nLndPts,k,jb)
                wisoim(nLndPts,k,jb)  = wisoim (nLndPts,k,jb)
                wisoi (nLndPts,k,jb)  = wisoim  (nLndPts,k,jb)
             END DO
             DO k=1,nsnolay
                tsno0    (nLndPts,k,jb) = tsno0  (nLndPts,k,jb)
                tsno     (nLndPts,k,jb) = tsno   (nLndPts,k,jb)
                tsnom    (nLndPts,k,jb) = tsnom  (nLndPts,k,jb)

                hsno      (nLndPts,k,jb) = hsno    (nLndPts,k,jb) 

             END DO

          END IF
       END DO
    END IF
    nLndPts=0
    DO i=1,nCols
       IF( iMask(i) >= 1_i8)THEN 
          nLndPts=nLndPts+1
          !          wliqmin  =0.00_r8         !local !
          !          wsnomin  =0.00_r8         !local !
          !          grunof   (nLndPts,jb)=0.00_r8 !local !
          !          ginvap   (nLndPts,jb)=0.00_r8 !local !
          !          gsuvap   (nLndPts,jb)=0.00_r8 !local !
          !          gtrans   (nLndPts,jb)=0.00_r8 !local !
          !          gtransu  (nLndPts,jb)=0.00_r8 !local !
          !          gtransl  (nLndPts,jb)=0.00_r8 !local !
          !          gdrain   (nLndPts,jb)=0.00_r8 !local !
          !          gadjust  (nLndPts,jb)=0.00_r8 !local !
          !          topparu  (nLndPts,jb)=0.00_r8 !local !
          !          topparl  (nLndPts,jb)=0.00_r8 !local !
          !          su       (nLndPts,jb)=0.00_r8 !local !
          !          ss       (nLndPts,jb)=0.00_r8 !local !
          !          sl       (nLndPts,jb)=0.00_r8 !local !
          !          agcub    (nLndPts,jb)=0.00_r8!local !
          !          agcuc    (nLndPts,jb)=0.00_r8!local !
          !          ancub    (nLndPts,jb)=0.00_r8!local !
          !          ancuc    (nLndPts,jb)=0.00_r8!local !
          !          totcondub(nLndPts,jb)=0.00_r8!local !
          !          totconduc(nLndPts,jb)=0.00_r8!local !
          !          agcls    (nLndPts,jb)=0.00_r8!local !
          !          agcl4    (nLndPts,jb)=0.00_r8!local !
          !          agcl3    (nLndPts,jb)=0.00_r8!local !
          !          ancls    (nLndPts,jb)=0.00_r8!local !
          !          ancl4    (nLndPts,jb)=0.00_r8!local !
          !          ancl3    (nLndPts,jb)=0.00_r8!local !
          !          totcondls(nLndPts,jb)=0.00_r8!local !
          !          totcondl3(nLndPts,jb)=0.00_r8!local !
          !         totcondl4(nLndPts,jb)=0.00_r8!local !
          !          consoi   (nLndPts,1:nsoilay,jb)=0.00_r8!local ! 
          !          qglif    (nLndPts,1:4,jb)=0.00_r8!local ! 
          !          hvasug   (nLndPts,jb)=0.00_r8!local ! 
          !          hvasui   (nLndPts,jb)=0.00_r8!local ! 
          !          stressl  (nLndPts,1:nsoilay,jb)=0.00_r8!local ! 
          !          stressu  (nLndPts,1:nsoilay,jb)=0.00_r8!local ! 
          !          stresstl (nLndPts,jb)=0.00_r8!local !  
          !          stresstu (nLndPts,jb)=0.00_r8!local !         
          !          upsoiu   (nLndPts,1:nsoilay,jb)=0.00_r8!local !         
          !          upsoil   (nLndPts,1:nsoilay,jb)=0.00_r8!local !         
          !          heatg    (nLndPts,jb)=0.00_r8!local !                  
          !          heati    (nLndPts,jb)=0.00_r8!local !                  
          !          asurd    (nLndPts,1:nband,jb)=0.00_r8!local !                 
          !          asuri    (nLndPts,1:nband,jb)=0.00_r8!local !                 
          !          tmm     (ncount)  = tmgm  (ncount,latco)

          tg (nLndPts,jb) = tgm     (nLndPts,jb) 
          ts (nLndPts,jb) = tsm     (nLndPts,jb)
          tl (nLndPts,jb) = tlm     (nLndPts,jb)
          ti (nLndPts,jb) = tim     (nLndPts,jb)

          wliqu(nLndPts,jb) = wliqum  (nLndPts,jb) 
          wliqs(nLndPts,jb) = wliqsm   (nLndPts,jb) 
          wliql(nLndPts,jb) = wliqlm   (nLndPts,jb) 
          wsnou(nLndPts,jb) = wsnoum   (nLndPts,jb)          
          wsnos(nLndPts,jb) = wsnosm   (nLndPts,jb)  
          wsnol(nLndPts,jb) = wsnolm   (nLndPts,jb) 

          !          tm0     (ncount)  = tmgp  (ncount,latco)
          DO k=1,nsoilay 
!             tsoi0(nLndPts,k,jb) = tsoi0 (nLndPts,k,jb)
!             tsoim(nLndPts,k,jb) = tsoim (nLndPts,k,jb)
             tsoi (nLndPts,k,jb) = tsoim  (nLndPts,k,jb)
             tsoi_landice(nLndPts,k)= tsoi (nLndPts,k,jb) 
          END DO
          DO k=1,nsoilay
!             wsoi0(nLndPts,k,jb)  = wsoi0 (nLndPts,k,jb)
!             wsoim(nLndPts,k,jb)  = wsoim (nLndPts,k,jb)
             wsoi (nLndPts,k,jb)  = wsoim  (nLndPts,k,jb)
          END DO
          DO k=1,nsoilay
!             wsoi0(nLndPts,k,jb)  = wsoi0 (nLndPts,k,jb)
!             wsoim(nLndPts,k,jb)  = wsoim (nLndPts,k,jb)
             wisoi (nLndPts,k,jb)  = wisoim  (nLndPts,k,jb)
          END DO
          DO k=1,nsnolay
             tsno     (nLndPts,k,jb) = tsnom   (nLndPts,k,jb)
          END DO

       END IF
    END DO


    IF(InitMod)THEN
       nint=2
       IntSib=5
    ELSE
       nint=1
       IntSib=1
    END IF
    IF(TRIM(iswrad).NE.'NON'.AND.TRIM(ilwrad).NE.'NON') THEN
       IF(InitMod .AND. nLndPts >= 1)THEN
          DO ind=1,nint
             nLndPts=0
             DO i=1,nCols
                IF (iMask(i) >= 1_i8) THEN
                   nLndPts=nLndPts+1
                   !
                   !     precipitation
                   !
                   raina  (nLndPts  ) =0.0e0_r8*1.0e-3_r8  !convert mm/s to m/s
                END IF
             END DO
             DO itr=1,IntSib
                CALL IbisDrv (tod,pi,stef,vonk,grav,tmelt,hfus,hvap,hsub,ch2o,cice,cair,cvap,rair,rvap,cappa, &
                     rhow,nLndPts,nband,nsoilay,nsnolay,npft,epsilon,dtc3x,doalb,&
                     ginvap       (1:nLndPts,jb),gsuvap      (1:nLndPts,jb),gtrans         (1:nLndPts,jb), &
                     gtransu      (1:nLndPts,jb),gtransl        (1:nLndPts,jb),grunof         (1:nLndPts,jb),gdrain           (1:nLndPts,jb), &
                     gadjust      (1:nLndPts,jb),a10scalparamu(1:nLndPts,jb),a10daylightu(1:nLndPts,jb),a10scalparaml(1:nLndPts,jb), &
                     a10daylightl (1:nLndPts,jb),vmax_pft        (1:npft)    ,tau15                     ,kc15                       , &
                     ko15,cimax,gammaub,alpha3,theta3,beta3,coefmub,coefbub,gsubmin,gammauc,coefmuc,coefbuc, &
                     gsucmin,gammals,coefmls,coefbls,gslsmin,gammal3,coefml3,coefbl3, &
                     gsl3min,gammal4,alpha4,theta4,beta4,coefml4,coefbl4,gsl4min,bps(1:nLndPts), &
                     wliqu        (1:nLndPts,jb),wliqumax                    ,wsnou         (1:nLndPts,jb),wsnoumax                 , &
                     tu              (1:nLndPts,jb),wliqs        (1:nLndPts,jb),wliqsmax                     ,wsnos           (1:nLndPts,jb), &
                     wsnosmax                  ,ts                (1:nLndPts,jb),wliql         (1:nLndPts,jb),wliqlmax                 , &  
                     wsnol        (1:nLndPts,jb),wsnolmax                    ,tl          (1:nLndPts,jb),topparu           (1:nLndPts,jb), &
                     topparl      (1:nLndPts,jb),fl                (1:nLndPts,jb),fu          (1:nLndPts,jb),lai      (1:nLndPts,1:2,jb), &
                     sai          (1:nLndPts,1:2,jb),rhoveg      (1:nband,1:2),tauveg        (1:nband,1:2),orieh    (1:2)               , &
                     oriev    (1:2)           ,wliqmin                    ,wsnomin                     ,t12           (1:nLndPts,jb), &
                     tdripu                   ,tblowu                    ,tdrips                     ,tblows                       , &
                     t34              (1:nLndPts,jb),tdripl                    ,tblowl,ztop     (1:nLndPts,1:2,jb),za     (1:nLndPts), & 
                     alaiml                   ,zbot     (1:nLndPts,1:2,jb),alaimu                     ,froot         (1:nLndPts,1:nsoilay,1:2,jb), &
                     q34              (1:nLndPts,jb),q12          (1:nLndPts,jb),su          (1:nLndPts,jb),cleaf                       , &
                     dleaf        (1:2)          ,ss                (1:nLndPts,jb),cstem                     ,dstem                       , &
                     sl              (1:nLndPts,jb),cgrass                    ,ciub         (1:nLndPts,jb),ciuc           (1:nLndPts,jb), &
                     exist (1:nLndPts,1:npft,jb),csub         (1:nLndPts,jb),gsub         (1:nLndPts,jb),csuc           (1:nLndPts,jb), &
                     gsuc              (1:nLndPts,jb),agcub        (1:nLndPts,jb),agcuc         (1:nLndPts,jb),ancub           (1:nLndPts,jb), &
                     ancuc        (1:nLndPts,jb),totcondub        (1:nLndPts,jb),totconduc   (1:nLndPts,jb),cils           (1:nLndPts,jb), &
                     cil3              (1:nLndPts,jb),cil4         (1:nLndPts,jb),csls         (1:nLndPts,jb),gsls           (1:nLndPts,jb), &
                     csl3              (1:nLndPts,jb),gsl3         (1:nLndPts,jb),csl4         (1:nLndPts,jb),gsl4           (1:nLndPts,jb), &
                     agcls        (1:nLndPts,jb),agcl4        (1:nLndPts,jb),agcl3         (1:nLndPts,jb),ancls           (1:nLndPts,jb), &
                     ancl4        (1:nLndPts,jb),ancl3        (1:nLndPts,jb),totcondls   (1:nLndPts,jb),totcondl3    (1:nLndPts,jb), &
                     totcondl4    (1:nLndPts,jb),chu(1:nVegClass),chs(1:nVegClass),chl(1:nVegClass),frac  (1:nLndPts,1:npft,jb),tlsub          (1:nLndPts,jb),z0sno,rhos, &
                     consno                   ,hsnotop                    ,hsnomin                     ,fimin                       , &
                     fimax                    ,fi                (1:nLndPts,jb),tsno(1:nLndPts,1:nsnolay,jb),hsno(1:nLndPts,1:nsnolay,jb),&
                     sand(1:nLndPts,1:nsoilay,jb),clay(1:nLndPts,1:nsoilay,jb),poros(1:nLndPts,1:nsoilay,jb),wsoi(1:nLndPts,1:nsoilay,jb), &
                     wisoi(1:nLndPts,1:nsoilay,jb),consoi(1:nLndPts,1:nsoilay,jb),zwpmax             ,wpud(1:nLndPts,jb)         , &
                     wipud        (1:nLndPts,jb),wpudmax                    ,qglif   (1:nLndPts,1:4,jb),tsoi(1:nLndPts,1:nsoilay,jb), &
                     hvasug       (1:nLndPts,jb),hvasui        (1:nLndPts,jb),albsav         (1:nLndPts,jb),albsan           (1:nLndPts,jb), &
                     tg              (1:nLndPts,jb),ti                (1:nLndPts,jb),z0soi(1:nLndPts,jb)       ,swilt(1:nLndPts,1:nsoilay,jb), &
                     sfield(1:nLndPts,1:nsoilay,jb),stressl(1:nLndPts,1:nsoilay,jb),stressu(1:nLndPts,1:nsoilay,jb),stresstl(1:nLndPts,jb), &
                     stresstu (1:nLndPts,jb)    ,csoi(1:nLndPts,1:nsoilay,jb),rhosoi(1:nLndPts,1:nsoilay,jb),hsoi(1:nLndPts,1:nsoilay+1,jb)   , &
                     suction(1:nLndPts,1:nsoilay,jb),bex (1:nLndPts,1:nsoilay,jb),upsoiu(1:nLndPts,1:nsoilay,jb),upsoil(1:nLndPts,1:nsoilay,jb), &
                     heatg   (1:nLndPts,jb),heati   (1:nLndPts,jb),hydraul(1:nLndPts,1:nsoilay,jb),porosflo(1:nLndPts,1:nsoilay,jb), &
                     ibex(1:nLndPts,1:nsoilay,jb),bperm(1:nLndPts,jb)          ,hflo      (1:nLndPts,1:nsoilay+1,jb),ta       (1:nLndPts)       , &
                     asurd (1:nLndPts,1:nBand,jb),asuri(1:nLndPts,1:nBand,jb),coszen(1:nLndPts)              ,solad (1:nLndPts,1:nBand)  , &
                     solai (1:nLndPts,1:nBand)  ,fira      (1:nLndPts)      ,raina (1:nLndPts)             ,qa    (1:nLndPts)               , &
                     psurf    (1:nLndPts)          ,snowa     (1:nLndPts)      ,ua    (1:nLndPts)             ,o2conc                       , &
                     co2conc  (1:nLndPts,jb) ,td         (1:nLndPts,jb)            ,vzero    (1:nLndPts,jb)   ,ndaypy                       , &
                     nppdummy(1:nLndPts,1:npft,jb) ,cbiow (1:nLndPts,1:npft,jb),sapfrac(1:nLndPts,jb)      ,cbior(1:nLndPts,1:npft,jb), &
                     tco2root(1:nLndPts,jb)          ,tneetot(1:nLndPts,jb)            ,tco2mic (1:nLndPts,jb)       ,a10td    (1:nLndPts,jb)    , &
                     a10ancub(1:nLndPts,jb)          ,a10ancuc(1:nLndPts,jb)     ,a10ancls(1:nLndPts,jb)    ,a10ancl3(1:nLndPts,jb)     , &
                     a10ancl4(1:nLndPts,jb)          ,ndtimes(1:nLndPts,jb), &
                     adrain  (1:nLndPts,jb)          ,adsnow  (1:nLndPts,jb),tnpptot(1:nLndPts,jb),adaet   (1:nLndPts,jb)    ,adtrunoff(1:nLndPts,jb)    , &
                     adsrunoff(1:nLndPts,jb)    ,addrainage(1:nLndPts,jb)   ,adrh    (1:nLndPts,jb)    ,adsnod   (1:nLndPts,jb)    , &
                     adsnof   (1:nLndPts,jb)    ,adwsoi    (1:nLndPts,jb)   ,adtsoi  (1:nLndPts,jb)    ,adwisoi  (1:nLndPts,jb)    , &
                     adtlaysoi(1:nLndPts,jb)    ,adwlaysoi (1:nLndPts,jb)   ,adwsoic (1:nLndPts,jb)    ,adtsoic  (1:nLndPts,jb)    , &
                     adco2mic (1:nLndPts,jb)    ,adco2root (1:nLndPts,jb)   ,adco2soi(1:nLndPts,jb)    ,adco2ratio(1:nLndPts,jb)   , &
                     adnmintot(1:nLndPts,jb)    ,decompl   (1:nLndPts,jb)   ,decomps (1:nLndPts,jb)    ,tnmin    (1:nLndPts,jb)    , &
                     ndaypm                     ,nmtimes(1:nLndPts,jb),amrain   (1:nLndPts,jb)    , &
                     amsnow    (1:nLndPts,jb)   ,amaet     (1:nLndPts,jb)   ,amtrunoff(1:nLndPts,jb)   ,amsrunoff(1:nLndPts,jb)    , &
                     amdrainage(1:nLndPts,jb)   ,amtemp    (1:nLndPts,jb)   ,amqa     (1:nLndPts,jb)    , &
                     amsolar   (1:nLndPts,jb)   ,amirup   (1:nLndPts,jb)   ,amirdown (1:nLndPts,jb)    , &
                     amsens    (1:nLndPts,jb)   ,amlatent  (1:nLndPts,jb)   ,amlaiu   (1:nLndPts,jb)   ,amlail   (1:nLndPts,jb)    , &
                     amtsoi    (1:nLndPts,jb)   ,amwsoi    (1:nLndPts,jb)   ,amwisoi  (1:nLndPts,jb)   ,amvwc    (1:nLndPts,jb)    , &
                     amawc     (1:nLndPts,jb)   ,amsnod    (1:nLndPts,jb)   ,amsnof   (1:nLndPts,jb)   ,amnpp(1:nLndPts,1:npft,jb) , &
                     amnpptot  (1:nLndPts,jb)   ,amco2mic  (1:nLndPts,jb)   ,amco2root(1:nLndPts,jb)   ,amco2soi (1:nLndPts,jb)    , &
                     amco2ratio(1:nLndPts,jb)   ,amneetot  (1:nLndPts,jb)   ,amnmintot(1:nLndPts,jb)   ,amalbedo (1:nLndPts,jb)    , &
                     amtsoil(1:nLndPts,1:nsoilay,jb) ,amwsoil(1:nLndPts,1:nsoilay,jb),amwisoil(1:nLndPts,1:nsoilay,jb), nytimes(1:nLndPts,jb), &
                     aysolar   (1:nLndPts,jb)   ,ayirup    (1:nLndPts,jb)   ,ayirdown (1:nLndPts,jb)   ,aysens        (1:nLndPts,jb)   , &
                     aylatent  (1:nLndPts,jb)   ,ayprcp    (1:nLndPts,jb)   ,ayaet    (1:nLndPts,jb)   ,aytrans        (1:nLndPts,jb)   , &
                     aytrunoff (1:nLndPts,jb)   ,aysrunoff (1:nLndPts,jb)   ,aydrainage(1:nLndPts,jb)  ,aydwtot        (1:nLndPts,jb)   , & 
                     aywsoi    (1:nLndPts,jb)   ,aywisoi   (1:nLndPts,jb)   ,aytsoi   (1:nLndPts,jb)   ,ayvwc        (1:nLndPts,jb)   , &
                     ayawc     (1:nLndPts,jb)   ,aystresstu(1:nLndPts,jb)   ,aystresstl(1:nLndPts,jb)  ,aygpp(1:nLndPts,1:npft,jb) , &
                     aygpptot  (1:nLndPts,jb)   ,aynpp(1:nLndPts,1:npft,jb) ,aynpptot (1:nLndPts,jb)   ,ayco2mic  (1:nLndPts,jb)   , &
                     ayco2root (1:nLndPts,jb)   ,ayco2soi  (1:nLndPts,jb)   ,ayneetot (1:nLndPts,jb)   ,ayrootbio (1:nLndPts,jb)   , &
                     aynmintot (1:nLndPts,jb)   ,ayalit    (1:nLndPts,jb)   ,ayblit   (1:nLndPts,jb)   ,aycsoi        (1:nLndPts,jb)   , & 
                     aycmic    (1:nLndPts,jb)   ,ayanlit   (1:nLndPts,jb)   ,aybnlit  (1:nLndPts,jb)   ,aynsoi        (1:nLndPts,jb)   , &
                     ayalbedo  (1:nLndPts,jb)   ,totalit   (1:nLndPts,jb)   ,totrlit  (1:nLndPts,jb)   ,totcsoi  (1:nLndPts,jb)    , &
                     totcmic   (1:nLndPts,jb)   , &
                     totanlit  (1:nLndPts,jb)   ,totrnlit  (1:nLndPts,jb)   ,totnsoi  (1:nLndPts,jb)   ,totnmic        (1:nLndPts,jb)   , &
                     totlit    (1:nLndPts,jb)   ,totfall   (1:nLndPts,jb)   ,totnlit  (1:nLndPts,jb)   ,firefac        (1:nLndPts,jb)   , &
                     wtot      (1:nLndPts,jb)   ,storedn   (1:nLndPts,jb)   ,yrleach  (1:nLndPts,jb)   ,ynleach        (1:nLndPts,jb)   , & 
                     falll     (1:nLndPts,jb)   ,fallr     (1:nLndPts,jb)   ,fallw    (1:nLndPts,jb)   ,clitlm        (1:nLndPts,jb)   , &
                     clitls    (1:nLndPts,jb)   ,clitrm    (1:nLndPts,jb)   ,clitrs   (1:nLndPts,jb)   ,clitwm        (1:nLndPts,jb)   , &
                     clitws    (1:nLndPts,jb)   ,csoislop  (1:nLndPts,jb)   ,csoislon (1:nLndPts,jb)   ,csoipas        (1:nLndPts,jb)   , &
                     clitll    (1:nLndPts,jb)   ,clitrl    (1:nLndPts,jb)   ,clitwl   (1:nLndPts,jb)   ,tc        (1:nLndPts,jb)   , &
                     agddu     (1:nLndPts,jb)   ,tempu     (1:nLndPts,jb)   ,agddl    (1:nLndPts,jb)   ,templ        (1:nLndPts,jb)   , &
                     dropu     (1:nLndPts,jb)   ,dropls    (1:nLndPts,jb)   ,dropl4   (1:nLndPts,jb)   ,dropl3        (1:nLndPts,jb)   , &
                     plai(1:nLndPts,1:npft,jb)  ,iday,imonth,iyear,iyear0,isimveg,spinmax, &
                     amts2     (1:nLndPts,jb)   ,amtransu  (1:nLndPts,jb)  ,amtransl (1:nLndPts,jb),   amsuvap  (1:nLndPts,jb)        ,&
                     aminvap   (1:nLndPts,jb)   ,ux             (1:nLndPts)      ,uy       (1:nLndPts)      ,taux        (1:nLndPts)      , &
                     tauy           (1:nLndPts),ts2(1:nLndPts),qs2(1:nLndPts),deltat   (1:nLndPts,jb)   ,gdd0        (1:nLndPts,jb)   , &  
                     gdd0this  (1:nLndPts,jb)   ,tcthis    (1:nLndPts,jb)   ,twthis   (1:nLndPts,jb)   ,tcmin        (1:nLndPts,jb)   , &
                     gdd5           (1:nLndPts,jb)   ,gdd5this  (1:nLndPts,jb)   ,TminL    (1:npft)       ,TminU        (1:npft)       , &
                     Twarm     (1:npft)          ,GDD       (1:npft)            ,aleaf    (1:npft)       ,awood        (1:npft)       , &
                     cbiol(1:nLndPts,1:npft,jb) ,aroot     (1:npft)            ,disturbf (1:nLndPts)      ,disturbo  (1:nLndPts)      , &
                     specla    (1:npft)          ,biomass(1:nLndPts,1:npft,jb),totlaiu (1:nLndPts,jb)   ,totlail        (1:nLndPts,jb)   , &  
                     totbiou   (1:nLndPts,jb)   ,totbiol   (1:nLndPts,jb)   ,woodnorm                     ,vegtype0  (1:nLndPts,jb)   , &
                     tauwood0  (1:npft)          ,tauwood   (1:nLndPts,1:npft,jb)            ,tauleaf  (1:npft)       ,tauroot        (1:npft)       , &
                     xminlai,cdisturb  (1:nLndPts,jb)   ,ayanpp(1:nLndPts,1:npft,jb),ayanpptot(1:nLndPts,jb),jdt,&
                     imonthprev,isimco2,isimfire, co2initm(1:nLndPts)     ,tw(1:nLndPts,jb),&
                     fvapa  (1:nLndPts),fsena (1:nLndPts),z0(1:nLndPts),ustar(1:nLndPts), hc(1:nLndPts),hg (1:nLndPts),ec (1:nLndPts),&
                     eg (1:nLndPts) ,d (1:nLndPts) ,cu (1:nLndPts),firb  (1:nLndPts),tgpptot (1:nLndPts),bstar1(1:nLndPts) ,&
                     ynleach_p (1:nLndPts,jb) ,tnmin_p   (1:nLndPts,jb) ,totnmic_p (1:nLndPts,jb) ,totnlit_p (1:nLndPts,jb) , &
                     totanlit_p(1:nLndPts,jb) ,totrnlit_p(1:nLndPts,jb),totnsoi_p (1:nLndPts,jb) ,storedn_p (1:nLndPts,jb) ,adnpp(1:nLndPts,1:npft,jb) ,&
                     adfalll     (1:nLndPts,jb)   ,adfallr     (1:nLndPts,jb)   ,adfallw    (1:nLndPts,jb),adcbiol(1:nLndPts,1:npft,jb),&
                     adcbior(1:nLndPts,1:npft,jb),adcbiow(1:nLndPts,1:npft,jb),adplai(1:nLndPts,1:npft,jb),beta1,beta2,&
                     stressfac,avmuir_factor, nVegClass ,rootmode)

                nLndPts=0
                DO i=1,nCols
                   IF (iMask(i) >= 1_i8) THEN
                      nLndPts=nLndPts+1
                      tu (nLndPts,jb) = tum     (nLndPts,jb  )
                      ts (nLndPts,jb) = tsm     (nLndPts,jb  )
                      tl (nLndPts,jb) = tlm     (nLndPts,jb  )
                      ti (nLndPts,jb) = tim     (nLndPts,jb  )
                      t12(nLndPts,jb) = ta     (nLndPts  )
                      t34(nLndPts,jb) = ta     (nLndPts  )        
                      tg(nLndPts,jb)  = tgm     (nLndPts,jb) 
                      q12(nLndPts,jb) = qa     (nLndPts  )        
                      q34(nLndPts,jb) = qa     (nLndPts  )        
                   END IF
                END DO

             END DO
             DO k=1,nsoilay
                DO i=1,nLndPts
!                   tsoim(i,k,jb) = tsoi (i,k,jb)
                   tsoi (i,k,jb) = tsoim(i,k,jb)
                END DO
             END DO
             DO k=1,nsoilay
                DO i=1,nLndPts
 !                  wsoim(i,k,jb) = wsoi(i,k,jb)
                   wsoi(i,k,jb)  = wsoim(i,k,jb)
                   wisoi(i,k,jb)  = wisoim(i,k,jb)

                END DO
             END DO
             DO k=1,nsnolay
                DO i=1,nLndPts
                   tsno     (i,k,jb) = tsnom   (i,k,jb)

                END DO
             END DO
             DO i=1,nLndPts
                !    capac(i,1)=capacm(i,1) wliqu
                !    capac(i,2)=capacm(i,2) wliqu
                !    tu   (i,jb)  =tum   (i,jb)
                IF(ind.EQ.1) THEN
                   tmin (i) =tg (i,jb)
                ELSE
                   tmax (i) =tg (i,jb)
                END IF
                tgm   (i,jb)=tg (i,jb)
                tg   (i,jb) =tgm(i,jb)
                wliqu(i,jb) = wliqum  (i,jb) 
                wliqs(i,jb) = wliqsm   (i,jb) 
                wliql(i,jb) = wliqlm   (i,jb) 
                wsnou(i,jb) = wsnoum   (i,jb)
                wsnos(i,jb) = wsnosm   (i,jb)  
                wsnol(i,jb) = wsnolm   (i,jb) 

             END DO
          END DO
          DO k=1,nsoilay
             DO i=1,nLndPts
                !          td   (i,k) =tdm   (i,k)
                tsoi (i,k,jb) = tsoim(i,k,jb)
                !    td   (i,k) =0.9_r8*0.5_r8*(tmax(i)+tmin(i))+0.1_r8*tdm(i,k)
                tsoi (i,k,jb) = 0.9_r8*0.5_r8*(tmax(i)+tmin(i))+0.1_r8*tsoim(i,k,jb)
                !    tdm  (i,k) =td(i,k)
                tsoim(i,k,jb) = tsoi (i,k,jb)
                !    td0  (i,k) =td(i,k)
                tsoi0(i,k,jb) = tsoi (i,k,jb)
             END DO
          END DO
          !
          !     this is a start of equilibrium tg,tc comp.
          !
          DO i=1,nLndPts
             IF(coszen(i).LT.0.0e0_r8) THEN
                tgm  (i,jb)  =tmin(i)
                tg0  (i,jb)  =tmin(i)
             END IF
          END DO
       END IF
    END IF
    !PPPPP
    IF(nLndPts.GE.1) THEN
       nLndPts=0
       DO i=1,nCols
          IF (iMask(i) >= 1_i8) THEN
             nLndPts=nLndPts+1
             !
             !     precipitation
             !
             raina  (nLndPts  ) = (ppli(i) + ppci(i) -snow   (i))*(1.0_r8/dtc3x) !1.0e-3_r8  !convert mm/s to m/s
          END IF
       END DO
       CALL IbisDrv (tod,pi,stef,vonk,grav,tmelt,hfus,hvap,hsub,ch2o,cice,cair,cvap,rair,rvap,cappa, &
            rhow,nLndPts,nband,nsoilay,nsnolay,npft,epsilon,dtc3x,doalb,&
            ginvap       (1:nLndPts,jb),gsuvap      (1:nLndPts,jb),gtrans         (1:nLndPts,jb), &
            gtransu      (1:nLndPts,jb),gtransl        (1:nLndPts,jb),grunof         (1:nLndPts,jb),gdrain           (1:nLndPts,jb), &
            gadjust      (1:nLndPts,jb),a10scalparamu(1:nLndPts,jb),a10daylightu(1:nLndPts,jb),a10scalparaml(1:nLndPts,jb), &
            a10daylightl (1:nLndPts,jb),vmax_pft        (1:npft)    ,tau15                     ,kc15                       , &
            ko15,cimax,gammaub,alpha3,theta3,beta3,coefmub,coefbub,gsubmin,gammauc,coefmuc,coefbuc, &
            gsucmin,gammals,coefmls,coefbls,gslsmin,gammal3,coefml3,coefbl3, &
            gsl3min,gammal4,alpha4,theta4,beta4,coefml4,coefbl4,gsl4min,bps(1:nLndPts), &
            wliqu        (1:nLndPts,jb),wliqumax                    ,wsnou         (1:nLndPts,jb),wsnoumax                 , &
            tu              (1:nLndPts,jb),wliqs        (1:nLndPts,jb),wliqsmax                     ,wsnos           (1:nLndPts,jb), &
            wsnosmax                  ,ts                (1:nLndPts,jb),wliql         (1:nLndPts,jb),wliqlmax                 , &  
            wsnol        (1:nLndPts,jb),wsnolmax                    ,tl          (1:nLndPts,jb),topparu           (1:nLndPts,jb), &
            topparl      (1:nLndPts,jb),fl                (1:nLndPts,jb),fu          (1:nLndPts,jb),lai      (1:nLndPts,1:2,jb), &
            sai          (1:nLndPts,1:2,jb),rhoveg      (1:nband,1:2),tauveg        (1:nband,1:2),orieh    (1:2)               , &
            oriev    (1:2)           ,wliqmin                    ,wsnomin                     ,t12           (1:nLndPts,jb), &
            tdripu                   ,tblowu                    ,tdrips                     ,tblows                       , &
            t34              (1:nLndPts,jb),tdripl                    ,tblowl,ztop     (1:nLndPts,1:2,jb),za     (1:nLndPts), & 
            alaiml                   ,zbot     (1:nLndPts,1:2,jb),alaimu                     ,froot        (1:nLndPts,1:nsoilay,1:2,jb), &
            q34              (1:nLndPts,jb),q12          (1:nLndPts,jb),su          (1:nLndPts,jb),cleaf                       , &
            dleaf        (1:2)          ,ss                (1:nLndPts,jb),cstem                     ,dstem                       , &
            sl              (1:nLndPts,jb),cgrass                    ,ciub         (1:nLndPts,jb),ciuc           (1:nLndPts,jb), &
            exist (1:nLndPts,1:npft,jb),csub         (1:nLndPts,jb),gsub         (1:nLndPts,jb),csuc           (1:nLndPts,jb), &
            gsuc              (1:nLndPts,jb),agcub        (1:nLndPts,jb),agcuc         (1:nLndPts,jb),ancub           (1:nLndPts,jb), &
            ancuc        (1:nLndPts,jb),totcondub        (1:nLndPts,jb),totconduc   (1:nLndPts,jb),cils           (1:nLndPts,jb), &
            cil3              (1:nLndPts,jb),cil4         (1:nLndPts,jb),csls         (1:nLndPts,jb),gsls           (1:nLndPts,jb), &
            csl3              (1:nLndPts,jb),gsl3         (1:nLndPts,jb),csl4         (1:nLndPts,jb),gsl4           (1:nLndPts,jb), &
            agcls        (1:nLndPts,jb),agcl4        (1:nLndPts,jb),agcl3         (1:nLndPts,jb),ancls           (1:nLndPts,jb), &
            ancl4        (1:nLndPts,jb),ancl3        (1:nLndPts,jb),totcondls   (1:nLndPts,jb),totcondl3    (1:nLndPts,jb), &
            totcondl4    (1:nLndPts,jb),chu(1:nVegClass),chs(1:nVegClass),chl(1:nVegClass),frac  (1:nLndPts,1:npft,jb),tlsub          (1:nLndPts,jb),z0sno,rhos, &
            consno                   ,hsnotop                    ,hsnomin                     ,fimin                       , &
            fimax                    ,fi                (1:nLndPts,jb),tsno(1:nLndPts,1:nsnolay,jb),hsno(1:nLndPts,1:nsnolay,jb),&
            sand(1:nLndPts,1:nsoilay,jb),clay(1:nLndPts,1:nsoilay,jb),poros(1:nLndPts,1:nsoilay,jb),wsoi(1:nLndPts,1:nsoilay,jb), &
            wisoi(1:nLndPts,1:nsoilay,jb),consoi(1:nLndPts,1:nsoilay,jb),zwpmax             ,wpud(1:nLndPts,jb)         , &
            wipud        (1:nLndPts,jb),wpudmax                    ,qglif   (1:nLndPts,1:4,jb),tsoi(1:nLndPts,1:nsoilay,jb), &
            hvasug       (1:nLndPts,jb),hvasui        (1:nLndPts,jb),albsav         (1:nLndPts,jb),albsan           (1:nLndPts,jb), &
            tg              (1:nLndPts,jb),ti                (1:nLndPts,jb),z0soi(1:nLndPts,jb)       ,swilt(1:nLndPts,1:nsoilay,jb), &
            sfield(1:nLndPts,1:nsoilay,jb),stressl(1:nLndPts,1:nsoilay,jb),stressu(1:nLndPts,1:nsoilay,jb),stresstl(1:nLndPts,jb), &
            stresstu (1:nLndPts,jb)    ,csoi(1:nLndPts,1:nsoilay,jb),rhosoi(1:nLndPts,1:nsoilay,jb),hsoi(1:nLndPts,1:nsoilay+1,jb)   , &
            suction(1:nLndPts,1:nsoilay,jb),bex (1:nLndPts,1:nsoilay,jb),upsoiu(1:nLndPts,1:nsoilay,jb),upsoil(1:nLndPts,1:nsoilay,jb), &
            heatg   (1:nLndPts,jb),heati   (1:nLndPts,jb),hydraul(1:nLndPts,1:nsoilay,jb),porosflo(1:nLndPts,1:nsoilay,jb), &
            ibex(1:nLndPts,1:nsoilay,jb),bperm(1:nLndPts,jb)          ,hflo      (1:nLndPts,1:nsoilay+1,jb),ta       (1:nLndPts)       , &
            asurd (1:nLndPts,1:nBand,jb),asuri(1:nLndPts,1:nBand,jb),coszen(1:nLndPts)              ,solad (1:nLndPts,1:nBand)  , &
            solai (1:nLndPts,1:nBand)  ,fira      (1:nLndPts)      ,raina (1:nLndPts)             ,qa    (1:nLndPts)               , &
            psurf    (1:nLndPts)          ,snowa     (1:nLndPts)      ,ua    (1:nLndPts)             ,o2conc                       , &
            co2conc  (1:nLndPts,jb),td         (1:nLndPts,jb)            ,vzero    (1:nLndPts,jb)   ,ndaypy                       , &
            nppdummy(1:nLndPts,1:npft,jb) ,cbiow (1:nLndPts,1:npft,jb),sapfrac(1:nLndPts,jb)      ,cbior(1:nLndPts,1:npft,jb), &
            tco2root(1:nLndPts,jb)          ,tneetot(1:nLndPts,jb)            ,tco2mic (1:nLndPts,jb)       ,a10td    (1:nLndPts,jb)    , &
            a10ancub(1:nLndPts,jb)          ,a10ancuc(1:nLndPts,jb)     ,a10ancls(1:nLndPts,jb)    ,a10ancl3(1:nLndPts,jb)     , &
            a10ancl4(1:nLndPts,jb)          ,ndtimes(1:nLndPts,jb), &
            adrain  (1:nLndPts,jb)          ,adsnow  (1:nLndPts,jb),tnpptot(1:nLndPts,jb),adaet   (1:nLndPts,jb)    ,adtrunoff(1:nLndPts,jb), &
            adsrunoff(1:nLndPts,jb)    ,addrainage(1:nLndPts,jb)   ,adrh    (1:nLndPts,jb)    ,adsnod   (1:nLndPts,jb)    , &
            adsnof   (1:nLndPts,jb)    ,adwsoi    (1:nLndPts,jb)   ,adtsoi  (1:nLndPts,jb)    ,adwisoi  (1:nLndPts,jb)    , &
            adtlaysoi(1:nLndPts,jb)    ,adwlaysoi (1:nLndPts,jb)   ,adwsoic (1:nLndPts,jb)    ,adtsoic  (1:nLndPts,jb)    , &
            adco2mic (1:nLndPts,jb)    ,adco2root (1:nLndPts,jb)   ,adco2soi(1:nLndPts,jb)    ,adco2ratio(1:nLndPts,jb)   , &
            adnmintot(1:nLndPts,jb)    ,decompl   (1:nLndPts,jb)   ,decomps (1:nLndPts,jb)    ,tnmin    (1:nLndPts,jb)    , &
            ndaypm                   ,nmtimes(1:nLndPts,jb),amrain   (1:nLndPts,jb)    , &
            amsnow    (1:nLndPts,jb)   ,amaet     (1:nLndPts,jb)   ,amtrunoff(1:nLndPts,jb)   ,amsrunoff(1:nLndPts,jb)    , &
            amdrainage(1:nLndPts,jb)   ,amtemp    (1:nLndPts,jb)   ,amqa     (1:nLndPts,jb)    , &
            amsolar   (1:nLndPts,jb)   ,amirup   (1:nLndPts,jb)   ,amirdown (1:nLndPts,jb)    , &
            amsens    (1:nLndPts,jb)   ,amlatent  (1:nLndPts,jb)   ,amlaiu   (1:nLndPts,jb)   ,amlail   (1:nLndPts,jb)    , &
            amtsoi    (1:nLndPts,jb)   ,amwsoi    (1:nLndPts,jb)   ,amwisoi  (1:nLndPts,jb)   ,amvwc    (1:nLndPts,jb)    , &
            amawc     (1:nLndPts,jb)   ,amsnod    (1:nLndPts,jb)   ,amsnof   (1:nLndPts,jb)   ,amnpp(1:nLndPts,1:npft,jb) , &
            amnpptot  (1:nLndPts,jb)   ,amco2mic  (1:nLndPts,jb)   ,amco2root(1:nLndPts,jb)   ,amco2soi (1:nLndPts,jb)    , &
            amco2ratio(1:nLndPts,jb)   ,amneetot  (1:nLndPts,jb)   ,amnmintot(1:nLndPts,jb)   ,amalbedo (1:nLndPts,jb)    , &
            amtsoil(1:nLndPts,1:nsoilay,jb) ,amwsoil(1:nLndPts,1:nsoilay,jb),amwisoil(1:nLndPts,1:nsoilay,jb), nytimes(1:nLndPts,jb), &
            aysolar   (1:nLndPts,jb)   ,ayirup    (1:nLndPts,jb)   ,ayirdown (1:nLndPts,jb)   ,aysens        (1:nLndPts,jb)   , &
            aylatent  (1:nLndPts,jb)   ,ayprcp    (1:nLndPts,jb)   ,ayaet    (1:nLndPts,jb)   ,aytrans        (1:nLndPts,jb)   , &
            aytrunoff (1:nLndPts,jb)   ,aysrunoff (1:nLndPts,jb)   ,aydrainage(1:nLndPts,jb)  ,aydwtot        (1:nLndPts,jb)   , & 
            aywsoi    (1:nLndPts,jb)   ,aywisoi   (1:nLndPts,jb)   ,aytsoi   (1:nLndPts,jb)   ,ayvwc        (1:nLndPts,jb)   , &
            ayawc     (1:nLndPts,jb)   ,aystresstu(1:nLndPts,jb)   ,aystresstl(1:nLndPts,jb)  ,aygpp(1:nLndPts,1:npft,jb) , &
            aygpptot  (1:nLndPts,jb)   ,aynpp(1:nLndPts,1:npft,jb) ,aynpptot (1:nLndPts,jb)   ,ayco2mic  (1:nLndPts,jb)   , &
            ayco2root (1:nLndPts,jb)   ,ayco2soi  (1:nLndPts,jb)   ,ayneetot (1:nLndPts,jb)   ,ayrootbio (1:nLndPts,jb)   , &
            aynmintot (1:nLndPts,jb)   ,ayalit    (1:nLndPts,jb)   ,ayblit   (1:nLndPts,jb)   ,aycsoi        (1:nLndPts,jb)   , & 
            aycmic    (1:nLndPts,jb)   ,ayanlit   (1:nLndPts,jb)   ,aybnlit  (1:nLndPts,jb)   ,aynsoi        (1:nLndPts,jb)   , &
            ayalbedo  (1:nLndPts,jb)   ,totalit   (1:nLndPts,jb)   ,totrlit  (1:nLndPts,jb)   ,totcsoi  (1:nLndPts,jb)    , &
            totcmic   (1:nLndPts,jb)   , &
            totanlit  (1:nLndPts,jb)   ,totrnlit  (1:nLndPts,jb)   ,totnsoi  (1:nLndPts,jb)   ,totnmic        (1:nLndPts,jb)   , &
            totlit    (1:nLndPts,jb)   ,totfall   (1:nLndPts,jb)   ,totnlit  (1:nLndPts,jb)   ,firefac        (1:nLndPts,jb)   , &
            wtot           (1:nLndPts,jb)   ,storedn   (1:nLndPts,jb)   ,yrleach  (1:nLndPts,jb)   ,ynleach        (1:nLndPts,jb)   , & 
            falll     (1:nLndPts,jb)   ,fallr     (1:nLndPts,jb)   ,fallw    (1:nLndPts,jb)   ,clitlm        (1:nLndPts,jb)   , &
            clitls    (1:nLndPts,jb)   ,clitrm    (1:nLndPts,jb)   ,clitrs   (1:nLndPts,jb)   ,clitwm        (1:nLndPts,jb)   , &
            clitws    (1:nLndPts,jb)   ,csoislop  (1:nLndPts,jb)   ,csoislon (1:nLndPts,jb)   ,csoipas        (1:nLndPts,jb)   , &
            clitll    (1:nLndPts,jb)   ,clitrl    (1:nLndPts,jb)   ,clitwl   (1:nLndPts,jb)   ,tc        (1:nLndPts,jb)   , &
            agddu     (1:nLndPts,jb)   ,tempu     (1:nLndPts,jb)   ,agddl    (1:nLndPts,jb)   ,templ        (1:nLndPts,jb)   , &
            dropu     (1:nLndPts,jb)   ,dropls    (1:nLndPts,jb)   ,dropl4   (1:nLndPts,jb)   ,dropl3        (1:nLndPts,jb)   , &
            plai(1:nLndPts,1:npft,jb)  ,iday,imonth,iyear,iyear0,isimveg,spinmax, &
            amts2     (1:nLndPts,jb)   ,amtransu  (1:nLndPts,jb)  ,amtransl (1:nLndPts,jb),   amsuvap  (1:nLndPts,jb)        ,&
            aminvap   (1:nLndPts,jb)   ,ux             (1:nLndPts)      ,uy       (1:nLndPts)      ,taux        (1:nLndPts)      , &
            tauy(1:nLndPts),ts2(1:nLndPts),qs2(1:nLndPts),deltat(1:nLndPts,jb),gdd0(1:nLndPts,jb), &  
            gdd0this  (1:nLndPts,jb)   ,tcthis    (1:nLndPts,jb)   ,twthis   (1:nLndPts,jb)   ,tcmin        (1:nLndPts,jb)   , &
            gdd5           (1:nLndPts,jb)   ,gdd5this  (1:nLndPts,jb)   ,TminL    (1:npft)       ,TminU        (1:npft)       , &
            Twarm     (1:npft)          ,GDD       (1:npft)            ,aleaf    (1:npft)       ,awood        (1:npft)       , &
            cbiol(1:nLndPts,1:npft,jb) ,aroot     (1:npft)            ,disturbf (1:nLndPts)      ,disturbo  (1:nLndPts)      , &
            specla    (1:npft)          ,biomass(1:nLndPts,1:npft,jb),totlaiu (1:nLndPts,jb)   ,totlail        (1:nLndPts,jb)   , &  
            totbiou   (1:nLndPts,jb)   ,totbiol   (1:nLndPts,jb)   ,woodnorm                     ,vegtype0  (1:nLndPts,jb)   , &
            tauwood0  (1:npft)          ,tauwood   (1:nLndPts,1:npft,jb)            ,tauleaf  (1:npft)       ,tauroot        (1:npft)       , &
            xminlai                  ,cdisturb  (1:nLndPts,jb)   ,ayanpp(1:nLndPts,1:npft,jb),ayanpptot(1:nLndPts,jb)   , &
            jdt,imonthprev,isimco2,isimfire,&
            co2initm (1:nLndPts)     ,tw(1:nLndPts,jb),fvapa  (1:nLndPts),fsena (1:nLndPts),z0(1:nLndPts),ustar(1:nLndPts) ,hc(1:nLndPts),&
            hg (1:nLndPts),ec (1:nLndPts),eg (1:nLndPts),d(1:nLndPts),cu(1:nLndPts),firb  (1:nLndPts),tgpptot (1:nLndPts),bstar1(1:nLndPts),&
            ynleach_p (1:nLndPts,jb) ,tnmin_p   (1:nLndPts,jb) ,totnmic_p (1:nLndPts,jb) ,totnlit_p (1:nLndPts,jb) ,&
            totanlit_p(1:nLndPts,jb) ,totrnlit_p(1:nLndPts,jb),totnsoi_p (1:nLndPts,jb) ,storedn_p (1:nLndPts,jb) , adnpp(1:nLndPts,1:npft,jb),&
            adfalll     (1:nLndPts,jb)   ,adfallr     (1:nLndPts,jb)   ,adfallw    (1:nLndPts,jb),adcbiol(1:nLndPts,1:npft,jb),&
            adcbior(1:nLndPts,1:npft,jb),adcbiow(1:nLndPts,1:npft,jb),adplai(1:nLndPts,1:npft,jb),beta1,beta2 ,&
            stressfac, avmuir_factor,nVegClass ,rootmode)

    END IF

    nLndPts=0
    DO i=1,nCols
       IF (iMask(i) >= 1_i8) THEN
          nLndPts=nLndPts+1
          !tgeff(i)=SQRT ( SQRT (( firb(i)*stef )))fu fi fl        tgeff(i)=SQRT ( SQRT (( firb(i)*stef )))
          tsea(i)            =SQRT ( SQRT (( MIN(MAX(firb(nLndPts),0.0_r8),1000.0_r8)*(1.0_r8  /stef) )))!(tu(nLndPts,jb) + tl(nLndPts,jb)+ti(nLndPts,jb) + tg(nLndPts,jb))/4.0_r8

          tsea(i)= MAX(218.00_r8,tsea(i))

          TSK (i)=SQRT ( SQRT (( MIN(MAX(firb(nLndPts),0.0_r8),1000.0_r8)*(1.0_r8  /stef) )))!(tu(nLndPts,jb) + tl(nLndPts,jb)+ti(nLndPts,jb) + tg(nLndPts,jb))/4.0_r8
       END IF
    END DO
    !     
    !     temperature and snow depths in Antarctica and Groenland
    !     
!    CALL sfc_sice1(nLndPts,nsoilay,dtc3x   ,&
!       psurf         (1:nLndPts  )   ,&
!       ux            (1:nLndPts  )   ,&
!       uy            (1:nLndPts  )   ,&
!       ta            (1:nLndPts  )   ,&
!       qa            (1:nLndPts  )   ,&
!       vegtype0      (1:nLndPts,jb)  ,&
!       dlwflx        (1:nLndPts  )   ,&
!       sfcnsw        (1:nLndPts  )   ,&
!       sfcdsw        (1:nLndPts  )   ,&
!       firb          (1:nLndPts  )   ,&
!       tprcp_landice (1:nLndPts  )   ,&
!       zorl_landice  (1:nLndPts  )   ,&
!       tsoi_landice  (1:nLndPts,1:nsoilay),&
!       sheleg_landice(1:nLndPts  )  ,&
!       tsea_landice  (1:nLndPts ),&
!       qsurf_landice (1:nLndPts ),&
!       bstar_landice (1:nLndPts ),&
!       evap_landice  (1:nLndPts ),&
!       hflx_landice  (1:nLndPts ),&
!       taux_landice  (1:nLndPts ),&
!       tauy_landice  (1:nLndPts ),&
!       ustar_landice (1:nLndPts) )

    DO k = 1, nsoilay
       nLndPts=0
       DO i=1,nCols
          IF (iMask(i) >= 1_i8)THEN 
             nLndPts=nLndPts+1
             IF (INT(vegtype0  (nLndPts,jb))  == 15_i8) THEN
                !zorl (i) = zorl_landice (nLndPts  ) 
                !tsea (i) = tsea_landice  (nLndPts )
                !taux  (nLndPts )=taux_landice  (nLndPts )
                !tauy  (nLndPts )=tauy_landice  (nLndPts )
                !ustar (nLndPts )=ustar_landice (nLndPts ) 
                !bstar (nLndPts )=bstar_landice (nLndPts )
                !tsoi  (nLndPts,k,jb)=tsoi_landice  (nLndPts,k)
                !fvapa (nLndPts) = -evap_landice  (nLndPts )
                !fsena (nLndPts) = -hflx_landice  (nLndPts )
                !PRINT*,evap_landice  (nLndPts )*hvap,hflx_landice  (nLndPts ),-fsena (nLndPts),-fvapa (nLndPts)*hvap
                !
                !                 Antarctica ice
                !
                tsoi (nLndPts,k,jb) = MIN(MAX(218.00_r8,tsoi  (nLndPts,k,jb)),273.16_r8)
                wisoi(nLndPts,k,jb) = MIN(MAX(0.94_r8   ,wisoi (nLndPts,k,jb)),1.00_r8  )
                wsoi (nLndPts,k,jb) = MIN(MAX(0.4_r8   ,wsoi  (nLndPts,k,jb)),1.0_r8  )
             END IF
          END IF
       END DO
    END DO
    DO k =1, nsnolay
       nLndPts=0
       DO i=1,nCols
          IF (iMask(i) >= 1_i8)THEN 
             nLndPts=nLndPts+1
             IF (INT(vegtype0  (nLndPts,jb))  == 15_i8) THEN
                !
                !                 Antarctica snow
                !
                fi   (nLndPts  ,jb) = 1.0_r8
                tg   (nLndPts  ,jb) = MIN(MAX(218.00_r8,tg   (nLndPts  ,jb)),273.16_r8)
                tsno (nLndPts,k,jb) = MIN(MAX(218.00_r8,tsno (nLndPts,k,jb)),273.16_r8)
                hsno (nLndPts,k,jb) = 2.0_r8
                IF(k==1)hsno(nLndPts,1,jb) = 0.050_r8
                IF(k==2)hsno(nLndPts,2,jb) = 2.0_r8
                IF(k==3)hsno(nLndPts,3,jb) = 2.0_r8
             END IF
          END IF
       END DO
    END DO

    !
    !     sib time integaration and time filter 
    !
    DO i=1,nLndPts
       !      qm(i)=MAX(1.0e-12_r8,qm(i))
    END DO
    CALL sextrp( &
       wsoi   (1:nLndPts,1:nsoilay,jb), &! ,    INTENT(in	)  Var_StepP    , &
       wsoi0  (1:nLndPts,1:nsoilay,jb), &! ,    INTENT(inout) Var_Step0    , &
       wsoim  (1:nLndPts,1:nsoilay,jb), &! ,    INTENT(inout) Var_StepM    , &
       nLndPts         , &     !, INTENT(in   ) :: nmax
       nsoilay         , &     !, INTENT(in   ) :: zmax
       istrt           , &     !, INTENT(in   ) :: istrt
       epsflt          , &     !, INTENT(in   ) :: epsflt
       intg              )     !, INTENT(in   ) :: intg
    CALL sextrp( &
       wisoi   (1:nLndPts,1:nsoilay,jb), &! ,    INTENT(in	)  Var_StepP    , &
       wisoi0  (1:nLndPts,1:nsoilay,jb), &! ,    INTENT(inout) Var_Step0    , &
       wisoim  (1:nLndPts,1:nsoilay,jb), &! ,    INTENT(inout) Var_StepM    , &
       nLndPts         , &     !, INTENT(in   ) :: nmax
       nsoilay         , &     !, INTENT(in   ) :: zmax
       istrt           , &     !, INTENT(in   ) :: istrt
       epsflt          , &     !, INTENT(in   ) :: epsflt
       intg              )     !, INTENT(in   ) :: intg

    CALL sextrp( &
       tsoi   (1:nLndPts,1:nsoilay,jb), &! ,    INTENT(in	)  Var_StepP    , &
       tsoi0  (1:nLndPts,1:nsoilay,jb), &! ,    INTENT(inout) Var_Step0    , &
       tsoim  (1:nLndPts,1:nsoilay,jb), &! ,    INTENT(inout) Var_StepM    , &
       nLndPts         , &     !, INTENT(in   ) :: nmax
       nsoilay         , &     !, INTENT(in   ) :: zmax
       istrt           , &     !, INTENT(in   ) :: istrt
       epsflt          , &     !, INTENT(in   ) :: epsflt
       intg              )     !, INTENT(in   ) :: intg

    CALL sextrp( &
       tg   (1:nLndPts,jb), &! ,    INTENT(in	)  Var_StepP    , &
       tg0  (1:nLndPts,jb), &! ,    INTENT(inout) Var_Step0    , &
       tgm  (1:nLndPts,jb), &! ,    INTENT(inout) Var_StepM    , &
       nLndPts         , &     !, INTENT(in   ) :: nmax
       1               , &     !, INTENT(in   ) :: zmax
       istrt           , &     !, INTENT(in   ) :: istrt
       epsflt          , &     !, INTENT(in   ) :: epsflt
       intg              )     !, INTENT(in   ) :: intg

    CALL sextrp( &
       tsno   (1:nLndPts,1:nsnolay,jb), &! ,    INTENT(in	)  Var_StepP    , &
       tsno0  (1:nLndPts,1:nsnolay,jb), &! ,    INTENT(inout) Var_Step0    , &
       tsnom  (1:nLndPts,1:nsnolay,jb), &! ,    INTENT(inout) Var_StepM    , &
       nLndPts         , &     !, INTENT(in   ) :: nmax
       nsnolay            , &     !, INTENT(in   ) :: zmax
       istrt           , &     !, INTENT(in   ) :: istrt
       epsflt          , &     !, INTENT(in   ) :: epsflt
       intg              )     !, INTENT(in   ) :: intg


    CALL sextrp( &
       wliqu   (1:nLndPts,jb), &! ,    INTENT(in	)  Var_StepP    , &
       wliqu0  (1:nLndPts,jb), &! ,    INTENT(inout) Var_Step0    , &
       wliqum  (1:nLndPts,jb), &! ,    INTENT(inout) Var_StepM    , &
       nLndPts         , &     !, INTENT(in   ) :: nmax
       1               , &     !, INTENT(in   ) :: zmax
       istrt           , &     !, INTENT(in   ) :: istrt
       epsflt          , &     !, INTENT(in   ) :: epsflt
       intg              )     !, INTENT(in   ) :: intg

    CALL sextrp( &
       wliqs   (1:nLndPts,jb), &! ,    INTENT(in	)  Var_StepP    , &
       wliqs0  (1:nLndPts,jb), &! ,    INTENT(inout) Var_Step0    , &
       wliqsm  (1:nLndPts,jb), &! ,    INTENT(inout) Var_StepM    , &
       nLndPts         , &     !, INTENT(in   ) :: nmax
       1               , &     !, INTENT(in   ) :: zmax
       istrt           , &     !, INTENT(in   ) :: istrt
       epsflt          , &     !, INTENT(in   ) :: epsflt
       intg              )     !, INTENT(in   ) :: intg

    CALL sextrp( &
       wliql   (1:nLndPts,jb), &! ,    INTENT(in	)  Var_StepP    , &
       wliql0  (1:nLndPts,jb), &! ,    INTENT(inout) Var_Step0    , &
       wliqlm  (1:nLndPts,jb), &! ,    INTENT(inout) Var_StepM    , &
       nLndPts         , &     !, INTENT(in   ) :: nmax
       1               , &     !, INTENT(in   ) :: zmax
       istrt           , &     !, INTENT(in   ) :: istrt
       epsflt          , &     !, INTENT(in   ) :: epsflt
       intg              )     !, INTENT(in   ) :: intg

    CALL sextrp( &
       wsnou   (1:nLndPts,jb), &! ,    INTENT(in	)  Var_StepP    , &
       wsnou0  (1:nLndPts,jb), &! ,    INTENT(inout) Var_Step0    , &
       wsnoum  (1:nLndPts,jb), &! ,    INTENT(inout) Var_StepM    , &
       nLndPts         , &     !, INTENT(in   ) :: nmax
       1               , &     !, INTENT(in   ) :: zmax
       istrt           , &     !, INTENT(in   ) :: istrt
       epsflt          , &     !, INTENT(in   ) :: epsflt
       intg              )     !, INTENT(in   ) :: intg

    CALL sextrp( &
       wsnos   (1:nLndPts,jb), &! ,    INTENT(in	)  Var_StepP    , &
       wsnos0  (1:nLndPts,jb), &! ,    INTENT(inout) Var_Step0    , &
       wsnosm  (1:nLndPts,jb), &! ,    INTENT(inout) Var_StepM    , &
       nLndPts         , &     !, INTENT(in   ) :: nmax
       1               , &     !, INTENT(in   ) :: zmax
       istrt           , &     !, INTENT(in   ) :: istrt
       epsflt          , &     !, INTENT(in   ) :: epsflt
       intg              )     !, INTENT(in   ) :: intg

    CALL sextrp( &
       wsnol   (1:nLndPts,jb), &! ,    INTENT(in	)  Var_StepP    , &
       wsnol0  (1:nLndPts,jb), &! ,    INTENT(inout) Var_Step0    , &
       wsnolm  (1:nLndPts,jb), &! ,    INTENT(inout) Var_StepM    , &
       nLndPts         , &     !, INTENT(in   ) :: nmax
       1               , &     !, INTENT(in   ) :: zmax
       istrt           , &     !, INTENT(in   ) :: istrt
       epsflt          , &     !, INTENT(in   ) :: epsflt
       intg              )     !, INTENT(in   ) :: intg

    CALL sextrp( &
       tu   (1:nLndPts,jb), &! ,    INTENT(in	)  Var_StepP    , &
       tu0  (1:nLndPts,jb), &! ,    INTENT(inout) Var_Step0    , &
       tum  (1:nLndPts,jb), &! ,    INTENT(inout) Var_StepM    , &
       nLndPts         , &     !, INTENT(in   ) :: nmax
       1               , &     !, INTENT(in   ) :: zmax
       istrt           , &     !, INTENT(in   ) :: istrt
       epsflt          , &     !, INTENT(in   ) :: epsflt
       intg              )     !, INTENT(in   ) :: intg

    CALL sextrp( &
       ts   (1:nLndPts,jb), &! ,    INTENT(in	)  Var_StepP    , &
       ts0  (1:nLndPts,jb), &! ,    INTENT(inout) Var_Step0    , &
       tsm  (1:nLndPts,jb), &! ,    INTENT(inout) Var_StepM    , &
       nLndPts         , &     !, INTENT(in   ) :: nmax
       1               , &     !, INTENT(in   ) :: zmax
       istrt           , &     !, INTENT(in   ) :: istrt
       epsflt          , &     !, INTENT(in   ) :: epsflt
       intg              )     !, INTENT(in   ) :: intg

    CALL sextrp( &
       tl   (1:nLndPts,jb), &! ,    INTENT(in	)  Var_StepP    , &
       tl0  (1:nLndPts,jb), &! ,    INTENT(inout) Var_Step0    , &
       tlm  (1:nLndPts,jb), &! ,    INTENT(inout) Var_StepM    , &
       nLndPts         , &     !, INTENT(in   ) :: nmax
       1               , &     !, INTENT(in   ) :: zmax
       istrt           , &     !, INTENT(in   ) :: istrt
       epsflt          , &     !, INTENT(in   ) :: epsflt
       intg              )     !, INTENT(in   ) :: intg

    CALL sextrp( &
       ti   (1:nLndPts,jb), &! ,    INTENT(in	)  Var_StepP    , &
       ti0  (1:nLndPts,jb), &! ,    INTENT(inout) Var_Step0    , &
       tim  (1:nLndPts,jb), &! ,    INTENT(inout) Var_StepM    , &
       nLndPts         , &     !, INTENT(in   ) :: nmax
       1               , &     !, INTENT(in   ) :: zmax
       istrt           , &     !, INTENT(in   ) :: istrt
       epsflt          , &     !, INTENT(in   ) :: epsflt
       intg              )     !, INTENT(in   ) :: intg

    nLndPts=0
    DO i=1,nCols
       IF (iMask(i) >= 1_i8)THEN 
          nLndPts=nLndPts+1
          sens    (i) = MIN(MAX(-fsena (nLndPts),MinFlux),MaxFlux)
          evap    (i) = MIN(MAX(-fvapa (nLndPts)*hvap,MinFlux),MaxFlux)

          tgm  (nLndPts,jb  ) = tgm     (nLndPts,jb )
          tg   (nLndPts,jb  ) = tg      (nLndPts,jb )
          tg0  (nLndPts,jb  ) = tg0     (nLndPts,jb )

          DO k=1,nsoilay 
             tsoi0(nLndPts,k,jb) = tsoi0 (nLndPts,k,jb)
             tsoim(nLndPts,k,jb) = tsoim (nLndPts,k,jb)
             tsoi (nLndPts,k,jb) = tsoi  (nLndPts,k,jb)
          END DO
          DO k=1,nsoilay
             wsoi0(nLndPts,k,jb)  = wsoi0 (nLndPts,k,jb)
             wsoim(nLndPts,k,jb)  = wsoim (nLndPts,k,jb)
             wsoi (nLndPts,k,jb)  = wsoi  (nLndPts,k,jb)
          END DO
       END IF
    END DO
    nLndPts=0
    DO i=1,nCols
       IF(iMask(i) >= 1_i8 ) THEN
          nLndPts=nLndPts+1
          IF ( iMask(i) == 15_i8 ) THEN
             sm0(nLndPts,1)    = poros(nLndPts,1,jb)
          ELSE
             sm0(nLndPts,1)      =wsoi(nLndPts,1,jb)*poros(nLndPts,1,jb)
             sm0(nLndPts,2)      =wsoi(nLndPts,nsoilay/2,jb)*poros(nLndPts,nsoilay/2,jb)
             sm0(nLndPts,3)      =wsoi(nLndPts,nsoilay,jb)*poros(nLndPts,nsoilay,jb)
          END IF
       END IF
    END DO

    !
    !     sea or sea ice
    ! gu gv gps colrad sens evap umom vmom rmi rhi cond stor zorl rnet ztn2 THETA_2M VELC_2m MIXQ_2M
    ! THETA_10M VELC_10M MIXQ_10M
    ! mmax=ncols-nmax+1
    ! including case 1D physics

    IF(initlz >= 0 .AND. kt == 0 .AND. jdt == 1)THEN
       Tsfc0(1:nCols,jb)=gt  (1:nCols,1)
       Qsfc0(1:nCols,jb)=gq  (1:nCols,1)
       Tsfcm(1:nCols,jb)=gt  (1:nCols,1)
       Qsfcm(1:nCols,jb)=gq  (1:nCols,1)
       tsfc (1:nCols)=gt  (1:nCols,1)
       qsfc (1:nCols)=gq  (1:nCols,1)
    END IF
    DO i=1,nCols
       IF(mskant(i) == 1_i8)THEN
          xsea (i) = tseam(i)
          tsfc (i) = Tsfcm(i,jb)
          qsfc (i) = Qsfcm(i,jb)
       END IF
    END DO

   CALL seasfc( &
           tmtx  (1:nCols,1:kMax,1:3)  ,umtx  (1:nCols,1:kMax,1:4),qmtx  (1:nCols,1:kMax,1:3)  ,&
           kmax                        ,kmax                      ,slrad (1:nCols)             ,&
           tsurf (1:nCols)             ,qsurf (1:nCols)           ,gu    (1:nCols,1:kMax)      ,&
           gv    (1:nCols,1:kMax)      ,gt    (1:nCols,1:kMax)    ,gq    (1:nCols,1:kMax)      ,&
           prsi(1:ncols,1:kMax+1),prsl(1:ncols,1:kMax),phii(1:nCols,1:kMax+1),phil(1:nCols,1:kMax),&
           xsea  (1:nCols)           ,dtc3x                       ,SIN(colrad(1:nCols))        ,&
           sens  (1:nCols)             ,evap  (1:nCols)           ,umom  (1:nCols)             ,&
           vmom  (1:nCols)             ,rmi   (1:nCols)           ,rhi   (1:nCols)             ,&
           cond  (1:nCols)             ,stor  (1:nCols)           ,zorl  (1:nCols)             ,&
           nCols                       ,speedm(1:nCols)           ,bstar (1:nCols)             ,&
           Ustarm(1:nCols)             ,z0sea (1:nCols)           ,rho   (1:nCols)             ,&
           qsfc  (1:nCols)             ,tsfc  (1:nCols)           ,MskAnt(1:nCols)             ,&
           iMask (1:nCols)             ,zenith (1:nCols)          ,ppli  (1:nCols)             ,&
           ppci  (1:nCols)             ,LwSfcDown(1:nCols)        ,xvisb (1:nCols)             ,&
           xvisd (1:nCols)             ,xnirb(1:nCols)            ,xnird (1:nCols)             ,&
           HML   (1:nCols)             ,HUML (1:nCols)            ,HVML (1:nCols)              ,&
           TSK   (1:nCols)             ,GSW(1:nCols)              ,GLW(1:nCols)                ,&
           cldtot(1:nCols,1:kMax)      ,sfcnsw(1:nCols)           ,month(1:nCols)             ,&
           sfcnlw(1:nCols)             ,pblh  (1:nCols)           ,QCF (1:nCols,1:kMax)        ,&
           QCL  (1:nCols,1:kMax)       ,mlsi  (1:nCols)           ,jb                          ,&
           Mmlen (1:nCols)             ,colrad(1:nCols)           ,idatec,dump(1:nCols,1:kMax ))


    DO i=1,nCols
       IF(mskant(i) == 1_i8 .AND. tsea(i) <= 0.0e0_r8 .AND. tsurf(i) < tice+0.01e0_r8 ) THEN
          IF(intg.EQ.2) THEN
             IF(istrt.EQ.0) THEN
                tseam(i)=filta*tsea (i) + epsflt*(tseam(i)+xsea(i))
                qsfc (i)=MAX(1.0e-12_r8,qsfc(i))
                Tsfcm(i,jb)=filta*Tsfc0 (i,jb) + epsflt*(Tsfcm(i,jb)+tsfc(i))
                Qsfcm(i,jb)=filta*Qsfc0 (i,jb) + epsflt*(Qsfcm(i,jb)+qsfc(i))
             END IF
             tsea (i) = xsea(i)
             qsfc (i) = MAX(1.0e-12_r8,qsfc(i))
             Tsfc0(i,jb) = tsfc(i)
             Qsfc0(i,jb) = qsfc(i)
          ELSE
             tsea (i) = xsea(i)
             tseam(i) = xsea(i)
             qsfc (i) = MAX(1.0e-12_r8,qsfc(i))
             Tsfc0(i,jb) = tsfc(i)
             Qsfc0(i,jb) = qsfc(i)
             Tsfcm(i,jb) = tsfc(i)
             Qsfcm(i,jb) = qsfc(i)
          END IF
       END IF
       IF(mskant(i) == 1_i8 .AND. tsea(i).LT.0.0e0_r8.AND.tsurf(i).GE.tice+0.01e0_r8) THEN
          tseam(i) = tsea (i)
          Tsfcm(i,jb) = Tsfc0(i,jb)
          Qsfcm(i,jb) = Qsfc0(i,jb)
       END IF
    END DO

    nLndPts=0
    DO i=1,nCols
       IF(iMask(i) >= 1_i8)THEN 
          nLndPts=nLndPts+1
          tgrd(nLndPts)= tg    (nLndPts,jb)
          roff(nLndPts)= grunof(nLndPts,jb)
          ect(nLndPts )= gtrans(nLndPts,jb) * hltm * dtc3x!Transpiracao no topo da copa (J/m*m)
          eci(nLndPts )= ginvap(nLndPts,jb) * hltm * dtc3x!...Evaporacao da agua interceptada no topo da copa (J/m*m) ! (kg m-2 s-1 * J/kg*dt)
          egt(nLndPts )= gtransu(nLndPts,jb)* hltm * dtc3x!Transpiracao na base da copa (J/m*m)
          egi(nLndPts )= gtransl(nLndPts,jb)* hltm * dtc3x !Evaporacao da neve (J/m*m)
          egs(nLndPts )= gsuvap (nLndPts,jb)* hltm * dtc3x !Evaporacao do solo arido (J/m*m)

          !WRITE(*,'(A,5F12.5)')'pkubota', ect(nLndPts ),eci(nLndPts ),egt(nLndPts ),egi(nLndPts ),egs(nLndPts )

          iMaskIBIS(nLndPts,jb) =  INT(vegtype0  (nLndPts,jb)) 
          umom     (i) =  taux  (nLndPts)
          vmom     (i) =  tauy  (nLndPts)
          bstar    (i) =  bstar1(nLndPts)
          sens     (i) =  MIN(MAX(-fsena (nLndPts),MinFlux),MaxFlux)
          evap     (i) =  MIN(MAX(-fvapa (nLndPts)*hvap,MinFlux),MaxFlux)
!          Tsfc0    (i,jb)       = tg    (nLndPts,jb)
          Tsfc0    (i,jb)       = SQRT ( SQRT (( MIN(MAX(firb(nLndPts),0.0_r8),1000.0_r8)*(1.0_r8  /stef) )))
          ! Tsfc0    (i,jb)       = ts2    (nLndPts)
          ! Qsfc0    (i,jb)       = MAX   (1.0e-12_r8,qa (nLndPts))
          Qsfc0    (i,jb)       = MAX   (1.0e-12_r8,qs2 (nLndPts))
          Tsfcm    (i,jb)       = tg    (nLndPts,jb)
          !Tsfcm    (i,jb)       = SQRT ( SQRT (( MIN(MAX(firb(nLndPts),0.0_r8),1000.0_r8)*(1.0_r8  /stef) )))
          ! Tsfcm    (i,jb)       = ts2    (nLndPts)
          ! Qsfcm    (i,jb)       = MAX   (1.0e-12_r8,qa (nLndPts))       
          Qsfcm    (i,jb)       = MAX   (1.0e-12_r8,qs2 (nLndPts))       
          w0       (nLndPts,1,jb) = wsoi0(nLndPts,1,jb)
          w0       (nLndPts,2,jb) = wsoi0(nLndPts,nsoilay/2,jb)
          w0       (nLndPts,3,jb) = wsoi0(nLndPts,nsoilay  ,jb)
          wm       (nLndPts,1,jb) = wsoim(nLndPts,1,jb)
          wm       (nLndPts,2,jb) = wsoim(nLndPts,nsoilay/2,jb)
          wm       (nLndPts,3,jb) = wsoim(nLndPts,nsoilay,jb)
          capac0   (nLndPts,1,jb) = wliqu(nLndPts,jb) + wliqs(nLndPts,jb)
          capac0   (nLndPts,2,jb) = wliql(nLndPts,jb)
          capacm   (nLndPts,1,jb) = wliqu(nLndPts,jb) + wliqs(nLndPts,jb)
          capacm   (nLndPts,2,jb) = wliql(nLndPts,jb)
          td0      (nLndPts,jb)   = tsoi0(nLndPts,nsoilay,jb)  
          tdm      (nLndPts,jb)   = tsoim(nLndPts,nsoilay,jb)  
          tc0      (nLndPts,jb)   = tu0    (nLndPts,jb)
          tcm      (nLndPts,jb)   = tum    (nLndPts,jb)
          tg0      (nLndPts,jb)   = tg0    (nLndPts,jb)
          tgm      (nLndPts,jb)   = tgm    (nLndPts,jb)
          tm0      (nLndPts,jb)   = ts2   (nLndPts)
          qm0      (nLndPts,jb)   = qs2   (nLndPts)
       END IF
    END DO
    cflx=0.0_r8
    kk=0
    DO itrac=nClass+1,nClass+nAeros
       IF(TRIM(typechem(itrac))=='CO2')THEN
          kk=kk+1
          nLndPts=0
          DO i=1,nCols
             IF(isimco2 == 0)THEN
                IF(iMask(i) >= 1_i8)THEN 
                   nLndPts=nLndPts+1
                   ! local ! instantaneous net ecosystem exchange of co2 per timestep (kg_C m-2/timestep)
                   ! (mol-CO2 / m-2 / second) to (kg-CO2 / m-2 / second)
                   cflx  (i,indexchem(itrac)) = 0.0_r8
                   !cflx  (i,indexchem(itrac)) = (tgpptot (nLndPts)/0.044_r8)
                ELSE
                   cflx  (i,indexchem(itrac)) = 0.0_r8
                END IF
             ELSE
                IF(iMask(i) >= 1_i8)THEN 
                   nLndPts=nLndPts+1
                   ! local ! instantaneous net ecosystem exchange of co2 per timestep (mol-CO2 / m-2 / second)
                   ! (mol-CO2 / m-2 / second) to (kg-CO2 / m-2 / second)
                   cflx  (i,indexchem(itrac)) = -(tneetot (nLndPts,jb)*(0.044_r8)) -gbaco2
                   co2diag(i) =cflxm  (i,indexchem(itrac),jb) - totflux(kk)
                   !cflx  (i,itrac) = (tgpptot (nLndPts)/0.044_r8)
                ELSE
                   cflx  (i,indexchem(itrac)) = gco2flx(i,jb)- gbaco2!*(0.044_r8)
                   co2diag(i) =cflxm  (i,indexchem(itrac),jb)  - totflux(kk)
               END IF
             END IF
          END DO
       END IF
    END DO

    !    asurd (1:nLndPts,1:nBand,jb),asuri(1:nLndPts,1:nBand,jb)
    IF(dodia(nDiag_biomau))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN 
             nLndPts=nLndPts+1
             diag(i) = totbiou (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_biomau,jb)
    ENDIF
    IF(dodia(nDiag_biomal))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN 
             nLndPts=nLndPts+1
             diag(i) = totbiol (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_biomal,jb)
    ENDIF
    IF(dodia(nDiag_tlaiup))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN 
             nLndPts=nLndPts+1
             diag(i) = totlaiu(nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_tlaiup,jb)
    ENDIF
    IF(dodia(nDiag_tlailw))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = totlail(nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_tlailw,jb)
    ENDIF
    IF(dodia(nDiag_tstnsp))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = storedn(nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_tstnsp,jb)
    ENDIF



    IF(dodia(nDiag_somdfa))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = decompl   (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_somdfa,jb)
    ENDIF
    IF(dodia(nDiag_lidecf))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = decomps (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_lidecf,jb)
    ENDIF
    IF(dodia(nDiag_wsttot))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = wtot (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_wsttot,jb)
    ENDIF

    IF(dodia(nDiag_facuca))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = fu (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_facuca,jb)
    ENDIF
    IF(dodia(nDiag_fsfclc))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = fl (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_fsfclc,jb)
    ENDIF
    IF(dodia(nDiag_frsnow))THEN !! fractional snow cover
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = fi (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_frsnow,jb)
    ENDIF

    IF(dodia(nDiag_insnpp))THEN ! instantaneous npp (mol-CO2 / m-2 / second)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = tnpptot (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_insnpp,jb)
    ENDIF

    IF(dodia(nDiag_insnee))THEN ! instantaneous net ecosystem exchange of co2 per timestep (kg_C m-2/timestep)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = co2diag(i)
          ELSE
             diag(i) = co2diag(i)
          END IF
       END DO
       CALL updia(diag,nDiag_insnee,jb)
    ENDIF

    IF(dodia(nDiag_grbdy0))THEN ! annual total growing degree days for current year > 0C
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = gdd0this (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_grbdy0,jb)
    ENDIF

    IF(dodia(nDiag_grbdy5))THEN ! annual total growing degree days for current year > 5C  
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = gdd5this (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_grbdy5,jb)
    ENDIF

    IF(dodia(nDiag_avet2m))THEN ! monthly average 2-m surface-air temperature 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = amts2 (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_avet2m,jb)
    ENDIF

    IF(dodia(nDiag_monnpp))THEN ! monthly total npp for ecosystem (kg-C/m**2/month)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = amnpptot (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_monnpp,jb)
    ENDIF

    IF(dodia(nDiag_monnee))THEN ! monthly total net ecosystem exchange of CO2 (kg-C/m**2/month)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = amneetot (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_monnee,jb)
    ENDIF

    IF(dodia(nDiag_yeanpp))THEN ! annual total npp for ecosystem (kg-c/m**2/yr)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = aynpptot (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_yeanpp,jb)
    ENDIF

    IF(dodia(nDiag_yeanee))THEN ! annual total npp for ecosystem (kg-c/m**2/yr)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = ayneetot (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_yeanee,jb)
    ENDIF

    IF(dodia(nDiag_upclai))THEN ! upper canopy single-sided leaf area index (area leaf/area veg)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = lai (nLndPts,2,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_upclai,jb)
    ENDIF

    IF(dodia(nDiag_lwclai))THEN ! lower canopy single-sided leaf area index (area leaf/area veg)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = lai (nLndPts,1,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_lwclai,jb)
    ENDIF

    IF(dodia(nDiag_pfts01))THEN  ! pft tropical broadleaf evergreen trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = plai (nLndPts,1,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_pfts01,jb)
    ENDIF

    IF(dodia(nDiag_pfts02))THEN  ! pft tropical broadleaf drought-deciduous trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = plai (nLndPts,2,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_pfts02,jb)
    ENDIF

    IF(dodia(nDiag_pfts03))THEN  ! pft warm-temperate broadleaf evergreen trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = plai (nLndPts,3,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_pfts03,jb)
    ENDIF

    IF(dodia(nDiag_pfts04))THEN  ! pft temperate conifer evergreen trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = plai (nLndPts,4,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_pfts04,jb)
    ENDIF

    IF(dodia(nDiag_pfts05))THEN  ! pft temperate broadleaf cold-deciduous trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = plai (nLndPts,5,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_pfts05,jb)
    ENDIF

    IF(dodia(nDiag_pfts06))THEN  ! pft boreal conifer evergreen trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = plai (nLndPts,6,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_pfts06,jb)
    ENDIF

    IF(dodia(nDiag_pfts07))THEN  ! pft boreal broadleaf cold-deciduous trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = plai (nLndPts,7,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_pfts07,jb)
    ENDIF

    IF(dodia(nDiag_pfts08))THEN  ! pft boreal conifer cold-deciduous trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = plai (nLndPts,8,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_pfts08,jb)
    ENDIF

    IF(dodia(nDiag_pfts09))THEN  ! pft evergreen shrubs
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = plai (nLndPts,9,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_pfts09,jb)
    ENDIF

    IF(dodia(nDiag_pfts10))THEN  ! pft cold-deciduous shrubs
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = plai (nLndPts,10,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_pfts10,jb)
    ENDIF

    IF(dodia(nDiag_pfts11))THEN  ! pft cool (c4) grasses   
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = plai (nLndPts,11,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_pfts11,jb)
    ENDIF

    IF(dodia(nDiag_pfts12))THEN   ! pft cool (c3) grasses   
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = plai (nLndPts,12,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_pfts12,jb)
    ENDIF

    IF(dodia(nDiag_biol01))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = cbiol (nLndPts,1,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_biol01,jb)
    ENDIF

    IF(dodia(nDiag_biol02))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = cbiol (nLndPts,2,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_biol02,jb)
    ENDIF

    IF(dodia(nDiag_biol03))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = cbiol (nLndPts,3,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_biol03,jb)
    ENDIF

    IF(dodia(nDiag_biol04))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = cbiol (nLndPts,4,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_biol04,jb)
    ENDIF

    IF(dodia(nDiag_biol05))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = cbiol (nLndPts,5,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_biol05,jb)
    ENDIF

    IF(dodia(nDiag_biol06))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = cbiol (nLndPts,6,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_biol06,jb)
    ENDIF

    IF(dodia(nDiag_biol07))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = cbiol (nLndPts,7,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_biol07,jb)
    ENDIF

    IF(dodia(nDiag_biol08))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = cbiol (nLndPts,8,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_biol08,jb)
    ENDIF

    IF(dodia(nDiag_biol09))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = cbiol (nLndPts,9,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_biol09,jb)
    ENDIF

    IF(dodia(nDiag_biol10))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = cbiol (nLndPts,10,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_biol10,jb)
    ENDIF

    IF(dodia(nDiag_biol11))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = cbiol (nLndPts,11,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_biol11,jb)
    ENDIF

    IF(dodia(nDiag_biol12))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = cbiol (nLndPts,12,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_biol12,jb)
    ENDIF

    IF(dodia(nDiag_ynpp01))THEN  ! ynpp tropical broadleaf evergreen trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = aynpp (nLndPts,1,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_ynpp01,jb)
    ENDIF

    IF(dodia(nDiag_ynpp02))THEN ! ynpp tropical broadleaf drought-deciduous trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = aynpp (nLndPts,2,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_ynpp02,jb)
    ENDIF

    IF(dodia(nDiag_ynpp03))THEN  ! ynpp warm-temperate broadleaf evergreen trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = aynpp (nLndPts,3,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_ynpp03,jb)
    ENDIF

    IF(dodia(nDiag_ynpp04))THEN  ! ynpp temperate conifer evergreen trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = aynpp (nLndPts,4,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_ynpp04,jb)
    ENDIF

    IF(dodia(nDiag_ynpp05))THEN  ! ynpp temperate broadleaf cold-deciduous trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = aynpp (nLndPts,5,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_ynpp05,jb)
    ENDIF

    IF(dodia(nDiag_ynpp06))THEN  ! ynpp boreal conifer evergreen trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = aynpp (nLndPts,6,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_ynpp06,jb)
    ENDIF

    IF(dodia(nDiag_ynpp07))THEN  !ynpp boreal broadleaf cold-deciduous trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = aynpp (nLndPts,7,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_ynpp07,jb)
    ENDIF

    IF(dodia(nDiag_ynpp08))THEN  ! ynpp boreal conifer cold-deciduous trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = aynpp (nLndPts,8,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_ynpp08,jb)
    ENDIF

    IF(dodia(nDiag_ynpp09))THEN  ! ynpp evergreen shrubs
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = aynpp (nLndPts,9,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_ynpp09,jb)
    ENDIF

    IF(dodia(nDiag_ynpp10))THEN  ! ynpp cold-deciduous shrubs
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = aynpp (nLndPts,10,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_ynpp10,jb)
    ENDIF

    IF(dodia(nDiag_ynpp11))THEN  !! ynpp warm (c4) grasses
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = aynpp (nLndPts,11,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_ynpp11,jb)
    ENDIF

    IF(dodia(nDiag_ynpp12))THEN  ! ynpp cool (c3) grasses
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = aynpp (nLndPts,12,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_ynpp12,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cmontp = 170 !coldest monthly temperature                             (C)
    IF(dodia(nDiag_cmontp))THEN  !coldest monthly temperature                             (C)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = tc (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cmontp,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_wmontp = 171 !warmest monthly temperature                             (C)
    IF(dodia(nDiag_wmontp))THEN   !warmest monthly temperature                             (C)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = tw (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_wmontp,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_atogpp = 172 !annual total gpp for ecosystem                               (kg-c/m**2/yr)
    IF(dodia(nDiag_atogpp))THEN  !annual total gpp for ecosystem                               (kg-c/m**2/yr)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = aygpptot (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_atogpp,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_toigpp = 173 !instantaneous gpp                                (mol-CO2 / m-2 / second)
    IF(dodia(nDiag_toigpp))THEN  !instantaneous gpp                                (mol-CO2 / m-2 / second)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = tgpptot (nLndPts)
          END IF
       END DO
       CALL updia(diag,nDiag_toigpp,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_fxcsol = 174 !instantaneous fine co2 flux from soil                       (mol-CO2 / m-2 / second)
    IF(dodia(nDiag_fxcsol))THEN  !instantaneous fine co2 flux from soil                       (mol-CO2 / m-2 / second)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = tco2root(nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_fxcsol,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_mcsoil = 175 !instantaneous microbial co2 flux from soil       (mol-CO2 / m-2 / second)
    IF(dodia(nDiag_mcsoil))THEN  !instantaneous microbial co2 flux from soil       (mol-CO2 / m-2 / second)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = tco2mic (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_mcsoil,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cagcub = 176 !canopy average gross photosynthesis rate - broadleaf  (mol_co2 m-2 s-1)
    IF(dodia(nDiag_cagcub))THEN !canopy average gross photosynthesis rate - broadleaf  (mol_co2 m-2 s-1)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = agcub (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cagcub,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cagcuc = 177 !canopy average gross photosynthesis rate - conifer    (mol_co2 m-2 s-1)
    IF(dodia(nDiag_cagcuc))THEN  !canopy average gross photosynthesis rate - conifer    (mol_co2 m-2 s-1)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = agcuc (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cagcuc,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cagcls = 178 !canopy average gross photosynthesis rate - shrubs     (mol_co2 m-2 s-1)
    IF(dodia(nDiag_cagcls))THEN  !canopy average gross photosynthesis rate - shrubs     (mol_co2 m-2 s-1)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = agcls (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cagcls,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cagcl4 = 179 !canopy average gross photosynthesis rate - c4 grasses (mol_co2 m-2 s-1)
    IF(dodia(nDiag_cagcl4))THEN  !canopy average gross photosynthesis rate - c4 grasses (mol_co2 m-2 s-1)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = agcl4 (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cagcl4,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cagcl3 = 180 !canopy average gross photosynthesis rate - c3 grasses (mol_co2 m-2 s-1)
    IF(dodia(nDiag_cagcl3))THEN  !canopy average gross photosynthesis rate - c3 grasses (mol_co2 m-2 s-1)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = agcl3 (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cagcl3,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cancub = 181 !canopy average net photosynthesis rate - broadleaf    (mol_co2 m-2 s-1)
    IF(dodia(nDiag_cancub))THEN  !canopy average net photosynthesis rate - broadleaf    (mol_co2 m-2 s-1)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = ancub (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cancub,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cancuc = 182 !canopy average net photosynthesis rate - conifer      (mol_co2 m-2 s-1)
    IF(dodia(nDiag_cancuc))THEN  !canopy average net photosynthesis rate - conifer      (mol_co2 m-2 s-1)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = ancuc (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cancuc,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cancls = 183 !canopy average net photosynthesis rate - shrubs            (mol_co2 m-2 s-1)
    IF(dodia(nDiag_cancls))THEN  !canopy average net photosynthesis rate - shrubs            (mol_co2 m-2 s-1)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = ancls (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cancls,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cancl4 = 184 !canopy average net photosynthesis rate - c4 grasses   (mol_co2 m-2 s-1)
    IF(dodia(nDiag_cancl4))THEN  !canopy average net photosynthesis rate - c4 grasses   (mol_co2 m-2 s-1)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = ancl4 (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cancl4,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cancl3 = 185 !canopy average net photosynthesis rate - c3 grasses   (mol_co2 m-2 s-1)
    IF(dodia(nDiag_cancl3))THEN  !canopy average net photosynthesis rate - c3 grasses   (mol_co2 m-2 s-1)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = ancl3 (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cancl3,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cicoub = 186 !intercellular co2 concentration - broadleaf            (mol_co2/mol_air)
    IF(dodia(nDiag_cicoub))THEN  !intercellular co2 concentration - broadleaf            (mol_co2/mol_air)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = ciub (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cicoub,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cicouc = 187 !intercellular co2 concentration - conifer             (mol_co2/mol_air)
    IF(dodia(nDiag_cicouc))THEN  !intercellular co2 concentration - conifer             (mol_co2/mol_air)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = ciuc (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cicouc,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cscoub = 188 !leaf boundary layer co2 concentration - broadleaf     (mol_co2/mol_air)
    IF(dodia(nDiag_cscoub))THEN  ! !leaf boundary layer co2 concentration - broadleaf     (mol_co2/mol_air)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = csub (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cscoub,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_gscoub = 189 !upper canopy stomatal conductance - broadleaf            (mol_co2 m-2 s-1)
    IF(dodia(nDiag_gscoub))THEN  !upper canopy stomatal conductance - broadleaf            (mol_co2 m-2 s-1)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = gsub (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_gscoub,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cscouc = 190 !leaf boundary layer co2 concentration - conifer            (mol_co2/mol_air)
    IF(dodia(nDiag_cscouc))THEN   !leaf boundary layer co2 concentration - conifer            (mol_co2/mol_air)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = csuc (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cscouc,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_gscouc = 191 !upper canopy stomatal conductance - conifer            (mol_co2 m-2 s-1)
    IF(dodia(nDiag_gscouc))THEN  !upper canopy stomatal conductance - conifer            (mol_co2 m-2 s-1)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = gsuc (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_gscouc,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cicols = 192 !intercellular co2 concentration - shrubs              (mol_co2/mol_air)
    IF(dodia(nDiag_cicols))THEN  !intercellular co2 concentration - shrubs              (mol_co2/mol_air)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = cils (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cicols,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cicol3 = 193 !intercellular co2 concentration - c3 plants            (mol_co2/mol_air)
    IF(dodia(nDiag_cicol3))THEN  !intercellular co2 concentration - c3 plants            (mol_co2/mol_air)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = cil3 (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cicol3,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cicol4 = 194 !intercellular co2 concentration - c4 plants            (mol_co2/mol_air)
    IF(dodia(nDiag_cicol4))THEN  !intercellular co2 concentration - c4 plants            (mol_co2/mol_air)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = cil4 (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cicol4,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cscols = 195 !leaf boundary layer co2 concentration - shrubs            (mol_co2/mol_air)
    IF(dodia(nDiag_cscols))THEN  !leaf boundary layer co2 concentration - shrubs            (mol_co2/mol_air)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = csls (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cscols,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_gscols = 196 !lower canopy stomatal conductance - shrubs            (mol_co2 m-2 s-1)
    IF(dodia(nDiag_gscols))THEN  !lower canopy stomatal conductance - shrubs            (mol_co2 m-2 s-1)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = gsls (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_gscols,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cscol3 = 197 !leaf boundary layer co2 concentration - c3 plants     (mol_co2/mol_air)
    IF(dodia(nDiag_cscol3))THEN  !leaf boundary layer co2 concentration - c3 plants     (mol_co2/mol_air)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = csl3 (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cscol3,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_gscol3 = 198 !lower canopy stomatal conductance - c3 grasses            (mol_co2 m-2 s-1)
    IF(dodia(nDiag_gscol3))THEN  !lower canopy stomatal conductance - c3 grasses            (mol_co2 m-2 s-1)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = gsl3 (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_gscol3,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cscol4 = 199 !leaf boundary layer co2 concentration - c4 plants     (mol_co2/mol_air)
    IF(dodia(nDiag_cscol4))THEN  !leaf boundary layer co2 concentration - c4 plants     (mol_co2/mol_air)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = csl4 (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cscol4,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_gscol4 = 200 !lower canopy stomatal conductance - c4 grasses            (mol_co2 m-2 s-1)
    IF(dodia(nDiag_gscol4))THEN  ! lower canopy stomatal conductance - c4 grasses            (mol_co2 m-2 s-1)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = gsl4 (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_gscol4,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_tcthis = 201 !coldest monthly temperature of current year              (C)
    IF(dodia(nDiag_tcthis))THEN  !coldest monthly temperature of current year              (C)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = tcthis (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_tcthis,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_twthis = 201 !warmest monthly temperature of current year              (C)
    IF(dodia(nDiag_twthis))THEN  !warmest monthly temperature of current year              (C)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = twthis (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_twthis,jb)
    ENDIF

!         ynleach_p (1:nLndPts,jb) ,tnmin_p   (1:nLndPts,jb) ,totnmic_p (1:nLndPts,jb) ,totnlit_p (1:nLndPts,jb) ,&
!         totanlit_p(1:nLndPts,jb) ,totrnlit_p(1:nLndPts,jb) ,totnsoi_p (1:nLndPts,jb) ,storedn_p (1:nLndPts,jb) )

    idateprev=idatec
  END SUBROUTINE Ibis_Interface



  SUBROUTINE IbisDrv(mcsec        ,pi          ,stef         ,vonk         ,grav        , &
       tmelt        ,hfus        ,hvap         ,hsub         ,ch2o        , &
       cice         ,cair        ,cvap         ,rair         ,rvap        , &
       cappa        ,rhow        ,npoi         ,nband        ,nsoilay     , & 
       nsnolay      ,npft        ,epsilon      ,dtime        ,doalb       , &
       ginvap       ,gsuvap      ,gtrans       ,gtransu      ,gtransl     , &
       grunof       ,gdrain      ,gadjust      ,a10scalparamu,a10daylightu, &
       a10scalparaml,a10daylightl,vmax_pft     ,tau15        ,kc15        , &
       ko15         ,cimax       ,gammaub      ,alpha3       ,theta3      , &
       beta3        ,coefmub     ,coefbub      ,gsubmin      ,gammauc     , &
       coefmuc      ,coefbuc     ,gsucmin      ,gammals      ,coefmls     , & 
       coefbls      ,gslsmin     ,gammal3      ,coefml3      ,coefbl3     , &
       gsl3min      ,gammal4     ,alpha4       ,theta4       ,beta4       , &
       coefml4      ,coefbl4     ,gsl4min      ,bps          , wliqu        ,wliqumax    , & 
       wsnou        ,wsnoumax    ,tu           ,wliqs        ,wliqsmax    , & 
       wsnos        ,wsnosmax    ,ts           ,wliql        ,wliqlmax    , &  
       wsnol        ,wsnolmax    ,tl           ,topparu      ,topparl     , &
       fl           ,fu          ,lai          ,sai          ,rhoveg      , &   
       tauveg       ,orieh       ,oriev        ,wliqmin      ,wsnomin     , & 
       t12          ,tdripu      ,tblowu       ,tdrips       ,tblows      , &
       t34          ,tdripl      ,tblowl       ,ztop         ,za          ,alaiml      , &
       zbot         ,alaimu      ,froot        ,q34          ,q12         , &
       su                 ,cleaf       ,dleaf        ,ss           ,cstem       , & 
       dstem        ,sl          ,cgrass       ,ciub         ,ciuc        , &
       exist        ,csub        ,gsub         ,csuc         ,gsuc        , &
       agcub        ,agcuc       ,ancub        ,ancuc        ,totcondub   , &
       totconduc    ,cils        ,cil3         ,cil4         ,csls        , &
       gsls         ,csl3        ,gsl3         ,csl4         ,gsl4        , &
       agcls        ,agcl4       ,agcl3        ,ancls        ,ancl4       , &
       ancl3       ,totcondls    ,totcondl3    ,totcondl4    ,chu         , &
       chs          ,chl         ,frac         ,tlsub          ,z0sno       , & 
       rhos         ,consno      ,hsnotop      ,hsnomin      ,fimin       , &
       fimax        ,fi          ,tsno         ,hsno         ,sand               , &
       clay         ,poros       ,wsoi         ,wisoi        ,consoi      , &  
       zwpmax       ,wpud        ,wipud        ,wpudmax      ,qglif       , &         
       tsoi         ,hvasug      ,hvasui       ,albsav       ,albsan      , &
       tg           ,ti          ,z0soi        ,swilt        ,sfield      , &
       stressl      ,stressu     ,stresstl     ,stresstu     ,csoi        , &         
       rhosoi       ,hsoi        ,suction      ,bex          ,upsoiu      , &  
       upsoil       ,heatg       ,heati        ,hydraul      ,porosflo    , &
       ibex         ,bperm       ,hflo            ,ta           ,asurd       , &  
       asuri        ,coszen      ,solad        ,solai        ,fira               , & 
       raina        ,qa          ,psurf        ,snowa        ,ua               , &   
       o2conc       ,co2conc     ,&
       td                 ,vzero       ,ndaypy       ,nppdummy     , & 
       cbiow        ,sapfrac      ,cbior       , & 
       tco2root    ,tneetot      ,tco2mic      ,a10td       , &
       a10ancub     ,a10ancuc    ,a10ancls     ,a10ancl3     ,a10ancl4    , & 
       ndtimes      ,adrain       ,adsnow      ,tnpptot, &
       adaet        ,adtrunoff   ,adsrunoff    ,addrainage   ,adrh        , &
       adsnod       ,adsnof      ,adwsoi       ,adtsoi       ,adwisoi     , &
       adtlaysoi    ,adwlaysoi   ,adwsoic      ,adtsoic      ,adco2mic    , &
       adco2root    ,adco2soi    ,adco2ratio   ,adnmintot    ,decompl     , &    
       decomps      ,tnmin       ,ndaypm       ,nmtimes     , & 
       amrain       ,amsnow      ,amaet        ,amtrunoff    ,amsrunoff   , &
       amdrainage   ,amtemp      ,amqa         , &
       amsolar      ,amirup      ,amirdown     ,amsens       ,amlatent    , &  
       amlaiu       ,amlail      ,amtsoi       ,amwsoi       ,amwisoi     , &   
       amvwc        ,amawc       ,amsnod       ,amsnof       ,amnpp       , &
       amnpptot     ,amco2mic    ,amco2root    ,amco2soi     ,amco2ratio  , &
       amneetot     ,amnmintot   ,amalbedo     ,amtsoil      ,amwsoil     , & 
       amwisoil     ,nytimes      ,aysolar      ,ayirup      , &
       ayirdown     ,aysens      ,aylatent     ,ayprcp       ,ayaet       , &  
       aytrans      ,aytrunoff   ,aysrunoff    ,aydrainage   ,aydwtot     , & 
       aywsoi       ,aywisoi     ,aytsoi       ,ayvwc        ,ayawc       , &  
       aystresstu   ,aystresstl  ,aygpp            ,aygpptot     ,aynpp       , & 
       aynpptot     ,ayco2mic    ,ayco2root    ,ayco2soi     ,ayneetot    , &
       ayrootbio    ,aynmintot   ,ayalit       ,ayblit       ,aycsoi      , & 
       aycmic       ,ayanlit     ,aybnlit      ,aynsoi       ,ayalbedo    ,&
       totalit     , &
       totrlit      ,totcsoi     ,totcmic      ,totanlit     ,totrnlit    , &
       totnsoi      ,totnmic     ,totlit       ,totfall      ,totnlit     , &
       firefac      ,wtot        ,storedn      ,yrleach      ,ynleach     , & 
       falll        ,fallr       ,fallw        ,clitlm       ,clitls      , &
       clitrm       ,clitrs      ,clitwm       ,clitws       ,csoislop    , &
       csoislon     ,csoipas     ,clitll       ,clitrl       ,clitwl      , &  
       tc           ,agddu              ,tempu            ,agddl        ,templ       , &
       dropu        ,dropls      ,dropl4       ,dropl3       ,plai               , &
       iday        ,imonth       ,iyear        ,iyear0       , &
       isimveg      ,spinmax     ,amts2        ,amtransu     ,amtransl    , &
       amsuvap      ,aminvap     ,ux            ,uy                  ,taux        , &
       tauy         ,ts2              ,qs2          ,deltat       ,gdd0        , &  
       gdd0this     ,tcthis      ,twthis       ,tcmin        ,gdd5        , & 
       gdd5this     ,TminL       ,TminU        ,Twarm        ,GDD         , & 
       aleaf        ,awood       ,cbiol        ,aroot        ,disturbf    , &
       disturbo     ,specla      ,biomass      ,totlaiu      ,totlail     , &  
       totbiou      ,totbiol     ,woodnorm     ,vegtype0     ,tauwood0    , & 
       tauwood      ,tauleaf     ,tauroot      ,xminlai      ,cdisturb    , & 
       ayanpp       ,ayanpptot   , &
       nstep         ,imonthprev   , &
       isimco2      ,isimfire    , &
       co2initm     ,tw           ,fvapa       ,fsena        , &
       z0           ,ustar       ,hc           ,hg          ,ec           , &
       eg           ,dispu       ,cu           ,firb        ,tgpptot      , &
       bstar        ,ynleach_p   ,tnmin_p      ,totnmic_p   ,totnlit_p    , &
       totanlit_p   ,totrnlit_p  ,totnsoi_p    ,storedn_p   ,adnpp        , &
       adfalll      ,adfallr     ,adfallw      ,adcbiol     ,adcbior      , &
       adcbiow      ,adplai      ,beta1        ,beta2       ,&
                     stressfac,avmuir_factor, nVegClass ,rootmode)

    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN        ) :: mcsec!global  ! current seconds in day (0 - (86400 - dtime))
    REAL(KIND=r8), INTENT(IN        ) :: pi   !global
    REAL(KIND=r8), INTENT(IN        ) :: stef !global  ! stefan-boltzmann constant (W m-2 K-4)
    REAL(KIND=r8), INTENT(IN        ) :: vonk !global  ! von karman constant (dimensionless)
    REAL(KIND=r8), INTENT(IN        ) :: grav !global  ! gravitational acceleration (m s-2)
    REAL(KIND=r8), INTENT(IN        ) :: tmelt!global  ! freezing point of water (K)
    REAL(KIND=r8), INTENT(IN        ) :: hfus !global  ! latent heat of fusion of water (J kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: hvap !global  ! latent heat of vaporization of water (J kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: hsub !global  ! latent heat of sublimation of ice (J kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: ch2o !global  ! specific heat of liquid water (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: cice !global  ! specific heat of ice (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: cair !global  ! specific heat of dry air at constant pressure (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: cvap !global  ! specific heat of water vapor at constant pressure (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: rair !global  ! gas constant for dry air (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: rvap !global  ! gas constant for water vapor (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: cappa!global  ! rair/cair
    REAL(KIND=r8), INTENT(IN        ) :: rhow !global  ! density of liquid water (all types) (kg m-3)

    ! 
    !
    INTEGER, INTENT(IN   ) :: nVegClass
    INTEGER, INTENT(IN   ) :: npoi   !global  
    INTEGER, INTENT(IN   ) :: nband  !global  
    INTEGER, INTENT(IN   ) :: nsoilay!global   ! number of soil layers
    INTEGER, INTENT(IN   ) :: nsnolay!global   ! number of snow layers
    INTEGER, INTENT(IN   ) :: npft   !global   ! number of plant functional types
    REAL(KIND=r8), INTENT(IN   ) :: epsilon!global   ! small quantity to avoid zero-divides and other
    ! truncation or machine-limit troubles with small
    ! values. should be slightly greater than o(1)
    ! machine precision
    REAL(KIND=r8), INTENT(IN   )  :: dtime !global   ! model timestep (seconds)
    LOGICAL, INTENT(IN   )  :: doalb !global    ! true if surface albedo calculation time step

    !      INCLUDE 'comhyd.h'
    REAL(KIND=r8), INTENT(OUT  ) :: ginvap (npoi)!local ! total evaporation rate from all intercepted h2o (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: gsuvap (npoi)!local ! total evaporation rate from surface (snow/soil) (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: gtrans (npoi)!local ! total transpiration rate from all vegetation canopies (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: gtransu(npoi)!local ! transpiration from upper canopy (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: gtransl(npoi)!local ! transpiration from lower canopy (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: grunof (npoi)!local ! surface runoff rate (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: gdrain (npoi)!local ! drainage rate out of bottom of lowest soil layer (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: gadjust(npoi)!local ! h2o flux due to adjustments in subroutine wadjust (kg_h2o m-2 s-1)
    !      INCLUDE 'comsum.h'
    REAL(KIND=r8), INTENT(INOUT) :: a10scalparamu(npoi)!global ! 10-day average day-time scaling parameter - upper canopy (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: a10daylightu (npoi)!global ! 10-day average day-time PAR - upper canopy (micro-Ein m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: a10scalparaml(npoi)!global ! 10-day average day-time scaling parameter - lower canopy (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: a10daylightl (npoi)!global ! 10-day average day-time PAR - lower canopy (micro-Ein m-2 s-1)
    !      INCLUDE 'compft.h'
    REAL(KIND=r8), INTENT(IN   ) :: vmax_pft(npft)!global ! nominal vmax of top leaf at 15 C (mol-co2/m**2/s) [not used]
    REAL(KIND=r8), INTENT(IN   ) :: tau15           !global ! co2/o2 specificity ratio at 15 degrees C (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: kc15           !global ! co2 kinetic parameter (mol/mol)
    REAL(KIND=r8), INTENT(IN   ) :: ko15           !global ! o2 kinetic parameter (mol/mol) 
    REAL(KIND=r8), INTENT(IN   ) :: cimax           !global ! maximum value for ci (needed for model stability)
    REAL(KIND=r8), INTENT(IN   ) :: gammaub           !global ! leaf respiration coefficient
    REAL(KIND=r8), INTENT(IN   ) :: alpha3           !global ! intrinsic quantum efficiency for C3 plants (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: theta3           !global ! photosynthesis coupling coefficient for C3 plants (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: beta3           !global ! photosynthesis coupling coefficient for C3 plants (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: coefmub           !global ! 'm' coefficient for stomatal conductance relationship
    REAL(KIND=r8), INTENT(IN   ) :: coefbub           !global ! 'b' coefficient for stomatal conductance relationship
    REAL(KIND=r8), INTENT(IN   ) :: gsubmin           !global ! absolute minimum stomatal conductance
    REAL(KIND=r8), INTENT(IN   ) :: gammauc           !global ! leaf respiration coefficient
    REAL(KIND=r8), INTENT(IN   ) :: coefmuc           !global ! 'm' coefficient for stomatal conductance relationship  
    REAL(KIND=r8), INTENT(IN   ) :: coefbuc           !global ! 'b' coefficient for stomatal conductance relationship  
    REAL(KIND=r8), INTENT(IN   ) :: gsucmin           !global ! absolute minimum stomatal conductance
    REAL(KIND=r8), INTENT(IN   ) :: gammals           !global ! leaf respiration coefficient
    REAL(KIND=r8), INTENT(IN   ) :: coefmls           !global ! 'm' coefficient for stomatal conductance relationship 
    REAL(KIND=r8), INTENT(IN   ) :: coefbls           !global ! 'b' coefficient for stomatal conductance relationship 
    REAL(KIND=r8), INTENT(IN   ) :: gslsmin           !global ! absolute minimum stomatal conductance
    REAL(KIND=r8), INTENT(IN   ) :: gammal3           !global ! leaf respiration coefficient
    REAL(KIND=r8), INTENT(IN   ) :: coefml3           !global ! 'm' coefficient for stomatal conductance relationship 
    REAL(KIND=r8), INTENT(IN   ) :: coefbl3           !global ! 'b' coefficient for stomatal conductance relationship 
    REAL(KIND=r8), INTENT(IN   ) :: gsl3min           !global ! absolute minimum stomatal conductance
    REAL(KIND=r8), INTENT(IN   ) :: gammal4           !global ! leaf respiration coefficient
    REAL(KIND=r8), INTENT(IN   ) :: alpha4           !global ! intrinsic quantum efficiency for C4 plants (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: theta4           !global ! photosynthesis coupling coefficient for C4 plants (dimensionless) 
    REAL(KIND=r8), INTENT(IN   ) :: beta4           !global ! photosynthesis coupling coefficient for C4 plants (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: coefml4           !global ! 'm' coefficient for stomatal conductance relationship
    REAL(KIND=r8), INTENT(IN   ) :: coefbl4           !global ! 'b' coefficient for stomatal conductance relationship
    REAL(KIND=r8), INTENT(IN   ) :: gsl4min           !global ! absolute minimum stomatal conductance
    REAL(KIND=r8), INTENT(INOUT) :: bps(npoi)
    !      include 'comveg.h'
    REAL(KIND=r8), INTENT(INOUT) :: wliqu    (npoi)  !global ! intercepted liquid h2o on upper canopy leaf area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wliqumax              !global ! maximum intercepted water on a unit upper canopy leaf area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: wsnou    (npoi)  !global ! intercepted frozen h2o (snow) on upper canopy leaf area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wsnoumax              !global ! intercepted snow capacity for upper canopy leaves (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: tu       (npoi)  !global ! temperature of upper canopy leaves (K)
    REAL(KIND=r8), INTENT(INOUT) :: wliqs    (npoi)  !global ! intercepted liquid h2o on upper canopy stem area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wliqsmax              !global ! maximum intercepted water on a unit upper canopy stem area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: wsnos    (npoi)  !global ! intercepted frozen h2o (snow) on upper canopy stem area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wsnosmax              !global ! intercepted snow capacity for upper canopy stems (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: ts       (npoi)  !global ! temperature of upper canopy stems (K)
    REAL(KIND=r8), INTENT(INOUT) :: wliql    (npoi)  !global ! intercepted liquid h2o on lower canopy leaf and stem area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wliqlmax              !global ! maximum intercepted water on a unit lower canopy stem & leaf area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: wsnol    (npoi)  !global ! intercepted frozen h2o (snow) on lower canopy leaf & stem area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wsnolmax              !global ! intercepted snow capacity for lower canopy leaves & stems (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: tl       (npoi)  !global ! temperature of lower canopy leaves & stems(K)
    REAL(KIND=r8), INTENT(INOUT) :: topparu  (npoi)  !local  ! total photosynthetically active raditaion absorbed 
    ! by top leaves of upper canopy (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: topparl  (npoi)  !local  ! total photosynthetically active raditaion absorbed
    ! by top leaves of lower canopy (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: fl       (npoi)   !global ! fraction of snow-free area covered by lower  canopy
    REAL(KIND=r8), INTENT(INOUT) :: fu       (npoi)   !global ! fraction of overall area covered by upper canopy
    REAL(KIND=r8), INTENT(INOUT) :: lai      (npoi,2) !global ! canopy single-sided leaf area index (area leaf/area veg)
    REAL(KIND=r8), INTENT(INOUT) :: sai      (npoi,2) !global ! current single-sided stem area index
    REAL(KIND=r8), INTENT(IN   ) :: rhoveg   (nband,2)!global ! reflectance of an average leaf/stem
    REAL(KIND=r8), INTENT(IN   ) :: tauveg   (nband,2)!global  ! transmittance of an average leaf/stem
    REAL(KIND=r8), INTENT(IN   ) :: orieh    (2)      !global! fraction of leaf/stems with horizontal orientation
    REAL(KIND=r8), INTENT(IN   ) :: oriev    (2)      !global! fraction of leaf/stems with vertical
    REAL(KIND=r8), INTENT(INOUT) :: wliqmin           !local ! minimum intercepted water on unit vegetated area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: wsnomin           !local ! minimum intercepted snow on unit vegetated area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: t12      (npoi)   !global ! air temperature at z12 (K)
    REAL(KIND=r8), INTENT(IN   ) :: tdripu            !global ! decay time for dripoff of liquid intercepted by upper canopy leaves (sec)
    REAL(KIND=r8), INTENT(IN   ) :: tblowu             !global ! decay time for blowoff of snow intercepted by upper canopy leaves (sec)
    REAL(KIND=r8), INTENT(IN   ) :: tdrips             !global ! decay time for dripoff of liquid intercepted by upper canopy stems (sec) 
    REAL(KIND=r8), INTENT(IN   ) :: tblows             !global ! decay time for blowoff of snow intercepted by upper canopy stems (sec)
    REAL(KIND=r8), INTENT(INOUT) :: t34      (npoi)   !global ! air temperature at z34 (K)
    REAL(KIND=r8), INTENT(IN   ) :: tdripl            !global ! decay time for dripoff of liquid intercepted
    ! by lower canopy leaves & stem (sec)
    REAL(KIND=r8), INTENT(IN   ) :: tblowl             ! global          ! decay time for blowoff of snow intercepted by lower canopy leaves & stems (sec)
    REAL(KIND=r8), INTENT(INOUT) :: za       (npoi)          ! local ! height above the surface of atmospheric forcing (m)

    REAL(KIND=r8), INTENT(INOUT) :: ztop     (npoi,2) ! global  ! height of plant top above ground (m)
    REAL(KIND=r8), INTENT(IN   ) :: alaiml             ! global ! lower canopy leaf & stem maximum area (2 sided) for
    ! normalization of drag coefficient (m2 m-2)
    REAL(KIND=r8), INTENT(INOUT) :: zbot     (npoi,2) ! global  ! height of lowest branches above ground (m)
    REAL(KIND=r8), INTENT(IN   ) :: alaimu             ! global  ! upper canopy leaf & stem area (2 sided) for 
    ! normalization of drag coefficient (m2 m-2)
    REAL(KIND=r8), INTENT(INOUT) :: froot    (npoi,nsoilay,2)! global! fraction of root in soil layer 
    REAL(KIND=r8), INTENT(INOUT) :: q34      (npoi)   ! global! specific humidity of air at z34
    REAL(KIND=r8), INTENT(INOUT) :: q12      (npoi)   ! global! specific humidity of air at z12
    REAL(KIND=r8), INTENT(INOUT) :: su       (npoi)   ! local ! air-vegetation transfer coefficients (*rhoa) for
    ! upper canopy leaves (m s-1 * kg m-3) (A39a Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(IN   ) :: cleaf             ! global! empirical constant in upper canopy leaf-air 
    ! aerodynamic transfer coefficient (m s-0.5) (A39a Pollard & Thompson 95)
    REAL(KIND=r8), INTENT(IN   ) :: dleaf    (2)             ! global ! typical linear leaf dimension in aerodynamic transfer coefficient (m)
    REAL(KIND=r8), INTENT(INOUT) :: ss       (npoi)   ! local! air-vegetation transfer coefficients (*rhoa) for 
    ! upper canopy stems (m s-1 * kg m-3) (A39a Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(IN   ) :: cstem             ! global ! empirical constant in upper canopy stem-air 
    ! aerodynamic transfer coefficient (m s-0.5) (A39a Pollard & Thompson 95)
    REAL(KIND=r8), INTENT(IN   ) :: dstem    (2)             ! global ! typical linear stem dimension in aerodynamic transfer coefficient (m)
    REAL(KIND=r8), INTENT(INOUT) :: sl       (npoi)   ! local ! air-vegetation transfer coefficients (*rhoa) for 
    ! lower canopy leaves & stems (m s-1*kg m-3) (A39a Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(IN   ) :: cgrass             ! global ! empirical constant in lower canopy-air aerodynamic 
    ! transfer coefficient (m s-0.5) (A39a Pollard & Thompson 95)
    REAL(KIND=r8), INTENT(INOUT) :: ciub     (npoi)         ! global ! intercellular co2 concentration - broadleaf (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: ciuc     (npoi)         ! global ! intercellular co2 concentration - conifer        (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: exist    (npoi,npft)  ! global ! probability of existence of each plant functional type in a gridcell
    REAL(KIND=r8), INTENT(INOUT) :: csub     (npoi)         ! global ! leaf boundary layer co2 concentration - broadleaf (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: gsub     (npoi)         ! global ! upper canopy stomatal conductance - broadleaf  (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: csuc     (npoi)         ! global ! leaf boundary layer co2 concentration - conifer   (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: gsuc     (npoi)         ! global ! upper canopy stomatal conductance - conifer    (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: agcub    (npoi)         ! local  ! canopy average gross photosynthesis rate - broadleaf  (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: agcuc    (npoi)         ! local  ! canopy average gross photosynthesis rate - conifer    (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: ancub    (npoi)         ! local  ! canopy average net photosynthesis rate - broadleaf    (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: ancuc    (npoi)         ! local  ! canopy average net photosynthesis rate - conifer          (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: totcondub(npoi)         ! local  ! 
    REAL(KIND=r8), INTENT(INOUT) :: totconduc(npoi)         ! local  !
    REAL(KIND=r8), INTENT(INOUT) :: cils     (npoi)         ! global ! intercellular co2 concentration - shrubs        (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: cil3     (npoi)         ! global ! intercellular co2 concentration - c3 plants (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: cil4     (npoi)         ! global ! intercellular co2 concentration - c4 plants (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: csls     (npoi)         ! global ! leaf boundary layer co2 concentration - shrubs   (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: gsls     (npoi)         ! global ! lower canopy stomatal conductance - shrubs     (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: csl3     (npoi)         ! global ! leaf boundary layer co2 concentration - c3 plants (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: gsl3     (npoi)         ! global ! lower canopy stomatal conductance - c3 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: csl4     (npoi)         ! global ! leaf boundary layer co2 concentration - c4 plants (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: gsl4     (npoi)         ! global ! lower canopy stomatal conductance - c4 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: agcls    (npoi)         ! local  ! canopy average gross photosynthesis rate - shrubs          (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: agcl4    (npoi)         ! local ! canopy average gross photosynthesis rate - c4 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: agcl3    (npoi)         ! local ! canopy average gross photosynthesis rate - c3 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: ancls    (npoi)         ! local ! canopy average net photosynthesis rate - shrubs          (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: ancl4    (npoi)         ! local ! canopy average net photosynthesis rate - c4 grasses   (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: ancl3    (npoi)         ! local ! canopy average net photosynthesis rate - c3 grasses   (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: totcondls(npoi)         ! local ! 
    REAL(KIND=r8), INTENT(INOUT) :: totcondl3(npoi)         ! local !
    REAL(KIND=r8), INTENT(INOUT) :: totcondl4(npoi)         ! local !
    REAL(KIND=r8), INTENT(IN   ) :: chu(1:nVegClass)             ! global ! heat capacity of upper canopy leaves per unit leaf area (J kg-1 m-2)
    REAL(KIND=r8), INTENT(IN   ) :: chs(1:nVegClass)             ! global ! heat capacity of upper canopy stems per unit stem area (J kg-1 m-2)
    REAL(KIND=r8), INTENT(IN   ) :: chl(1:nVegClass)             ! global ! heat capacity of lower canopy leaves & stems per unit leaf/stem area (J kg-1 m-2)
    REAL(KIND=r8), INTENT(INOUT) :: frac     (npoi,npft)  ! global ! fraction of canopy occupied by each plant functional type
    REAL(KIND=r8), INTENT(INOUT) :: tlsub    (npoi)         ! global ! temperature of lower canopy vegetation buried by snow (K)
    !      INCLUDE 'comsat.h'    
    !      include 'comsno.h'
    REAL(KIND=r8), INTENT(IN   ) :: z0sno  ! global ! roughness length of snow surface (m)
    REAL(KIND=r8), INTENT(IN   ) :: rhos   ! global ! density of snow (kg m-3)
    REAL(KIND=r8), INTENT(IN   ) :: consno ! global ! thermal conductivity of snow (W m-1 K-1)
    REAL(KIND=r8), INTENT(IN   ) :: hsnotop! global ! thickness of top snow layer (m)
    REAL(KIND=r8), INTENT(IN   ) :: hsnomin! global ! minimum total thickness of snow (m)
    REAL(KIND=r8), INTENT(IN   ) :: fimin  ! global ! minimum fractional snow cover
    REAL(KIND=r8), INTENT(IN   ) :: fimax  ! global ! maximum fractional snow cover
    REAL(KIND=r8), INTENT(INOUT) :: fi     (npoi)! global ! fractional snow cover
    REAL(KIND=r8), INTENT(INOUT) :: tsno   (npoi,nsnolay)! global ! temperature of snow layers (K)
    REAL(KIND=r8), INTENT(INOUT) :: hsno   (npoi,nsnolay)! global ! thickness of snow layers (m)

    !      INCLUDE 'comsoi.h'
    REAL(KIND=r8), INTENT(IN   ) :: sand    (npoi,nsoilay)! global ! percent sand of soil
    REAL(KIND=r8), INTENT(IN   ) :: clay    (npoi,nsoilay)! global ! percent clay of soil
    REAL(KIND=r8), INTENT(IN   ) :: poros   (npoi,nsoilay)! global ! porosity (mass of h2o per unit vol at sat / rhow)
    REAL(KIND=r8), INTENT(INOUT) :: wsoi    (npoi,nsoilay)! global ! fraction of soil pore space containing liquid water
    REAL(KIND=r8), INTENT(INOUT) :: wisoi   (npoi,nsoilay)! global ! fraction of soil pore space containing ice
    REAL(KIND=r8), INTENT(INOUT) :: consoi  (npoi,nsoilay)! local  ! thermal conductivity of each soil layer (W m-1 K-1)
    REAL(KIND=r8), INTENT(IN   ) :: zwpmax                 ! global ! assumed maximum fraction of soil surface 
    ! covered by puddles (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: wpud    (npoi)! global ! liquid content of puddles per soil area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: wipud   (npoi)! global ! ice content of puddles per soil area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wpudmax         ! global ! normalization constant for puddles (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: qglif   (npoi,4) ! local ! 1: fraction of soil evap (fvapg) from soil liquid
    ! 2: fraction of soil evap (fvapg) from soil ice
    ! 3: fraction of soil evap (fvapg) from puddle liquid
    ! 4: fraction of soil evap (fvapg) from puddle ice
    REAL(KIND=r8), INTENT(INOUT) :: tsoi    (npoi,nsoilay)! global        ! soil temperature for each layer (K)
    REAL(KIND=r8), INTENT(INOUT) :: hvasug  (npoi)        ! local ! latent heat of vap/subl, for soil surface (J kg-1)
    REAL(KIND=r8), INTENT(INOUT) :: hvasui  (npoi)        ! local ! latent heat of vap/subl, for snow surface (J kg-1)
    REAL(KIND=r8), INTENT(IN   ) :: albsav  (npoi)        ! global ! saturated soil surface albedo (visible waveband)
    REAL(KIND=r8), INTENT(IN   ) :: albsan  (npoi)        ! global ! saturated soil surface albedo (near-ir waveband)
    REAL(KIND=r8), INTENT(INOUT) :: tg      (npoi)        ! global ! soil skin temperature (K)
    REAL(KIND=r8), INTENT(INOUT) :: ti      (npoi)        ! global ! snow skin temperature (K)
    REAL(KIND=r8), INTENT(IN   ) :: z0soi   (npoi)        ! global ! roughness length of soil surface (m)
    REAL(KIND=r8), INTENT(IN   ) :: swilt   (npoi,nsoilay)! global ! wilting soil moisture value (fraction of pore space)
    REAL(KIND=r8), INTENT(IN   ) :: sfield  (npoi,nsoilay)! global ! field capacity soil moisture value (fraction of pore space)
    REAL(KIND=r8), INTENT(INOUT) :: stressl (npoi,nsoilay)! local ! soil moisture stress factor for the lower canopy (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: stressu (npoi,nsoilay)! local ! soil moisture stress factor for the upper canopy (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: stresstl(npoi)        ! local ! sum of stressl over all 6 soil layers (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: stresstu(npoi)        ! local ! sum of stressu over all 6 soil layers (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: csoi    (npoi,nsoilay)! global ! specific heat of soil, no pore spaces (J kg-1 deg-1)
    REAL(KIND=r8), INTENT(IN   ) :: rhosoi  (npoi,nsoilay)! global ! soil density (without pores, not bulk) (kg m-3)
    REAL(KIND=r8), INTENT(IN   ) :: hsoi    (npoi,nsoilay+1)   ! global ! soil layer thickness (m)
    REAL(KIND=r8), INTENT(IN   ) :: suction (npoi,nsoilay)! global ! saturated matric potential (m-h2o)
    REAL(KIND=r8), INTENT(IN   ) :: bex     (npoi,nsoilay)! global ! exponent "b" in soil water potential
    REAL(KIND=r8), INTENT(INOUT) :: upsoiu  (npoi,nsoilay)! local  ! soil water uptake from transpiration (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: upsoil  (npoi,nsoilay)! local  ! soil water uptake from transpiration (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: heatg   (npoi)         ! local  ! net heat flux into soil surface (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: heati   (npoi)         ! local  ! net heat flux into snow surface (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: hydraul (npoi,nsoilay)! global ! saturated hydraulic conductivity (m/s)
    REAL(KIND=r8), INTENT(INOUT) :: porosflo(npoi,nsoilay)! global ! porosity after reduction by ice content
    INTEGER      , INTENT(IN   ) :: ibex    (npoi,nsoilay)! global ! nint(bex), used for cpu speed
    REAL(KIND=r8), INTENT(IN   ) :: bperm   (npoi )        ! global ! lower b.c. for soil profile drainage 
    ! (0.0 = impermeable; 1.0 = fully permeable)
    REAL(KIND=r8), INTENT(INOUT) :: hflo    (npoi,nsoilay+1)  ! downward heat transport through soil layers (W m-2)


    !   INCLUDE 'comatm.h'
    REAL(KIND=r8), INTENT(IN   ) :: ta     (npoi)         ! global ! air temperature (K)
    REAL(KIND=r8), INTENT(INOUT) :: asurd  (npoi,nband)   ! local  ! direct albedo of surface system
    REAL(KIND=r8), INTENT(INOUT) :: asuri  (npoi,nband)   ! local  ! diffuse albedo of surface system 
    REAL(KIND=r8), INTENT(IN   ) :: coszen (npoi)         ! global ! cosine of solar zenith angle
    REAL(KIND=r8), INTENT(IN   ) :: solad  (npoi,nband)   ! global ! direct downward solar flux (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: solai  (npoi,nband)   ! global ! diffuse downward solar flux (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: fira   (npoi)         ! global ! incoming ir flux (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: raina  (npoi)         ! global ! rainfall rate (mm/s or kg m-2 s-1)
    REAL(KIND=r8), INTENT(IN   ) :: qa     (npoi)         ! global ! specific humidity (kg_h2o/kg_air)
    REAL(KIND=r8), INTENT(IN   ) :: psurf  (npoi)         ! global ! surface pressure (Pa)
    REAL(KIND=r8), INTENT(IN   ) :: snowa  (npoi)         ! global ! snowfall rate (mm/s or kg m-2 s-1 of water)
    REAL(KIND=r8), INTENT(IN   ) :: ua     (npoi)         ! global ! wind speed (m s-1)
    REAL(KIND=r8), INTENT(IN   ) :: o2conc                 ! global ! o2 concentration (mol/mol)
    REAL(KIND=r8), INTENT(INOUT) :: co2conc(npoi)                 ! global ! co2 concentration (mol/mol)
    REAL(KIND=r8), INTENT(INOUT) :: z0(npoi)
    REAL(KIND=r8), INTENT(INOUT) :: ustar(npoi)
    REAL(KIND=r8), INTENT(OUT) :: hc (npoi)
    REAL(KIND=r8), INTENT(OUT) :: hg (npoi)
    REAL(KIND=r8), INTENT(OUT) :: ec (npoi)
    REAL(KIND=r8), INTENT(OUT) :: eg (npoi)
    REAL(KIND=r8), INTENT(INOUT) :: dispu    (npoi)          ! local ! zero-plane displacement height for upper canopy (m)
    REAL(KIND=r8), INTENT(INOUT) :: cu       (npoi)          ! local ! air transfer coefficient (*rhoa) (m s-1 kg m-3) for
    !         upper air region (z12 --> za) (A35 Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(INOUT) :: firb     (npoi)          ! local ! net upward ir radiation at reference
    ! atmospheric level za (W m-2)
    REAL(KIND=r8), INTENT(OUT) :: bstar(npoi)

   REAL(KIND=r8), INTENT(INOUT) :: ynleach_p (npoi) ! annual total amount P leached from soil profile   (kg_P m-2/yr)
   REAL(KIND=r8), INTENT(OUT  ) :: tnmin_p   (npoi)   ! instantaneous phosphorus mineralization         (kg_P m-2/timestep)
   REAL(KIND=r8), INTENT(OUT  ) :: totnmic_p (npoi)   ! total phosphorus residing in microbial pool     (kg_P m-2)
   REAL(KIND=r8), INTENT(OUT  ) :: totnlit_p (npoi)   ! total phosphorus in all litter pools            (kg_P m-2)
   REAL(KIND=r8), INTENT(OUT  ) :: totanlit_p(npoi)   ! total standing aboveground phosphorus in litter (kg_P m-2)
   REAL(KIND=r8), INTENT(OUT  ) :: totrnlit_p(npoi)   ! total root litter phosphorus belowground        (kg_P m-2)
   REAL(KIND=r8), INTENT(OUT  ) :: totnsoi_p (npoi)   ! total phosphorus in soil                        (kg_P m-2)
   REAL(KIND=r8), INTENT(INOUT) :: storedn_p (npoi)   ! total storage of P in soil profile              (kg_P m-2) 

    !   INCLUDE 'com1d.h'
    REAL(KIND=r8) :: fwetu    (npoi)          ! local ! fraction of upper canopy leaf area wetted by intercepted liquid and/or snow
    REAL(KIND=r8) :: rliqu    (npoi)          ! local ! proportion of fwetu due to liquid
    REAL(KIND=r8) :: fwets    (npoi)          ! local ! fraction of upper canopy stem area wetted by intercepted liquid and/or snow
    REAL(KIND=r8) :: rliqs    (npoi)          ! local ! proportion of fwets due to liquid
    REAL(KIND=r8) :: fwetl    (npoi)          ! local ! fraction of lower canopy stem & leaf area wetted by
    !         intercepted liquid and/or snow
    REAL(KIND=r8) :: rliql    (npoi)          ! local ! proportion of fwetl due to liquid
    REAL(KIND=r8) :: solu     (npoi)          ! local ! solar flux (direct + diffuse) absorbed by upper 
    !         canopy leaves per unit canopy area (W m-2)
    REAL(KIND=r8) :: sols     (npoi)          ! local ! solar flux (direct + diffuse) absorbed by upper 
    !         canopy stems per unit canopy area (W m-2)
    REAL(KIND=r8) :: soll     (npoi)          ! local ! solar flux (direct + diffuse) absorbed by lower 
    !         canopy leaves and stems per unit canopy area (W m-2)
    REAL(KIND=r8) :: solg     (npoi)          ! local ! solar flux (direct + diffuse) absorbed by unit 
    !         snow-free soil (W m-2)
    REAL(KIND=r8) :: soli     (npoi)          ! local ! solar flux (direct + diffuse) absorbed by unit 
    ! snow surface (W m-2)
    REAL(KIND=r8) :: scalcoefl(npoi,4)     ! local ! term needed in lower canopy scaling
    REAL(KIND=r8) :: scalcoefu(npoi,4)     ! local ! term needed in upper canopy scaling
    INTEGER :: indsol   (npoi)                  ! local ! index of current strip for points with positive coszen
    REAL(KIND=r8) :: albsod   (npoi)          ! local ! direct  albedo for soil surface (visible or IR)
    REAL(KIND=r8) :: albsoi   (npoi)          ! local ! diffuse albedo for soil surface (visible or IR)
    REAL(KIND=r8) :: albsnd   (npoi)          ! local ! direct  albedo for snow surface (visible or IR)
    REAL(KIND=r8) :: albsni   (npoi)          ! local ! diffuse albedo for snow surface (visible or IR)
    REAL(KIND=r8) :: relod    (npoi)          ! local ! upward direct radiation per unit icident direct beam on lower canopy (W m-2)
    REAL(KIND=r8) :: reloi    (npoi)          ! local ! upward diffuse radiation per unit incident diffuse 
    ! radiation on lower canopy (W m-2)
    REAL(KIND=r8) :: reupd    (npoi)          ! local ! upward direct radiation per unit incident direct 
    ! radiation on upper canopy (W m-2)
    REAL(KIND=r8) :: reupi    (npoi)          ! local ! upward diffuse radiation per unit incident diffuse 
    ! radiation on upper canopy (W m-2)
    REAL(KIND=r8) :: ablod    (npoi)          ! local ! fraction of direct  radiation absorbed by lower canopy
    REAL(KIND=r8) :: abloi    (npoi)          ! local ! fraction of diffuse radiation absorbed by lower canopy
    REAL(KIND=r8) :: flodd    (npoi)          ! local ! downward direct radiation per unit incident direct
    ! radiation on lower canopy (W m-2)
    REAL(KIND=r8) :: dummy    (npoi)          ! local ! placeholder, always = 0: no direct flux produced for diffuse incident
    REAL(KIND=r8) :: flodi    (npoi)          ! local ! downward diffuse radiation per unit incident direct
    ! radiation on lower canopy (W m-2)
    REAL(KIND=r8) :: floii    (npoi)          ! local ! downward diffuse radiation per unit incident 
    ! diffuse radiation on lower canopy
    REAL(KIND=r8) :: terml    (npoi,7)     ! local ! term needed in lower canopy scaling
    REAL(KIND=r8) :: termu    (npoi,7)     ! local ! term needed in upper canopy scaling
    REAL(KIND=r8) :: abupd    (npoi)          ! local ! fraction of direct  radiation absorbed by upper canopy
    REAL(KIND=r8) :: abupi    (npoi)          ! local ! fraction of diffuse radiation absorbed by upper canopy
    REAL(KIND=r8) :: fupdd    (npoi)          ! local ! downward direct radiation per unit incident direct
    ! beam on upper canopy (W m-2)
    REAL(KIND=r8) :: fupdi    (npoi)          ! local ! downward diffuse radiation per unit icident direct
    ! radiation on upper canopy (W m-2)
    REAL(KIND=r8) :: fupii    (npoi)          ! local ! downward diffuse radiation per unit incident diffuse
    ! radiation on upper canopy (W m-2)
    REAL(KIND=r8) :: sol2d    (npoi)          ! local ! direct downward radiation  out of upper canopy 
    ! per unit vegetated (upper) area (W m-2)
    REAL(KIND=r8) :: sol2i    (npoi)          ! local ! diffuse downward radiation out of upper
    ! canopy per unit vegetated (upper) area(W m-2)
    REAL(KIND=r8) :: sol3d    (npoi)          ! local ! direct downward radiation  out of upper
    ! canopy + gaps per unit grid cell area (W m-2)
    REAL(KIND=r8) :: sol3i    (npoi)          ! local ! diffuse downward radiation out of upper
    ! canopy + gaps per unit grid cell area (W m-2)
    REAL(KIND=r8) :: firs     (npoi)          ! local ! ir radiation absorbed by upper canopy stems (W m-2)
    REAL(KIND=r8) :: firu     (npoi)          ! local ! ir raditaion absorbed by upper canopy leaves (W m-2)
    REAL(KIND=r8) :: firl     (npoi)          ! local ! ir radiation absorbed by lower canopy leaves and stems (W m-2)
    REAL(KIND=r8) :: firg     (npoi)          ! local ! ir radiation absorbed by soil/ice (W m-2)
    REAL(KIND=r8) :: firi     (npoi)          ! local ! ir radiation absorbed by snow (W m-2)
    REAL(KIND=r8) :: snowg    (npoi)          ! local ! snowfall rate at soil level (kg h2o m-2 s-1)
    REAL(KIND=r8) :: tsnowg   (npoi)          ! local ! snowfall temperature at soil level (K) 
    REAL(KIND=r8) :: tsnowl   (npoi)          ! local ! snowfall temperature below upper canopy (K)
    REAL(KIND=r8) :: pfluxl   (npoi)          ! local ! heat flux on lower canopy leaves & stems due to intercepted h2o (W m-2)
    REAL(KIND=r8) :: raing    (npoi)          ! local ! rainfall rate at soil level (kg m-2 s-1)
    REAL(KIND=r8) :: traing   (npoi)          ! local ! rainfall temperature at soil level (K)
    REAL(KIND=r8) :: trainl   (npoi)          ! local ! rainfall temperature below upper canopy (K)
    REAL(KIND=r8) :: snowl    (npoi)          ! local ! snowfall rate below upper canopy (kg h2o m-2 s-1)
    REAL(KIND=r8) :: tsnowu   (npoi)          ! local ! snowfall temperature above upper canopy (K)
    REAL(KIND=r8) :: pfluxu   (npoi)          ! local ! heat flux on upper canopy leaves due to intercepted h2o (W m-2)
    REAL(KIND=r8) :: rainu    (npoi)          ! local ! rainfall rate above upper canopy (kg m-2 s-1)
    REAL(KIND=r8) :: trainu   (npoi)          ! local ! rainfall temperature above upper canopy (K)
    REAL(KIND=r8) :: snowu    (npoi)          ! local ! snowfall rate above upper canopy (kg h2o m-2 s-1)
    REAL(KIND=r8) :: pfluxs   (npoi)          ! local ! heat flux on upper canopy stems due to intercepted h2o (W m-2)
    REAL(KIND=r8) :: rainl    (npoi)          ! local ! rainfall rate below upper canopy (kg m-2 s-1)
!    REAL(KIND=r8) :: bps                  ! local ! (ps/p) ** (rair/cair) for atmospheric level  (const)
    REAL(KIND=r8) :: cp       (npoi)          ! local ! specific heat of air at za (allowing for h2o vapor) (J kg-1 K-1)
    REAL(KIND=r8) :: bdl      (npoi)          ! local ! aerodynamic coefficient ([(tau/rho)/u**2] for
    ! laower canopy (A31/A30 Pollard & Thompson 1995)
    REAL(KIND=r8) :: dil      (npoi)          ! local ! inverse of momentum diffusion coefficient within lower canopy (m)
    REAL(KIND=r8) :: z3       (npoi)          ! local ! effective top of the lower canopy (for momentum) (m)
    REAL(KIND=r8) :: z4       (npoi)          ! local ! effective bottom of the lower canopy (for momentum) (m)
    REAL(KIND=r8) :: z34      (npoi)          ! local ! effective middle of the lower canopy (for momentum) (m)
    REAL(KIND=r8) :: exphl    (npoi)          ! local ! exp(lamda/2*(z3-z4)) for lower canopy (A30 Pollard & Thompson)
    REAL(KIND=r8) :: expl     (npoi)          ! local ! exphl**2
    REAL(KIND=r8) :: displ    (npoi)          ! local ! zero-plane displacement height for lower canopy (m)
    REAL(KIND=r8) :: bdu      (npoi)          ! local ! aerodynamic coefficient ([(tau/rho)/u**2] for upper
    ! canopy (A31/A30 Pollard & Thompson 1995)
    REAL(KIND=r8) :: diu      (npoi)          ! local ! inverse of momentum diffusion coefficient within upper canopy (m)
    REAL(KIND=r8) :: z1       (npoi)          ! local ! effective top of upper canopy (for momentum) (m)
    REAL(KIND=r8) :: z2       (npoi)          ! local ! effective bottom of the upper canopy (for momentum) (m)
    REAL(KIND=r8) :: z12      (npoi)          ! local ! effective middle of the upper canopy (for momentum) (m)
    REAL(KIND=r8) :: exphu    (npoi)          ! local ! exp(lamda/2*(z3-z4)) for upper canopy (A30 Pollard & Thompson)
    REAL(KIND=r8) :: expu     (npoi)          ! local ! exphu**2
    REAL(KIND=r8) :: alogg    (npoi)          ! local ! log of soil roughness
    REAL(KIND=r8) :: alogi    (npoi)          ! local ! log of snow roughness
    REAL(KIND=r8) :: alogav   (npoi)          ! local ! average of alogi and alogg 
    REAL(KIND=r8) :: alog4    (npoi)          ! local ! log (max(z4, 1.1*z0sno, 1.1*z0soi)) 
    REAL(KIND=r8) :: alog3    (npoi)          ! local ! log (z3 - displ)
    REAL(KIND=r8) :: alog2    (npoi)          ! local ! log (z2 - displ)
    REAL(KIND=r8) :: alog1    (npoi)          ! local ! log (z1 - dispu) 
    REAL(KIND=r8) :: aloga    (npoi)          ! local ! log (za - dispu) 
    REAL(KIND=r8) :: u2       (npoi)          ! local ! wind speed at level z2 (m s-1)
    REAL(KIND=r8) :: alogu    (npoi)          ! local ! log (roughness length of upper canopy)
    REAL(KIND=r8) :: alogl    (npoi)          ! local ! log (roughness length of lower canopy)
    REAL(KIND=r8) :: richl    (npoi)          ! local ! richardson number for air above upper canopy (z3 to z2)
    REAL(KIND=r8) :: straml   (npoi)          ! local ! momentum correction factor for stratif between
    ! upper & lower canopy (z3 to z2) (louis et al.)
    REAL(KIND=r8) :: strahl   (npoi)          ! local ! heat/vap correction factor for stratif between
    !         upper & lower canopy (z3 to z2) (louis et al.)
    REAL(KIND=r8) :: richu    (npoi)          ! local ! richardson number for air between upper & lower canopy (z1 to za)
    REAL(KIND=r8) :: stramu   (npoi)          ! local ! momentum correction factor for stratif above
    !         upper canopy (z1 to za) (louis et al.)
    REAL(KIND=r8) :: strahu   (npoi)          ! local ! heat/vap correction factor for stratif above
    !         upper canopy (z1 to za) (louis et al.)
    REAL(KIND=r8) :: u1       (npoi)          ! local ! wind speed at level z1 (m s-1)
    REAL(KIND=r8) :: u12      (npoi)          ! local ! wind speed at level z12 (m s-1)
    REAL(KIND=r8) :: u3       (npoi)          ! local ! wind speed at level z3 (m s-1)
    REAL(KIND=r8) :: u34      (npoi)          ! local ! wind speed at level z34 (m s-1)
    REAL(KIND=r8) :: u4       (npoi)          ! local ! wind speed at level z4 (m s-1)
    REAL(KIND=r8) :: cl       (npoi)          ! local ! air transfer coefficient (*rhoa) (m s-1 kg m-3)
    !         between the 2 canopies (z34 --> z12) (A36 Pollard & Thompson 1995)
    REAL(KIND=r8) :: sg       (npoi)          ! local ! air-soil transfer coefficient
    REAL(KIND=r8) :: si       (npoi)          ! local ! air-snow transfer coefficient
    REAL(KIND=r8) :: fwetux   (npoi)          ! local ! fraction of upper canopy leaf area wetted if dew forms
    REAL(KIND=r8) :: fwetsx   (npoi)          ! local ! fraction of upper canopy stem area wetted if dew forms
    REAL(KIND=r8) :: fwetlx   (npoi)          ! local ! fraction of lower canopy leaf and stem area wetted if dew forms
    REAL(KIND=r8) :: fsena    (npoi)          ! local ! downward sensible heat flux between za & z12 at za (W m-2)
    REAL(KIND=r8) :: fseng    (npoi)          ! local ! upward sensible heat flux between soil surface & air at z34 (W m-2)
    REAL(KIND=r8) :: fseni    (npoi)          ! local ! upward sensible heat flux between snow surface & air at z34 (W m-2)
    REAL(KIND=r8) :: fsenu    (npoi)          ! local ! sensible heat flux from upper canopy leaves to air (W m-2)
    REAL(KIND=r8) :: fsens    (npoi)          ! local ! sensible heat flux from upper canopy stems to air (W m-2)
    REAL(KIND=r8) :: fsenl    (npoi)          ! local ! sensible heat flux from lower canopy to air (W m-2)
    REAL(KIND=r8) :: fvapa    (npoi)          ! local ! downward h2o vapor flux between za & z12 at za (kg m-2 s-1)
    REAL(KIND=r8) :: fvaput   (npoi)          ! local ! h2o vapor flux (transpiration from dry parts) 
    ! between upper canopy leaves and air at z12 (kg m-2 s-1/ LAI upper canopy/ fu)
    REAL(KIND=r8) :: fvaps    (npoi)       ! local ! h2o vapor flux (evaporation from wet surface)
    !            between upper canopy stems and air at z12 (kg m-2 s-1 / SAI lower canopy / fu)
    REAL(KIND=r8) :: fvaplw   (npoi)       ! local ! h2o vapor flux (evaporation from wet surface) 
    !            between lower canopy leaves & stems and air at z34 (kg m-2 s-1/ LAI lower canopy/ fl)
    REAL(KIND=r8) :: fvaplt   (npoi)       ! local ! h2o vapor flux (transpiration) 
    !            between lower canopy & air at z34 (kg m-2 s-1 / LAI lower canopy / fl)
    REAL(KIND=r8) :: fvapg    (npoi)       ! local ! h2o vapor flux (evaporation) between soil & air 
    !         at z34 (kg m-2 s-1/bare ground fraction)
    REAL(KIND=r8) :: fvapi    (npoi)       ! local ! h2o vapor flux (evaporation) between snow & air at z34 (kg m-2 s-1 / fi )
    REAL(KIND=r8) :: fvapuw   (npoi)       ! local ! h2o vapor flux (evaporation from wet parts)
    ! between upper canopy leaves and air at z12 (kg m-2 s-1/ LAI upper canopy/ fu)
    REAL(KIND=r8), INTENT(INOUT) :: td (npoi)   ! global! daily average temperature (K)
    REAL(KIND=r8), INTENT(IN   ) :: vzero(npoi) ! global! a real array of zeros, of length npoi

    INTEGER, INTENT(IN   ) :: ndaypy               ! global! number of days per year
    REAL(KIND=r8), INTENT(OUT  ) :: nppdummy (npoi,npft)! local ! canopy NPP before accounting for stem and root respiration
    REAL(KIND=r8) :: tgpp     (npoi,npft)                 ! local ! instantaneous GPP for each pft (mol-CO2 / m-2 / second)
    REAL(KIND=r8), INTENT(OUT  ) :: tgpptot  (npoi)                         ! local ! instantaneous gpp (mol-CO2 / m-2 / second)
    REAL(KIND=r8) :: tnpp     (npoi,npft)                 ! local ! instantaneous NPP for each pft (mol-CO2 / m-2 / second)
    REAL(KIND=r8), INTENT(INOUT) :: cbiow    (npoi,npft)  ! global! carbon in woody biomass pool (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: sapfrac  (npoi)         ! global! fraction of woody biomass that is in sapwood
    REAL(KIND=r8), INTENT(INOUT) :: cbior    (npoi,npft)  ! global! carbon in fine root biomass pool (kg_C m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: tnpptot  (npoi)                         ! local ! instantaneous npp (mol-CO2 / m-2 / second)
    REAL(KIND=r8), INTENT(INOUT) :: tco2root (npoi)         ! local ! instantaneous fine co2 flux from soil (mol-CO2 / m-2 / second)
    REAL(KIND=r8), INTENT(OUT  ) :: tneetot  (npoi)         ! local ! instantaneous net ecosystem exchange of co2 per timestep(mol-CO2 / m-2 / second)
    REAL(KIND=r8), INTENT(INOUT) :: tco2mic  (npoi)         ! local ! instantaneous microbial co2 flux from soil (mol-CO2 / m-2 / second)
    REAL(KIND=r8), INTENT(INOUT) :: a10td    (npoi)       ! global! 10-day average daily air temperature (K)
    REAL(KIND=r8), INTENT(INOUT) :: a10ancub (npoi)       ! global! 10-day average canopy photosynthesis rate - broadleaf (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: a10ancuc (npoi)       ! global! 10-day average canopy photosynthesis rate - conifer (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: a10ancls (npoi)       ! global! 10-day average canopy photosynthesis rate - shrubs (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: a10ancl3 (npoi)       ! global! 10-day average canopy photosynthesis rate - c3 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: a10ancl4 (npoi)       ! global! 10-day average canopy photosynthesis rate - c4 grasses (mol_co2 m-2 s-1)

    INTEGER, INTENT(INOUT) :: ndtimes        (npoi)! global! counter for daily average calculations
    REAL(KIND=r8), INTENT(INOUT) :: adrain    (npoi)! global! daily average rainfall rate (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: adsnow    (npoi)! global! daily average snowfall rate (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: adaet     (npoi)! global! daily average aet (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: adtrunoff (npoi)! global! daily average total runoff (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: adsrunoff (npoi)! global! daily average surface runoff (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: addrainage(npoi)! global! daily average drainage (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: adrh      (npoi)! global! daily average rh (percent)
    REAL(KIND=r8), INTENT(INOUT) :: adsnod    (npoi)! global! daily average snow depth (m)
    REAL(KIND=r8), INTENT(INOUT) :: adsnof    (npoi)! global! daily average snow fraction (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: adwsoi    (npoi)! global! daily average soil moisture (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: adtsoi    (npoi)! global! daily average soil temperature (c)
    REAL(KIND=r8), INTENT(INOUT) :: adwisoi   (npoi)! global! daily average soil ice (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: adtlaysoi (npoi)! global! daily average soil temperature (c) of top layer
    REAL(KIND=r8), INTENT(INOUT) :: adwlaysoi (npoi)! global! daily average soil moisture of top layer(fraction)
    REAL(KIND=r8), INTENT(INOUT) :: adwsoic   (npoi)! global! daily average soil moisture using root profile weighting (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: adtsoic   (npoi)! global! daily average soil temperature (c) using profile weighting
    REAL(KIND=r8), INTENT(INOUT) :: adco2mic  (npoi)! global! daily accumulated co2 respiration from microbes (kg_C m-2 /day)
    REAL(KIND=r8), INTENT(INOUT) :: adco2root (npoi)! global! daily accumulated co2 respiration from roots (kg_C m-2 /day)
    REAL(KIND=r8), INTENT(INOUT) :: adco2soi  (npoi)! global! daily accumulated co2 respiration from soil(total) (kg_C m-2 /day)
    REAL(KIND=r8), INTENT(INOUT) :: adco2ratio(npoi)! global! ratio of root to total co2 respiration
    REAL(KIND=r8), INTENT(INOUT) :: adnmintot (npoi)! global! daily accumulated net nitrogen mineralization (kg_N m-2 /day)
    REAL(KIND=r8), INTENT(INOUT) :: decompl   (npoi)! global! litter decomposition factor              (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: decomps   (npoi)! global! soil organic matter decomposition factor      (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: tnmin     (npoi)! global! instantaneous nitrogen mineralization (kg_N m-2/timestep)

    INTEGER, INTENT(IN   ) :: ndaypm          (12)          ! global! number of days per month

    INTEGER, INTENT(INOUT) :: nmtimes        (npoi)                   ! global! counter for monthly average calculations
    REAL(KIND=r8), INTENT(INOUT) :: amrain        (npoi)     ! global! monthly average rainfall rate (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: amsnow        (npoi)     ! global! monthly average snowfall rate (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: amaet        (npoi)     ! global! monthly average aet (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: amtrunoff    (npoi)     ! global! monthly average total runoff (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: amsrunoff    (npoi)     ! global! monthly average surface runoff (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: amdrainage   (npoi)     ! global! monthly average drainage (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: amtemp        (npoi)     ! global! monthly average air temperature (C)
    REAL(KIND=r8), INTENT(INOUT) :: amqa         (npoi)     ! global! monthly average specific humidity (kg-h2o/kg-air)
    REAL(KIND=r8), INTENT(INOUT) :: amsolar        (npoi)     ! global! monthly average incident solar radiation (W/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: amirup        (npoi)     ! global! monthly average upward ir radiation (W/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: amirdown     (npoi)     ! global! monthly average downward ir radiation (W/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: amsens        (npoi)     ! global! monthly average sensible heat flux (W/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: amlatent     (npoi)     ! global! monthly average latent heat flux (W/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: amlaiu        (npoi)     ! global! monthly average lai for upper canopy (m**2/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: amlail        (npoi)     ! global! monthly average lai for lower canopy (m**2/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: amtsoi        (npoi)     ! global! monthly average 1m soil temperature (C)
    REAL(KIND=r8), INTENT(INOUT) :: amwsoi        (npoi)     ! global! monthly average 1m soil moisture (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: amwisoi        (npoi)     ! global! monthly average 1m soil ice (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: amvwc        (npoi)     ! global! monthly average 1m volumetric water content (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: amawc        (npoi)     ! global! monthly average 1m plant-available water content (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: amsnod        (npoi)     ! global! monthly average snow depth (m)
    REAL(KIND=r8), INTENT(INOUT) :: amsnof        (npoi)     ! global! monthly average snow fraction (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: amnpp        (npoi,npft)! global! monthly total npp for each plant type (kg-C/m**2/month)
    REAL(KIND=r8), INTENT(INOUT) :: adnpp        (npoi,npft)! global! monthly total npp for each plant type (kg-C/m**2/day)
    REAL(KIND=r8), INTENT(OUT  ) :: amnpptot     (npoi)     ! local ! monthly total npp for ecosystem (kg-C/m**2/month)
    REAL(KIND=r8), INTENT(INOUT) :: amco2mic     (npoi)     ! global! monthly total CO2 flux from microbial respiration (kg-C/m**2/month)
    REAL(KIND=r8), INTENT(INOUT) :: amco2root    (npoi)     ! global! monthly total CO2 flux from soil due to root
    ! respiration (kg-C/m**2/month)
    REAL(KIND=r8), INTENT(OUT  ) :: amco2soi     (npoi)     ! local ! monthly total soil CO2 flux from microbial
    ! and root respiration (kg-C/m**2/month)
    REAL(KIND=r8), INTENT(OUT  ) :: amco2ratio   (npoi)       ! local ! monthly ratio of root to total co2 flux
    REAL(KIND=r8), INTENT(OUT  ) :: amneetot     (npoi)       ! local ! monthly total net ecosystem exchange of CO2 (kg-C/m**2/month)
    REAL(KIND=r8), INTENT(INOUT) :: amnmintot    (npoi)       ! global! monthly total N mineralization from microbes (kg-N/m**2/month)
    REAL(KIND=r8), INTENT(INOUT) :: amts2        (npoi)       ! global
    REAL(KIND=r8), INTENT(INOUT) :: amtransu     (npoi)       ! global
    REAL(KIND=r8), INTENT(INOUT) :: amtransl     (npoi)       ! global
    REAL(KIND=r8), INTENT(INOUT) :: amsuvap      (npoi)       ! global
    REAL(KIND=r8), INTENT(INOUT) :: aminvap      (npoi)       ! global
    REAL(KIND=r8), INTENT(INOUT) :: amalbedo     (npoi)         
    REAL(KIND=r8), INTENT(INOUT) :: amtsoil    (npoi, nsoilay) 
    REAL(KIND=r8), INTENT(INOUT) :: amwsoil    (npoi, nsoilay) 
    REAL(KIND=r8), INTENT(INOUT) :: amwisoil   (npoi, nsoilay)

    INTEGER, INTENT(INOUT) :: nytimes (npoi)                ! global! counter for yearly average calculations
    REAL(KIND=r8), INTENT(INOUT) :: aysolar        (npoi)     ! global! annual average incident solar radiation (w/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: ayirup        (npoi)     ! global! annual average upward ir radiation (w/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: ayirdown   (npoi)       ! global! annual average downward ir radiation (w/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: aysens        (npoi)     ! global! annual average sensible heat flux (w/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: aylatent   (npoi)       ! global! annual average latent heat flux (w/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: ayprcp        (npoi)     ! global! annual average precipitation (mm/yr)
    REAL(KIND=r8), INTENT(INOUT) :: ayaet        (npoi)     ! global! annual average aet (mm/yr)
    REAL(KIND=r8), INTENT(INOUT) :: aytrans        (npoi)     ! global! annual average transpiration (mm/yr)
    REAL(KIND=r8), INTENT(INOUT) :: aytrunoff  (npoi)       ! global! annual average total runoff (mm/yr)
    REAL(KIND=r8), INTENT(INOUT) :: aysrunoff  (npoi)       ! global! annual average surface runoff (mm/yr)
    REAL(KIND=r8), INTENT(INOUT) :: aydrainage (npoi)       ! global! annual average drainage (mm/yr)
    REAL(KIND=r8), INTENT(INOUT) :: aydwtot        (npoi)     ! global! annual average soil+vegetation+snow water 
    ! recharge (mm/yr or kg_h2o/m**2/yr)
    REAL(KIND=r8), INTENT(INOUT) :: aywsoi        (npoi)     ! global! annual average 1m soil moisture (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: aywisoi        (npoi)     ! global! annual average 1m soil ice (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: aytsoi        (npoi)     ! global! annual average 1m soil temperature (C)
    REAL(KIND=r8), INTENT(INOUT) :: ayvwc        (npoi)     ! global! annual average 1m volumetric water content (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: ayawc        (npoi)     ! global! annual average 1m plant-available water content (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: aystresstu (npoi)       ! global! annual average soil moisture stress 
    ! parameter for upper canopy (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: aystresstl(npoi)        ! global! annual average soil moisture stress 
    ! parameter for lower canopy (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: aygpp     (npoi,npft)   ! global! annual gross npp for each plant type(kg-c/m**2/yr)
    REAL(KIND=r8), INTENT(OUT  ) :: aygpptot  (npoi)        ! local ! annual total gpp for ecosystem (kg-c/m**2/yr)
    REAL(KIND=r8), INTENT(INOUT) :: aynpp     (npoi,npft)   ! global! annual total npp for each plant type(kg-c/m**2/yr)
    REAL(KIND=r8), INTENT(OUT  ) :: aynpptot  (npoi)        ! local ! annual total npp for ecosystem (kg-c/m**2/yr)
    REAL(KIND=r8), INTENT(INOUT) :: ayco2mic  (npoi)        ! global! annual total CO2 flux from microbial respiration (kg-C/m**2/yr)
    REAL(KIND=r8), INTENT(INOUT) :: ayco2root (npoi)        ! global! annual total CO2 flux from soil due to root respiration (kg-C/m**2/yr)
    REAL(KIND=r8), INTENT(OUT  ) :: ayco2soi  (npoi)        ! local ! annual total soil CO2 flux from microbial and 
    ! root respiration (kg-C/m**2/yr)
    REAL(KIND=r8), INTENT(OUT  ) :: ayneetot  (npoi)        ! local ! annual total NEE for ecosystem (kg-C/m**2/yr)
    REAL(KIND=r8), INTENT(INOUT) :: ayrootbio (npoi)        ! global! annual average live root biomass (kg-C / m**2)
    REAL(KIND=r8), INTENT(INOUT) :: aynmintot (npoi)        ! global! annual total nitrogen mineralization (kg-N/m**2/yr)
    REAL(KIND=r8), INTENT(INOUT) :: ayalit    (npoi)        ! global! aboveground litter (kg-c/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: ayblit    (npoi)        ! global! belowground litter (kg-c/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: aycsoi    (npoi)        ! global! total soil carbon (kg-c/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: aycmic    (npoi)        ! global! total soil carbon in microbial biomass (kg-c/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: ayanlit   (npoi)        ! global! aboveground litter nitrogen (kg-N/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: aybnlit   (npoi)        ! global! belowground litter nitrogen (kg-N/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: aynsoi    (npoi)        ! global! total soil nitrogen (kg-N/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: ayalbedo  (npoi)
    REAL(KIND=r8), INTENT(INOUT) :: totalit   (npoi)           ! global! total standing aboveground litter (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: totrlit   (npoi)           ! global! total root litter carbon belowground (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: totcsoi   (npoi)           ! global! total carbon in all soil pools (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: totcmic   (npoi)           ! global! total carbon residing in microbial pools (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: totanlit  (npoi)           ! global! total standing aboveground nitrogen in litter (kg_N m-2)
    REAL(KIND=r8), INTENT(INOUT) :: totrnlit  (npoi)        ! global! total root litter nitrogen belowground (kg_N m-2)
    REAL(KIND=r8), INTENT(INOUT) :: totnsoi   (npoi)        ! global! total nitrogen in soil (kg_N m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: totnmic   (npoi)        ! local! total nitrogen residing in microbial pool (kg_N m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: totlit    (npoi)        ! local! total carbon in all litter pools (kg_C m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: totfall   (npoi)        ! local! total litterfall and root turnover (kg_C m-2/year)
    REAL(KIND=r8), INTENT(OUT  ) :: totnlit   (npoi)        ! local! total nitrogen in all litter pools (kg_N m-2)

    REAL(KIND=r8), INTENT(INOUT) :: firefac   (npoi)        ! global! factor that respresents the annual average
    REAL(KIND=r8), INTENT(INOUT) :: wtot      (npoi)        ! global! total amount of water stored in snow, soil,
    ! puddels, and on vegetation (kg_h2o)
    ! fuel dryness of a grid cell, and hence characterizes the readiness to burn

    REAL(KIND=r8), INTENT(INOUT) :: storedn (npoi)           ! global! total storage of N in soil profile (kg_N m-2) 
    REAL(KIND=r8), INTENT(INOUT) :: yrleach (npoi)           ! global! annual total amount C leached from soil profile (kg_C m-2/yr)
    REAL(KIND=r8), INTENT(INOUT) :: ynleach (npoi)
    REAL(KIND=r8), INTENT(INOUT) :: falll   (npoi)          ! global ! annual leaf litter fall (kg_C m-2/year)
    REAL(KIND=r8), INTENT(INOUT) :: fallr   (npoi)          ! global ! annual root litter input                    (kg_C m-2/year)
    REAL(KIND=r8), INTENT(INOUT) :: fallw   (npoi)          ! global ! annual wood litter fall                    (kg_C m-2/year)
    REAL(KIND=r8), INTENT(INOUT) :: clitlm  (npoi)          ! global! carbon in leaf litter pool - metabolic       (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitls  (npoi)          ! global! carbon in leaf litter pool - structural      (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitrm  (npoi)          ! global! carbon in fine root litter pool - metabolic  (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitrs  (npoi)          ! global! carbon in fine root litter pool - structural (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitwm  (npoi)          ! global! carbon in woody litter pool - metabolic      (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitws  (npoi)          ! global! carbon in woody litter pool - structural     (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: csoislop(npoi)          ! global! carbon in soil - slow protected humus           (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: csoislon(npoi)          ! global! carbon in soil - slow nonprotected humus     (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: csoipas (npoi)          ! global! carbon in soil - passive humus                   (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitll  (npoi)          ! global! carbon in leaf litter pool - lignin           (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitrl  (npoi)          ! global! carbon in fine root litter pool - lignin     (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitwl  (npoi)          ! global! carbon in woody litter pool - lignin           (kg_C m-2)


    REAL(KIND=r8), INTENT(INOUT) :: tc      (npoi)          ! global  ! coldest monthly temperature (C)
    REAL(KIND=r8), INTENT(INOUT) :: agddu   (npoi)          ! global  ! annual accumulated growing degree days for bud
    ! burst, upper canopy (day-degrees)
    REAL(KIND=r8), INTENT(INOUT) :: tempu   (npoi)          ! global  ! cold-phenology trigger for trees (non-dimensional)
    REAL(KIND=r8), INTENT(INOUT) :: agddl   (npoi)          ! global  ! annual accumulated growing degree days for bud burst,
    ! lower canopy (day-degrees)
    REAL(KIND=r8), INTENT(INOUT) :: templ   (npoi)           ! global  ! cold-phenology trigger for grasses/shrubs (non-dimensional)
    REAL(KIND=r8), INTENT(INOUT) :: dropu   (npoi)          ! global  ! drought-phenology trigger for trees (non-dimensional)
    REAL(KIND=r8), INTENT(INOUT) :: dropls  (npoi)          ! global  ! drought-phenology trigger for shrubs (non-dimensional)
    REAL(KIND=r8), INTENT(INOUT) :: dropl4  (npoi)          ! global  ! drought-phenology trigger for c4 grasses (non-dimensional)
    REAL(KIND=r8), INTENT(INOUT) :: dropl3  (npoi)          ! global  ! drought-phenology trigger for c3 grasses (non-dimensional)
    REAL(KIND=r8), INTENT(INOUT) :: plai    (npoi,npft)     ! global  ! total leaf area index of each plant functional type
    !
    ! Arguments (input)
    !
    INTEGER, INTENT(IN   ) :: iday         ! day number  (passed in)
    INTEGER, INTENT(IN   ) :: imonth         ! month number (passed in)
    INTEGER, INTENT(IN   ) :: iyear
    INTEGER, INTENT(IN   ) :: iyear0
    INTEGER, INTENT(IN   ) :: isimveg 
    INTEGER, INTENT(IN   ) :: spinmax 
    REAL(KIND=r8), INTENT(IN        ) :: ux  (npoi)
    REAL(KIND=r8), INTENT(IN        ) :: uy  (npoi)
    REAL(KIND=r8), INTENT(OUT  ) :: taux(npoi)
    REAL(KIND=r8), INTENT(OUT  ) :: tauy(npoi)
    REAL(KIND=r8), INTENT(INOUT) :: ts2 (npoi)
    REAL(KIND=r8), INTENT(INOUT) :: qs2 (npoi)
    REAL(KIND=r8), INTENT(IN        ) :: deltat   (npoi)      ! absolute minimum temperature -
    ! temp on average of coldest month (C)
    REAL(KIND=r8), INTENT(INOUT) :: gdd0     (npoi)         ! growing degree days > 0C 
    REAL(KIND=r8), INTENT(INOUT) :: gdd0this (npoi)         ! annual total growing degree days for current year
    REAL(KIND=r8), INTENT(INOUT) :: tcthis   (npoi)      ! coldest monthly temperature of current year (C)
    REAL(KIND=r8), INTENT(INOUT) :: twthis   (npoi)      ! warmest monthly temperature of current year (C)
    REAL(KIND=r8), INTENT(INOUT) :: tcmin    (npoi)      ! coldest daily temperature of current year (C)
    REAL(KIND=r8), INTENT(INOUT) :: gdd5     (npoi)      ! growing degree days > 5C
    REAL(KIND=r8), INTENT(INOUT) :: gdd5this (npoi)      ! annual total growing degree days for current year
    REAL(KIND=r8), INTENT(IN   ) :: TminL    (npft)      ! Absolute minimum temperature -- lower limit (upper canopy PFTs)
    REAL(KIND=r8), INTENT(IN   ) :: TminU    (npft)      ! Absolute minimum temperature -- upper limit (upper canopy PFTs)
    REAL(KIND=r8), INTENT(IN   ) :: Twarm    (npft)      ! Temperature of warmest month (lower canopy PFTs)
    REAL(KIND=r8), INTENT(IN   ) :: GDD      (npft)      ! minimum GDD needed (base 5 C for upper canopy PFTs, 
    REAL(KIND=r8), INTENT(IN   ) :: aleaf    (npft)          ! carbon allocation fraction to leaves
    REAL(KIND=r8), INTENT(IN   ) :: awood    (npft)          ! carbon allocation fraction to wood 
    REAL(KIND=r8), INTENT(INOUT) :: cbiol    (npoi,npft) ! carbon in leaf biomass pool (kg_C m-2)
    REAL(KIND=r8), INTENT(IN   ) :: aroot    (npft)         ! carbon allocation fraction to fine roots
    REAL(KIND=r8), INTENT(OUT  ) :: disturbf (npoi)         ! annual fire disturbance regime (m2/m2/yr)
    REAL(KIND=r8), INTENT(OUT  ) :: disturbo (npoi)         ! fraction of biomass pool lost every year to disturbances other than fire
    REAL(KIND=r8), INTENT(IN   ) :: specla   (npft)         ! specific leaf area (m**2/kg) 
    REAL(KIND=r8), INTENT(OUT  ) :: biomass  (npoi,npft) ! total biomass of each plant functional type  (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: totlaiu  (npoi)         ! total leaf area index for the upper canopy
    REAL(KIND=r8), INTENT(INOUT) :: totlail  (npoi)         ! total leaf area index for the lower canopy
    REAL(KIND=r8), INTENT(INOUT) :: totbiou  (npoi)         ! total biomass in the upper canopy (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: totbiol  (npoi)         ! total biomass in the lower canopy (kg_C m-2)
    REAL(KIND=r8), INTENT(IN   ) :: woodnorm                   ! value of woody biomass for upper canopy closure
    ! (ie when wood = woodnorm fu = 1.0) (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: vegtype0 (npoi)      ! annual vegetation type - ibis classification
    REAL(KIND=r8), INTENT(IN   ) :: tauwood0 (npft)      ! normal (unstressed) turnover time for wood biomass (years)
    REAL(KIND=r8), INTENT(OUT  ) :: tauwood  (npoi,npft)      ! wood biomass turnover time constant (years)
    REAL(KIND=r8), INTENT(IN   ) :: tauleaf  (npft)      ! foliar biomass turnover time constant (years)
    REAL(KIND=r8), INTENT(IN   ) :: tauroot  (npft)      ! fine root biomass turnover time constant (years)
    REAL(KIND=r8), INTENT(IN   ) :: xminlai                 ! Minimum LAI for each existing PFT
    REAL(KIND=r8), INTENT(OUT  ) :: cdisturb (npoi)         ! annual amount of vegetation carbon lost 
    ! to atmosphere due to fire  (biomass burning) (kg_C m-2/year)
    REAL(KIND=r8), INTENT(OUT  ) :: ayanpp   (npoi,npft)   ! annual above-ground npp for each plant type(kg-c/m**2/yr)
    REAL(KIND=r8), INTENT(OUT  ) :: ayanpptot(npoi)             ! annual above-ground npp for ecosystem (kg-c/m**2/yr)
    !REAL(KIND=r8), INTENT(IN   ) :: garea    (npoi)   ! area of each gridcell (m**2)
    REAL(KIND=r8), INTENT(INOUT) :: tw      (npoi)      ! warmest monthly temperature (C)
    REAL(KIND=r8), INTENT(INOUT) :: adfalll   (npoi)          ! global ! annual leaf litter fall (kg_C m-2/year)
    REAL(KIND=r8), INTENT(INOUT) :: adfallr   (npoi)          ! global ! annual root litter input    (kg_C m-2/year)
    REAL(KIND=r8), INTENT(INOUT) :: adfallw   (npoi)          ! global ! annual wood litter fall    (kg_C m-2/year)
    REAL(KIND=r8), INTENT(INOUT) :: adcbiol(npoi,npft)
    REAL(KIND=r8), INTENT(INOUT) :: adcbior(npoi,npft)
    REAL(KIND=r8), INTENT(INOUT) :: adcbiow(npoi,npft)
    REAL(KIND=r8), INTENT(INOUT) :: adplai    (npoi,npft)     ! global  ! total leaf area index of each plant functional type
   CHARACTER(len=*) , INTENT(IN   ) ::  rootmode
    !
    ! Arguments (input)
    !
    INTEGER, INTENT(IN   ) :: nstep      ! atm time step index
    !INTEGER, INTENT(IN   ) :: idayprev   ! day in month of previous timestep
    INTEGER, INTENT(IN   ) :: imonthprev ! month of previous timestep 
    !INTEGER, INTENT(IN   ) :: iyearprev  ! year of previous timestep 
    !INTEGER, INTENT(IN   ) :: idayout         ! write out daily output
    !INTEGER, INTENT(IN   ) :: imonthout  ! write out daily output
    !INTEGER, INTENT(IN   ) :: iyearout
    INTEGER, INTENT(IN   ) :: isimco2
    INTEGER, INTENT(IN   ) :: isimfire
    INTEGER, PARAMETER     :: lenMonth(12)=(/31,28,31,30,31,30,31,31,30,31,30,31/) 
    REAL(KIND=r8), INTENT(IN   ) :: co2initm(npoi)  !mol/mol
    REAL(KIND=r8), INTENT(IN   ) :: beta1(nVegClass)
    REAL(KIND=r8), INTENT(IN   ) :: beta2(nVegClass)
    REAL(KIND=r8), INTENT(IN   ) :: stressfac(nVegClass)
    REAL(KIND=r8), INTENT(IN   ) :: avmuir_factor(nVegClass,2)

    !REAL(KIND=r8), INTENT(IN   ) :: calday             ! current julian day (1-365.99)
    !  decimals=fraction of day

    !  INTEGER :: iday         ! current day in month (1-31 or 1-30, or 1-28)
    !   integer:: imonth         ! current month (1 - 12)
    !INTEGER :: mcdate
    !REAL(KIND=r8) :: co2vmrgcm            ! modified co2 volume mixing ratio
    !REAL(KIND=r8) :: dtibis             ! time step as passed from GCM, dtime from
    INTEGER :: lenMonthly
    INTEGER :: yi,mi,di,hi,LenYearbyDay,nday2y
    REAL(KIND=r8) :: ftod,xday2
    INTEGER :: j
    INTEGER :: i
    INTEGER :: ndyn
    terml=0.0_r8
    termu=0.0_r8
    abupd=0.0_r8
    abupi=0.0_r8
    fupdd=0.0_r8
    !
    !
      IF (mcsec == 0.0_r8) THEN
           CALL DailyDynaVeg(isimfire  ,&
                              npoi      ,&!  INTEGER      , INTENT(IN   ) :: npoi                 ! total number of land points
                              npft     ,&!  INTEGER      , INTENT(IN   ) :: npft                 ! number of plant functional types
                              woodnorm ,&!  REAL(KIND=r8), INTENT(IN   ) :: woodnorm                ! value of woody biomass for upper canopy closure! (ie when wood = woodnorm fu = 1.0) (kg_C m-2)
                              xminlai  ,&!  REAL(KIND=r8), INTENT(IN   ) :: xminlai                ! Minimum LAI for each existing PFT
                              specla   ,&!  REAL(KIND=r8), INTENT(IN   ) :: specla    (npft)         ! specific leaf area (m**2/kg) 
                              aleaf    ,&!  REAL(KIND=r8), INTENT(IN   ) :: aleaf     (npft)         ! carbon allocation fraction to leaves
                              awood    ,&!  REAL(KIND=r8), INTENT(IN   ) :: awood     (npft)         ! carbon allocation fraction to wood 
                              tauwood0 ,&!  REAL(KIND=r8), INTENT(IN   ) :: tauwood0  (npft)         ! normal (unstressed) turnover time for wood biomass (years)
                              tauleaf  ,&!  REAL(KIND=r8), INTENT(IN   ) :: tauleaf   (npft)         ! foliar biomass turnover time constant (years)
                              tauroot  ,&!  REAL(KIND=r8), INTENT(IN   ) :: tauroot   (npft)         ! fine root biomass turnover time constant (years)
                              aroot    ,&!  REAL(KIND=r8), INTENT(IN   ) :: aroot     (npft)         ! carbon allocation fraction to fine roots
                              exist    ,&!  REAL(KIND=r8), INTENT(IN   ) :: exist     (npoi,npft)  ! probability of existence of each plant functional type in a gridcell
                              adco2mic ,&!  REAL(KIND=r8), INTENT(IN   ) :: adco2mic  (npoi)         ! global! daily accumulated co2 respiration from microbes (kg_C m-2 /day)
                              adnpp    ,&!  REAL(KIND=r8), INTENT(INOUT) :: adnpp     (npoi,npft)  ! annual total npp for each plant type(kg-c/m**2/day)
                              adcbiol  ,&!  REAL(KIND=r8), INTENT(INOUT) :: adcbiol   (npoi,npft)  ! carbon in leaf biomass pool (kg_C m-2)
                              adcbior  ,&!  REAL(KIND=r8), INTENT(INOUT) :: adcbior   (npoi,npft)  ! carbon in fine root biomass pool (kg_C m-2)
                              adcbiow  ,&!  REAL(KIND=r8), INTENT(INOUT) :: adcbiow   (npoi,npft)  ! carbon in woody biomass pool (kg_C m-2)
                              adplai   ,&!  REAL(KIND=r8), INTENT(INOUT) :: adplai    (npoi,npft)  ! global  ! total leaf area index of each plant functional type
                              adfalll  ,&!  REAL(KIND=r8), INTENT(OUT  ) :: adfalll   (npoi)         ! global ! annual leaf litter fall (kg_C m-2/day)
                              adfallr  ,&!  REAL(KIND=r8), INTENT(OUT  ) :: adfallr   (npoi)         ! global ! annual root litter input(kg_C m-2/day)
                              adfallw  ,&!  REAL(KIND=r8), INTENT(OUT  ) :: adfallw   (npoi)         ! global ! annual wood litter fall (kg_C m-2/day)
                              fu       ,&!  REAL(KIND=r8), INTENT(OUT  ) :: fu        (npoi) ! fraction of overall area covered by upper canopy
                              fl       ,&!  REAL(KIND=r8), INTENT(OUT  ) :: fl        (npoi) ! fraction of snow-free area covered by lower  canopy
                              zbot     ,&!  REAL(KIND=r8), INTENT(OUT  ) :: zbot      (npoi,2) ! height of lowest branches above ground (m)
                              ztop     ,&!  REAL(KIND=r8), INTENT(OUT  ) :: ztop      (npoi,2) ! height of plant top above ground (m)
                              sai      ,&!  REAL(KIND=r8), INTENT(OUT  ) :: sai       (npoi,2) ! current single-sided stem area index
                              sapfrac  ,&!  REAL(KIND=r8), INTENT(OUT  ) :: sapfrac   (npoi) ! fraction of woody biomass that is in 
                              vegtype0 ,&!REAL(KIND=r8), INTENT(INOUT) :: vegtype0 (npoi)      ! annual vegetation type - ibis classification
                              totlit   ,&
                              firefac  ,& 
                              totlaiu  ,&
                              totlail  ,&
                              totbiou  ,&
                              totbiol   )

     
      END IF

!
!
    CALL Ibis    (mcsec         ,pi          ,stef         ,vonk         ,grav        , &
         tmelt        ,hfus        ,hvap         ,hsub         ,ch2o        , &
         cice         ,cair        ,cvap         ,rair         ,rvap        , &
         cappa        ,rhow        ,npoi         ,nband        ,nsoilay     , & 
         nsnolay      ,npft        ,epsilon      ,dtime        ,doalb       , &
         ginvap       ,gsuvap      ,gtrans       ,gtransu      ,gtransl     , &
         grunof       ,gdrain      ,gadjust      ,a10scalparamu,a10daylightu, &
         a10scalparaml,a10daylightl,vmax_pft     ,tau15        ,kc15        , &
         ko15         ,cimax       ,gammaub      ,alpha3       ,theta3      , &
         beta3        ,coefmub     ,coefbub      ,gsubmin      ,gammauc     , &
         coefmuc      ,coefbuc     ,gsucmin      ,gammals      ,coefmls     , & 
         coefbls      ,gslsmin     ,gammal3      ,coefml3      ,coefbl3     , &
         gsl3min      ,gammal4     ,alpha4       ,theta4       ,beta4       , &
         coefml4      ,coefbl4     ,gsl4min      ,wliqu        ,wliqumax    , & 
         wsnou        ,wsnoumax    ,tu           ,wliqs        ,wliqsmax    , & 
         wsnos        ,wsnosmax    ,ts           ,wliql        ,wliqlmax    , &  
         wsnol        ,wsnolmax    ,tl           ,topparu      ,topparl     , &
         fl           ,fu          ,lai          ,sai          ,rhoveg      , &   
         tauveg       ,orieh       ,oriev        ,wliqmin      ,wsnomin     , & 
         t12          ,tdripu      ,tblowu       ,tdrips       ,tblows      , &
         t34          ,tdripl      ,tblowl       ,ztop         ,alaiml      , &
         zbot         ,alaimu      ,froot        ,q34          ,q12         , &
         su           ,cleaf       ,dleaf        ,ss           ,cstem       , & 
         dstem        ,sl          ,cgrass       ,ciub         ,ciuc        , &
         exist        ,csub        ,gsub         ,csuc         ,gsuc        , &
         agcub        ,agcuc       ,ancub        ,ancuc        ,totcondub   , &
         totconduc    ,cils        ,cil3         ,cil4         ,csls        , &
         gsls         ,csl3        ,gsl3         ,csl4         ,gsl4        , &
         agcls        ,agcl4       ,agcl3        ,ancls        ,ancl4       , &
         ancl3        ,totcondls   ,totcondl3    ,totcondl4    ,chu         , &
         chs          ,chl         ,frac         ,tlsub        ,z0sno       , & 
         rhos         ,consno      ,hsnotop      ,hsnomin      ,fimin       , &
         fimax        ,fi          ,tsno         ,hsno         ,sand        , &
         clay         ,poros       ,wsoi         ,wisoi        ,consoi      , &  
         zwpmax       ,wpud        ,wipud        ,wpudmax      ,qglif       , &         
         tsoi         ,hvasug      ,hvasui       ,albsav       ,albsan      , &
         tg           ,ti          ,z0soi        ,swilt        ,sfield      , &
         stressl      ,stressu     ,stresstl     ,stresstu     ,csoi        , &         
         rhosoi       ,hsoi        ,suction      ,bex          ,upsoiu      , &  
         upsoil       ,heatg       ,heati        ,hydraul      ,porosflo    , &
         ibex         ,bperm       ,hflo         ,ta           ,asurd       , &  
         asuri        ,coszen      ,solad        ,solai        ,fira        , & 
         raina        ,qa          ,psurf        ,snowa        ,ua          , &   
         o2conc       ,co2conc     ,fwetu        ,rliqu        ,fwets       , &
         rliqs        ,fwetl       ,rliql        ,solu         , &!PK
         sols         ,soll        ,solg         ,soli         ,scalcoefl   , &
         scalcoefu    ,indsol      ,albsod       ,albsoi       ,albsnd      , &
         albsni       ,relod       ,reloi        ,reupd        ,reupi       , &
         ablod        ,abloi       ,flodd        ,dummy        ,flodi       , &
         floii         ,terml      ,termu        ,abupd        ,abupi       , & 
         fupdd        ,fupdi       ,fupii        ,sol2d        ,sol2i       , & 
         sol3d        ,sol3i       ,firb         ,firs         ,firu        , &
         firl         ,firg        ,firi         ,snowg        ,tsnowg      , &    
         tsnowl       ,pfluxl      ,raing        ,traing       ,trainl      , &
         snowl        ,tsnowu      ,pfluxu       ,rainu        ,trainu      , &   
         snowu        ,pfluxs      ,rainl        ,bps          ,cp          , &  
         za           ,bdl         ,dil          ,z3           ,z4          , & 
         z34          ,exphl       ,expl         ,displ        ,bdu         , &
         diu          ,z1          ,z2           ,z12          ,exphu       , &
         expu         ,dispu       ,alogg        ,alogi        ,alogav      , & 
         alog4        ,alog3       ,alog2        ,alog1        ,aloga       , &          
         u2           ,alogu       ,alogl        ,richl        ,straml      , &
         strahl       ,richu       ,stramu       ,strahu       ,u1          , & 
         u12          ,u3          ,u34          ,u4           ,cu          , &
         cl           ,sg          ,si           ,fwetux       ,fwetsx      , &  
         fwetlx       ,fsena       ,fseng        ,fseni        ,fsenu       , &
         fsens        ,fsenl       ,fvapa        ,fvaput       ,fvaps       , &
         fvaplw       ,fvaplt      ,fvapg        ,fvapi        ,fvapuw      , &
         td           ,vzero       ,ndaypy       ,nppdummy     ,tgpp        , & 
         tgpptot      ,tnpp        ,cbiow        ,sapfrac      ,cbior       , & 
         tnpptot      ,tco2root    ,tneetot      ,tco2mic      ,a10td       , &
         a10ancub     ,a10ancuc    ,a10ancls     ,a10ancl3     ,a10ancl4    , & 
         ndtimes      ,adrain      ,adsnow      , &
         adaet        ,adtrunoff   ,adsrunoff    ,addrainage   ,adrh        , &
         adsnod       ,adsnof      ,adwsoi       ,adtsoi       ,adwisoi     , &
         adtlaysoi    ,adwlaysoi   ,adwsoic      ,adtsoic      ,adco2mic    , &
         adco2root    ,adco2soi    ,adco2ratio   ,adnmintot    ,decompl     , &    
         decomps      ,tnmin       ,ndaypm       ,nmtimes     , & 
         amrain       ,amsnow      ,amaet        ,amtrunoff    ,amsrunoff   , &
         amdrainage   ,amtemp      ,amqa         , &
         amsolar      ,amirup      ,amirdown     ,amsens       ,amlatent    , &  
         amlaiu       ,amlail      ,amtsoi       ,amwsoi       ,amwisoi     , &   
         amvwc        ,amawc       ,amsnod       ,amsnof       ,amnpp       , &
         amnpptot     ,amco2mic    ,amco2root    ,amco2soi     ,amco2ratio  , &
         amneetot     ,amnmintot   ,nytimes      ,aysolar      ,ayirup      , &
         ayirdown     ,aysens      ,aylatent     ,ayprcp       ,ayaet       , &  
         aytrans      ,aytrunoff   ,aysrunoff    ,aydrainage   ,aydwtot     , & 
         aywsoi       ,aywisoi     ,aytsoi       ,ayvwc        ,ayawc       , &  
         aystresstu   ,aystresstl  ,aygpp        ,aygpptot     ,aynpp       , & 
         aynpptot     ,ayco2mic    ,ayco2root    ,ayco2soi     ,ayneetot    , &
         ayrootbio    ,aynmintot   ,ayalit       ,ayblit       ,aycsoi      , & 
         aycmic       ,ayanlit     ,aybnlit      ,aynsoi       ,ayalbedo    , &
         beta1        ,beta2       ,stressfac    ,avmuir_factor,totalit     , &
         totrlit      ,totcsoi     ,totcmic      ,totanlit     ,totrnlit    , &
         totnsoi      ,totnmic     ,totlit       ,totfall      ,totnlit     , &
         firefac      ,wtot        ,storedn      ,yrleach      ,ynleach     , & 
         falll        ,fallr       ,fallw        ,clitlm       ,clitls      , &
         clitrm       ,clitrs      ,clitwm       ,clitws       ,csoislop    , &
         csoislon     ,csoipas     ,clitll       ,clitrl       ,clitwl      , &  
         tc           ,agddu       ,tempu        ,agddl        ,templ       , &
         dropu        ,dropls      ,dropl4       ,dropl3       ,plai        , &
         iday         ,imonth      ,iyear        ,iyear0      , &
         isimveg      ,spinmax     ,amts2        ,amtransu    ,amtransl     , &
         amsuvap      ,aminvap     ,amalbedo     ,amtsoil     ,amwsoil      , & 
         amwisoil     ,ux          ,uy           ,taux        ,tauy         , &
         ts2          ,qs2         ,gdd0this     ,gdd5this    ,bstar        , &
         vegtype0     ,ynleach_p    ,tnmin_p     ,totnmic_p    ,totnlit_p   , &
         totanlit_p   ,totrnlit_p  ,totnsoi_p    ,storedn_p   ,adnpp       , &
         adfalll      ,adfallr     ,adfallw      ,adcbiol     ,adcbior     ,&
         adcbiow      ,adplai      ,nstep        ,nVegClass   ,rootmode)

    DO i = 1, npoi
       z0    (i) =  EXP(((alogu(i))))
       ustar (i) =  ua (i) * cu (i)
       !ustar(i) =  ua (i) * (cu (i)/rhoa(i))

       hc    (i) =  (fsena  (i) + fsenu    (i) + fsens     (i) + fsenl    (i) +fsenl(i) )*dtime
       hg    (i) =  (fseng  (i) + fseni    (i)                                          )*dtime

       !hltm   = 2.52e6_r8            !  latent heat of vaporization (J kg^-1)

       ec    (i) =  (fvaput (i) + fvaps    (i) + fvaplw    (i)  + fvaplt (i) +fvapuw (i))*hltm*dtime
       eg    (i) =  (fvapg  (i) + fvapi    (i)                                                 )*hltm*dtime
    END DO


    IF (nstep > 20) THEN
       ! IF (nstep /= 0) THEN
       !        
       lenMonthly=lenMonth(imonth)
       yi=iyear
       mi=imonth
       di=iday
       hi=0
       ftod=0.0_r8
       CALL jull(yi,mi,di,hi,ftod,xday2,nday2y,LenYearbyDay)
       IF(LenYearbyDay == 366 .AND. imonth == 2) lenMonthly=29
       IF ((mcsec == 86400.0_r8-(dtime/2.0_r8)) .AND. (iday == lenMonthly) ) THEN 
          !
          !     end of calculations done once a month
          !
          DO i = 1, npoi
             tcthis(i) = MIN (tcthis(i), (amts2(i) - 273.160_r8))
             twthis(i) = MAX (twthis(i), (amts2(i) - 273.160_r8))
          END DO
       END IF
       !
       ! calculations done once a year
       !
       IF ((nstep > 360) .AND.(mcsec == 86400.0_r8-(dtime/2.0_r8)) .AND. (iday == lenMonthly) .AND. (imonth == 12) ) THEN
          !         IF ((mcsec == 0.0_r8) .and. (iday == 1) .and. (imonth == 1) ) THEN

          !
          ! get new o2 and co2 concentrations for this year
          !
          IF (isimco2.EQ.1) CALL co2 (npoi,co2initm,  &! INTENT(IN   )
               co2conc,  &! INTENT(OUT  )
               iyear     )! INTENT(IN   )

          !
          ! perform vegetation dynamics
          !
          ndyn = 1
          DO j = 1, ndyn

             IF (isimveg /= 0) CALL dynaveg2 (isimfire , &! INTENT(IN   )        dynaveg1 (isimfire , &! INTENT(IN   )
                  tauwood0  , &! INTENT(IN   )                  tauwood0 , &! INTENT(IN   )
                  tauwood   , &! INTENT(OUT  )                  tauwood  , &! INTENT(OUT  )
                  tauleaf   , &! INTENT(IN   )                  tauleaf  , &! INTENT(IN   )
                  tauroot   , &! INTENT(IN   )                  tauroot  , &! INTENT(IN   )
                  xminlai   , &! INTENT(IN   )                  xminlai  , &! INTENT(IN   )
                  falll     , &! INTENT(OUT  )                  falll    , &! INTENT(OUT  )
                  fallr     , &! INTENT(OUT  )                  fallr    , &! INTENT(OUT  )
                  fallw     , &! INTENT(OUT  )                  fallw    , &! INTENT(OUT  )
                  cdisturb  , &! INTENT(OUT  )                  cdisturb , &! INTENT(OUT  )
                  exist     , &! INTENT(OUT  )                  exist    , &! INTENT(IN   )
                  aleaf     , &! INTENT(IN   )                  aleaf    , &! INTENT(IN   )
                  awood     , &! INTENT(IN   )                  awood    , &! INTENT(IN   )
                  cbiol     , &! INTENT(INOUT) global           cbiol    , &! INTENT(INOUT) global
                  cbior     , &! INTENT(INOUT) global           cbior    , &! INTENT(INOUT) global
                  cbiow     , &! INTENT(INOUT) global           cbiow    , &! INTENT(INOUT) global
                  aroot     , &! INTENT(IN   )                  aroot    , &! INTENT(IN   )
                  disturbf  , &! INTENT(OUT  )                  disturbf , &! INTENT(OUT  )
                  disturbo  , &! INTENT(OUT  )                  disturbo , &! INTENT(OUT  )
                  firefac   , &! INTENT(IN   )                  firefac  , &! INTENT(IN   )
                  totlit    , &! INTENT(IN   )                  totlit   , &! INTENT(IN   )
                  specla    , &! INTENT(IN   )                  specla   , &! INTENT(IN   )
                  plai      , &! INTENT(INOUT) local                 plai          , &! INTENT(INOUT) local
                  biomass   , &! INTENT(OUT  )                  biomass  , &! INTENT(OUT  )
                  totlaiu   , &! INTENT(INOUT) local                 totlaiu  , &! INTENT(INOUT) local
                  totlail   , &! INTENT(INOUT) local                 totlail  , &! INTENT(INOUT) local
                  totbiou   , &! INTENT(INOUT) local                 totbiou  , &! INTENT(INOUT) local
                  totbiol   , &! INTENT(OUT  )                  totbiol  , &! INTENT(OUT  )
                  fu             , &! INTENT(OUT  )                  fu          , &! INTENT(OUT  )
                  woodnorm  , &! INTENT(IN   )                  woodnorm , &! INTENT(IN   )
                  fl             , &! INTENT(OUT  )                  fl          , &! INTENT(OUT  )
                  zbot      , &! INTENT(OUT  )                  zbot          , &! INTENT(OUT  )
                  ztop      , &! INTENT(OUT  )                  ztop          , &! INTENT(OUT  )
                  sai       , &! INTENT(OUT  )                  sai          , &! INTENT(OUT  )
                  sapfrac   , &! INTENT(OUT  )                  sapfrac  , &! INTENT(OUT  )
                  vegtype0  , &! INTENT(OUT  )                  vegtype0 , &! INTENT(OUT  )
                  gdd5      , &! INTENT(IN   )                  gdd5          , &! INTENT(IN   )
                  gdd0      , &! INTENT(IN   )                  gdd0          , &! INTENT(IN   )
                  aynpp     , &! INTENT(INOUT) global           aynpp    , &! INTENT(INOUT) global
                  ayanpp    , &! INTENT(OUT  )                  ayanpp   , &! INTENT(OUT  )
                  ayneetot  , &! INTENT(INOUT) global           ayneetot , &! INTENT(INOUT) global
                  ayanpptot , &! INTENT(OUT  )                  ayanpptot, &! INTENT(OUT  )
                  aynpptot  , &! INTENT(OUT  )                  npoi          , &!
                  ayco2mic  , &! INTENT(IN   )                  npft            )! , isim_ac, year)
                  npoi      , &!
                  npft        )! , isim_ac, year)!

             !IF (isimveg /= 0) CALL dynaveg1 (isimfire , &! INTENT(IN        )
             !                             tauwood0 , &! INTENT(IN   )
             !                             tauwood  , &! INTENT(OUT  )
             !                             tauleaf  , &! INTENT(IN   )
             !                             tauroot  , &! INTENT(IN   )
             !                             xminlai  , &! INTENT(IN   )
             !                             falll    , &! INTENT(OUT  )
             !                             fallr    , &! INTENT(OUT  )
             !                             fallw    , &! INTENT(OUT  )
             !                             cdisturb , &! INTENT(OUT  )
             !                             exist    , &! INTENT(IN   )
             !                             aleaf    , &! INTENT(IN   )
             !                             awood    , &! INTENT(IN   )
             !                             cbiol    , &! INTENT(INOUT) global
             !                             cbior    , &! INTENT(INOUT) global
             !                             cbiow    , &! INTENT(INOUT) global
             !                             aroot    , &! INTENT(IN   )
             !                             disturbf , &! INTENT(OUT  )
             !                             disturbo , &! INTENT(OUT  )
             !                             firefac  , &! INTENT(IN   )
             !                             totlit   , &! INTENT(IN   )
             !                             specla   , &! INTENT(IN   )
             !                             plai     , &! INTENT(INOUT) local
             !                             biomass  , &! INTENT(OUT  )
             !                             totlaiu  , &! INTENT(INOUT) local
             !                             totlail  , &! INTENT(INOUT) local
             !                             totbiou  , &! INTENT(INOUT) local
             !                             totbiol  , &! INTENT(OUT  )
             !                             fu       , &! INTENT(OUT  )
             !                             woodnorm , &! INTENT(IN   )
             !                             fl       , &! INTENT(OUT  )
             !                             zbot     , &! INTENT(OUT  )
             !                             ztop     , &! INTENT(OUT  )
             !                             sai      , &! INTENT(OUT  )
             !                             sapfrac  , &! INTENT(OUT  )
             !                             vegtype0 , &! INTENT(OUT  )
             !                             gdd5     , &! INTENT(IN   )
             !                             gdd0     , &! INTENT(IN   )
             !                             aynpp    , &! INTENT(INOUT) global
             !                             ayanpp   , &! INTENT(OUT  )
             !                             ayneetot , &! INTENT(INOUT) global
             !                             ayanpptot, &! INTENT(OUT  )
             !                             npoi     , &!
             !                             npft        )! , isim_ac, year)
             !
             !
          END DO

          !
          !
          !     recalculate bioclimatic parameters (used in dynaveg, calculated
          !     even in fixed vegetation case when fixed vegetation is an
          !     initialisation of a dynamic run)
          !
          CALL climanl2(TminL    , &! INTENT(IN   )
               TminU    , &! INTENT(IN   )
               Twarm    , &! INTENT(IN   )
               GDD      , &! INTENT(IN   )
               gdd0     , &! INTENT(INOUT)
               gdd0this , &! INTENT(IN   )
               tc       , &! INTENT(INOUT)
               tw       , &! INTENT(INOUT)
               tcthis   , &! INTENT(IN   )
               twthis   , &! INTENT(IN   )
               tcmin    , &! INTENT(INOUT) local
               gdd5     , &! INTENT(INOUT) local
               gdd5this , &! INTENT(IN   )
               exist    , &! INTENT(INOUT)
               deltat   , &! INTENT(IN   )
               npoi     , &! INTENT(IN   )
               npft       )! INTENT(IN   )

       END IF

       IF (imonthprev .NE. imonth) THEN
          !
          ! write restart files
          !
          !CALL wrestart (mdcur, imonthprev, iyearprev, iyear0)
       END IF
       !
       !     End of test on 1st time step
       !
    END IF
  END SUBROUTINE IbisDrv






  SUBROUTINE Ibis   (mcsec        ,pi          ,stef         ,vonk         ,grav        , &
       tmelt        ,hfus        ,hvap         ,hsub         ,ch2o        , &
       cice         ,cair        ,cvap         ,rair         ,rvap        , &
       cappa        ,rhow        ,npoi         ,nband        ,nsoilay     , & 
       nsnolay      ,npft        ,epsilon      ,dtime        ,doalb       , &
       ginvap       ,gsuvap      ,gtrans       ,gtransu      ,gtransl     , &
       grunof       ,gdrain      ,gadjust      ,a10scalparamu,a10daylightu, &
       a10scalparaml,a10daylightl,vmax_pft     ,tau15        ,kc15        , &
       ko15         ,cimax       ,gammaub      ,alpha3       ,theta3      , &
       beta3        ,coefmub     ,coefbub      ,gsubmin      ,gammauc     , &
       coefmuc      ,coefbuc     ,gsucmin      ,gammals      ,coefmls     , & 
       coefbls      ,gslsmin     ,gammal3      ,coefml3      ,coefbl3     , &
       gsl3min      ,gammal4     ,alpha4       ,theta4       ,beta4       , &
       coefml4      ,coefbl4     ,gsl4min      ,wliqu        ,wliqumax    , & 
       wsnou        ,wsnoumax    ,tu           ,wliqs        ,wliqsmax    , & 
       wsnos        ,wsnosmax    ,ts           ,wliql        ,wliqlmax    , &  
       wsnol        ,wsnolmax    ,tl           ,topparu      ,topparl     , &
       fl           ,fu          ,lai          ,sai          ,rhoveg      , &   
       tauveg       ,orieh       ,oriev        ,wliqmin      ,wsnomin     , & 
       t12          ,tdripu      ,tblowu       ,tdrips       ,tblows      , &
       t34          ,tdripl      ,tblowl       ,ztop         ,alaiml      , &
       zbot         ,alaimu      ,froot        ,q34          ,q12         , &
       su           ,cleaf       ,dleaf        ,ss           ,cstem       , & 
       dstem        ,sl          ,cgrass       ,ciub         ,ciuc        , &
       exist        ,csub        ,gsub         ,csuc         ,gsuc        , &
       agcub        ,agcuc       ,ancub        ,ancuc        ,totcondub   , &
       totconduc    ,cils        ,cil3         ,cil4         ,csls        , &
       gsls         ,csl3        ,gsl3         ,csl4         ,gsl4        , &
       agcls        ,agcl4       ,agcl3        ,ancls        ,ancl4       , &
       ancl3       ,totcondls    ,totcondl3    ,totcondl4    ,chu         , &
       chs          ,chl         ,frac         ,tlsub        ,z0sno       , & 
       rhos         ,consno      ,hsnotop      ,hsnomin      ,fimin       , &
       fimax        ,fi          ,tsno         ,hsno         ,sand        , &
       clay         ,poros       ,wsoi         ,wisoi        ,consoi      , &  
       zwpmax       ,wpud        ,wipud        ,wpudmax      ,qglif       , &         
       tsoi         ,hvasug      ,hvasui       ,albsav       ,albsan      , &
       tg           ,ti          ,z0soi        ,swilt        ,sfield      , &
       stressl      ,stressu     ,stresstl     ,stresstu     ,csoi        , &         
       rhosoi       ,hsoi        ,suction      ,bex          ,upsoiu      , &  
       upsoil       ,heatg       ,heati        ,hydraul      ,porosflo    , &
       ibex         ,bperm       ,hflo         ,ta           ,asurd       , &  
       asuri        ,coszen      ,solad        ,solai        ,fira        , & 
       raina        ,qa          ,psurf        ,snowa        ,ua          , &   
       o2conc       ,co2conc     ,fwetu        ,rliqu        ,fwets       , &
       rliqs        ,fwetl       ,rliql        ,solu         ,&         
       sols         ,soll        ,solg         ,soli         ,scalcoefl   , &
       scalcoefu    ,indsol      ,albsod       ,albsoi       ,albsnd      , &
       albsni       ,relod       ,reloi        ,reupd        ,reupi       , &
       ablod        ,abloi       ,flodd        ,dummy        ,flodi       , &
       floii        ,terml       ,termu        ,abupd        ,abupi       , & 
       fupdd        ,fupdi       ,fupii        ,sol2d        ,sol2i       , & 
       sol3d        ,sol3i       ,firb         ,firs         ,firu        , &
       firl         ,firg        ,firi         ,snowg        ,tsnowg      , &    
       tsnowl       ,pfluxl      ,raing        ,traing       ,trainl      , &
       snowl        ,tsnowu      ,pfluxu       ,rainu        ,trainu      , &   
       snowu        ,pfluxs      ,rainl        ,bps          ,cp          , &  
       za           ,bdl         ,dil          ,z3           ,z4          , & 
       z34          ,exphl       ,expl         ,displ        ,bdu         , &
       diu          ,z1          ,z2           ,z12          ,exphu       , &
       expu         ,dispu       ,alogg        ,alogi        ,alogav      , & 
       alog4        ,alog3       ,alog2        ,alog1        ,aloga       , &          
       u2           ,alogu       ,alogl        ,richl        ,straml      , &
       strahl       ,richu       ,stramu       ,strahu       ,u1          , & 
       u12          ,u3          ,u34          ,u4           ,cu          , &
       cl           ,sg          ,si           ,fwetux       ,fwetsx      , &  
       fwetlx       ,fsena       ,fseng        ,fseni        ,fsenu       , &
       fsens        ,fsenl       ,fvapa        ,fvaput       ,fvaps       , &
       fvaplw       ,fvaplt      ,fvapg        ,fvapi        ,fvapuw      , &
       td           ,vzero       ,ndaypy       ,nppdummy     ,tgpp        , & 
       tgpptot      ,tnpp        ,cbiow        ,sapfrac      ,cbior       , & 
       tnpptot      ,tco2root    ,tneetot      ,tco2mic      ,a10td       , &
       a10ancub     ,a10ancuc    ,a10ancls     ,a10ancl3     ,a10ancl4    , & 
       ndtimes      ,adrain       ,adsnow      , &
       adaet        ,adtrunoff   ,adsrunoff    ,addrainage   ,adrh        , &
       adsnod       ,adsnof      ,adwsoi       ,adtsoi       ,adwisoi     , &
       adtlaysoi    ,adwlaysoi   ,adwsoic      ,adtsoic      ,adco2mic    , &
       adco2root    ,adco2soi    ,adco2ratio   ,adnmintot    ,decompl     , &    
       decomps      ,tnmin       ,ndaypm       ,nmtimes     , & 
       amrain       ,amsnow      ,amaet        ,amtrunoff    ,amsrunoff   , &
       amdrainage   ,amtemp      ,amqa         ,&
       amsolar      ,amirup      ,amirdown     ,amsens       ,amlatent    , &  
       amlaiu       ,amlail      ,amtsoi       ,amwsoi       ,amwisoi     , &   
       amvwc        ,amawc       ,amsnod       ,amsnof       ,amnpp       , &
       amnpptot     ,amco2mic    ,amco2root    ,amco2soi     ,amco2ratio  , &
       amneetot     ,amnmintot   ,nytimes      ,aysolar      ,ayirup      , &
       ayirdown     ,aysens      ,aylatent     ,ayprcp       ,ayaet       , &  
       aytrans      ,aytrunoff   ,aysrunoff    ,aydrainage   ,aydwtot     , & 
       aywsoi       ,aywisoi     ,aytsoi       ,ayvwc        ,ayawc       , &  
       aystresstu   ,aystresstl  ,aygpp        ,aygpptot     ,aynpp       , & 
       aynpptot     ,ayco2mic    ,ayco2root    ,ayco2soi     ,ayneetot    , &
       ayrootbio    ,aynmintot   ,ayalit       ,ayblit       ,aycsoi      , & 
       aycmic       ,ayanlit     ,aybnlit      ,aynsoi       ,ayalbedo    , &
       beta1        ,beta2       ,stressfac    ,avmuir_factor,totalit     , &
       totrlit      ,totcsoi     ,totcmic      ,totanlit     ,totrnlit    , &
       totnsoi      ,totnmic     ,totlit       ,totfall      ,totnlit     , &
       firefac      ,wtot        ,storedn      ,yrleach      ,ynleach     , & 
       falll        ,fallr       ,fallw        ,clitlm       ,clitls      , &
       clitrm       ,clitrs      ,clitwm       ,clitws       ,csoislop    , &
       csoislon     ,csoipas     ,clitll       ,clitrl       ,clitwl      , &  
       tc           ,agddu       ,tempu        ,agddl        ,templ       , &
       dropu        ,dropls      ,dropl4       ,dropl3       ,plai        , &
       iday        ,imonth       ,iyear        ,iyear0       , &
       isimveg      ,spinmax     ,amts2        ,amtransu     ,amtransl    , &
       amsuvap      ,aminvap     ,amalbedo     ,amtsoil      ,amwsoil     , & 
       amwisoil     , ux         ,uy           ,taux         ,tauy        , &
       ts2          ,qs2         ,gdd0this     ,gdd5this     ,bstar       , &
       vegtype0     ,ynleach_p    ,tnmin_p     ,totnmic_p    ,totnlit_p   , &
       totanlit_p   ,totrnlit_p  ,totnsoi_p    ,storedn_p    ,adnpp       ,&
       adfalll      ,adfallr     ,adfallw      ,adcbiol      ,adcbior     ,&
       adcbiow      ,adplai      ,nstep        ,nVegClass   ,rootmode)

    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN        ) :: mcsec!global  ! current seconds in day (0 - (86400 - dtime))
    REAL(KIND=r8), INTENT(IN        ) :: pi   !global
    REAL(KIND=r8), INTENT(IN        ) :: stef !global  ! stefan-boltzmann constant (W m-2 K-4)
    REAL(KIND=r8), INTENT(IN        ) :: vonk !global  ! von karman constant (dimensionless)
    REAL(KIND=r8), INTENT(IN        ) :: grav !global  ! gravitational acceleration (m s-2)
    REAL(KIND=r8), INTENT(IN        ) :: tmelt!global  ! freezing point of water (K)
    REAL(KIND=r8), INTENT(IN        ) :: hfus !global  ! latent heat of fusion of water (J kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: hvap !global  ! latent heat of vaporization of water (J kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: hsub !global  ! latent heat of sublimation of ice (J kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: ch2o !global  ! specific heat of liquid water (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: cice !global  ! specific heat of ice (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: cair !global  ! specific heat of dry air at constant pressure (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: cvap !global  ! specific heat of water vapor at constant pressure (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: rair !global  ! gas constant for dry air (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: rvap !global  ! gas constant for water vapor (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: cappa!global  ! rair/cair
    REAL(KIND=r8), INTENT(IN        ) :: rhow !global  ! density of liquid water (all types) (kg m-3)

    ! 
    !
    INTEGER, INTENT(IN   ) :: nVegClass
    INTEGER, INTENT(IN   ) :: npoi   !global  
    INTEGER, INTENT(IN   ) :: nband  !global  
    INTEGER, INTENT(IN   ) :: nsoilay!global   ! number of soil layers
    INTEGER, INTENT(IN   ) :: nsnolay!global   ! number of snow layers
    INTEGER, INTENT(IN   ) :: npft   !global   ! number of plant functional types
    REAL(KIND=r8), INTENT(IN   ) :: epsilon!global   ! small quantity to avoid zero-divides and other
    ! truncation or machine-limit troubles with small
    ! values. should be slightly greater than o(1)
    ! machine precision
    REAL(KIND=r8), INTENT(IN   )  :: dtime !global   ! model timestep (seconds)
    LOGICAL, INTENT(IN   )  :: doalb !global    ! true if surface albedo calculation time step

    !      INCLUDE 'comhyd.h'
    REAL(KIND=r8), INTENT(OUT  ) :: ginvap (npoi)!local ! total evaporation rate from all intercepted h2o (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: gsuvap (npoi)!local ! total evaporation rate from surface (snow/soil) (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: gtrans (npoi)!local ! total transpiration rate from all vegetation canopies (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: gtransu(npoi)!local ! transpiration from upper canopy (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: gtransl(npoi)!local ! transpiration from lower canopy (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: grunof (npoi)!local ! surface runoff rate (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: gdrain (npoi)!local ! drainage rate out of bottom of lowest soil layer (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: gadjust(npoi)!local ! h2o flux due to adjustments in subroutine wadjust (kg_h2o m-2 s-1)
    !      INCLUDE 'comsum.h'
    REAL(KIND=r8), INTENT(INOUT) :: a10scalparamu(npoi)!global ! 10-day average day-time scaling parameter - upper canopy (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: a10daylightu (npoi)!global ! 10-day average day-time PAR - upper canopy (micro-Ein m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: a10scalparaml(npoi)!global ! 10-day average day-time scaling parameter - lower canopy (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: a10daylightl (npoi)!global ! 10-day average day-time PAR - lower canopy (micro-Ein m-2 s-1)
    !      INCLUDE 'compft.h'
    REAL(KIND=r8), INTENT(IN   ) :: vmax_pft(npft)!global ! nominal vmax of top leaf at 15 C (mol-co2/m**2/s) [not used]
    REAL(KIND=r8), INTENT(IN   ) :: tau15           !global ! co2/o2 specificity ratio at 15 degrees C (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: kc15           !global ! co2 kinetic parameter (mol/mol)
    REAL(KIND=r8), INTENT(IN   ) :: ko15           !global ! o2 kinetic parameter (mol/mol) 
    REAL(KIND=r8), INTENT(IN   ) :: cimax           !global ! maximum value for ci (needed for model stability)
    REAL(KIND=r8), INTENT(IN   ) :: gammaub           !global ! leaf respiration coefficient
    REAL(KIND=r8), INTENT(IN   ) :: alpha3           !global ! intrinsic quantum efficiency for C3 plants (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: theta3           !global ! photosynthesis coupling coefficient for C3 plants (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: beta3           !global ! photosynthesis coupling coefficient for C3 plants (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: coefmub           !global ! 'm' coefficient for stomatal conductance relationship
    REAL(KIND=r8), INTENT(IN   ) :: coefbub           !global ! 'b' coefficient for stomatal conductance relationship
    REAL(KIND=r8), INTENT(IN   ) :: gsubmin           !global ! absolute minimum stomatal conductance
    REAL(KIND=r8), INTENT(IN   ) :: gammauc           !global ! leaf respiration coefficient
    REAL(KIND=r8), INTENT(IN   ) :: coefmuc           !global ! 'm' coefficient for stomatal conductance relationship  
    REAL(KIND=r8), INTENT(IN   ) :: coefbuc           !global ! 'b' coefficient for stomatal conductance relationship  
    REAL(KIND=r8), INTENT(IN   ) :: gsucmin           !global ! absolute minimum stomatal conductance
    REAL(KIND=r8), INTENT(IN   ) :: gammals           !global ! leaf respiration coefficient
    REAL(KIND=r8), INTENT(IN   ) :: coefmls           !global ! 'm' coefficient for stomatal conductance relationship 
    REAL(KIND=r8), INTENT(IN   ) :: coefbls           !global ! 'b' coefficient for stomatal conductance relationship 
    REAL(KIND=r8), INTENT(IN   ) :: gslsmin           !global ! absolute minimum stomatal conductance
    REAL(KIND=r8), INTENT(IN   ) :: gammal3           !global ! leaf respiration coefficient
    REAL(KIND=r8), INTENT(IN   ) :: coefml3           !global ! 'm' coefficient for stomatal conductance relationship 
    REAL(KIND=r8), INTENT(IN   ) :: coefbl3           !global ! 'b' coefficient for stomatal conductance relationship 
    REAL(KIND=r8), INTENT(IN   ) :: gsl3min           !global ! absolute minimum stomatal conductance
    REAL(KIND=r8), INTENT(IN   ) :: gammal4           !global ! leaf respiration coefficient
    REAL(KIND=r8), INTENT(IN   ) :: alpha4           !global ! intrinsic quantum efficiency for C4 plants (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: theta4           !global ! photosynthesis coupling coefficient for C4 plants (dimensionless) 
    REAL(KIND=r8), INTENT(IN   ) :: beta4           !global ! photosynthesis coupling coefficient for C4 plants (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: coefml4           !global ! 'm' coefficient for stomatal conductance relationship
    REAL(KIND=r8), INTENT(IN   ) :: coefbl4           !global ! 'b' coefficient for stomatal conductance relationship
    REAL(KIND=r8), INTENT(IN   ) :: gsl4min           !global ! absolute minimum stomatal conductance
    !      include 'comveg.h'
    REAL(KIND=r8), INTENT(INOUT) :: wliqu    (npoi)!global ! intercepted liquid h2o on upper canopy leaf area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wliqumax              !global ! maximum intercepted water on a unit upper canopy leaf area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: wsnou    (npoi)!global ! intercepted frozen h2o (snow) on upper canopy leaf area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wsnoumax              !global ! intercepted snow capacity for upper canopy leaves (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: tu       (npoi)!global ! temperature of upper canopy leaves (K)
    REAL(KIND=r8), INTENT(INOUT) :: wliqs    (npoi)!global ! intercepted liquid h2o on upper canopy stem area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wliqsmax              !global ! maximum intercepted water on a unit upper canopy stem area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: wsnos    (npoi)!global ! intercepted frozen h2o (snow) on upper canopy stem area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wsnosmax              !global ! intercepted snow capacity for upper canopy stems (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: ts       (npoi)!global ! temperature of upper canopy stems (K)
    REAL(KIND=r8), INTENT(INOUT) :: wliql    (npoi)!global ! intercepted liquid h2o on lower canopy leaf and stem area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wliqlmax              !global ! maximum intercepted water on a unit lower canopy stem & leaf area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: wsnol    (npoi)!global ! intercepted frozen h2o (snow) on lower canopy leaf & stem area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wsnolmax              !global ! intercepted snow capacity for lower canopy leaves & stems (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: tl       (npoi)!global ! temperature of lower canopy leaves & stems(K)
    REAL(KIND=r8), INTENT(INOUT) :: topparu  (npoi)!local  ! total photosynthetically active raditaion absorbed 
    ! by top leaves of upper canopy (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: topparl  (npoi)!local  ! total photosynthetically active raditaion absorbed
    ! by top leaves of lower canopy (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: fl       (npoi)!global ! fraction of snow-free area covered by lower  canopy
    REAL(KIND=r8), INTENT(IN   ) :: fu       (npoi)!global ! fraction of overall area covered by upper canopy
    REAL(KIND=r8), INTENT(INOUT) :: lai      (npoi,2)!global ! canopy single-sided leaf area index (area leaf/area veg)
    REAL(KIND=r8), INTENT(IN   ) :: sai      (npoi,2)!global ! current single-sided stem area index
    REAL(KIND=r8), INTENT(IN   ) :: rhoveg   (nband,2)!global ! reflectance of an average leaf/stem
    REAL(KIND=r8), INTENT(IN   ) :: tauveg   (nband,2)  ! transmittance of an average leaf/stem
    REAL(KIND=r8), INTENT(IN   ) :: orieh    (2) ! fraction of leaf/stems with horizontal orientation
    REAL(KIND=r8), INTENT(IN   ) :: oriev    (2) ! fraction of leaf/stems with vertical
    REAL(KIND=r8), INTENT(INOUT) :: wliqmin      ! local ! minimum intercepted water on unit vegetated area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: wsnomin      ! local ! minimum intercepted snow on unit vegetated area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: t12      (npoi) !global ! air temperature at z12 (K)
    REAL(KIND=r8), INTENT(IN   ) :: tdripu       ! global ! decay time for dripoff of liquid intercepted by upper canopy leaves (sec)
    REAL(KIND=r8), INTENT(IN   ) :: tblowu          ! global ! decay time for blowoff of snow intercepted by upper canopy leaves (sec)
    REAL(KIND=r8), INTENT(IN   ) :: tdrips          ! global ! decay time for dripoff of liquid intercepted by upper canopy stems (sec) 
    REAL(KIND=r8), INTENT(IN   ) :: tblows          ! global ! decay time for blowoff of snow intercepted by upper canopy stems (sec)
    REAL(KIND=r8), INTENT(INOUT) :: t34      (npoi)! global ! air temperature at z34 (K)
    REAL(KIND=r8), INTENT(IN   ) :: tdripl       ! global ! decay time for dripoff of liquid intercepted
    ! by lower canopy leaves & stem (sec)
    REAL(KIND=r8), INTENT(IN   ) :: tblowl          ! global          ! decay time for blowoff of snow intercepted by lower canopy leaves & stems (sec)
    REAL(KIND=r8), INTENT(INOUT) :: ztop     (npoi,2) ! global  ! height of plant top above ground (m)
    REAL(KIND=r8), INTENT(IN   ) :: alaiml          ! global ! lower canopy leaf & stem maximum area (2 sided) for
    ! normalization of drag coefficient (m2 m-2)
    REAL(KIND=r8), INTENT(INOUT) :: zbot     (npoi,2) ! global  ! height of lowest branches above ground (m)
    REAL(KIND=r8), INTENT(IN   ) :: alaimu                ! global  ! upper canopy leaf & stem area (2 sided) for 
    ! normalization of drag coefficient (m2 m-2)
    REAL(KIND=r8), INTENT(INOUT) :: froot    (npoi,nsoilay,2)! global! fraction of root in soil layer 
    REAL(KIND=r8), INTENT(INOUT) :: q34      (npoi)         ! global! specific humidity of air at z34
    REAL(KIND=r8), INTENT(INOUT) :: q12      (npoi)         ! global! specific humidity of air at z12
    REAL(KIND=r8), INTENT(INOUT) :: su       (npoi)         ! local ! air-vegetation transfer coefficients (*rhoa) for
    ! upper canopy leaves (m s-1 * kg m-3) (A39a Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(IN   ) :: cleaf                 ! global ! empirical constant in upper canopy leaf-air 
    ! aerodynamic transfer coefficient (m s-0.5) (A39a Pollard & Thompson 95)
    REAL(KIND=r8), INTENT(IN   ) :: dleaf    (2)         ! global ! typical linear leaf dimension in aerodynamic transfer coefficient (m)
    REAL(KIND=r8), INTENT(INOUT) :: ss       (npoi)         ! local! air-vegetation transfer coefficients (*rhoa) for 
    ! upper canopy stems (m s-1 * kg m-3) (A39a Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(IN   ) :: cstem                 ! global ! empirical constant in upper canopy stem-air 
    ! aerodynamic transfer coefficient (m s-0.5) (A39a Pollard & Thompson 95)
    REAL(KIND=r8), INTENT(IN   ) :: dstem    (2)         ! global ! typical linear stem dimension in aerodynamic transfer coefficient (m)
    REAL(KIND=r8), INTENT(INOUT) :: sl       (npoi)         ! local! air-vegetation transfer coefficients (*rhoa) for 
    ! lower canopy leaves & stems (m s-1*kg m-3) (A39a Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(IN   ) :: cgrass                 ! global ! empirical constant in lower canopy-air aerodynamic 
    ! transfer coefficient (m s-0.5) (A39a Pollard & Thompson 95)
    REAL(KIND=r8), INTENT(INOUT) :: ciub     (npoi)         ! global ! intercellular co2 concentration - broadleaf (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: ciuc     (npoi)         ! global ! intercellular co2 concentration - conifer        (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(IN   ) :: exist    (npoi,npft)! global ! probability of existence of each plant functional type in a gridcell
    REAL(KIND=r8), INTENT(INOUT) :: csub     (npoi)         ! global ! leaf boundary layer co2 concentration - broadleaf (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: gsub     (npoi)         ! global ! upper canopy stomatal conductance - broadleaf  (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: csuc     (npoi)         ! global ! leaf boundary layer co2 concentration - conifer   (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: gsuc     (npoi)         ! global ! upper canopy stomatal conductance - conifer    (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: agcub    (npoi)         ! local  ! canopy average gross photosynthesis rate - broadleaf  (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: agcuc    (npoi)         ! local  ! canopy average gross photosynthesis rate - conifer    (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: ancub    (npoi)         ! local ! canopy average net photosynthesis rate - broadleaf    (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: ancuc    (npoi)         ! local ! canopy average net photosynthesis rate - conifer          (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: totcondub(npoi)         ! local ! 
    REAL(KIND=r8), INTENT(INOUT) :: totconduc(npoi)         ! local !
    REAL(KIND=r8), INTENT(INOUT) :: cils     (npoi)         ! global ! intercellular co2 concentration - shrubs        (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: cil3     (npoi)         ! global ! intercellular co2 concentration - c3 plants (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: cil4     (npoi)         ! global ! intercellular co2 concentration - c4 plants (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: csls     (npoi)         ! global ! leaf boundary layer co2 concentration - shrubs   (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: gsls     (npoi)         ! global ! lower canopy stomatal conductance - shrubs     (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: csl3     (npoi)         ! global ! leaf boundary layer co2 concentration - c3 plants (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: gsl3     (npoi)         ! global ! lower canopy stomatal conductance - c3 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: csl4     (npoi)         ! global ! leaf boundary layer co2 concentration - c4 plants (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: gsl4     (npoi)         ! global ! lower canopy stomatal conductance - c4 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: agcls    (npoi)         ! local  ! canopy average gross photosynthesis rate - shrubs          (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: agcl4    (npoi)         ! local ! canopy average gross photosynthesis rate - c4 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: agcl3    (npoi)         ! local ! canopy average gross photosynthesis rate - c3 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: ancls    (npoi)         ! local ! canopy average net photosynthesis rate - shrubs          (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: ancl4    (npoi)         ! local ! canopy average net photosynthesis rate - c4 grasses   (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: ancl3    (npoi)         ! local ! canopy average net photosynthesis rate - c3 grasses   (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: totcondls(npoi)         ! local ! 
    REAL(KIND=r8), INTENT(INOUT) :: totcondl3(npoi)         ! local !
    REAL(KIND=r8), INTENT(INOUT) :: totcondl4(npoi)         ! local !
    REAL(KIND=r8), INTENT(IN   ) :: chu(1:nVegClass)             ! global ! heat capacity of upper canopy leaves per unit leaf area (J kg-1 m-2)
    REAL(KIND=r8), INTENT(IN   ) :: chs(1:nVegClass)             ! global ! heat capacity of upper canopy stems per unit stem area (J kg-1 m-2)
    REAL(KIND=r8), INTENT(IN   ) :: chl(1:nVegClass)             ! global ! heat capacity of lower canopy leaves & stems per unit leaf/stem area (J kg-1 m-2)
    REAL(KIND=r8), INTENT(INOUT) :: frac     (npoi,npft)! global ! fraction of canopy occupied by each plant functional type
    REAL(KIND=r8), INTENT(INOUT) :: tlsub    (npoi)         ! global ! temperature of lower canopy vegetation buried by snow (K)
    !      INCLUDE 'comsat.h'    
    !      include 'comsno.h'
    REAL(KIND=r8), INTENT(IN   ) :: z0sno  ! global ! roughness length of snow surface (m)
    REAL(KIND=r8), INTENT(IN   ) :: rhos   ! global ! density of snow (kg m-3)
    REAL(KIND=r8), INTENT(IN   ) :: consno ! global ! thermal conductivity of snow (W m-1 K-1)
    REAL(KIND=r8), INTENT(IN   ) :: hsnotop! global ! thickness of top snow layer (m)
    REAL(KIND=r8), INTENT(IN   ) :: hsnomin! global ! minimum total thickness of snow (m)
    REAL(KIND=r8), INTENT(IN   ) :: fimin  ! global ! minimum fractional snow cover
    REAL(KIND=r8), INTENT(IN   ) :: fimax  ! global ! maximum fractional snow cover
    REAL(KIND=r8), INTENT(INOUT) :: fi     (npoi)! global ! fractional snow cover
    REAL(KIND=r8), INTENT(INOUT) :: tsno   (npoi,nsnolay)! global ! temperature of snow layers (K)
    REAL(KIND=r8), INTENT(INOUT) :: hsno   (npoi,nsnolay)! global ! thickness of snow layers (m)

    !      INCLUDE 'comsoi.h'
    REAL(KIND=r8), INTENT(IN   ) :: sand    (npoi,nsoilay)! global ! percent sand of soil
    REAL(KIND=r8), INTENT(IN   ) :: clay    (npoi,nsoilay)! global ! percent clay of soil
    REAL(KIND=r8), INTENT(IN   ) :: poros   (npoi,nsoilay)! global ! porosity (mass of h2o per unit vol at sat / rhow)
    REAL(KIND=r8), INTENT(INOUT) :: wsoi    (npoi,nsoilay)! global ! fraction of soil pore space containing liquid water
    REAL(KIND=r8), INTENT(INOUT) :: wisoi   (npoi,nsoilay)! global ! fraction of soil pore space containing ice
    REAL(KIND=r8), INTENT(INOUT) :: consoi  (npoi,nsoilay)! local  ! thermal conductivity of each soil layer (W m-1 K-1)
    REAL(KIND=r8), INTENT(IN   ) :: zwpmax                   ! global ! assumed maximum fraction of soil surface 
    ! covered by puddles (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: wpud    (npoi)! global ! liquid content of puddles per soil area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: wipud   (npoi)! global ! ice content of puddles per soil area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wpudmax                        ! normalization constant for puddles (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: qglif   (npoi,4) ! local ! 1: fraction of soil evap (fvapg) from soil liquid
    ! 2: fraction of soil evap (fvapg) from soil ice
    ! 3: fraction of soil evap (fvapg) from puddle liquid
    ! 4: fraction of soil evap (fvapg) from puddle ice
    REAL(KIND=r8), INTENT(INOUT) :: tsoi    (npoi,nsoilay)! global        ! soil temperature for each layer (K)
    REAL(KIND=r8), INTENT(INOUT) :: hvasug  (npoi) ! local ! latent heat of vap/subl, for soil surface (J kg-1)
    REAL(KIND=r8), INTENT(INOUT) :: hvasui  (npoi)! local ! latent heat of vap/subl, for snow surface (J kg-1)
    REAL(KIND=r8), INTENT(IN   ) :: albsav  (npoi)! global ! saturated soil surface albedo (visible waveband)
    REAL(KIND=r8), INTENT(IN   ) :: albsan  (npoi)! global ! saturated soil surface albedo (near-ir waveband)
    REAL(KIND=r8), INTENT(INOUT) :: tg      (npoi)! global ! soil skin temperature (K)
    REAL(KIND=r8), INTENT(INOUT) :: ti      (npoi)! global ! snow skin temperature (K)
    REAL(KIND=r8), INTENT(IN   ) :: z0soi   (npoi)! global ! roughness length of soil surface (m)
    REAL(KIND=r8), INTENT(IN   ) :: swilt   (npoi,nsoilay)! global ! wilting soil moisture value (fraction of pore space)
    REAL(KIND=r8), INTENT(IN   ) :: sfield  (npoi,nsoilay)! global ! field capacity soil moisture value (fraction of pore space)
    REAL(KIND=r8), INTENT(INOUT) :: stressl (npoi,nsoilay)! local ! soil moisture stress factor for the lower canopy (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: stressu (npoi,nsoilay)! local ! soil moisture stress factor for the upper canopy (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: stresstl(npoi) ! local ! sum of stressl over all 6 soil layers (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: stresstu(npoi)! local ! sum of stressu over all 6 soil layers (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: csoi    (npoi,nsoilay)! global ! specific heat of soil, no pore spaces (J kg-1 deg-1)
    REAL(KIND=r8), INTENT(IN   ) :: rhosoi  (npoi,nsoilay)! global ! soil density (without pores, not bulk) (kg m-3)
    REAL(KIND=r8), INTENT(IN   ) :: hsoi    (npoi,nsoilay+1)   ! global ! soil layer thickness (m)
    REAL(KIND=r8), INTENT(IN   ) :: suction (npoi,nsoilay)! global ! saturated matric potential (m-h2o)
    REAL(KIND=r8), INTENT(IN   ) :: bex     (npoi,nsoilay)! global ! exponent "b" in soil water potential
    REAL(KIND=r8), INTENT(INOUT) :: upsoiu  (npoi,nsoilay)! local  ! soil water uptake from transpiration (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: upsoil  (npoi,nsoilay)! local  ! soil water uptake from transpiration (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: heatg   (npoi)        ! local        ! net heat flux into soil surface (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: heati   (npoi)        ! local        ! net heat flux into snow surface (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: hydraul (npoi,nsoilay)! global ! saturated hydraulic conductivity (m/s)
    REAL(KIND=r8), INTENT(INOUT) :: porosflo(npoi,nsoilay)! global ! porosity after reduction by ice content
    INTEGER, INTENT(IN   ) :: ibex    (npoi,nsoilay)! global ! nint(bex), used for cpu speed
    REAL(KIND=r8), INTENT(IN   ) :: bperm (npoi)  ! global! lower b.c. for soil profile drainage 
    ! (0.0 = impermeable; 1.0 = fully permeable)
    ! (0.0 = impermeable; 1.0 = fully permeable)
    REAL(KIND=r8), INTENT(INOUT) :: hflo    (npoi,nsoilay+1)  ! downward heat transport through soil layers (W m-2)


    !   INCLUDE 'comatm.h'
    REAL(KIND=r8), INTENT(IN   ) :: ta     (npoi)         ! global ! air temperature (K)
    REAL(KIND=r8), INTENT(INOUT) :: asurd  (npoi,nband) ! local  ! direct albedo of surface system
    REAL(KIND=r8), INTENT(INOUT) :: asuri  (npoi,nband) ! local  ! diffuse albedo of surface system 
    REAL(KIND=r8), INTENT(IN   ) :: coszen (npoi)         ! global ! cosine of solar zenith angle
    REAL(KIND=r8), INTENT(IN   ) :: solad  (npoi,nband) ! global ! direct downward solar flux (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: solai  (npoi,nband) ! global ! diffuse downward solar flux (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: fira   (npoi)         ! global ! incoming ir flux (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: raina  (npoi)         ! global ! rainfall rate (mm/s or kg m-2 s-1)
    REAL(KIND=r8), INTENT(IN   ) :: qa     (npoi)         ! global ! specific humidity (kg_h2o/kg_air)
    REAL(KIND=r8), INTENT(IN   ) :: psurf  (npoi)         ! global ! surface pressure (Pa)
    REAL(KIND=r8), INTENT(IN   ) :: snowa  (npoi)         ! global ! snowfall rate (mm/s or kg m-2 s-1 of water)
    REAL(KIND=r8), INTENT(IN   ) :: ua     (npoi)         ! global ! wind speed (m s-1)
    REAL(KIND=r8), INTENT(IN   ) :: o2conc                 ! global ! o2 concentration (mol/mol)
    REAL(KIND=r8), INTENT(IN   ) :: co2conc (npoi)                ! global ! co2 concentration (mol/mol)

    !   INCLUDE 'com1d.h'
    REAL(KIND=r8), INTENT(INOUT) :: fwetu    (npoi)         ! local ! fraction of upper canopy leaf area wetted by intercepted liquid and/or snow
    REAL(KIND=r8), INTENT(INOUT) :: rliqu    (npoi)         ! local ! proportion of fwetu due to liquid
    REAL(KIND=r8), INTENT(INOUT) :: fwets    (npoi)         ! local ! fraction of upper canopy stem area wetted by intercepted liquid and/or snow
    REAL(KIND=r8), INTENT(INOUT) :: rliqs    (npoi)         ! local ! proportion of fwets due to liquid
    REAL(KIND=r8), INTENT(INOUT) :: fwetl    (npoi)         ! local ! fraction of lower canopy stem & leaf area wetted by
    ! intercepted liquid and/or snow
    REAL(KIND=r8), INTENT(INOUT) :: rliql    (npoi)         ! local ! proportion of fwetl due to liquid
    REAL(KIND=r8), INTENT(INOUT) :: solu     (npoi)         ! local ! solar flux (direct + diffuse) absorbed by upper 
    ! canopy leaves per unit canopy area (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: sols     (npoi)         ! local ! solar flux (direct + diffuse) absorbed by upper 
    ! canopy stems per unit canopy area (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: soll     (npoi)         ! local ! solar flux (direct + diffuse) absorbed by lower 
    ! canopy leaves and stems per unit canopy area (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: solg     (npoi)         ! local ! solar flux (direct + diffuse) absorbed by unit 
    ! snow-free soil (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: soli     (npoi)         ! local ! solar flux (direct + diffuse) absorbed by unit 
    ! snow surface (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: scalcoefl(npoi,4)   ! local ! term needed in lower canopy scaling
    REAL(KIND=r8), INTENT(INOUT) :: scalcoefu(npoi,4)   ! local ! term needed in upper canopy scaling
    INTEGER, INTENT(INOUT) :: indsol   (npoi)         ! local ! index of current strip for points with positive coszen
    REAL(KIND=r8), INTENT(INOUT) :: albsod   (npoi)         ! local ! direct  albedo for soil surface (visible or IR)
    REAL(KIND=r8), INTENT(INOUT) :: albsoi   (npoi)         ! local ! diffuse albedo for soil surface (visible or IR)
    REAL(KIND=r8), INTENT(INOUT) :: albsnd   (npoi)         ! local ! direct  albedo for snow surface (visible or IR)
    REAL(KIND=r8), INTENT(INOUT) :: albsni   (npoi)         ! local ! diffuse albedo for snow surface (visible or IR)
    REAL(KIND=r8), INTENT(OUT  ) :: relod    (npoi)         ! local ! upward direct radiation per unit icident direct beam on lower canopy (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: reloi    (npoi)         ! local ! upward diffuse radiation per unit incident diffuse 
    ! radiation on lower canopy (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: reupd    (npoi)         ! local ! upward direct radiation per unit incident direct 
    ! radiation on upper canopy (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: reupi    (npoi)         ! local ! upward diffuse radiation per unit incident diffuse 
    ! radiation on upper canopy (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: ablod    (npoi)         ! local ! fraction of direct  radiation absorbed by lower canopy
    REAL(KIND=r8), INTENT(INOUT) :: abloi    (npoi)         ! local ! fraction of diffuse radiation absorbed by lower canopy
    REAL(KIND=r8), INTENT(INOUT) :: flodd    (npoi)         ! local ! downward direct radiation per unit incident direct
    ! radiation on lower canopy (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: dummy    (npoi)         ! local ! placeholder, always = 0: no direct flux produced for diffuse incident
    REAL(KIND=r8), INTENT(INOUT) :: flodi    (npoi)         ! local ! downward diffuse radiation per unit incident direct
    ! radiation on lower canopy (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: floii    (npoi)         ! local ! downward diffuse radiation per unit incident 
    ! diffuse radiation on lower canopy
    REAL(KIND=r8), INTENT(INOUT) :: terml    (npoi,7)         ! local ! term needed in lower canopy scaling
    REAL(KIND=r8), INTENT(INOUT) :: termu    (npoi,7)         ! local ! term needed in upper canopy scaling
    REAL(KIND=r8), INTENT(INOUT) :: abupd    (npoi)         ! local ! fraction of direct  radiation absorbed by upper canopy
    REAL(KIND=r8), INTENT(INOUT) :: abupi    (npoi)         ! local ! fraction of diffuse radiation absorbed by upper canopy
    REAL(KIND=r8), INTENT(INOUT) :: fupdd    (npoi)         ! local ! downward direct radiation per unit incident direct
    ! beam on upper canopy (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: fupdi    (npoi)         ! local ! downward diffuse radiation per unit icident direct
    ! radiation on upper canopy (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: fupii    (npoi)         ! local ! downward diffuse radiation per unit incident diffuse
    ! radiation on upper canopy (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: sol2d    (npoi)         ! local ! direct downward radiation  out of upper canopy 
    ! per unit vegetated (upper) area (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: sol2i    (npoi)         ! local ! diffuse downward radiation out of upper
    ! canopy per unit vegetated (upper) area(W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: sol3d    (npoi)         ! local ! direct downward radiation  out of upper
    ! canopy + gaps per unit grid cell area (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: sol3i    (npoi)         ! local ! diffuse downward radiation out of upper
    ! canopy + gaps per unit grid cell area (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: firb     (npoi)         ! local ! net upward ir radiation at reference
    ! atmospheric level za (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: firs     (npoi)         ! local ! ir radiation absorbed by upper canopy stems (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: firu     (npoi)         ! local ! ir raditaion absorbed by upper canopy leaves (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: firl     (npoi)         ! local ! ir radiation absorbed by lower canopy leaves and stems (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: firg     (npoi)         ! local ! ir radiation absorbed by soil/ice (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: firi     (npoi)         ! local ! ir radiation absorbed by snow (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: snowg    (npoi)         ! local ! snowfall rate at soil level (kg h2o m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: tsnowg   (npoi)         ! local ! snowfall temperature at soil level (K) 
    REAL(KIND=r8), INTENT(INOUT) :: tsnowl   (npoi)         ! local ! snowfall temperature below upper canopy (K)
    REAL(KIND=r8), INTENT(INOUT) :: pfluxl   (npoi)         ! local ! heat flux on lower canopy leaves & stems due to intercepted h2o (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: raing    (npoi)         ! local ! rainfall rate at soil level (kg m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: traing   (npoi)         ! local ! rainfall temperature at soil level (K)
    REAL(KIND=r8), INTENT(INOUT) :: trainl   (npoi)         ! local ! rainfall temperature below upper canopy (K)
    REAL(KIND=r8), INTENT(INOUT) :: snowl    (npoi)         ! local ! snowfall rate below upper canopy (kg h2o m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: tsnowu   (npoi)         ! local ! snowfall temperature above upper canopy (K)
    REAL(KIND=r8), INTENT(INOUT) :: pfluxu   (npoi)         ! local ! heat flux on upper canopy leaves due to intercepted h2o (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: rainu    (npoi)         ! local  ! rainfall rate above upper canopy (kg m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: trainu   (npoi)         ! local  ! rainfall temperature above upper canopy (K)
    REAL(KIND=r8), INTENT(INOUT) :: snowu    (npoi)         ! local ! snowfall rate above upper canopy (kg h2o m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: pfluxs   (npoi)         ! local ! heat flux on upper canopy stems due to intercepted h2o (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: rainl    (npoi)         ! local ! rainfall rate below upper canopy (kg m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: bps      (npoi)            ! local ! (ps/p) ** (rair/cair) for atmospheric level  (const)
    REAL(KIND=r8):: rhoa     (npoi)         ! local ! air density at za (allowing for h2o vapor) (kg m-3)
    REAL(KIND=r8), INTENT(INOUT) :: cp       (npoi)         ! local ! specific heat of air at za (allowing for h2o vapor) (J kg-1 K-1)
    REAL(KIND=r8), INTENT(INOUT) :: za       (npoi)         ! local ! height above the surface of atmospheric forcing (m)
    REAL(KIND=r8), INTENT(INOUT) :: bdl      (npoi)         ! local ! aerodynamic coefficient ([(tau/rho)/u**2] for
    ! laower canopy (A31/A30 Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(INOUT) :: dil      (npoi)         ! local ! inverse of momentum diffusion coefficient within lower canopy (m)
    REAL(KIND=r8), INTENT(INOUT) :: z3       (npoi)         ! local ! effective top of the lower canopy (for momentum) (m)
    REAL(KIND=r8), INTENT(INOUT) :: z4       (npoi)         ! local ! effective bottom of the lower canopy (for momentum) (m)
    REAL(KIND=r8), INTENT(INOUT) :: z34      (npoi)         ! local ! effective middle of the lower canopy (for momentum) (m)
    REAL(KIND=r8), INTENT(INOUT) :: exphl    (npoi)         ! local ! exp(lamda/2*(z3-z4)) for lower canopy (A30 Pollard & Thompson)
    REAL(KIND=r8), INTENT(INOUT) :: expl     (npoi)         ! local ! exphl**2
    REAL(KIND=r8), INTENT(INOUT) :: displ    (npoi)         ! local ! zero-plane displacement height for lower canopy (m)
    REAL(KIND=r8), INTENT(INOUT) :: bdu      (npoi)         ! local ! aerodynamic coefficient ([(tau/rho)/u**2] for upper
    ! canopy (A31/A30 Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(INOUT) :: diu      (npoi)         ! local ! inverse of momentum diffusion coefficient within upper canopy (m)
    REAL(KIND=r8), INTENT(INOUT) :: z1       (npoi)         ! local ! effective top of upper canopy (for momentum) (m)
    REAL(KIND=r8), INTENT(INOUT) :: z2       (npoi)         ! local ! effective bottom of the upper canopy (for momentum) (m)
    REAL(KIND=r8), INTENT(INOUT) :: z12      (npoi)         ! local ! effective middle of the upper canopy (for momentum) (m)
    REAL(KIND=r8), INTENT(INOUT) :: exphu    (npoi)         ! local ! exp(lamda/2*(z3-z4)) for upper canopy (A30 Pollard & Thompson)
    REAL(KIND=r8), INTENT(INOUT) :: expu     (npoi)         ! local ! exphu**2
    REAL(KIND=r8), INTENT(INOUT) :: dispu    (npoi)         ! local ! zero-plane displacement height for upper canopy (m)
    REAL(KIND=r8), INTENT(INOUT) :: alogg    (npoi)         ! local ! log of soil roughness
    REAL(KIND=r8), INTENT(INOUT) :: alogi    (npoi)         ! local ! log of snow roughness
    REAL(KIND=r8), INTENT(INOUT) :: alogav   (npoi)         ! local ! average of alogi and alogg 
    REAL(KIND=r8), INTENT(INOUT) :: alog4    (npoi)         ! local ! log (max(z4, 1.1*z0sno, 1.1*z0soi)) 
    REAL(KIND=r8), INTENT(INOUT) :: alog3    (npoi)         ! local ! log (z3 - displ)
    REAL(KIND=r8), INTENT(INOUT) :: alog2    (npoi)         ! local ! log (z2 - displ)
    REAL(KIND=r8), INTENT(INOUT) :: alog1    (npoi)         ! local ! log (z1 - dispu) 
    REAL(KIND=r8), INTENT(INOUT) :: aloga    (npoi)         ! local ! log (za - dispu) 
    REAL(KIND=r8), INTENT(INOUT) :: u2       (npoi)         ! local ! wind speed at level z2 (m s-1)
    REAL(KIND=r8), INTENT(INOUT) :: alogu    (npoi)         ! local ! log (roughness length of upper canopy)
    REAL(KIND=r8), INTENT(INOUT) :: alogl    (npoi)         ! local ! log (roughness length of lower canopy)
    REAL(KIND=r8), INTENT(INOUT) :: richl    (npoi)         ! local ! richardson number for air above upper canopy (z3 to z2)
    REAL(KIND=r8), INTENT(INOUT) :: straml   (npoi)         ! local ! momentum correction factor for stratif between
    ! upper & lower canopy (z3 to z2) (louis et al.)
    REAL(KIND=r8), INTENT(INOUT)  :: strahl   (npoi)         ! local ! heat/vap correction factor for stratif between
    ! upper & lower canopy (z3 to z2) (louis et al.)
    REAL(KIND=r8), INTENT(INOUT)  :: richu    (npoi)         ! local ! richardson number for air between upper & lower canopy (z1 to za)
    REAL(KIND=r8), INTENT(INOUT)  :: stramu   (npoi)         ! local ! momentum correction factor for stratif above
    ! upper canopy (z1 to za) (louis et al.)
    REAL(KIND=r8), INTENT(INOUT)  :: strahu   (npoi)         ! local ! heat/vap correction factor for stratif above
    ! upper canopy (z1 to za) (louis et al.)
    REAL(KIND=r8), INTENT(INOUT)  :: u1       (npoi)         ! local ! wind speed at level z1 (m s-1)
    REAL(KIND=r8), INTENT(INOUT)  :: u12      (npoi)         ! local ! wind speed at level z12 (m s-1)
    REAL(KIND=r8), INTENT(INOUT)  :: u3       (npoi)         ! local ! wind speed at level z3 (m s-1)
    REAL(KIND=r8), INTENT(INOUT)  :: u34      (npoi)         ! local ! wind speed at level z34 (m s-1)
    REAL(KIND=r8), INTENT(INOUT)  :: u4       (npoi)         ! local ! wind speed at level z4 (m s-1)
    REAL(KIND=r8), INTENT(INOUT)  :: cu       (npoi)         ! local ! air transfer coefficient (*rhoa) (m s-1 kg m-3) for
    ! upper air region (z12 --> za) (A35 Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(INOUT)  :: cl       (npoi)         ! local ! air transfer coefficient (*rhoa) (m s-1 kg m-3)
    ! between the 2 canopies (z34 --> z12) (A36 Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(INOUT)  :: sg       (npoi)         ! local ! air-soil transfer coefficient
    REAL(KIND=r8), INTENT(INOUT)  :: si       (npoi)         ! local ! air-snow transfer coefficient
    REAL(KIND=r8), INTENT(INOUT)  :: fwetux   (npoi)         ! local ! fraction of upper canopy leaf area wetted if dew forms
    REAL(KIND=r8), INTENT(INOUT)  :: fwetsx   (npoi)         ! local ! fraction of upper canopy stem area wetted if dew forms
    REAL(KIND=r8), INTENT(INOUT)  :: fwetlx   (npoi)         ! local ! fraction of lower canopy leaf and stem area wetted if dew forms
    REAL(KIND=r8), INTENT(INOUT)  :: fsena    (npoi)         ! local ! downward sensible heat flux between za & z12 at za (W m-2)
    REAL(KIND=r8), INTENT(INOUT)  :: fseng    (npoi)         ! local ! upward sensible heat flux between soil surface & air at z34 (W m-2)
    REAL(KIND=r8), INTENT(INOUT)  :: fseni    (npoi)         ! local ! upward sensible heat flux between snow surface & air at z34 (W m-2)
    REAL(KIND=r8), INTENT(INOUT)  :: fsenu    (npoi)         ! local ! sensible heat flux from upper canopy leaves to air (W m-2)
    REAL(KIND=r8), INTENT(INOUT)  :: fsens    (npoi)         ! local ! sensible heat flux from upper canopy stems to air (W m-2)
    REAL(KIND=r8), INTENT(INOUT)  :: fsenl    (npoi)         ! local ! sensible heat flux from lower canopy to air (W m-2)
    REAL(KIND=r8), INTENT(INOUT)  :: fvapa    (npoi)         ! local ! downward h2o vapor flux between za & z12 at za (kg m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT)  :: fvaput   (npoi)         ! local ! h2o vapor flux (transpiration from dry parts) 
    ! between upper canopy leaves and air at z12 (kg m-2 s-1/ LAI upper canopy/ fu)
    REAL(KIND=r8), INTENT(INOUT)  :: fvaps    (npoi)         ! local ! h2o vapor flux (evaporation from wet surface)
    ! between upper canopy stems and air at z12 (kg m-2 s-1 / SAI lower canopy / fu)
    REAL(KIND=r8), INTENT(INOUT)  :: fvaplw   (npoi)         ! local ! h2o vapor flux (evaporation from wet surface) 
    ! between lower canopy leaves & stems and air at z34 (kg m-2 s-1/ LAI lower canopy/ fl)
    REAL(KIND=r8), INTENT(INOUT)  :: fvaplt   (npoi)         ! local ! h2o vapor flux (transpiration) 
    ! between lower canopy & air at z34 (kg m-2 s-1 / LAI lower canopy / fl)
    REAL(KIND=r8), INTENT(INOUT)  :: fvapg    (npoi)         ! local ! h2o vapor flux (evaporation) between soil & air 
    ! at z34 (kg m-2 s-1/bare ground fraction)
    REAL(KIND=r8), INTENT(INOUT)  :: fvapi    (npoi)         ! local ! h2o vapor flux (evaporation) between snow & air at z34 (kg m-2 s-1 / fi )
    REAL(KIND=r8), INTENT(INOUT)  :: fvapuw   (npoi)         ! local ! h2o vapor flux (evaporation from wet parts)
    ! between upper canopy leaves and air at z12 (kg m-2 s-1/ LAI upper canopy/ fu)
    REAL(KIND=r8), INTENT(INOUT)  :: td       (npoi)      ! global! daily average temperature (K)

    REAL(KIND=r8), INTENT(IN   )  :: vzero   (npoi)         ! global! a real array of zeros, of length npoi


    INTEGER, INTENT(IN   ) :: ndaypy              ! global! number of days per year
    REAL(KIND=r8), INTENT(OUT  ) :: nppdummy (npoi,npft)! local ! canopy NPP before accounting for stem and root respiration
    REAL(KIND=r8), INTENT(INOUT) :: tgpp     (npoi,npft)! local ! instantaneous GPP for each pft (mol-CO2 / m-2 / second)
    REAL(KIND=r8), INTENT(OUT  ) :: tgpptot  (npoi)         ! local ! instantaneous gpp (mol-CO2 / m-2 / second)
    REAL(KIND=r8), INTENT(INOUT) :: tnpp     (npoi,npft)! local ! instantaneous NPP for each pft (mol-CO2 / m-2 / second)
    REAL(KIND=r8), INTENT(IN   ) :: cbiow    (npoi,npft)! global! carbon in woody biomass pool (kg_C m-2)
    REAL(KIND=r8), INTENT(IN   ) :: sapfrac  (npoi)         ! global! fraction of woody biomass that is in sapwood
    REAL(KIND=r8), INTENT(IN   ) :: cbior    (npoi,npft)! global! carbon in fine root biomass pool (kg_C m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: tnpptot  (npoi)         ! local ! instantaneous npp (mol-CO2 / m-2 / second)
    REAL(KIND=r8), INTENT(INOUT) :: tco2root (npoi)         ! local ! instantaneous fine co2 flux from soil (mol-CO2 / m-2 / second)
    REAL(KIND=r8), INTENT(OUT  ) :: tneetot  (npoi)         ! local ! instantaneous net ecosystem exchange of co2 per timestep (mol-CO2 / m-2 / second)
    REAL(KIND=r8), INTENT(INOUT) :: tco2mic  (npoi)         ! local ! instantaneous microbial co2 flux from soil (mol-CO2 / m-2 / second)
    REAL(KIND=r8), INTENT(INOUT) :: a10td    (npoi)     ! global! 10-day average daily air temperature (K)
    REAL(KIND=r8), INTENT(INOUT) :: a10ancub (npoi)     ! global! 10-day average canopy photosynthesis rate - broadleaf (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: a10ancuc (npoi)     ! global! 10-day average canopy photosynthesis rate - conifer (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: a10ancls (npoi)     ! global! 10-day average canopy photosynthesis rate - shrubs (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: a10ancl3 (npoi)     ! global! 10-day average canopy photosynthesis rate - c3 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: a10ancl4 (npoi)     ! global! 10-day average canopy photosynthesis rate - c4 grasses (mol_co2 m-2 s-1)

    INTEGER, INTENT(INOUT) :: ndtimes(npoi)             ! global! counter for daily average calculations
    REAL(KIND=r8), INTENT(INOUT) :: adrain    (npoi)! global! daily average rainfall rate (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: adsnow    (npoi)! global! daily average snowfall rate (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: adaet     (npoi)! global! daily average aet (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: adtrunoff (npoi)! global! daily average total runoff (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: adsrunoff (npoi)! global! daily average surface runoff (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: addrainage(npoi)! global! daily average drainage (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: adrh      (npoi)! global! daily average rh (percent)
    REAL(KIND=r8), INTENT(INOUT) :: adsnod    (npoi)! global! daily average snow depth (m)
    REAL(KIND=r8), INTENT(INOUT) :: adsnof    (npoi)! global! daily average snow fraction (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: adwsoi    (npoi)! global! daily average soil moisture (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: adtsoi    (npoi)! global! daily average soil temperature (c)
    REAL(KIND=r8), INTENT(INOUT) :: adwisoi   (npoi)! global! daily average soil ice (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: adtlaysoi (npoi)! global! daily average soil temperature (c) of top layer
    REAL(KIND=r8), INTENT(INOUT) :: adwlaysoi (npoi)! global! daily average soil moisture of top layer(fraction)
    REAL(KIND=r8), INTENT(INOUT) :: adwsoic   (npoi)! global! daily average soil moisture using root profile weighting (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: adtsoic   (npoi)! global! daily average soil temperature (c) using profile weighting
    REAL(KIND=r8), INTENT(INOUT) :: adco2mic  (npoi)! global! daily accumulated co2 respiration from microbes (kg_C m-2 /day)
    REAL(KIND=r8), INTENT(INOUT) :: adco2root (npoi)! global! daily accumulated co2 respiration from roots (kg_C m-2 /day)
    REAL(KIND=r8), INTENT(INOUT) :: adco2soi  (npoi)! global! daily accumulated co2 respiration from soil(total) (kg_C m-2 /day)
    REAL(KIND=r8), INTENT(INOUT) :: adco2ratio(npoi)! global! ratio of root to total co2 respiration
    REAL(KIND=r8), INTENT(INOUT) :: adnmintot (npoi)! global! daily accumulated net nitrogen mineralization (kg_N m-2 /day)
    REAL(KIND=r8), INTENT(INOUT) :: decompl   (npoi)! global! litter decomposition factor              (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: decomps   (npoi)! global! soil organic matter decomposition factor      (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: tnmin     (npoi)! global! instantaneous nitrogen mineralization (kg_N m-2/timestep)

    INTEGER, INTENT(IN   ) :: ndaypm    (12)  ! global! number of days per month


    INTEGER, INTENT(INOUT) :: nmtimes        (npoi)           ! global! counter for monthly average calculations
    REAL(KIND=r8), INTENT(INOUT) :: amrain        (npoi)     ! global! monthly average rainfall rate (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: amsnow        (npoi)     ! global! monthly average snowfall rate (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: amaet        (npoi)     ! global! monthly average aet (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: amtrunoff  (npoi)     ! global! monthly average total runoff (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: amsrunoff  (npoi)     ! global! monthly average surface runoff (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: amdrainage (npoi)     ! global! monthly average drainage (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: amtemp        (npoi)     ! global! monthly average air temperature (C)
    REAL(KIND=r8), INTENT(INOUT) :: amqa        (npoi)     ! global! monthly average specific humidity (kg-h2o/kg-air)
    REAL(KIND=r8), INTENT(INOUT) :: amsolar        (npoi)     ! global! monthly average incident solar radiation (W/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: amirup        (npoi)     ! global! monthly average upward ir radiation (W/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: amirdown   (npoi)     ! global! monthly average downward ir radiation (W/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: amsens        (npoi)     ! global! monthly average sensible heat flux (W/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: amlatent   (npoi)     ! global! monthly average latent heat flux (W/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: amlaiu        (npoi)     ! global! monthly average lai for upper canopy (m**2/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: amlail        (npoi)     ! global! monthly average lai for lower canopy (m**2/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: amtsoi        (npoi)     ! global! monthly average 1m soil temperature (C)
    REAL(KIND=r8), INTENT(INOUT) :: amwsoi        (npoi)     ! global! monthly average 1m soil moisture (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: amwisoi        (npoi)     ! global! monthly average 1m soil ice (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: amvwc        (npoi)     ! global! monthly average 1m volumetric water content (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: amawc        (npoi)     ! global! monthly average 1m plant-available water content (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: amsnod        (npoi)     ! global! monthly average snow depth (m)
    REAL(KIND=r8), INTENT(INOUT) :: amsnof        (npoi)     ! global! monthly average snow fraction (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: amnpp        (npoi,npft)! global! monthly total npp for each plant type (kg-C/m**2/month)
    REAL(KIND=r8), INTENT(INOUT) :: adnpp      (npoi,npft)! global! monthly total npp for each plant type (kg-C/m**2/day)
    REAL(KIND=r8), INTENT(OUT  ) :: amnpptot   (npoi)     ! local ! monthly total npp for ecosystem (kg-C/m**2/month)
    REAL(KIND=r8), INTENT(INOUT) :: amco2mic   (npoi)     ! global! monthly total CO2 flux from microbial respiration (kg-C/m**2/month)
    REAL(KIND=r8), INTENT(INOUT) :: amco2root  (npoi)     ! global! monthly total CO2 flux from soil due to root
    ! respiration (kg-C/m**2/month)
    REAL(KIND=r8), INTENT(OUT  ) :: amco2soi   (npoi)     ! local ! monthly total soil CO2 flux from microbial
    !          and root respiration (kg-C/m**2/month)
    REAL(KIND=r8), INTENT(OUT  ) :: amco2ratio (npoi)     ! local ! monthly ratio of root to total co2 flux
    REAL(KIND=r8), INTENT(OUT  ) :: amneetot   (npoi)     ! local ! monthly total net ecosystem exchange of CO2 (kg-C/m**2/month)
    REAL(KIND=r8), INTENT(INOUT) :: amnmintot  (npoi)     ! global! monthly total N mineralization from microbes (kg-N/m**2/month)
    REAL(KIND=r8), INTENT(INOUT) :: amalbedo   (npoi)         
    REAL(KIND=r8), INTENT(INOUT) :: amtsoil    (npoi, nsoilay) 
    REAL(KIND=r8), INTENT(INOUT) :: amwsoil    (npoi, nsoilay) 
    REAL(KIND=r8), INTENT(INOUT) :: amwisoil   (npoi, nsoilay)
    REAL(KIND=r8), INTENT(INOUT) :: amts2      (npoi)     ! global
    REAL(KIND=r8), INTENT(INOUT) :: amtransu   (npoi)     ! global
    REAL(KIND=r8), INTENT(INOUT) :: amtransl   (npoi)     ! global
    REAL(KIND=r8), INTENT(INOUT) :: amsuvap    (npoi)     ! global
    REAL(KIND=r8), INTENT(INOUT) :: aminvap    (npoi)     ! global
    INTEGER, INTENT(INOUT) :: nytimes         (npoi)           ! global! counter for yearly average calculations
    REAL(KIND=r8), INTENT(INOUT) :: aysolar        (npoi)     ! global! annual average incident solar radiation (w/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: ayirup        (npoi)     ! global! annual average upward ir radiation (w/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: ayirdown   (npoi)     ! global! annual average downward ir radiation (w/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: aysens        (npoi)     ! global! annual average sensible heat flux (w/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: aylatent   (npoi)     ! global! annual average latent heat flux (w/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: ayprcp        (npoi)     ! global! annual average precipitation (mm/yr)
    REAL(KIND=r8), INTENT(INOUT) :: ayaet        (npoi)     ! global! annual average aet (mm/yr)
    REAL(KIND=r8), INTENT(INOUT) :: aytrans        (npoi)     ! global! annual average transpiration (mm/yr)
    REAL(KIND=r8), INTENT(INOUT) :: aytrunoff  (npoi)     ! global! annual average total runoff (mm/yr)
    REAL(KIND=r8), INTENT(INOUT) :: aysrunoff  (npoi)     ! global! annual average surface runoff (mm/yr)
    REAL(KIND=r8), INTENT(INOUT) :: aydrainage (npoi)     ! global! annual average drainage (mm/yr)
    REAL(KIND=r8), INTENT(INOUT) :: aydwtot        (npoi)     ! global! annual average soil+vegetation+snow water 
    ! recharge (mm/yr or kg_h2o/m**2/yr)
    REAL(KIND=r8), INTENT(INOUT) :: aywsoi        (npoi)     ! global! annual average 1m soil moisture (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: aywisoi        (npoi)     ! global! annual average 1m soil ice (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: aytsoi        (npoi)     ! global! annual average 1m soil temperature (C)
    REAL(KIND=r8), INTENT(INOUT) :: ayvwc        (npoi)     ! global! annual average 1m volumetric water content (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: ayawc        (npoi)     ! global! annual average 1m plant-available water content (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: aystresstu (npoi)     ! global! annual average soil moisture stress 
    ! parameter for upper canopy (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: aystresstl(npoi)      ! global! annual average soil moisture stress 
    ! parameter for lower canopy (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: aygpp     (npoi,npft) ! global! annual gross npp for each plant type(kg-c/m**2/yr)
    REAL(KIND=r8), INTENT(OUT  ) :: aygpptot  (npoi)      ! local ! annual total gpp for ecosystem (kg-c/m**2/yr)
    REAL(KIND=r8), INTENT(INOUT) :: aynpp     (npoi,npft) ! global! annual total npp for each plant type(kg-c/m**2/yr)
    REAL(KIND=r8), INTENT(OUT  ) :: aynpptot  (npoi)      ! local ! annual total npp for ecosystem (kg-c/m**2/yr)
    REAL(KIND=r8), INTENT(INOUT) :: ayco2mic  (npoi)      ! global! annual total CO2 flux from microbial respiration (kg-C/m**2/yr)
    REAL(KIND=r8), INTENT(INOUT) :: ayco2root (npoi)      ! global! annual total CO2 flux from soil due to root respiration (kg-C/m**2/yr)
    REAL(KIND=r8), INTENT(OUT  ) :: ayco2soi  (npoi)      ! local ! annual total soil CO2 flux from microbial and 
    ! root respiration (kg-C/m**2/yr)
    REAL(KIND=r8), INTENT(OUT  ) :: ayneetot  (npoi)      ! local! annual total NEE for ecosystem (kg-C/m**2/yr)
    REAL(KIND=r8), INTENT(INOUT) :: ayrootbio (npoi)      ! global! annual average live root biomass (kg-C / m**2)
    REAL(KIND=r8), INTENT(INOUT) :: aynmintot (npoi)      ! global! annual total nitrogen mineralization (kg-N/m**2/yr)
    REAL(KIND=r8), INTENT(INOUT) :: ayalit    (npoi)      ! global! aboveground litter (kg-c/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: ayblit    (npoi)      ! global! belowground litter (kg-c/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: aycsoi    (npoi)      ! global! total soil carbon (kg-c/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: aycmic    (npoi)      ! global! total soil carbon in microbial biomass (kg-c/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: ayanlit   (npoi)      ! global! aboveground litter nitrogen (kg-N/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: aybnlit   (npoi)      ! global! belowground litter nitrogen (kg-N/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: aynsoi    (npoi)      ! global! total soil nitrogen (kg-N/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: ayalbedo  (npoi)  
    REAL(KIND=r8), INTENT(INOUT) :: totalit   (npoi)           ! global! total standing aboveground litter (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: totrlit   (npoi)           ! global! total root litter carbon belowground (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: totcsoi   (npoi)           ! global! total carbon in all soil pools (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: totcmic   (npoi)           ! global! total carbon residing in microbial pools (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: totanlit  (npoi)           ! global! total standing aboveground nitrogen in litter (kg_N m-2)
    REAL(KIND=r8), INTENT(INOUT) :: totrnlit  (npoi)      ! global! total root litter nitrogen belowground (kg_N m-2)
    REAL(KIND=r8), INTENT(INOUT) :: totnsoi   (npoi)      ! global! total nitrogen in soil (kg_N m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: totnmic   (npoi)      ! local! total nitrogen residing in microbial pool (kg_N m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: totlit    (npoi)      ! local! total carbon in all litter pools (kg_C m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: totfall   (npoi)      ! local! total litterfall and root turnover (kg_C m-2/year)
    REAL(KIND=r8), INTENT(OUT  ) :: totnlit   (npoi)      ! local! total nitrogen in all litter pools (kg_N m-2)

    REAL(KIND=r8), INTENT(INOUT) :: firefac   (npoi)     ! global! factor that respresents the annual average
    REAL(KIND=r8), INTENT(INOUT) :: wtot      (npoi)     ! global! total amount of water stored in snow, soil,
    ! puddels, and on vegetation (kg_h2o)
    ! fuel dryness of a grid cell, and hence characterizes the readiness to burn

    REAL(KIND=r8), INTENT(INOUT) :: storedn (npoi)        ! global ! total storage of N in soil profile (kg_N m-2) 
    REAL(KIND=r8), INTENT(INOUT) :: yrleach (npoi)        ! global ! annual total amount C leached from soil profile (kg_C m-2/yr)
    REAL(KIND=r8), INTENT(INOUT) :: ynleach (npoi)
    REAL(KIND=r8), INTENT(IN   ) :: falll   (npoi)     ! global ! annual leaf litter fall (kg_C m-2/year)
    REAL(KIND=r8), INTENT(IN   ) :: fallr   (npoi)     ! global ! annual root litter input                    (kg_C m-2/year)
    REAL(KIND=r8), INTENT(IN   ) :: fallw   (npoi)     ! global! annual wood litter fall                    (kg_C m-2/year)
    REAL(KIND=r8), INTENT(INOUT) :: clitlm  (npoi)     ! global! carbon in leaf litter pool - metabolic       (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitls  (npoi)     ! global! carbon in leaf litter pool - structural      (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitrm  (npoi)     ! global! carbon in fine root litter pool - metabolic  (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitrs  (npoi)     ! global! carbon in fine root litter pool - structural (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitwm  (npoi)     ! global! carbon in woody litter pool - metabolic      (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitws  (npoi)     ! global! carbon in woody litter pool - structural     (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: csoislop(npoi)     ! global! carbon in soil - slow protected humus           (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: csoislon(npoi)     ! global! carbon in soil - slow nonprotected humus     (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: csoipas (npoi)     ! global! carbon in soil - passive humus                   (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitll  (npoi)     ! global! carbon in leaf litter pool - lignin           (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitrl  (npoi)     ! global! carbon in fine root litter pool - lignin     (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitwl  (npoi)     ! global! carbon in woody litter pool - lignin           (kg_C m-2)


    REAL(KIND=r8), INTENT(IN   ) :: tc      (npoi)        ! global  ! coldest monthly temperature (C)
    REAL(KIND=r8), INTENT(INOUT) :: agddu   (npoi)        ! global  ! annual accumulated growing degree days for bud
    ! burst, upper canopy (day-degrees)
    REAL(KIND=r8), INTENT(INOUT) :: tempu   (npoi)        ! global  ! cold-phenology trigger for trees (non-dimensional)
    REAL(KIND=r8), INTENT(INOUT) :: agddl   (npoi)        ! global  ! annual accumulated growing degree days for bud burst,
    ! lower canopy (day-degrees)
    REAL(KIND=r8), INTENT(INOUT) :: templ   (npoi)        ! global  ! cold-phenology trigger for grasses/shrubs (non-dimensional)
    REAL(KIND=r8), INTENT(INOUT) :: dropu   (npoi)        ! global  ! drought-phenology trigger for trees (non-dimensional)
    REAL(KIND=r8), INTENT(INOUT) :: dropls  (npoi)        ! global  ! drought-phenology trigger for shrubs (non-dimensional)
    REAL(KIND=r8), INTENT(INOUT) :: dropl4  (npoi)        ! global  ! drought-phenology trigger for c4 grasses (non-dimensional)
    REAL(KIND=r8), INTENT(INOUT) :: dropl3  (npoi)        ! global  ! drought-phenology trigger for c3 grasses (non-dimensional)
    REAL(KIND=r8), INTENT(IN   ) :: plai    (npoi,npft)! global  ! total leaf area index of each plant functional type

    REAL(KIND=r8), INTENT(INOUT) :: ynleach_p (npoi)   ! annual total amount P leached from soil profile (kg_P m-2/yr)
    REAL(KIND=r8), INTENT(OUT  ) :: tnmin_p   (npoi)   ! instantaneous phosphorus mineralization         (kg_P m-2/timestep)
    REAL(KIND=r8), INTENT(OUT  ) :: totnmic_p (npoi)   ! total phosphorus residing in microbial pool     (kg_P m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: totnlit_p (npoi)   ! total phosphorus in all litter pools            (kg_P m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: totanlit_p(npoi)   ! total standing aboveground phosphorus in litter (kg_P m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: totrnlit_p(npoi)   ! total root litter phosphorus belowground        (kg_P m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: totnsoi_p (npoi)   ! total phosphorus in soil                        (kg_P m-2)
    REAL(KIND=r8), INTENT(INOUT) :: storedn_p (npoi)   ! total storage of P in soil profile              (kg_P m-2) 
    REAL(KIND=r8), INTENT(IN   ) :: adfalll   (npoi)     ! global ! annual leaf litter fall (kg_C m-2/year)
    REAL(KIND=r8), INTENT(IN   ) :: adfallr   (npoi)     ! global ! annual root litter input    (kg_C m-2/year)
    REAL(KIND=r8), INTENT(IN   ) :: adfallw   (npoi)     ! global! annual wood litter fall    (kg_C m-2/year)
    REAL(KIND=r8), INTENT(INOUT) :: adcbiol (npoi,npft)   
    REAL(KIND=r8), INTENT(INOUT) :: adcbior (npoi,npft)  
    REAL(KIND=r8), INTENT(INOUT) :: adcbiow (npoi,npft)
    REAL(KIND=r8), INTENT(IN   ) :: adplai    (npoi,npft)! global  ! total leaf area index of each plant functional type
    CHARACTER(len=*) , INTENT(IN   ) ::  rootmode
    INTEGER, INTENT(IN   ) :: nstep
    REAL(KIND=r8), INTENT(IN   ) :: beta1(nVegClass)
    REAL(KIND=r8), INTENT(IN   ) :: beta2(nVegClass)
    REAL(KIND=r8), INTENT(IN   ) :: stressfac(nVegClass)
    REAL(KIND=r8), INTENT(IN   ) :: avmuir_factor(nVegClass,2)

    !
    ! Arguments (input)
    !
    INTEGER, INTENT(IN   ) :: iday         ! day number  (passed in)
    INTEGER, INTENT(IN   ) :: imonth         ! month number (passed in)
    INTEGER, INTENT(IN   ) :: iyear
    INTEGER, INTENT(IN   ) :: iyear0
    INTEGER, INTENT(IN   ) :: isimveg 
    INTEGER, INTENT(IN   ) :: spinmax 
    REAL(KIND=r8), INTENT(IN        ) :: ux  (npoi)
    REAL(KIND=r8), INTENT(IN        ) :: uy  (npoi)
    REAL(KIND=r8), INTENT(OUT  ) :: taux(npoi)
    REAL(KIND=r8), INTENT(OUT  ) :: tauy(npoi)
    REAL(KIND=r8), INTENT(INOUT) :: ts2 (npoi)
    REAL(KIND=r8), INTENT(INOUT) :: qs2 (npoi)
    REAL(KIND=r8), INTENT(INOUT) :: gdd0this(npoi)         
    REAL(KIND=r8), INTENT(INOUT) :: gdd5this(npoi) 
    REAL(KIND=r8), INTENT(OUT  ) :: bstar(npoi)
    REAL(KIND=r8), INTENT(IN   ) :: vegtype0 (npoi)      ! annual vegetation type - ibis classification

    REAL(KIND=r8) :: tthreshold (npoi)! temperature threshold for budburst and senescence
    REAL(KIND=r8) :: gthreshold (npoi)! temperature threshold for budburst and senescence
    REAL(KIND=r8) :: avglaiu    (npoi)  ! average lai of upper canopy 
    REAL(KIND=r8) :: avglail    (npoi)! average lai of lower canopy 

    !
    INTEGER :: ib  
    INTEGER :: nsol                  ! number of points in indsol
    INTEGER :: spin

    !
    ! Compute current time step zenith angle
    !
    !      call ibiszen (calday  ,loni    ,lati    ,coszen, kpti ,kptj)


    IF (mcsec == 0.0_r8) THEN
       ! 
       ! Calculates phenology once a day (beginning of the day)
       !

       CALL pheno(tc      , &! INTENT(IN        )
            agddu   , &! INTENT(INOUT) global
            tempu   , &! INTENT(INOUT) global
            agddl   , &! INTENT(INOUT) global
            templ   , &! INTENT(INOUT) global
            dropu   , &! INTENT(INOUT) global
            dropls  , &! INTENT(INOUT) global
            dropl4  , &! INTENT(INOUT) global
            dropl3  , &! INTENT(INOUT) global
            vegtype0, &! INTENT(INOUT) global
            froot   , &! INTENT(INOUT) global
            hsoi    , &
            beta1   , &
            beta2   , &
            plai    , &! INTENT(IN        )
            adplai  , &! INTENT(IN   )
            frac    , &! INTENT(OUT  )
            lai     , &! INTENT(OUT  )
            fl      , &! INTENT(IN        )
            fu      , &! INTENT(IN        )
            zbot    , &! INTENT(OUT  )
            ztop    , &! INTENT(OUT  )
            a10td   , &! INTENT(IN        )
            a10ancub, &! INTENT(IN        )
            a10ancls, &! INTENT(IN        )
            a10ancl4, &! INTENT(IN        )
            a10ancl3, &! INTENT(IN        )
            td      , &! INTENT(IN        )
            tthreshold, &! INTENT(OUT  )
            gthreshold, &! INTENT(OUT  )
            avglaiu   , &! INTENT(OUT  )
            avglail   , &! INTENT(OUT  )
            adnpp     , &! INTENT(IN   )
            adtsoi    , &! INTENT(IN   )
            adwsoi    , &! INTENT(IN   )
            adwisoi   , &! INTENT(IN   )
            poros     , &! INTENT(IN   )
            rhow      , &! INTENT(IN   )
            npoi      , &! INTENT(IN        )
            npft      , &! INTENT(IN        )
            nsoilay    , &! INTENT(IN   )
            nVegClass, &! INTENT(IN   )
            rootmode, &! INTENT(IN   )
            epsilon   )! INTENT(IN        )

       !
       IF (isimveg .EQ. 1) THEN
          !
          ! call soil biogeochemistry model
          !

          !
          ! Soil carbon acceleration model deleted in favor of spinmax
          ! specification at each restart (AAM - 3/14/02)
          !

          !          if (soicspin .eq. 1) then
          !
          !             if ((iyear - iyear0) .le.
          !    >          (spinfrac * (nspinsoil - iyear0 - eqyears))) then
          !                spinmax = int(spincons)
          !
          !             else if ((iyear - iyear0) .lt.
          !    >              (nspinsoil - iyear0 -  eqyears)) then
          !
          !                slope   = spincons / ((nspinsoil - iyear0 - eqyears) -
          !    >                    (spinfrac * (nspinsoil - iyear0 - eqyears)))
          !
          !                spinmax = int (spincons - (slope * ((iyear - iyear0) -
          !    >                (spinfrac * (nspinsoil - iyear0 - eqyears)))))
          !
          !                spinmax = max(spinmax,1)
          !
          !             else
          !
          !                spinmax = 1
          !
          !             endif            ! if (iyear - iyear0) ....
          !
          !          else 
          !
          !             spinmax = 1
          !
          !          endif               ! if (soicspin = 1)

          DO  spin = 1, spinmax
             CALL soilbgc (iyear    , &! INTENT(IN   )
                  iyear0   , &! INTENT(IN   )
                  imonth   , &! INTENT(IN   )
                  iday     , &! INTENT(IN   )
                  spin     , &! INTENT(IN   )
                  spinmax  , &! INTENT(IN   )
                  ayprcp   , &! INTENT(IN   )
                  adfalll  , &! INTENT(IN   )
                  adfallr  , &! INTENT(IN   )
                  adfallw  , &! INTENT(IN   )
                  falll    , &! INTENT(IN   )
                  fallr    , &! INTENT(IN   )
                  fallw    , &! INTENT(IN   )
                  clitlm   , &! INTENT(INOUT)
                  clitls   , &! INTENT(INOUT)
                  clitrm   , &! INTENT(INOUT)
                  clitrs   , &! INTENT(INOUT)
                  clitwm   , &! INTENT(INOUT)
                  clitws   , &! INTENT(INOUT)
                  csoislop , &! INTENT(INOUT)
                  csoislon , &! INTENT(INOUT)
                  csoipas  , &! INTENT(INOUT)
                  totcmic  , &! INTENT(INOUT)
                  clitll   , &! INTENT(INOUT)
                  clitrl   , &! INTENT(INOUT)
                  clitwl   , &! INTENT(INOUT)
                  decomps  , &! INTENT(IN   )
                  decompl  , &! INTENT(IN   )
                  tnmin    , &! INTENT(OUT  )
                  totnmic  , &! INTENT(OUT  )
                  totlit   , &! INTENT(OUT  )
                  totalit  , &! INTENT(OUT  )
                  totrlit  , &! INTENT(OUT  )
                  totcsoi  , &! INTENT(OUT  )
                  totfall  , &! INTENT(OUT  )
                  totnlit  , &! INTENT(OUT  )
                  totanlit , &! INTENT(OUT  )
                  totrnlit , &! INTENT(OUT  )
                  totnsoi  , &! INTENT(OUT  )
                  tco2mic  , &! INTENT(OUT  )
                  storedn  , &! INTENT(INOUT)
                  yrleach  , &! INTENT(INOUT)
                  ynleach  , &! INTENT(INOUT)
                  ynleach_p ,&! INTENT(INOUT)
                  tnmin_p   ,&! INTENT(OUT  )
                  totnmic_p ,&! INTENT(OUT  )
                  totnlit_p ,&! INTENT(OUT  )
                  totanlit_p,&! INTENT(OUT  )
                  totrnlit_p,&! INTENT(OUT  )
                  totnsoi_p ,&! INTENT(OUT  )
                  storedn_p ,&! INTENT(INOUT)
                  csoi      ,&! global ! specific heat of soil, no pore spaces (J kg-1 deg-1)
                  ta        ,&! INTENT(IN   )
                  vegtype0  ,&! INTENT(IN   )
                  hsoi     , &! INTENT(IN   )
                  sand     , &! INTENT(IN   )
                  clay     , &! INTENT(IN   )
                  npoi     , &! INTENT(IN   )
                  nsoilay  , &! INTENT(IN   )
                  ndaypy     )! INTENT(IN   )

          END DO
          !
       END IF                  ! if (isimveg = 1)
       !
    END IF                    ! if (msec < dtime)

    !call the land surface model

    CALL lsxmain(ginvap       ,gsuvap       , & 
         gtrans       ,gtransu      , & 
         gtransl      ,grunof       , & 
         gdrain       ,gadjust      , & 
         a10scalparamu,a10daylightu , & 
         a10scalparaml,a10daylightl , & 
         vmax_pft     ,tau15        , & 
         kc15              ,ko15            , & 
         cimax        ,gammaub      , & 
         alpha3       ,theta3       , & 
         beta3        ,coefmub      , & 
         coefbub      ,gsubmin      , & 
         gammauc      ,coefmuc      , & 
         coefbuc      ,gsucmin      , & 
         gammals      ,coefmls      , & 
         coefbls      ,gslsmin      , & 
         gammal3      ,coefml3      , & 
         coefbl3      ,gsl3min      , & 
         gammal4      ,alpha4       , & 
         theta4       ,beta4        , & 
         coefml4      ,coefbl4      , & 
         gsl4min      ,wliqu        , & 
         wliqumax     ,wsnou        , & 
         wsnoumax     ,tu           , & 
         wliqs        ,wliqsmax     , & 
         wsnos        ,wsnosmax     , & 
         ts              ,wliql        , & 
         wliqlmax     ,wsnol        , & 
         wsnolmax     ,tl           , & 
         topparu      ,topparl      , & 
         fl              ,fu           , & 
         lai              ,sai          , & 
         rhoveg       ,tauveg       , & 
         orieh        ,oriev        , & 
         wliqmin      ,wsnomin      , & 
         t12              ,tdripu       , & 
         tblowu       ,tdrips       , & 
         tblows       ,t34            , & 
         tdripl       ,tblowl       , & 
         ztop              ,alaiml       , & 
         zbot              ,alaimu       , & 
         froot        ,q34            , & 
         q12              ,su            , & 
         cleaf        ,dleaf        , & 
         ss              ,cstem        , & 
         dstem        ,sl            , & 
         cgrass       ,ciub            , & 
         ciuc              ,exist        , & 
         csub              ,gsub            , & 
         csuc              ,gsuc            , & 
         agcub        ,agcuc        , & 
         ancub        ,ancuc        , & 
         totcondub    ,totconduc    , & 
         cils              ,cil3              , & 
         cil4              ,csls              , & 
         gsls              ,csl3              , & 
         gsl3              ,csl4              , & 
         gsl4              ,agcls        , & 
         agcl4        ,agcl3        , & 
         ancls        ,ancl4        , & 
         ancl3        ,totcondls    , & 
         totcondl3    ,totcondl4    , & 
         chu              ,chs              , & 
         chl              ,frac              , & 
         tlsub        , & 
         z0sno        , & 
         rhos              , & 
         consno       , & 
         hsnotop      , & 
         hsnomin      , & 
         fimin        , & 
         fimax        , & 
         fi              , & 
         tsno              , & 
         hsno              , & 
         sand              , & 
         clay              , & 
         poros        , & 
         wsoi              , & 
         wisoi        , & 
         consoi       , & 
         zwpmax       , & 
         wpud              , & 
         wipud        , & 
         wpudmax      , & 
         qglif        , & 
         tsoi              , & 
         hvasug       , & 
         hvasui       , & 
         albsav       , & 
         albsan       , & 
         tg              , & 
         ti              , & 
         z0soi        , & 
         swilt        , & 
         sfield       , & 
         stressl      , & 
         stressu      , & 
         stresstl     , & 
         stresstu     , & 
         csoi              , & 
         rhosoi       , & 
         hsoi              , & 
         suction      , & 
         bex              , & 
         upsoiu       , & 
         upsoil       , & 
         heatg        , & 
         heati        , & 
         hydraul      , & 
         porosflo     , & 
         ibex              , & 
         bperm        , & 
         hflo              , & 
         ta              , & 
         asurd        , & 
         asuri        , & 
         coszen       , & 
         solad        , & 
         solai        , & 
         fira              , & 
         raina        , & 
         qa              , & 
         psurf        , & 
         snowa        , & 
         ua              , & 
         o2conc       , & 
         co2conc      , & 
         fwetu        , & 
         rliqu        , & 
         fwets        , & 
         rliqs        , & 
         fwetl        , & 
         rliql        , & 
                                !nsol              , & 
         solu              , & 
         sols              , & 
         soll              , & 
         solg              , & 
         soli              , & 
         scalcoefl    , & 
         scalcoefu    , & 
         indsol       , & 
         albsod       , & 
         albsoi       , & 
         albsnd       , & 
         albsni       , & 
         relod        , & 
         reloi        , & 
         reupd        , & 
         reupi        , & 
         ablod        , & 
         abloi        , & 
         flodd        , & 
         dummy        , & 
         flodi        , & 
         floii        , & 
         terml        , & 
         termu        , & 
         abupd        , & 
         abupi        , & 
         fupdd        , & 
         fupdi        , & 
         fupii        , & 
         sol2d        , & 
         sol2i        , & 
         sol3d        , & 
         sol3i        , & 
         firb              , & 
         firs              , & 
         firu              , & 
         firl              , & 
         firg              , & 
         firi              , & 
         snowg        , & 
         tsnowg       , & 
         tsnowl       , & 
         pfluxl       , & 
         raing        , & 
         traing       , & 
         trainl       , & 
         snowl        , & 
         tsnowu       , & 
         pfluxu       , & 
         rainu        , & 
         trainu       , & 
         snowu        , & 
         pfluxs       , & 
         rainl        , & 
         bps              , & 
         rhoa              , & 
         cp              , & 
         za              , & 
         bdl              , & 
         dil              , & 
         z3              , & 
         z4              , & 
         z34              , & 
         exphl        , & 
         expl              , & 
         displ        , & 
         bdu              , & 
         diu              , & 
         z1              , & 
         z2              , & 
         z12              , & 
         exphu        , & 
         expu              , & 
         dispu        , & 
         alogg        , & 
         alogi        , & 
         alogav       , & 
         alog4        , & 
         alog3        , & 
         alog2        , & 
         alog1        , & 
         aloga        , & 
         u2              , & 
         alogu        , & 
         alogl        , & 
         richl        , & 
         straml       , & 
         strahl       , & 
         richu        , & 
         stramu       , & 
         strahu       , & 
         u1              , & 
         u12              , & 
         u3              , & 
         u34              , & 
         u4              , & 
         cu              , & 
         cl              , & 
         sg              , & 
         si              , & 
         fwetux       , & 
         fwetsx       , & 
         fwetlx       , & 
         fsena        , & 
         fseng        , & 
         fseni        , & 
         fsenu        , & 
         fsens        , & 
         fsenl        , & 
         fvapa        , & 
         fvaput       , & 
         fvaps        , & 
         fvaplw       , & 
         fvaplt       , & 
         fvapg        , & 
         fvapi        , & 
         fvapuw       , & 
         npoi              , & 
         nband        , & 
         nsoilay      , & 
         nsnolay      , & 
         npft              , & 
         epsilon      , & 
         dtime        , & 
         stef              , & 
         vonk              , & 
         grav              , & 
         tmelt        , & 
         hfus              , & 
         hvap              , & 
         hsub              , & 
         ch2o              , & 
         cice              , & 
         cair              , & 
         cvap              , & 
         rair              , & 
         rvap              , & 
         cappa        , & 
         rhow              , & 
         vzero        , & 
         pi              , &
         ux              , &! INTENT(IN        ) !global
         uy              , &! INTENT(IN        ) !global
         taux              , &! INTENT(OUT        ) !local
         tauy              , &! INTENT(OUT        ) !local
         bstar        , &! INTENT(OUT        ) !local
         ts2          , &! INTENT(OUT        ) !local
         qs2          , & ! INTENT(OUT        ) !local
         vegtype0     , &! INTENT(IN   ) !local
         stressfac    , &
	 avmuir_factor, &
         nVegClass     )
    !
    ! accumulate some variables every timestep
    !
    CALL  sumnow(a10td   , &! INTENT(INOUT) !global
         a10ancub, &! INTENT(INOUT) !global
         a10ancuc, &! INTENT(INOUT) !global
         a10ancls, &! INTENT(INOUT) !global
         a10ancl3, &! INTENT(INOUT) !global
         a10ancl4, &! INTENT(INOUT) !global
         nppdummy, &! INTENT(OUT  ) !local
         frac    , &! INTENT(IN   ) !global
         ancub   , &! INTENT(IN   ) !global
         lai           , &! INTENT(IN   ) !global
         fu           , &! INTENT(IN   ) !global
         ancuc   , &! INTENT(IN   ) !global
         ancls   , &! INTENT(IN   ) !global
         fl           , &! INTENT(IN   ) !global
         ancl4   , &! INTENT(IN   ) !global
         ancl3   , &! INTENT(IN   ) !global
         tgpp    , &! INTENT(OUT  ) !local
         agcub   , &! INTENT(IN   ) !global
         agcuc   , &! INTENT(IN   ) !global
         agcls   , &! INTENT(IN   ) !global
         agcl4   , &! INTENT(IN   ) !global
         agcl3   , &! INTENT(IN   ) !global
         tgpptot , &! INTENT(OUT  ) !local
         ts      , &! INTENT(IN   ) !global
         froot   , &! INTENT(IN   ) !global
         tnpp           , &! INTENT(OUT  ) !local
         cbiow   , &! INTENT(IN   ) !global
         sapfrac , &! INTENT(IN   ) !global
         cbior   , &! INTENT(IN   ) !global
         tnpptot , &! INTENT(OUT  ) !local
         tco2root, &! INTENT(OUT  ) !local
         tneetot , &! INTENT(OUT  ) !local
         tco2mic , &! INTENT(IN   ) !global
         tsoi           , &! INTENT(IN   ) !global
         fi           , &! INTENT(IN   ) !global
         td           , &! INTENT(IN   ) !global
         npoi    , &! INTENT(IN   ) !global
         nsoilay , &! INTENT(IN   ) !global
         npft           , &! INTENT(IN   ) !global
         ndaypy  , &! INTENT(IN   ) !global
         dtime     )! INTENT(IN   ) !global

    CALL sumday(adnpp     , &
         tnpp      , &! INTENT(IN  ) !local
         raina     , &! INTENT(IN   )
         snowa     , &! INTENT(IN   )
         fvapa     , &! INTENT(IN   )
         grunof    , &! INTENT(IN   )
         gdrain    , &! INTENT(IN   )
         hsno      , &! INTENT(IN   )
         fi        , &! INTENT(IN   )
         hsoi      , &! INTENT(IN   )
         tsoi      , &! INTENT(IN   )
         wsoi      , &! INTENT(IN   )
         wisoi     , &! INTENT(IN   )
         ndtimes   , &! INTENT(INOUT) global
         adrain    , &! INTENT(INOUT) global
         adsnow    , &! INTENT(INOUT) global
         adaet     , &! INTENT(INOUT) global
         adtrunoff , &! INTENT(INOUT) global
         adsrunoff , &! INTENT(INOUT) global
         addrainage, &! INTENT(INOUT) global
         adrh      , &! INTENT(INOUT) global
         adsnod    , &! INTENT(INOUT) global
         adsnof    , &! INTENT(INOUT) global
         adwsoi    , &! INTENT(INOUT) global
         adtsoi    , &! INTENT(INOUT) global
         adwisoi   , &! INTENT(INOUT) global
         adtlaysoi , &! INTENT(INOUT) global
         adwlaysoi , &! INTENT(INOUT) global
         adwsoic   , &! INTENT(INOUT) global
         adtsoic   , &! INTENT(INOUT) global
         adco2mic  , &! INTENT(INOUT) global
         adco2root , &! INTENT(INOUT) global
         adco2soi  , &! INTENT(INOUT) global
         adco2ratio, &! INTENT(INOUT) global
         adnmintot , &! INTENT(INOUT) global
         froot     , &! INTENT(IN   )
         tco2mic   , &! INTENT(IN   )
         tco2root  , &! INTENT(IN   )
         decompl   , &! INTENT(INOUT) global
         decomps   , &! INTENT(INOUT) global
         tnmin     , &! INTENT(IN   )
         npoi      , &! INTENT(IN   )
         npft      , &! INTENT(IN   ) !global
         nsoilay   , &! INTENT(IN   )
         nsnolay   , &! INTENT(IN   )
         dtime     , &! INTENT(IN   )
         td            , &! INTENT(IN   )
         gdd0this  , &! INTENT(INOUT) global
         gdd5this  , &! INTENT(INOUT) global
         ts2       , &! INTENT(INOUT) global
         mcsec             ) ! INTENT(INOUT) global

    CALL summonth (dtime     , &! INTENT(IN   )!global
         mcsec     , &! INTENT(IN   )!global
         iday      , &! INTENT(IN   )!global
         imonth    , &! INTENT(IN   )!global
         nmtimes   , &! INTENT(INOUT)!global
         amrain    , &! INTENT(INOUT)!global
         amsnow    , &! INTENT(INOUT)!global
         amaet     , &! INTENT(INOUT)!global
         amtrunoff , &! INTENT(INOUT)!global
         amsrunoff , &! INTENT(INOUT)!global
         amdrainage, &! INTENT(INOUT)!global
         amtemp    , &! INTENT(INOUT)!global
         amqa            , &! INTENT(INOUT)!global
         amsolar   , &! INTENT(INOUT)!global
         amirup    , &! INTENT(INOUT)!global
         amirdown  , &! INTENT(INOUT)!global
         amsens    , &! INTENT(INOUT)!global
         amlatent  , &! INTENT(INOUT)!global
         amlaiu    , &! INTENT(INOUT)!global
         amlail    , &! INTENT(INOUT)!global
         amtsoi    , &! INTENT(INOUT)!global
         amwsoi    , &! INTENT(INOUT)!global
         amwisoi   , &! INTENT(INOUT)!global
         amvwc     , &! INTENT(INOUT)!global
         amawc     , &! INTENT(INOUT)!global
         amsnod    , &! INTENT(INOUT)!global
         amsnof    , &! INTENT(INOUT)!global
         amnpp            , &! INTENT(INOUT)!global
         amnpptot  , &! INTENT(OUT  )!local
         amco2mic  , &! INTENT(INOUT)!global
         amco2root , &! INTENT(INOUT)!global
         amco2soi  , &! INTENT(OUT  )!local
         amco2ratio, &! INTENT(OUT  )!local
         amneetot  , &! INTENT(OUT  )!local
         amnmintot , &! INTENT(INOUT)!global
         amts2     , &! INTENT(INOUT)!global
         amtransu  , &! INTENT(INOUT)!global
         amtransl  , &! INTENT(INOUT)!global
         amsuvap   , &! INTENT(INOUT)!global
         aminvap   , &! INTENT(INOUT)!global
         amalbedo  , &! INTENT(INOUT)!global
         amtsoil   , &! INTENT(INOUT)!global
         amwsoil   , &! INTENT(INOUT)!global
         amwisoil  , &! INTENT(INOUT)!global
         ts2       , &! INTENT(INOUT)!global
         fu        , &! INTENT(IN   )!global
         lai       , &! INTENT(IN   )!global
         fl            , &! INTENT(IN   )!global
         tnpp      , &! INTENT(IN   )!global
         tco2mic   , &! INTENT(IN   )!global
         tco2root  , &! INTENT(IN   )!global
         tnmin     , &! INTENT(IN   )!global
         hsoi      , &! INTENT(IN   )!global
         tsoi      , &! INTENT(IN   )!global
         wsoi      , &! INTENT(IN   )!global
         wisoi            , &! INTENT(IN   )!global
         poros     , &! INTENT(IN   )!global
         swilt            , &! INTENT(IN   )!global
         hsno      , &! INTENT(IN   )!global
         fi        , &! INTENT(IN   )!global
         grunof    , &! INTENT(IN   )!global
         gdrain    , &! INTENT(IN   )!global
         gtransu   , &! INTENT(IN   )!global
         gtransl   , &! INTENT(IN   )!global
         gsuvap    , &! INTENT(IN   )!global
         ginvap    , &! INTENT(IN   )!global
         asurd     , &! INTENT(IN   )!global
         asuri     , &! INTENT(IN   )!global
         fvapa     , &! INTENT(IN   )!global
         firb            , &! INTENT(IN   )!global
         fsena     , &! INTENT(IN   )!global
         raina     , &! INTENT(IN   )!global
         snowa     , &! INTENT(IN   )!global
         ta        , &! INTENT(IN   )!global
         qa        , &! INTENT(IN   )!global
         solad            , &! INTENT(IN   )!global
         solai     , &! INTENT(IN   )!global
         fira      , &! INTENT(IN   )!global
         npoi      , &! INTENT(IN   )!global
         nband            , &! INTENT(IN   )!global
         nsoilay   , &! INTENT(IN   )!global
         nsnolay   , &! INTENT(IN   )!global
         npft      , &! INTENT(IN   )!global
         ndaypm    , &! INTENT(IN   )!global
         hvap        )! INTENT(IN   )!global

    CALL sumyear(dtime     , &! INTENT(IN   )
         mcsec     , &! INTENT(IN   )
         iday      , &! INTENT(IN   )
         imonth    , &! INTENT(IN   )
         wliqu     , &! INTENT(IN   )
         wsnou     , &! INTENT(IN   )
         fu            , &! INTENT(IN   )
         lai            , &! INTENT(IN   )
         wliqs     , &! INTENT(IN   )
         wsnos     , &! INTENT(IN   )
         sai       , &! INTENT(IN   )
         wliql            , &! INTENT(IN   )
         wsnol     , &! INTENT(IN   )
         fl        , &! INTENT(IN   )
         tgpp      , &! INTENT(IN   )
         tnpp      , &! INTENT(IN   )
         firefac   , &! INTENT(INOUT) global
         tco2mic   , &! INTENT(IN   )
         tco2root  , &! INTENT(IN   )
         cbior     , &! INTENT(IN   )
         tnmin     , &! INTENT(IN   )
         totalit   , &! INTENT(IN   )
         totrlit   , &! INTENT(IN   )
         totcsoi   , &! INTENT(IN   )
         totcmic   , &! INTENT(IN   )
         totanlit  , &! INTENT(IN   )
         totrnlit  , &! INTENT(IN   )
         totnsoi   , &! INTENT(IN   )
         nytimes   , &! INTENT(INOUT) global
         aysolar   , &! INTENT(INOUT) global
         ayirup    , &! INTENT(INOUT) global
         ayirdown  , &! INTENT(INOUT) global
         aysens    , &! INTENT(INOUT) global
         aylatent  , &! INTENT(INOUT) global
         ayprcp    , &! INTENT(INOUT) global
         ayaet     , &! INTENT(INOUT) global
         aytrans   , &! INTENT(INOUT) global
         aytrunoff , &! INTENT(INOUT) global
         aysrunoff , &! INTENT(INOUT) global
         aydrainage, &! INTENT(INOUT) global
         aydwtot   , &! INTENT(INOUT) global
         aywsoi    , &! INTENT(INOUT) global
         aywisoi   , &! INTENT(INOUT) global
         aytsoi    , &! INTENT(INOUT) global
         ayvwc     , &! INTENT(INOUT) global
         ayawc     , &! INTENT(INOUT) global
         aystresstu, &! INTENT(INOUT) global
         aystresstl, &! INTENT(INOUT) global
         aygpp     , &! INTENT(INOUT) global
         aygpptot  , &! INTENT(OUT  ) local
         aynpp            , &! INTENT(INOUT) global
         aynpptot  , &! INTENT(OUT  ) local
         ayco2mic  , &! INTENT(INOUT) global
         ayco2root , &! INTENT(INOUT) global
         ayco2soi  , &! INTENT(OUT  ) global
         ayneetot  , &! INTENT(OUT  ) global
         ayrootbio , &! INTENT(INOUT) global
         aynmintot , &! INTENT(INOUT) global
         ayalit    , &! INTENT(INOUT) global
         ayblit    , &! INTENT(INOUT) global
         aycsoi    , &! INTENT(INOUT) global
         aycmic    , &! INTENT(INOUT) global
         ayanlit   , &! INTENT(INOUT) global
         aybnlit   , &! INTENT(INOUT) global
         aynsoi    , &! INTENT(INOUT) global
         ayalbedo  , &! INTENT(INOUT) global
         hsoi            , &! INTENT(IN   ) global
         wpud            , &! INTENT(IN   ) global
         wipud     , &! INTENT(IN   ) global
         poros     , &! INTENT(IN   ) global
         wsoi            , &! INTENT(IN   ) global
         wisoi     , &! INTENT(IN   ) global
         tsoi            , &! INTENT(IN   ) global
         swilt     , &! INTENT(IN   ) global
         stresstu  , &! INTENT(IN   ) global
         stresstl  , &! INTENT(IN   ) global
         fi            , &! INTENT(IN   ) global
         rhos      , &! INTENT(IN   ) global
         hsno            , &! INTENT(IN   ) global
         gtrans    , &! INTENT(IN   ) global
         grunof    , &! INTENT(IN   ) global
         gdrain    , &! INTENT(IN   ) global
         wtot            , &! INTENT(INOUT) global
         firb            , &! INTENT(IN   ) global
         fsena     , &! INTENT(IN   ) global
         fvapa     , &! INTENT(IN   ) global
         solad            , &! INTENT(IN   ) global
         solai     , &! INTENT(IN   ) global
         fira      , &! INTENT(IN   ) global
         raina     , &! INTENT(IN   ) global
         snowa            , &! INTENT(IN   ) global
         asurd     , &! INTENT(IN   ) global
         asuri     , &! INTENT(IN   ) global
         npoi            , &! INTENT(IN   ) global
         nband     , &! INTENT(IN   ) global
         nsoilay   , &! INTENT(IN   ) global
         nsnolay   , &! INTENT(IN   ) global
         npft      , &! INTENT(IN   ) global
         ndaypy    , &! INTENT(IN   ) global
         hvap      , &! INTENT(IN   ) global
         rhow        )! INTENT(IN   ) global
    IF (doalb) THEN
       !
       ! Compute next time step zenith angle
       !
       !         CALL ibiszen (calday1, loni, lati, coszen, kpti ,kptj)
       !
       ! Compute albedos (used in next step radiation computations)
       ! set up for solar calculations
       !
       !         CALL solset(loopi, kpti, kptj)
       !
       ! set up for solar calculations
       !
       CALL solset(npoi     , &! INTENT(IN   )
            nsol     , &! INTENT(OUT  )
            nband    , &! INTENT(IN   )
            solu     , &! INTENT(OUT  )
            sols     , &! INTENT(OUT  )
            soll     , &! INTENT(OUT  )
            solg     , &! INTENT(OUT  )
            soli     , &! INTENT(OUT  )
            scalcoefl, &! INTENT(OUT  )
            scalcoefu, &! INTENT(OUT  )
            indsol   , &! INTENT(OUT  )
            topparu  , &! INTENT(OUT  )
            topparl  , &! INTENT(OUT  )
            asurd    , &! INTENT(OUT  )
            asuri    , &! INTENT(OUT  )
            coszen     )! INTENT(IN   )  

       !
       ! solar calculations for each waveband
       !
       DO  ib = 1, nband
          !
          ! solsur sets surface albedos for soil and snow
          ! solalb performs the albedo calculations
          ! solarf uses the unit-incident-flux results from solalb
          ! to obtain absorbed fluxes sol[u,s,l,g,i] and 
          ! incident pars sunp[u,l]
          !
          CALL solsur (ib       , &! INTENT(IN   )
               tmelt    , &! INTENT(IN   )
               nsol     , &! INTENT(IN   )
               albsod   , &! INTENT(OUt  )
               albsoi   , &! INTENT(OUt  )
               albsnd   , &! INTENT(OUt  )
               albsni   , &! INTENT(OUt  )
               indsol   , &! INTENT(IN   )
               wsoi     , &! INTENT(IN   )
               wisoi    , &! INTENT(IN   )
               albsav   , &! INTENT(IN   )
               albsan   , &! INTENT(IN   )
               tsno     , &! INTENT(IN   )
               coszen   , &! INTENT(IN   )
               npoi     , &! INTENT(IN   )
               nsoilay  , &! INTENT(IN   )
               nsnolay    )! INTENT(IN   )

          CALL solalb (ib       , &! INTENT(IN   )
               nVegClass , &! INTENT(IN   )
               vegtype0     , &! INTENT(IN   )
               avmuir_factor, &! INTENT(IN   )
               relod    , &! INTENT(OUT  )
               reloi    , &! INTENT(OUT  )
               indsol   , &! INTENT(IN   )
               reupd    , &! INTENT(OUT  )
               reupi    , &! INTENT(OUT  )
               albsnd   , &! INTENT(IN   )
               albsni   , &! INTENT(IN   )
               albsod   , &! INTENT(IN   )
               albsoi   , &! INTENT(IN   )
               fl       , &! INTENT(IN   )
               fu       , &! INTENT(IN   )
               fi       , &! INTENT(IN   )
               asurd    , &! INTENT(INOUT)! local
               asuri    , &! INTENT(INOUT)! local
               npoi     , &! INTENT(IN   )
               nband    , &! INTENT(IN   )
               nsol     , &! INTENT(IN   )
               ablod    , &! INTENT(OUT  )
               abloi    , &! INTENT(OUT  )
               flodd    , &! INTENT(OUT  )
               dummy    , &! INTENT(OUT  )
               flodi    , &! INTENT(OUT  )
               floii    , &! INTENT(OUT  )
               coszen   , &! INTENT(IN   )
               terml    , &! INTENT(OUT  )
               termu    , &! INTENT(OUT  )
               lai      , &! INTENT(IN   )
               sai      , &! INTENT(IN   )
               abupd    , &! INTENT(OUT  )
               abupi    , &! INTENT(OUT  )
               fupdd    , &! INTENT(OUT  )
               fupdi    , &! INTENT(OUT  )
               fupii    , &! INTENT(OUT  )
               fwetl    , &! INTENT(IN   )
               rliql    , &! INTENT(IN   )
               rliqu    , &! INTENT(IN   )
               rliqs    , &! INTENT(IN   )
               fwetu    , &! INTENT(IN   )
               fwets    , &! INTENT(IN   )
               rhoveg   , &! INTENT(IN   )
               tauveg   , &! INTENT(IN   )
               orieh    , &! INTENT(IN   )
               oriev    , &! INTENT(IN   )
               tl       , &! INTENT(IN   )
               ts       , &! INTENT(IN   )
               tu       , &! INTENT(IN   )
               pi       , &! INTENT(IN   )
               tmelt    , &! INTENT(IN   )
               epsilon    )! INTENT(IN   )

          CALL solarf (ib       , & ! INTENT(IN        ) 
               nsol     , & ! INTENT(IN        ) 
               solu     , & ! INTENT(INOUT) !global
               indsol   , & ! INTENT(IN        ) 
               abupd    , & ! INTENT(IN        ) 
               abupi    , & ! INTENT(IN        ) 
               sols     , & ! INTENT(INOUT) !global
               sol2d    , & ! INTENT(OUT  ) 
               fupdd    , & ! INTENT(IN        ) 
               sol2i    , & ! INTENT(OUT  ) 
               fupii    , & ! INTENT(IN        ) 
               fupdi    , & ! INTENT(IN        ) 
               sol3d    , & ! INTENT(OUT  ) 
               sol3i    , & ! INTENT(OUT  ) 
               soll     , & ! INTENT(INOUT) !global
               ablod    , & ! INTENT(IN        ) 
               abloi    , & ! INTENT(IN        ) 
               flodd    , & ! INTENT(IN        ) 
               flodi    , & ! INTENT(IN        ) 
               floii    , & ! INTENT(IN        ) 
               solg     , & ! INTENT(INOUT) !global
               albsod   , & ! INTENT(IN        ) 
               albsoi   , & ! INTENT(IN        ) 
               soli     , & ! INTENT(INOUT) !global
               albsnd   , & ! INTENT(IN        ) 
               albsni   , & ! INTENT(IN        ) 
               scalcoefu, & ! INTENT(OUT  ) 
               termu    , & ! INTENT(IN        ) 
               scalcoefl, & ! INTENT(OUT  ) 
               terml    , & ! INTENT(IN        ) 
               lai      , & ! INTENT(IN        ) 
               sai      , & ! INTENT(IN        ) 
               fu       , & ! INTENT(IN        )
               fl       , & ! INTENT(IN        )
               topparu  , & ! INTENT(OUT  ) 
               topparl  , & ! INTENT(OUT  ) 
               solad    , & ! INTENT(IN        )
               solai    , & ! INTENT(IN        )
               npoi     , & ! INTENT(IN        ) 
               nband    , & ! INTENT(IN        ) 
               epsilon    ) ! INTENT(IN        ) 
          !
       END DO
       !
    END IF

  END SUBROUTINE Ibis






  SUBROUTINE Albedo_IBIS( &
       ! Model information
       latco          ,nCols          ,kmax           , &
       imask          , &
       ! Model Geometry
       zenith         , &
       ! Time info
       month2         ,month          , &
       ! Atmospheric fields
       wind           , tsea          , &
       ! Microphysics
       taud           , &
       ! LW Radiation fields at last integer hour
       LwSfcDown      , &
       ! Radiation field (Interpolated) at time = tod
       xVisBeam       ,xVisDiff       ,xNirBeam       , &
       xNirDiff       , &
       ! Surface Albedo
       avisb          ,avisd          ,anirb          , &
       anird )

    IMPLICIT NONE
    ! Model information
    INTEGER         , INTENT(IN   ) :: latco
    INTEGER         , INTENT(IN   ) :: nCols
    INTEGER         , INTENT(IN   ) :: kmax
    INTEGER(KIND=i8), INTENT(IN   ) :: imask  (ncols)
    ! Model Geometry
    REAL   (KIND=r8), INTENT(IN   ) :: zenith (nCols)       ! cosine of solar zenith angle   
    ! Time info
    INTEGER         , INTENT(INOUT) :: month2 (ncols)
    INTEGER         , INTENT(IN   ) :: month  (ncols)
    ! Atmospheric fields
    REAL   (KIND=r8), INTENT(IN   ) :: wind   (ncols)!wind speed in m/s
    REAL   (KIND=r8), INTENT(IN   ) :: tsea   (nCols)       ! cosine of solar zenith angle   
    ! Microphysics
    REAL   (KIND=r8), INTENT(IN   ) :: taud(ncols,kMax)
    ! LW Radiation fields at last integer hour
    REAL   (KIND=r8), INTENT(IN   ) :: LwSfcDown(1:nCols)
    ! Radiation field (Interpolated) at time = tod
    REAL   (KIND=r8), INTENT(IN   ) :: xVisBeam (1:nCols)
    REAL   (KIND=r8), INTENT(IN   ) :: xVisDiff (1:nCols)
    REAL   (KIND=r8), INTENT(IN   ) :: xNirBeam (1:nCols)
    REAL   (KIND=r8), INTENT(IN   ) :: xNirDiff (1:nCols)
    ! Surface Albedo
    REAL   (KIND=r8), INTENT(OUT  ) :: avisb (ncols)
    REAL   (KIND=r8), INTENT(OUT  ) :: avisd (ncols)
    REAL   (KIND=r8), INTENT(OUT  ) :: anirb (ncols)
    REAL   (KIND=r8), INTENT(OUT  ) :: anird (ncols)

    INTEGER       :: nsol         ! number of points in indsol
    REAL(KIND=r8) :: solu   (nCols)! solar flux (direct + diffuse) absorbed by upper canopy leaves per unit canopy area (W m-2)
    REAL(KIND=r8) :: sols   (nCols)! solar flux (direct + diffuse) absorbed by upper canopy stems per unit canopy area (W m-2)
    REAL(KIND=r8) :: soll   (nCols)! solar flux (direct + diffuse) absorbed by lower canopy leaves and stems per unit canopy area (W m-2)
    REAL(KIND=r8) :: solg   (nCols)! solar flux (direct + diffuse) absorbed by unit snow-free soil (W m-2)
    REAL(KIND=r8) :: soli   (nCols)! solar flux (direct + diffuse) absorbed by unit snow surface (W m-2)
    REAL(KIND=r8) :: scalcoefl(nCols,4)   ! term needed in lower canopy scaling
    REAL(KIND=r8) :: scalcoefu(nCols,4)   ! term needed in upper canopy scaling
    INTEGER       :: indsol (nCols)         ! index of current strip for points with positive coszen
    REAL(KIND=r8) :: topparu(nCols)        ! total photosynthetically active raditaion absorbed by top leaves of upper canopy (W m-2)
    REAL(KIND=r8) :: topparl(nCols)        ! total photosynthetically active raditaion absorbed by top leaves of lower canopy (W m-2)
    REAL(KIND=r8) :: albsod (nCols)          ! direct  albedo for soil surface (visible or IR)
    REAL(KIND=r8) :: albsoi (nCols)          ! diffuse albedo for soil surface (visible or IR)
    REAL(KIND=r8) :: albsnd (nCols)          ! direct  albedo for snow surface (visible or IR)
    REAL(KIND=r8) :: albsni (nCols)          ! diffuse albedo for snow surface (visible or IR)
    REAL(KIND=r8) :: relod  (nCols)         ! upward direct radiation per unit icident direct beam on lower canopy (W m-2)
    REAL(KIND=r8) :: reloi  (nCols)         ! upward diffuse radiation per unit incident diffuse radiation on lower canopy (W m-2)
    REAL(KIND=r8) :: reupd  (nCols)         ! upward direct radiation per unit incident direct radiation on upper canopy (W m-2)
    REAL(KIND=r8) :: reupi  (nCols)         ! upward diffuse radiation per unit incident diffuse radiation on upper canopy (W m-2)
    REAL(KIND=r8) :: ablod  (nCols)          ! fraction of direct  radiation absorbed by lower canopy
    REAL(KIND=r8) :: abloi  (nCols)          ! fraction of diffuse radiation absorbed by lower canopy
    REAL(KIND=r8) :: flodd  (nCols)          ! downward direct radiation per unit incident direct radiation on lower canopy (W m-2)
    REAL(KIND=r8) :: dummy  (nCols)          ! placeholder, always = 0: no direct flux produced for diffuse incident
    REAL(KIND=r8) :: flodi  (nCols)         ! downward diffuse radiation per unit incident direct radiation on lower canopy (W m-2)
    REAL(KIND=r8) :: floii  (nCols)         ! downward diffuse radiation per unit incident diffuse radiation on lower canopy
    REAL(KIND=r8) :: terml  (nCols,7)          ! term needed in lower canopy scaling
    REAL(KIND=r8) :: termu  (nCols,7)          ! term needed in upper canopy scaling
    REAL(KIND=r8) :: abupd  (nCols)       ! fraction of direct  radiation absorbed by upper canopy
    REAL(KIND=r8) :: abupi  (nCols)         ! fraction of diffuse radiation absorbed by upper canopy
    REAL(KIND=r8) :: fupdd  (nCols)         ! downward direct radiation per unit incident direct beam on upper canopy (W m-2)
    REAL(KIND=r8) :: fupdi  (nCols)         ! downward diffuse radiation per unit icident direct radiation on upper canopy (W m-2)
    REAL(KIND=r8) :: fupii  (nCols)         ! downward diffuse radiation per unit incident diffuse radiation on upper canopy (W m-2)
    REAL(KIND=r8) :: fwetu  (nCols)         ! fraction of upper canopy leaf area wetted by intercepted liquid and/or snow
    REAL(KIND=r8) :: rliqu  (nCols)         ! proportion of fwetu due to liquid
    REAL(KIND=r8) :: fwets  (nCols)         ! fraction of upper canopy stem area wetted by intercepted liquid and/or snow
    REAL(KIND=r8) :: rliqs  (nCols)         ! proportion of fwets due to liquid
    REAL(KIND=r8) :: fwetl  (nCols)         ! fraction of lower canopy stem & leaf area wetted by intercepted liquid and/or snow
    REAL(KIND=r8) :: rliql  (nCols)         ! proportion of fwetl due to liquid
    REAL(KIND=r8) :: coszen (nCols)       ! cosine of solar zenith angle
    REAL(KIND=r8) :: f
    REAL(KIND=r8) :: ocealb
    REAL(KIND=r8) :: IceOceanAlb (nCols,2,2)  
    REAL(KIND=r8), PARAMETER :: tice=271.16_r8

!     --------------------------- INPUT ---------------------------------------
!   
!    specify the parameters for albedo here:

    REAL(KIND=r8) ::         tau (nCols)             !aerosol/cloud optical depth
    REAL(KIND=r8) ::         chl  (nCols)            !chlorophyll concentration in mg/m3
    
    INTEGER       :: ncount,xband
    INTEGER       :: ib  ,i, npoi,k
!     --------------------------- Initialization ---------------------------------------


!     --------------------------- LOCAL ---------------------------------------


    npoi=0
    DO i=1,nCols
       IF (iMask(i) >= 1_i8) THEN
          npoi=npoi+1 
          coszen(npoi) = zenith(i)
       END IF
    END DO

    !
    ! calculate areal fractions wetted by intercepted h2o
    !
    IF(npoi>=1)THEN
    CALL fwetcal(npoi               , &! INTENT(IN   )
         fwetu (1:npoi)        , &! INTENT(OUT  )
         rliqu (1:npoi)        , &! INTENT(OUT  )
         fwets (1:npoi)        , &! INTENT(OUT  )
         rliqs (1:npoi)        , &! INTENT(OUT  )
         fwetl (1:npoi)        , &! INTENT(OUT  )
         rliql (1:npoi)        , &! INTENT(OUT  )
         wliqu (1:npoi,latco)  , &! INTENT(IN   )
         wliqumax              , &! INTENT(IN   ) ::
         wsnou (1:npoi,latco)  , &! INTENT(IN   ) ::
         wsnoumax              , &! INTENT(IN   ) ::
         tu   (1:npoi,latco)   , &! INTENT(IN   )
         wliqs(1:npoi,latco)   , &! INTENT(IN   )
         wliqsmax              , &! INTENT(IN   )
         wsnos(1:npoi,latco)   , &! INTENT(IN   )
         wsnosmax              , &! INTENT(IN   )
         ts   (1:npoi,latco)   , &! INTENT(IN   )
         wliql(1:npoi,latco)   , &! INTENT(IN   )
         wliqlmax              , &! INTENT(IN   )
         wsnol(1:npoi,latco)   , &! INTENT(IN   )
         wsnolmax              , &! INTENT(IN   )
         tl   (1:npoi,latco)   , &! INTENT(IN   )
         epsilon               , &! INTENT(IN   )
         tmelt                 )  ! INTENT(IN   )
    !
    ! set up for solar calculations
    !
    CALL solset(npoi                 , &! INTENT(IN   )
         nsol                 , &! INTENT(OUT  )
         nband                , &! INTENT(IN   )
         solu     (1:npoi)    , &! INTENT(OUT  )
         sols     (1:npoi)    , &! INTENT(OUT  )
         soll     (1:npoi)    , &! INTENT(OUT  )
         solg     (1:npoi)    , &! INTENT(OUT  )
         soli     (1:npoi)    , &! INTENT(OUT  )
         scalcoefl(1:npoi,1:4), &! INTENT(OUT  )
         scalcoefu(1:npoi,1:4), &! INTENT(OUT  )
         indsol   (1:npoi)    , &! INTENT(OUT  )
         topparu  (1:npoi)    , &! INTENT(OUT  )
         topparl  (1:npoi)    , &! INTENT(OUT  )
         asurd    (1:npoi,1:nband,latco) , &! INTENT(OUT  )
         asuri    (1:npoi,1:nband,latco) , &! INTENT(OUT  )
         coszen   (1:npoi))      ! INTENT(IN   )  
    !
    ! solar calculations for each waveband
    !
    xband=nband
    DO  ib = 1, nband
       !
       ! solsur sets surface albedos for soil and snow
       ! solalb performs the albedo calculations
       ! solarf uses the unit-incident-flux results from solalb
       ! to obtain absorbed fluxes sol[u,s,l,g,i] and 
       ! incident pars sunp[u,l]
       !
       CALL solsur (ib                               , &! INTENT(IN   )
            tmelt                            , &! INTENT(IN   )
            nsol                             , &! INTENT(IN   )
            albsod (1:npoi)                  , &! INTENT(OUt  )
            albsoi (1:npoi)                  , &! INTENT(OUt  )
            albsnd (1:npoi)                  , &! INTENT(OUt  )
            albsni (1:npoi)                  , &! INTENT(OUt  )
            indsol (1:npoi)                  , &! INTENT(IN   )
            wsoi   (1:npoi,1:nsoilay,latco)     , &! INTENT(IN   )
            wisoi  (1:npoi,1:nsoilay,latco)     , &! INTENT(IN   )
            albsav (1:npoi,latco)               , &! INTENT(IN   )
            albsan (1:npoi,latco)               , &! INTENT(IN   )
            tsno   (1:npoi,1:nsnolay,latco)     , &! INTENT(IN   )
            coszen (1:npoi)                  , &! INTENT(IN   )
            npoi                             , &! INTENT(IN   )
            nsoilay                          , &! INTENT(IN   )
            nsnolay                            )! INTENT(IN   )

       CALL solalb (ib               , &! INTENT(IN   )
            nVegClass        , &! INTENT(IN   )
            vegtype0 (1:npoi,latco), &! INTENT(IN   )
            avmuir_factor(1:nVegClass,1:2), &! INTENT(IN   )
            relod (1:npoi)   , &! INTENT(OUT  )
            reloi (1:npoi)   , &! INTENT(OUT  )
            indsol(1:npoi)   , &! INTENT(IN   )
            reupd (1:npoi)   , &! INTENT(OUT  )
            reupi (1:npoi)   , &! INTENT(OUT  )
            albsnd(1:npoi)   , &! INTENT(IN   )
            albsni(1:npoi)   , &! INTENT(IN   )
            albsod(1:npoi)   , &! INTENT(IN   )
            albsoi(1:npoi)   , &! INTENT(IN   )
            fl    (1:npoi,latco), &! INTENT(IN   )
            fu    (1:npoi,latco), &! INTENT(IN   )
            fi    (1:npoi,latco), &! INTENT(IN   )
            asurd (1:npoi,1:xband,latco) , &! INTENT(INOUT)! local
            asuri (1:npoi,1:xband,latco) , &! INTENT(INOUT)! local
            npoi             , &! INTENT(IN   )
            nband            , &! INTENT(IN   )
            nsol             , &! INTENT(IN   )
            ablod (1:npoi)   , &! INTENT(OUT  )
            abloi (1:npoi)   , &! INTENT(OUT  )
            flodd (1:npoi)   , &! INTENT(OUT  )
            dummy (1:npoi)   , &! INTENT(OUT  )
            flodi (1:npoi)   , &! INTENT(OUT  )
            floii (1:npoi)   , &! INTENT(OUT  )
            coszen(1:npoi)   , &! INTENT(IN   )
            terml (1:npoi,1:7)  , &! INTENT(OUT  )
            termu (1:npoi,1:7)  , &! INTENT(OUT  )
            lai   (1:npoi,1:2,latco), &! INTENT(IN   )
            sai   (1:npoi,1:2,latco), &! INTENT(IN   )
            abupd (1:npoi)      , &! INTENT(OUT  )
            abupi (1:npoi)      , &! INTENT(OUT  )
            fupdd (1:npoi)      , &! INTENT(OUT  )
            fupdi (1:npoi)      , &! INTENT(OUT  )
            fupii (1:npoi)      , &! INTENT(OUT  )
            fwetl (1:npoi)      , &! INTENT(IN   )
            rliql (1:npoi)      , &! INTENT(IN   )
            rliqu (1:npoi)      , &! INTENT(IN   )
            rliqs (1:npoi)      , &! INTENT(IN   )
            fwetu (1:npoi)      , &! INTENT(IN   )
            fwets (1:npoi)      , &! INTENT(IN   )
            rhoveg(1:nband,1:2) , &! INTENT(IN   )
            tauveg(1:nband,1:2) , &! INTENT(IN   )
            orieh (1:2)         , &! INTENT(IN   )
            oriev (1:2)         , &! INTENT(IN   )
            tl    (1:npoi,latco)   , &! INTENT(IN   )
            ts    (1:npoi,latco)   , &! INTENT(IN   )
            tu    (1:npoi,latco)   , &! INTENT(IN   )
            pi                  , &! INTENT(IN   )
            tmelt               , &! INTENT(IN   )
            epsilon                )! INTENT(IN   )


    END DO
    END IF
!     --------------------------- INPUT ---------------------------------------
!   
!    specify the parameters for albedo here:
    IceOceanAlb=0.0_r8
    chl = 0.10_r8              !chlorophyll concentration in mg/m3
    ! Two spectral surface albedos for direct (dir) and diffuse (dif)
    ! incident radiation are calculated. The spectral intervals are:
    !   s (shortwave)  = 0.2-0.7 micro-meters
    !   l (longwave)   = 0.7-5.0 micro-meters
    !
    tau=0.0_r8
    DO k=1,kMax
       DO i=1,ncols
          tau(i)=tau(i)+ taud(i,k) ! tau = SUM(taud(i,1:kMax))              !aerosol/cloud optical depth
       END DO
    END DO   
    DO i=1,ncols
        tau(i)=MIN(MAX(tau(i),0.0_r8),25.0_r8)      !aerosol/cloud optical depth
    END DO
    IF(TRIM(SLABOCEAN) == 'SLAB')THEN
       ncount=0
       DO i=1,ncols
          IF(imask(i) >= 1_i8) THEN
             ncount=ncount+1
             avisb(i)=asuri(ncount,1,latco)                   !asurd  (npoi,nband)   ! local  ! direct albedo of surface system
             avisd(i)=asurd(ncount,1,latco)                   !asuri  (npoi,nband)   ! local  ! diffuse albedo of surface system 
             anirb(i)=asuri(ncount,2,latco)
             anird(i)=asurd(ncount,2,latco)
          ELSE 
          !ELSE IF(ABS(tsea(i)).GE.271.16e0_r8 +0.01e0_r8) THEN
             IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) > tice+0.01_r8) THEN
                f=MAX(zenith(i),0.0e0_r8 )
                avisb(i)=GetOceanAlb(i,tau(i),f,wind(i),chl(i),0.2_r8,0.7_r8) !   s (shortwave)  = 0.2-0.7 micro-meters
                avisd(i)=oceald
                anirb(i)=GetOceanAlb(i,tau(i),f,wind(i),chl(i),0.7_r8,5.0_r8) !   l (longwave)   = 0.7-5.0 micro-meters
                anird(i)=oceald
             ELSE IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) <= tice+0.01_r8) THEN
                IF (TRIM(ICEMODEL)=='SSIB')THEN
                   f=MAX(zenith(i),0.0e0_r8 )
                   IceOceanAlb(i,1:2,1:2)=GetIceOceanAlb(i,latco,month(i),xVisBeam(i),xVisDiff(i),&
                                       xNirBeam(i),xNirDiff(i),f,LwSfcDown(i))
                   avisb(i)=IceOceanAlb(i,1,1)
                   anirb(i)=IceOceanAlb(i,2,1)
                   avisd(i)=IceOceanAlb(i,1,2)
                   anird(i)=IceOceanAlb(i,2,2)
                ELSE IF (TRIM(ICEMODEL)=='COLA')THEN
                   avisb(i)=icealv!icealv =         0.8e0_r8! constant icealv
                   avisd(i)=icealv!icealv =         0.8e0_r8! constant icealv
                   anirb(i)=icealn!icealn =         0.4e0_r8! constant icealn
                   anird(i)=icealn!icealn =         0.4e0_r8! constant icealn
                ELSE
                   STOP "ICEMODEL ->OPTIONS"
                END IF
             END IF
          END IF
       END DO
    ELSE IF(TRIM(SLABOCEAN) == 'COLA')THEN
       ncount=0
       DO i=1,ncols
          IF(imask(i) >= 1_i8) THEN
             ncount=ncount+1
             avisb(i)=asuri(ncount,1,latco)                   !asurd  (npoi,nband)   ! local  ! direct albedo of surface system
             avisd(i)=asurd(ncount,1,latco)                   !asuri  (npoi,nband)   ! local  ! diffuse albedo of surface system 
             anirb(i)=asuri(ncount,2,latco)
             anird(i)=asurd(ncount,2,latco)
          ELSE
             IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) > tice+0.01_r8) THEN
        !     IF(ABS(tsea(i)).GE.271.16e0_r8 +0.01e0_r8) THEN
                f=MAX(zenith(i),0.0e0_r8 )
                ocealb=0.12347e0_r8 +f*(0.34667e0_r8+f*(-1.7485e0_r8 + &
                    f*(2.04630e0_r8 -0.74839e0_r8 *f)))
                avisb(i)=ocealb
                avisd(i)=oceald!oceald =      0.0419e0_r8! constant oceald
                anirb(i)=ocealb
                anird(i)=oceald!oceald =      0.0419e0_r8! constant oceald
             ELSE IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) <= tice+0.01_r8) THEN
                IF (TRIM(ICEMODEL)=='SSIB')THEN
                   f=MAX(zenith(i),0.0e0_r8 )
                   IceOceanAlb(i,1:2,1:2)=GetIceOceanAlb(i,latco,month(i),xVisBeam(i),xVisDiff(i),&
                                          xNirBeam(i),xNirDiff(i),f,LwSfcDown(i))
                   avisb(i)=IceOceanAlb(i,1,1)
                   anirb(i)=IceOceanAlb(i,2,1)
                   avisd(i)=IceOceanAlb(i,1,2)
                   anird(i)=IceOceanAlb(i,2,2)
                ELSE IF (TRIM(ICEMODEL)=='COLA')THEN
                   avisb(i)=icealv!icealv =         0.8e0_r8! constant icealv
                   avisd(i)=icealv!icealv =         0.8e0_r8! constant icealv
                   anirb(i)=icealn!icealn =         0.4e0_r8! constant icealn
                   anird(i)=icealn!icealn =         0.4e0_r8! constant icealn
                ELSE
                   STOP "ICEMODEL ->OPTIONS"
                END IF
             END IF
          END IF
       END DO
    ELSE
       WRITE(0,*)"ERRO SLABOCEAN",TRIM(SLABOCEAN)
       STOP 
    END IF

  END SUBROUTINE Albedo_IBIS
  !
  !***************************************************************************
  !                      (imonth,iday,iyear)
  REAL(KIND=r8) FUNCTION julday (imonth,iday,iyear,tod)
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: imonth
    INTEGER, INTENT(IN   ) :: iday
    INTEGER, INTENT(IN   ) :: iyear
    REAL(KIND=r8)   , INTENT(IN   ) :: tod
    !
    ! compute the julian day from a normal date
    !
    julday= iday  &
         + MIN(1,MAX(0,imonth-1))*31  &
         + MIN(1,MAX(0,imonth-2))*(28+(1-MIN(1,MOD(iyear,4))))  &
         + MIN(1,MAX(0,imonth-3))*31  &
         + MIN(1,MAX(0,imonth-4))*30  &
         + MIN(1,MAX(0,imonth-5))*31  &
         + MIN(1,MAX(0,imonth-6))*30  &
         + MIN(1,MAX(0,imonth-7))*31  &
         + MIN(1,MAX(0,imonth-8))*31  &
         + MIN(1,MAX(0,imonth-9))*30  &
         + MIN(1,MAX(0,imonth-10))*31  &
         + MIN(1,MAX(0,imonth-11))*30  &
         + MIN(1,MAX(0,imonth-12))*31  &
         + tod/86400.0

  END FUNCTION julday

  SUBROUTINE sextrp &
       (Var_StepP    , &
        Var_Step0    , &
        Var_StepM    , &
        nmax         , &
        zmax         , &
        istrt        , &
        epsflt       , &
        intg           )
!    INTEGER         , PARAMETER     :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers

    INTEGER         , INTENT(in   ) :: nmax
    INTEGER         , INTENT(in   ) :: zmax
    INTEGER         , INTENT(in   ) :: istrt
    REAL(KIND=r8)   , INTENT(in   ) :: epsflt
    INTEGER         , INTENT(in   ) :: intg
    REAL(KIND=r8),    INTENT(in   ) :: Var_StepP(nmax,zmax)
    REAL(KIND=r8),    INTENT(inout) :: Var_Step0(nmax,zmax)
    REAL(KIND=r8),    INTENT(inout) :: Var_StepM(nmax,zmax)

    INTEGER :: i, k

    IF (intg == 2) THEN
       IF (istrt >= 1) THEN
          DO k=1,zmax
             DO i = 1, nmax
                Var_Step0    (i,k)=Var_StepP    (i,k)
             END DO
          END DO
       ELSE
          DO k=1,zmax
             DO i = 1, nmax
                IF(Var_Step0    (i,k) > 0.0_r8 ) THEN
                   Var_Step0(i,k)=Var_Step0(i,k)+epsflt*(Var_StepP(i,k)+Var_StepM(i,k)-2.0_r8  *Var_Step0(i,k))
                END IF
             END DO
          END DO
          DO k=1,zmax
             DO i = 1, nmax
                Var_StepM    (i,k)=Var_Step0    (i,k)
             END DO
          END DO
          DO k=1,zmax
             DO i = 1, nmax
                Var_Step0    (i,k)=Var_StepP     (i,k)
             END DO
          END DO
       END IF
    ELSE
       DO k=1,zmax
          DO i = 1, nmax
             Var_StepM    (i,k)=Var_StepP    (i,k)
             Var_Step0    (i,k)=Var_StepP    (i,k)
          END DO
       END DO
    END IF
  END SUBROUTINE sextrp

END MODULE Sfc_Ibis_Interface
