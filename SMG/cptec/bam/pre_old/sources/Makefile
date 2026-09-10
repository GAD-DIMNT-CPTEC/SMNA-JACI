#
#  $Author: pkubota $
#  $Date: 2007/03/06 14:43:32 $
#  $Revision: 1.3 $
#
# $Id: Makefile,v 1.3 2007/03/06 14:43:32 pkubota Exp $
#
SHELL=/bin/sh
#
# Makefile for Pos-Proccessing output global model
#
# Targets are the suffixes of the system-specific makefiles in
# the makefiles/ directory.
# For example, to build PosGrib for Solaris, give the command
#
#   make solaris
#
# This builds an intermediate library in the util/ directory,
# then builds the nedit and nc executables in the source/ directory.
#
aptions1="alb            albc        chp         nmi        dst        dstc       lndmk                                    "
aptions2="rolen          rolenc      snowc       soimo      soimoc     sstclm     ocmclm         co2clm                    "
aptions3="sstwkl         sstwklncep  snowwklncep snowweekly topospec   topogt30   toponavy       vartopo                   "
aptions4="albsib         vegmsk      parsib      ndviclm    ndvi       ibismsk    ibisveg        sstdirseason              "
aptions5="sandmskclm     sandmsk     claymskclm  claymsk    clmtclm    clmt       deltatclm      deltat                    "
aptions6="sib2clm        sib2        soiltexclm  soiltex    soilmsclm  soilms     sib2claymskclm sib2claymsk               "
aptions7="sib2sandmskclm sib2sandmsk sstdirec    co2direc   topograd   chpp       sstdailydirec  soilmwklcptec soilmweekly "

opt=" "
MAKE=make
PATH2=/scratchin/grupos/bam/home/paulo.kubota/bam_hybrid-2.2/pre/exec
PATH2=../../exec


all:
	@echo "Please specify target:"
	@echo "(For example, type \"make linux\" for a Linux system.)"
	@(cd  makefiles && ls -C Makefile* | sed -e 's/Makefile.//g')
	@echo "(For example, type \"make linux opt=alb\" for albedo option)"
	@echo "alb            albc            chp          nmi        dst        dstc           lndmk                                    "
	@echo "rolen          rolenc          snowc        soimo      soimoc     sstclm         ocmclm         co2clm                    "
	@echo "sstwkl         sstwklncep      snowwklncep  snowweekly topospec   topogt30       toponavy       vartopo                   "
	@echo "albsib         vegmsk          parsib       ndviclm    ndvi       ibismsk        ibisveg        sstdirseason              "
	@echo "sandmskclm     sandmsk         claymskclm   claymsk    clmtclm    clmt           deltatclm      deltat                    "
	@echo "sib2clm        sib2            soiltexclm   soiltex    soilmsclm  soilms         sib2claymskclm sib2claymsk               "
	@echo "sib2sandmskclm sib2sandmsk     sstdirec     co2direc   topograd   chpp           sstdailydirec  soilmwklcptec soilmweekly "

.DEFAULT:

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f bacio/v2.0.2/src_32bits/Makefile.$@ ];\
	   then ln -fs ../../../makefiles/Makefile.$@ bacio/v2.0.2/src_32bits/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f bacio/v2.0.2/src_64bits/Makefile.$@ ];\
	   then ln -fs ../../../makefiles/Makefile.$@ bacio/v2.0.2/src_64bits/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f w3lib-1.4/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ w3lib-1.4/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f nemsio/v2.2.2/src/Makefile.$@ ];\
	   then ln -fs ../../../makefiles/Makefile.$@ nemsio/v2.2.2/src/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f Albedo/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ Albedo/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f AlbedoClima/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ AlbedoClima/Makefile.$@; fi)

#	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f Chopping_serial/Makefile.$@ ];\
#	   then ln -fs ../makefiles/Makefile.$@ Chopping_serial/Makefile.$@; fi)

#	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f Chopping_parallel/DummyMPI/Makefile.$@ ];\
#	   then ln -fs ../../makefiles/Makefile.$@ Chopping_parallel/DummyMPI/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f Chopping_parallel/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ Chopping_parallel/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f NormalModes/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ NormalModes/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f DeepSoilTemperature/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ DeepSoilTemperature/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f DeepSoilTemperatureClima/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ DeepSoilTemperatureClima/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f LandSeaMask/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ LandSeaMask/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f RoughnessLength/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ RoughnessLength/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f NDVIClima/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ NDVIClima/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f NDVI/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ NDVI/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f RoughnessLengthClima/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ RoughnessLengthClima/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f SnowClima/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ SnowClima/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f SoilMoisture/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ SoilMoisture/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f SoilMoistureClima/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ SoilMoistureClima/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f SSTClima/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ SSTClima/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f OCMClima/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ OCMClima/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f FLUXCO2Clima/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ FLUXCO2Clima/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f SSTWeekly/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ SSTWeekly/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f SNOWWeekly/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ SNOWWeekly/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f SoilMoistureWeekly/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ SoilMoistureWeekly/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f SSTWeeklyNCEP/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ SSTWeeklyNCEP/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f SNOWWeeklyNCEP/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ SNOWWeeklyNCEP/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f SoilMoistureWeeklyCPTEC/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ SoilMoistureWeeklyCPTEC/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f TopoSpectral/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ TopoSpectral/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f TopoWaterPercGT30/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ TopoWaterPercGT30/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f TopoWaterPercNavy/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ TopoWaterPercNavy/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f VarTopo/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ VarTopo/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f VegetationAlbedoSSiB/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ VegetationAlbedoSSiB/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f VegetationMask/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ VegetationMask/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f VegetationMaskIBIS/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ VegetationMaskIBIS/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f VegetationMaskSSiB/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ VegetationMaskSSiB/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f VegetationMaskIBISClima/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ VegetationMaskIBISClima/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f PorceSandMaskIBISClima/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ PorceSandMaskIBISClima/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f PorceSandMaskIBIS/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ PorceSandMaskIBIS/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f PorceClayMaskIBISClima/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ PorceClayMaskIBISClima/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f PorceClayMaskIBIS/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ PorceClayMaskIBIS/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f TemperatureClima/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ TemperatureClima/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f Temperature/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ Temperature/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f DeltaTempColdestClima/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ DeltaTempColdestClima/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f DeltaTempColdest/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ DeltaTempColdest/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f VegetationMaskSiB2Clima/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ VegetationMaskSiB2Clima/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f VegetationMaskSiB2/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ VegetationMaskSiB2/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f SoilTextureMaskSiB2Clima/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ SoilTextureMaskSiB2Clima/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f SoilTextureMaskSiB2/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ SoilTextureMaskSiB2/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f CLimaSoilMoistureClima/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ CLimaSoilMoistureClima/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f CLimaSoilMoisture/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ CLimaSoilMoisture/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f PorceClayMaskSiB2Clima/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ PorceClayMaskSiB2Clima/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f PorceClayMaskSiB2/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ PorceClayMaskSiB2/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f PorceSandMaskSiB2Clima/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ PorceSandMaskSiB2Clima/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f PorceSandMaskSiB2/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ PorceSandMaskSiB2/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f SSTMonthlyDirec/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ SSTMonthlyDirec/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f SSTSeasonDirec/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ SSTSeasonDirec/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f CO2MonthlyDirec/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ CO2MonthlyDirec/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f SSTDailyDirec/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ SSTDailyDirec/Makefile.$@; fi)

	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f TopographyGradient/Makefile.$@ ];\
	   then ln -fs ../makefiles/Makefile.$@ TopographyGradient/Makefile.$@; fi)

#
# builds the binaries
#

	@- (cd bacio/v2.0.2/src_32bits;$(MAKE) -f Makefile.$@;$(MAKE) -f Makefile.$@)
	@- (cd bacio/v2.0.2/src_64bits;$(MAKE) -f Makefile.$@;$(MAKE) -f Makefile.$@)
	@- (cd w3lib-1.4;$(MAKE) -f Makefile.$@;$(MAKE) -f Makefile.$@)
	@- (cd nemsio/v2.2.2/src;$(MAKE) -f Makefile.$@;$(MAKE) -f Makefile.$@)

	@ (if [ ${opt} = "alb" ];then cd Albedo;$(MAKE) -f Makefile.$@;mv Albedo $(PATH2)/Albedo; fi)
	@ (if [ ${opt} = " " ];then   cd Albedo;$(MAKE) -f Makefile.$@;mv Albedo $(PATH2)/Albedo; fi)

	@ (if [ ${opt} = "albc" ];then cd AlbedoClima;$(MAKE) -f Makefile.$@;mv AlbedoClima $(PATH2)/AlbedoClima; fi)
	@ (if [ ${opt} = " " ];then   cd AlbedoClima;$(MAKE) -f Makefile.$@;mv AlbedoClima $(PATH2)/AlbedoClima; fi)

#	@ (if [ ${opt} = "chp" ];then cd Chopping_serial;$(MAKE) -f Makefile.$@;mv Chopping_serial $(PATH2)/Chopping_serial; fi)
#	@ (if [ ${opt} = " " ];then   cd Chopping_serial;$(MAKE) -f Makefile.$@;mv Chopping_serial $(PATH2)/Chopping_serial; fi)

#	@ (if [ ${opt} = "mpid" ];then cd Chopping_parallel/DummyMPI;$(MAKE) -f Makefile.$@ ;cd ../;cp ummyMPI/include/mpif.h .; fi)
#	@ (if [ ${opt} = " " ];then    cd Chopping_parallel/DummyMPI;$(MAKE) -f Makefile.$@ ;cd ../;cp DummyMPI/include/mpif.h .; fi)

	@ (if [ ${opt} = "chpp" ];then cd Chopping_parallel;$(MAKE) -f Makefile.$@;mv Chopping_parallel $(PATH2)/Chopping_parallel; fi)
	@ (if [ ${opt} = " " ];then   cd Chopping_parallel;$(MAKE) -f Makefile.$@;mv Chopping_parallel $(PATH2)/Chopping_parallel; fi)

	@ (if [ ${opt} = "nmi" ];then cd NormalModes;$(MAKE) -f Makefile.$@;mv NormalModes $(PATH2)/NormalModes; fi)
	@ (if [ ${opt} = " " ];then   cd NormalModes;$(MAKE) -f Makefile.$@;mv NormalModes $(PATH2)/NormalModes; fi)

	@ (if [ ${opt} = "dst" ];then cd DeepSoilTemperature;$(MAKE) -f Makefile.$@;mv DeepSoilTemperature $(PATH2)/DeepSoilTemperature; fi)
	@ (if [ ${opt} = " " ];then   cd DeepSoilTemperature;$(MAKE) -f Makefile.$@;mv DeepSoilTemperature $(PATH2)/DeepSoilTemperature; fi)

	@ (if [ ${opt} = "dstc" ];then cd DeepSoilTemperatureClima;$(MAKE) -f Makefile.$@;mv DeepSoilTemperatureClima $(PATH2)/DeepSoilTemperatureClima; fi)
	@ (if [ ${opt} = " " ];then   cd DeepSoilTemperatureClima;$(MAKE) -f Makefile.$@;mv DeepSoilTemperatureClima $(PATH2)/DeepSoilTemperatureClima; fi)

	@ (if [ ${opt} = "lndmk" ];then cd LandSeaMask;$(MAKE) -f Makefile.$@;mv LandSeaMask $(PATH2)/LandSeaMask; fi)
	@ (if [ ${opt} = " " ];then   cd LandSeaMask;$(MAKE) -f Makefile.$@;mv LandSeaMask $(PATH2)/LandSeaMask; fi)

	@ (if [ ${opt} = "rolen" ];then cd RoughnessLength;$(MAKE) -f Makefile.$@;mv RoughnessLength $(PATH2)/RoughnessLength; fi)
	@ (if [ ${opt} = " " ];then   cd RoughnessLength;$(MAKE) -f Makefile.$@;mv RoughnessLength $(PATH2)/RoughnessLength; fi)

	@ (if [ ${opt} = "rolenc" ];then cd RoughnessLengthClima;$(MAKE) -f Makefile.$@;mv RoughnessLengthClima $(PATH2)/RoughnessLengthClima; fi)
	@ (if [ ${opt} = " " ];then   cd RoughnessLengthClima;$(MAKE) -f Makefile.$@;mv RoughnessLengthClima $(PATH2)/RoughnessLengthClima; fi)

	@ (if [ ${opt} = "ndviclm" ];then cd NDVIClima;$(MAKE) -f Makefile.$@;mv NDVIClima $(PATH2)/NDVIClima; fi)
	@ (if [ ${opt} = " " ];then   cd NDVIClima;$(MAKE) -f Makefile.$@;mv NDVIClima $(PATH2)/NDVIClima; fi)

	@ (if [ ${opt} = "ndvi" ];then cd NDVI;$(MAKE) -f Makefile.$@;mv NDVI $(PATH2)/NDVI; fi)
	@ (if [ ${opt} = " " ];then   cd NDVI;$(MAKE) -f Makefile.$@;mv NDVI $(PATH2)/NDVI; fi)

	@ (if [ ${opt} = "snowc" ];then cd SnowClima;$(MAKE) -f Makefile.$@;mv SnowClima $(PATH2)/SnowClima; fi)
	@ (if [ ${opt} = " " ];then   cd SnowClima;$(MAKE) -f Makefile.$@;mv SnowClima $(PATH2)/SnowClima; fi)

	@ (if [ ${opt} = "soimo" ];then cd SoilMoisture;$(MAKE) -f Makefile.$@;mv SoilMoisture $(PATH2)/SoilMoisture; fi)
	@ (if [ ${opt} = " " ];then   cd SoilMoisture;$(MAKE) -f Makefile.$@;mv SoilMoisture $(PATH2)/SoilMoisture; fi)

	@ (if [ ${opt} = "soimoc" ];then cd SoilMoistureClima;$(MAKE) -f Makefile.$@;mv SoilMoistureClima $(PATH2)/SoilMoistureClima; fi)
	@ (if [ ${opt} = " " ];then   cd SoilMoistureClima;$(MAKE) -f Makefile.$@;mv SoilMoistureClima $(PATH2)/SoilMoistureClima; fi)

	@ (if [ ${opt} = "sstclm" ];then cd SSTClima;$(MAKE) -f Makefile.$@;mv SSTClima $(PATH2)/SSTClima; fi)
	@ (if [ ${opt} = " " ];then   cd SSTClima;$(MAKE) -f Makefile.$@;mv SSTClima $(PATH2)/SSTClima; fi)

	@ (if [ ${opt} = "ocmclm" ];then cd OCMClima;$(MAKE) -f Makefile.$@;mv OCMClima $(PATH2)/OCMClima; fi)
	@ (if [ ${opt} = " " ];then   cd OCMClima;$(MAKE) -f Makefile.$@;mv OCMClima $(PATH2)/OCMClima; fi)

	@ (if [ ${opt} = "co2clm" ];then cd FLUXCO2Clima;$(MAKE) -f Makefile.$@;mv FLUXCO2Clima $(PATH2)/FLUXCO2Clima; fi)
	@ (if [ ${opt} = " " ];then   cd FLUXCO2Clima;$(MAKE) -f Makefile.$@;mv FLUXCO2Clima $(PATH2)/FLUXCO2Clima; fi)

	@ (if [ ${opt} = "sstwkl" ];then cd SSTWeekly;$(MAKE) -f Makefile.$@;mv SSTWeekly $(PATH2)/SSTWeekly; fi)
	@ (if [ ${opt} = " " ];then   cd SSTWeekly;$(MAKE) -f Makefile.$@;mv SSTWeekly $(PATH2)/SSTWeekly; fi)

	@ (if [ ${opt} = "snowweekly" ];then cd SNOWWeekly ;$(MAKE) -f Makefile.$@;mv SNOWWeekly $(PATH2)/SNOWWeekly; fi)
	@ (if [ ${opt} = " " ];then   cd SNOWWeekly;$(MAKE) -f Makefile.$@;mv SNOWWeekly $(PATH2)/SNOWWeekly; fi)

	@ (if [ ${opt} = "soilmweekly" ];then cd SoilMoistureWeekly ;$(MAKE) -f Makefile.$@;mv SoilMoistureWeekly $(PATH2)/SoilMoistureWeekly; fi)
	@ (if [ ${opt} = " " ];then   cd SoilMoistureWeekly;$(MAKE) -f Makefile.$@;mv SoilMoistureWeekly $(PATH2)/SoilMoistureWeekly; fi)

	@ (if [ ${opt} = "sstwklncep" ];then cd SSTWeeklyNCEP;$(MAKE) -f Makefile.$@;mv SSTWeeklyNCEP $(PATH2)/SSTWeeklyNCEP; fi)
	@ (if [ ${opt} = " " ];then   cd SSTWeeklyNCEP;$(MAKE) -f Makefile.$@;mv SSTWeeklyNCEP $(PATH2)/SSTWeeklyNCEP; fi)

	@ (if [ ${opt} = "snowwklncep" ];then cd SNOWWeeklyNCEP;$(MAKE) -f Makefile.$@;mv SNOWWeeklyNCEP $(PATH2)/SNOWWeeklyNCEP; fi)
	@ (if [ ${opt} = " " ];then   cd SNOWWeeklyNCEP;$(MAKE) -f Makefile.$@;mv SNOWWeeklyNCEP $(PATH2)/SNOWWeeklyNCEP; fi)

	@ (if [ ${opt} = "soilmwklcptec" ];then cd SoilMoistureWeeklyCPTEC;$(MAKE) -f Makefile.$@;mv SoilMoistureWeeklyCPTEC $(PATH2)/SoilMoistureWeeklyCPTEC; fi)
	@ (if [ ${opt} = " " ];then   cd SoilMoistureWeeklyCPTEC;$(MAKE) -f Makefile.$@;mv SoilMoistureWeeklyCPTEC $(PATH2)/SoilMoistureWeeklyCPTEC; fi)

	@ (if [ ${opt} = "topospec" ];then cd TopoSpectral;$(MAKE) -f Makefile.$@;mv TopoSpectral $(PATH2)/TopoSpectral; fi)
	@ (if [ ${opt} = " " ];then   cd TopoSpectral;$(MAKE) -f Makefile.$@;mv TopoSpectral $(PATH2)/TopoSpectral; fi)

	@ (if [ ${opt} = "topogt30" ];then cd TopoWaterPercGT30;$(MAKE) -f Makefile.$@;mv TopoWaterPercGT30 $(PATH2)/TopoWaterPercGT30; fi)
	@ (if [ ${opt} = " " ];then   cd TopoWaterPercGT30;$(MAKE) -f Makefile.$@;mv TopoWaterPercGT30 $(PATH2)/TopoWaterPercGT30; fi)

	@ (if [ ${opt} = "toponavy" ];then cd TopoWaterPercNavy;$(MAKE) -f Makefile.$@;mv TopoWaterPercNavy $(PATH2)/TopoWaterPercNavy; fi)
	@ (if [ ${opt} = " " ];then   cd TopoWaterPercNavy;$(MAKE) -f Makefile.$@;mv TopoWaterPercNavy $(PATH2)/TopoWaterPercNavy; fi)

	@ (if [ ${opt} = "vartopo" ];then cd VarTopo;$(MAKE) -f Makefile.$@;mv VarTopo $(PATH2)/VarTopo; fi)
	@ (if [ ${opt} = " " ];then   cd VarTopo;$(MAKE) -f Makefile.$@;mv VarTopo $(PATH2)/VarTopo; fi)

	@ (if [ ${opt} = "albsib" ];then cd VegetationAlbedoSSiB;$(MAKE) -f Makefile.$@;mv VegetationAlbedoSSiB $(PATH2)/VegetationAlbedoSSiB; fi)
	@ (if [ ${opt} = " " ];then   cd VegetationAlbedoSSiB;$(MAKE) -f Makefile.$@;mv VegetationAlbedoSSiB $(PATH2)/VegetationAlbedoSSiB; fi)

	@ (if [ ${opt} = "vegmsk" ];then cd VegetationMask;$(MAKE) -f Makefile.$@;mv VegetationMask $(PATH2)/VegetationMask; fi)
	@ (if [ ${opt} = " " ];then   cd VegetationMask;$(MAKE) -f Makefile.$@;mv VegetationMask $(PATH2)/VegetationMask; fi)

	@ (if [ ${opt} = "ibisveg" ];then cd VegetationMaskSSiB;$(MAKE) -f Makefile.$@;mv VegetationMaskSSiB $(PATH2)/VegetationMaskSSiB; fi)
	@ (if [ ${opt} = " " ];then   cd VegetationMaskSSiB;$(MAKE) -f Makefile.$@;mv VegetationMaskSSiB $(PATH2)/VegetationMaskSSiB; fi)

	@ (if [ ${opt} = "ibismsk" ];then cd VegetationMaskIBIS;$(MAKE) -f Makefile.$@;mv VegetationMaskIBIS $(PATH2)/VegetationMaskIBIS; fi)
	@ (if [ ${opt} = " " ];then   cd VegetationMaskIBIS;$(MAKE) -f Makefile.$@;mv VegetationMaskIBIS $(PATH2)/VegetationMaskIBIS; fi)

	@ (if [ ${opt} = "ibisveg" ];then cd VegetationMaskIBISClima;$(MAKE) -f Makefile.$@;mv VegetationMaskIBISClima $(PATH2)/VegetationMaskIBISClima; fi)
	@ (if [ ${opt} = " " ];then   cd VegetationMaskIBISClima;$(MAKE) -f Makefile.$@;mv VegetationMaskIBISClima $(PATH2)/VegetationMaskIBISClima; fi)

	@ (if [ ${opt} = "sandmskclm" ];then cd PorceSandMaskIBISClima;$(MAKE) -f Makefile.$@;mv PorceSandMaskIBISClima $(PATH2)/PorceSandMaskIBISClima; fi)
	@ (if [ ${opt} = " " ];then   cd PorceSandMaskIBISClima;$(MAKE) -f Makefile.$@;mv PorceSandMaskIBISClima $(PATH2)/PorceSandMaskIBISClima; fi)

	@ (if [ ${opt} = "sandmsk" ];then cd PorceSandMaskIBIS;$(MAKE) -f Makefile.$@;mv PorceSandMaskIBIS $(PATH2)/PorceSandMaskIBIS; fi)
	@ (if [ ${opt} = " " ];then   cd PorceSandMaskIBIS;$(MAKE) -f Makefile.$@;mv PorceSandMaskIBIS $(PATH2)/PorceSandMaskIBIS; fi)

	@ (if [ ${opt} = "claymskclm" ];then cd PorceClayMaskIBISClima;$(MAKE) -f Makefile.$@;mv PorceClayMaskIBISClima $(PATH2)/PorceClayMaskIBISClima; fi)
	@ (if [ ${opt} = " " ];then   cd PorceClayMaskIBISClima;$(MAKE) -f Makefile.$@;mv PorceClayMaskIBISClima $(PATH2)/PorceClayMaskIBISClima; fi)

	@ (if [ ${opt} = "claymsk" ];then cd PorceClayMaskIBIS;$(MAKE) -f Makefile.soilmsclm$@;mv PorceClayMaskIBIS $(PATH2)/PorceClayMaskIBIS; fi)
	@ (if [ ${opt} = " " ];then   cd PorceClayMaskIBIS;$(MAKE) -f Makefile.$@;mv PorceClayMaskIBIS $(PATH2)/PorceClayMaskIBIS; fi)

	@ (if [ ${opt} = "clmtclm" ];then cd TemperatureClima;$(MAKE) -f Makefile.$@;mv TemperatureClima $(PATH2)/TemperatureClima; fi)
	@ (if [ ${opt} = " " ];then   cd TemperatureClima;$(MAKE) -f Makefile.$@;mv TemperatureClima $(PATH2)/TemperatureClima; fi)

	@ (if [ ${opt} = "clmt" ];then cd Temperature;$(MAKE) -f Makefile.$@;mv Temperature $(PATH2)/Temperature; fi)
	@ (if [ ${opt} = " " ];then   cd Temperature;$(MAKE) -f Makefile.$@;mv Temperature $(PATH2)/Temperature; fi)

	@ (if [ ${opt} = "deltatclm" ];then cd DeltaTempColdestClima;$(MAKE) -f Makefile.$@;mv DeltaTempColdestClima $(PATH2)/DeltaTempColdestClima; fi)
	@ (if [ ${opt} = " " ];then   cd DeltaTempColdestClima;$(MAKE) -f Makefile.$@;mv DeltaTempColdestClima $(PATH2)/DeltaTempColdestClima; fi)

	@ (if [ ${opt} = "deltat" ];then cd DeltaTempColdest;$(MAKE) -f Makefile.$@;mv DeltaTempColdest $(PATH2)/DeltaTempColdest; fi)
	@ (if [ ${opt} = " " ];then   cd DeltaTempColdest;$(MAKE) -f Makefile.$@;mv DeltaTempColdest $(PATH2)/DeltaTempColdest; fi)

	@ (if [ ${opt} = "sib2clm" ];then cd VegetationMaskSiB2Clima;$(MAKE) -f Makefile.$@;mv VegetationMaskSiB2Clima $(PATH2)/VegetationMaskSiB2Clima; fi)
	@ (if [ ${opt} = " " ];then   cd VegetationMaskSiB2Clima;$(MAKE) -f Makefile.$@;mv VegetationMaskSiB2Clima $(PATH2)/VegetationMaskSiB2Clima; fi)

	@ (if [ ${opt} = "sib2" ];then cd VegetationMaskSiB2;$(MAKE) -f Makefile.$@;mv VegetationMaskSiB2 $(PATH2)/VegetationMaskSiB2; fi)
	@ (if [ ${opt} = " " ];then   cd VegetationMaskSiB2;$(MAKE) -f Makefile.$@;mv VegetationMaskSiB2 $(PATH2)/VegetationMaskSiB2; fi)

	@ (if [ ${opt} = "soiltexclm" ];then cd SoilTextureMaskSiB2Clima;$(MAKE) -f Makefile.$@;mv SoilTextureMaskSiB2Clima $(PATH2)/SoilTextureMaskSiB2Clima; fi)
	@ (if [ ${opt} = " " ];then   cd SoilTextureMaskSiB2Clima;$(MAKE) -f Makefile.$@;mv SoilTextureMaskSiB2Clima $(PATH2)/SoilTextureMaskSiB2Clima; fi)

	@ (if [ ${opt} = "soiltex" ];then cd SoilTextureMaskSiB2;$(MAKE) -f Makefile.$@;mv SoilTextureMaskSiB2 $(PATH2)/SoilTextureMaskSiB2; fi)
	@ (if [ ${opt} = " " ];then   cd SoilTextureMaskSiB2;$(MAKE) -f Makefile.$@;mv SoilTextureMaskSiB2 $(PATH2)/SoilTextureMaskSiB2; fi)

	@ (if [ ${opt} = "soilmsclm" ];then cd CLimaSoilMoistureClima;$(MAKE) -f Makefile.$@;mv CLimaSoilMoistureClima $(PATH2)/CLimaSoilMoistureClima; fi)
	@ (if [ ${opt} = " " ];then   cd CLimaSoilMoistureClima;$(MAKE) -f Makefile.$@;mv CLimaSoilMoistureClima $(PATH2)/CLimaSoilMoistureClima; fi)

	@ (if [ ${opt} = "soilms" ];then cd CLimaSoilMoisture;$(MAKE) -f Makefile.$@;mv CLimaSoilMoisture $(PATH2)/CLimaSoilMoisture; fi)
	@ (if [ ${opt} = " " ];then   cd CLimaSoilMoisture;$(MAKE) -f Makefile.$@;mv CLimaSoilMoisture $(PATH2)/CLimaSoilMoisture; fi)

	@ (if [ ${opt} = "sib2claymskclm" ];then cd PorceClayMaskSiB2Clima;$(MAKE) -f Makefile.$@;mv PorceClayMaskSiB2Clima $(PATH2)/PorceClayMaskSiB2Clima; fi)
	@ (if [ ${opt} = " " ];then   cd PorceClayMaskSiB2Clima;$(MAKE) -f Makefile.$@;mv PorceClayMaskSiB2Clima $(PATH2)/PorceClayMaskSiB2Clima; fi)

	@ (if [ ${opt} = "sib2claymsk" ];then cd PorceClayMaskSiB2;$(MAKE) -f Makefile.$@;mv PorceClayMaskSiB2 $(PATH2)/PorceClayMaskSiB2; fi)
	@ (if [ ${opt} = " " ];then   cd PorceClayMaskSiB2;$(MAKE) -f Makefile.$@;mv PorceClayMaskSiB2 $(PATH2)/PorceClayMaskSiB2; fi)

	@ (if [ ${opt} = "sib2sandmskclm" ];then cd PorceSandMaskSiB2Clima;$(MAKE) -f Makefile.$@;mv PorceSandMaskSiB2Clima $(PATH2)/PorceSandMaskSiB2Clima; fi)
	@ (if [ ${opt} = " " ];then   cd PorceSandMaskSiB2Clima;$(MAKE) -f Makefile.$@;mv PorceSandMaskSiB2Clima $(PATH2)/PorceSandMaskSiB2Clima; fi)

	@ (if [ ${opt} = "sib2sandmsk" ];then cd PorceSandMaskSiB2;$(MAKE) -f Makefile.$@;mv PorceSandMaskSiB2 $(PATH2)/PorceSandMaskSiB2; fi)
	@ (if [ ${opt} = " " ];then   cd PorceSandMaskSiB2;$(MAKE) -f Makefile.$@;mv PorceSandMaskSiB2 $(PATH2)/PorceSandMaskSiB2; fi)

	@ (if [ ${opt} = "sstdirec" ];then cd SSTMonthlyDirec;$(MAKE) -f Makefile.$@;mv SSTMonthlyDirec $(PATH2)/SSTMonthlyDirec; fi)
	@ (if [ ${opt} = " " ];then   cd SSTMonthlyDirec;$(MAKE) -f Makefile.$@;mv SSTMonthlyDirec $(PATH2)/SSTMonthlyDirec; fi)

	@ (if [ ${opt} = "sstdirseason" ];then cd SSTSeasonDirec;$(MAKE) -f Makefile.$@;mv SSTSeasonDirec $(PATH2)/SSTSeasonDirec; fi)
	@ (if [ ${opt} = " " ];then   cd SSTSeasonDirec;$(MAKE) -f Makefile.$@;mv SSTSeasonDirec $(PATH2)/SSTSeasonDirec; fi)

	@ (if [ ${opt} = "co2direc" ];then cd CO2MonthlyDirec;$(MAKE) -f Makefile.$@;mv CO2MonthlyDirec $(PATH2)/CO2MonthlyDirec; fi)
	@ (if [ ${opt} = " " ];then   cd CO2MonthlyDirec;$(MAKE) -f Makefile.$@;mv CO2MonthlyDirec $(PATH2)/CO2MonthlyDirec; fi)

	@ (if [ ${opt} = "sstdailydirec" ];then cd SSTDailyDirec;$(MAKE) -f Makefile.$@;mv SSTDailyDirec $(PATH2)/SSTDailyDirec; fi)
	@ (if [ ${opt} = " " ];then   cd SSTDailyDirec;$(MAKE) -f Makefile.$@;mv SSTDailyDirec $(PATH2)/SSTDailyDirec; fi)

	@ (if [ ${opt} = "topograd" ];then cd TopographyGradient;$(MAKE) -f Makefile.$@;mv TopographyGradient $(PATH2)/TopographyGradient; fi)
	@ (if [ ${opt} = " " ];then   cd TopographyGradient;$(MAKE) -f Makefile.$@;mv TopographyGradient $(PATH2)/TopographyGradient; fi)

# We need a "dev-all" target that builds the docs plus binaries, but
# that doesn't work since we require the user to specify the target.  More
# thought is needed

clean:
	@echo "$(MAKE)"

	(cd w3lib-1.4;$(MAKE) -f Makefile.common clean)
	(cd bacio/v2.0.2/src_32bits/;$(MAKE) -f Makefile.common clean)
	(cd bacio/v2.0.2/src_64bits/;$(MAKE) -f Makefile.common clean)
	(cd nemsio/v2.2.2/src/;$(MAKE) -f Makefile.common clean)
	(cd Albedo;$(MAKE) -f Makefile.common clean)
	(cd AlbedoClima;$(MAKE) -f Makefile.common clean)
#	(cd Chopping_serial;$(MAKE) -f Makefile.common clean)
	(cd Chopping_parallel;$(MAKE) -f Makefile.common clean)
	(cd NormalModes;$(MAKE) -f Makefile.common clean)
	(cd DeepSoilTemperature;$(MAKE) -f Makefile.common clean)
	(cd DeepSoilTemperatureClima;$(MAKE) -f Makefile.common clean)
	(cd LandSeaMask;$(MAKE) -f Makefile.common clean)
	(cd RoughnessLength;$(MAKE) -f Makefile.common clean)
	(cd RoughnessLengthClima;$(MAKE) -f Makefile.common clean)
	(cd NDVIClima;$(MAKE) -f Makefile.common clean)
	(cd NDVI;$(MAKE) -f Makefile.common clean)	
	(cd SnowClima;$(MAKE) -f Makefile.common clean)
	(cd SoilMoisture;$(MAKE) -f Makefile.common clean)
	(cd SoilMoistureClima;$(MAKE) -f Makefile.common clean)
	(cd SSTClima;$(MAKE) -f Makefile.common clean)
	(cd OCMClima;$(MAKE) -f Makefile.common clean)
	(cd FLUXCO2Clima;$(MAKE) -f Makefile.common clean)
	(cd SSTWeekly;$(MAKE) -f Makefile.common clean)
	(cd SNOWWeekly;$(MAKE) -f Makefile.common clean)
	(cd SoilMoistureWeekly;$(MAKE) -f Makefile.common clean)
	(cd SSTWeeklyNCEP;$(MAKE) -f Makefile.common clean)
	(cd SNOWWeeklyNCEP;$(MAKE) -f Makefile.common clean)
	(cd SoilMoistureWeeklyCPTEC;$(MAKE) -f Makefile.common clean)
	(cd TopoSpectral;$(MAKE) -f Makefile.common clean)
	(cd TopoWaterPercGT30;$(MAKE) -f Makefile.common clean)
	(cd TopoWaterPercNavy;$(MAKE) -f Makefile.common clean)
	(cd VarTopo;$(MAKE) -f Makefile.common clean)
	(cd VegetationAlbedoSSiB;$(MAKE) -f Makefile.common clean)
	(cd VegetationMask;$(MAKE) -f Makefile.common clean)
	(cd VegetationMaskIBIS;$(MAKE) -f Makefile.common clean)	
	(cd VegetationMaskSSiB;$(MAKE) -f Makefile.common clean)
	(cd VegetationMaskIBISClima;$(MAKE) -f Makefile.common clean)
	(cd PorceSandMaskIBISClima;$(MAKE) -f Makefile.common clean)
	(cd PorceSandMaskIBIS;$(MAKE) -f Makefile.common clean)
	(cd PorceClayMaskIBISClima;$(MAKE) -f Makefile.common clean)
	(cd PorceClayMaskIBIS;$(MAKE) -f Makefile.common clean)
	(cd TemperatureClima;$(MAKE) -f Makefile.common clean)
	(cd Temperature;$(MAKE) -f Makefile.common clean)
	(cd DeltaTempColdestClima;$(MAKE) -f Makefile.common clean)
	(cd DeltaTempColdest;$(MAKE) -f Makefile.common clean)
	(cd VegetationMaskSiB2Clima;$(MAKE) -f Makefile.common clean)
	(cd VegetationMaskSiB2;$(MAKE) -f Makefile.common clean)
	(cd SoilTextureMaskSiB2Clima;$(MAKE) -f Makefile.common clean)
	(cd SoilTextureMaskSiB2;$(MAKE) -f Makefile.common clean)
	(cd CLimaSoilMoistureClima;$(MAKE) -f Makefile.common clean)
	(cd CLimaSoilMoisture;$(MAKE) -f Makefile.common clean)
	(cd PorceClayMaskSiB2Clima;$(MAKE) -f Makefile.common clean)
	(cd PorceClayMaskSiB2;$(MAKE) -f Makefile.common clean)
	(cd PorceSandMaskSiB2Clima;$(MAKE) -f Makefile.common clean)
	(cd PorceSandMaskSiB2;$(MAKE) -f Makefile.common clean)
	(cd SSTMonthlyDirec;$(MAKE) -f Makefile.common clean)
	(cd SSTSeasonDirec;$(MAKE) -f Makefile.common clean)
	(cd CO2MonthlyDirec;$(MAKE) -f Makefile.common clean)
	(cd SSTDailyDirec;$(MAKE) -f Makefile.common clean)
	(cd TopographyGradient;$(MAKE) -f Makefile.common clean)
