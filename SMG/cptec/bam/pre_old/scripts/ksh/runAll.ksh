#! /bin/ksh
#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#

set -x
export dirhome=/gfs/dk12/pkubota/mcgaibis-4.0.0/pre
export dirdata=/gfs/dk12/pkubota/mcgaibis-4.0.0
export dirgrads=/usr/local/grads

# Machine options: SX6; Linux
export MAQUI=SX6

# Set  Output Res for Chopping
export RESOUT=62
export KMOUT=28

# Set  T170 Quadratic
export RESO=62
export IM=192
export JM=96

# Set  T170 Quadratic

export SetLinear=FALSE

if [ "$SetLinear" = "TRUE" ]; then
if [ ${RESOUT} -lt 10000 ]; then
export TRUNC=TL${RESOUT}
if [ ${RESOUT} -lt 1000 ]; then
export TRUNC=TL0${RESOUT}
if [ ${RESOUT} -lt 100 ]; then
export TRUNC=TL00${RESOUT}
fi
fi
fi
else
if [ ${RESOUT} -lt 10000 ]; then
export TRUNC=TQ${RESOUT}
if [ ${RESOUT} -lt 1000 ]; then
export TRUNC=TQ0${RESOUT}
if [ ${RESOUT} -lt 100 ]; then
export TRUNC=TQ00${RESOUT}
fi
fi
fi
fi

if [ ${KMOUT} -lt 1000 ]; then
export prefixVert=${KMOUT}
if [ ${KMOUT} -lt 100 ]; then
export prefixVert=0${KMOUT}
if [ ${KMOUT} -lt 10 ]; then
export prefixVert=00${KMOUT}
fi
fi
fi

if [ ${JM} -lt "10000" ]; then
export prefixGrid=0${JM}
if [ ${JM} -lt 1000 ]; then
export prefixGrid=00${JM}
if [ ${JM} -lt 100 ]; then
export prefixGrid=000${JM}
if [ ${JM} -lt 10 ]; then
export prefixGrid=0000${JM}
fi
fi
fi
fi
#
#set run date
export DATA=2004032600

echo " ***** Gerando arquivos do pre-processamento do Global *****" 
#
# Set <pre-scripts>=1 to execute or <pre-scripts>=0 not execute
# 
TopoWaterPercNavy=0
TopoWaterPercGT30=0 
LandSeaMask=1
Chopping=1
VarTopo=1
TopoSpectral=1
VegetationMaskSSiB=1
VegetationMask=1
VegetationMaskSiB2Clima=1
VegetationMaskSiB2=1
VegetationMaskIBISClima=1
VegetationMaskIBIS=1
VegetationAlbedoSSiB=1
DeepSoilTemperatureClima=1
DeepSoilTemperature=1
RoughnessLengthClima=1
RoughnessLength=1
SoilMoistureClima=1
SoilMoisture=1
AlbedoClima=1
Albedo=1
SnowClima=1
SSTClima=1
SSTWeeklyNCEP=1 
SSTWeekly=1
CLimaSoilMoistureClima=1
CLimaSoilMoisture=1
ClmtClima=1
Clmt=1
DeltaTempColdestClima=1
DeltaTempColdest=1
NDVIClima=1
NDVI=1
PorceClayMaskIBISClima=1
PorceClayMaskIBIS=1
PorceClayMaskSiB2Clima=1
PorceClayMaskSiB2=1
PorceSandMaskIBISClima=1
PorceSandMaskIBIS=1
PorceSandMaskSiB2Clima=1
PorceSandMaskSiB2=1
SoilTextureMaskSiB2Clima=1
SoilTextureMaskSiB2=1
SSTMonthlyDirec=1



if [ $TopoWaterPercNavy -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_TopoWaterPercNavy.ksh 1 1 TopoNavy hold
   ${dirhome}/scripts/run_TopoWaterPercNavy.ksh  1 1 TopoNavy hold
   file_out=${dirdata}/pre/dataout/TopoNavy.dat
   file_out2=${dirdata}/pre/dataout/WaterNavy.dat
   sleep 5
   if test ! -s ${file_out} ; then 
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else if test ! -s ${file_out2} ; then
      echo "Problema!!! Nao foi gerado ${file_out2} "
      exit 1
   else
      echo
      echo "Fim do TopoWaterPercNavy"
      echo
   fi
   fi
fi

if [ $TopoWaterPercGT30 -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_TopoWaterPercGT30.ksh  1 1 GT30 hold
   ${dirhome}/scripts/run_TopoWaterPercGT30.ksh 1 1 GT30 hold
   file_out=${dirdata}/pre/dataout/TopoGT30.dat
   file_out2=${dirdata}/pre/dataout/WaterGT30.dat
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else if test ! -s ${file_out2} ; then
      echo "Problema!!! Nao foi gerado ${file_out2} "
      exit 1
   else
      echo
      echo "Fim do TopoWaterPercGT30"
      echo
   fi
   fi
fi

if [ $LandSeaMask -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_LandSeaMask.ksh  1 1 SeaMask hold 
   ${dirhome}/scripts/run_LandSeaMask.ksh  1 1 SeaMask hold
   file_out=${dirdata}/pre/dataout/LandSeaMaskNavy.G${prefixGrid}.dat
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do LandSeaMask"
      echo
   fi
fi

if [ $Chopping -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_Chopping.ksh  1 1 Chopping hold
   ${dirhome}/scripts/run_Chopping.ksh  1 1 Chopping hold
   file_out=${dirdata}/model/datain/GANLNMC${DATA}S.unf.${TRUNC}L0${KMOUT}
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do Chopping"
      echo
   fi
fi 
if [ $VarTopo -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_VarTopo.ksh  1 1 VarTop hold
   ${dirhome}/scripts/run_VarTopo.ksh  1 1 VarTop  hold
   file_out=${dirdata}/pre/dataout/Topography.G${prefixGrid}
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do VarTopo"
      echo
   fi
fi

if [ $TopoSpectral -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_TopoSpectral.ksh   1 1 TopoSpectra hold
   ${dirhome}/scripts/run_TopoSpectral.ksh    1 1 TopoSpectra hold
   file_out=${dirdata}/model/datain/TopoVariance.G${prefixGrid}
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do TopoSpectral"
      echo
   fi
fi 

if [ $VegetationMaskSSiB -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_VegetationMaskSSiB.ksh 1 1 VegMasSSiB hold
   ${dirhome}/scripts/run_VegetationMaskSSiB.ksh  1 1 VegMasSSiB hold
   file_out=${dirdata}/pre/dataout/VegetationMaskClima.dat
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do VegetationMaskSSiB"
      echo
   fi
fi 

if [ $VegetationMask -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_VegetationMask.ksh  1 1 VegetationMask  hold
   ${dirhome}/scripts/run_VegetationMask.ksh   1 1 VegetationMask  hold
   file_out=${dirdata}/pre/dataout/VegetationMask.G${prefixGrid}
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do VegetationMask"
      echo
   fi
fi 

if [ $VegetationMaskSiB2Clima -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_VegetationMaskSiB2Clima.ksh 1 1 VegMasSiB2 hold
   ${dirhome}/scripts/run_VegetationMaskSiB2Clima.ksh  1 1 VegMasSiB2 hold
   file_out=${dirdata}/pre/dataout/VegetationMaskClima.dat
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do VegetationMaskSiB2Clima"
      echo
   fi
fi 

if [ $VegetationMaskSiB2 -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_VegetationMaskSiB2.ksh  1 1 VegetationMaskSiB2  hold
   ${dirhome}/scripts/run_VegetationMaskSiB2.ksh   1 1 VegetationMaskSiB2  hold
   file_out=${dirdata}/pre/dataout/VegetationMaskSiB2.G${prefixGrid}
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do VegetationMaskSiB2"
      echo
   fi
fi 

if [ $VegetationMaskIBISClima -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_VegetationMaskIBISClima.ksh 1 1 VegMasIBIS hold
   ${dirhome}/scripts/run_VegetationMaskIBISClima.ksh  1 1 VegMasIBIS hold
   file_out=${dirdata}/pre/dataout/VegetationMaskIBISClima.dat
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do VegetationMaskIBISClima"
      echo
   fi
fi 

if [ $VegetationMaskIBIS -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_VegetationMaskIBIS.ksh 1 1 VegMasIBIS hold
   ${dirhome}/scripts/run_VegetationMaskIBIS.ksh  1 1 VegMasIBIS hold
   file_out=${dirdata}/pre/dataout/VegetationMaskIBISClima.dat
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do VegetationMaskIBIS"
      echo
   fi
fi 

if [ $VegetationAlbedoSSiB -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_VegetationAlbedoSSiB.ksh  1 1 VegAlbSSiB  hold
   ${dirhome}/scripts/run_VegetationAlbedoSSiB.ksh   1 1 VegAlbSSiB  hold
   file_out=${dirdata}/model/datain/VegetationMask.G${prefixGrid}
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do VegetationAlbedoSSiB"
      echo
   fi
fi 

if [ $DeepSoilTemperatureClima -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_DeepSoilTemperatureClima.ksh 1 1 DeepSoilTemClim hold
   ${dirhome}/scripts/run_DeepSoilTemperatureClima.ksh 1 1 DeepSoilTemClim hold
   file_out=${dirdata}/pre/dataout/DeepSoilTemperatureClima.dat
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do DeepSoilTemperatureClima"
      echo
   fi
fi 

if [ $DeepSoilTemperature -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_DeepSoilTemperature.ksh  1 1 DeepSoilTemp hold
   ${dirhome}/scripts/run_DeepSoilTemperature.ksh 1 1 DeepSoilTemp hold
   file_out=${dirdata}/model/datain/DeepSoilTemperature.G${prefixGrid}
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do DeepSoilTemperature"
      echo
   fi
fi 

if [ $RoughnessLengthClima -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_RoughnessLengthClima.ksh 1 1 RouLenClm hold
   ${dirhome}/scripts/run_RoughnessLengthClima.ksh 1 1 RouLenClm hold 
   file_out=${dirdata}/pre/dataout/RoughnessLengthClima.dat
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do RoughnessLengthClima"
      echo
   fi
fi 

if [ $RoughnessLength -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_RoughnessLength.ksh  1 1 RougLeng hold
   ${dirhome}/scripts/run_RoughnessLength.ksh 1 1 RougLeng hold
   file_out=${dirdata}/model/datain/RoughnessLength.G${prefixGrid}
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do RoughnessLength"
      echo
   fi
fi 

if [ $SoilMoistureClima -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_SoilMoistureClima.ksh  1 1 SoilMoisClm hold
   ${dirhome}/scripts/run_SoilMoistureClima.ksh 1 1 SoilMoisClm hold
   file_out=${dirdata}/pre/dataout/SoilMoistureClima.dat
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do SoilMoistureClima"
      echo
   fi
fi 

if [ $SoilMoisture -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_SoilMoisture.ksh  1 1 SoilMoist hold
   ${dirhome}/scripts/run_SoilMoisture.ksh 1 1 SoilMoist hold
   file_out=${dirdata}/model/datain/SoilMoisture.G${prefixGrid}
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do SoilMoisture"
      echo
   fi
fi 

if [ $AlbedoClima -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_AlbedoClima.ksh 1 1 AlbClm hold
   ${dirhome}/scripts/run_AlbedoClima.ksh 1 1 AlbClm hold
   file_out=${dirdata}/pre/dataout/AlbedoClima.dat
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do AlbedoClima"
      echo
   fi
fi 

if [ $Albedo -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_Albedo.ksh 1 1 Alb hold
   ${dirhome}/scripts/run_Albedo.ksh 1 1 Alb hold
   file_out=${dirdata}/pre/dataout/Albedo.G${prefixGrid}
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do Albedo"
      echo
   fi
fi 

if [ $SnowClima -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_SnowClima.ksh  1 1 SnowClm hold 
   ${dirhome}/scripts/run_SnowClima.ksh 1 1 SnowClm hold 
   file_out=${dirdata}/model/datain/Snow${DATA}S.unf.G${prefixGrid}
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do SnowClima"
      echo
   fi
fi 

if [ $SSTClima -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_SSTClima.ksh   1 1 SSTClim hold
   ${dirhome}/scripts/run_SSTClima.ksh  1 1 SSTClim hold
   datt=`echo ${DATA} |cut -c 1-8`
   file_out=${dirdata}/model/datain/SSTClima$datt.G${prefixGrid}
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do SSTClima"
      echo
   fi
fi 

if [ $SSTWeeklyNCEP -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_SSTWeeklyNCEP.ksh 1 1 SSTWeeklyNCEP hold
   ${dirhome}/scripts/run_SSTWeeklyNCEP.ksh 1 1 SSTWeeklyNCEP hold
   datt=`echo ${DATA} |cut -c 1-8`
   file_out=${dirdata}/pre/dataout/SSTWeekly.$datt
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else

      echo
      echo "Fim do SSTWeeklyNCEP"
      echo
   fi
fi 

if [ $SSTWeekly -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_SSTWeekly.ksh  1 1 SSTWeekly hold
   ${dirhome}/scripts/run_SSTWeekly.ksh 1 1 SSTWeekly hold
   datt=`echo ${DATA} |cut -c 1-8`
   file_out=${dirdata}/model/datain/SSTWeekly$datt.G${prefixGrid}
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do SSTWeekly"
      echo
   fi
fi 

if [ $CLimaSoilMoistureClima -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_CLimaSoilMoistureClima.ksh 1 1 CLimaSoilMoistureClima hold
   ${dirhome}/scripts/run_CLimaSoilMoistureClima.ksh 1 1 CLimaSoilMoistureClima hold
   file_out=${dirdata}/pre/dataout/CLimaSoilMoistureClima.dat
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do CLimaSoilMoistureClima"
      echo
   fi
fi 

if [ $CLimaSoilMoisture -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_CLimaSoilMoisture.ksh 1 1 CLimaSoilMoisture hold
   ${dirhome}/scripts/run_CLimaSoilMoisture.ksh 1 1 CLimaSoilMoisture hold
   file_out=${dirdata}/model/datain/CLimaSoilMoisture.G${prefixGrid}
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do CLimaSoilMoisture"
      echo
   fi
fi 

if [ $ClmtClima -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_ClmtClima.ksh 1 1 ClmtClima hold
   ${dirhome}/scripts/run_ClmtClima.ksh 1 1 ClmtClima hold
   file_out=${dirdata}/pre/dataout/TemperatureClima.dat
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do ClmtClima"
      echo
   fi
fi 

if [ $Clmt -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_Clmt.ksh 1 1 Clmt hold
   ${dirhome}/scripts/run_Clmt.ksh 1 1 Clmt hold
   file_out=${dirdata}/model/datain/Temperature.G${prefixGrid}
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do Clmt"
      echo
   fi
fi 

if [ $DeltaTempColdestClima -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_DeltaTempColdestClima.ksh 1 1 DeltaTempColdestClima hold
   ${dirhome}/scripts/run_DeltaTempColdestClima.ksh 1 1 DeltaTempColdestClima hold
   file_out=${dirdata}/pre/dataout/DeltaTempColdestClima.dat
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do DeltaTempColdestClima"
      echo
   fi
fi 

if [ $DeltaTempColdest -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_DeltaTempColdest.ksh 1 1 DeltaTempColdest hold
   ${dirhome}/scripts/run_DeltaTempColdest.ksh 1 1 DeltaTempColdest hold
   file_out=${dirdata}/model/datain/DeltaTempColdes.G${prefixGrid}
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do DeltaTempColdest"
      echo
   fi
fi 

if [ $NDVIClima -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_NDVIClima.ksh 1 1 NDVIClima hold
   ${dirhome}/scripts/run_NDVIClima.ksh 1 1 NDVIClima hold
   file_out=${dirdata}/pre/dataout/NDVIClima.dat
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do NDVIClima"
      echo
   fi
fi 

if [ $NDVI -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_NDVI.ksh 1 1 NDVIClima hold
   ${dirhome}/scripts/run_NDVI.ksh 1 1 NDVI hold
   file_out=${dirdata}/model/datain/NDVI.G${prefixGrid}
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do NDVI"
      echo
   fi
fi 

if [ $PorceClayMaskIBISClima -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_PorceClayMaskIBISClima.ksh 1 1 PorceClayMaskIBISClima hold
   ${dirhome}/scripts/run_PorceClayMaskIBISClima.ksh 1 1 PorceClayMaskIBISClima hold
   file_out=${dirdata}/pre/dataout/PorceClayMaskIBISClima.dat
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do PorceClayMaskIBISClima"
      echo
   fi
fi 

if [ $PorceClayMaskIBIS -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_PorceClayMaskIBIS.ksh 1 1 PorceClayMaskIBIS hold
   ${dirhome}/scripts/run_PorceClayMaskIBIS.ksh 1 1 PorceClayMaskIBIS hold
   file_out=${dirdata}/model/datain/PorceClayMaskIBIS.G${prefixGrid}
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do PorceClayMaskIBIS"
      echo
   fi
fi 

if [ $PorceClayMaskSiB2Clima -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_PorceClayMaskSiB2Clima.ksh 1 1 PorceClayMaskSiB2Clima hold
   ${dirhome}/scripts/run_PorceClayMaskSiB2Clima.ksh 1 1 PorceClayMaskSiB2Clima hold
   file_out=${dirdata}/pre/dataout/PorceClayMaskSiB2Clima.dat
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do PorceClayMaskSiB2Clima"
      echo
   fi
fi 
if [ $PorceClayMaskSiB2 -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_PorceClayMaskSiB2.ksh 1 1 PorceClayMaskSiB2 hold
   ${dirhome}/scripts/run_PorceClayMaskSiB2.ksh 1 1 PorceClayMaskSiB2 hold
   file_out=${dirdata}/model/datain/PorceClayMaskSiB2.G${prefixGrid}
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do PorceClayMaskSiB2"
      echo
   fi
fi 
if [ $PorceSandMaskIBISClima -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_PorceSandMaskIBISClima.ksh 1 1 PorceSandMaskIBISClima hold
   ${dirhome}/scripts/run_PorceSandMaskIBISClima.ksh 1 1 PorceSandMaskIBISClima hold
   file_out=${dirdata}/pre/dataout/PorceSandMaskIBISClima.dat
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do PorceSandMaskIBISClima"
      echo
   fi
fi 

if [ $PorceSandMaskIBIS -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_PorceSandMaskIBIS.ksh 1 1 PorceSandMaskIBIS hold
   ${dirhome}/scripts/run_PorceSandMaskIBIS.ksh 1 1 PorceSandMaskIBIS hold
   file_out=${dirdata}/model/datain/PorceSandMaskIBIS.G${prefixGrid}
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do PorceSandMaskIBIS"
      echo
   fi
fi 
if [ $PorceSandMaskSiB2Clima -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_PorceSandMaskSiB2Clima.ksh 1 1 PorceSandMaskSiB2Clima hold
   ${dirhome}/scripts/run_PorceSandMaskSiB2Clima.ksh 1 1 PorceSandMaskSiB2Clima hold
   file_out=${dirdata}/pre/dataout/PorceSandMaskSiB2Clima.dat
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do PorceSandMaskSiB2Clima"
      echo
   fi
fi 
if [ $PorceSandMaskSiB2 -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_PorceSandMaskSiB2.ksh 1 1 PorceSandMaskSiB2 hold
   ${dirhome}/scripts/run_PorceSandMaskSiB2.ksh 1 1 PorceSandMaskSiB2 hold
   file_out=${dirdata}/model/datain/PorceSandMaskSiB2.G${prefixGrid}
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do PorceSandMaskSiB2"
      echo
   fi
fi 
if [ $SoilTextureMaskSiB2Clima -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_SoilTextureMaskSiB2Clima.ksh 1 1 SoilTextureMaskSiB2Clima hold
   ${dirhome}/scripts/run_SoilTextureMaskSiB2Clima.ksh 1 1 SoilTextureMaskSiB2Clima hold
   file_out=${dirdata}/pre/dataout/SoilTextureMaskClima.dat
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do SoilTextureMaskSiB2Clima"
      echo
   fi
fi 
if [ $SoilTextureMaskSiB2 -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_SoilTextureMaskSiB2.ksh 1 1 SoilTextureMaskSiB2 hold
   ${dirhome}/scripts/run_SoilTextureMaskSiB2.ksh 1 1 SoilTextureMaskSiB2 hold
   file_out=${dirdata}/model/datain/SoilTextureMaskSiB2.G${prefixGrid}
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do SoilTextureMaskSiB2"
      echo
   fi
fi 

if [ $SSTMonthlyDirec -eq 1 ]; then
   echo executando: ${dirhome}/scripts/run_SSTMonthlyDirec.ksh 1 1 SSTMonthlyDirec hold
   ${dirhome}/scripts/run_SSTMonthlyDirec.ksh 1 1 SSTMonthlyDirec hold
   DATA1=`echo ${DATA} | awk '{print substr($1,1,8)}'`
   file_out=${dirdata}/model/datain/SSTMonthlyDirec${DATA1}.G${prefixGrid}
   sleep 5
   if test ! -s ${file_out} ; then
      echo "Problema!!! Nao foi gerado ${file_out} "
      exit 1
   else
      echo
      echo "Fim do SSTMonthlyDirec"
      echo
   fi
fi 
