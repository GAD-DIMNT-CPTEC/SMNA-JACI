#! /bin/bash 
#-------------------------------------------------------------------------------------------------#
#                         Brazilian Global Atmospheric Model - BAM_V2.2.1                         #
#-------------------------------------------------------------------------------------------------#
# Descrição:                                                                                      #
#     Script para executar o ciclo (pré-processamento + modelo) do BAM                            #
#                                                                                                 #
# Uso:                                                                                            #
#     ./run_cycle.sh LABELI LABELF                                                                #
#                                                                                                 #
# Exemplo:                                                                                        #
#     ./run_cycle.sh 2020061300 2020062300                                                        #
#                                                                                                 #
#revisões:                                                                                        #
#     * 22-04-2021: Khamis, E. G.  - código inicial                                               #
#     * 30-04-2021: Khamis, E. G.  - inclusao de funcoes e correcoes nas dependencias             #
#                                    dos arquivos fixos gerados no pre em sua primeira execucao    #
#     * 13-07-2021: Khamis, E. G.  - alteracoes para a inclusao dos processos CO2MonthlyDirec,    #
#                                    SSTDailyDirec e SSTSeasonDirec                               #
#     * 22-11-2023: Aravequia, J. A. - Adicionada opção headnode (EGEON)                          #
#                                                                                                 #
# TODO:                                                                                           #
#     * adicionar o pos-processamento ao ciclo                                                    #
#                                                                                                 #
# DIMNT/CPTEC/INPE, 2021                                                                          #
#-------------------------------------------------------------------------------------------------#

#
# Funções comuns
#

function usageprincipal()
{
  echo ""
  echo -e "\033[34;1m Brazilian Global Atmospheric Model - BAM_V2.2.1\033[m"
  echo ""
  echo " Descrição:"
  echo "     Script para executar o ciclo do BAM"
  echo ""
  echo " Uso:"
  echo "     ./run_cycle.sh LABELI LABELF"
  echo ""
  echo " Exemplo:"
  echo "     ./run_cycle.sh 2020061300 2020062300"
  echo ""
  echo " Opções:"
  echo "     * LABELI....: data da análise inicial" 
  echo "     * LABELF....: data da análise final" 
  echo ""
  echo -e "\033[33;1m DIMNT/CPTEC/INPE, 2021 \033[m"
  echo ""    
}

function rmFiles()
{
  pre_out_aux=$1
  model_in_aux=$2

  rm -f ${pre_out_aux}/*.G0*
  rm -f ${pre_out_aux}/GANL*
  rm -f ${pre_out_aux}/OCMC*
  rm -f ${pre_out_aux}/CO2MonthlyDirec*
  rm -f ${pre_out_aux}/SSTDailyDirec*
  rm -f ${pre_out_aux}/SSTSeasonDirec*
  rm -f ${pre_out_aux}/SnowClima*
  rm -f ${pre_out_aux}/SSTWeekly*
  rm -f ${pre_out_aux}/SoilMoistureWeekly*
  rm -f ${pre_out_aux}/TopographyGradient*
  rm -f ${pre_out_aux}/AlbedoClima.dat
  rm -f ${pre_out_aux}/PorceSandMaskIBISClima.dat
  rm -f ${pre_out_aux}/SoilTextureMaskClima.dat
  rm -f ${pre_out_aux}/VegetationMaskClima2.dat
  rm -f ${pre_out_aux}/CLimaSoilMoistureClima.dat
  rm -f ${pre_out_aux}/NDVIClima.dat
  rm -f ${pre_out_aux}/PorceSandMaskIBISClimaG.dat
  rm -f ${pre_out_aux}/SoilTextureMaskClima2.dat
  rm -f ${pre_out_aux}/VegetationMaskClima2G.dat
  rm -f ${pre_out_aux}/DeepSoilTemperatureClima.dat
  rm -f ${pre_out_aux}/PorceClayMaskIBISClima.dat
  rm -f ${pre_out_aux}/PorceSandMaskSiB2Clima.dat
  rm -f ${pre_out_aux}/TemperatureClima.dat
  rm -f ${pre_out_aux}/VegetationMaskClimaG.dat
  rm -f ${pre_out_aux}/DeltaTempColdestClima.dat
  rm -f ${pre_out_aux}/PorceClayMaskIBISClimaG.dat
  rm -f ${pre_out_aux}/PorceSandMaskSiB2ClimaG.dat
  rm -f ${pre_out_aux}/VegetationMaskIBISClima.dat
  rm -f ${pre_out_aux}/PorceClayMaskSiB2Clima.dat
  rm -f ${pre_out_aux}/RoughnessLengthClima.dat
  rm -f ${pre_out_aux}/VegetationMaskIBISClimaG.dat
  rm -f ${pre_out_aux}/VegetationMaskIBISClimaG.dat
  rm -f ${pre_out_aux}/LandSeaMaskNavy.G00096.dat
  rm -f ${pre_out_aux}/PorceClayMaskSiB2ClimaG.dat
  rm -f ${pre_out_aux}/SoilMoistureClima.dat
  rm -f ${pre_out_aux}/VegetationMaskClima.dat

  rm -f ${model_in_aux}/GANL*
  rm -f ${model_in_aux}/CO2MonthlyDirec*
  rm -f ${model_in_aux}/SSTDailyDirec*
  rm -f ${model_in_aux}/SSTSeasonDirec*
  rm -f ${model_in_aux}/SSTWeekly*
  rm -f ${model_in_aux}/Snow*
  rm -f ${model_in_aux}/SoilMoistureWeekly*
  rm -f ${model_in_aux}/TopographyGradient*
  rm -f ${model_in_aux}/OZON*
  rm -f ${model_in_aux}/TRAC*
}

#
# Verifica o número de argumentos passados junto com o script:
# Pelo menos 2 argumentos devem ser passados;
# Caso contrário, será chamada a função para imprimir o cabeçalho e o script eh encerrado
#

if [ "$#" -lt 2 ]
then
  usageprincipal
  exit 0
fi

#
# Pegando os argumentos da linha de comando
#

if [ -z "${1}" ]
then
  echo "LABELI is not set" 
  exit 1
else
  export LABELI=${1}
fi

if [ -z "${2}" ]
then
  echo "LABELF is not set" 
  exit 2
else
  export LABELF=${2}
fi

#
# Ajuste das variáveis de ambiente
#

export inctime=${HOME}/bin/inctime
#export wgrib2=${HOME}/bin/wgrib2
#export oper_path=/lustre_xc50/ioper/data/external

export local_dir=$(dirname $(readlink -e ${0})) 
#export home_dir=/lustre_xc50/${USER}/up
export home_dir=${local_dir}
export dados=${home_dir}/dados
export bam_dir=${home_dir}/BAM_V2.2.1
export bam_run_dir=${home_dir}/BAM_V2.2.1/run
export pre_exe_dir=${home_dir}/bin
export model_exe_dir=${home_dir}/bin
#export pre_exe_dir=${bam_dir}/pre/build
#export model_exe_dir=${bam_dir}/model/source
export pre_exe_name=ParPre_MPI
export model_exe_name=ParModel_MPI

#datanow=`date  "+%Y%m%d%H"`
export datai=${LABELI} #2019112300
export dataf=${LABELF} #2019120100
export data=${datai}

#YYYYMMDD=`echo ${data} |cut -c1-8`
#hh=`echo ${data} |cut -c9-10`

#
#
# Início do ciclo
#
#

#  data=$(${inctime} ${data} +1d %y4%m2%d2%h2)

while [ ${data} -le ${dataf} ]
do

  # variaveis dependentes da data atual devem estar internas ao loop 
  export pre_in=${home_dir}/${data}/pre/datain
  export pre_out=${home_dir}/${data}/pre/dataout
  export pre_exec=${home_dir}/${data}/pre/exec
  export model_in=${home_dir}/${data}/model/datain
  export model_exec=${home_dir}/${data}/model/exec
  export data_run=${home_dir}/${data}/run
  
  mkdir -p ${pre_in}
  mkdir -p ${pre_out}
  mkdir -p ${pre_exec}
  mkdir -p ${model_in}
  mkdir -p ${model_exec}
  mkdir -p ${data_run}

  # remove saidas do pre-processamento caso tenha falhado anteriormente
  rmFiles ${pre_out} ${model_in}
  
  ln -s ${dados}/pre/datain/*           ${home_dir}/${data}/pre/datain/.
  ln -s ${dados}/pre/databcs            ${home_dir}/${data}/pre/databcs
  ln -s ${dados}/pre/datasst            ${home_dir}/${data}/pre/datasst
  ln -s ${dados}/pre/dataco2            ${home_dir}/${data}/pre/dataco2
  ln -s ${dados}/pre/dataTop            ${home_dir}/${data}/pre/dataTop
  ln -s ${dados}/model/datain/*         ${home_dir}/${data}/model/datain/.
  ln -s ${dados}/pre/dataout/HPRIME.dat ${home_dir}/${data}/pre/dataout/.
  
  # copia dados de entrada operacionais (trecho abaixo incluso no runPre)
  #ln -s ${oper_path}/${data}/dataout/NCEP/gblav.T00Z.atmanl.nemsio.${data} ${pre_in}
  #ln -s ${oper_path}/${data}/dataout/NCEP/rtgssthr_grb_0.083.grib2.${data:0:8} ${pre_in}
  #ln -s ${oper_path}/${data}/dataout/Umid_Solo/GL_SM.GPNR.${data}.vfm ${pre_in}
  #rec=`${wgrib2} -s -YY -d 1 -order we:ns ${pre_in}/rtgssthr_grb_0.083.grib2.${data:0:8} -ieee ${pre_in}/gdas1.T${data:8:10}Z.sstgrd.${data}`
  ##rec=`${wgrib2} -s -YY -d 1 -order we:ns ${pre_in}/gdas1.T${hh}Z.sstgrb2.${YYYYMMDD}${hh} -ieee  ${pre_in}/gdas1.T${hh}Z.sstgrd.${YYYYMMDD}${hh}`
  #echo ${rec}
  #sleep 5s
  #until [ -s ${pre_in}/gdas1.T00Z.sstgrd.${data} ]; do sleep 1s; done

  cp -pvf ${bam_run_dir}/* ${data_run}
  cp -pfv ${pre_exe_dir}/${pre_exe_name} ${home_dir}/${data}/pre/exec
  cp -pfv ${model_exe_dir}/${model_exe_name} ${home_dir}/${data}/model/exec
  sleep 5s

  sed -i "s;homebase_dir;${home_dir}/${data};g" ${data_run}/EnvironmentalVariables
  sed -i "s;subtbase_dir;${home_dir}/${data};g" ${data_run}/EnvironmentalVariables
  sed -i "s;workbase_dir;${home_dir}/${data};g" ${data_run}/EnvironmentalVariables

  echo "Submit to " $hpc_name
  if [ ${data} -eq ${datai} ] 
  then
    echo ""
    echo "Rodando o pre-processamento completo na resolucao TQ0299L64..." 
    # realiza o pre-processamento completo em sua primeira execucao
    # TODO: checar -O (gera Ozonio) e -T (gera tracadores)

    case ${hpc_name} in

      egeon)  /bin/bash ${data_run}/runPre -v -t 299 -l 64 -np 120 -N 120 -d 1 -I ${data} -p SMT -n 0 -s -G -ti 1534 -li 64 -Gt grid -Gp gblav -O -T
              ### EGEON has more cpus per node . -np Ntasks -N TasksPerNode -d ThreadsPerMPITask
              ## ajustes para rodar em 4 nos (Filas PESQ1 e PESQ2)
      ;;
      XC50)  /bin/bash ${data_run}/runPre -v -t 299 -l 64 -np 240 -N 40 -d 1 -I ${data} -p SMT -n 0 -s -G -ti 1534 -li 64 -Gt grid -Gp gblav -O -T
      ;;

    esac

  else
    echo ""
    echo "Rodando o pre-processamento parcial com os dados gerados na data inicial..." 
    # arquivos fixos gerados no pre' pela primeira vez (data inicial)
    ln -s ${home_dir}/${datai}/model/datain/AlbedoSSiB ${model_in} 
    ln -s ${home_dir}/${datai}/model/datain/DeepSoilTemperature.G* ${model_in}
    ln -s ${home_dir}/${datai}/model/datain/DeltaTempColdes.G* ${model_in} 
    ln -s ${home_dir}/${datai}/model/datain/PorceClayMaskIBIS.G* ${model_in} 
    ln -s ${home_dir}/${datai}/model/datain/PorceSandMaskIBIS.G* ${model_in} 
    ln -s ${home_dir}/${datai}/model/datain/RoughnessLength.G* ${model_in} 
    ln -s ${home_dir}/${datai}/model/datain/SoilMoisture.G* ${model_in} 
    ln -s ${home_dir}/${datai}/model/datain/Temperature.G* ${model_in} 
    ln -s ${home_dir}/${datai}/model/datain/temp_ltm_month.G* ${model_in} 
    ln -s ${home_dir}/${datai}/model/datain/TopographyRec.G* ${model_in} 
    ln -s ${home_dir}/${datai}/model/datain/TopoVariance.G* ${model_in} 
    ln -s ${home_dir}/${datai}/model/datain/VegetationMask.G* ${model_in} 
    ln -s ${home_dir}/${datai}/model/datain/VegetationMaskIBIS.G* ${model_in} 
    ln -s ${home_dir}/${datai}/model/datain/VegetationSSiB ${model_in} 
    ln -s ${home_dir}/${datai}/pre/dataout/Topography.TQ* ${pre_out} 
    ln -s ${home_dir}/${datai}/pre/dataout/Albedo.G* ${pre_out} # dependencia de SnowClima no pre'
    ln -s ${home_dir}/${datai}/pre/dataout/ModelLandSeaMask.G* ${pre_out} # dependencia de SSTWeekly no pre'
    ln -s ${home_dir}/${datai}/pre/dataout/LandSeaMask.G* ${pre_out} # dependencia de SnowClima no pre' (checar)
    # dependencias do pre' a ignorar no arquivo cfg a partir das rodadas seguintes
    # ignorar TopoSpectral no Chopping
    sed -i s,'modCfg(3)%dependencies','!modCfg(3)%dependencies',g   ${data_run}/PRE_cfg.nml_default
    # ignorar LandSeaMask e manter o Chopping no CO2MonthlyDirec
    sed -i   '/modCfg(4)%dep/s/,/, !/'                              ${data_run}/PRE_cfg.nml_default
    # ignorar LandSeaMask e manter o Chopping no OCMClima
    sed -i   '/modCfg(11)%dep/s/,/, !/'                             ${data_run}/PRE_cfg.nml_default
    # ignorar Albedo e LandSeaMask no SnowClima
    sed -i s,'modCfg(18)%dependencies','!modCfg(18)%dependencies',g ${data_run}/PRE_cfg.nml_default
    # ignorar LandSeaMask e manter o Chopping no SSTDailyDirec
    sed -i   '/modCfg(24)%dep/s/,/, !/'                             ${data_run}/PRE_cfg.nml_default
    # ignorar LandSeaMask e manter o Chopping no SSTSeasonDirec
    sed -i   '/modCfg(26)%dep/s/,/, !/'                             ${data_run}/PRE_cfg.nml_default
    # realiza o pre-processamento parcial - o modelo ira utilizar tambem as saidas fixas da primeira execucao

    case ${hpc_name} in

      egeon)  /bin/bash ${data_run}/runPre -v -t 299 -l 64 -np 120 -N 120 -d 1 -I ${data} -p SMT -n 2 -s -G -ti 1534 -li 64 -Gt grid -Gp gblav -O -T
              ### EGEON has more cpus per node . -np Ntasks -N TasksPerNode -d ThreadsPerMPITask
              ## ajustes para rodar em 4 nos (Filas PESQ1 e PESQ2)
      ;;
      XC50)  /bin/bash ${data_run}/runPre -v -t 299 -l 64 -np 240 -N 40 -d 1 -I ${data} -p SMT -n 2 -s -G -ti 1534 -li 64 -Gt grid -Gp gblav -O -T
      ;;

    esac
    /bin/bash ${data_run}/runPre -v -t 299 -l 64 -np 240 -N 40 -d 1 -I ${data} -p SMT -n 2 -s -G -ti 1534 -li 64 -Gt grid -Gp gblav -O -T
  fi

  # adequando os nomes das saidas do chopping para serem lidas pelo modelo
  ln -s ${model_in}/GANLSMT${data}S.unf.TQ0299L064 ${model_in}/GANLNMC${data}S.unf.TQ0299L064
  ln -s ${model_in}/GANLSMT${data}S.unf.TQ1024L064 ${model_in}/GANLNMC${data}S.unf.TQ1024L064
  ln -s ${model_in}/OZONSMT${data}S.grd.G00450L064 ${model_in}/OZONNMC${data}S.grd.G00450L064
  ln -s ${model_in}/OZONSMT${data}S.unf.TQ1024L064 ${model_in}/OZONNMC${data}S.unf.TQ1024L064
  ln -s ${model_in}/TRACSMT${data}S.grd.G00450L064 ${model_in}/TRACNMC${data}S.grd.G00450L064
  ln -s ${model_in}/TRACSMT${data}S.unf.TQ1024L064 ${model_in}/TRACNMC${data}S.unf.TQ1024L064
 
  # variavel auxiliar para calcular a data de uma previsao de 2 dias a partir da data corrente
  aux=$(${inctime} ${data} +2d %y4%m2%d2%h2)

  echo ""
  echo "Rodando o modelo na resolucao TQ0299L64..."
  # realiza a execucao do modelo
  case ${hpc_name} in

    egeon)  /bin/bash ${data_run}/runModel -v -np 64 -N 16 -d 8 -t 299 -l 64 -I ${data} -F ${aux} -W ${aux} -p NMC -s sstwkl -ts 24 -tr 6
             ### EGEON has more cpus per node . -np Ntasks -N TasksPerNode -d ThreadsPerMPITask
             ## ajustes para rodar em 4 nos (Filas PESQ1 e PESQ2)
    ;;
    XC50)  /bin/bash ${data_run}/runModel -v -np 800 -N 4 -d 10 -t 299 -l 64 -I ${data} -F ${aux} -W ${aux} -p NMC -s sstwkl -ts 24 -tr 6
    ;;

  esac
  data=$(${inctime} ${data} +1d %y4%m2%d2%h2)

done

exit 0
