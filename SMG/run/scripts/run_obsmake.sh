#! /bin/bash 
#-----------------------------------------------------------------------------#
#           Group on Data Assimilation Development - GDAD/CPTEC/INPE          #
#-----------------------------------------------------------------------------#
#BOP
#
# !SCRIPT: script utilizado para extrair as observações para os ciclos de
#          assimilação de dados
#
# !DESCRIPTION:
#
# !CALLING SEQUENCE:
#
#   ./run_obsmake.sh <opções>
#
#      As <opções> válidas são
#          * START_DATE : Data da condição inicial
#
#          exemplo:
#          ./run_obsmake.sh 2015043006
#
# !REVISION HISTORY:
#
# !REMARKS:
#
#EOP
#-----------------------------------------------------------------------------#
#BOC

# Carregando as variaveis do sistema
SCRIPT_PATH="$(realpath "${BASH_SOURCE[0]}")"
RootDir="$(dirname "$SCRIPT_PATH")"
export SMG_ROOT=${RootDir}
source ${SMG_ROOT}/../../config_smg.ksh vars_export

#-----------------------------------------------------------------------------#
# return usage from main program
#-----------------------------------------------------------------------------#
usage() {
   echo
   echo "Usage:"
   sed -n '/^#BOP/,/^#EOP/{/^#BOP/d;/^#EOP/d;p}' ${BASH_SOURCE}
}

# Data da condição inicial
if [ -z ${1} ]; then
   echo -e "\e[31;1m >> Erro: \e[m\e[33;1m Data da condição inicial do modelo não foi passada\e[m"
   usage
   exit -1
fi

export START_DATE=${1}

DATE=${START_DATE}

echo ""
echo -e "\033[34;1m >>  Descompactando Observacoes... \033[m"
echo ""

YMD=${START_DATE:0:8}
HH=${START_DATE:8:10}
YM=${START_DATE:0:6}
DH=${START_DATE:6:10}
Y=${START_DATE:0:4}
M=${START_DATE:4:2}
D=${START_DATE:6:2}
#DIRFILES=${ncep_ext}/ASSIMDADOS
DIRFILES=${ncep_ext}/${Y}/${M}/${D}
    echo ${subt_gsi_datain_obs}

count=0
#ls -1 ${DIRFILES}/*${YMD}*.gz | while read file; do
#    echo -e "\e[32;1m${file}\e[m"
#    tar -xvzf ${file} -C ${subt_gsi_datain_obs}
#    if [ $? -eq 0 ];then
#       count=$((count+1))
#    fi
#done

if [ ! -d ${subt_gsi_datain_obs} ]; then mkdir -p ${subt_gsi_datain_obs}; fi

for obsfile in $(find ${DIRFILES} -type f -size +0c -name "gdas.*")
do
  ln -sfv ${obsfile} ${subt_gsi_datain_obs}/      
  #cp -v ${obsfile} ${subt_gsi_datain_obs}/      
  count=$((count+1))
done        

if [ ${count} -gt 0 ];then
   echo -e ""
   echo -e "\e[34;1m Foram obtidos\e[m \e[37;1m${count}\e[m \e[34;1marquivos.\e[m"
   echo -e ""
   exit 0
else
   echo -e "\033[31;1m !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \033[m"
   echo -e "\033[31;1m !!! Nenhum Arquivo de Observacoes disponível !!! \033[m"
   echo -e "\033[31;1m !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \033[m"
   exit 1
fi

## Primeiro verifica se os arquivos existem no seguinte diretório
##
## Obs.: Trocar pelo diretorio da DMD
#DMD='/stornext/online6/das/bruna.silveira/dados_testecase_AD'
#DMD='/stornext/online6/das/gdad/OBS_prepbufr'
#
#DirFiles=${DMD}/${YM}/${DH}
#
#if [ "$(ls -A ${DirFiles})" ]; then
#   count=0
#   while read line; do
#      count=$((count+1))
#      echo -e "\e[32;1m$line\e[m"
#      cp -pfr ${DirFiles}/${line} ${subt_gsi_datain_obs}
#   done < <(ls -1 ${DirFiles})
#
#   echo -e ""
#   echo -e "\e[34;1m Foram obtidos\e[m \e[37;1m${count}\e[m \e[34;1marquivos.\e[m"
#   echo -e ""
#else
#
#   # Descompacta os arquivos de observacao que estao em $ncep_ext:
##   file=${ncep_ext}/${YM}/${DH}/gdas1.bufr_${START_DATE}.tgz
##   file=${ncep_ext}/ASSIMDADOS/gdas1.bufr_${START_DATE}.tgz
#
#   if [ -e ${file} ]; then
#     count=0
#     while read line; do
#        count=$((count + 1))
#        echo -e "\e[32;1m$line\e[m"
#     done  < <(tar -zxvf ${file} -C ${subt_gsi_datain_obs})
#
#     echo -e ""
#     echo -e "\e[34;1m Foram obtidos\e[m \e[37;1m${count}\e[m \e[34;1marquivos.\e[m"
#     echo -e ""
#   else
#      echo ""
#      echo -e "\033[31;1m !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \033[m"
#      echo -e "\033[31;1m !!! Arquivo de Observacoes não disponível !!! \033[m"
#      echo ""
#      echo -e "\033[32;1m ${file} \033[m"
#      echo ""
#      echo -e "\033[31;1m !!! Abortando .....                       !!! \033[m"
#      echo -e "\033[31;1m !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \033[m"
#      echo ""
#      exit 1
#   fi
#fi
#
#exit 0
