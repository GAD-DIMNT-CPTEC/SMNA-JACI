#! /bin/bash

# Descomentar para debugar
#set -o xtrace

### Define hpc_name below is needed because in some cases (ex.: First Run),
### this script may be called from the command line.
### 
lognode=`cat /proc/sys/kernel/hostname | cut  -b 1-6`

case $lognode in

  clogin)
    STR=`uname -a`
    SUB='cray'
    if [[ "$STR" == *"$SUB"* ]]; then
      echo -n "This will run on cray XC50 ..."
      export hpc_name="XC50"
    fi
    ;;

  headno)
    STR=`uname -a`
    SUB='egeon'
    if [[ "$STR" == *"$SUB"* ]]; then
      echo -n "This will run on EGEON Cluster ..."
      export hpc_name="egeon"
    fi
    ;;

  *)
    mach=`cat /proc/sys/kernel/hostname`
    echo -n "The configurations for "$mach" is not defined yet !"
    echo -n "1) Add the machine to the defined systems in etc/mach ; and"
    echo -n "2) add an option for it in the function copy_fixed_files in etc/functions"
    exit
    ;;
esac

# Carregando as variaveis do sistema
SCRIPT_PATH="$(realpath "${BASH_SOURCE[0]}")"
RootDir="$(dirname "$SCRIPT_PATH")"
export SMG_ROOT=${RootDir}
source ${SMG_ROOT}/../../config_smg.ksh vars_export
source ${SMG_ROOT}/run/smg_functions.sh

# carregando funcoes do pre-processamento

source ${home_run_bam}/runPre.func

# Verificando argumentos de entrada
if [ -z "${1}" ]
then
  echo "LABELANL is not set" 
  exit 3
else
  export LABELANL=${1}
fi
if [ -z "${2}" ]
then
  echo "LABELFCT is not set" 
  exit 3
else
  export LABELFCT=${2} 
fi
if [ -z "${3}" ]
then
  echo "PREFIX is not set" 
  exit 3
else
  export PREFIX=${3}
fi
if [ -z "${4}" ]
then
  echo "TRC is not set" 
  exit 3
else
  export TRC=${4}
fi
if [ -z "${5}" ]
then
  echo "NLV is not set" 
  exit 3
else
  export NLV=${5}
fi
if [ -z "${6}" ]
then
  echo "NPROC is not set"
  case ${hpc_name} in
     	egeon)  echo "setting to 64"
             	export NPROC=64  # ntasks
	;;
  	XC50)  	echo "setting to 480"
  		export NPROC=480             	
	;;
  esac
else
  export NPROC=${6}
fi
if [ "$#" == 7 ]
then 
  export RUNPOS=$(echo ${7} | tr [:upper:] [:lower:])   
else 
  export RUNPOS="yes"
fi

case ${hpc_name} in
   egeon) tasks_per_node=16
	        cpus_per_task=8
	;;
   XC50)  tasks_per_node=10
	        cpus_per_task=4
   ;;
esac

getBAMSize ${TRC}
export postfix=$(printf "G%5.5dL%3.3d \n" $JM $NLV)
export MRES=`echo ${TRC} ${NLV} | awk '{printf("TQ%4.4dL%3.3d\n",$1,$2)}'`
export ANL="GANL${PREFIX}${LABELANL}S.unf.${MRES}"
### export LABELFGS=`${inctime} ${LABELANL} +6h %y4%m2%d2%h2`
export LABELFGS=`date -u +%Y%m%d%H -d "${LABELANL:0:8} ${LABELANL:8:2} +6 hours" `
### export YESTERDAY=$(${inctime} ${LABELANL} -1d %y4%m2%d200)
export YESTERDAY=`date -u +%Y%m%d00 -d "${LABELANL:0:8} ${LABELANL:8:2} -1 days" `
export modelDataIn=${subt_model_bam}/datain
export gsiDataOut=${subt_gsi_dataout}/${LABELANL}

echo -e ""
echo -e "\033[34;1m >> Submetendo o MCGA:\033[m \033[31;1m${LABELANL}\033[m"
echo -e " "
echo -e "\033[34;1m > Resolucao (espectral) : \033[m \033[31;1m${MRES}\033[m"
echo -e "\033[34;1m > Resolucao (grade)     : \033[m \033[31;1m${postfix}\033[m"
echo -e "\033[34;1m > Condicao Inicial      : \033[m \033[31;1m${LABELANL}\033[m"
echo -e "\033[34;1m > Previsao ate          : \033[m \033[31;1m${LABELFCT}\033[m"
echo -e "\033[34;1m > Pos-proc. Ativado     : \033[m \033[31;1m${RUNPOS}\033[m"

#
# mudando para diretorio dos scripts do BAM
#

cd ${home_run_bam}

#
# rodando o somente o Chopping do pré para pegar o arquivo de ozônio
# saida gerada no bam/model/datain/ 

# only BAM /bin/bash runPre -v -t ${TRC} -l ${NLV} -I ${LABELANL} -s -n chp -O

# only BAM STATUS=$?
# only BAM echo "1st call to runPre. Status: "${STATUS}
# only BAM if [ ${STATUS} -ne 0 ];then
# only BAM    exit ${STATUS}
# only BAM fi

# only BAM cp -f ${modelDataIn}/OZONSMT${LABELANL}S.grd.${postfix} ${modelDataIn}/OZON${PREFIX}${LABELANL}S.grd.${postfix}

#
# remove arquivos desnecessarios
#

# only BAM rm -fr ${modelDataIn}/GANLSMT${LABELANL}S.unf.*
# only BAM rm -fr ${modelDataIn}/OZONSMT${LABELANL}S.unf.*

#
# Copiando Analise do GSI para o Model/DataIn
#

# only BAM cp -pfr ${gsiDataOut}/GANL${PREFIX}${LABELANL}S.unf.${MRES} ${modelDataIn}

#
# Rodando os demais processos do pré e usando a análise do GSI
#
# /bin/bash runPre -v -t 299 -l 64 -I ${LABELANL}  -n 0 -p SMT -s -O -T -G -Gp gblav -Gt Grid

# only BAM /bin/bash runPre -v -t ${TRC} -l ${NLV} -I ${LABELANL} -p CPT -n das
# only BAM STATUS=$?
# only BAM echo "2nd call to runPre. Status: "${STATUS}
# only BAM if [ ${STATUS} -ne 0 ];then
# only BAM    exit ${STATUS}
# only BAM fi

# Rodando o Modelo  
# Obs.: This script is for long run, e.g. 240 h, when data assimilition cycle was already done
#       then, here "runModel" is called with -f (forecast only mode)

/bin/bash runModel -f -das -v -np ${NPROC} -N ${tasks_per_node} -d ${cpus_per_task} \
                   -t ${TRC} -l ${NLV} -I ${LABELANL} -F ${LABELFCT} -W  ${LABELFCT} \
                   -p ${PREFIX} -s sstwkl -ts 6 -tr 6 -i 2 -s sstwkl

                                               ## -r  removed to avoid writing restart files

# Pos-processa as previsoes caso a variavel RUNPOS possua o valor Yes ou Y
if [ ${RUNPOS} == "yes" -o ${RUNPOS} == "y" ]
then

  # Verifica se o executavel se encontra em ${subt_pos_bam_run}
  if [ ! -e ${subt_pos_bam}/exec/PostGrib ]
  then

    cp -v ${home_pos_bam}/exec/PostGrib ${subt_pos_bam}/exec

  fi

  cd ${home_run_bam}
  
## NOTES FOR PostGrib 
## @egeon: ntasks-per-node  must be <=6 (6, 4), 8 or more creates some corrupted grib files
## @egeon: ThreadsPerMPITask do not seems to speed the processing time  

  ###       -np  (MPI Tasks) ; -N (tasks_per node)  ; -d (ThreadsperMPItasks)

  echo   "./runPos -np 6 -N 6 -d 8 -t ${TRC} -l ${NLV} -I ${LABELANL} -F ${LABELFCT} -p ${PREFIX} > /dev/null 2>&1"
  /bin/bash ./runPos -np 6 -N 6 -d 8 -t ${TRC} -l ${NLV} -I ${LABELANL} -F ${LABELFCT} -p ${PREFIX} 
  ## > /dev/null 2>&1
  # STATUS=$?
  # if [ ${STATUS} -ne 0 ];then
  #    exit ${STATUS}
  # fi

fi

exit 0
