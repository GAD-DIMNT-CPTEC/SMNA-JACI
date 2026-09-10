#!/bin/bash
#-----------------------------------------------------------------------------#
#           Group on Data Assimilation Development - GDAD/CPTEC/INPE          #
#-----------------------------------------------------------------------------#
#BOP
#
# !SCRIPT:run_blsdas
#
# !DESCRIPTION: script para rodar o blsdas no ciclo do GSI+BAM
#
#
#      ./runBldas <opções>
#
#         As <opções> válidas são:
#            * -t   <val> : truncamento [default: 62]
#            * -l   <val> : numero de niveis [default: 28]
#            * -p   <val> : prefixo dos arquivos do BAM (condição inicial e previsões) [default: CPT]
#            * -i   <val> : Data que deve gerar a condição inicial
#            * -Mnp <val> : numero de processadores [default: 72]
#            * -ts  <val> : TimeStep da previsão [default: 6]
#            * -h   <val> : mostra esta ajuda
#
#  example:
#
#     ./runBldas -t 299 -l 64 -d 2013010100 -Mnp 480

#
# !REVISION HISTORY:
# 07 May 2018 - J. G. de Mattos - Initial Version
#
# !REMARKS:
#
#
#EOP
#-----------------------------------------------------------------------------#
#BOC
# Carregando as variaveis do sistema
SCRIPT_PATH="$(realpath "${BASH_SOURCE[0]}")"
RootDir="$(dirname "$SCRIPT_PATH")"
export SMG_ROOT=${RootDir}
source ${SMG_ROOT}/../../config_smg.ksh vars_export

subwrd ( ) {
   str=$(echo "${@}" | awk '{ for (i=1; i<=NF-1; i++) printf("%s ",$i)}')
   n=$(echo "${@}" | awk '{ print $NF }')
   echo "${str}" | awk -v var=${n} '{print $var}'
}
# Local (path) onde está este script
LOCALDIR=$(dirname $(readlink -e ${0})) # Local (path) onde está este script

#
# Verificando o numero de argumentos
# caso nao existam argumentos, exibe help
#
if [ $# -eq 0 ];then
   cat < ${0} | sed -n '/^#BOP/,/^#EOP/p'
   exit 0
fi

#
# Pegando as opções que foram passadas pela linha de comando
#
i=1
flag=0
while [ 1 ]; do

   arg=$(echo "${@}" | awk -v var=${i} '{print $var}')
   i=$((i+1))

   if [ -z ${arg} ];then
      break;
   fi

   while [ 1 ];do

      if [ ${arg} = '-t' ]; then TRC=$(subwrd ${@} ${i}); i=$((i+1)); break; fi
      if [ ${arg} = '-l' ]; then LV=$(subwrd ${@} ${i}); i=$((i+1)); break; fi
      if [ ${arg} = '-p' ]; then PREFIX=$(subwrd ${@} ${i}); i=$((i+1)); break; fi
      if [ ${arg} = '-i' ]; then LABELI=$(subwrd ${@} ${i}); i=$((i+1)); break; fi
      if [ ${arg} = '-Mnp' ]; then MPITasks=$(subwrd ${@} ${i}); i=$((i+1)); break; fi
      if [ ${arg} = '-ts' ]; then DHFCT=$(subwrd ${@} ${i}); i=$((i+1)); break; fi
      if [ ${arg} = '-h' ]; then cat < ${0} | sed -n '/^#BOP/,/^#EOP/p' ; i=$((i+0)); exit 0; fi

      flag=1
      i=$((i-1))

      break;
   done
   if [ ${flag} -eq 1 ]; then break; fi
done

# truncamento
if [ -z ${TRC} ];then
   TRC=62
fi

# numero de niveis verticais
if [ -z ${LV} ];then
   LV=28
fi

# prefixo dos arquivos
if [ -z ${PREFIX} ];then
   PREFIX=CPT
fi

# Data da Condição Inicial (cold start)
if [ -z ${LABELI} ];then
   echo -e "\033[31;1m LABELI not set \033[m"
   exit 1
else
   ### LABELFGS=$(${inctime} ${LABELI} -6h %y4%m2%d2%h2)
   LABELFGS=`date -u +%Y%m%d%H -d "${LABELI:0:8} ${LABELI:8:2} -6 hours" `
fi

# Numero de processadores que foram utilizados pelo BAM
if [ -z ${MPITasks} ];then
   MPITasks=72
fi

# TimeStep da Previsão
if [ -z ${DHFCT} ];then
   DHFCT=6
fi

#
# SETTING THE APPROPRIATED ENVIRONMENT
#

case ${TRC} in
     21) IMax=64; JMax=32;;
     31) IMax=96; JMax=48;;
     42) IMax=128; JMax=64;;
     62) IMax=192; JMax=96;;
    106) IMax=320; JMax=160;;
    126) IMax=384; JMax=192;;
    133) IMax=400; JMax=200;;
    159) IMax=480; JMax=240;;
    170) IMax=512; JMax=256;;
    213) IMax=640; JMax=320;;
    254) IMax=768; JMax=384;;
    299) IMax=900; JMax=450;;
    319) IMax=960; JMax=480;;
    341) IMax=1024; JMax=512;;
    382) IMax=1152; JMax=576;;
    511) IMax=1536; JMax=768;;
    533) IMax=1600; JMax=800;;
    666) IMax=2000; JMax=1000;;
    863) IMax=2592; JMax=1296;;
   1279) IMax=3840; JMax=1920;;
   1332) IMax=4000; JMax=2000;;
   *)echo -e "\033[32;1m Truncamento desconhecido ${TRC} \033[m"
esac

MRES=$(printf "TQ%04dL%03d" ${TRC} ${LV})
LABL=$(printf "G%05d" ${JMax})

# PBS
walltime=00:45:00
queue=pesq
queue_name="BLSDAS${LABELI:6:10}"

# Configurando os Diretorios
# Considera-se que o sistema sera rodado no scratchin por isso tudo
# sera configurado a partir do home do usuario no scratch1 (${SUBMIT_HOME})
# Diretorio temporario para a rodada do SMG
if [ ! -z ${subt_blsdas_run} ];then
   RunBLSDAS=${subt_blsdas_run}
   if [ -e ${RunBLSDAS} ]; then
      rm -fr ${RunBLSDAS}/*
   else
      mkdir -p ${RunBLSDAS}
   fi
   cd ${RunBLSDAS}
else
   echo -e "\033[34;1m BLSDAS\033[m \033[31;1menvironment is not set ! \033[m"
   echo -e "\033[31;1m Please verify it in config_smg.sh \033[m"
   echo -  " "
   exit -1
fi

# Verifica se o Diretorio de saida do BLSDAS existe
if [ ! -e ${work_blsdas_dataout} ];then
   mkdir -p ${work_blsdas_dataout}
   ln -s ${work_blsdas_dataout} ${subt_blsdas_dataout}
fi
if [ ! -L ${subt_blsdas_dataout} ];then
   rm -fr ${subt_blsdas_dataout}
   ln -s ${work_blsdas_dataout} ${subt_blsdas_dataout}
fi

# Copiando executavel do blsdas para o diretorio de rodada
FILE=${home_blsdas_bin}/blsdas.x
if [ -e ${FILE} ]
then
  ExecBLSDAS=$(basename ${FILE})
  cp -pfr ${FILE} ${RunBLSDAS}/${ExecBLSDAS}
  chmod +x  ${RunBLSDAS}/${ExecBLSDAS}
else
  echo -e "\033[34;1m[\033[m\033[31;1m Falhou \033[m\033[34;1m]\033[m"
  echo -e "\033[31;1m !!! Arquivo Nao Encontrado !!! \033[m"
  echo -e "\033[31;1m ${FILE} \033[m"
  rm -fr ${RunBLSDAS}
  exit 1
fi

#---------------------------------------------------------------------#
# Configurando parametros para a rodada
#
# Modelo BAM
ModelDirOut=${work_model_bam_dataout}/${MRES}/${LABELFGS}
AtmFG=${ModelDirOut}/GFCT${PREFIX}%fy4%fm2%fd2%fh2%y4%m2%d2%h2F.%e.${MRES}

#ModelDirRst=${work_model_bam_dataout}/${MRES}/${LABELI}/RST
#MParRst=${ModelDirRst}/GFCT${PREFIX}%y4%m2%d2%h2%y4%m2%d2%h2F.unf.${MRES}.readrstP%e
#SfcRstI=${ModelDirRst}/GFCT${PREFIX}%y4%m2%d2%h2%y4%m2%d2%h2F.unf.${MRES}.sibprgP%e
#SfcRstO=${ModelDirRst}/GFCT${PREFIX}%y4%m2%d2%h2%y4%m2%d2%h2F.unf.${MRES}.sibprgP%e
#for file in ${AtmFG} ${MParRst} ${SfcRstI};do
for file in ${AtmFG};do
#   nfile=$(echo ${file} | sed -e 's/%e/000/g' -e "s/%y4%m2%d2%h2/${LABELI}/g" -e "s/%fy4%fm2%fd2%fh2/${LABELFGS}/g")
   nfile=$(echo ${file} | sed -e 's/%e/fct/g' -e "s/%y4%m2%d2%h2/${LABELI}/g" -e "s/%fy4%fm2%fd2%fh2/${LABELFGS}/g")
   if [ ! -e ${nfile} ];then
      echo -e "\033[31;1mFile not found !\033[m"
      echo -e "\033[33;1m${nfile}\033[m"
#      exit -1
   fi
done

#Observacoes
ObsFile=${work_gsi_datain_obs}/gdas1.t%h2z.prepbufr.nr.%y4%m2%d2

#BLSDAS
AtmANL=${subt_blsdas_dataout}/%y4%m2%d2%h2/BLSDAS%y4%m2%d2%h2.atm.${LABL}
SfcANL=${subt_blsdas_dataout}/%y4%m2%d2%h2/BLSDAS%y4%m2%d2%h2.sfc.${LABL}
Stats=${subt_blsdas_dataout}/blsdas.stats

DirOut=${subt_blsdas_dataout}/${LABELI}
if [ ! -e ${DirOut} ];then
   mkdir -p ${DirOut}
fi
#
#---------------------------------------------------------------------#

# Alterando o blsdas.conf
sed -e "s;#LABELI#;${LABELI};g" \
    -e "s;#DHFCT#;${DHFCT};g" \
    -e "s;#AtmFG#;${AtmFG};g" \
    -e "s;#ObsFile#;${ObsFile};g" \
    -e "s;#AtmANL#;${AtmANL};g" \
    -e "s;#SfcANL#;${SfcANL};g" \
    -e "s;#Stats#;${Stats};g" \
    -e "s;#MPITasks#;${MPITasks};g" \
    -e "s;#MParRst#;${MParRst};g" \
    -e "s;#SfcRstI#;${SfcRstI};g" \
    -e "s;#SfcRstO#;${SfcRstO};g" \
    ${LOCALDIR}/blsdas_scripts/blsdas.conf.template > ${RunBLSDAS}/blsdas.conf

# Criando qsub

tmstp=$(date +'%s')

cat << EOF > ${RunBLSDAS}/blsdas.qsub
#!/bin/bash
#PBS -o ${HSTMAQ}:${RunBLSDAS}/Out.blsdas.${LABELI}.${tmstp}.out
#PBS -j oe
#PBS -l walltime=${walltime}
#PBS -l mppwidth=24
#PBS -l mppnppn=24
#PBS -l mppdepth=1
#PBS -V
#PBS -S /bin/bash
#PBS -N ${queue_name}
#PBS -q ${queue}
#PBS -A CPTEC
###PBS -W block=true

export ATP_ENABLED=1
ulimit -s unlimited
ulimit -c unlimited


cd ${RunBLSDAS}

SECONDS=0

/usr/bin/time -v ${RunBLSDAS}/$(basename ${ExecBLSDAS}) > log.log 2>&1

duration=\${SECONDS}
echo "\$((\$duration / 60)) minutes and \$((\$duration % 60)) seconds elapsed."

#> ${RunBLSDAS}/monitor.blsdas

EOF

# Executando o modelo
cd ${RunBLSDAS}
qsub -W block=true blsdas.qsub

#until [ -e ${RunBLSDAS}/monitor.blsdas ]; do sleep 1s; done
#rm -fr ${BAMRUN}/monitor.blsdas

# Copy file to model datain
File=$(echo $SfcANL | sed -e "s/%y4%m2%d2%h2/${LABELI}/g")
cp -pvfr ${File} ${subt_model_bam_datain}

exit 0
#EOC
#-----------------------------------------------------------------------------#
