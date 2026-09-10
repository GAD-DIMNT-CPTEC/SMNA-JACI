#! /bin/bash

#set -e
# Descomente a linha abaixo para debugar
#set -o xtrace

# Carregando as variaveis do sistema
dir_now=`pwd`

cd $SMG_ROOT

SCRIPT_PATH="$(realpath "${BASH_SOURCE[0]}")"
RootDir="$(dirname "$SCRIPT_PATH")"
export SMG_ROOT=${RootDir}
source ${SMG_ROOT}/../../config_smg.ksh vars_export
source ${SMG_ROOT}/run/smg_functions.sh     ## has inctime function wrote in bash script

cd $dir_now

# Lendo parametros de entrada
if [ -z "${1}" ]
then
  echo "LABELANL is not set"
  exit 3
else
  export LABELANL=${1}
export LABELFGS=`${inctime} ${LABELANL} -6h %y4%m2%d2%h2`
  ###   export LABELFGS=`date -u +%Y%m%d%H -d "${LABELANL:0:8} ${LABELANL:8:2} -6 hours" `
export LABELFCT=`${inctime} ${LABELANL} +6h %y4%m2%d2%h2`
  ###   export LABELFCT=`date -u +%Y%m%d%H -d "${LABELANL:0:8} ${LABELANL:8:2} +6 hours" `
fi
if [ -z "${2}" ]
then
  echo "PREFIX is not set"
  exit 3
else
  export PREFIX=${2}
fi
if [ -z "${3}" ]
then
  echo "TRC is not set"
  exit 3
else
  export TRC=${3}
fi
if [ -z "${4}" ]
then
  echo "NLV is not set"
  exit 3
else
  export NLV=${4}
fi
if [ -z "${5}" ]
then
  echo "NPROC is not set"
  exit 3
else
  export NPROC=${5}
  if [ ${NPROC} -lt 24 ]
  then
    echo "NPROC less then 24"
    echo "setting NPROC to 24"
    export NPROC=24
  fi
fi
if [ -z "${6}" ]
then
  echo "TRC Model is not set"
  TRCB=${TRC}
else
  export TRCB=${6}
fi

export MRES=`echo ${TRC} ${NLV} | awk '{printf("TQ%4.4dL%3.3d\n",$1,$2)}'`
export MRESB=`echo ${TRCB} ${NLV} | awk '{printf("TQ%4.4dL%3.3d\n",$1,$2)}'`

#-------------------------------------------------------------------#
# Defininfo numero de pontos de grade
#

case ${TRC} in
   62) IMAX=192; JMAX=96;;
   126)IMAX=384; JMAX=192;;
   213)IMAX=640; JMAX=320;;
   254)IMAX=768; JMAX=384;;
   299)IMAX=900; JMAX=450;;
   666)IMAX=2000; JMAX=1000;;
   *)echo "Truncamento desconhecido ${MRES}"
esac
KMAX=${NLV}
JMAX=$((JMAX+2))

#
#-------------------------------------------------------------------#

#
# Defining byte order
#

BYTE_ORDER=Big_Endian


# Configurando os Diretorios
# Considera-se que o sistema sera rodado no scratchin por isso tudo
# sera configurado a partir do home do usuario no scratch1 (${SUBMIT_HOME})
# Diretorio temporario para a rodada do SMG
RunGSI=${subt_run_gsi}
if [ -e ${RunGSI} ]; then
   rm -fr ${RunGSI}/*
else
   mkdir -p ${RunGSI}
fi
cd ${RunGSI}

# Diretorio de saida do GSI
export OutGSI=${work_gsi_dataout}

# Preparing the environment for GSI run
export MPICH_UNEX_BUFFER_SIZE=100000000
export MPICH_MAX_SHORT_MSG_SIZE=4096
export MPICH_PTL_UNEX_EVENTS=50000
export MPICH_PTL_OTHER_EVENTS=2496

# Copiando executavel do gsi para o diretorio de rodada
FILE=${home_cptec}/bin/gsi.x
if [ -e ${FILE} ]
then
  ExecGSI=$(basename ${FILE})
  cp -pvfr ${FILE} ${RunGSI}/${ExecGSI}
  chmod +x  ${RunGSI}/${ExecGSI}
else
  echo -e "\033[34;1m[\033[m\033[31;1m Falhou \033[m\033[34;1m]\033[m"
  echo -e "\033[31;1m !!! Arquivo Nao Encontrado !!! \033[m"
  echo -e "\033[31;1m ${FILE} \033[m"
  rm -fr ${RunGSI}
  exit 1
fi

if [ ! -e ${RunGSI}/bin ]; then mkdir -p ${RunGSI}/bin; fi


# Configurando parametros para a rodada
mpicmd="aprun -n ${NPROC}"
yyyymmdd=${LABELANL:0:8}

#Acrescimo de variaveis para a geracao do satbias_angle
#yyyy2=`echo ${LABELAN}|cut -c 1-4`
yy1=${LABELANL:0:4}
mm1=`echo ${LABELANL}|cut -c 5-6`
dd1=`echo ${LABELANL}|cut -c 7-8`

yymmdd=${LABELANL:2:${#LABELANL}}
hm3=`${inctime} ${LABELANL} -3h %h2` #15
## hm3=`date -u +%H -d "${LABELANL:0:8} ${LABELANL:8:2} -3 hours" `
hh0=`${inctime} ${LABELANL} +0h %h2` #18
## hh0=${LABELANL:8:2} 
hp3=`${inctime} ${LABELANL} +3h %h2` #21
## hp3=`date -u +%H -d "${LABELANL:0:8} ${LABELANL:8:2} +3 hours" `
expid=cptec
bkg=bkg

save=${OutGSI}/${LABELANL}
if [ ! -d ${save} ]; then mkdir -p ${save}; fi

# Escolha das observacoes
ObsDir=${work_gsi_datain_obs}

cat << EOF > ${RunGSI}/nmlobs.gsi
#
#Simple namelist for GSI observations selection 
#

FLAG                  FILENAME                                     ALIAS

0  ${ObsDir}/gmaoairs.${yyyymmdd}.t${hh0}z.bufr               airsbufr
0  ${home_gsi}/rc/airs_bufr.table                             airs_bufr.table
0  ${ObsDir}/madeup.tcvitals                                  tcvitl
0  ${ObsDir}/gdas1.t${hh0}z.airsev.tm00.bufr_d.${yyyymmdd}    airsbufr
1  ${ObsDir}/gdas1.t${hh0}z.1bamua.tm00.bufr_d.${yyyymmdd}    amsuabufr
0  ${ObsDir}/gdas1.t${hh0}z.1bamub.tm00.bufr_d.${yyyymmdd}    amsubbufr
0  ${ObsDir}/gdas1.t${hh0}z.1bhrs3.tm00.bufr_d.${yyyymmdd}    hirs3bufr
0  ${ObsDir}/gdas1.t${hh0}z.1bhrs4.tm00.bufr_d.${yyyymmdd}    hirs4bufr
0  ${ObsDir}/gdas1.t${hh0}z.1bmhs.tm00.bufr_d.${yyyymmdd}     mhsbufr
0  ${ObsDir}/gdas1.t${hh0}z.amsre.tm00.bufr_d.${yyyymmdd}     amsrebufr
0  ${ObsDir}/gdas1.t${hh0}z.goesfv.tm00.bufr_d.${yyyymmdd}    gsnd1bufr
1  ${ObsDir}/gdas1.t${hh0}z.gpsro.tm00.bufr_d.${yyyymmdd}     gpsrobufr
0  ${ObsDir}/gdas1.t${hh0}z.mtiasi.tm00.bufr_d.${yyyymmdd}    iasibufr
0  ${ObsDir}/gdas1.t${hh0}z.ncepssmit.tm00.bufr_d.${yyyymmdd} ssmitbufr
0  ${ObsDir}/gdas1.t${hh0}z.osbuv.tm00.bufr_d.${yyyymmdd}     sbuvbufr
0  ${ObsDir}/gdas1.t${hh0}z.sptrmm.tm00.bufr_d.${yyyymmdd}    tmirrbufr
0  ${ObsDir}/gdas1.t${hh0}z.spssmi.tm00.bufr_d.${yyyymmdd}    ssmirrbufr
1  ${ObsDir}/gdas1.t${hh0}z.prepbufr.nr.${yyyymmdd}           prepbufr
0  ${ObsDir}/gdas1.t${hh0}z.osbuv8.tm00.bufr_d.${yyyymmdd}    sbuvbufr
0  ${ObsDir}/gdas1.t${hh0}z.avcs18.tm00.bufr_d.${yyyymmdd}    avh18bufr
0  ${ObsDir}/gdas1.t${hh0}z.avcsmeta.tm00.bufr_d.${yyyymmdd}  avhmabufr
0  ${ObsDir}/gdas1.t${hh0}z.esamua.tm00.bufr_d.${yyyymmdd}    amsuabufrears
0  ${ObsDir}/gdas1.t${hh0}z.esamub.tm00.bufr_d.${yyyymmdd}    amsubbufrears
0  ${ObsDir}/gdas1.t${hh0}z.eshrs3.tm00.bufr_d.${yyyymmdd}    hirs3bufrears
0  ${ObsDir}/gdas1.t${hh0}z.geoimr.tm00.bufr_d.${yyyymmdd}    geoimr
1  ${ObsDir}/gdas1.t${hh0}z.satwnd.tm00.bufr_d.${yyyymmdd}    satwnd

1  ${ObsDir}/gdas.t${hh0}z.atms.tm00.bufr_d                   atmsbufr
EOF

# Informa os tipos de observacoes selecionada e cria os links simbolicos
#cat ${subt_run_gsi}/nmlobs.gsi | grep -E '^ *1' > ${subt_run_gsi}/gsiobs.avail
#nobs=`cat ${subt_run_gsi}/gsiobs.avail | wc -l`
#while read tipo
#do
#  nome=`echo ${tipo} | awk -F " " '{print $3}'`
#  orig=`echo ${tipo} | awk -F " " '{print $2}'`
#  dest=`echo ${tipo} | awk -F " " '{print $3}'`
#  ln -sf ${orig} ${dest}
#  echo "Utilizando tipo de observacao: ${nome}"
#done < ${subt_run_gsi}/gsiobs.avail
#echo "Numero de tipos de observacoes selecionada: ${nobs}"
let cont_nobs=0

# Conta os tipos de observacoes utilizadas,imprime no stdout e faz o link para o GSI assimilar
echo ""
echo -e "\033[32;2m > Selecao da observacoes...\033[m"
cat ${subt_run_gsi}/nmlobs.gsi | grep -E '^ *1' > ${subt_run_gsi}/gsiobs.avail
nobs=`cat ${subt_run_gsi}/gsiobs.avail | wc -l`
echo -e "\033[32;2m > ${nobs} tipos de observacoes foram selecionadas:\033[m"
while read tipo
do
  nome=`echo ${tipo} | awk -F " " '{print $3}'`
  orig=`echo ${tipo} | awk -F " " '{print $2}'`
  dest=`echo ${tipo} | awk -F " " '{print $3}'`
#  ln -sf ${orig} ${dest}
  cp -pfr ${orig} ${dest}
  ls `find ${nome} -type l -printf "%p %l\n" | awk -F " " '{print $2}'` > /dev/null 2>&1
  if [ `echo $?` -ne "0" ]
  then
    echo -e "\033[31;1m > Falha ao criar link para observacao: ${nome}\033[m"
  else
    echo -e "\033[34;1m * ${nome}\033[m"
    cont_nobs=$((${cont_nobs}+1))
  fi
done < ${subt_run_gsi}/gsiobs.avail
echo -e "\033[32;2m > Tipos de observacoes que serao utilizadas: ${cont_nobs}\033[m"
echo ""


#-------------------------------------------------------------------------------#
# Linkando arquivos para background
#

for inc in $(seq -3 3 3); do
   TIME=$(printf "%02g" $((inc+6)))
   LABEL=$(${inctime} ${LABELANL} ${inc}h %y4%m2%d2%h2)
   ### LABEL=`date -u +%Y%m%d%H -d "${LABELANL:0:8} ${LABELANL:8:2} +${inc} hours" `
   FFCT=GFCT${PREFIX}${LABELFGS}${LABEL}F.fct.${MRESB}
   FDIR=GFCT${PREFIX}${LABELFGS}${LABEL}F.dir.${MRESB}

   cp -pfr ${work_model_bam_dataout}/${MRESB}/${LABELFGS}/${FFCT} BAM.fct.${TIME}
   cp -pfr ${work_model_bam_dataout}/${MRESB}/${LABELFGS}/${FDIR} BAM.dir.${TIME}
done

#
#-------------------------------------------------------------------------------#
#  link Fixed fields to working directory
#

#cp -pfr ${public_fix}/global_anavinfo.l${NLV}.txt        anavinfo
cp -pfr ${public_gsi}/Berror/SMG/gsir4.berror_stats.gcv.BAM.${MRES} berror_stats
#cp -pfr ${public_fix}/global_satangbias.txt              satbias_angle
cp -pfr ${public_fix}/atms_beamwidth.txt                 atms_beamwidth.txt
#cp -pfr ${public_fix}/global_scaninfo.txt                scaninfo
cp -pfr ${public_fix}/global_ozinfo.txt                  ozinfo
cp -pfr ${public_fix}/global_pcpinfo.txt                 pcpinfo
cp -pfr ${public_fix}/prepobs_errtable.global            errtable
cp -pfr ${public_fix}/prepobs_prep.bufrtable             prepobs_prep.bufrtable
cp -pfr ${public_fix}/bufrtab.012                        bftab_sstphr

# Copiando arquivos de configuracao do usuário
cp -pfr ${home_gsi_fix}/global_anavinfo.l${NLV}.txt       anavinfo
cp -pfr ${home_gsi_fix}/global_satinfo.txt                 satinfo
cp -pfr ${home_gsi_fix}/global_convinfo.txt                convinfo

file=${OutGSI}/${LABELFGS}/satbias_out
fileOut=satbias_in
sample=${public_fix}/sample.satbias

if [ -e ${file} ];then
   cp -pfr ${file} ${fileOut}
else
   echo -e " "
   echo -e "\e[31;1m Atenção:\e[m\e[36;1m O arquivo\e[m\e[37;1m $(basename ${file})\e[m\e[36;1m não existe:\e[m"
   echo -e "\e[31;1m         \e[36;1m criando\e[m\e[37;1m ${fileOut}\e[m\e[36;1m a partir do arquivo \e[m\e[37;1m ${sample}\e[m"
   echo -e " "
   echo -e "\e[35;1m         Este procedimento deve ocorrer somente no\e[m\e[31;1m 1° ciclo de assimilação\e[m\e[35;1m !\e[m"
   echo -e "\e[35;1m         Caso este não seja o 1° ciclo de assimilação verifique o motivo.\e[m"
   echo -e " "

  cp -pfr ${sample} ${fileOut}

fi

#Testa se os arquivos de correção de bias referente ao ângulo foram gerados no ciclo anterior, caso contrário copia arquivo fixo

file=${OutGSI}/${LABELFGS}/satbias_ang.out.${LABELFGS}
fileOut=satbias_angle
sample=${public_fix}/global_satangbias.txt

if [ -e ${file} ];then
   cp -pfr ${file} ${fileOut}
else
   echo -e " "
   echo -e "\e[31;1m Atenção:\e[m\e[36;1m O arquivo\e[m\e[37;1m $(basename ${file})\e[m\e[36;1m não existe:\e[m"
   echo -e "\e[31;1m         \e[36;1m criando\e[m\e[37;1m ${fileOut}\e[m\e[36;1m a partir do arquivo \e[m\e[37;1m ${sample}\e[m"
   echo -e " "
   echo -e "\e[35;1m         Este procedimento deve ocorrer somente no\e[m\e[31;1m 1° ciclo de assimilação\e[m\e[35;1m !\e[m"
   echo -e "\e[35;1m         Caso este não seja o 1° ciclo de assimilação verifique o motivo.\e[m"
   echo -e " "

   cp -pfr ${sample} ${fileOut}

fi


#

#
#-------------------------------------------------------------------------------#
# link CRTM Spectral and Transmittance coefficients
#

#ln -sf ${public_crtm}/${BYTE_ORDER}/Nalli.IRwater.EmisCoeff.bin    Nalli.IRwater.EmisCoeff.bin
#ln -sf ${public_crtm}/${BYTE_ORDER}/NPOESS.IRice.EmisCoeff.bin     NPOESS.IRice.EmisCoeff.bin
#ln -sf ${public_crtm}/${BYTE_ORDER}/NPOESS.IRland.EmisCoeff.bin    NPOESS.IRland.EmisCoeff.bin
#ln -sf ${public_crtm}/${BYTE_ORDER}/NPOESS.IRsnow.EmisCoeff.bin    NPOESS.IRsnow.EmisCoeff.bin
#ln -sf ${public_crtm}/${BYTE_ORDER}/NPOESS.VISice.EmisCoeff.bin    NPOESS.VISice.EmisCoeff.bin
#ln -sf ${public_crtm}/${BYTE_ORDER}/NPOESS.VISland.EmisCoeff.bin   NPOESS.VISland.EmisCoeff.bin
#ln -sf ${public_crtm}/${BYTE_ORDER}/NPOESS.VISsnow.EmisCoeff.bin   NPOESS.VISsnow.EmisCoeff.bin
#ln -sf ${public_crtm}/${BYTE_ORDER}/NPOESS.VISwater.EmisCoeff.bin  NPOESS.VISwater.EmisCoeff.bin
#ln -sf ${public_crtm}/${BYTE_ORDER}/FASTEM5.MWwater.EmisCoeff.bin  FASTEM5.MWwater.EmisCoeff.bin
#ln -sf ${public_crtm}/${BYTE_ORDER}/AerosolCoeff.bin               AerosolCoeff.bin
#ln -sf ${public_crtm}/${BYTE_ORDER}/CloudCoeff.bin                 CloudCoeff.bin

cp -v ${plus_crtm}/${BYTE_ORDER}/Nalli.IRwater.EmisCoeff.bin    Nalli.IRwater.EmisCoeff.bin
cp -v ${plus_crtm}/${BYTE_ORDER}/NPOESS.IRice.EmisCoeff.bin     NPOESS.IRice.EmisCoeff.bin
cp -v ${plus_crtm}/${BYTE_ORDER}/NPOESS.IRland.EmisCoeff.bin    NPOESS.IRland.EmisCoeff.bin
cp -v ${plus_crtm}/${BYTE_ORDER}/NPOESS.IRsnow.EmisCoeff.bin    NPOESS.IRsnow.EmisCoeff.bin
cp -v ${plus_crtm}/${BYTE_ORDER}/NPOESS.VISice.EmisCoeff.bin    NPOESS.VISice.EmisCoeff.bin
cp -v ${plus_crtm}/${BYTE_ORDER}/NPOESS.VISland.EmisCoeff.bin   NPOESS.VISland.EmisCoeff.bin
cp -v ${plus_crtm}/${BYTE_ORDER}/NPOESS.VISsnow.EmisCoeff.bin   NPOESS.VISsnow.EmisCoeff.bin
cp -v ${plus_crtm}/${BYTE_ORDER}/NPOESS.VISwater.EmisCoeff.bin  NPOESS.VISwater.EmisCoeff.bin
cp -v ${plus_crtm}/${BYTE_ORDER}/FASTEM5.MWwater.EmisCoeff.bin  FASTEM5.MWwater.EmisCoeff.bin
cp -v ${plus_crtm}/${BYTE_ORDER}/AerosolCoeff.bin               AerosolCoeff.bin
cp -v ${plus_crtm}/${BYTE_ORDER}/CloudCoeff.bin                 CloudCoeff.bin

#
#-------------------------------------------------------------------------------#
# Copy CRTM coefficient files based on entries in satinfo file
#

for file in `awk '{if($1!~"!"){print $1}}' ./satinfo | sort | uniq` ;do
#   ln -s ${public_crtm}/${BYTE_ORDER}/${file}.SpcCoeff.bin ./
#   ln -s ${public_crtm}/${BYTE_ORDER}/${file}.TauCoeff.bin ./
   cp -v ${plus_crtm}/${BYTE_ORDER}/${file}.SpcCoeff.bin ./
   cp -v ${plus_crtm}/${BYTE_ORDER}/${file}.TauCoeff.bin ./
   # Carlos (02/07/2025) - estou colocando aqui os dados para o cálculo da profundidade óptica do algorítmo ODPS
   # por ser mais adequado para dados do microondas (no caso do ATMS, os arquivos atms_*.TauCoeff.bin - com excessão do npp
   # só esão disponíveis pelo ODPS, pelo pacote crtm-2.4.0_emc.1 disponível no GitHub)
   cp -v ${plus_crtm}/TauCoeff/ODPS/${BYTE_ORDER}/${file}.TauCoeff.bin ${runDir}
done

#
#-------------------------------------------------------------------------------#

File=${save}/${exp}_${label}.log
if [ -e ${File} ]; then rm ${File}; fi

#-------------------------------------------------------------------#
# Configure gsiparam
#

sed "s/#CENTER#/cptec/g" ${home_gsi_fix}/gsiparm.anl > ${subt_run_gsi}/gsiparm.anl
sed -i -e "s/#IMAX#/${IMAX}/g" \
       -e "s/#JMAX#/${JMAX}/g" \
       -e "s/#KMAX#/${KMAX}/g" \
       -e "s/#TRUNC#/${TRC}/g" \
       -e "s/#TRUNCB#/${TRCB}/g" \
       -e "s/#NLVL#/$((${NLV}-1))/g" \
       -e "s/#NLV#/${NLV}/g" \
       -e "s/#LABELANL#/${LABELANL}/" ${subt_run_gsi}/gsiparm.anl
#-------------------------------------------------------------------#
TIME=`date '+%H:%M:%S'`

cat<< EOF > ${RunGSI}/qsub_gdad.qsb
#!/bin/bash
#PBS -o ${save}/gdad_anl.${LABELANL}_${TIME}.out
#PBS -e ${save}/gdad_anl.${LABELANL}_${TIME}.err
#PBS -l walltime=00:45:00
#PBS -l mppwidth=${NPROC}
#PBS -l mppnppn=12
#PBS -V
#PBS -S /bin/bash
#PBS -N gdad_anl
#PBS -q pesq
#PBS -j oe
#PBS -A CPTEC
###PBS -W block=true


# LGGG this section still needs investigation
#setenv NST 0
#setenv BZBERROR 1
#setenv BZQERROR 0
#setenv FDDA 0
#limit stacksize unlimited

export ATP_ENABLED=1
ulimit -c unlimited
ulimit -s unlimited

cd ${RunGSI}

time ${mpicmd} ${RunGSI}/${ExecGSI} > ${save}/stdout_${LABELANL}_$(date +"%Y.%m.%d:%H.%M.%S").log

#touch ${RunGSI}/monitor.t

EOF

chmod 755 ${RunGSI}/qsub_gdad.qsb
qsub -W block=true ${RunGSI}/qsub_gdad.qsb

#until [ -e ${RunGSI}/monitor.t ]; do sleep 1s; done
#rm -fr ${RunGSI}/monitor.t

AnlOut=${RunGSI}/BAM.anl

if [ -e $AnlOut ]; then

 echo -e "\033[32;1m #---------------------------------------#\033[m"
 echo -e "\033[32;1m #                                       #\033[m"
 echo -e "\033[32;2m # Arquivo de analise gerado com sucesso #\033[m"
 echo -e "\033[32;1m #                                       #\033[m"
 echo -e "\033[32;1m #---------------------------------------#\033[m"

 cp -pvfr ${AnlOut} ${save}/GANL${PREFIX}${LABELANL}S.unf.${MRESB}
 cp -pvfr  ${AnlOut} ${subt_pre_bam_datain}/GANL${PREFIX}${LABELANL}S.unf.${MRESB}

else

 echo -e "\033[31;1m #--------------------------------------#\033[m"
 echo -e "\033[31;1m #                                      #\033[m"
 echo -e "\033[31;1m #    Falha durante o processo de AD    # \033[m"
 echo -e "\033[31;1m #                                      #\033[m"
 echo -e "\033[31;1m #--------------------------------------#\033[m"
 exit -1

fi

# Loop over first and last outer loops to generate innovation
# diagnostic files for indicated observation types (groups)
#
# NOTE:  Since we set miter=2 in GSI namelist SETUP, outer
#        loop 03 will contain innovations with respect to
#        the analysis.  Creation of o-a innovation files
#        is triggered by write_diag(3)=.true.  The setting
#        write_diag(1)=.true. turns on creation of o-g
#        innovation files.
#

loops="01 03"
loops=$(ls *conv_* | awk -F"_" '{print $2}' | sort -u | xargs)
for loop in $loops; do

   case $loop in
     01) string=ges;;
     03) string=anl;;
      *) string=$loop;;
   esac

#  Collect diagnostic files for obs types (groups) below
   listall="conv hirs2_n14 msu_n14 sndr_g08 sndr_g11 sndr_g12 sndr_g13 sndr_g08_prep sndr_g11_prep sndr_g12_prep sndr_g13_prep
            sndrd1_g11 sndrd2_g11 sndrd3_g11 sndrd4_g11 sndrd1_g12 sndrd2_g12 sndrd3_g12 sndrd4_g12 sndrd1_g13 sndrd2_g13
            sndrd3_g13 sndrd4_g13 sndrd1_g14 sndrd2_g14 sndrd3_g14 sndrd4_g14 sndrd1_g15 sndrd2_g15 sndrd3_g15 sndrd4_g15
            hirs3_n15 hirs3_n16 hirs3_n17 amsua_n15 amsua_n16 amsua_n17 amsub_n15 amsub_n16 amsub_n17 hsb_aqua airs_aqua
            amsua_aqua imgr_g08 imgr_g11 imgr_g12 ssmi_f13 ssmi_f14 imgr_g14 imgr_g15 ssmi_f15 hirs4_n18 hirs4_metop-a amsua_n18
            amsua_metop-a mhs_n18 mhs_metop-a amsre_low_aqua amsre_mid_aqua amsre_hig_aqua ssmis_las_f16 ssmis_uas_f16 ssmis_img_f16
            ssmis_env_f16 ssmis_las_f17 ssmis_uas_f17 ssmis_img_f17 ssmis_env_f17 ssmis_las_f18 ssmis_uas_f18 ssmis_img_f18 ssmis_env_f18
            ssmis_las_f19 ssmis_uas_f19 ssmis_img_f19 ssmis_env_f19 ssmis_las_f20 ssmis_uas_f20 ssmis_img_f20 ssmis_env_f20 iasi_metop-a 
            hirs4_n19 amsua_n19 mhs_n19 seviri_m08 seviri_m09 seviri_m10 cris_npp atms_npp hirs4_metop-b amsua_metop-b mhs_metop-b iasi_metop-b"
   for type in $listall; do
      count=0
      if [[ -f pe0000.${type}_${loop} ]]; then
         count=`ls pe*${type}_${loop}* | wc -l`
      fi
      if [[ $count -gt 0 ]]; then
         cat pe*${type}_${loop}* > diag_${type}_${string}.${LABELANL}
      fi
   done
done


# Leva os arquivos de saida que contem a estatistica do SMG
cp -pfr diag_* ${save}/
cp -pfr fort.2* ${save}/
cp -pfr gsiparm.anl ${save}/
cp -pvfr satbias_out ${save}/
cp -pvfr satbias_in ${save}/
cp -pvfr satbias_angle ${save}/


mkdir -p ${save}/diag
cp -pfr pe* ${save}/diag

# Concatenando os arquivos pe* das observações de radiância de satélites usando o script gsidiags
cd ${save}/diag
ln -s ../gsiparm.anl .
echo "      Run gsidiags ${yyyymmdd} ${hh0}0000 cptec set"
${scripts_smg}/gsi_scripts/gsidiags ${yyyymmdd} ${hh0}0000 cptec set > gsidiags.log 2>&1 
cd ..

# Remover os arquivos de BKG do diretorio /scratchin/grupos/assim_dados/home/gdad/GSI/Ana/Bkg/cptec e remover as Observacoes 
#rm -fr ${RunGSI}/*

##################################################################################################################
# Rodando um utilitário para gerar a correção de bias de satélite
# baseada no angulo em um modo cíclico
#

# Create the ram work directory and cd into it
workdirSatAng=${RunGSI}/angupdate

echo -e "\033[32;2m > Criando a pasta\033[m"
echo -e "\033[32;2m > Para gerar correção do ângulo variando no tempo !!!\033[m"
echo ""

if [ -d "${workdirSatAng}" ]; then
  rm -rf ${workdirSatAng}
fi
mkdir -p ${workdirSatAng}
cd ${workdirSatAng}

# Loop over first and last outer loops to generate innovation
# diagnostic files for indicated observation types (groups)
#
ls -l ${RunGSI}/diag_* > listpe
#  Collect diagnostic files for obs types (groups) below
  listall="hirs2_n14 msu_n14 sndr_g08 sndr_g11 sndr_g12 sndr_g13 sndr_g08_prep sndr_g11_prep 
           sndr_g12_prep sndr_g13_prep sndrd1_g11 sndrd2_g11 sndrd3_g11 sndrd4_g11 sndrd1_g12
           sndrd2_g12 sndrd3_g12 sndrd4_g12 sndrd1_g13 sndrd2_g13 sndrd3_g13 sndrd4_g13 sndrd1_g14
           sndrd2_g14 sndrd3_g14 sndrd4_g14 sndrd1_g15 sndrd2_g15 sndrd3_g15 sndrd4_g15 hirs3_n15
           hirs3_n16 hirs3_n17 amsua_n15 amsua_n16 amsua_n17 amsub_n15 amsub_n16 amsub_n17 hsb_aqua
           airs_aqua amsua_aqua imgr_g08 imgr_g11 imgr_g12 ssmi_f13 ssmi_f14 imgr_g14 imgr_g15 ssmi_f15
           hirs4_n18 hirs4_metop-a amsua_n18 amsua_metop-a mhs_n18 mhs_metop-a amsre_low_aqua amsre_mid_aqua
           amsre_hig_aqua ssmis_las_f16 ssmis_uas_f16 ssmis_img_f16 ssmis_env_f16 ssmis_las_f17 ssmis_uas_f17
           ssmis_img_f17 ssmis_env_f17 ssmis_las_f18 ssmis_uas_f18 ssmis_img_f18 ssmis_env_f18 ssmis_las_f19
           ssmis_uas_f19 ssmis_img_f19 ssmis_env_f19 ssmis_las_f20 ssmis_uas_f20 ssmis_img_f20 ssmis_env_f20
           iasi_metop-a hirs4_n19 amsua_n19 mhs_n19 seviri_m08 seviri_m09 seviri_m10 cris_npp atms_npp hirs4_metop-b
           amsua_metop-b mhs_metop-b iasi_metop-b"

   for type in $listall; do
      count=`grep diag_${type}_ges.${LABELANL} listpe | wc -l`
      if [[ $count -gt 0 ]]; then
         ln -fs ${RunGSI}/diag_${type}_ges.${LABELANL} ./diag_${type}.${LABELANL}
      fi
   done



# Copy GSI executable, background file,
echo -e "\033[32;2m > Copiando o executável da atualização da correção do bias do ângulo \033[m"

# Copiando executavel do gsi para o diretorio de rodada
ANGEXE=${home_cptec}/bin/global_angupdate
if [ -e ${ANGEXE} ]
then
  Execang=$(basename ${ANGEXE})
  cp ${ANGEXE} ${workdirSatAng}/${Execang}
  chmod +x  ${workdirSatAng}
else
  echo -e "\033[34;1m[\033[m\033[31;1m Falhou \033[m\033[34;1m]\033[m"
  echo -e "\033[31;1m !!! Arquivo Nao Encontrado !!! \033[m"
  echo -e "\033[31;1m ${ANGEXE} \033[m"

  exit 1
fi


cp ${RunGSI}/satbias_angle ./satbias_ang.in

# Set some parameters for use by the GSI executable and to build the namelist
echo "       Build the namelist for gsi_angupdate "

# Build the GSI namelist on-the-fly
iy=$(echo ${LABELANL} |cut -c1-4)
im=$(echo ${LABELANL} |cut -c5-6)
id=$(echo ${LABELANL} |cut -c7-8)
ih=$(echo ${LABELANL} |cut -c9-10)
cat << EOF > global_angupdate.namelist
 &setup
  jpch=2680,nstep=90,nsize=20,wgtang=0.008333333,wgtlap=0.0,
  iuseqc=1,dtmax=1.0,
  iyy1=${iy},imm1=${im},idd1=${id},ihh1=${ih},
  iyy2=${iy},imm2=${im},idd2=${id},ihh2=${ih},
  dth=01,ndat=50
 /
 &obs_input
  dtype(01)='hirs3',     dplat(01)='n17',       dsis(01)='hirs3_n17',
  dtype(02)='hirs4',     dplat(02)='metop-a',   dsis(02)='hirs4_metop-a',
  dtype(03)='goes_img',  dplat(03)='g11',       dsis(03)='imgr_g11',
  dtype(04)='goes_img',  dplat(04)='g12',       dsis(04)='imgr_g12',
  dtype(05)='airs',      dplat(05)='aqua',      dsis(05)='airs281SUBSET_aqua',
  dtype(06)='amsua',     dplat(06)='n15',       dsis(06)='amsua_n15',
  dtype(07)='amsua',     dplat(07)='n18',       dsis(07)='amsua_n18',
  dtype(08)='amsua',     dplat(08)='metop-a',   dsis(08)='amsua_metop-a',
  dtype(09)='amsua',     dplat(09)='aqua',      dsis(09)='amsua_aqua',
  dtype(10)='mhs',       dplat(10)='n18',       dsis(10)='mhs_n18',
  dtype(11)='mhs',       dplat(11)='metop-a',   dsis(11)='mhs_metop-a',
  dtype(12)='ssmi',      dplat(12)='f15',       dsis(12)='ssmi_f15',
  dtype(13)='amsre_low', dplat(13)='aqua',      dsis(13)='amsre_aqua',
  dtype(14)='amsre_mid', dplat(14)='aqua',      dsis(14)='amsre_aqua',
  dtype(15)='amsre_hig', dplat(15)='aqua',      dsis(15)='amsre_aqua',
  dtype(16)='ssmis_las', dplat(16)='f16',       dsis(16)='ssmis_f16',
  dtype(17)='ssmis_uas', dplat(17)='f16',       dsis(17)='ssmis_f16',
  dtype(18)='ssmis_img', dplat(18)='f16',       dsis(18)='ssmis_f16',
  dtype(19)='ssmis_env', dplat(19)='f16',       dsis(19)='ssmis_f16',
  dtype(20)='sndrd1',    dplat(20)='g12',       dsis(20)='sndrD1_g12',
  dtype(21)='sndrd2',    dplat(21)='g12',       dsis(21)='sndrD2_g12',
  dtype(22)='sndrd3',    dplat(22)='g12',       dsis(22)='sndrD3_g12',
  dtype(23)='sndrd4',    dplat(23)='g12',       dsis(23)='sndrD4_g12',
  dtype(24)='sndrd1',    dplat(24)='g11',       dsis(24)='sndrD1_g11',
  dtype(25)='sndrd2',    dplat(25)='g11',       dsis(25)='sndrD2_g11',
  dtype(26)='sndrd3',    dplat(26)='g11',       dsis(26)='sndrD3_g11',
  dtype(27)='sndrd4',    dplat(27)='g11',       dsis(27)='sndrD4_g11',
  dtype(28)='sndrd1',    dplat(28)='g13',       dsis(28)='sndrD1_g13',
  dtype(29)='sndrd2',    dplat(29)='g13',       dsis(29)='sndrD2_g13',
  dtype(30)='sndrd3',    dplat(30)='g13',       dsis(30)='sndrD3_g13',
  dtype(31)='sndrd4',    dplat(31)='g13',       dsis(31)='sndrD4_g13',
  dtype(32)='iasi',      dplat(32)='metop-a',   dsis(32)='iasi616_metop-a',
  dtype(33)='hirs4',     dplat(33)='n19',       dsis(33)='hirs4_n19',
  dtype(34)='amsua',     dplat(34)='n19',       dsis(34)='amsua_n19',
  dtype(35)='mhs',       dplat(35)='n19',       dsis(35)='mhs_n19',
  dtype(36)='amsub',     dplat(36)='n17',       dsis(36)='amsub_n17',
  dtype(37)='hirs4',     dplat(37)='metop-b',   dsis(37)='hirs4_metop-b',
  dtype(38)='amsua',     dplat(38)='metop-b',   dsis(38)='amsua_metop-b',
  dtype(39)='mhs',       dplat(39)='metop-b',   dsis(39)='mhs_metop-b',
  dtype(40)='iasi',      dplat(40)='metop-b',   dsis(40)='iasi616_metop-b',
  dtype(41)='atms',      dplat(41)='npp',       dsis(41)='atms_npp',
  dtype(42)='cris',      dplat(42)='npp',       dsis(42)='cris_npp',
  dtype(43)='sndrd1',    dplat(43)='g14',       dsis(43)='sndrD1_g14',
  dtype(44)='sndrd2',    dplat(44)='g14',       dsis(44)='sndrD2_g14',
  dtype(45)='sndrd3',    dplat(45)='g14',       dsis(45)='sndrD3_g14',
  dtype(46)='sndrd4',    dplat(46)='g14',       dsis(46)='sndrD4_g14',
  dtype(47)='sndrd1',    dplat(47)='g15',       dsis(47)='sndrD1_g15',
  dtype(48)='sndrd2',    dplat(48)='g15',       dsis(48)='sndrD2_g15',
  dtype(49)='sndrd3',    dplat(49)='g15',       dsis(49)='sndrD3_g15',
  dtype(50)='sndrd4',    dplat(50)='g15',       dsis(50)='sndrD4_g15',
 /
EOF
#
###################################################
#fila=$(cat $home_mod_scp/fila)
# Build the PBS script on-the-fly and run

cat << EOF > ${workdirSatAng}/qsub.satbang.qsb
#!/bin/csh -x
#PBS -o ${workdirSatAng}/satbang.${LABELANL}_${TIME}.out
#PBS -e ${workdirSatAng}/satbang.${LABELANL}_${TIME}.err
#PBS -j oe
#PBS -l walltime=00:05:00
#PBS -l mppwidth=24
#PBS -l mppnppn=24
#PBS -V
#PBS -S /bin/bash
#PBS -N satbang
#PBS -q pesq
#PBS -A CPTEC

cd ${workdirSatAng}/

time aprun -n 1 ${workdirSatAng}/global_angupdate

#touch ${workdirSatAng}/monitor.t
EOF

# Submetendo o processo do GSI no Gray e testando o resultado e esperando o processo
chmod 755 ${workdirSatAng}/qsub.satbang.qsb
qsub -W block=true ${workdirSatAng}/qsub.satbang.qsb

error=$?
if [ ${error} -ne 0 ]; then
  echo "       ERROR: ${GSI} crashed  Exit status=${error}"
  exit ${error}
fi

#until [ -e ${workdirSatAng}/monitor.t ]; do sleep 1s;done
#rm -fr ${workdirSatAng}/monitor.t
bout=${workdirSatAng}/satbias_ang.out
cd ..
cp -pvfr ${bout} ${save}/satbias_ang.out.${LABELANL}



if [ -e $bout ]; then

 echo -e "\033[32;1m #------------------------------------------------#\033[m"
 echo -e "\033[32;1m #                                                #\033[m"
 echo -e "\033[32;2m # Arquivo de correção de bias gerado com sucesso #\033[m"
 echo -e "\033[32;1m #                                                #\033[m"
 echo -e "\033[32;1m #------------------------------------------------#\033[m"

else

 echo -e "\033[31;1m #----------------------------------------------------#\033[m"
 echo -e "\033[31;1m #                                                    #\033[m"
 echo -e "\033[31;1m #    Falha durante o processo de geracação do Bias   # \033[m"
 echo -e "\033[31;1m #                                                    #\033[m"
 echo -e "\033[31;1m #----------------------------------------------------#\033[m"
 exit -1

fi

exit 0
