#!/bin/bash
#-----------------------------------------------------------------------------#
#           Group on Data Assimilation Development - GDAD/CPTEC/INPE          #
#-----------------------------------------------------------------------------#
#BOP
#
# !DESCRIPTION: script auxiliar utilizado pelo runGSI. Este script possui todas
#               as funções utilizadas pelo script runGSI.
#
#
# !REVISION HISTORY:
# 20 Ago 2018 - J. G. de Mattos - Initial Version
#
# !REMARKS:
#
#    * para o correto funcionamento deve-se primeiro carregar este script
#      com o comando source
#      - source runGSI_functions.sh
#    * O passo seguinte deve ser executar a função constants!
#
# !BUGS:
#
#    Not yet!
#
#EOP
#-----------------------------------------------------------------------------#
#BOC
#-----------------------------------------------------------------------------#
# assign value $2 to variable $1
#-----------------------------------------------------------------------------#
assign(){
eval export $1=$2
}

#-----------------------------------------------------------------------------#
# return a subword of a string
#-----------------------------------------------------------------------------#

subwrd ( ) {
   str=$(echo "${@}" | awk '{ for (i=1; i<=NF-1; i++) printf("%s ",$i)}')
   n=$(echo "${@}" | awk '{ print $NF }')
   echo "${str}" | awk -v var=${n} '{print $var}'
}

#-----------------------------------------------------------------------------#
# return usage from main program
#-----------------------------------------------------------------------------#
usage ( ){
   echo
   echo "Usage:"
   sed -n '/^#BOP/,/^#EOP/{/^#BOP/d;/^#EOP/d;p}' ${runGSI}
}

#-----------------------------------------------------------------------------#
# Define some constants to GSI use
#-----------------------------------------------------------------------------#
constants ( ) {

   # Get local run of the runGSI
    for i in $(seq 0 $((${#FUNCNAME[@]}-1)) );do
       if [ ${FUNCNAME[$i]} == 'main' ];then
          export runGSIDir=$(dirname $(readlink -e ${BASH_SOURCE[$i]}))
          export runGSI=$(readlink -e ${BASH_SOURCE[$i]})
          break
       fi
    done

   # Define endianness
   export BYTE_ORDER=Big_Endian
   #export BYTE_ORDER=Little_Endian

   # caminho arquivos Little -> cd /pesq/share/das/public_tmp/GSI3.4/tmp/crtm-2.4.0_emc.1/fix

   ##### CAROL 12AGT2026 ---> arranjo para a JACI
   # 1. Mant�m a sua vari�vel BYTE_ORDER do jeito que est� nas v�rias linhas do script
   #export BYTE_ORDER=Big_Endian

   # 2. Converte para mai�sculas e remove sublinhados para o padr�o da Intel
   #    Transforma "Big_Endian" -> "BIG_ENDIAN"
   #INTEL_BYTE_ORDER=$(echo "$BYTE_ORDER" | tr '[:lower:]' '[:upper:]' | sed 's/_ENDIAN/_ENDIAN/')

   # 3. Exporta as duas vari�veis que o compilador Intel (ifort/ifx/ftn) reconhece
   #export FORT_CONVERT_BYTE_ORDER="$INTEL_BYTE_ORDER"
   #export FORT_CONVERT_BIG_ENDIAN=1

   # 1. Garante que o ambiente padr�o continue Little Endian (para o BAM e Sat�lites)
   #unset FORT_CONVERT_BYTE_ORDER
   #unset FORT_CONVERT_BIG_ENDIAN

   # 2. Ativa a convers�o por sintaxe com par�nteses/sublinhado que o Bash aceita
   #export FORT_CONVERT_BIG_ENDIAN="10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,40,41,42,43,44,45,46,47,48,49,50"


   # PBS default variables
   #
   #   * Some of these variables
   #     can be modified by command
   #     line
  case ${hpc_name} in

    egeon)
        export MaxCoresPerNode=60
        export MTasks=120                     # Number of Processors
        export ThreadsPerMPITask=1             # Number of cores hosting OpenMP threads
        export TasksPerNode=$((${MaxCoresPerNode}/${ThreadsPerMPITask})) # Number of Processors used by each MPI tasks
        export PEs=$((${MTasks}/${ThreadsPerMPITask}))
        export Nodes=$(((${MTasks}+${MaxCoresPerNode}-1)/${MaxCoresPerNode}))
        #export Queue=PESQ1
        #export Queue=PESQ2
	export Queue=batch
        export WallTime=01:00:00
        export BcCycles=0

        # MPI environmental variables
        export MPICH_UNEX_BUFFER_SIZE=100000000
        export MPICH_MAX_SHORT_MSG_SIZE=4096
        export MPICH_PTL_UNEX_EVENTS=50000
        export MPICH_PTL_OTHER_EVENTS=2496
    ;;
    jaci)
        export MaxCoresPerNode=128
        export MTasks=120                     # Number of Processors
        export ThreadsPerMPITask=1             # Number of cores hosting OpenMP threads
        export TasksPerNode=$((${MaxCoresPerNode}/${ThreadsPerMPITask})) # Number of Processors used by each MPI tasks
        export PEs=$((${MTasks}/${ThreadsPerMPITask}))
        export Nodes=$(((${MTasks}+${MaxCoresPerNode}-1)/${MaxCoresPerNode}))
        #export Queue=PESQ1
        #export Queue=PESQ2
	export Queue=pesqmidi
        export WallTime=01:00:00
        export BcCycles=0

        # MPI / Cray PALS environmental variables
        export MPICH_NO_BUFFER_ALIAS_CHECK=TRUE
        export OMP_STACKSIZE=1G
        export KMP_STACKSIZE=1G
    ;;
    XC50)
       export MaxCoresPerNode=40
       export MPITasks=120                    # Number of Processors
       export ThreadsPerMPITask=1             # Number of cores hosting OpenMP threads
       export TasksPerNode=$((${MaxCoresPerNode}/${ThreadsPerMPITask})) # Number of Processors used by each MPI tasks
       export PEs=$((${MPITasks}/${ThreadsPerMPITask}))
       export Nodes=$(((${MPITasks}+${MaxCoresPerNode}-1)/${MaxCoresPerNode}))
       export Queue=pesq
       export WallTime=01:00:00
       export BcCycles=0

       # MPI environmental variables
       export MPICH_UNEX_BUFFER_SIZE=100000000
       export MPICH_MAX_SHORT_MSG_SIZE=4096
       export MPICH_PTL_UNEX_EVENTS=50000
       export MPICH_PTL_OTHER_EVENTS=2496
   ;;

 esac

   # files used by this script
   export execGSI=${home_cptec}/bin/gsi.x
   export parmGSI=${home_gsi_fix}/gsiparm.anl
   export obsGSI=${home_gsi_fix}/obsfiles.rc
   export cldRadInfo=${home_gsi_fix}/cloudy_radiance_info.txt
#   export SatBiasSample=${public_fix}/comgsi_satbias_in
#   export SatBiasPCSample=${public_fix}/comgsi_satbias_pc_in
# Colocando arquivos de BC do SMNA ver issue 85 do readDiag
#   export SatBiasSample=${public_fix}/gdas1.t00z.abias
#   export SatBiasPCSample=${public_fix}/gdas1.t00z.abias_pc
   export SatBiasSample=${public_fix}/satbias_out_50ciclos
   export SatBiasPCSample=${public_fix}/satbias_pc.out_50ciclos
   export ScanInfo=${home_gsi_fix}/global_scaninfo.txt
   export SatBiasAngSample=${home_gsi_fix}/global_satangbias.txt
   export execBCAng=${home_cptec}/bin/global_angupdate
   export parmBCAng=${home_gsi_fix}/global_angupdate.namelist

   # Satbias files
   export satbiasIn=satbias_in
   export satbiasOu=satbias_out
   export satbiasPCIn=satbias_pc
   export satbiasPCOu=satbias_pc.out
   export satbiasAngIn=satbias_ang.in
   export satbiasAngOu=satbias_ang.out


}
#-----------------------------------------------------------------------------#
# Get BAM trunc
#-----------------------------------------------------------------------------#
BAM_CoordSize ( ){

   trunc=${1}

   case ${trunc} in
      21)TimeStep=3600; IMax=64; JMax=32;;
      31)TimeStep=1800; IMax=96; JMax=48;;
      42)TimeStep=1800; IMax=128; JMax=64;;
      62)TimeStep=900; IMax=192; JMax=96;;
      106)TimeStep=900; IMax=320; JMax=160;;
      126)TimeStep=600; IMax=384; JMax=192;;
      133)TimeStep=600; IMax=400; JMax=200;;
      159)TimeStep=600; IMax=480; JMax=240;;
      170)TimeStep=450; IMax=512; JMax=256;;
      213)TimeStep=300; IMax=640; JMax=320;;
      254)TimeStep=255; IMax=768; JMax=384;;
      299)TimeStep=200; IMax=900; JMax=450;;
      319)TimeStep=225; IMax=960; JMax=480;;
      341)TimeStep=200; IMax=1024; JMax=512;;
      382)TimeStep=180; IMax=1152; JMax=576;;
      511)TimeStep=150; IMax=1536; JMax=768;;
      533)TimeStep=150; IMax=1600; JMax=800;;
      666)TimeStep=240; IMax=2000; JMax=1000;;
      863)TimeStep=150; IMax=2592; JMax=1296;;
      1279)TimeStep=20; IMax=3840; JMax=1920;;
      1332)TimeStep=20; IMax=4000; JMax=2000;;
      *)echo "Truncamento desconhecido ${trunc}"
   esac
}

#-----------------------------------------------------------------------------#
# Get options from runGSI line command
#-----------------------------------------------------------------------------#
ParseOpts( ) {

   i=1
   flag=0
   while true; do

      arg=$(echo "${@}" | awk -v var=${i} '{print $var}')
      i=$((i+1))

      if [ -z ${arg} ];then
         break;
      fi

      while true ;do

         if [ ${arg} = '-t' ]; then BkgTrunc=$(subwrd ${@} ${i}); i=$((i+1)); break; fi
         if [ ${arg} = '-l' ]; then BkgNLevs=$(subwrd ${@} ${i}); i=$((i+1)); break; fi
         if [ ${arg} = '-p' ]; then BkgPrefix=$(subwrd ${@} ${i}); i=$((i+1)); break; fi

         if [ ${arg} = '-I' ]; then AnlDate=$(subwrd ${@} ${i}); i=$((i+1)); break; fi
         if [ ${arg} = '-T' ]; then AnlTrunc=$(subwrd ${@} ${i}); i=$((i+1)); break; fi
         if [ ${arg} = '-bc' ]; then BcCycles=$(subwrd ${@} ${i}); i=$((i+1)); break; fi

         if [ ${arg} = '-np' ]; then MPITasks=$(subwrd ${@} ${i}); i=$((i+1)); break; fi
         if [ ${arg} = '-N' ]; then TasksPerNode=$(subwrd ${@} ${i}); i=$((i+1)); break; fi
         if [ ${arg} = '-d' ]; then ThreadsPerMPITask=$(subwrd ${@} ${i}); i=$((i+1)); break; fi

         if [ ${arg} = '-om' ]; then ObsMod=$(subwrd ${@} ${i}); i=$((i+1)); break; fi

         if [ ${arg} = '-h' ]; then cat < ${0} | sed -n '/^#BOP/,/^#EOP/p' ; i=$((i+0)); exit 0; fi

         flag=1
         i=$((i-1))

         break;
      done
      if [ ${flag} -eq 1 ]; then break; fi
   done

   # Truncamento do background
   if [ -z ${BkgTrunc} ];then
      echo -e "\e[31;1m >> Erro: \e[m\e[33;1m Truncamento do background não foi passado\e[m"
      usage
      exit -1
   fi

   # Numero de niveis do background
   if [ -z ${BkgNLevs} ];then
      echo -e "\e[31;1m >> Erro: \e[m\e[33;1m Número de níveis do background não foi passado\e[m"
      usage
      exit -1
   fi

   # Prefixo do arquivo de background
   if [ -z ${BkgPrefix} ];then
      echo -e "\e[31;1m >> Erro: \e[m\e[33;1m Prefixo do arquivo de background não foi passado\e[m"
      usage
      exit -1
   fi

   # Data para a produção da Análise
   if [ -z ${AnlDate} ];then
      echo -e "\e[31;1m >> Erro: \e[m\e[33;1m A data da análise não foi passada\e[m"
      usage
      exit -1
   fi

   # Execução com 1 (oneobtest) ou mais (fullobs) observações
   if [ -z ${ObsMod} ]; then
     ObsMod="false"
   fi

   # Truncamento em que a análise será calculada
   # Obs: Este é o mesmo truncamento utilizado no
   #      cálculo prévio da matriz de covariância
   #      do erro. Então esta será a resolução em
   #      que o GSI irá processar internamente a
   #      análise. Pode ser diferente do Truncamento
   #      do modelo.
   if [ -z ${AnlTrunc} ];then
      AnlTrunc=${BkgTrunc}
   fi

   # Não haverá modificações na coordenada vertical
   export AnlNLevs=${BkgNLevs}

   # Configurando datas
   export BkgDate=$(${inctime} ${AnlDate} -6h %y4%m2%d2%h2)
   export FctDate=$(${inctime} ${AnlDate} +6h %y4%m2%d2%h2)

   # Configurando Sulfixo dos arquivos do Modelo
   export AnlMRES=$(printf 'TQ%4.4dL%3.3d' ${AnlTrunc} ${AnlNLevs})
   export BkgMRES=$(printf 'TQ%4.4dL%3.3d' ${BkgTrunc} ${BkgNLevs})
}

#-----------------------------------------------------------------------------#
# Choose and link/copy observations from ObsDir to GSI dir run
#-----------------------------------------------------------------------------#
linkObs ( ){
   local runDate=${1}
   local runDir=${2}

   local verbose=true
   
   # add bufr path in EGEON
   local obsDir=${ncep_ext}/${runDate:0:4}/${runDate:4:2}/${runDate:6:2}
# local para ler os dados de datas recentes na opera��o (ultimos 30 dias)   
#   local obsDir=${obsDir}:${ncep_ext}/${runDate:0:8}00/dataout/NCEP
# local para testes com o fomrato do arquivo prepbufr do INPE   
#   local obsDir=/mnt/beegfs/luiz.sapucci/obs_V1.3.3/dataout/${runDate:0:8} 
# local para ler dados salvos em experimentos   
#   local obsDir=/pesq/dados/das/poper/luiz.sapucci/NCEPdataSMNA/PrepBufr/${runDate}
   local obsDir=/oper/dados/dboper/raw/arch/mod/ncep/gdas/${runDate:0:4}/${runDate:4:2}/${runDate:6:2}
   
#   local obsDir=${obsDir}:/lustre_xc50/ioper/data/external/ASSIMDADOS
#   local obsDir=${obsDir}:/lustre_xc50/joao_gerd/data/${runDate}
#   local obsDir=${obsDir}:/lustre_xc50/joao_gerd/data/obs/${runDate:0:6}/${runDate:6:4}

   local IFS=":"; read -a obsPath < <(echo "${obsDir}")
   local nPaths=${#obsPath[@]}
   echo "@linkObs : obsPath : " $obsPath
   echo "obsGSI : "${obsGSI}
   local IFS=" "
   local names=$(sed -n '/OBS_INPUT::/,/::/{/OBS_INPUT/d;/::/d;/^!/d;p}' ${parmGSI} | awk '{print $1}' | sort -u | xargs)
   echo "parmGSI : names : " $names
   count=0
   for name in ${names};do
      i=0
      IsObsList=$(grep -iw ${name} ${obsGSI} | wc -l)
      if [ $IsObsList -eq 0 ]; then
         echo -e " " 
         echo $name" não está na lista de "$obsGSI;
         echo -e " "
         continue                                     # skip to the next name in $names
      fi
      while [ $i -le $((nPaths-1)) ];do
         filemask=${obsPath[$i]}/$(grep -iw ${name} ${obsGSI} | awk '{print $1}'; exit ${PIPESTATUS[0]} )
         if [ $? -eq 0 ];then
            file=$(${inctime} ${runDate} +0h ${filemask})
            ## change -e to -f avoid that cases when grep return empty and the if test the whole obsPath 
            #  which causes a entire directory copy to $rundir/$name ....
            if [ -f ${file} ];then
               cp -pfr  ${file} ${runDir}/${name} 2> /dev/null

               if [ $? -eq 0 ];then
                  echo -en "\033[34;1m[\033[m\033[32;1m OK\033[m"
                  count=$((count + 1))
               else
                  echo -en "\033[34;1m[\033[m\033[31;1m FAIL\033[m"
               fi

               if [ ${verbose} == 'true' ];then
                  echo -e "\033[34;1m ]\033[m\033[34;1m link\033[m\033[37;1m ${name}\033[m @ [ ${file} ]" 
               else
                  echo -e "\033[34;1m ]\033[m\033[34;1m link\033[m\033[37;1m ${name}\033[m"
               fi

               break

            elif [ $i -eq $((nPaths-1)) ];then
               echo -e " "
               echo -e "\033[31;1m File not found \033[m\033[34;1m $(basename ${file})\033[m\033[31;1m\033[m"
               echo -e " "
            fi
         else
            echo -e " "
            echo -e "\033[31;1m Dont seach for \033[m\033[34;1m ${name}\033[m\033[31;1m observation file \033[m"
            echo -e "\033[31;1m Observation not included in\033[m\033[34;1m ${obsGSI}\033[m\033[31;1m file ! \033[m"
            echo -e " "
         fi
         i=$((i+1))
      done

   done
   echo -e "\033[34;1m Found\033[m\033[32;1m ${count}\033[m\033[34;1m observation files to use\033[m"
}

#-----------------------------------------------------------------------------#
# Copy fixed files to use in GSI run
#-----------------------------------------------------------------------------#
FixedFiles ( ) {

   local Trunc=${1}
   local NLevs=${2}
   local runDir=${3}

   local mres=$(printf 'TQ%4.4dL%3.3d' ${Trunc} ${NLevs})
   echo "Copy of Public fixed files"
   # Public fixed files
   # (GSI)
   ### cp -pfr ${public_fix}/gsir4.berror_stats.gcv.BAM.${mres} ${runDir}/berror_stats  
   cp -pfr ${public_fix}/gsir4.berror_stats.gcv.BAM.TQ0254L064 ${runDir}/berror_stats        ### using TQ0254 for TQ0299 for now
   cp -pfr ${public_fix}/atms_beamwidth.txt                 ${runDir}/atms_beamwidth.txt
   cp -pfr ${public_fix}/global_ozinfo.txt                  ${runDir}/ozinfo
   cp -pfr ${public_fix}/global_pcpinfo.txt                 ${runDir}/pcpinfo
   cp -pfr ${public_fix}/prepobs_errtable.global            ${runDir}/errtable
   if [ ${ObsMod} == "true" ]; then
     cp -pfr ${public_fix}/prepobs_prep.bufrtable             ${runDir}/prepobs_prep.bufrtable
   fi
   cp -pfr ${public_fix}/bufrtab.012                        ${runDir}/bftab_sstphr

   # (CRTM)
#   ln -sf ${public_crtm}/${BYTE_ORDER}/Nalli.IRwater.EmisCoeff.bin    ${runDir}/Nalli.IRwater.EmisCoeff.bin
#   ln -sf ${public_crtm}/${BYTE_ORDER}/NPOESS.IRice.EmisCoeff.bin     ${runDir}/NPOESS.IRice.EmisCoeff.bin
#   ln -sf ${public_crtm}/${BYTE_ORDER}/NPOESS.IRland.EmisCoeff.bin    ${runDir}/NPOESS.IRland.EmisCoeff.bin
#   ln -sf ${public_crtm}/${BYTE_ORDER}/NPOESS.IRsnow.EmisCoeff.bin    ${runDir}/NPOESS.IRsnow.EmisCoeff.bin
#   ln -sf ${public_crtm}/${BYTE_ORDER}/NPOESS.VISice.EmisCoeff.bin    ${runDir}/NPOESS.VISice.EmisCoeff.bin
#   ln -sf ${public_crtm}/${BYTE_ORDER}/NPOESS.VISland.EmisCoeff.bin   ${runDir}/NPOESS.VISland.EmisCoeff.bin
#   ln -sf ${public_crtm}/${BYTE_ORDER}/NPOESS.VISsnow.EmisCoeff.bin   ${runDir}/NPOESS.VISsnow.EmisCoeff.bin
#   ln -sf ${public_crtm}/${BYTE_ORDER}/NPOESS.VISwater.EmisCoeff.bin  ${runDir}/NPOESS.VISwater.EmisCoeff.bin
#   ln -sf ${public_crtm}/${BYTE_ORDER}/FASTEM5.MWwater.EmisCoeff.bin  ${runDir}/FASTEM5.MWwater.EmisCoeff.bin
#   ln -sf ${public_crtm}/${BYTE_ORDER}/AerosolCoeff.bin               ${runDir}/AerosolCoeff.bin
#   ln -sf ${public_crtm}/${BYTE_ORDER}/CloudCoeff.bin                 ${runDir}/CloudCoeff.bin

   # como estava antes -> 12AGT2026
   cp  ${public_crtm}/${BYTE_ORDER}/Nalli.IRwater.EmisCoeff.bin    ${runDir}/Nalli.IRwater.EmisCoeff.bin
   cp  ${public_crtm}/${BYTE_ORDER}/NPOESS.IRice.EmisCoeff.bin     ${runDir}/NPOESS.IRice.EmisCoeff.bin
   cp  ${public_crtm}/${BYTE_ORDER}/NPOESS.IRland.EmisCoeff.bin    ${runDir}/NPOESS.IRland.EmisCoeff.bin
   cp  ${public_crtm}/${BYTE_ORDER}/NPOESS.IRsnow.EmisCoeff.bin    ${runDir}/NPOESS.IRsnow.EmisCoeff.bin
   cp  ${public_crtm}/${BYTE_ORDER}/NPOESS.VISice.EmisCoeff.bin    ${runDir}/NPOESS.VISice.EmisCoeff.bin
   cp  ${public_crtm}/${BYTE_ORDER}/NPOESS.VISland.EmisCoeff.bin   ${runDir}/NPOESS.VISland.EmisCoeff.bin
   cp  ${public_crtm}/${BYTE_ORDER}/NPOESS.VISsnow.EmisCoeff.bin   ${runDir}/NPOESS.VISsnow.EmisCoeff.bin
   cp  ${public_crtm}/${BYTE_ORDER}/NPOESS.VISwater.EmisCoeff.bin  ${runDir}/NPOESS.VISwater.EmisCoeff.bin
   cp  ${public_crtm}/${BYTE_ORDER}/FASTEM5.MWwater.EmisCoeff.bin  ${runDir}/FASTEM5.MWwater.EmisCoeff.bin
   cp  ${public_crtm}/${BYTE_ORDER}/AerosolCoeff.bin               ${runDir}/AerosolCoeff.bin
   cp  ${public_crtm}/${BYTE_ORDER}/CloudCoeff.bin                 ${runDir}/CloudCoeff.bin


   # Carol - test CRTM -> JACI 12AGT2026
   #cp  ${plus_crtm}/EmisCoeff/IR_Water/${BYTE_ORDER}/Nalli.IRwater.EmisCoeff.bin               ${runDir}/Nalli.IRwater.EmisCoeff.bin
   #cp  ${plus_crtm}/EmisCoeff/IR_Ice/SEcategory/${BYTE_ORDER}/NPOESS.IRice.EmisCoeff.bin       ${runDir}/NPOESS.IRice.EmisCoeff.bin
   #cp  ${plus_crtm}/EmisCoeff/IR_Land/SEcategory/${BYTE_ORDER}/NPOESS.IRland.EmisCoeff.bin     ${runDir}/NPOESS.IRland.EmisCoeff.bin
   #cp  ${plus_crtm}/EmisCoeff/IR_Snow/SEcategory/${BYTE_ORDER}/NPOESS.IRsnow.EmisCoeff.bin     ${runDir}/NPOESS.IRsnow.EmisCoeff.bin
   #cp  ${plus_crtm}/EmisCoeff/VIS_Ice/SEcategory/${BYTE_ORDER}/NPOESS.VISice.EmisCoeff.bin     ${runDir}/NPOESS.VISice.EmisCoeff.bin
   #cp  ${plus_crtm}/EmisCoeff/VIS_Land/SEcategory/${BYTE_ORDER}/NPOESS.VISland.EmisCoeff.bin   ${runDir}/NPOESS.VISland.EmisCoeff.bin
   #cp  ${plus_crtm}/EmisCoeff/VIS_Snow/SEcategory/${BYTE_ORDER}/NPOESS.VISsnow.EmisCoeff.bin   ${runDir}/NPOESS.VISsnow.EmisCoeff.bin
   #cp  ${plus_crtm}/EmisCoeff/VIS_Water/SEcategory/${BYTE_ORDER}/NPOESS.VISwater.EmisCoeff.bin ${runDir}/NPOESS.VISwater.EmisCoeff.bin
   #cp  ${plus_crtm}/EmisCoeff/MW_Water/${BYTE_ORDER}/FASTEM5.MWwater.EmisCoeff.bin  ${runDir}/FASTEM5.MWwater.EmisCoeff.bin
   #cp  ${plus_crtm}/AerosolCoeff/${BYTE_ORDER}/AerosolCoeff.bin             ${runDir}/AerosolCoeff.bin
   #cp  ${plus_crtm}/CloudCoeff/${BYTE_ORDER}/CloudCoeff.bin                 ${runDir}/CloudCoeff.bin

   # User fixed files
   cp -pfr ${home_gsi_fix}/global_anavinfo.l${NLevs}.txt   ${runDir}/anavinfo
   cp -pfr ${home_gsi_fix}/global_satinfo.txt              ${runDir}/satinfo
   cp -pfr ${home_gsi_fix}/global_convinfo.txt             ${runDir}/convinfo
   cp -pfr ${home_gsi_fix}/global_scaninfo.txt             ${runDir}/scaninfo

   #
   #-------------------------------------------------------------------------------#
   # Copy CRTM coefficient files based on entries in satinfo file
   #
   for file in $(awk '{if($1!~"!"){print $1}}' ${runDir}/satinfo | sort | uniq) ;do
      #ln -sf ${plus_crtm}/SpcCoeff/${BYTE_ORDER}/${file}.SpcCoeff.bin ${runDir}
      cp ${plus_crtm}/SpcCoeff/${BYTE_ORDER}/${file}.SpcCoeff.bin ${runDir}
      # ln -s ${public_crtm}/${BYTE_ORDER}/${file}.SpcCoeff.bin ${runDir}
      #ln -sf ${plus_crtm}/TauCoeff/ODAS/${BYTE_ORDER}/${file}.TauCoeff.bin ${runDir}
       
      # Carlos (02/07/2025) - estou colocando aqui os dados para o cálculo da profundidade óptica do algorítmo ODPS
      # por ser mais adequado para dados do microondas (no caso do ATMS, os arquivos atms_*.TauCoeff.bin - com excessão do npp
      # só esão disponíveis pelo ODPS, pelo pacote crtm-2.4.0_emc.1 disponível no GitHub)
      cp ${plus_crtm}/TauCoeff/ODPS/${BYTE_ORDER}/${file}.TauCoeff.bin ${runDir}
#      cp -v ${plus_crtm}/TauCoeff/${BYTE_ORDER}/${file}.TauCoeff.bin ${runDir}

      # ln -s ${public_crtm}/${BYTE_ORDER}/${file}.TauCoeff.bin ${runDir}
   done

   emiscoef_IRwater=${plus_crtm}/EmisCoeff/IR_Water/${BYTE_ORDER}/Nalli.IRwater.EmisCoeff.bin
   emiscoef_IRice=${plus_crtm}/EmisCoeff/IR_Ice/SEcategory/${BYTE_ORDER}/NPOESS.IRice.EmisCoeff.bin
   emiscoef_IRland=${plus_crtm}/EmisCoeff/IR_Land/SEcategory/${BYTE_ORDER}/NPOESS.IRland.EmisCoeff.bin
   emiscoef_IRsnow=${plus_crtm}/EmisCoeff/IR_Snow/SEcategory/${BYTE_ORDER}/NPOESS.IRsnow.EmisCoeff.bin
   emiscoef_VISice=${plus_crtm}/EmisCoeff/VIS_Ice/SEcategory/${BYTE_ORDER}/NPOESS.VISice.EmisCoeff.bin
   emiscoef_VISland=${plus_crtm}/EmisCoeff/VIS_Land/SEcategory/${BYTE_ORDER}/NPOESS.VISland.EmisCoeff.bin
   emiscoef_VISsnow=${plus_crtm}/EmisCoeff/VIS_Snow/SEcategory/${BYTE_ORDER}/NPOESS.VISsnow.EmisCoeff.bin
   emiscoef_VISwater=${plus_crtm}/EmisCoeff/VIS_Water/SEcategory/${BYTE_ORDER}/NPOESS.VISwater.EmisCoeff.bin
   emiscoef_MWwater=${plus_crtm}/EmisCoeff/MW_Water/${BYTE_ORDER}/FASTEM6.MWwater.EmisCoeff.bin
   aercoef=${plus_crtm}/AerosolCoeff/${BYTE_ORDER}/AerosolCoeff.bin
   cldcoef=${plus_crtm}/CloudCoeff/${BYTE_ORDER}/CloudCoeff.bin

#   ln -sf $emiscoef_IRwater ${runDir}/Nalli.IRwater.EmisCoeff.bin
#   ln -sf $emiscoef_IRice ${runDir}/NPOESS.IRice.EmisCoeff.bin
#   ln -sf $emiscoef_IRsnow ${runDir}/NPOESS.IRsnow.EmisCoeff.bin
#   ln -sf $emiscoef_IRland ${runDir}/NPOESS.IRland.EmisCoeff.bin
#   ln -sf $emiscoef_VISice ${runDir}/NPOESS.VISice.EmisCoeff.bin
#   ln -sf $emiscoef_VISland ${runDir}/NPOESS.VISland.EmisCoeff.bin
#   ln -sf $emiscoef_VISsnow ${runDir}/NPOESS.VISsnow.EmisCoeff.bin
#   ln -sf $emiscoef_VISwater ${runDir}/NPOESS.VISwater.EmisCoeff.bin
#   ln -sf $emiscoef_MWwater ${runDir}/FASTEM6.MWwater.EmisCoeff.bin
   cp $emiscoef_MWwater ${runDir}/FASTEM6.MWwater.EmisCoeff.bin
#   ln -sf $aercoef  ${runDir}/AerosolCoeff.bin
#   ln -sf $cldcoef  ${runDir}/CloudCoeff.bin

   echo

}
#-----------------------------------------------------------------------------#
# Link/Copy info Files (satinfo, convinfo, ozinfo). These files are dependent
# of the anl time.
#-----------------------------------------------------------------------------#
getInfoFiles (){
   anlDate=${1}
   infoFiles=${home_gsi_fix}/infoFiles
   for infoFile in satinfo convinfo ozinfo;do
      while read file;do
         infoDate=${file##*.}
         #echo -e $infoFile "${file##*.} -- $anlDate"
         if [ $anlDate -ge ${infoDate} ];then
            assign ${infoFile} ${file}
         fi
      done < <(ls -1 ${infoFiles}/global_*${infoFile}* )
   done

   cp -pvfr ${satinfo}  ${runDir}/satinfo
   cp -pvfr ${convinfo} ${runDir}/convinfo
   cp -pvfr ${ozinfo}   ${runDir}/ozinfo

}
#-----------------------------------------------------------------------------#
# Link/Copy coefficient files for both air bias correction and angle dependent
# bias into the GSI run directory before running the GSI executable.
#-----------------------------------------------------------------------------#
getSatBias ( ){

   local AnlDate=${1}
   local bkgDate=${2}
   local runDir=${3}
   local DataOut=${4}

   local cold=0

   # Find for a bias corrected file

   # if are running a bias correcion cycle
   # find in current GSI running dir by
   # satbias file
   FindDir01=${runDir}

   # but nomaly find from a previous cycle
   FindDir02=${DataOut}/${BkgDate}

   local FileSatbiasOu=$(find ${FindDir01} ${FindDir02} -iname "${satbiasOu}*" -print 2> /dev/null | sort -nr | head -n 1)

   if [ ${#FileSatbiasOu} -eq 0 ];then

     echo -e "\033[31;1m #--------------------------------------#\033[m"
     echo -e "\033[31;1m #     Usando o os coeficientes de BC     #\033[m"
     echo -e "\033[31;1m #        de 50 ciclos do SMNA            #\033[m"
     echo -e "\033[31;1m #     Caso não seja o 1° ciclo verificar  #\033[m"
     echo -e "\033[31;1m #    porque não copiou o arquivo anterior  #\033[m"
     echo -e "\033[31;1m #              ${satbiasOu}               #\033[m"
     echo -e "\033[31;1m #--------------------------------------#\033[m"

     if [ -e ${runDir}/${satbiasIn} ];then
        #
        # if file exist, may be corrupted!
        #
        rm -fr ${runDir}/${satbiasIn}
     fi

     echo cp -pfr ${SatBiasSample} ${runDir}/${satbiasIn}
     cp -pfr ${SatBiasSample} ${runDir}/${satbiasIn}

     cold=1

   else

     echo cp -pfr ${FileSatbiasOu} ${runDir}/${satbiasIn}
     cp -pfr ${FileSatbiasOu} ${runDir}/${satbiasIn}
     
   fi

#   local FileSatbiasAngOu=$(find ${FindDir01} ${FindDir02} -iname "${satbiasAngOu}*" -print 2> /dev/null | sort -nr | head -n 1)

##     See AdvancedGSIUserGuide_v3.5.0.0.pdf , Ming Hu
##     Using 8.4.4 Enhanced Radiance Bias Correction, so SatbiasAng is not used 
#
#  if [ ${#FileSatbiasAngOu} -eq 0 ];then
#
#     echo -e "\033[31;1m #--------------------------------------#\033[m"
#     echo -e "\033[31;1m #    1° ciclo de assimilação           #\033[m"
#     echo -e "\033[31;1m #  Caso não seja o 1° ciclo verificar  #\033[m"
#     echo -e "\033[31;1m #    porque não copiou o arquivo       #\033[m"
#     echo -e "\033[31;1m #          ${satbiasAngOu}         #\033[m"
#     echo -e "\033[31;1m #--------------------------------------#\033[m"
#
#     if [ -e ${runDir}/${satbiasAngIn} ];then
#        #
#        # if file exist, may be corrupted!
#        #
#        rm -fr ${runDir}/${satbiasAngIn}
#     fi
#
#     echo cp -pfr ${SatBiasAngSample} ${runDir}/${satbiasAngIn}
#     cp -pfr ${SatBiasAngSample} ${runDir}/${satbiasAngIn}
#
#   else
#
#     echo cp -pfr ${FileSatbiasAngOu} ${runDir}/${satbiasAngIn}
#     cp -pfr ${FileSatbiasAngOu} ${runDir}/${satbiasAngIn}
#
#   fi

   local FileSatbiasPCOu=$(find ${FindDir01} ${FindDir02} -iname "${satbiasPCOu}*" -print 2> /dev/null | sort -nr | head -n 1)

   if [ ${#FileSatbiasPCOu} -eq 0 ];then

     echo -e "\033[31;1m #--------------------------------------#\033[m"
     echo -e "\033[31;1m #    Usando o os coeficientes de BC_PC   #\033[m"
     echo -e "\033[31;1m #        de 50 ciclos do SMNA            #\033[m"
     echo -e "\033[31;1m #  Caso não seja o 1° ciclo verificar  #\033[m"
     echo -e "\033[31;1m #    porque não copiou o arquivo       #\033[m"
     echo -e "\033[31;1m #          ${satbiasPCOu}         #\033[m"
     echo -e "\033[31;1m #--------------------------------------#\033[m"

     if [ -e ${runDir}/${satbiasPCIn} ];then
        #
        # if file exist, may be corrupted!
        #
        rm -fr ${runDir}/${satbiasPCIn}
     fi

     echo cp -pfr ${SatBiasPCSample} ${runDir}/${satbiasPCIn}
     cp -pfr ${SatBiasPCSample} ${runDir}/${satbiasPCIn}

   else
     echo cp -pfr ${FileSatbiasPCOu} ${runDir}/${satbiasPCIn}
     cp -pfr ${FileSatbiasPCOu} ${runDir}/${satbiasPCIn}

   fi

   return ${cold}

}


#-----------------------------------------------------------------------------#
# Function to submit GSI job at tupa supercomputer
#-----------------------------------------------------------------------------#
subGSI() {
   runDir=${1}
   imax=${2}
   jmax=${3}
   kmax=${4}
   atrc=${5}
   btrc=${6}
   andt=${7}
   onet=${8}
   cold=${9}

   cp ${cldRadInfo} ${runDir}/$(basename ${cldRadInfo})

   if [ $cold == '.true.' ];then
     echo "  using "${parmGSI}.cold " to setup GSI run"
      sed "s/#CENTER#/cptec/g" ${parmGSI}.cold > ${runDir}/$(basename ${parmGSI})
   else
     echo "  using "${parmGSI} " to setup GSI run"
     sed "s/#CENTER#/cptec/g" ${parmGSI} > ${runDir}/$(basename ${parmGSI})
   fi

   sed -i -e "s/#IMAX#/${imax}/g" \
          -e "s/#JMAX#/${jmax}/g" \
          -e "s/#KMAX#/${kmax}/g" \
          -e "s/#TRUNC#/${atrc}/g" \
          -e "s/#TRUNCB#/${btrc}/g" \
          -e "s/#NLVL#/$((${kmax}-1))/g" \
          -e "s/#NLV#/${kmax}/g" \
          -e "s/#LABELANL#/${andt}/g" \
          -e "s/#ONEOBTEST#/${onet}/" ${runDir}/$(basename ${parmGSI})

  runTime=$(date '+runTime-%H:%M:%S')

  #
  # Sanity Check
  #
  sanity=$((${TasksPerNode}*${ThreadsPerMPITask}))
  if [ ${sanity} -ne ${MaxCoresPerNode} ];then
     echo -e "\e[31;1m >> Erro: \e[m\e[33;1m Redefina Numero de Processos MPI e openMP\e[m"
     exit -1
  fi
ana_date=${andt:4:6}
case ${hpc_name} in
  egeon)
     cat << EOF > ${runDir}/gsi.qsb
#!/bin/bash
#SBATCH --nodes=${Nodes}
#SBATCH --time=${WallTime}
#SBATCH --ntasks=${MTasks}
#SBATCH --job-name=GSI-SMNA
#SBATCH --mem=480G
#SBATCH --cpus-per-task=1
#SBATCH --partition=${Queue}

# Enable ro debug after run gsi
# must use Stack Trace Analysis Tool (STAT)
export ATP_ENABLED=1

cd ${runDir}
pwd

export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK

source ${home_gsi}/env.sh egeon ${compiler}

module list

echo  "STARTING AT `date` "

ulimit -c unlimited
ulimit -s unlimited

export I_MPI_DEBUG=15

mpirun -np \${SLURM_NTASKS} ./$(basename ${execGSI}) > gsiStdout_${andt}.${runTime}.log

EOF

    cd ${runDir}

    PID=$(sbatch -W  gsi.qsb; exit ${PIPESTATUS[0]})
  ;;
  jaci)
     cat << EOF > ${runDir}/gsi.qsb
#!/bin/bash
#PBS -o ${runDir}/Out.gsi.${andt}.${runTime}.out
#PBS -j oe
#PBS -l walltime=${WallTime}
#PBS -l select=${Nodes}:ncpus=${MaxCoresPerNode}:mpiprocs=${MTasks}
#PBS -V
#PBS -S /bin/bash
#PBS -N GSI-SMNA
#PBS -q ${Queue}
#PBS -l place=scatter:excl

set -euo pipefail

cd ${runDir}
pwd

# Configuracoes de Threads e Memoria para o GSI
export OMP_NUM_THREADS=${ThreadsPerMPITask}
#export OMP_STACKSIZE=1G
#export KMP_STACKSIZE=1G
#export KMP_AFFINITY=disabled
#export MPICH_NO_BUFFER_ALIAS_CHECK=TRUE
#export BYTE_ORDER=Little_Endian

ulimit -c unlimited
ulimit -s unlimited

# 1. Carrega PrgEnv-cray para o 'module swap' do env.sh n�o falhar no n� de login
module load PrgEnv-cray 2>/dev/null || true

# 2. Inicializa as vari�veis como vazias ("") apenas para o Bash n�o reclamar de 'unbound variable'
export FFLAGS="${FFLAGS:-}"
export FCFLAGS="${FCFLAGS:-}"

# Carrega ambiente e modulos
source ${home_gsi}/env.sh jaci ${compiler}

# Forca a carga do PALS para garantir que o mpiexec se comunique com o PBS Pro
module load cray-pals 2>/dev/null || true
module load cray-libpals 2>/dev/null || true

module list

echo "STARTING AT $(date)"

# Execucao do GSI via Cray PALS
time mpiexec -n ${MTasks} --depth=${ThreadsPerMPITask} ./$(basename ${execGSI}) > gsiStdout_${andt}.${runTime}.log 2>&1

EOF

     cd ${runDir}

     PID=$(qsub -W block=true gsi.qsb; exit ${PIPESTATUS[0]})
     #qstatus=\$?
  ;;
  XC50)
     cat << EOF > ${runDir}/gsi.qsb
#!/bin/bash
#PBS -o \${PBS_O_WORKDIR}/gsiAnl.${andt}.${runTime}.out
#PBS -S /bin/bash
#PBS -q ${Queue}
#PBS -l nodes=${Nodes}:ppn=${MaxCoresPerNode}
#PBS -N gsiAnl
#PBS -A CPTEC

ulimit -c unlimited
ulimit -s unlimited

# Enable ro debug after run gsi
# must use Stack Trace Analysis Tool (STAT)
export ATP_ENABLED=1

cd \${PBS_O_WORKDIR}

aprun -n ${PEs} -N ${TasksPerNode} -d ${ThreadsPerMPITask} $(basename ${execGSI}) > gsiStdout_${andt}.${runTime}.log

EOF

     cd ${runDir}

     PID=$(qsub -W block=true gsi.qsb; exit ${PIPESTATUS[0]})

   ;;
esac
return $?
exit 1
}

#-----------------------------------------------------------------------------#
# Function to submet radiance Bias Correction from satellite angle
#-----------------------------------------------------------------------------#
subBCAng() {

   runDir=${1}
   andt=${2}

   cp ${parmBCAng} ${runDir}/$(basename ${parmBCAng})

   sed -i -e "s/#YEAR#/${andt:0:4}/g" \
          -e "s/#MONTH#/${andt:4:2}/g" \
          -e "s/#DAY#/${andt:6:2}/g" \
          -e "s/#HOUR#/${andt:8:2}/g" ${runDir}/$(basename ${parmBCAng})

  runTime=$(date '+runTime-%H:%M:%S')

case ${hpc_name} in
  egeon)
     cat << EOF > ${runDir}/angupdate.qsb
#!/bin/bash
#SBATCH --nodes=${Nodes}
#SBATCH --time=00:15:00
#SBATCH --ntasks=${MTasks}
#SBATCH --job-name=AngUp${ana_date}
#SBATCH --mem=480G
#SBATCH --cpus-per-task=1
#SBATCH --partition=${Queue}

ulimit -c unlimited
ulimit -s unlimited

# Enable ro debug after run gsi
# must use Stack Trace Analysis Tool (STAT)
export ATP_ENABLED=1

cd ${runDir}
pwd

export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK

source ${home_gsi}/env.sh egeon ${compiler}

export LD_LIBRARY_PATH=${home_gsi}/libsrc/crtm/lib64:${LD_LIBRARY_PATH}
module load intel
module load impi
module load curl-7.85.0-gcc-9.4.0-qbney7y
module load cmake/3.21.3
module load openblas
module load netcdf
module load netcdf-fortran


ldd ./$(basename ${execBCAng})
module list

echo  "STARTING AT `date` "

ulimit -c unlimited
ulimit -s unlimited

export I_MPI_DEBUG=15

mpirun -np \${SLURM_NTASKS} ./$(basename ${execBCAng}) > gsiAngUpdateStdout_${andt}.${runTime}.log

EOF

    cd ${runDir}

    PID=$(sbatch -W  angupdate.qsb; exit ${PIPESTATUS[0]})
  ;;
  XC50)

cat << EOF > ${runDir}/angupdate.qsb
#!/bin/bash
#PBS -o \${PBS_O_WORKDIR}/gsiAngUpdate.${andt}.${runTime}.out
#PBS -S /bin/bash
#PBS -q ${Queue}
#PBS -l nodes=1:ppn=40
#PBS -l walltime=00:05:00
#PBS -N gsiAngUpdate
#PBS -A CPTEC

# Enable to debug after run gsi angupdate
# must use Stack Trace Analysis Tool (STAT)
export ATP_ENABLED=1

cd \${PBS_O_WORKDIR}

aprun -n 40 $(basename ${execBCAng}) > gsiAngUpdateStdout_${andt}.${runTime}.log

EOF

cd ${runDir}

PID=$(qsub -W block=true angupdate.qsb; exit ${PIPESTATUS[0]})
   ;;
esac

return $?

}
teste(){
echo 'hello'
}
#-----------------------------------------------------------------------------#
# Function to Merge diagnostic files from GSI run
#-----------------------------------------------------------------------------#
mergeDiagFiles ( ) {

   local runDir=${1}
   local AnDate=${2}

   local miter=$(printf "%2.2d" $(($(grep -i miter ${parmGSI} | awk -F"," '{print $1}' | awk -F"=" '{print $2}')+1)))
   local files=$(find -P -O3 ${runDir} -maxdepth 1 -iname "pe*.*" -printf '%f\n'| awk -F"." '{print $2}' | sort -u)

   for file in $(echo ${files});do

      #tmpv=${file##pe0000.}
      #loop=$(echo ${tmpv}|awk -F"_" '{print $NF}')
      #type=${tmpv%%_${loop}}
      loop=$(echo ${file}|awk -F"_" '{print $NF}')
      type=${file%%_${loop}}
      echo -en "\033[34;1mMerging diag files \033[m\033[32;1m${type}\033[m\033[34;1m, outer loop \033[m\033[32;1m${loop}\033[m \033[34;1m[\033[m"

      cat $(find -P -O3 ${runDir} -maxdepth 1 -iname "pe*${type}*${loop}" | sort -n) > ${runDir}/diag_${type}_${loop}.${AnDate} 2> /dev/null

      if [ $? -eq 0 ];then
         echo -e "\033[32;1m OK \033[m\033[34;1m]\033[m"

#         # remove pe's files
#         find ${runDir} -iname "pe*${type}*${loop}" -exec rm -f {} \;

         # copy pe's files to diag dir #
         #
         # This procedure will be removed later
         #
         if [ ! -e ${runDir}/diag ];then
            mkdir -p ${runDir}/diag
         fi

         find -P -O3 ${runDir} -maxdepth 1 -iname "pe*${type}*${loop}" -exec cp -f {} ${runDir}/diag \;

         # link to be used by global_angleupdate
         #if [ ${loop} = ${miter} ];then
         #   ln -sf ${runDir}/diag_${type}_${loop}.${AnDate} ${runDir}/diag_${type}.${AnDate}
         #fi
      else
         echo -e "\033[31;1m FAIL \033[m\033[34;1m]\033[m"
      fi

   done


}
#-----------------------------------------------------------------------------#
# Function to copy files generated by gsi and gsiAngUpdate
#-----------------------------------------------------------------------------#
copyFiles (){

   local fromDir=${1}
   local toDir=${2}

   local CP='cp -pf'
#   local MV='mv -f'
   local MV='cp -pf'

   # arquivos gerados pelo gsi
   find -P -O3 ${fromDir} -maxdepth 1 -type f -iname "diag_*" -exec mv -f {} ${toDir} \;
   find -P -O3 ${fromDir} -maxdepth 1 -type f -iname "fort.*" -exec mv -f {} ${toDir} \;

   mv -f ${fromDir}/diag ${toDir}

   ${MV} ${fromDir}/BAM.anl ${toDir}/GANL${BkgPrefix}${AnlDate}S.unf.${BkgMRES}
   # copy para manter 
   ${CP} ${fromDir}/${satbiasIn} ${toDir}/.
   echo ${CP} ${fromDir}/${satbiasOu} ${toDir}/.
   ${CP} ${fromDir}/${satbiasOu} ${toDir}/.
   ${CP} ${fromDir}/${satbiasPCIn} ${toDir}/.
   echo ${CP} ${fromDir}/${satbiasPCOu} ${toDir}/.
   ${CP} ${fromDir}/${satbiasPCOu} ${toDir}/.
#   ${CP} ${fromDir}/${satbiasAngIn} ${toDir}/.
#   echo ${CP} ${fromDir}/${satbiasAngOu} ${toDir}/.
#   ${CP} ${fromDir}/${satbiasAngOu} ${toDir}/.

   # arquivos de configuracao
   ${MV} ${fromDir}/gsiparm.anl ${toDir}
   ${CP} ${fromDir}/anavinfo ${toDir}
   ${CP} ${fromDir}/satinfo  ${toDir}
   ${CP} ${fromDir}/convinfo ${toDir}

   # arquivos de log
#   ${MV} ${fromDir}/gsiAnl.* ${toDir}
   ${MV} ${fromDir}/gsiStdout* ${toDir}
#   mv -f ${fromDir}/gsiAngUpdate* ${toDir}

}
#EOC
#-----------------------------------------------------------------------------#
function .log () {
  local LOG_LEVELS=([0]="emerg"  [1]="alert"   [2]="crit"   [3]="err"    [4]="warning" [5]="notice" [6]="info"   [7]="debug")
  local LOG_COLORS=([0]="[31;1m" [1]="[33;1m" [2]="[31;1m" [3]="[31;1m" [4]="[33;1m"  [5]="[32;1m" [6]="[34;1m" [7]="[32;1m")
  local LEVEL=${1}
  local __VERBOSE=${#LOG_LEVELS[@]}
  shift
  if [ ${__VERBOSE} -ge ${LEVEL} ]; then
    echo -e "\033${LOG_COLORS[$LEVEL]}[${LOG_LEVELS[$LEVEL]}]\033[m" "$@ ${LOG_COLORS}"
  fi
}
