#! /bin/ksh
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $ 8128-3823
#

if [[($MAQUI = 'Linux') || ($MAQUI = 'linux')|| ($MAQUI = 'LINUX')]];then
export queue=pesq
else
export queue=Longa
fi

####################################################################
function espera_qsub {

if [ $TIMEaprox != 'indef' ]; then TIMEaprox=" / $TIMEaprox"; else TIMEaprox=""; fi
if [ x$NJOB = 'x' ]; then smslabel Info "$info Nao definido a variavel NJOB" 2>/dev/null; echo "Nao definido a variavel NJOB"; break ; fi

sleep 1
stdo=`qstat -f ${NJOB} | grep Stdout | cut -d= -f2 | awk -F ":" '{print $2}'`
stdo=`echo ${stdo} | sed s:\%s:${NJOB}:g`
rm -f ${stdo}

cerr=0
while [ 1 ]
do

  if [ -s ${stdo} ]; then
        break
  fi

  CONTRL=`cat /adm/log/qstat.log | grep ^${NJOB} | wc -w | awk '{print $1}'`
  if [ ${CONTRL} -lt 3 ]; then CONTRL=`qstat ${NJOB} |wc -l | awk '{print $1}'`; fi

  if [ ${CONTRL} -lt 3 ]; then
    ERRNO=`qstat | grep "errno" | wc -l`

    if [ ${ERRNO} -ge 01 ]; then
      smslabel Info "$info Can't connect to BatchServer" 2>/dev/null
      echo "Can't connect to BatchServer"
    else
      cerr=$(($cerr+1))
      smslabel Info "$info Verifing Files ${cerr}/12" 2>/dev/null
      echo "Verifing Files ${cerr}/12"
      if [ ${cerr} -ge 1 ]; then
        break
      fi
    fi
    sleep 5
  else
    cerr=0
    EXS=`qstat | grep ${NJOB} | awk '{print $6}'`          # Captura o status da fila
    ERRNO=`qstat | grep "errno" | wc -l`

    if [ ${ERRNO} -ge 01 ]; then
      smslabel Info "$info Can't connect to BatchServer" 2>/dev/null
      echo "Can't connect to BatchServer"
    else
      case ${EXS} in

      QUE)		                   # Status de Fila
       smslabel Info "$info Queued ..." 2>/dev/null
       echo "QUEUED..."
      ;; 
      RUN)		                   # Status de Rodando
       TIME=`qstat | grep ${NJOB} | awk '{print $10}'` # Captura a quanto tempo estah rodando
       TIME=`echo ${TIME} | cut -d"." -f1`
       smslabel Info "$info Running: ${TIME}${TIMEaprox} Secs ..." 2>/dev/null
       echo "Running: ${TIME}${TIMEaprox} Secs ..."
      ;;
      STG)		                   # Status de Stage
       smslabel Info "$info Staged ..." 2>/dev/null
       echo "STAGED..."
      ;;
      EXT)		                   # Status de Stage
       smslabel Info "$info EXT ..." 2>/dev/null
       echo "Finalizando..."
       break
      ;;
      *)		                   # Default para outros status
       smslabel Info "$info Status Unknow: ${EXS}" 2>/dev/null
       echo "Status Unknow: ${EXS}"
       break
      ;;
      esac
    fi
    sleep 5
  fi

done
}
