#!/bin/ksh
#
#  $Author: bonatti $
#  $Date: 2007/09/18 18:07:15 $
#  $Revision: 1.2 $
#
. ../configenv_pre.ksh

varname=SSTWeekly
ieeefiles='10,30,50,60'
#
direxe=${dirdata}/pre/exec
dirsrc=${dirhome}/sources/${varname}
dirout=${dirhome}/scripts/output
dirrun=${dirhome}/scripts
echo " "
host=`hostname`
echo " ${host}"
RUNTM=`date +'%d_%H:%M'`
#
cat <<EOT0 > ${dirrun}/aux/${varname}.nml
 &InputDim
  Mend=${RESO},             ! Spectral Resolution Horizontal Truncation
  Kmax=${KMOUT},              ! Number of Layers of the Initial Condition for the Global Model
  Idim=360,             ! Number of Longitudes For Climatological SST Data
  Jdim=180,             ! Number of Latitudes For Climatological SST Data
  SSTSeaIce=-1.749,     ! SST Value in Celsius Degree Over Sea Ice (-1.749 NCEP, -1.799 CAC)
  LatClimSouth=-50.0,   ! Southern Latitude For Climatological SST Data
  LatClimNorth=60.0,    ! Northern Latitude For Climatological SST Data
  ClimWindow=.FALSE.,   ! Flag to Climatological SST Data Window
  Linear=.TRUE.,        ! Flag for Linear (T) or Area Weighted (F) Interpolation
  LinearGrid=.FALSE.,   ! Flag to Set Linear (T) or Quadratic Gaussian Grid (T)
  GrADS=.TRUE.,         ! Flag for GrADS Outputs
  DateICn='${DATA}', ! Date of th Initial Condition for the Global Model
  Preffix='GANLNMC',    ! Preffix of the Initial Condition for the Global Model
  Suffix='S.unf.',      ! Suffix of the Initial Condition for the Global Model
  DirMain='${dirdata}/ '             ! Main Datain/Dataout Directory
 /
EOT0
cat <<EOT1 > ${dirrun}/aux/set${varname}.ksh
#PBS -S /bin/ksh
#PBS -V
#PBS -l cpunum_prc=1
#PBS -l tasknum_prc=1
#PBS -l memsz_job=1gb
#PBS -q ${queue}
#PBS -o ${host}-e:${dirout}/${varname}.${RUNTM}.out
#PBS -e ${host}-e:${dirout}/${varname}.${RUNTM}.err
#
if [[ (${MAQUI} = "Linux") || (${MAQUI} = "linux") ]]; then
  export F_UFMTENDIAN=${ieeefiles}
  export GFORTRAN_CONVERT_UNIT=big_endian:${ieeefiles}
  echo ${F_UFMTENDIAN}
  echo "GFORTRAN_CONVERT_UNIT= " ${GFORTRAN_CONVERT_UNIT}
fi
#
cd ${dirrun}/aux
time ${direxe}/${varname} 
EOT1
#
#   Change mode to be executable
#
chmod +x ${dirrun}/aux/set${varname}.ksh
echo qsub -q ${queue} ${dirrun}/aux/set${varname}.ksh
submet=`/usr/bin/nqsII/qsub -q ${queue} -N ${varname} ${dirrun}/aux/set${varname}.ksh`
NJOB=`echo $submet |awk '{print $2}' |awk -F "." '{print $1}'`
export NJOB; export TIMEaprox="~2";  espera_qsub

