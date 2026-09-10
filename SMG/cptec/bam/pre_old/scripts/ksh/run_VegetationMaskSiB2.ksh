#!/bin/ksh
#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#
. ../configenv_pre.ksh

varname=VegetationMaskSiB2
ieeefiles='20,30,50'
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
  Imax=${IM},     ! Number of Longitudes at Output Gaussian Grid
  Jmax=${JM},     ! Number of Latitudes at Output Gaussian Grid
  Idim=360,    ! Number of Longitudes in Climatological SSiB Vegetation Mask Data
  Jdim=180,    ! Number of Latitudes in Climatological SSiB Vegetation Mask Data
  GrADS=.TRUE. ! Flag for GrADS Outputs
  DirMain='${dirdata}/ ' ! Main Datain/Dataout Directory
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

