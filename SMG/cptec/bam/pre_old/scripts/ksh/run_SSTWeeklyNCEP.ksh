#!/bin/ksh
#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#
. ../configenv_pre.ksh

varname=SSTWeeklyNCEP
ieeefiles='10,20'
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
  Idim=360,         ! Number of Longitudes at the NCEP SST Weekly Grid
  Jdim=180,         ! Number of Latitudes at the NCEP SST Weekly Grid
  Date='${DATA}' ! Date of Initial Conditions (yyyymmddhh)
  GrADS=.TRUE.      ! Flag for GrADS Outputs
  DirMain='${dirdata}/ ' ! Main Datain/Dataout Directory
 /
EOT0

#
#  Run GetSST.ksh
#  To deGRIB NCEP SST file with the script GetSST.ksh
#  and to format properly the file for SSTWeekly)
#
hh=00
cd ${dirdata}/pre/datain
rm -f dump
rec=`${dirgrads}/bin/wgrib -s -4yr -d 1 -ieee gdas1.T${hh}Z.sstgrb.${DATA}`
echo ${rec}
date=`awk 'BEGIN {print substr("'${rec}'",7,10)}'`
echo ${date}
#
mv dump gdas1.T${hh}Z.sstgrd.${date}


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
# Run SSTWeekly
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

