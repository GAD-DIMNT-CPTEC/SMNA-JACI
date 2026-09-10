#!/bin/ksh
#
#  $Author: bonatti $
#  $Date: 2008/11/13 16:44:29 $
#  $Revision: 1.3 $
#
. ../configenv_pre.ksh

varname=Chopping
ieeefiles='10,15,20,25,30,35,40,45,50,55,60,65,70,75'
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
 &ChopNML
  MendInp=${RESIN},  ! Spectral Horizontal Resolution of Input Data
  KmaxInp=${KMIN},   ! Number of Layers of Input Data
  MendOut=${RESOUT},        ! Spectral Horizontal Resolution of Output Data
  KmaxOut=${KMOUT},         ! Number of Layers of Output Data
  MendCut=${RESOUT},        ! Spectral Resolution Cut Off for Topography Smoothing
  Iter=10,            ! Number of Iteractions in Topography Smoothing
  SmthPerCut=0.12,    ! Percentage for Topography Smoothing
  GetOzone=T,         ! Flag to Produce Ozone Files
  GetTracers=F,       ! Flag to Produce Tracers Files
  GrADS=T,            ! Flag to Produce GrADS Files
  GrADSOnly=F,        ! Flag to Only Produce GrADS Files (Do Not Produce Inputs for Model)
  GDASOnly=F,         ! Flag to Only Produce Input CPTEC Analysis File
  SmoothTopo=F,       ! Flag to Performe Topography Smoothing
  RmGANL=T,           ! Flag to Remove GANL File if Desired
  LinearGrid=${SetLinear},       ! Flag to Set Linear (T) or Quadratic Gaussian Grid (T)
  DateLabel='${DATA}', ! Date Label: yyyymmddhh or DateLabel='          '
                          !       If Year (yyyy), Month (mm) and Day (dd) Are Unknown
  UTC='00',               ! UTC Hour: hh, Must Be Given if Label='          ', else UTC=' '
  NCEPName='gdas1 ',                       ! NCEP Analysis Preffix for Input File Name
  DirMain='${dirdata}/ '  ! Main User Data Directory
  DirHome='${dirhome}/ ' ! Home User Source Directory
 /
EOT0
cat <<EOT1 > ${dirrun}/aux/set${varname}.ksh
#PBS -S /bin/ksh
#PBS -V
#PBS -l cpunum_prc=1
#PBS -l tasknum_prc=1
#PBS -l memsz_job=10gb
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
export NJOB; export TIMEaprox="~126";  espera_qsub

