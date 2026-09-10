#!/bin/bash
# 
# script to run CPTEC Global Model on PC Clusters under MPI Scali
# and Sun Grid Engine without OpenMP
#
# assumptions: assume present at the same directory:
#              ParModel_MPI (Global Model Executable file)
#              MODELIN (Global Model input Namelist file)
#
# usage: run_multi_UNA cpu_mpi cpu_node name hold
# where:
# cpu_mpi: integer, the desired number of mpi processes
# cpu_node: integer, the desired number of mpi processes per shared memory node
# name: character, the job name (for SGE)
# hold: any, present or not;
#            if absent, script finishes after queueing job;
#            if present, script holds till job completion
if [ "$#" == 4 ]
then hold="-sync y"
else hold=
fi
export FEXE=`pwd`
export cpu_mpi=$1
export cpu_node=$2
export RES=$3
num=$(($cpu_mpi+$cpu_node-1))
fra=$(($num/$cpu_node))
cpu_tot=$(($fra*$cpu_node))
echo fila=mpi-npn${cpu_node} total cpus=${cpu_tot}

######################################################
#dirhome=/home/pkubota/mcgacptec/pre
#dirdata=/mpp/pkubota/mcgacptec
#dirgrads=/usr/local/grads
#
# Set  Res for Chopping
#RESIN=382
#KMIN=64
#RESOUT=042
#KMOUT=28
#
#SetLinear=F
#
#RESO=042
#IM=128
#JM=64
#prefix=${JM}
##set run date
#export DATA=2007090300
#
## Machine options: SX6; Linux
#MAQUI=Linux
#####################################################

varname=Chopping
ieeefiles='10,15,20,25,30,35,40,45,50,55,60,65,70,75'
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
  MendMin=127               ! Minimum Spectral Resolution For Doing Topography Smoothing
  MendCut=${RESOUT},        ! Spectral Resolution Cut Off for Topography Smoothing
  Iter=10,            ! Number of Iteractions in Topography Smoothing
  SmthPerCut=0.12,    ! Percentage for Topography Smoothing
  GetOzone=T,         ! Flag to Produce Ozone Files
  GetTracers=F,       ! Flag to Produce Tracers Files
  GrADS=T,            ! Flag to Produce GrADS Files
  GrADSOnly=F,        ! Flag to Only Produce GrADS Files (Do Not Produce Inputs for Model)
  GDASOnly=F,         ! Flag to Only Produce Input CPTEC Analysis File
  SmoothTopo=F,       ! Flag to Performe Topography Smoothing
  RmGANL=F,           ! Flag to Remove GANL File if Desired
  LinearGrid=${SetLinear},           ! Flag to Set Linear or Quadratic Gaussian Grid
  DateLabel='${DATA}', ! Date Label: yyyymmddhh or DateLabel='          '
                          !       If Year (yyyy), Month (mm) and Day (dd) Are Unknown
  UTC='12',               ! UTC Hour: hh, Must Be Given if Label='          ', else UTC=' '
  NCEPName='gdas1 ',                       ! NCEP Analysis Preffix for Input File Name
  DirMain='${dirdata}/ '  ! Main User Data Directory
  DirHome='${dirhome}/ ' ! Home User Source Directory
 /
EOT0
cat <<EOT1 > ${dirrun}/aux/set${varname}.bash
#!/bin/bash
#$ -pe mpi-npn${cpu_node} ${cpu_tot}
#$ -o una1:${FEXE}/Out.MPI${cpu_mpi}
#$ -j y
#$ -V
#$ -S /bin/bash
#$ -N $RES
#
if [[ (${MAQUI} = "Linux") || (${MAQUI} = "linux") ]]; then
  export F_UFMTENDIAN=${ieeefiles}
  export GFORTRAN_CONVERT_UNIT=big_endian:${ieeefiles}
  echo ${F_UFMTENDIAN}
  echo "GFORTRAN_CONVERT_UNIT= " ${GFORTRAN_CONVERT_UNIT}
fi
#
export KMP_STACKSIZE=128m
ulimit -s unlimited

cd ${dirrun}/aux
time ${direxe}/${varname} 
EOT1
#
#   Change mode to be executable
#
chmod +x ${dirrun}/aux/set${varname}.bash
qsub $hold ${dirrun}/aux/set${varname}.bash
