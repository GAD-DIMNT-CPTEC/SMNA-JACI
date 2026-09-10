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
then hold=""
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
#export DATA=2003010100
#
## Machine options: SX6; Linux
#MAQUI=Linux
#####################################################

varname=Chopping_parallel
ieeefiles='10,15,20,25,30,35,40,45,50,55,60,65,70,75'
direxe=${dirdata}/pre/exec
dirsrc=${dirhome}/sources/${varname}
dirout=${dirhome}/scripts/output
dirrun=${dirhome}/scripts
echo " "
host=`hostname`
echo " ${host}"
RUNTM=`date +'%d_%H:%M'`
cp -f ${dirsrc}/Delta* ${dirdata}/pre/datain/
#
cat <<EOT0 > ${direxe}/${varname}.nml
 &ChopNML
  MendInp=${RESIN},                  !: Spectral Horizontal Resolution of Input Data
  KmaxInp=${KMIN},                   !: Number of Layers of Input Data
  MendOut=${RESOUT},        ! Spectral Horizontal Resolution of Output Data
  KmaxOut=${KMOUT},         ! Number of Layers of Output Data
  MendMin=127               ! Minimum Spectral Resolution For Doing Topography Smoothing
  MendCut=${MendCut},        ! Spectral Resolution Cut Off for Topography Smoothing
  Iter=10,            ! Number of Iteractions in Topography Smoothing
  nproc_vert=1,       ! Number of processors to be used in the vertical
  ibdim_size=192,     ! Basic block-size to grid data structure
  tamBlock=512,       ! number of fft's allocated in each fft-block'
  SmthPerCut=${SmthPerCut}     !0.18,    ! Percentage for Topography Smoothing
  GetOzone=${GetOzone},         ! Flag to Produce Ozone Files
  GetTracers=${GetTracers},       ! Flag to Produce Tracers Files
  GrADS=${GrADS},            ! Flag to Produce GrADS Files
  GrADSOnly=${GrADSOnly},        ! Flag to Only Produce GrADS Files (Do Not Produce Inputs for Model)
  GDASOnly=${GDASOnly},         ! Flag to Only Produce Input CPTEC Analysis File
  SmoothTopo=${SmoothTopo},       ! Flag to Performe Topography Smoothing
  RmGANL=${RmGANL},           ! Flag to Remove GANL File if Desired
  LinearGrid=${SetLinear}, ! Flag to Set Linear or Quadratic Gaussian Grid
  givenfouriergroups=F,  ! True if nproc_vert is to be used. If false data partitioning is done automatically.
  DateLabel='${DATA}', ! Date Label: yyyymmddhh or DateLabel='          '
                          !       If Year (yyyy), Month (mm) and Day (dd) Are Unknown
  UTC='${UTC}',               ! UTC Hour: hh, Must Be Given if Label='          ', else UTC=' '
  NCEPName='${AnlPref} ',                       ! NCEP Analysis Preffix for Input File Name
  DirMain='${dirdata}/ '  ! Main User Data Directory
  DirHome='${dirhome}/ ' ! Home User Source Directory
 /
EOT0


export PBS_SERVER=${pbs_server1}
optserver=`printf "$PBS_SERVER \n" | cut -c1-3`
if [[ (${optserver} = "aux") ]]; then
export MPPBS="#"
else
export MPPBS="#PBS -l mppwidth=${cpu_mpi}"
fi

cat <<EOT1 > ${direxe}/set${varname}.bash
#!/bin/bash -x
#PBS -o ${host}:${direxe}/${varname}.MPI${cpu_mpi}
#PBS -j oe
#PBS -l walltime=${AUX_WALLTIME}
${MPPBS}
#PBS -l mppnppn=${cpu_node}
#PBS -A ${QUOTA}
#PBS -V
#PBS -S /bin/bash
#PBS -N $RES
#PBS -q ${QUEUE}
#

.   /opt/modules/default/etc/modules.sh
module load stat
module list
#
export PBS_SERVER=${pbs_server1}

echo $PBS_O_QUEUE
if [[ (${MAQUI} = "Linux") || (${MAQUI} = "linux") ]]; then
  export F_UFMTENDIAN=${ieeefiles}
  export GFORTRAN_CONVERT_UNIT=big_endian:${ieeefiles}
  echo ${F_UFMTENDIAN}
  echo "GFORTRAN_CONVERT_UNIT= " ${GFORTRAN_CONVERT_UNIT}
fi
#
export KMP_STACKSIZE=128m
ulimit -s unlimited

cd ${direxe}
date
optserver=`printf "$PBS_SERVER \n" | cut -c1-3`
if [[ (\${optserver} = "aux") ]]; then
${direxe}/${varname} -i ${direxe}/${varname}.nml
else
time aprun -n ${cpu_mpi} -N ${cpu_node} ${direxe}/${varname} -i ${direxe}/${varname}.nml
fi
date

EOT1
#
#   Change mode to be executable
#
#qsub   ${direxe}/set${varname}.bash
${QSUB} $hold ${direxe}/set${varname}.bash
if [[ ${it} -eq 1 ]];then
#FIRST=`qsub ${direxe}/set${varname}.bash`
export FIRST
echo $FIRST
else
#SECOND=`qsub -W depend=afterok:$FIRST ${direxe}/set${varname}.bash`
echo $SECOND
fi


#it=2
#while [ ${it} -gt 0 ];do
#it=`qstat @aux20 @eslogin13 | grep $USER | grep $RES | wc -l`
#done

