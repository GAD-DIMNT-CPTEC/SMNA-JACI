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

#################################################
#export dirhome=/stornext/home/paulo.kubota/agcm_cptec/pre
#export dirdata=/scratchin/grupos/pad/home/paulo.kubota/agcm_cptec
#export dirgrads=/usr/local/grads
#
## Machine options: SX6; Linux
#export MAQUI=Linux
#
## Set  Res for Chopping
#export RESIN=382
#export KMIN=64
#export RESOUT=42
#export KMOUT=28
#export SetLinear=FALSE
#export RESO=42
#export IM=128
#export JM=64
#export prefix=${JM}
##set run date
#export DATA=2003010100
####################################################

varname=VegetationMaskIBISClima
ieeefiles='10,20,30'
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
cat <<EOT0 > ${direxe}/${varname}.nml
 &InputDim
  Idim=720,     ! Number of Longitudes in SSiB Vegetation Mask Data
  Jdim=360,     ! Number of Latitudes in SSiB Vegetation Mask Data
  GrADS=.TRUE., ! Flag for GrADS Outputs
  DirMain='${dirdata}/ ' ! Main Datain/Dataout Directory
 /

EOT0

export PBS_SERVER=${pbs_server2}
optserver=`printf "$PBS_SERVER \n" | cut -c1-3`
if [[ (${optserver} = "aux") ]]; then
export MPPBS="#"
else
export MPPBS="#PBS -l mppwidth=${cpu_mpi}"
fi


cat <<EOT1 > ${direxe}/set${varname}.bash
#!/bin/bash
#PBS -o ${host}:${direxe}/Out.MPI${cpu_mpi}
#PBS -j oe
#PBS -l walltime=${AUX_WALLTIME}
#PBS -A ${QUOTA}
${MPPBS}
#PBS -l mppnppn=${cpu_node}
#PBS -V
#PBS -S /bin/bash
#PBS -N $RES
#PBS -q ${AUX_QUEUE}

#
export PBS_SERVER=${pbs_server2}


if [[ (${MAQUI} = "Linux") || (${MAQUI} = "linux") ]]; then
  export F_UFMTENDIAN=${ieeefiles}
  export GFORTRAN_CONVERT_UNIT=big_endian:${ieeefiles}
  echo ${F_UFMTENDIAN}
  echo "GFORTRAN_CONVERT_UNIT= " ${GFORTRAN_CONVERT_UNIT}
fi
export KMP_STACKSIZE=128m
ulimit -s unlimited


#
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
chmod +x ${direxe}/set${varname}.bash
${QSUB} $hold ${direxe}/set${varname}.bash
it=2
while [ ${it} -gt 0 ];do
it=`qstat @aux20 @eslogin13 | grep $USER | grep $RES | wc -l`
done

