#!/bin/bash -x
#PBS -o eslogin01:/scratchin/grupos/dmdcc/home/denis.eiras/repo/novo_pre/branches/minimal_prototype/Mod_ChoppingParallel/Out.MPI24_Chopping_parallel
#PBS -j oe
#PBS -l walltime=01:30:00
#PBS -A CPTEC
#PBS -l mppwidth=24
#PBS -l mppnppn=24
#PBS -V
#PBS -S /bin/bash
#PBS -N Chopping
#PBS -q pesq

#
auxvarname=Chopping_parallel
if [[ (${auxvarname} = "Chopping_parallel") ]]; then
.  /opt/modules/default/etc/modules.sh
   module load stat
   module list
#   
   export PBS_SERVER=eslogin13
#   echo ${PBS_O_QUEUE} # doing nothing
else
   export PBS_SERVER=eslogin13
fi

if [[ (Linux = "Linux") || (Linux = "linux") ]]; then
  export F_UFMTENDIAN=10,15,20,25,30,35,40,45,50,55,60,65,70,75
#  export GFORTRAN_CONVERT_UNIT=big_endian:10,15,20,25,30,35,40,45,50,55,60,65,70,75
  echo "F_UFMTENDIAN = ${F_UFMTENDIAN}"
  echo "GFORTRAN_CONVERT_UNIT = ${GFORTRAN_CONVERT_UNIT}"
fi
export KMP_STACKSIZE=128m

if [[ (${auxvarname} = "FLUXCO2Clima") ]]; then
   ulimit -s 65532
else
   ulimit -s unlimited
fi

#
cd /scratchin/grupos/dmdcc/home/denis.eiras/repo/novo_pre/branches/minimal_prototype/Mod_ChoppingParallel
date

#optserver=esl
#if [[ (${optserver} = "aux") ]]; then
#/scratchin/grupos/dmdcc/home/denis.eiras/SMG_pre/datainout/bam/pre/exec/Chopping_parallel -i /scratchin/grupos/dmdcc/home/denis.eiras/SMG_pre/datainout/bam/pre/exec/Chopping_parallel.nml
#else
#time aprun -n 24 -N 24 /scratchin/grupos/dmdcc/home/denis.eiras/SMG_pre/datainout/bam/pre/exec/Chopping_parallel -i /scratchin/grupos/dmdcc/home/denis.eiras/SMG_pre/datainout/bam/pre/exec/Chopping_parallel.nml
#fi

if [[ (${auxvarname} = "Chopping_parallel") ]]; then
time aprun -n 24 -N 24 /scratchin/grupos/dmdcc/home/denis.eiras/repo/novo_pre/branches/minimal_prototype/Mod_ChoppingParallel/Chopping_parallel -i /scratchin/grupos/dmdcc/home/denis.eiras/repo/novo_pre/branches/minimal_prototype/Mod_ChoppingParallel/Chopping_parallel.nml
else
time aprun -n 1 -N 1 /scratchin/grupos/dmdcc/home/denis.eiras/SMG_pre/datainout/bam/pre/exec/Chopping_parallel -i /scratchin/grupos/dmdcc/home/denis.eiras/SMG_pre/datainout/bam/pre/exec/Chopping_parallel.nml
fi

#wait
date
#wait
> /scratchin/grupos/dmdcc/home/denis.eiras/repo/novo_pre/branches/minimal_prototype/Mod_ChoppingParallel/.Chopping_parallel.ok

