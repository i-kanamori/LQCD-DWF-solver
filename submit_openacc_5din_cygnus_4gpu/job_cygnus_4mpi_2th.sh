#!/bin/bash
#PBS -A BRIDGE
#PBS -q gpu
#PBS -b 1
#PBS -l elapstim_req=0:30:00
#PBS -T openmpi
#PBS -v NQSV_MPI_VER=gdr/4.1.5/nvhpc22.11_gcc8.3.1-cuda11.8

module load openmpi/${NQSV_MPI_VER}

#export MV2_ENABLE_AFFINITY=0

cd ${PBS_O_WORKDIR}/
export OMP_NUM_THREADS=2

mpiexec ${NQSII_MPIOPTS} -np 4 -npernode 4 -x UCX_MAX_RNDV_RAILS=4 ${PBS_O_WORKDIR}/../benchmark/domainwall/bridge.elf alt_accel
