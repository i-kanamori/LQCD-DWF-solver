#!/bin/bash
#PJM -L node=12:torus
#PJM -L elapse=0:05:00
#PJM -L rscgrp=small
#PJM -g ra000001
#PJM --mpi proc=48
#PJM -j
#PJM -S

#export OMP_NUM_THREADS=12
export FLIB_BARRIER=HARD

BIN=../benchmark/domainwall/bridge.elf

#mpiexec $BIN corelib
mpiexec $BIN alt_qxs

