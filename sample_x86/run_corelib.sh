BIN=../benchmark/domainwall/bridge.elf

# general data layout: not for simd, RIRIRI... for complex numbers
mpiexec -np 4 $BIN corelib

