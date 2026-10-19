BIN=../benchmark/domainwall/bridge.elf

# data layout for x86 machine (one-dim simd tilying for lattice sites), RIRI... for complex numbers
mpiexec -np 4 $BIN alt_simd

