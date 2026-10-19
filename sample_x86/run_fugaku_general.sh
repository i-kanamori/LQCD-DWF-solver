BIN=../benchmark/domainwall/bridge.elf

# data layout for fugaku (two-dim simd tilying for lattice sites), RRRR...IIII... for complex numbers
mpiexec -np 4 $BIN alt_qxs

