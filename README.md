# LQCD-DWF-solver: solver benchmark


An iterative solver for lattice QCD with domain-wall-type fermions.
It is a kernelized version of the proxy application `LQCD-DWF-HMC`
(https://github.com/i-kanamori/LQCD-DWF-HMC), which is based on
Bridge++ (https://bridge.kek.jp/Lattice-code/).
Some of the architecture-specific implementations differ from those in the original Bridge++.


The benchmark solves an even-odd preconditioned linear equation for 12 different right-hand sides and calculates meson two-point functions, i.e., several partial norms of the solution vectors.

Except for the general version (`corelib` below), it uses a CG solver with a mixed-precision scheme. The general version uses double precision only. Both versions use even-odd (red-black) preconditioning.

This program uses hybrid OpenMP/MPI parallelization with `MPI_THREAD_FUNNELED`.
The nearest-neighbor communications use MPI persistent communication. The 4-dimensional structured lattice is divided among the MPI processes.

## License

GPLv3.

## Target-specific implementations

* Fugaku: The Fujitsu extension of MPI is available (specified at compile time). The SIMD degrees of freedom are assigned to the x-y directions of the structured lattice, and the implementation uses ACLE.
* SIMD version for x86 (`alt_simd`): uses AVX/AVX2/AVX-512, and
the SIMD degrees of freedom are assigned to the x-direction.
* GPU version (OpenACC or CUDA, `alt_accel`): assumes that the x-direction is used for the vector (or warp) degrees of freedom.
* Non-SIMD data layout (`corelib`): should run on any
architecture with C++11 support.


## Compilation

Edit the Makefile for your target. Several samples are provided:
 * `Makefile_fugaku_acle`: uses ACLE and FJMPI, with the clang mode of the Fujitsu compiler
 * `Makefile_fugaku_acle_fapp`: uses ACLE and FJMPI, with the clang mode of the Fujitsu compiler and the fapp profiler
 * `Makefile_fugaku_general`: uses the data layout for Fugaku without ACLE or FJMPI, with the clang mode of the Fujitsu compiler
 * `Makefile_corelib_only`: without the data layout for Fugaku
 * `Makefile_simd_avx512`: uses AVX-512 intrinsics on x86 machines
 * `Makefile_simd_avx2`: uses AVX2 intrinsics on x86 machines
 * `Makefile_openacc`: uses OpenACC on NVIDIA GPUs
 * `Makefile_cuda`: uses CUDA on NVIDIA GPUs

### Important switches in the Makefile

 * `target`: specifies the compiler. See `Makefile_target.inc` for the details of each target, including the compiler options.
 * `use_mpi`: set `fjmpi` to use the Fujitsu extension of persistent communication; set `gpu_aware` to use GPU-aware MPI.
 * `use_alternative`: set `yes` to use architecture-specific tuning
 * `use_alt_qxs`: set `yes` to use A64FX (Fugaku)-specific tuning
 * `use_qxs_arch`: set `general` to emulate the ACLE implementation
 * `use_alt_accel` (OpenACC/CUDA version): set `yes` to use the GPU version
 * `use_device_arch` (OpenACC/CUDA version): sets the target GPU generation
 * `use_alt_simd` (SIMD version): set `yes` to use the SIMD version for x86 machines
 * `simd_arch` (SIMD version): sets the target SIMD architecture

### How to build

```
# build the library
make
# build the benchmark executable
cd benchmark/domainwall
make
```


## Execution

See the `submit_*` and `sample_x86*` directories for examples. In short:

```
mpiexec bridge.elf alt_qxs   # optimized for A64FX (with or without ACLE)
mpiexec bridge.elf alt_accel # runs on GPUs (with OpenACC/CUDA)
mpiexec bridge.elf alt_simd  # optimized for x86 (with AVX/AVX2/AVX-512)
mpiexec bridge.elf corelib   # baseline before tuning
```

## Input description

The program reads 2 input YAML files and requires 1 configuration data file
(`conf_04040408.txt`), which is specified in one of the YAML files.

The YAML files can be edited if needed.

 * `main.yaml`: specifies the problem size (`lattice_size`), the MPI process grid (`grid_size`), and the number of OpenMP threads (`number_of_thread`). `verbose_level` controls the amount of logging. The other fields should be kept as they are.

 * `test_alt_Spectrum_Domainwall_5din_Hadron2ptFunction.yaml`: specifies the detailed parameters, which should be left unchanged. It also specifies the configuration data file (`conf_04040408.txt`); edit the path if the file is located elsewhere. To see the convergence behavior of the solver, set `verbose_level` in the `Solver_xxx` section to `Detailed`. If the run takes too long, increase `quark_mass` in the `Fopr` section (any value is allowed, but in practice it should be at most 1.0).


## Output example

All output is written to stdout.

### Results for the solution vectors

The lines following `Solving quark propagator:` show the label of the source vector (`color` and `spin`) and the residual norm.
```
Solving quark propagator:
  color spin   Nconv      diff           diff2
Richardson  nconv = 4, diff2 =  4.203788057279973e-32
    0    0        4     2.0503e-16     2.6184e-31
Richardson  nconv = 4, diff2 =  5.276297017670305e-32
    1    0        4     2.2970e-16     1.3779e-31
Richardson  nconv = 4, diff2 =  4.887325589330260e-32
    2    0        4     2.2107e-16     1.4508e-31
...
```
The `Richardson ...` lines do not appear in the `corelib` case.


Partial sums of several linear combinations of the solution vectors
are printed in the `2-point correlator:` section. These values
provide a non-trivial check of the solver and the matrix-vector multiplication
kernel.
```
2-point correlator:
PS <-- PS correlator:
     0    1.622108950386e-01    2.635501058632e-20
     1    4.982525189422e-02    2.095970631002e-20
...
```
The first value (`1.622108950386e-01` above) is used to check the result,
which is reported in the `test result:` line near the end of the output.
The value is compared with `expected_result` in the `Test_Spectrum` section of the input file; the value provided there is valid only for the lattice size 32x8x8x12.
