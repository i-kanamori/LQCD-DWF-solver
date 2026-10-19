#####  GNU C/C++  ##############################################

ifeq ($(strip $(target)),PC_GNU)

ifeq ($(use_mpi),yes)
  CXX = mpic++
else
  CXX = g++
endif

ifeq ($(use_thread),omp)
  CXXFLAGS += -fopenmp
  LDFLAGS  += -fopenmp
endif

ifeq ($(use_opt_level),low)
  OPTFLAGS = -O0
endif
ifeq ($(use_opt_level),std)
  OPTFLAGS = -O
endif
ifeq ($(use_opt_level),high)
  OPTFLAGS = -O3
endif

ifeq ($(use_static),yes)
  LDFLAGS  += -static
endif

ifeq ($(use_debug),yes)
  CXXFLAGS += -g -Wall -Wno-unused -Wno-sign-compare -Wno-unknown-pragmas -fbounds-check -fsanitize=undefined -fno-sanitize-recover
  CPP_DEFS += -DDEBUG
  LDFLAGS += -fsanitize=undefined -fno-sanitize-recover
else
  CPP_DEFS += -DNDEBUG
endif

ifeq ($(use_cpp11),yes)
  CXXFLAGS += -std=gnu++11
endif

#CXXFLAGS += -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64
CXXFLAGS += $(OPTFLAGS)

ifeq ($(use_alt_simd),yes)
  ifeq ($(simd_arch),AVX1)
    CXXFLAGS += -mavx
  endif
  ifeq ($(simd_arch),AVX2)
    CXXFLAGS += -mavx2 -mfma
  endif
  ifeq ($(simd_arch),AVX512)
    CXXFLAGS += -mavx512f
  endif
endif

LD = $(CXX)
LDLIBS = -lm

CXXDEP = $(CXX)
CXXDEPFLAGS = -MM

AR = ar r
RANLIB = ranlib

TARGET_MACHINE = PC_GNU

#####  CLANG C/C++  ##############################################

else ifeq ($(strip $(target)),PC_CLANG)

ifeq ($(use_mpi),yes)
  CXX = mpic++
else
  CXX = clang++
endif

ifeq ($(use_thread),omp)
  CXXFLAGS += -fopenmp
  LDFLAGS  += -fopenmp
endif

ifeq ($(use_opt_level),low)
  OPTFLAGS = -O0
endif
ifeq ($(use_opt_level),std)
  OPTFLAGS = -O
endif
ifeq ($(use_opt_level),high)
  OPTFLAGS = -O3
endif

ifeq ($(use_static),yes)
  LDFLAGS  += -static
endif

ifeq ($(use_debug),yes)
  CXXFLAGS += -g -Weverything
  CPP_DEFS += -DDEBUG
  LDFLAGS += -fsanitize=undefined -fno-sanitize-recover
else
  CPP_DEFS += -DNDEBUG
endif

ifeq ($(use_cpp11),yes)
  CXXFLAGS += -std=gnu++11
endif

ifeq ($(use_python),yes)
  CXXFLAGS += -fPIC
endif

#CXXFLAGS += -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64
CXXFLAGS += $(OPTFLAGS)

LD = $(CXX)
LDLIBS = -lm

CXXDEP = $(CXX)
CXXDEPFLAGS = -MM

AR = ar r
RANLIB = ranlib

TARGET_MACHINE = PC_CLANG

#####  Intel C++  ##############################################

else ifeq ($(strip $(target)),PC_INTEL)

ifeq ($(use_mpi),yes)
  #CXX = mpic++
  CXX = mpiicpc
  CXXFLAGS += -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX
else
  CXX = icpc
endif

ifeq ($(use_thread),omp)
  CXXFLAGS += -qopenmp  # -openmp
  LDFLAGS  += -qopenmp  # -openmp
endif

ifeq ($(use_opt_level),low)
  OPTFLAGS = -O0
endif
ifeq ($(use_opt_level),std)
  OPTFLAGS = -O2
endif
ifeq ($(use_opt_level),high)
  #OPTFLAGS = -fast
  OPTFLAGS = -O3 -no-prec-div -xHost  # -ipo
endif

ifeq ($(use_static),yes)
  LDFLAGS  += -static-intel  # -static
  ## for intel MPI library
  # LDFLAGS  += -static_mpi
else
  CXXFLAGS += -mcmodel=medium  # -shared-intel
  LDFLAGS  += -mcmodel=medium  # -shared-intel
endif

ifeq ($(use_debug),yes)
  CXXFLAGS += -g -debug all -traceback
  CPP_DEFS += -DDEBUG
else
  CPP_DEFS += -DNDEBUG
endif

ifeq ($(use_cpp11),yes)
  CXXFLAGS += -std=c++11
  ## NB. cpp11+c99 causes compile errors for icpc 15.0.3
  #  $(error Incompatible options: cpp11 conflicts c99 for INTEL)
  $(warning cpp11 may fail for PC_INTEL, try cpp11=no)
endif

CXXFLAGS += $(OPTFLAGS)

LD = $(CXX)
LDLIBS = -lm

CXXDEP = $(CXX)
CXXDEPFLAGS = -MM

AR = ar r
RANLIB = ranlib

TARGET_MACHINE = PC_INTEL

#####  NVIDIA HPC SDK C++  #####################################

else ifeq ($(strip $(target)),PC_NVIDIA)

CPP_DEFS += -DUSE_STRDUP

ifeq ($(use_mpi),yes)
  CXX = mpic++
else ifeq ($(use_mpi),gpu_aware)
  CXX = mpic++
else
  CXX = nvc++
endif

ifeq ($(use_thread),omp)
  CXXFLAGS += -mp
  LDFLAGS  += -mp
endif

ifeq ($(use_opt_level),low)
  OPTFLAGS = -O0
endif
ifeq ($(use_opt_level),std)
  OPTFLAGS = -O2
endif
ifeq ($(use_opt_level),high)
  OPTFLAGS = -fast -O3
endif

ifeq ($(use_static),yes)
  LDFLAGS  += -Wl,-Bstatic
endif

ifeq ($(use_debug),yes)
  CXXFLAGS += -g -Minform=inform -Mbounds -Ktrap=fp
  CPP_DEFS += -DDEBUG
else
  CPP_DEFS += -DNDEBUG
endif

ifeq ($(use_cpp11),yes)
  CXXFLAGS += -std=c++11
endif

CXXFLAGS += $(OPTFLAGS)

#for profiling
CXXFLAGS += -lineinfo -g -fno-omit-frame-pointer

LD = $(CXX)
LDLIBS = -lm

ifeq ($(use_fftw_library),yes)
  ifeq ($(use_thread),omp)
    LDLIBS += -lgomp
  endif
endif

ifeq ($(use_alt_accel),yes)

ifeq ($(use_accel_impl),openacc)
  CXXFLAGS += -acc
  LDFLAGS  += -acc
endif

# CUDA_ARCH is handed over to src/lib_alt_Accel/BridgeCUDA (nvcc), so that
# use_device_arch controls the device code of the CUDA kernels as well.
# Without this, BridgeCUDA used its own hard-coded architecture list and
# silently produced a library with no code for the actual GPU.
ifeq ($(use_device_arch),kepler)
#  CXXFLAGS += -ta=tesla:ptxinfo,cc35,maxregcount:
#  LDFLAGS  += -ta=tesla:ptxinfo,cc35,maxregcount:
  CXXFLAGS += -gpu=ptxinfo,cc35,maxregcount:
  LDFLAGS  += -gpu=ptxinfo,cc35,maxregcount:
  CUDA_ARCH = -arch=sm_35
endif
ifeq ($(use_device_arch),pascal)
  CXXFLAGS += -gpu=ptxinfo,cc60,maxregcount:
  LDFLAGS  += -gpu=ptxinfo,cc60,maxregcount:
  CUDA_ARCH = -arch=sm_60
endif
ifeq ($(use_device_arch),volta)
  CXXFLAGS += -gpu=ptxinfo,cc70,maxregcount:
  LDFLAGS  += -gpu=ptxinfo,cc70,maxregcount:
  CUDA_ARCH = -arch=sm_70
endif
ifeq ($(use_device_arch),ampere)
  CXXFLAGS += -gpu=ptxinfo,cc80,maxregcount:
  LDFLAGS  += -gpu=ptxinfo,cc80,maxregcount:
  CUDA_ARCH = -arch=sm_80
endif
ifeq ($(use_device_arch),hopper)
  CXXFLAGS += -gpu=ptxinfo,cc90,maxregcount:
  LDFLAGS  += -gpu=ptxinfo,cc90,maxregcount:
  CUDA_ARCH = -arch=sm_90
endif
ifeq ($(use_device_arch),blackwell)
  # Data center Blackwell (B100, B200), sm_100. Blackwell is split over
  # two compute capability families: 10.x here, and 12.x for gb10 below.
  # Code of one family does not run on the other one.
  CXXFLAGS += -gpu=ptxinfo,cc100,maxregcount:
  LDFLAGS  += -gpu=ptxinfo,cc100,maxregcount:
  CUDA_ARCH = -arch=sm_100
endif
ifeq ($(use_device_arch),gb10)
  # GB10 (Grace-Blackwell superchip as in NVIDIA DGX Spark), sm_121.
  # Note that sm_100 (blackwell) code does NOT run on it: cubins are not
  # compatible across major compute capability versions.
  CXXFLAGS += -gpu=ptxinfo,cc121,maxregcount:
  LDFLAGS  += -gpu=ptxinfo,cc121,maxregcount:
  CUDA_ARCH = -arch=sm_121
endif
ifeq ($(use_device_arch),auto)
  # Use nvidia-smi to detect the GPU architecture.
  GPU_CC := $(shell nvidia-smi --query-gpu=compute_cap --format=csv,noheader,nounits 2>/dev/null | head -n 1)
  # Without a visible GPU neither nvc++ nor nvcc can determine the target:
  # nvcc would silently fall back to its default architecture (-arch=native
  # only warns) and produce a library that does not run on this machine.
  # Fail here instead of shipping such a library.
  ifeq ($(strip $(GPU_CC)),)
    $(error use_device_arch = auto: no GPU detected. \
      nvidia-smi is unavailable or reports no device. \
      Set use_device_arch to an explicit architecture such as gb10 or hopper, \
      or to all-major when building on a node without a GPU.)
  endif
  GPU_ARCH := cc$(subst .,,$(GPU_CC))
  CXXFLAGS += -gpu=ptxinfo,$(GPU_ARCH),maxregcount:
  LDFLAGS  += -gpu=ptxinfo,$(GPU_ARCH),maxregcount:
  # nvcc detects the GPU by itself, no need to enumerate architectures.
  CUDA_ARCH = -arch=native
endif
ifeq ($(use_device_arch),all-major)
  # Fat binary covering every major architecture of the toolkit, without
  # PTX: only GPUs known at build time are supported, an unknown one is
  # reported by BridgeACC::check_device_code() instead of being served by
  # a JIT fallback. Combine with jit if that fallback is wanted.
  # The list is derived from nvcc itself, so it needs no maintenance when
  # a new GPU generation appears. For each major compute capability the
  # lowest supported minor is taken, because a cubin runs on the same
  # major version with an equal or higher minor (sm_120 runs on sm_121).
  CUDA_ARCH = $(shell nvcc --list-gpu-code | sed 's/sm_//' | \
    awk '{maj=substr($$0,1,length($$0)-1); min=substr($$0,length($$0))+0; \
          if(!(maj in m) || min<m[maj]) m[maj]=min} \
         END{for(k in m) printf "-gencode=arch=compute_%s%s,code=sm_%s%s ", \
             k,m[k],k,m[k]}')
  CXXFLAGS += -gpu=ptxinfo,ccall-major,maxregcount:
  LDFLAGS  += -gpu=ptxinfo,ccall-major,maxregcount:
endif
ifeq ($(use_device_arch),jit)
  # PTX only: no cubin at all, the driver JIT-compiles at start-up. Runs
  # on every GPU supported by the driver, including ones not yet known to
  # this toolkit, at the price of the JIT time on the first launch.
  # The baseline is the oldest architecture of the toolkit, since PTX is
  # forward compatible only.
  CUDA_JIT_ARCH := $(shell nvcc --list-gpu-arch | head -n 1)
  CUDA_ARCH = -arch=$(CUDA_JIT_ARCH)
  CXXFLAGS += -gpu=ptxinfo,$(subst compute_,cc,$(CUDA_JIT_ARCH)),maxregcount:
  LDFLAGS  += -gpu=ptxinfo,$(subst compute_,cc,$(CUDA_JIT_ARCH)),maxregcount:
endif

# Detect a value that matched none of the branches above. Checking the
# result instead of listing the accepted names keeps this check correct
# when architectures are added or renamed.
ifeq ($(strip $(CUDA_ARCH)),)
  $(error use_device_arch = '$(use_device_arch)' is not supported. \
    See the use_device_arch branches in Makefile_target.inc \
    for the accepted values.)
endif

ifeq ($(use_accel_impl),cuda)
  # Detect the CUDA_PATH
  #CUDA_PATH_DESU ?= $(CUDA_PATH)

  # Ensure CUDA_PATH is set
  #ifeq ($(CUDA_PATH_DESU),)
  #  $(error "CUDA_PATH is not detected, please set it manually. ")
  #endif

  # Set include and library paths
  #CXXFLAGS += -I$(CUDA_PATH_DESU)/include -cuda
  #LDFLAGS  += -L$(CUDA_PATH_DESU)/lib64 -L$(CUDA_PATH_DESU)/lib -cuda
  CXXFLAGS +=  -cuda
  LDFLAGS  +=  -cuda

endif

endif

CXXDEP = $(CXX)
CXXDEPFLAGS = -MM

AR = ar r
RANLIB = ranlib

TARGET_MACHINE = PC_NVIDIA


#####  Fugaku_CLANG (using MPI)  ################################

else ifeq ($(strip $(target)),Fugaku_CLANG)

CPP_DEFS += -DFUJITSU_FX

ifeq ($(use_mpi),yes)
CXX = mpiFCCpx
#  CXX = mpiFCC
else ifeq ($(use_mpi),fjmpi)
CXX = mpiFCCpx
#  CXX = mpiFCC
else
CXX = FCCpx
#  CXX = FCC
endif

CXXFLAGS = -Nclang

ifeq ($(use_thread),omp)
  CXXFLAGS += -fopenmp
endif

ifeq ($(use_opt_level),low)
  OPTFLAGS = -O0
endif
ifeq ($(use_opt_level),std)
  OPTFLAGS = -O2
endif
ifeq ($(use_opt_level),high)
  OPTFLAGS = -Ofast -mllvm -inline-threshold=1000
endif

ifeq ($(use_debug),yes)
  CXXFLAGS += -g -Knolargepage
  CPP_DEFS += -DDEBUG
else
  CPP_DEFS += -DNDEBUG
endif

ifeq ($(use_profiler),fapp)
  CPP_DEFS += -DUSE_FAPP
endif

ifeq ($(use_qxs_arch),acle)
  CPP_DEFS += -DUSE_QXS_ACLE
endif

## NB. Fugaku_CLANG employs c++14 as default
#ifeq ($(use_cpp11),yes)
#  CXXFLAGS += -std=c++11
#endif

CXXFLAGS += $(OPTFLAGS)

LD       = $(CXX)
LDFLAGS  = $(CXXFLAGS)
LDLIBS   = -lm

ifeq ($(use_qws_library),yes)
  LDLIBS += -ltofucom -L/opt/FJSVtcs/pwrm/aarch64/lib64 -lpwr
endif

CXXDEP = $(CXX)
CXXDEPFLAGS = -MM

AR = ar r
RANLIB = ranlib

TARGET_MACHINE = Fugaku_FCC

#####  NEC SX-AURORA @KEK  #############################################

else ifeq ($(strip $(target)),SX_AURORA)

CPP_DEFS += -DNECSX
CPP_DEFS += -DUSE_STRDUP
CPP_DEFS += -DNOEXCEPTION

ifeq ($(use_mpi),yes)
  CXX = mpinc++
else
  CXX = nc++
endif

ifeq ($(use_thread),omp)
  CXXFLAGS += -fopenmp
  LDFLAGS  += -fopenmp
endif

ifeq ($(use_opt_level),low)
  OPTFLAGS = -O0 -report-all
endif
ifeq ($(use_opt_level),std)
  OPTFLAGS = -O2 -report-all
endif
ifeq ($(use_opt_level),high)
  OPTFLAGS = -O4 -report-all
endif
ifeq ($(use_opt_level),trace)
  OPTFLAGS = -O4
  CXXFLAGS += -report-all -ftrace -proginf
  LDFLAGS  += -ftrace
endif

ifeq ($(use_static),yes)
  LDFLAGS  += -static
endif

ifeq ($(use_debug),yes)
  CXXFLAGS += -g -traceback
  CPP_DEFS += -DDEBUG
else
  CPP_DEFS += -DNDEBUG
endif

ifeq ($(use_cpp11),yes)
  CXXFLAGS +=
endif

# no exceptions
CXXFLAGS += -fno-exceptions

# use restrict keyword
CXXFLAGS += -frestrict


#CXXFLAGS += -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64                      
CXXFLAGS += $(OPTFLAGS)

LD = $(CXX)
LDLIBS = -lm

CXXDEP = $(CXX)
CXXDEPFLAGS = -M

AR = ar r
RANLIB = ranlib

TARGET_MACHINE = SX_AURORA

#####
else

$(error No such a target: $(target))

endif
