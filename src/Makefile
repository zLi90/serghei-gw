#####################################################
## Makefile to build SERGHEI
## Copyright SERGHEI Workgroup 2019-2022
#####################################################
## USAGE
## Invoke make with your prefered architecture, for example
## make arch=cpu
## make arch=gpu device=Volta70
## If you invoke simply "make", cpu is the default architecture
## For GPU compilation device is mandatory

.PHONY: main

## DEFINE MODEL SETUP
## In this section you can control the compilation of model components
## The defaults settings are hardcoded in SERGHEI
## By uncommenting and defining these flags you can manually control which components are compiled
## You can also pass these definitions from command line, e.g., make SERGHEI_PARTICLE_TRACKING=1
#SERGHEI_SUBSURFACE_MODEL=0
#SERGHEI_PARTICLE_TRACKING=0
#SERGHEI_VEGETATION_MODEL=0
#SERGHEI_TOOLS=0


## Define library paths
PNETCDF_PATH=${HOME}/PnetCDF
KOKKOS_PATH = ${HOME}/kokkos
KOKKOS_SRC_PATH = ${KOKKOS_PATH}

## by default build for cpu, alternative should be gpu
ifndef ($arch)
	arch = cpu
	archCheck = checkCPU
endif

KOKKOS_DEBUG=no
ifeq ($(debug),yes)
	KOKKOS_DEBUG=yes
endif



########### LEAVE THIS ALONE ##############
# Define binary folders and build executable name
BIN_DIR = ../bin
EXE = $(BIN_DIR)/serghei
dummy_build_folder := $(shell mkdir -p $(BIN_DIR))


CXXFLAGS = -O3 -I${PNETCDF_PATH}/include
MPICXX=`which mpic++`

ifdef SERGHEI_SUBSURFACE_MODEL
	SERGHEI_MODEL_FLAGS += -DSERGHEI_SUBSURFACE_MODEL=$(SERGHEI_SUBSURFACE_MODEL)
endif

ifdef SERGHEI_PARTICLE_TRACKING
	SERGHEI_MODEL_FLAGS += -DSERGHEI_PARTICLE_TRACKING=$(SERGHEI_PARTICLE_TRACKING)
endif
ifdef SERGHEI_VEGETATION_MODEL
	SERGHEI_MODEL_FLAGS += -DSERGHEI_VEGETATION_MODEL=$(SERGHEI_VEGETATION_MODEL)
endif
ifdef SERGHEI_TOOLS
	SERGHEI_MODEL_FLAGS += -DSERGHEI_TOOLS=$(SERGHEI_TOOLS)
endif
ifdef DEBUG_PARALLEL_DECOMPOSITION
	SERGHEI_MODEL_FLAGS += -DDEBUG_PARALLEL_DECOMPOSITION=$(DEBUG_PARALLEL_DECOMPOSITION)
endif
ifdef DEBUG_WORKFLOW
	SERGHEI_MODEL_FLAGS += -DDEBUG_WORKFLOW=$(DEBUG_WORKFLOW)
endif
ifdef DEBUG_KOKKOS_SETUP
	SERGHEI_MODEL_FLAGS += -DDEBUG_KOKKOS_SETUP=$(DEBUG_KOKKOS_SETUP)
endif
ifdef DEBUG_BOUNDARY
	SERGHEI_MODEL_FLAGS += -DDEBUG_BOUNDARY=$(DEBUG_BOUNDARY)
endif
ifdef FRICTION_MODEL
	SERGHEI_MODEL_FLAGS += -DFRICTION_MODEL=$(FRICTION_MODEL)
endif
ifdef DEBUG_DT
	SERGHEI_MODEL_FLAGS += -DDEBUG_DT=$(DEBUG_DT)
endif
ifdef SERGHEI_DEBUG_TOOLS
	SERGHEI_MODEL_FLAGS += -DSERGHEI_DEBUG_TOOLS=$(SERGHEI_DEBUG_TOOLS)
endif
ifeq ($(debug),yes)
	SERGHEI_MODEL_FLAGS += -g
endif

ifeq ($(arch),cpu)
	CXX = ${MPICXX}
	CXXFLAGS += -DKOKKOS_ENABLE_CXX17 -std=c++17
	KOKKOS_DEVICES = "OpenMP"
	KOKKOS_OPTIONS = "aggressive_vectorization,disable_profiling"
endif

ifeq ($(arch),gpu)
	archCheck=checkGPU
	CXX = ${KOKKOS_PATH}/bin/nvcc_wrapper -x cu -ccbin `which mpic++`
	KOKKOS_DEVICES = "Cuda"
	KOKKOS_OPTIONS += "disable_profiling"
	KOKKOS_CUDA_OPTIONS += "enable_lambda"
endif

CXXFLAGS += $(SERGHEI_MODEL_FLAGS)
LINK = ${CXX}

LDFLAGS = -L${PNETCDF_PATH}/lib -lpnetcdf

SRC = driver.cpp
OBJ = $(BIN_DIR)/$(notdir $(SRC:.cpp=.o))

# what does this do?
DEPFLAGS = -M

# exports all variables to be used by sub-make
export

main: $(archCheck) $(arch)

checkCPU:

checkGPU:
ifndef device
	$(error GPU architecture undefined. Please set e.g., device=Volta70. You can choose from Ampere80, Volta70, Pascal61, Maxwell50, etc... See for more options  https://github.com/kokkos/kokkos/wiki/Compiling#table-43-architecture-variables)
endif
KOKKOS_ARCH = ${device}
# For device keywords see https://github.com/kokkos/kokkos/wiki/Compiling#table-43-architecture-variables



cpu:
	@echo "\033[92mBuilding for CPU...\033[0m"
	@make -j4 -f Makefile.cpu

gpu:
	@echo "\033[92mBuilding for GPU...\033[0m"
	@make -j4 -f Makefile.gpu



clean:
	rm -f *.gch *.o *.dat serghei
	rm -rf $(BIN_DIR)

kclean: clean
	rm -f KokkosCore_config.h KokkosCore_config.tmp libkokkos.a
