#   This is the CMakeLists file for building and compiling SERGHEI-GW
#   To use, open a terminal, cd to the SERGHEI-GW directory and type "cmake ."
#   Then, type "make" to compile SERGHEI-GW
#   Note: If you changed any cmake options herein, before rebuilding, delete "CMakeCache.txt" to erase old cmake settings

cmake_minimum_required (VERSION 3.14)

#   If needed, specify the CXX compiler here, otherwise comment out these lines
set(CMAKE_C_COMPILER gcc-14)
set(CMAKE_CXX_COMPILER g++-14)

# Name the project
project (serghei LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

#   Specify the path to PnetCDF
set(PNETCDF_ROOT "/Users/zhili/Codes/pnetcdf")

#   Specify the path where KokkosKernels is built
set(KokkosKernels_ROOT "/Users/zhili/Codes/kokkoskernels-install")

#   Specify if using GPU or not
set(Enable_GPU "OFF")

#   Specify SERGHEI Macros (You can add more macros below)
add_compile_definitions(SERGHEI_FRICTION_MODEL=1)
execute_process(COMMAND git describe COMMAND tr -d '\n' OUTPUT_VARIABLE serghei_git_version)
add_compile_definitions(SERGHEI_GIT_VERSION="${serghei_git_version}")

option(SERGHEI_SUBSURFACE_MODEL ON "Enable subsurface model")
option(SERGHEI_SWE_MODEL OFF "Enable surface flow model")
option(SERGHEI_SUBSURFACE_TRANSPORT ON "Enable subsurface transport model")
option(SERGHEI_SURFACE_TRANSPORT ON "Enable surface transport model")
option(CROP_GROWTH_MODEL ON "Enable crop growth model")
option(SERGHEI_PARTICLE_TRACKING OFF "Enable particle tracking")
option(SERGHEI_VEGETATION_MODEL OFF "Enable vegetation model")
option(SERGHEI_WAVE_MODEL OFF "Enable wave-driven RE boundary module")
option(SERGHEI_USE_PNETCDF ON "Enable PNETCDF I/O support")
option(SERGHEI_TOOLS ON "Enable additional tools")
option(SERGHEI_FRICTION_MODEL ON "Enable friction model")

if(SERGHEI_SWE_MODEL AND SERGHEI_SUBSURFACE_MODEL)
  set(SERGHEI_SWE_RE_VALUE 1)
else()
  set(SERGHEI_SWE_RE_VALUE 0)
endif()

if(SERGHEI_SUBSURFACE_MODEL)
  set(SERGHEI_SUBSURFACE_TRANSPORT_VALUE $<BOOL:${SERGHEI_SUBSURFACE_TRANSPORT}>)
else()
  set(SERGHEI_SUBSURFACE_TRANSPORT_VALUE 0)
  if(SERGHEI_SUBSURFACE_TRANSPORT)
    message(WARNING "SERGHEI_SUBSURFACE_TRANSPORT requires SERGHEI_SUBSURFACE_MODEL=ON. Disabling subsurface transport.")
  endif()
endif()

if(SERGHEI_SWE_MODEL)
  set(SERGHEI_SURFACE_TRANSPORT_VALUE $<BOOL:${SERGHEI_SURFACE_TRANSPORT}>)
else()
  set(SERGHEI_SURFACE_TRANSPORT_VALUE 0)
  if(SERGHEI_SURFACE_TRANSPORT)
    message(WARNING "SERGHEI_SURFACE_TRANSPORT requires SERGHEI_SWE_MODEL=ON. Disabling surface transport.")
  endif()
endif()

add_compile_definitions(SERGHEI_RE_MODEL=$<BOOL:${SERGHEI_SUBSURFACE_MODEL}>)
add_compile_definitions(SERGHEI_SWE_MODEL=$<BOOL:${SERGHEI_SWE_MODEL}>)
add_compile_definitions(SERGHEI_SWE_RE=${SERGHEI_SWE_RE_VALUE})
add_compile_definitions(SERGHEI_SUBSURFACE_TRANSPORT=${SERGHEI_SUBSURFACE_TRANSPORT_VALUE})
add_compile_definitions(SERGHEI_SURFACE_TRANSPORT=${SERGHEI_SURFACE_TRANSPORT_VALUE})
add_compile_definitions(CROP_GROWTH_MODEL=$<BOOL:${CROP_GROWTH_MODEL}>)
add_compile_definitions(SERGHEI_LPT=$<BOOL:${SERGHEI_PARTICLE_TRACKING}>)
add_compile_definitions(SERGHEI_VEGETATION_MODEL=$<BOOL:${SERGHEI_VEGETATION_MODEL}>)
add_compile_definitions(SERGHEI_WAVE_MODEL=$<BOOL:${SERGHEI_WAVE_MODEL}>)
add_compile_definitions(SERGHEI_USE_PNETCDF=$<BOOL:${SERGHEI_USE_PNETCDF}>)

# debugging options
option(SERGHEI_DEBUG_PARALLEL_DECOMPOSITION OFF "Enable parallel decomposition debugging")
option(SERGHEI_DEBUG_WORKFLOW OFF "Enable workflow debugging")
option(SERGHEI_DEBUG_KOKKOS_SETUP OFF "Enable Kokkos setup debugging")
option(SERGHEI_DEBUG_BOUNDARY OFF "Enable boundary debugging")
option(SERGHEI_DEBUG_DT OFF "Enable timestep debugging")
option(SERGHEI_DEBUG_TOOLS OFF "Enable tools debugging")
option(SERGHEI_DEBUG_MASS_CONS OFF "Enable mass conservation debugging")
option(SERGHEI_DEBUG_INFILTRATION OFF "Enable infiltration debugging")
option(SERGHEI_DEBUG_MPI OFF "Enable MPI debugging")
option(SERGHEI_DEBUG_OUTPUT OFF "Enable output debugging")

# just to print out the compile defintions
get_directory_property( DirDefs DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_DEFINITIONS )
message(STATUS "SERGHEI Compilation definitions")
foreach( d ${DirDefs} )
    message( STATUS "\t" ${d})
endforeach()

#   Add the subdirectory where the source codes are stored
add_subdirectory(src)

