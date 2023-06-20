cmake_minimum_required (VERSION 3.10)

IF(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.12.0")
  MESSAGE(STATUS "Setting policy CMP0074 to use <Package>_ROOT variables")
  CMAKE_POLICY(SET CMP0074 NEW)
ENDIF()

project (serghei)

add_subdirectory(src)


#----------------------------------------------------
# Execution with OpenMP:
# cmake -DEnable_GPU=False -DKokkos_ROOT=${path-to-kokkos-build(OpenMP)}
#----------------------------------------------------
# Execution with CUDA:
# cmake -DEnable_GPU=True -DKokkos_ROOT=${path-to-kokkos-build(CUDA)}
#----------------------------------------------------
