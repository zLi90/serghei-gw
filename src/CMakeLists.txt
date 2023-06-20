
file(GLOB sources *.h)

add_executable(serghei ${sources} driver.cpp)

target_include_directories(serghei PUBLIC /Users/zli/Codes/PnetCDF/include)
target_link_libraries(serghei PUBLIC /Users/zli/Codes/PnetCDF/lib/libpnetcdf.a)

message(STATUS "set compiler")
set(CMAKE_C_COMPILER "/usr/local/bin/gcc-11")
set(CMAKE_CXX_COMPILER "/usr/local/bin/g++-11")

IF(Enable_GPU)
	if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  		set(CMAKE_CUDA_ARCHITECTURES 80)
	endif()
ELSE()
	set(OPENMP_LIBRARIES "/usr/local/Cellar/gcc/11.2.0/lib")
         set(OPENMP_INCLUDES "/usr/local/Cellar/gcc/11.2.0/include")
	find_package(OpenMP REQUIRED)
	include_directories("${OPENMP_INCLUDES}")
	link_directories("${OPENMP_LIBRARIES}")
ENDIF()

find_package(MPI REQUIRED)
target_link_libraries(serghei PUBLIC MPI::MPI_CXX)

find_package(Kokkos REQUIRED)
target_link_libraries(serghei PUBLIC Kokkos::kokkos)

add_custom_command(TARGET serghei POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:serghei> ../bin/serghei
    COMMENT "Executable in ../bin/serghei"
)
