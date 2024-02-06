
file(GLOB sources *.h)

add_executable(serghei ${sources} driver.cpp)

target_include_directories(serghei PUBLIC /Users/zli/Codes/PnetCDF/include)
target_link_libraries(serghei PUBLIC /Users/zli/Codes/PnetCDF/lib/libpnetcdf.a)

IF(Enable_GPU)
	if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  		set(CMAKE_CUDA_ARCHITECTURES 80)
	endif()
ELSE()
	find_package(OpenMP REQUIRED)
	target_link_libraries(serghei PUBLIC OpenMP::OpenMP_CXX)
ENDIF()

find_package(MPI REQUIRED)
target_link_libraries(serghei PUBLIC MPI::MPI_CXX)

find_package(KokkosKernels REQUIRED)
target_link_libraries(serghei PUBLIC Kokkos::kokkoskernels)

add_custom_command(TARGET serghei POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:serghei> ../bin/serghei
    COMMENT "Executable in ../bin/serghei"
)
