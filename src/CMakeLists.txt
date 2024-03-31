
file(GLOB LIB_HEADERS
	${CMAKE_CURRENT_SOURCE_DIR}/*.h
	${CMAKE_CURRENT_SOURCE_DIR}/*.hpp
)

add_executable(${CMAKE_PROJECT_NAME} ${LIB_HEADERS} driver.cpp)

#	If use GPU(CUDA), specify the CUDA architecture
#	Otherwise, link OpenMP
IF(Enable_GPU)
	if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  		set(CMAKE_CUDA_ARCHITECTURES 80)
	endif()
ELSE()
	find_package(OpenMP REQUIRED)
	target_link_libraries(serghei PUBLIC OpenMP::OpenMP_CXX)
ENDIF()

#	Link MPI
find_package(MPI REQUIRED)
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE MPI::MPI_CXX)

#	Link PNetCDF library
target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE ${PNETCDF_ROOT}/include)
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE     ${PNETCDF_ROOT}/lib/libpnetcdf.a)

#	Link KokkosKernels
find_package(KokkosKernels REQUIRED)
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE     Kokkos::kokkoskernels)

#	Move the executable into bin/
add_custom_command(TARGET serghei POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:serghei> ../bin/serghei
    COMMENT "SERGHEI executable in ../bin/serghei"
)
