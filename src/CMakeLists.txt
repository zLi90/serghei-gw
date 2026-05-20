file(GLOB LIB_HEADERS
	${CMAKE_CURRENT_SOURCE_DIR}/*.h
	${CMAKE_CURRENT_SOURCE_DIR}/*.hpp
	#${CMAKE_CURRENT_SOURCE_DIR}/driver.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp
)


add_executable(${CMAKE_PROJECT_NAME} ${LIB_HEADERS} driver.cpp)
#add_executable(dambreak ${LIB_HEADERS} tests/dambreak.cpp)

#	If use GPU(CUDA), specify the CUDA architecture
#	Otherwise, link OpenMP
IF(NOT Enable_GPU)
	find_package(OpenMP REQUIRED)
	target_link_libraries(serghei PUBLIC OpenMP::OpenMP_CXX)
ENDIF()

#	Link MPI
find_package(MPI REQUIRED)
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE MPI::MPI_CXX)

#	Link PNetCDF library
target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE ${PnetCDF_INCLUDE_DIRS})
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE     ${PnetCDF_LINK_LIBRARIES})

#	Link KokkosKernels
if(SERGHEI_SUBSURFACE_MODEL)
	find_package(KokkosKernels REQUIRED)
	target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE     Kokkos::kokkoskernels)
else()
	#	Link Kokkos
	find_package(Kokkos REQUIRED)
	target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE     Kokkos::kokkos)
endif()

#	Move the executable into bin/
add_custom_command(TARGET serghei POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:serghei> ../bin/serghei
    COMMENT "SERGHEI executable in ../bin/serghei"
)

if (Build_Tests)
	file( GLOB TEST_SOURCES tests/* )
	foreach( testsourcefile ${TEST_SOURCES} )
    	# Cut off the file extension and directory path
    	get_filename_component( testname ${testsourcefile} NAME_WE )
    	add_executable( ${testname} ${testsourcefile} )
		target_link_libraries(${testname} PRIVATE MPI::MPI_CXX)
		target_include_directories(${testname} PRIVATE ${PnetCDF_INCLUDE_DIRS})
		target_link_libraries(${testname} PRIVATE ${PnetCDF_LINK_LIBRARIES})
		target_link_libraries(${testname} PRIVATE     Kokkos::kokkos)
		if(SERGHEI_SUBSURFACE_MODEL)
			target_link_libraries(${testname} PRIVATE     Kokkos::kokkoskernels)
		endif()

		add_custom_command(TARGET ${testname} POST_BUILD
    		COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${testname}> ../bin/tests/${testname}
    		COMMENT "Test executable ../bin/tests/${testname}"
		)
	endforeach( testsourcefile ${TEST_SOURCES} )
endif()



