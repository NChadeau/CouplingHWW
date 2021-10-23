########################################################
# cmake file for building Marlin example Package
# @author Jan Engels, Desy IT
CMAKE_MINIMUM_REQUIRED(VERSION 3.6 FATAL_ERROR)
########################################################

# include directories
INCLUDE_DIRECTORIES( ./include )

# add library
AUX_SOURCE_DIRECTORY( ./src library_sources )


file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

file(GLOB headers ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp)

add_library( PlotsLib ${sources} ${headers} )


file( GLOB APP_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/exec/*.cxx )
foreach( testsourcefile ${APP_SOURCES} )
	get_filename_component(testname ${testsourcefile} NAME_WE)
	add_executable( ${testname} ${testsourcefile} )
	target_link_libraries( ${testname} PlotsLib ${ROOT_COMPONENT_LIBRARIES} )
	install(TARGETS ${testname} RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin)
endforeach( testsourcefile ${APP_SOURCES} )

