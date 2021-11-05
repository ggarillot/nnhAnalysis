########################################################
# cmake file for building Marlin example Package
# @author Jan Engels, DESY
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################

# project name
PROJECT( nnhAnalysis )

# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 1 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

LIST( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake )

### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES REQUIRED )
# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )

FIND_PACKAGE( Eigen3 REQUIRED )
INCLUDE_DIRECTORIES( SYSTEM ${EIGEN3_INCLUDE_DIR} )
LIST( APPEND ${PROJECT_NAME}_DEPENDS_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR} )


FIND_PACKAGE( Marlin 1.0 REQUIRED )
FIND_PACKAGE( FastJet REQUIRED )
FIND_PACKAGE( ROOT 6.18 REQUIRED COMPONENTS ROOTDataFrame ROOTVecOps)

FOREACH( pkg Marlin FastJet ROOT )
    IF( ${pkg}_FOUND )
		INCLUDE_DIRECTORIES( SYSTEM ${${pkg}_INCLUDE_DIRS} )
        LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
        ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
    ENDIF()
ENDFOREACH()

### LIBRARY #################################################################

# include directories
include_directories( ${PROJECT_SOURCE_DIR}/include )
#INSTALL_DIRECTORY( ./include DESTINATION . FILES_MATCHING PATTERN "*.h" )

# definitions to pass to the compiler
ADD_DEFINITIONS( "-Wall -ansi -pedantic -std=c++17 -Weffc++" )
ADD_DEFINITIONS( "-Wno-long-long" )
ADD_DEFINITIONS( "-O3" )

file(GLOB_RECURSE sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB_RECURSE headers ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/include/*.hpp)

#add library
add_library( ${PROJECT_NAME} SHARED ${sources} ${headers} )
target_link_libraries( ${PROJECT_NAME} ${ROOT_COMPONENT_LIBRARIES} )
install( TARGETS ${PROJECT_NAME} LIBRARY DESTINATION ${PROJECT_SOURCE_DIR}/lib )

file( GLOB APP_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/exec/*.cxx )
foreach( testsourcefile ${APP_SOURCES} )
	get_filename_component(testname ${testsourcefile} NAME_WE)
	add_executable( ${testname} ${testsourcefile} )
	target_link_libraries( ${testname} ${PROJECT_NAME} ${ROOT_LIBRARIES} )
	install(TARGETS ${testname} RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin)
endforeach( testsourcefile ${APP_SOURCES} )