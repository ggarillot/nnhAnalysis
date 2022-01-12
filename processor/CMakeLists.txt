include_directories( ${CMAKE_CURRENT_LIST_DIR}/include )

file(GLOB_RECURSE sources ${CMAKE_CURRENT_LIST_DIR}/src/*.cc)
file(GLOB_RECURSE headers ${CMAKE_CURRENT_LIST_DIR}/include/*.hh)

#add library
add_library(nnhProcessorLib SHARED ${sources} ${headers} )
target_link_libraries( nnhProcessorLib ${ROOT_COMPONENT_LIBRARIES} )
install( TARGETS nnhProcessorLib LIBRARY DESTINATION ${CMAKE_CURRENT_LIST_DIR}/lib )