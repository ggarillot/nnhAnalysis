include_directories( ${CMAKE_CURRENT_LIST_DIR}/include )

file(GLOB_RECURSE sources ${CMAKE_CURRENT_LIST_DIR}/src/*.cc)
file(GLOB_RECURSE headers ${CMAKE_CURRENT_LIST_DIR}/include/*.hh)

#add library
add_library(nnhProcessor SHARED ${sources} ${headers} )
target_link_libraries( nnhProcessor ${ROOT_COMPONENT_LIBRARIES} )
install( TARGETS nnhProcessor LIBRARY DESTINATION ${CMAKE_CURRENT_LIST_DIR}/lib )