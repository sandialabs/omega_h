add_executable(osh_castle main.cpp)
target_link_libraries(osh_castle PRIVATE omega_h)
set_property(TARGET osh_castle PROPERTY CXX_STANDARD "14")
set_property(TARGET osh_castle PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(TARGET osh_castle PROPERTY CXX_EXTENSIONS OFF)
add_test(NAME osh_castle COMMAND  COMMAND $<TARGET_FILE:osh_castle>)

add_custom_command(OUTPUT square.msh
	COMMAND ${GMSH_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/square.geo
	-2 -o square.msh
	MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/square.geo
)
add_custom_target(osh_castle_mesh ALL DEPENDS square.msh)
