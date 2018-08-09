add_executable(osh_castle main.cpp)
target_link_libraries(osh_castle PRIVATE omega_h)
add_test(NAME osh_castle COMMAND  COMMAND $<TARGET_FILE:osh_castle>)

add_custom_command(OUTPUT square.msh
	COMMAND ${GMSH_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/square.geo
	-2 -o square.msh
	MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/square.geo
)
add_custom_target(osh_castle_mesh ALL DEPENDS square.msh)
