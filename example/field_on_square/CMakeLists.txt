add_executable(osh_field_on_square main.cpp)
target_link_libraries(osh_field_on_square PRIVATE omega_h)
add_test(NAME osh_field_on_square COMMAND $<TARGET_FILE:osh_field_on_square>)

add_custom_command(OUTPUT square.msh
	COMMAND ${GMSH_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/square.geo
	-2 -o square.msh
	MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/square.geo
)
add_custom_target(osh_field_on_square_mesh ALL DEPENDS square.msh)
