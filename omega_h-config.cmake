get_filename_component(OSH_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
get_filename_component(OSH_LIB_DIR "${SELF_DIR}" PATH)
include(${OSH_LIB_DIR}/omega_h-targets.cmake)
get_filename_component(Omega_h_INCLUDE_DIRS "${OSH_LIB_DIR}/../include" ABSOLUTE)
