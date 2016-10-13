file(GLOB sources *.cc)
file(GLOB headers *.hh)
set(library_sources ${sources})
list(REMOVE_ITEM library_sources ${CMAKE_CURRENT_SOURCE_DIR}/main.cc)
add_library(quadprog STATIC ${library_sources} ${headers})
set_property(TARGET quadprog PROPERTY POSITION_INDEPENDENT_CODE 1)

add_executable(main main.cc ${headers})
target_link_libraries(main quadprog)