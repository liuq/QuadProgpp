# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

file(GLOB sources *.cc)
file(GLOB headers *.hh)
set(library_sources ${sources})
list(REMOVE_ITEM library_sources ${CMAKE_CURRENT_SOURCE_DIR}/main.cc)
add_library(quadprog STATIC ${library_sources} ${headers})
set_property(TARGET quadprog PROPERTY POSITION_INDEPENDENT_CODE 1)
set_property(TARGET quadprog PROPERTY PUBLIC_HEADER ${headers})

add_executable(main main.cc ${headers})
target_link_libraries(main quadprog)

install(TARGETS quadprog
  EXPORT quadprog-targets
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/QuadProg++
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install(EXPORT quadprog-targets DESTINATION cmake)
