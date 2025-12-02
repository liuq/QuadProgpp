# Examples
add_executable(simple_example simple_example.cc)
target_link_libraries(simple_example PRIVATE QuadProg::quadprogpp)

# Optional: semi-definite example if Boland extension is enabled
if(QUADPROGPP_ENABLE_SEMIDEFINITE AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/semidefinite_example.cc)
    add_executable(semidefinite_example semidefinite_example.cc)
    target_link_libraries(semidefinite_example PRIVATE QuadProg::quadprogpp)
endif()

# Install examples
install(TARGETS simple_example
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
