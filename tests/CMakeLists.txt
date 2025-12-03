# Fetch Catch2
include(FetchContent)
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.5.0
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(Catch2)

# List all test files
set(TEST_SOURCES
    test_matrices.cc
    test_triangular.cc
    test_solvers.cc
    test_quadprog.cc
    test_quadprog_stress.cc
)

# Create test executable
add_executable(quadprogpp_tests ${TEST_SOURCES})

# Link against the header-only library
target_link_libraries(quadprogpp_tests 
    PRIVATE 
        QuadProg::quadprogpp
        Catch2::Catch2WithMain
)

# If spdlog is enabled, we need to link it for the tests too
# because the headers will include spdlog headers
if(QUADPROGPP_USE_SPDLOG)
    target_link_libraries(quadprogpp_tests PRIVATE spdlog::spdlog)
endif()

# Set C++ standard (should inherit from quadprogpp, but be explicit)
target_compile_features(quadprogpp_tests PRIVATE cxx_std_20)

# Enable warnings
if(MSVC)
    target_compile_options(quadprogpp_tests PRIVATE /W4)
else()
    target_compile_options(quadprogpp_tests PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Discover and register tests with CTest
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(CTest)
include(Catch)
catch_discover_tests(quadprogpp_tests
    EXTRA_ARGS --order rand  # Randomize test order to catch dependencies
)

# Optional: Add custom test targets for convenience
add_custom_target(test_verbose
    COMMAND ${CMAKE_CTEST_COMMAND} --verbose --output-on-failure
    DEPENDS quadprogpp_tests
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_custom_target(test_fast
    COMMAND quadprogpp_tests --order rand --durations yes
    DEPENDS quadprogpp_tests
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# Separate target for stress tests (since they take longer)
add_custom_target(test_stress
    COMMAND quadprogpp_tests "[stress]" --order rand
    DEPENDS quadprogpp_tests
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# Benchmark tests (separate since they're slow)
add_custom_target(test_benchmark
    COMMAND quadprogpp_tests "[benchmark]" --benchmark-samples 100
    DEPENDS quadprogpp_tests
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)