cmake_minimum_required(VERSION 3.14)

project(QuadProgpp 
    VERSION 2.0.0
    DESCRIPTION "Modern C++ Quadratic Programming Solver (Goldfarb-Idnani Algorithm) header-only library"
    LANGUAGES CXX
)

# C++ Standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Options
option(QUADPROGPP_USE_EIGEN "Use Eigen for matrix operations" OFF)
option(QUADPROGPP_USE_ARMADILLO "Use Armadillo for matrix operations" OFF)
option(QUADPROGPP_BUILD_TESTS "Build test suite" ON)
option(QUADPROGPP_BUILD_EXAMPLES "Build examples" ON)
option(QUADPROGPP_USE_SPDLOG "Use spdlog for logging" OFF)
option(QUADPROGPP_ENABLE_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(QUADPROGPP_ENABLE_BOUNDS_CHECKING "Enable bounds checking" ON)
option(QUADPROGPP_USE_SPDLOG "Use spdlog for logging" OFF)
option(QUADPROGPP_ENABLE_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(QUADPROGPP_ENABLE_BOUNDS_CHECKING "Enable bounds checking" ON)

# FetchContent module
include(FetchContent)

# Configure matrix backend
set(QUADPROGPP_MATRIX_BACKEND "builtin" CACHE STRING "Matrix backend to use")
set_property(CACHE QUADPROGPP_MATRIX_BACKEND PROPERTY STRINGS "builtin" "eigen" "armadillo")

# Handle matrix backend selection
if(QUADPROGPP_USE_EIGEN OR QUADPROGPP_MATRIX_BACKEND STREQUAL "eigen")
    find_package(Eigen3 3.3 QUIET)
    if(NOT Eigen3_FOUND)
        message(STATUS "Eigen3 not found locally, fetching from GitHub...")
        FetchContent_Declare(
            Eigen
            GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
            GIT_TAG 3.4.0
            GIT_SHALLOW TRUE
        )
        FetchContent_MakeAvailable(Eigen)
    endif()
    set(QUADPROGPP_HAS_EIGEN TRUE)
    set(QUADPROGPP_MATRIX_BACKEND "eigen")
elseif(QUADPROGPP_USE_ARMADILLO OR QUADPROGPP_MATRIX_BACKEND STREQUAL "armadillo")
    find_package(Armadillo QUIET)
    if(Armadillo_FOUND)
        set(QUADPROGPP_HAS_ARMADILLO TRUE)
        set(QUADPROGPP_MATRIX_BACKEND "armadillo")
    else()
        message(WARNING "Armadillo requested but not found, falling back to builtin")
        set(QUADPROGPP_MATRIX_BACKEND "builtin")
    endif()
else()
    set(QUADPROGPP_MATRIX_BACKEND "builtin")
endif()

message(STATUS "QuadProgpp matrix backend: ${QUADPROGPP_MATRIX_BACKEND}")

# Configure header
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/include/quadprog/config.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/include/quadprog/config.h
)

# Main library - INTERFACE (header-only)
add_library(quadprogpp INTERFACE)
add_library(QuadProg::quadprogpp ALIAS quadprogpp)

target_include_directories(quadprogpp
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

# Set C++ standard requirement for users
target_compile_features(quadprogpp INTERFACE cxx_std_20)

# Link matrix backend
if(QUADPROGPP_HAS_EIGEN)
    target_link_libraries(quadprogpp INTERFACE Eigen3::Eigen)
    target_compile_definitions(quadprogpp INTERFACE QUADPROGPP_USE_EIGEN)
elseif(QUADPROGPP_HAS_ARMADILLO)
    target_link_libraries(quadprogpp INTERFACE ${ARMADILLO_LIBRARIES})
    target_include_directories(quadprogpp INTERFACE ${ARMADILLO_INCLUDE_DIRS})
    target_compile_definitions(quadprogpp INTERFACE QUADPROGPP_USE_ARMADILLO)
endif()

# Bounds checking
if(QUADPROGPP_ENABLE_BOUNDS_CHECKING)
    target_compile_definitions(quadprogpp INTERFACE QUADPROGPP_BOUNDS_CHECK)
endif()

# Examples
if(QUADPROGPP_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(QUADPROGPP_USE_SPDLOG)
    find_package(spdlog QUIET)
    if(NOT spdlog_FOUND)
        message(STATUS "spdlog not found locally, fetching from GitHub...")
        FetchContent_Declare(
            spdlog
            GIT_REPOSITORY https://github.com/gabime/spdlog.git
            GIT_TAG v1.16.0
            GIT_SHALLOW TRUE
        )
        FetchContent_MakeAvailable(spdlog)
    endif()
#    target_compile_definitions(quadprogpp PRIVATE QUADPROGPP_USE_SPDLOG)
    target_link_libraries(quadprogpp PRIVATE spdlog::spdlog)
endif()

if (QUADPROGPP_ENABLE_WARNINGS_AS_ERRORS)
    if(MSVC)
        target_compile_options(quadprogpp PRIVATE /WX)
    else()
        target_compile_options(quadprogpp PRIVATE -Werror)
    endif()
endif()

# if (QUADPROGPP_BUILD_TESTS)
#     FetchContent_Declare(
#         Catch2
#         GIT_REPOSITORY https://github.com/catchorg/Catch2.git
#         GIT_TAG v3.11.0
#         GIT_SHALLOW TRUE
#     )
#     FetchContent_MakeAvailable(Catch2)

#     # Test executable
#     add_executable(quadprogpp_tests
#         tests/test_main.cc        
#         tests/test_matrices.cc
#         tests/test_triangular.cc
#         tests/test_solvers.cc
#         tests/test_quadprog.cc
#         tests/test_quadprog_stress.cc
#     )

#     target_link_libraries(quadprogpp_tests 
#         PRIVATE 
#             QuadProg::quadprogpp
#             Catch2::Catch2WithMain
#     )

#     # If logging is enabled, link spdlog
#     if(QUADPROGPP_USE_SPDLOG)
#         target_link_libraries(quadprogpp_tests PRIVATE spdlog::spdlog)
#     endif()

#     # Discover tests
#     include(CTest)
#     include(Catch)
#     catch_discover_tests(quadprogpp_tests)
# endif()

# Compiler warnings
if(MSVC)
    target_compile_options(quadprogpp INTERFACE /W4 /WX-)
else()
    target_compile_options(quadprogpp INTERFACE -Wall -Wextra -Wpedantic)
endif()

# Tests
if(QUADPROGPP_BUILD_TESTS)
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/CMakeLists.txt)
        enable_testing()
        add_subdirectory(tests)
    else()
        message(STATUS "Tests directory not yet implemented, skipping")
    endif()
endif()

# Installation
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Install header files
install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING 
        PATTERN "*.h"
        PATTERN "*.hh"
        PATTERN "*.tpp"
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/include/quadprog/config.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/quadprog
)

# Export targets - this works for INTERFACE libraries
install(TARGETS quadprogpp
    EXPORT QuadProgppTargets
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT QuadProgppTargets
    FILE QuadProgppTargets.cmake
    NAMESPACE QuadProg::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/QuadProgpp
)

# Package config
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/QuadProgppConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/QuadProgppConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/QuadProgpp
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/QuadProgppConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/QuadProgppConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/QuadProgppConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/QuadProgpp
)