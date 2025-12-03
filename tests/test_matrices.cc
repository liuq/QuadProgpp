#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <quadprog/quadprog.h>

#if defined(QUADPROGPP_MATRIX_BACKEND_BUILTIN)

using namespace quadprog;
using Catch::Matchers::WithinAbs;

TEST_CASE("Matrix basic operations", "[matrix]") {
    SECTION("Construction and access") {
        Matrix<double> M(3, 2);
        REQUIRE(M.rows() == 3);
        REQUIRE(M.cols() == 2);
        
        M(1, 1) = 5.0;
        REQUIRE(M(1, 1) == 5.0);
    }
    
    SECTION("Resize") {
        Matrix<double> M(2, 2);
        M.resize(4, 3);
        REQUIRE(M.rows() == 4);
        REQUIRE(M.cols() == 3);
    }
    
    SECTION("Bounds checking in debug") {
#ifndef QUADPROGPP_BOUNDS_CHECK
    SKIP("Bounds checking not enabled");
#endif
        Matrix<double> M(2, 2);
        REQUIRE_THROWS_AS(M(3, 0), std::out_of_range);
        REQUIRE_THROWS_AS(M(0, 5), std::out_of_range);
    }
}

TEST_CASE("Vector operations", "[vector]") {
    SECTION("Construction and access") {
        Vector<double> v(5);
        REQUIRE(v.size() == 5);
        
        v(2) = 3.14;
        REQUIRE_THAT(v(2), WithinAbs(3.14, 1e-10));
    }
    
    SECTION("Initialization with value") {
        Vector<double> v(3, 2.5);
        REQUIRE_THAT(v(0), WithinAbs(2.5, 1e-10));
        REQUIRE_THAT(v(1), WithinAbs(2.5, 1e-10));
        REQUIRE_THAT(v(2), WithinAbs(2.5, 1e-10));
    }
}

#endif  // QUADPROGPP_MATRIX_BACKEND_BUILTIN