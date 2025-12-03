#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <quadprog/quadprog.h>

using namespace quadprog;
using Catch::Matchers::WithinAbs;

TEST_CASE("Lower triangular matrix", "[triangular]") {
    SECTION("Construction and access") {
        LowerTriangularMatrix<double> L(3);
        
        L(0, 0) = 1.0;
        L(1, 0) = 2.0; L(1, 1) = 3.0;
        L(2, 0) = 4.0; L(2, 1) = 5.0; L(2, 2) = 6.0;
        
        REQUIRE_THAT(L(0, 0), WithinAbs(1.0, 1e-10));
        REQUIRE_THAT(L(1, 0), WithinAbs(2.0, 1e-10));
        REQUIRE_THAT(L(2, 1), WithinAbs(5.0, 1e-10));
    }
    
    SECTION("Upper part access should throw") {
#ifndef QUADPROGPP_BOUNDS_CHECK
    SKIP("Bounds checking not enabled");
#endif
        LowerTriangularMatrix<double> L(3);
        REQUIRE_THROWS_AS(L(0, 1), std::out_of_range);
        REQUIRE_THROWS_AS(L(1, 2), std::out_of_range);
    }
}

TEST_CASE("Transpose creates zero-copy view", "[triangular][transpose]") {
    LowerTriangularMatrix<double> L(3);
    L(0, 0) = 1.0;
    L(1, 0) = 2.0; L(1, 1) = 3.0;
    L(2, 0) = 4.0; L(2, 1) = 5.0; L(2, 2) = 6.0;
    
    SECTION("Upper from lower transpose") {
        auto U = transpose(L);
        
        // U(i,j) should equal L(j,i)
        REQUIRE_THAT(U(0, 1), WithinAbs(2.0, 1e-10));  // L(1,0)
        REQUIRE_THAT(U(0, 2), WithinAbs(4.0, 1e-10));  // L(2,0)
        REQUIRE_THAT(U(1, 2), WithinAbs(5.0, 1e-10));  // L(2,1)
        REQUIRE_THAT(U(0, 0), WithinAbs(1.0, 1e-10));  // L(0,0)
    }
    
    SECTION("Modifications through view affect original") {
        auto U = transpose(L);
        U(0, 1) = 99.0;  // Should modify L(1,0)
        
        REQUIRE_THAT(L(1, 0), WithinAbs(99.0, 1e-10));
        REQUIRE_THAT(U(0, 1), WithinAbs(99.0, 1e-10));
    }
    
    SECTION("Double transpose returns to original indexing") {
        auto U = transpose(L);
        auto L2 = transpose(U);
        
        REQUIRE_THAT(L2(1, 0), WithinAbs(2.0, 1e-10));
        REQUIRE_THAT(L2(2, 1), WithinAbs(5.0, 1e-10));
        
        // Modify through double transpose
        L2(2, 1) = 42.0;
        REQUIRE_THAT(L(2, 1), WithinAbs(42.0, 1e-10));
    }
}