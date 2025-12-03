#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <quadprog/quadprog.h>
#include <cmath>

#if defined(QUADPROGPP_MATRIX_BACKEND_BUILTIN)

using namespace quadprog;
using Catch::Matchers::WithinAbs;

TEST_CASE("Forward elimination", "[solvers]") {
    LowerTriangularMatrix<double> L(3);
    L(0, 0) = 2.0;
    L(1, 0) = 1.0; L(1, 1) = 3.0;
    L(2, 0) = 4.0; L(2, 1) = 1.0; L(2, 2) = 2.0;
    
    Vector<double> b = {2.0, 8.0, 13.0};
    Vector<double> y(3);
    
    forward_elimination(L, y, b);
    
    // Verify L*y = b
    REQUIRE_THAT(L(0, 0) * y(0), 
                 WithinAbs(b(0), 1e-10));
    REQUIRE_THAT(L(1, 0) * y(0) + L(1, 1) * y(1), 
                 WithinAbs(b(1), 1e-10));
    REQUIRE_THAT(L(2, 0) * y(0) + L(2, 1) * y(1) + L(2, 2) * y(2), 
                 WithinAbs(b(2), 1e-10));
}

TEST_CASE("Backward substitution", "[solvers]") {
    UpperTriangularMatrix<double> U(3);
    // Note: This is created via transpose internally, but we set up the test
    LowerTriangularMatrix<double> L(3);
    L(0, 0) = 2.0;
    L(1, 0) = 0.0; L(1, 1) = 3.0;
    L(2, 0) = 0.0; L(2, 1) = 1.0; L(2, 2) = 2.0;
    
    auto U2 = transpose(L);
    
    Vector<double> b = {13.0, 8.0, 2.0};
    Vector<double> x(3);
    
    backward_substitution(U2, x, b);
    
    // Verify U*x = b
    double error = 0.0;
    error += std::abs(U2(0, 0) * x(0) + U2(0, 1) * x(1) + U2(0, 2) * x(2) - b(0));
    error += std::abs(U2(1, 1) * x(1) + U2(1, 2) * x(2) - b(1));
    error += std::abs(U2(2, 2) * x(2) - b(2));
    
    REQUIRE(error < 1e-10);
}

TEST_CASE("Cholesky solve", "[solvers]") {
    // Create SPD matrix A = L*L^T
    LowerTriangularMatrix<double> L(3);
    L(0, 0) = 2.0;
    L(1, 0) = 1.0; L(1, 1) = 2.0;
    L(2, 0) = 3.0; L(2, 1) = 1.0; L(2, 2) = 1.5;
    
    // Known solution
    Vector<double> x_true = {1.0, 2.0, 3.0};
    
    // Compute b = A*x = L*L^T*x
    // First compute L^T*x
    Vector<double> Ltx(3);
    auto Lt = transpose(L);
    for (size_t i = 0; i < 3; i++) {
        Ltx(i) = 0.0;
        for (size_t j = i; j < 3; j++) {
            Ltx(i) += Lt(i, j) * x_true(j);
        }
    }
    
    // Then compute L*(L^T*x)
    Vector<double> b(3);
    for (size_t i = 0; i < 3; i++) {
        b(i) = 0.0;
        for (size_t j = 0; j <= i; j++) {
            b(i) += L(i, j) * Ltx(j);
        }
    }
    
    // Now solve
    Vector<double> x_solved(3);
    cholesky_solve(L, x_solved, b);
    
    // Check solution
    for (size_t i = 0; i < 3; i++) {
        REQUIRE_THAT(x_solved(i), WithinAbs(x_true(i), 1e-8));
    }
}

#endif // defined(QUADPROGPP_MATRIX_BACKEND_BUILTIN)