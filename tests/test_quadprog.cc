#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <quadprog/quadprog.h>

using namespace quadprog;
using Catch::Matchers::WithinAbs;

TEST_CASE("Simple QP problem", "[quadprog]") {
    // min 0.5*x^T*G*x + g0^T*x
    // where x = [x1, x2]
    // G = [2, 0; 0, 2]
    // g0 = [-2, -5]
    // Solution should be [1, 2.5] (unconstrained)
    
    Matrix<double> G(2, 2);
    G(0, 0) = 2.0; G(0, 1) = 0.0;
    G(1, 0) = 0.0; G(1, 1) = 2.0;
    
    Vector<double> g0 = {-2.0, -5.0};
    
    Matrix<double> CE(2, 0);  // No equality constraints
    Vector<double> ce0(0);
    
    Matrix<double> CI(2, 0);  // No inequality constraints
    Vector<double> ci0(0);
    
    auto result = solve_quadprog(G, g0, CE, ce0, CI, ci0);

    REQUIRE(result.is_success());
    auto x = result.solution;
    REQUIRE(x.size() == 2u);
    
    REQUIRE_THAT(x(0), WithinAbs(1.0, 1e-6));
    REQUIRE_THAT(x(1), WithinAbs(2.5, 1e-6));
}

TEST_CASE("QP with equality constraint", "[quadprog]") {
    // min 0.5*x^T*G*x
    // subject to x1 + x2 = 1
    
    Matrix<double> G(2, 2);
    G(0, 0) = 2.0; G(0, 1) = 0.0;
    G(1, 0) = 0.0; G(1, 1) = 2.0;
    
    Vector<double> g0 = {0.0, 0.0};
    
    Matrix<double> CE(2, 1);
    CE(0, 0) = 1.0;
    CE(1, 0) = 1.0;
    
    Vector<double> ce0 = {-1.0};
    
    Matrix<double> CI(2, 0);
    Vector<double> ci0(0);
        
    auto result = solve_quadprog(G, g0, CE, ce0, CI, ci0);

    auto x = result.solution;
    REQUIRE(result.is_success());
    REQUIRE(x.size() == 2u);
    
    // Check constraint is satisfied
    REQUIRE_THAT(x(0) + x(1), WithinAbs(1.0, 1e-6));
}

TEST_CASE("Size mismatch throws exception", "[quadprog][error]") {
    Matrix<double> G(2, 2);
    Vector<double> g0 = {1.0, 2.0, 3.0};  // Wrong size!
    
    Matrix<double> CE(2, 0);
    Vector<double> ce0(0);
    Matrix<double> CI(2, 0);
    Vector<double> ci0(0);
    Vector<double> x(2);
    
    REQUIRE_THROWS_AS(
        solve_quadprog(G, g0, CE, ce0, CI, ci0),
        std::invalid_argument
    );
}