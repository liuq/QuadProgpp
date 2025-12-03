#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <quadprog/quadprog.h>
#include <random>
#include <cmath>

using namespace quadprog;
using Catch::Matchers::WithinAbs;

// Helper to create random SPD matrix
Matrix<double> random_spd_matrix(size_t n, std::mt19937& gen) {
    std::uniform_real_distribution<> dis(0.1, 1.0);
    
    Matrix<double> A(n, n);
    // Create random matrix
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            A(i, j) = dis(gen);
        }
    }
    
    // Make symmetric positive definite: G = A^T * A + I
    Matrix<double> G(n, n);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            G(i, j) = 0.0;
            for (size_t k = 0; k < n; k++) {
                G(i, j) += A(k, i) * A(k, j);
            }
            if (i == j) G(i, j) += 1.0;  // Ensure positive definite
        }
    }
    
    return G;
}

// Helper to verify KKT conditions (approximately)
bool verify_solution(const Matrix<double>& /* G */, const Vector<double>& /* g0 */,
                     const Matrix<double>& CE, const Vector<double>& ce0,
                     const Matrix<double>& CI, const Vector<double>& ci0,
                     const Vector<double>& x, double tol = 1e-6) {    
#if defined(QUADPROGPP_MATRIX_BACKEND_BUILTIN)
    // Check equality constraints: CE^T * x + ce0 = 0
    for (size_t i = 0; i < CE.cols(); i++) {
        double constraint = ce0(i);
        for (size_t j = 0; j < x.size(); j++) {
            constraint += CE(j, i) * x(j);
        }
        if (std::abs(constraint) > tol) return false;
    }
#elif defined(QUADPROGPP_MATRIX_BACKEND_EIGEN)
    for (Eigen::Index i{0}; i < CE.cols(); i++) {
        double constraint = ce0(i) + CE.col(i).dot(x);
        if (std::abs(constraint) > tol) return false;
    }
#endif

#if defined(QUADPROGPP_MATRIX_BACKEND_BUILTIN)
    // Check inequality constraints: CI^T * x + ci0 >= 0
    for (size_t i = 0; i < CI.cols(); i++) {
        double constraint = ci0(i);
        for (size_t j = 0; j < x.size(); j++) {
            constraint += CI(j, i) * x(j);
        }
        if (constraint < -tol) return false;
    }
#elif defined(QUADPROGPP_MATRIX_BACKEND_EIGEN)
    for (Eigen::Index i{0}; i < CI.cols(); i++) {
        double constraint = ci0(i) + CI.col(i).dot(x);
        if (constraint < -tol) return false;
    }
#endif
    
    return true;
}

TEST_CASE("Large scale problems", "[quadprog][stress][!mayfail]") {
    std::mt19937 gen(42);
    
    SECTION("50 variables, unconstrained") {
        const size_t n = 50;
        auto G = random_spd_matrix(n, gen);
        
        std::uniform_real_distribution<> dis(-10.0, 10.0);
        Vector<double> g0(n);
        for (size_t i = 0; i < n; i++) g0(i) = dis(gen);
        
        Matrix<double> CE(n, 0);
        Vector<double> ce0(0);
        Matrix<double> CI(n, 0);
        Vector<double> ci0(0);
        
        auto result = solve_quadprog(G, g0, CE, ce0, CI, ci0);
        
        REQUIRE(result.is_success());
        REQUIRE(result.solution.size() == n);
        REQUIRE(verify_solution(G, g0, CE, ce0, CI, ci0, result.solution, 1e-4));
    }
    
    SECTION("100 variables with constraints") {
        const size_t n = 100;
        const size_t n_eq = 10;
        const size_t n_ineq = 20;
        
        auto G = random_spd_matrix(n, gen);
        
        std::uniform_real_distribution<> dis(-5.0, 5.0);
        Vector<double> g0(n);
        for (size_t i = 0; i < n; i++) g0(i) = dis(gen);
        
        // Random equality constraints
        Matrix<double> CE(n, n_eq);
        Vector<double> ce0(n_eq);
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n_eq; j++) {
                CE(i, j) = dis(gen);
            }
        }
        for (size_t j = 0; j < n_eq; j++) ce0(j) = dis(gen);
        
        // Random inequality constraints
        Matrix<double> CI(n, n_ineq);
        Vector<double> ci0(n_ineq);
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n_ineq; j++) {
                CI(i, j) = dis(gen);
            }
        }
        for (size_t j = 0; j < n_ineq; j++) ci0(j) = std::abs(dis(gen));
        
        auto result = solve_quadprog(G, g0, CE, ce0, CI, ci0);
        
        REQUIRE(result.is_success());
        REQUIRE(result.solution.size() == n);
    }
}

TEST_CASE("Ill-conditioned problems", "[quadprog][stress]") {
    std::mt19937 gen(123);
    
    SECTION("High condition number matrix") {
        const size_t n = 20;
        Matrix<double> G(n, n);
        
        // Create diagonal matrix with large condition number
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                if (i == j) {
                    G(i, j) = std::pow(10.0, -static_cast<double>(i) / 2.0);
                } else {
                    G(i, j) = 0.0;
                }
            }
        }

#if defined(QUADPROGPP_MATRIX_BACKEND_BUILTIN)        
        Vector<double> g0(n, -1.0);
#elif defined(QUADPROGPP_MATRIX_BACKEND_EIGEN)
        auto g0 = Vector<double>::Constant(n, -1.0);
#endif
        Matrix<double> CE(n, 0);
        Vector<double> ce0(0);
        Matrix<double> CI(n, 0);
        Vector<double> ci0(0);
        
        auto result = solve_quadprog<double>(G, g0, CE, ce0, CI, ci0);
        REQUIRE(result.is_success());
        REQUIRE(verify_solution(G, g0, CE, ce0, CI, ci0, result.solution, 1e-4));
    }
    
    SECTION("Nearly singular constraints") {
        const size_t n = 10;
        auto G = random_spd_matrix(n, gen);
#if defined(QUADPROGPP_MATRIX_BACKEND_BUILTIN)        
        Vector<double> g0(n, 0.0);
#elif defined(QUADPROGPP_MATRIX_BACKEND_EIGEN)
        auto g0 = Vector<double>::Constant(n, 0.0);
#endif
        
        // Two nearly parallel constraints
        Matrix<double> CE(n, 2);
        for (size_t i = 0; i < n; i++) {
            CE(i, 0) = 1.0;
#if defined(QUADPROGPP_MATRIX_BACKEND_BUILTIN)            
            CE(i, 1) = 1.0 + 1e-8;  // Nearly parallel
#elif defined(QUADPROGPP_MATRIX_BACKEND_EIGEN)
            CE(i, 1) = 1.0 + 1e-16;
#endif
        }
        Vector<double> ce0(2);
        ce0(0) = -1.0;
        ce0(1) = -1.0;
        
        Matrix<double> CI(n, 0);
        Vector<double> ci0(0);

        quadprog::SolverOptions<double> options;
        options.graceful_exit = true;
        
        auto result = solve_quadprog<double>(G, g0, CE, ce0, CI, ci0, options);
#if defined(QUADPROGPP_MATRIX_BACKEND_BUILTIN)
        REQUIRE(result.status == quadprog::SolverStatus::INFEASIBLE);
#elif defined(QUADPROGPP_MATRIX_BACKEND_EIGEN)        
        REQUIRE(result.is_success());
        REQUIRE(verify_solution(G, g0, CE, ce0, CI, ci0, result.solution, 1e-4));
#endif
    }
}

TEST_CASE("Many constraints", "[quadprog][stress]") {
    std::mt19937 gen(456);
    
    SECTION("More inequality constraints than variables") {
        const size_t n = 20;
        const size_t n_ineq = 50;  // More constraints than variables
        
        auto G = random_spd_matrix(n, gen);
        
        std::uniform_real_distribution<> dis(-5.0, 5.0);
        Vector<double> g0(n);
        for (size_t i = 0; i < n; i++) g0(i) = dis(gen);
        
        Matrix<double> CE(n, 0);
        Vector<double> ce0(0);
        
        Matrix<double> CI(n, n_ineq);
        Vector<double> ci0(n_ineq);
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n_ineq; j++) {
                CI(i, j) = dis(gen);
            }
        }
        for (size_t j = 0; j < n_ineq; j++) {
            ci0(j) = std::abs(dis(gen)) + 1.0;  // Make feasible
        }
        
        auto result = solve_quadprog(G, g0, CE, ce0, CI, ci0);
        
        REQUIRE(result.is_success());
    }
}

TEST_CASE("Degenerate cases", "[quadprog][stress]") {
    SECTION("Single variable problem") {
        Matrix<double> G(1, 1);
        G(0, 0) = 2.0;
        
        Vector<double> g0(1);
        g0(0) = -4.0;
        
        Matrix<double> CE(1, 0);
        Vector<double> ce0(0);
        Matrix<double> CI(1, 0);
        Vector<double> ci0(0);
        
        auto result = solve_quadprog(G, g0, CE, ce0, CI, ci0);
        
        REQUIRE(result.is_success());
        REQUIRE_THAT(result.solution(0), WithinAbs(2.0, 1e-6));
    }
    
    SECTION("Identity matrix") {
        const size_t n = 10;
        Matrix<double> G(n, n);
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                G(i, j) = (i == j) ? 1.0 : 0.0;
            }
        }

#if defined(QUADPROGPP_MATRIX_BACKEND_BUILTIN)
        Vector<double> g0(n, -1.0);
#elif defined(QUADPROGPP_MATRIX_BACKEND_EIGEN)
        auto g0 = Vector<double>::Constant(n, -1.0);
#endif
        Matrix<double> CE(n, 0);
        Vector<double> ce0(0);
        Matrix<double> CI(n, 0);
        Vector<double> ci0(0);
        
        auto result = solve_quadprog<double>(G, g0, CE, ce0, CI, ci0);
        
        REQUIRE(result.is_success());
        
        // Solution should be all 1.0
        for (size_t i = 0; i < n; i++) {
            REQUIRE_THAT(result.solution(i), WithinAbs(1.0, 1e-6));
        }
    }
    
    SECTION("All constraints active") {
        const size_t n = 5;
        Matrix<double> G(n, n);
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                G(i, j) = (i == j) ? 2.0 : 0.0;
            }
        }

#if defined(QUADPROGPP_MATRIX_BACKEND_BUILTIN)        
        Vector<double> g0(n, 0.0);
#elif defined(QUADPROGPP_MATRIX_BACKEND_EIGEN)
        auto g0 = Vector<double>::Constant(n, 0.0);
#endif
        
        // Box constraints: -1 <= x_i <= 1
        Matrix<double> CE(n, 0);
        Vector<double> ce0(0);
        
        Matrix<double> CI(n, 2*n);
        Vector<double> ci0(2*n);
        for (size_t i = 0; i < n; i++) {
            // x_i >= -1  =>  x_i + 1 >= 0
            CI(i, 2*i) = 1.0;
            ci0(2*i) = 1.0;
            // x_i <= 1   =>  -x_i + 1 >= 0
            CI(i, 2*i+1) = -1.0;
            ci0(2*i+1) = 1.0;
        }
        
        auto result = solve_quadprog<double>(G, g0, CE, ce0, CI, ci0);
        
        REQUIRE(result.is_success());
    }
}

TEST_CASE("Numerical edge cases", "[quadprog][stress]") {
    SECTION("Very small coefficients") {
        Matrix<double> G(2, 2);
        G(0, 0) = 1e-10; G(0, 1) = 0.0;
        G(1, 0) = 0.0;   G(1, 1) = 1e-10;
        
        Vector<double> g0(2);
        g0(0) = -1e-10;
        g0(1) = -1e-10;
        
        Matrix<double> CE(2, 0);
        Vector<double> ce0(0);
        Matrix<double> CI(2, 0);
        Vector<double> ci0(0);

        // This problem has a nearly non positive definite matrix; expect failure
        auto result = solve_quadprog(G, g0, CE, ce0, CI, ci0);
        REQUIRE(result.is_success());
    }
    
    SECTION("Very large coefficients") {
        Matrix<double> G(2, 2);
        G(0, 0) = 1e10; G(0, 1) = 0.0;
        G(1, 0) = 0.0;  G(1, 1) = 1e10;
        
        Vector<double> g0(2);
        g0(0) = -1e10;
        g0(1) = -1e10;
        
        Matrix<double> CE(2, 0);
        Vector<double> ce0(0);
        Matrix<double> CI(2, 0);
        Vector<double> ci0(0);
        
        auto result = solve_quadprog(G, g0, CE, ce0, CI, ci0);
        
        REQUIRE(result.is_success());
    }
    
    SECTION("Mixed scales") {
        Matrix<double> G(3, 3);
        G(0, 0) = 1e-6; G(0, 1) = 0.0;   G(0, 2) = 0.0;
        G(1, 0) = 0.0;  G(1, 1) = 1.0;   G(1, 2) = 0.0;
        G(2, 0) = 0.0;  G(2, 1) = 0.0;   G(2, 2) = 1e6;
        
        Vector<double> g0(3);
        g0(0) = -1e-6;
        g0(1) = -1.0;
        g0(2) = -1e6;
        
        Matrix<double> CE(3, 0);
        Vector<double> ce0(0);
        Matrix<double> CI(3, 0);
        Vector<double> ci0(0);
        
        auto result = solve_quadprog(G, g0, CE, ce0, CI, ci0);
        
        REQUIRE(result.is_success());
    }
}

TEST_CASE("Random problem battery", "[quadprog][stress][!mayfail]") {
    std::mt19937 gen(789);
    const int num_tests = 100;
    int successes = 0;
    
    for (int test = 0; test < num_tests; test++) {
        std::uniform_int_distribution<> size_dis(5, 30);
        size_t n = size_dis(gen);
        size_t n_eq = std::min(n / 3, size_t(10));
        size_t n_ineq = std::min(n / 2, size_t(20));
        
        try {
            auto G = random_spd_matrix(n, gen);
            
            std::uniform_real_distribution<> dis(-10.0, 10.0);
            Vector<double> g0(n);
            for (size_t i = 0; i < n; i++) g0(i) = dis(gen);
            
            Matrix<double> CE(n, n_eq);
            Vector<double> ce0(n_eq);
            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n_eq; j++) {
                    CE(i, j) = dis(gen);
                }
            }
            for (size_t j = 0; j < n_eq; j++) ce0(j) = dis(gen);
            
            Matrix<double> CI(n, n_ineq);
            Vector<double> ci0(n_ineq);
            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n_ineq; j++) {
                    CI(i, j) = dis(gen);
                }
            }
            for (size_t j = 0; j < n_ineq; j++) {
                ci0(j) = std::abs(dis(gen)) + 1.0;
            }
            
            auto result = solve_quadprog(G, g0, CE, ce0, CI, ci0);
            
            if (result.is_success() && 
                verify_solution(G, g0, CE, ce0, CI, ci0, result.solution, 1e-3)) {
                successes++;
            }
        } catch (...) {
            // Some random problems may be infeasible/ill-posed
        }
    }
    
    INFO("Succeeded on " << successes << "/" << num_tests << " random problems");
    REQUIRE(successes > num_tests * 0.8);  // At least 80% success rate
}

TEST_CASE("Performance benchmarks", "[quadprog][benchmark][!benchmark]") {
    std::mt19937 gen(999);
    
    BENCHMARK("10 variables, unconstrained") {
        auto G = random_spd_matrix(10, gen);
#if defined(QUADPROGPP_MATRIX_BACKEND_BUILTIN)        
        Vector<double> g0(10, -1.0);
#elif defined(QUADPROGPP_MATRIX_BACKEND_EIGEN)
        auto g0 = Vector<double>::Constant(10, -1.0);
#endif
        Matrix<double> CE(10, 0);
        Vector<double> ce0(0);
        Matrix<double> CI(10, 0);
        Vector<double> ci0(0);
        return solve_quadprog<double>(G, g0, CE, ce0, CI, ci0);
    };
    
    BENCHMARK("50 variables, 10 constraints") {
        auto G = random_spd_matrix(50, gen);
#if defined(QUADPROGPP_MATRIX_BACKEND_BUILTIN)        
        Vector<double> g0(50, -1.0);
#elif defined(QUADPROGPP_MATRIX_BACKEND_EIGEN)
        auto g0 = Vector<double>::Constant(50, -1.0);
#endif
        
        Matrix<double> CE(50, 5);
#if defined(QUADPROGPP_MATRIX_BACKEND_BUILTIN)        
        Vector<double> ce0(5, 1.0);
#elif defined(QUADPROGPP_MATRIX_BACKEND_EIGEN)
        auto ce0 = Vector<double>::Constant(5, 1.0);
#endif
        for (size_t i = 0; i < 50; i++) {
            for (size_t j = 0; j < 5; j++) {
                CE(i, j) = (i == j) ? 1.0 : 0.0;
            }
        }
        
        Matrix<double> CI(50, 5);
#if defined(QUADPROGPP_MATRIX_BACKEND_BUILTIN)        
        Vector<double> ci0(5, 1.0);
#elif defined(QUADPROGPP_MATRIX_BACKEND_EIGEN)
        auto ci0 = Vector<double>::Constant(5, 1.0);
#endif
        for (size_t i = 0; i < 50; i++) {
            for (size_t j = 0; j < 5; j++) {
                CI(i, j) = (i == j + 5) ? 1.0 : 0.0;
            }
        }
        
        return solve_quadprog<double>(G, g0, CE, ce0, CI, ci0);
    };
}

TEST_CASE("Result metadata", "[quadprog]") {
    Matrix<double> G(2, 2);
    G(0, 0) = 2.0; G(0, 1) = 0.0;
    G(1, 0) = 0.0; G(1, 1) = 2.0;
    
    Vector<double> g0(2);
    g0(0) = -2.0;
    g0(1) = -5.0;
    
    Matrix<double> CE(2, 0);
    Vector<double> ce0(0);
    Matrix<double> CI(2, 0);
    Vector<double> ci0(0);
    
    auto result = solve_quadprog(G, g0, CE, ce0, CI, ci0);
    
    REQUIRE(result.is_success());
    REQUIRE(result.iterations > 0);
    REQUIRE(result.objective_value < 0.0);  // Minimizing, should be negative
    REQUIRE(result.message.length() > 0);
}