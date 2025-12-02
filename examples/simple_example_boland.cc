#include <quadprog/quadprog.h>
#include <iostream>
#include <iomanip>

/**
 * Simple example: minimize x1^2 + x2^2 - 2*x1 - 5*x2
 * subject to: x1 >= 0, x2 >= 0, x1 + x2 >= 1
 * 
 * This is Example 1 from Boland (1997) paper
 */

int main() {
    using namespace quadprog;
    
    std::cout << "QuadProg++ v" << QUADPROGPP_VERSION << std::endl;
    std::cout << "Simple QP Example" << std::endl;
    std::cout << std::string(50, '=') << std::endl;
    
    // Problem formulation
    std::cout << "\nProblem:" << std::endl;
    std::cout << "  min  x1^2 + x2^2 - 2*x1 - 5*x2" << std::endl;
    std::cout << "  s.t. x1 >= 0" << std::endl;
    std::cout << "       x2 >= 0" << std::endl;
    std::cout << "       x1 + x2 >= 1" << std::endl;
    
    // G matrix (Hessian): [[2, 0], [0, 2]]
    Matrix<double> G(2, 2);
    G(0, 0) = 2.0; G(0, 1) = 0.0;
    G(1, 0) = 0.0; G(1, 1) = 2.0;
    
    // g0 vector (linear term): [-2, -5]
    Vector<double> g0(2);
    g0(0) = -2.0;
    g0(1) = -5.0;
    
    // No equality constraints
    Matrix<double> CE(2, 0);
    Vector<double> ce0;
    
    // Inequality constraints: CI^T x >= ci0
    // x1 >= 0, x2 >= 0, x1 + x2 >= 1
    Matrix<double> CI(2, 3);
    CI(0, 0) = 1.0; CI(1, 0) = 0.0;  // x1 >= 0
    CI(0, 1) = 0.0; CI(1, 1) = 1.0;  // x2 >= 0
    CI(0, 2) = 1.0; CI(1, 2) = 1.0;  // x1 + x2 >= 1
    
    Vector<double> ci0(3);
    ci0(0) = 0.0;
    ci0(1) = 0.0;
    ci0[2] = 1.0;
    
    // Solution vector
    Vector<double> x;
    
    // Solve
    std::cout << "\nSolving..." << std::endl;
    auto result = solve_quadprog(G, g0, CE, ce0, CI, ci0);
    
    // Print results
    std::cout << std::string(50, '=') << std::endl;
    std::cout << "Results:" << std::endl;
    std::cout << std::string(50, '=') << std::endl;
    
    if (result.is_success()) {
        std::cout << "Status: SUCCESS" << std::endl;
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Solution: x = [" << x(0) << ", " << x(1) << "]" << std::endl;
        std::cout << "Objective value: " << result.objective_value << std::endl;
        std::cout << "Iterations: " << result.iterations << std::endl;
        
        std::cout << "\nActive constraints: ";
        if (result.active_set.empty()) {
            std::cout << "none";
        } else {
            for (size_t i : result.active_set) {
                std::cout << i << " ";
            }
        }
        std::cout << std::endl;
        
        // Verify solution
        double obj = 0.5 * (x(0)*x(0)*G(0,0) + 2*x(0)*x(1)*G(0,1) + x(1)*x(1)*G(1,1))
                   + g0(0)*x(0) + g0(1)*x(1);
        std::cout << "\nVerification:" << std::endl;
        std::cout << "  Computed objective: " << obj << std::endl;
        std::cout << "  Constraint x1 >= 0: " << (x(0) >= -1e-6 ? "✓" : "✗") << std::endl;
        std::cout << "  Constraint x2 >= 0: " << (x(1) >= -1e-6 ? "✓" : "✗") << std::endl;
        std::cout << "  Constraint x1+x2 >= 1: " << ((x(0)+x(1)) >= 1.0-1e-6 ? "✓" : "✗") << std::endl;
        
        // Expected solution for this problem: x = [0.25, 0.875], obj = 0.9375
        std::cout << "\nExpected solution: x = [0.25, 0.875], obj = 0.9375" << std::endl;
        
    } else {
        std::cout << "Status: " << static_cast<int>(result.status) << std::endl;
        std::cout << "Message: " << result.message << std::endl;
    }
    
    return result.is_success() ? 0 : 1;
}
