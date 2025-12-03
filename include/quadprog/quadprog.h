#pragma once

#include <quadprog/config.h>
#include <vector>
#include <stdexcept>
#include <memory>
#include <cmath>
#include <concepts>
#include <limits>

#if defined(QUADPROGPP_MATRIX_BACKEND_BUILTIN)
#include <quadprog/array_impl.hh>
#elif defined(QUADPROGPP_MATRIX_BACKEND_EIGEN)
// Eigen backend
#include <Eigen/Dense>
template <typename T>
    using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template <typename T>
    using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;
#endif
namespace quadprog
{
    // ============================================================================
    // Solver Configuration
    // ============================================================================

    template <std::floating_point T>
    struct SolverOptions
    {
        T tolerance = std::numeric_limits<T>::epsilon(); // Numerical tolerance
        bool enable_logging = false;                     // Enable debug logging
        bool use_semidefinite_extension = false;         // Use Boland method for semi-definite
        bool graceful_exit = false;                      // Exit without throwing on failure

        SolverOptions() = default;
    };

    // ============================================================================
    // Solver Result
    // ============================================================================

    enum class SolverStatus
    {
        UNSOLVED,
        SUCCESS,
        INFEASIBLE,
        UNBOUNDED,
        MAX_ITERATIONS,
        NUMERICAL_ERROR
    };

    template <std::floating_point T>
    struct SolverResult
    {
        SolverStatus status = SolverStatus::UNSOLVED;
        Vector<T> solution;
        T objective_value = 0.0;
        size_t iterations = 0;
        std::vector<size_t> active_set;
        std::string message = "Unsolved";

        bool is_success() const { return status == SolverStatus::SUCCESS; }
    };

    // ============================================================================
    // Main Solver Interface
    // ============================================================================

    /**
     * @brief Solve a quadratic programming problem
     *
     * Solves: min 0.5 * x^T G x + g0^T x
     *         subject to: CE^T x = ce0
     *                    CI^T x >= ci0
     *
     * @param G Positive definite or semi-definite matrix (n x n)
     * @param g0 Linear term (n)
     * @param CE Equality constraint matrix (n x p)
     * @param ce0 Equality constraint vector (p)
     * @param CI Inequality constraint matrix (n x m)
     * @param ci0 Inequality constraint vector (m)
     * @param x Solution vector (n) - will be resized
     * @param options Solver options
     * @return SolverResult with status and diagnostics
     */
    template <std::floating_point T>
    SolverResult<T> solve_quadprog(
        const Matrix<T> &G,
        const Vector<T> &g0,
        const Matrix<T> &CE,
        const Vector<T> &ce0,
        const Matrix<T> &CI,
        const Vector<T> &ci0,
        const SolverOptions<T> &options = SolverOptions<T>());

    /**
     * @brief Legacy interface matching original QuadProg++
     *
     * @return Objective value at solution, or infinity if infeasible
     */
    double solve_quadprog_legacy(
        Matrix<double> &G,
        Vector<double> &g0,
        const Matrix<double> &CE,
        const Vector<double> &ce0,
        const Matrix<double> &CI,
        const Vector<double> &ci0,
        Vector<double> &x);

#ifdef QUADPROGPP_ENABLE_SEMIDEFINITE
    /**
     * @brief Solve semi-definite QP using Boland's extension
     *
     * Handles positive semi-definite G matrices by partitioning variables
     * into quadratic and linear parts.
     */
    template <std::floating_point T>
    SolverResult<T> solve_quadprog_semidefinite(
        const Matrix<T> &G,
        const Vector<T> &g0,
        const Matrix<T> &CE,
        const Vector<T> &ce0,
        const Matrix<T> &CI,
        const Vector<T> &ci0,
        const SolverOptions<T> &options = SolverOptions());
#endif

} // namespace quadprog

// Include implementation

#include <quadprog/quadprog.tpp>
