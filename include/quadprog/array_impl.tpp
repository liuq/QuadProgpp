#include <quadprog/quadprog.h>
#include <quadprog/array_impl.hh>
#include <quadprog/internal/logging.h>

// This file contains the implementation of built-in matrix operations
// when not using Eigen or Armadillo

namespace quadprog
{

    // Built-in matrix operations are implemented inline in the header
    // for now. Future implementations can be added here.
    template <std::floating_point T>
    LowerTriangularMatrix<T> cholesky_decompose(const Matrix<T> &G, T tolerance)
    {
        const size_t n = G.rows();
        LowerTriangularMatrix<T> L(n, 0.0);

        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j <= i; ++j)
            {
                T sum = G(i, j);

                for (size_t k = 0; k < j; ++k)
                {
                    sum -= L(i, k) * L(j, k);
                }

                if (i == j)
                {
                    if (sum <= tolerance)
                    {
                        QUADPROG_TRACE_MATRIX("G", G);
                        QUADPROG_TRACE("Matrix is not positive definite, sum at index {} is {}", i, sum);
                        throw std::invalid_argument("Matrix G is not positive definite");
                    }
                    L(i, j) = std::sqrt(sum);
                }
                else
                {
                    L(i, j) = sum / L(j, j);
                }
            }
        }

        return L;
    }

    // Forward elimination: solves Ly = b where L is lower triangular
    template <std::floating_point T>
    void forward_elimination(const LowerTriangularMatrix<T> &L, Vector<T> &y, const Vector<T> &b)
    {
        size_t n = L.nrows();
        y(0) = b(0) / L(0, 0);
        for (size_t i = 1; i < n; i++)
        {
            y(i) = b(i);
            for (size_t j = 0; j < i; j++)
                y(i) -= L(i, j) * y(j);
            y(i) /= L(i, i);
        }
    }

    // Backward substitution: solves Ux = b where U is upper triangular
    template <std::floating_point T>
    void backward_substitution(const UpperTriangularMatrix<T> &U, Vector<T> &x, const Vector<T> &b)
    {
        size_t n = U.nrows();
        x(n - 1) = b(n - 1) / U(n - 1, n - 1);
        for (int i = n - 2; i >= 0; i--) // Note: going backwards
        {
            x(i) = b(i);
            for (size_t j = i + 1; j < n; j++) // Note: j > i (upper part)
                x(i) -= U(i, j) * x(j);
            x(i) /= U(i, i);
        }
    }

    // Solve L * L^T * x = b
    template <std::floating_point T>
    void cholesky_solve(const LowerTriangularMatrix<T> &L, Vector<T> &x, const Vector<T> &b)
    {
        size_t n = L.nrows(); // Use size_t instead of int for consistency
#ifndef QUADPROGPP_NO_BOUNDS_CHECK
        if (b.size() != n || x.size() != n)
        {
            throw std::invalid_argument("Dimension mismatch in cholesky_solve");
        }
#endif
        Vector<T> y(n);

        /* Solve L * y = b */
        forward_elimination(L, y, b);

        QUADPROG_TRACE_VECTOR("y after forward elimination", y);

        /* Solve L^T * x = y */
        backward_substitution(transpose(L), x, y);
        QUADPROG_TRACE_VECTOR("x after backward substitution", x);
    }

} // namespace quadprog
