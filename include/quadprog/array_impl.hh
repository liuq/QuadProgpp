#pragma once

#include <quadprog/config.h>
#include <quadprog/internal/logging.h>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <concepts>
#include <memory>

// This file contains the implementation of built-in matrix operations
// when not using Eigen or Armadillo

namespace quadprog
{
    // ============================================================================
    // Matrix/Vector type builtin implementations 
    // ============================================================================
    template <typename T>
    class MatrixImpl
    {
    protected:
        std::vector<T> data_;
        size_t rows_;
        size_t cols_;

    public:
        MatrixImpl() : rows_(0), cols_(0) {}

        MatrixImpl(size_t r, size_t c) : MatrixImpl(r, c, T(0)) {}

        MatrixImpl(size_t r, size_t c, const T &val)
            : data_(r * c, val), rows_(r), cols_(c) {}

        size_t size() const { return data_.size(); }

        void resize(size_t r, size_t c)
        {
            rows_ = r;
            cols_ = c;
            data_.resize(r * c);
            std::fill(data_.begin(), data_.end(), T(0));
        }

        T &operator()(size_t i, size_t j)
        {
#ifdef QUADPROGPP_BOUNDS_CHECK // In release mode, skip the check for performance
            if (i >= rows_ || j >= cols_)
                throw std::out_of_range("Matrix index out of range");
#endif
            return data_[i * cols_ + j];
        }

        const T &operator()(size_t i, size_t j) const
        {
#ifdef QUADPROGPP_BOUNDS_CHECK // In release mode, skip the check for performance
            if (i >= rows_ || j >= cols_)
                throw std::out_of_range("Matrix index out of range");
#endif
            return data_[i * cols_ + j];
        }

        T *data() { return data_.data(); }
        const T *data() const { return data_.data(); }

        // Row access
        std::vector<T> row(size_t i) const
        {
#ifdef QUADPROGPP_BOUNDS_CHECK // In release mode, skip the check for performance
            if (i >= rows_)
                throw std::out_of_range("Matrix row index out of range");
#endif
            std::vector<T> result(cols_);
            for (size_t j = 0; j < cols_; ++j)
            {
                result[j] = (*this)(i, j);
            }
            return result;
        }

        // Column access
        std::vector<T> col(size_t j) const
        {
#ifdef QUADPROGPP_BOUNDS_CHECK // In release mode, skip the check for performance
            if (j >= cols_)
                throw std::out_of_range("Matrix column index out of range");
#endif
            std::vector<T> result(rows_);
            for (size_t i = 0; i < rows_; ++i)
            {
                result[i] = (*this)(i, j);
            }
            return result;
        }

        static MatrixImpl<T> make_matrix(size_t rows, size_t cols)
        {
            return MatrixImpl<T>(rows, cols, T(0));
        }

        inline constexpr size_t rows() const { return rows_; }
        inline constexpr size_t cols() const { return cols_; }

        constexpr bool is_triangular() const { return false; }
        constexpr bool is_lower_triangular() const { return false; }
        constexpr bool is_upper_triangular() const { return false; }
    };

    // Forward declarations
    template <typename T>
    class LowerTriangularMatrix;

    template <typename T>
    class UpperTriangularMatrix;

    // Forward declare transpose functions
    template <typename T>
    UpperTriangularMatrix<T> transpose(const LowerTriangularMatrix<T> &L);

    template <typename T>
    LowerTriangularMatrix<T> transpose(const UpperTriangularMatrix<T> &U);

    template <typename T>
    class LowerTriangularMatrix
    {
    protected:
    protected:
        std::shared_ptr<std::vector<T>> data_;
        size_t n_;

        // Private constructor - only accessible via transpose()
        LowerTriangularMatrix(std::shared_ptr<std::vector<T>> data, size_t n)
            : data_(data), n_(n) {}

        friend UpperTriangularMatrix<T> transpose<>(const LowerTriangularMatrix<T> &L);
        friend LowerTriangularMatrix<T> transpose<>(const UpperTriangularMatrix<T> &U);

    public:
        LowerTriangularMatrix(size_t n)
            : data_(std::make_shared<std::vector<T>>(n * (n + 1) / 2, T(0))), n_(n) {}

        LowerTriangularMatrix(size_t n, const T &val)
            : data_(std::make_shared<std::vector<T>>(n * (n + 1) / 2, val)), n_(n) {}

        size_t nrows() const { return n_; }
        size_t ncols() const { return n_; }
        size_t size() const { return n_; }

        T &operator()(size_t i, size_t j)
        {
#ifdef QUADPROGPP_BOUNDS_CHECK // In release mode, skip the check for performance
            if (i < j)
                throw std::out_of_range("Accessing upper part of LowerTriangularMatrix");
#endif
            return (*data_)[i * (i + 1) / 2 + j];
        }

        const T &operator()(size_t i, size_t j) const
        {
#ifdef QUADPROGPP_BOUNDS_CHECK // In release mode, skip the check for performance
            if (i < j)
                throw std::out_of_range("Accessing upper part of LowerTriangularMatrix");
#endif
            return (*data_)[i * (i + 1) / 2 + j];
        }

        constexpr bool is_triangular() const { return true; }
        constexpr bool is_lower_triangular() const { return true; }
        constexpr bool is_upper_triangular() const { return false; }

    protected:
        friend class UpperTriangularMatrix<T>;

        // Access to shared data for creating views
        std::shared_ptr<std::vector<T>> data_ptr() const { return data_; }
    };

    template <typename T>
    class UpperTriangularMatrix
    {
    protected:
    protected:
        std::shared_ptr<std::vector<T>> data_; // Shared with LowerTriangularMatrix
        size_t n_;

        // Private constructor - only accessible via transpose()
        UpperTriangularMatrix(std::shared_ptr<std::vector<T>> data, size_t n)
            : data_(data), n_(n) {}

        friend UpperTriangularMatrix<T> transpose<>(const LowerTriangularMatrix<T> &L);
        friend LowerTriangularMatrix<T> transpose<>(const UpperTriangularMatrix<T> &U);

    public:
        size_t nrows() const { return n_; }
        size_t ncols() const { return n_; }
        size_t size() const { return n_; }

        UpperTriangularMatrix(size_t n)
            : data_(std::make_shared<std::vector<T>>(n * (n + 1) / 2, T(0))), n_(n) {}

        UpperTriangularMatrix(size_t n, const T &val)
            : data_(std::make_shared<std::vector<T>>(n * (n + 1) / 2, val)), n_(n) {}

        T &operator()(size_t i, size_t j)
        {
#ifdef QUADPROGPP_BOUNDS_CHECK // In release mode, skip the check for performance
            if (i > j)
                throw std::out_of_range("Accessing lower part of UpperTriangularMatrix");
#endif
            // U(i,j) = L(j,i) - transpose the indices and use lower triangular indexing
            return (*data_)[j * (j + 1) / 2 + i];
        }

        const T &operator()(size_t i, size_t j) const
        {
#ifndef QUADPROGPP_BOUNDS_CHECK // In release mode, skip the check for performance
            if (i > j)
                throw std::out_of_range("Accessing lower part of UpperTriangularMatrix");
#endif
            return (*data_)[j * (j + 1) / 2 + i];
        }

        constexpr bool is_triangular() const { return true; }
        constexpr bool is_lower_triangular() const { return false; }
        constexpr bool is_upper_triangular() const { return true; }

    protected:
        friend class LowerTriangularMatrix<T>;

        // Access to shared data for creating views
        std::shared_ptr<std::vector<T>> data_ptr() const { return data_; }
    };

    // Transpose operations - zero-copy views
    template <typename T>
    UpperTriangularMatrix<T> transpose(const LowerTriangularMatrix<T> &L)
    {
        return UpperTriangularMatrix<T>(L.data_ptr(), L.size());
    }

    template <typename T>
    LowerTriangularMatrix<T> transpose(const UpperTriangularMatrix<T> &U)
    {
        return LowerTriangularMatrix<T>(U.data_ptr(), U.size());
    }

    template <typename T>
    class VectorImpl : public std::vector<T>
    {
    public:
        using std::vector<T>::vector; // Inherit constructors

        VectorImpl(size_t size) : std::vector<T>(size, T(0)) {}

        void resize(size_t r, size_t c)
        {
            std::vector<T>::resize(r * c, T(0));
        }

        // Use the proper reference type from vector
        typename std::vector<T>::reference operator()(size_t i)
        {
#ifndef QUADPROGPP_BOUNDS_CHECK 
            if (i >= this->size())
                throw std::out_of_range("Vector index out of range");
#endif
            return std::vector<T>::operator[](i);
        }

        typename std::vector<T>::const_reference operator()(size_t i) const
        {
#ifndef QUADPROGPP_BOUNDS_CHECK
            if (i >= this->size())
                throw std::out_of_range("Vector index out of range");
#endif
            return std::vector<T>::operator[](i);
        }

        static VectorImpl<T> make_vector(size_t size)
        {
            return VectorImpl<T>(size, T(0));
        }

        friend T scalar_product(const VectorImpl<T> &x, const VectorImpl<T> &y);

        // Disable operator[] to avoid confusion
        T operator[](size_t i) = delete;
        const T &operator[](size_t i) const = delete;
    };

    template <typename T>
    T scalar_product(const VectorImpl<T> &x, const VectorImpl<T> &y)
    {
        if (x.size() != y.size())
        {
            throw std::invalid_argument("Vector size mismatch in scalar_product");
        }
        T result = 0.0;
        for (size_t i = 0; i < x.size(); ++i)
        {
            result += x(i) * y(i);
        }
        return result;
    }

    template <typename T>
    using Matrix = MatrixImpl<T>;

    template <typename T>
    using Vector = VectorImpl<T>;

    // ============================================================================
    // Utility Functions
    // ============================================================================

    /**
     * @brief Check if matrix is positive definite
     */
    // bool is_positive_definite(const Matrix<double> &G, double tolerance = DEFAULT_TOLERANCE);

    // /**
    //  * @brief Compute matrix rank
    //  */
    // size_t matrix_rank(const Matrix<double> &M, double tolerance = DEFAULT_TOLERANCE);

    /**
     * @brief Cholesky decomposition: G = L * L^T
     * @return L, lower triangular matrix
     * @throw std::invalid_argument if G is not positive definite
     */
    template <std::floating_point T>
    LowerTriangularMatrix<T> cholesky_decompose(const Matrix<T> &G, T tolerance);

    /**
     * @brief Solve L * y = b with forward elimination, where L is lower triangular
     */
    template <std::floating_point T>
    void forward_elimination(const LowerTriangularMatrix<T> &L, Vector<T> &y, const Vector<T> &b);

    /**
     * @brief Solve U * y = b with backward substitution, where U is upper triangular
     */
    template <std::floating_point T>
    void backward_substitution(const UpperTriangularMatrix<T> &U, Vector<T> &y, const Vector<T> &b);

    /**
     * @brief Solve L * L^T * x = b using the Cholesky factor L
     */
    template <std::floating_point T>
    void cholesky_solve(const LowerTriangularMatrix<T> &L, Vector<T> &x, const Vector<T> &b);

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



