#pragma once

#include <quadprog/config.h>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <concepts>
#include <memory>

// Matrix backend includes
#ifdef QUADPROGPP_MATRIX_BACKEND_EIGEN
#include <Eigen/Dense>
#elif defined(QUADPROGPP_MATRIX_BACKEND_ARMADILLO)
#include <armadillo>
#endif

namespace quadprog
{
    // ============================================================================
    // Matrix/Vector type abstractions
    // ============================================================================

#ifdef QUADPROGPP_MATRIX_BACKEND_EIGEN
    using Matrix = Eigen::MatrixXd;
    using Vector = Eigen::VectorXd;

    inline Matrix make_matrix(size_t rows, size_t cols)
    {
        return Matrix::Zero(rows, cols);
    }

    inline Vector make_vector(size_t size)
    {
        return Vector::Zero(size);
    }

#elif defined(QUADPROGPP_MATRIX_BACKEND_ARMADILLO)
    using Matrix = arma::mat;
    using Vector = arma::vec;

    inline Matrix make_matrix(size_t rows, size_t cols)
    {
        return arma::mat(rows, cols, arma::fill::zeros);
    }

    inline Vector make_vector(size_t size)
    {
        return arma::vec(size, arma::fill::zeros);
    }

#else                          // BUILTIN
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

        size_t rows() const { return rows_; }
        size_t cols() const { return cols_; }
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

        inline constexpr size_t nrows() const { return rows_; }
        inline constexpr size_t ncols() const { return cols_; }

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
#ifndef NDEBUG
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
#ifndef NDEBUG
            if (i >= this->size())
                throw std::out_of_range("Vector index out of range");
#endif
            return std::vector<T>::operator[](i);
        }

        typename std::vector<T>::const_reference operator()(size_t i) const
        {
#ifndef NDEBUG
            if (i >= this->size())
                throw std::out_of_range("Vector index out of range");
#endif
            return std::vector<T>::operator[](i);
        }

        static VectorImpl<T> make_vector(size_t size)
        {
            return VectorImpl<T>(size, T(0));
        }

        friend double scalar_product(const VectorImpl<T> &x, const VectorImpl<T> &y);

        // Disable operator[] to avoid confusion
        T operator[](size_t i) = delete;
        const T &operator[](size_t i) const = delete;
    };

    template <typename T>
    double scalar_product(const VectorImpl<T> &x, const VectorImpl<T> &y)
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

#endif

} // namespace quadprog

// include implementation
#include <quadprog/array_impl.tpp>