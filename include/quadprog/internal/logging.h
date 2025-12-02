#pragma once

#ifdef QUADPROGPP_ENABLE_LOGGING
#include <spdlog/logger.h>
#include <spdlog/sinks/stdout_color_sinks.h>

namespace quadprog
{
    namespace logging
    {
        inline spdlog::logger *get_logger()
        {
            static spdlog::logger *logger = []()
            {
                auto l = spdlog::stdout_color_mt("quadprog").get();
                l->set_level(spdlog::level::trace);
                return l;
            }();
            return logger;
        }

#define g_logger quadprog::logging::get_logger()
    }
}
#include <sstream>
#include <iomanip>

// Format a vector for logging
// Format a vector for logging with optional size limit
template <typename Vector>
std::string format_vector(const Vector &v, int precision = 4, size_t max_elements = 0)
{
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision);

    // If max_elements is 0, show all
    size_t elements_to_show = (max_elements == 0) ? v.size() : std::min(v.size(), max_elements);

    oss << "[";
    if (elements_to_show < v.size())
    {
        oss << "size=" << v.size() << ": ";
    }

    for (size_t i = 0; i < elements_to_show; ++i)
    {
        if (i > 0)
            oss << ", ";
        oss << v(i);
    }

    if (elements_to_show < v.size())
    {
        oss << ", ...";
    }
    oss << "]";

    return oss.str();
}

// Format a matrix for logging with optional size limits
template <typename Matrix>
std::string format_matrix(const Matrix &m, int precision = 4,
                          size_t max_rows = 0, size_t max_cols = 0)
{
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision);

    // If max_rows/max_cols is 0, show all
    size_t rows_to_show = (max_rows == 0) ? m.nrows() : std::min(m.nrows(), max_rows);
    size_t cols_to_show = (max_cols == 0) ? m.ncols() : std::min(m.ncols(), max_cols);

    oss << "\n  [" << m.nrows() << "x" << m.ncols() << " matrix";
    if (rows_to_show < m.nrows() || cols_to_show < m.ncols())
    {
        oss << ", showing " << rows_to_show << "x" << cols_to_show;
    }
    oss << "]\n";
    
    for (size_t i = 0; i < rows_to_show; ++i)
    {
        size_t start_col, end_col;
        if (m.is_upper_triangular()) {
            start_col = i;
            end_col = std::min(cols_to_show, m.ncols());
        }
        else if (m.is_lower_triangular()) {
            start_col = 0;
            end_col = std::min(cols_to_show, i + 1);
        }
        else
        {
            start_col = 0;
            end_col = std::min(cols_to_show, m.ncols());
        } 
        oss << "  [";
        for (size_t j = start_col; j < end_col; ++j)
        {
            if (j > 0)
                oss << ", ";
            oss << std::setw(precision + 4) << m(i, j);
        }
        if (cols_to_show < m.ncols())
        {
            oss << ", ...";
        }
        oss << "]\n";
    }
    if (rows_to_show < m.nrows())
    {
        oss << "  ...\n";
    }

    return oss.str();
}

#define QUADPROG_TRACE(...) \
    if (quadprog::logging::get_logger()->should_log(spdlog::level::trace)) \
        quadprog::logging::get_logger()->trace(__VA_ARGS__)

#define QUADPROG_TRACE_VECTOR(name, vec, ...) \
    if (quadprog::logging::get_logger()->should_log(spdlog::level::trace)) \
        quadprog::logging::get_logger()->trace("{}: {}", name, format_vector(vec __VA_OPT__(,) __VA_ARGS__))

#define QUADPROG_TRACE_MATRIX(name, mat, ...) \
    if (quadprog::logging::get_logger()->should_log(spdlog::level::trace)) \
        quadprog::logging::get_logger()->trace("{}: {}", name, format_matrix(mat __VA_OPT__(,) __VA_ARGS__))
#else
#define QUADPROG_TRACE(...) ((void)0)
#define QUADPROG_TRACE_VECTOR(name, vec, ...) ((void)0)
#define QUADPROG_TRACE_MATRIX(name, mat, ...) ((void)0)
#endif