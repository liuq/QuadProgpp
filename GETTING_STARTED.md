# QuadProg++ Modernization - Quick Start Guide

## What Has Been Created

We've set up a modern C++17 framework for QuadProg++ with:

### ✅ Completed Structure

```
quadprogpp-modern/
├── CMakeLists.txt                 # Modern CMake with optional backends
├── README.md                      # Complete documentation
├── cmake/
│   └── QuadProgppConfig.cmake.in  # Package configuration
├── include/quadprog/
│   ├── config.h.in                # Build configuration
│   └── quadprog.h                 # Modern C++17 API
├── examples/
│   ├── CMakeLists.txt
│   └── simple_example.cc          # Example from Boland paper
└── docs/
    └── ROADMAP.md                 # Complete implementation plan
```

### Key Features

1. **Multi-Backend Support**
   - Built-in matrix implementation (no dependencies)
   - Optional Eigen3 support (FetchContent if not found)
   - Optional Armadillo support
   - Switchable via CMake options

2. **Modern C++20 API**
   ```cpp
   SolverResult solve_quadprog(
       const Matrix& G, const Vector& g0,
       const Matrix& CE, const Vector& ce0,
       const Matrix& CI, const Vector& ci0,
       const SolverOptions& options = {}
   );
   ```

3. **Backward Compatibility**
   - Legacy API still available
   - Can be used as drop-in replacement

4. **Boland Extension Ready**
   - Conditional compilation
   - Separate implementation file
   - Clear extension points

## Next Steps

### Step 1: Implement Built-in Matrix Class

Create `src/array_impl.cc`:

```cpp
// Implement basic linear algebra operations:
// - Matrix-vector multiplication
// - Matrix-matrix multiplication  
// - Cholesky decomposition with rank detection
// - Basic BLAS-like operations
```

Key function to implement:
```cpp
bool cholesky_decompose(const Matrix& A, Matrix& L, double tol) {
    // Detect rank by checking diagonal elements
    // If A(i,i) < tol, A is not full rank
    // This is where Boland extension kicks in
}
```

### Step 2: Port Goldfarb-Idnani Algorithm

Create `src/quadprog.cc`:

1. **Copy logic from original QuadProg++.cc**
2. **Modernize**:
   - Replace `Array<T>` with `std::vector<T>` or `Matrix`
   - Use `auto` for type deduction
   - Range-based loops
   - Const correctness

3. **Add diagnostics**:
   - Track iterations
   - Record active set
   - Return `SolverResult` instead of just `double`

### Step 3: Add Rank Detection

Modify Cholesky decomposition phase:

```cpp
// Around line 200-240 in original code
for (size_t j = 0; j < n; j++) {
    // Compute R[j][j]
    double sum = /* ... */;
    
    // CRITICAL: Rank detection
    if (sum < tolerance) {
        // Matrix is singular/semi-definite
        // Store j as index of "linear" variable
        linear_vars.push_back(j);
        R[j][j] = 0.0;  // Mark as zero
        continue;  // Skip this variable
    }
    
    R[j][j] = std::sqrt(sum);
    // ... rest of Cholesky
}
```

### Step 4: Implement Boland Extension

Create `src/quadprog_semidefinite.cc`:

Key functions from paper:

```cpp
// Section 3.1: Transform to special form
SDQPPSpecialForm transform_problem(
    const Matrix& G_original, 
    const Vector& g0,
    ...
) {
    // 1. Find rank of G via eigendecomposition
    // 2. Partition variables: x (quadratic), y (linear)
    // 3. Ensure D has full row rank
    // 4. Return transformed problem
}

// Section 3.3: Find initial active set
std::vector<size_t> find_initial_active_set(
    const Matrix& D, 
    const Vector& p,
    size_t l  // number of linear vars
) {
    // Solve: find A with |A| = l, D_A invertible, D_A^{-1}*p >= 0
    // This is Phase I of simplex method
    // Can use existing simplex implementations or write simple version
}

// Section 4.1: Look-ahead deactivation  
LookaheadStep compute_lookahead_direction(
    const ActiveSetData& current,
    size_t j  // constraint to activate
) {
    // When D_A would lose rank after deactivation:
    // 1. Compute (x̄, ȳ) = solution of SDQEP(A ∪ {j})
    //    using formulas from Proposition 1
    // 2. Step direction = (x̄ - x, ȳ - y) / t₁
    // 3. where t₁ from Lemma 3
}
```

### Step 5: Testing

Create `tests/test_goldfarb_idnani.cc`:

```cpp
TEST(GoldfarbIdnani, SimpleExample) {
    // Example from Boland paper Section 4.2.2
    Matrix G = make_matrix(2, 2);
    G(0,0) = 2.0; G(1,1) = 2.0;  // Diagonal, positive definite
    
    Vector g0 = {-2.0, -5.0};
    // ... constraints ...
    
    auto result = solve_quadprog(G, g0, CE, ce0, CI, ci0);
    
    EXPECT_TRUE(result.is_success());
    EXPECT_NEAR(x[0], 0.25, 1e-6);
    EXPECT_NEAR(x[1], 0.875, 1e-6);
    EXPECT_NEAR(result.objective_value, 0.9375, 1e-6);
}

TEST(Boland, SemidefiniteExample) {
    // Example with rank-deficient G
    Matrix G = make_matrix(2, 2);
    G(0,0) = 1.0; G(0,1) = 0.0;
    G(1,0) = 0.0; G(1,1) = 0.0;  // Rank 1
    
    // Should automatically detect and use Boland extension
    auto result = solve_quadprog(G, g0, CE, ce0, CI, ci0);
    
    EXPECT_TRUE(result.is_success());
    // Verify solution...
}
```

## Build and Test

```bash
# Configure
cmake -B build \
    -DCMAKE_BUILD_TYPE=Debug \
    -DQUADPROGPP_BUILD_TESTS=ON \
    -DQUADPROGPP_BUILD_EXAMPLES=ON

# Build
cmake --build build

# Run tests
cd build && ctest --output-on-failure

# Run example
./build/examples/simple_example
```

## Implementation Priority

### High Priority (Week 1-2)
1. ✅ Project structure (DONE)
2. **Built-in matrix implementation**
3. **Port Goldfarb-Idnani to modern C++**
4. **Basic tests**

### Medium Priority (Week 3-4)
5. **Rank detection in Cholesky**
6. **Variable partitioning for Boland**
7. **Initial active set finder**
8. **Look-ahead deactivation**

### Lower Priority (Week 5+)
9. Eigen/Armadillo backends
10. Advanced examples
11. Performance benchmarks
12. Full documentation

## Where to Get Help

### Key References

1. **Original Paper**: Goldfarb & Idnani (1983)
   - Section 2: Algorithm description
   - Implementation details

2. **Boland Paper** (attached PDF):
   - Section 3.1: Problem transformation
   - Section 3.2-3.4: Conditions and matrices
   - Section 4.1: Look-ahead method
   - Section 4.2.2: Example walkthrough

3. **Original Code**: QuadProgpp from GitHub
   - `QuadProg++.cc`: Lines 200-500 (core algorithm)
   - Cholesky factorization
   - Active set updates

### Critical Code Sections

From original QuadProg++.cc to understand:

- **Lines ~200-240**: Cholesky factorization (where to add rank detection)
- **Lines ~250-350**: Main solver loop
- **Lines ~400-450**: Constraint activation/deactivation
- **Lines ~500+**: Dual updates and step length calculation

## Questions?

Key decisions to make:

1. **Matrix implementation**: Start with built-in or integrate Eigen immediately?
   → Recommendation: **Built-in first** for simplicity

2. **Testing framework**: Google Test, Catch2, or custom?
   → Recommendation: **Google Test** (industry standard)

3. **Boland integration**: Separate library or always included?
   → Recommendation: **Optional compile flag** (current setup)

4. **Performance vs. Clarity**: Optimize now or later?
   → Recommendation: **Clarity first**, optimize after working

---

**Ready to start?** Begin with implementing the built-in matrix class in `src/array_impl.cc`!
