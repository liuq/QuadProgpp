#include <cmath>
#include <limits>
#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <vector>
#include <quadprog/array_impl.hh>
#include <quadprog/quadprog.h>
#include <quadprog/internal/logging.h>

namespace quadprog
{

  // ============================================================================
  // Utility Functions
  // ============================================================================

  // Givens rotation and related helper functions
  void compute_d(Vector<double> &d, const Matrix<double> &J, const Vector<double> &np);
  void update_z(Vector<double> &z, const Matrix<double> &J, const Vector<double> &d, size_t iq);
  void update_r(const Matrix<double> &R, Vector<double> &r, const Vector<double> &d, size_t iq);
  bool add_constraint(Matrix<double> &R, Matrix<double> &J, Vector<double> &d, size_t &iq, double &rnorm);
  // TODO: check p and l types
  void delete_constraint(Matrix<double> &R, Matrix<double> &J, Vector<int> &A, Vector<double> &u, size_t n, size_t p, size_t &iq, int l);
  double distance(double a, double b);

  // ============================================================================
  // Main Solver (Stub Implementation)
  // ============================================================================

  SolverResult solve_quadprog(
      const Matrix<double> &G,
      const Vector<double> &g0,
      const Matrix<double> &CE,
      const Vector<double> &ce0,
      const Matrix<double> &CI,
      const Vector<double> &ci0,
      const SolverOptions &options)
  {
    SolverResult result;

    // Validate inputs
    const size_t n = G.rows();
    const size_t p = CE.cols();
    const size_t m = CI.cols();

    if (G.cols() != n)
    {
      if (!options.graceful_exit)
        throw std::invalid_argument("G matrix must be square, got non-square matrix " + std::to_string(G.rows()) + "x" + std::to_string(G.cols()));

      result.status = SolverStatus::NUMERICAL_ERROR;
      result.message = "G matrix must be square";
      return result;
    }

    if (g0.size() != n)
    {
      if (!options.graceful_exit)
        throw std::invalid_argument("g0 vector size mismatch with respect to G, expected " + std::to_string(n) + ", got " + std::to_string(g0.size()));
      result.status = SolverStatus::NUMERICAL_ERROR;
      result.message = "g0 vector size mismatch";
      return result;
    }

    if (CE.nrows() != n)
    {
      if (!options.graceful_exit)
        throw std::invalid_argument("The matrix CE is incompatible (incorrect number of rows " + std::to_string(CE.nrows()) + " , expecting " + std::to_string(n) + ")");
      result.status = SolverStatus::NUMERICAL_ERROR;
      result.message = "CE matrix row size mismatch";
      return result;
    }

    if (ce0.size() != p)
    {
      if (!options.graceful_exit)
        throw std::invalid_argument("The vector ce0 is incompatible (incorrect dimension " + std::to_string(ce0.size()) + ", expecting " + std::to_string(p) + ")");
      result.status = SolverStatus::NUMERICAL_ERROR;
      result.message = "ce0 vector size mismatch";
      return result;
    }

    if (CI.nrows() != n)
    {
      if (!options.graceful_exit)
        throw std::invalid_argument("The matrix CI is incompatible (incorrect number of rows " + std::to_string(CI.nrows()) + " , expecting " + std::to_string(n) + ")");
      result.status = SolverStatus::NUMERICAL_ERROR;
      result.message = "CI matrix row size mismatch";
      return result;
    }
    if (ci0.size() != m)
    {
      if (!options.graceful_exit)
        throw std::invalid_argument("The vector ci0 is incompatible (incorrect dimension " + std::to_string(ci0.size()) + ", expecting " + std::to_string(m) + ")");
      result.status = SolverStatus::NUMERICAL_ERROR;
      result.message = "ci0 vector size mismatch";
      return result;
    }

    // Check if G is positive definite
    // bool is_pd = is_positive_definite(G, options.tolerance);

    // if (!is_pd && !options.use_semidefinite_extension)
    // {
    //   if (!options.graceful_exit)
    //     throw std::invalid_argument("G is not positive definite. Enable semi-definite extension for semi-definite problems.");
    //   result.status = SolverStatus::NUMERICAL_ERROR;
    //   result.message = "G is not positive definite. Enable semi-definite extension for semi-definite problems.";
    //   return result;
    // }

    size_t i, j, k; /* indices */    
    int l; 
    int ip;            // this is the index of the constraint to be added to the active set
    Matrix<double> R(n, n), J(n, n);
    Vector<double> s(m + p), z(n), r(m + p), d(n), np(n), u(m + p), x(n), x_old(n), u_old(m + p);
    double f_value, psi, c1, c2, sum, ss, R_norm;
    constexpr double inf = std::numeric_limits<double>::infinity();
    double t, t1, t2; /* t is the step lenght, which is the minimum of the partial step length t1 and the full step length t2 */
    Vector<int> A(m + p, 0), A_old(m + p, 0), iai(m + p, 0);
    size_t iq, iter = 0;
    Vector<bool> iaexcl(m + p, false);

    QUADPROG_TRACE("Starting solve_quadprog");
    QUADPROG_TRACE_MATRIX("G", G);
    QUADPROG_TRACE_VECTOR("g0", g0);
    QUADPROG_TRACE_MATRIX("CE", CE);
    QUADPROG_TRACE_VECTOR("ce0", ce0);
    QUADPROG_TRACE_MATRIX("CI", CI);
    QUADPROG_TRACE_VECTOR("ci0", ci0);

    /* Preprocessing phase */

    /* compute the trace of the original matrix G */
    c1 = 0.0;
    for (i = 0; i < n; i++)
      c1 += G(i, i);
    /* decompose the matrix G in the form L^T L */
    LowerTriangularMatrix<double> G_ = cholesky_decompose(G, options.tolerance);
    QUADPROG_TRACE_MATRIX("L from Cholesky", G_);
    R_norm = 1.0; /* this variable will hold the norm of the matrix R */
    /* compute the inverse of the factorized matrix G^-1, this is the initial value for H */
    c2 = 0.0;
    for (i = 0; i < n; i++)
    {
      // Set d to the i-th unit vector
      d(i) = 1.0;
      forward_elimination(G_, z, d);
      for (j = 0; j < n; j++)
        J(i, j) = z(j);
      c2 += z(i);
      // Reset d
      d(i) = 0.0;
    }
    QUADPROG_TRACE_MATRIX("G inverse", J);
    /* c1 * c2 is an estimate for cond(G) */

    /*
     * Find the unconstrained minimizer of the quadratic form 1/2 * x G x + g0 x
     * this is a feasible point in the dual space
     * x = G^-1 * g0
     */
    cholesky_solve(G_, x, g0);
    // Negate x
    for (i = 0; i < n; i++)
      x(i) = -x(i);
    /* and compute the current solution value */
    f_value = 0.5 * scalar_product<double>(g0, x);
    QUADPROG_TRACE("Unconstrained solution: {}", f_value);
    QUADPROG_TRACE_VECTOR("x", x);

    /* Add equality constraints to the working set A */
    iq = 0;
    for (i = 0; i < p; i++)
    {
      for (j = 0; j < n; j++)
        np(j) = CE(j, i);
      compute_d(d, J, np);
      update_z(z, J, d, iq);
      update_r(R, r, d, iq);

      QUADPROG_TRACE_MATRIX("R before adding constraint", R, n, iq);
      QUADPROG_TRACE_VECTOR("d", d);
      QUADPROG_TRACE_VECTOR("z", z);
      QUADPROG_TRACE_VECTOR("r", r);

      /* compute full step length t2: i.e., the minimum step in primal space s.t. the contraint becomes feasible */
      t2 = 0.0;
      if (std::fabs(scalar_product<double>(z, z)) > DEFAULT_TOLERANCE) // i.e. z != 0
        t2 = (-scalar_product<double>(np, x) - ce0(i)) / scalar_product<double>(z, np);

      /* set x = x + t2 * z */
      for (k = 0; k < n; k++)
        x(k) += t2 * z(k);

      /* set u = u+ */
      u(iq) = t2;
      for (k = 0; k < iq; k++)
        u(k) -= t2 * r(k);

      /* compute the new solution value */
      f_value += 0.5 * (t2 * t2) * scalar_product<double>(z, np);
      A(iq) = -i - 1;

      if (!add_constraint(R, J, d, iq, R_norm))
      {
        // Equality constraints are linearly dependent
        if (!options.graceful_exit)
        {
          throw std::runtime_error("Equality constraints are linearly dependent");
        }
        result.status = SolverStatus::INFEASIBLE;
        result.message = "Equality constraints are linearly dependent";
        return result;
      }
    }

    /* set iai = K \ A */
    for (i = 0; i < m; i++)
      iai(i) = i;

  l1:
    iter++;
    QUADPROG_TRACE("Iteration {} from l1", iter);
    QUADPROG_TRACE_VECTOR("x", x);
    /* step 1: choose a violated constraint */
    for (i = p; i < iq; i++)
    {
      ip = A(i);
      iai(ip) = -1;
    }

    /* compute s[x] = ci^T * x + ci0 for all elements of K \ A */
    ss = 0.0;
    psi = 0.0; /* the sum of all infeasibilities */
    ip = 0;    /* ip will be the index of the chosen violated constraint */
    for (i = 0; i < m; i++)
    {
      iaexcl(i) = true;
      sum = 0.0;
      for (j = 0; j < n; j++)
        sum += CI(j, i) * x(j);
      sum += ci0(i);
      s(i) = sum;
      psi += std::min(0.0, sum);
    }
    QUADPROG_TRACE_VECTOR("s", s);

    // TODO: check tolerance scaling
    if (fabs(psi) <= m * options.tolerance * c1 * c2)
    {
      /* numerically there are not infeasibilities anymore */
      result.status = SolverStatus::SUCCESS;
      result.objective_value = f_value;
      result.solution = x;
      result.iterations = iter;
      result.message = "Optimization successful";
      QUADPROG_TRACE("Optimization successful in {} iterations, objective value {}", iter, f_value);
      return result;
    }

    /* save old values for u and A */
    for (i = 0; i < iq; i++)
    {
      u_old(i) = u(i);
      A_old(i) = A(i);
    }
    /* and for x */
    for (i = 0; i < n; i++)
      x_old(i) = x(i);

  l2: /* Step 2: check for feasibility and determine a new S-pair */
    for (i = 0; i < m; i++)
    {
      if (s(i) < ss && iai(i) != -1 && iaexcl(i))
      {
        ss = s(i);
        ip = i;
      }
    }
    if (ss >= 0.0)
    {
      /* all constraints are satisfied */
      result.status = SolverStatus::SUCCESS;
      result.objective_value = f_value;
      result.solution = x;
      result.iterations = iter;
      result.message = "Optimization successful";
      QUADPROG_TRACE("Optimization successful in {} iterations, objective value {}", iter, f_value);
      return result;
    }

    /* set np = n[ip] */
    for (i = 0; i < n; i++)
      np(i) = CI(i, ip);
    /* set u = [u 0]^T */
    u(iq) = 0.0;
    /* add ip to the active set A */
    A(iq) = ip;
    QUADPROG_TRACE("Trying with constraint {}", ip);
    QUADPROG_TRACE_VECTOR("np", np);

  l2a: /* Step 2a: determine step direction */
    /* compute z = H np: the step direction in the primal space (through J, see the paper) */
    compute_d(d, J, np);
    update_z(z, J, d, iq);
    /* compute N* np (if q > 0): the negative of the step direction in the dual space */
    update_r(R, r, d, iq);

    QUADPROG_TRACE_VECTOR("Step direction z", z);
    QUADPROG_TRACE_VECTOR("r", r, iq + 1);
    QUADPROG_TRACE_VECTOR("u", u, iq + 1);
    QUADPROG_TRACE_VECTOR("d", d);
    QUADPROG_TRACE_VECTOR("A", A, iq + 1);

    /* Step 2b: compute step lengths */
    l = 0;
    /* Compute t1: partial step length (maximum step in dual space without violating dual feasibility */
    t1 = inf;
    /* find the index l s.t. it reaches the minimum of u+[x] / r */
    for (k = p; k < iq; k++)
    {
      if (r(k) > 0.0)
      {
        if (u(k) / r(k) < t1)
        {
          t1 = u(k) / r(k);
          l = A(k);
        }
      }
    }
    /* Compute t2: full step length (minimum step in primal space such that the constraint ip becomes feasible */
    if (fabs(scalar_product<double>(z, z)) > std::numeric_limits<double>::epsilon()) // i.e. z != 0
    {
      t2 = -s(ip) / scalar_product<double>(z, np);
      if (t2 < 0) // patch suggested by Takano Akio for handling numerical inconsistencies
        t2 = inf;
    }
    else
      t2 = inf; /* +inf */

    /* the step is chosen as the minimum of t1 and t2 */
    t = std::min(t1, t2);
    QUADPROG_TRACE("Step sizes: {} (t1 = {}, t2 = {}) ", t, t1, t2);

    /* Step 2c: determine new S-pair and take step: */
    /* case (i): step in dual space only */
    if (t >= inf)
    {
      /* QPP is unbounded */
      result.status = SolverStatus::UNBOUNDED;
      result.message = "Quadratic program is unbounded";
      result.objective_value = std::numeric_limits<double>::infinity();
      QUADPROG_TRACE("Quadratic program is unbounded");
      return result;
    }
    /* case (ii): step in dual space */
    if (t2 >= inf)
    {
      /* set u = u +  t * [-r 1] and drop constraint l from the active set A */
      for (k = 0; k < iq; k++)
        u(k) -= t * r(k);
      u(iq) += t;
      iai(l) = l;
      delete_constraint(R, J, A, u, n, p, iq, l);
      // FIXME: check the logging
      //    QUADPROG_TRACE("Deleting constraint {}", l);
      QUADPROG_TRACE(" in dual space: {}", f_value);
      QUADPROG_TRACE_VECTOR("x", x);
      QUADPROG_TRACE_VECTOR("z", z);
      QUADPROG_TRACE_VECTOR("A", A, iq + 1);

      goto l2a;
    }

    /* case (iii): step in primal and dual space */

    /* set x = x + t * z */
    for (k = 0; k < n; k++)
      x(k) += t * z(k);
    /* update the solution value */
    f_value += t * scalar_product<double>(z, np) * (0.5 * t + u(iq));
    /* u = u + t * [-r 1] */
    for (k = 0; k < iq; k++)
      u(k) -= t * r(k);
    u(iq) += t;
    QUADPROG_TRACE(" in both spaces: {}", f_value);
    QUADPROG_TRACE_VECTOR("x", x);
    QUADPROG_TRACE_VECTOR("u", u, iq + 1);
    QUADPROG_TRACE_VECTOR("r", r, iq + 1);
    QUADPROG_TRACE_VECTOR("A", A, iq + 1);

    if (fabs(t - t2) < std::numeric_limits<double>::epsilon())
    {
      QUADPROG_TRACE("Full step has taken {}", t);
      QUADPROG_TRACE_VECTOR("x", x);
      /* full step has taken */
      /* add constraint ip to the active set*/
      if (!add_constraint(R, J, d, iq, R_norm))
      {
        iaexcl(ip) = false;
        delete_constraint(R, J, A, u, n, p, iq, ip);
        QUADPROG_TRACE_MATRIX("R", R);
        QUADPROG_TRACE_VECTOR("A", A, iq);
        QUADPROG_TRACE_VECTOR("iai", iai);
        for (i = 0; i < m; i++)
          iai(i) = i;
        for (i = p; i < iq; i++)
        {
          A(i) = A_old(i);
          u(i) = u_old(i);
          iai(A(i)) = -1;
        }
        for (i = 0; i < n; i++)
          x(i) = x_old(i);
        goto l2; /* go to step 2 */
      }
      else
        iai(ip) = -1;
      QUADPROG_TRACE_MATRIX("R", R);
      QUADPROG_TRACE_VECTOR("A", A, iq);
      QUADPROG_TRACE_VECTOR("iai", iai);
      goto l1;
    }

    /* a patial step has taken */
    QUADPROG_TRACE("Partial step has taken {}", t);
    QUADPROG_TRACE_VECTOR("x", x);
    /* drop constraint l */
    iai(l) = l;
    delete_constraint(R, J, A, u, n, p, iq, l);
    QUADPROG_TRACE_MATRIX("R", R);
    QUADPROG_TRACE_VECTOR("A", A, iq);

    /* update s[ip] = CI * x + ci0 */
    sum = 0.0;
    for (k = 0; k < n; k++)
      sum += CI(k, ip) * x(k);
    s(ip) = sum + ci0(ip);

    QUADPROG_TRACE_VECTOR("s", s, m);

    goto l2a;
  }

  double solve_quadprog_legacy(
      Matrix<double> &G,
      Vector<double> &g0,
      const Matrix<double> &CE,
      const Vector<double> &ce0,
      const Matrix<double> &CI,
      const Vector<double> &ci0,
      Vector<double> &x)
  {
    auto result = solve_quadprog(G, g0, CE, ce0, CI, ci0);

    if (result.is_success())
    {
      x = result.solution;
      return result.objective_value;
    }
    else
    {
      return std::numeric_limits<double>::infinity();
    }
  }

  inline void compute_d(Vector<double> &d, const Matrix<double> &J, const Vector<double> &np)
  {
    size_t n = d.size();
    double sum;

    /* compute d = H^T * np */
    for (size_t i = 0; i < n; i++)
    {
      sum = 0.0;
      for (size_t j = 0; j < n; j++)
        sum += J(j, i) * np(j);
      d(i) = sum;
    }
  }

  inline void update_z(Vector<double> &z, const Matrix<double> &J, const Vector<double> &d, size_t iq)
  {
    size_t n = z.size();

    /* setting of z = H * d */
    for (size_t i = 0; i < n; i++)
    {
      z(i) = 0.0;
      for (size_t j = iq; j < n; j++)
        z(i) += J(i, j) * d(j);
    }
  }

  inline void update_r(const Matrix<double> &R, Vector<double> &r, const Vector<double> &d, size_t iq)
  {
    double sum;

    /* setting of r = R^-1 d */
    for (size_t i = iq; i-- > 0;) // go backwards from iq-1 to 0
//    for (int i = iq - 1; i >= 0; i--)
    {
      sum = 0.0;
      for (size_t j = i + 1; j < iq; j++)
        sum += R(i, j) * r(j);
      r(i) = (d(i) - sum) / R(i, i);
    }
  }

  bool add_constraint(Matrix<double> &R, Matrix<double> &J, Vector<double> &d, size_t &iq, double &R_norm)
  {
    size_t n = d.size();
    QUADPROG_TRACE("Adding constraint at position {}", iq);
    double cc, ss, h, t1, t2, xny;

    /* we have to find the Givens rotation which will reduce the element
      d[j] to zero.
      if it is already zero we don't have to do anything, except of
      decreasing j */
    for (size_t j = n; j-- > iq + 1;) // go backwards from n-1 to iq+1
//    for (int j = n - 1; j >= static_cast<int>(iq + 1); j--)
    {
      /* The Givens rotation is done with the matrix (cc cs, cs -cc).
      If cc is one, then element (j) of d is zero compared with element
      (j - 1). Hence we don't have to do anything.
      If cc is zero, then we just have to switch column (j) and column (j - 1)
      of J. Since we only switch columns in J, we have to be careful how we
      update d depending on the sign of gs.
      Otherwise we have to apply the Givens rotation to these columns.
      The i - 1 element of d has to be updated to h. */
      cc = d(j - 1);
      ss = d(j);
      h = distance(cc, ss);
      if (std::fabs(h) < std::numeric_limits<double>::epsilon()) // h == 0
        continue;
      d(j) = 0.0;
      ss = ss / h;
      cc = cc / h;
      if (cc < 0.0)
      {
        cc = -cc;
        ss = -ss;
        d(j - 1) = -h;
      }
      else
        d(j - 1) = h;
      xny = ss / (1.0 + cc);
      for (size_t k = 0; k < n; k++)
      {
        t1 = J(k, j - 1);
        t2 = J(k, j);
        J(k, j - 1) = t1 * cc + t2 * ss;
        J(k, j) = xny * (t1 + J(k, j - 1)) - t2;
      }
    }
    /* update the number of constraints added*/
    iq++;
    /* To update R we have to put the iq components of the d vector
      into column iq - 1 of R
      */
    for (size_t i = 0; i < iq; i++)
      R(i, iq - 1) = d(i);
    QUADPROG_TRACE_MATRIX("R", R, iq, iq);
    QUADPROG_TRACE_MATRIX("J", J);
    QUADPROG_TRACE_VECTOR("d", d, iq);

    if (fabs(d(iq - 1)) <= std::numeric_limits<double>::epsilon() * R_norm)
    {
      // degenerate problem
      return false;
    }
    R_norm = std::max<double>(R_norm, fabs(d(iq - 1)));
    return true;
  }

  void delete_constraint(Matrix<double> &R, Matrix<double> &J, Vector<int> &A, Vector<double> &u, size_t n, size_t p, size_t &iq, int l)
  {
    QUADPROG_TRACE("Deleting constraint at position {}", l);
    size_t qq = 0; // just to prevent warnings from smart compilers
    double cc, ss, h, xny, t1, t2;

    bool found = false;
    /* Find the index qq for active constraint l to be removed */
    for (size_t i = p; i < iq; i++)
      if (A(i) == l)
      {
        qq = i;
        found = true;
        break;
      }

    if (!found)
    {
      QUADPROG_TRACE("Attempt to delete non existing constraint {}", l);
      throw std::invalid_argument("Attempt to delete non existing constraint " + std::to_string(l));
    }
    /* remove the constraint from the active set and the duals */
    for (size_t i = qq; i < iq - 1; i++)
    {
      A(i) = A(i + 1);
      u(i) = u(i + 1);
      for (size_t j = 0; j < n; j++)
        R(j, i) = R(j, i + 1);
    }

    A(iq - 1) = A(iq);
    u(iq - 1) = u(iq);
    A(iq) = 0;
    u(iq) = 0.0;
    for (size_t j = 0; j < iq; j++)
      R(j, iq - 1) = 0.0;
    /* constraint has been fully removed */
    iq--;
    QUADPROG_TRACE("{} constraints remain after deletion", iq);

    if (iq == 0) // no constraints left
      return;

    for (size_t j = qq; j < iq; j++)
    {
      cc = R(j, j);
      ss = R(j + 1, j);
      h = distance(cc, ss);
      if (fabs(h) < std::numeric_limits<double>::epsilon()) // h == 0
        continue;
      cc = cc / h;
      ss = ss / h;
      R(j + 1, j) = 0.0;
      if (cc < 0.0)
      {
        R(j, j) = -h;
        cc = -cc;
        ss = -ss;
      }
      else
        R(j, j) = h;

      xny = ss / (1.0 + cc);
      for (size_t k = j + 1; k < iq; k++)
      {
        t1 = R(j, k);
        t2 = R(j + 1, k);
        R(j, k) = t1 * cc + t2 * ss;
        R(j + 1, k) = xny * (t1 + R(j, k)) - t2;
      }
      for (size_t k = 0; k < n; k++)
      {
        t1 = J(k, j);
        t2 = J(k, j + 1);
        J(k, j) = t1 * cc + t2 * ss;
        J(k, j + 1) = xny * (J(k, j) + t1) - t2;
      }
    }
  }

  inline double distance(double a, double b)
  {
    double a1, b1, t;
    a1 = std::fabs(a);
    b1 = std::fabs(b);
    if (a1 > b1)
    {
      t = (b1 / a1);
      return a1 * std::sqrt(1.0 + t * t);
    }
    else if (b1 > a1)
    {
      t = (a1 / b1);
      return b1 * std::sqrt(1.0 + t * t);
    }
    return a1 * std::sqrt(2.0);
  }

} // namespace quadprog
