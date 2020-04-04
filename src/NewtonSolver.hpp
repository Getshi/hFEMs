#ifndef _NEWTONSOLVER_H_
#define _NEWTONSOLVER_H_

#include <memory>
#include "Solvable.hpp"

struct NewtonSolverSettings {
  // regularization
  scalar reg_start = 5e4;
  scalar reg_end = 5e0;
  int reg_steps = 400;
  scalar step_limit = 0;
  // stopping criterion
  scalar projgrad_epsilon = 1e-8;//1e-5; // |z| < eps * |z0|
  scalar lagrange_epsilon = 1e-8; // |Cq-d| < eps
  #ifndef NDEBUG
  scalar max_linear_error = 1e-2; // debug or iterative solver TODO
  #endif
};

// Solve min_q f(q) subject to q+=q- and CL q = dL
class NewtonSolver {
 public:
  NewtonSolver(std::shared_ptr<Solvable> solvable, NewtonSolverSettings settings = NewtonSolverSettings());

  void initialize();
  void step();
  bool isFinished() { return m_finished; }
  VectorXs getSolution() { return m_Ctilde * m_y; }

 private:
  NewtonSolverSettings m_settings;
  std::shared_ptr<Solvable> m_solvable;
  bool m_finished;
  SparseXXs m_Ctilde, m_CLCt;
  VectorXs m_y;
  VectorXs m_dL;
  int m_iterations;
  // stopping criterion
  Eigen::SparseQR<SparseXXs, Eigen::COLAMDOrdering<int>> m_CCCCsolver;
  scalar m_z0; // initial norm of proj grad
  scalar m_Eprev;  // initial/previous energy


  VectorXs reproject(const VectorXs& q);
  scalar computeProjectedGradientNorm(const VectorXs& dEdy);
  VectorXs computeSearchDirection();
  void backtrackingLinesearch(const VectorXs & dy);
};

// non-member functions
SparseXXs generateCtilde(const std::vector<std::pair<int, int>>& pairs,
                         int ndof);

#endif /* end of include guard: _NEWTONSOLVER_H_ */
