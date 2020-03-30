#ifndef _NEWTONSOLVER_H_
#define _NEWTONSOLVER_H_

#include <memory>
#include "Solvable.hpp"

// Solve min_q f(q) subject to q+=q- and CL q = dL
class NewtonSolver {
 public:
  NewtonSolver(std::shared_ptr<Solvable> solvable);

  void initialize();
  void step();
  bool isFinished() { return m_finished; }

 private:
  std::shared_ptr<Solvable> m_solvable;
  bool m_finished;
  SparseXXs m_Ctilde, m_CLCt;
  VectorXs m_y;
  VectorXs m_dL;
  scalar m_Eprev;  // initial/previous energy
  Eigen::SparseQR<SparseXXs, Eigen::COLAMDOrdering<int>> m_CCCCsolver;

  VectorXs reproject(const VectorXs& q);
};

// non-member functions
SparseXXs generateCtilde(const std::vector<std::pair<int, int>>& pairs,
                         int ndof);

#endif /* end of include guard: _NEWTONSOLVER_H_ */
