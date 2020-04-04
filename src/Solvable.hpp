#ifndef _SOLVABLE_H_
#define _SOLVABLE_H_

#include <utility>
#include <vector>
#include "EigenDefinitions.hpp"

class Solvable {
  // f(q): R^N -> R
  // for finding min_q f <--- grad_q f = 0 using Newton's method
  // with constraints q+ = q- and CL q = dL
 public:
  Solvable() {}
  virtual ~Solvable() {}

  virtual int getNumDof() = 0;
  virtual std::vector<std::pair<int, int>> getPeriodicPairs()     = 0;
  virtual std::pair<SparseXXs, VectorXs> getLagrangeConstraints() = 0;

  virtual scalar val(const VectorXs& q)     = 0;  // f(q)
  virtual VectorXs grad(const VectorXs& q)  = 0;  // df/dq (q)
  virtual SparseXXs hess(const VectorXs& q) = 0;  // d2f/dqdq (q)

  virtual scalar computeMaximumStep(const VectorXs& dq) { return 0; }
};

#endif /* end of include guard: _SOLVABLE_H_ */
