#ifndef _TESTSOLVABLE_H_
#define _TESTSOLVABLE_H_

#include "Solvable.hpp"

class TestSolvable : public Solvable {
 public:
  TestSolvable(int ndof, const VectorXs & lengths, const std::vector<std::pair<int, int>>& pairs,
               const MatrixXXs& CL, const VectorXs& dL) {
    assert(CL.rows() == dL.rows());
    assert(lengths.rows() == ndof-1);
    m_l = lengths;
    m_ndof = ndof;
    m_pairs = pairs;
    m_dL    = dL;
    // convert dense np->eigenmat to sparse CL
    std::vector<Triplet> triplets;
    triplets.reserve(CL.rows() * CL.cols());
    for (int i = 0; i < CL.rows(); i++) {
      for (int j = 0; j < CL.cols(); j++) {
        if (abs(CL(i, j)) > EPSILON) {
          triplets.emplace_back(i, j, CL(i, j));
        }
      }
    }
    m_CL.resize(CL.rows(), CL.cols());
    m_CL.setFromTriplets(triplets.begin(), triplets.end());
  }
  ~TestSolvable() {}

  std::vector<std::pair<int, int>> getPeriodicPairs() { return m_pairs; };
  std::pair<SparseXXs, VectorXs> getLagrangeConstraints() {
    return std::make_pair(m_CL, m_dL);
  }

  int getNumDof() { return m_ndof; }

  scalar val(const VectorXs& q) {  // f(q)
    scalar E = 0;
    // sum q_i^2
    // for (int i = 0; i < q.rows(); i++) {
    //     E += q(i) * q(i);
    // }
    // sum (q_i+1 - q_i)^2
    for (int i = 0; i < q.rows() - 1; i++) {
      E += (q(i + 1) - q(i) - m_l(i)) * (q(i + 1) - q(i) - m_l(i));
    }
    return E;
  }

  VectorXs grad(const VectorXs& q) {  // df/dq (q)
    VectorXs gradE = VectorXs::Zero(q.rows());
    // 2 q_j
    // for (int i = 0; i < q.rows(); i++) {
    //     gradE(i) += 2 * q(i);
    // }
    // gradE_i += -2 (q_i+1 - q_i), gradE_i+1 += 2 (q_i+1 - q_i)
    for (int i = 0; i < q.rows() - 1; i++) {
      gradE(i) += -2 * (q(i + 1) - q(i) - m_l(i));
      gradE(i + 1) += 2 * (q(i + 1) - q(i) - m_l(i));
    }
    return gradE;
  }

  SparseXXs hess(const VectorXs& q) {  // d2f/dqdq (q)
    std::vector<Triplet> triplets;
    triplets.reserve(4 * q.rows());
    // 2 dij
    // for (int i = 0; i < q.rows(); i++) {
    //     for (int j = 0; j < q.rows(); j++) {
    //         if (i == j)
    //             triplets.emplace_back(i,j, 2);
    //     }
    // }
    // hessE_ii += 2, hessE_i+1i+1 += 2, hessE_ii+1 += -2
    for (int i = 0; i < q.rows() - 1; i++) {
      triplets.emplace_back(i, i, 2);
      triplets.emplace_back(i, i + 1, -2);
      triplets.emplace_back(i + 1, i, -2);
      triplets.emplace_back(i + 1, i + 1, 2);
    }
    SparseXXs hessE(q.rows(), q.rows());
    hessE.setFromTriplets(triplets.begin(), triplets.end());
    return hessE;
  }

 private:
  std::vector<std::pair<int, int>> m_pairs;
  SparseXXs m_CL;
  VectorXs m_dL;
  VectorXs m_l;
  int m_ndof;
};

#endif /* end of include guard: _TESTSOLVABLE_H_ */