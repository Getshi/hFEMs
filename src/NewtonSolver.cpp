#include "NewtonSolver.hpp"
#include <iostream>
#include "debug/debug_tools.hpp"

NewtonSolver::NewtonSolver(std::shared_ptr<Solvable> solvable)
    : m_solvable(solvable) {
  initialize();
}

void NewtonSolver::initialize() {
  m_finished = false;

  // constraints CL q = dL solved using Lagrange multipliers
  SparseXXs CL;
  std::tie(CL, m_dL) = m_solvable->getLagrangeConstraints();

  int ndof = m_solvable->getNumDof();  // get ndof

  // periodicity
  auto pairs = m_solvable->getPeriodicPairs();  // {(i,j)}  -> q_i = q_j
  m_Ctilde   = generateCtilde(
      pairs, ndof);  // NOTE: dtilde is 0 for simplified periodicity

  // precompute CL Ctilde and factorize (CL Ct Ct.T CL.T)
  m_CLCt         = (CL * m_Ctilde).pruned();
  SparseXXs CCCC = m_CLCt * m_CLCt.transpose();
  m_CCCCsolver.compute(CCCC);
  if (m_CCCCsolver.info() != Eigen::Success) {
    Debug::error("Warning! CCCC-QR factorization failed!");
  }

  // initial guess q=0
  VectorXs q0 = VectorXs::Zero(ndof);

  // project initial guess onto free variables s.t. CL Ctilde y = dL
  // NOTE: done for generality. for the special case of initial q=0 and dL=0, it
  // would suffice to set y=0
  m_y = reproject(q0);

  scalar m_Eprev = m_solvable->val(m_Ctilde * m_y);

  // TODO
  // z0 = stoppingconditionvalue()
}

void NewtonSolver::step() {
  // TODO

  // compute dy:
  //   dEdy = Ct.T @ solvable.gradE
  //   rhs = [-dEdy; CL Ct y - dL]
  //   check stopping criterion:
  //     TODO
  //     finished = true
  //   ddEdydy = Ct.T @ solvable.hessE @ Ct
  //   LHS = [ddEdydy + alpha I, Ct.T CL.T; CL Ct, 0]
  //   solve dy
  // linesearch:
  //   y0 = y
  //   for i = 0 .. 10
  //     y = y0 + 0.1^i dy
  //     solvable.set(Ct y)
  //     E = solvable
  //     if E < E0:
  //       E0 = E
  //       break
  //   if E > E0:
  //     abort
  //     finished = true
}

VectorXs NewtonSolver::reproject(const VectorXs& q) {
  // min y |Ctilde y + dtilde - q| st CL (Ctilde y + dtilde) = dL

  // | Ctilde.T Ctilde    Ctilde.T CL.T | | y   |   | Ctilde.T (q - dtilde) |
  // |     CL Ctilde         0          | | lam | = |     dL - CL dtilde    |

  int ndof  = m_solvable->getNumDof();
  int nfree = m_Ctilde.cols();
  int neqs  = m_dL.rows();

  SparseXXs CtTCt = m_Ctilde.transpose() * m_Ctilde;

  std::vector<Triplet> triplets;
  triplets.reserve(CtTCt.nonZeros() + 2 * m_CLCt.nonZeros());

  // add top left
  for (int k = 0; k < CtTCt.outerSize(); ++k)
    for (SparseXXs::InnerIterator it(CtTCt, k); it; ++it) {
      triplets.emplace_back(it.row(), it.col(), it.value());
    }
  // add top right, bottom left
  for (int k = 0; k < CtTCt.outerSize(); ++k)
    for (SparseXXs::InnerIterator it(CtTCt, k); it; ++it) {
      triplets.emplace_back(nfree + it.row(), it.col(), it.value()); // BL
      triplets.emplace_back(it.col(), nfree + it.row(), it.value()); // TR
    }

  SparseXXs LHS(nfree+neqs, nfree+neqs);
  LHS.setFromTriplets(triplets.begin(), triplets.end());

  VectorXs rhs(nfree + neqs);
  rhs.head(nfree) = m_Ctilde.transpose() * q;  // Ctilde.T (q-dtilde); NOTE: dtilde = 0
  rhs.tail(neqs) = m_dL; // dL - CL dtilde

  Eigen::LeastSquaresConjugateGradient<SparseXXs> solver;
  // Eigen::SimplicialLDLT<SparseXXs> solver;
  solver.compute(LHS);
  VectorXs y         = solver.solve(rhs);
  bool system_solved = solver.info() == Eigen::Success;
  if (!system_solved)
    Debug::error("Couldnt solve projection system!");

  return y.head(nfree);
}

// non-member functions

// generate Ctilde for simple periodicity q_i = q_j, given pairs (i,j)
SparseXXs generateCtilde(const std::vector<std::pair<int, int>>& pairs,
                         int ndof) {
  // q = Ctilde y; thus Ctilde_ij takes reduced indices j to full indices i
  // define the reduced order with pairs p as [p0, p1, ..., pn, other
  // non-paired... ]

  // def i -> ix if i in pair[ix] else -1
  int npairs        = (int)pairs.size();
  VectorXs i2pairix = VectorXs::Constant(ndof, -1);
  for (int ix = 0; ix < npairs; ix++) {
    auto ij = pairs[ix];
    assert(ij.first != ij.second);  // not self-periodic
    assert(i2pairix[ij.first] == -1 &&
           i2pairix[ij.second] == -1);  // non-duplicate, i.e. not pre-assigned
    i2pairix[ij.first]  = ix;
    i2pairix[ij.second] = ix;
  }

  int nreduced = ndof - npairs;
  SparseXXs Ctilde(ndof, nreduced);
  std::vector<Triplet> triplets;
  triplets.reserve(ndof);

  int nextfree = npairs;
  for (int i = 0; i < ndof; i++) {
    // map j to i
    int j, ix;
    ix = i2pairix[i];
    if (ix < 0) {
      j = nextfree;
      nextfree++;
    } else {
      j = ix;
    }
    triplets.emplace_back(i, j, 1);
  }

  Ctilde.setFromTriplets(triplets.begin(), triplets.end());
}
