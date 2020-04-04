#include "NewtonSolver.hpp"
#include <iostream>
#include "debug/debug_tools.hpp"

NewtonSolver::NewtonSolver(std::shared_ptr<Solvable> solvable, NewtonSolverSettings settings)
    : m_settings(settings), m_solvable(solvable) {
  initialize();
}

void NewtonSolver::initialize() {
  m_finished = false;
  m_iterations = 0;

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

  m_Eprev = m_solvable->val(m_Ctilde * m_y);

  VectorXs dEdy = m_Ctilde.transpose() * m_solvable->grad(m_Ctilde * m_y);
  m_z0 = computeProjectedGradientNorm(dEdy);
}


scalar NewtonSolver::computeProjectedGradientNorm(const VectorXs & dEdy) {
  // gradient norm of free variables y in tangent space of lagrange constraint CL q = dL:
  // A A.T lam = A z
  // |P(z)| = sqrt(z.T z - lam.T A z)
  // where A = CL Ctilde, z = grad_y E

  auto & z = dEdy;
  auto & A = m_CLCt;
  VectorXs Az = A * z;
  VectorXs lam = m_CCCCsolver.solve(Az);
  if (m_CCCCsolver.info() != Eigen::Success) {
    Debug::error("Could not solve lagrange projection!");
    return false;
  }

  scalar znorm = std::sqrt(z.squaredNorm() - lam.dot(Az));
  return znorm;
}


VectorXs NewtonSolver::computeSearchDirection() {
  // | Ctilde.T H Ctilde + aI    Ctilde.T CL.T | | dy  |   | - Ctilde.T dEdq |
  // |     CL Ctilde                   0       | | lam | = |  dL - CL dtilde |

  VectorXs q = m_Ctilde * m_y;
  VectorXs dEdy = m_Ctilde.transpose() * m_solvable->grad(q);
  VectorXs residual = m_dL - m_CLCt * m_y; // CL (Ctilde y + dtilde) - dL = 0

  int nfree = m_Ctilde.cols();
  int neqs = m_dL.rows();

  // rhs = -[dEdy; CL Ct y - dL]
  VectorXs rhs(nfree + neqs);
  rhs.head(nfree) = -dEdy;
  rhs.tail(neqs) = residual;

  // stopping criterion
  scalar znorm = computeProjectedGradientNorm(dEdy);
  bool crit_projgrad = znorm < m_z0 * m_settings.projgrad_epsilon;
  bool crit_lagrange = residual.norm() < m_settings.lagrange_epsilon;
  assert(crit_lagrange); // given initial projection and consistent Newton updates 
  // Debug::logf("Stop crit (proj. grad): %.2e %.2e %s\n",
  // znorm, znorm / m_z0, crit_projgrad ? "True" : "False");
  if (crit_projgrad && crit_lagrange) {
    m_finished = true;
    Debug::log("Stopping criterion fulfilled.");
    return VectorXs::Zero(m_y.rows());
  }

  // construct LHS
  SparseXXs ClTHCt = m_Ctilde.transpose() * m_solvable->hess(q) * m_Ctilde;
  std::vector<Triplet> triplets;
  triplets.reserve(ClTHCt.nonZeros() + 2 * m_CLCt.nonZeros() + nfree);

  // add top left
  for (int k = 0; k < ClTHCt.outerSize(); ++k)
    for (SparseXXs::InnerIterator it(ClTHCt, k); it; ++it) {
      if (std::abs(it.value()) > EPSILON)
        triplets.emplace_back(it.row(), it.col(), it.value());
    }
  // add top left regularizer: alpha I
  // logarithmically decaying regularizer
  scalar alpha_hat = m_settings.reg_start * std::pow(m_settings.reg_end/m_settings.reg_start,std::min(std::max(0.0,m_iterations * 1.0 / m_settings.reg_steps),1.0));
  scalar alpha = std::sqrt(rhs.cwiseAbs2().maxCoeff()) * alpha_hat;
  // alpha=1; alpha = std::max(alpha, 1e-5);
  assert(!(alpha < 0));
  if (alpha > EPSILON) {
    for (int i = 0; i < nfree; i++)
      triplets.emplace_back(i, i, alpha);
  }
  // add top right, bottom left
  for (int k = 0; k < m_CLCt.outerSize(); ++k)
    for (SparseXXs::InnerIterator it(m_CLCt, k); it; ++it) {
      triplets.emplace_back(nfree + it.row(), it.col(), it.value()); // BL
      triplets.emplace_back(it.col(), nfree + it.row(), it.value()); // TR
    }

  SparseXXs LHS(nfree+neqs, nfree+neqs);
  LHS.setFromTriplets(triplets.begin(), triplets.end());

  // TODO note that could iterate increasing regularization for linear solve
  // if so then wrap this..
  // NOTE also could fail when we are at a saddlepoint and grad == 0 -> no reg 
  // e.g. at perfect solution!

  Eigen::SimplicialLDLT<SparseXXs> solver;
  // factorize
  solver.compute(LHS);
  if (solver.info() != Eigen::Success) {
    Debug::errorf("Linear solver failed to factorize!\n");

    Debug::log("\n",LHS.toDense());
    // continue; // NOTE could incrementally add more regularization
    return VectorXs::Zero(nfree);
  }
  VectorXs dy = solver.solve(rhs);
  if (solver.info() != Eigen::Success) {
    Debug::errorf("Linear solver failed to solve!\n");
    // continue; // NOTE could incrementally add more regularization
    return VectorXs::Zero(nfree);
  }


  #ifndef NDEBUG
  // r.x = x.A.x should be > 0
  scalar dot = dy.dot(rhs); 
  if (dot <= 0.0) {
    if (dot == 0.0) {
      Debug::warningf("System was semidefinite!\n");
    } else {
      Debug::errorf("System was indefinite!\n");
      // continue; // NOTE could incrementally add more regularization
      return VectorXs::Zero(nfree);
    }
  }
  // For safety and mistrust of Eigen and measure error
  // |Ax - b|^2 / |b|^2 < max_error^2
  if ((LHS.selfadjointView<Eigen::Lower>() * dy - rhs).squaredNorm() >
      rhs.squaredNorm() * m_settings.max_linear_error *
          m_settings.max_linear_error * 1.2) { // *1.2 to be conservative wrt iterative solver error threshold
    Debug::errorf("Linear solve error above threshold!\n");
    // continue; // NOTE could incrementally add more regularization
    return VectorXs::Zero(nfree);
  }
  #endif

  return dy.head(nfree);
}

void NewtonSolver::backtrackingLinesearch(const VectorXs & dy) {
  VectorXs y0 = m_y;
  for (int i = 0; i < 10; i++) {
    m_y = y0 + std::pow(0.1,i) * dy; // step
    scalar E = m_solvable->val(m_Ctilde * m_y); // check energy
    if (E < m_Eprev) {
      m_Eprev = E;
      return; // success
    }
  }

  // simple linesearch failed
  // for now defaulting to aborting the solve
  m_finished = true;
  m_y = y0;
  Debug::log("Aborting because linesearch did not improve.");
}

void NewtonSolver::step() {
  // solve linear system for update to free variables
  // and also check stopping criterion
  VectorXs dy = computeSearchDirection();
  if (m_finished)
    return;

  if (m_settings.step_limit > 0) {
    scalar maxstep = m_solvable->computeMaximumStep(m_Ctilde * dy);
    if (maxstep > m_settings.step_limit)
      dy *= m_settings.step_limit / maxstep;
  }

  // simple linesearch for multiples of dy
  backtrackingLinesearch(dy);

  m_iterations++;
}

VectorXs NewtonSolver::reproject(const VectorXs& q) {
  // min y |Ctilde y + dtilde - q| st CL (Ctilde y + dtilde) = dL

  // | Ctilde.T Ctilde    Ctilde.T CL.T | | y   |   | Ctilde.T (q - dtilde) |
  // |     CL Ctilde         0          | | lam | = |     dL - CL dtilde    |

  // int ndof  = m_solvable->getNumDof();
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
  for (int k = 0; k < m_CLCt.outerSize(); ++k)
    for (SparseXXs::InnerIterator it(m_CLCt, k); it; ++it) {
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

  // periodicity is allowed as 1->Many mapping
  //   i.e. one master can have multiple copies
  //   but we assume that no copy occurs twice (i.e copy cant be a copy of multiple masters, and copy cannot itself be a master)
  // For simplicity, we assume that pairs are given as {master, copy}

  // def i -> ix if i in pair[ix] else -1
  int npairs        = (int)pairs.size();
  VectorXs i2reducedi = VectorXs::Constant(ndof, -1);
  int reducedi = 0;
  for (int ix = 0; ix < npairs; ix++) {
    auto ij = pairs[ix];
    assert(ij.first != ij.second);  // not self-periodic
    assert(i2reducedi[ij.second] == -1);  // non-duplicate, i.e. copies are not pre-assigned
    if (i2reducedi[ij.first] == -1) { // unassigned master
      i2reducedi[ij.first]  = reducedi;
      i2reducedi[ij.second] = reducedi;
      reducedi++;
    } else { // master already exists with an assigned reduced index
      int preassigned_reducedi = i2reducedi[ij.first];
      i2reducedi[ij.second] = preassigned_reducedi;
    }
  }
  // now append the rest (non-periodic verts)
  for (int i = 0; i < ndof; i++) {
    if (i2reducedi[i] != -1) // already assigned periodic verts
      continue;
    i2reducedi[i] = reducedi;
    reducedi++;
  }

  int nreduced = reducedi;
  SparseXXs Ctilde(ndof, nreduced);
  assert(ndof >= nreduced);

  std::vector<Triplet> triplets;
  triplets.reserve(ndof);

  for (int i = 0; i < ndof; i++) {
    int j = i2reducedi[i];
    assert(j != -1);
    triplets.emplace_back(i, j, 1);
  }

  Ctilde.setFromTriplets(triplets.begin(), triplets.end());
  return Ctilde;
}
