#include <memory>
#include "../src/TestSolvable.hpp"
#include "../src/VegaFEMSolvable.hpp"
#include "../src/NewtonSolver.hpp"
#include "../src/debug/debug_tools.hpp"


/*
// TODO
// expose PY
// check
//  grads and hess
//  check initial guess for various defo
//  check solved state for various defo
// --- try adding step limit
*/

int main (int argc, char *argv[]) {
  Debug::log("Test start.");

  // // DEBUG Expansion
  // Vector6s strain;
  // strain << 0.3,0.6,0, -0.6,0.2,0.3;
  // Vector2s min, max;
  // min << -1.0, -2.0;
  // max << 3.0, 3.0;
  // scalar dx = 0.1;
  // Expansion expansion(strain, min, max, dx);
  // int nrows,ncols;
  // const MatrixXXs* phi;
  // std::tie(nrows,ncols,phi) = expansion. getGrid();
  // std::cout<<"nrows, ncols = " << nrows <<","<<ncols<<"\n";
  // Eigen::IOFormat eigenFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "[", "]", "[", "]");
  // std::cout<<"V = np.array("<<phi->format(eigenFmt)<<",dtype=np.float)\n";

  // DEBUG VegaFEMSolvable

  std::string filename = "/home/gsperl/Projects/hFEMs/meshes/testmesh.veg";
  Vector6s strain;
  // strain << 0,0,0,0,0,0;  // sx sa sy II00 II01 II11
  strain << 0.3, 0.0, 0.0, 0.2, 0.0, 0.0;  // sx sa sy II00 II01 II11

  auto solvable = std::make_shared<VegaFEMSolvable>(filename, strain);
  int ndof = solvable->getNumDof();
  Debug::log("Energy at start:", solvable->val(VectorXs::Zero(ndof)));

  NewtonSolverSettings settings;
  settings.projgrad_epsilon=0.1; // TODO absolute stop. crit?
  auto solver = std::make_shared<NewtonSolver>(solvable, settings);


  // // DEBUG TestSolvable OLD
  // int ndof = 10;
  // VectorXs l = VectorXs::Zero(ndof-1);
  // l(0) = 2;
  // l(1) = 2;
  // l(2) = 2;
  // l(5) = -2;
  // l(6) = 1;
  // l(8) = -2;
  // MatrixXXs CL = MatrixXXs::Ones(1,ndof);
  // VectorXs dL = VectorXs::Zero(1);
  // std::vector<std::pair<int, int>> pairs = {{0,ndof-1}};

  // auto solvable = std::make_shared<TestSolvable>(ndof, l, pairs, CL, dL);

  // Debug::log(solvable->val(VectorXs::Ones(ndof)));

  // NewtonSolverSettings settings;
  // settings.projgrad_epsilon=0.99; // TODO this relative stop cond doesnt seem super perfect, not sure what to do. pick something absolute and problem specific?
  // auto solver = std::make_shared<NewtonSolver>(solvable, settings);
  
  // for (size_t i = 0; i < 1000; i++) {
  //   if (solver->isFinished())
  //     break;
  //   auto q = solver->getSolution();
  //   Debug::log("Iteration", i, solvable->val(q), q.transpose());
  //   solver->step();
  // }

  Debug::log("Test end.");
  return 0;
}