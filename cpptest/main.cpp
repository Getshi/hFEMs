#include <memory>
#include "../src/TestSolvable.hpp"
#include "../src/VegaFEMSolvable.hpp"
#include "../src/NewtonSolver.hpp"
#include "../src/debug/debug_tools.hpp"


int main (int argc, char *argv[]) {
  Debug::log("Test start.");

  auto tmp = std::make_shared<VegaFEMSolvable>();


  Debug::log("Test init.");

  int ndof = 10;
  VectorXs l = VectorXs::Zero(ndof-1);
  l(0) = 2;
  l(1) = 2;
  l(2) = 2;
  l(5) = -2;
  l(6) = 1;
  l(8) = -2;
  MatrixXXs CL = MatrixXXs::Ones(1,ndof);
  VectorXs dL = VectorXs::Zero(1);
  std::vector<std::pair<int, int>> pairs = {{0,ndof-1}};

  auto solvable = std::make_shared<TestSolvable>(ndof, l, pairs, CL, dL);

  Debug::log(solvable->val(VectorXs::Ones(ndof)));

  NewtonSolverSettings settings;
  settings.projgrad_epsilon=0.99; // TODO this relative stop cond doesnt seem super perfect, not sure what to do. pick something absolute and problem specific?
  auto solver = std::make_shared<NewtonSolver>(solvable, settings);
  
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