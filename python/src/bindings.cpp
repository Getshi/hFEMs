#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;
using namespace pybind11::literals;  // for shorter kwargs

#ifndef VERSION_INFO  // should get defined in setup.py
#define VERSION_INFO "dev"
#endif

#include "../../src/NewtonSolver.hpp"
#include "../../src/TestSolvable.hpp"

// since Solvable is virtual, we need a trampoline python class for
// using/passing obj of that type
// class PySolvable : public Solvable {
// public:
//     using Solvable::Solvable; // Inherit the constructors

//     /* Trampoline (need one for each virtual function) */
//     // PYBIND11_OVERLOAD_PURE(return type, parent class, function name
//     matching cpp and py, arguments)
//     std::vector<std::pair<int, int>> getPeriodicPairs() override {
//         PYBIND11_OVERLOAD_PURE(
//             std::vector<std::pair<int, int>>, Solvable, getPeriodicPairs,
//         );
//     }
//     std::pair<SparseXXs, VectorXs> getLagrangeConstraints() override {
//         PYBIND11_OVERLOAD_PURE(
//             std::pair<SparseXXs, VectorXs>, Solvable, getLagrangeConstraints,
//         );
//     }
//     scalar val(const VectorXs & q) override {
//         PYBIND11_OVERLOAD_PURE(
//             scalar, Solvable, val, q
//         );
//     }
//     VectorXs grad(const VectorXs & q) override {
//         PYBIND11_OVERLOAD_PURE(
//             VectorXs, Solvable, grad, q
//         );
//     }
//     SparseXXs hess(const VectorXs & q) override {
//         PYBIND11_OVERLOAD_PURE(
//             SparseXXs, Solvable, hess, q
//         );
//     }
// };

PYBIND11_MODULE(_solver, m) {
  // general module metadata
  m.doc()               = "......";  // optional module docstring
  m.attr("__version__") = VERSION_INFO;

  // m.def("pyfun", &cppfun, "docstring",
  //       "kwvarname"_a, "kwvarname"_a = default_value);
  // m.def("add",[](int i, int j) {return i+j;});

  // py::class_<Solvable, PySolvable /* <--- trampoline*/>(m, "Solvable")
  //     .def(py::init<>());

  // py::class_<TestSolvable, std::shared_ptr<TestSolvable> /* <- holder type
  // */, Solvable /* <- specify C++ parent type */>(m, "TestSolvable")
  py::class_<TestSolvable, std::shared_ptr<TestSolvable>>(m, "TestSolvable")
      .def(py::init<const std::vector<std::pair<int, int>> &, const MatrixXXs &,
                    const VectorXs &>())
      .def("val", &TestSolvable::val)
      .def("grad", &TestSolvable::grad)
      .def("hess", &TestSolvable::hess);

  // NOTE: using lambda constructor avoids pybind unexposed inherited virtual
  // parameter class passing...

  py::class_<NewtonSolver>(m, "NewtonSolver")
      .def(py::init([](std::shared_ptr<TestSolvable> solvable) {
        return new NewtonSolver(solvable);
      }))
      .def("test", [](NewtonSolver &ns, int n) { return ns.test(n); });
  // .def(py::init<std::shared_ptr<Solvable>>());
}

/*


#include <iostream>
#include <memory>
#include <tuple>

typedef Eigen::Matrix<double, 11, 1> Vec11;
typedef Eigen::Matrix<double, 6, 1> Vec6;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VecX;

using namespace Mathematica;


template <typename T>
void setDebug(const std::string &str, T value) {
  Debug::get<T>(str) = value;
}

template <typename T>
const T &getDebug(const std::string &str) {
  return Debug::get<T>(str);
}

PYBIND11_MODULE(PatternSolver, m) {
  // general module metadata
  m.doc() = "homogenized yarn level cloth";  // optional module docstring
  m.attr("__version__") = VERSION_INFO;

  // m.def("pyfun", &cppfun, "docstring",
  //       "kwvarname"_a, "kwvarname"_a = default_value);

  // general methods
  m.def("setDebug", setDebug<bool>)
      .def("setDebug", setDebug<double>)
      .def("setDebug", setDebug<int>)
      .def("getDebug", getDebug<bool>)
      .def("getDebug", getDebug<double>)
      .def("getDebug", getDebug<int>);

  // // partially exposed YarnSolverSettings
  // py::class_<YarnSolverSettings>(m, "YarnSolverSettings")
  //     .def(py::init<>())
  //     .def_readwrite("stopValue", &YarnSolverSettings::stopValue);
  // // partially exposed HeightSolverSettings
  // py::class_<HeightExpansionSettings>(m, "HeightExpansionSettings")
  //     .def(py::init<const Vec11 &C>())
  //     .def_readwrite("stopValue", &HeightExpansionSettings::stopValue);
  //
  // // simulation settings
  // py::class_<SimulationSettings>(m, "SimulationSettings")
  //     .def(py::init<const std::string &>())
  //     .def_readwrite("protofile", &SimulationSettings::protofile);
  //
  // m.def("test_structarg", [](SimulationSettings &settings) {
  //   return std::make_tuple(settings.protofile, settings.solver.solverType,
  //                          settings.solver.curveatureForceLambda);
  // });

  // microrod simulation interface
  py::class_<Simulation>(m, "Simulation")
      .def(py::init([](const std::string &protofile, const Vec6 &ek, double
thetaref,
                       int gridRes, bool direct_solver, double step_factor,
                       double reg_start, double reg_end, double reg_steps,
                       double yE, double G, double gamma, double density,
                       double kc, double solverscale_x, double stop_tangent,
                       double stop_energy, int stop_energy_history,
                       double thickness_penalty) {
             auto sim = std::make_unique<Simulation>();

             // Simulation *sim = new Simulation();
             SimulationSettings settings;
             settings.protofile = protofile;

             settings.solver = YarnSolverSettings{
                 direct_solver,    step_factor, reg_start,
                 reg_end,          reg_steps,   solverscale_x,
                 stop_tangent,     stop_energy, stop_energy_history,
                 thickness_penalty
                 //  lbfgs,
                 //  LBFGS<scalar>::Params{}
             };
             settings.ek             = ek;
             settings.theta          = thetaref;
             settings.gridResolution = gridRes;
             settings.yE             = yE;
             settings.G              = G;
             settings.gamma          = gamma;
             settings.density        = density;
             settings.kc             = kc;
             sim->reset(settings);
             return sim;
           }),
           "protofile"_a, "ek"_a, "thetaref"_a = 0, "gridRes"_a = 101,
"direct"_a = true,
           "step_factor"_a = 0.2, "reg_start"_a = 50.0, "reg_end"_a = 0.05,
           "reg_steps"_a = 200, "yE"_a = 1e4, "G"_a = 4e3, "gamma"_a = 0.1,
           "density"_a = 1.2e2, "kc"_a = 5e-1, "solverscale_x"_a = 6e3,
           "stop_tangent"_a = 1e-5, "stop_energy"_a = 1e-5,
           "stop_energy_history"_a = 50,
           "thickness_penalty"_a   = 0.0)  //, "lbfgs"_a = false)
      .def("step", &Simulation::update)
      .def("isFinished", &Simulation::isFinished)
      .def("stopType", &Simulation::stopType)
      .def("getResults", &Simulation::getResults)
      .def("serializeState",
           [](Simulation &sim, const std::string &filename) {
             sim.serialize(filename);
           })
      .def("serializeGrid",
           [](Simulation &sim, const std::string &filename) {
             sim.serialize_grid(filename);
           })
      .def("dbgSaveFrameState",
           [](Simulation &sim) {
             sim.getYarnSolver().dbgSaveFrameState(sim.getYarnPatch());
           })
      .def("dbgPerturb", &Simulation::debugPerturb)
      .def("dbgSaveFrameState",
           [](Simulation &sim) {
             sim.getYarnSolver().dbgSaveFrameState(sim.getYarnPatch());
           })
      .def("dbgLoadFrameState",
           [](Simulation &sim) {
             sim.getYarnSolver().dbgLoadFrameState(sim.getYarnPatch());
           })
      .def("setX",
           [](Simulation &sim, const VectorXs &newX, bool do_transport) {
             sim.getYarnPatch().getX() = newX;
             if (do_transport) {
               threadutils::for_each(
                   sim.getYarnPatch().getYarns(),
                   [&](auto yarn) { yarn->updateStrandState(); });
             }
           })
      .def("getX", [](Simulation &sim) { return sim.getYarnPatch().getX(); })
      .def("setY",
           [](Simulation &sim, const VectorXs &newY) {
             sim.getYarnSolver().setState(sim.getYarnPatch(), newY);
           })
      .def("getY",
           [](Simulation &sim) {
             return sim.getYarnSolver().getState(sim.getYarnPatch());
           })
      .def("twistDofMask",
           [](Simulation &sim) {
             auto &yp      = sim.getYarnPatch();
             VectorXs mask = VectorXs::Zero(yp.getNumDof());
             for (auto yarn : yp.getYarns()) {
               for (int i = 0; i < yarn->getNumEdges(); ++i) {
                 mask[yarn->getVertexIndex(i) * 4 + 3] = 1.0;
               }
             }
             return mask;
           })
      .def("actualDofMask",  // ie no tip twists
           [](Simulation &sim) {
             auto &yp = sim.getYarnPatch();
             VectorXs mask(yp.getNumDof());
             mask.setOnes();
             for (auto yarn : yp.getYarns())
               mask[yarn->getVertexIndex(yarn->getNumVertices() - 1) * 4 + 3] =
                   0.0;
             return mask;
           })
      .def("dbgE",
           [](Simulation &sim) {
             return sim.getYarnSolver().getEnergy(sim.getYarnPatch());
           })
      .def("dbggradE",
           [](Simulation &sim) {
             return sim.getYarnSolver().getGradient(sim.getYarnPatch());
           })
      .def("dbggradyE",
           [](Simulation &sim) {
             return sim.getYarnSolver().dbggetgradyE(sim.getYarnPatch());
           })
      .def("dbghessE",
           [](Simulation &sim) {
             return sim.getYarnSolver().getHessian(sim.getYarnPatch());
           })
      .def("dbgupdateCollisions", [](Simulation &sim) {
        return sim.getYarnSolver().updateCollisions(sim.getYarnPatch());
      });
}

*/