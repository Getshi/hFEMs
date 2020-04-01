#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/iostream.h>
#include <pybind11/stl.h>
namespace py = pybind11;
using namespace pybind11::literals;  // for shorter kwargs

#ifndef VERSION_INFO  // should get defined in setup.py
#define VERSION_INFO "dev"
#endif

#include "../../src/NewtonSolver.hpp"
#include "../../src/TestSolvable.hpp"


// NOTE .def(name,func,redirect()) to every method that should redirect stdout/stderr
typedef typename py::call_guard<py::scoped_ostream_redirect, py::scoped_estream_redirect> redirect;


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

  // m.def("create_ostream_capsule", []() {
  //   py::print("constructing o");
  //   return py::capsule((void *) new py::scoped_ostream_redirect(), [](void *ptr) {
  //     delete static_cast<py::scoped_ostream_redirect *>(ptr);
  //     py::print("destructing o");
  //   });
  // });
  // m.def("create_estream_capsule", []() {
  //   return py::capsule((void *) new py::scoped_estream_redirect(), [](void *ptr) {
  //     delete static_cast<py::scoped_estream_redirect *>(ptr);
  //   });
  // });

  // py::class_<Solvable, PySolvable /* <--- trampoline*/>(m, "Solvable")
  //     .def(py::init<>());

  // py::class_<TestSolvable, std::shared_ptr<TestSolvable> /* <- holder type
  // */, Solvable /* <- specify C++ parent type */>(m, "TestSolvable")
  py::class_<TestSolvable, std::shared_ptr<TestSolvable>>(m, "TestSolvable")
      .def(py::init<int, const VectorXs &, const std::vector<std::pair<int, int>> &, const MatrixXXs &,
                    const VectorXs &>())
      .def("val", &TestSolvable::val, redirect())
      .def("grad", &TestSolvable::grad, redirect())
      .def("hess", &TestSolvable::hess, redirect());

  py::class_<NewtonSolver>(m, "NewtonSolver")
      .def(py::init([](std::shared_ptr<TestSolvable> solvable) {
        // NOTE: using lambda constructor avoids pybind unexposed inherited virtual parameter class passing...
        return new NewtonSolver(solvable);
      }));
      // NOTE add .def(name,func,redirect()) to every method that should redirect stdout/stderr 
      // .def("test", [](NewtonSolver &ns, int n) { return ns.test(n); });
  // .def(py::init<std::shared_ptr<Solvable>>());
}