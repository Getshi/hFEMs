#ifndef _VEGAFEMSOLVABLE_H_
#define _VEGAFEMSOLVABLE_H_

#include <memory>
#include "Solvable.hpp"
#include "Expansion.hpp"
#include <isotropicHyperelasticFEM.h>
#include <StVKIsotropicMaterial.h> 

class VegaFEMSolvable : public Solvable {
 public:
  VegaFEMSolvable(/* VEGA INITS */);
  ~VegaFEMSolvable();

  std::vector<std::pair<int, int>> getPeriodicPairs();
  std::pair<SparseXXs, VectorXs> getLagrangeConstraints();

  int getNumDof();
  scalar val(const VectorXs& q);      // f(q)
  VectorXs grad(const VectorXs& q);   // df/dq (q)
  SparseXXs hess(const VectorXs& q);  // d2f/dqdq (q)

 private:
  std::unique_ptr<Expansion> m_expansion;
  std::unique_ptr<TetMesh> m_mesh;
  std::unique_ptr<StVKIsotropicMaterial> m_material;
  std::unique_ptr<IsotropicHyperelasticFEM> m_deformableModel;

  VectorXs convert_q_to_u(const VectorXs & q);
};

Vector3s convert_vec(const Vec3d & v);

#endif /* end of include guard: _VEGAFEMSOLVABLE_H_ */
