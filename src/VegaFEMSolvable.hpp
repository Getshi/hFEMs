#ifndef _VEGAFEMSOLVABLE_H_
#define _VEGAFEMSOLVABLE_H_

#include <memory>
#include "Solvable.hpp"
#include "Expansion.hpp"
#include <isotropicHyperelasticFEM.h>
#include <StVKIsotropicMaterial.h> 

class VegaFEMSolvable : public Solvable {
 public:
  VegaFEMSolvable(const std::string & filename, const Vector6s & strain);
  ~VegaFEMSolvable();

  std::vector<std::pair<int, int>> getPeriodicPairs();
  std::pair<SparseXXs, VectorXs> getLagrangeConstraints();

  int getNumDof();
  scalar val(const VectorXs& q);      // f(q)
  VectorXs grad(const VectorXs& q);   // df/dq (q)
  SparseXXs hess(const VectorXs& q);  // d2f/dqdq (q)
  scalar computeMaximumStep(const VectorXs& dq);
  
  VectorXs getWorldPositions(const VectorXs& q); // x = xshell + Rq
  std::tuple<int, int, const MatrixXXs *> getGrid() {
    return m_expansion->getGrid();
  }

 private:
  std::unique_ptr<TetMesh> m_mesh;
  std::unique_ptr<StVKIsotropicMaterial> m_material;
  std::unique_ptr<IsotropicHyperelasticFEM> m_deformableModel;
  std::unique_ptr<SparseMatrix> m_stiffnessMatrix;
  std::unique_ptr<Expansion> m_expansion;
  VectorXs m_xshell;
  SparseXXs m_R; // maybe faster if we store R^T instead
  std::vector<std::pair<int, int>> m_pairs;

  VectorXs convert_q_to_u(const VectorXs & q);
};


bool load_pairs_from_file(const std::string & filename, std::vector<std::pair<int, int>> & pairs);
Vector3s convert_vec(const Vec3d & v);

#endif /* end of include guard: _VEGAFEMSOLVABLE_H_ */
