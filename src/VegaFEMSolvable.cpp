#include "VegaFEMSolvable.hpp"
#include <volumetricMeshLoader.h>
#include <fstream>
#include "debug/debug_tools.hpp"

VegaFEMSolvable::VegaFEMSolvable(const std::string & filename, const Vector6s & strain) {

  bool pairs_loaded = load_pairs_from_file(filename, m_pairs);
  DECLARE_UNUSED(pairs_loaded);
  assert(pairs_loaded);

  m_mesh = std::unique_ptr<TetMesh>(
      static_cast<TetMesh*>(VolumetricMeshLoader::load(filename.c_str())));
  if (!m_mesh)
    printf("Error: failed to load mesh.\n");
  else
    printf("Success. Number of vertices: %d . Number of elements: %d .\n",
           m_mesh->getNumVertices(), m_mesh->getNumElements());

  int enableCompressionResistance = 1;
  double compressionResistance    = 0.1;

  m_material = std::make_unique<StVKIsotropicMaterial>(
      m_mesh.get(), enableCompressionResistance, compressionResistance);
  m_deformableModel = std::make_unique<IsotropicHyperelasticFEM>(
      m_mesh.get(), m_material.get());

  // precompute (constant) stiffnessmatrix topology
  SparseMatrix* raw_stiffnessMatrix;
  m_deformableModel->GetStiffnessMatrixTopology(&raw_stiffnessMatrix);
  m_stiffnessMatrix = std::unique_ptr<SparseMatrix>(raw_stiffnessMatrix);

  // Expansion
  auto bb       = m_mesh->getBoundingBox();
  Vector2s bmin = convert_vec(bb.bmin()).head<2>();
  Vector2s bmax = convert_vec(bb.bmax()).head<2>();
  scalar dx     = std::min(bmax(0) - bmin(0), bmax(1) - bmin(1)) * 1e-1;
  m_expansion   = std::make_unique<Expansion>(strain, bmin, bmax, dx);

  // precompute xshell and R at each vertex
  int nverts = m_mesh->getNumVertices();
  m_xshell.resize(nverts * 3);
  m_R.resize(nverts * 3, nverts * 3);
  std::vector<Triplet> triplets;
  triplets.reserve(nverts * 9);
  for (int i = 0; i < nverts; i++) {
    Vector3s X                 = convert_vec(m_mesh->getVertex(i));
    m_xshell.segment<3>(3 * i) = m_expansion->shell(X);
    Matrix3s R                 = m_expansion->R(X.head<2>());
    for (int ki = 0; ki < 3; ki++) {
      for (int kj = 0; kj < 3; kj++) {
        if (std::abs(R(ki, kj)) > EPSILON)
          triplets.emplace_back(3 * i + ki, 3 * i + kj, R(ki, kj));
      }
    }
  }
  m_R.setFromTriplets(triplets.begin(), triplets.end());
}
VegaFEMSolvable::~VegaFEMSolvable() {}

std::vector<std::pair<int, int>> VegaFEMSolvable::getPeriodicPairs() {
  // // find matching vertices on opposite boundaries
  // // assuming mesh has vertices on an axis-aligned hexahedral grid, perfectly
  // contained in its aabb, no holes on the boundaries

  // auto bb = m_mesh->getBoundingBox();
  // Vector2s bmin = convert_vec(bb.bmin()).head<2>();
  // Vector2s bmax = convert_vec(bb.bmax()).head<2>();
  // int nverts = m_mesh->getNumVertices();
  // Vector2s dx = bmax - bmin;

  // auto near = [](scalar a, scalar b){ return std::abs(a-b) < EPSILON; };

  // Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic> locs =
  // Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic>::Zero(nverts,4); //
  // vertex i is near boundary j

  // // assign boundaries
  // int UP=0,LEFT=1,DOWN=2,RIGHT=3; // wrt view of coords (X[0],X[1])
  // for (int i = 0; i < nverts; i++) {
  //   const auto & X = m_mesh->getVertex(i);
  //   if (near(X[0],bmin[0])) {
  //     locs(i, LEFT) = true;
  //     assert(!(near(X[0],bmax[0]))); // cant be on both sides at once
  //   } else if (near(X[0],bmax[0])) {
  //     locs(i, RIGHT) = true;
  //   }

  //   if (near(X[1],bmin[1])) {
  //     locs(i, DOWN) = true;
  //     assert(!(near(X[1],bmax[1]))); // cant be on both sides at once
  //   } else if (near(X[1],bmax[1])) {
  //     locs(i, UP) = true;
  //   }
  // }
  // find all left (inc top/bottom)
  // find all right (inc top/bottom)
  // find pairs between the two

  // find all top (exc left/right)
  // find all bottom (exc left/right)
  // find pairs between the two

  // TODO check how newtonsolver deals with empty periodicity

  // because m_pairs stores vertex periodicity, need to convert to DOF periodicity
  std::vector<std::pair<int, int>> dof_pairs;
  dof_pairs.reserve(m_pairs.size()*3);
  for (int i = 0; i < (int)m_pairs.size(); i++) {
    dof_pairs.emplace_back(3 * m_pairs[i].first + 0, 3 * m_pairs[i].second + 0);
    dof_pairs.emplace_back(3 * m_pairs[i].first + 1, 3 * m_pairs[i].second + 1);
    dof_pairs.emplace_back(3 * m_pairs[i].first + 2, 3 * m_pairs[i].second + 2);
  }

  return dof_pairs;
}

std::pair<SparseXXs, VectorXs> VegaFEMSolvable::getLagrangeConstraints() {
  assert(m_mesh);
  int nverts = m_mesh->getNumVertices();
  int ndof   = nverts * 3;
  std::vector<Triplet> triplets;

  // TODO nice correct quadrature
  // remember that constraint is of form CL q = dL, i.e. using q = R^T utilde

  // actually not sure.. energy force and hess are also using linear tets, does
  // that imply R evaluated on vertices only?
  // kinda i think. e.g. hessian is -R^T K R, where K depends linearly on u, and
  // R we could define linear interpolated within tet or actually nonlinear, but
  // then evaluation of hessian would require quadrature
  // i.e. H = sum_tets int_tetX -R^T(X) K(u) R(X). complicated, not worth.

  // or maybe try correct volume integrals but assuming R etc per vertex and
  // linear interpol
  // so as to not need gaussian quadrature

  // updated understanding: for per vertex q_i, R_i, x_i etc. (i.e. R or xshell are also discretized with same linear tet. basis functions)
  // then:
  //   int_mesh q dV = 0
  //   sum_tet int_tet sum_i q_i bary_i(X) dV  = 0
  //   sum_tet sum_i q_i int_tet bary_i(X) dV  = 0
  //   sum_tet sum_i q_i |V|/4  = 0
  //   CL q = dL
  //   -> dL = 0, CL = for each vertex sum |V|/4 of each tet
  //   can normalize CL also by some constant, e.g. total volume

  scalar Vtotal = 0;
  for (int el = 0; el < m_mesh->getNumElements(); el++) {
    Vtotal += m_mesh->getElementVolume(el);
  }
  scalar fac = 1/Vtotal;

  SparseXXs CL(3, ndof);
  CL.reserve(VectorXi::Constant(ndof, 1)); // each column will have 1 nonzero
  for (int el = 0; el < m_mesh->getNumElements(); el++) {
    scalar V = m_mesh->getElementVolume(el);
    scalar localfac = V * 0.25 * fac; // |V|/4 scaled by 1/Vtotal
    for (int i = 0; i < 4; i++) {
      int vix = m_mesh->getVertexIndex(el, i);
      CL.coeffRef(0, 3 * vix + 0) += localfac;
      CL.coeffRef(1, 3 * vix + 1) += localfac;
      CL.coeffRef(2, 3 * vix + 2) += localfac;
    }
  }
  CL.makeCompressed();

  VectorXs dL = VectorXs::Zero(CL.rows()) * fac;
  return std::make_pair(CL, dL);
}

int VegaFEMSolvable::getNumDof() {
  assert(m_mesh);
  return m_mesh->getNumVertices() * 3;
}

VectorXs VegaFEMSolvable::getWorldPositions(const VectorXs& q) {
  // x = xshell + Rq
  return m_xshell + m_R * q;
}

VectorXs VegaFEMSolvable::convert_q_to_u(const VectorXs& q) {
  // X + u = xshell + Rq
  VectorXs u = m_xshell + m_R * q;  // x
  int nverts = m_mesh->getNumVertices();
  for (int i = 0; i < nverts; i++) {
    u.segment<3>(3 * i) -= convert_vec(m_mesh->getVertex(i));  // x - X
  }
  return u;
}

scalar VegaFEMSolvable::val(const VectorXs& q) {
  // u = xshell + Rq - X
  VectorXs u = convert_q_to_u(q);
  // VectorXs u = q; //DEBUG
  scalar E   = m_deformableModel->ComputeEnergy(u.data());
  return E;
}

VectorXs VegaFEMSolvable::grad(const VectorXs& q) {
  // u = xshell + Rq - X
  VectorXs u = convert_q_to_u(q);
  // VectorXs u = q; //DEBUG
  VectorXs G(u.rows());
  m_deformableModel->ComputeForces(u.data(), G.data());

  // dE/dq = dE/du du/dq;
  // dE/du = -f
  // u = xshell + Rq - X; -> du/dq = R
  // -> dE/dq = -R^T f
  return (m_R.transpose() * G);
  // return (f);//DEBUG
  // NOTE for some reason it seems like in vega dE/du = -G
}

SparseXXs VegaFEMSolvable::hess(const VectorXs& q) {
  // u = xshell + Rq - X
  VectorXs u = convert_q_to_u(q);
  // VectorXs u = q; //DEBUG

  // K = -d2E/dudu /// well actually same as with forces because vega gives "G" and here "dG/du" somehow means minus sign is wrong
  m_deformableModel->GetTangentStiffnessMatrix(u.data(),
                                               m_stiffnessMatrix.get());

  std::vector<Triplet> triplets;
  triplets.reserve(m_stiffnessMatrix->GetNumEntries());
  int nrows = m_stiffnessMatrix->GetNumRows();
  // int ncols = nrows;  // assume symmetry, vega's GetNumColumns() can be wrong
                      // when columns are 0
  for (int i = 0; i < nrows; i++)
    for (int sj = 0; sj < m_stiffnessMatrix->GetRowLength(i); sj++) {
      int val = m_stiffnessMatrix->GetEntries()[i][sj];
      int col = m_stiffnessMatrix->GetColumnIndices()[i][sj];
      triplets.emplace_back(i, col, val);
    }
  SparseXXs dGdu(u.rows(), u.rows());
  dGdu.setFromTriplets(triplets.begin(), triplets.end());
  // d2E/dqdq = R^T (dGdu) R
  SparseXXs hessE = m_R.transpose() * dGdu * m_R;
  // SparseXXs hessE = dGdu; // DEBUG
  return hessE;
}

scalar VegaFEMSolvable::computeMaximumStep(const VectorXs& dq) {
  scalar maxdq2 = 0;
  int nverts = m_mesh->getNumVertices();
  for (int i = 0; i < nverts; i++) {
    maxdq2 = std::max(maxdq2, dq.segment<3>(3*i).squaredNorm());
  }
  return std::sqrt(maxdq2);
}


Vector3s convert_vec(const Vec3d& v) {
  Vector3s v2;
  v2[0] = v[0];
  v2[1] = v[1];
  v2[2] = v[2];
  return v2;
}

bool load_pairs_from_file(const std::string& filename,
                          std::vector<std::pair<int, int>>& pairs) {
  // *PERIODICPAIRS NPAIRS
  // i0 j0
  // i1 j1
  // ...

  std::ifstream in(filename);
  if (!in)
    return false;

  pairs.clear();

  bool accumulating = false;
  std::string str_line;
  while (std::getline(in, str_line)) {
    if (!accumulating) {
      if (str_line.rfind("*PERIODICPAIRS", 0) == 0) {
        accumulating = true;
        std::stringstream ss(str_line);
        std::string tmp;
        int n;
        ss >> tmp >> n;
        pairs.reserve(n);
      }
    } else {  // accumulating
      if (str_line.empty())
        break;
      else {
        std::stringstream ss(str_line);
        int i, j;
        ss >> i >> j;
        pairs.emplace_back(i, j);
      }
    }
  }

  return true;
}