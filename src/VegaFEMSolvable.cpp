#include "VegaFEMSolvable.hpp"
#include <volumetricMeshLoader.h>
#include "debug/debug_tools.hpp"

/*

// TODO
// get hextet meshing, write model to file
// initialize deformable model here from file
// wrap energy functions
//  convert q=g to u
    // X+u = xshell + R g
    // u = xshell - X + R g; precompute xshell and R at each vert into big
matrices
    // call stuff
    // convert sparse to eigensparse

// new function positions from g (for rendering)

// wrap periodic pairs & lagrange
// curved expansion

*/

VegaFEMSolvable::VegaFEMSolvable(/* VEGA INITS */) {
  Debug::log("0000");
  Debug::log("1111");


  std::string filename = "/home/gsperl/Projects/hFEMs/meshes/testmesh.veg";

  m_mesh = std::unique_ptr<TetMesh>(
      static_cast<TetMesh*>(VolumetricMeshLoader::load(filename.c_str())));
  if (!m_mesh)
    printf("Error: failed to load mesh.\n");
  else
    printf("Success. Number of vertices: %d . Number of elements: %d .\n",
           m_mesh->getNumVertices(),
           m_mesh->getNumElements());

  int enableCompressionResistance = 1;
  double compressionResistance    = 0.1;

  m_material = std::make_unique<StVKIsotropicMaterial>(
      m_mesh.get(), enableCompressionResistance, compressionResistance);
  m_deformableModel = std::make_unique<IsotropicHyperelasticFEM>(
      m_mesh.get(), m_material.get());


  // Expansion
  scalar sx, sa, sy;
  sx = 0.0;
  sa = 0.0;
  sy = 0.0;
  Matrix2s II = Matrix2s::Identity();
  auto bb = m_mesh->getBoundingBox();
  Vector2s bmin = convert_vec(bb.bmin).head<2>();
  Vector2s bmax = convert_vec(bb.bmax).head<2>();
  scalar dx = std::min(bmax(0)-bmin(0),bmax(1)-bmin(1)) * 1e-2;
  m_expansion = std::make_unique<Expansion>(sx, sa, sy, II, bmin, bmax, dx);

  // TODO precompute xshell per vertex, and R
}
VegaFEMSolvable::~VegaFEMSolvable() {}

std::vector<std::pair<int, int>> VegaFEMSolvable::getPeriodicPairs() {
  // TODO assume we have/get coords of left/right/top/bottom as aabb
  // find all left (inc top/bottom)
  // find all right (inc top/bottom)
  // find pairs between the two

  // find all top (exc left/right)
  // find all bottom (exc left/right)
  // find pairs between the two
  std::vector<std::pair<int, int>> pairs = {}; // TODO check how newtonsolver deals with empty periodicity
  return pairs;
}

std::pair<SparseXXs, VectorXs> VegaFEMSolvable::getLagrangeConstraints() {
  assert(m_mesh);
  int nverts = m_mesh->getNumVertices();
  int ndof   = nverts * 3;
  std::vector<Triplet> triplets;

  // TODO nice correct quadrature
  // remember that constraint is of form CL q = dL, i.e. using q = R^T utilde

  // stupid simple:
  scalar Ninv = 1 / nverts;
  SparseXXs CL(3, ndof);
  CL.reserve(VectorXi::Constant(ndof, 1));
  for (int i = 0; i < nverts; i++) {
    CL.insert(0, 3 * i + 0) = Ninv;
    CL.insert(1, 3 * i + 1) = Ninv;
    CL.insert(2, 3 * i + 2) = Ninv;
  }
  CL.makeCompressed();

  VectorXs dL = VectorXs::Zero(CL.rows());
  return std::make_pair(CL, dL);
}

int VegaFEMSolvable::getNumDof() {
  assert(m_mesh);
  return m_mesh->getNumVertices() * 3;
}

VectorXs VegaFEMSolvable::convert_q_to_u(const VectorXs & q) {
  // u = xshell + Rq - X
  VectorXs u(q.rows());

  // int nverts = m_mesh->getNumVertices();
  // for (int i = 0; i < nverts; i++) {
  //   //
  // }

  // TODO if we prcomp xshell and Rglobal into flattened dof form then
  // VectorXs u = xshell + R * q;
  // u_i -= X_i

  // m_mesh->getVertex(vix) // this is undeformed position X

}


scalar VegaFEMSolvable::val(const VectorXs& q) {
  // u = xshell + Rq - X
  VectorXs u = convert_q_to_u(q);
  scalar E = m_deformableModel->ComputeEnergy(u.data());
  return E;
}
VectorXs VegaFEMSolvable::grad(const VectorXs& q) {
  // u = xshell + Rq - X
  VectorXs u = convert_q_to_u(q);
  VectorXs f(u.rows());
  m_deformableModel->ComputeForces(u.data(), f.data());
  return -f;
}

SparseXXs VegaFEMSolvable::hess(const VectorXs& q) {
  // u = xshell + Rq - X
  VectorXs u = convert_q_to_u(q);
  // TODO
  // m_deformableModel->GetTangentStiffnessMatrix
  // call hessian (opt precomputed topology)
  // void GetStiffnessMatrixTopology(SparseMatrix ** tangentStiffnessMatrix);
  // void GetTangentStiffnessMatrix(const double * u, SparseMatrix *
  // tangentStiffnessMatrix); // use both methods for speed and check how vega
  // solver uses it ?
  // multiply -1
  // convert to eigen
  SparseXXs hessE;
  return hessE;
}

Vector3s convert_vec(const Vec3d & v) {
  Vector3s v2;
  v2[0] = v[0];
  v2[1] = v[1];
  v2[2] = v[2];
  return v2;
}