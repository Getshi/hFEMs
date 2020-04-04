#include "Expansion.hpp"
#include "debug/debug_tools.hpp"

Vector2s Grid::nodal_xi(int i, int j) const {
  // xi(i,j) = offset + {j*dx, i*dx}
  assert(i >= 0 && j >= 0 && i < nrows && j < ncols && "in bounds");
  Vector2s local;
  local << j * dx, i * dx;
  return offset + local;
}

Vector3s Grid::interpolate_phi(const Vector2s &xi) const {
  // interpolate: grid values
  assert(nodal_phi.rows() == nrows * ncols && nodal_phi.cols() == 3);
  
  // phi(xi) = (phiBL, phiBR, phiTL, phiTR) . (1-b 1-a; 1-b a; b 1-a; b a)

  // find corresponding cell (BL vertex)
  int i = int((xi(1) - offset(1)) / dx);
  int j = int((xi(0) - offset(0)) / dx);
  assert(i >= 0 && j >= 0 && i < nrows - 1 && j < ncols - 1 &&
         "accessing outside grid!");
  i = std::min(std::max(i, 0), nrows - 1);  // clamp anyway
  j = std::min(std::max(j, 0), ncols - 1);

  Vector2s xiBL = nodal_xi(i,j);
  scalar a   = (xi(0) - xiBL(0)) / dx;  // interpolation factors
  scalar b   = (xi(1) - xiBL(1)) / dx;

  auto phiBL        = nodal_phi.row(i * ncols + j);
  auto phiBR        = nodal_phi.row(i * ncols + j + 1);
  auto phiTL        = nodal_phi.row((i + 1) * ncols + j);
  auto phiTR        = nodal_phi.row((i + 1) * ncols + j + 1);
  Vector3s phivalue = (1 - b) * (1 - a) * phiBL + (1 - b) * a * phiBR +
                      b * (1 - a) * phiTL + b * a * phiTR;
  return phivalue;
}

Matrix3s Grid::interpolate_FN(const Vector2s &xi) const {
  // interpolate: FN(xi) = [phi_,1 phi_,2 n]
  assert(nodal_phi.rows() == nrows * ncols && nodal_phi.cols() == 3);

  // find corresponding cell (BL vertex)
  int i = int((xi(1) - offset(1)) / dx);
  int j = int((xi(0) - offset(0)) / dx);
  assert(i >= 0 && j >= 0 && i < nrows - 1 && j < ncols - 1 &&
         "accessing outside grid!");
  i = std::min(std::max(i, 0), nrows - 1);  // clamp anyway
  j = std::min(std::max(j, 0), ncols - 1);

  Vector2s xiBL = nodal_xi(i,j);
  scalar a   = (xi(0) - xiBL(0)) / dx;  // interpolation factors
  scalar b   = (xi(1) - xiBL(1)) / dx;
  scalar invdx = 1.0 / dx;

  auto phiBL        = nodal_phi.row(i * ncols + j);
  auto phiBR        = nodal_phi.row(i * ncols + j + 1);
  auto phiTL        = nodal_phi.row((i + 1) * ncols + j);
  auto phiTR        = nodal_phi.row((i + 1) * ncols + j + 1);
  // compute tangents from derivative of bilinear interpolation
  Vector3s a1 = (1 - b) * (-invdx) * phiBL + (1 - b) * invdx * phiBR +
                b * (-invdx) * phiTL + b * invdx * phiTR;
  Vector3s a2 = (-invdx) * (1 - a) * phiBL + (-invdx) * a * phiBR +
                invdx * (1 - a) * phiTL + invdx * a * phiTR;

  // normal
  Vector3s nn = (a1.cross(a2)).normalized();

  Matrix3s FN;
  FN.col(0) = a1;
  FN.col(1) = a2;
  FN.col(2) = nn;
  return FN;
}

void Grid::recenter_phi(const Vector3s &c) {
  // translate such that phi(xi1=0,xi2=0) = c
  Vector3s dphi = this->interpolate_phi(Vector2s::Zero());
  int ndof = nrows * ncols;
  assert(nodal_phi.rows() == ndof && nodal_phi.cols() == 3);
  for (int ix = 0; ix < ndof; ix++) {
    nodal_phi.row(ix) -= dphi;
  }
}

Expansion::Expansion(const Vector6s &strain, Vector2s min, Vector2s max,
                     scalar dx) {
  // prepare grid
  // get grid size and placement to cover the range [min,max]
  Vector2s extents = max - min;
  int rows         = int(extents(1) / dx + 0.5) + 2;  // round up and add two
  int cols         = int(extents(0) / dx + 0.5) + 2;  // round up and add two
  // align center of (bigger) grid with center of [min,max]
  // i.e. define (xi1=0,xi2=0) at center of [min,max]
  Vector2s c_region = 0.5 * (min + max);
  Vector2s c_grid;
  c_grid << (cols-1) * dx * 0.5, (rows-1) * dx * 0.5;
  m_grid.offset = c_region - c_grid;
  m_grid.nrows  = rows;
  m_grid.ncols  = cols;
  m_grid.dx     = dx;

  // prepare required funcs
  Matrix3x2s S;  // [phi_,1 phi_,2]
  Matrix2s dN;   // [n_,1 n_,2]
  std::tie(S, dN) = convert_strain(strain);
  std::function<Matrix3x2s(const Vector2s &)> eval_poisson_F =
      [&](const Vector2s &xi) {
        Vector2s dn  = dN * xi;
        Matrix3s R   = sincRotationMatrix(dn(0), dn(1));
        Matrix3x2s F = R * S;
        return F;
      };

  // solve on grid:
  //  ∇.∇ phi_i = ∇ . F_(row i)
  //  N . ∇ phi_i = N . F_(row i) on boundary
  solve_poisson(m_grid, eval_poisson_F);

  // recenter solution s.t. phi(0,0) = c
  Vector3s c  = Vector3s::Zero();
  c.head<2>() = c_region;
  m_grid.recenter_phi(c);
}

Matrix3s Expansion::R(Vector2s xi) {
  // irrotational ( grid.FN(xi) )
  Matrix3s F = m_grid.interpolate_FN(xi);
  Eigen::JacobiSVD<Matrix3s> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Matrix3s R =
      svd.matrixU() * svd.matrixV().transpose();  // NOTE: S ~~ V Sigma V^T
  return R;
}

Vector3s Expansion::shell(Vector3s xih) {
  // grid.phi() + h grid.F()col3
  Vector3s phi = m_grid.interpolate_phi(xih.head<2>());
  Matrix3s F   = m_grid.interpolate_FN(xih.head<2>());
  return phi + xih(2) * F.col(2);  // phi(xi1,xi2) + h * n(xi1,xi2)
}

void solve_poisson(Grid &grid,
                   std::function<Matrix3x2s(const Vector2s &)> compute_F) {
  // solve on grid:
  //  ∇.∇ phi_i = ∇ . F_(row i)
  //  N . ∇ phi_i = N . F_(row i) on boundary

  int ndof = grid.nrows * grid.ncols;  // scalar for each grid node

  MatrixXXs rhs(ndof, 3);  // all three right hand sides as columns

  // build LHS and rhs using discrete laplacian and discrete divergence
  // each row can have up to 5 entries (if not used, simply set to 0)
  std::vector<Triplet> triplets(ndof * 5);
  for (int ix = 0; ix < ndof; ix++) {
    int i = int(ix / grid.ncols);
    int j = ix % grid.ncols;

    addLaplaceNeumannRow(i, j, grid, compute_F, triplets, rhs);
  }

  // remove nullspace from rhs by projecting with (I - ee^T), e is vec of ones
  {
    Vector3s sum = Vector3s::Zero();
    for (int ix = 0; ix < ndof; ix++) {
      sum += rhs.row(ix);
    }
    sum *= 1.0 / ndof;
    for (int ix = 0; ix < ndof; ix++) {
      rhs.row(ix) -= sum;
    }
  }

  // assemble and factorize LHS
  SparseXXs LHS(ndof, ndof);
  LHS.setFromTriplets(triplets.begin(), triplets.end());
  Eigen::SimplicialLDLT<SparseXXs> solver(LHS);

  if (solver.info() != Eigen::Success)
    Debug::errorf("Linear solver failed to factorize!\n");



  // solve (column-wise)
  grid.nodal_phi.resize(ndof, 3);
  for (int i = 0; i < 3; ++i) {
    grid.nodal_phi.col(i) = solver.solve(rhs.col(i));

    if (solver.info() != Eigen::Success) {
      Debug::errorf("Linear solver failed to solve!\n");
    }
  }
}

// i think this creates L and r for L phi = div r
// but scaled, i.e. actually Lhat = L*-dx^2, div rhat = div r * -dx^2
// and using finite differencing for div r
void addLaplaceNeumannRow(int i, int j, const Grid &grid,
                          std::function<Matrix3x2s(const Vector2s &)> compute_F,
                          std::vector<Triplet> &triplets, MatrixXXs &rhs) {
  bool bdryL = j == 0;
  bool bdryR = j == grid.ncols - 1;
  bool bdryB = i == 0;
  bool bdryT = i == grid.nrows - 1;
  int ix     = i * grid.ncols + j;

  Vector2s xi = grid.nodal_xi(i, j);

  int n_bounds = int(bdryL) + int(bdryR) + int(bdryT) + int(bdryB);
  DECLARE_UNUSED(n_bounds);
  assert(!(bdryL && bdryR) && !(bdryB && bdryT) && n_bounds < 3);

  Triplet empty(ix, ix, 0.0);

  rhs.row(ix).setZero();
  int diag = 0.0;

  Vector2s dx1, dx2;
  dx1 << grid.dx, 0;
  dx2 << 0, grid.dx;

  if (!bdryL) {
    rhs.row(ix) += compute_F(xi - 0.5 * dx1).col(0) * grid.dx;
    diag++;
    triplets[ix * 5 + 4] = Triplet(ix, i * grid.ncols + j - 1, -1);
  }
  if (!bdryR) {
    rhs.row(ix) -= compute_F(xi + 0.5 * dx1).col(0) * grid.dx;
    diag++;
    triplets[ix * 5 + 3] = Triplet(ix, i * grid.ncols + j + 1, -1);
  }
  if (!bdryB) {
    rhs.row(ix) += compute_F(xi - 0.5 * dx2).col(1) * grid.dx;
    diag++;
    triplets[ix * 5 + 2] = Triplet(ix, (i - 1) * grid.ncols + j, -1);
  }
  if (!bdryT) {
    rhs.row(ix) -= compute_F(xi + 0.5 * dx2).col(1) * grid.dx;
    diag++;
    triplets[ix * 5 + 1] = Triplet(ix, (i + 1) * grid.ncols + j, -1);
  }

  triplets[ix * 5 + 0] = Triplet(ix, ix, diag);
}

std::pair<Matrix3x2s, Matrix2s> convert_strain(const Vector6s &strain) {
  // compute irrotational in-plane deformation using polar/svd
  Matrix2s S2x2;
  scalar sx1 = strain(0)+1, sa = strain(1), sy1 = strain(2)+1;
  S2x2 << sx1, sa * sy1, 0, std::sqrt(1 - sa * sa) * sy1;
  Eigen::JacobiSVD<Matrix2s> svd(S2x2,
                                 Eigen::ComputeFullU | Eigen::ComputeFullV);
  S2x2 = svd.matrixV() * svd.singularValues().asDiagonal() *
         svd.matrixV().transpose();
  Matrix3x2s S = Matrix3x2s::Zero();
  S.block<2, 2>(0, 0) = S2x2;

  // second fundamental form
  Matrix2s II;
  II << strain(3), strain(4), strain(4), strain(5);

  // compute normal derivatives
  Matrix2s dN = -S2x2.transpose().inverse() * II;

  return std::make_pair(S, dN);
}

scalar sinc(scalar x) {
  if (std::abs(x) < 1e-10) {
    // 4th order taylor approximation
    scalar x2                = x * x;
    static const scalar i6   = 1.0 / 6.0;
    static const scalar i120 = 1.0 / 120.0;
    return 1 - x2 * i6 + x2 * x2 * i120;
  } else {
    return std::sin(x) / x;
  }
}

Matrix3s sincRotationMatrix(const scalar &a, const scalar &b) {
  // this is a numerically safe (because using sinc) formula of constructing
  // a rodrigues rotation matrix with u={-b,a,0}
  // or equivalently e^([0,0,a; 0,0,b; -a,-b,0])
  scalar theta = std::sqrt(a * a + b * b);
  scalar c0    = sinc(theta);
  scalar c1    = sinc(theta / 2);
  c1           = 0.5 * c1 * c1;

  Matrix3s R;
  R << 1 - a * a * c1, -a * b * c1, a * c0, -a * b * c1, 1 - b * b * c1, b * c0,
      -a * c0, -b * c0, std::cos(theta);
  return R;
}