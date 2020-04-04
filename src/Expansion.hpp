#ifndef _EXPANSION_H_
#define _EXPANSION_H_

#include "EigenDefinitions.hpp"

struct Grid {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;  // storing fixed size Eigen matrix

  int nrows, ncols;
  scalar dx;
  Vector2s offset;
  MatrixXXs nodal_phi;  // (nrows*ncols) x 3 Matrix of surface values at nodes

  Vector2s nodal_xi(int i, int j) const;  // xi(i,j) = offset + {j*dx, i*dx}
  Vector3s interpolate_phi(
      const Vector2s &xi) const;  // bilinearly interpolate phi
  Matrix3s interpolate_FN(
      const Vector2s &xi) const;  // FN(xi) = [phi_,1 phi_,2 n]
  void recenter_phi(
      const Vector3s &c);  // translate such that phi(xi1=0,xi2=0) = c
};

class Expansion {  // xshell(xi1,xi2,h) = phi(xi1,xi2) + h n(xi1,xi2)
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;  // storing fixed size Eigen matrix
  Expansion(const Vector6s &strain, Vector2s min, Vector2s max, scalar dx);

  Matrix3s R(Vector2s xi);
  Vector3s shell(Vector3s xih);

  std::tuple<int, int, const MatrixXXs *> getGrid() const {
    return std::make_tuple(m_grid.nrows, m_grid.ncols, &m_grid.nodal_phi);
  }

 private:
  Grid m_grid;
};

std::pair<Matrix3x2s, Matrix2s> convert_strain(const Vector6s &strain);
scalar sinc(
    scalar x);  // sin(x)/x, with taylor around 0 for numerical stability
Matrix3s sincRotationMatrix(const scalar &a,
                            const scalar &b);  // numerically stable rotation
                                               // e^([0,0,a; 0,0,b; -a,-b,0])

//  ∇.∇ phi_i = ∇ . F_(row i); Neumann: N . ∇ phi_i = N . F_(row i)
void solve_poisson(Grid &grid,
                   std::function<Matrix3x2s(const Vector2s &)> compute_F);
void addLaplaceNeumannRow(int i, int j, const Grid &grid,
                          std::function<Matrix3x2s(const Vector2s &)> compute_F,
                          std::vector<Triplet> &triplets, MatrixXXs &rhs);

#endif /* end of include guard: _EXPANSION_H_ */
