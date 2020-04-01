#ifndef _EXPANSION_H_
#define _EXPANSION_H_

#include "EigenDefinitions.h"

class Expansion {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW; // storing fixed size Eigen matrix
  Expansion(scalar sx, scalar sa, scalar sy, Matrix2s II, Vector2s min, Vector2s max, scalar dx) {
    // TODO
  }

  Matrix3s R(Vector2s xi) {return Matrix3s::Identity();}
  Vector3s shell(Vector3s xih) {return xih;}
};
/*

Expansion

  init(sx,sa,sy,II, bounds, dx):
    prep grid
    prep S, dN
    solve phi

  def: R = rotation(xi1,xi2)
  def: xyz = shell(xi1,xi2,h) 
    phi(xi1,xi2)
    + h * n(xi1.xi2)


*/


#endif /* end of include guard: _EXPANSION_H_ */
