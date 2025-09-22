// Copyright (C) 2025 PISM Authors
//
// Local melange rheology (skeleton).

#ifndef PISM_RHEOLOGY_MELANGE_LOCAL_HH
#define PISM_RHEOLOGY_MELANGE_LOCAL_HH

#include <memory>

#include "pism/util/Vector2d.hh"
#include "pism/util/Grid.hh"

namespace pism {

namespace array {
class Scalar;
class Array3D;
}
class Config;

// To DO: maybe move this to a common header file? It is copied from FlowLaw.hh
/*!
    * This uses the definition of squared second invariant from Hutter and several others, namely the
    * output is @f$ D^2 = \frac 1 2 D_{ij} D_{ij} @f$ where incompressibility is used to compute
    * @f$ D_{zz}. @f$
    *
    * This is the approximation of the full second invariant corresponding to the shallow shelf
    * approximation. In particular, we assume that @f$ u @f$ and @f$ v @f$ are depth-independent (@f$
    * u_z = v_z = 0 @f$) and neglect horizontal derivatives of the vertical velocity (@f$ w_x = w_y = 0
    * @f$).
    */
static inline double secondInvariant_2D(const Vector2d &U_x, const Vector2d &U_y) {
    const double
    u_x = U_x.u,
    u_y = U_y.u,
    v_x = U_x.v,
    v_y = U_y.v,
    w_z = -(u_x + v_y);         // w_z is computed assuming incompressibility of ice
    return 0.5 * (u_x * u_x + v_y * v_y + w_z * w_z + 0.5*(u_y + v_x)*(u_y + v_x));
}


namespace melange {

//TO DO: make this implement the common MelangeRheology interface once defined
class MelangeLocal {
public:
  explicit MelangeLocal(std::shared_ptr<const Grid> grid);
  ~MelangeLocal() = default;

  //TO DO: add the actual public API for the melange model
};

} // namespace melange
} // namespace pism

#endif // PISM_RHEOLOGY_MELANGE_LOCAL_HH


