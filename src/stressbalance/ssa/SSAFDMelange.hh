// Copyright (C) 2004--2024 Jed Brown, Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef _SSAFDMelange_H_
#define _SSAFDMelange_H_

#include <array>
#include <cassert>
#include <memory>
#include <stdexcept>

#include "pism/stressbalance/ssa/SSAFDBase.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/petscwrappers/Viewer.hh"
#include "pism/util/petscwrappers/KSP.hh"
#include "pism/util/petscwrappers/Mat.hh"
#include "pism/util/petscwrappers/DM.hh"
#include "pism/util/petscwrappers/Vec.hh"
#include "pism/util/array/Staggered.hh"

namespace pism {
namespace stressbalance {


//! Melange SSA solver: the finite difference implementation.
class SSAFDMelange : public SSAFDBase {
public:
  SSAFDMelange(std::shared_ptr<const Grid> g, bool regional_mode);
  virtual ~SSAFDMelange() = default;

protected:

  void init_impl();

  void pc_setup_bjacobi();

  void pc_setup_asm();

  void solve(const Inputs &inputs);

  void picard_iteration(const Inputs &inputs, double nuH_regularization,
                        double nuH_iter_failure_underrelax);

  void picard_manager(const Inputs &inputs, double nuH_regularization,
                      double nuH_iter_failure_underrelax);

  void picard_strategy_regularization(const Inputs &inputs);

  std::array<double, 2> compute_nuH_norm(const array::Staggered &nuH,
                                         array::Staggered &nuH_old);

  void assemble_matrix(const Inputs &inputs, const array::Vector1 &velocity,
                       const array::Staggered1 &nuH, const array::CellType1 &cell_type, Mat A);

  void write_system_petsc(const std::string &namepart);

  void update_nuH_viewers(const array::Staggered &nuH);

  // below are the implementations of the methods from SSAFDBase that are overridden
  // for Melange

  void initialize_iterations(const Inputs &inputs) override;

  void compute_nuH(const array::Scalar1 &ice_thickness, const array::CellType2 &cell_type,
    const pism::Vector2d *const *velocity, const array::Staggered &hardness,
    double nuH_regularization, array::Staggered1 &result) override;

  void compute_nuH_everywhere(const array::Scalar1 &ice_thickness,
                const pism::Vector2d *const *velocity,
                const array::Staggered &hardness, double nuH_regularization,
                array::Staggered &result) override;

  void compute_nuH_cfbc(const array::Scalar1 &ice_thickness,
          const array::CellType2 &cell_type,
          const pism::Vector2d* const* velocity,
          const array::Staggered &hardness, double nuH_regularization,
          array::Staggered &result) override;

  void compute_driving_stress(const array::Scalar &ice_thickness,
                const array::Scalar1 &surface_elevation,
                const array::CellType1 &cell_type,
                const array::Scalar1 *no_model_mask, const EnthalpyConverter &EC,
                array::Vector &result) const override;

  void adjust_driving_stress(const array::Scalar &ice_thickness,
                const array::Scalar1 &surface_elevation,
                const array::CellType1 &cell_type, const array::Scalar1 *no_model_mask,
                array::Vector &driving_stress) const override;

                void assemble_rhs(const Inputs &inputs, const array::CellType1 &cell_type,
                  const array::Vector &driving_stress, double bc_scaling, array::Vector &result) const override;

  void fd_operator(const Geometry &geometry, const array::Scalar *bc_mask, double bc_scaling,
                 const array::Scalar &basal_yield_stress,
                 IceBasalResistancePlasticLaw *basal_sliding_law,
                 const pism::Vector2d *const *velocity, const array::Staggered1 &nuH,
                 const array::CellType1 &cell_type, Mat *A, Vector2d **Ax) const override;

  



  array::Staggered1 m_nuH_old;

  petsc::KSP m_KSP;
  petsc::Mat m_A;

  array::Vector1 m_velocity_old;

  unsigned int m_default_pc_failure_count;
  unsigned int m_default_pc_failure_max_count;
  
  std::shared_ptr<petsc::Viewer> m_nuh_viewer;



  class KSPFailure : public RuntimeError {
  public:
    KSPFailure(const char* reason);
  };

  class PicardFailure : public RuntimeError {
  public:
    PicardFailure(const std::string &message);
  };
};

} // end of namespace stressbalance
} // end of namespace pism

#endif /* _SSAFDMelange_H_ */