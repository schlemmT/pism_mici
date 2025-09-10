//TO DO
#include "pism/stressbalance/ssa/SSAFDMelange.hh"

using namespace pism;
using namespace pism::stressbalance;

//TO DO
void SSAFDMelange::assemble_matrix(const Inputs &inputs, const array::Vector1 &velocity,
                       const array::Staggered1 &nuH, const array::CellType1 &cell_type, Mat A) {
}

//TO DO
void SSAFDMelange::write_system_petsc(const std::string &namepart) {
}

//TO DO
void SSAFDMelange::update_nuH_viewers(const array::Staggered &nuH) {
}

//TO DO
void SSAFDMelange::compute_nuH(const array::Scalar1 &ice_thickness, const array::CellType2 &cell_type,
    const pism::Vector2d *const *velocity, const array::Staggered &hardness,
    double nuH_regularization, array::Staggered1 &result) {
}

//TO DO
void SSAFDMelange::compute_nuH_everywhere(const array::Scalar1 &ice_thickness,
                const pism::Vector2d *const *velocity,
                const array::Staggered &hardness, double nuH_regularization,
                array::Staggered &result) {
}

//TO DO
void SSAFDMelange::compute_nuH_cfbc(const array::Scalar1 &ice_thickness,
          const array::CellType2 &cell_type,
          const pism::Vector2d* const* velocity,
          const array::Staggered &hardness, double nuH_regularization,
          array::Staggered &result) {
}

//TO DO
void SSAFDMelange::compute_driving_stress(const array::Scalar &ice_thickness,
                const array::Scalar1 &surface_elevation,
                const array::CellType1 &cell_type,
                const array::Scalar1 *no_model_mask, const EnthalpyConverter &EC,
                array::Vector &result) const {
}

//TO DO
void SSAFDMelange::adjust_driving_stress(const array::Scalar &ice_thickness,
                const array::Scalar1 &surface_elevation,
                const array::CellType1 &cell_type, const array::Scalar1 *no_model_mask,
                array::Vector &driving_stress) const {
}

//TO DO
void SSAFDMelange::assemble_rhs(const Inputs &inputs, const array::CellType1 &cell_type,
                  const array::Vector &driving_stress, double bc_scaling, array::Vector &result) const {
}

//TO DO
void SSAFDMelange::fd_operator(const Geometry &geometry, const array::Scalar *bc_mask, double bc_scaling,
                 const array::Scalar &basal_yield_stress,
                 IceBasalResistancePlasticLaw *basal_sliding_law,
                 const pism::Vector2d *const *velocity, const array::Staggered1 &nuH,
                 const array::CellType1 &cell_type, Mat *A, Vector2d **Ax) const {
}


