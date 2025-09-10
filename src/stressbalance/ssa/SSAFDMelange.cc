#include "pism/stressbalance/ssa/SSAFDMelange.hh"

using namespace pism;
using namespace pism::stressbalance;

/*!
Because the FD implementation of the SSA uses Picard iteration, a PETSc KSP
and Mat are used directly.  In particular we set up \f$A\f$
(Mat m_A) and a \f$b\f$ (= Vec m_rhs) and iteratively solve
linear systems
  \f[ A x = b \f]
where \f$x\f$ (= Vec m_velocity_global).  A PETSc SNES object is never created.
copied from SSAFD.cc
 */
SSAFDMelange::SSAFDMelange(std::shared_ptr<const Grid> grid, bool regional_mode)
    : SSAFDBase(grid, regional_mode),
      m_nuH_old(grid, "nuH_old"),
      m_velocity_old(grid, "velocity_old") {

  m_velocity_old.metadata(0)
      .long_name("old SSA velocity field; used for re-trying with a different epsilon")
      .units("m s^-1");

  m_nuH_old.metadata(0)
      .long_name("ice thickness times effective viscosity (before an update)")
      .units("Pa s m");

  // The nuH viewer:
  m_nuh_viewer = nullptr;

  // PETSc objects and settings
  {
    auto dm = m_velocity_global.dm();

    PetscErrorCode ierr;
    ierr = DMSetMatType(*dm, MATAIJ);
    PISM_CHK(ierr, "DMSetMatType");

    ierr = DMCreateMatrix(*dm, m_A.rawptr());
    PISM_CHK(ierr, "DMCreateMatrix");

    ierr = KSPCreate(m_grid->com, m_KSP.rawptr());
    PISM_CHK(ierr, "KSPCreate");

    ierr = KSPSetOptionsPrefix(m_KSP, "ssafd_");
    PISM_CHK(ierr, "KSPSetOptionsPrefix");

    // Use non-zero initial guess (i.e. SSA velocities from the last
    // solve() call).
    ierr = KSPSetInitialGuessNonzero(m_KSP, PETSC_TRUE);
    PISM_CHK(ierr, "KSPSetInitialGuessNonzero");

    // Use the initial residual norm.
    ierr = KSPConvergedDefaultSetUIRNorm(m_KSP);
    PISM_CHK(ierr, "KSPConvergedDefaultSetUIRNorm");
  }
}

//TO DO
void SSAFDMelange::init_impl() {
}

//TO DO
void SSAFDMelange::pc_setup_bjacobi() {
}

//TO DO
void SSAFDMelange::pc_setup_asm() {
}

//TO DO
void SSAFDMelange::solve(const Inputs &inputs) {
}

//TO DO
void SSAFDMelange::picard_iteration(const Inputs &inputs, double nuH_regularization,
                        double nuH_iter_failure_underrelax) {
}

//TO DO
void SSAFDMelange::picard_manager(const Inputs &inputs, double nuH_regularization,
                      double nuH_iter_failure_underrelax) {
}

//TO DO
void SSAFDMelange::picard_strategy_regularization(const Inputs &inputs) {
}

// copied from SSAFD.cc
SSAFDMelange::KSPFailure::KSPFailure(const char* reason)
  : RuntimeError(ErrorLocation(), std::string("SSAFDMelange KSP (linear solver) failed: ") + reason){
  // empty
}

// copied from SSAFD.cc
SSAFDMelange::PicardFailure::PicardFailure(const std::string &message)
  : RuntimeError(ErrorLocation(), "SSAFDMelange Picard iterations failed: " + message) {
  // empty
}
//TO DO
std::array<double, 2> SSAFDMelange::compute_nuH_norm(const array::Staggered &nuH,
                                         array::Staggered &nuH_old) {
  return {};
}

//TO DO
void SSAFDMelange::initialize_iterations(const Inputs &inputs) {
}

 



