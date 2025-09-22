// Copyright (C) 2025 PISM Authors
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

#ifndef _MELANGELOCAL_H_
#define _MELANGELOCAL_H_

#include "pism/melange/Melange.hh"
#include "pism/stressbalance/ssa/SSAFDMelange.hh"
#include "pism/rheology/melange/MelangeRheologyFactory.hh"

namespace pism {

namespace melange {

//! \brief A melange model that uses local rheology and SSA solver.
/*!
  This is the main implementation of the melange model, similar to how
  NullTransport implements the Hydrology interface.

  This model:
  1. Wraps the existing SSAFDMelange solver for melange flow
  2. Uses MelangeRheologyFactory to create appropriate rheology
  3. Handles melange formation from calving
  4. Computes back pressure on ice shelf
  5. Manages melange disintegration

  The physics are implemented in the wrapped components:
  - SSAFDMelange: handles the stress balance and flow
  - MelangeRheology: handles the rheological behavior
*/
class MelangeLocal : public Melange {
public:
  MelangeLocal(std::shared_ptr<const Grid> g);
  virtual ~MelangeLocal() = default;

protected:
  //! Virtual implementations - similar to NullTransport pattern
  virtual void restart_impl(const File &input_file, int record);
  virtual void bootstrap_impl(const File &input_file,
                              const array::Scalar &ice_thickness);
  virtual void init_impl(const array::Scalar &melange_thickness,
                         const array::Scalar &melange_pressure);
  virtual MaxTimestep max_timestep_impl(double t) const;
  virtual void update_impl(double t, double dt, const Inputs& inputs);

  //! Diagnostics implementation
  virtual std::map<std::string, Diagnostic::Ptr> diagnostics_impl() const;

private:
  void initialization_message() const;
  
  //! Wraps your existing components
  std::shared_ptr<stressbalance::SSAFDMelange> m_ssa_solver;
  std::shared_ptr<rheology::melange::MelangeRheology> m_rheology;
  
  //! Configuration
  bool m_use_ssa_solver;
  bool m_use_rheology;
  double m_formation_rate_factor;
  double m_disintegration_rate_factor;
  
  //! Internal state
  array::Scalar m_melange_thickness_old;
  array::Vector m_melange_velocity_old;
  
  //! Helper methods
  void update_melange_formation(const Inputs& inputs, double dt);
  void update_melange_flow(const Inputs& inputs, double dt);
  void update_melange_disintegration(const Inputs& inputs, double dt);
  void update_melange_state(const Inputs& inputs, double dt);
  
  //! SSA solver interface
  void solve_melange_flow(const Inputs& inputs);
  
  //! Rheology interface
  void update_melange_rheology(const Inputs& inputs);
};

} // end of namespace melange
} // end of namespace pism

#endif /* _MELANGELOCAL_H_ */
