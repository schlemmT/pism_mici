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

#include "pism/melange/MelangeNull.hh"
#include "pism/util/MaxTimestep.hh"

namespace pism {
namespace melange {

//! Constructor - similar to NullTransport constructor
MelangeNull::MelangeNull(std::shared_ptr<const Grid> g)
  : Melange(g) {
  // Initialize everything to zero
  m_melange_thickness.set(0.0);
  m_melange_pressure.set(0.0);
  m_melange_velocity.set(0.0);
  m_melange_mass_change.set(0.0);
  m_melange_back_pressure.set(0.0);
}

//! Restart implementation - do nothing
void MelangeNull::restart_impl(const File &input_file, int record) {
  // Null model - do nothing
  (void)input_file;
  (void)record;
}

//! Bootstrap implementation - do nothing
void MelangeNull::bootstrap_impl(const File &input_file, const array::Scalar &ice_thickness) {
  // Null model - do nothing
  (void)input_file;
  (void)ice_thickness;
}

//! Init implementation - do nothing
void MelangeNull::init_impl(const array::Scalar &melange_thickness,
                            const array::Scalar &melange_pressure) {
  // Null model - do nothing
  (void)melange_thickness;
  (void)melange_pressure;
}

//! Max timestep - no restrictions
MaxTimestep MelangeNull::max_timestep_impl(double t) const {
  // Null model - no timestep restrictions
  (void)t;
  return MaxTimestep();
}

//! Update implementation - do nothing
void MelangeNull::update_impl(double t, double dt, const Inputs& inputs) {
  // Null model - do nothing
  (void)t;
  (void)dt;
  (void)inputs;
}

//! Initialization message
void MelangeNull::initialization_message() const {
  m_log->message(2, "  Using null melange model (no melange physics)\n");
}

} // end of namespace melange
} // end of namespace pism
