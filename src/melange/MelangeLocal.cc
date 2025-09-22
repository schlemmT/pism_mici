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

#include "pism/melange/MelangeLocal.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/File.hh"
#include "pism/util/array/CellType.hh"
#include "pism/geometry/Geometry.hh"

namespace pism {
namespace melange {

//! Constructor - similar to NullTransport constructor
MelangeLocal::MelangeLocal(std::shared_ptr<const Grid> g)
  : Melange(g),
    m_melange_thickness_old(g, "melange_thickness_old"),
    m_melange_velocity_old(g, "melange_velocity_old") {

  // Configuration - similar to NullTransport's m_diffuse_tillwat, etc.
  m_use_ssa_solver = m_config->get_flag("melange.use_ssa_solver", true);
  m_use_rheology = m_config->get_flag("melange.use_rheology", true);
  // Density is now a constant from constants.sea_water.density
  m_formation_rate_factor = m_config->get_number("melange.formation_rate_factor", 1.0);
  m_disintegration_rate_factor = m_config->get_number("melange.disintegration_rate_factor", 1.0);

  // Initialize wrapped components
  if (m_use_ssa_solver) {
    m_ssa_solver = std::make_shared<stressbalance::SSAFDMelange>(g, false);
  }
  
  if (m_use_rheology) {
    m_rheology = rheology::melange::MelangeRheologyFactory::create(g, m_config);
  }

  // Initialize old state
  m_melange_thickness_old.set(0.0);
  m_melange_velocity_old.set(0.0);
}

//! Restart implementation - similar to NullTransport::restart_impl()
void MelangeLocal::restart_impl(const File &input_file, int record) {
  // Read melange state from file
  m_melange_thickness.read(input_file, record);
  m_melange_pressure.read(input_file, record);
  m_melange_velocity.read(input_file, record);
  
  // Initialize old state
  m_melange_thickness_old.copy_from(m_melange_thickness);
  m_melange_velocity_old.copy_from(m_melange_velocity);
  
  // Restart wrapped components if they exist
  if (m_ssa_solver) {
    // SSA solver restart would go here
    // m_ssa_solver->restart(input_file, record);
  }
  
  if (m_rheology) {
    // Rheology restart would go here
    // m_rheology->restart(input_file, record);
  }
}

//! Bootstrap implementation - similar to NullTransport::bootstrap_impl()
void MelangeLocal::bootstrap_impl(const File &input_file, const array::Scalar &ice_thickness) {
  // Initialize melange state from ice thickness
  m_melange_thickness.set(0.0);  // Start with no melange
  
  // Initialize pressure
  m_melange_pressure.set(0.0);
  m_melange_velocity.set(0.0);
  
  // Initialize old state
  m_melange_thickness_old.copy_from(m_melange_thickness);
  m_melange_velocity_old.copy_from(m_melange_velocity);
}

//! Init implementation - similar to NullTransport::init_impl()
void MelangeLocal::init_impl(const array::Scalar &melange_thickness,
                             const array::Scalar &melange_pressure) {
  // Copy initial state
  m_melange_thickness.copy_from(melange_thickness);
  m_melange_pressure.copy_from(melange_pressure);
  
  // Initialize old state
  m_melange_thickness_old.copy_from(m_melange_thickness);
  m_melange_velocity_old.copy_from(m_melange_velocity);
}

//! Max timestep - similar to NullTransport::max_timestep_impl()
MaxTimestep MelangeLocal::max_timestep_impl(double t) const {
  // Melange physics might need smaller timesteps
  double dt_max = m_config->get_number("time_stepping.maximum_time_step");
  
  if (m_use_ssa_solver) {
    // SSA solver might restrict timestep
    double dt_ssa = m_config->get_number("stress_balance.ssa.max_timestep", dt_max);
    dt_max = std::min(dt_max, dt_ssa);
  }
  
  return MaxTimestep(dt_max);
}

//! Update implementation - similar to NullTransport::update_impl()
void MelangeLocal::update_impl(double t, double dt, const Inputs& inputs) {
  // Store old state
  m_melange_thickness_old.copy_from(m_melange_thickness);
  m_melange_velocity_old.copy_from(m_melange_velocity);
  
  // Reset mass change tracking
  m_melange_mass_change.set(0.0);
  
  // Update melange physics in order:
  // 1. Formation from calving
  update_melange_formation(inputs, dt);
  
  // 2. Flow using SSA solver
  if (m_use_ssa_solver) {
    update_melange_flow(inputs, dt);
  }
  
  // 3. Disintegration
  update_melange_disintegration(inputs, dt);
  
  // 4. Update state variables
  update_melange_state(inputs, dt);
  
  // Compute back pressure
  if (inputs.water_column_pressure) {
    compute_back_pressure(m_melange_thickness, 
                         *inputs.ice_thickness, *inputs.water_column_pressure, 
                         m_melange_back_pressure);
  }
}

//! Update melange formation - similar to NullTransport::diffuse_till_water()
void MelangeLocal::update_melange_formation(const Inputs& inputs, double dt) {
  // Simple formation model: calving rate -> melange thickness
  if (inputs.calving_rate) {
    // Only calving contributes to melange formation (not frontal melt)
    array::Scalar formation_rate(*inputs.calving_rate);
    formation_rate.scale(m_formation_rate_factor);
    m_melange_thickness.add(formation_rate, dt);
    
    // Track mass change (positive for formation)
    m_melange_mass_change.add(formation_rate, 1.0);
  } else if (inputs.retreat_rate) {
    // Fallback: use retreat rate if calving rate not available
    array::Scalar formation_rate(*inputs.retreat_rate);
    formation_rate.scale(m_formation_rate_factor);
    m_melange_thickness.add(formation_rate, dt);
    
    // Track mass change (positive for formation)
    m_melange_mass_change.add(formation_rate, 1.0);
  }
}

//! Update melange flow using SSA solver
void MelangeLocal::update_melange_flow(const Inputs& inputs, double dt) {
  if (not m_ssa_solver) {
    return;
  }
  
  // Set up SSA solver inputs
  stressbalance::Inputs ssa_inputs;
  ssa_inputs.geometry = inputs.geometry;
  ssa_inputs.basal_yield_stress = nullptr;  // No basal sliding for melange
  ssa_inputs.basal_sliding_law = nullptr;
  
  // Solve for melange velocity
  solve_melange_flow(ssa_inputs);
  
  // Update melange velocity
  // (This would copy from SSA solver output)
  // m_melange_velocity.copy_from(m_ssa_solver->velocity());
}

//! Update melange disintegration
void MelangeLocal::update_melange_disintegration(const Inputs& inputs, double dt) {
  // Simple disintegration model: frontal melt -> melange loss
  if (inputs.frontal_melt_rate) {
    // Only frontal melt contributes to melange disintegration (not calving)
    array::Scalar disintegration_rate(*inputs.frontal_melt_rate);
    disintegration_rate.scale(m_disintegration_rate_factor);
    m_melange_thickness.add(disintegration_rate, -dt);
    
    // Track mass change (negative for disintegration)
    m_melange_mass_change.add(disintegration_rate, -1.0);
  } else if (inputs.retreat_rate and inputs.calving_rate) {
    // Fallback: compute from retreat rate - calving rate if both available
    array::Scalar disintegration_rate(*inputs.retreat_rate);
    disintegration_rate.add(-1.0, *inputs.calving_rate);
    disintegration_rate.scale(m_disintegration_rate_factor);
    m_melange_thickness.add(disintegration_rate, -dt);
    
    // Track mass change (negative for disintegration)
    m_melange_mass_change.add(disintegration_rate, -1.0);
  }
}

//! Update melange state variables
void MelangeLocal::update_melange_state(const Inputs& inputs, double dt) {
  // Update pressure based on new thickness
  compute_melange_pressure(m_melange_thickness, 
                          *inputs.ice_thickness, m_melange_pressure);
  
  // Density is now a constant, no need to update
  
  // Update rheology if available
  if (m_rheology) {
    update_melange_rheology(inputs);
  }
}

//! Solve melange flow using SSA solver
void MelangeLocal::solve_melange_flow(const stressbalance::Inputs& inputs) {
  if (not m_ssa_solver) {
    return;
  }
  
  // This is where we would call the SSA solver
  // m_ssa_solver->solve(inputs);
  
  // For now, just set velocity to zero
  m_melange_velocity.set(0.0);
}

//! Update melange rheology
void MelangeLocal::update_melange_rheology(const Inputs& inputs) {
  if (not m_rheology) {
    return;
  }
  
  // This is where we would update the rheology
  // m_rheology->update(m_melange_thickness, m_melange_velocity, ...);
}

//! Diagnostics implementation
std::map<std::string, Diagnostic::Ptr> MelangeLocal::diagnostics_impl() const {
  // Similar to Hydrology diagnostics
  std::map<std::string, Diagnostic::Ptr> result;
  
  // Add diagnostics here
  // result["melange_thickness"] = Diagnostic::Ptr(new MelangeThickness(this));
  // result["melange_pressure"] = Diagnostic::Ptr(new MelangePressure(this));
  // etc.
  
  return result;
}

//! Initialization message
void MelangeLocal::initialization_message() const {
  m_log->message(2, 
    "  Using local melange model with SSA solver: %s, Rheology: %s\n",
    m_use_ssa_solver ? "yes" : "no",
    m_use_rheology ? "yes" : "no");
}


} // end of namespace melange
} // end of namespace pism
