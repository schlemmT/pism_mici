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

#include "pism/melange/Melange.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/File.hh"
#include "pism/util/array/CellType.hh"
#include "pism/geometry/Geometry.hh"

namespace pism {
namespace melange {

//! Inputs constructor
Inputs::Inputs() {
  no_model_mask = nullptr;  
  geometry = nullptr;
  ice_velocity = nullptr;
  basal_melt_rate = nullptr;
  calving_rate = nullptr;
  frontal_melt_rate = nullptr;
  water_column_pressure = nullptr;
}

//! Melange constructor
Melange::Melange(std::shared_ptr<const Grid> g)
  : Component(g),
    m_melange_thickness(g, "melange_thickness"),
    m_melange_velocity(g, "melange_velocity"),
    m_hydrostatic_pressure(g, "hydrostatic_pressure"),
    m_granular_pressure(g, "granular_pressure"),
    m_melange_mass_change(g, "melange_mass_change"),
    m_melange_back_pressure(g, "melange_back_pressure") {

  // Set up metadata - similar to Hydrology
  m_melange_thickness.metadata()
    .long_name("thickness of melange layer")
    .units("m")
    .standard_name("land_ice_melange_thickness");

  m_hydrostatic_pressure.metadata()
    .long_name("hydrostatic pressure in melange layer")
    .units("Pa")
    .standard_name("melange_hydrostatic_pressure");

  m_granular_pressure.metadata()
    .long_name("granular pressure in melange layer")
    .units("Pa")
    .standard_name("melange_granular_pressure");

  m_melange_velocity.metadata(0)
    .long_name("melange velocity in x direction")
    .units("m s-1")
    .standard_name("land_ice_melange_velocity_x");

  m_melange_velocity.metadata(1)
    .long_name("melange velocity in y direction")
    .units("m s-1")
    .standard_name("land_ice_melange_velocity_y");

  m_melange_mass_change.metadata()
    .long_name("rate of change of melange mass")
    .units("kg s-1")
    .standard_name("melange_mass_change_rate");

  m_melange_back_pressure.metadata()
    .long_name("back pressure from melange on ice shelf")
    .units("Pa")
    .standard_name("melange_back_pressure");


  // Initialize to zero
  m_melange_thickness.set(0.0);
  m_melange_velocity.set(0.0);
  m_hydrostatic_pressure.set(0.0);
  m_granular_pressure.set(0.0);
  m_melange_mass_change.set(0.0);
  m_melange_back_pressure.set(0.0);
}

//! Public interface methods 
void Melange::restart(const File &input_file, int record) {
  restart_impl(input_file, record);
}

void Melange::bootstrap(const File &input_file) {
  bootstrap_impl(input_file);
}

void Melange::init(const array::Scalar &melange_thickness,
                   const array::Vector &melange_velocity) {
  init_impl(melange_thickness, melange_velocity);
}

void Melange::update(double t, double dt, const Inputs& inputs) {
  update_impl(t, dt, inputs);
}

//! Accessor methods - similar to Hydrology's till_water_thickness(), etc.
const array::Scalar& Melange::melange_thickness() const {
  return m_melange_thickness;
}

const array::Vector& Melange::melange_velocity() const {
  return m_melange_velocity;
}

const array::Scalar& Melange::hydrostatic_pressure() const {
  return m_hydrostatic_pressure;
}

const array::Scalar& Melange::granular_pressure() const {
  return m_granular_pressure;
}

const array::Scalar& Melange::melange_mass_change() const {
  return m_melange_mass_change;
}

const array::Scalar& Melange::melange_back_pressure() const {
  return m_melange_back_pressure;
}


//! Helper methods
void Melange::compute_hydrostatic_pressure(const array::Scalar &melange_thickness) {
  const double g = m_config->get_number("constants.standard_gravity");
  const double rho_ice = m_config->get_number("constants.ice.density");
  const double packing_density = m_config->get_number("melange.packing_density");
  const double rho_melange = rho_ice * packing_density;
  
  // Hydrostatic pressure:
  // P = rho_melange * g * h_melange
  // where rho_melange = rho_ice * packing_density
  m_hydrostatic_pressure.copy_from(melange_thickness);
  m_hydrostatic_pressure.scale(g * rho_melange);
}

void Melange::compute_granular_pressure(const array::Scalar &melange_thickness) {
  const double g = m_config->get_number("constants.standard_gravity");
  const double rho_ice = m_config->get_number("constants.ice.density");
  const double rho_water = m_config->get_number("constants.water.density");
  const double packing_density = m_config->get_number("melange.packing_density");
  const double rho_melange = rho_ice * packing_density;
  
  // Granular static pressure: 
  // P = 0.5 rho_melange g ( 1 - rho_melange / rho_water ) h_melange
  // where rho_melange = rho_ice * packing_density
  const double granular_factor = 0.5 * (1 - rho_melange / rho_water);
  m_granular_pressure.copy_from(melange_thickness);
  m_granular_pressure.scale(g * rho_melange * granular_factor);
}

void Melange::compute_back_pressure(const array::Scalar &melange_thickness,
                                    const Geometry &geometry,
                                    const array::Scalar &water_column_pressure,
                                    array::Scalar &result) const {
  // Compute back pressure on ice shelf from melange
  // This is where melange affects ice dynamics
  // TODO: Implement back pressure calculation
  
  // For now, set result to zero
  result.set(0.0);
}

//! Default implementations for virtual methods
void Melange::restart_impl(const File &input_file, int record) {
  // Default implementation - read from file
  m_melange_thickness.read(input_file, record);
  m_melange_velocity.read(input_file, record);
  
  // Compute pressure from melange thickness
  compute_hydrostatic_pressure(m_melange_thickness);
  compute_granular_pressure(m_melange_thickness);
}

void Melange::bootstrap_impl(const File &input_file) {
  // Default implementation - initialize to zero
  m_melange_thickness.set(0.0);
  m_melange_velocity.set(0.0);
  
  // Compute pressure from thickness 
  compute_hydrostatic_pressure(m_melange_thickness);
  compute_granular_pressure(m_melange_thickness);
}

void Melange::init_impl(const array::Scalar &melange_thickness,
                        const array::Vector &melange_velocity) {
  // Default implementation - copy input
  m_melange_thickness.copy_from(melange_thickness);
  m_melange_velocity.copy_from(melange_velocity);
  
  // Compute pressure from thickness 
  compute_hydrostatic_pressure(m_melange_thickness);
  compute_granular_pressure(m_melange_thickness);
}

std::map<std::string, Diagnostic::Ptr> Melange::diagnostics_impl() const {
  // Default implementation - no diagnostics
  return {};
}

void Melange::define_model_state_impl(const File &output) const {
  // Define output variables
  m_melange_thickness.define(output, io::PISM_DOUBLE);
  m_melange_velocity.define(output, io::PISM_DOUBLE);
  m_hydrostatic_pressure.define(output, io::PISM_DOUBLE);
  m_granular_pressure.define(output, io::PISM_DOUBLE);
  m_melange_mass_change.define(output, io::PISM_DOUBLE);
  m_melange_back_pressure.define(output, io::PISM_DOUBLE);
}

void Melange::write_model_state_impl(const File &output) const {
  // Write state variables
  m_melange_thickness.write(output);
  m_melange_velocity.write(output);
  m_hydrostatic_pressure.write(output);
  m_granular_pressure.write(output);
  m_melange_mass_change.write(output);
  m_melange_back_pressure.write(output);
}

} // end of namespace melange
} // end of namespace pism
