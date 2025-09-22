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

#ifndef _PISMMELANGE_H_
#define _PISMMELANGE_H_

#include "pism/util/array/Vector.hh"
#include "pism/util/Component.hh"

namespace pism {

//! @brief Melange models and related diagnostics.
namespace melange {

//! \brief Input data structure for melange models.
/*!
  This structure contains all the input data that melange models need
  to compute their state. Similar to hydrology::Inputs.
*/
class Inputs {
public:
  Inputs();

  //! a mask array used in  regional modeling mode to define which parts of the computational domain should be excluded from certain physics calculations
  const array::Scalar1      *no_model_mask;
  
  //! ice sheet geometry
  const Geometry *geometry;
  
  //! ice thickness [m]
  const array::Scalar        *ice_thickness;
  
  //! surface elevation [m]
  const array::Scalar        *surface_elevation;
  
  //! cell type mask (grounded, floating, ice-free ocean, etc.) - from geometry
  //! Available mask values: MASK_ICE_FREE_BEDROCK, MASK_GROUNDED, MASK_FLOATING, MASK_ICE_FREE_OCEAN
  //! Use for: identifying where melange can form (ice-free ocean), grounding line detection
  const array::CellType      *cell_type;
  
  //! bed topography [m] - for determining fjord geometry
  const array::Scalar        *bed_elevation;
  
  //! sea level elevation [m]
  const array::Scalar        *sea_level;
  
  //! ice velocity [m/s]
  const array::Vector        *ice_velocity;
  
  //! basal melt rate [m/s] - from ocean model (shelf_base_mass_flux converted to m/s)
  const array::Scalar        *basal_melt_rate;
  
  //! calving rate [m/s] - sum of all calving mechanisms (eigen, hayhurst, vonmises, cliff, linear)
  const array::Scalar        *calving_rate;
  
  //! frontal melt rate [m/s] - from frontal melt model
  const array::Scalar        *frontal_melt_rate;
  
  //! combined retreat rate [m/s] - calving_rate + frontal_melt_rate
  const array::Scalar        *retreat_rate;
};

//! \brief The PISM melange model interface.
/*!
  This is a virtual base class.

  The purpose of this class and its derived classes is to provide
  \code
  melange_thickness()      // thickness of melange layer [m]
  melange_pressure()       // pressure in melange [Pa]
  melange_velocity()       // velocity of melange [m/s]
  \endcode

  These correspond to state variables in the melange model.

  Additional modeled fields, for diagnostic purposes, are
  \code
  melange_mass_change()           // rate of mass change [kg/s]
  melange_formation_rate()        // rate of melange formation [m/s]
  melange_back_pressure()         // back pressure on ice shelf [Pa]
  \endcode

  This interface is used to solve an SSA equation for the melange velocity
  using different melange rheologies.

  The melange class uses the ice geometry, ocean state,
  calving rates, and ice velocity in determining the evolution
  of the melange state variables.
*/
class Melange : public Component {
public:
  Melange(std::shared_ptr<const Grid> g);
  virtual ~Melange() = default;

  //! Restart from a PISM output file
  void restart(const File &input_file, int record);

  //! Bootstrap from a file with ice thickness and no melange state variables
  void bootstrap(const File &input_file,
                 const array::Scalar &ice_thickness);

  //! Initialize with given state
  void init(const array::Scalar &melange_thickness,
            const array::Scalar &melange_pressure);

  //! Update the melange state
  void update(double t, double dt, const Inputs& inputs);

  //! State variables
  const array::Scalar& melange_thickness() const;
  const array::Scalar& melange_pressure() const;
  const array::Vector& melange_velocity() const;

  //! Essential diagnostic quantities
  const array::Scalar& melange_mass_change() const;
  const array::Scalar& melange_back_pressure() const;

protected:
  //! Virtual implementations 
  virtual void restart_impl(const File &input_file, int record);
  virtual void bootstrap_impl(const File &input_file,
                              const array::Scalar &ice_thickness);
  virtual void init_impl(const array::Scalar &melange_thickness,
                         const array::Scalar &melange_pressure);
  virtual void update_impl(double t, double dt, const Inputs& inputs) = 0;
  virtual std::map<std::string, Diagnostic::Ptr> diagnostics_impl() const;

  virtual void define_model_state_impl(const File &output) const;
  virtual void write_model_state_impl(const File &output) const;

  //! Helper methods - similar to Hydrology's compute_overburden_pressure()
  void compute_melange_pressure(const array::Scalar &melange_thickness,
                                const array::Scalar &ice_thickness,
                                array::Scalar &result) const;

  void compute_back_pressure(const array::Scalar &melange_thickness,
                             const array::Scalar &ice_thickness,
                             const array::Scalar &water_column_pressure,
                             array::Scalar &result) const;

protected:
  //! State variables - similar to Hydrology's m_Wtill, m_W, etc.
  
  //! thickness of melange layer [m]
  array::Scalar m_melange_thickness;
  
  //! pressure in melange [Pa]
  array::Scalar m_melange_pressure;
  
  //! velocity of melange [m/s]
  array::Vector m_melange_velocity;

  //! Essential diagnostic quantities
  //! rate of mass change [kg/s]
  array::Scalar m_melange_mass_change;
  
  //! back pressure on ice shelf [Pa]
  array::Scalar m_melange_back_pressure;


private:
  virtual void initialization_message() const = 0;
};

} // end of namespace melange
} // end of namespace pism

#endif /* _PISMMELANGE_H_ */
