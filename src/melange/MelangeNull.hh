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

#ifndef _MELANGENULL_H_
#define _MELANGENULL_H_

#include "pism/melange/Melange.hh"

namespace pism {

namespace melange {

//! \brief A null melange model that does nothing.
/*!
  This is a minimal implementation that provides the melange interface
  but performs no physics. Similar to how NullTransport works for hydrology.
  
  This model:
  - Always returns zero melange thickness, pressure, density, velocity
  - Performs no formation, flow, or disintegration
  - Provides no back pressure on ice shelf
  - Useful for testing or when melange physics are not needed
*/
class MelangeNull : public Melange {
public:
  MelangeNull(std::shared_ptr<const Grid> g);
  virtual ~MelangeNull() = default;

protected:
  //! Virtual implementations - all do nothing
  virtual void restart_impl(const File &input_file, int record);
  virtual void bootstrap_impl(const File &input_file);
  virtual void init_impl(const array::Scalar &melange_thickness,
                         const array::Vector &melange_velocity);
  virtual MaxTimestep max_timestep_impl(double t) const;
  virtual void update_impl(double t, double dt, const Inputs& inputs);

private:
  void initialization_message() const;
};

} // end of namespace melange
} // end of namespace pism

#endif /* _MELANGENULL_H_ */
