/* Copyright (C) 2022 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "IceModelVec2Int.hh"
#include "IceModelVec_impl.hh"

namespace pism {

IceModelVec2Int::IceModelVec2Int(IceGrid::ConstPtr grid, const std::string &name)
  : IceModelVec2S(grid, name) {
  m_impl->interpolation_type = NEAREST;
}

IceModelVec2Int::IceModelVec2Int(IceGrid::ConstPtr grid, const std::string &name,
                                 IceModelVecKind ghostedp, int width)
  : IceModelVec2S(grid, name, ghostedp, width) {
  m_impl->interpolation_type = NEAREST;
}

} // end of namespace pism
