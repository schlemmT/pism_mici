// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2019, 2020, 2023, 2024 PISM Authors
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

#include "pism/util/io/NC4_Par.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/IO_Flags.hh"

// netcdf_par.h has to be included *after* mpi.h and after netcdf.h
//
// note that we don't need to define MPI_INCLUDED because this code is built *only* if we
// have a parallel NetCDF library.
extern "C" {
#include <netcdf.h>
#include <netcdf_par.h>
}

namespace pism {
namespace io {

//! \brief Prints an error message; for debugging.
static void check(const ErrorLocation &where, int return_code) {
  if (return_code != NC_NOERR) {
    throw RuntimeError(where, nc_strerror(return_code));
  }
}

void NC4_Par::open_impl(const std::string &fname, io::Mode mode) {
  MPI_Info info = MPI_INFO_NULL;
  int stat;

  int open_mode = mode == io::PISM_READONLY ? NC_NOWRITE : NC_WRITE;
  open_mode = open_mode | NC_MPIIO;

  stat = nc_open_par(fname.c_str(), open_mode, m_com, info, &m_file_id);

  check(PISM_ERROR_LOCATION, stat);
}

void NC4_Par::create_impl(const std::string &fname) {
  MPI_Info info = MPI_INFO_NULL;
  int stat;

  stat = nc_create_par(fname.c_str(),
                       NC_NETCDF4 | NC_MPIIO,
                       m_com, info, &m_file_id);

  check(PISM_ERROR_LOCATION, stat);
}

void NC4_Par::set_access_mode(int varid) const {
  int stat;

  // Use collective parallel access mode because it is faster.
  stat = nc_var_par_access(m_file_id, varid, NC_COLLECTIVE);
  check(PISM_ERROR_LOCATION, stat);
}

void NC4_Par::set_compression_level_impl(int level) const {
  m_compression_level = level;
}


} // end of namespace io
} // end of namespace pism
