// Copyright (C) 2009--2020 Constantine Khroulev
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

#include <petsc.h>
#include <cassert>

#include "iceModelVec2T.hh"
#include "pism/util/io/File.hh"
#include "pism_utilities.hh"
#include "Time.hh"
#include "IceGrid.hh"
#include "ConfigInterface.hh"

#include "error_handling.hh"
#include "io/io_helpers.hh"
#include "pism/util/Logger.hh"
#include "pism/util/interpolation.hh"
#include "pism/util/Context.hh"
#include "pism/util/IceModelVec_impl.hh"
#include "pism/util/VariableMetadata.hh"

namespace pism {

struct IceModelVec2T::Data {
  Data()
    : array(nullptr),
      N(0),
      first(-1),
      period(0),
      reference_time(0.0) {
    // empty
  }
  //! all the times available in filename
  std::vector<double> time;

  //! time bounds
  std::vector<double> time_bounds;

  //! file to read (regrid) from
  std::string filename;

  // DM with dof equal to the number of records kept in memory
  std::shared_ptr<petsc::DM> da;

  //! a 3D Vec used to store records
  petsc::Vec v;

  void ***array;

  //! maximum number of records to store in memory
  unsigned int n_records;

  //! number of records kept in memory
  unsigned int N;

  //! number of evaluations per year used to compute temporal averages
  unsigned int n_evaluations_per_year;

  //! in-file index of the first record stored in memory ("int" to allow first==-1 as an
  //! "invalid" first value)
  int first;

  InterpolationType interp_type;
  std::shared_ptr<Interpolation> interp;

  // forcing period, in years
  unsigned int period;

  // reference time, in seconds
  double reference_time;
};

/*!
 * Allocate an instance that will be used to load and use a forcing field from a file.
 *
 * Checks the number of records in a file and allocates storage accordingly.
 *
 * If `periodic` is true, allocate enough storage to hold all the records, otherwise
 * allocate storage for at most `max_buffer_size` records.
 *
 * @param[in] grid computational grid
 * @param[in] file input file
 * @param[in] short_name variable name in `file`
 * @param[in] standard_name standard name (if available); leave blank to ignore
 * @param[in] max_buffer_size maximum buffer size for non-periodic fields
 * @param[in] evaluations_per_year number of evaluations per year to use when averaging
 * @param[in] periodic true if this forcing field should be interpreted as periodic
 */
IceModelVec2T::Ptr IceModelVec2T::ForcingField(IceGrid::ConstPtr grid,
                                               const File &file,
                                               const std::string &short_name,
                                               const std::string &standard_name,
                                               int max_buffer_size,
                                               int evaluations_per_year,
                                               bool periodic,
                                               InterpolationType interpolation_type) {

  int n_records = file.nrecords(short_name, standard_name,
                                    grid->ctx()->unit_system());

  if (not periodic) {
    n_records = std::min(n_records, max_buffer_size);
  }
  // In the periodic case we try to keep all the records in RAM.

  // Allocate storage for one record if the variable was not found. This is needed to be
  // able to cheaply allocate and then discard an "-atmosphere given" model
  // (atmosphere::Given) when "-surface given" (Given) is selected.
  n_records = std::max(n_records, 1);

  // LCOV_EXCL_START
  if (n_records > IceGrid::max_dm_dof) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "cannot allocate storage for %d records of %s (%s)"
                                  " (exceeds the maximum of %d)",
                                  n_records, short_name.c_str(), standard_name.c_str(),
                                  IceGrid::max_dm_dof);
  }
  // LCOV_EXCL_STOP

  if (periodic and interpolation_type == LINEAR) {
    interpolation_type = LINEAR_PERIODIC;
  }

  return IceModelVec2T::Ptr(new IceModelVec2T(grid, short_name, n_records,
                                              evaluations_per_year, interpolation_type));
}


IceModelVec2T::IceModelVec2T(IceGrid::ConstPtr grid, const std::string &short_name,
                             unsigned int n_records,
                             unsigned int n_evaluations_per_year,
                             InterpolationType interpolation_type)
  : IceModelVec2S(grid, short_name, WITHOUT_GHOSTS, 1),
    m_data(new Data())
{
  m_impl->report_range = false;

  m_data->interp_type = interpolation_type;
  m_data->n_records = n_records;
  m_data->n_evaluations_per_year = n_evaluations_per_year;

  if (not (m_data->interp_type == PIECEWISE_CONSTANT or
           m_data->interp_type == LINEAR or
           m_data->interp_type == LINEAR_PERIODIC)) {
    throw RuntimeError(PISM_ERROR_LOCATION, "unsupported interpolation type");
  }

  // LCOV_EXCL_START
  if (n_records > IceGrid::max_dm_dof) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "cannot allocate storage for %d records of %s"
                                  " (exceeds the maximum of %d)",
                                  n_records, short_name.c_str(), IceGrid::max_dm_dof);
  }
  // LCOV_EXCL_STOP

  // initialize the m_data->da member:
  m_data->da = m_impl->grid->get_dm(n_records, this->m_impl->da_stencil_width);

  // allocate the 3D Vec:
  PetscErrorCode ierr = DMCreateGlobalVector(*m_data->da, m_data->v.rawptr());
  PISM_CHK(ierr, "DMCreateGlobalVector");
}

IceModelVec2T::~IceModelVec2T() {
  delete m_data;
}

unsigned int IceModelVec2T::n_records() {
  return m_data->n_records;
}

double*** IceModelVec2T::array3() {
  return reinterpret_cast<double***>(m_data->array);
}

void IceModelVec2T::begin_access() const {
  if (m_impl->access_counter == 0) {
    PetscErrorCode ierr = DMDAVecGetArrayDOF(*m_data->da, m_data->v, &m_data->array);
    PISM_CHK(ierr, "DMDAVecGetArrayDOF");
  }

  // this call will increment the m_access_counter
  IceModelVec2S::begin_access();
}

void IceModelVec2T::end_access() const {
  // this call will decrement the m_access_counter
  IceModelVec2S::end_access();

  if (m_impl->access_counter == 0) {
    PetscErrorCode ierr = DMDAVecRestoreArrayDOF(*m_data->da, m_data->v, &m_data->array);
    PISM_CHK(ierr, "DMDAVecRestoreArrayDOF");
    m_data->array = NULL;
  }
}

void IceModelVec2T::init(const std::string &fname, unsigned int period, double reference_time) {

  auto ctx = m_impl->grid->ctx();
  const Logger &log = *ctx->log();

  m_data->filename       = fname;
  m_data->period         = period;
  m_data->reference_time = reference_time;

  // We find the variable in the input file and
  // try to find the corresponding time dimension.

  File file(m_impl->grid->com, m_data->filename, PISM_GUESS, PISM_READONLY);
  auto var = file.find_variable(m_impl->metadata[0].get_name(),
                                m_impl->metadata[0].get_string("standard_name"));
  if (not var.exists) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "can't find %s (%s) in %s.",
                                  m_impl->metadata[0].get_string("long_name").c_str(),
                                  m_impl->metadata[0].get_name().c_str(),
                                  m_data->filename.c_str());
  }

  auto time_name = io::time_dimension(ctx->unit_system(),
                                      file, var.name);

  if (not time_name.empty()) {
    // we're found the time dimension
    TimeseriesMetadata time_dimension(time_name, time_name, ctx->unit_system());

    auto time_units = ctx->time()->units_string();
    time_dimension.set_string("units", time_units);

    io::read_timeseries(file, time_dimension,
                        *ctx->time(), log, m_data->time);

    std::string bounds_name = file.read_text_attribute(time_name, "bounds");

    if (m_data->time.size() > 1) {

      if (m_data->interp_type == PIECEWISE_CONSTANT) {
        if (bounds_name.empty()) {
          // no time bounds attribute
          throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                        "Variable '%s' does not have the time_bounds attribute.\n"
                                        "Cannot use time-dependent forcing data '%s' (%s) without time bounds.",
                                        time_name.c_str(),
                                        m_impl->metadata[0].get_string("long_name").c_str(),
                                        m_impl->metadata[0].get_name().c_str());
        }

        // read time bounds data from a file
        TimeBoundsMetadata tb(bounds_name, time_name, ctx->unit_system());
        tb.set_string("units", time_units);

        io::read_time_bounds(file, tb, *ctx->time(), log, m_data->time_bounds);

        // time bounds data overrides the time variable: we make t[j] be the
        // left end-point of the j-th interval
        for (unsigned int k = 0; k < m_data->time.size(); ++k) {
          m_data->time[k] = m_data->time_bounds[2*k + 0];
        }
      } else {
        // fake time step length used to generate the right end point of the last interval
        // TODO: figure out if there is a better way to do this.
        double dt = 1.0;
        size_t N = m_data->time.size();
        m_data->time_bounds.resize(2 * N);
        for (size_t k = 0; k < N; ++k) {
          m_data->time_bounds[2 * k + 0] = m_data->time[k];
          m_data->time_bounds[2 * k + 1] = k + 1 < N ? m_data->time[k + 1] : m_data->time[k] + dt;
        }
      }

    } else {
      // only one time record; set fake time bounds:
      m_data->time_bounds = {m_data->time[0] - 1.0, m_data->time[0] + 1};
    }

  } else {
    // no time dimension; assume that we have only one record and set the time
    // to 0
    m_data->time = {0.0};

    // set fake time bounds:
    m_data->time_bounds = {-1.0, 1.0};
  }

  if (not is_increasing(m_data->time)) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "times have to be strictly increasing (read from '%s').",
                                  m_data->filename.c_str());
  }

  if (m_data->period != 0) {
    if ((size_t)m_data->n_records < m_data->time.size()) {
      throw RuntimeError(PISM_ERROR_LOCATION,
                         "buffer has to be big enough to hold all records of periodic data");
    }

    // read periodic data right away (we need to hold it all in memory anyway)
    update(0);
  }
}

//! Initialize as constant in time and space
void IceModelVec2T::init_constant(double value) {

  // set constant value everywhere
  set(value);
  set_record(0);

  // set the time to zero
  m_data->time = {0.0};
  m_data->N = 1;
  m_data->first = 0;

  // set fake time bounds:
  m_data->time_bounds = {-1.0, 1.0};
}

//! Read some data to make sure that the interval (t, t + dt) is covered.
void IceModelVec2T::update(double t, double dt) {

  if (m_data->filename.empty()) {
    // We are not reading data from a file.
    return;
  }

  if (m_data->time_bounds.size() == 0) {
    update(0);
    return;
  }

  if (m_data->period != 0) {
    // we read all data in IceModelVec2T::init() (see above)
    return;
  }

  if (m_data->N > 0) {
    unsigned int last = m_data->first + (m_data->N - 1);

    // find the interval covered by data held in memory:
    double
      t0 = m_data->time_bounds[m_data->first * 2],
      t1 = m_data->time_bounds[last * 2 + 1];

    // just return if we have all the data we need:
    if (t >= t0 and t + dt <= t1) {
      return;
    }
  }

  Interpolation I(m_data->interp_type, m_data->time, {t, t + dt});

  unsigned int
    first = I.left(0),
    last  = I.right(1),
    N     = last - first + 1;

  // check if all the records necessary to cover this interval fit in the
  // buffer:
  if (N > m_data->n_records) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "cannot read %d records of %s (buffer size: %d)",
                                  N, m_impl->name.c_str(), m_data->n_records);
  }

  update(first);
}

//! Update by reading at most n_records records from the file.
void IceModelVec2T::update(unsigned int start) {

  unsigned int time_size = (int)m_data->time.size();

  if (start >= time_size) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "IceModelVec2T::update(int start): start = %d is invalid", start);
  }

  unsigned int missing = std::min(m_data->n_records, time_size - start);

  if (start == static_cast<unsigned int>(m_data->first)) {
    // nothing to do
    return;
  }

  int kept = 0;
  if (m_data->first >= 0) {
    unsigned int last = m_data->first + (m_data->N - 1);
    if ((m_data->N > 0) && (start >= (unsigned int)m_data->first) && (start <= last)) {
      int discarded = start - m_data->first;
      kept = last - start + 1;
      discard(discarded);
      missing -= kept;
      start += kept;
      m_data->first += discarded;
    } else {
      m_data->first = start;
    }
  } else {
    m_data->first = start;
  }

  if (missing <= 0) {
    return;
  }

  m_data->N = kept + missing;

  auto t = m_impl->grid->ctx()->time();

  auto log = m_impl->grid->ctx()->log();
  if (this->n_records() > 1) {
    log->message(4,
               "  reading \"%s\" into buffer\n"
               "          (short_name = %s): %d records, time intervals (%s, %s) through (%s, %s)...\n",
               metadata().get_string("long_name").c_str(), m_impl->name.c_str(), missing,
               t->date(m_data->time_bounds[start*2]).c_str(),
               t->date(m_data->time_bounds[start*2 + 1]).c_str(),
               t->date(m_data->time_bounds[(start + missing - 1)*2]).c_str(),
               t->date(m_data->time_bounds[(start + missing - 1)*2 + 1]).c_str());
    m_impl->report_range = false;
  } else {
    m_impl->report_range = true;
  }

  File file(m_impl->grid->com, m_data->filename, PISM_GUESS, PISM_READONLY);

  const bool allow_extrapolation = m_impl->grid->ctx()->config()->get_flag("grid.allow_extrapolation");

  for (unsigned int j = 0; j < missing; ++j) {
    {
      petsc::VecArray tmp_array(m_impl->v);
      io::regrid_spatial_variable(m_impl->metadata[0], *m_impl->grid, file, start + j, CRITICAL,
                                  m_impl->report_range, allow_extrapolation,
                                  0.0, m_impl->interpolation_type, tmp_array.get());
    }

    log->message(5, " %s: reading entry #%02d, year %s...\n",
                 m_impl->name.c_str(),
                 start + j,
                 t->date(m_data->time[start + j]).c_str());

    set_record(kept + j);
  }
}

//! Discard the first N records, shifting the rest of them towards the "beginning".
void IceModelVec2T::discard(int number) {

  if (number == 0) {
    return;
  }

  m_data->N -= number;

  AccessList l{this};

  double ***a3 = array3();
  for (Points p(*m_impl->grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    for (unsigned int k = 0; k < m_data->N; ++k) {
      a3[j][i][k] = a3[j][i][k + number];
    }
  }
}

//! Sets the record number n to the contents of the (internal) Vec v.
void IceModelVec2T::set_record(int n) {

  AccessList l{this};

  double  **a2 = array();
  double ***a3 = array3();
  for (Points p(*m_impl->grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    a3[j][i][n] = a2[j][i];
  }
}

//! Sets the (internal) Vec v to the contents of the nth record.
void IceModelVec2T::get_record(int n) {

  AccessList l{this};

  double  **a2 = array();
  double ***a3 = array3();
  for (Points p(*m_impl->grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    a2[j][i] = a3[j][i][n];
  }
}

//! @brief Given the time t determines the maximum possible time-step this IceModelVec2T
//! allows.
MaxTimestep IceModelVec2T::max_timestep(double t) const {
  // only allow going to the next record

  // find the index k such that m_data->time[k] <= x < m_data->time[k + 1]
  size_t k = gsl_interp_bsearch(m_data->time.data(), t, 0, m_data->time.size());

  // end of the corresponding interval
  double
    t_next = m_data->time_bounds[2 * k + 1],
    dt     = std::max(t_next - t, 0.0);

  if (dt > 1.0) {               // never take time-steps shorter than 1 second
    return MaxTimestep(dt);
  } else if (k + 1 < m_data->time.size()) {
    dt = m_data->time_bounds[2 * (k + 1) + 1] - m_data->time_bounds[2 * (k + 1)];
    return MaxTimestep(dt);
  } else {
    return MaxTimestep();
  }
}

/*
 * \brief Use piecewise-constant interpolation to initialize
 * IceModelVec2T with the value at time `t`.
 *
 * \note This method does not check if an update() call is necessary!
 *
 * @param[in] t requested time
 *
 */
void IceModelVec2T::interp(double t) {

  init_interpolation({t});

  get_record(m_data->interp->left(0));
}


/**
 * Compute the average value over the time interval `[t, t + dt]`.
 *
 * @param t  start of the time interval, in seconds
 * @param dt length of the time interval, in seconds
 *
 */
void IceModelVec2T::average(double t, double dt) {

  double dt_years = units::convert(m_impl->grid->ctx()->unit_system(),
                                   dt, "seconds", "years"); // *not* time->year(dt)

  // if only one record, nothing to do
  if (m_data->time.size() == 1) {
    return;
  }

  // Determine the number of small time-steps to use for averaging:
  int M = (int) ceil(m_data->n_evaluations_per_year * (dt_years));
  if (M < 1) {
    M = 1;
  }

  std::vector<double> ts(M);
  double ts_dt = dt / M;
  for (int k = 0; k < M; k++) {
    ts[k] = t + k * ts_dt;
  }

  init_interpolation(ts);

  AccessList l{this};

  double **a2 = array();
  for (Points p(*m_impl->grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    a2[j][i] = average(i, j);
  }
}

/**
 * \brief Compute weights for the piecewise-constant interpolation.
 * This is used *both* for time-series and "snapshots".
 *
 * @param ts requested times, in seconds
 *
 */
void IceModelVec2T::init_interpolation(const std::vector<double> &ts) {

  assert(m_data->first >= 0);

  auto time = m_impl->grid->ctx()->time();

  // Compute "periodized" times if necessary.
  std::vector<double> times_requested(ts.size());
  if (m_data->period != 0) {
    for (unsigned int k = 0; k < ts.size(); ++k) {
      times_requested[k] = time->mod(ts[k] - m_data->reference_time, m_data->period);
    }
  } else {
    times_requested = ts;
  }

  m_data->interp.reset(new Interpolation(m_data->interp_type, &m_data->time[m_data->first], m_data->N,
                                   times_requested.data(), times_requested.size(),
                                   time->years_to_seconds(m_data->period)));
}

/**
 * \brief Compute values of the time-series using precomputed indices
 * (and piecewise-constant interpolation).
 *
 * @param i,j map-plane grid point
 * @param result pointer to an allocated array of `weights.size()` `double`
 *
 */
void IceModelVec2T::interp(int i, int j, std::vector<double> &result) {
  double ***a3 = array3();

  result.resize(m_data->interp->alpha().size());

  m_data->interp->interpolate(a3[j][i], result.data());
}

//! \brief Finds the average value at i,j over the interval (t, t +
//! dt) using the rectangle rule.
/*!
  Can (and should) be optimized. Later, though.
 */
double IceModelVec2T::average(int i, int j) {
  unsigned int M = m_data->interp->alpha().size();
  double result = 0.0;

  if (m_data->N == 1) {
    double ***a3 = array3();
    result = a3[j][i][0];
  } else {
    std::vector<double> values(M);

    interp(i, j, values);

    // rectangular rule (uses the fact that points are equally-spaced
    // in time)
    result = 0;
    for (unsigned int k = 0; k < M; ++k) {
      result += values[k];
    }
    result /= (double)M;
  }
  return result;
}



} // end of namespace pism
