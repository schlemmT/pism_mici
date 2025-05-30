// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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

#include <gsl/gsl_poly.h>
#include <cassert>

#include "pism/coupler/ocean/GivenTH.hh"
#include "pism/coupler/util/options.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Grid.hh"
#include "pism/util/Time.hh"
#include "pism/util/array/Forcing.hh"

namespace pism {
namespace ocean {

GivenTH::Constants::Constants(const Config &config) {
  // coefficients of the in situ melting point temperature
  // parameterization:
  a[0] = -0.0575;
  a[1] =  0.0901;
  a[2] = -7.61e-4;
  // coefficients of the in situ melting point potential temperature
  // parameterization:
  b[0] = -0.0575;
  b[1] =  0.0921;
  b[2] = -7.85e-4;

  // FIXME: this should not be hard-wired. Eventually we should be able
  // to use the spatially-variable top-of-the-ice temperature.
  shelf_top_surface_temperature    = -20.0; // degrees Celsius

  gamma_T                          = config.get_number("ocean.th.gamma_T");  
  gamma_S                          = config.get_number("ocean.th.gamma_S");  
  water_latent_heat_fusion         = config.get_number("constants.fresh_water.latent_heat_of_fusion");
  sea_water_density                = config.get_number("constants.sea_water.density");
  sea_water_specific_heat_capacity = config.get_number("constants.sea_water.specific_heat_capacity");
  ice_density                      = config.get_number("constants.ice.density");
  ice_specific_heat_capacity       = config.get_number("constants.ice.specific_heat_capacity");
  ice_thermal_diffusivity          = config.get_number("constants.ice.thermal_conductivity") / (ice_density * ice_specific_heat_capacity);
  limit_salinity_range             = config.get_flag("ocean.th.clip_salinity");
}

GivenTH::GivenTH(std::shared_ptr<const Grid> g)
  : CompleteOceanModel(g, std::shared_ptr<OceanModel>()) {

  ForcingOptions opt(*m_grid->ctx(), "ocean.th");

  {
    unsigned int buffer_size = m_config->get_number("input.forcing.buffer_size");

    File file(m_grid->com, opt.filename, io::PISM_NETCDF3, io::PISM_READONLY);

    m_theta_ocean = std::make_shared<array::Forcing>(m_grid,
                                                file,
                                                "theta_ocean",
                                                "", // no standard name
                                                buffer_size,
                                                opt.periodic,
                                                LINEAR);

    m_salinity_ocean = std::make_shared<array::Forcing>(m_grid,
                                                   file,
                                                   "salinity_ocean",
                                                   "", // no standard name
                                                   buffer_size,
                                                   opt.periodic,
                                                   LINEAR);
  }

  m_theta_ocean->metadata(0)
      .long_name("potential temperature of the adjacent ocean")
      .units("kelvin");

  m_salinity_ocean->metadata(0)
      .long_name("salinity of the adjacent ocean")
      .units("g/kg");
}

void GivenTH::init_impl(const Geometry &geometry) {

  m_log->message(2,
             "* Initializing the 3eqn melting parameterization ocean model\n"
             "  reading ocean temperature and salinity from a file...\n");

  ForcingOptions opt(*m_grid->ctx(), "ocean.th");

  // potential temperature is required
  m_theta_ocean->init(opt.filename, opt.periodic);

  // read ocean salinity from a file if present, otherwise use a constant
  {
    File input(m_grid->com, opt.filename, io::PISM_GUESS, io::PISM_READONLY);

    auto variable_name = m_salinity_ocean->metadata().get_name();

    if (input.variable_exists(variable_name)) {
      m_salinity_ocean->init(opt.filename, opt.periodic);
    } else {
      double salinity = m_config->get_number("constants.sea_water.salinity", "g / kg");

      m_salinity_ocean = array::Forcing::Constant(m_grid, variable_name, salinity);

      m_log->message(2, "  Variable '%s' not found; using constant salinity: %f (g / kg).\n",
                     variable_name.c_str(), salinity);
    }
  }

  // read time-independent data right away:
  if (m_theta_ocean->buffer_size() == 1 && m_salinity_ocean->buffer_size() == 1) {
    update(geometry, time().current(), 0); // dt is irrelevant
  }

  const double
    ice_density   = m_config->get_number("constants.ice.density"),
    water_density = m_config->get_number("constants.sea_water.density"),
    g             = m_config->get_number("constants.standard_gravity");

  compute_average_water_column_pressure(geometry, ice_density, water_density, g,
                                           *m_water_column_pressure);
}

void GivenTH::update_impl(const Geometry &geometry, double t, double dt) {
  m_theta_ocean->update(t, dt);
  m_salinity_ocean->update(t, dt);

  m_theta_ocean->average(t, dt);
  m_salinity_ocean->average(t, dt);

  Constants c(*m_config);

  const array::Scalar &ice_thickness = geometry.ice_thickness;

  array::Scalar &temperature = *m_shelf_base_temperature;
  array::Scalar &mass_flux = *m_shelf_base_mass_flux;

  array::AccessScope list{ &ice_thickness, m_theta_ocean.get(), m_salinity_ocean.get(),
      &temperature, &mass_flux};

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    double potential_temperature_celsius = (*m_theta_ocean)(i,j) - 273.15;

    double
      shelf_base_temp_celsius = 0.0,
      shelf_base_massflux     = 0.0;

    pointwise_update(c,
                     (*m_salinity_ocean)(i,j),
                     potential_temperature_celsius,
                     ice_thickness(i,j),
                     &shelf_base_temp_celsius,
                     &shelf_base_massflux);

    // Convert from Celsius to kelvin:
    temperature(i,j) = shelf_base_temp_celsius + 273.15;
    mass_flux(i,j)   = shelf_base_massflux;
  }

  // convert mass flux from [m s-1] to [kg m-2 s-1]:
  m_shelf_base_mass_flux->scale(m_config->get_number("constants.ice.density"));

  const double
    ice_density   = m_config->get_number("constants.ice.density"),
    water_density = m_config->get_number("constants.sea_water.density"),
    g             = m_config->get_number("constants.standard_gravity");

  compute_average_water_column_pressure(geometry, ice_density, water_density, g,
                                           *m_water_column_pressure);
}

MaxTimestep GivenTH::max_timestep_impl(double t) const {
  (void) t;

  return MaxTimestep("ocean th");
}

//* Evaluate the parameterization of the melting point temperature.
/** The value returned is in degrees Celsius.
 */
static double melting_point_temperature(GivenTH::Constants c,
                                        double salinity, double ice_thickness) {
  return c.a[0] * salinity + c.a[1] + c.a[2] * ice_thickness;
}

/** Melt rate, obtained by solving the salt flux balance equation.
 *
 * @param c model constants
 * @param sea_water_salinity sea water salinity
 * @param basal_salinity shelf base salinity
 *
 * @return shelf base melt rate, in [m/s]
 */
static double shelf_base_melt_rate(GivenTH::Constants c,
                                   double sea_water_salinity, double basal_salinity) {

  return c.gamma_S * c.sea_water_density * (sea_water_salinity - basal_salinity) / (c.ice_density * basal_salinity);
}

/** @brief Compute temperature and melt rate at the base of the shelf.
 * Based on [@ref HellmerOlbers1989] and [@ref HollandJenkins1999].
 *
 * See the manual for details.
 *
 * @param[in] constants model constants
 * @param[in] sea_water_salinity sea water salinity
 * @param[in] sea_water_potential_temperature sea water potential temperature
 * @param[in] thickness ice shelf thickness
 * @param[out] shelf_base_temperature_out resulting basal temperature
 * @param[out] shelf_base_melt_rate_out resulting basal melt rate
 *
 * @return 0 on success
 */
void GivenTH::pointwise_update(const Constants &constants, double sea_water_salinity,
                               double sea_water_potential_temperature, double thickness,
                               double *shelf_base_temperature_out,
                               double *shelf_base_melt_rate_out) {

  assert(thickness >= 0.0);

  // This model works for sea water salinity in the range of [4, 40]
  // psu. Ensure that input salinity is in this range.
  const double
    min_salinity = 4.0,
    max_salinity = 40.0;

  if (constants.limit_salinity_range) {
    if (sea_water_salinity < min_salinity) {
      sea_water_salinity = min_salinity;
    } else if (sea_water_salinity > max_salinity) {
      sea_water_salinity = max_salinity;
    }
  }

  double basal_salinity = sea_water_salinity;
  subshelf_salinity(constants, sea_water_salinity, sea_water_potential_temperature,
                    thickness, &basal_salinity);

  // Clip basal salinity so that we can use the freezing point
  // temperature parameterization to recover shelf base temperature.
  if (constants.limit_salinity_range) {
    if (basal_salinity <= min_salinity) {
      basal_salinity = min_salinity;
    } else if (basal_salinity >= max_salinity) {
      basal_salinity = max_salinity;
    }
  }

  *shelf_base_temperature_out = melting_point_temperature(constants, basal_salinity, thickness);

  *shelf_base_melt_rate_out = shelf_base_melt_rate(constants, sea_water_salinity, basal_salinity);

  // no melt if there is no ice
  if (thickness == 0.0) {
    *shelf_base_melt_rate_out = 0.0;
  }
}


/** @brief Compute the basal salinity and make sure that it is
 * consistent with the basal melt rate.
 *
 * @param[in] c constants
 * @param[in] sea_water_salinity sea water salinity
 * @param[in] sea_water_potential_temperature sea water potential temperature
 * @param[in] thickness ice shelf thickness
 * @param[out] shelf_base_salinity resulting shelf base salinity
 *
 * @return 0 on success
 */
void GivenTH::subshelf_salinity(const Constants &c, double sea_water_salinity,
                                double sea_water_potential_temperature, double thickness,
                                double *shelf_base_salinity) {

  double basal_salinity = sea_water_salinity;

  // first, assume that there is melt at the shelf base:
  {
    subshelf_salinity_melt(c, sea_water_salinity, sea_water_potential_temperature,
                           thickness, &basal_salinity);

    double basal_melt_rate = shelf_base_melt_rate(c, sea_water_salinity, basal_salinity);

    if (basal_melt_rate > 0.0) {
      // computed basal melt rate is consistent with the assumption used
      // to compute basal salinity
      *shelf_base_salinity = basal_salinity;
      return;
    }
  }

  // Assuming that there is melt resulted in an inconsistent
  // (salinity, melt_rate) pair. Assume that there is freeze-on at the base.
  {
    subshelf_salinity_freeze_on(c, sea_water_salinity, sea_water_potential_temperature,
                                thickness, &basal_salinity);

    double basal_melt_rate = shelf_base_melt_rate(c, sea_water_salinity, basal_salinity);

    if (basal_melt_rate < 0.0) {
      // computed basal melt rate is consistent with the assumption
      // used to compute basal salinity
      *shelf_base_salinity = basal_salinity;
      return;
    }
  }

  // Both assumptions (above) resulted in inconsistencies. Revert to
  // the "diffusion-only" case, which may be less accurate, but is
  // generic and is always consistent.
  {
    subshelf_salinity_diffusion_only(c, sea_water_salinity, sea_water_potential_temperature,
                                     thickness, &basal_salinity);

    *shelf_base_salinity = basal_salinity;
  }
}

/** Compute basal salinity in the basal melt case.
 *
 * We use the parameterization of the temperature gradient from [@ref
 * Hellmeretal1998], equation 13:
 *
 * @f[ T_{\text{grad}} = -\Delta T\, \frac{\frac{\partial h}{\partial t}}{\kappa}, @f]
 *
 * where @f$ \Delta T @f$ is the difference between the ice
 * temperature at the top of the ice column and its bottom:
 * @f$ \Delta T = T^S - T^B. @f$ With this parameterization, we have
 *
 * @f[ Q_T^I = \rho_I\, c_{pI}\, {\frac{\partial h}{\partial t}}\, (T^S - T^B). @f]
 *
 * Then the coefficients of the quadratic equation for basal salinity
 * (see pointwise_update()) are
 *
 * @f{align*}{
 * A &= a_{0}\,\gamma_S\,c_{pI}-b_{0}\,\gamma_T\,c_{pW}\\
 * B &= \gamma_S\,\left(L-c_{pI}\,\left(T^S+a_{0}\,S^W-a_{2}\,h-a_{1}\right)\right)+
 *      \gamma_T\,c_{pW}\,\left(\Theta^W-b_{2}\,h-b_{1}\right)\\
 * C &= -\gamma_S\,S^W\,\left(L-c_{pI}\,\left(T^S-a_{2}\,h-a_{1}\right)\right)
 * @f}
 *
 * @param[in] c physical constants, stored here to avoid looking them up in a double for loop
 * @param[in] sea_water_salinity salinity of the ocean immediately adjacent to the shelf, [g/kg]
 * @param[in] sea_water_potential_temperature potential temperature of the sea water, [degrees Celsius]
 * @param[in] thickness thickness of the ice shelf, [meters]
 * @param[out] shelf_base_salinity resulting shelf base salinity
 *
 * @return 0 on success
 */
void GivenTH::subshelf_salinity_melt(const Constants &c, double sea_water_salinity,
                                     double sea_water_potential_temperature, double thickness,
                                     double *shelf_base_salinity) {

  const double
    c_pI    = c.ice_specific_heat_capacity,
    c_pW    = c.sea_water_specific_heat_capacity,
    L       = c.water_latent_heat_fusion,
    T_S     = c.shelf_top_surface_temperature,
    S_W     = sea_water_salinity,
    Theta_W = sea_water_potential_temperature;

  // We solve a quadratic equation for Sb, the salinity at the shelf
  // base.
  //
  // A*Sb^2 + B*Sb + C = 0
  const double A = c.a[0] * c.gamma_S * c_pI - c.b[0] * c.gamma_T * c_pW;
  const double B = (c.gamma_S * (L - c_pI * (T_S + c.a[0] * S_W - c.a[2] * thickness - c.a[1])) +
                    c.gamma_T * c_pW * (Theta_W - c.b[2] * thickness - c.b[1]));
  const double C = -c.gamma_S * S_W * (L - c_pI * (T_S - c.a[2] * thickness - c.a[1]));

  double S1 = 0.0, S2 = 0.0;
  const int n_roots = gsl_poly_solve_quadratic(A, B, C, &S1, &S2);

  assert(n_roots > 0);
  assert(S2 > 0.0);             // The bigger root should be positive.

  *shelf_base_salinity = S2;
}

/** Compute basal salinity in the basal freeze-on case.
 *
 * In this case we assume that the temperature gradient at the shelf base is zero:
 *
 * @f[ T_{\text{grad}} = 0. @f]
 *
 * Please see pointwise_update() for details.
 *
 * In this case the coefficients of the quadratic equation for the
 * basal salinity are:
 *
 * @f{align*}{
 * A &= -b_{0}\,\gamma_T\,c_{pW} \\
 * B &= \gamma_S\,L+\gamma_T\,c_{pW}\,\left(\Theta^W-b_{2}\,h-b_{1}\right) \\
 * C &= -\gamma_S\,S^W\,L\\
 * @f}
 *
 * @param[in] c model constants
 * @param[in] sea_water_salinity sea water salinity
 * @param[in] sea_water_potential_temperature sea water temperature
 * @param[in] thickness ice shelf thickness
 * @param[out] shelf_base_salinity resulting basal salinity
 *
 * @return 0 on success
 */
void GivenTH::subshelf_salinity_freeze_on(const Constants &c, double sea_water_salinity,
                                          double sea_water_potential_temperature, double thickness,
                                          double *shelf_base_salinity) {

  const double
    c_pW    = c.sea_water_specific_heat_capacity,
    L       = c.water_latent_heat_fusion,
    S_W     = sea_water_salinity,
    Theta_W = sea_water_potential_temperature,
    h       = thickness;

  // We solve a quadratic equation for Sb, the salinity at the shelf
  // base.
  //
  // A*Sb^2 + B*Sb + C = 0
  const double A = -c.b[0] * c.gamma_T * c_pW;
  const double B = c.gamma_S * L + c.gamma_T * c_pW * (Theta_W - c.b[2] * h - c.b[1]);
  const double C = -c.gamma_S * S_W * L;

  double S1 = 0.0, S2 = 0.0;
  const int n_roots = gsl_poly_solve_quadratic(A, B, C, &S1, &S2);

  assert(n_roots > 0);
  assert(S2 > 0.0);             // The bigger root should be positive.

  *shelf_base_salinity = S2;
}

/** @brief Compute basal salinity in the case of no basal melt and no
 * freeze-on, with the diffusion-only temperature distribution in the
 * ice column.
 *
 * In this case the temperature gradient at the base ([@ref
 * HollandJenkins1999], equation 21) is
 *
 * @f[ T_{\text{grad}} = \frac{\Delta T}{h}, @f]
 *
 * where @f$ h @f$ is the ice shelf thickness and @f$ \Delta T = T^S -
 * T^B @f$ is the difference between the temperature at the top and
 * the bottom of the shelf.
 *
 * In this case the coefficients of the quadratic equation for the basal salinity are:
 *
 * @f{align*}{
 * A &= - \frac{b_{0}\,\gamma_T\,h\,\rho_W\,c_{pW}-a_{0}\,\rho_I\,c_{pI}\,\kappa}{h\,\rho_W}\\
 * B &= \frac{\rho_I\,c_{pI}\,\kappa\,\left(T^S-a_{2}\,h-a_{1}\right)}{h\,\rho_W}
 +\gamma_S\,L+\gamma_T\,c_{pW}\,\left(\Theta^W-b_{2}\,h-b_{1}\right)\\
 * C &= -\gamma_S\,S^W\,L\\
 * @f}
 *
 * @param[in] c model constants
 * @param[in] sea_water_salinity sea water salinity
 * @param[in] sea_water_potential_temperature sea water potential temperature
 * @param[in] thickness ice shelf thickness
 * @param[out] shelf_base_salinity resulting basal salinity
 *
 * @return 0 on success
 */
void GivenTH::subshelf_salinity_diffusion_only(const Constants &c, double sea_water_salinity,
                                               double sea_water_potential_temperature,
                                               double thickness, double *shelf_base_salinity) {

  const double
    c_pI    = c.ice_specific_heat_capacity,
    c_pW    = c.sea_water_specific_heat_capacity,
    L       = c.water_latent_heat_fusion,
    T_S     = c.shelf_top_surface_temperature,
    S_W     = sea_water_salinity,
    Theta_W = sea_water_potential_temperature,
    h       = thickness,
    rho_W   = c.sea_water_density,
    rho_I   = c.ice_density,
    kappa   = c.ice_thermal_diffusivity;

  // We solve a quadratic equation for Sb, the salinity at the shelf
  // base.
  //
  // A*Sb^2 + B*Sb + C = 0
  const double A = -(c.b[0] * c.gamma_T * h * rho_W * c_pW - c.a[0] * rho_I * c_pI * kappa) / (h * rho_W);
  const double B = ((rho_I * c_pI * kappa * (T_S - c.a[2] * h - c.a[1])) / (h * rho_W) +
                    c.gamma_S * L + c.gamma_T * c_pW * (Theta_W - c.b[2] * h - c.b[1]));
  const double C = -c.gamma_S * S_W * L;

  double S1 = 0.0, S2 = 0.0;
  const int n_roots = gsl_poly_solve_quadratic(A, B, C, &S1, &S2);

  assert(n_roots > 0);
  assert(S2 > 0.0);             // The bigger root should be positive.

  *shelf_base_salinity = S2;
}

} // end of namespace ocean
} // end of namespace pism
