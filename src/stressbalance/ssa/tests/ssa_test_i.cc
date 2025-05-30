// Copyright (C) 2010--2018, 2021, 2022, 2023, 2024 Ed Bueler, Constantine Khroulev, and David Maxwell
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

#include <memory>
static char help[] =
  "\nSSA_TESTI\n"
  "  Testing program for the finite element implementation of the SSA.\n"
  "  Does a time-independent calculation.  Does not run IceModel or a derived\n"
  "  class thereof. Uses verification test I. Also may be used in a PISM\n"
  "  software (regression) test.\n\n";

#include "pism/stressbalance/ssa/SSAFD.hh"
#include "pism/stressbalance/ssa/SSAFEM.hh"
#include "pism/stressbalance/ssa/SSATestCase.hh"
#include "pism/util/Context.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/petscwrappers/PetscInitializer.hh"
#include "pism/util/pism_options.hh"
#include "pism/verification/tests/exactTestsIJ.h"

namespace pism {
namespace stressbalance {

const double m_schoof = 10; // (pure number)
const double L_schoof = 40e3; // meters
const double aspect_schoof = 0.05; // (pure)
const double H0_schoof = aspect_schoof * L_schoof;
                                       // = 2000 m THICKNESS
const double B_schoof = 3.7e8; // Pa s^{1/3}; hardness
                                     // given on p. 239 of Schoof; why so big?

std::shared_ptr<Grid> ssa_test_i_grid(std::shared_ptr<Context> ctx, int Mx, int My) {
  return SSATestCase::grid(ctx, Mx, My,
                           std::max(60.0e3, ((Mx - 1) / 2) * (2.0 * (3.0 * L_schoof) / (My - 1))),
                           3.0 * L_schoof, grid::CELL_CORNER, grid::NOT_PERIODIC);
}

class SSATestCaseI: public SSATestCase {
public:
  SSATestCaseI(std::shared_ptr<SSA> ssa) : SSATestCase(ssa) {
    EnthalpyConverter EC(*m_config);
    // 0.01 water fraction
    m_ice_enthalpy.set(EC.enthalpy(273.15, 0.01, 0.0));
  }

protected:
  virtual void initializeSSACoefficients();

  virtual void exactSolution(int i, int j,
    double x, double y, double *u, double *v);

};

void SSATestCaseI::initializeSSACoefficients() {

  m_bc_mask.set(0);
  m_geometry.ice_thickness.set(H0_schoof);

  array::AccessScope list{&m_tauc, &m_bc_values, &m_bc_mask, &m_geometry.ice_surface_elevation, &m_geometry.bed_elevation};

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    // Evaluate the exact solution and yield stress. Exact u, v will only be used at the
    // grid edge.
    auto I = exactI(m_schoof, m_grid->x(i), m_grid->y(j));

    m_geometry.bed_elevation(i, j) = I.bed;

    m_tauc(i,j) = I.tauc;

    if (grid::domain_edge(*m_grid, i, j)) {
      m_bc_mask(i, j)   = 1;
      m_bc_values(i, j) = { I.u, I.v };
    }
  }

  m_geometry.ensure_consistency(0.0);

  m_tauc.update_ghosts();
  m_bc_mask.update_ghosts();
  m_bc_values.update_ghosts();
}


void SSATestCaseI::exactSolution(int /*i*/, int /*j*/,
                                 double x, double y,
                                 double *u, double *v) {
  auto I = exactI(m_schoof, x, y);
  *u = I.u;
  *v = I.v;
}

} // end of namespace stressbalance
} // end of namespace pism

int main(int argc, char *argv[]) {

  using namespace pism;
  using namespace pism::stressbalance;

  MPI_Comm com = MPI_COMM_WORLD;
  petsc::Initializer petsc(argc, argv, help);

  com = PETSC_COMM_WORLD;

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  try {
    std::shared_ptr<Context> ctx = context_from_options(com, "ssa_testi");
    Config::Ptr config = ctx->config();

    std::string usage = "\n"
      "usage of SSA_TESTi:\n"
      "  run ssa_testi -Mx <number> -My <number> -ssa_method <fd|fem>\n"
      "\n";

    bool stop = show_usage_check_req_opts(*ctx->log(), "ssa_testi", {}, usage);

    if (stop) {
      return 0;
    }

    // Parameters that can be overridden by command line options
    unsigned int Mx = config->get_number("grid.Mx");
    unsigned int My = config->get_number("grid.My");

    auto method = config->get_string("stress_balance.ssa.method");
    auto output_file = config->get_string("output.file");

    bool write_output = config->get_string("output.size") != "none";

    // we have to set parameters *before* `ssa` is allocated
    config->set_flag("basal_resistance.pseudo_plastic.enabled", false);

    config->set_string("stress_balance.ssa.flow_law", "isothermal_glen");
    config->set_number("flow_law.isothermal_Glen.ice_softness",
                       pow(B_schoof, -config->get_number("stress_balance.ssa.Glen_exponent")));

    // The finite difference code uses the following flag to treat the non-periodic grid correctly.
    config->set_flag("stress_balance.ssa.compute_surface_gradient_inward", true);
    config->set_number("stress_balance.ssa.epsilon", 0.0); // don't use this lower bound

    auto grid = ssa_test_i_grid(ctx, Mx, My);
    SSATestCaseI testcase(SSATestCase::solver(grid, method));
    testcase.init();
    testcase.run();
    testcase.report("I");
    if (write_output) {
      testcase.write(output_file);
    }
  }
  catch (...) {
    handle_fatal_errors(com);
    return 1;
  }

  return 0;
}
