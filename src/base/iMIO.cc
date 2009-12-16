// Copyright (C) 2004-2009 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <cstring>
#include <cstdio>
#include <petscda.h>
#include "iceModel.hh"
#include <algorithm>
#include <sstream>
#include <set>

//! Save model state in NetCDF format.
/*! 
Optionally allows saving of full velocity field.

Calls dumpToFile() and writeMatlabVars() to do the actual work.
 */
PetscErrorCode  IceModel::writeFiles(const char* default_filename) {
  PetscErrorCode ierr;
  char filename[PETSC_MAX_PATH_LEN];

  ierr = stampHistoryEnd(); CHKERRQ(ierr);

  PetscTruth o_set;
  ierr = PetscOptionsGetString(PETSC_NULL, "-o", filename, PETSC_MAX_PATH_LEN, &o_set); CHKERRQ(ierr);

  // Use the default if the output file name was not given:
  if (!o_set)
    strncpy(filename, default_filename, PETSC_MAX_PATH_LEN);

  if (!ends_with(filename, ".nc")) {
    ierr = verbPrintf(2, grid.com,
		      "PISM WARNING: output file name does not have the '.nc' suffix!\n");
    CHKERRQ(ierr);
  }

  ierr = verbPrintf(2, grid.com, "Writing model state to file `%s'\n", filename); CHKERRQ(ierr);
  ierr = dumpToFile(filename); CHKERRQ(ierr);

  // save the config file 
  char config_out[PETSC_MAX_PATH_LEN];
  PetscTruth dump_config;
  ierr = PetscOptionsGetString(PETSC_NULL, "-dump_config", config_out, PETSC_MAX_PATH_LEN, &dump_config);
  if (dump_config) {
    ierr = config.write(config_out); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode IceModel::dumpToFile(const char *filename) {
  PetscErrorCode ierr;
  NCTool nc(&grid);

  // Prepare the file
  ierr = nc.open_for_writing(filename, false, true); CHKERRQ(ierr);
  // append == false, check_dims == true
  ierr = nc.append_time(grid.year); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  ierr = mapping.write(filename); CHKERRQ(ierr);
  ierr = global_attributes.write(filename); CHKERRQ(ierr);

  PetscTruth override_used;
  ierr = check_option("-config_override", override_used); CHKERRQ(ierr);
  if (override_used) {
    overrides.update_from(config);
    ierr = overrides.write(filename); CHKERRQ(ierr);
  }

  ierr = write_model_state(filename);  CHKERRQ(ierr);

  if (atmosPCC != PETSC_NULL) {
    ierr = atmosPCC->writeCouplingFieldsToFile(grid.year,filename); CHKERRQ(ierr);
  } else {
    SETERRQ(1,"PISM ERROR: atmosPCC == PETSC_NULL");
  }

  if (oceanPCC != PETSC_NULL) {
    ierr = oceanPCC->writeCouplingFieldsToFile(grid.year,filename); CHKERRQ(ierr);
  } else {
    SETERRQ(1,"PISM ERROR: oceanPCC == PETSC_NULL");
  }

  ierr = write_extra_fields(filename); CHKERRQ(ierr); // chance for derived classes to do more

  return 0;
}

//! \brief Writes variables listed in vars to filename, using nctype to write
//! fields stored in dedicated IceModelVecs.
PetscErrorCode IceModel::write_variables(const char *filename, set<string> vars,
					 nc_type nctype) {
  PetscErrorCode ierr;
  IceModelVec *v;

  set<string>::iterator i = vars.begin();
  while (i != vars.end()) {
    v = variables.get(*i);

    if (v == NULL) {
      ++i;
    } else {
      if (*i == "mask") {
	ierr = v->write(filename); CHKERRQ(ierr); // use the default data type
      } else {
	ierr = v->write(filename, nctype); CHKERRQ(ierr);
      }

      vars.erase(i++);		// note that it only erases variables that were
				// found (and saved)
    }
  }

  // All the remaining names in vars must be of diagnostic quantities.
  i = vars.begin();
  while (i != vars.end()) {
    ierr = compute_by_name(*i, v); CHKERRQ(ierr);

    if (v == NULL)
      ++i;
    else {
      ierr = v->write(filename, NC_FLOAT); CHKERRQ(ierr); // diagnostic quantities are always written in float
      vars.erase(i++);
    }
  }

  // check if we have any variables we didn't write
  if (!vars.empty()) {
    ierr = verbPrintf(2, grid.com,
		      "PISM WARNING: skipping the following variables: "); CHKERRQ(ierr);
    for (i = vars.begin(); i != vars.end(); ++i) {
      ierr = verbPrintf(2, grid.com, "%s, ", (*i).c_str()); CHKERRQ(ierr);
    }
    ierr = verbPrintf(2, grid.com, "\b\b\n"); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode IceModel::write_model_state(const char* filename) {
  PetscErrorCode ierr;

  string tmp = config.get_string("output_variables");
  istringstream list(tmp);
  set<string> vars;
  
  // split the list; note that this also removes any duplicate entries
  while (getline(list, tmp, ' ')) {
    if (!tmp.empty())		// this ignores multiple spaces separating variable names
      vars.insert(tmp);
  }

  // add more variables (if needed)
  if (config.get_flag("use_ssa_velocity")) {
    vars.insert("vubarSSA");
    vars.insert("vvbarSSA");
  }

  if (config.get_flag("force_full_diagnostics")) {
    ierr = verbPrintf(2, grid.com, "Writing full 3D velocities...\n"); CHKERRQ(ierr);
    vars.insert("uvel");
    vars.insert("vvel");
    vars.insert("wvel");
    vars.insert("uvelsurf");
    vars.insert("vvelsurf");
  }

  if (config.get_flag("do_age"))
    vars.insert("age");

  // FIXME: temporarily, so that we can compare to IceEnthalpyModel results;
  //   what to do with pressure-adjusted temp in longer term?
  PetscTruth write_temp_pa;
  ierr = check_option("-temp_pa", write_temp_pa); CHKERRQ(ierr);  
  if (write_temp_pa == PETSC_TRUE) {
    // write temp_pa = pressure-adjusted temp in Celcius
    //   use Tnew3 (global) as temporary, allocated space for this purpose
    ierr = verbPrintf(2, grid.com,
      "  writing pressure-adjusted ice temperature (deg C) 'temp_pa' ...\n"); CHKERRQ(ierr);
    vars.insert("temp_pa");
  }

  ierr = write_variables(filename, vars, NC_DOUBLE);

  return 0;
}


//! Writes extra fields to the output file \c filename. Does nothing in the base class.
PetscErrorCode IceModel::write_extra_fields(const char* /*filename*/) {
  // Do nothing.
  return 0;
}

//! Read a saved PISM model state in NetCDF format, for complete initialization of an evolution or diagnostic run.
/*! When initializing from a NetCDF input file, the input file determines the
  number of grid points (Mx,My,Mz,Mbz) and the dimensions (Lx,Ly,Lz) of the
  computational box.
 */
PetscErrorCode IceModel::initFromFile(const char *filename) {
  PetscErrorCode  ierr;
  NCTool nc(&grid);

  ierr = verbPrintf(2, grid.com, "initializing from NetCDF file '%s'...\n",
		    filename); CHKERRQ(ierr);

  ierr = nc.open_for_reading(filename); CHKERRQ(ierr);

  // Find the index of the last record in the file:
  int last_record;
  ierr = nc.get_dim_length("t", &last_record); CHKERRQ(ierr);
  last_record -= 1;

  // Read the model state, mapping and climate_steady variables:
  set<IceModelVec*> vars = variables.get_variables();

  set<IceModelVec*>::iterator i = vars.begin();
  while (i != vars.end()) {

    string intent = (*i)->string_attr("pism_intent");
    if ((intent == "model_state") || (intent == "mapping") ||
	(intent == "climate_steady")) {
      ierr = (*i)->read(filename, last_record); CHKERRQ(ierr);
    }

    ++i;
  }
 
  // Read vubarSSA and vvbarSSA if SSA is on, if not asked to ignore them and
  // if they are present in the input file.
  bool have_ssa_velocities = false;
  if (config.get_flag("use_ssa_velocity")) {
    string word;
    ierr = nc.get_att_text(NC_GLOBAL, "pism_ssa_velocities_are_valid", word); CHKERRQ(ierr);

    have_ssa_velocities = (word == "true") || (word == "yes") || (word == "on");
  }

  PetscTruth dontreadSSAvels = PETSC_FALSE;
  ierr = check_option("-dontreadSSAvels", dontreadSSAvels); CHKERRQ(ierr);
  
  if (have_ssa_velocities && (!dontreadSSAvels)) {
    ierr = verbPrintf(3,grid.com,"Reading vubarSSA and vvbarSSA...\n"); CHKERRQ(ierr);

    ierr = vubarSSA.read(filename, last_record); CHKERRQ(ierr);
    ierr = vvbarSSA.read(filename, last_record); CHKERRQ(ierr);
  }

  // read mapping parameters if present
  bool mapping_exists;
  ierr = nc.find_variable("mapping", NULL, mapping_exists); CHKERRQ(ierr);
  if (mapping_exists) {
    ierr = mapping.read(filename); CHKERRQ(ierr);
    ierr = mapping.print(); CHKERRQ(ierr);
  }

  string history;
  ierr = nc.get_att_text(NC_GLOBAL, "history", history); CHKERRQ(ierr);
  global_attributes.prepend_history(history);

  ierr = nc.close(); CHKERRQ(ierr);

  return 0;
}

//! Manage regridding based on user options.  Call IceModelVec::regrid() to do each selected variable.
/*!
For each variable selected by option <tt>-regrid_vars</tt>, we regrid it onto the current grid from 
the NetCDF file specified by <tt>-regrid_from</tt>.

The default, if <tt>-regrid_vars</tt> is not given, is to regrid the 3 dimensional 
quantities \c tau3, \c T3, \c Tb3.  This is consistent with one standard purpose of 
regridding, which is to stick with current geometry through the downscaling procedure.  
Most of the time the user should carefully specify which variables to regrid.
 */
PetscErrorCode IceModel::regrid() {
  PetscErrorCode ierr;
  char filename[PETSC_MAX_PATH_LEN], tmp[TEMPORARY_STRING_LENGTH];
  PetscTruth regridVarsSet, regrid_from_set;
  NCTool nc(&grid);

  ierr = check_old_option_and_stop(grid.com, "-regrid", "-regrid_from"); CHKERRQ(ierr);

  ierr = PetscOptionsBegin(grid.com, PETSC_NULL,
			   "Options controlling regridding",
			   PETSC_NULL); CHKERRQ(ierr);

  // Get the regridding file name:
  ierr = PetscOptionsString("-regrid_from", "Specifies the file to regrid from", "", "",
			    filename, PETSC_MAX_PATH_LEN,
			    &regrid_from_set); CHKERRQ(ierr);

  ierr = PetscOptionsString("-regrid_vars", "Specifies the list of variable to regrid", "",
			    "age,temp,litho_temp",
			    tmp, TEMPORARY_STRING_LENGTH,
			    &regridVarsSet); CHKERRQ(ierr);

  // Done with the options.
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  // Return if no regridding is requested:
  if (!regrid_from_set) return 0;

  ierr = verbPrintf(2, grid.com, "regridding from file %s ...\n",filename); CHKERRQ(ierr);
  
  string var_name;
  set<string> vars;
  istringstream list(tmp);
  if (regridVarsSet) {
    // split the list; note that this also removes any duplicate entries
    while (getline(list, var_name, ','))
      vars.insert(var_name);
  } else {
    vars.insert("age");
    vars.insert("temp");
    vars.insert("litho_temp");
  }

  // create "local interpolation context" from dimensions, limits, and lengths
  //   extracted from regridFile, and from information about the part of the
  //   grid owned by this processor

  ierr = nc.open_for_reading(filename);
  
  grid_info g;
  // Note that after this call g.z_len and g.zb_len are zero if the
  // corresponding dimension does not exist.
  ierr = nc.get_grid_info(g); CHKERRQ(ierr);

  double *zlevs = NULL, *zblevs = NULL; // NULLs correspond to 2D-only regridding
  if ((g.z_len != 0) && (g.zb_len != 0)) {
    ierr = nc.get_vertical_dims(zlevs, zblevs); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com,
		      "PISM WARNING: at least one of 'z' and 'zb' is absent in '%s'.\n"
		      "              3D regridding is disabled.\n",
		      filename);
    CHKERRQ(ierr);
  }
  ierr = nc.close(); CHKERRQ(ierr);

  LocalInterpCtx lic(g, zlevs, zblevs, grid); // will be de-allocated at 'return 0' below.

  set<string>::iterator i;
  for (i = vars.begin(); i != vars.end(); ++i) {
    IceModelVec *v = variables.get(*i);

    if (v == NULL) {
      ierr = PetscPrintf(grid.com, "PISM ERROR: unknown variable name: %s\n",
			 (*i).c_str()); CHKERRQ(ierr);
      PetscEnd();
    }

    string pism_intent = v->string_attr("pism_intent");
    if (pism_intent != "model_state") {
      ierr = verbPrintf(2, grid.com, "  WARNING: skipping '%s' (only model_state variables can be regridded)...\n",
			(*i).c_str()); CHKERRQ(ierr);
      continue;
    }

    if ( ((v->grid_type() == GRID_3D) && lic.regrid_2d_only) ||
	 ((v->grid_type() == GRID_3D_BEDROCK) && lic.no_regrid_bedrock) )
      {
      ierr = verbPrintf(2, grid.com, "  WARNING: skipping '%s'...\n",
			(*i).c_str()); CHKERRQ(ierr);
      continue;
    }

    ierr = v->regrid(filename, lic, true); CHKERRQ(ierr);

  }

  // Note that deleting a NULL pointer is safe.
  delete [] zlevs;  delete [] zblevs;
  return 0;
}

//! Initializes the snapshot-saving mechanism.
PetscErrorCode IceModel::init_snapshots() {
  PetscErrorCode ierr;
  PetscTruth save_at_set = PETSC_FALSE, save_to_set = PETSC_FALSE;
  char tmp[TEMPORARY_STRING_LENGTH] = "\0";
  current_snapshot = 0;

  ierr = check_old_option_and_stop(grid.com, "-save_to", "-save_file"); CHKERRQ(ierr);
  ierr = check_old_option_and_stop(grid.com, "-save_at", "-save_times"); CHKERRQ(ierr);

  ierr = PetscOptionsGetString(PETSC_NULL, "-save_file", tmp,
			       PETSC_MAX_PATH_LEN, &save_to_set); CHKERRQ(ierr);
  snapshots_filename = tmp;

  ierr = PetscOptionsGetString(PETSC_NULL, "-save_times", tmp,
			       TEMPORARY_STRING_LENGTH, &save_at_set); CHKERRQ(ierr);

  if (save_to_set ^ save_at_set) {
    ierr = PetscPrintf(grid.com,
		       "PISM ERROR: you need to specify both -save_file and -save_times to save snapshots.\n");
    CHKERRQ(ierr);
    PetscEnd();
  }

  if (!save_to_set && !save_at_set) {
    save_snapshots = false;
    return 0;
  }

  ierr = parse_times(grid.com, tmp, snapshot_times);
  if (ierr != 0) {
    ierr = PetscPrintf(grid.com, "PISM ERROR: parsing the -save_times argument failed.\n"); CHKERRQ(ierr);
    PetscEnd();
  }

  save_snapshots = true;
  snapshots_file_is_ready = false;
  split_snapshots = false;

  PetscTruth split;
  ierr = check_option("-split_snapshots", split); CHKERRQ(ierr);
  if (split) {
    split_snapshots = true;
  } else if (!ends_with(snapshots_filename, ".nc")) {
    ierr = verbPrintf(2, grid.com,
		      "PISM WARNING: snapshots file name does not have the '.nc' suffix!\n");
    CHKERRQ(ierr);
  }
  
  if (split) {
    ierr = verbPrintf(2, grid.com, "saving snapshots to '%s+year.nc'; ",
		      snapshots_filename.c_str()); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com, "saving snapshots to '%s'; ",
		      snapshots_filename.c_str()); CHKERRQ(ierr);
  }

  ierr = verbPrintf(2, grid.com, "times requested: %s\n", tmp); CHKERRQ(ierr);

  return 0;
}

//! Writes a snapshot of the model state (if necessary)
PetscErrorCode IceModel::write_snapshot() {
  PetscErrorCode ierr;
  NCTool nc(&grid);
  double saving_after = -1.0e30; // initialize to avoid compiler warning; this
				 // value is never used, because saving_after
				 // is only used if save_now == true, and in
				 // this case saving_after is guaranteed to be
				 // initialized. See the code below.
  char filename[PETSC_MAX_PATH_LEN];

  // determine if the user set the -save_times and -save_file options
  if (!save_snapshots)
    return 0;

  // do we need to save *now*?
  if ( (grid.year >= snapshot_times[current_snapshot]) && (current_snapshot < snapshot_times.size()) ) {
    saving_after = snapshot_times[current_snapshot];

    while (snapshot_times[current_snapshot] <= grid.year)
      current_snapshot++;
  } else {
    // we don't need to save now, so just return
    return 0;
  }

  if (split_snapshots) {
    snapshots_file_is_ready = false;	// each snapshot is written to a separate file
    snprintf(filename, PETSC_MAX_PATH_LEN, "%s-%06.0f.nc",
	     snapshots_filename.c_str(), grid.year);
  } else {
    strncpy(filename, snapshots_filename.c_str(), PETSC_MAX_PATH_LEN);
  }

  ierr = verbPrintf(2, grid.com, 
		    "\nsaving snapshot to %s at %.5f a, for time-step goal %.5f a\n\n",
		    filename, grid.year,saving_after);
  CHKERRQ(ierr);

  // create line for history in .nc file, including time of write

  string date_str = timestamp();
  char tmp[TEMPORARY_STRING_LENGTH];
  snprintf(tmp, TEMPORARY_STRING_LENGTH,
	   "%s: %s snapshot at %10.5f a, for time-step goal %10.5f a\n",
	   date_str.c_str(), executable_short_name.c_str(), grid.year, saving_after);

  if (!snapshots_file_is_ready) {

    // Prepare the snapshots file:
    ierr = nc.open_for_writing(filename, false, true); CHKERRQ(ierr);
    // append == false, check_dims == true
    ierr = nc.close(); CHKERRQ(ierr);

    ierr = global_attributes.write(filename); CHKERRQ(ierr);
    ierr = mapping.write(filename); CHKERRQ(ierr);
    snapshots_file_is_ready = true;
  }
    
  ierr = nc.open_for_writing(filename, true, true); CHKERRQ(ierr);
  // append == true, check_dims == true
  ierr = nc.append_time(grid.year); CHKERRQ(ierr);
  ierr = nc.write_history(tmp); CHKERRQ(ierr); // append the history
  ierr = nc.close(); CHKERRQ(ierr);

  ierr = write_model_state(filename);  CHKERRQ(ierr);

  // Let couplers write their fields:
  if (atmosPCC != PETSC_NULL) {
    ierr = atmosPCC->writeCouplingFieldsToFile(grid.year,filename); CHKERRQ(ierr);
  } else {
    SETERRQ(1,"PISM ERROR: atmosPCC == PETSC_NULL");
  }

  if (oceanPCC != PETSC_NULL) {
    ierr = oceanPCC->writeCouplingFieldsToFile(grid.year,filename); CHKERRQ(ierr);
  } else {
    SETERRQ(1,"PISM ERROR: oceanPCC == PETSC_NULL");
  }
    
  ierr = write_extra_fields(filename); CHKERRQ(ierr);

  return 0;
}
