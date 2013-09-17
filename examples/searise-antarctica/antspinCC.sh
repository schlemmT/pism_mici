#!/bin/bash

# Copyright (C) 2009-2013  Maria Martin and Ed Bueler
##################################################################################
# Spinup of Antarctic ice sheet model using data from Anne Le Brocq (from SeaRISE wiki).
# Uses PIK physics and enthalpy model (see Publications at www.pism-docs.org) 
# and modified configuration parameters with constant climate.  Uses constant
# and precip and a parameterization for artm as in Martin et al (2012).
# WARNING: especially for the last part, with finer resolutions, output is large!
##################################################################################


SCRIPTNAME="#(antspinCC.sh)"


echo "$SCRIPTNAME   run preprocess.sh before this..."
#inputs generated by preprocess.sh:
#  present day           = pism_Antarctica_5km.nc 					
#  glacial cycle forcing = pism_dT.nc pism_dSL.nc
#  future climate data   = ar4_ant_acab_anomaly_scalefactor_{1.0,1.5,2.0}.nc
#                          ar4_ant_artm_anomaly_scalefactor_{1.0,1.5,2.0}.nc

set -e  # exit on error

echo "$SCRIPTNAME   Constant-climate spinup-script using SeaRISE-Antarctica data"
echo "$SCRIPTNAME      and -ssa_sliding and -pik"
echo "$SCRIPTNAME   Run as './antspinCC.sh NN' for NN procs and 30km grid"


# naming files, directories, executables
RESDIR=
BOOTDIR=
PISM_EXEC=pismr
PISM_MPIDO="mpiexec -n "

# input data:
PISM_INDATANAME=${BOOTDIR}pism_Antarctica_5km.nc
PISM_TEMPSERIES=${BOOTDIR}pism_dT.nc
PISM_SLSERIES=${BOOTDIR}pism_dSL.nc

NN=4  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "antspinup.sh 8" then NN = 8
  NN="$1"
fi
echo "$SCRIPTNAME              NN = $NN"
set -e  # exit on error

# check if env var PISM_DO was set (i.e. PISM_DO=echo for a 'dry' run)
if [ -n "${PISM_DO:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME         PISM_DO = $PISM_DO  (already set)"
else
  PISM_DO="" 
fi
DO=$PISM_DO

# grids
THIRTYKMGRID="-Mx 200 -My 200 -Lz 5000 -Lbz 2000 -Mz 41 -Mbz 16"
TWENTYKMGRID="-Mx 300 -My 300 -Lz 5000 -Lbz 2000 -Mz 81 -Mbz 21"
FIFTEENKMGRID="-Mx 400 -My 400 -Lz 5000 -Lbz 2000 -Mz 81 -Mbz 21"
TWELVEKMGRID="-Mx 500 -My 500 -Lz 5000 -Lbz 2000 -Mz 101 -Mbz 31"
TENKMGRID="-Mx 600 -My 600 -Lz 5000 -Lbz 2000 -Mz 101 -Mbz 31"
SEVENKMGRID="-Mx 900 -My 900 -Lz 5000 -Lbz 2000 -Mz 151 -Mbz 31"
FIVEKMGRID="-Mx 1200 -My 1200 -Lz 5000 -Lbz 2000 -Mz 201 -Mbz 51"

# skips:  
SKIPTHIRTYKM=10
SKIPTWENTYKM=10
SKIPFIFTEENKM=10
SKIPTWELVEKM=50
SKIPTENKM=100
SKIPSEVENKM=100
SKIPFIVEKM=200

# these coarse grid defaults are for development/regression runs, not
# "production" or science
GRID=$THIRTYKMGRID
SKIP=$SKIPTHIRTYKM
GRIDNAME=30km

SIA_ENHANCEMENT="-sia_e 5.6"

#PIK-stuff; notes:
# 1)   '-pik' = '-cfbc -part_grid -part_redist -kill_icebergs'
# 2)   -meltfactor_pik 5e-3 is default when using -ocean pik
PIKPHYS="-ssa_method fd -ssa_e 0.6 -pik -eigen_calving -eigen_calving_K 2.0e18 -thickness_calving -calving_at_thickness 50.0"
PIKPHYS_COUPLING="-atmosphere given -atmosphere_given_file $PISM_INDATANAME -surface simple -ocean pik -meltfactor_pik 1.5e-2"

# sliding related options:
PARAMS="-pseudo_plastic -pseudo_plastic_q 0.25 -till_effective_fraction_overburden 0.01"
#TILLPHI="-topg_to_phi 5.0,20.0,-1000.0,0.0,10.0" # as in Martin et al 2012
TILLPHI="-topg_to_phi 15.0,40.0,-300.0,700.0"
FULLPHYS="-ssa_sliding -hydrology null $PARAMS $TILLPHI"

# use these if KSP "diverged" errors occur
STRONGKSP="-ksp_type gmres -ksp_norm_type unpreconditioned -ksp_pc_side right -pc_type asm -sub_pc_type lu"


echo "$SCRIPTNAME             PISM = $PISM_EXEC"
echo "$SCRIPTNAME         FULLPHYS = $FULLPHYS"
echo "$SCRIPTNAME          PIKPHYS = $PIKPHYS"
echo "$SCRIPTNAME PIKPHYS_COUPLING = $PIKPHYS_COUPLING"


# #######################################
# bootstrap and SHORT smoothing run to 100 years
# #######################################
stage=earlyone
INNAME=$PISM_INDATANAME
RESNAMEONE=${RESDIR}${stage}_${GRIDNAME}.nc
RUNTIME=1
echo
echo "$SCRIPTNAME  bootstrapping plus SIA run for $RUNTIME a"
cmd="$PISM_MPIDO $NN $PISM_EXEC -skip -skip_max $SKIP -boot_file ${INNAME} $GRID \
	$SIA_ENHANCEMENT $PIKPHYS_COUPLING -ocean_kill \
	-y $RUNTIME -o $RESNAMEONE"
$DO $cmd
#exit # <-- uncomment to stop here

# replace topg with topgsmooth (does nothing if $RESNAMEONE is absent)
echo
echo "$SCRIPTNAME  replacing topg in $RESNAMEONE with topgsmooth"
cmd="ncatted -O -a standard_name,topg,d,, $RESNAMEONE"
$DO $cmd
cmd="ncrename -O -v topg,topgoriginal $RESNAMEONE"
$DO $cmd
cmd="ncrename -O -v topgsmooth,topg $RESNAMEONE"
$DO $cmd
cmd="ncatted -O -a standard_name,topg,o,c,bedrock_altitude $RESNAMEONE"
$DO $cmd
cmd="ncatted -O -a long_name,topg,d,, $RESNAMEONE"
$DO $cmd
cmd="ncatted -O -a pism_intent,topg,d,, $RESNAMEONE"
$DO $cmd


stage=smoothing
RESNAME=${RESDIR}${stage}_${GRIDNAME}.nc
RUNTIME=100
echo
echo "$SCRIPTNAME  short SIA run for $RUNTIME a"
cmd="$PISM_MPIDO $NN $PISM_EXEC -skip -skip_max $SKIP -i $RESNAMEONE \
	$SIA_ENHANCEMENT $PIKPHYS_COUPLING -ocean_kill \
	-y $RUNTIME -o $RESNAME"
$DO $cmd

# #######################################
# run with -no_mass (no surface change) on 15km for 200ka
# #######################################
stage=nomass
INNAME=$RESNAME
RESNAME=${RESDIR}${stage}_${GRIDNAME}.nc
TSNAME=${RESDIR}ts_${stage}_${GRIDNAME}.nc
RUNTIME=200000 
EXTRANAME=${RESDIR}extra_${stage}_${GRIDNAME}.nc
expackage="-extra_times 0:1000:$RUNTIME -extra_vars bmelt,tillwat,csurf,temppabase,diffusivity,hardav"
echo
echo "$SCRIPTNAME  -no_mass (no surface change) SIA for $RUNTIME a"
cmd="$PISM_MPIDO $NN $PISM_EXEC -i $INNAME $PIKPHYS_COUPLING  \
	$SIA_ENHANCEMENT -no_mass \
  	-ys 0 -y $RUNTIME \
	-extra_file $EXTRANAME $expackage \
  	-o $RESNAME"
$DO $cmd
#exit # <-- uncomment to stop here


# #######################################
# run into steady state with constant climate forcing
# #######################################
stage=run
INNAME=$RESNAME
RESNAME=${RESDIR}${stage}_${GRIDNAME}.nc
TSNAME=${RESDIR}ts_${stage}_${GRIDNAME}.nc
RUNTIME=100000 
EXTRANAME=${RESDIR}extra_${stage}_${GRIDNAME}.nc
exvars="thk,usurf,cbase,cbar,mask,diffusivity,tauc,bmelt,tillwat,temppabase,hardav"
expackage="-extra_times 0:1000:$RUNTIME -extra_vars $exvars"

echo
echo "$SCRIPTNAME  run into steady state with constant climate forcing for $RUNTIME a"
cmd="$PISM_MPIDO $NN $PISM_EXEC -skip -skip_max $SKIP -i $INNAME \
	$SIA_ENHANCEMENT $PIKPHYS_COUPLING $PIKPHYS $FULLPHYS \
	-ys 0 -y $RUNTIME \
	-ts_file $TSNAME -ts_times 0:1:$RUNTIME \
	-extra_file $EXTRANAME $expackage \
	$STRONGKSP \
	-o $RESNAME -o_size big"
$DO $cmd


exit # <-- uncomment to continue

# #######################################
## do a regridding to 10km and reset year to 0 (demonstrate for only 100 model years):
# #######################################
stage=run_regrid_10km
INNAME=$RESNAME
RESNAME=${RESDIR}$stage.nc
TSNAME=${RESDIR}ts_$stage.nc
RUNTIME=100
EXTRANAME=${RESDIR}extra_$stage.nc
expackage="-extra_times 0:10:$RUNTIME -extra_vars $exvars"

echo
echo "$SCRIPTNAME  continue but regrid to 10km"
cmd="$PISM_MPIDO $NN $PISM_EXEC -skip -skip_max $SKIPTENKM \
    -boot_file $PISM_INDATANAME $TENKMGRID \
    -regrid_file $INNAME -regrid_vars litho_temp,thk,enthalpy,tillwat,bmelt \
    $SIA_ENHANCEMENT $PIKPHYS_COUPLING $PIKPHYS $FULLPHYS \
    -ys 0 -y $RUNTIME \
    -ts_file $TSNAME -ts_times 0:1:$RUNTIME \
    -extra_file $EXTRANAME $expackage \
    $STRONGKSP \
    -o $RESNAME -o_size big"

$DO $cmd
#exit # <-- uncomment to stop here

# one can do a regridding to 6.7km, 5km and so on, if the number of model years
# (RUNTIME) is appropriately shortened, and sufficient memory is available

echo
echo "$SCRIPTNAME  spinup done"

