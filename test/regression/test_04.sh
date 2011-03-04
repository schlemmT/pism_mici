#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

test="Test # 4: regridding during bootstrapping: coarse -> fine -> coarse."
files="foo.nc bar.nc baz.nc"

OPTS="-surface constant -Lz 4000 -y 0 -o_size small "
COARSE="-Mx 11 -My 11 -Mz 11"
FINE=" -Mx 21 -My 21 -Mz 21"

set -e

# Create a file to bootstrap from:
$PISM_PATH/pismv -test G $COARSE -y 0 -o foo.nc 

# Coarse -> fine:
$PISM_PATH/pismr -boot_file foo.nc $FINE   $OPTS -o bar.nc

# Fine -> coarse:
$PISM_PATH/pismr -boot_file bar.nc $COARSE $OPTS -o baz.nc

set +e

# Compare:
$PISM_PATH/nccmp.py -v topg foo.nc baz.nc
if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0
