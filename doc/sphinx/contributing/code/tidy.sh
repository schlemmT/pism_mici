#!/bin/bash

set -e
set -u

# PISM's build directory
build_dir=~/local/build/pism

# Extract compiler options specifying MPI's include directories
mpi=$(mpicc -show | grep -E -o -e "-I[^ ]+")

# Run clang-tidy
clang-tidy --extra-arg=${mpi} -p ${build_dir} $@
