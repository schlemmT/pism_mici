/* VariableMetadata and SpatialVariableMetadata don't have default constructors. */
/* This should go before array::Array so that array::Array::metadata()
   is wrapped properly. */
%feature("valuewrapper") pism::VariableMetadata;
%feature("valuewrapper") pism::SpatialVariableMetadata;

%ignore pism::Attribute::operator=;
%ignore pism::VariableMetadata::operator[];

%include "util/VariableMetadata.hh"
