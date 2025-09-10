// Copyright (C) 2025 PISM Authors

//TO DO: implement selection based on configuration

#include "pism/rheology/melange/MelangeRheologyFactory.hh"

// Forward includes for concrete types once available
// #include "pism/rheology/melange/MelangeRheology.hh"
// #include "pism/rheology/melange/LocalRheology.hh"
// #include "pism/rheology/melange/GranularFluidity.hh"

namespace pism {
namespace rheology {
namespace melange {

class MelangeRheology; // forward declaration (remove when real header exists)

std::shared_ptr<MelangeRheology>
MelangeRheologyFactory::create(std::shared_ptr<const Grid> grid,
                               std::shared_ptr<const Config> config) {
  (void)grid; (void)config;
  //TO DO: inspect config (e.g., config->get_string("melange.rheology")) and construct the right class
  return nullptr;
}

} // namespace melange
} // namespace rheology
} // namespace pism


