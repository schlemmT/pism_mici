// Copyright (C) 2025 PISM Authors
//
// Factory for selecting a melange rheology implementation.

#ifndef PISM_RHEOLOGY_MELANGE_FACTORY_HH
#define PISM_RHEOLOGY_MELANGE_FACTORY_HH

#include <memory>

namespace pism {
class Grid; // forward declaration
class Config; // forward declaration

namespace rheology {
namespace melange {

class MelangeRheology; // forward declaration of the common interface

class MelangeRheologyFactory {
public:
  // Create a melange rheology instance based on configuration (e.g., melange.rheology)
  static std::shared_ptr<MelangeRheology>
  create(std::shared_ptr<const Grid> grid,
         std::shared_ptr<const Config> config);
};

} // namespace melange
} // namespace rheology
} // namespace pism

#endif // PISM_RHEOLOGY_MELANGE_FACTORY_HH


