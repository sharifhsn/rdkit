#pragma once

#include "rust/cxx.h"
#include <GraphMol/GraphMol.h>

namespace RDKit {
rust::String mol_to_inchi(const std::shared_ptr<ROMol> &mol);
rust::String mol_to_inchi_key(const std::shared_ptr<ROMol> &mol);
} // namespace RDKit
