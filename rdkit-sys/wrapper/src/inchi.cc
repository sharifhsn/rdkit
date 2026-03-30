#include "rust/cxx.h"
#include <GraphMol/GraphMol.h>
#include <GraphMol/inchi.h>

namespace RDKit {

rust::String mol_to_inchi(const std::shared_ptr<ROMol> &mol) {
	ExtraInchiReturnValues rv;
	return MolToInchi(*mol, rv);
}

rust::String mol_to_inchi_key(const std::shared_ptr<ROMol> &mol) {
	return MolToInchiKey(*mol);
}

} // namespace RDKit
