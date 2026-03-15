#include "rust/cxx.h"
#include <GraphMol/Descriptors/Crippen.h>
#include <GraphMol/Descriptors/Lipinski.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Descriptors/MolSurf.h>
#include <GraphMol/Descriptors/Property.h>
#include <GraphMol/ROMol.h>

namespace RDKit {
using Descriptors::Properties;

std::shared_ptr<Properties> new_properties() { return std::shared_ptr<Properties>(new Properties()); }

std::unique_ptr<std::vector<std::string>> get_property_names(const std::shared_ptr<Properties> &props) {
	std::vector<std::string> names       = props->getPropertyNames();
	std::vector<std::string> *names_heap = new std::vector<std::string>(names);
	return std::unique_ptr<std::vector<std::string>>(names_heap);
}

std::unique_ptr<std::vector<double>> compute_properties(const std::shared_ptr<Properties> &props,
                                                        const std::shared_ptr<ROMol> &mol) {
	std::vector<double> computed = props->computeProperties(*mol);
	auto computed_heap           = new std::vector<double>(computed);
	return std::unique_ptr<std::vector<double>>(computed_heap);
}

double calc_amw(const std::shared_ptr<ROMol> &mol) { return Descriptors::calcAMW(*mol); }

double calc_clogp(const std::shared_ptr<ROMol> &mol) {
	double logp, mr;
	Descriptors::calcCrippenDescriptors(*mol, logp, mr);
	return logp;
}

unsigned int calc_num_hbd(const std::shared_ptr<ROMol> &mol) { return Descriptors::calcNumHBD(*mol); }

unsigned int calc_num_hba(const std::shared_ptr<ROMol> &mol) { return Descriptors::calcNumHBA(*mol); }

double calc_tpsa(const std::shared_ptr<ROMol> &mol) { return Descriptors::calcTPSA(*mol); }

unsigned int calc_num_rotatable_bonds(const std::shared_ptr<ROMol> &mol) {
	return Descriptors::calcNumRotatableBonds(*mol);
}

} // namespace RDKit