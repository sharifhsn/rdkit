use std::collections::HashMap;

use cxx::SharedPtr;

use crate::ROMol;

/// Lipinski/Veber drug-likeness descriptors for a single molecule.
pub struct LipinskiDescriptors {
    pub mw: f64,
    pub clogp: f64,
    pub num_hbd: u32,
    pub num_hba: u32,
    pub tpsa: f64,
    pub num_rotatable_bonds: u32,
}

/// Calculate individual molecular descriptors without the overhead
/// of computing all 40+ descriptors via [`Properties::compute_properties`].
impl ROMol {
    pub fn calc_amw(&self) -> f64 {
        rdkit_sys::descriptors_ffi::calc_amw(&self.ptr)
    }

    pub fn calc_clogp(&self) -> f64 {
        rdkit_sys::descriptors_ffi::calc_clogp(&self.ptr)
    }

    pub fn calc_num_hbd(&self) -> u32 {
        rdkit_sys::descriptors_ffi::calc_num_hbd(&self.ptr)
    }

    pub fn calc_num_hba(&self) -> u32 {
        rdkit_sys::descriptors_ffi::calc_num_hba(&self.ptr)
    }

    pub fn calc_tpsa(&self) -> f64 {
        rdkit_sys::descriptors_ffi::calc_tpsa(&self.ptr)
    }

    pub fn calc_num_rotatable_bonds(&self) -> u32 {
        rdkit_sys::descriptors_ffi::calc_num_rotatable_bonds(&self.ptr)
    }

    pub fn lipinski_descriptors(&self) -> LipinskiDescriptors {
        LipinskiDescriptors {
            mw: self.calc_amw(),
            clogp: self.calc_clogp(),
            num_hbd: self.calc_num_hbd(),
            num_hba: self.calc_num_hba(),
            tpsa: self.calc_tpsa(),
            num_rotatable_bonds: self.calc_num_rotatable_bonds(),
        }
    }
}

pub struct Properties {
    ptr: SharedPtr<rdkit_sys::descriptors_ffi::Properties>,
}

impl Default for Properties {
    fn default() -> Self {
        Properties::new()
    }
}

impl Properties {
    pub fn new() -> Self {
        Properties {
            ptr: rdkit_sys::descriptors_ffi::new_properties(),
        }
    }

    pub fn compute_properties(&self, ro_mol: &ROMol) -> HashMap<String, f64> {
        let names = rdkit_sys::descriptors_ffi::get_property_names(&self.ptr);
        let computed = rdkit_sys::descriptors_ffi::compute_properties(&self.ptr, &ro_mol.ptr);

        assert!(!names.is_empty());
        assert!(computed.len() == names.len());

        names
            .into_iter()
            .zip(computed.as_slice())
            .map(|(k, v)| (k.to_string(), *v))
            .collect()
    }
}
