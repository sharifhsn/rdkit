use crate::ROMol;

impl ROMol {
    pub fn to_inchi(&self) -> String {
        rdkit_sys::inchi_ffi::mol_to_inchi(&self.ptr)
    }

    pub fn to_inchi_key(&self) -> String {
        rdkit_sys::inchi_ffi::mol_to_inchi_key(&self.ptr)
    }
}
