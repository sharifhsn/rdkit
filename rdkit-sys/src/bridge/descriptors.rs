#[cxx::bridge(namespace = "RDKit")]
pub mod ffi {
    unsafe extern "C++" {
        include!("wrapper/include/ro_mol.h");
        include!("wrapper/include/descriptors.h");

        pub type ROMol = crate::ro_mol_ffi::ROMol;
        pub type Properties;

        pub fn new_properties() -> SharedPtr<Properties>;
        pub fn get_property_names(
            properties: &SharedPtr<Properties>,
        ) -> UniquePtr<CxxVector<CxxString>>;
        pub fn compute_properties(
            properties: &SharedPtr<Properties>,
            mol: &SharedPtr<ROMol>,
        ) -> UniquePtr<CxxVector<f64>>;

        pub fn calc_amw(mol: &SharedPtr<ROMol>) -> f64;
        pub fn calc_clogp(mol: &SharedPtr<ROMol>) -> f64;
        pub fn calc_num_hbd(mol: &SharedPtr<ROMol>) -> u32;
        pub fn calc_num_hba(mol: &SharedPtr<ROMol>) -> u32;
        pub fn calc_tpsa(mol: &SharedPtr<ROMol>) -> f64;
        pub fn calc_num_rotatable_bonds(mol: &SharedPtr<ROMol>) -> u32;
    }
}
