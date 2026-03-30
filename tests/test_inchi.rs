use rdkit::ROMol;

#[test]
fn test_ethanol_inchi() {
    let mol = ROMol::from_smiles("CCO").unwrap();
    let inchi = mol.to_inchi();
    assert!(inchi.starts_with("InChI="), "got: {inchi}");
}

#[test]
fn test_ethanol_inchi_key() {
    let mol = ROMol::from_smiles("CCO").unwrap();
    let key = mol.to_inchi_key();
    assert_eq!(key, "LFQSCWFLJHTTHZ-UHFFFAOYSA-N");
}

#[test]
fn test_quercetin_inchi_key() {
    // Quercetin — common phytochemical in Dr. Duke's database
    let mol = ROMol::from_smiles("O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12").unwrap();
    let key = mol.to_inchi_key();
    assert_eq!(key, "REFJWTPEDVJJIY-UHFFFAOYSA-N");
}
