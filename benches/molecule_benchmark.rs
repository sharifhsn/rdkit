#![allow(soft_unstable)]
#![feature(test)]
extern crate test;

use rdkit::{has_substruct_match, ROMol, RWMol};

#[bench]
fn bench_molecules(bencher: &mut test::bench::Bencher) {
    let smiles1 = "c1ccccc1CCCCCCCC";

    bencher.iter(|| {
        ROMol::from_smiles(smiles1).unwrap();
    })
}

#[bench]
fn bench_fingerprint(bencher: &mut test::bench::Bencher) {
    let smiles1 = "c1ccccc1CCCCCCCC";
    let mol1 = ROMol::from_smiles(smiles1).unwrap();

    bencher.iter(|| mol1.rdk_fingerprint())
}

#[bench]
fn bench_tanimoto(bencher: &mut test::bench::Bencher) {
    let smiles1 = "c1ccccc1CCCCCCCC";
    let mol1 = ROMol::from_smiles(smiles1).unwrap();
    let smiles2 = "c1ccccc1CCCCCC";
    let mol2 = ROMol::from_smiles(smiles2).unwrap();

    bencher.iter(|| {
        let mol1_fingerprint = mol1.rdk_fingerprint();

        let mol2_fingerprint = mol2.rdk_fingerprint();

        mol1_fingerprint.tanimoto_distance(&mol2_fingerprint)
    });
}

#[bench]
fn bench_tanimoto_precomputed(bencher: &mut test::bench::Bencher) {
    let mol1 = ROMol::from_smiles("c1ccccc1CCCCCCCC").unwrap();
    let mol2 = ROMol::from_smiles("c1ccccc1CCCCCC").unwrap();
    let fp1 = mol1.rdk_fingerprint();
    let fp2 = mol2.rdk_fingerprint();

    bencher.iter(|| fp1.tanimoto_distance(&fp2));
}

#[bench]
fn bench_rwmol_smiles(bencher: &mut test::bench::Bencher) {
    let mol = ROMol::from_smiles("c1ccccc1CCCCCCCC").unwrap();
    let rw_mol = mol.as_rw_mol(false, -1);

    bencher.iter(|| rw_mol.as_smiles());
}

#[bench]
fn bench_lipinski_descriptors(bencher: &mut test::bench::Bencher) {
    let mol = ROMol::from_smiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O").unwrap();

    bencher.iter(|| mol.lipinski_descriptors());
}

#[bench]
fn bench_has_substruct_match(bencher: &mut test::bench::Bencher) {
    let mol = ROMol::from_smiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O").unwrap();
    let query = RWMol::from_smarts("[Br,Cl,I][CX4;CH,CH2]")
        .unwrap()
        .to_ro_mol();

    bencher.iter(|| has_substruct_match(&mol, &query));
}
