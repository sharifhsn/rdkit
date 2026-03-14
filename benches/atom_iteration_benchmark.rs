#![allow(soft_unstable)]
#![feature(test)]
extern crate test;

use rdkit::ROMol;

/// Drug-like molecules of varying size for realistic benchmarking.
/// These cover common pharmaceutical scaffolds and natural products.
const SMILES_SET: &[&str] = &[
    // aspirin (21 atoms with H)
    "CC(=O)Oc1ccccc1C(=O)O",
    // ibuprofen
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    // caffeine
    "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
    // diazepam
    "O=C1CN=C(c2ccccc2)c2cc(Cl)ccc2N1C",
    // atorvastatin (lipitor) — large drug molecule
    "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O",
    // taxol core — complex natural product
    "CC1=C2C(OC(=O)c3ccccc3)C(O)C4(OC(=O)C(O)(CC(OC(=O)c5ccccc5)C1O)C24C)C(=O)c1ccc(OC)cc1",
    // vancomycin fragment — large, many atoms
    "OC1C(O)C(OC2C(O)C(O)C(O)C(CO)O2)OC(CO)C1NC(=O)C1CC(O)CN1C(=O)C(NC(=O)C1CC(=O)NC(=O)C1O)C(O)c1ccc(O)cc1",
];

/// Benchmark: parse SMILES into molecules.
/// Baseline cost — everything else is on top of this.
#[bench]
fn bench_parse_smiles(b: &mut test::bench::Bencher) {
    b.iter(|| {
        for smiles in SMILES_SET {
            test::black_box(ROMol::from_smiles(smiles).unwrap());
        }
    });
}

/// Benchmark: iterate all atoms and read one cheap property (atomic number).
/// Each iteration: N atoms * 2 FFI calls (get_atom + get_atomic_num).
#[bench]
fn bench_atom_one_property(b: &mut test::bench::Bencher) {
    let mut mols: Vec<ROMol> = SMILES_SET
        .iter()
        .map(|s| ROMol::from_smiles(s).unwrap())
        .collect();

    b.iter(|| {
        for mol in &mut mols {
            let n = mol.num_atoms(true);
            for i in 0..n {
                let atom = mol.atom_with_idx(i);
                test::black_box(atom.get_atomic_num());
            }
        }
    });
}

/// Benchmark: iterate all atoms, read 7 properties per atom.
/// This is the realistic "featurization" workload — the hot path
/// in ML pipelines, QSAR descriptor computation, etc.
/// Each iteration: N atoms * 8 FFI calls (1 get_atom + 7 properties).
#[bench]
fn bench_atom_all_properties(b: &mut test::bench::Bencher) {
    let mut mols: Vec<ROMol> = SMILES_SET
        .iter()
        .map(|s| ROMol::from_smiles(s).unwrap())
        .collect();

    b.iter(|| {
        for mol in &mut mols {
            let n = mol.num_atoms(true);
            for i in 0..n {
                let atom = mol.atom_with_idx(i);
                test::black_box(atom.symbol());
                test::black_box(atom.get_atomic_num());
                test::black_box(atom.get_formal_charge());
                test::black_box(atom.get_is_aromatic());
                test::black_box(atom.get_hybridization_type());
                test::black_box(atom.get_degree());
                test::black_box(atom.get_total_num_hs());
            }
        }
    });
}

/// Benchmark: full pipeline — parse + featurize.
/// Measures end-to-end cost of the most common workflow.
#[bench]
fn bench_parse_and_featurize(b: &mut test::bench::Bencher) {
    b.iter(|| {
        for smiles in SMILES_SET {
            let mut mol = ROMol::from_smiles(smiles).unwrap();
            let n = mol.num_atoms(true);
            for i in 0..n {
                let atom = mol.atom_with_idx(i);
                test::black_box(atom.symbol());
                test::black_box(atom.get_atomic_num());
                test::black_box(atom.get_formal_charge());
                test::black_box(atom.get_is_aromatic());
                test::black_box(atom.get_hybridization_type());
                test::black_box(atom.get_degree());
                test::black_box(atom.get_total_num_hs());
            }
        }
    });
}

/// Benchmark: atom_ref with one property (const path, no &mut needed).
#[bench]
fn bench_atom_ref_one_property(b: &mut test::bench::Bencher) {
    let mols: Vec<ROMol> = SMILES_SET
        .iter()
        .map(|s| ROMol::from_smiles(s).unwrap())
        .collect();

    b.iter(|| {
        for mol in &mols {
            let n = mol.num_atoms(true);
            for i in 0..n {
                let atom = mol.atom_ref(i);
                test::black_box(atom.get_atomic_num());
            }
        }
    });
}

/// Benchmark: atom_ref with all 7 properties (const path).
/// Compare directly against bench_atom_all_properties.
#[bench]
fn bench_atom_ref_all_properties(b: &mut test::bench::Bencher) {
    let mols: Vec<ROMol> = SMILES_SET
        .iter()
        .map(|s| ROMol::from_smiles(s).unwrap())
        .collect();

    b.iter(|| {
        for mol in &mols {
            let n = mol.num_atoms(true);
            for i in 0..n {
                let atom = mol.atom_ref(i);
                test::black_box(atom.symbol());
                test::black_box(atom.get_atomic_num());
                test::black_box(atom.get_formal_charge());
                test::black_box(atom.get_is_aromatic());
                test::black_box(atom.get_hybridization_type());
                test::black_box(atom.get_degree());
                test::black_box(atom.get_total_num_hs());
            }
        }
    });
}

/// Benchmark: clone + featurize via atom_with_idx (the old &ROMol workflow).
/// This is the cost users pay when they only have &ROMol.
#[bench]
fn bench_clone_and_featurize(b: &mut test::bench::Bencher) {
    let mols: Vec<ROMol> = SMILES_SET
        .iter()
        .map(|s| ROMol::from_smiles(s).unwrap())
        .collect();

    b.iter(|| {
        for mol in &mols {
            let mut mol = mol.clone();
            let n = mol.num_atoms(true);
            for i in 0..n {
                let atom = mol.atom_with_idx(i);
                test::black_box(atom.symbol());
                test::black_box(atom.get_atomic_num());
                test::black_box(atom.get_formal_charge());
                test::black_box(atom.get_is_aromatic());
                test::black_box(atom.get_hybridization_type());
                test::black_box(atom.get_degree());
                test::black_box(atom.get_total_num_hs());
            }
        }
    });
}

/// Benchmark: clone cost alone.
/// atom_with_idx requires &mut self, so callers who have &ROMol must clone.
/// This measures that tax.
#[bench]
fn bench_clone_molecules(b: &mut test::bench::Bencher) {
    let mols: Vec<ROMol> = SMILES_SET
        .iter()
        .map(|s| ROMol::from_smiles(s).unwrap())
        .collect();

    b.iter(|| {
        for mol in &mols {
            test::black_box(mol.clone());
        }
    });
}
