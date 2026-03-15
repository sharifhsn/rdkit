//! Virtual Screening Filter Benchmark (Rust)
//!
//! Rust equivalent of benches/vs_filter_benchmark.py — the same pipeline
//! from Pat Walters' rd_filters (Relay Therapeutics).
//!
//! Pipeline per molecule:
//!   1. SMILES parsing
//!   2. Descriptor calculation (MW, LogP, HBD, HBA, TPSA, RotBonds)
//!   3. Substructure alert matching (SMARTS patterns)
//!   4. Property range filtering
//!
//! Usage:
//!   cargo run --release --example vs_filter_benchmark [SMILES_FILE]

use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

use rayon::prelude::*;
use rdkit::{ROMol, RWMol, has_substruct_match};

/// Real drugs from DrugBank (same as Python benchmark)
const BUILTIN_SMILES: &[&str] = &[
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "CC(=O)Oc1ccccc1C(=O)O",
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C",
    "CC(C)NCC(O)c1ccc(O)c(O)c1",
    "c1ccc2c(c1)cc(cc2)C(=O)O",
    "CC(=O)NC1=CC=C(C=C1)O",
    "OC(=O)c1ccccc1O",
    "c1ccc(cc1)C(=O)NC",
    "c1ccc2[nH]c(-c3ccccc3)nc2c1",
    "Clc1ccc(cc1)C(c2ccccc2)n3ccnc3",
    "OC(=O)CCc1ccccc1",
    "CC(C)(C)NCC(O)c1ccc(O)c(CO)c1",
    "c1ccc(c(c1)C(=O)O)NC2=CC=CC=C2",
    "c1ccnc(c1)C(=O)N",
    "CC(C)OC(=O)C(C)(C)Oc1ccc(Cl)cc1",
    "c1ccc2c(c1)c(=O)n(c(=O)[nH]2)C",
    "CC1=CC=C(C=C1)S(=O)(=O)NC(=O)NN2CCCCCC2",
    "COc1ccc(cc1OC)C(=O)CC(=O)c2ccc(OC)c(OC)c2",
    "c1ccc(cc1)c2cc(nn2c3ccccc3)C(F)(F)F",
];

/// Structural alerts (same subset as Python benchmark)
const ALERT_SMARTS: &[&str] = &[
    "[Br,Cl,I][CX4;CH,CH2]",
    "[S,C](=[O,S])[F,Br,Cl,I]",
    "C(=O)OC(=O)",
    "OO",
    "N=C=[S,O]",
    "OS(=O)(=O)C(F)(F)F",
    "[N;R0][N;R0]C(=O)",
    "[S-]",
    "[SX2]C(=O)",
    "C#C",
    "[N;R0]=[N;R0]",
    "SC#N",
    "C(=O)Oc1c(F)c(F)c(F)c(F)c1(F)",
    "N=C=N",
    "[NH1,NH2][CX4][Cl,Br,I]",
    "C1(=O)OCC1",
    "[CH1](=O)",
    "[O,S][CH2]O",
    "C(=O)N(C=O)",
    "[NX3H][NX3H]",
    "c1ccc2c(c1)ccc(=O)o2",
    "c1ccc(c(c1)O)N=Nc2ccccc2",
    "[#6]S(=O)(=O)c1ccc(N)cc1",
    "c1cc(ccc1N=Nc2ccccc2)O",
    "c1ccc2c(c1)[nH]c(=O)[nH]2",
    "c1ccc2c(c1)c(=O)c3ccccc3o2",
    "O=C1C=CC(=O)C=C1",
    "c1cc2c(cc1)oc(=O)c(c2=O)",
    "c1ccc(cc1)S(=O)(=O)N",
    "c1cc(oc1)C=NNC(=O)",
];

fn evaluate_molecule(smi: &str, alert_patterns: &[ROMol]) -> Option<&'static str> {
    let mol = ROMol::from_smiles(smi).ok()?;

    let desc = mol.lipinski_descriptors();

    // Structural alert check
    for pat in alert_patterns {
        if has_substruct_match(&mol, pat) {
            return Some("ALERT");
        }
    }

    // Property range filter (Lipinski + Veber)
    if desc.mw > 500.0
        || desc.clogp > 5.0
        || desc.clogp < -5.0
        || desc.num_hbd > 5
        || desc.num_hba > 10
        || desc.tpsa > 200.0
        || desc.num_rotatable_bonds > 10
    {
        return Some("PROPERTY_FAIL");
    }

    Some("PASS")
}

fn generate_builtin_dataset(n: usize) -> Vec<String> {
    (0..n)
        .map(|i| BUILTIN_SMILES[i % BUILTIN_SMILES.len()].to_string())
        .collect()
}

fn load_smiles_file(path: &str) -> Vec<String> {
    std::fs::read_to_string(path)
        .expect("could not read file")
        .lines()
        .filter_map(|line| line.split_whitespace().next().map(String::from))
        .collect()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let smiles_list = if args.len() > 1 {
        load_smiles_file(&args[1])
    } else {
        eprintln!("No input file provided, using built-in drug-like dataset (50K molecules)");
        eprintln!();
        generate_builtin_dataset(50_000)
    };

    let n = smiles_list.len();

    // Pre-compile SMARTS patterns
    let alert_patterns: Vec<ROMol> = ALERT_SMARTS
        .iter()
        .filter_map(|sma| RWMol::from_smarts(sma).ok().map(|rw| rw.to_ro_mol()))
        .collect();

    eprintln!("Molecules: {n}");
    eprintln!("Alert patterns: {}", alert_patterns.len());
    eprintln!("Rayon threads: {}", rayon::current_num_threads());
    eprintln!();

    // --- Single-threaded pipeline ---
    let t0 = Instant::now();
    let mut st_pass = 0u64;
    let mut st_alert = 0u64;
    let mut st_prop_fail = 0u64;
    let mut st_invalid = 0u64;

    for smi in &smiles_list {
        match evaluate_molecule(smi, &alert_patterns) {
            Some("PASS") => st_pass += 1,
            Some("ALERT") => st_alert += 1,
            Some("PROPERTY_FAIL") => st_prop_fail += 1,
            Some(_) => {}
            None => st_invalid += 1,
        }
    }
    let t_single = t0.elapsed();

    // --- Parallel pipeline (Rayon) ---
    let t0 = Instant::now();
    let par_pass = AtomicU64::new(0);
    let par_alert = AtomicU64::new(0);
    let par_prop_fail = AtomicU64::new(0);
    let par_invalid = AtomicU64::new(0);

    smiles_list.par_iter().for_each(|smi| {
        match evaluate_molecule(smi, &alert_patterns) {
            Some("PASS") => {
                par_pass.fetch_add(1, Ordering::Relaxed);
            }
            Some("ALERT") => {
                par_alert.fetch_add(1, Ordering::Relaxed);
            }
            Some("PROPERTY_FAIL") => {
                par_prop_fail.fetch_add(1, Ordering::Relaxed);
            }
            Some(_) => {}
            None => {
                par_invalid.fetch_add(1, Ordering::Relaxed);
            }
        }
    });
    let t_parallel = t0.elapsed();

    // Verify results match
    assert_eq!(st_pass, par_pass.load(Ordering::Relaxed));
    assert_eq!(st_alert, par_alert.load(Ordering::Relaxed));
    assert_eq!(st_prop_fail, par_prop_fail.load(Ordering::Relaxed));

    println!(
        "{:<30} {:>10} {:>10} {:>12}",
        "Phase", "Time (s)", "us/mol", "mol/s"
    );
    println!("{}", "-".repeat(65));
    println!(
        "{:<30} {:>10.3} {:>10.1} {:>12.0}",
        "Single-threaded pipeline",
        t_single.as_secs_f64(),
        t_single.as_secs_f64() / n as f64 * 1e6,
        n as f64 / t_single.as_secs_f64()
    );
    println!(
        "{:<30} {:>10.3} {:>10.1} {:>12.0}",
        "Parallel pipeline (Rayon)",
        t_parallel.as_secs_f64(),
        t_parallel.as_secs_f64() / n as f64 * 1e6,
        n as f64 / t_parallel.as_secs_f64()
    );
    println!(
        "{:<30} {:>10.1}x",
        "Speedup",
        t_single.as_secs_f64() / t_parallel.as_secs_f64()
    );
    eprintln!();
    eprintln!(
        "Results: {} pass / {} alerts / {} property fail / {} invalid",
        st_pass, st_alert, st_prop_fail, st_invalid
    );
}
