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

use std::time::Instant;

use rdkit::{ROMol, RWMol, SubstructMatchParameters, substruct_match};

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
        .filter_map(|sma| {
            RWMol::from_smarts(sma)
                .ok()
                .map(|rw| rw.to_ro_mol())
        })
        .collect();

    let match_params = SubstructMatchParameters::default();

    eprintln!("Molecules: {n}");
    eprintln!("Alert patterns: {}", alert_patterns.len());
    eprintln!();

    // Phase 1: SMILES parsing only
    let t0 = Instant::now();
    let mols: Vec<Option<ROMol>> = smiles_list
        .iter()
        .map(|smi| ROMol::from_smiles(smi).ok())
        .collect();
    let t_parse = t0.elapsed();
    let n_invalid = mols.iter().filter(|m| m.is_none()).count();

    // Phase 2: Descriptor calculation only
    let t0 = Instant::now();
    for mol in &mols {
        if let Some(mol) = mol {
            let _ = mol.lipinski_descriptors();
        }
    }
    let t_desc = t0.elapsed();

    // Phase 3: Substructure alert matching only
    let t0 = Instant::now();
    for mol in &mols {
        if let Some(mol) = mol {
            for pat in &alert_patterns {
                let _ = substruct_match(mol, pat, &match_params);
            }
        }
    }
    let t_alerts = t0.elapsed();

    // Phase 4: Full pipeline (end-to-end)
    let t0 = Instant::now();
    let mut n_pass = 0u64;
    let mut n_alert = 0u64;
    let mut n_prop_fail = 0u64;

    for smi in &smiles_list {
        let mol = match ROMol::from_smiles(smi) {
            Ok(m) => m,
            Err(_) => continue,
        };

        let desc = mol.lipinski_descriptors();

        // Structural alert check
        let mut alerted = false;
        for pat in &alert_patterns {
            if !substruct_match(&mol, pat, &match_params).is_empty() {
                alerted = true;
                break;
            }
        }

        if alerted {
            n_alert += 1;
        } else if desc.mw > 500.0
            || desc.clogp > 5.0
            || desc.clogp < -5.0
            || (desc.num_hbd as f64) > 5.0
            || (desc.num_hba as f64) > 10.0
            || desc.tpsa > 200.0
            || (desc.num_rotatable_bonds as f64) > 10.0
        {
            n_prop_fail += 1;
        } else {
            n_pass += 1;
        }
    }
    let t_total = t0.elapsed();

    let n_valid = (n - n_invalid) as f64;

    println!(
        "{:<30} {:>10} {:>10} {:>12}",
        "Phase", "Time (s)", "us/mol", "mol/s"
    );
    println!("{}", "-".repeat(65));
    println!(
        "{:<30} {:>10.3} {:>10.1} {:>12.0}",
        "SMILES parsing",
        t_parse.as_secs_f64(),
        t_parse.as_secs_f64() / n as f64 * 1e6,
        n as f64 / t_parse.as_secs_f64()
    );
    println!(
        "{:<30} {:>10.3} {:>10.1} {:>12.0}",
        "Descriptor calculation",
        t_desc.as_secs_f64(),
        t_desc.as_secs_f64() / n_valid * 1e6,
        n_valid / t_desc.as_secs_f64()
    );
    println!(
        "{:<30} {:>10.3} {:>10.1} {:>12.0}",
        "Substructure alerts",
        t_alerts.as_secs_f64(),
        t_alerts.as_secs_f64() / n_valid * 1e6,
        n_valid / t_alerts.as_secs_f64()
    );
    println!(
        "{:<30} {:>10.3} {:>10.1} {:>12.0}",
        "Full pipeline (end-to-end)",
        t_total.as_secs_f64(),
        t_total.as_secs_f64() / n as f64 * 1e6,
        n as f64 / t_total.as_secs_f64()
    );
    eprintln!();
    eprintln!(
        "Results: {} pass / {} alerts / {} property fail / {} invalid",
        n_pass, n_alert, n_prop_fail, n_invalid
    );
}
