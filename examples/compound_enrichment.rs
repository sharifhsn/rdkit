//! Compound Enrichment Tool
//!
//! Pre-computes InChIKeys, molecular formulas, and descriptors locally,
//! eliminating PubChem API calls for ChEMBL enrichment pipelines.
//!
//! Usage:
//!   cargo run --release --example compound_enrichment [TSV_FILE]
//!
//! Input TSV: compound_name<TAB>SMILES (one per line, no header)
//! Output CSV to stdout: compound_name,canonical_smiles,inchi_key,formula,mw,status

use std::collections::HashSet;
use std::time::Instant;

use rayon::prelude::*;
use rdkit::ROMol;

/// Built-in test compounds when no file is provided
const BUILTIN_COMPOUNDS: &[(&str, &str)] = &[
    ("QUERCETIN", "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"),
    ("KAEMPFEROL", "O=c1c(O)c(-c2ccc(O)cc2)oc2cc(O)cc(O)c12"),
    ("LUTEOLIN", "O=c1cc(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"),
    ("APIGENIN", "O=c1cc(-c2ccc(O)cc2)oc2cc(O)cc(O)c12"),
    ("CURCUMIN", "O=C(/C=C/c1ccc(O)c(OC)c1)CC(=O)/C=C/c1ccc(O)c(OC)c1"),
    ("RESVERATROL", "Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1"),
    ("CATECHIN", "Oc1cc(O)c2c(c1)O[C@@H](c1ccc(O)c(O)c1)[C@H](O)C2"),
    ("EPICATECHIN", "Oc1cc(O)c2c(c1)O[C@@H](c1ccc(O)c(O)c1)[C@@H](O)C2"),
    ("NARINGENIN", "O=C1C[C@H](c2ccc(O)cc2)Oc2cc(O)cc(O)c21"),
    ("GENISTEIN", "O=c1c(-c2ccc(O)cc2)coc2cc(O)cc(O)c12"),
    // Noise compounds — should be flagged
    ("WATER", "O"),
    ("GLUCOSE", "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"),
    ("ETHANOL", "CCO"),
    ("ACETIC ACID", "CC(O)=O"),
    ("METHANE", "C"),
];

/// Macronutrients and common solvents to flag as noise
const NOISE_NAMES: &[&str] = &[
    "WATER",
    "GLUCOSE",
    "FRUCTOSE",
    "SUCROSE",
    "GALACTOSE",
    "MALTOSE",
    "LACTOSE",
    "STARCH",
    "CELLULOSE",
    "ETHANOL",
    "METHANOL",
    "ACETIC ACID",
    "CITRIC ACID",
    "LACTIC ACID",
    "CARBON DIOXIDE",
    "OXYGEN",
    "NITROGEN",
    "HYDROGEN",
    "AMMONIA",
    "UREA",
    "GLYCEROL",
    "PALMITIC ACID",
    "STEARIC ACID",
    "OLEIC ACID",
    "LINOLEIC ACID",
    "GLUTAMIC ACID",
    "ASPARTIC ACID",
    "ALANINE",
    "GLYCINE",
];

const MIN_MW: f64 = 50.0;

struct EnrichedRow {
    name: String,
    canonical_smiles: String,
    inchi_key: String,
    formula: String,
    mw: f64,
    status: String,
}

fn enrich_compound(name: &str, smiles: &str, noise_set: &HashSet<&str>) -> Option<EnrichedRow> {
    let mol = ROMol::from_smiles(smiles).ok()?;

    let canonical_smiles = mol.as_smiles();
    let inchi_key = mol.to_inchi_key();
    let formula = mol.calc_mol_formula();
    let mw = mol.calc_amw();

    let status = if noise_set.contains(name.to_uppercase().as_str()) {
        "NOISE:macronutrient".to_string()
    } else if mw < MIN_MW {
        "NOISE:low_mw".to_string()
    } else {
        "OK".to_string()
    };

    Some(EnrichedRow {
        name: name.to_string(),
        canonical_smiles,
        inchi_key,
        formula,
        mw,
        status,
    })
}

fn load_tsv(path: &str) -> Vec<(String, String)> {
    std::fs::read_to_string(path)
        .expect("could not read file")
        .lines()
        .filter_map(|line| {
            let mut parts = line.splitn(2, '\t');
            let name = parts.next()?.trim().to_string();
            let smiles = parts.next()?.trim().to_string();
            if name.is_empty() || smiles.is_empty() {
                return None;
            }
            Some((name, smiles))
        })
        .collect()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let compounds: Vec<(String, String)> = if args.len() > 1 {
        load_tsv(&args[1])
    } else {
        eprintln!("No input file provided, using built-in phytochemical dataset");
        eprintln!();
        BUILTIN_COMPOUNDS
            .iter()
            .map(|(n, s)| (n.to_string(), s.to_string()))
            .collect()
    };

    let noise_set: HashSet<&str> = NOISE_NAMES.iter().copied().collect();
    let n = compounds.len();

    eprintln!("Compounds: {n}");
    eprintln!("Rayon threads: {}", rayon::current_num_threads());

    let t0 = Instant::now();
    let rows: Vec<EnrichedRow> = compounds
        .par_iter()
        .filter_map(|(name, smiles)| enrich_compound(name, smiles, &noise_set))
        .collect();
    let elapsed = t0.elapsed();

    // CSV header
    println!("compound_name,canonical_smiles,inchi_key,formula,mw,status");
    let mut ok_count = 0usize;
    let mut noise_count = 0usize;
    let invalid_count = n - rows.len();
    for row in &rows {
        // Quote fields that might contain commas
        println!(
            "\"{}\",{},{},{},{:.2},{}",
            row.name, row.canonical_smiles, row.inchi_key, row.formula, row.mw, row.status
        );
        if row.status == "OK" {
            ok_count += 1;
        } else {
            noise_count += 1;
        }
    }

    eprintln!();
    eprintln!(
        "Enriched {n} compounds in {:.3}s ({:.1} compounds/sec)",
        elapsed.as_secs_f64(),
        n as f64 / elapsed.as_secs_f64()
    );
    eprintln!(
        "Results: {ok_count} OK / {noise_count} noise / {invalid_count} invalid"
    );
}
