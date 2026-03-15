#!/usr/bin/env python3
"""
Virtual Screening Filter Benchmark
====================================
Reproduces the core bottleneck from Pat Walters' rd_filters
(https://github.com/PatWalters/rd_filters), a production tool used at
Relay Therapeutics and in published drug discovery workflows.

Pipeline per molecule:
  1. SMILES parsing (MolFromSmiles)
  2. Lipinski/Veber descriptor calculation (MW, LogP, HBD, HBA, TPSA, RotBonds)
  3. Substructure alert matching (SMARTS patterns from Glaxo/PAINS rule sets)
  4. Property range filtering

This is the #1 bottleneck pattern in computational chemistry: taking a vendor
catalog (millions of SMILES) and filtering for drug-likeness + structural alerts.

Usage:
  python vs_filter_benchmark.py [SMILES_FILE]

If no file is provided, generates 50,000 drug-like SMILES from a built-in set.
For large-scale testing, download ChEMBL SMILES:
  curl -L 'https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_35_chemreps.txt.gz' | \
    gunzip | tail -n +2 | cut -f2 | head -1000000 > chembl_1M.smi
"""

import sys
import time
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt, MolLogP, NumHDonors, NumHAcceptors, TPSA
from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds

# --- Structural alerts (subset of Glaxo + PAINS from rd_filters) ---
# These are real medicinal chemistry alerts used in drug discovery.
ALERT_SMARTS = [
    # Glaxo reactive group filters (most common)
    "[Br,Cl,I][CX4;CH,CH2]",          # Reactive alkyl halides
    "[S,C](=[O,S])[F,Br,Cl,I]",       # Acid halides
    "C(=O)OC(=O)",                     # Acid anhydrides
    "OO",                              # Peroxides
    "N=C=[S,O]",                       # Isocyanates/Isothiocyanates
    "OS(=O)(=O)C(F)(F)F",             # Triflates
    "[N;R0][N;R0]C(=O)",              # Acylhydrazide
    "[S-]",                            # Thiol anion
    "[SX2]C(=O)",                      # Thioester
    "C#C",                             # Alkynes (terminal)
    "[N;R0]=[N;R0]",                   # Azo compounds
    "SC#N",                            # Thiocyanate
    "C(=O)Oc1c(F)c(F)c(F)c(F)c1(F)", # Pentafluorophenyl esters
    "N=C=N",                           # Carbodiimides
    "[NH1,NH2][CX4][Cl,Br,I]",        # Beta haloamines
    "C1(=O)OCC1",                      # Epoxides
    "[CH1](=O)",                       # Aldehyde
    "[O,S][CH2]O",                     # Hemiacetals
    "C(=O)N(C=O)",                     # Imides
    "[NX3H][NX3H]",                    # Hydrazine
    # PAINS filters (Pan-Assay INterference compounds)
    "c1ccc2c(c1)ccc(=O)o2",           # Coumarin
    "c1ccc(c(c1)O)N=Nc2ccccc2",       # Azo_A
    "[#6]S(=O)(=O)c1ccc(N)cc1",       # Sulfonamide_A
    "c1cc(ccc1N=Nc2ccccc2)O",         # Azo_B
    "c1ccc2c(c1)[nH]c(=O)[nH]2",     # Hydantoin
    "c1ccc2c(c1)c(=O)c3ccccc3o2",     # Flavone
    "O=C1C=CC(=O)C=C1",               # Quinone_A
    "c1cc2c(cc1)oc(=O)c(c2=O)",       # Chromone
    "c1ccc(cc1)S(=O)(=O)N",           # Sulfonamide_B
    "c1cc(oc1)C=NNC(=O)",             # Furan_imine
]

# Drug-like SMILES for self-contained testing (real drugs from DrugBank)
BUILTIN_SMILES = [
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",                    # Ibuprofen
    "CC(=O)Oc1ccccc1C(=O)O",                          # Aspirin
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",                   # Caffeine
    "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C",            # Testosterone
    "CC(C)NCC(O)c1ccc(O)c(O)c1",                       # Isoproterenol
    "c1ccc2c(c1)cc(cc2)C(=O)O",                        # Naproxen-like
    "CC(=O)NC1=CC=C(C=C1)O",                           # Acetaminophen
    "OC(=O)c1ccccc1O",                                 # Salicylic acid
    "c1ccc(cc1)C(=O)NC",                               # Benzamide
    "c1ccc2[nH]c(-c3ccccc3)nc2c1",                     # Benzimidazole
    "Clc1ccc(cc1)C(c2ccccc2)n3ccnc3",                  # Clotrimazole
    "OC(=O)CCc1ccccc1",                                # Hydrocinnamic acid
    "CC(C)(C)NCC(O)c1ccc(O)c(CO)c1",                   # Salbutamol
    "c1ccc(c(c1)C(=O)O)NC2=CC=CC=C2",                 # Fenamate
    "c1ccnc(c1)C(=O)N",                                # Nicotinamide
    "CC(C)OC(=O)C(C)(C)Oc1ccc(Cl)cc1",                # Clofibrate
    "c1ccc2c(c1)c(=O)n(c(=O)[nH]2)C",                 # Barbiturate
    "CC1=CC=C(C=C1)S(=O)(=O)NC(=O)NN2CCCCCC2",        # Tolazamide
    "COc1ccc(cc1OC)C(=O)CC(=O)c2ccc(OC)c(OC)c2",     # Curcumin-like
    "c1ccc(cc1)c2cc(nn2c3ccccc3)C(F)(F)F",            # Celecoxib-like
]


def build_alert_patterns(smarts_list):
    """Pre-compile SMARTS patterns (done once, not per-molecule)."""
    patterns = []
    for sma in smarts_list:
        pat = Chem.MolFromSmarts(sma)
        if pat is not None:
            patterns.append((pat, sma))
    return patterns


def evaluate_molecule(smiles, alert_patterns):
    """
    Core per-molecule evaluation — this is the hot loop.
    Returns (smiles, status, descriptors) or None if invalid.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    # Lipinski + Veber descriptors
    mw = MolWt(mol)
    logp = MolLogP(mol)
    hbd = NumHDonors(mol)
    hba = NumHAcceptors(mol)
    tpsa = TPSA(mol)
    rot = CalcNumRotatableBonds(mol)

    # Structural alert matching
    for pat, desc in alert_patterns:
        if mol.HasSubstructMatch(pat):
            return (smiles, "ALERT:" + desc, mw, logp, hbd, hba, tpsa, rot)

    # Property range filter (Lipinski's Rule of Five + Veber)
    if mw > 500 or logp > 5 or logp < -5 or hbd > 5 or hba > 10 or tpsa > 200 or rot > 10:
        return (smiles, "PROPERTY_FAIL", mw, logp, hbd, hba, tpsa, rot)

    return (smiles, "PASS", mw, logp, hbd, hba, tpsa, rot)


def load_smiles(filename):
    """Load SMILES from file (one per line, optionally with name after whitespace)."""
    smiles_list = []
    with open(filename) as f:
        for line in f:
            parts = line.strip().split()
            if parts:
                smiles_list.append(parts[0])
    return smiles_list


def generate_builtin_dataset(n=50000):
    """Generate n SMILES by cycling through built-in drug molecules."""
    return [BUILTIN_SMILES[i % len(BUILTIN_SMILES)] for i in range(n)]


def benchmark(smiles_list):
    """Run the full pipeline and time each phase."""
    alert_patterns = build_alert_patterns(ALERT_SMARTS)
    n = len(smiles_list)

    print(f"Molecules: {n:,}")
    print(f"Alert patterns: {len(alert_patterns)}")
    print()

    # Phase 1: SMILES parsing only
    t0 = time.perf_counter()
    mols = []
    n_invalid = 0
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            n_invalid += 1
        mols.append(mol)
    t_parse = time.perf_counter() - t0

    # Phase 2: Descriptor calculation only
    t0 = time.perf_counter()
    for mol in mols:
        if mol is None:
            continue
        MolWt(mol)
        MolLogP(mol)
        NumHDonors(mol)
        NumHAcceptors(mol)
        TPSA(mol)
        CalcNumRotatableBonds(mol)
    t_desc = time.perf_counter() - t0

    # Phase 3: Substructure alert matching only
    t0 = time.perf_counter()
    for mol in mols:
        if mol is None:
            continue
        for pat, _ in alert_patterns:
            mol.HasSubstructMatch(pat)
    t_alerts = time.perf_counter() - t0

    # Phase 4: Full pipeline (end-to-end)
    t0 = time.perf_counter()
    n_pass = 0
    n_alert = 0
    n_prop_fail = 0
    for smi in smiles_list:
        result = evaluate_molecule(smi, alert_patterns)
        if result is None:
            continue
        status = result[1]
        if status == "PASS":
            n_pass += 1
        elif status.startswith("ALERT"):
            n_alert += 1
        else:
            n_prop_fail += 1
    t_total = time.perf_counter() - t0

    n_valid = n - n_invalid

    print(f"{'Phase':<30} {'Time (s)':>10} {'us/mol':>10} {'mol/s':>12}")
    print("-" * 65)
    print(f"{'SMILES parsing':<30} {t_parse:>10.3f} {t_parse/n*1e6:>10.1f} {n/t_parse:>12,.0f}")
    print(f"{'Descriptor calculation':<30} {t_desc:>10.3f} {t_desc/n_valid*1e6:>10.1f} {n_valid/t_desc:>12,.0f}")
    print(f"{'Substructure alerts':<30} {t_alerts:>10.3f} {t_alerts/n_valid*1e6:>10.1f} {n_valid/t_alerts:>12,.0f}")
    print(f"{'Full pipeline (end-to-end)':<30} {t_total:>10.3f} {t_total/n*1e6:>10.1f} {n/t_total:>12,.0f}")
    print()
    print(f"Results: {n_pass:,} pass / {n_alert:,} alerts / {n_prop_fail:,} property fail / {n_invalid:,} invalid")
    print()

    # Breakdown as percentage
    t_overhead = t_total - t_parse - t_desc - t_alerts
    print("Time breakdown (approximate):")
    print(f"  Parsing:     {t_parse/t_total*100:5.1f}%")
    print(f"  Descriptors: {t_desc/t_total*100:5.1f}%")
    print(f"  Alerts:      {t_alerts/t_total*100:5.1f}%")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        smiles_list = load_smiles(sys.argv[1])
    else:
        print("No input file provided, using built-in drug-like dataset (50K molecules)")
        print()
        smiles_list = generate_builtin_dataset(50000)

    benchmark(smiles_list)
