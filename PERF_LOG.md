# Workflow Optimization Ideas

Post-merge roadmap for the Rust RDKit wrapper. The theme: **the Rust layer should be a workflow engine, not just a function-by-function FFI mirror.** Compute in C++ once, query/filter/search in Rust many times.

Inspired by [cheminee](https://github.com/rdkit-rs/cheminee), which pre-computes fingerprints + 43 descriptors into a Tantivy index at index time, then runs all queries in pure Rust with zero C++ overhead per query.

## Priority 1: SubstructLibrary

RDKit has `SubstructLibrary` (`GraphMol/SubstructLibrary/SubstructLibrary.h`) — a purpose-built collection that stores molecules + pattern fingerprints and does automatic FP pre-screening before graph isomorphism. Turns O(N*M) brute-force substruct matching into O(N*candidates).

Multiple storage backends: `MolHolder` (in-memory, fastest), `CachedMolHolder` (binary pickled, balanced), `CachedSmilesMolHolder` (SMILES, smallest memory).

This is the single biggest win for the VS filter use case.

## Priority 2: FingerprintArena + bulk Tanimoto

Contiguous packed fingerprint storage (not individual `BitVec`s). Enables:
- Bulk Tanimoto: one query against N targets
- Pairwise similarity matrix
- Diversity selection (MaxMinPicker)

RDKit has `FPBReader` / `MultiFPBReader` (`DataStructs/FPBReader.h`) — on-disk memory-mapped fingerprint databases with built-in popcount screening and `getTanimotoNeighbors(threshold)`. Searches millions of FPs without loading molecules.

Also: `CalcBitmapTanimoto` / `CalcBitmapAllProbeBitsMatch` in `BitOps.h` operate directly on byte arrays. 12 similarity metrics beyond Tanimoto (Dice, Tversky, Cosine, etc.).

`SimDivPickers/MaxMinPicker` does diversity selection without computing full distance matrix (`lazyPick` computes distances on demand).

## Priority 3: MolPickler binary serialization

`MolPickler` (`GraphMol/MolPickler.h`) — binary molecule serialization, ~7x faster than SMILES round-tripping. Configurable via `PropertyPickleOptions` flags (NoProps, MolProps, AtomProps, NoConformers, etc.). Critical for workflows that reload the same molecules repeatedly.

## Priority 4: Pipeline builder

Pure Rust orchestration layer. Scientists write declarative pipelines, internals handle parallelism, error recovery, batching:

```rust
MolPipeline::from_smiles(smiles_iter)
    .filter_parse_errors()
    .standardize()
    .par_chunks(1000)
    .compute(|mol| mol.lipinski_descriptors())
    .collect()
```

## Priority 5: Batch RDKit APIs (already exist in C++)

- `FingerprintGenerator::getFingerprints(vector<ROMol*>, numThreads)` — batch FP generation with built-in threading
- `cleanupInPlace(vector<RWMol*>, numThreads)` — batch standardization
- `RWMol::beginBatchEdit()` / `commitBatchEdit()` / `rollbackBatchEdit()` — atomic batch molecule editing
- `EmbedMultipleConfs()` with `numThreads` — parallel conformer generation

## Priority 6: RGroupDecomposition

`RGroupDecomposition` (`GraphMol/RGroupDecomposition/RGroupDecomp.h`) — critical for medicinal chemistry SAR. Decomposes a set of analogs against core scaffold(s), returns row-oriented or column-oriented R-group tables. Batch API: `RGroupDecompose(cores, mols)`.

## Lower priority

- **FMCS** (`GraphMol/FMCS/FMCS.h`) — maximum common substructure with timeout and progress callbacks
- **ScaffoldNetwork** (`GraphMol/ScaffoldNetwork/ScaffoldNetwork.h`) — hierarchical scaffold analysis with custom bond-breaking rules
- **Chemical reactions** (`GraphMol/ChemReactions/Reaction.h`) — `runReactants()`, plus `ChemTransforms.h` for `deleteSubstructs` / `replaceSubstructs`
- **MolDraw2D** — SVG/PNG molecule depiction with atom highlighting
- **Advanced 3D descriptors** — AUTOCORR3D, RDF, MORSE, WHIM, GETAWAY
- **Reduced graphs** — ErG fingerprints for pharmacophore-level similarity

## Key design principle

From cheminee: separate "RDKit invocation" (expensive, C++, done once) from "search/filter execution" (cheap, pure Rust, done many times). The Rust layer's job is to minimize how often you cross into C++.
