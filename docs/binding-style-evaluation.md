# Binding Style Evaluation: Issue #39

Evaluation of [#39 — Proposal: Simplify Rust Bindings](https://github.com/rdkit-rs/rdkit/issues/39),
which proposes switching from free-function FFI declarations to CXX method syntax.

## Background

The proposal comes from [@mvisani](https://github.com/mvisani), who built a
[proof-of-concept repo](https://github.com/mvisani/rdkit-rust-sys) demonstrating
the approach with `Atom` (~25 methods) and `ROMol` (~7 methods). The maintainer
(@xrl) was positive; Marco said he'd open a PR but one was never submitted
(Oct 2024).

## Current State

The codebase has **115 FFI declarations across 9 bridge modules** and
**93 C++ wrapper functions**. Two distinct patterns coexist:

### Pattern A: Free functions with SharedPtr (7 of 9 modules)

```rust
// rdkit-sys/src/bridge/ro_mol.rs
pub fn get_atomic_num(atom: Pin<&Atom>) -> i32;
pub fn mol_to_smiles(mol: &SharedPtr<ROMol>) -> String;
```

```cpp
// rdkit-sys/wrapper/src/ro_mol.cc
int get_atomic_num(const Atom &atom) { return atom.getAtomicNum(); }
rust::String mol_to_smiles(const std::shared_ptr<ROMol> &mol) {
    return MolToSmiles(*mol);
}
```

### Pattern B: Method syntax with self receiver (periodic_table module)

```rust
// rdkit-sys/src/bridge/periodic_table.rs
pub fn getAtomicWeight(self: &PeriodicTable, atomic_number: u32) -> f64;
pub fn getRvdw(self: &PeriodicTable, atomic_number: u32) -> f64;
```

No C++ wrapper needed — CXX calls the C++ methods directly.

**The codebase is already a hybrid.** The periodic_table module uses exactly the
approach Marco proposes. It works, compiles, and all 17 tests pass.

## What Marco's Approach Changes

The core difference is in the CXX bridge declaration:

```rust
// Current (Pattern A):
pub fn get_atomic_num(atom: Pin<&Atom>) -> i32;

// Proposed (Pattern B):
pub fn getAtomicNum(self: &Atom) -> i32;
```

Both approaches still need C++ wrapper functions for anything beyond trivial
getters (constructors, type conversions, template instantiation). The difference
is purely at the bridge declaration layer.

## Detailed Analysis

### Advantages of method syntax

1. **Eliminates C++ wrappers for simple getters.** When the Rust function name
   matches the C++ method name, CXX calls it directly — no `.cc` wrapper needed.
   This could remove ~40% of current wrapper code (the simple
   getter/setter pass-throughs like `get_atomic_num`, `get_is_aromatic`,
   `get_formal_charge`, etc.).

2. **More readable bridge code.** `self: &Atom` reads better than
   `atom: Pin<&Atom>` as a free function parameter. The intent (this is a
   method on Atom) is immediately clear.

3. **Proven in this codebase.** The periodic_table module already uses method
   syntax and works well. This isn't a theoretical proposal — it's extending an
   existing pattern.

4. **Enables `impl` blocks in bridge files.** Marco adds convenience methods in
   `impl Atom { ... }` directly in the bridge file, which can reduce the need
   for a separate safe wrapper layer for simple types.

5. **Better error handling.** Marco marks potentially-throwing methods with
   `-> Result<T>`, making error handling explicit. The current codebase has gaps
   here — issue #18 (segfault on invalid SMARTS) was caused by a missing null
   check that `Result<T>` would have caught at compile time.

### Disadvantages and risks

1. **Not actually less total code for complex methods.** Constructors, template
   instantiations, type conversions, and anything that doesn't map 1:1 to a C++
   method still needs a C++ wrapper function. For these, method syntax saves
   nothing.

2. **Naming convention tension.** CXX method syntax requires the Rust function
   name to match the C++ method name by default. RDKit uses camelCase
   (`getAtomicNum`), Rust uses snake_case (`get_atomic_num`). Options:
   - Accept camelCase in the bridge (periodic_table does this today)
   - Use `#[rust_name = "get_atomic_num"]` on each method (verbose but correct)
   - The safe wrapper layer re-exports as snake_case either way, so users of
     the `rdkit` crate never see this

3. **Ownership model mismatch.** Marco's PoC uses `UniquePtr<Atom>` for
   standalone atom construction. The current codebase uses `SharedPtr` everywhere,
   and atoms are always borrowed from their parent molecule. Mixing these would
   create two incompatible Atom types. This is a design decision, not a
   limitation of method syntax itself.

4. **Migration risk.** A wholesale migration would touch all 115 FFI declarations
   and 93 C++ wrappers — every file in both crates. High regression risk for
   zero capability gain (the public API surface doesn't change).

5. **CXX limitations apply to both approaches.**
   - `&mut T` is forbidden on opaque C++ types; must use `Pin<&mut T>` regardless
   - `Vec<SharedPtr<T>>` doesn't work across the FFI boundary regardless
   - Template functions must be explicitly instantiated in C++ regardless

### Side-by-side: What changes per module

| Module | FFI functions | C++ wrappers eliminable | Migration effort |
|--------|:---:|:---:|---|
| ro_mol | 29 | ~12 (simple getters) | High — largest module, atom access is complex |
| mol_ops | 19 | ~13 (getter/setter pairs) | Medium — repetitive but mechanical |
| periodic_table | 18 | Already done | N/A — already uses method syntax |
| substruct_match | 15 | ~7 (getter/setter pairs) | Medium |
| mol_standardize | 11 | ~3 | Low |
| scaffold_network | 10 | ~2 | Low |
| rw_mol | 5 | ~1 | Low |
| fingerprint | 5 | ~1 | Low |
| descriptors | 3 | ~1 | Low |
| **Total** | **115** | **~40** | |

## Recommendation

**Don't do a wholesale migration. Adopt method syntax for new modules going
forward.**

Rationale:
- The codebase already proves method syntax works (periodic_table has used it
  since the module was written)
- A full migration is high-effort, high-risk, zero-capability-gain
- New modules can use method syntax naturally — it's a per-module choice, not
  all-or-nothing
- The safe Rust wrapper layer (`rdkit` crate) hides the FFI style from end
  users anyway, so the public API is unaffected

### Guidelines for new bridge modules

1. **Use method syntax** (`self: &T`, `self: Pin<&mut T>`) when the C++ method
   name maps cleanly
2. **Use `#[rust_name = "snake_case"]`** to keep Rust-idiomatic names where
   it matters
3. **Use `-> Result<T>`** for any C++ function that can throw
4. **Keep SharedPtr ownership** — don't introduce UniquePtr-based standalone
   objects (atoms, bonds) since the molecule-owns-atoms model is safer and
   consistent with the existing codebase
5. **Skip C++ wrappers** for trivial getters where CXX can call the C++ method
   directly

### Opportunistic migration

If a module is being substantially modified for other reasons, it's reasonable
to switch it to method syntax at the same time. The highest-value targets would
be `mol_ops` (13 repetitive getter/setter pairs) and `substruct_match` (7
pairs), where the boilerplate reduction is most noticeable.

## References

- Issue #39: https://github.com/rdkit-rs/rdkit/issues/39
- Marco's PoC: https://github.com/mvisani/rdkit-rust-sys
- CXX method syntax: https://cxx.rs/extern-c++.html
- CXX SharedPtr: https://cxx.rs/binding/sharedptr.html
- BACKLOG.md Vec<SharedPtr> workaround: https://github.com/dtolnay/cxx/issues/741
- Periodic table bridge (existing method syntax): `rdkit-sys/src/bridge/periodic_table.rs`
