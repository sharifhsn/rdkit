# RDKit

Rust bindings for the [RDKit](https://www.rdkit.org/) cheminformatics C++ library. Provides a high-level, idiomatic Rust API for SMILES parsing, molecule normalization, fingerprinting, substructure search, and more.

This workspace contains two crates:
- **`rdkit`** — Safe, high-level Rust API
- **`rdkit-sys`** — Low-level FFI bindings via [CXX](https://cxx.rs/)

Requires RDKit 2023.09.1 or higher.

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
rdkit = "0.4"
```

You must have the RDKit C++ libraries and headers installed on your system before building. See the platform-specific instructions below.

### macOS (Homebrew)

```bash
brew install rdkit
```

This installs headers to `/opt/homebrew/include/rdkit` (Apple Silicon) or `/usr/local/include/rdkit` (Intel) and libraries are linked automatically.

### Linux (Build from Source)

Install build dependencies:

```bash
# Debian/Ubuntu
sudo apt-get install -y build-essential cmake libboost-iostreams-dev libboost-serialization-dev

# RHEL/Rocky/Fedora
sudo dnf install -y gcc-c++ cmake boost-devel boost-iostreams boost-serialization
```

Download and build RDKit:

```bash
# Download a release (see https://github.com/rdkit/rdkit/releases)
curl -L -O https://github.com/rdkit/rdkit/archive/refs/tags/Release_2024_09_6.tar.gz
tar xf Release_2024_09_6.tar.gz
cd rdkit-Release_2024_09_6

mkdir build && cd build
cmake .. \
  -D RDK_BUILD_PYTHON_WRAPPERS=OFF \
  -D RDK_OPTIMIZE_POPCNT=OFF \
  -D RDK_INSTALL_COMIC_FONTS=OFF \
  -D RDK_BUILD_FREETYPE_SUPPORT=OFF \
  -D RDK_INSTALL_STATIC_LIBS=ON \
  -D RDK_INSTALL_INTREE=OFF \
  -D RDK_BUILD_SWIG_JAVA_WRAPPER=OFF \
  -D RDK_BUILD_CPP_TESTS=OFF

sudo make install -j$(nproc)
```

This installs headers to `/usr/local/include/rdkit` and libraries to `/usr/local/lib` by default.

If you get shared library errors at runtime (`cannot open shared object file`), add the library path:

```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
```

Add that line to your `~/.bashrc` or `~/.zshrc` to make it permanent.

### Linux (Pre-compiled Tarballs)

Pre-compiled tarballs are available for Ubuntu 22.04 (used in CI):

```bash
# AMD64
curl -O https://rdkit-rs-debian.s3.eu-central-1.amazonaws.com/rdkit_2024_09_1_ubuntu_22_04_amd64.tar.gz

# ARM64
curl -O https://rdkit-rs-debian.s3.eu-central-1.amazonaws.com/rdkit_2024_09_1_ubuntu_22_04_arm64.tar.gz
```

Extract and install:

```bash
sudo tar xf rdkit_2024_09_1_ubuntu_22_04_amd64.tar.gz -C /tmp
sudo mv /tmp/rdkit-Release_2024_09_1/Code /usr/local/include/rdkit
sudo mv /tmp/rdkit-Release_2024_09_1/build/lib/* /usr/lib/
```

Note: These are built against Ubuntu 22.04's Boost version. If you're on a different distro or version, building from source (above) is more reliable.

### Conda

Enable the `dynamic-linking-from-conda` feature flag:

```toml
[dependencies]
rdkit = { version = "0.4", features = ["dynamic-linking-from-conda"] }
```

Then install RDKit in your conda environment:

```bash
conda install -c conda-forge rdkit
```

The build script reads `CONDA_PREFIX` to find headers and libraries. Make sure your conda environment is activated when running `cargo build`.

### Custom Install Location

If RDKit is installed somewhere non-standard, the build script checks these paths in order:

| Platform | Header paths | Library paths |
|----------|-------------|---------------|
| macOS ARM | `/opt/homebrew/include`, `/opt/homebrew/include/rdkit` | `/opt/homebrew/lib` |
| macOS x86 | `/usr/local/include`, `/usr/local/include/rdkit` | (system default) |
| Linux | `/usr/local/include`, `/usr/local/include/rdkit`, `/usr/include`, `/usr/include/rdkit` | (system default) |

## Usage

```rust
use rdkit::ROMol;

let mol = ROMol::from_smiles("c1ccccc1").unwrap();
println!("{}", mol.as_smiles());
```

## Building & Testing

```bash
cargo build --workspace
cargo test --workspace
cargo clippy --workspace
cargo fmt --check
```

## Formatting

### C++ formatting

C++ wrapper code uses clang-format (LLVM style, 80 columns). To format:

```bash
clang-format -i rdkit-sys/wrapper/src/*.cc rdkit-sys/wrapper/include/*.h
```

### Rust formatting

```bash
cargo fmt
```

## Releasing

The `rdkit-sys` crate is a member of the `rdkit` workspace. Both crates move in lockstep versions:

```bash
cargo workspaces version patch
cd rdkit-sys && cargo publish && cd ..
cargo publish
```
