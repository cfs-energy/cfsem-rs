# cfsem

[Docs - Rust](https://docs.rs/cfsem) | [Docs - Python](https://cfsem.readthedocs.io/)

Quasi-steady electromagnetics including filamentized approximations, Biot-Savart, and Grad-Shafranov.

To avoid duplication, most tests and example applications are found with the [Python bindings](https://github.com/cfs-energy/cfsem-py).

This library makes use of FMA (fused multiply-add) instructions; as a result, its performance benefits significantly from enabling the `+fma` flag for x86 targets. See [this project's rustc config](.cargo/config.toml) for an example configuration enabling FMA for x86 targets. aarch64 targets (such as Apple and Raspberry Pi processors) have FMA enabled by default.

## Installation

To include this library in a Rust project, add an entry to your Cargo.toml's `[dependencies]` section:

```toml
cfsem = "*"
```

For Python installation, see the docs for the Python library.

## Benchmarking

Benchmarks are configured in Cargo.toml, and can be run via cargo:

```bash
cargo bench
```

## Contributing

Contributions consistent with the goals and anti-goals of the package are welcome.

Please make an issue ticket to discuss changes before investing significant time into a branch.

Goals

* Library-level functions and formulas
* Comprehensive documentation including literature references, assumptions, and units-of-measure
* Quantitative unit-testing of formulas
* Performance (both speed and memory-efficiency)
  * Guide development of performance-sensitive functions with structured benchmarking
* Cross-platform compatibility
* Minimization of long-term maintenance overhead (both for the library, and for users of the library)
  * Semantic versioning
  * Automated linting and formatting tools
  * Centralized CI and toolchain configuration in as few files as possible

Anti-Goals

* Fanciness that increases environment complexity, obfuscates reasoning, or introduces platform restrictions
* Brittle CI or toolchain processes that drive increased maintenance overhead
* Application-level functionality (graphical interfaces, simulation frameworks, etc)

## License

Licensed under either of

* Apache License, Version 2.0, ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
* MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.
