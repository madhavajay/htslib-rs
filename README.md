# htslib-rs

Pure Rust HTSlib-compatible building blocks.

The current target is a Rust-only API with test parity against selected HTSlib tests. Format-specific behavior should use the local `noodles` submodule wherever possible, and this crate should add only the HTSlib-compatible behavior that noodles does not already provide.

See [docs/compatibility.md](docs/compatibility.md) for the current compatibility
contract and [docs/api-coverage.md](docs/api-coverage.md) for the API coverage
map.
