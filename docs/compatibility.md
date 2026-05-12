# Compatibility Contract

This port targets a pure Rust API with test parity against selected HTSlib
library behavior. It is not a C ABI replacement and does not currently provide
HTSlib command-line tools.

## Compatible With HTSlib C

The Rust helpers preserve HTSlib-observable behavior where the ported tests
depend on it, including:

- Local format detection for selected SAM, BAM, CRAM, VCF, FASTA, FASTQ, BED,
  and index fixtures.
- BGZF reading, writing, virtual offsets, EOF marker handling, GZI indexes,
  write modes, and worker-count-backed threading behavior.
- Local BAI, CSI, TBI, CRAI, FAI, and GZI index read/build/write/query behavior
  covered by the acceptance tests.
- HTSlib-style associated-index lookup order, including replaced-extension
  fallback and `##idx##` delimiter handling.
- Region parsing for the ported `hts_parse_reg` modes and edge cases.
- SAM/BAM/CRAM view/query behavior covered by upstream fixtures, including raw
  header preservation, indexed region output, multi-region de-duplication, and
  CRAM external-reference decoding through noodles' reference repository API.
- VCF/BCF header, typed INFO/FORMAT, span, allele-removal, serialization,
  translation, iterator, classification, sweep, and selected synced-reader
  behavior covered by the Rust acceptance tests.
- FASTA/FASTQ indexing, retrieval, and conversion golden-output behavior covered
  by upstream fixtures.
- HTSlib utility behavior currently covered by tests: expression evaluation,
  logging source-message style, region indexes, integer parsing/string printing,
  endian helpers, time helpers, Fisher exact tests, dynamic string behavior, hash
  behavior, and SAM base-modification reporting.

## Rust-Native By Design

The public target is a Rust API, so these pieces intentionally differ from C:

- Production code uses noodles crates and Rust standard-library abstractions
  instead of linking to HTSlib C.
- Errors are Rust `Result` values and typed Rust errors where available, rather
  than `errno`-driven return codes.
- Ownership and lifetime behavior follows Rust types rather than C allocation
  and destroy functions.
- Shared mutable C structures such as global plugin registries, C ABI
  introspection, and C-style thread-pool objects are not modeled unless a Rust
  compatibility test needs them.
- Deprecated HTSlib APIs remain tracked for compatibility decisions, but they
  are not automatically exposed as C-shaped Rust functions.

## Out Of Scope For Now

These are intentionally deferred by project decision:

- Command-line tool parity for `bgzip`, `htsfile`, `tabix`, `annot-tsv`, and
  other HTSlib executables.
- Remote I/O and plugin-backed hFILE backends.
- C ABI or shared-library export compatibility.
- Full item-level public API coverage; the current inventory is in
  [`public-api-inventory.md`](public-api-inventory.md), and final per-item
  classification is still open.

## Remaining Compatibility Work

Known incomplete areas are tracked in [`TODO.md`](../TODO.md), including
realignment, broader `test_view.c` and `test-bcf-sr.c` coverage, mutation
adapters, differential-test decisions, and final item-level API classification.
