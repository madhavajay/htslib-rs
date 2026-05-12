# htslib-rs

`htslib-rs` is a pure Rust HTSlib compatibility project. The current target is
a Rust-only API with test parity against selected HTSlib library behavior. It is
not currently a C ABI replacement and does not provide HTSlib command-line
tools.

The project vendors two upstream codebases as submodules:

- `htslib/`: upstream HTSlib C source and test fixtures.
- `noodles/`: the Rust bioinformatics format implementation this project
  builds on.

## Design

The core rule is: use noodles for file-format semantics wherever possible, and
write Rust compatibility adapters only for HTSlib-observable behavior that
noodles does not already provide.

Production Rust code does not link to HTSlib C. The C source and fixtures are
used as the compatibility reference and test corpus.

## What Uses Noodles

Most format parsing, writing, and indexing is delegated to noodles crates:

| Area | Noodles crate(s) | htslib-rs role |
| --- | --- | --- |
| SAM headers and records | `noodles-sam` | Header/record adapters, view output, aux handling, filters, pileup-related helpers |
| BAM and BAI | `noodles-bam` | BAM read/write/query helpers and HTSlib-style associated-index lookup |
| CRAM and CRAI | `noodles-cram`, `noodles-fasta` | CRAM read/write/query helpers and local FASTA reference repositories |
| VCF | `noodles-vcf` | Header/value adapters, canonical output, indexed queries, mutation helpers |
| BCF | `noodles-bcf` | BCF read/write/query helpers, typed FORMAT/INFO behavior, synced-reader outputs |
| BGZF and GZI | `noodles-bgzf` | HTSlib virtual offsets, GZI load/dump/query, EOF/truncation behavior, worker-count helpers |
| CSI | `noodles-csi` | BAM/SAM/VCF/BCF CSI build and query adapters |
| tabix and TBI | `noodles-tabix` | BED/GFF/VCF text indexing and region-output compatibility |
| FASTA and FAI | `noodles-fasta` | FASTA indexing and HTSlib-compatible retrieval formatting |
| FASTQ and FQI | `noodles-fastq` | FASTQ indexing and FASTQ/FASTA conversion fixtures |

## What Is Implemented From Scratch

The pure Rust compatibility layer implements HTSlib-specific behavior that is
not directly provided by noodles or that needs HTSlib-compatible edge semantics:

- HTSlib-style format detection for local files.
- Associated-index lookup order, including replaced-extension fallback and
  `##idx##` delimiter handling.
- HTSlib region parsing and region-index behavior.
- hFILE utility behavior currently needed by tests, including `haddextension`
  and `data:` URL decoding. Remote hFILE backends are not included.
- Logging level state and `hts_log_*` source-message style validation.
- Utility APIs and test-observable behavior for endian helpers, integer
  parsing, string printing, time helpers, `kstring`, `khash`, `khash_str2int`,
  `kbitset`, `klist`, `kroundup`, `ksort`, and selected `kfunc` math helpers.
- HTS expression parsing/evaluation.
- SAM filter expression behavior, base modification reporting, pileup fixtures,
  and BAQ/probabilistic alignment cases covered by the ported tests.
- VCF/BCF header IDs, typed INFO/FORMAT adapters, allele removal, variant
  classification, sweep behavior, and synced-reader pairing logic.

Rust ownership, error handling, and iterator APIs are intentionally Rust-native.
The project does not try to recreate C allocation/free APIs unless a specific
Rust compatibility test needs the observable behavior.

## Noodles Fork And HTSlib Compatibility Fixes

`htslib-rs` currently points the `noodles` submodule at:

```text
git@github.com:madhavajay/noodles.git
branch: htslib-rs-compat
commit: dca218cd3 Support large HTSlib coordinates and genotypes
```

That fork carries small compatibility fixes needed by HTSlib fixtures.

### BCF GT Integer Width

The upstream HTSlib fixture
`htslib/test/tabix/vcf_file.vcf` includes a large multiallelic record with
genotypes:

```text
GT    0/300    240/260
```

In BCF, genotype alleles are encoded as `(allele_index + 1) << 1`, with a
phasing bit. These values become `602`, `482`, and `522`, which do not fit in
BCF `Int8`.

The noodles fork fixes this by:

- keeping normal genotypes encoded as `Int8` when they fit;
- promoting `GT` values to `Int16` or `Int32` only when required;
- decoding `GT` values from `Int8`, `Int16`, or `Int32`;
- preserving vector-end sentinel handling across integer widths.

The focused noodles test is
`noodles-bcf::record::codec::encoder::samples::values::tests::test_write_genotype_values_with_large_allele_indexes`.
The `htslib-rs` acceptance test that exercises the upstream fixture is
`writes_bcf_round_trip_records_matching_htslib_vcf_fixture` in
`crates/htslib-test-harness/tests/variant_io.rs`.

### Large SAM Coordinates

HTSlib tests also cover large reference coordinates. The noodles fork updates
SAM writing so:

- `@SQ LN` is written without narrowing through signed 32-bit integers;
- record `POS` is written without enforcing the previous writer-side
  `2^31 - 1` limit.

This supports HTSlib large-coordinate fixtures such as
`htslib/test/longrefs/longref.sam`.

## Test Coverage

The Rust acceptance suite lives in `crates/htslib-test-harness/tests`. It ports
HTSlib C tests and scripted library behavior into Rust tests over upstream
fixtures.

The manifest in `crates/htslib-test-harness/src/lib.rs` currently reports:

- 41 passing HTSlib test groups.
- 0 failing groups.
- 0 unported selected groups.
- 1 explicitly out-of-scope group: C ABI/shared-library introspection.

Major covered areas include:

- `test_bgzf.c`: BGZF read/write, EOF marker, embedded EOF, GZI, virtual
  positions, seeking, getline, truncation, compression modes, worker counts.
- `test_faidx.c` and faidx scripts: FASTA/FASTQ indexing and retrieval.
- `test_index.c`: BAI, CSI, TBI, CRAI, explicit CSI min-shift, associated-index
  lookup.
- `test_view.c`: SAM/BAM/CRAM/VCF/BCF read/write/view behavior, record limits,
  indexed regions, multi-region de-duplication, associated index output,
  compressed writes, FASTQ/FASTA output, long BAM records, and benchmark reads.
- `test-vcf-api.c`: VCF header mutation, typed INFO/FORMAT access, rlen
  behavior, allele removal, BCF iterators, serialization.
- `test-bcf-sr.c` and deterministic replacements for `test-bcf-sr.pl`: synced
  reader summaries, BCF/VCF output, weird chromosome region/target queries, and
  allele-pairing modes.
- `test_realn.c`: BAQ tag apply/revert/no-op, MAPQ capping, `probaln_glocal`,
  non-extended and extended BAQ recalculation.
- `test_mod.c`, base-mod scripts, mpileup scripts, SAM filter scripts, FASTQ
  conversion scripts, tabix scripts, TLEN scripts, and VCF various/VCF 4.4
  scripted cases.
- Utility tests for expressions, region parsing, region indexes, endian helpers,
  time helpers, `kfunc`, `khash`, `kstring`, `kroundup`, `ksort`, `klist`,
  `kbitset`, integer parsing, logging, and nibbles.

The CI workflow runs two independent jobs:

- Rust: `cargo fmt --all --check`, `cargo clippy --workspace --all-targets -- -D warnings`,
  and `cargo test --workspace`.
- HTSlib C: configure and build the `htslib` submodule, then run `make check`.

## Missing Or Intentionally Deferred

These are out of scope for the current Rust-only target:

- C ABI compatibility and shared-library export behavior.
- HTSlib command-line tool parity for `bgzip`, `htsfile`, `tabix`, `annot-tsv`,
  and other executables.
- Remote I/O and plugin-backed hFILE backends.
- Generic C `hts_opt` option-string plumbing and shared C thread-pool object
  assignment, beyond format-level worker-count behavior exposed by noodles.
- Low-level CRAM internals such as block/container/codec/slice/transcode APIs
  that are not needed by the selected Rust library tests.
- C preprocessor macros, compiler attributes, allocator hooks, and pointer-level
  container internals that have no useful Rust API equivalent.

Deprecated HTSlib APIs remain tracked in the coverage map because downstream C
libraries may rely on them later, but they are not all exposed in the current
Rust-only API.

## Repository Map

- `crates/htslib-rs`: pure Rust compatibility library.
- `crates/htslib-test-harness`: HTSlib acceptance-test manifest and ported
  integration tests.
- `docs/api-coverage.md`: header and item-level coverage classification.
- `docs/compatibility.md`: compatibility contract and explicit scope decisions.
- `docs/public-api-inventory.md`: generated inventory of HTSlib public headers.
- `TODO.md`: completed porting plan and scope decisions.

## Development

Clone with submodules:

```sh
git clone --recurse-submodules git@github.com:madhavajay/htslib-rs.git
```

Run the Rust gate:

```sh
cargo fmt --all
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
```

Run upstream HTSlib C tests:

```sh
cd htslib
autoreconf -i
./configure --enable-plugins --with-libdeflate
make -j"$(nproc)"
make check
```
