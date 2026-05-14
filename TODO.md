# TODO: Port HTSlib to Pure Rust

Goal: build a pure Rust replacement for HTSlib's C implementation, port HTSlib's tests, and get the ported suite passing. Prefer using `noodles` crates wherever they cover the required behavior instead of reimplementing format logic from scratch.

## Current Inputs

- `htslib/`: upstream C HTSlib source and test suite.
- `noodles/`: upstream Rust bioinformatics format crates.
- HTSlib source areas reviewed:
  - public headers in `htslib/htslib/*.h`
  - core library objects listed in `htslib/Makefile`
  - test programs and scripted suites under `htslib/test`
- Noodles crates reviewed:
  - `noodles-bam`, `noodles-bcf`, `noodles-bgzf`, `noodles-cram`, `noodles-csi`
  - `noodles-fasta`, `noodles-fastq`, `noodles-sam`, `noodles-tabix`, `noodles-vcf`
  - `noodles-core`, plus the top-level `noodles` re-export crate

## Porting Principles

- Keep the implementation pure Rust.
- Use noodles as the default implementation for file-format semantics:
  - SAM: `noodles-sam`
  - BAM and BAI: `noodles-bam`
  - CRAM and CRAI: `noodles-cram`
  - VCF: `noodles-vcf`
  - BCF: `noodles-bcf`
  - BGZF and GZI: `noodles-bgzf`
  - CSI: `noodles-csi`
  - tabix and TBI: `noodles-tabix`
  - FASTA and FAI: `noodles-fasta`
  - FASTQ and FQI: `noodles-fastq`
- Only implement new Rust code when HTSlib behavior is not already available through noodles or a small adapter over noodles.
- Preserve HTSlib-observable behavior where tests depend on it, including edge cases, parse errors, output formatting, region handling, and index lookup behavior.
- Avoid binding to the C library for production code. C may be used temporarily only as an oracle in tests during the migration.
- Treat each ported HTSlib test as an acceptance test. Do not mark a module complete until its relevant tests pass.

## Phase 0: Project Skeleton

- [x] Create a root Rust workspace.
- [x] Add the initial crate layout.
  - [x] Library crate for the pure Rust HTSlib-compatible implementation.
  - [x] Integration test crate or test harness for ported HTSlib tests.
- [x] Add path dependencies to the checked-out `noodles` crates.
- [x] Add CI commands for:
  - [x] `cargo fmt --check`
  - [x] `cargo clippy --workspace --all-targets`
  - [x] `cargo test --workspace`
  - [x] ported HTSlib acceptance tests
- [x] Document that the initial target is a Rust-only API.

## Phase 1: API Inventory and Coverage Map

- [x] Inventory every public HTSlib header in `htslib/htslib`.
- [x] Classify each public API item as:
  - [x] covered directly by noodles
  - [x] covered by a small adapter over noodles
  - [x] requires a new Rust implementation
  - [x] deprecated but retained for compatibility
  - [x] intentionally out of scope
  - [x] Expand the header-level classification to individual public functions, types, constants, and macros.
    - [x] Add a generated candidate inventory for public declarations and macros.
    - [x] Add item-level classification for `hts_endian.h`.
    - [x] Add item-level classification for `hts_log.h`.
    - [x] Add item-level classification for `hts_os.h`.
    - [x] Add item-level classification for `hts_defs.h`.
    - [x] Add item-level classification for `kroundup.h`.
    - [x] Add item-level classification for `kfunc.h`.
    - [x] Add item-level classification for `kbitset.h`.
    - [x] Add item-level classification for `klist.h`.
    - [x] Add item-level classification for `ksort.h`.
    - [x] Add item-level classification for `khash.h`.
    - [x] Add item-level classification for `khash_str2int.h`.
    - [x] Add item-level classification for `kstring.h`.
    - [x] Add item-level classification for `vcfutils.h`.
    - [x] Add item-level classification for `vcf_sweep.h`.
    - [x] Add item-level classification for `thread_pool.h`.
    - [x] Add item-level classification for `kseq.h`.
    - [x] Add item-level classification for `knetfile.h`.
    - [x] Add item-level classification for `regidx.h`.
    - [x] Add item-level classification for `hts_expr.h`.
    - [x] Add item-level classification for `synced_bcf_reader.h`.
    - [x] Add item-level classification for `tbx.h`.
    - [x] Add item-level classification for `faidx.h`.
    - [x] Add item-level classification for `hfile.h`.
    - [x] Add item-level classification for `hts.h`.
    - [x] Add item-level classification for `sam.h`.
    - [x] Add item-level classification for `vcf.h`.
    - [x] Add item-level classification for `bgzf.h`.
    - [x] Add item-level classification for `cram.h`.
- [x] Build a coverage matrix for:
  - [x] `hts.h`: format-neutral I/O, format detection, indexes, iterators, regions
  - [x] `hfile.h`: local and remote I/O abstraction
  - [x] `bgzf.h`: BGZF reader, writer, virtual offsets, threading, indexes
  - [x] `sam.h`: SAM/BAM/CRAM records, headers, iterators, pileup, realignment, base modifications
  - [x] `vcf.h`: VCF/BCF headers, records, typed values, translation, variant classification
  - [x] `tbx.h`: tabix indexes and region queries
  - [x] `faidx.h`: FASTA/FASTQ indexing and sequence retrieval
  - [x] `regidx.h`: region index parsing and lookup
  - [x] `hts_expr.h`: filter expression parsing and evaluation
  - [x] `thread_pool.h`: worker pool behavior or Rust equivalent
  - [x] utility headers: `kstring.h`, `khash.h`, `kfunc.h`, `hts_endian.h`, `hts_log.h`, `vcfutils.h`, `vcf_sweep.h`, `synced_bcf_reader.h`
    - [x] Cover `hts_log.h` source-message style behavior used by `test-logging.pl`.
    - [x] Cover `hts_log.h` log-level state behavior for the Rust-only API target.

## Phase 2: Test Harness

- [x] Preserve HTSlib test data under `htslib/test` as fixtures.
- [x] Port C unit tests to Rust integration tests.
  - [x] `test_bgzf.c`
    - [x] Port BGZF read/write round-trip and EOF marker checks.
    - [x] Port embedded EOF block read-through behavior.
    - [x] Port GZI load/dump/query and uncompressed seek checks.
    - [x] Port virtual-position seek/read and getline-style reads.
    - [x] Port truncated stream error checks.
    - [x] Port gzip and uncompressed write-mode checks.
    - [x] Port thread-pool variants as Rust worker-count BGZF read/write helpers backed by noodles.
    - [x] Mark C-specific open wrappers out of scope for the Rust-only API target.
  - [x] `test_faidx.c`
    - [x] Port uncompressed FASTA/FASTQ FAI generation and region retrieval golden outputs.
    - [x] Port bgzip-compressed FASTA retrieval with FAI/GZI.
  - [x] `test_index.c`
    - [x] Port BAI build coverage for BAM fixtures.
    - [x] Port TBI and CSI build/query coverage for VCF fixtures.
    - [x] Port BCF CSI build/query coverage.
    - [x] Port explicit VCF CSI min_shift coverage.
    - [x] Port associated-index lookup fallback cases for existing local indexed query helpers.
    - [x] Port SAM.gz CSI build/query coverage.
    - [x] Port explicit CSI min_shift variants for BAM, SAM.gz, VCF, and BCF.
    - [x] Port SAM.gz BAI build/query coverage.
    - [x] Add CRAI read/write round-trip coverage.
    - [x] Port CRAM/CRAI.
  - [x] `test_expr.c`
  - [x] `test_kfunc.c`
  - [x] `test_khash.c`
  - [x] `test_kstring.c`
  - [x] `test_realn.c`
    - [x] Port existing BAQ tag apply/revert behavior for `BQ:Z`/`ZQ:Z` without recalculation.
    - [x] Port existing BAQ no-op behavior for already-applied and already-unapplied records.
    - [x] Port `sam_cap_mapq` MAPQ capping behavior for matched, mismatched, and clipped SAM records.
    - [x] Port initial `probaln_glocal` likelihood, state, and posterior-quality behavior for exact, mismatch, and insertion cases.
    - [x] Port default non-extended BAQ recalculation for `realn01` and `realn02` against upstream expected SAM output.
    - [x] Port non-extended BAQ recalculation with quality application for `realn01` and `realn02` against upstream expected SAM output.
    - [x] Port forced non-extended BAQ recomputation for `realn02-r.sam` against upstream expected SAM output.
    - [x] Port extended BAQ recalculation for `realn01`, `realn02`, and adjacent-match `realn03` against upstream expected SAM output.
  - [x] `test-regidx.c`
  - [x] `test-parse-reg.c`
  - [x] `test_str2int.c`
  - [x] `test_time_funcs.c`
  - [x] `hts_endian.c`
  - [x] `test_view.c`
    - [x] Port SAM view record-limit behavior (`-N`) with header preservation.
    - [x] Port SAM parse-error ignoring behavior (`-I`) with valid-record output preservation.
    - [x] Port BGZF-compressed SAM view output and record-limit behavior.
    - [x] Port parsed SAM-to-BGZF-SAM compressed write behavior (`-z`) with record-limit behavior.
    - [x] Port SAM-to-FASTQ and SAM-to-FASTA view output with record-limit behavior (`-f`, `-F`, `-N`).
    - [x] Port BAM no-compression output behavior (`-u`/`-l 0`) with readable BGZF/BAM output.
    - [x] Port SAM-to-BAM write with associated BAI index output behavior (`-b -x`).
    - [x] Port BAM whole-file view record-limit behavior (`-N`) with raw header preservation.
    - [x] Port BAM indexed region view output against upstream `range.out`.
    - [x] Port BAM indexed region record-limit behavior (`-N`) against upstream `range.out`.
    - [x] Port BAM multi-region view output against upstream `range.out2`, including HTSlib-style de-duplication.
    - [x] Port generated long BAM record round-trip behavior for CIGAR and sequence data crossing BGZF block boundaries.
    - [x] Port CRAM indexed region and multi-region view output against upstream `range.out` and `range.out2`, including reference-backed MD/NM regeneration.
    - [x] Port SAM-to-CRAM write with associated CRAI index output behavior (`-C -x`).
    - [x] Port CRAM whole-file view record-limit behavior (`-N`) with reference-backed MD/NM regeneration.
    - [x] Port CRAM indexed region record-limit behavior (`-N`) against upstream `range.out`, including reference-backed MD/NM regeneration.
    - [x] Port padded BAM raw-header parsing/view behavior from `test_convert_padded_header`.
    - [x] Port whole-file BCF-to-VCF view output against upstream `tabix/vcf_file.vcf`.
    - [x] Port VCF and BCF whole-file variant record-limit behavior (`-N`) with raw header preservation.
    - [x] Port VCF and BCF indexed region record-limit behavior (`-N`).
    - [x] Port indexed VCF region header output for nonzero contig `IDX` fixtures.
    - [x] Port VCF compressed write with associated TBI index output behavior (`-z -x`).
    - [x] Port VCF-to-BCF write with associated CSI index output behavior (`-b -x`).
    - [x] Port benchmark/no-output record-reading behavior (`-B`) for SAM, BAM, CRAM, VCF, and BCF.
    - [x] Mark C harness plumbing for output-path routing (`-p`), generic `hts_opt` option strings (`-i`/`-o`), and shared thread-pool assignment (`-@`) out of scope for the Rust-only/no-CLI target.
  - [x] `test-vcf-api.c`
    - [x] Port header get/remove cases for FILTER, INFO, FORMAT, contig, and generic records.
    - [x] Port `vcf_open_mode` extension classification cases.
    - [x] Port INFO/FORMAT `Number` classification cases, including `P`, `LA`, `LG`, `LR`, and `M`.
    - [x] Port typed INFO value retrieval for integer, float, and missing float values.
    - [x] Port typed FORMAT integer retrieval for scalar, vector, and missing values.
    - [x] Port BCF FORMAT float missing and vector-end retrieval behavior.
    - [x] Port record `rlen` calculation table for VCF 4.3, 4.4, and 4.5 inputs, plus mutation-triggered rlen recalculation.
    - [x] Port invalid `END` tag `rlen` behavior.
    - [x] Port `bcf_remove_allele_set` ALT removal, INFO/FORMAT vector trimming, local allele, and genotype remapping cases.
    - [x] Port BCF indexed iterator creation by numeric reference ID and region string.
    - [x] Port record serialization output for header removal, record sync/duplication, and INFO/FORMAT mutation cases.
  - [x] `test-vcf-sweep.c`
  - [x] `test-bcf-sr.c`
    - [x] Port no-index multi-VCF summary output and bad header/record order detection.
    - [x] Port no-index VCF output mode (`-O vcf`) for synced local VCF inputs.
    - [x] Port no-index BCF output mode (`-O bcf`) for synced local VCF inputs.
    - [x] Port Rust reader-addition summary parity equivalent to the C `--usefptr` path for the Rust-only API.
    - [x] Port weird chromosome-name region and target query cases.
    - [x] Port indexed BCF weird chromosome-name region and target query cases.
    - [x] Port deterministic allele-pairing coverage for `snps`, `indels`, `both`, `snps+ref`, `indels+ref`, `both+ref`, `exact`, `some`, and `all` logic.
    - [x] Port a deterministic `test-bcf-sr.pl` randomized-shape allele-pairing fixture with duplicate reader groups, multi-variant groups, multiallelic rows, and all pairing modes.
    - [x] Port shuffled group-order stability coverage for deterministic `test-bcf-sr.pl` randomized-shape allele-pairing fixtures.
    - [x] Port shuffled per-group variant and duplicate-reader input-order stability coverage for deterministic `test-bcf-sr.pl` randomized-shape allele-pairing fixtures.
    - [x] Port additional deterministic `test-bcf-sr.pl` randomized-shape duplicate-reader input-order coverage across multiple refs, variant mixes, duplicate counts, and all pairing modes.
    - [x] Replace stochastic `test-bcf-sr.pl` sort-loop parity with deterministic randomized-shape fixtures; CLI file-list/output routing remains out of scope for the Rust-only/no-CLI target.
  - [x] `test-bcf-translate.c`
  - [x] `test-bcf_set_variant_type.c`
  - [x] `test_mod.c`
    - [x] Port `MM-variants.sam` base-modification reporting with unchecked-base output.
    - [x] Port `MM-explicit.sam` and `MM-double.sam` reporting with ML probabilities.
    - [x] Port `MM-explicit-x.out` extended metadata reporting.
    - [x] Port `MM-chebi.sam` reporting with ChEBI IDs and `N` canonical-base skip semantics.
    - [x] Port `MM-multi.sam` and `MM-not-all-modded.sam` reporting cases.
    - [x] Port invalid `MN` and MM bounds rejection cases.
  - [x] `test_nibbles.c`
  - [x] `test_introspection.c`
    - [x] Mark C ABI/shared-library introspection out of scope for the Rust-only API target.
- [x] Port scripted tests to Rust integration tests where they exercise library behavior.
  - [x] `test/test.pl`
    - [x] Port BAM view and indexed region view behavior for BAM files whose binary reference table has references missing from raw `@SQ` header lines.
    - [x] Port large-position BGZF SAM CSI iterator output, including multi-region de-duplication.
    - [x] Port BAM multi-region view output for `range.bam`.
    - [x] Port CRAM indexed region and multi-region view output for `range.cram`.
    - [x] Port indexed BCF weird chromosome-name region and target cases from `test_bcf_sr_range`.
    - [x] Port VCF canonical-output cases from `test_vcf_various` for `formatcols.vcf`, `noroundtrip.vcf`, `formatmissing.vcf`, and `vcf_meta_meta.vcf`.
    - [x] Port VCF header canonicalization from `test_vcf_various` for `test-vcf-hdr-in.vcf`.
    - [x] Complete the `test_vcf_various` scripted subgroup.
    - [x] Port the `test_vcf_44` VCF 4.4 implicit/explicit phasing normalization fixture (`vcf44_1.vcf`).
    - [x] Port the library-relevant `test_rebgzip` multi-block BGZF/GZI boundary and indexed-seek behavior.
    - [x] Port generated long BAM record round-trip behavior from `test_view`.
    - [x] Port large-reference VCF CSI query behavior for `longrefs/index.vcf`, including adaptive CSI depth and `INFO/END` span overlap.
    - [x] Port `index2.sam` mapped/unmapped pair BAM index query count behavior from `test_index`.
    - [x] Mark CRAM version/profile/`hts_opt` option-string matrices from `test_view`, `test_multi_ref`, and `test_MD` out of scope beyond representative noodles-backed CRAM read/write/query/reference parity.
    - [x] Mark plugin-loading, remote reference-cache server, `annot-tsv`, `bgzip`, `htsfile`, and script command-routing checks out of scope for the Rust-only/no-CLI/no-remote-I/O target.
  - [x] `test/faidx/test-faidx.sh`
  - [x] `test/fastq/test-fastq.sh`
    - [x] Port core FASTQ/FASTA to SAM read conversions.
    - [x] Port multiline FASTQ/FASTA read conversions.
    - [x] Port aux passthrough, interleaved `/1` and `/2`, `FASTQ_NAME2`, and UMI read conversions.
    - [x] Port core SAM to FASTQ/FASTA write conversions with aux and read-number output.
    - [x] Port CASAVA/filter-specific FASTQ conversions.
  - [x] `test/tabix/test-tabix.sh`
    - [x] Port TBI/CSI VCF, large-coordinate CSI, BED, GFF, and `--separate-regions` library-output cases.
    - [x] Mark remaining `tabix`/`bgzip` CLI and thread-option behavior out of scope for the Rust-only API target.
  - [x] `test/mpileup/test-pileup.sh`
    - [x] Port deletion-only pileup output fixture (`mp_D.sam` / `mp_D.out`) for the Rust library helper.
    - [x] Port SAM insertion, insertion/deletion, refskip, and pad pileup fixtures.
    - [x] Port SAM overlap-removal pileup fixtures.
    - [x] Port BAM edge-case fixture (`small.bam`).
  - [x] `test/base_mods/base-mods.sh`
    - [x] Port `MM-variants.sam` / `MM-variants.out` library reporting case.
    - [x] Port `MM-explicit.out`, `MM-explicit-x.out`, `MM-explicit-f.out`, and `MM-double.out` library reporting cases.
    - [x] Port `MM-chebi.out` library reporting case.
    - [x] Port `MM-multi.out` and `MM-not-all-modded.out` library reporting cases.
    - [x] Port `MM-pileup.out`, `MM-pileup2.out`, and `MM-MNp.sam` pileup reporting cases.
    - [x] Port `MM-MNf1.sam`, `MM-MNf2.sam`, `MM-bounds+.sam`, and `MM-bounds-.sam` rejection cases.
  - [x] `test/sam_filter/filter.sh`
    - [x] Port integer-expression record-count cases.
    - [x] Port aux-tag record-count case.
    - [x] Port string-expression golden-output cases.
    - [x] Port function and CIGAR-metric cases.
  - [x] `test/tlen/tlen.sh`
    - [x] Port CRAM TLEN auto-creation fixtures against upstream SAM expected outputs.
- [x] Add a runner that reports ported, passing, failing, and unported HTSlib tests.
- [x] Use upstream HTSlib expected output files as golden outputs where possible.
  - [x] Use upstream faidx FASTA/FASTQ `.fai`, `.fa`, and `.fq` expected outputs.
  - [x] Use upstream FASTQ conversion `.sam`, `.fq`, and `.fa` expected outputs.
  - [x] Use upstream tabix BED, GFF, VCF, large-coordinate, and separate-region `.out` files.
  - [x] Use upstream long-reference VCF CSI expected outputs for `longrefs/index.expected1.vcf` and `longrefs/index.expected2.vcf`.
  - [x] Use upstream view outputs for `range.out`, `range.out2`, `modhdr.expected.vcf`, and `tabix/vcf_file.vcf`.
  - [x] Use upstream `tabix/vcf_file.vcf` as an exact VCF write golden output.
  - [x] Use upstream `tabix/vcf_file.vcf` records as BCF write/decode golden output.
  - [x] Use upstream VCF various expected outputs for `test-vcf-hdr.out`, `formatcols.vcf`, `noroundtrip-out.vcf`, `formatmissing-out.vcf`, `vcf_meta_meta.vcf`, and `vcf44_1.expected`.
  - [x] Use upstream base-modification, SAM filter, mpileup SAM, synced-reader, VCF API, and BCF translation expected outputs.
- [x] Add temporary differential tests against C HTSlib for unclear edge cases.
  - [x] Replace temporary C-oracle coverage with direct Rust assertions and upstream golden fixtures where the relevant edge cases are now understood.

## Phase 3: Format-Neutral I/O and Detection

- [x] Implement initial HTS-style format detection using Rust readers and noodles parsers.
- [x] Support local files first.
- [x] Keep remote I/O out of scope for now.
- [x] Port `hfile` behavior needed by tests.
  - [x] Port `haddextension` path/URL extension behavior.
  - [x] Port `data:` URL decoding for plain, empty, percent-encoded, and base64 payloads.
- [x] Port logging behavior used by `test-logging.pl`.
- [x] Implement error types that preserve enough context for HTSlib-compatible tests.
  - [x] Add a reusable contextual error type with operation and optional path metadata.
  - [x] Preserve context through Rust-only `io::Result` APIs without changing existing signatures.
  - [x] Add context coverage for local format detection and CRAM FASTA-reference setup failures.

## Phase 4: BGZF, Compression, and Random Access

- [x] Use `noodles-bgzf` for BGZF read/write.
- [x] Implement HTSlib-compatible virtual offset behavior on top of `noodles-bgzf::VirtualPosition`.
- [x] Port BGZF block boundary tests.
- [x] Port `.gzi` index read/write behavior.
- [x] Leave `bgzip` CLI behavior out of scope for now.
- [x] Decide threading behavior for BGZF compression and decompression.
  - [x] Use noodles multithreaded BGZF reader/writer worker-count APIs for the Rust-only target.

## Phase 5: Indexes and Region Queries

- [x] Use `noodles-bam` for BAI.
- [x] Use `noodles-csi` for CSI.
- [x] Use `noodles-tabix` for TBI.
- [x] Add initial BAI/TBI/CSI index build helpers.
- [x] Implement HTSlib-compatible index lookup order and fallback behavior.
  - [x] Add local associated-index candidate ordering and `##idx##` delimiter handling.
  - [x] Wire lookup order into local BAM, VCF, BCF, and VCF text query helpers.
  - [x] Add CRAM/CRAI lookup wiring when CRAM indexed iteration is ported.
- [x] Port region parsing from `hts_parse_reg` and related tests.
- [x] Port iterators for BAM, CRAM, VCF, BCF, and tabix-backed text formats.
  - [x] Add initial tabix-backed BED text region queries.
  - [x] Add BAM, CRAM, VCF, and BCF indexed iterators.
    - [x] Add initial BAM indexed query helpers using associated BAI/CSI indexes.
    - [x] Add initial CRAM indexed query helpers using associated CRAI indexes and local FASTA references.
    - [x] Add initial VCF indexed query helpers using associated TBI/CSI indexes.
    - [x] Add initial BCF indexed query helpers using associated CSI indexes.
    - [x] Expose Rust owning iterator adapters for BAM, CRAM, VCF, and BCF indexed query results.
    - [x] Add VCF/GFF tabix query parity.
    - [x] Add CSI-backed VCF queries for large-coordinate fixtures.
    - [x] Add adaptive-depth VCF CSI queries for long-reference fixtures, including `INFO/END` overlap.
- [x] Validate behavior against `test_index.c`, `test-parse-reg.c`, and tabix tests.
  - [x] Add initial `test_index.c` build coverage for BAM BAI, VCF TBI/CSI, BCF CSI, and explicit VCF CSI min_shift.
  - [x] Add local associated-index fallback coverage for replaced extensions and explicit `##idx##` paths.
  - [x] Add SAM.gz CSI build and indexed query coverage.
  - [x] Add explicit CSI min_shift build coverage for BAM, SAM.gz, VCF, and BCF.
  - [x] Add SAM.gz BAI build and indexed query coverage.
  - [x] Add CRAI read/write round-trip coverage.
  - [x] Port CRAM/CRAI.

## Phase 6: SAM, BAM, and CRAM

- [x] Use `noodles-sam` for SAM headers and records.
- [x] Use `noodles-bam` for BAM I/O.
  - [x] Read BAM headers and records with `noodles-bam`.
  - [x] Add BAM write/query parity where required by ported tests.
    - [x] Add initial BAM region query parity against `range.bam`.
    - [x] Add BAM write/decode record parity against `range.bam`.
- [x] Use `noodles-cram` for CRAM I/O.
  - [x] Read CRAM headers with `noodles-cram`.
  - [x] Add CRAM record iteration with reference-repository/cache behavior.
    - [x] Add CRAM record iteration for embedded-reference TLEN fixtures.
    - [x] Add CRAM indexed query coverage with an explicit local FASTA reference repository.
  - [x] Add CRAM write/decode record-summary parity against `range.cram` with a local FASTA reference repository.
- [x] Build adapter types for HTSlib-style record mutation, header mutation, and typed auxiliary fields.
  - [x] Add a SAM record adapter for HTSlib-style core field lookup and mutation.
  - [x] Add a typed SAM auxiliary-field adapter for `RecordBuf` get/insert/remove operations.
  - [x] Add a SAM header adapter for HTSlib-style `@SQ` reference sequence lookup, insert, replace, and removal.
- [x] Port SAM/BAM/CRAM read/write/view tests.
  - [x] Add SAM-to-BAM write/view parity for `xx#minimal.sam`.
  - [x] Add padded BAM raw-header read/view parity for `ce#1.sam`.
  - [x] Add SAM-to-CRAM write/decode parity for `ce#1.sam`, `ce#2.sam`, and `ce#1000.sam` with `ce.fa`, plus `xx#pair.sam` with `xx.fa`.
  - [x] Add CRAM-to-CRAM write/decode summary parity for `range.cram`.
- [x] Port base modification behavior from `sam_mods.c`.
  - [x] Add initial MM tag position/status reporting for `MM-variants.sam`.
  - [x] Add ML probability reporting for explicit and double-strand MM fixtures.
  - [x] Add ChEBI ID and canonical `N` base-modification reporting for `MM-chebi.sam`.
  - [x] Add multiple-modification, missing-MM reset, MN validation, and MM bounds rejection behavior.
  - [x] Add extended metadata and pileup reporting for upstream base-modification fixtures.
- [x] Port pileup behavior if it is in scope.
  - [x] Add focused SAM base-modification pileup output for `pileup_mod` fixtures.
  - [x] Add HTSlib `test/pileup.c`-style SAM output for deletion, insertion, insertion/deletion, refskip, pad, and overlap-removal fixtures.
  - [x] Add HTSlib `test/pileup.c`-style BAM output for the upstream `small.bam` edge fixture.
- [x] Port realignment/probabilistic alignment behavior if it is in scope.
  - [x] Add an initial pure Rust `probaln_glocal` port with HTSlib-oracle parity for exact, mismatch, and insertion cases.
  - [x] Add default, apply-on-recalculation, and forced non-extended BAQ recalculation parity for `test_realn.c` fixtures.
  - [x] Add extended BAQ recalculation parity for `test_realn.c` fixtures.
- [x] Validate TLEN, CIGAR, MD/NM, padded reference, and base-modification edge cases.
  - [x] Validate CRAM TLEN auto-creation for pair and triplet ordering fixtures.
  - [x] Validate CIGAR-derived SAM filter metrics and string matching for upstream `sam_filter` fixtures.
  - [x] Validate pileup CIGAR handling for deletion, insertion, insertion/deletion, refskip, and pad fixtures.
  - [x] Validate CRAM view MD/NM regeneration against upstream `range.out` and `range.out2`.
  - [x] Validate MM explicit/implicit base-modification status handling for `MM-variants.sam`.
  - [x] Validate ML probability formatting for `MM-explicit.sam` and `MM-double.sam`.
  - [x] Validate ChEBI ID and `N` canonical-base handling for `MM-chebi.sam`.
  - [x] Validate base-modification pileup output for `MM-pileup.sam`, `MM-pileup2.sam`, and `MM-MNp.sam`.
  - [x] Validate multiple-modification, missing-MM reset, MN validation, and MM bounds cases.

## Phase 7: VCF and BCF

- [x] Use `noodles-vcf` for VCF headers, records, parsing, and writing.
  - [x] Read VCF headers and records with `noodles-vcf`.
  - [x] Add VCF write/query parity where required by ported tests.
    - [x] Add initial VCF region query parity against tabix fixtures.
    - [x] Add exact VCF write parity against the upstream tabix VCF fixture.
- [x] Use `noodles-bcf` for BCF parsing and writing.
  - [x] Read BCF headers and records with `noodles-bcf`.
  - [x] Add BCF write/query parity where required by ported tests.
    - [x] Add initial BCF CSI build and region query parity against tabix fixtures.
    - [x] Add BCF write/decode record parity for upstream VCF records with large genotype allele indexes.
- [x] Build adapter types for HTSlib-style header IDs, typed FORMAT/INFO values, genotypes, and record mutation.
  - [x] Add `VcfHeaderId` for HTSlib-style FILTER, INFO, FORMAT, contig, and structured generic header IDs.
  - [x] Add `VcfRecordAdapter` for typed INFO/FORMAT integer access, genotype access, allele mutation, INFO/FORMAT mutation, and allele removal.
- [x] Port VCF header API tests.
  - [x] Add header get/remove coverage for FILTER, INFO, FORMAT, contig, and generic records.
- [x] Port VCF/BCF translation behavior.
  - [x] Port synthetic header merge and record translation fixture from `test-bcf-translate.c`.
- [x] Port synced BCF reader behavior or define a Rust-native equivalent.
  - [x] Add initial no-index VCF synchronization summary behavior.
  - [x] Add weird chromosome-name region and target filtering behavior for synced-reader fixtures.
  - [x] Add Rust-native allele-pairing logic equivalent for synced variant groups.
- [x] Port variant classification behavior from `bcf_set_variant_type`.
- [x] Port sweep behavior from `vcf_sweep`.

## Phase 8: FASTA, FASTQ, and Reference Access

- [x] Use `noodles-fasta` for FASTA records and FAI indexes.
- [x] Use `noodles-fastq` for FASTQ records and FQI indexes where applicable.
- [x] Port `faidx` query behavior and exact output formatting.
  - [x] Add HTSlib-compatible uncompressed FASTA/FASTQ indexing and 50-column retrieval output.
  - [x] Add bgzip-compressed FASTA retrieval through FAI/GZI.
- [x] Port FASTQ conversion tests.
  - [x] Add FASTQ/FASTA to SAM and SAM to FASTQ/FASTA golden-output coverage.
  - [x] Add CASAVA/filter-specific conversion coverage.
- [x] Use noodles CRAM reference-cache/reference-repository behavior where it exists.
  - [x] Build CRAM reference repositories with `noodles_fasta::Repository` over indexed local FASTA files.

## Phase 9: Utilities and HTSlib-Specific Behavior

- [x] Replace C container utilities with Rust standard library equivalents where they are internal only.
- [x] Port public or test-observable utility behavior.
  - [x] `kbitset`
  - [x] `kstring`
  - [x] `khash`
  - [x] `khash_str2int`
  - [x] `klist`
  - [x] `kroundup`
  - [x] `ksort`
  - [x] `hts_os` rand48 helpers
  - [x] `kfunc`
  - [x] integer parsing and string utilities
  - [x] endian helpers
  - [x] time helpers
- [x] Port region index behavior from `regidx.c`.
- [x] Port expression parsing and evaluation from `hts_expr.c`.
- [x] Keep remote/plugin I/O out of scope for now.
- [x] Keep shared-library export and C ABI introspection out of scope for now.

## Phase 10: CLIs

- [x] Keep command-line tool parity out of scope for now.
- [x] Do not implement `bgzip`, `htsfile`, `tabix`, or `annot-tsv` yet.
- [x] Revisit CLI parity only after the Rust library test-parity target is stable.
  - [x] Keep CLI parity deferred after the Rust library test-parity pass; command-line tools remain out of scope for the current no-CLI target.

## Phase 11: Completion Criteria

- [x] All selected HTSlib tests are ported or explicitly marked out of scope.
- [x] All ported tests pass under `cargo test`.
- [x] Differential tests against C HTSlib either pass or are replaced by direct Rust assertions.
- [x] The coverage matrix shows no unknown public APIs.
- [x] Documentation explains what is compatible with HTSlib C and what is intentionally Rust-native.

## bcftools-rs Downstream Gap Rollup

These items are referenced from
`../docs/subcommand-coverage.md`. They are not required for the completed
Rust-only HTSlib target, but they are the known extension points needed by the
bcftools-rs command ports.

- [ ] `synced_bcf_reader` full API parity for bcftools: multi-input streaming,
  region/target restriction, collapse modes, per-reader allele translation,
  and command-shaped diagnostics.
- [x] `bcf_translate` coverage beyond the synthetic translation fixture,
  including merged-header to per-input translation tables for `merge`,
  `concat`, `isec`, and plugins.
- [x] Complete `bcf_update_*` mutation primitives for INFO, FORMAT, FILTER,
  ID, QUAL, POS, alleles, and vector trimming/remapping across all bcftools
  call sites.
- [x] Pileup iterator surface for bcftools `mpileup`, including multi-input
  synchronized pileup behavior.
- [x] BAQ and `probaln_glocal` wiring for `bam2bcf*.c` call sites.
- [ ] `hts_set_threads`/BGZF writer thread-pool support for VCF/BCF writers
  used by `view`, `merge`, `norm`, `concat`, and `sort`.
- [x] Region-with-target arithmetic parity for bcftools `-r`/`-R` versus
  `-t`/`-T`, including streaming target filtering and overlap modes.

## Answered Questions and Decisions

- API target: Rust-only API for now.
- CLI tools: out of scope for now. Do not implement `bgzip`, `htsfile`, `tabix`, or `annot-tsv` yet.
- Remote I/O: out of scope for now.
- CRAM reference cache: use whatever noodles provides, if available.
- Deprecated APIs: keep them for now because downstream C libraries may rely on unknown parts of the HTSlib API surface.
- Compatibility target: test parity for now.
