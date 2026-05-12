# HTSlib API Coverage

This is the coverage map for the pure Rust port. Public headers now have item-level classification, while implementation work remains for APIs marked as requiring new Rust code or deferred adapters.

Status meanings:

- `covered by noodles`: noodles already has the main format implementation.
- `adapter needed`: noodles covers core behavior, but HTSlib-compatible API semantics need wrapper code.
- `new Rust needed`: no direct noodles equivalent was identified.
- `out of scope now`: explicitly deferred by project decision.

## Classification Snapshot

This is a header-level summary of the public HTSlib API surface. Function,
type, constant, and macro decisions are tracked in the item-level section below.

The extraction aid for that follow-up work is
[`public-api-inventory.md`](public-api-inventory.md), generated from the public
headers under `htslib/htslib`.

| Classification | Public API Areas |
| --- | --- |
| Covered directly by noodles | Core file-format parsers, writers, and indexes from `bgzf.h`, `cram.h`, `faidx.h`, `kseq.h`, `sam.h`, `tbx.h`, and `vcf.h`, via the `noodles-bgzf`, `noodles-cram`, `noodles-fasta`, `noodles-fastq`, `noodles-sam`, `noodles-bam`, `noodles-tabix`, `noodles-csi`, `noodles-vcf`, and `noodles-bcf` crates. |
| Covered by small adapters over noodles | HTSlib-compatible virtual offsets, GZI, index lookup order, region parsing, tabix query output, CRAM reference repositories, SAM/BAM/CRAM view helpers, VCF/BCF header/value helpers, BCF iterators, and sweep helpers. |
| Requires new Rust implementation | HTSlib-specific utility behavior for `hfile.h`, `hts_defs.h`, `hts_endian.h`, `hts_expr.h`, `hts_log.h`, `hts_os.h`, `kbitset.h`, `kfunc.h`, `khash.h`, `khash_str2int.h`, `klist.h`, `kroundup.h`, `ksort.h`, `kstring.h`, `regidx.h`, `synced_bcf_reader.h`, and `vcfutils.h`. |
| Deprecated but retained for compatibility | Deprecated public APIs remain tracked instead of being dropped, including `bgzf_is_bgzf`, `faidx_fetch_nseq`, `hts_detect_format2`, `knetfile.h`, `bcf_hdr_combine`, `bcf_hdr_fmt_text`, `bcf_is_snp`, `bcf_is_mnp`, `bcf_remove_alleles`, and legacy collapse constants. |
| Intentionally out of scope now | Command-line tools, remote I/O, plugin-backed hFILE backends, C ABI/shared-library introspection, and C-style shared thread-pool objects. |

| Header | HTSlib Area | Initial Status | Rust Direction |
| --- | --- | --- | --- |
| `bgzf.h` | BGZF I/O, virtual offsets, GZI, threading | adapter needed | Use `noodles-bgzf`; HTSlib-compatible virtual-offset, local detection, GZI, and Rust worker-count threading adapters exist. |
| `cram.h` | CRAM I/O and reference handling | adapter needed | Use `noodles-cram`; use noodles reference/cache behavior where available. |
| `faidx.h` | FASTA/FASTQ indexing and sequence retrieval | adapter needed | Use `noodles-fasta` and `noodles-fastq`; add HTSlib-compatible query/output behavior. |
| `hfile.h` | HTSlib I/O abstraction | new Rust needed | Implement local-file behavior only for now; remote I/O is out of scope. |
| `hts.h` | Format-neutral I/O, detection, indexes, iterators, regions | adapter needed | Initial local format detection exists; indexes and iterators should wrap noodles where possible. |
| `hts_defs.h` | Common exported definitions | new Rust needed | Preserve compatibility concepts only where visible in Rust API/tests. |
| `hts_endian.h` | Endian helpers | new Rust needed | Use Rust endian primitives; add compatibility functions if tests require them. |
| `hts_expr.h` | Filter expression parser/evaluator | implemented | Rust parser/evaluator with acceptance coverage from `test_expr.c`. |
| `hts_log.h` | Logging | implemented | Source-message style validation behavior from `test-logging.pl`. |
| `hts_os.h` | OS compatibility helpers | new Rust needed | Implement only behavior needed by Rust tests. |
| `kbitset.h` | Bitset utility | new Rust needed | Use Rust containers internally; preserve public/test-observable behavior. |
| `kfunc.h` | Numeric/statistical helpers | new Rust needed | Port behavior covered by `test_kfunc.c`. |
| `khash.h` | Hash table macros | new Rust needed | Use Rust collections internally; retain compatibility behavior where downstream tests need it. |
| `khash_str2int.h` | String-to-integer hash helpers | new Rust needed | Use Rust maps; add compatibility behavior when test coverage requires it. |
| `klist.h` | List macros | new Rust needed | Use Rust collections internally. |
| `knetfile.h` | Legacy network I/O | out of scope now | Remote I/O is out of scope. Deprecated API remains tracked for future compatibility. |
| `kroundup.h` | Allocation/rounding helpers | new Rust needed | Use Rust allocation patterns; add helpers only if externally visible. |
| `kseq.h` | FASTA/FASTQ parser macros | adapter needed | Prefer `noodles-fasta` and `noodles-fastq`. |
| `ksort.h` | Sorting macros | new Rust needed | Use Rust sorting internally; expose compatibility behavior only if needed. |
| `kstring.h` | Dynamic string utility | new Rust needed | Port test-observable behavior from `test_kstring.c`. |
| `regidx.h` | Region indexes | new Rust needed | Port region index behavior and `test-regidx.c`. |
| `sam.h` | SAM/BAM/CRAM headers, records, pileup, iterators | adapter needed | Use `noodles-sam`, `noodles-bam`, and `noodles-cram`; add HTSlib-compatible mutation and query semantics. |
| `synced_bcf_reader.h` | Multi-file VCF/BCF reader | adapter needed | Use noodles VCF/BCF readers; implement HTSlib-compatible synchronization semantics. |
| `tbx.h` | Tabix index/query | adapter needed | Use `noodles-tabix` and `noodles-csi`. |
| `thread_pool.h` | HTSlib thread pool | out of scope now | C-style shared thread-pool objects are out of scope for the Rust-only API; BGZF uses noodles worker-count APIs. |
| `vcf.h` | VCF/BCF headers, records, typed values | adapter needed | Use `noodles-vcf` and `noodles-bcf`; add HTSlib-compatible header IDs, typed values, and record mutation. |
| `vcf_sweep.h` | VCF sweep helper | implemented | Rust forward/backward sweep adapter covered by `test-vcf-sweep.c` parity cases. |
| `vcfutils.h` | VCF utilities | new Rust needed | Port behavior covered by VCF utility tests. |

## Item-Level Classification

This section tracks headers whose public items have been refined beyond the
header-level map. It is intentionally incomplete until every public header has
an item-level decision.

### `hts_endian.h`

| Public item(s) | Classification | Rust evidence |
| --- | --- | --- |
| `HTS_ENDIAN_H` | out of scope now | C include guard; no Rust equivalent needed. |
| `HTS_x86`, `HTS_LITTLE_ENDIAN`, `HTS_BIG_ENDIAN`, `HTS_ENDIAN_NEUTRAL`, `HTS_ALLOW_UNALIGNED` | out of scope now | Compile-time C optimization switches; Rust helpers use portable endian primitives and slice reads/writes. |
| `uint16_u`, `uint32_u`, `uint64_u` | out of scope now | C-only unaligned-load typedefs; Rust uses safe byte-slice conversion. |
| `le_to_u8`, `le_to_u16`, `le_to_u32`, `le_to_u64` | new Rust implementation | `htslib_rs::endian::{le_to_u8, le_to_u16, le_to_u32, le_to_u64}` with `endian.rs` acceptance coverage. |
| `u16_to_le`, `u32_to_le`, `u64_to_le` | new Rust implementation | `htslib_rs::endian::{u16_to_le, u32_to_le, u64_to_le}` with aligned and unaligned-slice tests. |
| `le_to_i8`, `le_to_i16`, `le_to_i32`, `le_to_i64` | new Rust implementation | `htslib_rs::endian::{le_to_i8, le_to_i16, le_to_i32, le_to_i64}` with two's-complement boundary tests. |
| `i16_to_le`, `i32_to_le`, `i64_to_le` | new Rust implementation | `htslib_rs::endian::{i16_to_le, i32_to_le, i64_to_le}` with boundary tests. |
| `le_to_float`, `le_to_double`, `float_to_le`, `double_to_le` | new Rust implementation | C-compatible aliases over `le_to_f32`, `le_to_f64`, `f32_to_le`, and `f64_to_le`, covered by float and double endian tests. |

### `hts_defs.h`

| Public item(s) | Classification | Rust evidence |
| --- | --- | --- |
| `HTSLIB_HTS_DEFS_H` | out of scope now | C include guard; no Rust equivalent needed. |
| `HTS_COMPILER_HAS`, `HTS_GCC_AT_LEAST` | out of scope now | C preprocessor feature probes; Rust uses `cfg` and the compiler's native feature model. |
| `HTS_NONSTRING`, `HTS_NORETURN`, `HTS_ACCESS`, `HTS_OPT3`, `HTS_ALIGN32`, `HTS_RESULT_USED`, `HTS_UNUSED` | out of scope now | C compiler attributes; Rust equivalents are native language attributes or lint behavior, not HTSlib API. |
| `HTS_DEPRECATED`, `HTS_DEPRECATED_ENUM` | deprecated but retained for compatibility | Deprecated C APIs remain tracked in the inventory and per-header coverage tables; the C attribute macro itself has no Rust API equivalent. |
| `HTS_PRINTF_FMT`, `HTS_FORMAT` | out of scope now | C `printf` format-checking attributes; Rust formatting is type checked by macros and does not use these declarations. |
| `HTSLIB_EXPORT` | out of scope now | C symbol-visibility/export macro; the current target is Rust-only API, not a C ABI. |
| `hts_prefetch` | new Rust implementation | `htslib_rs::hts_defs::hts_prefetch` preserves the optimization-hint call shape with no semantic side effects, covered by `hts_defs.rs` acceptance coverage. |

### `hts_log.h`

| Public item(s) | Classification | Rust evidence |
| --- | --- | --- |
| `HTS_LOG_H` | out of scope now | C include guard; no Rust equivalent needed. |
| `enum htsLogLevel`, `HTS_LOG_OFF`, `HTS_LOG_ERROR`, `HTS_LOG_WARNING`, `HTS_LOG_INFO`, `HTS_LOG_DEBUG`, `HTS_LOG_TRACE` | new Rust implementation | `htslib_rs::log_compat::LogLevel` preserves HTSlib integer values, including the gap before warning. |
| `hts_set_log_level`, `hts_get_log_level` | new Rust implementation | `htslib_rs::log_compat::{hts_set_log_level, hts_get_log_level}` backed by process-wide atomic state and covered by logging tests. |
| `hts_verbose` | new Rust implementation | `htslib_rs::log_compat::{hts_verbose, set_hts_verbose}` expose the integer state needed by the Rust-only API target. |
| `hts_log` | out of scope now | C variadic formatting API is not exposed directly in the Rust-only API; Rust logging call sites should use Rust formatting. |
| `hts_log_error`, `hts_log_warning`, `hts_log_info`, `hts_log_debug`, `hts_log_trace` | out of scope now | C `__func__` variadic convenience macros are not exposed directly; source-message style is validated by `check_htslib_log_message_style`. |

### `hts_os.h`

| Public item(s) | Classification | Rust evidence |
| --- | --- | --- |
| `HTSLIB_HTS_OS_H` | out of scope now | C include guard; no Rust equivalent needed. |
| `hts_srand48` | new Rust implementation | `htslib_rs::hts_os::hts_srand48` seeds the process-wide deterministic POSIX-compatible `rand48` state. |
| `hts_erand48` | new Rust implementation | `htslib_rs::hts_os::hts_erand48` advances caller-provided seed words and returns the deterministic floating-point output. |
| `hts_drand48` | new Rust implementation | `htslib_rs::hts_os::hts_drand48` advances the global state and returns deterministic floating-point output. |
| `hts_lrand48` | new Rust implementation | `htslib_rs::hts_os::hts_lrand48` advances the global state and returns the high 31 bits. |
| `srand48`, `erand48`, `drand48`, `lrand48` Windows macro redirects | out of scope now | C platform compatibility macros; Rust callers use the explicit `hts_*rand48` helpers. |
| disabled `is_cygpty` declaration | out of scope now | Disabled C-only compatibility declaration, not part of the Rust target. |
| MinGW `mkdir` macro, Windows `srandom`/`random` aliases, MSVC `ssize_t` fallback | out of scope now | C compiler and libc compatibility shims; no Rust equivalent needed. |

### `kroundup.h`

| Public item(s) | Classification | Rust evidence |
| --- | --- | --- |
| `KROUNDUP_H` | out of scope now | C include guard; no Rust equivalent needed. |
| `k_signed_type`, `k_high_bit_set` | out of scope now | C macro internals for type-generic overflow detection; the Rust implementation uses typed checked power-of-two operations. |
| `kroundup64` | new Rust implementation | `htslib_rs::kroundup::kroundup64` and `roundup_u64` round to the next power of two and saturate to `u64::MAX`, covered by boundary tests. |
| `kroundup32` | new Rust implementation | `htslib_rs::kroundup::kroundup32` and `roundup_u32` preserve the historical 32-bit macro behavior, covered by boundary tests. |
| `kroundup_size_t` | new Rust implementation | `htslib_rs::kroundup::kroundup_size_t` and `roundup_size_t` preserve `size_t`-style saturation for Rust `usize`, covered by boundary tests. |

### `kfunc.h`

| Public item(s) | Classification | Rust evidence |
| --- | --- | --- |
| `HTSLIB_KFUNC_H` | out of scope now | C include guard; no Rust equivalent needed. |
| `kf_lgamma` | new Rust implementation | `htslib_rs::math::kf_lgamma` exposes the log-gamma helper used internally by the Fisher exact implementation. |
| `kf_erfc` | requires new Rust implementation | Classified as public `kfunc.h` API; no Rust adapter is exposed yet because current ported HTSlib tests do not exercise it. |
| `kf_gammap`, `kf_gammaq` | requires new Rust implementation | Classified as public incomplete-gamma API; no Rust adapter is exposed yet because current ported HTSlib tests do not exercise it. |
| `kf_betai` | requires new Rust implementation | Classified as public incomplete-beta API; no Rust adapter is exposed yet because current ported HTSlib tests do not exercise it. |
| `kt_fisher_exact` | new Rust implementation | `htslib_rs::math::kt_fisher_exact` and `fisher_exact` preserve HTSlib probability, left-tail, right-tail, two-tail, and underflow behavior covered by `test_kfunc.c` cases. |

### `kbitset.h`

| Public item(s) | Classification | Rust evidence |
| --- | --- | --- |
| `KBITSET_H` | out of scope now | C include guard; no Rust equivalent needed. |
| `KBS_ELTBITS`, `KBS_ELT`, `KBS_MASK` | new Rust implementation | `htslib_rs::kbitset::{KBS_ELTBITS, kbs_elt, kbs_mask}` preserve the public bit-index arithmetic in Rust-callable form. |
| `struct kbitset_t` | new Rust implementation | `htslib_rs::kbitset::BitSet` is the Rust-owned replacement over safe storage. |
| `kbs_last_mask` | new Rust implementation | `htslib_rs::kbitset::kbs_last_mask` preserves the final-word mask calculation, including word-boundary behavior. |
| `kbs_init2`, `kbs_init` | new Rust implementation | `htslib_rs::kbitset::{kbs_init2, kbs_init}` expose empty and filled initialization over `BitSet`. |
| `kbs_resize2`, `kbs_resize` | new Rust implementation | `htslib_rs::kbitset::{kbs_resize2, kbs_resize}` preserve resize behavior and optional fill of newly added indexes. |
| `kbs_destroy` | out of scope now | Explicit C heap free has no Rust API equivalent; `BitSet` is dropped by ownership. |
| `kbs_clear`, `kbs_insert_all`, `kbs_insert`, `kbs_delete`, `kbs_exists` | new Rust implementation | `htslib_rs::kbitset` exposes C-shaped operation names backed by safe `BitSet` methods and acceptance coverage. |
| `struct kbitset_iter_t`, `kbs_start`, `kbs_next` | new Rust implementation | `htslib_rs::kbitset::BitSetIter` with `kbs_start` and `kbs_next` preserves ascending iteration state for Rust callers. |

### `klist.h`

| Public item(s) | Classification | Rust evidence |
| --- | --- | --- |
| `_AC_KLIST_H` | out of scope now | C include guard; no Rust equivalent needed. |
| `klib_unused` | out of scope now | C compiler attribute shim; Rust uses native lint and attribute behavior. |
| `KMEMPOOL_INIT2`, `KMEMPOOL_INIT`, `kmempool_t`, `kmp_init`, `kmp_destroy`, `kmp_alloc`, `kmp_free` | out of scope now | C node memory-pool implementation detail; the Rust replacement uses owned `VecDeque` storage and normal drop semantics. |
| `KLIST_INIT2`, `KLIST_INIT`, `klist_t` | new Rust implementation | `htslib_rs::klist::KList<T>` is the Rust-owned generic list replacement. |
| generated `kl_init_*`, `kl_init` | new Rust implementation | `htslib_rs::klist::kl_init` initializes an empty `KList<T>`. |
| generated `kl_destroy_*`, `kl_destroy` | new Rust implementation | `htslib_rs::klist::kl_destroy` consumes the owned list, relying on Rust drop semantics. |
| generated `kl_pushp_*`, `kl_pushp` | new Rust implementation | `htslib_rs::klist::kl_pushp` appends a default value and returns a mutable reference for fill-in behavior. |
| generated `kl_shift_*`, `kl_shift` | new Rust implementation | `htslib_rs::klist::kl_shift` preserves FIFO removal semantics and returns `None` for empty lists. |
| `kliter_t`, `kl_val`, `kl_next`, `kl_begin`, `kl_end` | new Rust implementation | `htslib_rs::klist::kl_begin` returns a Rust iterator over `KList<T>`; direct pointer sentinel macros are represented by standard iterator termination. |

### `ksort.h`

| Public item(s) | Classification | Rust evidence |
| --- | --- | --- |
| `AC_KSORT_H` | out of scope now | C include guard; no Rust equivalent needed. |
| `klib_unused`, `extern "C"` guards, `HTSLIB_EXPORT extern double hts_drand48(void)` | out of scope now | C compiler/linkage plumbing; Rust shuffle helpers accept an explicit random source and `hts_os::hts_drand48` exists separately. |
| `ks_isort_stack_t` | out of scope now | C introsort implementation detail; Rust delegates to standard library slice sorting. |
| `KSORT_SWAP` | new Rust implementation | `htslib_rs::ksort::ksort_swap` exposes checked slice-index swapping for Rust callers. |
| `KSORT_INIT`, `KSORT_INIT_STATIC`, `KSORT_INIT2`, `KSORT_INIT_`, `KSORT_INIT_GENERIC`, `KSORT_INIT_STR`, `KSORT_INIT_STATIC_GENERIC`, `KSORT_INIT_STATIC_STR`, `KSORT_INIT2_GENERIC`, `KSORT_INIT2_STR` | new Rust implementation | C code-generation macros are represented by generic Rust functions over slices and less-than closures, with integer and string acceptance coverage. |
| generated `ks_mergesort_*`, `ks_mergesort` | new Rust implementation | `htslib_rs::ksort::{mergesort_by, ks_mergesort_by}` provide stable sorting over Rust slices. |
| generated `ks_introsort_*`, `ks_introsort` | new Rust implementation | `htslib_rs::ksort::{introsort_by, ks_introsort_by}` provide unstable sorting over Rust slices. |
| generated `ks_combsort_*`, `ks_combsort` | new Rust implementation | `htslib_rs::ksort::{combsort_by, ks_combsort_by}` preserve the observable sorted-output contract. |
| generated `ks_heapsort_*`, `ks_heapsort` | new Rust implementation | `htslib_rs::ksort::{heapsort_by, ks_heapsort_by}` provide heap-sort behavior over Rust slices. |
| generated `ks_heapmake_*`, `ks_heapmake` | new Rust implementation | `htslib_rs::ksort::{heapmake_by, ks_heapmake_by}` build a max heap according to a caller-provided less-than predicate. |
| generated `ks_heapadjust_*`, `ks_heapadjust` | new Rust implementation | `htslib_rs::ksort::{heapadjust_by, ks_heapadjust_by}` repair a heap prefix and no-op for invalid Rust bounds. |
| generated `ks_ksmall_*`, `ks_ksmall` | new Rust implementation | `htslib_rs::ksort::{ksmall_by, ks_ksmall_by}` return the kth smallest value while permitting input reordering; out-of-range Rust calls return `None`. |
| generated `ks_shuffle_*`, `ks_shuffle` | new Rust implementation | `htslib_rs::ksort::{shuffle_by, ks_shuffle_by}` preserve Fisher-Yates swap behavior with an explicit unit-random callback. |
| `ks_lt_generic`, `ks_lt_str`, `ksstr_t` | new Rust implementation | `htslib_rs::ksort::{ks_lt_generic, ks_lt_str}` cover numeric and string comparator behavior; Rust `&str` replaces the C string typedef. |

### `khash.h`

| Public item(s) | Classification | Rust evidence |
| --- | --- | --- |
| `__AC_KHASH_H` | out of scope now | C include guard; no Rust equivalent needed. |
| `AC_VERSION_KHASH_H` | out of scope now | C header-version macro; compatibility is tracked by this coverage table and tests. |
| `khint32_t`, `khint64_t`, `khint_t`, `khiter_t` | out of scope now | C storage/index typedefs; Rust exposes safe key lookup and iteration instead of bucket indexes. |
| `kh_inline`, `klib_unused`, allocation macros `kcalloc`, `kmalloc`, `krealloc`, `kfree`, and internal flag macros `__ac_*` | out of scope now | C implementation and allocator details; Rust uses `HashMap` storage and ownership. |
| `__ac_HASH_UPPER`, `__KHASH_TYPE`, `__KHASH_PROTOTYPES`, `__KHASH_IMPL`, `KHASH_DECLARE`, `KHASH_INIT2`, `KHASH_INIT` | new Rust implementation | C macro-generated hash tables are represented by typed Rust wrappers over `HashMap`; the port currently exposes `StrIntMap` for tested string-to-integer behavior. |
| hash/equality helpers `kh_int_hash_func`, `kh_int_hash_equal`, `kh_int64_hash_func`, `kh_int64_hash_equal`, `kh_str_hash_func`, `kh_str_hash_equal`, `kh_kstr_hash_func`, `kh_kstr_hash_equal`, `kh_int_hash_func2` | out of scope now | C bucket hash implementation details; Rust `HashMap` supplies hashing and equality. |
| `khash_t(name)` | new Rust implementation | `htslib_rs::khash::StrIntMap` is the Rust-owned map type used by ported `test_khash.c` behavior. |
| `kh_init`, `kh_destroy`, `kh_clear` | new Rust implementation | `htslib_rs::khash::{kh_init_str_int, kh_destroy_str_int, kh_clear_str_int}` expose Rust-owned initialization, drop, and clear behavior. |
| `kh_put`, `kh_get`, `kh_del`, `kh_exist` | new Rust implementation | `htslib_rs::khash::{kh_put_str_int, kh_get_str_int, kh_del_str_int, kh_exist_str_int}` cover insertion, lookup, deletion, replacement, and existence checks. |
| `kh_size`, `kh_foreach`, `kh_foreach_value` | new Rust implementation | `htslib_rs::khash::{kh_size_str_int, kh_foreach_str_int}` expose size and key-value iteration; value-only iteration is represented by standard Rust iterator mapping. |
| `kh_resize`, `kh_grow_to_fit`, `kh_begin`, `kh_end`, `kh_key`, `kh_val`, `kh_value`, `kh_n_buckets`, `kh_stats` | out of scope now | Bucket-index, capacity, and probe-stat APIs are C table internals; the Rust API intentionally exposes safe map operations only for the Rust-only target. |
| `KHASH_SET_INIT_INT`, `KHASH_MAP_INIT_INT`, `KHASH_SET_INIT_INT64`, `KHASH_MAP_INIT_INT64`, `KHASH_SET_INIT_STR`, `KHASH_MAP_INIT_STR`, `KHASH_SET_INIT_KSTR`, `KHASH_MAP_INIT_KSTR`, `kh_cstr_t` | new Rust implementation | Macro families are represented by Rust generic collections and the tested `StrIntMap`; additional typed aliases should be added only when downstream tests require them. |

### `khash_str2int.h`

| Public item(s) | Classification | Rust evidence |
| --- | --- | --- |
| `HTSLIB_KHASH_STR2INT_H` | out of scope now | C include guard; no Rust equivalent needed. |
| `KHASH_MAP_INIT_STR(str2int, int)` | new Rust implementation | `htslib_rs::khash::Str2Int` is the Rust-owned string-to-`i32` replacement. |
| `khash_str2int_init` | new Rust implementation | `htslib_rs::khash::khash_str2int_init` creates an empty map. |
| `khash_str2int_destroy`, `khash_str2int_destroy_free` | new Rust implementation | `htslib_rs::khash::{khash_str2int_destroy, khash_str2int_destroy_free}` consume the owned map; Rust ownership drops keys and values. |
| `khash_str2int_has_key` | new Rust implementation | `htslib_rs::khash::khash_str2int_has_key` preserves existence-check behavior. |
| `khash_str2int_get` | new Rust implementation | `htslib_rs::khash::khash_str2int_get` returns `Some(value)` for present keys and `None` for absent keys. |
| `khash_str2int_inc` | new Rust implementation | `htslib_rs::khash::khash_str2int_inc` auto-assigns the next sequential integer and returns existing values unchanged. |
| `khash_str2int_set` | new Rust implementation | `htslib_rs::khash::khash_str2int_set` inserts or replaces explicit values and reports whether the key was new. |
| `khash_str2int_size` | new Rust implementation | `htslib_rs::khash::khash_str2int_size` returns the number of keys. |

### `kstring.h`

| Public item(s) | Classification | Rust evidence |
| --- | --- | --- |
| `KSTRING_H`, `KSTRING_T`, `KS_ATTR_PRINTF`, `HAVE___BUILTIN_CLZ`, `HTSLIB_SSIZE_T`, `HTSLIB_EOVERFLOW` | out of scope now | C include guards, typedef guards, compiler attributes, and platform shims; Rust uses native types and formatting. |
| `struct kstring_t`, `KS_INITIALIZE`, `ks_initialize` | new Rust implementation | `htslib_rs::kstring::KString` and `ks_initialize` provide an owned Rust dynamic byte string and reset behavior. |
| `ks_resize`, `ks_expand` | new Rust implementation | `htslib_rs::kstring::{ks_resize, ks_expand}` expose capacity reservation behavior over `Vec<u8>`. |
| `ks_str`, `ks_c_str`, `ks_len`, `ks_clear`, `ks_release`, `ks_free` | new Rust implementation | `htslib_rs::kstring` exposes C-shaped access, length, clear, release, and owned-drop helpers. |
| `kputsn`, `kputs`, `kputc`, `kputc_`, `kputsn_` | new Rust implementation | C-shaped append helpers are backed by `KString::{push_bytes,push_byte}` and covered by acceptance tests. |
| `kputuw`, `kputw`, `kputll`, `kputl` | new Rust implementation | Decimal integer append helpers are backed by Rust formatting and covered over boundary ranges from `test_kstring.c`. |
| `kinsert_char`, `kinsert_str` | new Rust implementation | `htslib_rs::kstring::{kinsert_char, kinsert_str}` preserve insertion-at-position behavior and out-of-bounds errors. |
| `kmemmem`, `kstrstr`, `kstrnstr` | new Rust implementation | `htslib_rs::kstring::{kmemmem, kstrstr, kstrnstr}` expose HTSlib-style memory, string, and bounded C-string search behavior. |
| `kgetline`, `kfgetline`, `kgetline2`, `kgets_func`, `kgets_func2` | adapter needed | `htslib_rs::kstring::{getline_from_chunks, kgetline_from_chunks}` covers line-termination semantics for Rust chunk providers; C callback and `FILE *` forms remain out of scope for the Rust-only API. |
| `kvsprintf`, `ksprintf` | requires new Rust implementation | C variadic formatting API is not exposed; Rust call sites should use Rust formatting or dedicated append helpers. |
| `kputd` | requires new Rust implementation | Public custom double formatter is classified but not yet exposed as a separate Rust helper. |
| `ksplit_core`, `ksplit` | requires new Rust implementation | Split-offset API is classified but not yet implemented because current ported tests do not exercise it. |
| `struct ks_tokaux_t`, `kstrtok` | requires new Rust implementation | Tokenizer state and behavior are classified but not yet implemented. |

### `vcfutils.h`

| Public item(s) | Classification | Rust evidence |
| --- | --- | --- |
| `HTSLIB_VCFUTILS_H` | out of scope now | C include guard; no Rust equivalent needed. |
| forward `struct kbitset_t` | out of scope now | C forward declaration; Rust uses `htslib_rs::kbitset::BitSet` where a bit-set abstraction is needed. |
| `bcf_trim_alleles` | requires new Rust implementation | Classified public utility for removing unused ALT alleles from genotype fields; no standalone Rust helper is exposed yet. |
| deprecated `bcf_remove_alleles` | deprecated but retained for compatibility | Deprecated mask API remains tracked; current Rust acceptance coverage focuses on the replacement allele-set behavior. |
| `bcf_remove_allele_set` | adapter needed | `variant_io_compat` contains HTSlib-style ALT allele removal logic for INFO/FORMAT `Number=A/R/G/LA/LR/LG` trimming and genotype remapping, covered by `test-vcf-api.c` parity cases; a direct Rust API alias remains to be exposed. |
| `bcf_calc_ac` | requires new Rust implementation | Allele-count calculation from INFO/AN, INFO/AC, and genotype fields is classified but not yet exposed as a Rust helper. |
| `GT_HOM_RR`, `GT_HOM_AA`, `GT_HET_RA`, `GT_HET_AA`, `GT_HAPL_R`, `GT_HAPL_A`, `GT_UNKN` | requires new Rust implementation | Genotype-type constants are classified; no public Rust enum/constant set is exposed yet. |
| `bcf_gt_type` | requires new Rust implementation | Genotype classification from BCF FORMAT/GT is classified but not yet exposed as a Rust helper. |
| `bcf_acgt2int` | new Rust implementation | `htslib_rs::variant::bcf_acgt2int` preserves A/C/G/T to integer conversion with case-insensitive input. |
| `bcf_int2acgt` | new Rust implementation | `htslib_rs::variant::bcf_int2acgt` maps valid integer codes back to A/C/G/T and rejects out-of-range indexes. |
| `bcf_ij2G` | new Rust implementation | `htslib_rs::variant::bcf_ij2g` preserves the diploid `Number=G` triangular-index mapping. |

### `vcf_sweep.h`

| Public item(s) | Classification | Rust evidence |
| --- | --- | --- |
| `HTSLIB_VCF_SWEEP_H` | out of scope now | C include guard; no Rust equivalent needed. |
| `bcf_sweep_t` | adapter needed | `htslib_rs::variant_io_compat::VariantSweep` is the Rust-owned in-memory sweep state over noodles-parsed VCF/BCF records. |
| `bcf_sweep_init` | adapter needed | `read_vcf_sweep_from_path` and `read_bcf_sweep_from_path` initialize sweep state from local files, matching the Rust-only and no-remote-I/O target. |
| `bcf_sweep_destroy` | out of scope now | Explicit C heap destruction is represented by Rust ownership and normal drop semantics. |
| `bcf_sweep_hdr` | adapter needed | `VariantSweep::header` returns the parsed noodles VCF header; acceptance coverage checks sample count from the swept header. |
| `bcf_sweep_fwd` | adapter needed | `VariantSweep::forward` preserves forward traversal over the loaded records, covered by `test-vcf-sweep.c` checksum parity. |
| `bcf_sweep_bwd` | adapter needed | `VariantSweep::backward` preserves reverse traversal from the current cursor, covered by `test-vcf-sweep.c` checksum parity. |

### `thread_pool.h`

| Public item(s) | Classification | Rust evidence |
| --- | --- | --- |
| `HTSLIB_THREAD_POOL_H` | out of scope now | C include guard; no Rust equivalent needed. |
| `hts_tpool`, `hts_tpool_process`, `hts_tpool_result` | out of scope now | C callback-driven shared thread-pool objects are outside the Rust-only API target. Rust format adapters expose worker-count choices directly where noodles supports them. |
| `hts_tpool_init`, `hts_tpool_size`, `hts_tpool_worker_id`, `hts_tpool_destroy`, `hts_tpool_kill` | out of scope now | Pool lifecycle and worker identity are C threading primitives; the Rust API does not expose a reusable HTSlib-style scheduler. |
| `hts_tpool_dispatch`, `hts_tpool_dispatch2`, `hts_tpool_dispatch3`, `hts_tpool_wake_dispatch` | out of scope now | C function-pointer job submission and nonblocking queue behavior are not part of the current Rust API surface. |
| `hts_tpool_process_init`, `hts_tpool_process_destroy`, `hts_tpool_process_attach`, `hts_tpool_process_detach`, `hts_tpool_process_ref_incr`, `hts_tpool_process_ref_decr` | out of scope now | C process-queue lifetime, scheduler attachment, and reference-count management are replaced by Rust ownership and per-adapter worker-count configuration. |
| `hts_tpool_process_flush`, `hts_tpool_process_reset`, `hts_tpool_process_qsize`, `hts_tpool_process_empty`, `hts_tpool_process_len`, `hts_tpool_process_sz`, `hts_tpool_process_shutdown`, `hts_tpool_process_is_shutdown` | out of scope now | Process-queue inspection and shutdown controls are C scheduler details; current acceptance coverage targets observable format behavior instead. |
| `hts_tpool_next_result`, `hts_tpool_next_result_wait`, `hts_tpool_delete_result`, `hts_tpool_result_data` | out of scope now | C result queue ownership and raw `void *` result access are not exposed in the Rust-only API. |
| `bgzf_thread_pool`, `fai_thread_pool`, `hts_set_thread_pool` integration points | adapter needed | Shared C pool assignment is out of scope, but equivalent BGZF read/write concurrency is covered by `read_all_with_worker_count` and `write_all_with_worker_count` over `noodles-bgzf`; FAI and generic htsFile shared-pool hooks remain deferred until a Rust API need appears. |

### `kseq.h`

| Public item(s) | Classification | Rust evidence |
| --- | --- | --- |
| `AC_KSEQ_H` | out of scope now | C include guard; no Rust equivalent needed. |
| `klib_unused` | out of scope now | C compiler attribute shim; Rust uses native lint and attribute behavior. |
| `KS_SEP_SPACE`, `KS_SEP_TAB`, `KS_SEP_LINE`, `KS_SEP_MAX` | adapter needed | Generic token-delimiter constants are C `kstream` details; line-oriented behavior used by the port is represented by Rust `BufRead` and `kstring::getline_from_chunks` coverage. |
| `kstream_t`, `__KS_TYPE`, `__KS_BASIC`, `KSTREAM_INIT2`, `KSTREAM_INIT`, `KSTREAM_DECLARE` | out of scope now | C macro-generated stream types and constructors are replaced by Rust readers and noodles parser/indexer types. |
| `ks_init`, `ks_destroy` | out of scope now | C heap lifecycle for macro-generated streams is replaced by Rust reader ownership and drop semantics. |
| `ks_err`, `ks_eof`, `ks_rewind` | adapter needed | Error, EOF, and rewind behavior is handled by Rust reader results and explicit seek/reader reconstruction in the relevant adapters. No standalone `kstream_t` API is exposed. |
| `ks_getc`, `ks_getuntil`, `ks_getuntil2` | adapter needed | Token and line reading is represented by Rust `BufRead` helpers; `kstring::kgetline_from_chunks` covers HTSlib-style line termination semantics exercised by BGZF getline tests. |
| `kseq_t`, `__KSEQ_TYPE`, `KSEQ_INIT2`, `KSEQ_INIT`, `KSEQ_DECLARE` | covered by noodles | FASTA/FASTQ parser state is represented by `noodles-fasta` and `noodles-fastq` reader/indexer types, with Rust helper APIs wrapping observable HTSlib behavior. |
| `kseq_init`, `kseq_destroy` | covered by noodles | Rust FASTA/FASTQ readers and indexers are initialized through noodles constructors and dropped by ownership. |
| `kseq_rewind` | adapter needed | Rewind behavior is represented by seeking or rebuilding Rust readers; no mutable C parser cursor is exposed. |
| `kseq_read` | adapter needed | FASTA/FASTQ reading for indexed retrieval and conversion is covered by `faidx_compat` and `fastq_compat`, using noodles where available and scoped Rust parsing for multiline conversion parity. Acceptance tests cover `test/faidx/test-faidx.sh` and `test/fastq/test-fastq.sh` golden outputs. |

### `knetfile.h`

| Public item(s) | Classification | Rust evidence |
| --- | --- | --- |
| `KNETFILE_H` | out of scope now | C include guard; no Rust equivalent needed. |
| `netread`, `netwrite`, `netclose` | out of scope now | C platform socket/file-descriptor macros; Rust uses standard `Read`, `Write`, and file APIs for local I/O. |
| `HTSLIB_SSIZE_T`, `ssize_t` fallback | out of scope now | C platform typedef shim; Rust uses native integer and `usize`/`isize` types. |
| `KNF_TYPE_LOCAL`, `KNF_TYPE_FTP`, `KNF_TYPE_HTTP` | deprecated but retained for compatibility | Legacy knet type tags remain tracked, but FTP/HTTP remote I/O is explicitly out of scope for the current Rust-only target. |
| `knetFile` | deprecated but retained for compatibility | Deprecated C struct and its FTP/HTTP fields remain inventoried; no Rust-owned equivalent is exposed because callers should use local Rust readers or future `hfile`-style adapters. |
| `knet_tell`, `knet_fileno` | deprecated but retained for compatibility | C raw-offset and file-descriptor accessors are tracked but not exposed in the Rust-only API. |
| `knet_open` | deprecated but retained for compatibility | Tracked as a legacy predecessor to `hopen`; remote FTP/HTTP opening is out of scope, and local open behavior is represented by Rust file APIs and higher-level format adapters. |
| `knet_dopen` | deprecated but retained for compatibility | C file-descriptor wrapping is tracked but not exposed; Rust local I/O uses owned readers/writers. |
| `knet_read`, `knet_seek`, `knet_close` | deprecated but retained for compatibility | Legacy read/seek/close APIs are tracked, with local behavior covered through Rust `Read`/`Seek`/drop and higher-level adapters rather than a knet compatibility type. |

### `regidx.h`

| Public item(s) | Classification | Rust evidence |
| --- | --- | --- |
| `HTSLIB_REGIDX_H` | out of scope now | C include guard; no Rust equivalent needed. |
| `REGIDX_MAX` | new Rust implementation | `htslib_rs::regidx::REGIDX_MAX` preserves the HTSlib coordinate clamp limit. |
| `regidx_t` | new Rust implementation | `htslib_rs::regidx::RegionIndex` is the Rust-owned in-memory region index. |
| `regitr_t` | new Rust implementation | Overlap and full-index iteration are represented by Rust iterators and borrowed `RegionRecord` values instead of a mutable C iterator struct. |
| `regitr_payload`, `REGITR_START`, `REGITR_END`, `REGITR_PAYLOAD`, `REGITR_OVERLAP` | adapter needed | C pointer/field-access macros are represented by `RegionRecord::{start,end,payload}` and `RegionIndex::{overlaps,has_overlap}`. |
| `regidx_parse_f`, `regidx_free_f` | out of scope now | C callback signatures and payload-free hooks are replaced by typed Rust parser functions and owned optional payload strings. |
| `regidx_parse_bed`, `regidx_parse_tab`, `regidx_parse_reg` | new Rust implementation | `parse_bed`, `parse_tab`, and `parse_region_line` preserve built-in parser behavior used by `test-regidx.c`, including comments, leading whitespace, coordinate conversion, and zero-coordinate errors. |
| `regidx_parse_vcf` | requires new Rust implementation | Public VCF-line parser is classified but no standalone Rust helper is exposed yet; current ported region-index tests do not exercise it. |
| `regidx_init`, `regidx_init_string` | adapter needed | `RegionIndex::new` plus `insert_line` and `parse_line` cover in-memory and string-fed construction; file-loading/autodetection wrappers are not exposed as C-shaped APIs yet. |
| `regidx_destroy` | out of scope now | Explicit C destruction is represented by Rust ownership and drop semantics. |
| `regidx_overlap`, `regitr_overlap` | new Rust implementation | `RegionIndex::overlaps` and `RegionIndex::has_overlap` preserve overlap detection and iteration behavior, covered by deterministic and randomized `test-regidx.c` parity cases. |
| `regidx_insert`, `regidx_insert_list`, `regidx_push` | adapter needed | `RegionIndex::insert_line` and `RegionIndex::push` cover single-record insertion and sorted storage; delimiter-list insertion is classified but not exposed as a separate Rust helper. |
| `regidx_seq_names`, `regidx_seq_nregs`, `regidx_nregs` | requires new Rust implementation | Sequence-name and count accessors are classified but not exposed as C-shaped helpers; current tests use `RegionIndex::iter` and overlap queries instead. |
| `regitr_init`, `regitr_destroy`, `regitr_reset` | out of scope now | C iterator allocation/reset is replaced by normal Rust iterator construction and ownership. |
| `regitr_loop` | adapter needed | Full-index traversal is represented by `RegionIndex::iter`, covered by sequential-access parity tests. |

### `hts_expr.h`

| Public item(s) | Classification | Rust evidence |
| --- | --- | --- |
| `HTS_EXPR_H` | out of scope now | C include guard; no Rust equivalent needed. |
| `hts_expr_val_t` | new Rust implementation | `htslib_rs::expr::Value` preserves the numeric, string, undefined, and explicitly-true states used by HTSlib expressions. |
| `hts_expr_val_exists` | new Rust implementation | `Value::exists` preserves defined-value checks for numeric NaN and null string cases. |
| `hts_expr_val_existsT` | new Rust implementation | `Value::exists_true` preserves HTSlib's defined-or-explicitly-true behavior. |
| `hts_expr_val_undef` | new Rust implementation | `Value::undefined` represents the public undefined-value state; internal evaluation also resets values to this state where HTSlib returns null/NaN. |
| `hts_expr_val_free` | out of scope now | Explicit C string-memory cleanup is represented by Rust ownership and drop semantics. |
| `hts_filter_t` | new Rust implementation | `htslib_rs::expr::Filter` is the Rust-owned parsed filter expression. |
| `HTS_EXPR_VAL_INIT` | new Rust implementation | `Value::number(0.0)` and `Value::undefined` cover explicit initialization states; Rust construction avoids uninitialized C storage. |
| `hts_filter_init` | new Rust implementation | `Filter::new` constructs a filter from expression text; parse errors are reported during evaluation like the current Rust API design. |
| `hts_filter_free` | out of scope now | Explicit C destruction is represented by Rust ownership and drop semantics. |
| `hts_expr_sym_func` | adapter needed | The C callback type is represented by the `Filter::eval_with` closure `Fn(&str) -> Option<(Value, usize)>`, allowing symbol lookup with consumed input length. |
| `hts_filter_eval2` | new Rust implementation | `Filter::eval_with` evaluates expressions and returns `Value`; acceptance tests port `test_expr.c` arithmetic, boolean, regex, null, string, math-function, and symbol-lookup cases. |
| deprecated `hts_filter_eval` | deprecated but retained for compatibility | The deprecated API is tracked but not exposed separately; Rust callers use `Filter::eval_with`, which owns result cleanup semantics. |

### `synced_bcf_reader.h`

| Public item(s) | Classification | Rust evidence |
| --- | --- | --- |
| `HTSLIB_SYNCED_BCF_READER_H` | out of scope now | C include guard; no Rust equivalent needed. |
| `COLLAPSE_NONE`, `COLLAPSE_SNPS`, `COLLAPSE_INDELS`, `COLLAPSE_ANY`, `COLLAPSE_SOME`, `COLLAPSE_BOTH` | deprecated but retained for compatibility | Legacy collapse constants are tracked; Rust pairing logic uses `SyncedPairLogic` variants instead. |
| `BCF_SR_PAIR_SNPS`, `BCF_SR_PAIR_INDELS`, `BCF_SR_PAIR_ANY`, `BCF_SR_PAIR_SOME`, `BCF_SR_PAIR_SNP_REF`, `BCF_SR_PAIR_INDEL_REF`, `BCF_SR_PAIR_EXACT`, `BCF_SR_PAIR_ID`, `BCF_SR_PAIR_BOTH`, `BCF_SR_PAIR_BOTH_REF` | new Rust implementation | `htslib_rs::variant_io_compat::SyncedPairLogic` covers SNP, indel, both, reference-compatible, exact, subset, and all-compatible pairing modes used by `test-bcf-sr`; ID-only pairing remains classified but not exposed as a standalone mode. |
| `bcf_sr_opt_t`, `bcf_sr_set_opt` | adapter needed | Reader options are represented by purpose-built Rust helper arguments and `SyncedPairLogic`; a variadic C-shaped option API is not exposed. |
| `bcf_sr_regions_t`, `bcf_sr_region_t`, `bcf_sr_regions_init`, `bcf_sr_regions_destroy`, `bcf_sr_regions_seek`, `bcf_sr_regions_next`, `bcf_sr_regions_overlap`, `bcf_sr_regions_flush` | adapter needed | Region and target filtering is represented by Rust parsing/filter helpers for weird chromosome-name VCF fixtures and indexed BCF queries; a standalone mutable synced-region iterator API remains to be exposed if needed. |
| `bcf_sr_t`, `bcf_srs_t` | adapter needed | The C reader and synced-reader state structs are replaced by Rust-owned summary and pairing helpers backed by noodles VCF/BCF readers. |
| `bcf_sr_error`, `bcf_sr_strerror` | requires new Rust implementation | Error categories are classified, but the Rust API currently uses `io::Result`/typed parser errors rather than an HTSlib-style synced-reader error enum and string table. |
| `bcf_sr_init`, `bcf_sr_destroy` | adapter needed | Rust synced-reader state is constructed by `synced_vcf_summary_no_index_from_paths`, `synced_vcf_summary_no_index_from_readers`, and `pair_synced_variant_groups`, with destruction handled by ownership/drop. |
| `bcf_sr_set_threads`, `bcf_sr_destroy_threads` | out of scope now | C thread-pool ownership is outside the Rust-only API; format-level concurrency is handled through noodles worker-count APIs where available. |
| `bcf_sr_add_reader` | adapter needed | `synced_vcf_summary_no_index_from_paths` covers local path reader addition for no-index summary parity. |
| `bcf_sr_add_hreader` | adapter needed | `synced_vcf_summary_no_index_from_readers` covers Rust reader-addition parity equivalent to the C `--usefptr` path without taking ownership of an `htsFile *`. |
| `bcf_sr_remove_reader` | requires new Rust implementation | Dynamic reader removal is classified but not exposed in the current Rust helper surface. |
| `bcf_sr_next_line`, `bcf_sr_has_line`, `bcf_sr_get_line`, `bcf_sr_swap_line`, `bcf_sr_region_done`, `bcf_sr_get_header`, `bcf_sr_get_reader` | adapter needed | No-index synced summaries, VCF-output helpers, and allele-pairing helpers cover the currently ported observable next-line behavior; direct C buffer/header/reader pointer macros are not exposed. |
| `bcf_sr_seek` | adapter needed | Indexed BCF weird chromosome-name region queries exercise seek/query behavior through noodles-backed BCF index helpers; a direct seek API over synced-reader state remains deferred. |
| `bcf_sr_set_samples` | requires new Rust implementation | Sample-subsetting API is classified but not yet exposed as synced-reader behavior. |
| `bcf_sr_set_targets`, `bcf_sr_set_regions` | adapter needed | Region and target query semantics are partially covered by weird chromosome-name fixture tests; broader file-backed targets, complements, allele-target matching, and duplicate-target behavior remain follow-up work. |

### `tbx.h`

| Public item(s) | Classification | Rust evidence |
| --- | --- | --- |
| `HTSLIB_TBX_H` | out of scope now | C include guard; no Rust equivalent needed. |
| `TBX_MAX_SHIFT` | covered by noodles | CSI/TBI binning limits are handled by `noodles-tabix` and `noodles-csi`; custom CSI `min_shift` adapters are covered by index tests. |
| `TBX_GENERIC`, `TBX_SAM`, `TBX_VCF`, `TBX_GAF`, `TBX_UCSC` | adapter needed | `htslib_rs::tabix_compat::TextFormat` covers the tested BED, GFF, and VCF presets; SAM/GAF/UCSC preset aliases are classified but not exposed as standalone Rust variants. |
| `tbx_conf_t`, `tbx_conf_gff`, `tbx_conf_bed`, `tbx_conf_psltbl`, `tbx_conf_sam`, `tbx_conf_vcf`, `tbx_conf_gaf` | adapter needed | Tabix configuration is represented by `TextFormat` and noodles indexer configuration for the tested text formats. |
| `tbx_t` | covered by noodles | Tabix index state is represented by `noodles_tabix::Index`, re-exported through `htslib_rs::tabix`. |
| `tbx_itr_destroy`, `tbx_itr_queryi`, `tbx_itr_querys`, `tbx_itr_querys1`, `tbx_itr_next`, `tbx_bgzf_itr_next` | adapter needed | Region query helpers such as `query_records_from_path`, `query_records_from_path_separate_regions`, and `query_csi_records_from_path` wrap noodles indexed readers instead of exposing C iterator handles. |
| `tbx_name2id` | requires new Rust implementation | Sequence-name-to-ID lookup is classified but not exposed as a C-shaped helper; current tests query by parsed regions. |
| `hts_get_bgzfp`, `tbx_readrec` | out of scope now | C internal BGZF pointer and record-decoding hooks are not part of the Rust-only API. |
| `tbx_index` | adapter needed | In-memory BGZF index building is represented by `build_bgzf_and_index` and `build_bgzf_and_csi`, which return the compressed bytes/writer plus noodles index. |
| `tbx_index_build`, `tbx_index_build2`, `tbx_index_build3` | adapter needed | Local-file index creation is represented by `write_bgzf_and_index`, `write_bgzf_and_csi`, and explicit min-shift CSI helpers; C thread-option behavior is out of scope while CLI/shared thread-pool APIs are deferred. |
| `tbx_index_load`, `tbx_index_load2`, `tbx_index_load3` | adapter needed | Local TBI/CSI loading uses `noodles_tabix::fs::read`, `noodles_csi::fs::read`, and associated-index lookup helpers; remote index streaming/saving flags are out of scope with remote I/O. |
| `tbx_seqnames` | requires new Rust implementation | Sequence-name extraction is classified but not exposed as a C-shaped helper; current acceptance coverage focuses on build/query output. |
| `tbx_destroy` | out of scope now | Explicit C destruction is represented by Rust ownership and drop semantics. |

### `faidx.h`

| Public item(s) | Classification | Rust evidence |
| --- | --- | --- |
| `HTSLIB_FAIDX_H` | out of scope now | C include guard; no Rust equivalent needed. |
| `faidx_t` | covered by noodles | FASTA and FASTQ index state is represented by `noodles_fasta::fai::Index` and `noodles_fastq::fai::Index`, exposed through `faidx_compat::Index` and `FastqIndex`. |
| `struct hts_tpool` forward declaration | out of scope now | Shared C thread-pool objects are outside the Rust-only API target. |
| `enum fai_format_options`, `FAI_NONE`, `FAI_FASTA`, `FAI_FASTQ` | adapter needed | Format selection is represented by separate FASTA/FASTQ helper functions and Rust types rather than a C enum. |
| `fai_build3`, `fai_build` | adapter needed | `build_index` and `build_fastq_index` build FASTA/FASTQ FAI records, with golden-output coverage from `test_faidx.c` and `test/faidx/test-faidx.sh`; bgzip FASTA retrieval also uses GZI support. |
| `fai_destroy` | out of scope now | Explicit C destruction is represented by Rust ownership and drop semantics. |
| `enum fai_load_options`, `FAI_CREATE` | adapter needed | Index creation-on-load is represented by explicit Rust build/read calls; no C-shaped flag API is exposed. |
| `fai_load3`, `fai_load`, `fai_load3_format`, `fai_load_format` | adapter needed | `read_index`, `read_fastq_index`, and callers' explicit file handling cover loaded index behavior for local files. |
| `fai_fetch`, `fai_fetch64`, `faidx_fetch_seq`, `faidx_fetch_seq64` | adapter needed | `fetch_region_sequence` and `fetch_sequence` cover FASTA sequence retrieval by region or explicit coordinates, including HTSlib-style output tests. |
| `fai_line_length` | adapter needed | `line_length` and `fastq_line_length` expose indexed line-base lookup for sequence names; region-string wrapper semantics are classified but not separately exposed. |
| `fai_fetchqual`, `fai_fetchqual64`, `faidx_fetch_qual`, `faidx_fetch_qual64` | adapter needed | `fetch_fastq_region_quality` covers FASTQ quality retrieval for indexed local files and is validated by golden-output tests. |
| deprecated `faidx_fetch_nseq` | deprecated but retained for compatibility | Deprecated sequence-count API remains tracked; use `faidx_nseq` semantics through index length where needed. |
| `faidx_has_seq` | adapter needed | `has_sequence` and `has_fastq_sequence` cover FASTA and FASTQ presence checks. |
| `faidx_nseq` | adapter needed | Sequence count is represented by the length of the noodles index; no C-shaped helper is exposed yet. |
| `faidx_iseq` | adapter needed | `sequence_name` and `fastq_sequence_name` cover sequence-name lookup by index. |
| `faidx_seq_len64`, deprecated `faidx_seq_len` | adapter needed | `sequence_len` and `fastq_sequence_len` expose indexed sequence lengths; the deprecated narrow return form is tracked but not exposed separately. |
| `fai_parse_region` | adapter needed | Region parsing is embedded in retrieval helpers and shared HTSlib-style region parsing; a direct `fai_parse_region` wrapper remains deferred. |
| `fai_adjust_region` | requires new Rust implementation | Boundary-adjustment helper is classified but not exposed as a standalone Rust API. |
| `fai_set_cache_size` | out of scope now | BGZF cache tuning is a C hFILE/BGZF implementation detail not exposed in the Rust-only API. |
| `fai_thread_pool` | out of scope now | Shared C thread-pool assignment is out of scope; format-level concurrency uses noodles worker-count APIs where available. |
| `fai_path` | adapter needed | Associated index path handling is partially covered by local index lookup and `##idx##` parsing in index/reference helpers; a direct `fai_path` API remains deferred. |

### `hfile.h`

| Public item(s) | Classification | Rust evidence |
| --- | --- | --- |
| `HTSLIB_HFILE_H` | out of scope now | C include guard; no Rust equivalent needed. |
| `HTSLIB_SSIZE_T`, `ssize_t` fallback | out of scope now | C platform typedef shim; Rust uses native integer and `usize`/`isize` types. |
| `hFILE_backend`, `hFILE` | out of scope now | The C buffered stream internals and backend vtable are not exposed in the Rust-only API; callers use standard Rust readers/writers and noodles format readers. |
| `kstring_t` forward declaration | out of scope now | C forward declaration only; Rust uses `kstring::KString` where a compatibility string is needed. |
| `hopen`, `hdopen` | adapter needed | Local file opening is represented by `std::fs::File` and higher-level format adapters; remote/plugin-backed opening is explicitly out of scope. |
| `hisremote` | out of scope now | Remote I/O is out of scope for now, so no Rust API exposes HTSlib's remote-storage classifier. |
| `haddextension` | new Rust implementation | `htslib_rs::hfile_compat::add_extension` preserves local path and URL query/fragment extension behavior, covered by `hfile_utils.rs`. |
| `hclose`, `hclose_abruptly` | out of scope now | Explicit C stream close semantics are represented by Rust ownership and drop; flushable writers use Rust `Write` APIs where needed. |
| `herrno`, `hclearerr` | out of scope now | C stream error flags are replaced by Rust `Result` values. |
| `hseek`, `htell` | adapter needed | Seeking/telling is handled through Rust `Seek` on local readers and through format-specific indexed query helpers, not a standalone `hFILE` stream. |
| `hgetc`, `hgetc2`, `hgetdelim`, `hgetln`, `hgets`, `khgetline`, `hpeek`, `hread`, `hread2` | adapter needed | Byte and line reading is represented by Rust `Read`/`BufRead` and `kstring` chunked-line helpers; no mutable C `hFILE` read API is exposed. |
| `hputc`, `hputc2`, `hputs`, `hputs2`, `hwrite`, `hwrite2`, `hflush`, `hfile_set_blksize` | adapter needed | Writing and flushing are represented by Rust `Write` and format-specific writer adapters; C buffer-size tuning is not exposed. |
| `hfile_mem_get_buffer`, `hfile_mem_steal_buffer` | requires new Rust implementation | In-memory hFILE buffer ownership APIs are classified but not implemented as standalone Rust compatibility helpers. |
| `hfile_list_schemes`, `hfile_list_plugins`, `hfile_has_plugin` | out of scope now | Plugin-backed hFILE backends are explicitly out of scope for the Rust-only target. |
| `data:` hFILE plugin behavior | new Rust implementation | `htslib_rs::hfile_compat::decode_data_url` covers plain, empty, percent-encoded, and base64 `data:` URL payloads used by hFILE tests. |

### `hts.h`

| Public item(s) | Classification | Rust evidence |
| --- | --- | --- |
| `HTSLIB_HTS_H`, `HTS_BGZF_TYPEDEF`, `BGZF`, `cram_fd`, `hFILE`, `hts_tpool`, `sam_hdr_t`, `hts_idx_t`, `hts_filter_t` forward declarations | out of scope now | C include guards, typedef guards, and forward declarations are represented by Rust modules and owned types rather than public incomplete C structs. |
| `HTS_PATH_SEPARATOR_CHAR`, `HTS_PATH_SEPARATOR_STR` | adapter needed | Path-list separators are tracked for REF_PATH/plugin compatibility, but remote/plugin behavior is out of scope and no Rust API exposes a global HTS_PATH parser yet. |
| deprecated `hts_expand`, `hts_expand0`, `hts_resize`, `HTS_RESIZE_CLEAR`, `hts_resize_array_` | deprecated but retained for compatibility | C allocation macros are tracked because downstream C-compatibility may need them later; Rust code uses `Vec`/ownership and does not expose realloc-style APIs. |
| `hts_lib_shutdown`, `hts_free` | out of scope now | Dynamic-library cleanup and cross-DLL allocation boundaries are C ABI concerns; current target is Rust-only API. |
| `enum htsFormatCategory`, `enum htsExactFormat`, `enum htsCompression`, `htsFormat` | adapter needed | `htslib_rs::format::{Category, Exact, Compression, Format}` covers local format classification for selected HTSlib fixtures; unsupported exact formats and compression backends remain classified. |
| `htsFile` | adapter needed | C's format-neutral mutable handle is represented by purpose-built local Rust helpers for SAM/BAM/CRAM, VCF/BCF, BGZF, tabix, and FASTA/FASTQ rather than a single public handle. |
| `htsThreadPool`, `hts_set_threads`, `hts_set_thread_pool`, `hts_set_cache_size` | adapter needed | C shared pools are out of scope, but observable BGZF worker-count behavior is covered by noodles multithreaded helpers; generic htsFile and cache APIs remain deferred. |
| `enum sam_fields` and `CRAM_OPT_REQUIRED_FIELDS` | adapter needed | SAM/CRAM record projection flags are classified; current CRAM adapters decode full records through noodles and do not expose required-field selection. |
| `enum hts_fmt_option`, `enum hts_profile_option`, `cram_option`, `hts_opt`, `HTS_FILE_OPTS_INIT`, `hts_opt_add`, `hts_opt_apply`, `hts_opt_free`, `hts_parse_format`, `hts_parse_opt_list`, `hts_set_opt` | adapter needed | Format-specific options are represented by explicit Rust helper arguments such as local FASTA reference repositories, worker counts, compression kind, and FASTQ conversion options; C variadic/list option parsing is not exposed. |
| `HTS_IDX_DELIM` | adapter needed | `index_compat` implements `##idx##` associated-index parsing and candidate ordering for local BAM, CRAM, VCF, BCF, and tabix-backed queries. |
| `seq_nt16_table`, `seq_nt16_str`, `seq_nt16_int` | new Rust implementation | `sequence::SEQ_NT16_STR`, `bam_seqi`, and `nibble_to_bases` cover nt16 decoding used by nibbles tests; full exported table aliases remain to be exposed if callers need them. |
| `hts_version`, `HTS_VERSION`, `hts_features`, `hts_test_feature`, `hts_feature_string`, `HTS_FEATURE_*` | out of scope now | C build/version feature introspection and compile-time feature bit masks are not part of the Rust-only API; shared-library introspection tests are out of scope. |
| `hts_detect_format`, deprecated `hts_detect_format2`, `hts_format_description`, `hts_format_file_extension`, deprecated `hts_file_type`, `FT_*` | adapter needed | `io::detect_format` and `format::detect_path` provide local fixture detection for SAM/BAM/CRAM/VCF/BCF/index/FASTA/FASTQ/BED cases; C hFILE and description/string helper forms remain unexposed. |
| `hts_open`, `hts_open_format`, `hts_hopen`, `hts_flush`, `hts_close`, `hts_get_format` | adapter needed | Opening, flushing, closing, and format access are represented by standard Rust I/O ownership and format-specific helper functions; no monolithic `htsFile` wrapper is exposed yet. |
| `hts_getline`, `hts_readlines`, `hts_readlist` | adapter needed | Line reading is covered through Rust `BufRead` and `kstring::getline_from_chunks`; file/list splitting helpers are classified but not exposed as C-shaped APIs. |
| `hts_set_fai_filename` | adapter needed | CRAM and FASTA reference handling is represented by explicit local FASTA/FAI repository parameters and error-context helpers. |
| `hts_set_filter_expression` | adapter needed | SAM filter expression behavior is covered by Rust expression parsing/evaluation over SAM records; a mutable `htsFile` filter setter is not exposed. |
| `hts_check_EOF` | adapter needed | BGZF EOF marker checks are covered by `bgzf_compat::has_eof_marker`; CRAM EOF behavior is handled through noodles reader outcomes rather than a generic htsFile API. |
| `HTS_IDX_NOCOOR`, `HTS_IDX_START`, `HTS_IDX_REST`, `HTS_IDX_NONE`, `HTS_FMT_CSI`, `HTS_FMT_BAI`, `HTS_FMT_TBI`, `HTS_FMT_CRAI`, `HTS_FMT_FAI`, `HTS_IDX_SAVE_REMOTE`, `HTS_IDX_SILENT_FAIL` | adapter needed | `IndexFormat` and associated-index helpers cover local BAI/CSI/TBI/CRAI selection; special iterator tids and remote-index save flags remain future compatibility work. |
| `HTS_POS_MAX`, `HTS_POS_MIN`, `PRIhts_pos`, `hts_pos_t`, `hts_pair_pos_t`, `hts_pair32_t`, `hts_pair64_t`, `hts_pair64_max_t`, `hts_reglist_t` | adapter needed | Rust uses `i64` coordinates and typed region structs in `region`/`regidx`; C pair/reglist memory layouts are not exposed. |
| `hts_readrec_func`, `hts_seek_func`, `hts_tell_func`, `hts_itr_t`, `hts_itr_multi_t`, `hts_itr_query_func`, `hts_itr_multi_query_func` | adapter needed | Callback-driven C iterators are represented by owning Rust iterators and query helpers for BAM, CRAM, VCF, BCF, and tabix-backed text formats. |
| `hts_bin_first`, `hts_bin_parent`, `hts_reg2bin`, `hts_bin_level`, `hts_bin_bot`, `hts_bin_maxpos` | covered by noodles | Binning math is used through `noodles-csi`, `noodles-bam`, `noodles-tabix`, and related index builders; standalone C-shaped bin helpers are classified but not exposed. |
| `hts_idx_init`, `hts_idx_destroy`, `hts_idx_push`, `hts_idx_finish`, `hts_idx_fmt`, `hts_idx_tbi_name`, `hts_idx_save`, `hts_idx_save_as`, `hts_idx_load`, `hts_idx_load2`, `hts_idx_load3` | adapter needed | Index read/build/write helpers exist for BAI, CSI, TBI, GZI, and CRAI using noodles; C's generic mutable index builder/loader object remains unexposed. |
| `hts_idx_get_meta`, `hts_idx_set_meta`, `hts_idx_get_stat`, `hts_idx_get_n_no_coor`, `hts_idx_seqnames`, `hts_idx_nseq`, `hts_id2name_f` | adapter needed | Current tests validate index build/query/read/write behavior and BAI reference-sequence counts; generic metadata/stat/name C accessors remain follow-up work. |
| `HTS_PARSE_THOUSANDS_SEP`, `HTS_PARSE_ONE_COORD`, `HTS_PARSE_LIST`, `hts_parse_decimal`, `hts_name2id_f`, `hts_parse_reg64`, `hts_parse_reg`, `hts_parse_region` | new Rust implementation | `region` parsing ports `test-parse-reg.c` cases, including thousands separators, one-coordinate mode, comma-list mode, quoted names, ambiguity, and invalid-region behavior. |
| `hts_itr_query`, `hts_itr_destroy`, `hts_itr_querys`, `hts_itr_next`, `hts_itr_multi_bam`, `hts_itr_multi_cram`, `hts_itr_regions`, `hts_itr_multi_next`, `hts_reglist_create`, `hts_reglist_free`, `hts_itr_multi_destroy` | adapter needed | Format-specific Rust query/view helpers cover single- and multi-region behavior, including HTSlib-style de-duplication for BAM/CRAM; generic callback iterator APIs remain unexposed. |
| deprecated `FT_*`, `hts_file_type` | deprecated but retained for compatibility | Legacy file-type bits remain tracked; local detection uses the Rust `Format` model. |
| `errmod_t`, `errmod_init`, `errmod_destroy`, `errmod_cal` | requires new Rust implementation | Revised MAQ error model APIs are classified but not implemented; no current ported tests cover them directly. |
| `probaln_par_t`, `probaln_glocal` | new Rust implementation | `probaln::ProbalnParams` and `probaln::probaln_glocal` cover initial HTSlib-oracle likelihood, state, and posterior-quality parity for exact, mismatch, and insertion cases and are exercised through BAQ recalculation parity for `test_realn.c` fixtures. |
| `hts_md5_context`, `hts_md5_init`, `hts_md5_update`, `hts_md5_final`, `hts_md5_reset`, `hts_md5_hex`, `hts_md5_destroy` | requires new Rust implementation | MD5 helpers are classified but not exposed as standalone Rust compatibility APIs. |
| `hts_crc32` | requires new Rust implementation | CRC32 public helper is classified; BGZF validation uses noodles/codec behavior and no C-shaped Rust helper is exposed. |
| `ed_is_big`, `ed_swap_2`, `ed_swap_2p`, `ed_swap_4`, `ed_swap_4p`, `ed_swap_8`, `ed_swap_8p` | new Rust implementation | Endian behavior is covered through `hts_endian` Rust helpers and acceptance tests; these older inline names remain classified but not separately exposed. |

### `sam.h`

| Public item(s) | Classification | Rust evidence |
| --- | --- | --- |
| `HTSLIB_SAM_H`, `HTSLIB_SSIZE_T`, `ssize_t` fallback | out of scope now | C include guards and platform typedef shims have no Rust API equivalent. |
| `SAM_FORMAT_VERSION` | covered by noodles | SAM version parsing/writing is handled by `noodles-sam`; the Rust adapters target upstream SAM/BAM/CRAM fixture behavior rather than exposing this macro. |
| `sam_hrecs_t`, `sam_hdr_t`, deprecated `bam_hdr_t`, `samFile` | adapter needed | SAM headers and file handles are represented by noodles SAM/BAM/CRAM readers and the Rust `SamHeaderAdapter`; C layout compatibility remains future work. |
| `BAM_CMATCH`, `BAM_CINS`, `BAM_CDEL`, `BAM_CREF_SKIP`, `BAM_CSOFT_CLIP`, `BAM_CHARD_CLIP`, `BAM_CPAD`, `BAM_CEQUAL`, `BAM_CDIFF`, `BAM_CBACK`, `BAM_CIGAR_STR`, `BAM_CIGAR_SHIFT`, `BAM_CIGAR_MASK`, `BAM_CIGAR_TYPE`, `bam_cigar_table`, `bam_cigar_op`, `bam_cigar_oplen`, `bam_cigar_opchr`, `bam_cigar_gen`, `bam_cigar_type` | adapter needed | CIGAR semantics are covered through noodles record CIGAR APIs and Rust SAM filter/pileup helpers; C-shaped constants/macros are classified but not all exposed as public aliases. |
| `BAM_FPAIRED`, `BAM_FPROPER_PAIR`, `BAM_FUNMAP`, `BAM_FMUNMAP`, `BAM_FREVERSE`, `BAM_FMREVERSE`, `BAM_FREAD1`, `BAM_FREAD2`, `BAM_FSECONDARY`, `BAM_FQCFAIL`, `BAM_FDUP`, `BAM_FSUPPLEMENTARY`, `bam_is_rev`, `bam_is_mrev` | adapter needed | Flag behavior is represented by noodles SAM flags and `SamRecordAdapter`; TLEN, view, filter, and FASTQ conversion tests exercise key flag semantics. |
| `bam1_core_t`, `bam1_t`, `bam_get_qname`, `bam_get_cigar`, `bam_get_seq`, `bam_get_qual`, `bam_get_aux`, `bam_get_l_aux`, `bam_seqi`, `bam_set_seqi` | adapter needed | Raw BAM record layout is replaced by noodles record types and safe adapters; `sequence::bam_seqi` and nibbles tests cover packed sequence lookup, while direct byte-layout access remains unexposed. |
| `sam_hdr_init`, `sam_hdr_destroy`, `sam_hdr_dup`, deprecated `bam_hdr_init`, `bam_hdr_destroy`, `bam_hdr_dup`, `sam_hdr_incr_ref` | adapter needed | Header lifecycle is represented by Rust ownership/clone/drop and `SamHeaderAdapter`; explicit C reference-count APIs are not exposed. |
| `bam_hdr_read`, `bam_hdr_write`, `sam_hdr_parse`, `sam_hdr_read`, `sam_hdr_write`, `sam_hdr_length`, `sam_hdr_str`, `sam_hdr_nref` | adapter needed | Header reading/writing and raw-header preservation are covered for SAM, BAM, and CRAM through noodles-backed helpers and view/write parity tests. |
| line-level header APIs `sam_hdr_add_lines`, `sam_hdr_add_line`, `sam_hdr_find_line_id`, `sam_hdr_find_line_pos`, `sam_hdr_remove_line_id`, `sam_hdr_remove_line_pos`, `sam_hdr_update_line`, `sam_hdr_remove_except`, `sam_hdr_remove_lines`, `sam_hdr_count_lines`, `sam_hdr_line_index`, `sam_hdr_line_name` | adapter needed | `SamHeaderAdapter` covers reference-sequence lookup, insert, replace, and removal for tested behavior; broader line/tag editing APIs are classified for future expansion. |
| tag-level header APIs `sam_hdr_find_tag_id`, `sam_hdr_find_tag_pos`, `sam_hdr_remove_tag_id`, `sam_hdr_find_hd`, `sam_hdr_find_tag_hd`, `sam_hdr_update_hd`, `sam_hdr_remove_tag_hd`, `sam_hdr_change_HD` | adapter needed | Header tag mutation is partially represented by `SamHeaderAdapter` and raw-header view preservation; C variadic tag APIs are not exposed directly. |
| reference/header helpers `sam_hdr_name2tid`, `sam_hdr_tid2name`, `sam_hdr_tid2len`, deprecated `bam_name2id`, `sam_hdr_pg_id`, `sam_hdr_add_pg`, `stringify_argv`, `sam_hdr_set`, `sam_hdr_get` | adapter needed | Reference-sequence lookup and mutation are covered by the SAM header adapter; PG-line generation, command-line stringification, and htsFile header setters/getters remain future compatibility work. |
| `bam_init1`, `bam_destroy1`, `BAM_USER_OWNS_STRUCT`, `BAM_USER_OWNS_DATA`, `bam_set_mempolicy`, `bam_get_mempolicy`, `bam_copy1`, `bam_dup1`, `bam_set1`, `bam_set_qname` | adapter needed | Alignment record lifecycle and mutation are represented by noodles-owned records and `SamRecordAdapter`; raw memory-policy controls are C-only and not exposed. |
| `bam_read1`, `bam_write1`, `sam_parse1`, `sam_format1`, `sam_read1`, `sam_write1` | adapter needed | SAM/BAM/CRAM read, write, parse, format, and view behavior is covered by noodles-backed helpers and upstream fixture parity tests, including SAM-to-BAM, SAM-to-CRAM, BAM/CRAM view, and record-limit cases. |
| `bam_cigar2qlen`, `bam_cigar2rlen`, `bam_endpos`, `bam_str2flag`, `bam_flag2str`, `sam_parse_cigar`, `bam_parse_cigar` | adapter needed | CIGAR-derived metrics and flag behavior are exercised through SAM filter and pileup parity tests; standalone C-shaped helpers remain to be exposed if needed. |
| BAM/SAM/CRAM index macros and APIs `bam_itr_destroy`, `bam_itr_queryi`, `bam_itr_querys`, `bam_itr_next`, `bam_index_load`, `bam_index_build`, `sam_idx_init`, `sam_idx_save`, `sam_index_load`, `sam_index_load2`, `sam_index_load3`, `sam_index_build`, `sam_index_build2`, `sam_index_build3`, `sam_itr_destroy`, `sam_itr_queryi`, `sam_itr_querys`, `sam_itr_regions`, `sam_itr_regarray`, `sam_itr_next`, `sam_itr_multi_next`, `sam_parse_region` | adapter needed | BAI/CSI/CRAI build/load/query helpers and owning region iterators exist for BAM, BGZF SAM, and CRAM, including associated-index lookup and multi-region de-duplication; C iterator handles and writer-side on-handle indexing remain unexposed. |
| `sam_open`, `sam_open_format`, `sam_flush`, `sam_close`, `sam_open_mode`, `sam_open_mode_opts` | adapter needed | Opening/closing and mode handling are represented by Rust path/reader/writer helpers and format-specific write options; C mode-string construction is only partially covered by VCF open-mode tests, not SAM-specific aliases. |
| `sam_passes_filter` | adapter needed | SAM filter expressions are implemented as Rust expression helpers over SAM records and covered by upstream `sam_filter` fixtures; the C `hts_filter_t *` entry point is not exposed. |
| `sam_format_aux1`, `bam_aux_first`, `bam_aux_next`, `bam_aux_get`, `bam_aux_tag`, `bam_aux_type`, `bam_aux_get_str`, `bam_aux2i`, `bam_aux2f`, `bam_aux2A`, `bam_aux2Z`, `bam_auxB_len`, `bam_auxB2i`, `bam_auxB2f`, `bam_aux_append`, `bam_aux_del`, `bam_aux_remove`, `bam_aux_update_str`, `bam_aux_update_int`, `bam_aux_update_float`, `bam_aux_update_array` | adapter needed | `SamAuxAdapter` covers typed auxiliary get/insert/remove operations for Rust records; full raw BAM aux pointer iteration, formatting, deletion-order semantics, and array updates remain future compatibility work. |
| `bam_pileup_cd`, `bam_pileup1_t`, `bam_plp_auto_f`, `bam_plp_t`, `bam_mplp_t`, `bam_plp_init`, `bam_plp_destroy`, `bam_plp_push`, `bam_plp_next`, `bam_plp_auto`, `bam_plp64_next`, `bam_plp64_auto`, `bam_plp_set_maxcnt`, `bam_plp_reset`, `bam_plp_constructor`, `bam_plp_destructor`, `bam_plp_insertion`, `bam_plp_insertion_mod`, `bam_mplp_init`, `bam_mplp_init_overlaps`, `bam_mplp_destroy`, `bam_mplp_set_maxcnt`, `bam_mplp_auto`, `bam_mplp64_auto`, `bam_mplp_reset`, `bam_mplp_constructor`, `bam_mplp_destructor` | adapter needed | Pileup output parity is implemented for selected SAM/BAM/mpileup and base-modification fixtures, including overlap-removal and insertion/deletion rendering; generic C pileup iterator state and callbacks remain unexposed. |
| `sam_cap_mapq`, `enum htsRealnFlags`, `sam_prob_realn` | requires new Rust implementation | `sam_cap_mapq` MAPQ capping is covered for matched, mismatched, and clipped SAM records. Existing BAQ `BQ:Z`/`ZQ:Z` apply/revert/no-op behavior plus default, apply-on-recalculation, forced, and extended BAQ recalculation are covered against upstream `test_realn.c` expected output; C-shaped flags remain represented by Rust helper functions rather than a raw C API. |
| `hts_base_mod`, `HTS_MOD_UNKNOWN`, `HTS_MOD_UNCHECKED`, `HTS_MOD_REPORT_UNCHECKED`, `hts_base_mod_state`, `hts_base_mod_state_alloc`, `hts_base_mod_state_free`, `bam_parse_basemod`, `bam_parse_basemod2`, `bam_mods_at_next_pos`, `bam_next_basemod`, `bam_mods_at_qpos`, `bam_mods_query_type`, `bam_mods_queryi`, `bam_mods_recorded` | adapter needed | Base-modification parsing and reporting is implemented for upstream MM/ML fixtures, ChEBI IDs, explicit/implicit status, MN validation, bounds rejection, and pileup output; the C stateful iterator API is represented by Rust reporting helpers rather than exposed directly. |

### `vcf.h`

| Public item(s) | Classification | Rust evidence |
| --- | --- | --- |
| `HTSLIB_VCF_H` | out of scope now | C include guard; no Rust equivalent needed. |
| header line/type/number/dictionary constants `BCF_HL_*`, `BCF_HT_*`, `BCF_VL_*`, `BCF_DT_*` | adapter needed | `VcfHeaderId` and header number-classification helpers cover FILTER, INFO, FORMAT, contig, generic records, and VCF 4.4/4.5 `P`, `LA`, `LG`, `LR`, and `M` cases. |
| `bcf_hrec_t`, `bcf_idinfo_t`, `bcf_idpair_t`, `bcf_hdr_t`, `bcf_type_shift` | adapter needed | Header state is represented by noodles VCF headers plus Rust adapters for HTSlib-style IDs, record lookup/removal, typed metadata, and raw-header preservation. |
| typed-value constants `BCF_BT_*`, `bcf_int*_missing`, `bcf_int*_vector_end`, `bcf_str_missing`, `BCF_MAX_BT_*`, `BCF_MIN_BT_*`, `bcf_float_missing`, `bcf_float_vector_end`, float missing/vector-end helpers | adapter needed | Typed INFO/FORMAT integer and float retrieval covers HTSlib missing and vector-end behavior, including BCF FORMAT float normalization and large genotype allele indexes. |
| variant constants and structs `VCF_REF`, `VCF_SNP`, `VCF_MNP`, `VCF_INDEL`, `VCF_OTHER`, `VCF_BND`, `VCF_OVERLAP`, `VCF_INS`, `VCF_DEL`, `VCF_ANY`, `bcf_variant_t` | new Rust implementation | `htslib_rs::variant::classify_variant` and `test-bcf_set_variant_type.c` parity cover variant classification behavior. |
| `bcf_fmt_t`, `bcf_info_t`, `bcf_dec_t`, `bcf1_t`, `BCF1_DIRTY_*`, `BCF_UN_*` | adapter needed | Record state is represented by noodles VCF/BCF records plus `VcfRecordAdapter`; dirty/unpack internals are not exposed as C-shaped mutable state. |
| `BCF_ERR_*`, `bcf_strerror` | requires new Rust implementation | BCF parser error-bit reporting is classified but current Rust APIs return typed `Result` values and test-observable parse errors rather than an HTSlib error bitmask formatter. |
| compatibility macros `bcf_init1`, `bcf_read1`, `vcf_read1`, `bcf_write1`, `vcf_write1`, `bcf_destroy1`, `bcf_empty1`, `vcf_parse1`, `bcf_clear1`, `vcf_format1` | adapter needed | The equivalent operations are represented by Rust constructors, readers, writers, parsers, and serializers; macro aliases are tracked but not exposed. |
| `bcf_hdr_init`, `bcf_hdr_destroy`, `bcf_init`, `bcf_destroy`, `bcf_empty`, `bcf_clear`, `vcfFile`, `bcf_open`, `vcf_open`, `bcf_flush`, `bcf_close`, `vcf_close` | adapter needed | Rust ownership and format-specific helpers replace C allocation and `htsFile` handle lifecycles. |
| `bcf_hdr_read`, `bcf_hdr_write`, `vcf_hdr_read`, `vcf_hdr_write`, `vcf_parse`, `vcf_format`, `bcf_read`, `bcf_write`, `vcf_read`, `vcf_write`, `bcf_readrec`, `vcf_write_line`, `vcf_open_mode` | adapter needed | VCF/BCF read, write, view, record-limit, exact VCF output, BCF-to-VCF, and open-mode extension cases are covered by noodles-backed helpers and upstream fixtures. |
| `bcf_hdr_set_samples`, `bcf_subset_format`, `bcf_hdr_subset`, `bcf_subset` | requires new Rust implementation | Sample subsetting is classified but not implemented as a standalone Rust API; current selected tests focus on typed value access, querying, and record mutation. |
| `bcf_hdr_dup`, deprecated `bcf_hdr_combine`, `bcf_hdr_merge`, `bcf_hdr_add_sample`, `bcf_hdr_set`, `bcf_hdr_format`, deprecated `bcf_hdr_fmt_text`, `bcf_hdr_append`, `bcf_hdr_printf`, `bcf_hdr_get_version`, `bcf_hdr_set_version`, `bcf_hdr_remove`, `bcf_hdr_seqnames`, `bcf_hdr_nsamples`, `bcf_hdr_parse`, `bcf_hdr_sync` | adapter needed | Header duplicate/merge/mutation and serialization behavior is partially covered by BCF translation, header get/remove cases, raw-header preservation, and exact output tests; deprecated APIs remain tracked. |
| header-record APIs `bcf_hdr_parse_line`, `bcf_hrec_format`, `bcf_hdr_add_hrec`, `bcf_hdr_get_hrec`, `bcf_hrec_dup`, `bcf_hrec_add_key`, `bcf_hrec_set_val`, `bcf_hrec_find_key`, `hrec_add_idx`, `bcf_hrec_destroy` | adapter needed | Structured header lookup/removal is covered for FILTER, INFO, FORMAT, contig, and generic records; direct mutable `bcf_hrec_t` construction APIs remain future adapter work. |
| `bcf_translate` | adapter needed | Synthetic BCF translation parity from `test-bcf-translate.c` is covered using noodles VCF record serialization and Rust header merging helpers. |
| deprecated `bcf_get_variant_types`, deprecated `bcf_get_variant_type`, `enum bcf_variant_match`, `bcf_has_variant_types`, `bcf_has_variant_type`, `bcf_variant_length`, deprecated `bcf_is_snp` | new Rust implementation | Variant type, length, and match-mode behavior is covered by `test-bcf_set_variant_type.c`; deprecated names remain tracked. |
| FILTER/allele/ID mutation APIs `bcf_update_filter`, `bcf_add_filter`, `bcf_remove_filter`, `bcf_has_filter`, `bcf_update_alleles`, `bcf_update_alleles_str`, `bcf_update_id`, `bcf_add_id` | adapter needed | `VcfRecordAdapter` and VCF API tests cover header-aware record mutation, allele changes, record serialization, invalid `END`, and `rlen` recalculation. |
| INFO update APIs `bcf_update_info_int32`, `bcf_update_info_float`, `bcf_update_info_flag`, `bcf_update_info_string`, `bcf_update_info`, `bcf_update_info_int64` | adapter needed | INFO integer/float/missing retrieval and record serialization after INFO mutation are covered; full C-shaped typed update surface remains partially exposed through adapters. |
| FORMAT/genotype update APIs `bcf_update_format_int32`, `bcf_update_format_float`, `bcf_update_format_char`, `bcf_update_genotypes`, `bcf_update_format_string`, `bcf_update_format` | adapter needed | FORMAT integer, BCF float missing/vector-end, genotype access, and mutation/serialization parity are covered by `test-vcf-api.c` cases. |
| genotype helpers `bcf_gt_phased`, `bcf_gt_unphased`, `bcf_gt_missing`, `bcf_gt_is_missing`, `bcf_gt_is_phased`, `bcf_gt_allele`, `bcf_alleles2gt`, `bcf_gt2alleles` | adapter needed | Genotype encoding/decoding is handled through noodles and Rust record adapters; diploid `Number=G` triangular mapping is exposed as `variant::bcf_ij2g`. |
| lookup/get APIs `bcf_get_fmt`, `bcf_get_info`, `bcf_get_fmt_id`, `bcf_get_info_id`, `bcf_get_info_*`, `bcf_get_info_values`, `bcf_get_format_*`, `bcf_get_format_string`, `bcf_get_format_values`, `bcf_get_genotypes` | adapter needed | Typed INFO/FORMAT retrieval is covered for integer, float, missing values, scalar/vector values, BCF float vector-end, and genotypes. |
| header ID helpers `bcf_hdr_id2int`, `bcf_hdr_int2id`, `bcf_hdr_name2id`, `bcf_hdr_id2name`, `bcf_seqname`, `bcf_seqname_safe`, `bcf_hdr_id2length`, `bcf_hdr_id2number`, `bcf_hdr_id2type`, `bcf_hdr_id2coltype`, `bcf_hdr_idinfo_exists`, `bcf_hdr_id2hrec` | adapter needed | `VcfHeaderId` and header adapters cover HTSlib-style ID lookup for FILTER, INFO, FORMAT, contig, and generic records, plus Number/type classification. |
| low-level BCF formatting/encoding helpers `bcf_fmt_array`, `bcf_fmt_sized_array`, `bcf_enc_vchar`, `bcf_enc_vint`, `bcf_enc_vfloat`, `bcf_format_gt_v2`, `bcf_format_gt`, `bcf_enc_size`, `bcf_enc_inttype`, `bcf_enc_int1`, `bcf_dec_int1`, `bcf_dec_typed_int1`, `bcf_dec_size` | covered by noodles | BCF encoding/decoding is delegated to `noodles-bcf`, including local patches for large GT encodings and vector-end padding; standalone helper aliases are classified but not exposed. |
| BCF index APIs `bcf_itr_destroy`, `bcf_itr_queryi`, `bcf_itr_querys`, `bcf_itr_querys1`, `bcf_itr_next`, `bcf_index_load`, `bcf_index_seqnames`, `bcf_index_load2`, `bcf_index_load3`, `bcf_index_build`, `bcf_index_build2`, `bcf_index_build3`, `bcf_idx_init`, `bcf_idx_save` | adapter needed | BCF CSI build/read/query, explicit min_shift, indexed iterator creation by numeric ID and region string, and associated-index lookup are covered; on-handle index initialization/save remains future work. |

### `bgzf.h`

| Public item(s) | Classification | Rust evidence |
| --- | --- | --- |
| `HTSLIB_BGZF_H` | out of scope now | C include guard; no Rust equivalent needed. |
| `HTSLIB_SSIZE_T`, `ssize_t` fallback | out of scope now | C platform typedef shim; Rust uses native integer and `usize`/`isize` types. |
| `BGZF_BLOCK_SIZE`, `BGZF_MAX_BLOCK_SIZE` | covered by noodles | BGZF block sizing is handled by `noodles-bgzf`; the Rust compatibility layer validates block headers when building GZI indexes. |
| `BGZF_ERR_ZLIB`, `BGZF_ERR_HEADER`, `BGZF_ERR_IO`, `BGZF_ERR_MISUSE`, `BGZF_ERR_MT`, `BGZF_ERR_CRC` | requires new Rust implementation | Error bit constants are classified but not exposed; Rust helpers currently return `io::Result` and typed virtual-offset errors. |
| `hFILE`, `hts_tpool`, `kstring_t`, `bgzf_mtaux_t`, `bgzidx_t`, `bgzf_cache_t`, `z_stream_s` forward/internal types | out of scope now | C stream, thread-pool, cache, zlib, and internal index types are represented by Rust ownership, noodles types, and standard I/O traits where needed. |
| `BGZF` | covered by noodles | BGZF reader/writer state is represented by `noodles_bgzf::io::{Reader,Writer,IndexedReader,MultithreadedReader,MultithreadedWriter}` and Rust helpers in `bgzf_compat`. |
| `HTS_BGZF_TYPEDEF` | out of scope now | C typedef guard; no Rust equivalent needed. |
| `bgzf_dopen`, `bgzf_fdopen`, `bgzf_open`, `bgzf_hopen`, `bgzf_close` | adapter needed | C handle/open/close wrappers are not exposed directly; Rust uses standard readers/writers and `read_all`, `write_all`, `write_all_with_kind`, and file-specific callers. C open-wrapper tests are out of scope for the Rust-only target. |
| `bgzf_read`, `bgzf_read_small` | covered by noodles | `read_all` and direct noodles readers cover BGZF decompression, including upstream fixture and truncated-stream acceptance tests. |
| `bgzf_write`, `bgzf_write_small`, `bgzf_flush`, `bgzf_flush_try` | covered by noodles | `write_all` and noodles writers cover BGZF writing, finishing, EOF marker emission, and round-trip tests. |
| `bgzf_block_write` | requires new Rust implementation | Index-guided BGZF block layout writing is classified but not exposed as a separate Rust helper. |
| `bgzf_peek`, `bgzf_getc`, `bgzf_getline` | adapter needed | Line-oriented behavior is covered by BGZF getline-style tests using noodles readers and Rust `BufRead`; one-byte peek/get APIs are not exposed directly. |
| `bgzf_raw_read`, `bgzf_raw_write`, `bgzf_read_block` | out of scope now | Raw block/underlying-stream access is a C low-level escape hatch and not exposed in the current Rust-only API. |
| `bgzf_tell` | adapter needed | `virtual_offset`, `virtual_offset_parts`, `virtual_offset_from_position`, and noodles `VirtualPosition` preserve HTSlib virtual-offset layout. |
| `bgzf_seek` | adapter needed | Virtual-position seek behavior is covered by tests using noodles BGZF readers and recorded virtual positions. |
| `bgzf_check_EOF` | adapter needed | `has_eof_marker` checks the canonical empty BGZF EOF block and is covered by EOF and embedded-EOF tests. |
| `bgzf_compression` | adapter needed | `detect_compression_kind` and `read_auto` classify and read plain, gzip, and BGZF streams for local byte inputs. |
| deprecated `bgzf_is_bgzf` | deprecated but retained for compatibility | Deprecated BGZF detection is tracked; callers should use compression detection or format detection helpers. |
| `bgzf_set_cache_size` | out of scope now | C BGZF cache tuning is not exposed in the Rust-only API. |
| `bgzf_thread_pool`, `bgzf_mt` | adapter needed | Shared C pools are out of scope, but equivalent observable concurrency is covered by `read_all_with_worker_count` and `write_all_with_worker_count` using noodles multithreaded BGZF APIs. |
| `bgzf_compress` | requires new Rust implementation | Single-block compression helper is classified but not exposed as a Rust compatibility function. |
| `bgzf_useek`, `bgzf_utell` | adapter needed | GZI-backed uncompressed-offset seeking is covered by `read_gzi`, `query_gzi`, `build_gzi`, and noodles `IndexedReader` tests. |
| `bgzf_index_build_init`, `bgzf_index_load`, `bgzf_index_load_hfile`, `bgzf_index_dump`, `bgzf_index_dump_hfile` | adapter needed | GZI index read/write/build/query behavior is covered by `read_gzi`, `write_gzi`, `build_gzi`, and path helpers; C on-handle and hFILE index attachment/dump APIs are not exposed directly. |

### `cram.h`

| Public item(s) | Classification | Rust evidence |
| --- | --- | --- |
| `HTSLIB_CRAM_H` | out of scope now | C include guard; no Rust equivalent needed. |
| `enum cram_block_method`, `CRAM_COMP_*`, legacy `GZIP`, `BZIP2`, `LZMA`, `RANS4x8`, `RANSNx16`, `ARITH`, `FQZ`, `TOK3` aliases | adapter needed | Compression methods are parsed and written through `noodles-cram`; C enum values and legacy aliases remain tracked for compatibility but are not exposed as a Rust enum yet. |
| `struct cram_method_details` | requires new Rust implementation | HTSlib codec-detail metadata is classified but not exposed; noodles hides most codec selection behind CRAM reader/writer configuration. |
| `enum cram_content_type`, `CT_ERROR`, `FILE_HEADER`, `COMPRESSION_HEADER`, `MAPPED_SLICE`, `UNMAPPED_SLICE`, `EXTERNAL`, `CORE` | covered by noodles | CRAM container, slice, core, and external block structure is decoded and encoded by `noodles-cram`; no standalone HTSlib content-type API is exposed. |
| `cram_fd`, `cram_container`, `cram_block`, `cram_slice`, `cram_codec`, `cram_metrics`, `cram_block_slice_hdr`, `cram_block_compression_hdr`, `cram_file_def`, `refs_t` | adapter needed | High-level CRAM file, header, record, reference, and CRAI behavior uses noodles-owned readers, writers, repositories, and index types. Low-level mutable C structs remain future compatibility work. |
| `cram_fd_get_header`, `cram_fd_set_header`, `cram_fd_get_version`, `cram_fd_set_version`, `cram_major_vers`, `cram_minor_vers` | adapter needed | `read_cram_header` and `read_cram_header_from_path` expose decoded SAM header behavior through noodles; direct CRAM version mutation/accessors remain unexposed. |
| `cram_fd_get_fp`, `cram_fd_set_fp` | out of scope now | Raw `hFILE *` access is a C stream escape hatch; the Rust-only target uses owned local readers and writers. |
| `cram_container_*` accessors for length, blocks, landmarks, records, bases, empty state, and coordinates | requires new Rust implementation | Container metadata is classified for future compatibility, but current Rust APIs expose record/header/query behavior rather than mutable container internals. |
| `cram_block_*` allocation, read/write, compress/uncompress, getters, setters, append/update, offset, and size helpers | requires new Rust implementation | Low-level block manipulation is an HTSlib internal-style API; current coverage relies on noodles CRAM read/write behavior and golden decode summaries. |
| `cram_codec_id2name`, `cram_codec_id2str`, `cram_codec_name2id`, `cram_codec_describe` | requires new Rust implementation | Codec name/description lookup is classified but not exposed separately from noodles CRAM encoding/decoding. |
| compression-header and data-series helpers `cram_get_block_by_id`, `cram_block_compression_hdr_*`, `cram_describe_encodings`, `cram_cid2ds_*`, `cram_stats_*`, `cram_get_encoding_by_id` | requires new Rust implementation | HTSlib's compression-header and encoding-stat inspection APIs are not currently surfaced by the Rust adapters. |
| slice-header helpers `cram_block_slice_hdr_*`, `cram_slice_hdr_get_num_blocks`, `cram_slice_hdr_get_embed_ref_id`, `cram_slice_hdr_set_embed_ref_id`, `cram_slice_hdr_get_ref_base_id`, `cram_slice_hdr_set_ref_base_id` | requires new Rust implementation | Slice internals are classified but current Rust tests validate observable record output, reference-backed MD/NM regeneration, and CRAI queries instead. |
| `cram_open`, `cram_dopen`, `cram_close`, `cram_seek`, `cram_flush`, `cram_eof`, `cram_check_EOF` | adapter needed | Local CRAM open/read/write/seek-style behavior is represented by path-based helpers, owning iterators, writer helpers, CRAI indexed queries, and Rust drop/finish semantics. File-descriptor wrappers are out of scope. |
| `cram_set_option`, `cram_set_voption`, `cram_set_header` | adapter needed | CRAM reference selection is represented by `cram_reference_repository_from_fasta_path` and writer/query helpers that accept explicit local FASTA repositories; variadic HTSlib options remain unexposed. |
| `cram_get_refs`, `cram_set_refs` | adapter needed | External-reference decoding uses `noodles_fasta::Repository` over indexed local FASTA files. Shared `refs_t` ownership and cache tuning are not exposed; noodles reference behavior is used where available. |
| `cram_index_extents`, `cram_num_containers`, `cram_num_containers_between`, `cram_container_num2offset`, `cram_container_offset2num` | adapter needed | `read_cram_crai`, `build_cram_crai`, `write_cram_crai`, and CRAM indexed query helpers cover local CRAI read/write/query behavior; container-number and byte-offset inspection helpers remain future API work. |
| `cram_transcode_rg`, `cram_filter`, `cram_copy_slice` | requires new Rust implementation | HTSlib-specific slice transformation/filtering APIs are classified but not covered by current Rust-only CRAM helpers or selected tests. |
| `int32_put_blk` | requires new Rust implementation | Block-level integer encoding is classified but not exposed because raw block mutation is not part of the current Rust API. |
| `typedef sam_hdr_t SAM_hdr`, `sam_hdr_parse_`, `sam_hdr_free`, `sam_hdr_add_PG` | deprecated but retained for compatibility | Legacy CRAM-header SAM aliases are tracked; Rust header behavior is implemented through `noodles-sam` and SAM header adapters rather than deprecated C names. |

## Implemented So Far

- Root Cargo workspace.
- `htslib-rs` library crate.
- `htslib-test-harness` acceptance-test crate.
- Local path dependency on the checked-out noodles crates.
- Initial Rust-only API documentation.
- Initial format/category/compression model.
- Local BGZF read/write helpers backed by `noodles-bgzf`.
- Local BGZF multithreaded read/write helpers backed by `noodles-bgzf` worker-count readers and writers.
- HTSlib-style BGZF write-mode adapters for plain (`wu`), gzip (`wg`), and BGZF streams.
- Local GZI read/write/build/query helpers backed by `noodles-bgzf`.
- BGZF acceptance coverage for upstream read fixtures, write/read round trips, worker-count read/write, EOF marker detection, embedded EOF blocks, GZI seek, virtual-position seek, getline-style reads, truncated-stream errors, and the `test_rebgzip` multi-block BGZF/GZI boundary fixture.
- C-specific BGZF open-wrapper tests are out of scope for the Rust-only API target.
- Local FASTA/FAI index read, build, write, lookup, and uncompressed region retrieval helpers using `noodles-fasta` FAI records.
- Synthetic BCF/VCF header merge and record translation fixture coverage using noodles VCF record serialization.
- Local FASTQ FAI index read, build, write, lookup, and uncompressed sequence/quality region retrieval helpers using `noodles-fastq` FAI records.
- FASTQ/FASTA to SAM and SAM to FASTQ/FASTA conversion helpers for HTSlib fastq scripted fixtures, including aux, paired, name2, UMI, CASAVA, barcode, and filtered-read cases.
- Local BAI, CSI, and TBI index loading helpers backed by `noodles-bam`, `noodles-csi`, and `noodles-tabix`.
- Local CRAI index loading, building, and writing helpers backed by `noodles-cram`.
- Local BAM BAI/CSI, SAM.gz BAI/CSI, CRAM CRAI, VCF TBI/CSI, and BCF CSI index build/write helpers backed by `noodles-bam`, `noodles-sam`, `noodles-cram`, `noodles-tabix`, `noodles-csi`, and `noodles-bcf`, including explicit CSI min_shift adapters.
- HTSlib-style local associated-index candidate ordering, replaced-extension fallback, and `##idx##` delimiter parsing.
- Local BED, GFF, and VCF BGZF+TBI index building and region query helpers backed by `noodles-tabix` and `noodles-csi`.
- Tabix-style separate-region output formatting for multi-region text queries.
- Local VCF BGZF+CSI index building and region query helpers for large-coordinate and long-reference fixtures, including adaptive CSI depth and `INFO/END` overlap behavior.
- `test/tabix/test-tabix.sh` command-line and thread-option checks are out of scope while CLI tools are deferred.
- Local SAM and BAM header reading and record-count helpers, plus CRAM header, embedded-reference record iteration, and CRAI indexed query helpers with local FASTA reference repositories, backed by `noodles-sam`, `noodles-bam`, and `noodles-cram`.
- CRAM external-reference decoding uses noodles' FASTA reference repository API over indexed local FASTA files.
- HTSlib `test_view.c`-style SAM text view record limiting with header preservation.
- HTSlib `test_view.c`-style SAM parse-error ignoring with valid-record output preservation.
- HTSlib `test_view.c`-style BGZF-compressed SAM view and record limiting.
- HTSlib `test_view.c`-style parsed SAM-to-BGZF-SAM compressed write output with record limiting.
- HTSlib `test_view.c`-style SAM-to-FASTQ and SAM-to-FASTA view output with record limiting.
- HTSlib `test_view.c`-style BAM no-compression output using noodles BGZF compression-level selection.
- HTSlib `test_view.c`-style SAM-to-BAM write with associated BAI index output.
- HTSlib `test_view.c`-style BAM whole-file, indexed region, and multi-region view output with raw BAM header order preservation, whole-file and region record limiting, and HTSlib-style multi-region de-duplication.
- HTSlib `test_view.c`-style generated long BAM record round trip for CIGAR and sequence data crossing BGZF block boundaries.
- HTSlib `test_convert_padded_header`-style BAM raw-header NUL padding read/view behavior.
- HTSlib `test_view.c`-style CRAM whole-file, indexed region, and multi-region view output with raw CRAM header preservation, whole-file and region record limiting, reference-backed MD/NM regeneration, and HTSlib-style multi-region de-duplication.
- HTSlib `test_view.c`-style SAM-to-CRAM write with associated CRAI index output.
- SAM-to-CRAM write/decode summary parity for `ce#1.sam`, `ce#2.sam`, `ce#1000.sam`, and `xx#pair.sam` with local FASTA reference repositories.
- HTSlib `test_view.c`-style whole-file BCF-to-VCF view output with raw BCF VCF header order preservation.
- HTSlib `test_view.c`-style VCF and BCF whole-file and indexed-region variant record limiting with header preservation.
- HTSlib `test_view.c`-style indexed VCF region header output for nonzero contig `IDX` fixtures.
- HTSlib `test_view.c`-style VCF compressed write with associated TBI index output.
- HTSlib `test_view.c`-style VCF-to-BCF write with associated CSI index output.
- HTSlib `test_view.c`-style benchmark/no-output record-reading counts for SAM, BAM, CRAM, VCF, and BCF.
- `test_view.c` C harness output routing, generic `hts_opt` option-string plumbing, shared thread-pool assignment, and CLI option parsing are out of scope for the Rust-only/no-CLI target.
- HTSlib `test/test.pl` VCF canonical-output coverage for `test-vcf-hdr-in.vcf`, `formatcols.vcf`, `noroundtrip.vcf`, `formatmissing.vcf`, `vcf_meta_meta.vcf`, and `vcf44_1.vcf`, including lenient structured-header spacing normalization, duplicate FORMAT key normalization, VCF 4.4 implicit/explicit phasing normalization, and HTSlib-style placeholder FORMAT/sample columns for sample-bearing records with no explicit FORMAT values.
- Local BGZF-compressed SAM indexed-query helpers backed by `noodles-sam` and HTSlib-style associated BAI/CSI index lookup.
- Local BAM region-query helpers backed by `noodles-bam` and HTSlib-style associated BAI/CSI index lookup.
- HTSlib-style CRAM TLEN auto-creation tie-break behavior for paired and triplet read fixtures from `test/tlen/tlen.sh`.
- HTSlib-style SAM MM/ML base-modification reporting for explicit/implicit status, extended metadata output, unchecked-base output, probabilities, double-strand calls, ChEBI IDs, canonical `N` skip semantics, MN validation, bounds rejection, missing-MM reset behavior, and the base-modification pileup golden outputs.
- HTSlib-style SAM filter expression evaluation over SAM records for upstream `filter.sh` integer, aux-tag, string, function, and CIGAR-metric cases.
- HTSlib `test/pileup.c`-style pileup output for deletion, insertion, deletion/insertion, refskip, pad CIGAR, overlap-removal SAM fixtures, and the BAM `small.bam` edge fixture.
- Local VCF and BCF header reading and record-count helpers backed by `noodles-vcf` and `noodles-bcf`.
- HTSlib-style VCF header get/remove helpers for FILTER, INFO, FORMAT, contig, structured generic, and unstructured generic records, plus header `Number` classification including VCF 4.5 FORMAT `P`, `LA`, `LG`, `LR`, and `M`.
- HTSlib-style typed VCF INFO value retrieval for integer, float, and missing float cases from `test-vcf-api.c`.
- HTSlib-style typed VCF FORMAT integer retrieval for scalar, vector, and missing values from `test-vcf-api.c`.
- HTSlib-style BCF FORMAT float missing and vector-end normalization from `test-vcf-api.c`.
- HTSlib-style VCF record span (`bcf1_t::rlen`) calculation over the upstream VCF 4.3, 4.4, and 4.5 table cases, including invalid `END` fallback behavior and mutation-triggered recalculation after allele, INFO, FORMAT, and copy operations.
- HTSlib-style VCF alternate allele removal for the upstream `bcf_remove_allele_set` cases, including INFO/FORMAT `Number=A/R/G/LA/LR/LG` trimming and genotype remapping.
- HTSlib-style VCF record serialization output for `test-vcf-api.c` header-removal, record sync/duplication, and INFO/FORMAT mutation cases.
- HTSlib-style VCF/BCF forward/backward sweep helpers for position and FORMAT/PL checksum behavior.
- Initial HTSlib-style synced VCF no-index summary, VCF-output, and BCF-output helpers with contig-order, record-order, and weird chromosome-name region/target filtering validation.
- Rust reader-addition synced VCF summary helper equivalent to the C `test-bcf-sr --usefptr` path for the Rust-only API.
- HTSlib-style indexed BCF weird chromosome-name region/target filtering validation for `test-bcf-sr` range fixtures.
- Rust-native synced variant group allele-pairing helper covering deterministic `test-bcf-sr` pairing modes for SNP, indel, combined SNP/indel, reference-compatible, exact, subset, and all-compatible rows, plus deterministic `test-bcf-sr.pl` randomized-shape fixtures with duplicate reader groups, multiallelic rows, shuffled group ordering, shuffled per-group variant/input ordering, and additional duplicate-reader input-order cases across multiple refs and variant mixes.
- `test-bcf-sr.c` CLI file-list/output routing and stochastic `test-bcf-sr.pl` execution are represented by Rust path/reader helpers and deterministic randomized-shape fixtures for the Rust-only/no-CLI target.
- Local VCF region-query helpers backed by `noodles-vcf` and HTSlib-style associated TBI/CSI index lookup.
- Local BCF CSI build and region-query helpers backed by `noodles-bcf` and HTSlib-style associated CSI index lookup.
- HTSlib-style BCF iterator creation by numeric reference ID and region string backed by `noodles-bcf` indexed queries.
- HTSlib-style BGZF virtual-offset helpers over `noodles-bgzf::VirtualPosition`.
- HTSlib-style region parsing with `HTS_PARSE_THOUSANDS_SEP`, `HTS_PARSE_ONE_COORD`, and `HTS_PARSE_LIST` equivalents.
- HTSlib-style filter expression parser/evaluator with arithmetic, bitwise, boolean, comparison, regex, null/default/exists, and math-function behavior from `test_expr.c`.
- HTSlib-style in-memory region index parsing, insertion, iteration, and overlap queries from `test-regidx.c`.
- HTSlib-style integer parsing and safe string printing from `test_str2int.c`.
- Item-level `hts_endian.h` classification with Rust helper coverage for unsigned, signed, float, double, and C-compatible floating-point alias names.
- Item-level `hts_defs.h` classification with C attribute/export macro scope decisions and `hts_prefetch` call-shape coverage.
- Item-level `hts_log.h` classification with log-level state coverage and `test-logging.pl` source-message style validation.
- Item-level `hts_os.h` classification with deterministic POSIX-compatible `rand48` helper coverage and C platform-shim decisions.
- Item-level `kfunc.h` classification with C-shaped `kf_lgamma` and `kt_fisher_exact` aliases and remaining statistical helpers marked for future implementation.
- Item-level `kbitset.h` classification with C-shaped safe Rust aliases for initialization, resize, mutation, lookup, and iteration.
- Item-level `klist.h` classification with Rust-owned list aliases for initialization, push-and-fill, FIFO shift, destruction, and iteration; C memory-pool internals are out of scope.
- Item-level `ksort.h` classification with C-shaped Rust aliases for sorting, heap, kth-element, string comparison, swap, and shuffle helpers.
- Item-level `khash.h` and `khash_str2int.h` classification with Rust `HashMap` wrappers and C-shaped aliases for the tested string-to-integer behavior.
- Item-level `kstring.h` classification with C-shaped Rust aliases for dynamic-string lifecycle, append, integer formatting, insertion, search, and chunked line-read behavior; tokenizer, split, and variadic formatting APIs remain classified for future work.
- Item-level `vcfutils.h` classification with small public A/C/G/T and genotype-index helpers implemented and larger allele/genotype utilities classified for follow-up.
- Item-level `vcf_sweep.h` classification with Rust-owned sweep state, header access, and forward/backward traversal mapped to the existing noodles-backed sweep adapter.
- Item-level `thread_pool.h` classification with C scheduler/process/result objects out of scope and noodles worker-count adapters retained for observable format-level threading behavior.
- Item-level `kseq.h` classification with C stream/parser macros mapped to Rust reader ownership, noodles FASTA/FASTQ APIs, and scoped FASTX conversion helpers.
- Item-level `knetfile.h` classification with deprecated legacy network APIs retained in the map while FTP/HTTP remote I/O and raw descriptor wrappers remain out of scope.
- Item-level `regidx.h` classification with Rust-owned region indexes, built-in parser coverage, overlap iteration, and remaining C-shaped sequence/count helpers marked for follow-up.
- Item-level `hts_expr.h` classification with Rust expression values, filter evaluation, closure-based symbol lookup, and deprecated C eval API tracking.
- Item-level `synced_bcf_reader.h` classification with Rust no-index summaries, reader-addition parity, weird chromosome region/target filtering, indexed BCF query coverage, and allele-pairing logic mapped against remaining C-shaped synced-reader APIs.
- Item-level `tbx.h` classification with noodles-backed TBI/CSI state, text-format build/query adapters, local index-loading behavior, and remaining C iterator/name-list helpers marked for follow-up.
- Item-level `faidx.h` classification with noodles-backed FASTA/FASTQ FAI state, build/load/query adapters, golden retrieval coverage, and C cache/thread/path helpers marked for follow-up or out of scope.
- Item-level `hfile.h` classification with local path extension and `data:` URL helpers implemented while raw C stream, remote, memory-buffer, and plugin APIs remain tracked or out of scope.
- Item-level `hts.h` classification with format detection, options, indexes, regions, iterators, binning, legacy helpers, realignment/error-model, MD5, CRC32, and endian APIs mapped to existing Rust coverage or follow-up work.
- Item-level `sam.h` classification with SAM/BAM/CRAM header, record, CIGAR, flag, I/O, index, iterator, aux, pileup, base-modification, and BAQ APIs mapped to noodles adapters or remaining Rust work.
- Item-level `vcf.h` classification with VCF/BCF header, record, typed value, genotype, mutation, variant classification, low-level BCF encoding, and index APIs mapped to noodles adapters or remaining Rust work.
- Item-level `bgzf.h` classification with noodles-backed BGZF I/O, virtual offsets, GZI, EOF/compression detection, worker-count concurrency, and low-level C handle/block APIs marked for follow-up or out of scope.
- Item-level `cram.h` classification with noodles-backed CRAM read/write/query/reference/index behavior mapped separately from low-level HTSlib block, container, codec, slice, and transcode APIs.
- Item-level `kroundup.h` classification with C-shaped `kroundup64`, `kroundup32`, and `kroundup_size_t` aliases.
- HTSlib-style `hts_os` deterministic POSIX `rand48` helpers, including global seeding, explicit seed-state updates, floating-point output, and 31-bit integer output.
- HTSlib-style `kbitset` operations over Rust storage, including init, resize, clear, insert-all, insert/delete, existence checks, and ascending iteration.
- HTSlib-style `khash_str2int` wrapper semantics over a Rust `HashMap`, including key existence, lookup, auto-incrementing insertion, explicit set, and size behavior.
- HTSlib-style `klist` FIFO behavior over Rust `VecDeque`, including push-by-value, push-and-fill, shift, size, and iteration behavior.
- HTSlib-style `kroundup` helpers for unsigned 32-bit, unsigned 64-bit, `usize`, and signed 32-bit allocation rounding with overflow saturation.
- HTSlib-style `ksort` behavior over Rust slices, including stable and unstable sorting adapters, heap helpers, kth-element selection, string sorting, and shuffle swaps.
- HTSlib-style VCF/BCF variant classification from `test-bcf_set_variant_type.c`.
- Test-status manifest for ported, failing, unported, and out-of-scope HTSlib tests.
- Temporary C-oracle checks have been replaced by direct Rust assertions and upstream golden fixtures for the selected test-parity scope.
- Local file format detection for selected HTSlib fixtures:
  - SAM
  - BGZF-compressed BAM
  - CRAM
  - VCF
  - BAI
  - BGZF-compressed CSI
  - FASTA
  - FASTQ
  - BED
- HTSlib-compatible `haddextension` path/URL extension behavior and `data:` URL decoding from `hfile.c`.
- HTSlib `test-logging.pl` source-message style validation for `hts_log_*` literals.
- Ported `test-parse-reg.c` built-in cases to Rust acceptance tests.
- Ported `test-regidx.c` deterministic and randomized overlap cases to Rust acceptance tests.
- Ported `test_str2int.c` integer-boundary and `hts_strprint` cases to Rust acceptance tests.
- Ported `test-bcf_set_variant_type.c` cases to Rust acceptance tests.
- Ported `test_realn.c` existing BAQ `BQ:Z`/`ZQ:Z` apply, revert, and no-op behavior for precomputed BAQ tags.
- Ported `sam_cap_mapq` MAPQ capping behavior for matched, mismatched, and clipped SAM records.
- Ported initial `probaln_glocal` behavior with HTSlib-oracle parity for exact, mismatch, and insertion cases.
- Ported default non-extended `test_realn.c` BAQ recalculation for `realn01` and `realn02` against upstream expected SAM output.
- Ported non-extended `test_realn.c` BAQ recalculation with quality application for `realn01` and `realn02` against upstream expected SAM output.
- Ported forced non-extended `test_realn.c` BAQ recomputation for `realn02-r.sam` against upstream expected SAM output.
- Ported extended `test_realn.c` BAQ recalculation for `realn01`, `realn02`, and adjacent-match `realn03` against upstream expected SAM output.
- Ported BAM records spanning BGZF block boundary fixture checks from `test/test.pl`.
- Ported generated long BAM record round-trip behavior from `test_view`.
- Ported the library-relevant `test_rebgzip` multi-block BGZF/GZI boundary and indexed-seek fixture checks from `test/test.pl`; byte-for-byte `bgzip -g` CLI reconstruction remains out of scope while CLI parity is deferred.
- Ported BAM view and indexed region view behavior for `test/test.pl` fixtures with references in the BAM binary reference table but no raw SAM `@SQ` header lines.
- Ported large-position BGZF SAM CSI iterator output from `test/test.pl`, including adaptive CSI depth for long references and multi-region de-duplication.
- Ported long-reference VCF CSI query fixtures from `test/test.pl` comments for `longrefs/index.vcf`, including adaptive CSI depth and `INFO/END` span overlap against `longrefs/index.expected1.vcf` and `longrefs/index.expected2.vcf`.
- Ported VCF canonical-output cases from `test/test.pl` `test_vcf_various` and `test_vcf_44` for `test-vcf-hdr-in.vcf`, `formatcols.vcf`, `noroundtrip.vcf`, `formatmissing.vcf`, `vcf_meta_meta.vcf`, and `vcf44_1.vcf`.
- Ported `test/test.pl` `index2.sam` mapped/unmapped pair BAM index query count behavior from `test_index`.
- `test/test.pl` plugin-loading, remote reference-cache server, `annot-tsv`, `bgzip`, `htsfile`, command-routing, and exhaustive CRAM `hts_opt` version/profile matrix checks are out of scope for the Rust-only/no-CLI/no-remote-I/O target; representative CRAM read/write/query/reference behavior is covered through noodles-backed helpers.
- Ported uncompressed FASTA/FASTQ and BGZF-compressed FASTA faidx golden-output checks from `test/faidx/faidx.tst`.
- Ported FASTQ/FASTA conversion golden-output checks from `test/fastq/fastq.tst`.

## Current Verification

- `cargo fmt --check`
- `cargo clippy --workspace --all-targets -- -D warnings`
- `cargo test --workspace`
