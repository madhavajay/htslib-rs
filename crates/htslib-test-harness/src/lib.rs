//! HTSlib acceptance-test tracking.

/// Port status for an upstream HTSlib test.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Status {
    /// The test has a Rust equivalent and currently passes.
    Passing,
    /// The test has a Rust equivalent but currently fails.
    Failing,
    /// The test has not been ported yet.
    Unported,
    /// The test is explicitly out of scope for the current project target.
    OutOfScope,
}

/// A tracked HTSlib test or acceptance-test group.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Test {
    /// Upstream test name or local acceptance-test group.
    pub name: &'static str,
    /// Current port status.
    pub status: Status,
    /// Short note describing the scope.
    pub note: &'static str,
}

/// The current HTSlib test-port manifest.
pub const TESTS: &[Test] = &[
    Test {
        name: "format_detection.rs",
        status: Status::Passing,
        note: "Initial Rust acceptance tests over HTSlib fixtures.",
    },
    Test {
        name: "alignment_io.rs",
        status: Status::Passing,
        note: "SAM/BAM read, padded BAM header handling, SAM-to-CRAM write/decode across multiple local-reference fixtures, SAM.gz BAI/CSI and BAM/CRAM indexed query, associated-index fallback, CRAM header, CRAM TLEN, SAM filter integer counts, MM/ML base-modification, and pileup coverage over HTSlib fixtures.",
    },
    Test {
        name: "variant_io.rs",
        status: Status::Passing,
        note: "Initial VCF/BCF read, VCF/BCF indexed query, VCF header get/remove, VCF API serialization, VCF/BCF sweep, synced no-index and reader-addition summary, VCF various and VCF 4.4 header/output canonicalization, BCF weird chromosome-name range/target, and associated-index fallback coverage over HTSlib fixtures.",
    },
    Test {
        name: "bgzf_boundaries.rs",
        status: Status::Passing,
        note: "BAM records split across BGZF blocks are readable.",
    },
    Test {
        name: "bgzf.rs",
        status: Status::Passing,
        note: "BGZF read/write, worker-count read/write, plain/gzip/BGZF write modes, EOF marker, embedded EOF, GZI, multi-block GZI boundaries, seek, getline, and truncated-stream coverage.",
    },
    Test {
        name: "tabix_query.rs",
        status: Status::Passing,
        note: "BED, GFF, VCF TBI/CSI, separate-region output, and associated-index fallback query coverage over HTSlib fixtures.",
    },
    Test {
        name: "index_build.rs",
        status: Status::Passing,
        note: "BAM BAI/CSI, SAM.gz BAI/CSI, CRAM CRAI, VCF TBI/CSI, BCF CSI, and explicit CSI min_shift coverage for test_index.c library behavior.",
    },
    Test {
        name: "hfile_utils.rs",
        status: Status::Passing,
        note: "haddextension path/URL extension and data: URL decoding compatibility coverage.",
    },
    Test {
        name: "logging.rs",
        status: Status::Passing,
        note: "HTSlib log-message source style coverage from test-logging.pl.",
    },
    Test {
        name: "faidx.rs",
        status: Status::Passing,
        note: "FASTA/FASTQ FAI generation, uncompressed retrieval, and BGZF FASTA retrieval golden-output coverage.",
    },
    Test {
        name: "fastq_conversion.rs",
        status: Status::Passing,
        note: "Core FASTQ/FASTA to SAM and SAM to FASTQ/FASTA golden-output conversion coverage.",
    },
    Test {
        name: "test_bgzf.c",
        status: Status::Passing,
        note: "BGZF read/write, EOF, embedded EOF, GZI, virtual-position, getline, truncation, compression-mode, and Rust worker-count variants are covered by bgzf.rs; C open wrappers are out of scope for the Rust-only API target.",
    },
    Test {
        name: "test_faidx.c",
        status: Status::Passing,
        note: "FASTA/FASTQ faidx helper behavior covered by faidx.rs.",
    },
    Test {
        name: "test_index.c",
        status: Status::Passing,
        note: "Covered by index_build.rs and alignment_io.rs for BAM/SAM/CRAM/VCF/BCF library index behavior.",
    },
    Test {
        name: "test_expr.c",
        status: Status::Passing,
        note: "HTS expression parser/evaluator behavior covered by expr.rs.",
    },
    Test {
        name: "test_kfunc.c",
        status: Status::Passing,
        note: "Numeric helper coverage.",
    },
    Test {
        name: "test_khash.c",
        status: Status::Passing,
        note: "Hash utility compatibility coverage.",
    },
    Test {
        name: "test_kstring.c",
        status: Status::Passing,
        note: "Dynamic string compatibility coverage.",
    },
    Test {
        name: "test_realn.c",
        status: Status::Passing,
        note: "Covered by alignment_io.rs existing BAQ BQ/ZQ apply/revert/no-op behavior, sam_cap_mapq cases, probaln_glocal unit parity, and default/apply/forced/extended BAQ recalculation parity.",
    },
    Test {
        name: "test-regidx.c",
        status: Status::Passing,
        note: "Region index coverage.",
    },
    Test {
        name: "test-parse-reg.c",
        status: Status::Passing,
        note: "HTSlib region parser coverage.",
    },
    Test {
        name: "test_str2int.c",
        status: Status::Passing,
        note: "String-to-integer utility coverage.",
    },
    Test {
        name: "test_time_funcs.c",
        status: Status::Passing,
        note: "Time helper coverage.",
    },
    Test {
        name: "hts_endian.c",
        status: Status::Passing,
        note: "Endian conversion helper coverage.",
    },
    Test {
        name: "test_view.c",
        status: Status::Passing,
        note: "Covered by alignment_io.rs SAM/BGZF-SAM header-preserving record-limit, SAM parse-error ignoring, parsed SAM-to-BGZF-SAM compressed write, SAM-to-FASTQ/FASTA output, BAM no-compression output, SAM-to-BAM write with BAI output, padded BAM raw-header handling, BAM whole-file and indexed region/multi-region view with record limits, generated long BAM record round trip, CRAM whole-file/indexed region/multi-region view with record limits, SAM-to-CRAM write with CRAI output, and SAM/BAM/CRAM benchmark counts, plus variant_io.rs BCF-to-VCF, VCF/BCF whole-file and indexed-region record-limit, indexed VCF region header output, VCF compressed write with TBI output, VCF-to-BCF write with CSI output, and VCF/BCF benchmark counts; C harness output routing, generic hts_opt strings, shared thread-pool assignment, and CLI option plumbing are out of scope for the Rust-only/no-CLI target.",
    },
    Test {
        name: "test-vcf-api.c",
        status: Status::Passing,
        note: "Covered by variant_io.rs header get/remove, header Number classification, typed INFO values, FORMAT integer and BCF float vector-end values, vcf_open_mode, rlen table, rlen mutation recalculation, invalid END, bcf_remove_allele_set, BCF iterator creation, and record serialization output cases.",
    },
    Test {
        name: "test-vcf-sweep.c",
        status: Status::Passing,
        note: "VCF/BCF forward/backward sweep position and FORMAT/PL checksum behavior covered by variant_io.rs.",
    },
    Test {
        name: "test-bcf-sr.c",
        status: Status::Passing,
        note: "Covered by variant_io.rs no-index VCF summary, VCF output, BCF output, reader-addition summary, order-error, VCF weird chromosome-name range/target, and indexed BCF weird chromosome-name range/target cases, plus bcf_sr_pairing.rs deterministic allele-pairing modes, randomized-script-shape fixture, shuffled group-order stability, shuffled per-group variant/input-order stability, and additional duplicate-reader input-order cases across multiple refs and variant mixes. Stochastic test-bcf-sr.pl sort-loop parity is represented by deterministic randomized-shape fixtures; CLI file-list/output routing is out of scope for the Rust-only/no-CLI target.",
    },
    Test {
        name: "test-bcf-translate.c",
        status: Status::Passing,
        note: "Synthetic header merge and record translation fixture covered by variant_io.rs.",
    },
    Test {
        name: "test-bcf_set_variant_type.c",
        status: Status::Passing,
        note: "Variant classification coverage.",
    },
    Test {
        name: "test_mod.c",
        status: Status::Passing,
        note: "Covered by alignment_io.rs MM/ML reporting, extended metadata, ChEBI IDs, MN validation, bounds cases, and pileup fixture coverage.",
    },
    Test {
        name: "test_nibbles.c",
        status: Status::Passing,
        note: "SAM/BAM nibble encoding coverage.",
    },
    Test {
        name: "test_introspection.c",
        status: Status::OutOfScope,
        note: "C ABI/shared-library introspection is out of scope for the Rust-only API target.",
    },
    Test {
        name: "test/test.pl",
        status: Status::Passing,
        note: "Library-relevant scripted behavior is covered by bgzf.rs multi-block BGZF/GZI boundary checks, alignment_io.rs BAM block-boundary checks, no-raw-@SQ BAM view/query behavior, BAM/CRAM range view output, large-position BGZF SAM CSI iterator output, index2 mapped/unmapped pair BAM query counts, tabix_query.rs long-reference VCF CSI query behavior, and variant_io.rs VCF various/VCF 4.4 header/output canonicalization cases. Remaining plugin, remote ref-cache server, annot-tsv, bgzip, htsfile, command-routing, and CRAM hts_opt matrix checks are out of scope for the Rust-only/no-CLI/no-remote-I/O target.",
    },
    Test {
        name: "test-logging.pl",
        status: Status::Passing,
        note: "Source-message style checks for hts_log_* literals covered by logging.rs.",
    },
    Test {
        name: "test/faidx/test-faidx.sh",
        status: Status::Passing,
        note: "Scripted FAI cases covered by faidx.rs golden-output tests.",
    },
    Test {
        name: "test/fastq/test-fastq.sh",
        status: Status::Passing,
        note: "Scripted FASTQ/FASTA conversion cases covered by fastq_conversion.rs golden-output tests.",
    },
    Test {
        name: "test/tabix/test-tabix.sh",
        status: Status::Passing,
        note: "Library TBI/CSI VCF, large-coordinate CSI, BED, GFF, associated-index lookup, and separate-region output cases are covered by tabix_query.rs; tabix/bgzip CLI and thread-option behavior is out of scope for the Rust-only API target.",
    },
    Test {
        name: "test/mpileup/test-pileup.sh",
        status: Status::Passing,
        note: "Covered by alignment_io.rs SAM pileup fixtures for deletion, insertion, insertion/deletion, refskip, pads, overlap-removal, and the BAM small.bam edge fixture.",
    },
    Test {
        name: "test/base_mods/base-mods.sh",
        status: Status::Passing,
        note: "Covered by alignment_io.rs base-modification golden-output and rejection tests.",
    },
    Test {
        name: "test/sam_filter/filter.sh",
        status: Status::Passing,
        note: "Covered by alignment_io.rs integer, aux-tag, string golden-output, function, and CIGAR-metric filter cases.",
    },
    Test {
        name: "test/tlen/tlen.sh",
        status: Status::Passing,
        note: "CRAM TLEN auto-creation fixtures covered by alignment_io.rs against upstream SAM expected outputs.",
    },
];

/// Summary counts for the HTSlib test-port manifest.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct Summary {
    /// Passing tests.
    pub passing: usize,
    /// Ported but failing tests.
    pub failing: usize,
    /// Tests not ported yet.
    pub unported: usize,
    /// Tests intentionally out of scope for the current target.
    pub out_of_scope: usize,
}

impl Summary {
    /// Returns the total number of tracked tests.
    pub const fn total(self) -> usize {
        self.passing + self.failing + self.unported + self.out_of_scope
    }
}

/// Summarizes the current test-port manifest.
pub fn summarize(tests: &[Test]) -> Summary {
    let mut summary = Summary::default();

    for test in tests {
        match test.status {
            Status::Passing => summary.passing += 1,
            Status::Failing => summary.failing += 1,
            Status::Unported => summary.unported += 1,
            Status::OutOfScope => summary.out_of_scope += 1,
        }
    }

    summary
}
