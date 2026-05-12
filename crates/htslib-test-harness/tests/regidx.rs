use htslib_rs::regidx::{Parser, RegionIndex, parse_region_line};

#[test]
fn ports_test_regidx_sequential_access() {
    let mut index = RegionIndex::new();

    for i in 0..10 {
        let pos = 10 * (i + 1);
        let line = format!("1\t{pos}\t{pos}\t{pos}");
        assert!(index.insert_line(&line, Parser::Tab).unwrap());
    }

    for (i, record) in index.iter().enumerate() {
        let expected = 10 * (i as i64 + 1);
        assert_eq!(record.start, expected - 1);
        assert_eq!(record.end, expected - 1);
        assert_eq!(
            record.payload.as_deref(),
            Some(expected.to_string().as_str())
        );
    }
}

#[test]
fn ports_test_regidx_tab_reg_and_bed_queries() {
    check_parser(Parser::Tab, |chr, start, end| {
        format!("{chr}\t{start}\t{end}\n")
    });
    check_parser(Parser::Region, |chr, start, end| {
        format!("{chr}:{start}-{end}\n")
    });
    check_parser(Parser::Bed, |chr, start, end| {
        format!("{chr}\t{}\t{end}\n", start - 1)
    });
}

fn check_parser(make_parser: Parser, line: impl Fn(&str, i64, i64) -> String) {
    let mut index = RegionIndex::new();

    for i in 1..10 {
        for (start, end) in [(10 * i, 10 * i), (10 * i + 1, 10 * i + 1)] {
            index
                .insert_line(&line("1", start, end), make_parser)
                .unwrap();
        }

        let start = 20000 * i;
        let end = start + 2000;
        index
            .insert_line(&line("1", start, end), make_parser)
            .unwrap();
    }

    for i in 1..10 {
        let start = 10 * i - 1;
        assert!(!index.has_overlap("1", start - 1, start - 1));

        let start = 10 * i;
        assert_eq!(index.overlaps("1", start - 1, start - 1).len(), 1);

        let start = 10 * i + 1;
        assert_eq!(index.overlaps("1", start - 1, start - 1).len(), 1);

        let start = 10 * i;
        let end = start + 1;
        assert_eq!(index.overlaps("1", start - 1, end - 1).len(), 2);

        let start = 20000 * i - 5000;
        let end = 20000 * i + 3000;
        assert_eq!(index.overlaps("1", start - 1, end - 1).len(), 1);
    }
}

#[test]
fn ports_test_regidx_custom_payload_cases() {
    let mut index = RegionIndex::new();

    for line in [
        "1 10000000 10000000 1:10000000-10000000",
        "1 20000000 20000001 1:20000000-20000001",
        "1 20000002 20000002 1:20000002-20000002",
        "1 30000000 30000000 1:30000000-30000000",
        "1 8000000000 8000000000 1:8000000000-8000000000",
    ] {
        index.insert_line(line, Parser::Tab).unwrap();
    }

    let hits = index.overlaps("1", 9999999, 9999999);
    assert_eq!(hits[0].payload.as_deref(), Some("1:10000000-10000000"));
    assert!(index.has_overlap("1", 9999998, 9999999));
    assert!(index.has_overlap("1", 9999998, 10000002));
    assert!(!index.has_overlap("1", 9999998, 9999998));

    for pos in [20000000, 20000002, 30000000, 8000000000] {
        assert!(index.has_overlap("1", pos - 1, pos - 1));
    }

    let wrapped = 8000000000_i64 & 0xffffffff;
    assert!(!index.has_overlap("1", wrapped - 1, wrapped - 1));
}

#[test]
fn ports_test_regidx_explicit_past_case() {
    let mut index = RegionIndex::new();
    index
        .insert_line("12:2064519-2064763", Parser::Region)
        .unwrap();

    let query = parse_region_line("12:2064488-2067434").unwrap().unwrap();
    assert!(index.has_overlap(&query.seq, query.start, query.end));
}

#[test]
fn ports_test_regidx_random_overlap_cases() {
    let mut rng = Lcg::new(1);

    for _ in 0..100 {
        let mut index = RegionIndex::new();
        let (query_start, query_end) = random_region(&mut rng, 0, 999);
        let mut expected = 0;

        for _ in 0..50 {
            let (start, end) = random_region(&mut rng, 0, 999);
            let line = format!("1\t{}\t{}\t1:{}-{}", start + 1, end + 1, start + 1, end + 1);
            index.insert_line(&line, Parser::Tab).unwrap();

            if end >= query_start && start <= query_end {
                expected += 1;
            }
        }

        let hits = index.overlaps("1", query_start, query_end);
        assert_eq!(hits.len(), expected);

        for hit in hits {
            assert!(hit.end >= query_start && hit.start <= query_end);
            assert_eq!(
                hit.payload.as_deref(),
                Some(format!("1:{}-{}", hit.start + 1, hit.end + 1).as_str())
            );
        }
    }
}

struct Lcg(u64);

impl Lcg {
    fn new(seed: u64) -> Self {
        Self(seed)
    }

    fn next(&mut self) -> u32 {
        self.0 = self.0.wrapping_mul(1103515245).wrapping_add(12345);
        ((self.0 / 65536) % 32768) as u32
    }
}

fn random_region(rng: &mut Lcg, min: u32, max: u32) -> (i64, i64) {
    let start = min + (u64::from(rng.next()) * u64::from(max - min) / 32767) as u32;
    let end = start + (u64::from(rng.next()) * u64::from(max - start) / 32767) as u32;
    (i64::from(start), i64::from(end))
}
