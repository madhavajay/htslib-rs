use std::collections::HashSet;

use htslib_test_harness::{Status, TESTS, summarize};

#[test]
fn manifest_has_unique_test_names() {
    let mut names = HashSet::new();

    for test in TESTS {
        assert!(
            names.insert(test.name),
            "duplicate test name: {}",
            test.name
        );
    }
}

#[test]
fn manifest_reports_current_port_status() {
    let summary = summarize(TESTS);

    assert_eq!(summary.total(), TESTS.len());
    assert_eq!(summary.passing, 41);
    assert_eq!(summary.failing, 0);
    assert_eq!(summary.unported, 0);
    assert!(
        TESTS
            .iter()
            .any(|test| test.status == Status::OutOfScope && test.name == "test_introspection.c")
    );
}
