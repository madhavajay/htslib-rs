use htslib_rs::math::{fisher_exact, kf_lgamma, kt_fisher_exact};

fn assert_close(actual: f64, expected: f64, label: &str, table: (i32, i32, i32, i32)) {
    assert!(
        (actual - expected).abs() <= 1e-8,
        "{label} mismatch for [{}, {} | {}, {}]: {actual} != {expected}",
        table.0,
        table.1,
        table.2,
        table.3
    );
}

fn check_fisher(
    table: (i32, i32, i32, i32),
    expected_left: f64,
    expected_right: f64,
    expected_two_tail: f64,
    expected_probability: f64,
) {
    let actual = fisher_exact(table.0, table.1, table.2, table.3);
    let alias_actual = kt_fisher_exact(table.0, table.1, table.2, table.3);

    assert_eq!(alias_actual, actual);

    assert_close(actual.left, expected_left, "left", table);
    assert_close(actual.right, expected_right, "right", table);
    assert_close(actual.two_tail, expected_two_tail, "two-tail", table);
    assert_close(
        actual.probability,
        expected_probability,
        "probability",
        table,
    );
}

#[test]
fn ports_kfunc_lgamma_alias() {
    assert_close(kf_lgamma(1.0), 0.0, "lgamma(1)", (0, 0, 0, 0));
    assert_close(kf_lgamma(5.0), 24.0_f64.ln(), "lgamma(5)", (0, 0, 0, 0));
}

#[test]
fn ports_test_kfunc_fisher_exact_cases() {
    check_fisher(
        (2, 1, 0, 31),
        1.0,
        0.005347593583,
        0.005347593583,
        0.005347593583,
    );
    check_fisher((2, 1, 0, 1), 1.0, 0.5, 1.0, 0.5);
    check_fisher((3, 1, 0, 0), 1.0, 1.0, 1.0, 1.0);
    check_fisher(
        (3, 15, 37, 45),
        0.021479750169,
        0.995659202564,
        0.033161943699,
        0.017138952733,
    );
    check_fisher(
        (12, 5, 29, 2),
        0.044554737835,
        0.994525206022,
        0.080268552074,
        0.039079943857,
    );

    check_fisher((781, 23171, 4963, 2455001), 1.0, 0.0, 0.0, 0.0);
    check_fisher((333, 381, 801722, 7664285), 1.0, 0.0, 0.0, 0.0);
    check_fisher((4155, 4903, 805463, 8507517), 1.0, 0.0, 0.0, 0.0);
    check_fisher((4455, 4903, 805463, 8507517), 1.0, 0.0, 0.0, 0.0);
    check_fisher((5455, 4903, 805463, 8507517), 1.0, 0.0, 0.0, 0.0);

    check_fisher(
        (1, 1, 100000, 1000000),
        0.991735477166,
        0.173555146661,
        0.173555146661,
        0.165290623827,
    );
    check_fisher((1000, 1000, 100000, 1000000), 1.0, 0.0, 0.0, 0.0);
    check_fisher((1000, 1000, 1000000, 100000), 0.0, 1.0, 0.0, 0.0);

    check_fisher((49999, 10001, 90001, 49999), 1.0, 0.0, 0.0, 0.0);
    check_fisher((50000, 10000, 90000, 50000), 1.0, 0.0, 0.0, 0.0);
    check_fisher((50001, 9999, 89999, 50001), 1.0, 0.0, 0.0, 0.0);
    check_fisher((10000, 50000, 130000, 10000), 0.0, 1.0, 0.0, 0.0);
}
