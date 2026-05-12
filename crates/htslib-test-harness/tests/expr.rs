use htslib_rs::expr::{Filter, Value};

struct Case {
    truth: bool,
    number: f64,
    string: Option<&'static str>,
    expr: &'static str,
}

#[test]
fn ports_test_expr_cases() -> Result<(), Box<dyn std::error::Error>> {
    let cases = [
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "1",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "+1",
        },
        Case {
            truth: true,
            number: -1.0,
            string: None,
            expr: "-1",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: "!7",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "!0",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "!(!7)",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "!!7",
        },
        Case {
            truth: true,
            number: 5.0,
            string: None,
            expr: "2+3",
        },
        Case {
            truth: true,
            number: -1.0,
            string: None,
            expr: "2+-3",
        },
        Case {
            truth: true,
            number: 6.0,
            string: None,
            expr: "1+2+3",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "-2+3",
        },
        Case {
            truth: false,
            number: f64::NAN,
            string: None,
            expr: "1+null",
        },
        Case {
            truth: false,
            number: f64::NAN,
            string: None,
            expr: "null-1",
        },
        Case {
            truth: false,
            number: f64::NAN,
            string: None,
            expr: "-null",
        },
        Case {
            truth: true,
            number: 6.0,
            string: None,
            expr: "2*3",
        },
        Case {
            truth: true,
            number: 6.0,
            string: None,
            expr: "1*2*3",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: "2*0",
        },
        Case {
            truth: true,
            number: 7.0,
            string: None,
            expr: "(7)",
        },
        Case {
            truth: true,
            number: 7.0,
            string: None,
            expr: "((7))",
        },
        Case {
            truth: true,
            number: 21.0,
            string: None,
            expr: "(1+2)*(3+4)",
        },
        Case {
            truth: true,
            number: 14.0,
            string: None,
            expr: "(4*5)-(-2*-3)",
        },
        Case {
            truth: false,
            number: f64::NAN,
            string: None,
            expr: "2*null",
        },
        Case {
            truth: false,
            number: f64::NAN,
            string: None,
            expr: "null/2",
        },
        Case {
            truth: false,
            number: f64::NAN,
            string: None,
            expr: "0/0",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "(1+2)*3==9",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "(1+2)*3!=8",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: "(1+2)*3!=9",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: "(1+2)*3==8",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: "1>2",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "1<2",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: "3<3",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: "3>3",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "9<=9",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "9>=9",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "2*4==8",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "16==0x10",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "15<0x10",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "17>0x10",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: "2*4!=8",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "4+2<3+4",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: "4*2<3+4",
        },
        Case {
            truth: true,
            number: 8.0,
            string: None,
            expr: "4*(2<3)+4",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "(1<2) == (3>2)",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "1<2 == 3>2",
        },
        Case {
            truth: false,
            number: f64::NAN,
            string: None,
            expr: "null <= 0",
        },
        Case {
            truth: false,
            number: f64::NAN,
            string: None,
            expr: "null >= 0",
        },
        Case {
            truth: false,
            number: f64::NAN,
            string: None,
            expr: "null < 0",
        },
        Case {
            truth: false,
            number: f64::NAN,
            string: None,
            expr: "null > 0",
        },
        Case {
            truth: false,
            number: f64::NAN,
            string: None,
            expr: "null == null",
        },
        Case {
            truth: false,
            number: f64::NAN,
            string: None,
            expr: "null != null",
        },
        Case {
            truth: false,
            number: f64::NAN,
            string: None,
            expr: "null < 10",
        },
        Case {
            truth: false,
            number: f64::NAN,
            string: None,
            expr: "10 > null",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "2 && 1",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: "2 && 0",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: "0 && 2",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "2 || 1",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "2 || 0",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "0 || 2",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "1 || 2 && 3",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "2 && 3 || 1",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "0 && 3 || 2",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: "0 && 3 || 0",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: " 5 - 5 && 1",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: "+5 - 5 && 1",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: "null && 1",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: "1 && null",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "!null && 1",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "1 && !null",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "1 && null-but-true",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: "null || 0",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: "0 || null",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "!null || 0",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "0 || !null",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "0 || null-but-true",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "null || 1",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "1 || null",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "3 & 1",
        },
        Case {
            truth: true,
            number: 2.0,
            string: None,
            expr: "3 & 2",
        },
        Case {
            truth: true,
            number: 3.0,
            string: None,
            expr: "1 | 2",
        },
        Case {
            truth: true,
            number: 3.0,
            string: None,
            expr: "1 | 3",
        },
        Case {
            truth: true,
            number: 7.0,
            string: None,
            expr: "1 | 6",
        },
        Case {
            truth: true,
            number: 2.0,
            string: None,
            expr: "1 ^ 3",
        },
        Case {
            truth: false,
            number: f64::NAN,
            string: None,
            expr: "1 | null",
        },
        Case {
            truth: false,
            number: f64::NAN,
            string: None,
            expr: "null | 1",
        },
        Case {
            truth: false,
            number: f64::NAN,
            string: None,
            expr: "1 & null",
        },
        Case {
            truth: false,
            number: f64::NAN,
            string: None,
            expr: "null & 1",
        },
        Case {
            truth: false,
            number: f64::NAN,
            string: None,
            expr: "0 ^ null",
        },
        Case {
            truth: false,
            number: f64::NAN,
            string: None,
            expr: "null ^ 0",
        },
        Case {
            truth: false,
            number: f64::NAN,
            string: None,
            expr: "1 ^ null",
        },
        Case {
            truth: false,
            number: f64::NAN,
            string: None,
            expr: "null ^ 1",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "(1^0)&(4^3)",
        },
        Case {
            truth: true,
            number: 2.0,
            string: None,
            expr: "1 ^(0&4)^ 3",
        },
        Case {
            truth: true,
            number: 2.0,
            string: None,
            expr: "1 ^ 0&4 ^ 3",
        },
        Case {
            truth: true,
            number: 6.0,
            string: None,
            expr: "(1|0)^(4|3)",
        },
        Case {
            truth: true,
            number: 7.0,
            string: None,
            expr: "1 |(0^4)| 3",
        },
        Case {
            truth: true,
            number: 7.0,
            string: None,
            expr: "1 | 0^4 | 3",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "4 & 2 || 1",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "(4 & 2) || 1",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: "4 & (2 || 1)",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "1 || 4 & 2",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "1 || (4 & 2)",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: "(1 || 4) & 2",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: " (2*3)&7  > 4",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: " (2*3)&(7 > 4)",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "((2*3)&7) > 4",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "((2*3)&7) > 4 && 2*2 <= 4",
        },
        Case {
            truth: true,
            number: 1.0,
            string: Some("plugh"),
            expr: "magic",
        },
        Case {
            truth: true,
            number: 1.0,
            string: Some(""),
            expr: "empty",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "magic == \"plugh\"",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "magic != \"xyzzy\"",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "\"abc\" < \"def\"",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "\"abc\" <= \"abc\"",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: "\"abc\" < \"ab\"",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: "\"abc\" <= \"ab\"",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: "\"abc\" > \"def\"",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "\"abc\" >= \"abc\"",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "\"abc\" > \"ab\"",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "\"abc\" >= \"ab\"",
        },
        Case {
            truth: false,
            number: f64::NAN,
            string: None,
            expr: "null == \"x\"",
        },
        Case {
            truth: false,
            number: f64::NAN,
            string: None,
            expr: "null != \"x\"",
        },
        Case {
            truth: false,
            number: f64::NAN,
            string: None,
            expr: "null < \"x\"",
        },
        Case {
            truth: false,
            number: f64::NAN,
            string: None,
            expr: "null > \"x\"",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "\"abbc\" =~ \"^a+b+c+$\"",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: "\"aBBc\" =~ \"^a+b+c+$\"",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "\"aBBc\" !~ \"^a+b+c+$\"",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "\"xyzzy plugh abracadabra\" =~ magic",
        },
        Case {
            truth: true,
            number: 1.0,
            string: Some(""),
            expr: "empty-but-true",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: "!empty-but-true",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "!!empty-but-true",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "1 && empty-but-true && 1",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: "1 && empty-but-true && 0",
        },
        Case {
            truth: false,
            number: f64::NAN,
            string: None,
            expr: "null",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "!null",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: "!!null",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: "!\"foo\"",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "!!\"foo\"",
        },
        Case {
            truth: true,
            number: f64::NAN,
            string: None,
            expr: "null-but-true",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: "!null-but-true",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "!!null-but-true",
        },
        Case {
            truth: true,
            number: 0.0,
            string: None,
            expr: "zero-but-true",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: "!zero-but-true",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "!!zero-but-true",
        },
        Case {
            truth: true,
            number: 2.0_f64.ln(),
            string: None,
            expr: "log(2)",
        },
        Case {
            truth: true,
            number: 9.0_f64.exp(),
            string: None,
            expr: "exp(9)",
        },
        Case {
            truth: true,
            number: 9.0,
            string: None,
            expr: "log(exp(9))",
        },
        Case {
            truth: true,
            number: 8.0,
            string: None,
            expr: "pow(2,3)",
        },
        Case {
            truth: true,
            number: 3.0,
            string: None,
            expr: "sqrt(9)",
        },
        Case {
            truth: false,
            number: f64::NAN,
            string: None,
            expr: "sqrt(-9)",
        },
        Case {
            truth: true,
            number: 2.0,
            string: None,
            expr: "default(2,3)",
        },
        Case {
            truth: true,
            number: 3.0,
            string: None,
            expr: "default(null,3)",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: "default(null,0)",
        },
        Case {
            truth: true,
            number: f64::NAN,
            string: None,
            expr: "default(null-but-true,0)",
        },
        Case {
            truth: true,
            number: f64::NAN,
            string: None,
            expr: "default(null-but-true,null)",
        },
        Case {
            truth: true,
            number: f64::NAN,
            string: None,
            expr: "default(null,null-but-true)",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "exists(\"foo\")",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "exists(12)",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "exists(\"\")",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "exists(0)",
        },
        Case {
            truth: false,
            number: 0.0,
            string: None,
            expr: "exists(null)",
        },
        Case {
            truth: true,
            number: 1.0,
            string: None,
            expr: "exists(null-but-true)",
        },
    ];

    for case in cases {
        let actual = Filter::new(case.expr).eval_with(test_lookup)?;

        assert_eq!(actual.is_true(), case.truth, "truth for {}", case.expr);
        assert_float_eq(actual.number_value(), case.number, case.expr);
        assert_eq!(actual.as_str(), case.string, "string for {}", case.expr);
    }

    Ok(())
}

fn test_lookup(src: &str) -> Option<(Value, usize)> {
    if src.starts_with("foo") {
        Some((Value::number(15551.0), 3))
    } else if src.starts_with('a') {
        Some((Value::number(1.0), 1))
    } else if src.starts_with('b') {
        Some((Value::number(2.0), 1))
    } else if src.starts_with('c') {
        Some((Value::number(3.0), 1))
    } else if src.starts_with("magic") {
        Some((Value::string("plugh"), 5))
    } else if src.starts_with("empty-but-true") {
        Some((Value::string("").with_true(), 14))
    } else if src.starts_with("empty") {
        Some((Value::string(""), 5))
    } else if src.starts_with("zero-but-true") {
        Some((Value::number(0.0).with_true(), 13))
    } else if src.starts_with("null-but-true") {
        Some((Value::undefined().with_true(), 13))
    } else if src.starts_with("null") || src.starts_with("nan") {
        Some((
            Value::undefined(),
            if src.starts_with("null") { 4 } else { 3 },
        ))
    } else {
        None
    }
}

fn assert_float_eq(actual: f64, expected: f64, expr: &str) {
    if actual.is_nan() && expected.is_nan() {
        return;
    }

    assert!(
        (actual - expected).abs() <= f64::EPSILON * expected.abs().max(1.0),
        "number for {expr}: expected {expected}, got {actual}"
    );
}
