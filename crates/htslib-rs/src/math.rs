//! Numeric helpers compatible with selected HTSlib `kfunc` behavior.

/// Result of a Fisher exact test.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FisherExact {
    /// Probability of the observed contingency table.
    pub probability: f64,
    /// Left-tail probability.
    pub left: f64,
    /// Right-tail probability.
    pub right: f64,
    /// Two-tail probability.
    pub two_tail: f64,
}

#[derive(Clone, Copy, Debug)]
struct HypergeoAccumulator {
    n11: i32,
    n1_: i32,
    n_1: i32,
    n: i32,
    p: f64,
}

/// Calculates Fisher's exact test for a 2x2 table.
///
/// This mirrors HTSlib's `kt_fisher_exact`, including its underflow shortcut
/// for very large tables.
pub fn fisher_exact(n11: i32, n12: i32, n21: i32, n22: i32) -> FisherExact {
    let n1_ = n11 + n12;
    let n_1 = n11 + n21;
    let n = n11 + n12 + n21 + n22;
    let max = n_1.min(n1_);
    let min = 0.max(n1_ + n_1 - n);

    if min == max {
        return FisherExact {
            probability: 1.0,
            left: 1.0,
            right: 1.0,
            two_tail: 1.0,
        };
    }

    let mut aux = HypergeoAccumulator {
        n11: 0,
        n1_: 0,
        n_1: 0,
        n: 0,
        p: 0.0,
    };
    let q = hypergeo_acc(n11, Some((n1_, n_1, n)), &mut aux);

    if q == 0.0 {
        let peak_is_right =
            i64::from(n11) * (i64::from(n) + 2) < (i64::from(n_1) + 1) * (i64::from(n1_) + 1);

        return if peak_is_right {
            FisherExact {
                probability: 0.0,
                left: 0.0,
                right: 1.0,
                two_tail: 0.0,
            }
        } else {
            FisherExact {
                probability: 0.0,
                left: 1.0,
                right: 0.0,
                two_tail: 0.0,
            }
        };
    }

    let mut p = hypergeo_acc(min, None, &mut aux);
    let mut left = 0.0;
    let mut i = min + 1;

    while p < 0.99999999 * q && i <= max {
        left += p;
        p = hypergeo_acc(i, None, &mut aux);
        i += 1;
    }

    i -= 1;
    if p < 1.00000001 * q {
        left += p;
    } else {
        i -= 1;
    }

    p = hypergeo_acc(max, None, &mut aux);
    let mut right = 0.0;
    let mut j = max - 1;

    while p < 0.99999999 * q && j >= 0 {
        right += p;
        p = hypergeo_acc(j, None, &mut aux);
        j -= 1;
    }

    j += 1;
    if p < 1.00000001 * q {
        right += p;
    } else {
        j += 1;
    }

    let mut two_tail = left + right;
    if two_tail > 1.0 {
        two_tail = 1.0;
    }

    if (i - n11).abs() < (j - n11).abs() {
        right = 1.0 - left + q;
    } else {
        left = 1.0 - right + q;
    }

    FisherExact {
        probability: q,
        left,
        right,
        two_tail,
    }
}

/// C-shaped alias for HTSlib's `kt_fisher_exact`.
pub fn kt_fisher_exact(n11: i32, n12: i32, n21: i32, n22: i32) -> FisherExact {
    fisher_exact(n11, n12, n21, n22)
}

fn hypergeo_acc(n11: i32, params: Option<(i32, i32, i32)>, aux: &mut HypergeoAccumulator) -> f64 {
    if let Some((n1_, n_1, n)) = params {
        aux.n11 = n11;
        aux.n1_ = n1_;
        aux.n_1 = n_1;
        aux.n = n;
    } else {
        if n11 % 11 != 0 && n11 + aux.n - aux.n1_ - aux.n_1 != 0 {
            if n11 == aux.n11 + 1 {
                aux.p *= f64::from(aux.n1_ - aux.n11) / f64::from(n11)
                    * f64::from(aux.n_1 - aux.n11)
                    / f64::from(n11 + aux.n - aux.n1_ - aux.n_1);
                aux.n11 = n11;
                return aux.p;
            }

            if n11 == aux.n11 - 1 {
                aux.p *= f64::from(aux.n11) / f64::from(aux.n1_ - n11)
                    * f64::from(aux.n11 + aux.n - aux.n1_ - aux.n_1)
                    / f64::from(aux.n_1 - n11);
                aux.n11 = n11;
                return aux.p;
            }
        }

        aux.n11 = n11;
    }

    aux.p = hypergeo(aux.n11, aux.n1_, aux.n_1, aux.n);
    aux.p
}

fn hypergeo(n11: i32, n1_: i32, n_1: i32, n: i32) -> f64 {
    (lbinom(n1_, n11) + lbinom(n - n1_, n_1 - n11) - lbinom(n, n_1)).exp()
}

fn lbinom(n: i32, k: i32) -> f64 {
    if k == 0 || n == k {
        0.0
    } else {
        ln_gamma(f64::from(n + 1)) - ln_gamma(f64::from(k + 1)) - ln_gamma(f64::from(n - k + 1))
    }
}

fn ln_gamma(z: f64) -> f64 {
    const COEFFS: [f64; 9] = [
        0.9999999999998099,
        676.5203681218851,
        -1259.1392167224028,
        771.3234287776531,
        -176.6150291621406,
        12.507343278686905,
        -0.13857109526572012,
        9.984369578019572e-6,
        1.5056327351493116e-7,
    ];

    let z = z - 1.0;
    let mut x = COEFFS[0];

    for (i, coefficient) in COEFFS.iter().enumerate().skip(1) {
        x += coefficient / (z + i as f64);
    }

    let t = z + 7.5;
    0.5 * (2.0 * std::f64::consts::PI).ln() + (z + 0.5) * t.ln() - t + x.ln()
}

/// C-shaped alias for HTSlib's `kf_lgamma`.
pub fn kf_lgamma(z: f64) -> f64 {
    ln_gamma(z)
}

#[cfg(test)]
mod tests {
    use super::{fisher_exact, kf_lgamma, kt_fisher_exact};

    #[test]
    fn test_fisher_exact() {
        let actual = fisher_exact(2, 1, 0, 31);

        assert!((actual.probability - 0.005347593583).abs() <= 1e-8);
        assert!((actual.left - 1.0).abs() <= 1e-8);
        assert!((actual.right - 0.005347593583).abs() <= 1e-8);
        assert!((actual.two_tail - 0.005347593583).abs() <= 1e-8);

        assert_eq!(kt_fisher_exact(2, 1, 0, 31), actual);
    }

    #[test]
    fn test_lgamma_alias() {
        assert!((kf_lgamma(5.0) - 24.0_f64.ln()).abs() <= 1e-12);
    }
}
