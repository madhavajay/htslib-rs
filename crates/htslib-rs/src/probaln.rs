//! Probabilistic banded glocal alignment compatible with HTSlib.

use std::io;

const EI: f64 = 0.25;
const EM: f64 = 0.33333333333;

/// Alignment parameters for `probaln_glocal`.
#[derive(Clone, Copy, Debug)]
pub struct ProbalnParams {
    /// Gap-open probability.
    pub d: f64,
    /// Gap-extension probability.
    pub e: f64,
    /// Band width.
    pub bw: usize,
}

impl ProbalnParams {
    /// Returns the HTSlib short-read defaults.
    pub fn illumina() -> Self {
        Self {
            d: 0.001,
            e: 0.1,
            bw: 10,
        }
    }
}

/// Result of a probabilistic glocal alignment.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ProbalnResult {
    /// Phred-scaled likelihood score.
    pub likelihood: i32,
    /// MAP state for each query base.
    pub state: Option<Vec<i32>>,
    /// Phred-scaled posterior probability that each MAP state is wrong.
    pub posterior_quality: Option<Vec<u8>>,
}

/// Performs probabilistic banded glocal alignment.
///
/// `reference` and `query` are encoded as `0, 1, 2, 3, 4` for
/// `A, C, G, T, N`. If `base_quality` is `None`, quality 30 is used.
pub fn probaln_glocal(
    reference: &[u8],
    query: &[u8],
    base_quality: Option<&[u8]>,
    params: ProbalnParams,
    with_map: bool,
) -> io::Result<ProbalnResult> {
    let l_ref = reference.len();
    let l_query = query.len();

    if base_quality.is_some_and(|q| q.len() != l_query) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "base quality length differs from query length",
        ));
    }

    if l_ref == 0 || l_query == 0 {
        return Ok(ProbalnResult {
            likelihood: 0,
            state: with_map.then(Vec::new),
            posterior_quality: with_map.then(Vec::new),
        });
    }

    let mut bw = l_ref.max(l_query).min(params.bw);
    bw = bw.max(l_ref.abs_diff(l_query));
    let bw2 = bw * 2 + 1;
    let i_dim = if bw2 < l_ref {
        bw2 * 3 + 6
    } else {
        l_ref * 3 + 6
    };

    let mut f = vec![0.0; (l_query + 1) * i_dim];
    let mut b = with_map.then(|| vec![0.0; (l_query + 1) * i_dim]);
    let mut s = vec![0.0; l_query + 2];
    let mut qual = Vec::with_capacity(l_query);

    for i in 0..l_query {
        let q = base_quality.map_or(30, |qualities| qualities[i]);
        qual.push(10f64.powf(-(f64::from(q)) / 10.0));
    }

    let s_m = 1.0 / (2.0 * l_query as f64 + 2.0);
    let s_i = s_m;
    let mut m = [0.0; 9];

    m[0] = (1.0 - params.d - params.d) * (1.0 - s_m);
    m[1] = params.d * (1.0 - s_m);
    m[2] = params.d * (1.0 - s_m);
    m[3] = (1.0 - params.e) * (1.0 - s_i);
    m[4] = params.e * (1.0 - s_i);
    m[5] = 0.0;
    m[6] = 1.0 - params.e;
    m[7] = 0.0;
    m[8] = params.e;

    let b_m = (1.0 - params.d) / l_ref as f64;
    let b_i = params.d / l_ref as f64;

    let k0 = band_index(bw, 0, 0);
    f[k0] = 1.0;
    s[0] = 1.0;

    {
        let end = l_ref.min(bw + 1);
        let mut sum = 0.0;

        for k in 1..=end {
            let u = band_index(bw, 1, k);
            let e = emission(reference[k - 1], query[0], qual[0]);
            f[i_dim + u] = e * b_m;
            f[i_dim + u + 1] = EI * b_i;
            sum += f[i_dim + u] + f[i_dim + u + 1];
        }

        s[1] = sum;
    }

    for i in 2..=l_query {
        let beg = 1.max(i.saturating_sub(bw));
        let end = l_ref.min(i + bw);
        let prev_scale = 1.0 / s[i - 1];
        let qli = qual[i - 1];
        let qyi = query[i - 1];
        let mut sum = 0.0;

        for k in beg..=end {
            let u = band_index(bw, i, k);
            let v11 = band_index(bw, i - 1, k - 1);
            let v10 = band_index(bw, i - 1, k);
            let v01 = band_index(bw, i, k - 1);
            let e = emission(reference[k - 1], qyi, qli);
            let row = i * i_dim;
            let prev = (i - 1) * i_dim;

            f[row + u] = e
                * (m[0] * prev_scale * f[prev + v11]
                    + m[3] * prev_scale * f[prev + v11 + 1]
                    + m[6] * prev_scale * f[prev + v11 + 2]);
            f[row + u + 1] =
                EI * (m[1] * prev_scale * f[prev + v10] + m[4] * prev_scale * f[prev + v10 + 1]);
            f[row + u + 2] = m[2] * f[row + v01] + m[8] * f[row + v01 + 2];

            sum += f[row + u] + f[row + u + 1] + f[row + u + 2];
        }

        s[i] = sum;
    }

    {
        let prev_scale = 1.0 / s[l_query];
        let mut sum = 0.0;

        for k in 1..=l_ref {
            let u = band_index(bw, l_query, k);

            if u < 3 || u >= i_dim {
                continue;
            }

            let row = l_query * i_dim;
            sum += prev_scale * f[row + u] * s_m + prev_scale * f[row + u + 1] * s_i;
        }

        s[l_query + 1] = sum;
    }

    let likelihood = likelihood_from_scales(&s, l_ref, l_query);

    let Some(ref mut b) = b else {
        return Ok(ProbalnResult {
            likelihood,
            state: None,
            posterior_quality: None,
        });
    };

    for k in 1..=l_ref {
        let u = band_index(bw, l_query, k);

        if u < 3 || u >= i_dim {
            continue;
        }

        let row = l_query * i_dim;
        b[row + u] = s_m / s[l_query] / s[l_query + 1];
        b[row + u + 1] = s_i / s[l_query] / s[l_query + 1];
    }

    for i in (1..l_query).rev() {
        let beg = 1.max(i.saturating_sub(bw));
        let end = l_ref.min(i + bw);
        let y = if i > 1 { 1.0 } else { 0.0 };
        let qli1 = qual[i];
        let qyi1 = query[i];

        for k in (beg..=end).rev() {
            let u = band_index(bw, i, k);
            let v11 = band_index(bw, i + 1, k + 1);
            let v10 = band_index(bw, i + 1, k);
            let v01 = band_index(bw, i, k + 1);
            let row = i * i_dim;
            let next = (i + 1) * i_dim;
            let e = if k >= l_ref {
                0.0
            } else {
                emission(reference[k], qyi1, qli1) * b[next + v11]
            };

            b[row + u] = e * m[0] + EI * m[1] * b[next + v10 + 1] + m[2] * b[row + v01 + 2];
            b[row + u + 1] = e * m[3] + EI * m[4] * b[next + v10 + 1];
            b[row + u + 2] = (e * m[6] + m[8] * b[row + v01 + 2]) * y;
        }

        let beg_u = band_index(bw, i, beg);
        let end_u = band_index(bw, i, end) + 2;
        let scale = 1.0 / s[i];

        for value in &mut b[i * i_dim + beg_u..=i * i_dim + end_u] {
            *value *= scale;
        }
    }

    {
        let end = l_ref.min(bw + 1);
        let mut sum = 0.0;

        for k in (1..=end).rev() {
            let u = band_index(bw, 1, k);

            if u < 3 || u >= i_dim {
                continue;
            }

            let e = emission(reference[k - 1], query[0], qual[0]);
            sum += e * b[i_dim + u] * b_m + EI * b[i_dim + u + 1] * b_i;
        }

        let u = band_index(bw, 0, 0);
        b[u] = sum / s[0];
    }

    let mut state = Vec::with_capacity(l_query);
    let mut posterior_quality = Vec::with_capacity(l_query);

    for (i, scale_i) in s.iter().copied().enumerate().take(l_query + 1).skip(1) {
        let beg = 1.max(i.saturating_sub(bw));
        let end = l_ref.min(i + bw);
        let row = i * i_dim;
        let scale = 1.0 / scale_i;
        let mut sum = 0.0;
        let mut max = 0.0;
        let mut max_k = -1;

        for k in beg..=end {
            let u = band_index(bw, i, k);
            let z1 = scale * f[row + u] * b[row + u];
            let z2 = scale * f[row + u + 1] * b[row + u + 1];

            if z1 > max {
                max = z1;
                max_k = (k as i32 - 1) << 2;
            }

            if z2 > max {
                max = z2;
                max_k = ((k as i32 - 1) << 2) | 1;
            }

            sum += z1 + z2;
        }

        let posterior = max / sum;
        let q = (-4.343 * (1.0 - posterior).ln() + 0.499) as i32;

        state.push(max_k);
        posterior_quality.push(q.min(99) as u8);
    }

    Ok(ProbalnResult {
        likelihood,
        state: Some(state),
        posterior_quality: Some(posterior_quality),
    })
}

fn band_index(bw: usize, i: usize, k: usize) -> usize {
    let x = i.saturating_sub(bw);

    k.saturating_add(1).saturating_sub(x) * 3
}

fn emission(reference: u8, query: u8, quality_error_probability: f64) -> f64 {
    if reference > 3 || query > 3 {
        1.0
    } else if reference == query {
        1.0 - quality_error_probability
    } else {
        quality_error_probability * EM
    }
}

fn likelihood_from_scales(scales: &[f64], l_ref: usize, l_query: usize) -> i32 {
    let mut p = 1.0;
    let mut phred = 0.0;

    for scale in scales {
        p *= scale;

        if p < 1e-100 {
            phred += -4.343 * p.ln();
            p = 1.0;
        }
    }

    phred += -4.343 * (p * l_ref as f64 * l_query as f64).ln();

    (phred + 0.499) as i32
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bases(s: &str) -> Vec<u8> {
        s.bytes()
            .map(|b| match b.to_ascii_uppercase() {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 4,
            })
            .collect()
    }

    #[test]
    fn aligns_exact_match_without_map() -> io::Result<()> {
        let result = probaln_glocal(
            &bases("ACGT"),
            &bases("ACGT"),
            None,
            ProbalnParams::illumina(),
            false,
        )?;

        assert_eq!(result.likelihood, 5);
        assert_eq!(result.state, None);
        assert_eq!(result.posterior_quality, None);

        Ok(())
    }

    #[test]
    fn aligns_with_map_outputs() -> io::Result<()> {
        let result = probaln_glocal(
            &bases("ACGT"),
            &bases("ACGT"),
            None,
            ProbalnParams::illumina(),
            true,
        )?;

        assert_eq!(result.likelihood, 5);
        assert_eq!(result.state, Some(vec![0, 4, 8, 12]));
        assert_eq!(result.posterior_quality, Some(vec![36, 52, 52, 36]));

        Ok(())
    }

    #[test]
    fn aligns_mismatch_like_htslib() -> io::Result<()> {
        let result = probaln_glocal(
            &bases("ACGT"),
            &bases("AGGT"),
            None,
            ProbalnParams::illumina(),
            true,
        )?;

        assert_eq!(result.likelihood, 40);
        assert_eq!(result.state, Some(vec![0, 4, 8, 12]));
        assert_eq!(result.posterior_quality, Some(vec![17, 17, 29, 31]));

        Ok(())
    }

    #[test]
    fn aligns_insertion_like_htslib() -> io::Result<()> {
        let result = probaln_glocal(
            &bases("ACGT"),
            &bases("ACGTT"),
            None,
            ProbalnParams::illumina(),
            true,
        )?;

        assert_eq!(result.likelihood, 39);
        assert_eq!(result.state, Some(vec![0, 4, 8, 12, 13]));
        assert_eq!(result.posterior_quality, Some(vec![36, 46, 33, 3, 3]));

        Ok(())
    }
}
