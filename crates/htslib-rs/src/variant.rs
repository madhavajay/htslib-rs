//! HTSlib-compatible VCF/BCF variant classification helpers.

/// Variant type bitmask used by HTSlib's VCF/BCF APIs.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct VariantType(u32);

impl VariantType {
    /// Reference allele or non-variant symbolic allele.
    pub const REF: Self = Self(0);
    /// Single nucleotide polymorphism.
    pub const SNP: Self = Self(1 << 0);
    /// Multi-nucleotide polymorphism.
    pub const MNP: Self = Self(1 << 1);
    /// Insertion or deletion.
    pub const INDEL: Self = Self(1 << 2);
    /// Other non-SNP/MNP/INDEL variant.
    pub const OTHER: Self = Self(1 << 3);
    /// Breakend.
    pub const BND: Self = Self(1 << 4);
    /// Overlapping deletion allele, `ALT=*`.
    pub const OVERLAP: Self = Self(1 << 5);
    /// Insertion, always accompanied by [`Self::INDEL`].
    pub const INS: Self = Self(1 << 6);
    /// Deletion, always accompanied by [`Self::INDEL`].
    pub const DEL: Self = Self(1 << 7);
    /// Any variant type, excluding [`Self::REF`].
    pub const ANY: Self = Self(
        Self::SNP.0
            | Self::MNP.0
            | Self::INDEL.0
            | Self::OTHER.0
            | Self::BND.0
            | Self::OVERLAP.0
            | Self::INS.0
            | Self::DEL.0,
    );

    /// Returns the raw HTSlib-compatible bitmask.
    pub const fn bits(self) -> u32 {
        self.0
    }

    /// Returns whether all bits in `other` are present.
    pub const fn contains(self, other: Self) -> bool {
        (self.0 & other.0) == other.0
    }
}

impl std::ops::BitOr for VariantType {
    type Output = Self;

    fn bitor(self, rhs: Self) -> Self::Output {
        Self(self.0 | rhs.0)
    }
}

impl std::ops::BitOrAssign for VariantType {
    fn bitor_assign(&mut self, rhs: Self) {
        self.0 |= rhs.0;
    }
}

impl std::ops::BitAnd for VariantType {
    type Output = Self;

    fn bitand(self, rhs: Self) -> Self::Output {
        Self(self.0 & rhs.0)
    }
}

/// Classification for a single reference and alternate allele pair.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Variant {
    /// Variant type.
    pub variant_type: VariantType,
    /// Number of affected bases, negative for deletions.
    ///
    /// HTSlib does not initialize this field for all classifications, e.g.,
    /// breakends and symbolic alleles. Those cases are represented as `None`.
    pub len: Option<i32>,
}

impl Variant {
    const fn new(variant_type: VariantType, len: Option<i32>) -> Self {
        Self { variant_type, len }
    }
}

/// Classifies a VCF reference and alternate allele using HTSlib rules.
///
/// This ports `bcf_set_variant_type` from HTSlib and intentionally keeps its
/// allele-string behavior instead of using a normalized variant model.
pub fn classify_variant(ref_allele: &str, alt_allele: &str) -> Variant {
    let ref_allele = ref_allele.as_bytes();
    let alt_allele = alt_allele.as_bytes();

    if alt_allele == b"*" {
        return Variant::new(VariantType::OVERLAP, Some(0));
    }

    if ref_allele.len() == 1 && alt_allele.len() == 1 {
        if alt_allele == b"." || ref_allele[0] == alt_allele[0] || alt_allele == b"X" {
            return Variant::new(VariantType::REF, Some(0));
        }

        return Variant::new(VariantType::SNP, Some(1));
    }

    if alt_allele.first() == Some(&b'<') {
        if alt_allele == b"<X>" || alt_allele == b"<*>" || alt_allele == b"<NON_REF>" {
            return Variant::new(VariantType::REF, Some(0));
        }

        return Variant::new(VariantType::OTHER, None);
    }

    if matches!(alt_allele.first(), Some(b']' | b'[')) {
        return Variant::new(VariantType::BND, None);
    }

    let mut r = 0;
    let mut a = 0;

    while r < ref_allele.len()
        && a < alt_allele.len()
        && ref_allele[r].eq_ignore_ascii_case(&alt_allele[a])
    {
        r += 1;
        a += 1;
    }

    if a < alt_allele.len() && r == ref_allele.len() {
        if matches!(alt_allele.last(), Some(b']' | b'[')) {
            return Variant::new(VariantType::BND, None);
        }

        let len = alt_allele.len() as i32 - a as i32;
        return Variant::new(VariantType::INDEL | VariantType::INS, Some(len));
    } else if r < ref_allele.len() && a == alt_allele.len() {
        let len = a as i32 - ref_allele.len() as i32;
        return Variant::new(VariantType::INDEL | VariantType::DEL, Some(len));
    } else if r == ref_allele.len() && a == alt_allele.len() {
        return Variant::new(VariantType::REF, Some(0));
    }

    let mut re = ref_allele.len() - 1;
    let mut ae = alt_allele.len() - 1;

    if matches!(alt_allele[ae], b']' | b'[') {
        return Variant::new(VariantType::BND, None);
    }

    while re > r && ae > a && ref_allele[re].eq_ignore_ascii_case(&alt_allele[ae]) {
        re -= 1;
        ae -= 1;
    }

    if ae == a {
        if re == r {
            return Variant::new(VariantType::SNP, Some(1));
        }

        let len = -((re - r) as i32);
        let variant_type = if ref_allele[re].eq_ignore_ascii_case(&alt_allele[ae]) {
            VariantType::INDEL | VariantType::DEL
        } else {
            VariantType::OTHER
        };

        return Variant::new(variant_type, Some(len));
    } else if re == r {
        let len = (ae - a) as i32;
        let variant_type = if ref_allele[re].eq_ignore_ascii_case(&alt_allele[ae]) {
            VariantType::INDEL | VariantType::INS
        } else {
            VariantType::OTHER
        };

        return Variant::new(variant_type, Some(len));
    }

    let variant_type = if re - r == ae - a {
        VariantType::MNP
    } else {
        VariantType::OTHER
    };

    let len = if re - r > ae - a {
        -((re - r + 1) as i32)
    } else {
        (ae - a + 1) as i32
    };

    Variant::new(variant_type, Some(len))
}

/// Converts an A/C/G/T base to its HTSlib `bcf_acgt2int` code.
pub fn bcf_acgt2int(base: char) -> Option<u8> {
    match base.to_ascii_uppercase() {
        'A' => Some(0),
        'C' => Some(1),
        'G' => Some(2),
        'T' => Some(3),
        _ => None,
    }
}

/// Converts an HTSlib A/C/G/T code to a base.
pub fn bcf_int2acgt(index: usize) -> Option<char> {
    [Some('A'), Some('C'), Some('G'), Some('T')]
        .get(index)
        .copied()
        .flatten()
}

/// Converts diploid allele indexes to a VCF `Number=G` index.
pub const fn bcf_ij2g(i: usize, j: usize) -> usize {
    j * (j + 1) / 2 + i
}

#[cfg(test)]
mod tests {
    use super::{VariantType, bcf_acgt2int, bcf_ij2g, bcf_int2acgt, classify_variant};

    #[test]
    fn test_classify_variant() {
        assert_eq!(classify_variant("A", "T").variant_type, VariantType::SNP);
        assert_eq!(
            classify_variant("A", "AA").variant_type,
            VariantType::INDEL | VariantType::INS
        );
        assert_eq!(
            classify_variant("AA", "A").variant_type,
            VariantType::INDEL | VariantType::DEL
        );
        assert_eq!(classify_variant("AA", "TT").variant_type, VariantType::MNP);
        assert_eq!(
            classify_variant("A", "*").variant_type,
            VariantType::OVERLAP
        );
        assert_eq!(classify_variant("A", ".").variant_type, VariantType::REF);
    }

    #[test]
    fn test_vcfutils_small_helpers() {
        assert_eq!(bcf_acgt2int('A'), Some(0));
        assert_eq!(bcf_acgt2int('c'), Some(1));
        assert_eq!(bcf_acgt2int('G'), Some(2));
        assert_eq!(bcf_acgt2int('t'), Some(3));
        assert_eq!(bcf_acgt2int('N'), None);

        assert_eq!(bcf_int2acgt(0), Some('A'));
        assert_eq!(bcf_int2acgt(1), Some('C'));
        assert_eq!(bcf_int2acgt(2), Some('G'));
        assert_eq!(bcf_int2acgt(3), Some('T'));
        assert_eq!(bcf_int2acgt(4), None);

        assert_eq!(bcf_ij2g(0, 0), 0);
        assert_eq!(bcf_ij2g(0, 1), 1);
        assert_eq!(bcf_ij2g(1, 1), 2);
        assert_eq!(bcf_ij2g(0, 2), 3);
        assert_eq!(bcf_ij2g(1, 2), 4);
        assert_eq!(bcf_ij2g(2, 2), 5);
    }
}
