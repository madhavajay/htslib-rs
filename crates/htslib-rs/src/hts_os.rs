//! Rust replacements for selected HTSlib OS compatibility helpers.

use std::sync::Mutex;

const RAND48_SEED_0: u16 = 0x330e;
const RAND48_SEED_1: u16 = 0xabcd;
const RAND48_SEED_2: u16 = 0x1234;
const RAND48_MULT: u64 = 0x0005_deec_e66d;
const RAND48_ADD: u64 = 0x000b;
const RAND48_MASK: u64 = (1 << 48) - 1;

static RAND48: Mutex<Rand48> =
    Mutex::new(Rand48::new([RAND48_SEED_0, RAND48_SEED_1, RAND48_SEED_2]));

/// A POSIX `rand48` generator compatible with HTSlib's fallback implementation.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Rand48 {
    seed: [u16; 3],
}

impl Rand48 {
    /// Creates a generator from the three 16-bit seed words used by `erand48`.
    pub const fn new(seed: [u16; 3]) -> Self {
        Self { seed }
    }

    /// Creates a generator using `srand48` seed expansion.
    pub fn seeded(seed: i64) -> Self {
        Self {
            seed: [RAND48_SEED_0, seed as u16, ((seed as u64) >> 16) as u16],
        }
    }

    /// Returns the current seed words.
    pub fn seed(&self) -> [u16; 3] {
        self.seed
    }

    /// Advances this generator and returns a floating-point value in `[0.0, 1.0)`.
    pub fn erand48(&mut self) -> f64 {
        dorand48(&mut self.seed);

        f64::from(self.seed[0]) * 2f64.powi(-48)
            + f64::from(self.seed[1]) * 2f64.powi(-32)
            + f64::from(self.seed[2]) * 2f64.powi(-16)
    }

    /// Advances this generator and returns the high 31 bits.
    pub fn lrand48(&mut self) -> i64 {
        dorand48(&mut self.seed);

        (i64::from(self.seed[2]) << 15) + (i64::from(self.seed[1]) >> 1)
    }
}

/// Seeds the global `rand48` generator.
pub fn hts_srand48(seed: i64) {
    *RAND48.lock().expect("rand48 mutex poisoned") = Rand48::seeded(seed);
}

/// Advances the supplied seed words and returns a value in `[0.0, 1.0)`.
pub fn hts_erand48(seed: &mut [u16; 3]) -> f64 {
    let mut rand = Rand48::new(*seed);
    let value = rand.erand48();
    *seed = rand.seed();
    value
}

/// Advances the global generator and returns a value in `[0.0, 1.0)`.
pub fn hts_drand48() -> f64 {
    RAND48.lock().expect("rand48 mutex poisoned").erand48()
}

/// Advances the global generator and returns the high 31 bits.
pub fn hts_lrand48() -> i64 {
    RAND48.lock().expect("rand48 mutex poisoned").lrand48()
}

fn dorand48(seed: &mut [u16; 3]) {
    let state = u64::from(seed[0]) | (u64::from(seed[1]) << 16) | (u64::from(seed[2]) << 32);
    let next = state.wrapping_mul(RAND48_MULT).wrapping_add(RAND48_ADD) & RAND48_MASK;

    seed[0] = next as u16;
    seed[1] = (next >> 16) as u16;
    seed[2] = (next >> 32) as u16;
}

#[cfg(test)]
mod tests {
    use super::{Rand48, hts_drand48, hts_erand48, hts_lrand48, hts_srand48};

    #[test]
    fn test_rand48_seeded_sequence() {
        let mut rand = Rand48::seeded(0);

        assert_eq!(rand.erand48(), 0.17082803610628972);
        assert_eq!(rand.lrand48(), 1_610_402_240);

        let mut rand = Rand48::seeded(1);

        assert_eq!(rand.erand48(), 0.041630344771878214);
        assert_eq!(rand.lrand48(), 976_015_093);
    }

    #[test]
    fn test_global_and_explicit_seed() {
        hts_srand48(0);
        assert_eq!(hts_drand48(), 0.17082803610628972);
        assert_eq!(hts_lrand48(), 1_610_402_240);

        let mut seed = [0x330e, 0, 0];
        assert_eq!(hts_erand48(&mut seed), 0.17082803610628972);
        assert_eq!(seed, [0x5101, 0x62dc, 0x2bbb]);
    }
}
