//! Rust equivalents for selected `hts_defs.h` public helpers.

/// Compatibility wrapper for HTSlib's `hts_prefetch`.
///
/// The C helper is an optimization hint with no required semantic effect. The
/// Rust API keeps the call site expressible while leaving prefetch decisions to
/// the optimizer and target-specific code.
pub fn hts_prefetch<T: ?Sized>(_value: &T) {}
