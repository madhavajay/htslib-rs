//! HTSlib-compatible UTC time helpers.

/// A broken-down UTC time using the same field conventions as C `struct tm`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BrokenDownTime {
    /// Seconds after the minute.
    pub sec: i32,
    /// Minutes after the hour.
    pub min: i32,
    /// Hours after midnight.
    pub hour: i32,
    /// Day of the month, normally 1 through 31.
    pub mday: i32,
    /// Months since January, normally 0 through 11.
    pub mon: i32,
    /// Years since 1900.
    pub year: i32,
}

/// Error returned when a broken-down time cannot be represented.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TimeError {
    /// Normalization or conversion would overflow, or the result is before the Unix epoch.
    Overflow,
}

/// Returns whether a full calendar year is a leap year.
pub fn year_is_leap(year: i64) -> bool {
    (year % 4 == 0 && year % 100 != 0) || year % 400 == 0
}

/// Returns the number of leap years before the start of `year`.
///
/// This mirrors HTSlib and is defined for `year >= 1`.
pub fn leaps_to_year_start(year: i64) -> i64 {
    let year = year - 1;

    year / 4 - year / 100 + year / 400
}

/// Converts a broken-down UTC time to seconds since 1970-01-01T00:00:00Z.
///
/// The input is normalized in place before conversion, matching HTSlib's
/// `hts_time_gm` behavior.
pub fn time_gm(target: &mut BrokenDownTime) -> Result<i64, TimeError> {
    const MONTH_START: [[i32; 12]; 2] = [
        [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334],
        [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335],
    ];

    normalise_tm(target)?;

    if target.year < 70 {
        return Err(TimeError::Overflow);
    }

    let years_from_epoch = i64::from(target.year - 70);
    let leaps = leaps_to_year_start(i64::from(target.year) + 1900) - leaps_to_year_start(1970);
    let year_days = 365 * (years_from_epoch - leaps) + 366 * leaps;
    let leap = usize::from(year_is_leap(i64::from(target.year) + 1900));
    let days =
        year_days + i64::from(MONTH_START[leap][target.mon as usize]) + i64::from(target.mday - 1);

    Ok(days * 86400
        + i64::from(target.hour) * 3600
        + i64::from(target.min) * 60
        + i64::from(target.sec))
}

fn normalise_tm(t: &mut BrokenDownTime) -> Result<(), TimeError> {
    const DAYS_PER_MON: [[i32; 12]; 2] = [
        [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
        [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
    ];
    const YEAR_DAYS: [i32; 2] = [365, 366];

    if t.sec > 62 {
        normalise(&mut t.min, &mut t.sec, 60)?;
    }

    normalise(&mut t.hour, &mut t.min, 60)?;
    normalise(&mut t.mday, &mut t.hour, 24)?;
    normalise(&mut t.year, &mut t.mon, 12)?;

    let mut year = i64::from(t.year) + 1900;

    while t.mday <= 0 {
        year -= 1;
        let leap = usize::from(year_is_leap(year + i64::from(t.mon > 1)));
        t.mday += YEAR_DAYS[leap];
    }

    while t.mday > 366 {
        let leap = usize::from(year_is_leap(year + i64::from(t.mon > 1)));
        t.mday -= YEAR_DAYS[leap];
        year += 1;
    }

    loop {
        let leap = usize::from(year_is_leap(year));
        let mdays = DAYS_PER_MON[leap][t.mon as usize];

        if t.mday <= mdays {
            break;
        }

        t.mday -= mdays;
        t.mon += 1;

        if t.mon >= 12 {
            year += 1;
            t.mon = 0;
        }
    }

    let normalized_year = year - 1900;

    if normalized_year != i64::from(t.year) {
        t.year = i32::try_from(normalized_year).map_err(|_| TimeError::Overflow)?;
    }

    Ok(())
}

fn normalise(tens: &mut i32, units: &mut i32, base: i32) -> Result<(), TimeError> {
    if *units < 0 || *units >= base {
        let delta = units.div_euclid(base);
        let normalized_tens = i64::from(*tens) + i64::from(delta);

        *tens = i32::try_from(normalized_tens).map_err(|_| TimeError::Overflow)?;
        *units -= delta * base;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{BrokenDownTime, time_gm, year_is_leap};

    #[test]
    fn test_year_is_leap() {
        assert!(year_is_leap(2020));
        assert!(!year_is_leap(2022));
        assert!(!year_is_leap(1900));
        assert!(year_is_leap(2000));
    }

    #[test]
    fn test_time_gm_normalized() {
        let mut utc = BrokenDownTime {
            sec: 10,
            min: 32,
            hour: 12,
            mday: 14,
            mon: 5,
            year: 122,
        };

        assert_eq!(time_gm(&mut utc), Ok(1655209930));
    }
}
