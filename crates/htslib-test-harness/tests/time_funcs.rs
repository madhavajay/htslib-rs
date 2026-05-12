use htslib_rs::time::{BrokenDownTime, TimeError, time_gm};

fn check_time_gm(
    year: i32,
    mon: i32,
    mday: i32,
    hour: i32,
    min: i32,
    sec: i32,
    expected: Result<i64, TimeError>,
) {
    let mut utc = BrokenDownTime {
        sec,
        min,
        hour,
        mday,
        mon: mon - 1,
        year: year - 1900,
    };

    assert_eq!(time_gm(&mut utc), expected);
}

#[test]
fn ports_test_time_funcs_specific_cases() {
    check_time_gm(2022, 6, 14, 12, 32, 10, Ok(1655209930));
    check_time_gm(1993, 9, 10514, 12, 32, 10, Ok(1655209930));
    check_time_gm(2020, 2, 28, 12, 0, 0, Ok(1582891200));
    check_time_gm(2020, 2, 29, 12, 0, 0, Ok(1582977600));
    check_time_gm(2020, 2, 30, 12, 0, 0, Ok(1583064000));
    check_time_gm(2020, 3, 0, 12, 0, 0, Ok(1582977600));
    check_time_gm(2019, 14, 1, 12, 0, 0, Ok(1580558400));
    check_time_gm(2019, 15, 1, 12, 0, 0, Ok(1583064000));
    check_time_gm(2019, 27, 1, 12, 0, 0, Ok(1614600000));
    check_time_gm(2019, 62, 1, 12, 0, 0, Ok(1706788800));
    check_time_gm(2019, 63, 1, 12, 0, 0, Ok(1709294400));
    check_time_gm(2021, 0, 31, 23, 59, 59, Ok(1609459199));
    check_time_gm(2021, -9, 1, 12, 0, 0, Ok(1583064000));
    check_time_gm(2021, -10, 1, 12, 0, 0, Ok(1580558400));
    check_time_gm(2021, -22, 1, 12, 0, 0, Ok(1549022400));
    check_time_gm(1970, 1, 1, 0, 0, 0, Ok(0));
    check_time_gm(1970, 1, 1, 0, 0, i32::MAX, Ok(i64::from(i32::MAX)));
    check_time_gm(2038, 1, 19, 3, 14, 7, Ok(i64::from(i32::MAX)));
    check_time_gm(2038, 1, 19, 3, 14, 8, Ok(i64::from(i32::MAX) + 1));
}

#[test]
fn rejects_dates_before_epoch() {
    check_time_gm(1969, 12, 31, 23, 59, 59, Err(TimeError::Overflow));
}
