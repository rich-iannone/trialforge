//! Calendar dates for the study clock.
//!
//! A [`Date`] is a count of days since 1970-01-01 (the Arrow `Date32` and R `Date` epoch), so it
//! crosses into both bindings without conversion. Civil-date arithmetic uses Howard Hinnant's
//! algorithms, which are exact over the whole `i32` range and need no time-zone data.

use core::fmt;

/// A calendar date, stored as days since 1970-01-01.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Date(i32);

/// Days from the SAS epoch (1960-01-01) to the Unix epoch (1970-01-01).
const SAS_EPOCH_OFFSET: i32 = 3653;

const MONTHS: [&str; 12] = [
    "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC",
];

impl Date {
    /// The date `days` after 1970-01-01.
    pub const fn from_days(days: i32) -> Self {
        Self(days)
    }

    /// Days since 1970-01-01.
    pub const fn days(self) -> i32 {
        self.0
    }

    /// Days since 1960-01-01, the value SAS stores for a date.
    pub const fn sas_days(self) -> i32 {
        self.0 + SAS_EPOCH_OFFSET
    }

    /// The date for a year, month (1–12) and day (1–31). Returns `None` for invalid dates.
    pub fn from_ymd(year: i32, month: u32, day: u32) -> Option<Self> {
        if !(1..=12).contains(&month) || day == 0 || day > days_in_month(year, month) {
            return None;
        }
        let y = if month <= 2 { year - 1 } else { year };
        let era = y.div_euclid(400);
        let yoe = y - era * 400;
        let mp = (month + 9) % 12;
        let doy = (153 * mp + 2) / 5 + day - 1;
        let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy as i32;
        Some(Self(era * 146_097 + doe - 719_468))
    }

    /// The (year, month, day) for this date.
    pub fn ymd(self) -> (i32, u32, u32) {
        let z = self.0 + 719_468;
        let era = z.div_euclid(146_097);
        let doe = z - era * 146_097;
        let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
        let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
        let mp = (5 * doy + 2) / 153;
        let day = (doy - (153 * mp + 2) / 5 + 1) as u32;
        let month = if mp < 10 { mp + 3 } else { mp - 9 } as u32;
        let year = yoe + era * 400 + i32::from(month <= 2);
        (year, month, day)
    }

    /// Parses a complete ISO 8601 date (`YYYY-MM-DD`). Partial dates are not dates.
    pub fn parse_iso(s: &str) -> Option<Self> {
        let b = s.as_bytes();
        if b.len() != 10 || b[4] != b'-' || b[7] != b'-' {
            return None;
        }
        let year = s[0..4].parse().ok()?;
        let month = s[5..7].parse().ok()?;
        let day = s[8..10].parse().ok()?;
        Self::from_ymd(year, month, day)
    }

    /// The date `n` days later (or earlier, for negative `n`).
    ///
    /// # Panics
    /// If the result is outside the `i32` day range (a bug in the caller, never valid data).
    #[must_use]
    pub const fn add_days(self, n: i32) -> Self {
        match self.0.checked_add(n) {
            Some(d) => Self(d),
            None => panic!("date arithmetic overflowed"),
        }
    }

    /// Days from `other` to `self`.
    pub const fn days_since(self, other: Self) -> i32 {
        self.0 - other.0
    }

    /// `ddMMMyy:hh:mm:ss` at midnight, as used in SAS transport headers.
    pub fn sas_header_datetime(self) -> String {
        let (y, m, d) = self.ymd();
        format!(
            "{d:02}{}{:02}:00:00:00",
            MONTHS[(m - 1) as usize],
            y.rem_euclid(100)
        )
    }
}

impl fmt::Display for Date {
    /// ISO 8601 `YYYY-MM-DD`.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (y, m, d) = self.ymd();
        write!(f, "{y:04}-{m:02}-{d:02}")
    }
}

/// A date and time of day, stored as seconds since 1970-01-01T00:00:00 (no time zone: the
/// study clock is local time throughout, as in most SDTM data).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DateTime(i64);

/// Seconds per day.
const DAY_SECONDS: i64 = 86_400;

impl DateTime {
    /// `date` at `hour:minute`.
    ///
    /// # Panics
    /// If `hour` > 23 or `minute` > 59.
    pub fn new(date: Date, hour: u32, minute: u32) -> Self {
        assert!(hour < 24 && minute < 60, "invalid time of day");
        Self(i64::from(date.days()) * DAY_SECONDS + i64::from(hour * 3600 + minute * 60))
    }

    /// The date-time `seconds` after 1970-01-01T00:00:00.
    pub const fn from_seconds(seconds: i64) -> Self {
        Self(seconds)
    }

    /// Seconds since 1970-01-01T00:00:00 (the Arrow timestamp and R `POSIXct` epoch).
    pub const fn seconds(self) -> i64 {
        self.0
    }

    /// Seconds since 1960-01-01T00:00:00, the value SAS stores for a datetime.
    pub const fn sas_seconds(self) -> i64 {
        self.0 + SAS_EPOCH_OFFSET as i64 * DAY_SECONDS
    }

    /// The calendar date.
    pub fn date(self) -> Date {
        Date::from_days(self.0.div_euclid(DAY_SECONDS) as i32)
    }

    /// Seconds since midnight.
    pub fn time_seconds(self) -> i64 {
        self.0.rem_euclid(DAY_SECONDS)
    }

    /// The date-time `minutes` later (or earlier, for negative `minutes`).
    #[must_use]
    pub const fn add_minutes(self, minutes: i64) -> Self {
        Self(self.0 + minutes * 60)
    }

    /// Hours from `other` to `self`.
    pub fn hours_since(self, other: Self) -> f64 {
        (self.0 - other.0) as f64 / 3600.0
    }

    /// ISO 8601 with seconds, `YYYY-MM-DDThh:mm:ss`.
    pub fn iso_seconds(self) -> String {
        let t = self.time_seconds();
        format!(
            "{}T{:02}:{:02}:{:02}",
            self.date(),
            t / 3600,
            t % 3600 / 60,
            t % 60
        )
    }

    /// Parses ISO 8601 `YYYY-MM-DDThh:mm` or `YYYY-MM-DDThh:mm:ss`.
    pub fn parse_iso(s: &str) -> Option<Self> {
        let (d, t) = s.split_once('T')?;
        let date = Date::parse_iso(d)?;
        let mut parts = t.split(':');
        let h: u32 = parts.next()?.parse().ok()?;
        let m: u32 = parts.next()?.parse().ok()?;
        let sec: u32 = parts.next().map_or(Some(0), |x| x.parse().ok())?;
        if parts.next().is_some() || h > 23 || m > 59 || sec > 59 || t.len() < 5 {
            return None;
        }
        Some(Self(
            i64::from(date.days()) * DAY_SECONDS + i64::from(h * 3600 + m * 60 + sec),
        ))
    }
}

impl fmt::Display for DateTime {
    /// ISO 8601 `YYYY-MM-DDThh:mm`, with `:ss` only when the seconds are not zero.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let t = self.time_seconds();
        let (h, m, s) = (t / 3600, t % 3600 / 60, t % 60);
        write!(f, "{}T{h:02}:{m:02}", self.date())?;
        if s != 0 {
            write!(f, ":{s:02}")?;
        }
        Ok(())
    }
}

/// SDTM study day (`--DY`) of `date` relative to `reference` (usually RFSTDTC). There is no
/// day 0: the reference date is day 1 and the day before it is day -1.
pub fn study_day(date: Date, reference: Date) -> i32 {
    let diff = date.days_since(reference);
    if diff >= 0 { diff + 1 } else { diff }
}

fn is_leap(year: i32) -> bool {
    (year % 4 == 0 && year % 100 != 0) || year % 400 == 0
}

fn days_in_month(year: i32, month: u32) -> u32 {
    match month {
        2 if is_leap(year) => 29,
        2 => 28,
        4 | 6 | 9 | 11 => 30,
        _ => 31,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn datetimes() {
        let d = Date::from_ymd(2024, 3, 5).unwrap();
        let t = DateTime::new(d, 8, 30);
        assert_eq!(t.to_string(), "2024-03-05T08:30");
        assert_eq!(DateTime::parse_iso("2024-03-05T08:30"), Some(t));
        assert_eq!(
            DateTime::parse_iso("2024-03-05T08:30:15")
                .unwrap()
                .to_string(),
            "2024-03-05T08:30:15"
        );
        assert_eq!(DateTime::parse_iso("2024-03-05"), None);
        assert_eq!(t.add_minutes(16 * 60).date(), d.add_days(1));
        assert_eq!(
            DateTime::new(Date::from_ymd(1960, 1, 1).unwrap(), 0, 1).sas_seconds(),
            60
        );
        assert_eq!(t.add_minutes(90).hours_since(t), 1.5);
        let before = DateTime::from_seconds(-1);
        assert_eq!(before.to_string(), "1969-12-31T23:59:59");
    }

    #[test]
    fn epochs() {
        assert_eq!(Date::from_ymd(1970, 1, 1), Some(Date(0)));
        assert_eq!(Date::from_ymd(1960, 1, 1).unwrap().sas_days(), 0);
        assert_eq!(Date::from_ymd(2000, 3, 1).unwrap().days(), 11_017);
    }

    #[test]
    fn round_trips_across_centuries() {
        for days in (-200_000..200_000).step_by(37) {
            let d = Date(days);
            let (y, m, dd) = d.ymd();
            assert_eq!(Date::from_ymd(y, m, dd), Some(d));
            if (0..=9999).contains(&y) {
                assert_eq!(Date::parse_iso(&d.to_string()), Some(d));
            }
        }
    }

    #[test]
    fn invalid_dates_are_rejected() {
        assert_eq!(Date::from_ymd(2023, 2, 29), None);
        assert!(Date::from_ymd(2024, 2, 29).is_some());
        assert_eq!(Date::from_ymd(1900, 2, 29), None);
        assert_eq!(Date::parse_iso("2023-13-01"), None);
        assert_eq!(Date::parse_iso("2023-01"), None);
    }

    #[test]
    fn study_day_has_no_day_zero() {
        let rf = Date::from_ymd(2023, 3, 10).unwrap();
        assert_eq!(study_day(rf, rf), 1);
        assert_eq!(study_day(rf.add_days(1), rf), 2);
        assert_eq!(study_day(rf.add_days(-1), rf), -1);
    }

    #[test]
    fn sas_header_datetime() {
        let d = Date::from_ymd(2025, 8, 29).unwrap();
        assert_eq!(d.sas_header_datetime(), "29AUG25:00:00:00");
    }
}
