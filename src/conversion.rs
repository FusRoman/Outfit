use crate::constants::{ArcSec, Degree};

/// Estimate the accuracy from the fractional seconds field of a time string.
///
/// Arguments
/// ---------
/// * `field`: A string representing the fractional seconds (e.g., "23.37").
///
/// Returns
/// -------
/// * `Option<f64>`: The estimated accuracy in arcseconds. Returns `None` if the input string
///   does not contain a valid fractional part.
fn estimate_accuracy(field: &str) -> Option<ArcSec> {
    let dot_pos = field.find('.')?;
    let total_len = field.len();
    let decimal_digits = total_len - dot_pos - 1;
    let factor = 10f64.powi(-(decimal_digits as i32));
    Some(factor)
}

/// Parse a right ascension string to degrees
///
/// Arguments
/// ---------
/// * `ra`: a string representing the right ascension in the format `HH MM SS.SS`
///
/// Returns
/// -------
/// * `Option<(Degree, ArcSec)>`: A tuple where the first element is the right ascension in degrees
///   and the second element is the accuracy in arcseconds. Returns `None` if the input format is invalid.
pub(crate) fn parse_ra_to_deg(ra: &str) -> Option<(Degree, ArcSec)> {
    let parts: Vec<&str> = ra.split_whitespace().collect();
    if parts.len() != 3 {
        return None;
    }

    let h: f64 = parts[0].parse().ok()?;
    let m: f64 = parts[1].parse().ok()?;
    let s_raw = parts[2];
    let s: f64 = s_raw.parse().ok()?;

    let ra_deg = (h + m / 60.0 + s / 3600.0) * 15.0;
    let acc_arcsec = estimate_accuracy(s_raw)? * 15.0; // convert seconds to arcsec in RA
    Some((ra_deg, acc_arcsec))
}

/// Parse a declination string to degrees
///
/// Arguments
/// ---------
/// * `dec`: a string representing the declination in the format `±DD MM SS.SS`
///
/// Returns
/// -------
/// * `Option<(Degree, ArcSec)>`: A tuple where the first element is the declination in degrees
///   and the second element is the accuracy in arcseconds. Returns `None` if the input format is invalid.
pub fn parse_dec_to_deg(dec: &str) -> Option<(Degree, ArcSec)> {
    let parts: Vec<&str> = dec.split_whitespace().collect();
    if parts.len() != 3 {
        return None;
    }

    let sign = if parts[0].starts_with('-') { -1.0 } else { 1.0 };
    let d: f64 = parts[0].trim_start_matches(&['-', '+'][..]).parse().ok()?;
    let m: f64 = parts[1].parse().ok()?;
    let s_raw = parts[2];
    let s: f64 = s_raw.parse().ok()?;

    let dec_deg = sign * (d + m / 60.0 + s / 3600.0);
    let acc_arcsec = estimate_accuracy(s_raw)?;
    Some((dec_deg, acc_arcsec))
}

#[cfg(test)]
mod observations_test {
    use super::*;

    #[test]
    fn test_ra_to_deg() {
        assert_eq!(parse_ra_to_deg("22 52 23.37"), Some((343.097375, 0.15)));
        assert_eq!(
            parse_ra_to_deg("23 58 57.68"),
            Some((359.7403333333333, 0.15))
        );
        assert_eq!(parse_ra_to_deg("04 41 04.77"), Some((70.269875, 0.15)));
        assert_eq!(parse_ra_to_deg("1 2 3.4.5"), None);
        assert_eq!(parse_ra_to_deg("1 2"), None);
    }

    #[test]
    fn test_dec_to_deg() {
        assert_eq!(
            parse_dec_to_deg("-00 30 14.2"),
            Some((-0.5039444444444444, 0.1))
        );
        assert_eq!(
            parse_dec_to_deg("+13 55 42.7"),
            Some((13.928527777777777, 0.1))
        );
        assert_eq!(
            parse_dec_to_deg("89 15 50.2"),
            Some((89.26394444444445, 0.1))
        );
        assert_eq!(parse_dec_to_deg("89 15 50.2.3"), None);
        assert_eq!(parse_dec_to_deg("89 15"), None);
    }
}
