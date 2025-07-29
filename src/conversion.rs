use crate::constants::{ArcSec, Degree};

/// Estimate the accuracy of a numeric string based on its decimal precision.
///
/// Arguments
/// ---------------
/// * `field`: a string slice containing the numeric value (e.g., `"56.78"`), typically the last component of an angle
/// * `factor`: a scale factor to apply to the accuracy (e.g., `1.0 / 60.0` for arcminutes, `1.0 / 3600.0` for arcseconds)
///
/// Return
/// ----------
/// * `Option<f64>`: the estimated accuracy scaled by `factor`, or `None` if the input is malformed
fn compute_accuracy(field: &str, factor: f64) -> Option<f64> {
    if let Some(dot_pos) = field.find('.') {
        let length = field.trim().len();
        let digits_after_dot = length - dot_pos - 1;
        let acc = 10f64.powi(-(digits_after_dot as i32)) * factor;
        Some(acc)
    } else {
        Some(1.0 * factor)
    }
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
    let acc_arcsec = compute_accuracy(s_raw, 1.0 / 3600.0)?;
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
    let acc_arcsec = compute_accuracy(s_raw, 1. / 3600.)?;
    Some((dec_deg, acc_arcsec))
}

#[cfg(test)]
mod observations_test {
    use super::*;

    #[test]
    fn test_ra_to_deg() {
        assert_eq!(
            parse_ra_to_deg("22 52 23.37"),
            Some((343.097375, 2.777777777777778e-6))
        );
        assert_eq!(
            parse_ra_to_deg("23 58 57.68"),
            Some((359.7403333333333, 2.777777777777778e-6))
        );
        assert_eq!(
            parse_ra_to_deg("04 41 04.77"),
            Some((70.269875, 2.777777777777778e-6))
        );
        assert_eq!(parse_ra_to_deg("1 2 3.4.5"), None);
        assert_eq!(parse_ra_to_deg("1 2"), None);
        assert_eq!(
            parse_ra_to_deg("06 50 13.370"),
            Some((102.55570833333333, 2.7777777777777776e-7))
        );
    }

    #[test]
    fn test_dec_to_deg() {
        assert_eq!(
            parse_dec_to_deg("-00 30 14.2"),
            Some((-0.5039444444444444, 2.777777777777778e-5))
        );
        assert_eq!(
            parse_dec_to_deg("+13 55 42.7"),
            Some((13.928527777777777, 2.777777777777778e-5))
        );
        assert_eq!(
            parse_dec_to_deg("89 15 50.2"),
            Some((89.26394444444445, 2.777777777777778e-5))
        );
        assert_eq!(parse_dec_to_deg("89 15 50.2.3"), None);
        assert_eq!(parse_dec_to_deg("89 15"), None);

        assert_eq!(
            parse_dec_to_deg("-14 47 05.4"),
            Some((-14.784833333333333, 2.777777777777778e-5))
        );
    }

    #[test]
    fn test_estimate_accuracy() {
        assert_eq!(
            compute_accuracy("23.3", 1. / 3600.),
            Some(2.777777777777778e-5)
        );
        assert_eq!(
            compute_accuracy("23", 1. / 3600.),
            Some(0.0002777777777777778)
        );
        assert_eq!(
            compute_accuracy("23.370", 1. / 3600.),
            Some(2.7777777777777776e-7)
        );
        assert_eq!(
            compute_accuracy("23.37", 1. / 3600.),
            Some(2.777777777777778e-6)
        );
    }
}
