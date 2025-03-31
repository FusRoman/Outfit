use crate::constants::Degree32;

/// Parse a right ascension string to degrees
///
/// Arguments
/// ---------
/// * `ra`: a string representing the right ascension in the format HH MM SS.SS
///
/// Return
/// ------
/// * a float representing the right ascension in degrees
pub (crate) fn parse_ra_to_deg(ra: &str) -> Option<Degree32> {
    let parts: Vec<&str> = ra.split_whitespace().collect();
    if parts.len() != 3 {
        return None;
    }
    let h: f32 = parts[0].parse().ok()?;
    let m: f32 = parts[1].parse().ok()?;
    let s: f32 = parts[2].parse().ok()?;
    Some((h + m / 60.0 + s / 3600.0) * 15.0)
}

/// Parse a declination string to degrees
///
/// Arguments
/// ---------
/// * `dec`: a string representing the declination in the format HH MM SS.SS
///
/// Return
/// ------
/// * a float representing the declination in degrees
pub (crate) fn parse_dec_to_deg(dec: &str) -> Option<Degree32> {
    let parts: Vec<&str> = dec.split_whitespace().collect();
    if parts.len() != 3 {
        return None;
    }
    let sign = if parts[0].starts_with('-') { -1.0 } else { 1.0 };
    let d: f32 = parts[0]
        .trim_start_matches(['-', '+'].as_ref())
        .parse()
        .ok()?;
    let m: f32 = parts[1].parse().ok()?;
    let s: f32 = parts[2].parse().ok()?;
    Some(sign * (d + m / 60.0 + s / 3600.0))
}


#[cfg(test)]
mod observations_test {
    use super::*;

    #[test]
    fn test_ra_to_deg() {
        assert_eq!(parse_ra_to_deg("23 58 57.68"), Some(359.7403333333333));
        assert_eq!(parse_ra_to_deg("04 41 04.77"), Some(70.269875));
        assert_eq!(parse_ra_to_deg("1 2 3.4.5"), None);
        assert_eq!(parse_ra_to_deg("1 2"), None);
    }

    #[test]
    fn test_dec_to_deg() {
        assert_eq!(parse_dec_to_deg("-00 30 14.2"), Some(-0.5039444444444444));
        assert_eq!(parse_dec_to_deg("+13 55 42.7"), Some(13.928527777777777));
        assert_eq!(parse_dec_to_deg("89 15 50.2"), Some(89.26394444444445));
        assert_eq!(parse_dec_to_deg("89 15 50.2.3"), None);
        assert_eq!(parse_dec_to_deg("89 15"), None);
    }
}
