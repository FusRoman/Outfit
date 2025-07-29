use hifitime::{Epoch, TimeScale};
use std::str::FromStr;

/// Transformation from date in the format YYYY-MM-ddTHH:mm:ss to modified julian date (MJD)
///
/// Argument
/// --------
/// * `date`: a vector of date in the format YYYY-MM-ddTHH:mm:ss
///
/// Return
/// ------
/// * a vector of float representing the input date in modified julian date (MJD)
pub fn date_to_mjd(date: &Vec<&str>) -> Vec<f64> {
    date.iter()
        .map(|x| Epoch::from_str(x).unwrap().to_mjd_utc_days())
        .collect::<Vec<f64>>()
}

/// Transformation from modified julian date (MJD) in julian date (JD)
///
/// Argument
/// --------
/// * `mjd`: a vector of MJD
///
/// Return
/// ------
/// * a vector of jd
pub fn mjd_to_jd(mjd: &[f64]) -> Vec<f64> {
    mjd.iter()
        .map(|x| Epoch::from_mjd_utc(*x).to_jde_utc_days())
        .collect()
}

/// Transformation from julian date (JD) in modified julian date (MJD)
///
/// Argument
/// --------
/// * `jd`: a vector of JD
///
/// Return
/// ------
/// * a vector of MJD
pub fn jd_to_mjd(jd: &[f64]) -> Vec<f64> {
    jd.iter()
        .map(|x| Epoch::from_jde_utc(*x).to_mjd_utc_days())
        .collect()
}

/// Transformation from date in the format YYYY MM DD.FFFFF UTC frame to modified julian date (MJD) TT frame
///
/// Argument
/// --------
/// * `date_str`: a string representing the date in the format YYYY MM DD.FFFFF in the UTC frame
///
/// Return
/// ------
/// * a float representing the input date in modified julian date (MJD) in the TT frame
pub fn frac_date_to_mjd(date_str: &str) -> Result<f64, String> {
    let parts: Vec<&str> = date_str.split_whitespace().collect();
    if parts.len() != 3 {
        return Err("Invalid format, expected: YYYY MM DD.FFFFF".to_string());
    }

    // Extract values
    let year = i32::from_str(parts[0]).map_err(|_| "invalid year")?;
    let month = u8::from_str(parts[1]).map_err(|_| "invalid month")?;
    let day_fraction = f64::from_str(parts[2]).map_err(|_| "invalid frac day")?;

    // Separation of day and fraction day
    let day = day_fraction.trunc() as u8;
    let fraction = day_fraction - day as f64;

    let hour = (fraction * 24.0).trunc() as u8;
    let minute = ((fraction * 24.0 - hour as f64) * 60.0).trunc() as u8;
    let second = (((fraction * 24.0 - hour as f64) * 60.0 - minute as f64) * 60.0) as u8;
    let nano = ((((fraction * 24.0 - hour as f64) * 60.0 - minute as f64) * 60.0 - second as f64)
        * 1e9) as u32;

    // Creation of epoch
    let epoch = Epoch::from_gregorian(year, month, day, hour, minute, second, nano, TimeScale::UTC);

    // Convert to MJD
    Ok(epoch.to_mjd_tt_days())
}

#[cfg(test)]
mod time_test {
    use super::*;

    #[test]
    fn test_date_to_mjd() {
        let date = vec!["2021-01-01T00:00:00", "2021-01-02T00:00:00"];
        let mjd = date_to_mjd(&date);
        assert_eq!(mjd, vec![59215.0, 59216.0]);
    }

    #[test]
    fn test_mjd_to_jd() {
        let mjd = vec![59215.0, 59216.0];
        let jd = mjd_to_jd(&mjd);
        assert_eq!(jd, vec![2459215.5, 2459216.5]);
    }

    #[test]
    fn test_jd_to_mjd() {
        let jd = vec![2459215.5, 2459216.5];
        let mjd = jd_to_mjd(&jd);
        assert_eq!(mjd, vec![59215.0, 59216.0]);
    }

    #[test]
    fn test_frac_date_to_mjd() {
        let mjd = frac_date_to_mjd("2021 1 1.0").unwrap();
        assert_eq!(mjd, 59215.00080074074);

        let mjd = frac_date_to_mjd("2021 1 1.5").unwrap();
        assert_eq!(mjd, 59215.50080074074);

        let mjd = frac_date_to_mjd("2021 1 1.75").unwrap();
        assert_eq!(mjd, 59215.75080074074);

        let mjd = frac_date_to_mjd("2021 1 1.875").unwrap();
        assert_eq!(mjd, 59215.87580074074);

        let mjd = frac_date_to_mjd("2021 1 1.999").unwrap();
        assert_eq!(mjd, 59215.99980074074);

        let mjd = frac_date_to_mjd("2021 1 1.9999").unwrap();
        assert_eq!(mjd, 59216.00070074073);

        let mjd = frac_date_to_mjd("1976 09 20.93878").unwrap();
        assert_eq!(mjd, 43041.93932611111);
    }
}
