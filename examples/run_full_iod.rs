use std::collections::HashMap;

use hifitime::ut1::Ut1Provider;
use outfit::{
    jpl_ephem::download_jpl_file::EphemFileSource, FitIOD, IODParams, JPLEphem, OutfitError,
};
use photom::{
    io::polars::{ContiguousChoice, FromPolarsArgs},
    observation_dataset::ObsDataset,
    observer::error_model::ObsErrorModel,
};
use polars::lazy::{
    dsl::{col, lit},
    frame::LazyFrame,
};
use rand::{rngs::StdRng, SeedableRng};

fn outfit_error_label(err: &OutfitError) -> &'static str {
    match err {
        OutfitError::InvalidJPLStringFormat(_) => "InvalidJPLStringFormat",
        OutfitError::InvalidJPLEphemFileSource(_) => "InvalidJPLEphemFileSource",
        OutfitError::InvalidJPLEphemFileVersion(_) => "InvalidJPLEphemFileVersion",
        OutfitError::InvalidUrl(_) => "InvalidUrl",
        OutfitError::UreqHttpError(_) => "UreqHttpError",
        OutfitError::IoError(_) => "IoError",
        OutfitError::ReqwestError(_) => "ReqwestError",
        OutfitError::UnableToCreateBaseDir(_) => "UnableToCreateBaseDir",
        OutfitError::Utf8PathError(_) => "Utf8PathError",
        OutfitError::JPLFileNotFound(_) => "JPLFileNotFound",
        OutfitError::RootFindingError(_) => "RootFindingError",
        OutfitError::ObservationNotFound(_) => "ObservationNotFound",
        OutfitError::InvalidErrorModel(_) => "InvalidErrorModel",
        OutfitError::InvalidErrorModelFilePath(_) => "InvalidErrorModelFilePath",
        OutfitError::NomParsingError(_) => "NomParsingError",
        OutfitError::NoiseInjectionError(_) => "NoiseInjectionError",
        OutfitError::SingularDirectionMatrix => "SingularDirectionMatrix",
        OutfitError::PolynomialRootFindingFailed => "PolynomialRootFindingFailed",
        OutfitError::SpuriousRootDetected => "SpuriousRootDetected",
        OutfitError::GaussNoRootsFound => "GaussNoRootsFound",
        OutfitError::InvalidSpkDataType(_) => "InvalidSpkDataType",
        OutfitError::InvalidIODParameter(_) => "InvalidIODParameter",
        OutfitError::InvalidRefSystem(_) => "InvalidRefSystem",
        OutfitError::VelocityCorrectionError(_) => "VelocityCorrectionError",
        OutfitError::InvalidOrbit(_) => "InvalidOrbit",
        OutfitError::InvalidConversion(_) => "InvalidConversion",
        OutfitError::InvalidFloatValue(_) => "InvalidFloatValue",
        OutfitError::RmsComputationFailed(_) => "RmsComputationFailed",
        OutfitError::GaussPrelimOrbitFailed(_) => "GaussPrelimOrbitFailed",
        OutfitError::NoViableOrbit { .. } => "NoViableOrbit",
        OutfitError::NoFeasibleTriplets { .. } => "NoFeasibleTriplets",
        OutfitError::NonFiniteScore(_) => "NonFiniteScore",
        OutfitError::ObserverIdIsNone(_) => "ObserverIdIsNone",
        OutfitError::ObsDatasetError(_) => "ObsDatasetError",
        OutfitError::ObsDatasetErrorRef(_) => "ObsDatasetErrorRef",
        OutfitError::TrajectoryIdNotFound(_) => "TrajectoryIdNotFound",
        OutfitError::NoTrajectoryIndex => "NoTrajectoryIndex",
    }
}

fn main() -> Result<(), OutfitError> {
    let path_data = "tests/data/test_data_traj_str.parquet";

    let lf = LazyFrame::scan_parquet(path_data.into(), Default::default())
        .expect("scan_parquet must succeed")
        .filter(col("traj_id").is_not_null())
        .filter(
            col("traj_id")
                .count()
                .over([col("traj_id")])
                .gt_eq(lit(3u32)),
        );

    let polars_args = FromPolarsArgs {
        error_model: Some(ObsErrorModel::FCCT14),
        do_rechunk: Some(false),
        contiguous_choice: Some(ContiguousChoice::ContiguousTraj),
    };

    let obs_dataset = ObsDataset::from_lazy(lf, polars_args).unwrap();
    println!("\n ==== Observation Dataset ==== \n");
    println!("{obs_dataset}");
    println!("\n ---- \n");

    let max_traj_size = obs_dataset
        .iter_traj_id()
        .unwrap()
        .map(|traj_id| obs_dataset.len_trajectory(traj_id).unwrap())
        .max()
        .unwrap();

    let default = IODParams::builder()
        .n_noise_realizations(10)
        .noise_scale(1.1)
        .max_obs_for_triplets(max_traj_size)
        .max_triplets(30)
        .build()
        .unwrap();

    println!("\n ==== IOD Parameters ==== \n");
    println!("{default:#?}");
    println!("\n ---- \n");

    let ut1_provider = Ut1Provider::download_from_jpl("latest_eop2.long")
        .expect("Download of the JPL short time scale UT1 data failed");

    let jpl_file: EphemFileSource = "horizon:DE440"
        .try_into()
        .expect("Failed to parse JPL ephemeris source");
    let jpl_ephem = JPLEphem::new(&jpl_file).expect("Failed to load JPL ephemeris from Horizon");

    // fully sequential version of the full IOD fitting, for comparison with the parallel version in run_full_iod_parallel.rs
    let full_orbit = obs_dataset
        .fit_full_iod(
            &jpl_ephem,
            &ut1_provider,
            &default,
            ObsErrorModel::FCCT14,
            &mut StdRng::seed_from_u64(42),
        )
        .unwrap();

    let number_total_fit = full_orbit.len();
    let number_success_fit = full_orbit
        .iter()
        .filter(|(_, orbit_res)| orbit_res.is_ok())
        .count();
    let number_failed_fit = number_total_fit - number_success_fit;

    println!("\n ==== Full IOD Results ==== \n");
    println!("Total number of fits: {number_total_fit}");
    println!("Number of successful fits: {number_success_fit}");
    println!("Number of failed fits: {number_failed_fit}");

    let success_rate = (number_success_fit as f64) / (number_total_fit as f64) * 100.0;
    println!("Success rate: {success_rate:.2}%");

    println!("\n === Successful Fits Details === \n");
    println!("Orbit quality :");

    let all_rms = full_orbit
        .iter()
        .filter_map(|(_, res)| res.as_ref().ok().map(|orbit| orbit.orbit_quality()))
        .collect::<Vec<_>>();

    let mean_rms = all_rms.iter().sum::<f64>() / (all_rms.len() as f64);
    let median_rms = {
        let mut sorted = all_rms.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mid = sorted.len() / 2;
        if sorted.len() % 2 == 0 {
            (sorted[mid - 1] + sorted[mid]) / 2.0
        } else {
            sorted[mid]
        }
    };
    let min_rms = all_rms.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_rms = all_rms.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    println!("  Mean RMS: {mean_rms:.4}");
    println!("  Median RMS: {median_rms:.4}");
    println!("  Min RMS: {min_rms:.4}");
    println!("  Max RMS: {max_rms:.4}");

    println!("RMS distribution:");
    let mut rms_bins: HashMap<String, usize> = HashMap::new();
    for rms in all_rms {
        let bin = if rms < 1.0 {
            "< 1.0"
        } else if rms < 5.0 {
            "1.0 - 5.0"
        } else if rms < 10.0 {
            "5.0 - 10.0"
        } else {
            ">= 10.0"
        };
        *rms_bins.entry(bin.to_string()).or_insert(0) += 1;
    }
    let mut sorted_bins: Vec<_> = rms_bins.iter().collect();
    sorted_bins.sort_by_key(|(_, count)| std::cmp::Reverse(*count));
    for (bin, count) in sorted_bins {
        println!("  {bin:<10} : {count}");
    }
    println!("\n ---- \n");

    println!("\n === Failed Fits Details === \n");
    let mut failed_fits_stats: HashMap<&str, usize> = HashMap::new();

    for (_, orbit_res) in full_orbit.iter().filter(|(_, res)| res.is_err()) {
        let error = orbit_res.as_ref().err().unwrap();

        *failed_fits_stats
            .entry(outfit_error_label(error))
            .or_insert(0) += 1;
    }

    let mut sorted: Vec<_> = failed_fits_stats.iter().collect();
    sorted.sort_by_key(|(_, count)| std::cmp::Reverse(*count));

    for (label, count) in sorted {
        println!("  {label:<35} : {count}");
    }

    Ok(())
}
