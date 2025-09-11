#![cfg(feature = "jpl-download")]
#![allow(non_snake_case)]
use camino::Utf8Path;
use rand::rngs::StdRng;
use rand::SeedableRng;

use outfit::constants::ObjectNumber;
use outfit::initial_orbit_determination::gauss_result::GaussResult;
use std::collections::{BTreeMap, HashMap};
use std::hash::BuildHasher;

use outfit::prelude::*;

// ======================= orbit classification =======================

/// A coarse taxonomy for small-body orbits in the Solar System.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum OrbitClass {
    /// Interior-Earth Objects: a < 1 AU and Q < 0.983 AU (entirely inside Earth's orbit)
    Atira,
    /// a < 1 AU and Q > 0.983 AU (crosses Earth's orbit from the inside)
    Aten,
    /// a ≥ 1 AU and q ≤ 1.017 AU (Earth-crossing from the outside)
    Apollo,
    /// 1.017 < q < 1.3 AU (near-Earth but non-crossing)
    Amor,
    /// Near-Earth Object fallback: q < 1.3 AU but not Atira/Aten/Apollo/Amor
    Neo,
    Hungaria,
    MainBelt,
    Hilda,
    Centaur,
    /// Trans-Neptunian Object
    Tno,
    /// a < 1 AU (non-Atira/Aten) — inner resonant / Vulcanoid-like
    Inner,
    Other,
}

/// Classify an orbit using simple rules on a, e (→ q, Q) and inclination i (rad).
///
/// Rules (coarse):
/// -----------------
/// * q = a(1-e), Q = a(1+e)
/// * **Atira**: a < 1.0 AU and **Q < 0.983 AU** (entirely interior to Earth's orbit)
/// * **Aten**:  a < 1.0 AU and **Q > 0.983 AU** (Earth-crossing from the inside)
/// * **Apollo**: a ≥ 1.0 AU and q ≤ 1.017 AU
/// * **Amor**: 1.017 < q < 1.3 AU
/// * **NEO**: q < 1.3 AU (and not Atira/Aten/Apollo/Amor)
/// * **Hungaria**: 1.78 ≤ a ≤ 2.0 AU and i > 16°
/// * **Hilda**: a ≈ 3.9 AU (here 3.7–4.1 AU)
/// * **Main belt**: 2.0 ≤ a ≤ 3.5 AU and e < 0.35
/// * **Centaur**: 5 < a < 30 AU
/// * **TNO**: a ≥ 30 AU
/// * **Inner**: a < 1 AU (non-Atira/Aten)
/// * **Other**: fallback
fn classify_orbit(a: f64, e: f64, i_rad: f64) -> OrbitClass {
    let q = a * (1.0 - e);
    let Q = a * (1.0 + e);
    let i_deg = i_rad.to_degrees();

    // Use aphelion Q to split Atira vs Aten when a<1 AU
    if a < 1.0 && Q < 0.983 {
        return OrbitClass::Atira;
    }
    if a < 1.0 && Q > 0.983 {
        return OrbitClass::Aten;
    }
    if a >= 1.0 && q <= 1.017 {
        return OrbitClass::Apollo;
    }
    if q > 1.017 && q < 1.3 {
        return OrbitClass::Amor;
    }
    if q < 1.3 {
        return OrbitClass::Neo;
    }

    if (1.78..=2.0).contains(&a) && i_deg > 16.0 {
        return OrbitClass::Hungaria;
    }
    if (3.7..=4.1).contains(&a) {
        return OrbitClass::Hilda;
    }
    if (2.0..=3.5).contains(&a) && e < 0.35 {
        return OrbitClass::MainBelt;
    }
    if (5.0..30.0).contains(&a) {
        return OrbitClass::Centaur;
    }
    if a >= 30.0 {
        return OrbitClass::Tno;
    }
    if a < 1.0 {
        return OrbitClass::Inner;
    }
    OrbitClass::Other
}

/// Extract a reference to `KeplerianElements` from a `GaussResult` if available.
fn kepler_of(res: &GaussResult) -> Option<&outfit::KeplerianElements> {
    match res {
        GaussResult::CorrectedOrbit(k) => k.as_keplerian(),
        GaussResult::PrelimOrbit(k) => k.as_keplerian(),
    }
}

// ======================= reporting structs =======================

#[derive(Debug, Clone)]
struct ErrorStats {
    count: usize,
    samples: Vec<String>,
    attempts_sum: usize,
    attempts_n: usize,
}

#[derive(Debug, Clone)]
struct OrbitElementStats {
    // distribution statistics for successes
    a_min: f64,
    a_max: f64,
    a_mean: f64,
    a_median: f64,
    a_p95: f64,
    e_min: f64,
    e_max: f64,
    e_mean: f64,
    e_median: f64,
    e_p95: f64,
    i_min: f64,
    i_max: f64,
    i_mean: f64,
    i_median: f64,
    i_p95: f64, // radians
    q_min: f64,
    q_max: f64,
    q_mean: f64,
    q_median: f64,
    q_p95: f64,
    Q_min: f64,
    Q_max: f64,
    Q_mean: f64,
    Q_median: f64,
    Q_p95: f64,
    // class counts
    class_counts: BTreeMap<OrbitClass, usize>,
}

#[derive(Debug, Clone)]
struct IodBatchSummary {
    total_objects: usize,
    succeeded_total: usize,
    corrected_count: usize,
    prelim_count: usize,
    failed_total: usize,

    // RMS stats on successes
    rms_min: f64,
    rms_max: f64,
    rms_mean: f64,
    rms_median: f64,
    rms_p95: f64,

    // Attempts aggregated from NoViableOrbit
    attempts_sum: usize,
    attempts_mean: f64,
    attempts_p95: f64,

    // Errors (stable kind → stats)
    error_stats: BTreeMap<&'static str, ErrorStats>,

    // Orbit elements & classes (only successes)
    elements: Option<OrbitElementStats>,

    // Top/bottom objects by RMS
    best_k: Vec<(ObjectNumber, f64)>,
    worst_k: Vec<(ObjectNumber, f64)>,
}

// ======================= summarizer =======================
#[allow(clippy::type_complexity)]
fn summarize_estimates<S>(
    results: &HashMap<ObjectNumber, Result<(GaussResult, f64), OutfitError>, S>,
    k: usize,
) -> IodBatchSummary
where
    S: BuildHasher,
{
    const MAX_SAMPLES_PER_KIND: usize = 3;

    let mut rms_values: Vec<(ObjectNumber, f64, bool)> = Vec::new();
    let mut corrected_count = 0usize;
    let mut prelim_count = 0usize;

    let mut buckets: BTreeMap<&'static str, ErrorStats> = BTreeMap::new();
    let mut attempts_all: Vec<usize> = Vec::new();

    // collect successful orbit elements
    let mut a: Vec<f64> = Vec::new();
    let mut e: Vec<f64> = Vec::new();
    let mut i: Vec<f64> = Vec::new();
    let mut q: Vec<f64> = Vec::new();
    let mut Q: Vec<f64> = Vec::new();
    let mut class_counts: BTreeMap<OrbitClass, usize> = BTreeMap::new();

    for (obj, res) in results {
        match res {
            Ok((orbit_res, rms)) => {
                let corrected = matches!(orbit_res, GaussResult::CorrectedOrbit(_));
                if corrected {
                    corrected_count += 1;
                } else {
                    prelim_count += 1;
                }
                rms_values.push((obj.clone(), *rms, corrected));

                if let Some(kep) = kepler_of(orbit_res) {
                    let a_au = kep.semi_major_axis;
                    let e_ = kep.eccentricity;
                    let i_rad = kep.inclination;
                    let q_au = a_au * (1.0 - e_);
                    let Q_au = a_au * (1.0 + e_);
                    a.push(a_au);
                    e.push(e_);
                    i.push(i_rad);
                    q.push(q_au);
                    Q.push(Q_au);

                    let cls = classify_orbit(a_au, e_, i_rad);
                    *class_counts.entry(cls).or_default() += 1;
                }
            }
            Err(e) => {
                let (stable_key, sample_msg, attempts) = flatten_error_for_bucket(e);
                let entry = buckets.entry(stable_key).or_insert_with(|| ErrorStats {
                    count: 0,
                    samples: Vec::new(),
                    attempts_sum: 0,
                    attempts_n: 0,
                });
                entry.count += 1;
                if let Some(n) = attempts {
                    entry.attempts_sum += n;
                    entry.attempts_n += 1;
                    attempts_all.push(n);
                }
                if entry.samples.len() < MAX_SAMPLES_PER_KIND {
                    entry.samples.push(sample_msg);
                }
            }
        }
    }

    let total_objects = results.len();
    let succeeded_total = corrected_count + prelim_count;
    let failed_total = total_objects.saturating_sub(succeeded_total);

    // ----- RMS stats
    let mut rms_only: Vec<f64> = rms_values.iter().map(|(_, r, _)| *r).collect();
    rms_only.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let (rms_min, rms_max, rms_mean, rms_median, rms_p95) = if rms_only.is_empty() {
        (f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN)
    } else {
        let min = *rms_only.first().unwrap();
        let max = *rms_only.last().unwrap();
        let mean = rms_only.iter().sum::<f64>() / (rms_only.len() as f64);
        let median = percentile_sorted(&rms_only, 50.0);
        let p95 = percentile_sorted(&rms_only, 95.0);
        (min, max, mean, median, p95)
    };

    // ----- attempts stats
    attempts_all.sort_unstable();
    let attempts_sum = attempts_all.iter().copied().sum::<usize>();
    let attempts_mean = if attempts_all.is_empty() {
        f64::NAN
    } else {
        attempts_sum as f64 / attempts_all.len() as f64
    };
    let attempts_p95 = if attempts_all.is_empty() {
        f64::NAN
    } else {
        percentile_sorted_usize(&attempts_all, 95.0) as f64
    };

    // ----- orbit element stats
    let elements = if a.is_empty() {
        None
    } else {
        a.sort_by(|x, y| x.partial_cmp(y).unwrap());
        e.sort_by(|x, y| x.partial_cmp(y).unwrap());
        i.sort_by(|x, y| x.partial_cmp(y).unwrap());
        q.sort_by(|x, y| x.partial_cmp(y).unwrap());
        Q.sort_by(|x, y| x.partial_cmp(y).unwrap());

        let a_stats = (
            a[0],
            *a.last().unwrap(),
            mean(&a),
            percentile_sorted(&a, 50.0),
            percentile_sorted(&a, 95.0),
        );
        let e_stats = (
            e[0],
            *e.last().unwrap(),
            mean(&e),
            percentile_sorted(&e, 50.0),
            percentile_sorted(&e, 95.0),
        );
        let i_stats = (
            i[0],
            *i.last().unwrap(),
            mean(&i),
            percentile_sorted(&i, 50.0),
            percentile_sorted(&i, 95.0),
        );
        let q_stats = (
            q[0],
            *q.last().unwrap(),
            mean(&q),
            percentile_sorted(&q, 50.0),
            percentile_sorted(&q, 95.0),
        );
        let Q_stats = (
            Q[0],
            *Q.last().unwrap(),
            mean(&Q),
            percentile_sorted(&Q, 50.0),
            percentile_sorted(&Q, 95.0),
        );

        Some(OrbitElementStats {
            a_min: a_stats.0,
            a_max: a_stats.1,
            a_mean: a_stats.2,
            a_median: a_stats.3,
            a_p95: a_stats.4,
            e_min: e_stats.0,
            e_max: e_stats.1,
            e_mean: e_stats.2,
            e_median: e_stats.3,
            e_p95: e_stats.4,
            i_min: i_stats.0,
            i_max: i_stats.1,
            i_mean: i_stats.2,
            i_median: i_stats.3,
            i_p95: i_stats.4,
            q_min: q_stats.0,
            q_max: q_stats.1,
            q_mean: q_stats.2,
            q_median: q_stats.3,
            q_p95: q_stats.4,
            Q_min: Q_stats.0,
            Q_max: Q_stats.1,
            Q_mean: Q_stats.2,
            Q_median: Q_stats.3,
            Q_p95: Q_stats.4,
            class_counts,
        })
    };

    // ----- Top/Bottom K by RMS
    rms_values.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    let best_k = rms_values
        .iter()
        .take(k)
        .map(|(o, r, _)| (o.clone(), *r))
        .collect();
    let worst_k = rms_values
        .iter()
        .rev()
        .take(k)
        .map(|(o, r, _)| (o.clone(), *r))
        .collect();

    IodBatchSummary {
        total_objects,
        succeeded_total,
        corrected_count,
        prelim_count,
        failed_total,
        rms_min,
        rms_max,
        rms_mean,
        rms_median,
        rms_p95,
        attempts_sum,
        attempts_mean,
        attempts_p95,
        error_stats: buckets,
        elements,
        best_k,
        worst_k,
    }
}

// ======================= printer =======================

fn print_summary(sum: &IodBatchSummary) {
    println!("\n=== IOD Batch Summary ===");
    println!("Objects total ........ : {}", sum.total_objects);
    println!(
        "  Success total ...... : {}  (Corrected: {}, Prelim: {})",
        sum.succeeded_total, sum.corrected_count, sum.prelim_count
    );
    println!("  Failures total ..... : {}", sum.failed_total);

    println!("\n-- RMS on successes --");
    println!(
        "  min / median / mean / p95 / max : {:.6} / {:.6} / {:.6} / {:.6} / {:.6}",
        sum.rms_min, sum.rms_median, sum.rms_mean, sum.rms_p95, sum.rms_max
    );

    if let Some(el) = &sum.elements {
        println!("\n-- Orbit elements (successes only) --");
        println!(
            "  a [AU]  : min={:.6}  median={:.6}  mean={:.6}  p95={:.6}  max={:.6}",
            el.a_min, el.a_median, el.a_mean, el.a_p95, el.a_max
        );
        println!(
            "  e [-]   : min={:.6}  median={:.6}  mean={:.6}  p95={:.6}  max={:.6}",
            el.e_min, el.e_median, el.e_mean, el.e_p95, el.e_max
        );
        println!(
            "  i [deg] : min={:.6}  median={:.6}  mean={:.6}  p95={:.6}  max={:.6}",
            el.i_min.to_degrees(),
            el.i_median.to_degrees(),
            el.i_mean.to_degrees(),
            el.i_p95.to_degrees(),
            el.i_max.to_degrees()
        );
        println!(
            "  q [AU]  : min={:.6}  median={:.6}  mean={:.6}  p95={:.6}  max={:.6}",
            el.q_min, el.q_median, el.q_mean, el.q_p95, el.q_max
        );
        println!(
            "  Q [AU]  : min={:.6}  median={:.6}  mean={:.6}  p95={:.6}  max={:.6}",
            el.Q_min, el.Q_median, el.Q_mean, el.Q_p95, el.Q_max
        );

        println!("\n-- Orbit classes (counts) --");
        if el.class_counts.is_empty() {
            println!("  (none)");
        } else {
            for (cls, cnt) in &el.class_counts {
                println!("  {:<10} : {}", format!("{cls:?}"), cnt);
            }
        }
    } else {
        println!("\n-- Orbit elements --\n  (no successful orbits)");
    }

    println!("\n-- Failure breakdown (by error kind) --");
    if sum.error_stats.is_empty() {
        println!("  (none)");
    } else {
        for (kind, stats) in &sum.error_stats {
            let attempts_avg = if stats.attempts_n > 0 {
                stats.attempts_sum as f64 / stats.attempts_n as f64
            } else {
                f64::NAN
            };
            println!(
                "  {:<28} : {:>4}   attempts(sum/avg on carrying) = {}/{}",
                kind,
                stats.count,
                stats.attempts_sum,
                fmt_opt(attempts_avg)
            );
            for (i, sample) in stats.samples.iter().take(3).enumerate() {
                println!("      · example[{i}] {sample}");
            }
        }
    }

    println!("\n-- Attempts (global from NoViableOrbit) --");
    if sum.attempts_sum == 0 && sum.failed_total > 0 {
        println!("  (no attempts recorded in errors)");
    } else if sum.failed_total == 0 {
        println!("  (no failures)");
    } else {
        println!(
            "  sum / mean / p95 : {} / {:.2} / {:.0}",
            sum.attempts_sum, sum.attempts_mean, sum.attempts_p95
        );
    }

    println!("\n-- Best by RMS --");
    if sum.best_k.is_empty() {
        println!("  (no successes)");
    } else {
        for (obj, rms) in &sum.best_k {
            println!("  {:>8}  rms={:.6}", display_obj(obj), rms);
        }
    }

    println!("\n-- Worst by RMS --");
    if sum.worst_k.is_empty() {
        println!("  (no successes)");
    } else {
        for (obj, rms) in &sum.worst_k {
            println!("  {:>8}  rms={:.6}", display_obj(obj), rms);
        }
    }
    println!();
}

// ======================= small helpers =======================

fn flatten_error_for_bucket(e: &OutfitError) -> (&'static str, String, Option<usize>) {
    use OutfitError::*;
    match e {
        NoViableOrbit { cause, attempts } => {
            let (k, _, _) = flatten_error_for_bucket(cause);
            (k, format!("{e}"), Some(*attempts))
        }
        GaussNoRootsFound => ("GaussNoRootsFound", format!("{e}"), None),
        PolynomialRootFindingFailed => ("PolynomialRootFindingFailed", format!("{e}"), None),
        SpuriousRootDetected => ("SpuriousRootDetected", format!("{e}"), None),
        SingularDirectionMatrix => ("SingularDirectionMatrix", format!("{e}"), None),
        RmsComputationFailed(_) => ("RmsComputationFailed", format!("{e}"), None),
        GaussPrelimOrbitFailed(_) => ("GaussPrelimOrbitFailed", format!("{e}"), None),
        InvalidOrbit(_) => ("InvalidOrbit", format!("{e}"), None),
        InvalidIODParameter(_) => ("InvalidIODParameter", format!("{e}"), None),
        InvalidConversion(_) => ("InvalidConversion", format!("{e}"), None),
        RootFindingError(_) => ("RootFindingError", format!("{e}"), None),
        NoiseInjectionError(_) => ("NoiseInjectionError", format!("{e}"), None),
        InvalidRefSystem(_) => ("InvalidRefSystem", format!("{e}"), None),
        ObservationNotFound(_) => ("ObservationNotFound", format!("{e}"), None),
        Parquet(_) => ("Parquet", format!("{e}"), None),
        IoError(_) => ("IoError", format!("{e}"), None),
        UreqHttpError(_) => ("UreqHttpError", format!("{e}"), None),
        #[cfg(feature = "jpl-download")]
        ReqwestError(_) => ("ReqwestError", format!("{e}"), None),
        InvalidJPLStringFormat(_) => ("InvalidJPLStringFormat", format!("{e}"), None),
        InvalidJPLEphemFileSource(_) => ("InvalidJPLEphemFileSource", format!("{e}"), None),
        InvalidJPLEphemFileVersion(_) => ("InvalidJPLEphemFileVersion", format!("{e}"), None),
        JPLFileNotFound(_) => ("JPLFileNotFound", format!("{e}"), None),
        Utf8PathError(_) => ("Utf8PathError", format!("{e}"), None),
        UnableToCreateBaseDir(_) => ("UnableToCreateBaseDir", format!("{e}"), None),
        NomParsingError(_) => ("NomParsingError", format!("{e}"), None),
        Parsing80ColumnFileError(_) => ("Parsing80ColumnFileError", format!("{e}"), None),
        InvalidErrorModel(_) => ("InvalidErrorModel", format!("{e}"), None),
        InvalidErrorModelFilePath(_) => ("InvalidErrorModelFilePath", format!("{e}"), None),
        InvalidSpkDataType(_) => ("InvalidSpkDataType", format!("{e}"), None),
        InvalidFloatValue(_) => ("InvalidFloatValue", format!("{e}"), None),
        _ => ("Other", format!("{e}"), None),
    }
}

fn percentile_sorted(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    let q = p.clamp(0.0, 100.0);
    let idx = ((q / 100.0) * ((sorted.len() - 1) as f64)).round() as usize;
    sorted[idx]
}

fn percentile_sorted_usize(sorted: &[usize], p: f64) -> usize {
    if sorted.is_empty() {
        return 0;
    }
    let q = p.clamp(0.0, 100.0);
    let idx = ((q / 100.0) * ((sorted.len() - 1) as f64)).round() as usize;
    sorted[idx]
}

fn mean(v: &[f64]) -> f64 {
    if v.is_empty() {
        f64::NAN
    } else {
        v.iter().sum::<f64>() / v.len() as f64
    }
}

fn display_obj(id: &ObjectNumber) -> String {
    match id {
        ObjectNumber::Int(n) => format!("{n}"),
        ObjectNumber::String(s) => s.clone(),
    }
}

fn fmt_opt(x: f64) -> String {
    if x.is_finite() {
        format!("{x:.2}")
    } else {
        "-".to_string()
    }
}

fn main() -> Result<(), OutfitError> {
    let mut env_state = Outfit::new("horizon:DE440", ErrorModel::FCCT14).unwrap();

    let test_data = "tests/data/test_from_fink.parquet";

    println!("Loading observations from {test_data}");
    let path_file = Utf8Path::new(test_data);

    let ztf_observer = env_state.get_observer_from_mpc_code(&"I41".into());

    let mut traj_set =
        TrajectorySet::new_from_parquet(&mut env_state, path_file, ztf_observer, 0.5, 0.5, None)?;

    println!(
        "Loading done: {} trajectories / {} observations",
        traj_set.len(),
        traj_set.total_observations()
    );

    println!("Trajectory Set statistics:");
    println!(
        "{:#}",
        traj_set.obs_count_stats().expect("TrajSet is empty")
    );
    println!("--------------------------\n");

    let mut rng = StdRng::from_os_rng();

    let params = IODParams::builder()
        .n_noise_realizations(10)
        .noise_scale(1.1)
        .max_obs_for_triplets(12)
        .max_triplets(30)
        .build()?;

    println!("Estimating orbits with params: {params:#}");
    println!("Orbit estimation in progress... (this may take a while)");

    let orbits = traj_set.estimate_all_orbits(&env_state, &mut rng, &params);

    println!("Orbit estimation done.");

    let summary = summarize_estimates(&orbits, 5);
    print_summary(&summary);

    Ok(())
}
