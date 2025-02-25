use std::str::FromStr;

enum StepUnit {
    Days,
    Hours,
    Minutes,
    Years,
    Months,
}

struct Step {
    value: u32,
    unit: StepUnit,
}

impl Step {
    pub fn new(value: u32, unit: StepUnit) -> Self {
        Step {
            value: value,
            unit: unit,
        }
    }

    pub fn to_string(&self) -> String {
        match self.unit {
            StepUnit::Days => return format!("{}d", self.value),
            StepUnit::Hours => return format!("{}h", self.value),
            StepUnit::Minutes => return format!("{}m", self.value),
            StepUnit::Months => return format!("{}mo", self.value),
            StepUnit::Years => return format!("{}y", self.value),
        }
    }
}

struct IntervalTime {
    start_time: f64,
    stop_time: f64,
    step: Step,
}

impl IntervalTime {
    pub fn new(start_time: f64, stop_time: f64, step: Step) -> IntervalTime {
        IntervalTime {
            start_time: start_time,
            stop_time: stop_time,
            step: step,
        }
    }
}

fn jpl_params(
    observatory: &str,
    target_body: &str,
    interval_time: &IntervalTime,
) -> [(String, String); 11] {
    let start_time = format!("JD{}", interval_time.start_time);
    let stop_time = format!("JD{}", interval_time.stop_time);
    let step = interval_time.step.to_string();
    [
        ("format".into(), "text".into()),
        ("CENTER".into(), observatory.into()),
        ("COMMAND".into(), target_body.into()),
        ("OBJ_DATA".into(), "NO".into()),
        ("MAKE_EPHEM".into(), "YES".into()),
        ("EPHEM_TYPE".into(), "VECTORS".into()),
        ("START_TIME".into(), start_time),
        ("STOP_TIME".into(), stop_time),
        ("STEP_SIZE".into(), step),
        ("CSV_FORMAT".into(), "YES".into()),
        ("TLIST_TYPE".into(), "JD".into()),
    ]
}

#[cfg(test)]
mod jpl_ephem_tests {
    use super::*;

    #[test]
    fn test_step() {
        let step_day = Step::new(1, StepUnit::Days);
        assert_eq!(step_day.to_string(), "1d");
        let step_hours = Step::new(50, StepUnit::Hours);
        assert_eq!(step_hours.to_string(), "50h");
        let step_minute = Step::new(30, StepUnit::Minutes);
        assert_eq!(step_minute.to_string(), "30m");
        let step_year = Step::new(5, StepUnit::Years);
        assert_eq!(step_year.to_string(), "5y");
        let step_month = Step::new(6, StepUnit::Months);
        assert_eq!(step_month.to_string(), "6mo")
    }

    #[test]
    fn test_interval_time() {
        let interval1 = IntervalTime::new(0.0, 10.0, Step::new(1, StepUnit::Days));
        assert_eq!(interval1.step.to_string(), "1d");
        assert_eq!(interval1.start_time, 0.0);
        assert_eq!(interval1.stop_time, 10.0);

    }
}