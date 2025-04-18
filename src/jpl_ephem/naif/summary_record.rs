use std::fmt;

use hifitime::Epoch;
use nom::{
    number::complete::{le_f64, le_i32},
    IResult,
};

use crate::jpl_ephem::naif::naif_ids::{naif_type::SpkDataType, NaifIds};

#[derive(Debug, PartialEq)]
pub struct Summary {
    pub start_epoch: f64,
    pub end_epoch: f64,
    pub target: i32,
    pub center: i32,
    pub frame_id: i32,
    pub data_type: i32,
    pub initial_addr: i32,
    pub final_addr: i32,
}

impl Summary {
    pub fn parse(input: &[u8]) -> IResult<&[u8], Self> {
        let (input, start_epoch) = le_f64(input)?;
        let (input, end_epoch) = le_f64(input)?;

        let (input, target) = le_i32(input)?;
        let (input, center) = le_i32(input)?;
        let (input, frame_id) = le_i32(input)?;
        let (input, data_type) = le_i32(input)?;
        let (input, initial_addr) = le_i32(input)?;
        let (input, final_addr) = le_i32(input)?;
        Ok((
            input,
            Summary {
                start_epoch,
                end_epoch,
                target: target,
                center: center,
                frame_id: frame_id,
                data_type: data_type,
                initial_addr: initial_addr,
                final_addr: final_addr,
            },
        ))
    }
}

impl fmt::Display for Summary {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let start = Epoch::from_et_seconds(self.start_epoch);
        let end = Epoch::from_et_seconds(self.end_epoch);

        let naif_target = match NaifIds::from_id(self.target) {
            Ok(target) => target,
            Err(_) => return Err(fmt::Error),
        };
        let naif_center = match NaifIds::from_id(self.center) {
            Ok(center) => center,
            Err(_) => return Err(fmt::Error),
        };

        let naif_type = match SpkDataType::from_i32(self.data_type) {
            Ok(naif_type) => naif_type,
            Err(_) => return Err(fmt::Error),
        };

        let fields = vec![
            ("start_epoch", format!("{}", start)),
            ("end_epoch", format!("{}", end)),
            ("target", format!("{}", naif_target)),
            ("center", format!("{}", naif_center)),
            ("frame_id", self.frame_id.to_string()),
            ("data_type", format!("{}", naif_type.to_string())),
            ("initial_addr", self.initial_addr.to_string()),
            ("final_addr", self.final_addr.to_string()),
        ];

        let label_width = fields.iter().map(|(k, _)| k.len()).max().unwrap_or(10);
        let value_width = fields.iter().map(|(_, v)| v.len()).max().unwrap_or(10);

        let border = format!(
            "+{:-<label$}+{:-<value$}+",
            "",
            "",
            label = label_width + 2,
            value = value_width + 2
        );

        writeln!(f, "{border}")?;
        writeln!(
            f,
            "| {:<label_width$} | {:<value_width$} |",
            "Field", "Value",
        )?;
        writeln!(f, "{border}")?;

        for (label, value) in fields {
            writeln!(f, "| {:<label_width$} | {:<value_width$} |", label, value,)?;
        }

        writeln!(f, "{border}")
    }
}

#[cfg(test)]
mod test_summary {
    use super::*;

    #[test]
    fn test_summary_display() {
        let summary = Summary {
            start_epoch: -14200747200.0,
            end_epoch: 20514081600.0,
            target: 3,
            center: 0,
            frame_id: 1,
            data_type: 2,
            initial_addr: 3021513,
            final_addr: 4051108,
        };

        let expected = r#"+--------------+-------------------------+
| Field        | Value                   |
+--------------+-------------------------+
| start_epoch  | 1549-12-31T00:00:00 ET  |
| end_epoch    | 2650-01-25T00:00:00 ET  |
| target       | EarthMoonBarycenter     |
| center       | Solar System Barycenter |
| frame_id     | 1                       |
| data_type    | Chebyshev Position Only |
| initial_addr | 3021513                 |
| final_addr   | 4051108                 |
+--------------+-------------------------+
"#;

        let output = format!("{}", summary);
        assert_eq!(output, expected);
    }
}
