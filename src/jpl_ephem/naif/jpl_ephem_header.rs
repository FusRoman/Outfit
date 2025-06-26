use nom::{
    bytes::complete::{tag, take_until},
    character::complete::{line_ending, not_line_ending},
    number::complete::double,
    IResult, Parser,
};

#[derive(Debug, PartialEq, Clone)]
pub struct JPLEphemHeader {
    pub version: String,
    pub creation_date: String,
    pub start_ephem: String,
    pub end_ephem: String,
    pub start_jd: f64,
    pub end_jd: f64,
}

impl JPLEphemHeader {
    fn parse_version(input: &str) -> IResult<&str, &str> {
        let (input, _) = take_until("JPL planetary and lunar ephemeris")(input)?;
        let (input, _) = tag("JPL planetary and lunar ephemeris ")(input)?;
        let (input, version) = not_line_ending(input)?;
        let (input, _) = line_ending(input)?;

        Ok((input, version.trim()))
    }

    fn parse_creation_date(input: &str) -> IResult<&str, &str> {
        let (input, _) = take_until("Integrated ")(input)?;
        let (input, _) = tag("Integrated ")(input)?;
        let (input, creation_date) = not_line_ending(input)?;
        let (input, _) = line_ending(input)?;

        Ok((input, creation_date.trim()))
    }

    fn parse_date_range(input: &str) -> IResult<&str, (&str, &str)> {
        let (input, _) = take_until("Time span covered by ephemeris:")(input)?;
        let (input, _) = tag("Time span covered by ephemeris:")(input)?;
        let (input, _) = (line_ending, line_ending).parse(input)?;
        let (input, (start_ephem, _, end_ephem)) =
            (take_until("to"), tag("to   "), not_line_ending).parse(input)?;

        Ok((input, (start_ephem.trim(), end_ephem.trim())))
    }

    fn parse_jd_range(input: &str) -> IResult<&str, (f64, f64)> {
        let (input, (_, _, start_jd, _, end_jd)) = (
            line_ending,
            tag("JD   "),
            |s| double(s),
            tag("   to   JD   "),
            |s| double(s),
        )
            .parse(input)?;

        Ok((input, (start_jd, end_jd)))
    }

    pub fn parse(input: &str) -> IResult<&str, Self> {
        let (input, version) = Self::parse_version(input)?;
        let (input, creation_date) = Self::parse_creation_date(input)?;

        let (input, (start_ephem, end_ephem)) = Self::parse_date_range(input)?;
        let (input, (start_jd, end_jd)) = Self::parse_jd_range(input)?;
        Ok((
            input,
            JPLEphemHeader {
                version: version.to_string(),
                creation_date: creation_date.trim().to_string(),
                start_ephem: start_ephem.trim().to_string(),
                end_ephem: end_ephem.trim().to_string(),
                start_jd,
                end_jd,
            },
        ))
    }
}

use std::fmt;

impl fmt::Display for JPLEphemHeader {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        const LABEL_WIDTH: usize = 20;
        const VALUE_WIDTH: usize = 30;

        let border = format!(
            "+{:-<label$}+{:-<value$}+",
            "",
            "",
            label = LABEL_WIDTH + 2,
            value = VALUE_WIDTH + 2
        );

        writeln!(
            f,
            "+{:^label$}+{:^value$}+",
            "JPL Ephemeris Header",
            "",
            label = LABEL_WIDTH + 2,
            value = VALUE_WIDTH + 2
        )?;
        writeln!(f, "{}", border)?;
        writeln!(
            f,
            "| {:<label$} | {:<value$} |",
            "Version",
            self.version,
            label = LABEL_WIDTH,
            value = VALUE_WIDTH
        )?;
        writeln!(
            f,
            "| {:<label$} | {:<value$} |",
            "Creation Date",
            self.creation_date,
            label = LABEL_WIDTH,
            value = VALUE_WIDTH
        )?;
        writeln!(
            f,
            "| {:<label$} | {:<value$} |",
            "Start Ephem",
            self.start_ephem,
            label = LABEL_WIDTH,
            value = VALUE_WIDTH
        )?;
        writeln!(
            f,
            "| {:<label$} | {:<value$} |",
            "End Ephem",
            self.end_ephem,
            label = LABEL_WIDTH,
            value = VALUE_WIDTH
        )?;
        writeln!(
            f,
            "| {:<label$} | {:<value$} |",
            "Start JD",
            format!("{:.6}", self.start_jd),
            label = LABEL_WIDTH,
            value = VALUE_WIDTH
        )?;
        writeln!(
            f,
            "| {:<label$} | {:<value$} |",
            "End JD",
            format!("{:.6}", self.end_jd),
            label = LABEL_WIDTH,
            value = VALUE_WIDTH
        )?;
        writeln!(f, "{}", border)?;

        Ok(())
    }
}


#[cfg(test)]
mod test_jpl_header {
    use super::*;

    #[test]
    fn test_jpl_header_display() {
        let header = JPLEphemHeader {
            version: "DE440".to_string(),
            creation_date: "25 June 2020".to_string(),
            start_ephem: "31-DEC-1549 00:00".to_string(),
            end_ephem: "25-JAN-2650 00:00".to_string(),
            start_jd: 2287184.5,
            end_jd: 2688976.5
        };

        let expected = r#"+ JPL Ephemeris Header +                                +
+----------------------+--------------------------------+
| Version              | DE440                          |
| Creation Date        | 25 June 2020                   |
| Start Ephem          | 31-DEC-1549 00:00              |
| End Ephem            | 25-JAN-2650 00:00              |
| Start JD             | 2287184.500000                 |
| End JD               | 2688976.500000                 |
+----------------------+--------------------------------+
"#;
        let output = format!("{}", header);
        assert_eq!(output, expected);
    }
}