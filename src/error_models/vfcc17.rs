use std::collections::HashMap;

use nom::{
    branch::alt,
    bytes::complete::{tag, take_until, take_while, take_while1},
    character::complete::{char, multispace0, space1},
    combinator::{map, opt},
    number::complete::float,
    sequence::{preceded, separated_pair, terminated},
    IResult, Parser,
};

use crate::{constants::MpcCode, outfit_errors::OutfitError};

use super::CatalogCode;

fn is_word_char(c: char) -> bool {
    c.is_alphanumeric() || "*".contains(c)
}

fn parse_word(input: &str) -> IResult<&str, &str> {
    take_while1(is_word_char)(input)
}

fn parse_station(input: &str) -> IResult<&str, &str> {
    parse_word(input)
}

fn parse_catalog_codes(input: &str) -> IResult<&str, Vec<String>> {
    preceded(
        tag("c="),
        map(take_while1(is_word_char), |s: &str| {
            s.chars().map(|c| c.to_string()).collect()
        }),
    )
    .parse(input)
}

fn parse_rms_values(input: &str) -> IResult<&str, (f32, f32)> {
    preceded(
        take_until("@"),
        preceded(
            tag("@"),
            separated_pair(
                preceded(multispace0, float),
                preceded(multispace0, char(',')),
                preceded(multispace0, float),
            ),
        ),
    )
    .parse(input)
}

pub fn parse_vfcc17_line(input: &str) -> IResult<&str, Vec<((MpcCode, CatalogCode), (f32, f32))>> {
    let (input, remain) = opt(take_until("!")).parse(input)?; // Ignore comments
    let input = remain.unwrap_or(input).trim();

    map(
        (
            parse_station,
            take_until("c="),
            parse_catalog_codes,
            parse_rms_values,
        ),
        |(station, _, catalogs, (rmsa, rmsd))| {
            catalogs
                .into_iter()
                .map(|cat| ((station.to_string(), cat), (rmsa, rmsd)))
                .collect()
        },
    )
    .parse(input)
}

#[cfg(test)]
mod test_vfcc17_parser {
    use super::*;

    #[test]
    fn test_vfcc17_parser() {
        let input = "ALL t=cBCVn c=*          p=    >            <            @  1.00,  1.00 ! Unknown catalog";
        let result = parse_vfcc17_line(input);
        assert!(result.is_ok());
        let (_, parsed) = result.unwrap();
        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0].0 .0, "ALL");
        assert_eq!(parsed[0].0 .1, "*");
        assert_eq!(parsed[0].1 .0, 1.);
        assert_eq!(parsed[0].1 .1, 1.);

        let input = "568 t=cC    c=t          p=_   >            <            @  0.20,  0.20  ! Micheli updated ";
        let result = parse_vfcc17_line(input);
        assert!(result.is_ok());
        let (_, parsed) = result.unwrap();
        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0].0 .0, "568");
        assert_eq!(parsed[0].0 .1, "t");
        assert_eq!(parsed[0].1 .0, 0.2);
        assert_eq!(parsed[0].1 .1, 0.2);
    }
}
