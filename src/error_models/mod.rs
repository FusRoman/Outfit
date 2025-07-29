mod vfcc17;

use std::{collections::HashMap, str::FromStr};

use nom::{
    branch::alt,
    bytes::complete::{tag, take_until, take_while},
    character::complete::{char, multispace0},
    combinator::{map, opt},
    number::complete::float,
    sequence::{preceded, separated_pair, terminated},
    IResult, Parser,
};

use crate::{constants::MpcCode, outfit_errors::OutfitError};
use vfcc17::parse_vfcc17_line;

type CatalogCode = String;
pub type ErrorModelData = HashMap<(MpcCode, CatalogCode), (f32, f32)>;

pub enum ErrorModel {
    FCCT14,
    CBM10,
    VFCC17,
}

static FCCT14_RULES: &str = include_str!("data_models/fcct14.rules");
static CBM10_RULES: &str = include_str!("data_models/cbm10.rules");
static VFCC17_RULES: &str = include_str!("data_models/vfcc17.rules");

pub(in crate::error_models) type ParseResult<'a> =
    IResult<&'a str, Vec<((MpcCode, CatalogCode), (f32, f32))>>;

fn is_alphanum(c: char) -> bool {
    c.is_alphanumeric()
}

fn parse_station(input: &str) -> IResult<&str, &str> {
    terminated(take_while(is_alphanum), tag(":")).parse(input)
}

fn parse_rms_values(input: &str) -> IResult<&str, (f32, f32)> {
    preceded(
        multispace0,
        preceded(
            char('@'),
            separated_pair(
                preceded(multispace0, float),
                char(','),
                preceded(multispace0, float),
            ),
        ),
    )
    .parse(input)
}

fn parse_catalog_codes(input: &str) -> IResult<&str, Vec<String>> {
    preceded(
        multispace0,
        preceded(
            // accepte soit "c=" soit rien du tout
            alt((tag("c="), tag(""))),
            map(
                nom::bytes::complete::take_while(|c: char| c.is_alphabetic() || c == '*'),
                |s: &str| s.chars().map(|c| c.to_string()).collect(),
            ),
        ),
    )
    .parse(input)
}

fn parse_full_line(input: &str) -> ParseResult {
    let (input, remain) = opt(take_until("!")).parse(input)?; //ignore comments
    let input = input.trim();

    let input = remain.unwrap_or(input);

    map(
        (parse_station, parse_catalog_codes, parse_rms_values),
        |(station, catalogs, (rmsa, rmsd))| {
            catalogs
                .into_iter()
                .map(|cat| ((station.to_string(), cat), (rmsa, rmsd)))
                .collect()
        },
    )
    .parse(input)
}

fn parse_full_file<F>(file: &str, parse_line: F) -> Result<ErrorModelData, OutfitError>
where
    F: Fn(&str) -> ParseResult,
{
    let error_map: ErrorModelData = file
        .lines()
        .filter(|line| !line.trim().is_empty() && !line.trim_start().starts_with('!'))
        .map(|line| {
            parse_line(line)
                .map_err(|_e| OutfitError::NomParsingError(line.to_string()))
                .map(|(_, pairs)| pairs)
        })
        .collect::<Result<Vec<_>, OutfitError>>()?
        .into_iter()
        .flatten()
        .collect();

    Ok(error_map)
}

impl ErrorModel {
    pub(crate) fn read_error_model_file(&self) -> Result<ErrorModelData, OutfitError> {
        let error_map: ErrorModelData = match self {
            ErrorModel::FCCT14 => parse_full_file(FCCT14_RULES, parse_full_line)?,
            ErrorModel::CBM10 => parse_full_file(CBM10_RULES, parse_full_line)?,
            ErrorModel::VFCC17 => {
                // Implement parsing logic for VFCC17
                parse_full_file(VFCC17_RULES, parse_vfcc17_line)?
            }
        };

        Ok(error_map)
    }
}

impl FromStr for ErrorModel {
    type Err = OutfitError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "FCCT14" => Ok(ErrorModel::FCCT14),
            "CBM10" => Ok(ErrorModel::CBM10),
            "VFCC17" => Ok(ErrorModel::VFCC17),
            _ => Err(OutfitError::InvalidErrorModel(format!(
                "Invalid error model: {s}"
            ))),
        }
    }
}

pub(crate) fn get_bias_rms(
    error_model: &ErrorModelData,
    mpc_code: MpcCode,
    catalog_code: CatalogCode,
) -> Option<(f32, f32)> {
    error_model
        .get(&(mpc_code.clone(), catalog_code.clone()))
        .cloned()
        .or_else(|| {
            error_model
                .get(&(mpc_code.clone(), "e".to_string()))
                .cloned()
        })
        .or_else(|| error_model.get(&(mpc_code, "c".to_string())).cloned())
        .or_else(|| error_model.get(&("ALL".to_string(), catalog_code)).cloned())
        .or_else(|| {
            error_model
                .get(&("ALL".to_string(), "e".to_string()))
                .cloned()
        })
        .or_else(|| {
            error_model
                .get(&("ALL".to_string(), "c".to_string()))
                .cloned()
        })
}

impl TryFrom<&str> for ErrorModel {
    type Error = OutfitError;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        value.parse()
    }
}

#[cfg(test)]
mod test_error_model {
    use super::*;

    #[test]
    fn test_parse_fcct14_line() {
        let line = "ALL:  c=eqru @ 0.33, 0.30";
        let result = parse_full_line(line);
        assert!(result.is_ok());
        let ((mpc_code, catalog_code), (rmsa, rmsd)) = &result.unwrap().1[0];
        assert_eq!(mpc_code, "ALL");
        assert_eq!(catalog_code, "e");
        assert_eq!(*rmsa, 0.33);
        assert_eq!(*rmsd, 0.3);

        let line = "ALL:  c=cd   @ 0.51, 0.40 ! CBM Generic Catalog weights";
        let result = parse_full_line(line);
        assert!(result.is_ok());
        let ((mpc_code, catalog_code), (rmsa, rmsd)) = &result.unwrap().1[0];
        assert_eq!(mpc_code, "ALL");
        assert_eq!(catalog_code, "c");
        assert_eq!(*rmsa, 0.51);
        assert_eq!(*rmsd, 0.4);

        let line = "699:c  @ 0.93, 0.78";
        let result = parse_full_line(line);
        assert!(result.is_ok());
        let ((mpc_code, catalog_code), (rmsa, rmsd)) = &result.unwrap().1[0];
        assert_eq!(mpc_code, "699");
        assert_eq!(catalog_code, "c");
        assert_eq!(*rmsa, 0.93);
        assert_eq!(*rmsd, 0.78);
    }

    #[test]
    fn test_read_error_model_file() {
        let error_model = ErrorModel::FCCT14;
        let result = error_model.read_error_model_file();
        assert!(result.is_ok());
        let data = result.unwrap();
        assert!(!data.is_empty());

        let error_model = ErrorModel::CBM10;
        let result = error_model.read_error_model_file();
        assert!(result.is_ok());
        let data = result.unwrap();
        assert!(!data.is_empty());

        let error_model = ErrorModel::VFCC17;
        let result = error_model.read_error_model_file();

        assert!(result.is_ok());
        let data = result.unwrap();
        assert!(!data.is_empty());
    }

    #[test]
    fn test_get_bias_rms() {
        let error_model = ErrorModel::FCCT14.read_error_model_file().unwrap();
        let bias_rms = get_bias_rms(&error_model, "ALL".to_string(), "c".to_string());
        assert!(bias_rms.is_some());
        let (rmsa, rmsd) = bias_rms.unwrap();
        assert_eq!(rmsa, 0.51);
        assert_eq!(rmsd, 0.4);

        let bias_rms = get_bias_rms(&error_model, "699".to_string(), "c".to_string());
        assert!(bias_rms.is_some());
        let (rmsa, rmsd) = bias_rms.unwrap();
        assert_eq!(rmsa, 0.47);
        assert_eq!(rmsd, 0.39);

        let error_model = ErrorModel::CBM10.read_error_model_file().unwrap();
        let bias_rms = get_bias_rms(&error_model, "ALL".to_string(), "c".to_string());
        assert!(bias_rms.is_some());
        let (rmsa, rmsd) = bias_rms.unwrap();
        assert_eq!(rmsa, 0.5);
        assert_eq!(rmsd, 0.5);

        let bias_rms = get_bias_rms(&error_model, "699".to_string(), "c".to_string());
        assert!(bias_rms.is_some());
        let (rmsa, rmsd) = bias_rms.unwrap();
        assert_eq!(rmsa, 0.84);
        assert_eq!(rmsd, 0.81);

        let error_model = ErrorModel::VFCC17.read_error_model_file().unwrap();
        let bias_rms = get_bias_rms(&error_model, "ALL".to_string(), "U".to_string());
        assert!(bias_rms.is_some());
        let (rmsa, rmsd) = bias_rms.unwrap();
        assert_eq!(rmsa, 0.6);
        assert_eq!(rmsd, 0.6);
        let bias_rms = get_bias_rms(&error_model, "699".to_string(), "*".to_string());
        assert!(bias_rms.is_some());
        let (rmsa, rmsd) = bias_rms.unwrap();
        assert_eq!(rmsa, 0.8);
        assert_eq!(rmsd, 0.8);
    }
}
