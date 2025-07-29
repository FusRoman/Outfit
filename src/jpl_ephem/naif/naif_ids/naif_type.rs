use std::{convert::TryFrom, fmt};

use crate::outfit_errors::OutfitError;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum SpkDataType {
    ModifiedDifferenceArray = 1,
    ChebyshevPositionOnly = 2,
    ChebyshevPositionVelocity = 3,
    Reserved4 = 4,
    TwoBodyDiscreteStates = 5,
    Reserved6 = 6,
    Reserved7 = 7,
    EquallySpacedLagrange = 8,
    UnequallySpacedLagrange = 9,
    TwoLineElements = 10,
    Reserved11 = 11,
    HermiteUniform = 12,
    HermiteNonUniform = 13,
    ChebyshevNonUniform = 14,
    PrecessingConic = 15,
    Reserved16 = 16,
    EquinoctialElements = 17,
    ESAHermiteLagrange = 18,
    ESAPiecewiseInterpolation = 19,
    ChebyshevVelocityOnly = 20,
    ExtendedModifiedDifferenceArray = 21,
}

impl SpkDataType {
    pub fn from_i32(value: i32) -> Result<Self, OutfitError> {
        SpkDataType::try_from(value)
    }

    pub fn to_i32(self) -> i32 {
        self as i32
    }
}

impl From<SpkDataType> for i32 {
    fn from(data_type: SpkDataType) -> Self {
        data_type as i32
    }
}

impl TryFrom<i32> for SpkDataType {
    type Error = OutfitError;

    fn try_from(value: i32) -> Result<Self, Self::Error> {
        use SpkDataType::*;
        match value {
            1 => Ok(ModifiedDifferenceArray),
            2 => Ok(ChebyshevPositionOnly),
            3 => Ok(ChebyshevPositionVelocity),
            4 => Ok(Reserved4),
            5 => Ok(TwoBodyDiscreteStates),
            6 => Ok(Reserved6),
            7 => Ok(Reserved7),
            8 => Ok(EquallySpacedLagrange),
            9 => Ok(UnequallySpacedLagrange),
            10 => Ok(TwoLineElements),
            11 => Ok(Reserved11),
            12 => Ok(HermiteUniform),
            13 => Ok(HermiteNonUniform),
            14 => Ok(ChebyshevNonUniform),
            15 => Ok(PrecessingConic),
            16 => Ok(Reserved16),
            17 => Ok(EquinoctialElements),
            18 => Ok(ESAHermiteLagrange),
            19 => Ok(ESAPiecewiseInterpolation),
            20 => Ok(ChebyshevVelocityOnly),
            21 => Ok(ExtendedModifiedDifferenceArray),
            _ => Err(OutfitError::InvalidSpkDataType(value)),
        }
    }
}

impl fmt::Display for SpkDataType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            SpkDataType::ModifiedDifferenceArray => "Modified Difference Array",
            SpkDataType::ChebyshevPositionOnly => "Chebyshev Position Only",
            SpkDataType::ChebyshevPositionVelocity => "Chebyshev Position Velocity",
            SpkDataType::Reserved4 => "Reserved 4",
            SpkDataType::TwoBodyDiscreteStates => "Two Body Discrete States",
            SpkDataType::Reserved6 => "Reserved 6",
            SpkDataType::Reserved7 => "Reserved 7",
            SpkDataType::EquallySpacedLagrange => "Equally Spaced Lagrange",
            SpkDataType::UnequallySpacedLagrange => "Unequally Spaced Lagrange",
            SpkDataType::TwoLineElements => "Two Line Elements",
            SpkDataType::Reserved11 => "Reserved 11",
            SpkDataType::HermiteUniform => "Hermite Uniform",
            SpkDataType::HermiteNonUniform => "Hermite Non Uniform",
            SpkDataType::ChebyshevNonUniform => "Chebyshev Non Uniform",
            SpkDataType::PrecessingConic => "Precessing Conic",
            SpkDataType::Reserved16 => "Reserved 16",
            SpkDataType::EquinoctialElements => "Equinoctial Elements",
            SpkDataType::ESAHermiteLagrange => "ESA Hermite Lagrange",
            SpkDataType::ESAPiecewiseInterpolation => "ESA Piecewise Interpolation",
            SpkDataType::ChebyshevVelocityOnly => "Chebyshev Velocity Only",
            SpkDataType::ExtendedModifiedDifferenceArray => "Extended Modified Difference Array",
        };
        write!(f, "{s}")
    }
}

#[cfg(test)]
mod test_spk_type {
    use super::*;

    #[test]
    fn test_from() {
        assert_eq!(
            SpkDataType::from_i32(1),
            Ok(SpkDataType::ModifiedDifferenceArray)
        );
        assert_eq!(
            SpkDataType::from_i32(2),
            Ok(SpkDataType::ChebyshevPositionOnly)
        );
        assert_eq!(
            SpkDataType::from_i32(3),
            Ok(SpkDataType::ChebyshevPositionVelocity)
        );
        assert_eq!(SpkDataType::from_i32(4), Ok(SpkDataType::Reserved4));
        assert_eq!(
            SpkDataType::from_i32(5),
            Ok(SpkDataType::TwoBodyDiscreteStates)
        );
        assert_eq!(SpkDataType::from_i32(6), Ok(SpkDataType::Reserved6));
        assert_eq!(SpkDataType::from_i32(7), Ok(SpkDataType::Reserved7));
        assert_eq!(
            SpkDataType::from_i32(8),
            Ok(SpkDataType::EquallySpacedLagrange)
        );
        assert_eq!(
            SpkDataType::from_i32(9),
            Ok(SpkDataType::UnequallySpacedLagrange)
        );
        assert_eq!(SpkDataType::from_i32(10), Ok(SpkDataType::TwoLineElements));
        assert_eq!(SpkDataType::from_i32(11), Ok(SpkDataType::Reserved11));
        assert_eq!(SpkDataType::from_i32(12), Ok(SpkDataType::HermiteUniform));
        assert_eq!(
            SpkDataType::from_i32(13),
            Ok(SpkDataType::HermiteNonUniform)
        );
        assert_eq!(
            SpkDataType::from_i32(14),
            Ok(SpkDataType::ChebyshevNonUniform)
        );
        assert_eq!(SpkDataType::from_i32(15), Ok(SpkDataType::PrecessingConic));
        assert_eq!(SpkDataType::from_i32(16), Ok(SpkDataType::Reserved16));
        assert_eq!(
            SpkDataType::from_i32(17),
            Ok(SpkDataType::EquinoctialElements)
        );
        assert_eq!(
            SpkDataType::from_i32(18),
            Ok(SpkDataType::ESAHermiteLagrange)
        );
        assert_eq!(
            SpkDataType::from_i32(19),
            Ok(SpkDataType::ESAPiecewiseInterpolation)
        );
        assert_eq!(
            SpkDataType::from_i32(20),
            Ok(SpkDataType::ChebyshevVelocityOnly)
        );
        assert_eq!(
            SpkDataType::from_i32(21),
            Ok(SpkDataType::ExtendedModifiedDifferenceArray)
        );
        assert_eq!(
            SpkDataType::from_i32(22),
            Err(OutfitError::InvalidSpkDataType(22))
        );
    }

    #[test]
    fn test_into() {
        let spk_type = SpkDataType::ChebyshevPositionOnly;
        let type_id: i32 = spk_type.into();
        assert_eq!(type_id, 2);

        let spk_type = SpkDataType::ChebyshevPositionVelocity;
        let type_id: i32 = spk_type.into();
        assert_eq!(type_id, 3);
    }

    #[test]
    fn test_to_string() {
        let spk_type = SpkDataType::ChebyshevPositionOnly;
        assert_eq!(spk_type.to_string(), "Chebyshev Position Only");

        let spk_type = SpkDataType::ChebyshevPositionVelocity;
        assert_eq!(spk_type.to_string(), "Chebyshev Position Velocity");
    }
}
