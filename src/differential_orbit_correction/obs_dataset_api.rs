use crate::{FullOrbitResult, OrbitalElements, OutfitError};

pub trait FitLSQ {
    fn fit_lsq(
        &self,
        initial_orbit: Option<OrbitalElements>,
    ) -> Result<FullOrbitResult, OutfitError>;
}
