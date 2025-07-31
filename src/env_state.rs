//! # Outfit environment state
//!
//! This module defines [`crate::env_state::OutfitEnv`], the **shared environment object** used across
//! the `Outfit` library. It provides access to:
//!
//! - A persistent **HTTP client** (for downloading ephemerides, observatory lists, etc.).
//! - A **UT1 provider** from [hifitime](https://docs.rs/hifitime) to handle Earth rotation
//!   parameters from JPL.
//!
//! This object is designed to be **cheaply cloneable** and passed to algorithms
//! that require access to external data sources or Earth orientation models.
//!
//! ## Overview
//!
//! The main responsibilities of `OutfitEnv` are:
//!
//! 1. Manage a global [`ureq::Agent`] HTTP client with sensible default settings.
//! 2. Download and initialize an [`hifitime::ut1::Ut1Provider`] from JPL’s `latest_eop2.long` file
//!    (Earth orientation parameters) at startup.
//! 3. Provide simple utilities for performing HTTP GET requests.
//!
//! ## Structure
//!
//! ```text
//! OutfitEnv
//! ├── http_client  (ureq::Agent)
//! └── ut1_provider (hifitime::Ut1Provider)
//! ```
//!
//! ## Usage
//!
//! ```rust,ignore
//! use outfit::env_state::OutfitEnv;
//!
//! // Create a new environment (downloads UT1 data from JPL)
//! let env = OutfitEnv::new();
//!
//! // Access the UT1 provider
//! let ut1 = &env.ut1_provider;
//!
//! // Make a GET request using the built-in HTTP client
//! let response = env.get_from_url("https://ssd.jpl.nasa.gov/api/horizons.api");
//! println!("Response: {}", &response[..100.min(response.len())]);
//! ```
//!
//! ## Notes
//!
//! - The [`crate::env_state::OutfitEnv`] struct is meant to be reused and shared between different
//!   parts of the crate to avoid redundant downloads and HTTP session creation.
//! - The UT1 provider is initialized once at startup; if fresh data is needed,
//!   the library must be restarted or the provider re-downloaded manually.
//!
//! ## See also
//!
//! - [`hifitime::ut1::Ut1Provider`] – Manages Earth orientation and UT1 corrections.
//! - [`ureq::Agent`] – Minimal HTTP client used internally.
use hifitime::ut1::Ut1Provider;
use std::convert::TryFrom;
use std::{fmt::Debug, time::Duration};
use ureq::{
    http::{self, Uri},
    Agent,
};

/// This object is passed to the various functions in the library
/// to provide access to the state of the library
///
/// # Fields
///
/// * `http_client` - A reqwest client used to make HTTP requests
/// * `ut1_provider` - A provider used to get the current UT1 time
/// * `observatories` - A lazy map of observatories from the Minor Planet Center.
///   The key is the MPC code and the value is the observer
#[derive(Debug, Clone)]
pub struct OutfitEnv {
    pub http_client: Agent,
    pub ut1_provider: Ut1Provider,
}

impl Default for OutfitEnv {
    fn default() -> Self {
        Self::new()
    }
}

impl OutfitEnv {
    /// Create a new Outfit object
    ///
    /// Return
    /// ------
    /// * A new Outfit object
    ///     - The UT1 provider is downloaded from the JPL
    ///     - The HTTP client is created with default settings
    ///     - The observatories are lazily loaded from the Minor Planet Center
    pub fn new() -> Self {
        let ut1_provider = OutfitEnv::initialize_ut1_provider();

        let config = Agent::config_builder()
            .timeout_global(Some(Duration::from_secs(10)))
            .build();
        let agent: Agent = config.into();

        OutfitEnv {
            http_client: agent,
            ut1_provider,
        }
    }

    fn initialize_ut1_provider() -> Ut1Provider {
        Ut1Provider::download_from_jpl("latest_eop2.long")
            .expect("Download of the JPL short time scale UT1 data failed")
    }

    pub(crate) fn get_from_url<U>(&self, url: U) -> String
    where
        Uri: TryFrom<U>,
        <Uri as TryFrom<U>>::Error: Into<http::Error>,
    {
        self.http_client
            .get(url)
            .call()
            .expect("Get request failed")
            .body_mut()
            .read_to_string()
            .expect("Failed to read response body")
    }
}
