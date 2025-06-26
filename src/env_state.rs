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
///     The key is the MPC code and the value is the observer
#[derive(Debug, Clone)]
pub struct OutfitEnv {
    pub http_client: Agent,
    pub ut1_provider: Ut1Provider,
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
            ut1_provider: ut1_provider,
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

    pub(crate) fn post_from_url<T, I, K, V>(&self, url: T, form: I) -> String
    where
        Uri: TryFrom<T>,
        <Uri as TryFrom<T>>::Error: Into<http::Error>,
        I: IntoIterator<Item = (K, V)>,
        K: AsRef<str>,
        V: AsRef<str>,
    {
        self.http_client
            .post(url)
            .send_form(form)
            .expect("Post request failed")
            .body_mut()
            .read_to_string()
            .expect("Failed to read response body")
    }
}
