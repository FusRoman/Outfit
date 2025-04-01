use std::fmt::Debug;
use hifitime::ut1::Ut1Provider;
use reqwest::{Client, IntoUrl};
use serde::Serialize;
use tokio::task;

/// This object is passed to the various functions in the library
/// to provide access to the state of the library
///
/// # Fields
///
/// * `http_client` - A reqwest client used to make HTTP requests
/// * `ut1_provider` - A provider used to get the current UT1 time
/// * `observatories` - A lazy map of observatories from the Minor Planet Center.
///     The key is the MPC code and the value is the observer
#[derive(Debug)]
pub struct OutfitEnv {
    pub http_client: Client,
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
        OutfitEnv {
            http_client: Client::new(),
            ut1_provider: ut1_provider,
        }
    }

    fn initialize_ut1_provider() -> Ut1Provider {
        task::block_in_place(|| {
            Ut1Provider::download_short_from_jpl()
                .expect("Download of the JPL short time scale UT1 data failed")
        })
    }

    pub(crate) fn get_from_url<U>(&self, url: U) -> String
    where
        U: IntoUrl,
    {
        task::block_in_place(|| {
            tokio::runtime::Handle::current().block_on(async {
                self.http_client
                    .get(url)
                    .send()
                    .await
                    .expect("Get request failed")
                    .text()
                    .await
                    .expect("Failed to get text from get request")
            })
        })
    }

    pub(crate) fn post_from_url<U, T: Serialize + ?Sized>(&self, url: U, form: &T) -> String
    where
        U: IntoUrl,
    {
        task::block_in_place(|| {
            tokio::runtime::Handle::current().block_on(async {
                self.http_client
                    .post(url)
                    .form(form)
                    .send()
                    .await
                    .expect("Post request failed")
                    .text()
                    .await
                    .expect("Failed to get text from post request")
            })
        })
    }
}
