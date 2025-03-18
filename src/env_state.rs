use crate::observers::observers::init_observatories;
use hifitime::ut1::Ut1Provider;
use once_cell::sync::OnceCell;
use reqwest::Client;
use tokio::task;
use crate::constants::MpcCodeObs;

/// This object is passed to the various functions in the library
/// to provide access to the state of the library
/// 
/// # Fields
///
/// * `http_client` - A reqwest client used to make HTTP requests
/// * `ut1_provider` - A provider used to get the current UT1 time
/// * `observatories` - A lazy map of observatories from the Minor Planet Center.
///     The key is the MPC code and the value is the observer
pub struct OutfitState {
    pub http_client: Client,
    pub ut1_provider: Ut1Provider,
    observatories: OnceCell<MpcCodeObs>,
}

impl OutfitState {
    /// Create a new OutfitState object
    ///
    /// Return
    /// ------
    /// * A new OutfitState object
    ///     - The UT1 provider is downloaded from the JPL
    ///     - The HTTP client is created with default settings
    ///     - The observatories are lazily loaded from the Minor Planet Center
    pub async fn new() -> Self {
        let ut1_provider = task::spawn_blocking(|| Ut1Provider::download_short_from_jpl().unwrap())
            .await
            .expect("Download of the JPL short time scale UT1 data failed");
        OutfitState {
            http_client: Client::new(),
            ut1_provider: ut1_provider,
            observatories: OnceCell::new(),
        }
    }

    /// Get the observatories from the Minor Planet Center
    /// 
    /// Return
    /// ------
    /// * A map of observatories from the Minor Planet Center
    ///    The key is the MPC code and the value is the observer
    pub fn get_observatories(&self) -> &MpcCodeObs {
        self.observatories.get_or_init(|| {
            task::block_in_place(|| {
                tokio::runtime::Handle::current().block_on(async {
                    init_observatories(self).await
                })
            })
        })
    }
}
