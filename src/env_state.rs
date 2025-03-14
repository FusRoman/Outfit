use crate::observers::observers::init_observatories;
use hifitime::ut1::Ut1Provider;
use once_cell::sync::OnceCell;
use reqwest::Client;
use tokio::task;
use crate::constants::MpcCodeObs;

/// Object used to store internal state for the orbit estimation
pub struct OutfitState {
    pub http_client: Client,
    pub ut1_provider: Ut1Provider,
    observatories: OnceCell<MpcCodeObs>,
}

impl OutfitState {
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
