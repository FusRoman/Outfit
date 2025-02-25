use hifitime::ut1::Ut1Provider;
use reqwest::Client;
use tokio::task;

/// Object used to store internal state for the orbit estimation
pub struct OutfitState {
    pub http_client: Client,
    pub ut1_provider: Ut1Provider,
}

impl OutfitState {
    pub async fn new() -> Self {
        let ut1_provider = task::spawn_blocking(|| Ut1Provider::download_short_from_jpl().unwrap())
            .await
            .expect("Download of the JPL short time scale UT1 data failed");
        OutfitState {
            http_client: Client::new(),
            ut1_provider: ut1_provider,
        }
    }
}
