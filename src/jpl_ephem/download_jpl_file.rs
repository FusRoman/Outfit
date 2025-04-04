use camino::Utf8Path;
use tokio::{fs::File, io::AsyncWriteExt};

use tokio_stream::StreamExt;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error + Send + Sync>>;

pub(super) async fn download_big_file(url: &str, path: &Utf8Path) -> Result<()> {
    let mut file = File::create(path).await?;
    println!("Downloading {}...", url);

    let mut stream = reqwest::get(url).await?.bytes_stream();

    while let Some(chunk_result) = stream.next().await {
        let chunk = chunk_result?;
        file.write_all(&chunk).await?;
    }

    file.flush().await?;

    println!("Downloaded {}", url);
    Ok(())
}
