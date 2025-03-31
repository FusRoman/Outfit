use camino::Utf8Path;
use outfit::{observations::ades_reader::parse_ades, outfit::Outfit};

#[tokio::test(flavor = "multi_thread", worker_threads = 2)]
async fn test_read_ades() {
    let mut outfit = Outfit::new();

    let ades = parse_ades(&mut outfit, &Utf8Path::new("tests/data/example_ades.xml"));

    println!("{:?}", ades);
}
