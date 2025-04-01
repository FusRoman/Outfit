use outfit::outfit::Outfit;

#[test]
fn test_outfit_observer_management() {
    let mut outfit = Outfit::new();

    let obs1 = outfit.new_observer(51.58206, -73.06644, 100., Some("Test Observer 1".into()));
    assert_eq!(obs1.name, Some("Test Observer 1".to_string()));
    assert_eq!(obs1.longitude, 51.58206);
    assert_eq!(obs1.rho_cos_phi, 0.29216347396649495);
    assert_eq!(obs1.rho_sin_phi, -0.9531782585730825);

    let obs2 = outfit.new_observer(52.58206, 23.4587, 1423., Some("Test Observer 2".into()));
    assert_eq!(obs2.name, Some("Test Observer 2".to_string()));
    assert_eq!(obs2.longitude, 52.58206);
    assert_eq!(obs2.rho_cos_phi, 0.9180389162887692);
    assert_eq!(obs2.rho_sin_phi, 0.39572170773696747);
}
