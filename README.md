# Outfit - An Orbit Solver Package for Sun-Orbiting Bodies

Outfit is a Rust package designed to estimate the Keplerian orbital parameters for bodies orbiting around the Sun using the Gauss method.

## Features

- **Gauss Method**: Implements the Gauss method to estimate orbital elements.
- **Keplerian Elements**: Outputs the six Keplerian orbital elements.
- **Asynchronous Requests**: Utilizes asynchronous requests for data fetching.

## Inputs

The inputs required are:
- **Right Ascension (RA)**: In degrees.
- **Declination (DEC)**: In degrees.
- **Observation Time**: In Modified Julian Date (MJD).

## Outputs

The output consists of the following Keplerian orbital elements:
- **Semi-Major Axis**: The longest diameter of the elliptical orbit.
- **Eccentricity**: The deviation of the orbit from a perfect circle.
- **Inclination**: The tilt of the orbit's plane with respect to the reference plane.
- **Longitude of the Ascending Node**: The angle from the reference direction to the ascending node.
- **Argument of Periapsis**: The angle from the ascending node to the periapsis.
- **Mean Anomaly**: The fraction of the orbital period that has elapsed since the last periapsis.

## Installation

To include Outfit in your project, add the following to your `Cargo.toml`:

```toml
[dependencies]
outfit = "0.1.0"
```


### Usage
Here is a basic example of how to use Outfit:

```rust
use outfit::gauss::GaussObs;
use nalgebra::Vector3;

#[tokio::main]
async fn main() {
    let ra = Vector3::new(1.6893715963476696, 1.6898894500811472, 1.7527345385664372);
    let dec = Vector3::new(1.0824680373855251, 0.94358050479462163, 0.82737624078999861);
    let time = Vector3::new(57028.479297592596, 57049.245147592592, 57063.977117592593);
    let gauss_obs = GaussObs::new(ra, dec, time, 203.744090000, 20.707233557, 3067.694).await;
    let orbit = gauss_obs.prelim_orbit().unwrap();
    println!("{:?}", orbit);
}
```

### License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

### Acknowledgements
Special thanks to the developers and contributors of the dependencies used in this project. Additional special thanks to the developers of the OrbFit software, which directly inspired the Outfit package. Link of the original OrbFit software: [OrbFit](http://adams.dm.unipi.it/orbfit/)
