# Outfit

Outfit is a Rust library for managing, analyzing, and determining the orbits of celestial objects from astrometric observations. It provides efficient tools to read, manipulate, and process observations from various formats (80-column MPC, ADES XML, Parquet), perform initial orbit determination (Gauss method), manage observers, and interface with JPL ephemerides.

## Main Features

- **Observation reading**:
  - 80-column MPC files
  - ADES (XML) files
  - Parquet files (batch observations)
- **Observer management**: creation, identification by MPC code, handling coordinates and altitudes
- **Initial orbit determination**:
  - Gauss method on observation triplets
  - Calculations in different reference systems
- **Error handling and uncertainty models**
- **JPL ephemerides interface (HORIZON, NAIF)**
- **Batch processing and benchmarks**
- **Modular and extensible**

## Installation

Add Outfit to your Rust project via Cargo:

```toml
[dependencies]
outfit = "1.0.0"
```

Enable the `jpl-download` feature for automatic ephemeris download:

```toml
[dependencies]
outfit = { version = "1.0.0", features = ["jpl-download"] }
```

## Tests and Benchmarks

- Run tests:
  ```bash
  cargo test
  ```
- Run benchmarks:
  ```bash
  cargo bench
  ```

## Main Dependencies

- [nalgebra] for linear algebra
- [parquet] for Parquet file reading
- [quick-xml] for ADES parsing
- [criterion] for benchmarks

## Roadmap

Planned and upcoming features for Outfit include:

### v2.0
- **Full least squares orbit fitting**  
  Add a complete least squares optimization for orbit fitting, allowing robust adjustment of orbital elements to all available observations.
- **Support for hyperbolic and parabolic orbits (e ≥ 1)**  
  Enable orbit determination and propagation for interstellar objects and other bodies with e ≥ 1. ([#29](https://github.com/FusRoman/Outfit/issues/29))
- **ANISE support for JPL ephemeris loading**  
  Reduce codebase complexity and add support for loading JPL ephemerides via ANISE/SPICE. ([#28](https://github.com/FusRoman/Outfit/issues/28))
- **RMS estimation for all file formats**  
  Extend RMS (root mean square) error estimation to all supported observation file formats. ([#22](https://github.com/FusRoman/Outfit/issues/22))
- **Radar observation support**  
  Add the ability to ingest and process radar observations for orbit determination. ([#13](https://github.com/FusRoman/Outfit/issues/13))

### v3.0
- **Vaisala method for initial orbit determination**  
  Implement the Vaisala method as an alternative approach for initial orbit determination.

For more details and the latest updates, see the [GitHub issues page](https://github.com/FusRoman/Outfit/issues).

## License

This project is licensed under the CeCILL-C Free Software License Agreement. See the `LICENSE` file for more information.

## Author

Developed by FusRoman.

## Acknowledgements

Special thanks to the developers of [OrbFit](https://adams.dm.unipi.it/orbfit/) for their outstanding work and for providing much of the inspiration and reference algorithms for this project.

---

For any questions or contributions, feel free to open an issue or pull request on the GitHub repository.
