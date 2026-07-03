# wind-profile

[![CI](https://github.com/NewEarthTech/wind-profile/actions/workflows/ci.yml/badge.svg)](https://github.com/NewEarthTech/wind-profile/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/wind-profile.svg)](https://crates.io/crates/wind-profile)
[![docs.rs](https://docs.rs/wind-profile/badge.svg)](https://docs.rs/wind-profile)
[![license](https://img.shields.io/crates/l/wind-profile.svg)](#license)
[![MSRV](https://img.shields.io/badge/MSRV-1.83-blue.svg)](https://www.rust-lang.org)

Vertical wind-profile extrapolation for wind energy in Rust: power law, log law, shear exponent from two model heights, and the IEC 61400-12-1 air density and density-corrected wind speed.

## Why this exists

Weather models publish wind at a handful of fixed heights: 10 m always, and depending on the model an 80 m or 100 m level (GFS carries 80 m, the ECMWF open data carries 100 m). Turbine hubs sit at 80-160 m and climbing. Power goes with the cube of wind speed below rated, so a 5% error in the extrapolated speed is roughly a 15% error in energy; the vertical profile is where the money is. This crate is that step done correctly: derive the local shear exponent from two model heights when the model gives them (which beats any tabulated default), extrapolate by the power law or the neutral log law, and put the result on a power-curve footing with the IEC density correction.

## Example

```rust
use wind_profile as wp;

// A GFS-style pair: 6.2 m/s at 10 m, 8.9 m/s at 80 m. Derive the local shear...
let alpha = wp::shear_exponent_from_two_levels(6.2, 10.0, 8.9, 80.0);

// ...carry the 80 m wind up to a 120 m hub...
let u_hub = wp::wind_speed_m_s_power_law(8.9, 80.0, 120.0, alpha);
assert!((u_hub - 9.55).abs() < 0.01);

// ...and normalize it against a power curve quoted at 1.225 kg/m^3.
let rho = wp::air_density_kg_m3(283.15, 99_000.0, 0.8);
let u_pc = wp::density_corrected_wind_speed_m_s(u_hub, rho, wp::REFERENCE_AIR_DENSITY_KG_M3);
assert!(u_pc < u_hub); // warm-sector low pressure: thinner air, lower effective speed
```

## What it covers

- Shear exponent from wind speeds at two heights, the local alternative to a climatological default.
- Power-law extrapolation, with the classic 1/7 neutral onshore exponent as a named constant.
- Neutral log-law extrapolation, with the Davenport-Wieringa roughness classification (`Z0_SEA_M` through `Z0_CHAOTIC_M`) as named constants.
- Air density from temperature, pressure, and relative humidity by the IEC 61400-12-1 equation (the standard's own vapor-pressure fit), plus the standard's `(rho / rho_ref)^(1/3)` density-corrected wind speed for power-curve application.

Both profile laws are neutral-stability forms, the honest baseline when all you have is model output; diabatic (Monin-Obukhov) corrections need surface-flux inputs and are deliberately out of scope. Units are explicit in every function name: `_m_s`, `_m`, `_k`, `_pa`, `_kg_m3`. No dependencies, no `unsafe`.

## References

- IEC 61400-12-1, *Power performance measurements of electricity producing wind turbines*: the air-density equation and the wind-speed density normalization.
- Wieringa (1992) and the Davenport classification for the roughness lengths.
- The 1/7 power law and the neutral logarithmic profile as given in any boundary-layer meteorology text (e.g. Stull, *An Introduction to Boundary Layer Meteorology*).

## License

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
- MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.

### Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted for inclusion in the work by you, as defined in the Apache-2.0 license, shall be dual licensed as above, without any additional terms or conditions.
