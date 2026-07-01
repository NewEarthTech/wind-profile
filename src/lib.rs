//! # wind-profile
//!
//! Vertical wind-profile extrapolation for wind energy: take the wind a weather model
//! gives you at its fixed output heights and move it to the height a turbine actually
//! harvests at, then correct for what the air weighs when it gets there.
//!
//! ## Why hub height is a first-class problem
//!
//! Numerical weather models publish wind at a handful of fixed heights: 10 m always,
//! and depending on the model an 80 m or 100 m level (GFS carries 80 m, the ECMWF open
//! data carries 100 m). Modern turbine hubs sit at 80-160 m and climbing. Turbine power
//! goes with the cube of wind speed below rated, so a relative error in extrapolated
//! speed shows up three times as large in energy: 5% short on wind is roughly 15% short
//! on generation. The vertical profile is not a rounding step, it is where the money is.
//!
//! ## Which law, when
//!
//! - **Power law**, [`wind_speed_m_s_power_law`]: the wind-energy engineering standard,
//!   `u(z) = u_ref * (z / z_ref)^alpha`. Use [`shear_exponent_from_two_levels`] to derive
//!   `alpha` from two model heights when the model provides them; that locality beats any
//!   tabulated exponent. With only one height, [`SHEAR_EXPONENT_ONE_SEVENTH`] is the
//!   classic neutral-onshore default (offshore shear is typically nearer 0.10).
//! - **Log law**, [`wind_speed_m_s_log_law`]: the surface-layer similarity form,
//!   `u(z) = u_ref * ln(z / z0) / ln(z_ref / z0)`. Use it when you know the terrain and
//!   can pick a roughness length; the `Z0_*` constants carry the Davenport-Wieringa
//!   classification from open sea to city centre.
//!
//! Both are neutral-stability forms, the honest baseline when all you have is model
//! output. Diabatic (Monin-Obukhov) corrections need surface-flux inputs the target user
//! rarely has and are deliberately out of scope for now.
//!
//! ```
//! use wind_profile as wp;
//!
//! // A GFS-style pair: 6.2 m/s at 10 m, 8.9 m/s at 80 m. Derive the local shear...
//! let alpha = wp::shear_exponent_from_two_levels(6.2, 10.0, 8.9, 80.0);
//! assert!((alpha - 0.1738).abs() < 1e-3);
//!
//! // ...and carry the 80 m wind up to a 120 m hub.
//! let u_hub = wp::wind_speed_m_s_power_law(8.9, 80.0, 120.0, alpha);
//! assert!((u_hub - 9.55).abs() < 0.01);
//! ```
//!
//! ## Air density, the other half of the power calculation
//!
//! Power also scales linearly with air density, and density at a warm low-pressure site
//! is easily 10% off the 1.225 kg/m^3 reference that power curves are quoted at.
//! [`air_density_kg_m3`] is the IEC 61400-12-1 equation (humidity-corrected, via the
//! standard's own vapor-pressure fit), and [`density_corrected_wind_speed_m_s`] is the
//! standard's `(rho / rho_ref)^(1/3)` wind-speed normalization that makes a measured
//! speed comparable against a reference-density power curve.
//!
//! ```
//! use wind_profile as wp;
//!
//! // Dry air at ISA sea level is the IEC reference density exactly.
//! let rho = wp::air_density_kg_m3(288.15, 101_325.0, 0.0);
//! assert!((rho - wp::REFERENCE_AIR_DENSITY_KG_M3).abs() < 1e-3);
//! ```
//!
//! Units are explicit in every function name: `_m_s` metres per second, `_m` metres,
//! `_k` kelvin, `_pa` pascals, `_kg_m3` kilograms per cubic metre. Inputs are physical
//! quantities (heights above ground positive and above the roughness length, speeds
//! positive where a ratio is taken); outside that domain the arithmetic yields NaN or
//! infinity rather than a guess. No dependencies, no `unsafe`.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

/// Specific gas constant of dry air (J/(kg*K)), the value fixed by IEC 61400-12-1.
pub const R_DRY_AIR: f64 = 287.05;
/// Specific gas constant of water vapor (J/(kg*K)), per IEC 61400-12-1.
pub const R_WATER_VAPOR: f64 = 461.5;
/// The reference air density power curves are quoted at (kg/m^3): ISA sea level.
pub const REFERENCE_AIR_DENSITY_KG_M3: f64 = 1.225;
/// The classic neutral-stability onshore shear exponent, 1/7. A default for when the
/// model gives only one wind height; prefer [`shear_exponent_from_two_levels`] whenever
/// two heights exist. Offshore shear is typically nearer 0.10.
pub const SHEAR_EXPONENT_ONE_SEVENTH: f64 = 1.0 / 7.0;

/// Roughness length (m), Davenport-Wieringa class 1 "sea": open water, tidal flat.
pub const Z0_SEA_M: f64 = 0.0002;
/// Roughness length (m), Davenport-Wieringa class 2 "smooth": featureless land, ice.
pub const Z0_SMOOTH_M: f64 = 0.005;
/// Roughness length (m), Davenport-Wieringa class 3 "open": flat grassland, few obstacles.
pub const Z0_OPEN_M: f64 = 0.03;
/// Roughness length (m), Davenport-Wieringa class 4 "roughly open": low crops, scattered obstacles.
pub const Z0_ROUGHLY_OPEN_M: f64 = 0.10;
/// Roughness length (m), Davenport-Wieringa class 5 "rough": high crops, obstacle rows.
pub const Z0_ROUGH_M: f64 = 0.25;
/// Roughness length (m), Davenport-Wieringa class 6 "very rough": orchards, bushland.
pub const Z0_VERY_ROUGH_M: f64 = 0.5;
/// Roughness length (m), Davenport-Wieringa class 7 "closed": forest, suburb.
pub const Z0_CLOSED_M: f64 = 1.0;
/// Roughness length (m), Davenport-Wieringa class 8 "chaotic": city centre, high-rise.
pub const Z0_CHAOTIC_M: f64 = 2.0;

/// The shear exponent `alpha` implied by wind speeds observed at two heights:
/// `alpha = ln(u_hi / u_lo) / ln(z_hi / z_lo)`.
///
/// This is the local, in-the-moment exponent, and with modern model output (10 m plus
/// 80 m or 100 m) it beats any climatological default: it carries the actual stability
/// and terrain of the hour. Feed it straight into [`wind_speed_m_s_power_law`].
///
/// Both speeds and both heights must be positive, with the heights distinct; the
/// exponent is meaningful (and the profile monotone) when both speeds sit on the same
/// side of the profile, the overwhelmingly common case.
pub fn shear_exponent_from_two_levels(
    u_lo_m_s: f64,
    z_lo_m: f64,
    u_hi_m_s: f64,
    z_hi_m: f64,
) -> f64 {
    (u_hi_m_s / u_lo_m_s).ln() / (z_hi_m / z_lo_m).ln()
}

/// Wind speed at height `z_m` by the power law: `u_ref * (z / z_ref)^alpha`.
///
/// The wind-energy engineering standard for vertical extrapolation. `shear_exponent`
/// comes from [`shear_exponent_from_two_levels`] when the model provides two heights,
/// or [`SHEAR_EXPONENT_ONE_SEVENTH`] as the neutral onshore default. Heights are above
/// ground and must be positive; at `z_m == z_ref_m` the input speed comes back exactly.
pub fn wind_speed_m_s_power_law(
    u_ref_m_s: f64,
    z_ref_m: f64,
    z_m: f64,
    shear_exponent: f64,
) -> f64 {
    u_ref_m_s * (z_m / z_ref_m).powf(shear_exponent)
}

/// Wind speed at height `z_m` by the neutral log law:
/// `u_ref * ln(z / z0) / ln(z_ref / z0)`.
///
/// The surface-layer similarity profile over terrain of roughness length
/// `roughness_length_m` (pick a `Z0_*` constant, or supply a site value). Valid for
/// heights well above the roughness length, under near-neutral stratification. Both
/// heights must exceed `roughness_length_m`, which must be positive; at
/// `z_m == z_ref_m` the input speed comes back exactly.
pub fn wind_speed_m_s_log_law(
    u_ref_m_s: f64,
    z_ref_m: f64,
    z_m: f64,
    roughness_length_m: f64,
) -> f64 {
    u_ref_m_s * (z_m / roughness_length_m).ln() / (z_ref_m / roughness_length_m).ln()
}

/// Air density (kg/m^3) from temperature (K), pressure (Pa), and relative humidity
/// (0 to 1), by the IEC 61400-12-1 equation.
///
/// `rho = (1/T) * (B/R0 - phi * Pw * (1/R0 - 1/Rw))`, with `Pw` the standard's own
/// vapor-pressure fit `0.0000205 * exp(0.0631846 * T)` (Pa). Humid air is lighter than
/// dry air at the same state, so the correction always lowers density. With
/// `relative_humidity = 0` this reduces to the dry ideal-gas law.
///
/// Feed the temperature and pressure at hub height. The fit targets ordinary surface
/// meteorology (roughly 258-313 K); it is not a general psychrometric model.
pub fn air_density_kg_m3(temperature_k: f64, pressure_pa: f64, relative_humidity: f64) -> f64 {
    let vapor_pressure_pa = 0.000_020_5 * (0.063_184_6 * temperature_k).exp();
    (pressure_pa / R_DRY_AIR
        - relative_humidity * vapor_pressure_pa * (1.0 / R_DRY_AIR - 1.0 / R_WATER_VAPOR))
        / temperature_k
}

/// Wind speed normalized to a reference air density, per IEC 61400-12-1:
/// `u * (rho / rho_ref)^(1/3)`.
///
/// Turbine power scales with `rho * u^3`, so a measured speed in air of density
/// `rho_kg_m3` delivers the power that a speed of this value would deliver at
/// `rho_ref_kg_m3`. Normalize site wind with this before reading a power curve quoted
/// at [`REFERENCE_AIR_DENSITY_KG_M3`]. Thin air (lower density) corrects the speed
/// downward.
pub fn density_corrected_wind_speed_m_s(u_m_s: f64, rho_kg_m3: f64, rho_ref_kg_m3: f64) -> f64 {
    u_m_s * (rho_kg_m3 / rho_ref_kg_m3).powf(1.0 / 3.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) {
        assert!((a - b).abs() < tol, "{a} vs {b} (tol {tol})");
    }

    #[test]
    fn shear_exponent_recovers_the_constructed_profile() {
        let (u_lo, z_lo, z_hi): (f64, f64, f64) = (5.0, 10.0, 120.0);
        let u_hi = u_lo * (z_hi / z_lo).powf(0.2);
        approx(
            shear_exponent_from_two_levels(u_lo, z_lo, u_hi, z_hi),
            0.2,
            1e-12,
        );
    }

    #[test]
    fn one_seventh_law_matches_the_textbook_number() {
        // 8 m/s at 10 m under the 1/7 law is 8 * 8^(1/7) = 10.767 m/s at 80 m.
        let u = wind_speed_m_s_power_law(8.0, 10.0, 80.0, SHEAR_EXPONENT_ONE_SEVENTH);
        approx(u, 10.767, 0.001);
    }

    #[test]
    fn log_law_matches_the_hand_computation() {
        // 8 m/s at 10 m over open grassland to 100 m: 8 * ln(100/0.03) / ln(10/0.03).
        let u = wind_speed_m_s_log_law(8.0, 10.0, 100.0, Z0_OPEN_M);
        approx(u, 11.171, 0.001);
    }

    #[test]
    fn both_laws_return_the_input_at_the_reference_height() {
        approx(wind_speed_m_s_power_law(7.3, 80.0, 80.0, 0.19), 7.3, 1e-12);
        approx(
            wind_speed_m_s_log_law(7.3, 80.0, 80.0, Z0_ROUGH_M),
            7.3,
            1e-12,
        );
    }

    #[test]
    fn dry_air_at_isa_sea_level_is_the_iec_reference_density() {
        approx(
            air_density_kg_m3(288.15, 101_325.0, 0.0),
            REFERENCE_AIR_DENSITY_KG_M3,
            1e-3,
        );
    }

    #[test]
    fn humid_air_is_lighter_than_dry_air() {
        let dry = air_density_kg_m3(293.15, 101_325.0, 0.0);
        let humid = air_density_kg_m3(293.15, 101_325.0, 1.0);
        approx(dry, 1.2041, 5e-4); // ideal-gas value at 20 C
        approx(humid, 1.1939, 5e-4); // saturated, per the IEC vapor-pressure fit
        assert!(humid < dry);
    }

    #[test]
    fn density_correction_is_the_cube_root_ratio() {
        // At reference density the speed is untouched.
        approx(
            density_corrected_wind_speed_m_s(10.0, 1.225, 1.225),
            10.0,
            1e-12,
        );
        // 10% thin air corrects the speed down by the cube root: 0.9^(1/3) = 0.96549.
        let u = density_corrected_wind_speed_m_s(10.0, 0.9 * 1.225, 1.225);
        approx(u, 9.6549, 1e-3);
    }

    #[test]
    fn derived_shear_carries_a_model_pair_to_hub_height() {
        // The crate-docs example, end to end at full precision.
        let alpha = shear_exponent_from_two_levels(6.2, 10.0, 8.9, 80.0);
        let u_hub = wind_speed_m_s_power_law(8.9, 80.0, 120.0, alpha);
        approx(u_hub, 9.5499, 1e-3);
        // Extrapolating from the 10 m level with the same alpha lands on the same
        // profile, so the two answers must agree.
        approx(
            wind_speed_m_s_power_law(6.2, 10.0, 120.0, alpha),
            u_hub,
            1e-9,
        );
    }
}
