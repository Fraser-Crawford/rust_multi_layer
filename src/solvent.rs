pub struct Solvent{
    pub(crate) molar_mass: f64,
    pub(crate) density: fn(f64) -> f64,
    pub(crate) specific_heat_capacity: fn(f64) -> f64,
    pub(crate) specific_latent_heat_vaporisation: fn(f64) -> f64,
    pub(crate) equilibrium_vapour_pressure: fn(f64) -> f64,
    pub(crate) vapour_binary_diffusion_coefficient: fn(f64) -> f64,
    pub(crate) surface_tension: fn(f64) -> f64,
    pub(crate) refractive_index: f64,
    pub(crate) thermal_conductivity: fn(f64) -> f64,
    pub(crate) viscosity: fn(f64) -> f64,
}