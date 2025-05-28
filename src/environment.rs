use std::f64::consts::PI;
use crate::solvent::Solvent;
use crate::water::water;

pub const K: f64 = 1.380649e-23;
pub const GAS_CONSTANT: f64 = 8.314;

pub struct Environment {
    pub(crate) molar_mass:f64,
    pub(crate) pressure: f64,
    pub(crate) temperature: f64,
    pub relative_humidity: f64,
    pub(crate) thermal_conductivity: fn(f64) -> f64,
    pub dynamic_viscosity: f64,
    pub speed: f64,
    pub specific_heat_capacity:f64
}

impl Environment {
    pub fn density(&self) -> f64 {
        (1e-3*self.molar_mass) * self.pressure / (GAS_CONSTANT * self.temperature)
    }
    pub fn vapour_pressure(&self,solvent: &Solvent) -> f64 {
        self.relative_humidity * (solvent.equilibrium_vapour_pressure)(self.temperature)
    }
    pub fn mean_free_path(&self) -> f64 {
        self.dynamic_viscosity / self.density() * (PI * 1e-3*self.molar_mass / (2.0*GAS_CONSTANT * self.temperature)).sqrt()
    }
}

pub fn atmosphere(temperature:f64,relative_humidity:f64,speed:f64)->Environment {
    let vapour_pressure_water = relative_humidity * (water().equilibrium_vapour_pressure)(temperature);
    let mole_fraction_water = vapour_pressure_water / 101325.0;
    let molar_mass = (1.0-mole_fraction_water) * 28.9647 + mole_fraction_water * water().molar_mass;
    Environment{
        molar_mass,
        pressure: 101325.0,
        temperature,
        relative_humidity,
        specific_heat_capacity:0.7175,
        thermal_conductivity: |t| 7.29955694e-05*t + 4.47229506e-03,
        dynamic_viscosity: 18.13e-6,
        speed
    }
}