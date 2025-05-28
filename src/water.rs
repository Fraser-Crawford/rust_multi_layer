use crate::fit::polynomial;
use crate::solvent::Solvent;

fn thermal_conductivity_water(t:f64)->f64{
    -0.432+5.725e-3*t-8.078e-6*t.powi(2)+1.861e-9*t.powi(3)
}
fn surface_tension(A:f64, B:f64, C:f64, D:f64, E:f64, T_crit:f64, T:f64) ->f64{
    let T_r = T / T_crit;
    let power = B + (C * T_r) + (D * T_r.powi(2)) + (E * T_r.powi(3));
    (A * (1.0 - T_r).powf(power)) / 1000.0
}
fn surface_tension_water(temperature:f64)->f64{
    surface_tension(134.15, 1.6146, -2.035, 1.5598, 0.0, 647.3, temperature)
}

pub fn water_viscosity(temperature:f64)->f64{
    let T_R = temperature/273.1;
    0.09607e-6*(2.9/T_R.powi(3)).exp()*density_water(temperature)
}
fn specific_heat_capacity_water(t:f64)->f64{
    polynomial(&[4.92521605e-10, -9.80939941e-07,  8.14780736e-04, -3.61290267e-01,
               9.02139176e+01, -1.20295166e+04,  6.73532071e+05,],t)
}

fn equilibrium_vapour_pressure_water(temperature:f64)->f64{
    let T_C = temperature - 273.15;
    1e3 * 0.61161 * ((18.678 - (T_C / 234.5)) * (T_C / (257.14 + T_C))).exp()
}

pub fn density_water(temperature:f64)->f64{
    let ref_T = 647.096;
    let ref_density = 322.0;
    let (b1, b2, b3, b4, b5, b6) = (1.99274064, 1.09965342, -0.510839303, -1.75493479, -45.5170352, -674694.45);

    let theta = temperature / ref_T;
    let tau = 1.0 - theta;
    ref_density * (1.0 + b1 * tau.powf(1.0/3.0) + b2 * tau.powf(2.0/3.0) + b3 * tau.powf(5.0/3.0) + b4 * tau.powf(16.0/3.0) + b5 * tau.powf(45.0/3.0) + b6 * tau.powf(110.0/3.0))
}

pub fn water()->Solvent{
    Solvent{
        molar_mass:18.02,
        density:density_water,
        specific_heat_capacity:specific_heat_capacity_water,
        specific_latent_heat_vaporisation:|T| 3.14566e6 - 2361.64 * T,
        equilibrium_vapour_pressure:equilibrium_vapour_pressure_water,
        vapour_binary_diffusion_coefficient: |t| 0.2190e-4*(t/273.15).powf(1.81),
        surface_tension: surface_tension_water,
        refractive_index: 1.335,
        thermal_conductivity: thermal_conductivity_water,
        viscosity: water_viscosity,
    }
}