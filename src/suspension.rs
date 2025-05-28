use std::f64::consts::PI;
use crate::environment::K;

pub struct Suspension {
    pub(crate) specific_heat_capacity: f64,
    pub particle_radius: f64,
    pub(crate) particle_density: f64,
    pub(crate) critical_volume_fraction: f64,
    pub(crate) maximum_volume_fraction: f64,
}

impl Suspension {
    pub fn diffusion(&self, viscosity:f64, temperature:f64)->f64{
        K*temperature/(6.0*PI*viscosity*self.particle_radius)
    }
    pub fn get(name:String,radius:f64)->Suspension{
        match name.as_str(){
            "silica"=>silica(radius),
            bad_suspension => {panic!("{} IS NOT A KNOWN SUSPENSION",bad_suspension)},
        }
    }
}

pub fn silica(particle_radius:f64)->Suspension{
    Suspension{
        specific_heat_capacity: 703.0,
        particle_radius,
        particle_density: 2200.0,
        critical_volume_fraction: 0.56,
        maximum_volume_fraction: PI/(3.0*2.0f64.sqrt()),
    }
}