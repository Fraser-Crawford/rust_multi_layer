use std::f64::consts::PI;
use crate::binary_definitions::Binary;
use crate::environment::Environment;
use crate::layer::{Droplet,LIMIT_THRESHOLD};
use crate::suspension::Suspension;

pub struct CrystalDroplet{
    droplet: Droplet,
    log_crystal_mass: f64,
    crystal_radius: f64,
    n:usize,
}

impl CrystalDroplet {
    pub fn new(log_solvent_masses:Vec<f64>,temperatures:Vec<f64>,log_solute_masses:Vec<f64>,log_particle_masses:Vec<f64>,environment: Environment,solution:Binary,suspension:Suspension,log_crystal_mass:f64,n:usize, crystal_density:f64,position:&[f64],velocity:&[f64],stationary:bool)->Self{
        let crystal_mass = log_crystal_mass.exp();
        let crystal_volume = crystal_mass / crystal_density;
        let crystal_radius = (3.0/(4.0*PI)*crystal_volume).powf(1.0/3.0);
        CrystalDroplet{droplet:Droplet::new(log_solvent_masses,temperatures,log_solute_masses,log_particle_masses,environment,solution,suspension,n,crystal_radius,crystal_mass,velocity,position,stationary),log_crystal_mass,n,crystal_radius}
    }
    pub fn initial(solute_concentration:f64, particle_concentration:f64, radius:f64, environment:Environment, solution:Binary, suspension:Suspension, n:usize, velocity: &Vec<f64>, stationary:bool) ->Self{
        CrystalDroplet{droplet:Droplet::initial(solute_concentration,particle_concentration,radius,environment,solution,suspension,n,velocity,stationary),log_crystal_mass:-100.0,n,crystal_radius:0.0}
    }
    pub fn get_initial_state(&self)->Vec<f64>{
        let mut state = self.droplet.get_initial_state();
        state.push(self.log_crystal_mass);
        state
    }
    pub fn get_derivative(&mut self, convective:bool, saturation:f64, growth_rate:f64, enthalpy: f64) ->Vec<f64> {
        let inner_concentration = self.droplet.get_inner_concentration();
        let mut crystal_mass_derivative = growth_rate*4.0*PI*self.crystal_radius.powi(2)*(inner_concentration - saturation);
        if self.log_crystal_mass.exp() < LIMIT_THRESHOLD{
            crystal_mass_derivative = crystal_mass_derivative.max(0.0)
        }
        let heat = crystal_mass_derivative*enthalpy;
        let mut droplet_derivative = self.droplet.get_derivative(convective);
        droplet_derivative[self.n] += heat/self.droplet.get_inner_heat_capacity();
        droplet_derivative.push(crystal_mass_derivative/self.log_crystal_mass.exp());
        droplet_derivative[self.n*2] -= crystal_mass_derivative/self.droplet.get_inner_mass();
        droplet_derivative
    }
}