mod binary_definitions;
mod environment;
mod fit;
mod water;
mod suspension;
mod solvent;
mod layer;
mod crystal_drop;

use pyo3::prelude::*;
use crate::binary_definitions::Binary;
use crate::crystal_drop::CrystalDroplet;
use crate::environment::{atmosphere, Environment};
use crate::layer::Droplet;
use crate::suspension::Suspension;

fn state_to_droplet(state:Vec<f64>, solution_string: String, environment:(f64, f64, f64), suspension_string: String, suspension_radius:f64, layers:usize) ->Droplet{
    let solution = Binary::get(solution_string);
    let suspension = Suspension::get(suspension_string,suspension_radius);
    let environment = atmosphere(environment.0,environment.1,environment.2);
    let log_solvent_masses = &state[0..layers];
    let temperatures = &state[layers..2*layers];
    let log_solute_mass = &state[2*layers..3*layers];
    let log_particle_mass = &state[3*layers..4*layers];
    Droplet::new(Vec::from(log_solvent_masses), Vec::from(temperatures), Vec::from(log_solute_mass), Vec::from(log_particle_mass), environment, solution, suspension, layers,0.0)
}

#[pyfunction]
pub fn get_initial_state(solution_string: String, environment:(f64,f64,f64), suspension_string: String, suspension_radius:f64, radius:f64,solute_concentration:f64,particle_concentration:f64,layers:usize)->Vec<f64>{
    let solution = Binary::get(solution_string);
    let suspension = Suspension::get(suspension_string,suspension_radius);
    let environment = atmosphere(environment.0,environment.1,environment.2);
    let mut droplet = Droplet::initial(solute_concentration,particle_concentration,radius,environment,solution,suspension,layers);
    droplet.get_initial_state()
}

#[pyfunction]
pub fn y_prime(state:Vec<f64>, solution_string: String, environment:(f64,f64,f64), suspension_string: String, suspension_radius:f64,layers:usize,convective:bool)->Vec<f64>{
    state_to_droplet(state, solution_string, environment, suspension_string, suspension_radius,layers).get_derivative(convective)
}

#[pyfunction]
pub fn efflorescence(state:Vec<f64>, solution_string: String, environment:(f64,f64,f64), suspension_string: String, suspension_radius:f64,layers:usize)->f64{
    let droplet = state_to_droplet(state.clone(), solution_string, environment, suspension_string, suspension_radius,layers);
    droplet.surface_activity()
}

#[pyfunction]
pub fn locking(state:Vec<f64>, solution_string: String, environment:(f64,f64,f64), suspension_string: String, suspension_radius:f64,layers:usize,critical_layer_thickness:f64)->f64{
    let droplet = state_to_droplet(state, solution_string, environment, suspension_string, suspension_radius,layers);
    droplet.locking(critical_layer_thickness)
}
#[pyfunction]
pub fn volumes(state:Vec<f64>, solution_string: String, suspension_string: String, suspension_radius:f64, layers:usize)->Vec<f64>{
    let droplet = state_to_droplet(state,solution_string,(293.0,0.0,0.0),suspension_string,suspension_radius,layers);
    droplet.layer_volumes()
}

pub fn state_to_crystal_droplet(state:Vec<f64>, solution_string:String, environment:(f64,f64,f64), suspension_string:String, suspension_radius:f64,layers:usize,crystal_density:f64)->CrystalDroplet{
    let solution = Binary::get(solution_string);
    let suspension = Suspension::get(suspension_string,suspension_radius);
    let environment = atmosphere(environment.0,environment.1,environment.2);
    let log_solvent_masses = &state[0..layers];
    let temperatures = &state[layers..2*layers];
    let log_solute_masses = &state[2*layers..3*layers];
    let log_particle_masses = &state[3*layers..4*layers];
    let log_crystal_mass = state[4*layers];
    CrystalDroplet::new(Vec::from(log_solvent_masses), Vec::from(temperatures), Vec::from(log_solute_masses), Vec::from(log_particle_masses), environment, solution, suspension, log_crystal_mass, layers,crystal_density)
}
#[pyfunction]
pub fn crystal_y_prime(state:Vec<f64>, solution_string: String, environment:(f64,f64,f64), suspension_string: String, suspension_radius:f64,layers:usize,convective:bool,saturation:f64,growth_rate:f64,enthalpy:f64,crystal_density:f64)->Vec<f64>{
    state_to_crystal_droplet(state, solution_string, environment, suspension_string, suspension_radius,layers,crystal_density).get_derivative(convective,saturation, growth_rate, enthalpy)
}

#[pyfunction]
pub fn get_initial_crystal_state(solute_concentration:f64,particle_concentration:f64,radius:f64,environment:(f64,f64,f64),solution_string:String,suspension_string:String,suspension_radius:f64,n:usize)->Vec<f64>{
    let solution = Binary::get(solution_string);
    let suspension = Suspension::get(suspension_string,suspension_radius);
    let environment = atmosphere(environment.0,environment.1,environment.2);
    CrystalDroplet::initial(solute_concentration,particle_concentration,radius,environment,solution,suspension,n).get_initial_state()
}



#[pymodule]
fn rust_model(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(y_prime, m)?)?;
    m.add_function(wrap_pyfunction!(get_initial_state, m)?)?;
    m.add_function(wrap_pyfunction!(efflorescence, m)?)?;
    m.add_function(wrap_pyfunction!(locking, m)?)?;
    m.add_function(wrap_pyfunction!(volumes, m)?)?;

    m.add_function(wrap_pyfunction!(crystal_y_prime, m)?)?;
    m.add_function(wrap_pyfunction!(get_initial_crystal_state, m)?)?;
    Ok(())
}

