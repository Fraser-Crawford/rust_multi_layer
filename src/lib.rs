mod binary_definitions;
mod environment;
mod fit;
mod general_droplet;
mod water;
mod suspension;
mod solvent;
mod layer;

use pyo3::prelude::*;
use crate::binary_definitions::Binary;
use crate::environment::atmosphere;
use crate::general_droplet::GeneralDroplet;
use crate::layer::Droplet;
use crate::suspension::Suspension;

fn state_to_general_droplet(state:Vec<f64>, solution_string: String, environment:(f64, f64, f64), suspension_string: String, suspension_radius:f64, layers:usize) ->GeneralDroplet{
    let solution = Binary::get(solution_string);
    let suspension = Suspension::get(suspension_string,suspension_radius);
    let environment = atmosphere(environment.0,environment.1,environment.2);
    let log_solvent_masses = &state[0..layers];
    let temperatures = &state[layers..2*layers];
    let log_solute_mass = &state[2*layers..3*layers];
    let log_particle_mass = &state[3*layers..4*layers];
    GeneralDroplet::new_from_state(solution, environment, suspension, Vec::from(log_solvent_masses),
                                   Vec::from(temperatures), Vec::from(log_solute_mass), Vec::from(log_particle_mass))
}

fn state_to_droplet(state:Vec<f64>, solution_string: String, environment:(f64, f64, f64), suspension_string: String, suspension_radius:f64, layers:usize) ->Droplet{
    let solution = Binary::get(solution_string);
    let suspension = Suspension::get(suspension_string,suspension_radius);
    let environment = atmosphere(environment.0,environment.1,environment.2);
    let log_solvent_masses = &state[0..layers];
    let temperatures = &state[layers..2*layers];
    let log_solute_mass = &state[2*layers..3*layers];
    let log_particle_mass = &state[3*layers..4*layers];
    Droplet::new(Vec::from(log_solvent_masses), Vec::from(temperatures), Vec::from(log_solute_mass), Vec::from(log_particle_mass), environment, solution, suspension, layers)
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
pub fn get_initial_general_state(solution_string: String, environment:(f64,f64,f64), suspension_string: String, suspension_radius:f64, radius:f64,solute_concentration:f64,particle_concentration:f64,layers:usize)->Vec<f64>{
    let solution = Binary::get(solution_string);
    let suspension = Suspension::get(suspension_string,suspension_radius);
    let environment = atmosphere(environment.0,environment.1,environment.2);
    let mut droplet = GeneralDroplet::new(solution,suspension,environment,radius,solute_concentration,particle_concentration,layers);
    droplet.get_state()
}

#[pyfunction]
pub fn y_prime(state:Vec<f64>, solution_string: String, environment:(f64,f64,f64), suspension_string: String, suspension_radius:f64,layers:usize,convective:bool)->Vec<f64>{
    state_to_droplet(state, solution_string, environment, suspension_string, suspension_radius,layers).get_derivative(convective)
}

#[pyfunction]
pub fn general_y_prime(state:Vec<f64>, solution_string: String, environment:(f64,f64,f64), suspension_string: String, suspension_radius:f64,layers:usize,convective:bool)->Vec<f64>{
    state_to_general_droplet(state, solution_string, environment, suspension_string, suspension_radius,layers).dxdt(convective)
}

#[pyfunction]
pub fn efflorescence(state:Vec<f64>, solution_string: String, environment:(f64,f64,f64), suspension_string: String, suspension_radius:f64,layers:usize)->f64{
    let droplet = state_to_droplet(state.clone(), solution_string, environment, suspension_string, suspension_radius,layers);
    droplet.surface_activity()
}

#[pyfunction]
pub fn general_efflorescence(state:Vec<f64>, solution_string: String, environment:(f64,f64,f64), suspension_string: String, suspension_radius:f64,layers:usize)->f64{
    let droplet = state_to_general_droplet(state.clone(), solution_string, environment, suspension_string, suspension_radius,layers);
    droplet.vapour_pressure/(droplet.solution.volatile.equilibrium_vapour_pressure)(state[layers])
}

#[pyfunction]
pub fn locking(state:Vec<f64>, solution_string: String, environment:(f64,f64,f64), suspension_string: String, suspension_radius:f64,layers:usize,critical_layer_thickness:f64)->f64{
    let droplet = state_to_droplet(state, solution_string, environment, suspension_string, suspension_radius,layers);
    droplet.locking(critical_layer_thickness)
}

#[pyfunction]
pub fn general_locking(state:Vec<f64>, solution_string: String, environment:(f64,f64,f64), suspension_string: String, suspension_radius:f64,layers:usize,critical_layer_thickness:f64)->f64{
    let suspension = Suspension::get(suspension_string.clone(),suspension_radius);
    let droplet = state_to_general_droplet(state, solution_string, environment, suspension_string, suspension_radius,layers);
    let positions = droplet.get_positions();
    let radial_position = droplet.radius-critical_layer_thickness;
    if radial_position <= 0.0{
        1.0
    } else {
        let index = positions.iter().enumerate().filter_map(
            |(index,&r)| if r >= radial_position {
                Some(index)
            } else {None}
        ).collect::<Vec<usize>>()[0];
        let m = (droplet.particle_concentrations[index] - droplet.particle_concentrations[index - 1]) / (
            positions[index] - positions[index - 1]);
        suspension.critical_volume_fraction - (droplet.particle_concentrations[index] + (radial_position - positions[index]) * m)/suspension.particle_density
    }
}
#[pyfunction]
pub fn volumes(state:Vec<f64>, solution_string: String, suspension_string: String, suspension_radius:f64, layers:usize)->Vec<f64>{
    let droplet = state_to_droplet(state,solution_string,(293.0,0.0,0.0),suspension_string,suspension_radius,layers);
    droplet.layer_volumes()
}
#[pyfunction]
pub fn layer_volumes(state:Vec<f64>, solution_string: String, suspension_string: String, suspension_radius:f64, layers:usize)->Vec<f64>{
    let droplet = state_to_general_droplet(state,solution_string,(293.0,0.0,0.0),suspension_string,suspension_radius,layers);
    droplet.layer_volumes
}

#[pymodule]
fn rust_model(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(y_prime, m)?)?;
    m.add_function(wrap_pyfunction!(get_initial_state, m)?)?;
    m.add_function(wrap_pyfunction!(efflorescence, m)?)?;
    m.add_function(wrap_pyfunction!(locking, m)?)?;
    m.add_function(wrap_pyfunction!(volumes, m)?)?;
    
    m.add_function(wrap_pyfunction!(general_y_prime, m)?)?;
    m.add_function(wrap_pyfunction!(get_initial_general_state, m)?)?;
    m.add_function(wrap_pyfunction!(general_efflorescence, m)?)?;
    m.add_function(wrap_pyfunction!(general_locking, m)?)?;
    m.add_function(wrap_pyfunction!(layer_volumes, m)?)?;
    Ok(())
}

