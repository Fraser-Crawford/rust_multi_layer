use std::f64::consts::PI;
use crate::binary_definitions::Binary;
use crate::suspension::Suspension;
use crate::environment::Environment;
use crate::fit::{asymmetric_gaussian, convection_coefficients};
use crate::general_droplet::{LIMIT_THRESHOLD, REDISTRIBUTION, SIGMA};

struct Layer{
    temperature: f64,
    solvent_mass: f64,
    particle_mass: f64,
    solute_mass: f64,
    volume: f64,
    solute_concentration:f64,
    particle_concentration: f64,
    solvent_concentration: f64,
    solute_diffusion:f64,
    solvent_diffusion: f64,
    particle_diffusion:f64,
    density:f64,
    solute_mass_fraction:f64,
    viscosity:f64,
    start:f64,
    heat_capacity:f64
}

impl Layer {
    pub fn initial(solute_concentration:f64,particle_concentration:f64,volume:f64,temperature:f64,solution:&Binary,suspension:&Suspension, start:&mut f64) -> Self {
        let particle_mass = particle_concentration*volume;
        let dry_volume = particle_mass/suspension.particle_density;
        let solute_mass_fraction = (solution.mfs)(solute_concentration);
        let solvent_mass_fraction = 1.0-solute_mass_fraction;
        let density = (solution.density)(solute_mass_fraction);
        let wet_volume = volume - dry_volume;
        let solvent_concentration = solvent_mass_fraction*density;
        let viscosity = (solution.viscosity)(solute_mass_fraction,temperature);
        let old_start = *start;
        *start = (3.0*volume/(4.0*PI)+start.powi(3)).powf(1.0/3.0);
        let solvent_mass = solvent_mass_fraction * density * wet_volume;
        Self {
            temperature,
            solvent_mass,
            particle_mass,
            solute_mass: solute_concentration * wet_volume,
            volume,
            solute_concentration,
            particle_concentration,
            solvent_concentration,
            solute_diffusion: (solution.non_volatile_diffusion)(solute_mass_fraction, temperature),
            solvent_diffusion: (solution.volatile_diffusion)(solvent_mass_fraction, temperature),
            particle_diffusion: suspension.diffusion(viscosity, temperature),
            density,
            solute_mass_fraction,
            viscosity: (solution.viscosity)(solute_mass_fraction, temperature),
            start: old_start,
            heat_capacity: (solution.volatile.specific_heat_capacity)(temperature)*solvent_mass + suspension.specific_heat_capacity*particle_mass
        }
    }
    pub fn new(log_solute_mass:f64, log_solvent_mass:f64, log_particle_mass:f64, temperature:f64, solution:&Binary, suspension:&Suspension, start: &mut f64) -> Self {
        let solute_mass = if log_solute_mass <= -80.0 || log_solute_mass > 0.0{
            0.0
        } else {
            log_solute_mass.exp()
        };
        let solvent_mass = if log_solvent_mass <= -80.0 || log_solvent_mass > 0.0{
            0.0
        } else {
            log_solvent_mass.exp()
        };
        let particle_mass = if log_particle_mass <= -80.0 || log_particle_mass > 0.0{
            0.0
        } else {
            log_particle_mass.exp()
        };
        let wet_mass = solute_mass+solvent_mass;
        let solute_mass_fraction = solute_mass/wet_mass;
        let density = (solution.density)(solute_mass_fraction);
        let wet_volume = wet_mass/density;
        let volume = wet_mass/density + particle_mass/suspension.particle_density;
        Layer::initial(solute_mass/wet_volume, particle_mass/volume, volume, temperature, solution, suspension, start)
    }
}

fn evaporation(radius:f64,surface_layer:&Layer,solution:&Binary,environment: &Environment)->(f64,f64){
    let d_eff = (solution.volatile.vapour_binary_diffusion_coefficient)(environment.temperature);
    let vapour_pressure = (solution.activity)(surface_layer.solute_mass_fraction,surface_layer.temperature)*(solution.volatile.equilibrium_vapour_pressure)(surface_layer.temperature);
    let vapour_ratio = (environment.pressure-vapour_pressure)/(environment.pressure-environment.vapour_pressure(&solution.volatile));
    let prandtl = environment.specific_heat_capacity*environment.dynamic_viscosity/(environment.thermal_conductivity)(environment.temperature);
    let reynolds = environment.density()*2.0*radius*environment.speed/environment.dynamic_viscosity;
    let knudsen = environment.mean_free_path()/radius;
    let beta = (1.0 + knudsen) / (1.0 + (4.0 / 3.0 * (1.0 + knudsen) + 0.377) * knudsen);
    let schmidt = environment.dynamic_viscosity/(environment.density()*(solution.volatile.vapour_binary_diffusion_coefficient)(environment.temperature));
    let sherwood = 1.0+0.3*reynolds.sqrt()*prandtl.powf(1.0/3.0);
    let nusselt = 1.0+0.3*reynolds.sqrt()*schmidt.powf(1.0/3.0);
    let dmdt = 4.0*PI*radius*environment.density()*(solution.volatile.molar_mass/environment.molar_mass)*d_eff*sherwood*vapour_ratio.ln()*beta;

    let conduction = nusselt * 4.0*PI*radius.powi(2) * (environment.thermal_conductivity)(environment.temperature) * (environment.temperature-surface_layer.temperature)/radius;

    let heat = (solution.volatile.specific_latent_heat_vaporisation)(surface_layer.temperature)*dmdt;

    let radiation = 4.0*PI*radius.powi(2)*SIGMA*(environment.temperature.powi(4)-surface_layer.temperature.powi(4));
    
    let dtdt = (conduction + heat - radiation)/surface_layer.heat_capacity;
    (dmdt,dtdt)
}
pub struct Droplet{
    radius:f64,
    layers:Vec<Layer>,
    solution:Binary,
    suspension:Suspension,
    n:usize,
    dmdt:f64,
    convection_volumes:Vec<f64>,
    redistribution_volumes:Vec<f64>,
    displacement_volumes:Vec<f64>,
    dtdt:f64,
}

impl Droplet {
    pub fn initial(solute_concentration:f64,particle_concentration:f64,radius:f64,environment:Environment,solution:Binary,suspension:Suspension,n:usize)->Self{
        let volume = 4.0/3.0*PI*radius.powi(3);
        let layer_volume = volume/n as f64;
        let mut start = 0.0;
        let coefficients = convection_coefficients(radius);
        let layers: Vec<Layer> = (0..n).map(|_i|Layer::initial(solute_concentration,particle_concentration,layer_volume,environment.temperature,&solution,&suspension,&mut start)).collect();
        let surface = layers.last().unwrap();
        let (dmdt,dtdt) = evaporation(radius,surface,&solution,&environment);
        let surface_speed = -dmdt/(4.0*PI*radius.powi(2)*layers.last().unwrap().density);
        Self{
            radius,
            convection_volumes: layers.iter().map(|layer|{
                0.5*asymmetric_gaussian((layer.start/radius).powi(2),coefficients)*environment.speed/0.02*(1.0+1e-3/1.81e-5)/(1.0+layer.viscosity/1.81e-5)
            }).collect(),
            redistribution_volumes: layers.windows(2).map(|pair|{
                (pair[1].volume-pair[0].volume)*REDISTRIBUTION*(n as f64).powi(2)
            }).collect(),
            displacement_volumes: layers.windows(2).map(|pair|{
                let fullness = (pair[1].particle_concentration-suspension.critical_volume_fraction*suspension.particle_density).max(0.0).powi(2)/
                    (suspension.maximum_volume_fraction*suspension.particle_density - suspension.critical_volume_fraction * suspension.particle_density);
                surface_speed*4.0*PI*pair[0].start.powi(3)/radius*fullness/suspension.particle_density
            }).collect(),
            layers,
            solution,
            suspension,
            n,
            dtdt,
            dmdt,
        }
    }
    pub fn new(log_solvent_masses:Vec<f64>,temperatures:Vec<f64>,log_solute_masses:Vec<f64>,log_particle_masses:Vec<f64>,environment: Environment,solution:Binary,suspension:Suspension,n:usize)->Self{
        let mut start = 0.0;
        let layers:Vec<Layer> = (0..n).map(|i|Layer::new(log_solute_masses[i],log_solvent_masses[i],log_particle_masses[i],temperatures[i],&solution,&suspension,&mut start)).collect();
        let volume = layers.iter().enumerate().fold(0.0,|acc,(index,layer)| {
            acc + layer.volume});
        
        let radius = (3.0*volume/(4.0*PI)).powf(1.0/3.0);
        let surface = layers.last().unwrap();
        let (dmdt,dtdt) = evaporation(radius,surface,&solution,&environment);
        let coefficients = convection_coefficients(radius);
        let surface_speed = -dmdt/(4.0*PI*radius.powi(2)*layers.last().unwrap().density);
        Self{
            radius,
            n,
            dmdt,
            convection_volumes: layers.iter().map(|layer|{
                0.5e-18*asymmetric_gaussian((layer.start/radius).powi(2),coefficients)*environment.speed/0.02*(1.0+1e-3/1.81e-5)/(1.0+layer.viscosity/1.81e-5)
            }).collect(),
            redistribution_volumes: layers.windows(2).map(|pair|{
                (pair[1].volume-pair[0].volume)*REDISTRIBUTION*(n as f64).powi(2)
            }).collect(),
            displacement_volumes: layers.windows(2).map(|pair|{
                let fullness = (pair[1].particle_concentration-suspension.critical_volume_fraction*suspension.particle_density).max(0.0)*pair[1].particle_concentration/
                    (suspension.maximum_volume_fraction*suspension.particle_density - suspension.critical_volume_fraction * suspension.particle_density);
                surface_speed*4.0*PI*pair[0].start.powi(3)/radius*fullness/suspension.particle_density
            }).collect(),
            layers,
            dtdt,
            solution,
            suspension,
        }
    }
    pub fn get_initial_state(&self)->Vec<f64>{
        let mut log_solvent_masses = Vec::new();
        let mut log_solute_masses = Vec::new();
        let mut log_particle_masses = Vec::new();
        let mut temperatures = Vec::new();
        for layer in self.layers.iter(){
            log_solvent_masses.push(layer.solvent_mass.ln().max(-100.0));
            log_solute_masses.push(layer.solute_mass.ln().max(-100.0));
            log_particle_masses.push(layer.particle_mass.ln().max(-100.0));
            temperatures.push(layer.temperature);
        }
        [log_solvent_masses,temperatures,log_solute_masses,log_particle_masses].concat()
    }
    pub fn surface_activity(&self)->f64{
        let surface = self.layers.last().unwrap();
        (self.solution.activity)(surface.solute_mass_fraction,surface.temperature)
    }
    pub fn locking(&self, critical_layer_thickness:f64)->f64{
        let radial_position = self.radius - critical_layer_thickness;
        if radial_position <= 0.0{
            -1.0
        } else {
            for layer in self.layers.iter().rev(){
                let discriminant = self.suspension.critical_volume_fraction - layer.particle_concentration/self.suspension.particle_density;
                if discriminant > 0.0{
                    return discriminant;
                } else {
                    if layer.start < radial_position {
                        return discriminant;
                    }
                }
            }
            1.0
        }
    }
    pub fn layer_volumes(&self)->Vec<f64>{
        self.layers.iter().map(|layer|{layer.volume}).collect()
    }
    pub fn solute_diffusion(&self, index:usize) ->f64{
        let layer0 = &self.layers[index];
        let layer1 = &self.layers[index+1];
        let gradient = (layer1.solute_concentration - layer0.solute_concentration)/(layer1.start-layer0.start);
        let diffusion = (layer0.solute_diffusion+layer1.solute_diffusion)/2.0;
        4.0*PI*diffusion*gradient*layer1.start.powi(2)
    }
    pub fn solute_convect(&self, index:usize) ->f64{
        self.convection_volumes[index]*(&self.layers[index+1].solute_concentration-&self.layers[index].solute_concentration)
    }
    pub fn solute_redistribution(&self, index:usize) ->f64{
        let dv = self.redistribution_volumes[index];
        let layer = if dv > 0.0 {
            &self.layers[index+1]
        } else {
            &self.layers[index]
        };
        dv*layer.solute_mass/layer.volume
    }
    pub fn solute_displace(&self, index:usize) ->f64{
        -self.layers[index].solute_concentration*self.displacement_volumes[index]
    }
    pub fn solvent_diffusion(&self, index:usize) ->f64{
        let layer0 = &self.layers[index];
        let layer1 = &self.layers[index+1];
        let gradient = (layer1.solvent_concentration - layer0.solvent_concentration)/(layer1.start-layer0.start);
        let diffusion = (layer0.solvent_diffusion+layer1.solvent_diffusion)/2.0;
        4.0*PI*diffusion*gradient*layer1.start.powi(2)
    }
    pub fn solvent_convect(&self, index:usize) ->f64{
        self.convection_volumes[index]*(&self.layers[index+1].solvent_concentration-&self.layers[index].solvent_concentration)
    }
    pub fn solvent_redistribution(&self, index:usize) ->f64{
        let dv = self.redistribution_volumes[index];
        let layer = if dv > 0.0 {
            &self.layers[index+1]
        } else {
            &self.layers[index]
        };
        dv*layer.solvent_mass/layer.volume
    }
    pub fn solvent_displace(&self, index:usize) ->f64{
        -self.layers[index].solvent_concentration*self.displacement_volumes[index]
    }
    pub fn particle_diffusion(&self, index:usize) ->f64{
        let layer0 = &self.layers[index];
        let layer1 = &self.layers[index+1];
        let gradient = (layer1.particle_concentration - layer0.particle_concentration)/(layer1.start-layer0.start);
        let diffusion = (layer0.particle_diffusion+layer1.particle_diffusion)/2.0;
        4.0*PI*diffusion*gradient*layer1.start.powi(2)
    }
    pub fn particle_convect(&self, index:usize) ->f64{
        self.convection_volumes[index]*(&self.layers[index+1].particle_concentration-&self.layers[index].particle_concentration)
    }
    pub fn particle_redistribution(&self, index:usize) ->f64{
        let dv = self.redistribution_volumes[index];
        let layer = if dv > 0.0 {
            &self.layers[index+1]
        } else {
            &self.layers[index]
        };
        dv*layer.particle_concentration
    }
    pub fn particle_displace(&self, index:usize) ->f64{
        self.layers[index].particle_concentration*self.displacement_volumes[index]
    }
    pub fn temperature_diffusion(&self, index:usize) ->f64{
        let layer0 = &self.layers[index];
        let layer1 = &self.layers[index+1];
        let kappa = (self.solution.volatile.thermal_conductivity)((layer1.temperature + layer0.temperature)/2.0)/(layer0.heat_capacity);
        kappa * 4.0*PI*layer1.start.powi(2) * (layer1.temperature - layer0.temperature) / (layer1.start-layer0.start)
    }
    pub fn get_derivative(&mut self, convective:bool)->Vec<f64>{
        let mut solute_derivative:Vec<f64> = vec![0.0;self.n];
        let mut solvent_derivative:Vec<f64> = vec![0.0;self.n];
        let mut particle_derivative:Vec<f64> = vec![0.0;self.n];
        let mut temperature_derivative:Vec<f64> = vec![0.0;self.n];
        if !convective{
            self.convection_volumes = vec![0.0;self.n]
        }
        for index in 0..self.n{
            if self.layers[index].solute_mass > 0.0{
                if index < self.n-1{
                    let derivative = self.solute_diffusion(index) + self.solute_convect(index) + self.solute_redistribution(index) + self.solute_displace(index);
                    solute_derivative[index] += derivative;
                    solute_derivative[index + 1] -= derivative;
                }
                solute_derivative[index] /= self.layers[index].solute_mass;
            }
            if self.layers[index].solvent_mass > 0.0 {
                if index < self.n-1{
                    let derivative = self.solvent_diffusion(index) + self.solvent_convect(index) + self.solvent_redistribution(index) + self.solvent_displace(index);
                    solvent_derivative[index] += derivative;
                    solvent_derivative[index + 1] -= derivative;
                } else {
                    solvent_derivative[index] += self.dmdt;
                    if self.layers[index].solvent_mass < LIMIT_THRESHOLD{
                        solvent_derivative[index] = solvent_derivative[index].max(0.0);
                    }
                }
                solvent_derivative[index] /= self.layers[index].solvent_mass;
            }
            if self.layers[index].particle_mass > 0.0 {
                if index < self.n-1{
                    let derivative = self.particle_diffusion(index) + self.particle_convect(index) + self.particle_redistribution(index) + self.particle_displace(index);
                    particle_derivative[index] += derivative;
                    particle_derivative[index + 1] -= derivative;
                }
                particle_derivative[index] /= self.layers[index].particle_mass;
            }
            if index < self.n-1{
                let derivative = self.temperature_diffusion(index);
                temperature_derivative[index] += derivative;
                temperature_derivative[index + 1] -= derivative;
            } else {
                temperature_derivative[index] += self.dtdt;
            }
        }

        let result = [solvent_derivative,temperature_derivative,solute_derivative,particle_derivative].concat();
        result
    }
}
