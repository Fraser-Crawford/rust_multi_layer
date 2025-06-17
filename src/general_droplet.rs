use std::f64::consts::PI;
use crate::environment::{Environment};
use crate::fit::{asymmetric_gaussian, convection_coefficients};
use crate::suspension::Suspension;
use crate::binary_definitions::Binary;
pub const SIGMA:f64 = 5.670374419e-8;
const REDISTRIBUTION:f64 = 1e1;
const LIMIT_THRESHOLD:f64 = 1e-22;
pub struct GeneralDroplet{
    //IMMUTABLE STATE
    pub solution: Binary,
    environment: Environment,
    suspension: Suspension,
    layers: usize,
    
    //MUTABLE STATE
    temperatures: Vec<f64>,
    solvent_masses: Vec<f64>,
    particle_masses:Vec<f64>,
    solute_masses:Vec<f64>,

    //DERIVED
    pub radius:f64,
    pub layer_volumes:Vec<f64>,
    boundaries: Vec<f64>,
    nusselt:f64,
    sherwood:f64,
    beta:f64,
    pub vapour_pressure:f64,
    solute_concentrations:Vec<f64>,
    pub particle_concentrations:Vec<f64>,
    solute_diffusion:Vec<f64>,
    viscosities:Vec<f64>,
    density:f64,
    specific_heat_capacity:f64,
    solvent_diffusion:Vec<f64>,
}

impl GeneralDroplet{
    pub fn new(solution: Binary,suspension: Suspension,environment: Environment,radius:f64,solute_concentration:f64,particle_concentration:f64,layers:usize)->Self{
        let mfs = (solution.mfs)(solute_concentration);
        let mut log_solute_masses = Vec::new();
        let mut log_particle_masses = Vec::new();
        let mut log_solvent_masses = Vec::new();
        let volume = 4.0/3.0*PI*radius.powi(3);
        let layer_volumes = vec![volume/layers as f64;layers];
        for layer_volume in layer_volumes{
            let dry_volume = layer_volume*particle_concentration/suspension.particle_density;
            let layer_solute_mass = (layer_volume-dry_volume)*solute_concentration;
            let layer_solvent_mass = if mfs == 0.0{
                (layer_volume-dry_volume)*(solution.volatile.density)(environment.temperature)
            } else {
                layer_solute_mass/mfs - layer_solute_mass
            };
            if solute_concentration <= 0.0{
                log_solute_masses.push(-100.0)
            } else {
                log_solute_masses.push(layer_solute_mass.ln());
            }
            if particle_concentration <= 0.0{
                log_particle_masses.push(-100.0)
            } else {
                log_particle_masses.push((layer_volume*particle_concentration).ln());
            }
            log_solvent_masses.push(layer_solvent_mass.ln())
        }
        let T = environment.temperature;
        let N = layers;
        Self::new_from_state(solution,environment,suspension,log_solvent_masses,vec![T;N],log_solute_masses,log_particle_masses)
    }
    pub fn get_layer_densities(&self)->Vec<f64>{
        (0..self.layers).map(|i|(self.solvent_masses[i]+self.particle_masses[i]+self.solute_masses[i])/self.layer_volumes[i]).collect()
    }
    pub fn get_positions(&self)->Vec<f64>{
        [&[0.0],&self.boundaries[..],&[self.radius]].concat()
    }
    pub fn redistribute(&self)->Vec<(f64,f64,f64)>{
        let mut result = vec![(0.0,0.0,0.0);self.layers];
        for i in 0..self.layers-1{
            let v0 = self.layer_volumes[i];
            let v1 = self.layer_volumes[i+1];
            let dv = (v1-v0)*REDISTRIBUTION*(self.layers as f64).powi(2);
            let (d_solute, d_particle, d_solvent) = if dv > 0.0{
                (dv*self.solute_masses[i+1]/(v1-self.particle_masses[i+1]/self.suspension.particle_density),
                 dv*self.particle_masses[i+1]/v1,
                 dv*self.solvent_masses[i+1]/(v1-self.particle_masses[i+1]/self.suspension.particle_density))
            } else {
                (dv*self.solute_masses[i]/(v0-self.particle_masses[i]/self.suspension.particle_density),
                 dv*self.particle_masses[i]/v0,
                 dv*self.solvent_masses[i]/(v0-self.particle_masses[i]/self.suspension.particle_density))
            };
            result[i].0 += d_solute;
            result[i].1 += d_particle;
            result[i].2 += d_solvent;
            result[i+1].0 -= d_solute;
            result[i+1].1 -= d_particle;
            result[i+1].2 -= d_solvent;
        }
        result
    }
    pub fn new_from_state(solution:Binary, environment: Environment, suspension:Suspension, log_solvent_masses:Vec<f64>, temperatures:Vec<f64>, log_solute_masses:Vec<f64>, log_particle_masses:Vec<f64>) ->Self{
        let layers = log_solute_masses.len();
        let solvent_masses: Vec<f64> = log_solvent_masses.iter().map(|m|m.exp()).collect();
        let solute_masses:Vec<f64> = log_solute_masses.iter().map(|m| if *m <= -100.0 {
            0.0
        } else {
            m.exp()
        }).collect();
        let particle_masses:Vec<f64> = log_particle_masses.iter().map(|m| if *m <= -100.0 || *m > 10.0{
            0.0
        } else {
            m.exp()
        }).collect();
        let layer_volumes: Vec<f64> = solvent_masses.iter().zip(solute_masses.iter()).zip(particle_masses.iter())
            .map(|((solvent_mass,solute_mass),particle_mass)|{
                let mfs = solute_mass/(solute_mass+solvent_mass);
                let density = (solution.density)(mfs);
                let volume = (solvent_mass + solute_mass) / density + particle_mass / suspension.particle_density;
                if volume < 0.0{
                    panic!("volume {} is less than zero", volume);
                }
                volume
            }).collect();
        let mut old_radius = 0.0f64;
        let boundaries = layer_volumes[..layer_volumes.len() - 1].iter().map(|v|{
            if *v < 0.0{
                panic!("volume {} is less than zero", v);
            }
            let r3 = 3.0*v/(4.0*PI) + old_radius.powi(3);
            old_radius = r3.powf(1.0/3.0);
            old_radius
        }).collect::<Vec<f64>>();
        let volume = layer_volumes.iter().sum::<f64>();
        let radius= (3.0f64*volume/(4.0f64*PI)).powf(1.0/3.0);
        let knudsen = environment.mean_free_path()/radius;
        let beta = (1.0 + knudsen) / (1.0 + (4.0 / 3.0 * (1.0 + knudsen) + 0.377) * knudsen);
        let reynolds = environment.density()*2.0*radius*environment.speed/environment.dynamic_viscosity;
        let prandtl = environment.specific_heat_capacity*environment.dynamic_viscosity/(environment.thermal_conductivity)(environment.temperature);
        let schmidt = environment.dynamic_viscosity/(environment.density()*(solution.volatile.vapour_binary_diffusion_coefficient)(environment.temperature));
        let sherwood = 1.0+0.3*reynolds.sqrt()*prandtl.powf(1.0/3.0);
        let nusselt = 1.0+0.3*reynolds.sqrt()*schmidt.powf(1.0/3.0);

        let mut solute_concentrations = vec![solute_masses[0]/layer_volumes[0]-particle_masses[0]/suspension.particle_density;layers+1];
        let mut particle_concentrations = vec![particle_masses[0]/layer_volumes[0];layers+1];

        for i in 1..layers+1{
            solute_concentrations[i] = solute_masses[i-1]/(layer_volumes[i-1] - particle_masses[i-1]/suspension.particle_density);
            particle_concentrations[i] = particle_masses[i-1]/layer_volumes[i-1];
        }
        let viscosities: Vec<f64> = solute_concentrations.iter().zip(&temperatures).zip(particle_concentrations.iter()).map(|((solute_concentration,T),particle_concentration)|{
            let volume_fraction = particle_concentration/suspension.particle_density;
            let mfs = (solution.mfs)(*solute_concentration);
            (solution.viscosity)(mfs,*T)*(1.0+2.5*volume_fraction+25.0/4.0*volume_fraction.powi(2))
        }).collect();
        let solute_diffusion: Vec<f64> = solute_concentrations.iter().zip(&temperatures).zip(particle_concentrations.iter()).map(|((solute_concentration,T),particle_concentration)|{
            let mfs = (solution.mfs)(*solute_concentration);
            (solution.non_volatile_diffusion)(mfs,*T)
        }).collect();
        let solvent_diffusion: Vec<f64> = solute_concentrations.iter().zip(&temperatures).zip(particle_concentrations.iter()).map(|((solute_concentration,T),particle_concentration)|{
            let mfw = 1.0-(solution.mfs)(*solute_concentration);
            (solution.volatile_diffusion)(mfw,*T)
        }).collect();
        let vapour_pressure = (solution.activity)((solution.mfs)(*solute_concentrations.last().unwrap()),*temperatures.last().unwrap())*(solution.volatile.equilibrium_vapour_pressure)(*temperatures.last().unwrap());
        let total_mass = solvent_masses.iter().sum::<f64>() + particle_masses.iter().sum::<f64>() + solute_masses.iter().sum::<f64>();
        let density = total_mass/volume;
        let specific_heat_capacity = ((solution.volatile.specific_heat_capacity)(*temperatures.last().unwrap())*(solute_masses.iter().sum::<f64>() + solvent_masses.iter().sum::<f64>())+(suspension.specific_heat_capacity*particle_masses.iter().sum::<f64>()))/total_mass;
        Self{solute_diffusion,solvent_diffusion,viscosities, particle_masses, vapour_pressure, solute_concentrations,
            environment,solution,suspension,
            solvent_masses, radius, nusselt, sherwood,
            temperatures, beta, particle_concentrations, solute_masses,density,specific_heat_capacity,layers,boundaries,layer_volumes}
    }
    pub fn get_state(&self)->Vec<f64>{
        let log_mass_solute: Vec<f64> = self.solute_masses.iter().map(|m|
            if *m <= 0.0{
                -100.0
            } else {
                m.ln()
            }
            ).collect();
        let log_mass_particles: Vec<f64> = self.particle_masses.iter().map(|m|
            if *m <= 0.0{
                -100.0
            } else {
                m.ln()
            }).collect();
        let log_mass_solvent:Vec<f64> = self.solvent_masses.iter().map(|m|
            if *m <= 0.0{
                -100.0
            } else {
                m.ln()
            }).collect();
        [&log_mass_solvent[..],
            &self.temperatures[..],
            &log_mass_solute[..],
            &log_mass_particles[..],].concat()
    }
    pub fn temperature_derivative(&self, mass_derivative:f64) -> Vec<f64> {
        let mut result = vec![0.0;self.layers];
        if self.layers>1{
            let boundaries = [&[0.0f64],&self.boundaries[..]].concat();
            let normal_boundaries = boundaries.iter().map(|r|r/self.radius).collect::<Vec<f64>>();
            for i in 0..self.layers-1{
                let t0 = self.temperatures[i];
                let t1 = self.temperatures[i+1];
                let r0 = normal_boundaries[i];
                let r1 = normal_boundaries[i+1];
                let kappa = (self.solution.volatile.thermal_conductivity)((t1 + t0)/2.0)/(self.density*self.specific_heat_capacity);
                let denominator = (r1.powi(3)-r0.powi(3))*(r1-r0)*self.radius.powi(2);
                let numerator = 3.0*kappa*(t1-t0)*r1.powi(2);
                result[i] += numerator/denominator;
                result[i+1] -= numerator/denominator;
            }
        }
        let all_positions = [&[0.0f64],&self.boundaries[..]].concat();
        let r1 = self.radius;
        let r0 = all_positions[self.layers-1];
        let volume = 4.0/3.0*PI*(r1.powi(3)-r0.powi(3));
        let mass = volume*self.density;
        let heat_capacity = mass*self.specific_heat_capacity;

        let conduction = 4.0*PI*r1.powi(2)*(self.environment.thermal_conductivity)(self.environment.temperature)*
            (self.environment.temperature-self.temperatures.last().unwrap())/r1*self.nusselt;
        let heat = (self.solution.volatile.specific_latent_heat_vaporisation)(*self.temperatures.last().unwrap())*mass_derivative;
        let radiation = 4.0*PI*r1.powi(2)*SIGMA*(self.environment.temperature.powi(4)-self.temperatures.last().unwrap().powi(4));
        result[self.layers-1] += (conduction+heat-radiation)/heat_capacity;
        result
    }
    pub fn solvent_mass_derivative(&self) -> f64 {
        let d_eff = (self.solution.volatile.vapour_binary_diffusion_coefficient)(self.environment.temperature);
        let vapour_ratio = (self.environment.pressure-self.vapour_pressure)/(self.environment.pressure-self.environment.vapour_pressure(&self.solution.volatile));
        4.0*PI*self.radius*self.environment.density()*(self.solution.volatile.molar_mass/self.environment.molar_mass)*d_eff*self.sherwood*vapour_ratio.ln()*self.beta
    }
    fn get_gradients(&self, normalised_boundaries:&[f64])->Vec<(f64,f64,f64)>{
        let solute_iter = self.solute_concentrations.iter().skip(2).zip(self.solute_concentrations.iter().skip(1));
        let particle_iter = self.particle_concentrations.iter().skip(2).zip(self.particle_concentrations.iter().skip(1));
        let boundary_iter = normalised_boundaries.iter().skip(2).zip(normalised_boundaries.iter().skip(1));
        solute_iter.zip(particle_iter).zip(boundary_iter).map(|(((s1,s0),(p1,p0)),(b1,b0))| {
            let mfs0 = (self.solution.mfs)(*s0);
            let mfs1 = (self.solution.mfs)(*s1);

            let solvent_concentration0 = (self.solution.density)(mfs0)*(1.0-mfs0);
            let solvent_concentration1 = (self.solution.density)(mfs1)*(1.0-mfs1);
            ((s1-s0)/(b1-b0),(p1-p0)/(b1-b0),(solvent_concentration1-solvent_concentration0)/(b1-b0))
        }).collect()
    }
    pub fn diffuse(&self)->Vec<(f64,f64,f64)>{
        if self.layers>1{
            let boundaries = [&[0.0f64],&self.boundaries[..],&[self.radius]].concat();
            let normal_boundaries = boundaries.iter().map(|r|r/self.radius).collect::<Vec<f64>>();
            let gradients = self.get_gradients(&normal_boundaries);
            let mut result = vec![(0.0,0.0,0.0);self.layers];
            (0..self.layers-1).for_each(|i|{
                let solute_diffusion = self.solute_diffusion[i];
                let particle_diffusion = self.suspension.diffusion(self.viscosities[i],self.temperatures[i]);
                let solvent_diffusion = self.solvent_diffusion[i];
                let solute_value = 4.0*PI*self.radius*solute_diffusion*gradients[i].0*normal_boundaries[i+1].powi(2);
                let particle_value = 4.0*PI*self.radius*particle_diffusion*gradients[i].1*normal_boundaries[i+1].powi(2);
                let solvent_value = 4.0*PI*self.radius*solvent_diffusion*gradients[i].2*normal_boundaries[i+1].powi(2);
                result[i].0 += solute_value;
                result[i+1].0 -= solute_value;
                result[i].1 += particle_value;
                result[i+1].1 -= particle_value;
                result[i].2 += solvent_value;
                result[i+1].2 -= solvent_value;
            });
            result
        } else {
            vec![(0.0,0.0,0.0)]
        }
    }
    pub fn convection(&self)->Vec<(f64,f64)>{
        let mut result = vec![(0.0,0.0);self.layers];
        let coefficients = convection_coefficients(self.radius);
        let layer_volumes = self.layer_volumes.clone();
        let solute_concentrations = layer_volumes.iter().zip(self.solute_masses.iter()).map(|(v,m)|{
            m/(v*1e18)
        }).collect::<Vec<f64>>();
        let particle_concentrations = layer_volumes.iter().zip(self.particle_masses.iter()).map(|(v,m)|{
            m/(v*1e18)
        }).collect::<Vec<f64>>();
        for (i,boundary) in self.boundaries.iter().enumerate(){
            let rate = 0.5*asymmetric_gaussian((boundary/self.radius).powi(2),coefficients)*self.environment.speed/0.02*(1.0+1e-3/1.81e-5)/(1.0+self.viscosities.last().unwrap()/1.81e-5);
            result[i].0 += rate*(solute_concentrations[i+1]-solute_concentrations[i]);
            result[i].1 += rate*(particle_concentrations[i+1]-particle_concentrations[i]);
            result[i+1].0 -= rate*(solute_concentrations[i+1]-solute_concentrations[i]);
            result[i+1].1 -= rate*(particle_concentrations[i+1]-particle_concentrations[i]);
        }
        result
    }
    pub fn displace(&self, dmdt:f64)->Vec<(f64,f64,f64)>{
        if self.layers>1{
            let mut result = vec![(0.0,0.0,0.0);self.layers];
            if dmdt > 0.0{
                result
            } else {
                let fullness:Vec<f64> = self.particle_concentrations.iter().skip(2).map(|c|{
                    (c-self.suspension.critical_volume_fraction*self.suspension.particle_density).max(0.0).powi(2)/
                        (self.suspension.maximum_volume_fraction*self.suspension.particle_density - self.suspension.critical_volume_fraction * self.suspension.particle_density)
                }).collect();
                let positions = &self.get_positions()[1..];
                let surface_speed = -dmdt/(4.0*PI*self.radius.powi(2)*self.density);
                (0..self.layers-1).for_each(|i|{
                    let volume = surface_speed*4.0*PI*positions[i].powi(3)/self.radius;
                    let particle_rate = volume*fullness[i];
                    let volume = particle_rate / self.suspension.particle_density;
                    let solute_rate = self.solute_concentrations[i]*volume;
                    let mfs = self.solute_masses[i]/(self.solute_masses[i]+self.solvent_masses[i]);
                    let solvent_rate = (self.solution.density)(mfs)*(1.0-mfs)*volume;
                    result[i].0 -= solute_rate;
                    result[i+1].0 += solute_rate;
                    result[i].2 -= solvent_rate;
                    result[i+1].2 += solvent_rate;
                    result[i].1 += particle_rate;
                    result[i+1].1 -= particle_rate;
                });
                result
            }
        } else {
            vec![(0.0,0.0,0.0)]
        }

    }
    pub fn dxdt(&self, convective:bool)->Vec<f64>{
        let mut solute_derivative = vec![0.0; self.layers];
        let mut particle_derivative = vec![0.0; self.layers];
        let mut solvent_derivative = vec![0.0; self.layers];
        let diffusion = self.diffuse();
        let convect = if convective{
            self.convection()
        } else {
            vec![(0.0,0.0);self.layers]
        };
        let mass_derivative = self.solvent_mass_derivative();
        let temperature_derivative = self.temperature_derivative(mass_derivative);
        let redistribution = self.redistribute();
        let displacement = self.displace(mass_derivative);
        for i in 0..self.layers{
            if self.solute_masses[i] > 0.0{
                solute_derivative[i] += (diffusion[i].0+redistribution[i].0+displacement[i].0+convect[i].0)/self.solute_masses[i];
            }
            if self.particle_masses[i] > 0.0{
                particle_derivative[i] += (diffusion[i].1+redistribution[i].1+displacement[i].1+convect[i].1)/self.particle_masses[i];
            }
            if self.solvent_masses[i] > 0.0{
                solvent_derivative[i] += (diffusion[i].2+redistribution[i].2+displacement[i].2)/self.solvent_masses[i];
            }
            if i == self.layers -1{
                let dmdt = mass_derivative/self.solvent_masses[self.layers-1];
                solvent_derivative[self.layers-1] += dmdt;
            }
            if self.solvent_masses[i] < LIMIT_THRESHOLD{
                solvent_derivative[i] = solvent_derivative[i].max(0.0);
            }
        }
        [solvent_derivative,
            temperature_derivative,
            solute_derivative,
            particle_derivative].concat()
    }
}

