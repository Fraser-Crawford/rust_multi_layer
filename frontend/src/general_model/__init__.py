import numpy
import numpy as np
import pandas as pd
from rust_model import get_initial_state, y_prime, efflorescence,locking,volumes,crystal_y_prime,get_initial_crystal_state
from scipy.integrate import solve_ivp
from dataclasses import dataclass, field
from typing import Callable
from scipy.sparse import coo_matrix

def simulate_measurement(df, soluteD, soluteRI, solventD=1000.0,solventRI=1.33):
    angle = np.radians(45)
    n2 = (1 + 2 * df.solution_density * (((soluteRI ** 2 - 1) * df.mfs) / ((soluteRI ** 2 + 2) * soluteD) + (
            (1 - df.mfs) * (solventRI ** 2 - 1)) / (solventD * (solventRI ** 2 + 2)))) / (
                    1 - df.solution_density * (((soluteRI ** 2 - 1) * df.mfs) / ((soluteRI ** 2 + 2) * soluteD) + (
                    (1 - df.mfs) * (solventRI ** 2 - 1)) / (solventD * (solventRI ** 2 + 2))))
    n = np.sqrt(np.array(n2,dtype=float))
    water_factor = np.cos(angle/2)+1.33*np.sin(angle/2)/np.sqrt(1+1.33**2-2*1.33*np.cos(angle/2))
    new_factor = np.cos(angle/2)+n*np.sin(angle/2)/np.sqrt(1+n**2-2*n*np.cos(angle/2))
    return df.radius*water_factor/new_factor

class Timer:
    def __init__(self, thresholds):
        self.next_threshold = 0
        self.thresholds = thresholds
    def reset(self):
        self.next_threshold = 0
    def check_time(self, time):
        if self.next_threshold < len(self.thresholds):
            if time > self.thresholds[self.next_threshold]:
                print(f"Reached {self.thresholds[self.next_threshold]:.2f} s!")
                self.next_threshold += 1


@dataclass
class Droplet:
    temperature:float | Callable[[float],float] #Temeperature of the droplet core during evaporation
    relative_humidity:float | Callable[[float],float]
    layers: int #Layers are the number of shells used during simulation
    solution: str = field(default="water")
    suspension: str = field(default="silica") # The suspension definition used by the simulation
    particle_radius: float = field(default=200e-9)
    convective:bool = field(default=True) #Whether or not the system is allowed to ciculate
    timer: Timer = field(default=None) #Timer class for allowing print outs of the current time
    gravity: tuple = field(default=(0,0,-9.81))
    air_speed: tuple = field(default=(0.0,0.0,0.0))
    velocity: tuple = field(default=(0.0,0.0,0.0))
    stationary: bool = field(default=True)
    def starting_state(self,radius:float,solute_concentration:float,particle_concentration:float):
        if type(self.relative_humidity) is float or type(self.relative_humidity) is int:
            rh = self.relative_humidity
        else:
            rh = self.relative_humidity(0.0)
        if type(self.temperature) is float or type(self.temperature) is int:
            temperature = self.temperature
        else:
            temperature = self.temperature(0.0)
        if type(self.air_speed) is tuple:
            air_speed = self.air_speed
        else:
            air_speed = self.air_speed(0.0)
        return np.array(get_initial_state(self.solution,temperature,rh,air_speed, self.gravity,
                                 self.suspension,self.particle_radius,radius,solute_concentration,
                                 particle_concentration,self.layers,self.velocity,self.stationary))

    def update_state(self,time, state, verbose):
        self.timer.check_time(time)
        if type(self.relative_humidity) is float or type(self.relative_humidity) is int:
            rh = self.relative_humidity
        else:
            rh = self.relative_humidity(time)
        if type(self.temperature) is float or type(self.temperature) is int:
            temperature = self.temperature
        else:
            temperature = self.temperature(time)
        if type(self.air_speed) is tuple:
            air_speed = self.air_speed
        else:
            air_speed = self.air_speed(time)
        derivative = np.array(y_prime(state,self.solution,temperature,rh,air_speed, self.gravity,
                       self.suspension,self.particle_radius,self.layers,self.convective,self.stationary))
        if verbose:
            print(f"STATE: {state}")
            print(f"DERIVATIVE: {derivative}")
        return derivative

    def efflorescence(self,time,state):
        if type(self.relative_humidity) is float or type(self.relative_humidity) is int:
            rh = self.relative_humidity
        else:
            rh = self.relative_humidity(time)
        if type(self.temperature) is float or type(self.temperature) is int:
            temperature = self.temperature
        else:
            temperature = self.temperature(time)
        if type(self.air_speed) is tuple:
            air_speed = self.air_speed
        else:
            air_speed = self.air_speed(time)
        activity = efflorescence(state,self.solution,temperature,rh,air_speed, self.gravity,
                       self.suspension,self.particle_radius,self.layers)
        return activity

    def locking(self,time,state,locking_threshold):
        if type(self.relative_humidity) is float or type(self.relative_humidity) is int:
            rh = self.relative_humidity
        else:
            rh = self.relative_humidity(time)
        if type(self.temperature) is float or type(self.temperature) is int:
            temperature = self.temperature
        else:
            temperature = self.temperature(time)
        if type(self.air_speed) is tuple:
            air_speed = self.air_speed
        else:
            air_speed = self.air_speed(time)
        value= locking(state,self.solution,temperature,rh,air_speed, self.gravity,
                       self.suspension,self.particle_radius,self.layers,locking_threshold)
        return value

    def equilibrate(self,time,state,threshold,verbose):
        dmdt = np.abs(self.update_state(time,state,verbose)[self.layers-1])*np.exp(state[self.layers-1])
        return dmdt - threshold

    def integrate(self,time:float,radius:float,solute_concentration=0.0,particle_concentration=0.0,
                  terminate_on_equilibration=False, equ_threshold=1e-4,
                  terminate_on_efflorescence=False, eff_threshold=0.5,
                  terminate_on_locking=False, locking_threshold=400e-9, timer=None, verbose=False,dense=False,rtol=1e-4):
        if timer is None:
            self.timer = Timer(np.linspace(0.0, time, 10))
        else:
            self.timer = timer
        x0 = self.starting_state(radius,solute_concentration,particle_concentration)
        events = []

        if terminate_on_equilibration:
            dmdt = np.abs(self.update_state(0.0,x0,verbose)[self.layers-1])*np.exp(x0[self.layers-1])
            equilibrated = lambda time, x:  self.equilibrate(time, x, equ_threshold*dmdt,verbose)
            equilibrated.terminal = True
            events += [equilibrated]

        if terminate_on_efflorescence:
            efflorescing = lambda time, x: self.efflorescence(time,x) - eff_threshold
            efflorescing.terminal = True
            events += [efflorescing]

        if terminate_on_locking:
            shell_formation = lambda time, x: self.locking(time,x,locking_threshold)
            shell_formation.terminal = True
            events += [shell_formation]

        dxdt = lambda time, x: self.update_state(time,x,verbose)
        trajectory = solve_ivp(dxdt, (0,time), x0, atol=1e-6,rtol=rtol, events=events, method="Radau",dense_output=dense)
        return trajectory

    def solve_at_time(self,time,trajectory):
        if trajectory.sol is not None:
            state = trajectory.sol(time)
            labels = ["radius","solution_density", "surface_temperature", "solvent_mass", "layer_mfs", "temperatures",
                      "mfs", "layer_positions", "layer_solute_concentrations",
                      "wet_layer_volumes", "solute_masses", "true_boundaries", "particle_masses",
                      "layer_particle_concentrations", "particle_volume_fraction", "solvent_masses",
                      "layer_solvent_concentrations"]
            variables = {key:  np.empty(1, dtype=object) for key in labels}
            earlier_droplet = DataDroplet(state, self.solution, self.suspension, self.particle_radius,
                                          self.layers)
            earlier_state = earlier_droplet.complete_state()
            for label, value in earlier_state.items():
                variables[label][0] = value
            variables['time'] = time
            return pd.DataFrame(variables)
        else:
            print("Trajectory was not densely calculated!")
            return None

    def complete_trajectory(self, trajectory):
        """Get the trajectory of all variables (including dependent ones) from a simulation (i.e.
        the output of UniformDroplet.integrate).

        Args:
            trajectory: the output of UniformDroplet.integrate, which gives the trajectory of independent
                        variables only.
        Returns:
            A pandas dataframe detailing the complete droplet history.
        """
        labels = ["radius","solution_density","surface_temperature","solvent_mass","layer_mfs","temperatures",
                  "mfs","layer_positions","layer_solute_concentrations",
                  "wet_layer_volumes","solute_masses","true_boundaries","particle_masses",
                  "layer_particle_concentrations","particle_volume_fraction","solvent_masses","layer_solvent_concentrations","position","velocity"]
        variables = {key: np.empty(trajectory.t.size, dtype=object) for key in labels}
        for i, state in enumerate(trajectory.y.T):
            earlier_droplet = DataDroplet(state, self.solution, self.suspension, self.particle_radius, self.layers)
            earlier_state = earlier_droplet.complete_state()
            for label, value in earlier_state.items():
                variables[label][i] = value

        variables['time'] = trajectory.t
        return pd.DataFrame(variables)

class DataDroplet:
    def __init__(self, state, solution, suspension, particle_radius, layers, particle_density=2200):
        state = np.array(state)
        self.solvent_masses = np.exp(state[:layers])
        self.solvent_mass = np.sum(self.solvent_masses)
        self.temperatures = state[layers:2*layers]
        self.solute_masses = np.exp(state[2 * layers:3 * layers])
        self.solute_mass = np.sum(self.solute_masses)
        self.layer_particle_mass = np.exp(state[3*layers:4*layers])

        self.velocity = np.array(state[4 * layers:3 + 4 * layers])
        self.position = np.array(state[3 + 4 * layers:6 + 4 * layers])

        self.particle_mass = np.sum(self.layer_particle_mass)

        radius = 0
        self.layer_positions = np.zeros(layers)
        self.layer_volumes = volumes(state, solution, suspension, particle_radius, layers)
        self.volume = np.sum(self.layer_volumes)
        for layer in range(layers):
            self.layer_positions[layer] = np.cbrt(self.layer_volumes[layer]*3/(4*np.pi)+radius**3)
            radius = self.layer_positions[layer]

        self.mfs = self.solute_mass/(self.solvent_mass+self.solute_mass)
        self.wet_layer_volumes = self.layer_volumes - self.layer_particle_mass / particle_density
        self.densities = (self.solute_masses+self.solvent_masses)/self.wet_layer_volumes
        self.radius = np.cbrt(self.volume*3/(4*np.pi))
        self.true_boundaries = np.concatenate(([0], self.layer_positions))
        self.layer_particle_concentration = self.layer_particle_mass/self.layer_volumes

        self.layer_concentrations = self.solute_masses /self.wet_layer_volumes
        self.layer_solvent_concentrations = self.solvent_masses /self.wet_layer_volumes
        self.particle_volume_fraction = self.particle_mass/(self.volume * particle_density)
        self.layer_mfs = self.solute_masses / (self.solvent_masses + self.solute_masses)
        self.density = (self.solvent_mass + self.solute_mass)/self.volume

    def complete_state(self):
        return dict(radius=self.radius,solution_density=self.density,surface_temperature=self.temperatures[-1], temperatures=self.temperatures,
                    solvent_mass=self.solvent_mass, mfs=self.mfs, layer_mfs = self.layer_mfs,
                    layer_positions=self.layer_positions, layer_solute_concentrations=self.layer_concentrations,
                    wet_layer_volumes=self.wet_layer_volumes, solute_masses=self.solute_masses, solvent_masses=self.solvent_masses, particle_masses=self.layer_particle_mass,
                    true_boundaries=self.true_boundaries, layer_particle_concentrations=self.layer_particle_concentration, particle_volume_fraction=self.particle_volume_fraction, layer_solvent_concentrations=self.layer_solvent_concentrations
                    , velocity=self.velocity,position = self.position)

@dataclass(kw_only=True)
class CrystalDroplet(Droplet):
    saturation: float
    growth_rate: float
    enthalpy: float
    crystal_density: float
    def starting_state(self, radius: float, solute_concentration: float, particle_concentration: float):
        if type(self.relative_humidity) is float or type(self.relative_humidity) is int:
            rh = self.relative_humidity
        else:
            rh = self.relative_humidity(0.0)
        if type(self.temperature) is float or type(self.temperature) is int:
            temperature = self.temperature
        else:
            temperature = self.temperature(0.0)
        if type(self.air_speed) is tuple:
            air_speed = self.air_speed
        else:
            air_speed = self.air_speed(0.0)
        return np.array(get_initial_crystal_state(solute_concentration,particle_concentration,radius,temperature,rh,air_speed, self.gravity,self.solution,self.suspension,self.particle_radius,self.layers,self.velocity,self.stationary))
    def update_state(self, time, state, verbose):
        self.timer.check_time(time)
        if type(self.relative_humidity) is float or type(self.relative_humidity) is int:
            rh = self.relative_humidity
        else:
            rh = self.relative_humidity(time)
        if type(self.temperature) is float or type(self.temperature) is int:
            temperature = self.temperature
        else:
            temperature = self.temperature(time)
        if type(self.air_speed) is tuple:
            air_speed = self.air_speed
        else:
            air_speed = self.air_speed(time)
        derivative = np.array(crystal_y_prime(state, self.solution,temperature,rh,air_speed, self.gravity,
                                      self.suspension, self.particle_radius, self.layers, self.convective, self.saturation, self.growth_rate, self.enthalpy, self.crystal_density,self.stationary))
        if verbose:
            print(f"STATE: {state}")
            print(f"DERIVATIVE: {derivative}")
        return derivative
    def complete_trajectory(self, trajectory):
        """Get the trajectory of all variables (including dependent ones) from a simulation (i.e.
        the output of UniformDroplet.integrate).

        Args:
            trajectory: the output of UniformDroplet.integrate, which gives the trajectory of independent
                        variables only.
        Returns:
            A pandas dataframe detailing the complete droplet history.
        """
        labels = ["radius","solution_density","surface_temperature","solvent_mass","layer_mfs","temperatures",
                  "mfs","layer_positions","layer_solute_concentrations",
                  "wet_layer_volumes","solute_masses","true_boundaries","particle_masses",
                  "layer_particle_concentrations","particle_volume_fraction","solvent_masses","layer_solvent_concentrations","crystal_mass","position","velocity"]
        variables = {key: np.empty(trajectory.t.size, dtype=object) for key in labels}
        for i, state in enumerate(trajectory.y.T):
            earlier_droplet = DataCrystalDroplet(state, self.solution, self.suspension, self.particle_radius, self.crystal_density,self.layers)
            earlier_state = earlier_droplet.complete_state()
            for label, value in earlier_state.items():
                variables[label][i] = value

        variables['time'] = trajectory.t
        return pd.DataFrame(variables)

class DataCrystalDroplet:
    def __init__(self, state, solution, suspension, particle_radius, crystal_density, layers, particle_density=2200):
        state = np.array(state)
        self.solvent_masses = np.exp(state[:layers])
        self.solvent_mass = np.sum(self.solvent_masses)
        self.temperatures = state[layers:2*layers]
        self.solute_masses = np.exp(state[2 * layers:3 * layers])
        self.solute_mass = np.sum(self.solute_masses)
        self.layer_particle_mass = np.exp(state[3*layers:4*layers])

        self.velocity = np.array(state[4*layers:3+4*layers])
        self.position = np.array(state[3+4 * layers:6 + 4 * layers])

        self.crystal_mass = np.exp(state[6 + 4 * layers])
        self.particle_mass = np.sum(self.layer_particle_mass)
        crystal_volume = self.crystal_mass/crystal_density
        radius = np.cbrt(3/(4*np.pi)*crystal_volume)
        self.layer_positions = np.zeros(layers)
        self.layer_volumes = volumes(state, solution, suspension, particle_radius, layers)
        self.volume = np.sum(self.layer_volumes)
        for layer in range(layers):
            self.layer_positions[layer] = np.cbrt(self.layer_volumes[layer]*3/(4*np.pi)+radius**3)
            radius = self.layer_positions[layer]

        self.mfs = self.solute_mass/(self.solvent_mass+self.solute_mass)
        self.wet_layer_volumes = self.layer_volumes - self.layer_particle_mass / particle_density
        self.densities = (self.solute_masses+self.solvent_masses)/self.wet_layer_volumes
        self.radius = np.cbrt(self.volume*3/(4*np.pi))
        self.true_boundaries = np.concatenate(([0], self.layer_positions))
        self.layer_particle_concentration = self.layer_particle_mass/self.layer_volumes

        self.layer_concentrations = self.solute_masses /self.wet_layer_volumes
        self.layer_solvent_concentrations = self.solvent_masses /self.wet_layer_volumes
        self.particle_volume_fraction = self.particle_mass/(self.volume * particle_density)
        self.layer_mfs = self.solute_masses / (self.solvent_masses + self.solute_masses)
        self.density = (self.solvent_mass + self.solute_mass)/self.volume

    def complete_state(self):
        return dict(radius=self.radius,solution_density=self.density,surface_temperature=self.temperatures[-1], temperatures=self.temperatures,
                    solvent_mass=self.solvent_mass, mfs=self.mfs, layer_mfs = self.layer_mfs,
                    layer_positions=self.layer_positions, layer_solute_concentrations=self.layer_concentrations,
                    wet_layer_volumes=self.wet_layer_volumes, solute_masses=self.solute_masses, solvent_masses=self.solvent_masses, particle_masses=self.layer_particle_mass,
                    true_boundaries=self.true_boundaries, layer_particle_concentrations=self.layer_particle_concentration, particle_volume_fraction=self.particle_volume_fraction,
                    layer_solvent_concentrations=self.layer_solvent_concentrations, crystal_mass = self.crystal_mass, velocity=self.velocity,position = self.position)