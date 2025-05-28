import numpy as np
import pandas as pd
from rust_model import get_initial_general_state, general_y_prime, general_efflorescence,general_locking,layer_volumes
from scipy.integrate import solve_ivp
from dataclasses import dataclass, field
from typing import Callable

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
class GeneralDroplet:
    temperature:float | Callable[[float],float]
    relative_humidity:float | Callable[[float],float]
    air_speed:float | Callable[[float],float]
    layers: int
    solution: str = field(default="water")
    suspension: str = field(default="silica")
    particle_radius: float = field(default=200e-9)
    timer: Timer = field(default=None)
    def starting_state(self,radius:float,solute_concentration:float,particle_concentration:float):
        if type(self.relative_humidity) is float or type(self.relative_humidity) is int:
            rh = self.relative_humidity
        else:
            rh = self.relative_humidity(0.0)
        if type(self.temperature) is float or type(self.temperature) is int:
            temperature = self.temperature
        else:
            temperature = self.temperature(0.0)
        if type(self.air_speed) is float or type(self.air_speed) is int:
            air_speed = self.air_speed
        else:
            air_speed = self.air_speed(0.0)
        return np.array(get_initial_general_state(self.solution,(temperature,rh,air_speed),
                                 self.suspension,self.particle_radius,radius,solute_concentration,
                                 particle_concentration,self.layers))

    def update_state(self,time, state):
        self.timer.check_time(time)
        if type(self.relative_humidity) is float or type(self.relative_humidity) is int:
            rh = self.relative_humidity
        else:
            rh = self.relative_humidity(time)
        if type(self.temperature) is float or type(self.temperature) is int:
            temperature = self.temperature
        else:
            temperature = self.temperature(time)
        if type(self.air_speed) is float or type(self.air_speed) is int:
            air_speed = self.air_speed
        else:
            air_speed = self.air_speed(time)
        derivative = np.array(general_y_prime(state,self.solution,(temperature,rh,air_speed),
                       self.suspension,self.particle_radius,self.layers))
        return derivative

    def efflorescence(self,state):
        return general_efflorescence(state,self.solution,(self.temperature,self.relative_humidity,self.air_speed),
                       self.suspension,self.particle_radius,self.layers)

    def locking(self,state,locking_threshold):
        value= general_locking(state,self.solution,(self.temperature,self.relative_humidity,self.air_speed),
                       self.suspension,self.particle_radius,self.layers,locking_threshold)
        return value

    def equilibrate(self,time,state,threshold):
        dmdt = np.abs(self.update_state(time,state)[self.layers-1])*np.exp(state[self.layers-1])
        return dmdt - threshold

    def integrate(self,time:float,radius:float,solute_concentration=0.0,particle_concentration=0.0,rtol=1e-6,
                  terminate_on_equilibration=False, equ_threshold=1e-4,
                  terminate_on_efflorescence=False, eff_threshold=0.5,
                  terminate_on_locking=False, locking_threshold=400e-9, timer=None):
        if timer is None:
            self.timer = Timer(np.linspace(0.0, time, 10))
        else:
            self.timer = timer
        x0 = self.starting_state(radius,solute_concentration,particle_concentration)
        events = []

        if terminate_on_equilibration:
            dmdt = np.abs(self.update_state(0.0,x0)[self.layers-1])*np.exp(x0[self.layers-1])
            equilibrated = lambda time, x:  self.equilibrate(time, x, equ_threshold*dmdt)
            equilibrated.terminal = True
            events += [equilibrated]

        if terminate_on_efflorescence:
            efflorescing = lambda time, x: self.efflorescence(x) - eff_threshold
            efflorescing.terminal = True
            events += [efflorescing]

        if terminate_on_locking:
            shell_formation = lambda time, x: self.locking(x,locking_threshold)
            shell_formation.terminal = True
            events += [shell_formation]

        dxdt = lambda time, x: self.update_state(time,x)
        trajectory = solve_ivp(dxdt, (0,time), x0, rtol=rtol, events=events, method="Radau")
        return trajectory

    def complete_trajectory(self, trajectory):
        """Get the trajectory of all variables (including dependent ones) from a simulation (i.e.
        the output of UniformDroplet.integrate).

        Args:
            trajectory: the output of UniformDroplet.integrate, which gives the trajectory of independent
                        variables only.
        Returns:
            A pandas dataframe detailing the complete droplet history.
        """
        labels = ["radius","surface_temperature","solvent_mass","layer_mfs",
                  "mfs","layer_positions","layer_solute_concentrations",
                  "wet_layer_volumes","solute_masses","true_boundaries",
                  "layer_particle_concentrations","particle_volume_fraction","solvent_masses","layer_solvent_concentrations"]
        variables = {key: np.empty(trajectory.t.size, dtype=object) for key in labels}
        for i, state in enumerate(trajectory.y.T):
            earlier_droplet = GeneralDataDroplet(state, self.solution,self.suspension,self.particle_radius,self.layers)
            earlier_state = earlier_droplet.complete_state()
            for label, value in earlier_state.items():
                variables[label][i] = value

        variables['time'] = trajectory.t
        return pd.DataFrame(variables)

class GeneralDataDroplet:
    def __init__(self, state, solution, suspension, particle_radius, layers, particle_density=2200):
        state = np.array(state)
        self.solvent_masses = np.exp(state[:layers])
        self.solvent_mass = np.sum(self.solvent_masses)
        self.temperatures = state[layers:2*layers]
        self.solute_masses = np.exp(state[2 * layers:3 * layers])
        self.solute_mass = np.sum(self.solute_masses)
        self.layer_particle_mass = np.exp(state[3*layers:4*layers])
        self.particle_mass = np.sum(self.layer_particle_mass)

        radius = 0
        self.layer_positions = np.zeros(layers - 1)
        self.layer_volumes = layer_volumes(state, solution, suspension, particle_radius, layers)
        self.volume = np.sum(self.layer_volumes)
        for layer in range(layers-1):
            self.layer_positions[layer] = np.cbrt(3/(4*np.pi)*(self.layer_volumes[layer]-4/3*np.pi*radius**3))
            radius = self.layer_volumes[layer]

        self.mfs = self.solute_mass/(self.solvent_mass+self.solute_mass)
        self.wet_layer_volumes = self.layer_volumes - self.layer_particle_mass / particle_density
        self.densities = (self.solute_masses+self.solvent_masses)/(self.layer_volumes-self.wet_layer_volumes)
        self.radius = np.cbrt(self.volume*3/(4*np.pi))
        self.true_boundaries = np.concatenate(([0], self.layer_positions, [self.radius]))
        self.layer_particle_concentration = self.layer_particle_mass/self.layer_volumes

        self.layer_concentrations = self.solute_masses / (self.layer_volumes - self.wet_layer_volumes)
        self.layer_solvent_concentrations = self.solvent_masses /(self.layer_volumes-self.wet_layer_volumes)
        self.particle_volume_fraction = self.particle_mass/(self.volume * particle_density)
        self.layer_mfs = self.solute_masses / (self.solvent_masses + self.solute_masses)

    def complete_state(self):
        return dict(radius=self.radius, surface_temperature=self.temperatures[-1],
                    solvent_mass=self.solvent_mass, mfs=self.mfs, layer_mfs = self.layer_mfs,
                    layer_positions=self.layer_positions, layer_solute_concentrations=self.layer_concentrations,
                    wet_layer_volumes=self.wet_layer_volumes, solute_masses=self.solute_masses, solvent_masses=self.solvent_masses,
                    true_boundaries=self.true_boundaries, layer_particle_concentrations=self.layer_particle_concentration, particle_volume_fraction=self.particle_volume_fraction, layer_solvent_concentrations=self.layer_solvent_concentrations)
