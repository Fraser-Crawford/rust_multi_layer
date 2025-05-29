# Installation:
To get the model running locally you'll need to install two wheels, one for the rust model backend and one for the python interface.
1. Download the repo (or clone it to your system to keep up to date with changes).
2. In target/wheels, copy the path of the .whl file and in your python environment command line, type `pip install (PATH OF .WHL FILE)`
3. In frontend/dist, copy the path of the other .whl file and likewise `pip install (PATH OF OTHER .WHL FILE)` \
***If you already have a previous version installed, add `--force-reinstall` to these commands to make sure you overwrite them***
4. You should now have the model installed!

# Usage:
The frontend can be imported into your program using the name `general_model`
-> `from general_model import *`
general_model has 3 components: `Timer`, which allows you to model simulation progress, `GeneralDroplet`, which contains the environment for the simulation and `GeneralDataDroplet` which deals with the history of a droplet in a simulation.

Here is an example simulation:
```
from general_model import *
droplet = GeneralDroplet(293, 0.5, 0.0, 10, solution="nacl", suspension="silica", particle_radius=200e-9)
trajectory = droplet.integrate(100, 23e-6, solute_concentration=25, particle_concentration=20, timer = None)
dataframe = droplet.complete_trajectory(trajectory)
```
Line one sets the environment:
1. Temperature is set to 293 K
2. Relative humidity is set to 50% (0.5)
3. Air flow speed is set to 0.0
4. The droplet is made of 10 layers
5. The solution is `"nacl"`, NaCl in water. If not specified, pure water is assumed.
6. There are silica nanoparticles in solution. If not specified then a dummy system of silica NP is assumed. 
7. These particles have a radius of 200e-9, if not specified then a dummy value of 200e-9 is used.

***All units are given in SI***

Line two sets the simulation running:
1. Simulate 100 seconds
2. Set the initial droplet radius to 23e-6 = 23 microns
3. Set the concentration of solute to 25 kg / m3 = 25 g / L
4. Set the concentration of nanoparticles to 20 kg / m3 = 20 g / L
5. Use the default timer (prints when 10% of the simulation time has passed)

Line three creates a dataframe from the simulation data. `dataframe.time` is the series of times, `dataframe.layer_solute_concentrations` gives the solute concentration in each layer for each timestep etc. 

# Advanced Usage:
Sometimes you may want the droplet environment to change over time such as the temperature or relative humidity changing. To do this, provide a function that takes in the time as the only argument to the respective environment variable.
For instance:
```
def temperature_curve(time, t_min):
    if time < 0.1:
        return 293 - time*10*(293-t_min)
    else:
        return t_min
droplet = GeneralDroplet(lambda time: temperature_curve(time,253), 0.5, 0.0, 10, solution="nacl", suspension="silica", particle_radius=200e-9)
```
This allows the simulation to change the temperature with time, here the temperature decreases over the first 0.1 seconds from 293 K to t_min. A `lambda` expression can be used to turn a many argument function into a single argument function.
***This works for temperature, relative humidity and air speed***. 

There are 3 enablable terminal events for the system as well. When these are detected, the simulation will stop early. When calling droplet.integrate, you can set the flags to enable possible termination:
1. `terminate_on_equilibration` will terminate the simulation if the rate of mass loss goes below the `equ_threshold` times the initial mass loss rate
2. `terminate_on_efflorescence` will terminate the simulation if the activity of the droplet's surface becomes less than `eff_threshold`
3. `terminate_on_locking` will terminate the simulation if a shell forms of nanoparticles above their critical_volume_fraction with thickness equal to `locking_threshold`

droplet.integrate can also take a `timer` keyword. `Timer` takes an array of values and will print whenever the simulation crosses these thresholds. By default, the `timer` is np.linspace(0,max_time,10), giving 10 equal spaced thresholds over the simulation time.

# Current solutions and suspension definitions
`"nacl"` is salt in water
`"water"` is pure water
`"sucrose_zobrist"` and `"sucrose"` are two different sucrose solutions. The zobrist parametrisation is less diffusive.


