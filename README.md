# Molecular-Dynamics
This is a Julia module to generate a simulation of a set of particles with periodic boundaries, simulating a molecular system using the Langevin thermostat and the Lennard-Jones Potential. The repository includes the module "moldynamics", which contains both the internal and the functions exported, for which the user can generate a simulation as stated before. Most of the code was written in order to be intuitive to understand, and has the capacity to be mnodified to add other types of interactions.

To simulate this molecular system, we use the Langevin equation as:

$$
m \ddot{\vec{x}} = \vec{F} - \gamma \dot{\vec{x}} + \vec{\eta}
$$

In which the case of this system, the force is given as the sum of all the interactions with neighboor particles. The simulation was created using the Verlet Algorithm.
