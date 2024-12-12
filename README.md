# Molecular-Dynamics
This is a Julia module to generate a simulation of a set of particles with periodic boundaries, simulating a molecular system using the Langevin thermostat and the Lennard-Jones Potential. The repository includes the module "moldynamics", which contains both the internal and the exported functions, for which the user can generate a simulation type as stated before. Most of the code was written in order to be intuitive to understand, and has the capacity to be mnodified to add other types of interactions. The user can implement their own interactive force, or external force, and the use the internal function "total_interactions_forces" to add it, or add it directly in the Verlet Algorithm

To simulate this molecular system, we use the Langevin equation, which is:

$$
m \ddot{\vec{v}} = \vec{F} - \gamma \dot{\vec{v}} + \vec{\eta}
$$

In which the case of this system, we have a bidimensional system $\vec{v} = (x, y)$, and the force is given as the sum of all the interactions with neighboor particles, and an external force if present. The simulation was created using the Verlet Algorithm, as stated before, with the following equations for the steps:

$$
x(t + \Delta t) = x(t) + \Delta t v(t) + \frac{b \Delta t^2}{2m} F(t) + \frac{b \Delta t}{2m} g\Delta W(t)
$$

$$
v(t + \Delta t) = a v(t) + \frac{\Delta t}{2m} \big(a F(t) + F(t + \Delta t)\big) + \frac{b}{m} g\Delta W(t)
$$

Where we create the following constants just to reduce the amount of terms in the equations:

$$
a \equiv \left(1 - \frac{\gamma \Delta t}{2m}\right)\left(1 + \frac{\gamma \Delta t}{2m}\right)^{-1}
$$

$$
b \equiv \left(1 + \frac{\gamma \Delta t}{2m}\right)^{-1}
$$

Theres also a "tutorial" in which the operation is more comprehensible, and the user can modify the parameters to obtain different results as of the simulation. Here we show the result of this tutorial, an animation of the particles' system:
