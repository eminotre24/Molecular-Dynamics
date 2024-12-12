#Test of the Module "moldynamics"
#Created by: AE Novales T
#Date: 11/12/2024

using .moldynamics
using Plots, LaTeXStrings

#Definitions of Parameters for Langevin Equation
#Default - Not considering units
viscosity = 0.5
mass = 1.0
temperature = 10.0
boltzmann_constant = 1.0
parameters_langevin = (viscosity, mass, temperature, boltzmann_constant)

#Simulation Parameters
box_size = 30
number_of_particles = 50

#Time Parameters
initial_time = 0.0
final_time = 50.0
time_resolution = 0.01
save_step = 0.1

time_data = (initial_time, final_time, time_resolution, save_step)

x, y, vx, vy, t, s = verlet_langevin_system(time_data, parameters_langevin, number_of_particles, box_size)

#Proccesing Trajectories
#One should "unwrap" the positions, as the periodic boundaries make the particles do "jumps" from one boundary to another

#Example with the x axis
time_vector = real_time(time_data, s)
plot(dpi=300)
for i in 1:number_of_particles
    pos = unwrap_positions(x[i, :], box_size)
    plot!(time_vector, pos, label="")
end

title!(L"Trajectory of the Particles in $x$")
xlabel!(L"t")
ylabel!(L"x")
savefig("example_trayectories_x.png")

#Create an Animation of the Particle's System
#anim_system_(x, y, s, box_size, "test_animation") #Better resolution in mp4 format
anim_system_gif(x, y, s, box_size, "test_animation_git")