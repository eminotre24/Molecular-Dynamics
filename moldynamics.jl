#Molecular Dynamics Using the Langevin Equation - moldynamics Module
#Created by: AE Novales T
#Date: 11/12/2024
#This is a julia module created to store all the functions utilized for simulating a molecular system,
# with periodic boundaries and interactions given by the Lennard-Jones Potential.

module moldynamics
    #Modules Used for this Module
    using Plots, LaTeXStrings, Statistics, Distributions

    export verlet_langevin_system, anim_system, anim_system_gif, unwrap_positions, real_time

    #Internal Functions
    function time_gen(time_data::NTuple{4, Float64})
        """
            Function to generate both the saving lenght and the resolution lenght. The saving lenght indicates how much steps will be saved, that is the "real"
            simulation steps, while the resolution steps is used to give a better "resolution" of the simulation, that is a better output of the data with a much
            smaller timestep
            Args:
                time_data: a tuple with the the time parameters for the simulation

            Returns:
                time_vector: time vector in resolution_steps 
                save_vector: time vector in save_steps 
                time_length: lenght of the time in resolution_steps 
                save_length: lenght of the time in save_steps 
        """
        time_vector = Float64.(time_data[1]:time_data[3]:time_data[2])
        time_length = length(time_vector)
        save_vector = zeros(Bool, time_length)
        save_length = 0
        for i ∈ 1:Int(time_data[4]/time_data[3]):time_length
            save_vector[i] = true
            save_length += 1
        end
        return time_vector, save_vector, time_length, save_length
    end

    function periodic_bound(x::Vector{Float64}, L::Int)
        #Periodic Boundaries, map the new position in the box
        return mod.(x, L)
    end

    function periodic_distance(x1::Float64, x2::Float64, box_size::Int)
        distances = [x1 - x2, x1 - (x2 - box_size), x1 - (x2 + box_size)]
        return distances[argmin(abs.(distances))]
    end

    function lenjon_forces(x::Float64, y::Float64, x_interaction::Float64, y_interaction::Float64, interaction_parameter::Float64, diameter::Float64, box_size::Int, cutoff::Float64)
        x_dif = periodic_distance(x, x_interaction, box_size)
        y_dif = periodic_distance(y, y_interaction, box_size)
        r = sqrt(x_dif^2 + y_dif^2)
        if r > cutoff
            return 0, 0
        else
            factor6 = (diameter/r)^6
            radforce = 24*interaction_parameter*(2 * factor6^2 + factor6)/r^2 #The extra 1/r is because of the conversion to cartesian coordinates
            xforce = radforce*x_dif
            yforce = radforce*y_dif
            return xforce, yforce
        end 
    end

    function total_interactions_forces(x::Vector{Float64}, y::Vector{Float64}, num_particles::Int, box_size::Int)
        forcesx = zeros(Float64, num_particles)
        forcesy = zeros(Float64, num_particles)
        for i ∈ 1:num_particles
            #Positions for particle analyzing
            xi = x[i]
            yi = y[i]
            for j ∈ 1:num_particles
                if i == j #Skip itself
                    continue
                end
                #Here go all forces of interaction
                fx, fy = lenjon_forces(xi, yi, x[j], y[j], 1.0, 0.05, box_size, 1.0)
                forcesx[i] += fx
                forcesy[i] += fy
            end
        end
        return forcesx, forcesy
    end    
    

    #Functions to Export
    function verlet_langevin_system(time_data::NTuple{4, Float64}, parameters_langevin::NTuple{4, Float64}, number_of_particles::Int, box_size::Int)
        """
            Core function, which generates the simulation of the system
            Args:
                time_data: a tuple with the the time parameters for the simulation
                parameters_langevin: a tuple with the constant parameters of the system
                number_of_particles: Number of particles of the system
                box_size: The lenght of the box in which particles will be confined

            Returns:
                x: Matrix with the x positions of the number_of_particles' particles ([particle, time]), given as the points saved using the save_step
                y: Matrix with the y positions of the number_of_particles' particles ([particle, time]), given as the points saved using the save_step
                vx: Matrix with the vx values of the number_of_particles' particles ([particle, time]), given as the points saved using the save_step
                vy: Matrix with the vy values of the number_of_particles' particles ([particle, time]), given as the points saved using the save_step
                time: time vector in resolution step
                save_length: lenght of the time in save_steps 
        """
        time, save, time_length, save_length = time_gen(time_data)
        dt = time_data[3]
    
        x = Matrix{Float64}(undef, number_of_particles, save_length)
        y = Matrix{Float64}(undef, number_of_particles, save_length)
        vx = Matrix{Float64}(undef, number_of_particles, save_length)
        vy = Matrix{Float64}(undef, number_of_particles, save_length)
    
        x[:, 1] = rand(number_of_particles) .* box_size
        y[:, 1] = rand(number_of_particles) .* box_size
        vx[:, 1] .= zeros(number_of_particles)
        vy[:, 1] .= zeros(number_of_particles)
    
        γ = parameters_langevin[1]
        m = parameters_langevin[2]
        t = parameters_langevin[3]
        kb = parameters_langevin[4]
        a = (1-γ*dt/2m)/(1+γ*dt/2m)
        b = (1+γ*dt/2m)^-1
        g = sqrt(2*γ*kb*t)
    
        save_counter = 2
        x_zero = x[:, 1]
        y_zero = y[:, 1]
        vx_zero = vx[:, 1]
        vy_zero = vy[:, 1]
    
        for i ∈ 2:time_length
            #Weiner Factor dη = gdW
            dWx = rand(Normal(0, dt), number_of_particles)
            dWy = rand(Normal(0, dt), number_of_particles)
    
            fxpast, fypast = total_interactions_forces(x_zero, y_zero, number_of_particles, box_size)
    
            #Positions
            px = x_zero + dt*b.*vx_zero + b*(dt^2)/(2m).*fxpast + b*dt/(2m)*g.*dWx 
            py = y_zero + dt*b.*vy_zero + b*(dt^2)/(2m).*fypast + b*dt/(2m)*g.*dWy 
            x_one = periodic_bound(px, box_size)
            y_one = periodic_bound(py, box_size)
            
            fxnew, fynew = total_interactions_forces(x_one, y_one, number_of_particles, box_size)
        
            #Velocities
            vx_one = a.*vx_zero + dt/(2m).*(a.*fxpast .+ fxnew) + b/m*g.*dWx
            vy_one = a.*vy_zero + dt/(2m).*(a.*fypast .+ fynew) + b/m*g.*dWy
    
            if save[i]
                x[:, save_counter] = x_one[:]
                y[:, save_counter] = y_one[:]
                vx[:, save_counter] = vx_one[:]
                vy[:, save_counter] = vy_one[:]
                save_counter += 1
            end
    
            x_zero[:] = x_one[:]
            y_zero[:] = y_one[:]
            vx_zero[:] = vx_one[:]
            vy_zero[:] = vy_one[:]
        end
        return x, y, vx, vy, time, save_length
    end

    function anim_system(x::Matrix{Float64}, y::Matrix{Float64}, save_length::Int, box_size::Int, filename::String)
        """
            Animate the particle system in mp4 format
            Args:
                x: Matrix with the x positions of the number_of_particles' particles ([particle, time]), given as the points saved using the save_step
                y: Matrix with the y positions of the number_of_particles' particles ([particle, time]), given as the points saved using the save_step
                save_length: lenght of the time in save_steps 
                filename: String of the name of the file omiting the format
        """
        plot()
        number_of_particles = size(x, 1)
        anim = @animate for i ∈ 1:save_length
            scatter(x[:, i], y[:, i], title="Trayectories of the Particles", xlabel=L"x", ylabel=L"y", color=1:number_of_particles, xlims=(0, box_size), ylims=(0, box_size), label=false, dpi=300)
        end
        gif(anim, filename * ".mp4", fps = 30)
    end

    function anim_system_gif(x::Matrix{Float64}, y::Matrix{Float64}, save_length::Int, box_size::Int, filename::String)
        """
            Animate the particle system in gif format. Used for the GitHub repository description.
            Args:
                x: Matrix with the x positions of the number_of_particles' particles ([particle, time]), given as the points saved using the save_step
                y: Matrix with the y positions of the number_of_particles' particles ([particle, time]), given as the points saved using the save_step
                save_length: lenght of the time in save_steps 
                filename: String of the name of the file omiting the format
        """
        plot()
        number_of_particles = size(x, 1)
        anim = @animate for i ∈ 1:save_length
            scatter(x[:, i], y[:, i], title="Trayectories of the Particles", xlabel=L"x", ylabel=L"y", color=1:number_of_particles, xlims=(0, box_size), ylims=(0, box_size), label=false, dpi=150)
        end
        gif(anim, filename * ".gif", fps = 15)
    end

    function unwrap_positions(positions::Vector{Float64}, box_size::Int)
        """
            Unwrapp a vector of positions for a single particle of the system, used to get the "real" trayectory of the particle in that axis. We
            take advantaje of the symetry of the box to generalize this function for both x and y
            Args:
                positions: Vector of a single particle's position in whichever axis, both x or y
                box_size: The lenght of the box in which particles will be confined
            Returns:
                unwrapped: the position vector of a single particle unwrapped, without jumps considering it passed to another zone
        """
        unwrapped = copy(positions)
        for i in 2:length(positions)
            Δ = unwrapped[i] - unwrapped[i-1]
            if Δ > box_size / 2
                # Crossed the upper boundary
                unwrapped[i:end] .-= box_size
            elseif Δ < -box_size / 2
                # Crossed the lower boundary
                unwrapped[i:end] .+= box_size
            end
        end
        return unwrapped
    end
    
    function real_time(time_data::NTuple{4, Float64}, save_length::Int)
        """
            Return the time vector as of the save_step, instead of the resolution step
            Args:
                time_data: a tuple with the the time parameters for the simulation
                save_length: lenght of the time in save_steps 
            Returns:
                Vector with the time given in the save_steps
        """
        return Float64.(range(time_data[1], time_data[2], length=save_length))
    end
end