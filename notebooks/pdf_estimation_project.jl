# DON'T DO THIS IN PRODUCTION CODE
# This is just to ensure that the required packages are installed (for the notebook to run)
# We had some issues on JupyterLab where the environment was not properly activated

# These are packages we think you will need
using Pkg
Pkg.activate(joinpath(@__DIR__, "."))
Pkg.instantiate()
Pkg.add("Plots")
Pkg.add("Distributions")
Pkg.add("Parameters")
Pkg.add("Statistics")
Pkg.add("Random")

using Parameters
using Distributions

import Random: rand, Dims, AbstractRNG
using Plots

@with_kw struct ProjectileParameters{HeightType, SpeedType, AngleType, MassType, DragType}
    initial_height::HeightType
    initial_speed::SpeedType
    angle_in_degrees::AngleType
    mass::MassType
    drag::DragType
end

function rand(uqParam::ProjectileParameters)
    h0 = uqParam.initial_height isa Real ? uqParam.initial_height : rand(uqParam.initial_height)
    v0 = uqParam.initial_speed isa Real ? uqParam.initial_speed : rand(uqParam.initial_speed)
    angle = uqParam.angle_in_degrees isa Real ? uqParam.angle_in_degrees : rand(uqParam.angle_in_degrees)
    m = uqParam.mass isa Real ? uqParam.mass : rand(uqParam.mass)
    ρ = uqParam.drag.ρ isa Real ? uqParam.drag.ρ : rand(uqParam.drag.ρ)
    C_d = uqParam.drag.C_d isa Real ? uqParam.drag.C_d : rand(uqParam.drag.C_d)
    A = uqParam.drag.A isa Real ? uqParam.drag.A : rand(uqParam.drag.A)
    return ProjectileParameters(
        initial_height = h0,
        initial_speed = v0,
        mass = m,
        angle_in_degrees = angle,
        drag = Drag(ρ=ρ, C_d=C_d, A=A)
    )
end

function rand(rng::AbstractRNG, uqParam::ProjectileParameters)
    h0 = uqParam.initial_height isa Real ? uqParam.initial_height : rand(rng, uqParam.initial_height)
    v0 = uqParam.initial_speed isa Real ? uqParam.initial_speed : rand(rng, uqParam.initial_speed)
    angle = uqParam.angle_in_degrees isa Real ? uqParam.angle_in_degrees : rand(rng, uqParam.angle_in_degrees)
    m = uqParam.mass isa Real ? uqParam.mass : rand(rng, uqParam.mass)
    ρ = uqParam.drag.ρ isa Real ? uqParam.drag.ρ : rand(rng, uqParam.drag.ρ)
    C_d = uqParam.drag.C_d isa Real ? uqParam.drag.C_d : rand(rng, uqParam.drag.C_d)
    A = uqParam.drag.A isa Real ? uqParam.drag.A : rand(rng, uqParam.drag.A)
    return ProjectileParameters(
        initial_height = h0,
        initial_speed = v0,
        mass = m,
        angle_in_degrees = angle,
        drag = Drag(ρ=ρ, C_d=C_d, A=A)
    )
end

function rand(rng::AbstractRNG, uqParam::ProjectileParameters, dims::Dims)
    ElementType = typeof(rand(rng, uqParam.initial_height))
    ParameterType = ProjectileParameters{ElementType, ElementType, ElementType, Drag{ElementType, ElementType, ElementType}}
    results = Array{ParameterType}(undef, dims)
    for idx in CartesianIndices(results)
        results[idx] = rand(rng, uqParam)
    end
    return results
end

function rand(uqParam::ProjectileParameters, dims::Dims)
    ElementType = typeof(rand(uqParam.initial_height))
    ParameterType = ProjectileParameters{ElementType, ElementType, ElementType, ElementType, Drag{ElementType, ElementType, ElementType}}
    results = Array{ParameterType}(undef, dims)
    for idx in CartesianIndices(results)
        results[idx] = rand(uqParam)
    end
    return results
end

"""
    projectile_motion(parameters::ProjectileParameters; endtime=50.0, resolution=0.01)
Simulate projectile motion given a `ProjectileParameters` instance.
# Arguments
- `parameters`: Instance of `ProjectileParameters` containing simulation parameters.
# Keyword Arguments
- `endtime`: Maximum simulation time (s)
- `resolution`: Time step for numerical integration (s)
# Returns
- `t`: Time vector
- `x`: Horizontal position vector (m)
- `y`: Vertical position vector (m)
# Description
This function extracts parameters from the `ProjectileParameters` instance
and calls the more general `projectile_motion` function to perform the simulation.
"""
function projectile_motion(parameters::ProjectileParameters; endtime=50.0, resolution=0.01)
    return projectile_motion(
        parameters.initial_height,
        endtime,
        parameters.initial_speed,
        parameters.angle_in_degrees,
        parameters.mass;
        drag = parameters.drag, dt = resolution
    )
end

"""
    projectile_motion(initial_height, endtime, initial_speed, launch_angle, mass; drag = v -> 0.0, dt = 0.01)

Simulate projectile motion with optional drag force.

# Arguments
- `initial_height`: Initial height of the projectile (m)
- `endtime`: Maximum simulation time (s)
- `initial_speed`: Initial speed of the projectile (m/s)
- `launch_angle`: Launch angle in degrees
- `mass`: Mass of the projectile (kg)

# Keyword Arguments
- `drag`: Drag force function that takes velocity as input (default: no drag)
- `dt`: Time step for numerical integration (default: 0.01 s)

# Returns
- `t`: Time vector
- `x`: Horizontal position vector (m)
- `y`: Vertical position vector (m)

# Description
Numerically integrates the equations of motion for a projectile under gravity
and optional drag force using Euler's method. The simulation stops when the
projectile hits the ground (y < 0) or when `endtime` is reached.
"""
function projectile_motion(initial_height, endtime, initial_speed, launch_angle, mass; drag = v -> 0.0, dt = 0.01)
    g = 9.81  # Acceleration due to gravity (m/s^2)
    θ = deg2rad(launch_angle)  # Convert angle to radians

    vx0 = initial_speed * cos(θ)
    vy0 = initial_speed * sin(θ)
    
    # Time vector
    t = 0:dt:endtime
    position = zeros(length(t), 2)
    velocity = zeros(length(t), 2)
    position[1, 1] = 0.0  # Initial x position
    position[1, 2] = initial_height  # Initial y position
    velocity[1, 1] = vx0  # Initial x velocity
    velocity[1, 2] = vy0  # Initial y velocity
    for i in 1:length(t)-1
        v = sqrt(velocity[i, 1]^2 + velocity[i, 2]^2)
        drag_force = drag(v)
        ax = -drag_force * (velocity[i, 1] / v) / mass
        ay = -g - drag_force * (velocity[i, 2] / v) / mass

        velocity[i+1, 1] = velocity[i, 1] + ax * dt
        velocity[i+1, 2] = velocity[i, 2] + ay * dt

        position[i+1, 1] = position[i, 1] + velocity[i, 1] * dt
        position[i+1, 2] = position[i, 2] + velocity[i, 2] * dt

        if position[i+1, 2] < 0
            position = position[1:i+1, :]
            t = t[1:i+1]
            break
        end
    end
    x = position[:, 1]
    y = position[:, 2]

    return t, x, y
end

function projectile_distance(parameters::ProjectileParameters; kwargs...)
    t, x, y = projectile_motion(parameters; kwargs...)
    return maximum(x)
end

function projectile_distance(args; kwargs...)
    t, x, y = projectile_motion(args...; kwargs...)
    return maximum(x)
end


@with_kw struct Drag{RhoType, CdType, AreaType}
    ρ::RhoType = 1.225  # Air density (kg/m^3)
    C_d::CdType = 0.47  # Drag coefficient
    A::AreaType  = 0.01  # Cross-sectional area (m^2)
end

function (d::Drag)(v)
    return 0.5 * d.ρ * d.C_d * d.A * v^2
end


function optimal_mlmc_levels(tolerance, costs, variances)
    total_cost = sum(sqrt.(variances .* costs))
    required_samples = ceil.( (2 / tolerance^2) * total_cost * sqrt.(variances ./ costs) )
    return required_samples
end

function estimate_mlmc_variances(model, base_parameters, timestepsizes = [2^(-i) for i in 1:6], M = 1000)
    variances = []
    for (level, dt) in enumerate(timestepsizes)
        samples = []
        for sample in 1:M
            parameters = rand(base_parameters)
            fine_solution = model(parameters; resolution = dt)
            if level == 1
                push!(samples, fine_solution)
            else
                coarse_solution = model(parameters; resolution = timestepsizes[level-1])
                push!(samples, fine_solution - coarse_solution)
            end
        end
        push!(variances, var(samples))
    end
    return variances
end

function simulate_mlmc_samples(model, base_parameters, tolerance, variances, timestepsizes = [2^(-i) for i in 1:6])
    required_samples = optimal_mlmc_levels(tolerance, timestepsizes .^ (-1), variances)
    samples = []
    for (level, M) in enumerate(required_samples)
        push!(samples, [])
        for sample in 1:M
            parameters = rand(base_parameters)
            fine_solution = model(parameters; resolution = timestepsizes[level])
            if level == 1
                push!(samples[level], fine_solution)
            else
                coarse_solution = model(parameters; resolution = timestepsizes[level-1])
                push!(samples[level], (fine_solution, coarse_solution))
            end
        end
    end
    return samples
end

function mc_samples(model, base_parameters, M, timestepsize)
    samples = []
    for sample in 1:M
        parameters = rand(base_parameters)
        solution = model(parameters; resolution = timestepsize)
        push!(samples, solution)
    end
    return samples
end

function mc_estimate(mc_samples, g = x -> x)
    return mean(g.(mc_samples))
end

function mlmc_estimate(mlmc_samples, g = x -> x)
    mean_estimate = zero(typeof(g(mlmc_samples[1][1])))

    for level in 1:length(mlmc_samples)
        if level == 1
            mean_estimate += mean(g.(mlmc_samples[level]))
        else
            # last element is coarse, first is fine (only two elements per sample)
            mean_estimate += mean(g.(first.(mlmc_samples[level])) .- g.(last.(mlmc_samples[level])))
        end
    end
    return mean_estimate
end


resolutions = [2.0^(-i) for i in 1:6]

# Define the uncertain parameters
parameters = ProjectileParameters(
    initial_height = Uniform(8.0, 20.0),  # Initial height in meters
    initial_speed = Uniform(15.0, 25.0),   # Initial speed in m/s
    angle_in_degrees = Uniform(30.0, 60.0), # Launch angle in degrees
    mass = Uniform(0.5, 2.0),             # Mass in kg
    drag = Drag(
        ρ = Uniform(1.0, 1.3),            # Air density
        C_d = Uniform(0.4, 0.6),          # Drag coefficient
        A = Uniform(0.005, 0.02)          # Cross-sectional area
    )
)

variances = estimate_mlmc_variances(projectile_distance, parameters, resolutions, 5000)
variances_single_level = [var([projectile_distance(rand(parameters); resolution=res) for _ in 1:5000]) for res in resolutions]

tolerance = 0.01
monte_carlo_samples = mc_samples(projectile_distance, parameters, 10000, resolutions[end])
mlmc_samples = simulate_mlmc_samples(projectile_distance, parameters, tolerance, variances, resolutions)

mc_mean = mc_estimate(monte_carlo_samples)
mlmc_mean = mlmc_estimate(mlmc_samples)

println("MC Estimate: ", mc_mean)
println("MLMC Estimate: ", mlmc_mean)


