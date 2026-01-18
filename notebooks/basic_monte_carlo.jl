# Make sure the environment is set up
import Pkg
Pkg.activate("..")
Pkg.instantiate()
Pkg.add("Distributions")
Pkg.add("Statistics")
Pkg.add("Random")
Pkg.add("Plots")

using Random
using Statistics

f(x) = sin(pi * x)

M = 10000  # Number of samples
approximate_integral = 0.0
for i in 1:M
    x = rand()  # Sample uniformly in [0, 1]
    approximate_integral += f(x)
end
approximate_integral /= M  # Average over samples
println("Estimated integral of sin(πx) over [0, 1]: ", approximate_integral)
println("Exact integral: ", 2 / pi)



using Random
using Statistics
using Distributions

# Define the function to integrate
I(x, y) = x^2 + y^2 <= 1 ? 1.0 : 0.0

M = 10000  # Number of samples
approximate_pi = 0.0
for i in 1:M
    x = rand(Uniform(-1, 1))  # Sample x uniformly in [-1, 1]
    y = rand(Uniform(-1, 1))  # Sample y uniformly in [-1, 1]
    approximate_pi += I(x, y)
end
approximate_pi = (approximate_pi / M) * 4  # Scale by area
println("Estimated value of π: ", approximate_pi)
println("Exact value of π: ", Float64(pi))


using Random
using Statistics
using Distributions
using Plots
M_values = collect(2 .^ (4:15))
reruns = 10
errors = zeros(length(M_values), reruns)
for (j, M) in enumerate(M_values)
    for r in 1:reruns
        approximate_pi = 0.0
        for i in 1:M
            x = rand(Uniform(-1, 1))  # Sample x uniformly in [-1, 1]
            y = rand(Uniform(-1, 1))  # Sample y uniformly in [-1, 1]
            approximate_pi += x.^2 + y.^2 <= 1 ? 1.0 : 0.0
        end
        approximate_pi = (approximate_pi / M)  * 4 # Scale by area
        errors[j, r] = abs(approximate_pi - pi)
    end
end

mean_errors = mean(errors, dims=2)
std_errors = std(errors, dims=2)

plot(M_values, mean_errors[:], yerror=std_errors[:], xscale=:log10, yscale=:log10,
     xlabel="Number of samples M", ylabel="Absolute error",
     title="Convergence of Monte Carlo estimate of π",
     label="Mean absolute error with std dev",
     legend=:topright)
variance = pi/4 - pi^2 / 16
theoretical_std = 4 * sqrt.(variance ./ M_values)
ylims!(minimum(mean_errors) / 4, maximum(mean_errors) * 2)
plot!(M_values, theoretical_std, lw=2, ls=:dash, label="Theoretical std dev")
