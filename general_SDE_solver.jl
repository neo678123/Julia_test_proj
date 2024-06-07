using Distributions
using Random
using Plots
using LinearAlgebra
using BenchmarkTools

function getRandomGaussian(mean, variance, num)
    d = Normal(mean, variance)
    return rand(d, num)
 end
 

function solveSDE(mu, sigma, times, initialCondition)
    N = length(times)
    dim = length(initialCondition)
    normals = reshape(getRandomGaussian(0, 1, (N-1) * dim), (dim, N-1))

    solution = hcat(initialCondition, zeros((dim, N-1)))
    for i = 2:N
        dt = times[i] - times[i-1]
        solution[1:dim, i] = solution[1:dim, i-1] + mu(times[i-1], solution[1:dim, i-1]) * dt +
                                sigma(times[i-1], solution[1:dim, i-1]) * normals[1:dim, i-1] * sqrt(dt)
    end

    return transpose(solution)
end

dt = 0.01
T = 3
times = [0:dt:T;]
dim = 500

X0 = 300
mu = 0.75
sigma = 0.30

solution = solveSDE((t,X) -> mu * X, (t, X) -> sigma * Diagonal(X), times, X0 * ones(dim))
plot(times, solution)