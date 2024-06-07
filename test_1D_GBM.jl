using Distributions
using Random
using Plots

function mu(t, X)
   return 0.75 * X
end

function sigma(t, X)
   return 0.30 * X 
end

function getRandomGaussian(mean, variance, num)
   d = Normal(mean, variance)
   return rand(d, num)
end

function getSDESolution(X0, T, dt)
   N = floor(Int64, T/dt)
   dtSqrt = sqrt(dt)
   solution = [[X0]; zeros((N-1))]
   dW = dtSqrt * getRandomGaussian(0, 1, N)
   for i = 1:N-1
      t = i * dt
      X = solution[i]
      solution[i+1] = X + mu(t, X) * dt + sigma(t, X) * dW[i+1]
   end

   return solution
end

dt = 0.01
T = 3
X0 = 300


times = [0:dt:T-dt;]
plot(times, getSDESolution(X0, T, dt))