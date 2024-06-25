using Plots
using BenchmarkTools

function refractiveIndex(q::Float64, n0::Float64)::Float64
    if abs(q) <= 1
        kSquared = Float64(n0^2 - 1)
        return sqrt(n0^2 - kSquared * q^2)
    end
    
    return Float64(1)
end

function hamiltonian(q::Float64, p::Float64, n0::Float64)::Float64
    return -sqrt(refractiveIndex(q, n0)^2 - p^2)
end

function symplecticEulerScheme(q0::Float64, p0::Float64, n0::Float64, h::Float64, N::Int32)
    q = vcat(q0, zeros(N-1))
    p = vcat(p0, zeros(N-1))

    kSquared = Float64(n0^2 - 1)
    H = hamiltonian(q0, p0, n0)

    for n = 2:N
        #println("H = $(hamiltonian(q[n-1],p[n-1], n0))")
        q[n] = q[n-1] - h / H * p[n-1]
        #H = hamiltonian(q[n],p[n-1], n0)
        if q[n] <= 1
            p[n] = p[n-1] + h * kSquared / H * q[n] 
        else 
            p[n] = p[n-1]
        end
    end

    return [q,p]
end

h = Float64(0.01)
zMax = Float64(6pi)

n0 = Float64(1.4)
q0 = Float64(1.2)
p0 = Float64(0.03)

N = Int32(ceil(zMax/h))

phaseSpace = symplecticEulerScheme(q0, p0, n0, h, N)
qValues = [phaseSpace[1]]
pValues = [phaseSpace[2]]
for i = 1:9
    sol = symplecticEulerScheme(q0 - Float64(i * 0.1), p0, n0, h, N)
    push!(qValues, sol[1])
    push!(pValues, sol[2])    
end


times = h * [0:1:N-1]

#println(phaseSpace)
plot(qValues, pValues)
#plot(times, map((q, p) -> hamiltonian(q, p, n0), qValues[1], pValues[1]))
global labels = ["q₀ = $q0, p₀ = $p0"]
for i = 1:9
    global labels = hcat(labels, 
    ["q₀ = $(Float64(round(1000 * (q0 - Float64(i * 0.1)))/1000)), p₀ = $p0"])
end
println(labels)

plot(times, qValues, label = labels, xlabel = "t", ylabel = "q(t)")
# plot(times, pValues, label = labels, xlabel = "t", ylabel = "p(t)")
# plot(qValues, pValues, label = labels, xlabel = "q(t)", ylabel = "p(t)")
# plot(times, map((q, p) -> hamiltonian(q, p, n0), qValues[10], pValues[10]), label = "...",
#        xlabel = "t", ylabel = "H(q(t), p(t))")