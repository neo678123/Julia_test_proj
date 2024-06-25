using Plots
using BenchmarkTools

function refractiveIndex(q::Float64, n0::Float64)::Float64
    if abs(q) <= 1
        kSquared = Float64(n0^2 - 1)
        return sqrt(n0 - kSquared * q^2)
    end
    
    return Float64(1)
end

function hamiltonian(q::Float64, p::Float64, n0::Float64)::Float64
    return -sqrt(refractiveIndex(q, n0)^2 + p^2)
end

function symplecticEulerScheme(q0::Float64, p0::Float64, n0::Float64, h::Float64, N::Int32)
    q = vcat(q0, zeros(N-1))
    p = vcat(p0, zeros(N-1))

    kSquared = Float64(n0^2 - 1)
    H = hamiltonian(q0, p0, n0)

    for n = 2:N
        q[n] = q[n-1] - h / H * p[n-1]
        p[n] = p[n-1] + h * kSquared / H * q[n] 
    end

    return [q,p]
end

h = Float64(0.01)
zMax = Float64(6 * pi)

n0 = Float64(1.4)
q0 = Float64(1.1)
p0 = Float64(0)

N = Int32(ceil(zMax/h))

phaseSpace = symplecticEulerScheme(q0, p0, n0, h, N)
qValues = [phaseSpace[1]]
pValues = [phaseSpace[2]]
for i = 1:7
    phaseSpace = symplecticEulerScheme(q0 - Float64(i * 0.1), p0, n0, h, N)
    push!(qValues, phaseSpace[1])
    push!(pValues, phaseSpace[2])    
end


times = h * [0:1:N-1]

#println(phaseSpace)
plot(qValues, pValues)
plot(times, map((q, p) -> hamiltonian(q, p, n0), qValues[3], pValues[3]))
