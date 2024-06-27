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

function symplecticRK2Scheme(q0::Float64, p0::Float64, n0::Float64, h::Float64, N::Int32,
    A::Matrix{Float64}, b::Vector{Float64})
    q = vcat(q0, zeros(N-1))
    p = vcat(p0, zeros(N-1))

    kSquared = Float64(n0^2 - 1)
    H = hamiltonian(q0, p0, n0)

    # This should be read transposed!
    M = hcat([1, -kSquared/H * h * A[1,1], 0, -kSquared/H * h * A[2,1]],
             [h/H * A[1,1], 1, h/H * A[2,1], 0],
             [0, -kSquared/H * h * A[1,2], 1, -kSquared/H * h * A[2,2]],
             [h/H * A[1,2], 0, h/H * A[2,2], 1]
        )

    # This matrix is used when q > 1, as dH/dq = 0 then
    L = hcat([1, -kSquared/H * h * A[1,1], 0, -kSquared/H * h * A[2,1]],
             [0, 1, 0, 0],
             [0, -kSquared/H * h * A[1,2], 1, -kSquared/H * h * A[2,2]],
             [0, 0, 0, 1]
        )

    # Expensive, but still O(1)
    # TODO: there are initial conditions where L is not needed
    Minv = inv(M)
    Linv = inv(L)

    for n = 2:N
        if abs(q[n-1]) <= 1
            v = [-p[n-1]/H, kSquared/H * q[n-1]]
            k = Minv * vcat(v,v)
        else
            v = [-p[n-1]/H, 0]
            k = Minv * vcat(v,v)
        end

        q[n] = q[n-1] + h * (b[1] * k[1] + b[2] * k[3])
        p[n] = p[n-1] + h * (b[1] * k[2] + b[2] * k[4])
    end

    return [q,p]
end

h = Float64(0.01)
zMax = Float64(6pi)

n0 = Float64(1.4)
q0 = Float64(1.2)
p0 = Float64(0.03)

N = Int32(ceil(zMax/h))
kSquared = Float64(n0^2 - 1)

A = hcat([1/4, 1/4 + 1/(2 * sqrt(3))], [1/4 - 1/(2 * sqrt(3)), 1/4])
b = [1/2, 1/2]

phaseSpace = symplecticRK2Scheme(q0, p0, n0, h, N, A, b)
qValues = [phaseSpace[1]]
pValues = [phaseSpace[2]]
for i = 1:9
    sol = symplecticRK2Scheme(q0 - Float64(i * 0.1), p0, n0, h, N, A, b)
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

plot(times, qValues, label = labels, xlabel = "t", ylabel = "q(t)", legend=:bottomright)
# plot(times, pValues, label = labels, xlabel = "t", ylabel = "p(t)")
# plot(qValues, pValues, label = labels, xlabel = "q(t)", ylabel = "p(t)")
# plot(times, map((q, p) -> hamiltonian(q, p, n0), qValues[10], pValues[10]), label = "...",
#        xlabel = "t", ylabel = "H(q(t), p(t))")