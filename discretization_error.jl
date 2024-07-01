using Plots
using LinearAlgebra
using LaTeXStrings
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
        q[n] = q[n-1] - h / H * p[n-1]

        if abs(q[n]) < 1
            p[n] = p[n-1] + h * kSquared / H * q[n] 
        else 
            p[n] = p[n-1]
        end
    end

    return [q,p]
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
        if abs(q[n-1]) < 1
            v = [-p[n-1]/H, kSquared/H * q[n-1]]
            k = Minv * vcat(v,v)
        else
            v = [-p[n-1]/H, 0]
            k = Linv * vcat(v,v)
        end

        q[n] = q[n-1] + h * (b[1] * k[1] + b[2] * k[3])
        p[n] = p[n-1] + h * (b[1] * k[2] + b[2] * k[4])
    end

    return [q,p]
end

function actualSolution(q0::Float64, p0::Float64, n0::Float64, h::Float64, N::Int32) 
    q = vcat(q0, zeros(N-1))
    p = vcat(p0, zeros(N-1))

    k = sqrt(Float64(n0^2 - 1))
    H = hamiltonian(q0, p0, n0)

    for n = 2:N
        x = cos(k/H * n * h)
        y = sin(k/H * n * h)
        q[n] = q0 * x - p0 * y
        p[n] = q0 * y + p0 * x
    end

    return [q,p]
end

function getL2ErrorSE(q0::Float64, p0::Float64, n0::Float64, h::Float64, N::Int32)
    actualSol = actualSolution(q0, p0, n0, h, N)
    SESol = symplecticEulerScheme(q0, p0, n0, h, N)

    error = reduce(vcat, actualSol .- SESol)
    return norm(error, 2)
end

function getL2ErrorRK2(q0::Float64, p0::Float64, n0::Float64, h::Float64, N::Int32,
    A::Matrix{Float64}, b::Vector{Float64})
    actualSol = actualSolution(q0, p0, n0, h, N)
    RK2Sol = symplecticRK2Scheme(q0, p0, n0, h, N, A, b)
    
    error = reduce(vcat, actualSol .- RK2Sol)
    return norm(error, 2)
end

zMax = Float64(pi)

n0 = Float64(1.4)
q0 = Float64(0.8)
p0 = Float64(0)

kSquared = Float64(n0^2 - 1)

A = hcat([1/4, 1/4 + 1/(2 * sqrt(3))], [1/4 - 1/(2 * sqrt(3)), 1/4])
b = [1/2, 1/2]

hValues = [0.4:-0.005:0.001;]
SEerrors = map(h -> getL2ErrorSE(q0, p0, n0, h, Int32(ceil(zMax/h))), hValues)
RK2errors = map(h -> getL2ErrorRK2(q0, p0, n0, h, Int32(ceil(zMax/h)), A, b), hValues)
println(hValues)
println(SEerrors)

labels = permutedims(["h = $h" for h in hValues])
print(labels)
plot(hValues, SEerrors, seriestype=:scatter, label = "Symplectic Euler", ylabel = L"|| (e_n) ||_{\ell^2}", xlabel = L"h")
plot!(hValues, RK2errors, seriestype=:scatter, label = "RK2")
png("L2_discretisation_errors_q0=0.8")