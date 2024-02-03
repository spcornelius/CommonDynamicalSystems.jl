export MagneticPendulum, magnet_dist, magnetic_potential, energy, closest_magnet

using Distances
using StaticArrays

struct MagneticPendulum{T, N} <: AutonomousODESys
    ω::T
    α::T
    h::T
    r_mag::NTuple{N, SVector{2, T}}

    function MagneticPendulum{T}(ω, α, h, r_mag...) where {T}
        N = length(r_mag)
        N > 0 || error("Must have at least one magnet.")
        all(rₘ -> length(rₘ) == 2, r_mag) || 
            error("Each magnet position must be a vector or tuple of length 2.")
    
        # convert each element of r_mag to a static vector of length 2
        r_mag = map(SVector{2}, r_mag)
        new{T, N}(ω, α, h, r_mag)
    end
end

function MagneticPendulum(ω, α, h, r_mag...)
    T = Base.promote_eltype(ω, α, h, r_mag...)
    MagneticPendulum{T}(ω, α, h, r_mag...)
end

# type conversion
MagneticPendulum{T}(sys::MagneticPendulum) where {T} =
    MagneticPendulum{T}(sys.ω, sys.α, sys.h, sys.r_mag...)

@inline magnet_dist(r, rₘ, h) = sqrt(sqeuclidean(r, rₘ) + h^2)

function rhs(du, u, sys::MagneticPendulum)
    @unpack ω, α, h, r_mag = sys
    @views r, v, dr, dv = u[1:2], u[3:4], du[1:2], du[3:4]
    @.. dr = v
    @.. dv = -ω^2 * r - α * v
    for rₘ in r_mag
        d³ = magnet_dist(r, rₘ, h)^3
        @.. dv += (rₘ - r)/d³
    end
    nothing
end

function jac(J::AbstractMatrix{T}, u, sys::MagneticPendulum) where {T}
    @unpack ω, α,  h, r_mag = sys
    r = @view u[1:2]

    fill!(J, zero(T))

    J[1, 3] = one(T)
    J[2, 4] = one(T)

    J[3, 1] = -ω^2
    J[4, 2] = -ω^2

    J[3, 3] = -α
    J[4, 4] = -α

    x, y = r
    for rₘ in r_mag
        xₘ, yₘ = rₘ
        d = magnet_dist(r, rₘ, h)
        c1 = inv(d^3)
        c2 = d^5
        c3 = 3 * (xₘ - x) * (yₘ - y)/c2
        c4 = 3/c2
        J[3, 1] += c4 * (xₘ - x)^2 - c1
        J[3, 2] += c3
        J[4, 1] += c3
        J[4, 2] += c4 * (yₘ - y)^2 - c1
    end
    nothing
end

function magnetic_potential(r, rₘ, h)
    x, y = r
    xₘ, yₘ = rₘ

    return -1/sqrt((x - xₘ)^2 + (y - yₘ)^2 + h^2)
end

function energy(u, p::MagneticPendulum)
    @unpack ω, α, h, r_mag = p
    @views r, v = u[1:2], u[3:4]

    # kinetic and simple-harmonic potential energies
    E = (v[1]^2 + v[2]^2)/2 + ω^2 * (r[1]^2 + r[2]^2)/2
    for rₘ in r_mag
        E += magnetic_potential(r, rₘ, h)
    end

    return E 
end

function closest_magnet(r, p::MagneticPendulum{<:Any, m}) where {m}
    @boundscheck begin 
        n = length(r)
        n == 2 || error("Supplied r has length $n. Expected r = (x, y).")
    end
    all(isfinite, r) ? argmin(euclidean(r, rₘ) for rₘ in p.r_mag) : 0
end