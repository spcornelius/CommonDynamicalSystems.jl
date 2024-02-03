export LiSprott

struct LiSprott{T <: Real} <: AutonomousODESys
    a::T
    b::T

    function LiSprott(a::aType, b::bType) where {aType <: Real, bType <: Real}
        T = Base.promote_type(aType, bType)
        new{T}(a, b)
    end
end

function rhs(dr, r::AbstractVector, p::LiSprott)
    x, y, z, u = r
    @unpack a, b = p
    dr[1] = -x + y
    dr[2] = -x * z + u
    dr[3] = x * y - a
    dr[4] = -b * y
    nothing
end

function jac(J::AbstractMatrix{T}, r::AbstractVector, p::LiSprott) where {T}
    x, y, z, _ = r
    @unpack b = p
    J[1, 1] = -one(T)
    J[1, 2] = one(T)
    J[2, 1] = -z
    J[2, 3] = -x
    J[2, 4] = one(T)
    J[3, 1] = y
    J[3, 2] = x
    J[4, 2] = -b
    nothing
end