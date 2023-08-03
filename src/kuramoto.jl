export Kuramoto

using Graphs
using SparseArrays

struct Kuramoto{AType, ωType} <: AutonomousODESys
    A::AType
    ω::ωType

    function Kuramoto(A::AType, ω::ωType) where {AType <: AbstractSparseMatrix,
                                                 ωType}
        nr, nc = size(A)
        nr == nc || error("Interaction matrix A must be square.")
        if length(ω) > 1
            length(ω) == nr || error("Dimension of ω must be compatible with A.")
        end
        A = copy(A)
        ω = copy(ω)
        @views fill!(A[diagind(A)], zero(eltype(A)))
        new{AType, ωType}(A, ω)
    end
end

# arbitrary (not-necessarily sparse) input for A
Kuramoto(A::AbstractMatrix, ω) = Kuramoto(sparse(A), ω)

# zero-frequency case
Kuramoto(A::AbstractMatrix{T}) where {T} = Kuramoto(A, zero(T))

# from a graph
Kuramoto(g::AbstractGraph, ω) = Kuramoto(adjacency_matrix(g), ω)

# zero-frequency case
Kuramoto(g::AbstractGraph) = Kuramoto(adjacency_matrix(g))

function rhs(dθdt, θ, sys::Kuramoto)
    @unpack A, ω = sys
    @.. dθdt = ω
    rows = rowvals(A)
    nzv = nonzeros(A)
    for j in 1:size(A, 2), r in nzrange(A, j)
        i = rows[r]
        Aij = nzv[r]
        dθdt[j] += Aij * sin(θ[i] - θ[j])
    end
    nothing
end

function jac(J::AbstractMatrix{T}, θ, sys::Kuramoto) where {T <: Real}
    @unpack A = sys
    rows = rowvals(A)
    nzv = nonzeros(A)

    fill!(J, zero(T))
    for j in 1:size(A, 2), r in nzrange(A, j)
        i = rows[r]
        i != j || continue
        Aij = nzv[r]
        v = Aij * cos(θ[i] - θ[j])
        J[j, i] = v
        J[j, j] -= v
    end
    nothing
end
