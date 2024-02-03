export KuramotoNetwork

using Graphs
using SparseArrays

struct KuramotoNetwork{T, ωType} <: AutonomousODESys
    A::SparseMatrixCSC{T, Int}
    ω::ωType

    function KuramotoNetwork(A::AbstractMatrix{T}, ω::ωType = zero(T)) where {T, ωType}
        A = sparse(A)
        nr, nc = size(A)
        nr == nc || 
            throw(DimensionMismatch("Interaction matrix A must be square."))

        length(ω) == 1 || length(ω) == nr || 
            throw(DimensionMismatch("Dimension of ω must be compatible with A."))
        
        num_selfloops = sum(A[i, i] != 0 for i=1:nr)
        num_selfloops == 0 ||
            @warn "Supplied interaction matrix has self-loops; ignoring."
        @views fill!(A[diagind(A)], zero(eltype(A)))
        new{T, ωType}(sparse(A), copy(ω))
    end
end

# from a graph
KuramotoNetwork(g::AbstractGraph, args...) = 
    KuramotoNetwork(adjacency_matrix(g), args...)

function rhs(dθdt, θ, sys::KuramotoNetwork)
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

function jac(J::AbstractMatrix{T}, θ, sys::KuramotoNetwork) where {T}
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
