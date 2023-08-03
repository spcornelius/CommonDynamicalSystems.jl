export find_fixed_point, is_stable

using ForwardDiff
using LinearAlgebra
using NLsolve
using SciMLBase

function find_fixed_point(rhs::Function, jac::Function, x₀::AbstractVector{T}, 
                          p = SciMLBase.NullParameters(); kw...) where {T}
    f! = (F, x) -> rhs(F, x, p, nothing)
    j! = (J, x) -> jac(J, x, p, nothing)
    n = length(x₀)
    J = zeros(T, n, n)
    result = nlsolve(f!, j!, x₀; kw...)
    result.f_converged || 
        error("Numerical optimization failed.")
    return result.zero
end

function find_fixed_point(rhs::Function, x₀::AbstractVector, 
                          p = SciMLBase.NullParameters(); kw...)
    f! = (F, x) -> rhs(F, x, p, nothing)
    result = nlsolve(f!, x₀; kw...)
    result.f_converged || 
        error("Numerical optimization failed.")
    return result.zero
end

function find_fixed_point(sys::sysType, x₀::x₀Type;
                          use_jac = true, kw...) where 
                          {sysType <: AutonomousODESys, x₀Type <: AbstractVector}
    T = eltype(x₀Type)
    if hasmethod(jac, (AbstractMatrix{T}, sysType, x₀Type)) && use_jac
        return find_fixed_point(rhs, jac, x₀, sys; kw...)
    else
        return find_fixed_point(rhs, x₀, sys; kw...)
    end
end

function is_stable(jac::Function, x::AbstractVector{T}, p = SciMLBase.NullParameters();
                   tol = zero(T)) where {T}
    n = length(x)
    J = zeros(T, n, n)
    jac(J, x, p, nothing)
    λs = eigvals!(J)
    maximum(real, λs) < -tol
end

function is_stable(sys::sysType, x::x₀Type; kw...) where
    {sysType <: AutonomousODESys, x₀Type <: AbstractVector}

    T = eltype(x₀Type)
    jac_ = if hasmethod(jac, (AbstractMatrix{T}, sysType, x₀Type)) && use_jac
        jac
    else
        F = similar(x)
        f! = (dx, x) -> rhs(dx, x, sys, nothing)
        (J, x, _, _) -> ForwardDiff.jacobian!(J, f!, F, x)
    end
    is_stable(jac_, x, sys; kw...)
end
