export find_fixed_point, is_stable

function find_fixed_point(sys::sysType, x₀::AbstractVector{T}; jac = true, kw...) where 
    {sysType <: AutonomousODESys, T}
    f! = (F, x) -> rhs(F, x, sys)
    result = if hasmethod(jac, (AbstractMatrix, sysType, AbstractVector{T})) && jac
        n = length(x₀)
        J = zeros(T, n, n)
        j! = (J, x) -> jac(J, x, sys)
        nlsolve(f!, j!, x₀; kw...)
    else
        nlsolve(f!, x₀; kw...)
    end
    result.f_converged || 
        error("Numerical optimization failed.")
    return result.zero
end

function is_stable(sys::AutonomousODESys, x::AbstractVector{T}; 
                   tol = zero(T)) where {T}
    n = length(x)
    J = zeros(T, n, n)
    jac(J, x, sys)
    λs = eigvals!(J)
    return maximum(real, λs) < -tol
end
