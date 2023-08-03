export test_jac

using SciMLBase: isinplace
using ForwardDiff

function print_diffs(J_provided, J)
    @info "\t\tProvided\t\tNumerical"
    c = 0
    for idx in eachindex(IndexCartesian(), J)
        v1 = J_provided[idx]
        v2 = J[idx]
        if !(v1 ≈ v2)
            @info "$(idx.I)\t\t$v1\t$v2"
            c += 1
        end
    end
    @info "Detected $(c) errors in provided Jacobian."
    nothing
end

function test_jac(f, jac, x₀::AbstractVector{T}, p, t;
                  verbose = false) where {T}
    n = length(x₀)
    jac_iip = isinplace(jac, 4)
    J_provided = if jac_iip
        J = zeros(T, n, n)
        jac(J, x₀, p, t)
        J
    else
        jac(x₀, args...)
    end

    J = zeros(n, n)
    f_iip = isinplace(f, 4)
    if f_iip
        F = similar(x₀)
        ForwardDiff.jacobian!(J, (dx, x) -> f(dx, x, p, t), F, x₀)
    else
        ForwardDiff.jacobian!(J, x -> f(x, p, t), x₀)
    end

    if J_provided ≈ J
        verbose && @info "No errors detected in provided Jacobian."
        return true
    else
        verbose && print_diffs(J_provided, J)
        return false
    end
end

test_jac(sys::AbstractODESys, x₀, t; kw...) =
    test_jac(rhs, jac, x₀, sys, t; kw...)

test_jac(sys::AutonomousODESys, x₀; kw...) =
    test_jac(rhs, jac, x₀, sys, 0.0; kw...)