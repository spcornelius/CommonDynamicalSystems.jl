export AbstractODESys, rhs, jac

# supertype for a system of ODEs. Different ODE systems can be defined
# by creating subtypes, which should implement in-place versions of
# rhs and (optionally) jac, representing the differential equations and
# jacobian respectively. They should also implement "state_dims", giving the
# phase space dimension for a particular instance
abstract type AbstractODESys end;

# dummy functions to define the ODESys interface
function rhs end
function jac end

# don't attempt broadcast over an ODESys
# facilitates broadcasting with "." over methods that accept ODESys alongisde
# other arguments you DO want to boradcast over
Base.broadcastable(sys::AbstractODESys) = Ref(sys)

# generic out-of-place versions of rhs & jac, for convience
function rhs(x, sys::AbstractODESys, t = 0.0)
    dxdt = similar(x)
    rhs(dxdt, x, sys, t)
    return dxdt
end

function jac(x, sys::AbstractODESys, t = 0.0)
    n = length(x)
    J = similar(x, n, n)
    jac(J, x, sys, t)
    return J
end