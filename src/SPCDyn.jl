module SPCDyn

using Distances
using FastBroadcast
using Graphs
using LinearAlgebra
using NLsolve
using SparseArrays
using StaticArrays
using UnPack

include("./base.jl")
include("./fixed_points.jl")
include("./kuramoto.jl")
include("./magnetic_pendulum.jl")
include("./util.jl")

end 
