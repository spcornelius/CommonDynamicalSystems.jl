module SPCDyn

using Distances
using FastBroadcast
using Graphs
using LinearAlgebra
using SparseArrays
using StaticArrays
using UnPack

include("./base.jl")
include("./kuramoto.jl")
include("./magnetic_pendulum.jl")

end 
