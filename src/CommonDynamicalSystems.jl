module CommonDynamicalSystems

using FastBroadcast
using UnPack

include("./base.jl")
include("./fixed_points.jl")
include("./kuramoto.jl")
include("./magnetic_pendulum.jl")
include("./util.jl")

end 
