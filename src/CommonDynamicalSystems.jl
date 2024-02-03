module CommonDynamicalSystems

using FastBroadcast
using UnPack

include("./base.jl")

include("./util/fixed_points.jl")
include("./util/jacobian.jl")

include("./continuous/kuramoto.jl")
include("./continuous/magnetic_pendulum.jl")

end 
