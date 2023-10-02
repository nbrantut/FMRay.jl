module FMM

using DataStructures
using Roots #this is for velocity stuff
import Statistics: middle

export Grid
export march, traceray
export refine, refine!, coarsen_derivative, coarsen_derivative!
export precomputeT, locatelookup
export domainsize, getcartposition, getcartindex
export vgroup, vphase, tan_group_angle, phase_angle

include("main.jl")
include("velocity.jl")
include("ray.jl")
include("lookup.jl")
include("refine.jl")
include("common.jl")

# precompile march
march(CartesianIndex(1,1,1), Grid(0.1, fill(1.0, (10,10,1))))

end
