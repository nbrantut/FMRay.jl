module FMM

using DataStructures
using Roots #this is for velocity stuff
import Statistics: middle

export march, traceray
export refine, refine!, coarsen_derivative, coarsen_derivative!
export precomputeT, locatelookup
export domainsize, getcartposition, getcartindex

include("main.jl")
include("velocity.jl")
include("ray.jl")
include("lookup.jl")
include("refine.jl")
include("common.jl")

end
