module FMM

using DataStructures
using Roots #this is for velocity stuff

export march

include("main.jl")
include("velocity.jl")

end
