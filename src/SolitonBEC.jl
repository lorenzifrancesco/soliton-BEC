module SolitonBEC

export Simulation, Apparatus, InitialState, run, mem_estimate

using Printf
using PhysicalConstants

include("physics.jl")
include("solver.jl")
include("utils.jl")

end #module