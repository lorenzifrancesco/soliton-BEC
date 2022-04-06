module SolitonBEC

export Simulation, Apparatus, InitialState, run_simulation, mem_estimate, graphics_config

using Printf
using PhysicalConstants
using FFTW
using Plots
using ProgressMeter

include("physics.jl")
include("solver.jl")
include("utils.jl")

end #module