module SolitonBEC

export Simulation, Apparatus, InitialState, run_simulation, mem_estimate, graphics_config

using Printf
using PhysicalConstants
using FFTW
using Plots; pyplot()
using ProgressMeter
using QuadGK

include("physics.jl")
include("solver.jl")
include("utils.jl")

end #module