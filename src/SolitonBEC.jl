module SolitonBEC

export Numerics, Simulation, Apparatus, Potential, InitialState, run_dynamics, run_ground_state, mem_estimate, graphics_config

using Printf
using PhysicalConstants
using FFTW
using Plots
pyplot()
using ProgressMeter
using QuadGK
using Elliptic

include("physics.jl")
include("solver.jl")
include("utils.jl")

end #module