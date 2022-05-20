module SolitonBEC

export Numerics, Simulation, Apparatus, Potential, InitialState, run_dynamics, run_ground_state, ssfm_solve, ground_state_solve, get_coefficients, potential, mem_estimate, graphics_config, plot_dynamics, plot_ground_state, ssfm_propagate

using Printf
using FFTW
using QuadGK
using Elliptic

include("physics.jl")
include("solver.jl")
include("utils.jl")

end #module