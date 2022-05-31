module SolitonBEC

export Numerics, Simulation, Apparatus, Potential, InitialState, run_dynamics, run_ground_state, ssfm_solve, ssfm_propagate, ground_state_solve, get_coefficients, potential, mem_estimate, graphics_config, plot_dynamics, plot_ground_state, pseudospectral_solve, propagation_design!

using Printf
using PhysicalConstants
using FFTW
using ProgressMeter
using QuadGK
using Elliptic
using DifferentialEquations

include("physics.jl")
include("solver.jl")
include("utils.jl")

end #module