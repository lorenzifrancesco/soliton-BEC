module SolitonBEC

export Numerics, Numerics_3D, Simulation, Apparatus, Potential, InitialState, InitialState_3D, run_dynamics, run_ground_state, ssfm_solve, ssfm_propagate, ground_state_solve, get_coefficients, potential, mem_estimate, graphics_config, plot_dynamics, plot_ground_state, pseudospectral_solve, propagation_design!, ssfm_propagate_3d, get_coefficients_3d, ssfm_solve_3d

using Printf
using PhysicalConstants
using FFTW
using ProgressMeter
using QuadGK
using Elliptic
using UnPack
using Plots
using CUDA

include("physics.jl")
include("solver.jl")
include("utils.jl")

end #module