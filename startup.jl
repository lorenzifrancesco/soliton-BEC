using Pkg
pkg"add Distributed"
pkg"add Elliptic"
pkg"add FFTW"
pkg"add FileIO"
pkg"add PhysicalConstants"
pkg"add PlotlyBase"
pkg"add Plots"
pkg"add Printf"
pkg"add ProgressMeter"
pkg"add QuadGK"
pkg"add LinearAlgebra"
pkg"add Random"
pkg"add SparseArrays"
pkg"add Test"
pkg"precompile"
include("./src/SolitonBEC.jl")
include("./test/test_transmission.jl")
