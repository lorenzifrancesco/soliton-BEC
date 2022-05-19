using Pkg
Pkg.activate(SolitonBEC)
using SolitonBEC
@time include("./test/test_transmission.jl")
