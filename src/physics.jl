struct Simulation
    T::Float64
    dt::Float64
    S::Float64
    ds::Float64
end


struct Apparatus
    m::Float64 # kg
    as::Float64 # m
    omega_perp::Float64 # rad/s
    N::Int64
    l_perp::Float64 # m
  end
  

  struct InitialState
    P0::Float64
    S0::Float64
    sech_flag::Float64
  end


 struct Coefficients
    α::ComplexF64
    β::Function
    γ::Function
 end


function get_coefficients(app::Apparatus, state::InitialState)
    hbar = PhysicalConstants.CODATA2018.h / (2 * pi) 
    α = im * hbar/(2*app.m)
    
    # modify for GPE, NPSE...
    β(s) = im * hbar / (8*app.m)
    σ(ψ) = app.l_perp^2 * sqrt(1+2*app.as(app.N - 1) * abs(ψ)^2)
    γ(ψ) = - im * hbar/ (2*app.m * σ(ψ)^2) - im * app.m * app.omega_perp^2 / (2*hbar) * σ(ψ)^2 - im * 2 * app.as * (app.N - 1)/(app.m * σ(ψ)^2) * abs(ψ)^2 
    return Coefficients(α, β, γ)
end


function run(sim::Simulation, app::Apparatus, state::InitialState)
    coeffs = get_coefficients(sim, app, state)
    solve(sim, coeffs)
end