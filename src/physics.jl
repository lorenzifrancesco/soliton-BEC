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
    # implement with PhysicalConstants.CODATA2018.h
    hbar = 6.62607015e-34 / (2 * pi) 
       

    ## NPSE
    # α = im * hbar/(2*app.m)
    # β(s) = im * hbar / (8*app.m)
    # σ(ψ::ComplexF64) = app.l_perp * sqrt(1+2*app.as*(app.N - 1) * abs(ψ)^2) NPSE
    # γ(ψ::ComplexF64) = - im * hbar/ (2*app.m * σ(ψ)^2) - im * app.m * app.omega_perp^2 / (2*hbar) * σ(ψ)^2 - im * 2 * app.as * (app.N - 1)/(app.m * σ(ψ)^2) * abs(ψ)^2 NPSE
    
    ## GPE
    α = im * hbar/(2*app.m) - im * app.omega_perp
    # β(s) = im * hbar / (8*app.m)
    β(s) = 0
    γ(ψ::ComplexF64) = - im * hbar * app.as * (app.N - 1) / (app.m * app.l_perp) * abs.(ψ)^2
    return Coefficients(α, β, γ)
end


function run_simulation(sim::Simulation, app::Apparatus, state::InitialState)
    coeffs = get_coefficients(app, state)
    ssfm_solve(sim, coeffs, state)
end