struct Simulation
    T::Float64
    dt::Float64
    S::Float64
    ds::Float64
    equation::String #("GPE", "NPSE")
end


struct Apparatus
    m::Float64 # kg ω
    as::Float64 # m
    ω_perp::Float64 # rad/s
    N::Int64
    γ::Float64
    ω_z::Float64
  end
  

  struct InitialState
    width::Float64 # m
    sech_flag::Float64 # 1 or 0
  end


 struct Coefficients
    α::ComplexF64
    β::Function
    γ::Function
 end


function get_coefficients(sim::Simulation, app::Apparatus, state::InitialState)
    # implement with PhysicalConstants.CODATA2018.h
    hbar = 6.62607015e-34 / (2 * pi) 
    l_perp = sqrt(hbar / (app.m * app.ω_perp))
    l_z = sqrt(hbar / (app.m * app.ω_z))
    npse_gamma = app.γ * hbar * app.ω_perp / l_z^6 * 1e-4
    @assert(sim.equation in ["NPSE", "GPE"])
    if (sim.equation == "NPSE")
        α = -im * hbar/(2*app.m)
        β1(s) = im * hbar / (8*app.m)
        σ(ψ::ComplexF64) = l_perp * sqrt(sqrt(1+2*app.as*(app.N - 1) * abs.(ψ)^2))
        γ1(ψ::ComplexF64) = - im * hbar/ (2*app.m * σ(ψ)^2) - im * app.m * app.ω_perp^2 / (2*hbar) * σ(ψ)^2 - im * 2 * hbar * app.as * (app.N - 1)/(app.m * σ(ψ)^2) * abs.(ψ)^2 - npse_gamma * app.N^2 * abs.(ψ)^4
        β=β1
        γ=γ1
    elseif (sim.equation == "GPE")
        α = -im * hbar/(2*app.m) 
        β2(s) = - im * app.ω_perp - im * hbar * exp(-(s-30e-6)^2/(5e-6)^2) *9e12 / (8*app.m)
        #display( 2 * hbar * app.as * (app.N - 1) / (app.m * l_perp^2))
        γ2(ψ::ComplexF64) = - im * 2 * hbar * app.as * (app.N - 1) / (app.m * l_perp^2) * abs.(ψ)^2
        β=β2
        γ=γ2
    end

    return Coefficients(α, β, γ)
end


function run_simulation(sim::Simulation, app::Apparatus, state::InitialState)
    coeffs = get_coefficients(sim, app, state)
    describe_simulation(sim, app, state)
    ssfm_solve(sim, coeffs, state, app)
end