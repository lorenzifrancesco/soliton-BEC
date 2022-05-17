struct Numerics
    T::Float64
    dt::Float64
    S::Float64
    ds::Float64
end

struct Simulation
    equation::String #("GPE", "NPSE")
    potential::String #("barrier", "ellipse")
end


struct Potential
    type::String
    width::Float64
    position::Float64
    energy::Float64
    ϵ::Float64
    a::Float64
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
    type::String
    width::Float64
    v0::Float64
    x0::Float64
end

struct Coefficients
    α::ComplexF64
    β::Function
    γ::Function
    initial::Function
end

function potential(sim::Simulation, app::Apparatus, pot::Potential, s::Float64)
    hbar = 6.62607015e-34 / (2 * pi)
    l_perp = sqrt(hbar / (app.m * app.ω_perp))
    l_z = sqrt(hbar / (app.m * app.ω_z))
    if pot.type == "barrier"
        value = pot.energy / (sqrt(2 * pi) * pot.width) * exp(-(s - pot.position)^2 / (2 * (pot.width))^2)
    elseif pot.type == "ellipse"
        target_error = 0.001 * pot.a
        lower_ϕ = -pi
        upper_ϕ = 3 * pi
        let error
            let midpoint
                error = target_error + 1
                while error > target_error
                    midpoint = (lower_ϕ + upper_ϕ) / 2
                    if s - pot.a * E(midpoint, pot.ϵ) > 0
                        lower_ϕ = midpoint
                    else
                        upper_ϕ = midpoint
                    end
                    error = abs(s - pot.a * E(midpoint, pot.ϵ))
                end
                value = -hbar^2 * (1 / pot.a * (sqrt(1 - pot.ϵ^2)) / (sin(midpoint)^2 + sqrt(1 - pot.ϵ^2) * cos(midpoint)^2)^(3 / 2))^2 / (8 * app.m)
            end
        end
    end
    return value
end


function β(sim::Simulation, app::Apparatus, pot::Potential, s::Float64)
    hbar = 6.62607015e-34 / (2 * pi)
    l_perp = sqrt(hbar / (app.m * app.ω_perp))
    l_z = sqrt(hbar / (app.m * app.ω_z))
    if (sim.equation == "NPSE")
        value = im * hbar / (8 * app.m)
    elseif (sim.equation == "GPE")
        value = -im * app.ω_perp
    end
    value += (-im / hbar) * potential(sim::Simulation, app::Apparatus, pot::Potential, s::Float64)
    return value
end


function γ(sim::Simulation, app::Apparatus, state::InitialState, ψ::ComplexF64)
    hbar = 6.62607015e-34 / (2 * pi)
    l_perp = sqrt(hbar / (app.m * app.ω_perp))
    l_z = sqrt(hbar / (app.m * app.ω_z))
    npse_gamma = app.γ * hbar * app.ω_perp / l_z^6 * 1e-4
    σ(ψ::ComplexF64) = l_perp * sqrt(sqrt(1 + 2 * app.as * (app.N - 1) * abs.(ψ)^2))
    if sim.equation == "NPSE"
        value = -im * hbar / (2 * app.m * σ(ψ)^2) - im * app.m * app.ω_perp^2 / (2 * hbar) * σ(ψ)^2 - im * 2 * hbar * app.as * (app.N - 1) / (app.m * σ(ψ)^2) * abs.(ψ)^2 - npse_gamma * app.N^2 * abs.(ψ)^4
    elseif sim.equation == "GPE"
        value = -im * 2 * hbar * app.as * (app.N - 1) / (app.m * l_perp^2) * abs.(ψ)^2
    end
    return value
end


function wave(sim::Simulation, app::Apparatus, state::InitialState, s::Float64)
    hbar = 6.62607015e-34 / (2 * pi)
    l_perp = sqrt(hbar / (app.m * app.ω_perp))
    interaction_g = abs(2*hbar^2 * app.as / app.m / l_perp^2)
    ggg = 2*hbar*app.ω_perp*abs(app.as)
    GSEnergy = -app.N^3/24*ggg^2 * app.m/hbar^2 
    Energy = GSEnergy
    k0 = app.m/hbar * state.v0 * 2

    if state.type == "gaussian" # Gaussian Pulse
        return sqrt(1 / (sqrt(2 * pi) * state.width)) * exp.(-(s-state.x0) .^ 2 / (4 * state.width^2)) * exp(im * (s-state.x0) * k0)
    elseif state.type == "sech"
        return sqrt(1 / (2 * state.width)) * 2 ./ (exp.(-((s-state.x0) / (state.width))) .+ exp.((s-state.x0) / (state.width))) * exp(im * (s-state.x0) * k0)
    end
end

function get_coefficients(sim::Simulation, app::Apparatus, pot::Potential, state::InitialState)
    # implement with PhysicalConstants.CODATA2018.h
    hbar = 6.62607015e-34 / (2 * pi)
    l_perp = sqrt(hbar / (app.m * app.ω_perp))
    l_z = sqrt(hbar / (app.m * app.ω_z))
    npse_gamma = app.γ * hbar * app.ω_perp / l_z^6 * 1e-4
    @assert(sim.equation in ["NPSE", "GPE"])
    @assert(pot.type in ["barrier", "ellipse", ""])

    α = -im * hbar / (2 * app.m)
    beta(s::Float64) = β(sim::Simulation, app::Apparatus, pot::Potential, s)
    gamma(ψ::ComplexF64) = γ(sim::Simulation, app::Apparatus, state::InitialState, ψ)
    initial_state(s::Float64) = wave(sim::Simulation, app::Apparatus, state::InitialState, s)
    return Coefficients(α, beta, gamma, initial_state)
end


function run_ground_state(num::Numerics, sim::Simulation, app::Apparatus, pot::Potential, state::InitialState)

    coeffs = get_coefficients(sim, app, pot, state)

    display("Running ground-state simulation")
    results = ground_state_solve(num, coeffs)
    plot_ground_state(results...)
end


function run_dynamics(num::Numerics, sim::Simulation, app::Apparatus, pot::Potential, state::InitialState)

    coeffs = get_coefficients(sim, app, pot, state)

    display("Running Dynamic simulation")
    results = ssfm_solve(num, coeffs)
    plot_dynamics(results...)
end